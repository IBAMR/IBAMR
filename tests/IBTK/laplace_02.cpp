// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Config files

#include <SAMRAI_config.h>

// Headers for basic PETSc objects
#include <petscsys.h>

// Headers for major SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <FaceData.h>
#include <GriddingAlgorithm.h>
#include <HierarchyDataOpsManager.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/SCLaplaceOperator.h>
#include <ibtk/muParserCartGridFunction.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <vector>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

// A test program to check that the side-centered Laplace operator
// discretization yields the expected order of accuracy.
//
// If the input database sets synchronization_test = TRUE, the program instead
// checks that HierarchyMathOps::synchronizeCoarseFineBoundary() gives the same
// result as synchronizing each side-centered field alone.

namespace
{
// Synchronize two side-centered fields together with
// synchronizeCoarseFineBoundary(), and separately one at a time as the source of
// div() with synchronization; print norms of the results and of their
// differences.
void
run_synchronization_test(Pointer<AppInitializer> app_initializer,
                         Pointer<CartesianGridGeometry<NDIM>> grid_geometry,
                         Pointer<PatchHierarchy<NDIM>> patch_hierarchy)
{
    constexpr int NUM_FIELDS = 2;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("synchronization_context");

    std::vector<Pointer<SideVariable<NDIM, double>>> together_vars, alone_vars;
    std::vector<int> together_idxs, alone_idxs;
    for (int k = 0; k < NUM_FIELDS; ++k)
    {
        together_vars.push_back(new SideVariable<NDIM, double>("together_" + std::to_string(k)));
        together_idxs.push_back(var_db->registerVariableAndContext(together_vars.back(), ctx, IntVector<NDIM>(1)));
        alone_vars.push_back(new SideVariable<NDIM, double>("alone_" + std::to_string(k)));
        alone_idxs.push_back(var_db->registerVariableAndContext(alone_vars.back(), ctx, IntVector<NDIM>(1)));
    }
    Pointer<CellVariable<NDIM, double>> div_var = new CellVariable<NDIM, double>("div");
    const int div_idx = var_db->registerVariableAndContext(div_var, ctx, IntVector<NDIM>(0));

    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (int k = 0; k < NUM_FIELDS; ++k)
        {
            level->allocatePatchData(together_idxs[k], 0.0);
            level->allocatePatchData(alone_idxs[k], 0.0);
        }
        level->allocatePatchData(div_idx, 0.0);
    }

    // Fill the fields with different smooth data on every level.
    Pointer<HierarchyDataOpsReal<NDIM, double>> hier_sc_data_ops =
        HierarchyDataOpsManager<NDIM>::getManager()->getOperationsDouble(together_vars[0], patch_hierarchy, true);
    for (int k = 0; k < NUM_FIELDS; ++k)
    {
        muParserCartGridFunction fcn("field_" + std::to_string(k),
                                     app_initializer->getComponentDatabase("field_" + std::to_string(k)),
                                     grid_geometry);
        fcn.setDataOnPatchHierarchy(together_idxs[k], together_vars[k], patch_hierarchy, 0.0);
        hier_sc_data_ops->copyData(alone_idxs[k], together_idxs[k]);
    }

    HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
    hier_math_ops.synchronizeCoarseFineBoundary(together_idxs);
    for (int k = 0; k < NUM_FIELDS; ++k)
    {
        hier_math_ops.div(div_idx, div_var, 1.0, alone_idxs[k], alone_vars[k], nullptr, 0.0, true);
    }

    // The norms are not weighted by control volume, so that they include the
    // coarse-grid data that is covered by finer levels.
    plog << std::setprecision(12);
    plog << "number of levels = " << patch_hierarchy->getNumberOfLevels() << "\n";
    for (int k = 0; k < NUM_FIELDS; ++k)
    {
        plog << "field " << k
             << " synchronized together: max norm = " << hier_sc_data_ops->maxNorm(together_idxs[k], -1) << "\n";
        plog << "field " << k << " synchronized together: L2 norm = " << hier_sc_data_ops->L2Norm(together_idxs[k], -1)
             << "\n";
        hier_sc_data_ops->subtract(alone_idxs[k], together_idxs[k], alone_idxs[k]);
        plog << "field " << k << " max norm of difference from synchronized alone = "
             << std::abs(hier_sc_data_ops->maxNorm(alone_idxs[k], -1)) << "\n";
    }
}
} // namespace

/*******************************************************************************
 * For each run, the input filename must be given on the command line.  In all *
 * cases, the command line is:                                                 *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "sc_laplace.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Initialize the AMR patch hierarchy.
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        int tag_buffer = 1;
        int level_number = 0;
        bool done = false;
        while (!done && (gridding_algorithm->levelCanBeRefined(level_number)))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            done = !patch_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }

        if (input_db->getBoolWithDefault("synchronization_test", false))
        {
            run_synchronization_test(app_initializer, grid_geometry, patch_hierarchy);
            return EXIT_SUCCESS;
        }

        // Create variables and register them with the variable database.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");

        Pointer<SideVariable<NDIM, double>> u_sc_var = new SideVariable<NDIM, double>("u_sc");
        Pointer<SideVariable<NDIM, double>> f_sc_var = new SideVariable<NDIM, double>("f_sc");
        Pointer<SideVariable<NDIM, double>> e_sc_var = new SideVariable<NDIM, double>("e_sc");

        const int u_sc_idx = var_db->registerVariableAndContext(u_sc_var, ctx, IntVector<NDIM>(1));
        const int f_sc_idx = var_db->registerVariableAndContext(f_sc_var, ctx, IntVector<NDIM>(1));
        const int e_sc_idx = var_db->registerVariableAndContext(e_sc_var, ctx, IntVector<NDIM>(1));

        Pointer<CellVariable<NDIM, double>> u_cc_var = new CellVariable<NDIM, double>("u_cc", NDIM);
        Pointer<CellVariable<NDIM, double>> f_cc_var = new CellVariable<NDIM, double>("f_cc", NDIM);
        Pointer<CellVariable<NDIM, double>> e_cc_var = new CellVariable<NDIM, double>("e_cc", NDIM);

        const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx, IntVector<NDIM>(0));
        const int f_cc_idx = var_db->registerVariableAndContext(f_cc_var, ctx, IntVector<NDIM>(0));
        const int e_cc_idx = var_db->registerVariableAndContext(e_cc_var, ctx, IntVector<NDIM>(0));

        // Register variables for plotting.
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);

        visit_data_writer->registerPlotQuantity(u_cc_var->getName(), "VECTOR", u_cc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity(u_cc_var->getName() + std::to_string(d), "SCALAR", u_cc_idx, d);
        }

        visit_data_writer->registerPlotQuantity(f_cc_var->getName(), "VECTOR", f_cc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity(f_cc_var->getName() + std::to_string(d), "SCALAR", f_cc_idx, d);
        }

        visit_data_writer->registerPlotQuantity(e_cc_var->getName(), "VECTOR", e_cc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity(e_cc_var->getName() + std::to_string(d), "SCALAR", e_cc_idx, d);
        }

        // Allocate data on each level of the patch hierarchy.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_sc_idx, 0.0);
            level->allocatePatchData(f_sc_idx, 0.0);
            level->allocatePatchData(e_sc_idx, 0.0);
            level->allocatePatchData(u_cc_idx, 0.0);
            level->allocatePatchData(f_cc_idx, 0.0);
            level->allocatePatchData(e_cc_idx, 0.0);
        }

        // Setup vector objects.
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int h_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        SAMRAIVectorReal<NDIM, double> u_vec("u", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> f_vec("f", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> e_vec("e", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        u_vec.addComponent(u_sc_var, u_sc_idx, h_sc_idx);
        f_vec.addComponent(f_sc_var, f_sc_idx, h_sc_idx);
        e_vec.addComponent(e_sc_var, e_sc_idx, h_sc_idx);

        u_vec.setToScalar(0.0);
        f_vec.setToScalar(0.0);
        e_vec.setToScalar(0.0);

        // Setup exact solutions.
        muParserCartGridFunction u_fcn("u", app_initializer->getComponentDatabase("u"), grid_geometry);
        muParserCartGridFunction f_fcn("f", app_initializer->getComponentDatabase("f"), grid_geometry);

        u_fcn.setDataOnPatchHierarchy(u_sc_idx, u_sc_var, patch_hierarchy, 0.0);
        f_fcn.setDataOnPatchHierarchy(e_sc_idx, e_sc_var, patch_hierarchy, 0.0);

        // Compute -L*u = f.
        PoissonSpecifications poisson_spec("poisson_spec");
        poisson_spec.setCConstant(0.0);
        poisson_spec.setDConstant(-1.0);
        std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs(NDIM, nullptr);
        SCLaplaceOperator laplace_op("laplace op");
        laplace_op.setPoissonSpecifications(poisson_spec);
        laplace_op.setPhysicalBcCoefs(bc_coefs);
        laplace_op.initializeOperatorState(u_vec, f_vec);
        laplace_op.apply(u_vec, f_vec);

        // Compute error and print error norms.
        e_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double>>(&e_vec, false),
                       Pointer<SAMRAIVectorReal<NDIM, double>>(&f_vec, false));
        const double max_norm = e_vec.maxNorm();
        const double l2_norm = e_vec.L2Norm();
        const double l1_norm = e_vec.L1Norm();

        if (IBTK_MPI::getRank() == 0)
        {
            std::ofstream out("output");
            out << "|e|_oo = " << max_norm << "\n";
            out << "|e|_2  = " << l2_norm << "\n";
            out << "|e|_1  = " << l1_norm << "\n";
        }

        // Optionally check that the face weights sum to the volume of the
        // physical domain in each coordinate direction.
        if (input_db->getBoolWithDefault("test_face_weights", false))
        {
            const int h_fc_idx = hier_math_ops.getFaceWeightPatchDescriptorIndex();
            std::array<double, NDIM> h_fc_sum{};
            for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM>> patch = level->getPatch(p());
                    Pointer<FaceData<NDIM, double>> h_fc_data = patch->getPatchData(h_fc_idx);
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        for (FaceIterator<NDIM> fi(patch->getBox(), axis); fi; fi++)
                        {
                            h_fc_sum[axis] += (*h_fc_data)(fi());
                        }
                    }
                }
            }
            IBTK_MPI::sumReduction(h_fc_sum.data(), NDIM);

            if (IBTK_MPI::getRank() == 0)
            {
                std::ofstream out("output", std::ios_base::app);
                out << "volume = " << hier_math_ops.getVolumeOfPhysicalDomain() << "\n";
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    out << "sum(h_fc[" << axis << "]) = " << h_fc_sum[axis] << "\n";
                }
            }
        }

        // Interpolate the side-centered data to cell centers for output.
        static const bool synch_cf_interface = true;
        hier_math_ops.interp(u_cc_idx, u_cc_var, u_sc_idx, u_sc_var, nullptr, 0.0, synch_cf_interface);
        hier_math_ops.interp(f_cc_idx, f_cc_var, f_sc_idx, f_sc_var, nullptr, 0.0, synch_cf_interface);
        hier_math_ops.interp(e_cc_idx, e_cc_var, e_sc_idx, e_sc_var, nullptr, 0.0, synch_cf_interface);

        // Set invalid values on coarse levels (i.e., coarse-grid values that
        // are covered by finer grid patches) to equal zero.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber() - 1; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            BoxArray<NDIM> refined_region_boxes;
            Pointer<PatchLevel<NDIM>> next_finer_level = patch_hierarchy->getPatchLevel(ln + 1);
            refined_region_boxes = next_finer_level->getBoxes();
            refined_region_boxes.coarsen(next_finer_level->getRatioToCoarserLevel());
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<CellData<NDIM, double>> e_cc_data = patch->getPatchData(e_cc_idx);
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                {
                    const Box<NDIM> refined_box = refined_region_boxes[i];
                    const Box<NDIM> intersection = Box<NDIM>::grow(patch_box, 1) * refined_box;
                    if (!intersection.empty())
                    {
                        e_cc_data->fillAll(0.0, intersection);
                    }
                }
            }
        }

        // Output data for plotting.
        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);

    } // cleanup dynamically allocated objects prior to shutdown
} // run_example
