// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2022 by the IBAMR developers
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
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc objects
#include <petscsys.h>

// Headers for major SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/SCLaplaceOperator.h>
#include <ibtk/muParserCartGridFunction.h>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

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
        Pointer<CartesianGridGeometryNd> grid_geometry = new CartesianGridGeometryNd(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchyNd> patch_hierarchy = new PatchHierarchyNd("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitializeNd> error_detector = new StandardTagAndInitializeNd(
            "StandardTagAndInitialize", NULL, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsosNd> box_generator = new BergerRigoutsosNd();
        Pointer<LoadBalancerNd> load_balancer =
            new LoadBalancerNd("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithmNd> gridding_algorithm =
            new GriddingAlgorithmNd("GriddingAlgorithm",
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    error_detector,
                                    box_generator,
                                    load_balancer);

        // Create variables and register them with the variable database.
        VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");

        Pointer<SideVariableNd<double> > u_sc_var = new SideVariableNd<double>("u_sc");
        Pointer<SideVariableNd<double> > f_sc_var = new SideVariableNd<double>("f_sc");
        Pointer<SideVariableNd<double> > e_sc_var = new SideVariableNd<double>("e_sc");

        const int u_sc_idx = var_db->registerVariableAndContext(u_sc_var, ctx, IntVectorNd(1));
        const int f_sc_idx = var_db->registerVariableAndContext(f_sc_var, ctx, IntVectorNd(1));
        const int e_sc_idx = var_db->registerVariableAndContext(e_sc_var, ctx, IntVectorNd(1));

        Pointer<CellVariableNd<double> > u_cc_var = new CellVariableNd<double>("u_cc", NDIM);
        Pointer<CellVariableNd<double> > f_cc_var = new CellVariableNd<double>("f_cc", NDIM);
        Pointer<CellVariableNd<double> > e_cc_var = new CellVariableNd<double>("e_cc", NDIM);

        const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx, IntVectorNd(0));
        const int f_cc_idx = var_db->registerVariableAndContext(f_cc_var, ctx, IntVectorNd(0));
        const int e_cc_idx = var_db->registerVariableAndContext(e_cc_var, ctx, IntVectorNd(0));

        // Register variables for plotting.
        Pointer<VisItDataWriterNd> visit_data_writer = app_initializer->getVisItDataWriter();
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

        // Allocate data on each level of the patch hierarchy.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
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

        SAMRAIVectorRealNd<double> u_vec("u", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorRealNd<double> f_vec("f", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorRealNd<double> e_vec("e", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

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
        std::vector<RobinBcCoefStrategyNd*> bc_coefs(NDIM, static_cast<RobinBcCoefStrategyNd*>(NULL));
        SCLaplaceOperator laplace_op("laplace op");
        laplace_op.setPoissonSpecifications(poisson_spec);
        laplace_op.setPhysicalBcCoefs(bc_coefs);
        laplace_op.initializeOperatorState(u_vec, f_vec);
        laplace_op.apply(u_vec, f_vec);

        // Compute error and print error norms.
        e_vec.subtract(Pointer<SAMRAIVectorRealNd<double> >(&e_vec, false),
                       Pointer<SAMRAIVectorRealNd<double> >(&f_vec, false));
        pout << "|e|_oo = " << e_vec.maxNorm() << "\n";
        pout << "|e|_2  = " << e_vec.L2Norm() << "\n";
        pout << "|e|_1  = " << e_vec.L1Norm() << "\n";

        // Interpolate the side-centered data to cell centers for output.
        static const bool synch_cf_interface = true;
        hier_math_ops.interp(u_cc_idx, u_cc_var, u_sc_idx, u_sc_var, NULL, 0.0, synch_cf_interface);
        hier_math_ops.interp(f_cc_idx, f_cc_var, f_sc_idx, f_sc_var, NULL, 0.0, synch_cf_interface);
        hier_math_ops.interp(e_cc_idx, e_cc_var, e_sc_idx, e_sc_var, NULL, 0.0, synch_cf_interface);

        // Set invalid values on coarse levels (i.e., coarse-grid values that
        // are covered by finer grid patches) to equal zero.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber() - 1; ++ln)
        {
            Pointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            BoxArrayNd refined_region_boxes;
            Pointer<PatchLevelNd> next_finer_level = patch_hierarchy->getPatchLevel(ln + 1);
            refined_region_boxes = next_finer_level->getBoxes();
            refined_region_boxes.coarsen(next_finer_level->getRatioToCoarserLevel());
            for (PatchLevelNd::Iterator p(level); p; p++)
            {
                Pointer<PatchNd> patch = level->getPatch(p());
                const BoxNd& patch_box = patch->getBox();
                Pointer<CellDataNd<double> > e_cc_data = patch->getPatchData(e_cc_idx);
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                {
                    const BoxNd refined_box = refined_region_boxes[i];
                    const BoxNd intersection = BoxNd::grow(patch_box, 1) * refined_box;
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
} // main
