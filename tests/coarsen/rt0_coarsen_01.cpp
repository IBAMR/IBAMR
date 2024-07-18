// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Headers for major SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/CartSideDoubleRT0Coarsen.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/muParserCartGridFunction.h>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

// Verify that we can correctly coarsen a piecewise linear solution with the
// RT0 coarsening class. Since the RT0 element is a vector-valued element that
// is, on Cartesian grids,
//
//     RT0_K = (P^1(x) * P^0(y), P^0(x) * P^1(y))
//
// we expect coarsening a vector field which is in that space to have zero error.

int
main(int argc, char* argv[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);
    {
        auto app_initializer = make_samrai_shared<AppInitializer>(argc, argv, "rt0.log");
        SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.
        auto grid_geometry = make_samrai_shared<CartesianGridGeometryNd>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        auto patch_hierarchy = make_samrai_shared<PatchHierarchyNd>("PatchHierarchy", grid_geometry);
        auto error_detector = make_samrai_shared<StandardTagAndInitializeNd>(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        auto box_generator = make_samrai_shared<BergerRigoutsosNd>();
        auto load_balancer =
            make_samrai_shared<LoadBalancerNd>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        auto gridding_algorithm =
            make_samrai_shared<GriddingAlgorithmNd>("GriddingAlgorithm",
                                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                                    error_detector,
                                                    box_generator,
                                                    load_balancer);

        SAMRAIPointer<VisItDataWriterNd> visit_data_writer = app_initializer->getVisItDataWriter();

        // Create variables and register them with the variable database.
        VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();
        SAMRAIPointer<VariableContext> ctx = var_db->getContext("context");
        SAMRAIPointer<SideVariableNd<double> > u_sc_var = make_samrai_shared<SideVariableNd<double> >("u_sc");
        const int u_sc_idx = var_db->registerVariableAndContext(u_sc_var, ctx, IntVectorNd(4));
        SAMRAIPointer<SideVariableNd<double> > exact_sc_var = make_samrai_shared<SideVariableNd<double> >("exact_sc");
        const int exact_sc_idx = var_db->registerVariableAndContext(exact_sc_var, ctx, IntVectorNd(4));
        // TODO u_cc_var is only for plotting (and testing): remove later
        // #define DO_PLOT
#ifdef DO_PLOT
        SAMRAIPointer<CellVariableNd<double> > u_cc_var = make_samrai_shared<CellVariableNd<double> >("u_cc", NDIM);
        const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx);
        SAMRAIPointer<CellVariableNd<double> > exact_cc_var =
            make_samrai_shared<CellVariableNd<double> >("exact_cc", NDIM);
        const int exact_cc_idx = var_db->registerVariableAndContext(exact_cc_var, ctx);
#endif

        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        const int tag_buffer = std::numeric_limits<int>::max();
        int level_number = 0;
        while ((gridding_algorithm->levelCanBeRefined(level_number)))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            ++level_number;
        }

        const int finest_level = patch_hierarchy->getFinestLevelNumber();
        for (int ln = 0; ln <= finest_level; ++ln)
        {
            SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_sc_idx, 0.0);
            level->allocatePatchData(exact_sc_idx, 0.0);
#ifdef DO_PLOT
            level->allocatePatchData(u_cc_idx, 0.0);
            level->allocatePatchData(exact_cc_idx, 0.0);
#endif
        }

        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);

        {
            muParserCartGridFunction u_fcn("u", app_initializer->getComponentDatabase("u"), grid_geometry);
            u_fcn.setDataOnPatchHierarchy(u_sc_idx, u_sc_var, patch_hierarchy, 0.0);
            u_fcn.setDataOnPatchHierarchy(exact_sc_idx, exact_sc_var, patch_hierarchy, 0.0);
        }
        SAMRAIPointer<VisItDataWriterNd> visit_writer = app_initializer->getVisItDataWriter();

        // Fill in ghost cells. Since this is a no coarse fine interface, the actual refine operator doesn't matter. We
        // will just be copying values.
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_cell_comps = { ITC(
            u_sc_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE", "NONE", false, nullptr, nullptr, "NONE") };
        HierarchyGhostCellInterpolation hier_ghost_cell;
        hier_ghost_cell.initializeOperatorState(ghost_cell_comps, patch_hierarchy);
        hier_ghost_cell.fillData(0.0);
        // Now we need to fill physical ghost cells to be NaNs to ensure that the coarsen operator only touches physical
        // cells.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevelNd::Iterator p(level); p; p++)
            {
                SAMRAIPointer<PatchNd> patch = level->getPatch(p());
                SAMRAIPointer<SideDataNd<double> > u_data = patch->getPatchData(u_sc_idx);
                SAMRAIPointer<PatchGeometryNd> pgeom = patch->getPatchGeometry();

                tbox::Array<BoundaryBoxNd> bdry_boxes =
                    PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
                for (int i = 0; i < bdry_boxes.size(); ++i)
                {
                    const BoundaryBoxNd& bdry_box = bdry_boxes[i];
                    const int bdry_axis = bdry_box.getLocationIndex() / 2;
                    const int upper_lower = bdry_box.getLocationIndex() % 2;
                    if (!pgeom->getTouchesRegularBoundary(bdry_axis, upper_lower)) continue;
                    const BoxNd& fill_box =
                        pgeom->getBoundaryFillBox(bdry_box, patch->getBox(), u_data->getGhostCellWidth());
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for (SideIteratorNd si(fill_box, axis); si; si++)
                        {
                            const SideIndexNd& idx = si();
                            if (!patch->getBox().contains(idx.toCell(0)) && !patch->getBox().contains(idx.toCell(1)))
                                (*u_data)(idx) = std::numeric_limits<double>::quiet_NaN();
                        }
                    }
                }
            }
        }

        // The rest is just book-keeping, this is the actual test:
        {
            IntVectorNd ratio;
            for (int d = 0; d < NDIM; ++d) ratio(d) = 4;
            auto coarsen_alg = make_samrai_shared<CoarsenAlgorithmNd>();
            SAMRAIPointer<CoarsenOperatorNd> coarsen_op = make_samrai_shared<IBTK::CartSideDoubleRT0Coarsen>(ratio);
            coarsen_alg->registerCoarsen(u_sc_idx, u_sc_idx, coarsen_op);
            SAMRAIPointer<CoarsenScheduleNd> coarsen_sched =
                coarsen_alg->createSchedule(patch_hierarchy->getPatchLevel(0), patch_hierarchy->getPatchLevel(1));
            coarsen_sched->coarsenData();
        }

        HierarchySideDataOpsRealNd<double> hier_data_ops(patch_hierarchy);
        hier_data_ops.subtract(exact_sc_idx, exact_sc_idx, u_sc_idx, false);

        pout << "max norm error = " << hier_data_ops.maxNorm(exact_sc_idx) << '\n';

#ifdef DO_PLOT
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity(u_cc_var->getName() + std::to_string(d), "SCALAR", u_cc_idx, d);
            visit_data_writer->registerPlotQuantity(
                exact_cc_var->getName() + std::to_string(d), "SCALAR", exact_cc_idx, d);
        }
#endif

        // Check the values on the coarsest level.
        SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(0);
        for (PatchLevelNd::Iterator p(level); p; p++)
        {
            pout << "patch number " << p() << '\n';
            SAMRAIPointer<PatchNd> patch = level->getPatch(p());
            SAMRAIPointer<SideDataNd<double> > e_data = patch->getPatchData(exact_sc_idx);
            const BoxNd patch_box = patch->getBox();

            // same as SideData::print, but elides zero values. We don't
            // print any information about the patch when no values are
            // above the cutoff.
            for (int axis = 0; axis < NDIM; ++axis)
            {
                pout << "Array side normal = " << axis << std::endl;
                for (int d = 0; d < e_data->getDepth(); ++d)
                {
                    pout << "Array depth = " << d << std::endl;
                    const ArrayDataNd<double>& data = e_data->getArrayData(axis);
                    for (SideIteratorNd i(patch_box, axis); i; i++)
                    {
                        const double value = data(i(), d);
                        if (std::abs(value) > 1e-12)
                        {
                            pout << "array" << i() << " = " << value << '\n';
                        }
                    }
                }
            }
        }

#ifdef DO_PLOT
        visit_writer->writePlotData(patch_hierarchy, 0, 0.0);
#endif
    }

    // At this point all SAMRAI, PETSc, and IBAMR objects have been cleaned
    // up, so we shut things down in the opposite order of initialization:
    SAMRAIManager::shutdown();
    PetscFinalize();
}
