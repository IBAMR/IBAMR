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
    PetscInitialize(&argc, &argv, NULL, NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);
    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "rt0.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", NULL, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();

        // Create variables and register them with the variable database.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");
        Pointer<SideVariable<NDIM, double> > u_sc_var = new SideVariable<NDIM, double>("u_sc");
        const int u_sc_idx = var_db->registerVariableAndContext(u_sc_var, ctx, IntVector<NDIM>(4));
        Pointer<SideVariable<NDIM, double> > exact_sc_var = new SideVariable<NDIM, double>("exact_sc");
        const int exact_sc_idx = var_db->registerVariableAndContext(exact_sc_var, ctx, IntVector<NDIM>(4));
        // TODO u_cc_var is only for plotting (and testing): remove later
        // #define DO_PLOT
#ifdef DO_PLOT
        Pointer<CellVariable<NDIM, double> > u_cc_var = new CellVariable<NDIM, double>("u_cc", NDIM);
        const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx);
        Pointer<CellVariable<NDIM, double> > exact_cc_var = new CellVariable<NDIM, double>("exact_cc", NDIM);
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
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
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
        Pointer<VisItDataWriter<NDIM> > visit_writer = app_initializer->getVisItDataWriter();

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
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_sc_idx);
                Pointer<PatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

                tbox::Array<BoundaryBox<NDIM> > bdry_boxes =
                    PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
                for (int i = 0; i < bdry_boxes.size(); ++i)
                {
                    const BoundaryBox<NDIM>& bdry_box = bdry_boxes[i];
                    const int bdry_axis = bdry_box.getLocationIndex() / 2;
                    const int upper_lower = bdry_box.getLocationIndex() % 2;
                    if (!pgeom->getTouchesRegularBoundary(bdry_axis, upper_lower)) continue;
                    const Box<NDIM>& fill_box =
                        pgeom->getBoundaryFillBox(bdry_box, patch->getBox(), u_data->getGhostCellWidth());
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for (SideIterator<NDIM> si(fill_box, axis); si; si++)
                        {
                            const SideIndex<NDIM>& idx = si();
                            if (!patch->getBox().contains(idx.toCell(0)) && !patch->getBox().contains(idx.toCell(1)))
                                (*u_data)(idx) = std::numeric_limits<double>::quiet_NaN();
                        }
                    }
                }
            }
        }

        // The rest is just book-keeping, this is the actual test:
        {
            IntVector<NDIM> ratio;
            for (int d = 0; d < NDIM; ++d) ratio(d) = 4;
            Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg = new CoarsenAlgorithm<NDIM>();
            Pointer<CoarsenOperator<NDIM> > coarsen_op = new IBTK::CartSideDoubleRT0Coarsen(ratio);
            coarsen_alg->registerCoarsen(u_sc_idx, u_sc_idx, coarsen_op);
            Pointer<CoarsenSchedule<NDIM> > coarsen_sched =
                coarsen_alg->createSchedule(patch_hierarchy->getPatchLevel(0), patch_hierarchy->getPatchLevel(1));
            coarsen_sched->coarsenData();
        }

        HierarchySideDataOpsReal<NDIM, double> hier_data_ops(patch_hierarchy);
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
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(0);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            pout << "patch number " << p() << '\n';
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM, double> > e_data = patch->getPatchData(exact_sc_idx);
            const Box<NDIM> patch_box = patch->getBox();

            // same as SideData::print, but elides zero values. We don't
            // print any information about the patch when no values are
            // above the cutoff.
            for (int axis = 0; axis < NDIM; ++axis)
            {
                pout << "Array side normal = " << axis << std::endl;
                for (int d = 0; d < e_data->getDepth(); ++d)
                {
                    pout << "Array depth = " << d << std::endl;
                    const ArrayData<NDIM, double>& data = e_data->getArrayData(axis);
                    for (SideIterator<NDIM> i(patch_box, axis); i; i++)
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
