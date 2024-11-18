// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2024 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/AppInitializer.h>
#include <ibtk/CartCellDoubleQuadraticCFInterpolation.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/muParserCartGridFunction.h>

#include <petscsys.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CoarsenAlgorithm.h>
#include <CoarsenSchedule.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <StandardTagAndInitialize.h>

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
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "ghost_cells.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);
        // Create cell-centered and node-centered quantities, and initialize them with a function read from the input
        // file
        Pointer<Variable<NDIM> > Q_var;
        if (input_db->getString("CENTERING").compare("CELL") == 0)
            Q_var = new CellVariable<NDIM, double>("Q");
        else if (input_db->getString("CENTERING").compare("SIDE") == 0)
            Q_var = new SideVariable<NDIM, double>("Q");

        auto Q_fcn = [](const VectorNd& x) -> double
        {
            double val = 1.0;
            for (int d = 0; d < NDIM; ++d) val += x[d] + x[d] * x[d];
            return val;
        };

        auto var_db = VariableDatabase<NDIM>::getDatabase();
        // Total amount patch indices. Note we need a single ghost cell width on the cell centered quantity.
        const int Q_idx = var_db->registerVariableAndContext(Q_var, var_db->getContext("CTX"), IntVector<NDIM>(1));
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

        // Allocate patch data
        int coarsest_ln = 0, finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(Q_idx);
        }

        // Fill in initial data
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                const double* const xlow = pgeom->getXLower();
                const hier::Index<NDIM>& idx_low = patch->getBox().lower();
                Pointer<CellData<NDIM, double> > Q_cc_data = patch->getPatchData(Q_idx);
                Pointer<SideData<NDIM, double> > Q_sc_data = patch->getPatchData(Q_idx);
                if (Q_cc_data)
                {
                    for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
                    {
                        const CellIndex<NDIM>& idx = ci();
                        VectorNd x;
                        for (int d = 0; d < NDIM; ++d)
                            x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
                        (*Q_cc_data)(idx) = Q_fcn(x);
                    }
                }
                else if (Q_sc_data)
                {
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for (SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
                        {
                            const SideIndex<NDIM>& idx = si();
                            VectorNd x;
                            for (int d = 0; d < NDIM; ++d)
                                x[d] = xlow[d] +
                                       dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + (d == axis ? 0.0 : 0.5));
                            (*Q_sc_data)(idx) = Q_fcn(x);
                        }
                    }
                }
            }
        }

        // Test out the filling of coarse fine interfaces. We only use the coarse fine interface operator, which fills
        // one layer of ghost cells using quadratic interpolation.
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> ghost_cell_comps(1);
        ghost_cell_comps[0] =
            InterpolationTransactionComponent(Q_idx, "NONE", true, "NONE", "NONE", false, nullptr /*bdry_conds*/);
        HierarchyGhostCellInterpolation hier_ghost_cell;
        hier_ghost_cell.initializeOperatorState(ghost_cell_comps, patch_hierarchy);
        hier_ghost_cell.fillData(0.0);

        {
            int ln = patch_hierarchy->getFinestLevelNumber();
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                const double* const xlow = pgeom->getXLower();
                const hier::Index<NDIM>& idx_low = patch->getBox().lower();
                Pointer<CellData<NDIM, double> > Q_cc_data = patch->getPatchData(Q_idx);
                Pointer<SideData<NDIM, double> > Q_sc_data = patch->getPatchData(Q_idx);
                IntVector<NDIM> ghost_cells =
                    Q_cc_data ? Q_cc_data->getGhostCellWidth() : Q_sc_data->getGhostCellWidth();
                // Only print ghost cells on the left side.
                Box<NDIM> ghost_box = patch->getBox();
                ghost_box.upper(0) = ghost_box.lower(0) + (ghost_cells(0) - 1);
                ghost_box.shift(0, -ghost_cells(0));

                // Shrink the box in other directions. CF-interfaces are only filled away from corners (where enough
                // cells exist to perform interpolation)
                for (int d = 1; d < NDIM; ++d)
                {
                    ghost_box.shorten(d, 2);
                    ghost_box.shorten(d, -2);
                }

                if (Q_cc_data)
                {
                    for (CellIterator<NDIM> ci(ghost_box); ci; ci++)
                    {
                        const CellIndex<NDIM>& idx = ci();
                        VectorNd x;
                        for (int d = 0; d < NDIM; ++d)
                            x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
                        if (!IBTK::abs_equal_eps((*Q_cc_data)(idx)-Q_fcn(x), 0.0))
                        {
                            pout << "Incorrect ghost value!\n";
                            pout << "On ghost cell " << idx << "\n";
                            pout << "Computed value " << (*Q_cc_data)(idx) << "\n";
                            pout << "Exact value:   " << Q_fcn(x) << "\n";
                        }
                    }
                }
                else if (Q_sc_data)
                {
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for (SideIterator<NDIM> si(ghost_box, axis); si; si++)
                        {
                            const SideIndex<NDIM>& idx = si();
                            VectorNd x;
                            for (int d = 0; d < NDIM; ++d)
                                x[d] = xlow[d] +
                                       dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + (d == axis ? 0.0 : 0.5));
                            if (!IBTK::abs_equal_eps((*Q_sc_data)(idx)-Q_fcn(x), 0.0))
                            {
                                pout << "Incorrect ghost value!\n";
                                pout << "On ghost cell " << idx << " and axis " << axis << "\n";
                                pout << "Computed value " << (*Q_sc_data)(idx) << "\n";
                                pout << "Exact value:   " << Q_fcn(x) << "\n";
                            }
                        }
                    }
                }
            }
        }

        // Deallocate patch data
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(Q_idx);
        }

    } // cleanup dynamically allocated objects prior to shutdown
} // main
