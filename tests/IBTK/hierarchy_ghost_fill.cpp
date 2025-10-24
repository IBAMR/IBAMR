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
        Pointer<SideVariable<NDIM, double> > Q_var = new SideVariable<NDIM, double>("Q_var", 1, true);
        Pointer<SideVariable<NDIM, double> > Q_var_ip = new SideVariable<NDIM, double>("Q_var_ip", 1, true);
        Pointer<SideVariable<NDIM, double> > Q_var_oop = new SideVariable<NDIM, double>("Q_var_oop", 1, true);
        Pointer<CellVariable<NDIM, double> > Q_draw_var = new CellVariable<NDIM, double>("Q_DRAW");

        auto Q_fcn = [](const VectorNd& x) -> double
        {
            double val = 1.0;
            // for (int d = 0; d < NDIM; ++d) val += double(d + 1) * x[d];
            for (int d = 0; d < NDIM; ++d) val += x[d] + x[d] * x[d];
            return std::sin(x[0]) * std::cos(x[1]);
            return val;
        };

        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();

        auto var_db = VariableDatabase<NDIM>::getDatabase();
        int GCW = input_db->getInteger("GCW");
        // Total amount patch indices. Note we need a single ghost cell width on the cell centered quantity.
        const int Q_idx = var_db->registerVariableAndContext(Q_var, var_db->getContext("CTX"), IntVector<NDIM>(0));
        const int Q_ip_idx =
            var_db->registerVariableAndContext(Q_var_ip, var_db->getContext("CTX"), IntVector<NDIM>(GCW));
        const int Q_oop_idx =
            var_db->registerVariableAndContext(Q_var_oop, var_db->getContext("CTX"), IntVector<NDIM>(GCW));
        const int Q_draw_idx = var_db->registerVariableAndContext(Q_draw_var, var_db->getContext("CTX"));
        visit_data_writer->registerPlotQuantity("Q", "SCALAR", Q_draw_idx);
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
            level->allocatePatchData(Q_ip_idx);
            level->allocatePatchData(Q_oop_idx);
            level->allocatePatchData(Q_draw_idx);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CellData<NDIM, double> > Q_draw_data = patch->getPatchData(Q_draw_idx);
                Q_draw_data->fillAll(0.0);
            }
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
                Pointer<SideData<NDIM, double> > Q_sc_data = patch->getPatchData(Q_idx);
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
                    {
                        const SideIndex<NDIM>& idx = si();
                        VectorNd x;
                        for (int d = 0; d < NDIM; ++d)
                            x[d] =
                                xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + (d == axis ? 0.0 : 0.5));
                        (*Q_sc_data)(idx) = Q_fcn(x);
                    }
                }
            }
        }

        // Test out the filling of coarse fine interfaces. We only use the coarse fine interface operator, which fills
        // one layer of ghost cells using quadratic interpolation.
        std::string refine_type = input_db->getString("REFINE_TYPE");
        std::string coarsen_type = input_db->getString("COARSEN_TYPE");
        std::string bdry_interp_type = input_db->getString("BDRY_INTERP_TYPE");
        bool use_cf_interpolation = input_db->getBool("USE_CF_INTERPOLATION");
        bool consistent_type_2_bdry = input_db->getBool("CONSISTENT_TYPE_2_BDRY");
        bool verbose_output = input_db->getBool("VERBOSE_OUTPUT");
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_cell_comps(2);
        HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy);
        hier_sc_data_ops.copyData(Q_ip_idx, Q_idx);
        ghost_cell_comps[0] = ITC(Q_ip_idx,
                                  refine_type,
                                  use_cf_interpolation,
                                  coarsen_type,
                                  bdry_interp_type,
                                  consistent_type_2_bdry,
                                  nullptr /*bdry_conds*/);
        ghost_cell_comps[1] = ITC(Q_oop_idx,
                                  Q_idx,
                                  refine_type,
                                  use_cf_interpolation,
                                  coarsen_type,
                                  bdry_interp_type,
                                  consistent_type_2_bdry,
                                  nullptr /*bdry_conds*/);
        HierarchyGhostCellInterpolation hier_ghost_cell;
        hier_ghost_cell.initializeOperatorState(ghost_cell_comps, patch_hierarchy);
        hier_ghost_cell.fillData(0.0);

        for (int ln = 1; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                const double* const xlow = pgeom->getXLower();
                const hier::Index<NDIM>& idx_low = patch->getBox().lower();
                Pointer<SideData<NDIM, double> > Q_ip_data = patch->getPatchData(Q_ip_idx);
                Pointer<SideData<NDIM, double> > Q_oop_data = patch->getPatchData(Q_oop_idx);
                Box<NDIM> ghost_box = Q_ip_data->getGhostBox();
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    pout << "Ghost box = " << ghost_box << "\n";
                    pout << "Reg box =   " << patch->getBox() << "\n";
                    for (SideIterator<NDIM> si(ghost_box, axis); si; si++)
                    {
                        const SideIndex<NDIM>& idx = si();
                        VectorNd x;
                        for (int d = 0; d < NDIM; ++d)
                            x[d] =
                                xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + (d == axis ? 0.0 : 0.5));

                        if (!IBTK::abs_equal_eps((*Q_ip_data)(idx), (*Q_oop_data)(idx)))
                        {
                            pout << "On index " << idx << ", axis " << axis << ", and level " << ln << "\n";
                            pout << "Values are different!\n";
                            pout << "in place:     " << (*Q_ip_data)(idx) << "\n";
                            pout << "out of place: " << (*Q_oop_data)(idx) << "\n";
                            pout << "exact:        " << Q_fcn(x) << "\n";
                        }
                        else if (verbose_output)
                        {
                            pout << "On index " << idx << ", axis " << axis << ", and level " << ln << "\n";
                            pout << "Values are the same!\n";
                            pout << "in place:     " << (*Q_ip_data)(idx) << "\n";
                            pout << "out of place: " << (*Q_oop_data)(idx) << "\n";
                            pout << "exact:        " << Q_fcn(x) << "\n";
                        }
                    }
                }
            }
        }

        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);

        // Deallocate patch data
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(Q_idx);
            level->deallocatePatchData(Q_ip_idx);
            level->deallocatePatchData(Q_oop_idx);
        }

    } // cleanup dynamically allocated objects prior to shutdown
} // main
