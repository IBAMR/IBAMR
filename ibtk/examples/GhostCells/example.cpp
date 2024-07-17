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
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/muParserCartGridFunction.h>

#include "TotalAmountRefineAndCoarsen.h"

#include <CoarsenAlgorithm.h>
#include <CoarsenSchedule.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>

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
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "ghost_cells.log");
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

        // Create cell-centered and node-centered quantities, and initialize them with a function read from the input
        // file
        Pointer<CellVariableNd<double> > Q_var = new CellVariableNd<double>("Q");
        Pointer<NodeVariableNd<double> > N_var = new NodeVariableNd<double>("Q Node");

        auto var_db = VariableDatabaseNd::getDatabase();
        // Total amount patch indices. Note we need a single ghost cell width on the cell centered quantity.
        const int Q_tot_idx =
            var_db->registerVariableAndContext(Q_var, var_db->getContext("Total Amount"), IntVectorNd(1));
        const int N_1_idx = var_db->registerVariableAndContext(N_var, var_db->getContext("Method 1"));
        const int N_2_idx = var_db->registerVariableAndContext(N_var, var_db->getContext("Method 2"));
        // Average amount patch index. Note we need a single ghost cell width.
        const int Q_avg_idx =
            var_db->registerVariableAndContext(Q_var, var_db->getContext("Average Amount"), IntVectorNd(1));
        // Output some components
        Pointer<VisItDataWriterNd> visit_data_writer = new VisItDataWriterNd(
            "data writer", app_initializer->getComponentDatabase("Main")->getString("viz_dump_dirname"));
        visit_data_writer->registerPlotQuantity("Q_total", "SCALAR", Q_tot_idx);
        visit_data_writer->registerPlotQuantity("Q_avg", "SCALAR", Q_avg_idx);
        visit_data_writer->registerPlotQuantity("N method 1", "SCALAR", N_1_idx);
        visit_data_writer->registerPlotQuantity("N method 2", "SCALAR", N_2_idx);

        // Allocate patch data
        int coarsest_ln = 0, finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(Q_tot_idx);
            level->allocatePatchData(N_1_idx);
            level->allocatePatchData(N_2_idx);
            level->allocatePatchData(Q_avg_idx);
        }

        // Fill in initial data
        Pointer<CartGridFunction> init_fcn = new muParserCartGridFunction(
            "Initial data", app_initializer->getComponentDatabase("InitialFcn"), grid_geometry);
        init_fcn->setDataOnPatchHierarchy(Q_avg_idx, Q_var, patch_hierarchy, 0.0);
        // Note that CartGridFunction will fill in cell concentrations, we must manually convert this to cell amounts
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevelNd::Iterator p(level); p; p++)
            {
                Pointer<PatchNd> patch = level->getPatch(p());
                Pointer<CartesianPatchGeometryNd> pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                double cell_volume = 1.0;
                for (int d = 0; d < NDIM; ++d) cell_volume *= dx[d];
                Pointer<CellDataNd<double> > Q_avg_data = patch->getPatchData(Q_avg_idx);
                Pointer<CellDataNd<double> > Q_tot_data = patch->getPatchData(Q_tot_idx);
                for (CellIteratorNd ci(patch->getBox()); ci; ci++)
                {
                    (*Q_tot_data)(ci()) = (*Q_avg_data)(ci()) * cell_volume;
                }
            }
        }

        // Method 1: Fill in ghost cell values of the cell average. Convert the data to cell totals (including ghost
        // cell values). Interpolate to nodes. We have the cell average. We can use HierarchyGhostCellInterpolation to
        // fill in ghost cells.
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> ghost_cell_comps(1);
        ghost_cell_comps[0] = InterpolationTransactionComponent(Q_avg_idx,
                                                                "CONSERVATIVE_LINEAR_REFINE",
                                                                true,
                                                                "CONSERVATIVE_COARSEN",
                                                                "LINEAR",
                                                                false,
                                                                nullptr /*bdry_conds*/);
        HierarchyGhostCellInterpolation hier_ghost_cell;
        hier_ghost_cell.initializeOperatorState(ghost_cell_comps, patch_hierarchy);
        hier_ghost_cell.fillData(0.0);

        // We have ghost cells filled in. We need to convert back to cell totals. We then interpoate to the nodes
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevelNd::Iterator p(level); p; p++)
            {
                Pointer<PatchNd> patch = level->getPatch(p());
                Pointer<CartesianPatchGeometryNd> pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                double cell_volume = 1.0;
                for (int d = 0; d < NDIM; ++d) cell_volume *= dx[d];
                Pointer<CellDataNd<double> > Q_tot_data = patch->getPatchData(Q_tot_idx);
                Pointer<CellDataNd<double> > Q_avg_data = patch->getPatchData(Q_avg_idx);
                // Note we loop over the ghost box to convert the ghost cells as well.
                for (CellIteratorNd ci(Q_tot_data->getGhostBox()); ci; ci++)
                {
                    const CellIndexNd& idx = ci();
                    (*Q_tot_data)(idx) = (*Q_avg_data)(idx)*cell_volume;
                }

                // Now interpolate to the nodes
                Pointer<NodeDataNd<double> > Q_node_data = patch->getPatchData(N_1_idx);
                for (NodeIteratorNd ni(patch->getBox()); ni; ni++)
                {
                    const NodeIndexNd& idx = ni();
                    (*Q_node_data)(idx) = 0.0;
                    const CellIndexNd c_idx(idx);
                    // We need neighboring cell indices
                    for (int i = 0; i < 2; ++i)
                        for (int j = 0; j < 2; ++j)
#if (NDIM == 3)
                            for (int k = 0; k < 2; ++k)
#endif
                                (*Q_node_data)(idx) += 0.25 * (*Q_tot_data)(c_idx - IntVectorNd(i,
                                                                                                j
#if (NDIM == 3)
                                                                                                ,
                                                                                                k
#endif
                                                                                                ));
                }
            }
        }

        // Method 2: Fill in ghost cells directly to the total amounts. Then we can interpolate to nodes
        // Register coarsen and refine operators
        grid_geometry->addSpatialCoarsenOperator(new TotalAmountCoarsen());
        grid_geometry->addSpatialRefineOperator(new TotalAmountRefine());
        /* Set up the coasen operations. These COARSEN the data
         * General process:
         * 1) Set up a coarsen algorithm
         * 2) Register a coarsen operator with the algorithm
         * 3) Fill a coarsen schedule with the coarsen algorithm
         * 4) To actually coarsen data, use coarsen schedule -> coarsen data()
         */
        Pointer<CoarsenOperatorNd> coarsen_op = grid_geometry->lookupCoarsenOperator(Q_var, "AMOUNT_CONSTANT_COARSEN");
        Pointer<CoarsenAlgorithmNd> coarsen_alg = new CoarsenAlgorithmNd();
        coarsen_alg->registerCoarsen(Q_tot_idx, Q_tot_idx, coarsen_op);
        std::vector<Pointer<CoarsenScheduleNd> > coarsen_scheds(finest_ln + 1);
        for (int ln = coarsest_ln + 1; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            Pointer<PatchLevelNd> coarser_level = patch_hierarchy->getPatchLevel(ln - 1);
            coarsen_scheds[ln] = coarsen_alg->createSchedule(coarser_level, level);
        }
        /* Set Refine Algorithms. This interpolates data onto finer grid
         * General process:
         * 1) Set up a refine algorithm
         * 2) Register a refine operation with the algorithm
         * 3) Fill a refine schedule with the refine algorithm
         * 4) Invoke fill data() inside refine schedule
         */
        Pointer<RefineOperatorNd> refine_op = grid_geometry->lookupRefineOperator(Q_var, "AMOUNT_CONSTANT_REFINE");
        Pointer<RefineAlgorithmNd> refine_alg = new RefineAlgorithmNd();
        refine_alg->registerRefine(Q_tot_idx, Q_tot_idx, Q_tot_idx, refine_op);
        std::vector<Pointer<RefineScheduleNd> > refine_scheds(finest_ln + 1);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            refine_scheds[ln] = refine_alg->createSchedule(level, ln - 1, patch_hierarchy, nullptr);
        }
        // We refine data first, then coarsen
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln) refine_scheds[ln]->fillData(0.0);
        for (int ln = finest_ln; ln > coarsest_ln; --ln) coarsen_scheds[ln]->coarsenData();

        // Now we can interpolate to cell nodes
        // We have ghost cells filled in. We need to convert back to cell totals. We then interpoate to the nodes
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevelNd::Iterator p(level); p; p++)
            {
                Pointer<PatchNd> patch = level->getPatch(p());
                Pointer<CellDataNd<double> > Q_tot_data = patch->getPatchData(Q_tot_idx);
                // Now interpolate to the nodes
                Pointer<NodeDataNd<double> > Q_node_data = patch->getPatchData(N_2_idx);
                for (NodeIteratorNd ni(patch->getBox()); ni; ni++)
                {
                    const NodeIndexNd& idx = ni();
                    const CellIndexNd c_idx(idx);
                    (*Q_node_data)(idx) = 0.0;
                    // We need neighboring cell indices
                    for (int i = 0; i < 2; ++i)
                        for (int j = 0; j < 2; ++j)
#if (NDIM == 3)
                            for (int k = 0; k < 2; ++k)
#endif
                                (*Q_node_data)(idx) += 0.25 * (*Q_tot_data)(c_idx - IntVectorNd(i,
                                                                                                j
#if (NDIM == 3)
                                                                                                ,
                                                                                                k
#endif
                                                                                                ));
                }
            }
        }

        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);

        // Deallocate patch data
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(Q_tot_idx);
            level->deallocatePatchData(N_1_idx);
            level->deallocatePatchData(N_2_idx);
            level->deallocatePatchData(Q_avg_idx);
        }

    } // cleanup dynamically allocated objects prior to shutdown
} // main
