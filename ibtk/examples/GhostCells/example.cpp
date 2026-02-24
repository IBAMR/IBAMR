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

// Config files
#include <ibtk/samrai_compatibility_names.h>

#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc objects
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <SAMRAIBergerRigoutsos.h>
#include <SAMRAICartesianGridGeometry.h>
#include <SAMRAICartesianPatchGeometry.h>
#include <SAMRAICellData.h>
#include <SAMRAICellIndex.h>
#include <SAMRAICellIterator.h>
#include <SAMRAICellVariable.h>
#include <SAMRAICoarsenAlgorithm.h>
#include <SAMRAICoarsenOperator.h>
#include <SAMRAICoarsenSchedule.h>
#include <SAMRAIGriddingAlgorithm.h>
#include <SAMRAIIntVector.h>
#include <SAMRAILoadBalancer.h>
#include <SAMRAINodeData.h>
#include <SAMRAINodeIndex.h>
#include <SAMRAINodeIterator.h>
#include <SAMRAINodeVariable.h>
#include <SAMRAIPatch.h>
#include <SAMRAIPatchHierarchy.h>
#include <SAMRAIPatchLevel.h>
#include <SAMRAIPointer.h>
#include <SAMRAIRefineAlgorithm.h>
#include <SAMRAIRefineOperator.h>
#include <SAMRAIRefineSchedule.h>
#include <SAMRAIStandardTagAndInitialize.h>
#include <SAMRAIVariableDatabase.h>
#include <SAMRAIVisItDataWriter.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/muParserCartGridFunction.h>

#include "TotalAmountRefineAndCoarsen.h"

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
        SAMRAIPointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "ghost_cells.log");
        SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        SAMRAIPointer<SAMRAICartesianGridGeometry> grid_geometry = new SAMRAICartesianGridGeometry(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        SAMRAIPointer<SAMRAIPatchHierarchy> patch_hierarchy = new SAMRAIPatchHierarchy("PatchHierarchy", grid_geometry);
        SAMRAIPointer<SAMRAIStandardTagAndInitialize> error_detector = new SAMRAIStandardTagAndInitialize(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        SAMRAIPointer<SAMRAIBergerRigoutsos> box_generator = new SAMRAIBergerRigoutsos();
        SAMRAIPointer<SAMRAILoadBalancer> load_balancer =
            new SAMRAILoadBalancer("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        SAMRAIPointer<SAMRAIGriddingAlgorithm> gridding_algorithm =
            new SAMRAIGriddingAlgorithm("GriddingAlgorithm",
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
        SAMRAIPointer<SAMRAICellVariable<double>> Q_var = new SAMRAICellVariable<double>("Q");
        SAMRAIPointer<SAMRAINodeVariable<double>> N_var = new SAMRAINodeVariable<double>("Q Node");

        auto var_db = SAMRAIVariableDatabase::getDatabase();
        // Total amount patch indices. Note we need a single ghost cell width on the cell centered quantity.
        const int Q_tot_idx =
            var_db->registerVariableAndContext(Q_var, var_db->getContext("Total Amount"), SAMRAIIntVector(1));
        const int N_1_idx = var_db->registerVariableAndContext(N_var, var_db->getContext("Method 1"));
        const int N_2_idx = var_db->registerVariableAndContext(N_var, var_db->getContext("Method 2"));
        // Average amount patch index. Note we need a single ghost cell width.
        const int Q_avg_idx =
            var_db->registerVariableAndContext(Q_var, var_db->getContext("Average Amount"), SAMRAIIntVector(1));
        // Output some components
        SAMRAIPointer<SAMRAIVisItDataWriter> visit_data_writer = new SAMRAIVisItDataWriter(
            "data writer", app_initializer->getComponentDatabase("Main")->getString("viz_dump_dirname"));
        visit_data_writer->registerPlotQuantity("Q_total", "SCALAR", Q_tot_idx);
        visit_data_writer->registerPlotQuantity("Q_avg", "SCALAR", Q_avg_idx);
        visit_data_writer->registerPlotQuantity("N method 1", "SCALAR", N_1_idx);
        visit_data_writer->registerPlotQuantity("N method 2", "SCALAR", N_2_idx);

        // Allocate patch data
        int coarsest_ln = 0, finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(Q_tot_idx);
            level->allocatePatchData(N_1_idx);
            level->allocatePatchData(N_2_idx);
            level->allocatePatchData(Q_avg_idx);
        }

        // Fill in initial data
        SAMRAIPointer<CartGridFunction> init_fcn = new muParserCartGridFunction(
            "Initial data", app_initializer->getComponentDatabase("InitialFcn"), grid_geometry);
        init_fcn->setDataOnPatchHierarchy(Q_avg_idx, Q_var, patch_hierarchy, 0.0);
        // Note that CartGridFunction will fill in cell concentrations, we must manually convert this to cell amounts
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            for (SAMRAIPatchLevel::Iterator p(level); p; p++)
            {
                SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
                SAMRAIPointer<SAMRAICartesianPatchGeometry> pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                double cell_volume = 1.0;
                for (int d = 0; d < NDIM; ++d) cell_volume *= dx[d];
                SAMRAIPointer<SAMRAICellData<double>> Q_avg_data = patch->getPatchData(Q_avg_idx);
                SAMRAIPointer<SAMRAICellData<double>> Q_tot_data = patch->getPatchData(Q_tot_idx);
                for (SAMRAICellIterator ci(patch->getBox()); ci; ci++)
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
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            for (SAMRAIPatchLevel::Iterator p(level); p; p++)
            {
                SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
                SAMRAIPointer<SAMRAICartesianPatchGeometry> pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                double cell_volume = 1.0;
                for (int d = 0; d < NDIM; ++d) cell_volume *= dx[d];
                SAMRAIPointer<SAMRAICellData<double>> Q_tot_data = patch->getPatchData(Q_tot_idx);
                SAMRAIPointer<SAMRAICellData<double>> Q_avg_data = patch->getPatchData(Q_avg_idx);
                // Note we loop over the ghost box to convert the ghost cells as well.
                for (SAMRAICellIterator ci(Q_tot_data->getGhostBox()); ci; ci++)
                {
                    const SAMRAICellIndex& idx = ci();
                    (*Q_tot_data)(idx) = (*Q_avg_data)(idx)*cell_volume;
                }

                // Now interpolate to the nodes
                SAMRAIPointer<SAMRAINodeData<double>> Q_node_data = patch->getPatchData(N_1_idx);
                for (SAMRAINodeIterator ni(patch->getBox()); ni; ni++)
                {
                    const SAMRAINodeIndex& idx = ni();
                    (*Q_node_data)(idx) = 0.0;
                    const SAMRAICellIndex c_idx(idx);
                    // We need neighboring cell indices
                    for (int i = 0; i < 2; ++i)
                        for (int j = 0; j < 2; ++j)
#if (NDIM == 3)
                            for (int k = 0; k < 2; ++k)
#endif
                                (*Q_node_data)(idx) += 0.25 * (*Q_tot_data)(c_idx - SAMRAIIntVector(i,
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
        SAMRAIPointer<SAMRAICoarsenOperator> coarsen_op =
            grid_geometry->lookupCoarsenOperator(Q_var, "AMOUNT_CONSTANT_COARSEN");
        SAMRAIPointer<SAMRAICoarsenAlgorithm> coarsen_alg = new SAMRAICoarsenAlgorithm();
        coarsen_alg->registerCoarsen(Q_tot_idx, Q_tot_idx, coarsen_op);
        std::vector<SAMRAIPointer<SAMRAICoarsenSchedule>> coarsen_scheds(finest_ln + 1);
        for (int ln = coarsest_ln + 1; ln <= finest_ln; ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            SAMRAIPointer<SAMRAIPatchLevel> coarser_level = patch_hierarchy->getPatchLevel(ln - 1);
            coarsen_scheds[ln] = coarsen_alg->createSchedule(coarser_level, level);
        }
        /* Set Refine Algorithms. This interpolates data onto finer grid
         * General process:
         * 1) Set up a refine algorithm
         * 2) Register a refine operation with the algorithm
         * 3) Fill a refine schedule with the refine algorithm
         * 4) Invoke fill data() inside refine schedule
         */
        SAMRAIPointer<SAMRAIRefineOperator> refine_op =
            grid_geometry->lookupRefineOperator(Q_var, "AMOUNT_CONSTANT_REFINE");
        SAMRAIPointer<SAMRAIRefineAlgorithm> refine_alg = new SAMRAIRefineAlgorithm();
        refine_alg->registerRefine(Q_tot_idx, Q_tot_idx, Q_tot_idx, refine_op);
        std::vector<SAMRAIPointer<SAMRAIRefineSchedule>> refine_scheds(finest_ln + 1);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            refine_scheds[ln] = refine_alg->createSchedule(level, ln - 1, patch_hierarchy, nullptr);
        }
        // We refine data first, then coarsen
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln) refine_scheds[ln]->fillData(0.0);
        for (int ln = finest_ln; ln > coarsest_ln; --ln) coarsen_scheds[ln]->coarsenData();

        // Now we can interpolate to cell nodes
        // We have ghost cells filled in. We need to convert back to cell totals. We then interpoate to the nodes
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            for (SAMRAIPatchLevel::Iterator p(level); p; p++)
            {
                SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
                SAMRAIPointer<SAMRAICellData<double>> Q_tot_data = patch->getPatchData(Q_tot_idx);
                // Now interpolate to the nodes
                SAMRAIPointer<SAMRAINodeData<double>> Q_node_data = patch->getPatchData(N_2_idx);
                for (SAMRAINodeIterator ni(patch->getBox()); ni; ni++)
                {
                    const SAMRAINodeIndex& idx = ni();
                    const SAMRAICellIndex c_idx(idx);
                    (*Q_node_data)(idx) = 0.0;
                    // We need neighboring cell indices
                    for (int i = 0; i < 2; ++i)
                        for (int j = 0; j < 2; ++j)
#if (NDIM == 3)
                            for (int k = 0; k < 2; ++k)
#endif
                                (*Q_node_data)(idx) += 0.25 * (*Q_tot_data)(c_idx - SAMRAIIntVector(i,
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
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(Q_tot_idx);
            level->deallocatePatchData(N_1_idx);
            level->deallocatePatchData(N_2_idx);
            level->deallocatePatchData(Q_avg_idx);
        }

    } // cleanup dynamically allocated objects prior to shutdown
} // main
