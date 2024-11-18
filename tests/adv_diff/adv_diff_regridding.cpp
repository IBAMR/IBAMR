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
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include "TagInterfaceRefinementCells.h"

#include <ibamr/app_namespaces.h>

#include "TagInterfaceRefinementCells.cpp"
/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "adv_diff.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<AdvDiffHierarchyIntegrator> time_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
            "AdvDiffSemiImplicitHierarchyIntegrator",
            app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Create an initial condition specification object.
        Pointer<CartGridFunction> Q_init = new muParserCartGridFunction(
            "Q_init", app_initializer->getComponentDatabase("QInitialConditions"), grid_geometry);

        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);

        // Create boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        std::vector<RobinBcCoefStrategy<NDIM>*> Q_bc_coefs;
        if (periodic_shift.min() > 0)
        {
            Q_bc_coefs.push_back(nullptr);
        }
        else
        {
            Q_bc_coefs.push_back(new muParserRobinBcCoefs(
                "Q_bc_coefs", app_initializer->getComponentDatabase("Q_bc_coefs"), grid_geometry));
        }

        // Set up the advected and diffused quantity.
        Pointer<CellVariable<NDIM, double> > Q_var = new CellVariable<NDIM, double>("Q");
        time_integrator->registerTransportedQuantity(Q_var);
        time_integrator->setDiffusionCoefficient(Q_var, input_db->getDouble("DIFF_Q"));
        time_integrator->setInitialConditions(Q_var, Q_init);
        time_integrator->setPhysicalBcCoefs(Q_var, Q_bc_coefs);

        Pointer<FaceVariable<NDIM, double> > u_adv_var = new FaceVariable<NDIM, double>("u_adv");
        time_integrator->registerAdvectionVelocity(u_adv_var);
        time_integrator->setAdvectionVelocityFunction(u_adv_var, u_init);
        time_integrator->setAdvectionVelocity(Q_var, u_adv_var);

        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CellVariable<NDIM, double> > F_var = new CellVariable<NDIM, double>("F", 1);
            Pointer<CartGridFunction> F_fcn = new muParserCartGridFunction(
                "F_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerSourceTerm(F_var);
            time_integrator->setSourceTermFunction(F_var, F_fcn);
            time_integrator->setSourceTerm(Q_var, F_var);
        }

        // Tag cells for refinement
        const double min_tag_val = input_db->getDouble("MIN_TAG_VAL");
        const double max_tag_val = input_db->getDouble("MAX_TAG_VAL");
        TagInterfaceRefinementCells tagger(time_integrator, Q_var, min_tag_val, max_tag_val);
        time_integrator->registerApplyGradientDetectorCallback(&callTagInterfaceRefinementCellsCallbackFunction,
                                                               static_cast<void*>(&tagger));

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && uses_visit && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
        }

        // Determine the accuracy of the computed solution.
        pout << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n"
             << "Computing error norms.\n\n";

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

        const Pointer<VariableContext> Q_ctx = time_integrator->getCurrentContext();
        const int Q_idx = var_db->mapVariableAndContextToIndex(Q_var, Q_ctx);
        const int Q_cloned_idx = var_db->registerClonedPatchDataIndex(Q_var, Q_idx);

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(Q_cloned_idx, loop_time);
        }

        Q_init->setDataOnPatchHierarchy(Q_cloned_idx, Q_var, patch_hierarchy, loop_time);

        HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(coarsest_ln, finest_ln);
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        hier_cc_data_ops.subtract(Q_cloned_idx, Q_idx, Q_cloned_idx);
        const double L1_norm = hier_cc_data_ops.L1Norm(Q_cloned_idx, wgt_cc_idx);
        const double L2_norm = hier_cc_data_ops.L2Norm(Q_cloned_idx, wgt_cc_idx);
        const double max_norm = hier_cc_data_ops.maxNorm(Q_cloned_idx, wgt_cc_idx);
        if (IBTK_MPI::getRank() == 0)
        {
            std::ofstream out("output");
            out << "Error in U at time " << loop_time << ":\n"
                << "  L1-norm:  " << std::setprecision(10) << L1_norm << "\n"
                << "  L2-norm:  " << L2_norm << "\n"
                << "  max-norm: " << max_norm << "\n";
        }

        if (dump_viz_data && uses_visit)
        {
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num + 1, loop_time);
        }

        // Cleanup boundary condition specification objects (when necessary).
        delete Q_bc_coefs[0];

    } // cleanup dynamically allocated objects prior to shutdown
} // main
