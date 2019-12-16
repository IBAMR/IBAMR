// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2019 by the IBAMR developers
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
#include <IBAMR_config.h>
#include <IBTK_config.h>

#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <HyperbolicLevelIntegrator.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/AdvectorExplicitPredictorPatchOps.h>
#include <ibamr/AdvectorPredictorCorrectorHyperbolicPatchOps.h>
#include <ibamr/FastSweepingLSMethod.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/muParserCartGridFunction.h>

#include <LocationIndexRobinBcCoefs.h>
#include <TimeRefinementIntegrator.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Application specific includes.
#include "LSLocateInterface.h"

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
    // Initialize MPI and SAMRAI.
    SAMRAI_MPI::init(&argc, &argv);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "advect.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Get solver configuration options.
        bool using_refined_timestepping = false;
        if (main_db->keyExists("timestepping"))
        {
            string timestepping_method = main_db->getString("timestepping");
            if (timestepping_method == "SYNCHRONIZED")
            {
                using_refined_timestepping = false;
            }
            else
            {
                using_refined_timestepping = true;
            }
        }
        if (using_refined_timestepping)
        {
            pout << "using subcycled timestepping.\n";
        }
        else
        {
            pout << "NOT using subcycled timestepping.\n";
        }

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<AdvectorExplicitPredictorPatchOps> explicit_predictor = new AdvectorExplicitPredictorPatchOps(
            "AdvectorExplicitPredictorPatchOps",
            app_initializer->getComponentDatabase("AdvectorExplicitPredictorPatchOps"));
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<AdvectorPredictorCorrectorHyperbolicPatchOps> hyp_patch_ops =
            new AdvectorPredictorCorrectorHyperbolicPatchOps(
                "AdvectorPredictorCorrectorHyperbolicPatchOps",
                app_initializer->getComponentDatabase("AdvectorPredictorCorrectorHyperbolicPatchOps"),
                explicit_predictor,
                grid_geometry);
        Pointer<HyperbolicLevelIntegrator<NDIM> > hyp_level_integrator =
            new HyperbolicLevelIntegrator<NDIM>("HyperbolicLevelIntegrator",
                                                app_initializer->getComponentDatabase("HyperbolicLevelIntegrator"),
                                                hyp_patch_ops,
                                                true,
                                                using_refined_timestepping);
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               hyp_level_integrator,
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
        Pointer<TimeRefinementIntegrator<NDIM> > time_integrator =
            new TimeRefinementIntegrator<NDIM>("TimeRefinementIntegrator",
                                               app_initializer->getComponentDatabase("TimeRefinementIntegrator"),
                                               patch_hierarchy,
                                               hyp_level_integrator,
                                               gridding_algorithm);

        // Setup the advection velocity.
        const bool u_is_div_free = main_db->getBoolWithDefault("u_is_div_free", false);
        if (u_is_div_free)
        {
            pout << "advection velocity u is discretely divergence free.\n";
        }
        else
        {
            pout << "advection velocity u is NOT discretely divergence free.\n";
        }
        Pointer<FaceVariable<NDIM, double> > u_var = new FaceVariable<NDIM, double>("u");
        Pointer<CartGridFunction> u_fcn = new muParserCartGridFunction(
            "u_fcn", app_initializer->getComponentDatabase("AdvectionVelocityFunction"), grid_geometry);
        hyp_patch_ops->registerAdvectionVelocity(u_var);
        hyp_patch_ops->setAdvectionVelocityIsDivergenceFree(u_var, u_is_div_free);
        hyp_patch_ops->setAdvectionVelocityFunction(u_var, u_fcn);

        // Setup the advected quantity.
        const ConvectiveDifferencingType difference_form =
            IBAMR::string_to_enum<ConvectiveDifferencingType>(main_db->getStringWithDefault(
                "difference_form", IBAMR::enum_to_string<ConvectiveDifferencingType>(ADVECTIVE)));
        pout << "solving the advection equation in "
             << IBAMR::enum_to_string<ConvectiveDifferencingType>(difference_form) << " form.\n";
        Pointer<CellVariable<NDIM, double> > Q_var = new CellVariable<NDIM, double>("Q");
        LocationIndexRobinBcCoefs<NDIM> physical_bc_coef(
            "physical_bc_coef", app_initializer->getComponentDatabase("LocationIndexRobinBcCoefs"));
        hyp_patch_ops->registerTransportedQuantity(Q_var);
        hyp_patch_ops->setAdvectionVelocity(Q_var, u_var);
        hyp_patch_ops->setConvectiveDifferencingType(Q_var, difference_form);
        hyp_patch_ops->setPhysicalBcCoefs(Q_var, &physical_bc_coef);

        // Set up visualization plot file writer.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit) hyp_patch_ops->registerVisItDataWriter(visit_data_writer);

        // Initialize hierarchy configuration and data on all patches.
        double dt_now = time_integrator->initializeHierarchy();

        // Create inital level set
        CircularInterface circle;
        circle.R = 0.2;
        circle.X0(0) = 0.5;
        circle.X0(1) = 0.5;
#if (NDIM == 3)
        circle.X0(2) = 0.5;
#endif
        input_db->getDoubleWithDefault("R", circle.R);
        if (input_db->keyExists("X0"))
        {
            input_db->getDoubleArray("X0", circle.X0.data(), NDIM);
        }

        Pointer<VariableContext> current_ctx = hyp_level_integrator->getCurrentContext();
        Pointer<VariableContext> scratch_ctx = hyp_level_integrator->getScratchContext();
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, current_ctx);
        const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, scratch_ctx);
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            if (level->checkAllocated(Q_scratch_idx)) continue;
            level->allocatePatchData(Q_scratch_idx, time_integrator->getIntegratorTime());
        }

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        Pointer<HierarchyMathOps> hier_math_ops =
            new HierarchyMathOps("HierarchyMathOps", patch_hierarchy, coarsest_ln, finest_ln);
        Pointer<FastSweepingLSMethod> level_set_ops =
            new FastSweepingLSMethod("FastSweepingLSMethod", app_initializer->getComponentDatabase("LevelSet"));
        level_set_ops->registerInterfaceNeighborhoodLocatingFcn(&circular_interface_neighborhood, (void*)&circle);
        level_set_ops->initializeLSData(Q_scratch_idx,
                                        hier_math_ops,
                                        time_integrator->getIntegratorStep(),
                                        time_integrator->getIntegratorTime(),
                                        /*initial_time*/ true);

        HierarchyCellDataOpsReal<NDIM, double> cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        cc_data_ops.copyData(Q_current_idx, Q_scratch_idx);
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            if (!level->checkAllocated(Q_scratch_idx)) continue;
            level->deallocatePatchData(Q_scratch_idx);
        }

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
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            double dt_new = time_integrator->advanceHierarchy(dt_now);
            loop_time += dt_now;
            dt_now = dt_new;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // At specified intervals, write visualization and restart files,
            // and print out timer data.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && uses_visit && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\nn";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
        }

        if (dump_viz_data && uses_visit)
        {
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num + 1, loop_time);
        }

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    SAMRAI_MPI::finalize();
} // main
