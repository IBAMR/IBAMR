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
#include <ibtk/samrai_compatibility_names.h>

#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <SAMRAIBergerRigoutsos.h>
#include <SAMRAICartesianGridGeometry.h>
#include <SAMRAICellVariable.h>
#include <SAMRAIFaceVariable.h>
#include <SAMRAIGriddingAlgorithm.h>
#include <SAMRAIHierarchyCellDataOpsReal.h>
#include <SAMRAIHyperbolicLevelIntegrator.h>
#include <SAMRAILoadBalancer.h>
#include <SAMRAILocationIndexRobinBcCoefs.h>
#include <SAMRAIPatchHierarchy.h>
#include <SAMRAIPointer.h>
#include <SAMRAIStandardTagAndInitialize.h>
#include <SAMRAITimeRefinementIntegrator.h>
#include <SAMRAIVariableDatabase.h>
#include <SAMRAIVisItDataWriter.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/AdvectorExplicitPredictorPatchOps.h>
#include <ibamr/AdvectorPredictorCorrectorHyperbolicPatchOps.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>

// Set up application namespace declarations
#include "QInit.h"
#include "UFunction.h"

#include <ibamr/app_namespaces.h>

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
        SAMRAIPointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "advect.log");
        SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();
        SAMRAIPointer<Database> main_db = app_initializer->getComponentDatabase("Main");

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
        SAMRAIPointer<AdvectorExplicitPredictorPatchOps> explicit_predictor = new AdvectorExplicitPredictorPatchOps(
            "AdvectorExplicitPredictorPatchOps",
            app_initializer->getComponentDatabase("AdvectorExplicitPredictorPatchOps"));
        SAMRAIPointer<SAMRAICartesianGridGeometry> grid_geometry = new SAMRAICartesianGridGeometry(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        SAMRAIPointer<AdvectorPredictorCorrectorHyperbolicPatchOps> hyp_patch_ops =
            new AdvectorPredictorCorrectorHyperbolicPatchOps(
                "AdvectorPredictorCorrectorHyperbolicPatchOps",
                app_initializer->getComponentDatabase("AdvectorPredictorCorrectorHyperbolicPatchOps"),
                explicit_predictor,
                grid_geometry);
        SAMRAIPointer<SAMRAIHyperbolicLevelIntegrator> hyp_level_integrator =
            new SAMRAIHyperbolicLevelIntegrator("HyperbolicLevelIntegrator",
                                                app_initializer->getComponentDatabase("HyperbolicLevelIntegrator"),
                                                hyp_patch_ops,
                                                true,
                                                using_refined_timestepping);
        SAMRAIPointer<SAMRAIPatchHierarchy> patch_hierarchy = new SAMRAIPatchHierarchy("PatchHierarchy", grid_geometry);
        SAMRAIPointer<SAMRAIStandardTagAndInitialize> error_detector =
            new SAMRAIStandardTagAndInitialize("StandardTagAndInitialize",
                                               hyp_level_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        SAMRAIPointer<SAMRAIBergerRigoutsos> box_generator = new SAMRAIBergerRigoutsos();
        SAMRAIPointer<SAMRAILoadBalancer> load_balancer =
            new SAMRAILoadBalancer("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        SAMRAIPointer<SAMRAIGriddingAlgorithm> gridding_algorithm =
            new SAMRAIGriddingAlgorithm("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);
        SAMRAIPointer<SAMRAITimeRefinementIntegrator> time_integrator =
            new SAMRAITimeRefinementIntegrator("TimeRefinementIntegrator",
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
        SAMRAIPointer<SAMRAIFaceVariable<double>> u_var = new SAMRAIFaceVariable<double>("u");
        UFunction u_fcn("UFunction", grid_geometry, app_initializer->getComponentDatabase("UFunction"));
        hyp_patch_ops->registerAdvectionVelocity(u_var);
        hyp_patch_ops->setAdvectionVelocityIsDivergenceFree(u_var, u_is_div_free);
        hyp_patch_ops->setAdvectionVelocityFunction(u_var, SAMRAIPointer<CartGridFunction>(&u_fcn, false));

        // Setup the advected quantity.
        const ConvectiveDifferencingType difference_form =
            IBAMR::string_to_enum<ConvectiveDifferencingType>(main_db->getStringWithDefault(
                "difference_form", IBAMR::enum_to_string<ConvectiveDifferencingType>(ADVECTIVE)));
        pout << "solving the advection equation in "
             << IBAMR::enum_to_string<ConvectiveDifferencingType>(difference_form) << " form.\n";
        SAMRAIPointer<SAMRAICellVariable<double>> Q_var = new SAMRAICellVariable<double>("Q");
        QInit Q_init("QInit", grid_geometry, app_initializer->getComponentDatabase("QInit"));
        SAMRAILocationIndexRobinBcCoefs physical_bc_coef(
            "physical_bc_coef", app_initializer->getComponentDatabase("LocationIndexRobinBcCoefs"));
        hyp_patch_ops->registerTransportedQuantity(Q_var);
        hyp_patch_ops->setAdvectionVelocity(Q_var, u_var);
        hyp_patch_ops->setConvectiveDifferencingType(Q_var, difference_form);
        hyp_patch_ops->setInitialConditions(Q_var, SAMRAIPointer<CartGridFunction>(&Q_init, false));
        hyp_patch_ops->setPhysicalBcCoefs(Q_var, &physical_bc_coef);

        // Set up visualization plot file writer.
        SAMRAIPointer<SAMRAIVisItDataWriter> visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit) hyp_patch_ops->registerVisItDataWriter(visit_data_writer);

        // Initialize hierarchy configuration and data on all patches.
        double dt_now = time_integrator->initializeHierarchy();

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
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
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

        // Determine the accuracy of the computed solution.
        pout << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n"
             << "Computing error norms.\n\n";

        SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
        const SAMRAIPointer<VariableContext> Q_ctx = hyp_level_integrator->getCurrentContext();
        const int Q_idx = var_db->mapVariableAndContextToIndex(Q_var, Q_ctx);
        const int Q_cloned_idx = var_db->registerClonedPatchDataIndex(Q_var, Q_idx);

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(Q_cloned_idx, loop_time);
        }
        Q_init.setDataOnPatchHierarchy(Q_cloned_idx, Q_var, patch_hierarchy, loop_time);

        HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(coarsest_ln, finest_ln);
        const int wgt_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        SAMRAIHierarchyCellDataOpsReal<double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        hier_cc_data_ops.subtract(Q_idx, Q_idx, Q_cloned_idx);
        pout << "Error in " << Q_var->getName() << " at time " << loop_time << ":\n"
             << "  L1-norm:  " << std::setprecision(10) << hier_cc_data_ops.L1Norm(Q_idx, wgt_idx) << "\n"
             << "  L2-norm:  " << hier_cc_data_ops.L2Norm(Q_idx, wgt_idx) << "\n"
             << "  max-norm: " << hier_cc_data_ops.maxNorm(Q_idx, wgt_idx) << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        if (dump_viz_data && uses_visit)
        {
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num + 1, loop_time);
        }

    } // cleanup dynamically allocated objects prior to shutdown
} // main
