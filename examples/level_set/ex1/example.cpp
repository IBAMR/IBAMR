// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2023 by the IBAMR developers
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
#include <HyperbolicLevelIntegrator.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/AdvectorExplicitPredictorPatchOps.h>
#include <ibamr/AdvectorPredictorCorrectorHyperbolicPatchOps.h>
#include <ibamr/RelaxationLSMethod.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
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
        SAMRAIPointer<CartesianGridGeometryNd> grid_geometry = new CartesianGridGeometryNd(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        SAMRAIPointer<AdvectorPredictorCorrectorHyperbolicPatchOps> hyp_patch_ops =
            new AdvectorPredictorCorrectorHyperbolicPatchOps(
                "AdvectorPredictorCorrectorHyperbolicPatchOps",
                app_initializer->getComponentDatabase("AdvectorPredictorCorrectorHyperbolicPatchOps"),
                explicit_predictor,
                grid_geometry);
        SAMRAIPointer<HyperbolicLevelIntegratorNd> hyp_level_integrator =
            new HyperbolicLevelIntegratorNd("HyperbolicLevelIntegrator",
                                            app_initializer->getComponentDatabase("HyperbolicLevelIntegrator"),
                                            hyp_patch_ops,
                                            true,
                                            using_refined_timestepping);
        SAMRAIPointer<PatchHierarchyNd> patch_hierarchy = new PatchHierarchyNd("PatchHierarchy", grid_geometry);
        SAMRAIPointer<StandardTagAndInitializeNd> error_detector =
            new StandardTagAndInitializeNd("StandardTagAndInitialize",
                                           hyp_level_integrator,
                                           app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        SAMRAIPointer<BergerRigoutsosNd> box_generator = new BergerRigoutsosNd();
        SAMRAIPointer<LoadBalancerNd> load_balancer =
            new LoadBalancerNd("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        SAMRAIPointer<GriddingAlgorithmNd> gridding_algorithm =
            new GriddingAlgorithmNd("GriddingAlgorithm",
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    error_detector,
                                    box_generator,
                                    load_balancer);
        SAMRAIPointer<TimeRefinementIntegratorNd> time_integrator =
            new TimeRefinementIntegratorNd("TimeRefinementIntegrator",
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
        SAMRAIPointer<FaceVariableNd<double> > u_var = new FaceVariableNd<double>("u");
        SAMRAIPointer<CartGridFunction> u_fcn = new muParserCartGridFunction(
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
        SAMRAIPointer<CellVariableNd<double> > Q_var = new CellVariableNd<double>("Q");
        LocationIndexRobinBcCoefsNd physical_bc_coef(
            "physical_bc_coef", app_initializer->getComponentDatabase("LocationIndexRobinBcCoefs"));
        hyp_patch_ops->registerTransportedQuantity(Q_var);
        hyp_patch_ops->setAdvectionVelocity(Q_var, u_var);
        hyp_patch_ops->setConvectiveDifferencingType(Q_var, difference_form);
        hyp_patch_ops->setPhysicalBcCoefs(Q_var, &physical_bc_coef);

        // Level set initial conditions
        SAMRAIPointer<CartGridFunction> Q_init = new muParserCartGridFunction(
            "Q_init", app_initializer->getComponentDatabase("QInitFunction"), grid_geometry);
        hyp_patch_ops->setInitialConditions(Q_var, Q_init);

        // Set up visualization plot file writer.
        SAMRAIPointer<VisItDataWriterNd> visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit) hyp_patch_ops->registerVisItDataWriter(visit_data_writer);

        // Initialize hierarchy configuration and data on all patches.
        double dt_now = time_integrator->initializeHierarchy();

        // Create inital level set
        CircularInterface circle;
        circle.R = input_db->getDouble("R");
        input_db->getDoubleArray("X0", circle.X0.data(), NDIM);

        SAMRAIPointer<VariableContext> current_ctx = hyp_level_integrator->getCurrentContext();
        SAMRAIPointer<VariableContext> scratch_ctx = hyp_level_integrator->getScratchContext();
        VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();
        const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, current_ctx);
        const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, scratch_ctx);

        SAMRAIPointer<CellVariableNd<double> > E_var = new CellVariableNd<double>("E");
        const int E_idx = var_db->registerVariableAndContext(E_var, scratch_ctx);

        // Heaviside
        SAMRAIPointer<CellVariableNd<double> > H_var = new CellVariableNd<double>("H");
        const int H_idx = var_db->registerVariableAndContext(H_var, scratch_ctx);

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            if (!level->checkAllocated(Q_scratch_idx))
                level->allocatePatchData(Q_scratch_idx, time_integrator->getIntegratorTime());
            if (!level->checkAllocated(E_idx)) level->allocatePatchData(E_idx, time_integrator->getIntegratorTime());
            if (!level->checkAllocated(H_idx)) level->allocatePatchData(H_idx, time_integrator->getIntegratorTime());
        }

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
        }

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        SAMRAIPointer<HierarchyMathOps> hier_math_ops =
            new HierarchyMathOps("HierarchyMathOps", patch_hierarchy, coarsest_ln, finest_ln);
        SAMRAIPointer<RelaxationLSMethod> level_set_ops =
            new RelaxationLSMethod("RelaxationLSMethod", app_initializer->getComponentDatabase("LevelSet"));
        level_set_ops->registerInterfaceNeighborhoodLocatingFcn(&circular_interface_neighborhood, (void*)&circle);
        level_set_ops->registerPhysicalBoundaryCondition(&physical_bc_coef);
        level_set_ops->initializeLSData(Q_scratch_idx,
                                        hier_math_ops,
                                        time_integrator->getIntegratorStep(),
                                        time_integrator->getIntegratorTime(),
                                        /*initial_time*/ true);

        // Compute L1 error from analytical solution
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevelNd::Iterator p(level); p; p++)
            {
                SAMRAIPointer<PatchNd> patch = level->getPatch(p());
                const BoxNd& patch_box = patch->getBox();
                SAMRAIPointer<CellDataNd<double> > E_data = patch->getPatchData(E_idx);
                SAMRAIPointer<CellDataNd<double> > Q_scratch_data = patch->getPatchData(Q_scratch_idx);
                SAMRAIPointer<CellDataNd<double> > H_data = patch->getPatchData(H_idx);
                for (BoxNd::Iterator it(patch_box); it; it++)
                {
                    CellIndexNd ci(it());

                    // Get physical coordinates
                    IBTK::Vector coord = IBTK::Vector::Zero();
                    SAMRAIPointer<CartesianPatchGeometryNd> patch_geom = patch->getPatchGeometry();
                    const double* patch_X_lower = patch_geom->getXLower();
                    const hier::IndexNd& patch_lower_idx = patch_box.lower();
                    const double* const patch_dx = patch_geom->getDx();
                    for (int d = 0; d < NDIM; ++d)
                    {
                        coord[d] =
                            patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
                    }
                    const double distance =
                        std::sqrt(std::pow((coord[0] - circle.X0(0)), 2.0) + std::pow((coord[1] - circle.X0(1)), 2.0)
#if (NDIM == 3)
                                  + std::pow((coord[2] - circle.X0(2)), 2.0)
#endif
                        );

                    (*E_data)(ci) = distance - circle.R;

                    const double phi =
                        -(*Q_scratch_data)(ci); // This make sure phi is positive inide the circle so that H = 1.
                    (*H_data)(ci) = IBTK::discontinuous_heaviside(phi);
                }
            }
        }

        HierarchyCellDataOpsRealNd<double> cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        cc_data_ops.subtract(E_idx, E_idx, Q_scratch_idx);
        const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();
        pout << "Error in Q after level set initialization:" << std::endl
             << "L1-norm:  " << std::setprecision(10) << cc_data_ops.L1Norm(E_idx, wgt_cc_idx) << std::endl;

        const double Numerical_volume = cc_data_ops.integral(H_idx, wgt_cc_idx);
        const double Exact_volume = M_PI * std::pow(circle.R, 2.0);
        const double error = std::abs(Numerical_volume - Exact_volume) / Exact_volume;
        pout << "Volume error of a circle:" << std::setprecision(10) << error << std::endl;

        double E_domain = 0.0;
        double E_interface = 0.0;
        int num_interface_pts = 0;
        // Compute L1 Norm for specific regions
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevelNd::Iterator p(level); p; p++)
            {
                SAMRAIPointer<PatchNd> patch = level->getPatch(p());
                const BoxNd& patch_box = patch->getBox();
                SAMRAIPointer<CellDataNd<double> > D_data = patch->getPatchData(Q_scratch_idx);
                SAMRAIPointer<CellDataNd<double> > E_data = patch->getPatchData(E_idx);
                SAMRAIPointer<CellDataNd<double> > W_data = patch->getPatchData(wgt_cc_idx);
                for (BoxNd::Iterator it(patch_box); it; it++)
                {
                    CellIndexNd ci(it());
                    const double phi = (*D_data)(ci);
                    const double err = (*E_data)(ci);
                    const double dV = (*W_data)(ci);
                    SAMRAIPointer<CartesianPatchGeometryNd> patch_geom = patch->getPatchGeometry();
                    const double* const patch_dx = patch_geom->getDx();

                    if (std::abs(phi) < 1.2 * patch_dx[0])
                    {
                        E_interface += std::abs(err) * dV;
                        num_interface_pts++;
                    }
                    if (phi > -0.8) E_domain += std::abs(err) * dV;
                }
            }
        }
        // Perform sum reduction
        num_interface_pts = IBTK_MPI::sumReduction(num_interface_pts);
        E_interface = IBTK_MPI::sumReduction(E_interface);
        E_domain = IBTK_MPI::sumReduction(E_domain);

        pout << "Error in Q near interface after level set initialization:" << std::endl
             << "L1-norm:  " << std::setprecision(10) << E_interface << std::endl;
        pout << "Error in Q in entire domain (minus center) after level set initialization:" << std::endl
             << "L1-norm:  " << std::setprecision(10) << E_domain << std::endl;
        pout << "Number of points within the interface (used to compute interface error):" << std::endl
             << num_interface_pts << std::endl;

        // Register for plotting
        visit_data_writer->registerPlotQuantity("Error", "SCALAR", E_idx);

        cc_data_ops.copyData(Q_current_idx, Q_scratch_idx);
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            if (level->checkAllocated(Q_scratch_idx)) level->deallocatePatchData(Q_scratch_idx);
        }

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

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

        if (dump_viz_data && uses_visit)
        {
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num + 1, loop_time);
        }

    } // cleanup dynamically allocated objects prior to shutdown
} // main
