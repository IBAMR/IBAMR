// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2022 by the IBAMR developers
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

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>

#include <LocationIndexRobinBcCoefs.h>
#include <TimeRefinementIntegrator.h>

// Set up application namespace declarations
#include "QInit.h"
#include "UFunction.h"

#include <ibamr/app_namespaces.h>

// A test based on the first advection example. Simply prints solution norms
// after integrating in time a few steps.

int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    std::ofstream output("output");
    { // cleanup dynamically allocated objects prior to shutdown
        // prevent a warning about timer initializations
        TimerManager::createManager(nullptr);

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        SAMRAIPointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "advect.log");
        SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();
        SAMRAIPointer<Database> main_db = app_initializer->getComponentDatabase("Main");

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
            output << "using subcycled timestepping.\n";
        }
        else
        {
            output << "NOT using subcycled timestepping.\n";
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
            output << "advection velocity u is discretely divergence free.\n";
        }
        else
        {
            output << "advection velocity u is NOT discretely divergence free.\n";
        }
        SAMRAIPointer<FaceVariableNd<double> > u_var = new FaceVariableNd<double>("u");
        UFunction u_fcn("UFunction", grid_geometry, app_initializer->getComponentDatabase("UFunction"));
        hyp_patch_ops->registerAdvectionVelocity(u_var);
        hyp_patch_ops->setAdvectionVelocityIsDivergenceFree(u_var, u_is_div_free);
        hyp_patch_ops->setAdvectionVelocityFunction(u_var, SAMRAIPointer<CartGridFunction>(&u_fcn, false));

        // Setup the advected quantity.
        const ConvectiveDifferencingType difference_form =
            IBAMR::string_to_enum<ConvectiveDifferencingType>(main_db->getStringWithDefault(
                "difference_form", IBAMR::enum_to_string<ConvectiveDifferencingType>(ADVECTIVE)));
        output << "solving the advection equation in "
               << IBAMR::enum_to_string<ConvectiveDifferencingType>(difference_form) << " form.\n";
        SAMRAIPointer<CellVariableNd<double> > Q_var = new CellVariableNd<double>("Q");
        QInit Q_init("QInit", grid_geometry, app_initializer->getComponentDatabase("QInit"));
        LocationIndexRobinBcCoefsNd physical_bc_coef(
            "physical_bc_coef", app_initializer->getComponentDatabase("LocationIndexRobinBcCoefs"));
        hyp_patch_ops->registerTransportedQuantity(Q_var);
        hyp_patch_ops->setAdvectionVelocity(Q_var, u_var);
        hyp_patch_ops->setConvectiveDifferencingType(Q_var, difference_form);
        hyp_patch_ops->setInitialConditions(Q_var, SAMRAIPointer<CartGridFunction>(&Q_init, false));
        hyp_patch_ops->setPhysicalBcCoefs(Q_var, &physical_bc_coef);

        // Initialize hierarchy configuration and data on all patches.
        double dt_now = time_integrator->initializeHierarchy();

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            plog << "\n";
            plog << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            plog << "At beginning of timestep # " << iteration_num << "\n";
            plog << "Simulation time is " << loop_time << "\n";

            double dt_new = time_integrator->advanceHierarchy(dt_now);
            loop_time += dt_now;
            dt_now = dt_new;

            plog << "\n";
            plog << "At end       of timestep # " << iteration_num << "\n";
            plog << "Simulation time is " << loop_time << "\n";
            plog << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            plog << "\n";

            // At specified intervals, write visualization and restart files,
            // and print out timer data.
            iteration_num += 1;
        }

        // Determine the accuracy of the computed solution.
        output << "\n"
               << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n"
               << "Computing error norms.\n\n";

        VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();
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

        HierarchyCellDataOpsRealNd<double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        hier_cc_data_ops.subtract(Q_idx, Q_idx, Q_cloned_idx);
        output << "Error in " << Q_var->getName() << " at time " << loop_time << ":\n"
               << "  L1-norm:  " << std::setprecision(10) << hier_cc_data_ops.L1Norm(Q_idx, wgt_idx) << "\n"
               << "  L2-norm:  " << hier_cc_data_ops.L2Norm(Q_idx, wgt_idx) << "\n"
               << "  max-norm: " << hier_cc_data_ops.maxNorm(Q_idx, wgt_idx) << "\n"
               << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    } // cleanup dynamically allocated objects prior to shutdown
} // main
