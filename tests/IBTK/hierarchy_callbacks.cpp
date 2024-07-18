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
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>

#include <LocationIndexRobinBcCoefs.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

struct DummyObject
{
public:
    DummyObject(std::string object_name) : d_object_name(std::move(object_name))
    {
    }
    std::string d_object_name;
};

void
preprocessCallback(double /*current_time*/, double /*new_time*/, int /*num_cycles*/, void* ctx)
{
    auto test = static_cast<DummyObject*>(ctx);
    if (test) pout << "Executing preprocess callback\n";
}

void
postprocessCallback(double /*current_time*/,
                    double /*new_time*/,
                    bool /*skip_synchronization_new_state_data*/,
                    int /*num_cycles*/,
                    void* ctx)
{
    auto test = static_cast<DummyObject*>(ctx);
    if (test) pout << "Executing postprocess callback\n";
}

void
integrateHierarchyCallback(double /*current_time*/, double /*new_time*/, int /*cycle_num*/, void* ctx)
{
    auto test = static_cast<DummyObject*>(ctx);
    if (test) pout << "Executing integrate hierarchy callback\n";
}

void
gradientDetectorCallback(SAMRAIPointer<BasePatchHierarchyNd> /*hierarchy*/,
                         int /*level_num*/,
                         double /*error_data_time*/,
                         int /*tag_index*/,
                         bool /*initial_time*/,
                         bool /*uses_richardson_extrapolation_too*/,
                         void* ctx)
{
    auto test = static_cast<DummyObject*>(ctx);
    if (test) pout << "Executing gradient detector callback\n";
}
void
regridHierarchyCallback(SAMRAIPointer<BasePatchHierarchyNd> /*hierarchy*/,
                        double /*data_time*/,
                        bool /*initial_time*/,
                        void* ctx)
{
    auto test = static_cast<DummyObject*>(ctx);
    if (test) pout << "Executing regrid hierarchy callback\n";
}

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
      // prevent a warning about timer initialization
        TimerManager::createManager(nullptr);

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        auto app_initializer = make_samrai_shared<AppInitializer>(argc, argv, "adv_diff.log");
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

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        SAMRAIPointer<AdvDiffHierarchyIntegrator> time_integrator =
            make_samrai_shared<AdvDiffSemiImplicitHierarchyIntegrator>(
                "AdvDiffSemiImplicitHierarchyIntegrator",
                app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
        auto grid_geometry = make_samrai_shared<CartesianGridGeometryNd>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        auto patch_hierarchy = make_samrai_shared<PatchHierarchyNd>("PatchHierarchy", grid_geometry);
        auto error_detector = make_samrai_shared<StandardTagAndInitializeNd>(
            "StandardTagAndInitialize",
            time_integrator,
            app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        auto box_generator = make_samrai_shared<BergerRigoutsosNd>();
        auto load_balancer =
            make_samrai_shared<LoadBalancerNd>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        auto gridding_algorithm =
            make_samrai_shared<GriddingAlgorithmNd>("GriddingAlgorithm",
                                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                                    error_detector,
                                                    box_generator,
                                                    load_balancer);

        // Setup the advection velocity.
        SAMRAIPointer<FaceVariableNd<double> > u_var = make_samrai_shared<FaceVariableNd<double> >("u");
        SAMRAIPointer<CartGridFunction> u_fcn = make_samrai_shared<muParserCartGridFunction>(
            "UFunction", app_initializer->getComponentDatabase("UFunction"), grid_geometry);
        const bool u_is_div_free = true;
        time_integrator->registerAdvectionVelocity(u_var);
        time_integrator->setAdvectionVelocityIsDivergenceFree(u_var, u_is_div_free);
        time_integrator->setAdvectionVelocityFunction(u_var, u_fcn);

        // Setup the advected and diffused quantity.
        const ConvectiveDifferencingType difference_form =
            IBAMR::string_to_enum<ConvectiveDifferencingType>(main_db->getStringWithDefault(
                "difference_form", IBAMR::enum_to_string<ConvectiveDifferencingType>(ADVECTIVE)));
        SAMRAIPointer<CellVariableNd<double> > Q_var = make_samrai_shared<CellVariableNd<double> >("Q");
        SAMRAIPointer<CartGridFunction> Q_init = make_samrai_shared<muParserCartGridFunction>(
            "QInit", app_initializer->getComponentDatabase("QInit"), grid_geometry);
        LocationIndexRobinBcCoefsNd physical_bc_coef(
            "physical_bc_coef", app_initializer->getComponentDatabase("LocationIndexRobinBcCoefs"));
        const double kappa = app_initializer->getComponentDatabase("QInit")->getDouble("kappa");
        time_integrator->registerTransportedQuantity(Q_var);
        time_integrator->setAdvectionVelocity(Q_var, u_var);
        time_integrator->setDiffusionCoefficient(Q_var, kappa);
        time_integrator->setConvectiveDifferencingType(Q_var, difference_form);
        time_integrator->setInitialConditions(Q_var, Q_init);
        time_integrator->setPhysicalBcCoef(Q_var, &physical_bc_coef);

        DummyObject callbackObject("temp_object");
        time_integrator->registerPreprocessIntegrateHierarchyCallback(preprocessCallback,
                                                                      static_cast<void*>(&callbackObject));
        time_integrator->registerIntegrateHierarchyCallback(integrateHierarchyCallback,
                                                            static_cast<void*>(&callbackObject));
        time_integrator->registerPostprocessIntegrateHierarchyCallback(postprocessCallback,
                                                                       static_cast<void*>(&callbackObject));
        time_integrator->registerApplyGradientDetectorCallback(gradientDetectorCallback,
                                                               static_cast<void*>(&callbackObject));
        time_integrator->registerRegridHierarchyCallback(regridHierarchyCallback, static_cast<void*>(&callbackObject));

        // Set up visualization plot file writer.
        SAMRAIPointer<VisItDataWriterNd> visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Close the restart manager.
        RestartManager::getManager()->closeRestartFile();

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
            // and print out timer data.
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
                pout << "\nWriting restart files...\n\nn";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
        }

    } // cleanup dynamically allocated objects prior to shutdown

} // main
