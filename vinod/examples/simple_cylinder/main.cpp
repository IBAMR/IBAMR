// Filename: main.cpp
// Created by: Vinod
//
// Simple IBAMR example: Flow past a circular cylinder
//
// This example demonstrates:
// - Basic IB method setup
// - Cylinder geometry creation
// - INS solver with IB coupling
// - Visualization output

#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>

#include <petscsys.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Application includes
#include <fstream>
#include <iostream>

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must  *
 * be given on the command line. For non-restarted case, command line is:     *
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
    // Initialize IBAMR and libraries
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // Cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        SAMRAI::tbox::Pointer<IBTK::AppInitializer> app_initializer = new IBTK::AppInitializer(argc, argv, "IB.log");
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const std::string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const std::string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            SAMRAI::tbox::Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create major algorithm and data objects
        SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geometry =
            new SAMRAI::geom::CartesianGridGeometry<NDIM>(
                "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));

        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy =
            new SAMRAI::hier::PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

        SAMRAI::tbox::Pointer<SAMRAI::mesh::StandardTagAndInitialize<NDIM> > error_detector =
            new SAMRAI::mesh::StandardTagAndInitialize<NDIM>(
                "StandardTagAndInitialize",
                NULL,
                app_initializer->getComponentDatabase("StandardTagAndInitialize"));

        SAMRAI::tbox::Pointer<SAMRAI::mesh::BergerRigoutsos<NDIM> > box_generator =
            new SAMRAI::mesh::BergerRigoutsos<NDIM>();

        SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > load_balancer =
            new SAMRAI::mesh::LoadBalancer<NDIM>(
                "LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));

        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_algorithm =
            new SAMRAI::mesh::GriddingAlgorithm<NDIM>(
                "GriddingAlgorithm",
                app_initializer->getComponentDatabase("GriddingAlgorithm"),
                error_detector,
                box_generator,
                load_balancer);

        // Create Eulerian integrator (fluid solver)
        SAMRAI::tbox::Pointer<IBAMR::INSHierarchyIntegrator> navier_stokes_integrator =
            new IBAMR::INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));

        // Create Lagrangian integrator (IB method)
        SAMRAI::tbox::Pointer<IBAMR::IBHierarchyIntegrator> time_integrator =
            new IBAMR::IBExplicitHierarchyIntegrator(
                "IBHierarchyIntegrator",
                app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                new IBAMR::IBMethod("IBMethod", app_initializer->getComponentDatabase("IBMethod")),
                navier_stokes_integrator);

        // Register visualization data writers
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(app_initializer->getVisItDataWriter());
        }

        // Initialize hierarchy configuration and data
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Print input database to log file
        SAMRAI::tbox::plog << "\n\nInput database:\n";
        input_db->printClassData(SAMRAI::tbox::plog);

        // Main time step loop
        double loop_time = time_integrator->getIntegratorTime();
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;

        int iteration_num = time_integrator->getIntegratorStep();

        while (!IBTK::IBTK_MPI::equalReduction(loop_time >= loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            SAMRAI::tbox::pout << "\n";
            SAMRAI::tbox::pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            SAMRAI::tbox::pout << "At beginning of timestep # " << iteration_num << "\n";
            SAMRAI::tbox::pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            SAMRAI::tbox::pout << "At end of timestep # " << iteration_num << "\n";
            SAMRAI::tbox::pout << "Simulation time is " << loop_time << "\n";
            SAMRAI::tbox::pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            SAMRAI::tbox::pout << "\n";

            // Determine whether to write visualization data
            if (dump_viz_data && (iteration_num % viz_dump_interval == 0))
            {
                SAMRAI::tbox::pout << "\nWriting visualization files...\n\n";
                time_integrator->setupPlotData();
                app_initializer->getVisItDataWriter()->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }

            // Determine whether to write restart data
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0))
            {
                SAMRAI::tbox::pout << "\nWriting restart files...\n\n";
                SAMRAI::tbox::RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }

            // Determine whether to write timer data
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0))
            {
                SAMRAI::tbox::pout << "\nWriting timer data...\n\n";
                SAMRAI::tbox::TimerManager::getManager()->print(SAMRAI::tbox::plog);
            }
        }

        // Print final statistics
        SAMRAI::tbox::pout << "\n";
        SAMRAI::tbox::pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        SAMRAI::tbox::pout << "Computing final statistics\n";
        SAMRAI::tbox::pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        SAMRAI::tbox::pout << "\n";

    } // Cleanup dynamically allocated objects

    return 0;
} // main
