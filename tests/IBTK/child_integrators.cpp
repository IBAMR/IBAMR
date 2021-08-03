// ---------------------------------------------------------------------
//
// Copyright (c) 2021 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>

#include <petscsys.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <LocationIndexRobinBcCoefs.h>
#include <SAMRAI_config.h>
#include <StandardTagAndInitialize.h>
#include <VisItDataWriter.h>

#include <ibamr/app_namespaces.h>

class DummyAdvDiffIntegrator : public AdvDiffHierarchyIntegrator
{
public:
    DummyAdvDiffIntegrator(std::string object_name, Pointer<Database> input_db)
        : AdvDiffHierarchyIntegrator(std::move(object_name), input_db, true)
    {
    }

    ~DummyAdvDiffIntegrator() = default;

    void integrateHierarchy(double current_time, double new_time, int cycle_num = 0) override
    {
        HierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);
        return;
    }

protected:
    void regridHierarchyBeginSpecialized() override
    {
        plog << d_object_name << ": regridHierarchyBeginSpecialized()\n";
    }
    void regridHierarchyEndSpecialized() override
    {
        plog << d_object_name << ": regridHierarchyEndSpecialized()\n";
    }
    double getMinimumTimeStepSizeSpecialized() override
    {
        plog << d_object_name << ": getMinimumTimeStepSizeSpecialized()\n";
        return std::numeric_limits<double>::min();
    }
    double getMaximumTimeStepSizeSpecialized() override
    {
        plog << d_object_name << ": getMaximumTimeStepSizeSpecialized()\n";
        return std::numeric_limits<double>::max();
    }
    void synchronizeHierarchyDataSpecialized(VariableContextType /*ctx_type*/) override
    {
        plog << d_object_name << ": synchronizeHierarchyDataSpecialized()\n";
    }
    void resetTimeDependentHierarchyDataSpecialized(const double /*new_time*/) override
    {
        plog << d_object_name << ": resetTimeDependentHierarchyDataSpecialized()\n";
    }
    void resetIntegratorToPreadvanceStateSpecialized() override
    {
        plog << d_object_name << ": resetIntegratorToPreadvanceStateSpecialized()\n";
    }
    bool atRegridPointSpecialized() const override
    {
        plog << d_object_name << ": atRegridPointSpecialized()\n";
        return false;
    }
    void setupPlotDataSpecialized() override
    {
        plog << d_object_name << ": setupPlotDataSpecialized()\n";
    }
    void initializeCompositeHierarchyDataSpecialized(double /*init_data_time*/, bool /*initial_time*/) override
    {
        plog << d_object_name << ": initializeCompositeHierarchyDataSpecialized()\n";
    }
    void initializeLevelDataSpecialized(Pointer<BasePatchHierarchy<NDIM> > /*hierarchy*/,
                                        int /*level_number*/,
                                        double /*init_data_time*/,
                                        bool /*can_be_refined*/,
                                        bool /*initial_time*/,
                                        Pointer<BasePatchLevel<NDIM> > /*old_level*/,
                                        bool /*allocate_data*/) override
    {
        plog << d_object_name << ": initializeLevelDataSpecialized()\n";
    }
    void resetHierarchyConfigurationSpecialized(Pointer<BasePatchHierarchy<NDIM> > /*hierarchy*/,
                                                int /*coarsest_level*/,
                                                int /*finest_level*/) override
    {
        plog << d_object_name << ": resetHierarchyConfigurationSpecialized()\n";
    }
    void putToDatabaseSpecialized(Pointer<Database> /*db*/) override
    {
        plog << d_object_name << ": putToDatabaseSpecialized()\n";
    }
    void addWorkloadEstimate(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/, const int /*workload_data_idx*/) override
    {
        plog << d_object_name << ": addWorkloadEstimate()\n";
    }
    void executePreprocessIntegrateHierarchyCallbackFcns(double /*current_time*/,
                                                         double /*new_time*/,
                                                         int /*num_cycles*/) override
    {
        plog << d_object_name << ": executePreprocessIntegrateHierarchyCallbackFcns()\n";
    }
};

void
generate_structure(const unsigned int& /*strct_num*/,
                   const int& /*ln*/,
                   int& num_vertices,
                   std::vector<IBTK::Point>& vertex_posn)
{
    num_vertices = 1;
    vertex_posn.resize(num_vertices);
    vertex_posn[0] = Point(0.5, 0.5);
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

#ifndef IBTK_HAVE_SILO
    // Suppress warnings caused by running without silo
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);
#endif

    { // cleanup dynamically allocated objects prior to shutdown
      // prevent a warning about timer initialization
        TimerManager::createManager(nullptr);

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "adv_diff.log");
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

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator =
            new DummyAdvDiffIntegrator("DummyAdvDiff", app_initializer->getComponentDatabase("DummyAdvDiff"));
        Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);
        Pointer<IBMethod> ib_method_ops = new IBMethod("IBMethod", app_initializer->getComponentDatabase("IBMethod"));
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);
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
        // Configure the IB solver.
        Pointer<IBRedundantInitializer> ib_initializer = new IBRedundantInitializer(
            "IBRedundantInitializer", app_initializer->getComponentDatabase("IBRedundantInitializer"));
        std::vector<std::string> struct_list_vec = { "dummy" };
        ib_initializer->setStructureNamesOnLevel(input_db->getInteger("MAX_LEVELS") - 1, struct_list_vec);
        ib_initializer->registerInitStructureFunction(generate_structure);
        ib_method_ops->registerLInitStrategy(ib_initializer);

        // Set up visualization plot file writer.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
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
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
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

    } // cleanup dynamically allocated objects prior to shutdown

} // main
