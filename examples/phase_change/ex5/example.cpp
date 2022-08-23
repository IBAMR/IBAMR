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
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/AllenCahnHierarchyIntegrator.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/PhaseChangeHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <ibamr/app_namespaces.h>

// Application
#include "LiquidFractionInitialCondition.h"
#include "SetFluidProperties.h"
#include "TemperatureInitialCondition.h"

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
    // Initialize IBAMR and libraries. Deinitialization is handled by this object
    // as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // Increase maximum patch data component indices
    SAMRAIManager::setMaxNumberPatchDataEntries(2500);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && !app_initializer->getVisItDataWriter().isNull();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.

        Pointer<AdvDiffHierarchyIntegrator> time_integrator;
        time_integrator = new AllenCahnHierarchyIntegrator(
            "AllenCahnHierarchyIntegrator", app_initializer->getComponentDatabase("AllenCahnHierarchyIntegrator"));

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

        // register liquid fraction
        Pointer<CellVariable<NDIM, double> > lf_var = new CellVariable<NDIM, double>("lf_var");
        Pointer<AllenCahnHierarchyIntegrator> ac_hier_integrator = time_integrator;
        ac_hier_integrator->registerLiquidFractionVariable(lf_var, true);

        // register Heaviside
        Pointer<CellVariable<NDIM, double> > H_var = new CellVariable<NDIM, double>("heaviside_var");
        time_integrator->registerTransportedQuantity(H_var, true);
        time_integrator->setDiffusionCoefficient(H_var, 0.0);

        // set Heaviside
        ac_hier_integrator->setHeavisideVariable(H_var);

        // register temperature
        Pointer<CellVariable<NDIM, double> > T_var = new CellVariable<NDIM, double>("Temperature");
        ac_hier_integrator->registerTemperatureVariable(T_var, true);

        Pointer<CartGridFunction> H_init = new muParserCartGridFunction(
            "H_init", app_initializer->getComponentDatabase("HeavisideInitialConditions"), grid_geometry);
        time_integrator->setInitialConditions(H_var, H_init);

        const double init_liquid_solid_interface_position = input_db->getDouble("INITIAL_INTERFACE_POSITION");
        const double init_liquid_temperature = input_db->getDouble("LIQUID_TEMPERATURE");
        const double init_solid_temperature = input_db->getDouble("SOLID_TEMPERATURE");

        Pointer<CartGridFunction> T_init = new TemperatureInitialCondition(
            "T_init", init_liquid_solid_interface_position, init_liquid_temperature, init_solid_temperature);
        ac_hier_integrator->setTemperatureInitialCondition(T_var, T_init);

        Pointer<CartGridFunction> lf_init =
            new LiquidFractionInitialCondition("lf_init", init_liquid_solid_interface_position);
        ac_hier_integrator->setLiquidFractionInitialCondition(lf_var, lf_init);

        Pointer<CellVariable<NDIM, double> > rho_cc_var = new CellVariable<NDIM, double>("rho_cc_var");
        ac_hier_integrator->registerDensityVariable(rho_cc_var, true);

        Pointer<CellVariable<NDIM, double> > Cp_var = new CellVariable<NDIM, double>("Cp");
        ac_hier_integrator->registerSpecificHeatVariable(Cp_var, true);

        // Create Eulerian boundary condition specification objects (when
        // necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();

        RobinBcCoefStrategy<NDIM>* H_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("HeavisideBcCoefs"))
        {
            H_bc_coef = new muParserRobinBcCoefs(
                "H_bc_coef", app_initializer->getComponentDatabase("HeavisideBcCoefs"), grid_geometry);
            time_integrator->setPhysicalBcCoef(H_var, H_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* T_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("TemperatureBcCoefs"))
        {
            T_bc_coef = new muParserRobinBcCoefs(
                "T_bc_coef", app_initializer->getComponentDatabase("TemperatureBcCoefs"), grid_geometry);
            ac_hier_integrator->setTemperaturePhysicalBcCoef(T_var, T_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* lf_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("LiquidFractionBcCoefs"))
        {
            lf_bc_coef = new muParserRobinBcCoefs(
                "lf_bc_coef", app_initializer->getComponentDatabase("LiquidFractionBcCoefs"), grid_geometry);
            ac_hier_integrator->setLiquidFractionPhysicalBcCoef(lf_var, lf_bc_coef);
        }

        // Array for input into callback function
        const double kappa_liquid = input_db->getDouble("KAPPA_L");
        const double kappa_solid = input_db->getDouble("KAPPA_S");
        const double Cp_liquid = input_db->getDouble("CP_L");
        const double Cp_solid = input_db->getDouble("CP_S");
        const double rho_liquid = input_db->getDouble("RHO_L");
        const double rho_solid = input_db->getDouble("RHO_S");

        // Callback functions can either be registered with the NS integrator, or
        // the advection-diffusion integrator
        SetFluidProperties* ptr_SetFluidProperties = new SetFluidProperties("SetFluidProperties",
                                                                            time_integrator,
                                                                            lf_var,
                                                                            lf_bc_coef,
                                                                            rho_liquid,
                                                                            rho_solid,
                                                                            kappa_liquid,
                                                                            kappa_solid,
                                                                            Cp_liquid,
                                                                            Cp_solid);

        ac_hier_integrator->registerResetDiffusionCoefficientFcn(&callSetLiquidSolidConductivityCallbackFunction,
                                                                 static_cast<void*>(ptr_SetFluidProperties));

        ac_hier_integrator->registerResetSpecificHeatFcn(&callSetLiquidSolidSpecificHeatCallbackFunction,
                                                         static_cast<void*>(ptr_SetFluidProperties));

        ac_hier_integrator->registerResetFluidDensityFcn(&callSetLiquidSolidDensityCallbackFunction,
                                                         static_cast<void*>(ptr_SetFluidProperties));

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Remove the AppInitializer
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

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        delete ptr_SetFluidProperties;

    } // cleanup dynamically allocated objects prior to shutdown
} // main
