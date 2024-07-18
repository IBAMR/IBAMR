// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2024 by the IBAMR developers
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
#include <ibamr/FifthOrderStokesWaveGenerator.h>
#include <ibamr/FirstOrderStokesWaveGenerator.h>
#include <ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h>
#include <ibamr/IrregularWaveBcCoef.h>
#include <ibamr/IrregularWaveGenerator.h>
#include <ibamr/LevelSetUtilities.h>
#include <ibamr/RelaxationLSMethod.h>
#include <ibamr/StokesFifthOrderWaveBcCoef.h>
#include <ibamr/StokesFirstOrderWaveBcCoef.h>
#include <ibamr/StokesSecondOrderWaveBcCoef.h>
#include <ibamr/SurfaceTensionForceFunction.h>
#include <ibamr/WaveDampingFunctions.h>
#include <ibamr/WaveGenerationFunctions.h>
#include <ibamr/vc_ins_utilities.h>

#include "ibtk/IndexUtilities.h"
#include <ibtk/AppInitializer.h>
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <ibamr/app_namespaces.h>

// Application
#include "LSLocateGasInterface.h"

#include "LSLocateGasInterface.cpp"

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 Pointer<INSVCStaggeredHierarchyIntegrator> time_integrator,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

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

    // Increase maximum patch data component indices
    SAMRAIManager::setMaxNumberPatchDataEntries(2500);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Whether or not it is from a restarted run
        const bool is_from_restart = app_initializer->isFromRestart();

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
        Pointer<INSVCStaggeredHierarchyIntegrator> time_integrator;
        const string discretization_form =
            app_initializer->getComponentDatabase("Main")->getString("discretization_form");
        const bool conservative_form = (discretization_form == "CONSERVATIVE");
        if (conservative_form)
        {
            time_integrator = new INSVCStaggeredConservativeHierarchyIntegrator(
                "INSVCStaggeredConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredConservativeHierarchyIntegrator"));
        }
        else if (!conservative_form)
        {
            time_integrator = new INSVCStaggeredNonConservativeHierarchyIntegrator(
                "INSVCStaggeredNonConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredNonConservativeHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << discretization_form << "\n"
                                                   << "Valid options are: CONSERVATIVE, NON_CONSERVATIVE");
        }

        // Set up the advection diffusion hierarchy integrator
        Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator;
        const string adv_diff_solver_type = app_initializer->getComponentDatabase("Main")->getStringWithDefault(
            "adv_diff_solver_type", "SEMI_IMPLICIT");
        if (adv_diff_solver_type == "SEMI_IMPLICIT")
        {
            adv_diff_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
                "AdvDiffSemiImplicitHierarchyIntegrator",
                app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << adv_diff_solver_type << "\n"
                                                   << "Valid options are: SEMI_IMPLICIT");
        }
        time_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);

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

        // Setup level set information
        const double fluid_height = input_db->getDouble("DEPTH");
        const string& ls_name = "level_set";
        Pointer<CellVariable<NDIM, double> > phi_var = new CellVariable<NDIM, double>(ls_name);
        adv_diff_integrator->registerTransportedQuantity(phi_var);
        adv_diff_integrator->setDiffusionCoefficient(phi_var, 0.0);
        adv_diff_integrator->setAdvectionVelocity(phi_var, time_integrator->getAdvectionVelocityVariable());

        Pointer<RelaxationLSMethod> level_set_ops =
            new RelaxationLSMethod("RelaxationLSMethod", app_initializer->getComponentDatabase("RelaxationLSMethod"));

        LSLocateGasInterface setLSLocateGasInterface(
            "LSLocateGasInterface", adv_diff_integrator, phi_var, fluid_height);
        level_set_ops->registerInterfaceNeighborhoodLocatingFcn(&callLSLocateGasInterfaceCallbackFunction,
                                                                static_cast<void*>(&setLSLocateGasInterface));

        IBAMR::LevelSetUtilities::SetLSProperties setSetLSProperties("SetLSProperties", level_set_ops);
        adv_diff_integrator->registerResetFunction(
            phi_var, &IBAMR::LevelSetUtilities::setLSDataPatchHierarchy, static_cast<void*>(&setSetLSProperties));

        // Setup the INS maintained material properties.
        Pointer<Variable<NDIM> > rho_var;
        if (conservative_form)
        {
            rho_var = new SideVariable<NDIM, double>("rho");
        }
        else
        {
            rho_var = new CellVariable<NDIM, double>("rho");
        }

        time_integrator->registerMassDensityVariable(rho_var);

        Pointer<CellVariable<NDIM, double> > mu_var = new CellVariable<NDIM, double>("mu");
        time_integrator->registerViscosityVariable(mu_var);

        // Array for input into callback function
        const double rho_inside = input_db->getDouble("RHO_I");
        const double rho_outside = input_db->getDouble("RHO_O");
        const double mu_inside = input_db->getDouble("MU_I");
        const double mu_outside = input_db->getDouble("MU_O");
        const double num_interface_cells = input_db->getDouble("NUM_INTERFACE_CELLS");

        // Callback functions can either be registered with the NS integrator, or the advection-diffusion integrator
        // Note that these will set the initial conditions for density and viscosity, based on level set information
        IBAMR::VCINSUtilities::SetFluidProperties setSetFluidProperties("SetFluidProperties",
                                                                        adv_diff_integrator,
                                                                        phi_var,
                                                                        rho_outside,
                                                                        rho_inside,
                                                                        mu_outside,
                                                                        mu_inside,
                                                                        num_interface_cells);
        time_integrator->registerResetFluidDensityFcn(&IBAMR::VCINSUtilities::callSetDensityCallbackFunction,
                                                      static_cast<void*>(&setSetFluidProperties));
        time_integrator->registerResetFluidViscosityFcn(&IBAMR::VCINSUtilities::callSetViscosityCallbackFunction,
                                                        static_cast<void*>(&setSetFluidProperties));

        // Register callback function for tagging refined cells for level set data
        const double tag_thresh = input_db->getDouble("LS_TAG_ABS_THRESH");
        const double tag_min_value = -tag_thresh;
        const double tag_max_value = tag_thresh;
        IBAMR::LevelSetUtilities::TagLSRefinementCells ls_tagger(
            adv_diff_integrator, phi_var, tag_min_value, tag_max_value);
        time_integrator->registerApplyGradientDetectorCallback(&IBAMR::LevelSetUtilities::tagLSCells,
                                                               static_cast<void*>(&ls_tagger));

        // Create Eulerian initial condition specification objects.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        time_integrator->registerVelocityInitialConditions(u_init);

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            time_integrator->registerPressureInitialConditions(p_init);
        }

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        const string& wave_type = input_db->getString("WAVE_TYPE");
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = nullptr;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();

                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                if (wave_type == "FIRST_ORDER_STOKES")
                {
                    u_bc_coefs[d] = new StokesFirstOrderWaveBcCoef(
                        bc_coefs_name, d, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
                }
                else if (wave_type == "SECOND_ORDER_STOKES")
                {
                    u_bc_coefs[d] = new StokesSecondOrderWaveBcCoef(
                        bc_coefs_name, d, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
                }
                else if (wave_type == "FIFTH_ORDER_STOKES")
                {
                    u_bc_coefs[d] = new StokesFifthOrderWaveBcCoef(
                        bc_coefs_name, d, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
                }
                else if (wave_type == "IRREGULAR")
                {
                    u_bc_coefs[d] = new IrregularWaveBcCoef(
                        bc_coefs_name, d, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
                }
                else
                {
                    TBOX_ERROR("Unknown WAVE_TYPE = " << wave_type << " specified in the input file" << std::endl);
                }
            }
            time_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Create a damping zone near the channel outlet to absorb water waves via time-splitting approach.
        // This method uses a time splitting approach and modifies the fluid momentum in the
        // post processing step.
        // Register callback function for tagging refined cells for level set data
        const double x_zone_start = input_db->getDouble("X_ZONE_START");
        const double x_zone_end = input_db->getDouble("X_ZONE_END");
        const double depth = input_db->getDouble("DEPTH");
        const double alpha = input_db->getDouble("ALPHA");
        WaveDampingData wave_damper;
        wave_damper.d_x_zone_start = x_zone_start;
        wave_damper.d_x_zone_end = x_zone_end;
        wave_damper.d_depth = depth;
        wave_damper.d_alpha = alpha;
        wave_damper.d_sign_gas_phase = -1;
        wave_damper.d_ins_hier_integrator = time_integrator;
        wave_damper.d_adv_diff_hier_integrator = adv_diff_integrator;
        wave_damper.d_phi_var = phi_var;

        const string wave_damping_method = input_db->getStringWithDefault("WAVE_DAMPING_METHOD", "RELAXATION");
        if (wave_damping_method == "RELAXATION")
        {
            time_integrator->registerPostprocessIntegrateHierarchyCallback(
                &WaveDampingFunctions::callRelaxationZoneCallbackFunction, static_cast<void*>(&wave_damper));
        }
        else if (wave_damping_method == "CONSERVING")
        {
            time_integrator->registerPostprocessIntegrateHierarchyCallback(
                &WaveDampingFunctions::callConservedWaveAbsorbingCallbackFunction, static_cast<void*>(&wave_damper));
        }
        else
        {
            TBOX_ERROR("Unknown WAVE_DAMPING_METHOD = " << wave_damping_method << " specified in the input file"
                                                        << std::endl);
        }

        // Create a generating zone near the channel inlet to absorb water waves reflected by the
        // WEC towards the channel inlet. This method also uses a time splitting approach and modifies
        // the fluid momentum in the post processing step.
        // Register callback function for tagging refined cells for level set data
        const string& wave_generator_type = input_db->getStringWithDefault("WAVE_GENERATOR_TYPE", "");
        Pointer<Database> wave_db =
            app_initializer->getComponentDatabase("VelocityBcCoefs_0")->getDatabase("wave_parameters_db");
        StokesWaveGeneratorStrategy* wave_generator = nullptr;
        if (wave_generator_type == "FIRST_ORDER_STOKES")
        {
            wave_generator = new FirstOrderStokesWaveGenerator("FIRST_ORDER_STOKES", wave_db);
        }
        if (wave_generator_type == "FIFTH_ORDER_STOKES")
        {
            wave_generator = new FifthOrderStokesWaveGenerator("FIFTH_ORDER_STOKES", wave_db);
        }
        if (wave_generator_type == "IRREGULAR")
        {
            wave_generator = new IrregularWaveGenerator("IRREGULAR", wave_db);
        }
        if (wave_generator)
        {
            const double inlet_zone_start = input_db->getDouble("INLET_ZONE_START");
            const double inlet_zone_end = input_db->getDouble("INLET_ZONE_END");
            wave_generator->d_wave_gen_data.d_x_zone_start = inlet_zone_start;
            wave_generator->d_wave_gen_data.d_x_zone_end = inlet_zone_end;
            wave_generator->d_wave_gen_data.d_alpha = alpha;
            wave_generator->d_wave_gen_data.d_sign_gas_phase = -1;
            wave_generator->d_wave_gen_data.d_ins_hier_integrator = time_integrator;
            wave_generator->d_wave_gen_data.d_adv_diff_hier_integrator = adv_diff_integrator;
            wave_generator->d_wave_gen_data.d_phi_var = phi_var;

            time_integrator->registerPostprocessIntegrateHierarchyCallback(
                &WaveGenerationFunctions::callStokesWaveRelaxationCallbackFunction, static_cast<void*>(wave_generator));
        }

        RobinBcCoefStrategy<NDIM>* phi_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("PhiBcCoefs"))
        {
            phi_bc_coef = new muParserRobinBcCoefs(
                "phi_bc_coef", app_initializer->getComponentDatabase("PhiBcCoefs"), grid_geometry);
            adv_diff_integrator->setPhysicalBcCoef(phi_var, phi_bc_coef);
        }
        level_set_ops->registerPhysicalBoundaryCondition(phi_bc_coef);

        // LS initial conditions
        if (input_db->keyExists("LevelSetInitialConditions"))
        {
            Pointer<CartGridFunction> phi_init = new muParserCartGridFunction(
                "phi_init", app_initializer->getComponentDatabase("LevelSetInitialConditions"), grid_geometry);
            adv_diff_integrator->setInitialConditions(phi_var, phi_init);
        }

        RobinBcCoefStrategy<NDIM>* rho_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("RhoBcCoefs"))
        {
            rho_bc_coef = new muParserRobinBcCoefs(
                "rho_bc_coef", app_initializer->getComponentDatabase("RhoBcCoefs"), grid_geometry);
            time_integrator->registerMassDensityBoundaryConditions(rho_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* mu_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("MuBcCoefs"))
        {
            mu_bc_coef = new muParserRobinBcCoefs(
                "mu_bc_coef", app_initializer->getComponentDatabase("MuBcCoefs"), grid_geometry);
            time_integrator->registerViscosityBoundaryConditions(mu_bc_coef);
        }

        // Forcing terms
        std::vector<double> grav_const(NDIM);
        input_db->getDoubleArray("GRAV_CONST", &grav_const[0], NDIM);
        Pointer<CartGridFunction> grav_force =
            new IBAMR::VCINSUtilities::GravityForcing("GravityForcing", time_integrator, grav_const);

        Pointer<SurfaceTensionForceFunction> surface_tension_force =
            new SurfaceTensionForceFunction("SurfaceTensionForceFunction",
                                            app_initializer->getComponentDatabase("SurfaceTensionForceFunction"),
                                            adv_diff_integrator,
                                            phi_var);

        Pointer<CartGridFunctionSet> eul_forces = new CartGridFunctionSet("eulerian_forces");
        eul_forces->addFunction(surface_tension_force);
        eul_forces->addFunction(grav_force);
        time_integrator->registerBodyForceFunction(eul_forces);

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

        // Get the probe points from the input file
        Pointer<Database> probe_db = input_db->getDatabase("ProbePoints");
        const int num_probes = (probe_db->getAllKeys()).getSize();
        std::vector<std::vector<double> > probe_points;
        std::vector<std::ofstream*> probe_streams;
        probe_points.resize(num_probes);
        probe_streams.resize(num_probes);
        for (int i = 0; i < num_probes; ++i)
        {
            std::string probe_name = "output";
            probe_points[i].resize(NDIM);
            probe_db->getDoubleArray(probe_name, &probe_points[i][0], NDIM);

            if (!IBTK_MPI::getRank())
            {
                if (is_from_restart)
                {
                    probe_streams[i] = new std::ofstream(probe_name.c_str(), std::fstream::app);
                    probe_streams[i]->precision(10);
                }
                else
                {
                    probe_streams[i] = new std::ofstream(probe_name.c_str(), std::fstream::out);

                    *probe_streams[i] << "Printing level set at cell center closest to point (" << probe_points[i][0]
                                      << ", " << probe_points[i][1]
#if (NDIM == 3)
                                      << ", " << probe_points[i][2]
#endif
                                      << ") " << std::endl;
                    probe_streams[i]->precision(10);
                }
            }
        }

        // File to write to for fluid mass data and irregular wave data.
        ofstream mass_file, wave_stream;

        if (!IBTK_MPI::getRank())
        {
            mass_file.open("mass_fluid.txt");

            if (wave_generator_type == "IRREGULAR")
            {
                IrregularWaveGenerator* irregular_wave_generator =
                    dynamic_cast<IrregularWaveGenerator*>(wave_generator);

                wave_stream.open("gen_output", fstream::out);
                wave_stream.precision(10);
                irregular_wave_generator->printWaveData(wave_stream);
                wave_stream.close();
            }
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

            // Compute the fluid mass in the domain from interpolated density
            const int rho_ins_idx = time_integrator->getLinearOperatorRhoPatchDataIndex();
#if !defined(NDEBUG)
            TBOX_ASSERT(rho_ins_idx >= 0);
#endif
            const int coarsest_ln = 0;
            const int finest_ln = patch_hierarchy->getFinestLevelNumber();
            HierarchySideDataOpsReal<NDIM, double> hier_rho_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
            HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
            hier_math_ops.setPatchHierarchy(patch_hierarchy);
            hier_math_ops.resetLevels(coarsest_ln, finest_ln);
            const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
            const double mass_fluid = hier_rho_data_ops.integral(rho_ins_idx, wgt_sc_idx);

            // Write to file
            if (!IBTK_MPI::getRank())
            {
                mass_file << std::setprecision(13) << loop_time << "\t" << mass_fluid << std::endl;
            }

            // Print out the level set values at probe locations
            // Note that it will print the nearest cell center
            // Max reduction over ls_val array to ensure that only processor 0 prints
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            int phi_idx = var_db->mapVariableAndContextToIndex(phi_var, adv_diff_integrator->getCurrentContext());
            std::vector<double> ls_val(num_probes);
            for (int i = 0; i < num_probes; ++i)
            {
                ls_val[i] = -std::numeric_limits<double>::max();
                bool found_point_in_patch = false;
                for (int ln = finest_ln; ln >= coarsest_ln; --ln)
                {
                    // Get the cell index for this point
                    Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                    const CellIndex<NDIM> cell_idx =
                        IndexUtilities::getCellIndex(&probe_points[i][0], level->getGridGeometry(), level->getRatio());
                    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        Pointer<Patch<NDIM> > patch = level->getPatch(p());
                        const Box<NDIM>& patch_box = patch->getBox();
                        const bool contains_probe = patch_box.contains(cell_idx);
                        if (!contains_probe) continue;

                        // Get the level set value at this particular cell and print to stream
                        Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(phi_idx);
                        ls_val[i] = (*phi_data)(cell_idx);
                        found_point_in_patch = true;
                        if (found_point_in_patch) break;
                    }
                    if (found_point_in_patch) break;
                }
            }
            IBTK_MPI::maxReduction(&ls_val[0], num_probes);
            if (!IBTK_MPI::getRank())
            {
                for (int i = 0; i < num_probes; ++i)
                {
                    *probe_streams[i] << loop_time << '\t' << ls_val[i] << std::endl;
                }
            }

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
            if (dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
            {
                output_data(patch_hierarchy, time_integrator, iteration_num, loop_time, postproc_data_dump_dirname);
            }
        }

        // Close file
        if (!IBTK_MPI::getRank())
        {
            mass_file.close();
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

        // Cleanup other dumb pointers
        delete wave_generator;
        for (int i = 0; i < num_probes; ++i) delete probe_streams[i];

    } // cleanup dynamically allocated objects prior to shutdown
} // main

void
output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
            Pointer<INSVCStaggeredHierarchyIntegrator> time_integrator,
            const int iteration_num,
            const double loop_time,
            const string& data_dump_dirname)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;

    // Write Cartesian data.
    string file_name = data_dump_dirname + "/" + "hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, IBTK_MPI::getRank());
    file_name += temp_buf;
    Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
    hier_db->create(file_name);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    ComponentSelector hier_data;
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(time_integrator->getVelocityVariable(),
                                                           time_integrator->getCurrentContext()));
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(time_integrator->getPressureVariable(),
                                                           time_integrator->getCurrentContext()));
    patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putInteger("iteration_num", iteration_num);
    hier_db->close();
    return;
} // output_data
