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
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <sstream>

#include "PhaseChangeTestUtilities.cpp"

#include <ibamr/app_namespaces.h>

// Application
#include <ibamr/PhaseChangeUtilities.h>

#include "LiquidFractionInitialCondition.cpp"
#include "LiquidFractionInitialCondition.h"
#include "TemperatureInitialCondition.cpp"
#include "TemperatureInitialCondition.h"

// Exercise the interpolation profiles through the integrator's public source
// evaluation on its live patch data. The test supplies fractions on both sides
// of independently computed transition points and temperatures around melting.
class InterpolationCheckingIntegrator : public RegridCountingIntegrator<AllenCahnHierarchyIntegrator>
{
public:
    using RegridCountingIntegrator<AllenCahnHierarchyIntegrator>::RegridCountingIntegrator;

    void setProfileCheck(const std::string& profile)
    {
        d_profile = profile;
    }

    int getProfileCheckCount() const
    {
        return d_profile_check_count;
    }

    double getProfileMaxError() const
    {
        return d_profile_max_error;
    }

    void computeDivergenceVelocitySourceTerm(int source_idx, double new_time) override
    {
        if (d_profile.empty())
        {
            AllenCahnHierarchyIntegrator::computeDivergenceVelocitySourceTerm(source_idx, new_time);
            return;
        }
        const bool hybrid = d_profile.compare(0, 7, "LINEAR_") == 0;
        const double slope = hybrid ? std::stod(d_profile.substr(7)) : 0.0;
        double lo = 0.0, hi = 1.0 / 3.0;
        for (int iteration = 0; iteration < 60; ++iteration)
        {
            const double f = 0.5 * (lo + hi);
            if (30.0 * f * (1.0 - f) * (1.0 - f) < slope)
            {
                lo = f;
            }
            else
            {
                hi = f;
            }
        }
        const double transition = hybrid ? 0.5 * (lo + hi) : 0.1;
        const double samples[] = { 0.0,         5e-11,       2e-10, transition - 1e-7,       transition + 1e-7,
                                   0.25,        0.5,         0.75,  1.0 - transition - 1e-7, 1.0 - transition + 1e-7,
                                   1.0 - 2e-10, 1.0 - 5e-11, 1.0 };
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int lf_idx = var_db->mapVariableAndContextToIndex(d_lf_var, getNewContext());
        const int T_idx = var_db->mapVariableAndContextToIndex(d_T_var, getNewContext());
        const int q_idx =
            var_db->mapVariableAndContextToIndex(var_db->getVariable("q_firstder_var"), getCurrentContext());
        for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> fraction = patch->getPatchData(lf_idx);
                Pointer<CellData<NDIM, double>> temperature = patch->getPatchData(T_idx);
                for (Box<NDIM>::Iterator it(patch->getBox()); it; it++)
                {
                    const CellIndex<NDIM> ci(it());
                    (*fraction)(ci) = samples[ci(1) % 13];
                    (*temperature)(ci) = d_T_melt + static_cast<double>(ci(0) % 3 - 1);
                }
            }
        }
        AllenCahnHierarchyIntegrator::computeDivergenceVelocitySourceTerm(source_idx, new_time);
        for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> fraction = patch->getPatchData(lf_idx);
                Pointer<CellData<NDIM, double>> temperature = patch->getPatchData(T_idx);
                Pointer<CellData<NDIM, double>> derivative = patch->getPatchData(q_idx);
                for (Box<NDIM>::Iterator it(patch->getBox()); it; it++)
                {
                    const CellIndex<NDIM> ci(it());
                    const double f = (*fraction)(ci), T = (*temperature)(ci);
                    double expected = 30.0 * f * f * (1.0 - f) * (1.0 - f);
                    if (d_profile == "QUADRATIC")
                    {
                        expected = 6.0 * f * (1.0 - f);
                    }
                    else if (hybrid && (f < transition || f > 1.0 - transition))
                    {
                        expected = slope * std::min(f, 1.0 - f);
                    }
                    if ((f >= 1.0 - 1e-10 && T <= d_T_melt) || (f <= 1e-10 && T >= d_T_melt))
                    {
                        expected = 1.0;
                    }
                    if (!std::isfinite((*derivative)(ci)))
                    {
                        TBOX_ERROR("Incorrect interpolation derivative for " << d_profile << " at f = " << f << "\n");
                    }
                    d_profile_max_error = std::max(d_profile_max_error, std::abs((*derivative)(ci)-expected));
                    ++d_profile_check_count;
                }
            }
        }
    }

private:
    std::string d_profile;
    int d_profile_check_count = 0;
    double d_profile_max_error = 0.0;
};

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
        const bool check_restart = input_db->getBoolWithDefault("check_restart", false);
        const bool check_amr = input_db->getBoolWithDefault("check_amr", false);
        const bool check_source_transfer = input_db->getBoolWithDefault("check_source_transfer", false);
        const bool from_restart = RestartManager::getManager()->isFromRestart();
        const bool check_profiles = input_db->getBoolWithDefault("check_profiles", false);

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
        time_integrator = new InterpolationCheckingIntegrator(
            "AllenCahnHierarchyIntegrator", app_initializer->getComponentDatabase("AllenCahnHierarchyIntegrator"));

        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // register liquid fraction
        Pointer<CellVariable<NDIM, double>> lf_var = new CellVariable<NDIM, double>("lf_var");
        Pointer<AllenCahnHierarchyIntegrator> ac_hier_integrator = time_integrator;
        ac_hier_integrator->registerLiquidFractionVariable(lf_var, true);
        Pointer<CellVariable<NDIM, double>> lf_gradient_var;
        if (check_amr || check_restart)
        {
            lf_gradient_var = new CellVariable<NDIM, double>("lf_gradient", NDIM);
            ac_hier_integrator->registerLiquidFractionGradientVariable(lf_gradient_var, true);
        }

        Pointer<InterpolationCheckingIntegrator> profile_integrator = time_integrator;
        if (check_profiles)
        {
            profile_integrator->setProfileCheck(input_db->getString("INTERPOLATION_FUNCTION_PROFILE"));
        }

        // register Heaviside
        Pointer<CellVariable<NDIM, double>> H_var = new CellVariable<NDIM, double>("heaviside_var");
        time_integrator->registerTransportedQuantity(H_var, true);
        time_integrator->setDiffusionCoefficient(H_var, 0.0);

        // set Heaviside
        ac_hier_integrator->registerHeavisideVariable(H_var);

        // register temperature
        Pointer<CellVariable<NDIM, double>> T_var = new CellVariable<NDIM, double>("Temperature");
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

        Pointer<CellVariable<NDIM, double>> rho_cc_var = new CellVariable<NDIM, double>("rho_cc_var");
        ac_hier_integrator->registerDensityVariable(rho_cc_var, true);

        Pointer<CellVariable<NDIM, double>> Cp_var = new CellVariable<NDIM, double>("Cp");
        ac_hier_integrator->registerSpecificHeatVariable(Cp_var, true);

        // Create Eulerian boundary condition specification objects (when
        // necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();

        RobinBcCoefStrategy<NDIM>* H_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("HeavisideBcCoefs"))
        {
            H_bc_coef = new muParserRobinBcCoefs(
                "H_bc_coef", app_initializer->getComponentDatabase("HeavisideBcCoefs"), grid_geometry);
            time_integrator->setPhysicalBcCoef(H_var, H_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* T_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("TemperatureBcCoefs"))
        {
            T_bc_coef = new muParserRobinBcCoefs(
                "T_bc_coef", app_initializer->getComponentDatabase("TemperatureBcCoefs"), grid_geometry);
            ac_hier_integrator->setTemperaturePhysicalBcCoef(T_var, T_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* lf_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("LiquidFractionBcCoefs"))
        {
            lf_bc_coef = new muParserRobinBcCoefs(
                "lf_bc_coef", app_initializer->getComponentDatabase("LiquidFractionBcCoefs"), grid_geometry);
            ac_hier_integrator->setLiquidFractionPhysicalBcCoef(lf_var, lf_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* k_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("ThermalConductivityBcCoefs"))
        {
            k_bc_coef = new muParserRobinBcCoefs(
                "k_bc_coef", app_initializer->getComponentDatabase("ThermalConductivityBcCoefs"), grid_geometry);
            ac_hier_integrator->registerThermalConductivityBoundaryConditions(k_bc_coef);
        }

        // Array for input into callback function
        const double kappa_liquid = input_db->getDouble("KAPPA_L");
        const double kappa_solid = input_db->getDouble("KAPPA_S");
        const double kappa_gas = input_db->getDouble("KAPPA_G");
        const double Cp_liquid = input_db->getDouble("CP_L");
        const double Cp_solid = input_db->getDouble("CP_S");
        const double Cp_gas = input_db->getDouble("CP_G");
        const double rho_liquid = input_db->getDouble("RHO_L");
        const double rho_solid = input_db->getDouble("RHO_S");
        const double rho_gas = input_db->getDouble("RHO_G");
        const double mu_liquid = input_db->getDouble("MU_L");
        const double mu_solid = input_db->getDouble("MU_S");
        const double mu_gas = input_db->getDouble("MU_G");

        // Callback functions can either be registered with the NS integrator, or
        // the advection-diffusion integrator
        IBAMR::PhaseChangeUtilities::SetFluidProperties set_fluid_properties("SetFluidProperties",
                                                                             time_integrator,
                                                                             H_var,

                                                                             H_bc_coef,

                                                                             lf_var,

                                                                             lf_bc_coef,
                                                                             rho_liquid,
                                                                             rho_solid,
                                                                             rho_gas,
                                                                             kappa_liquid,
                                                                             kappa_solid,
                                                                             kappa_gas,
                                                                             Cp_liquid,
                                                                             Cp_solid,
                                                                             Cp_gas,
                                                                             mu_liquid,
                                                                             mu_solid,
                                                                             mu_gas);

        ac_hier_integrator->registerResetDiffusionCoefficientFcn(
            &IBAMR::PhaseChangeUtilities::call_set_thermal_conductivity_callback,
            static_cast<void*>(&set_fluid_properties));

        ac_hier_integrator->registerResetSpecificHeatFcn(&IBAMR::PhaseChangeUtilities::call_set_specific_heat_callback,
                                                         static_cast<void*>(&set_fluid_properties));

        ac_hier_integrator->registerResetDensityFcn(&IBAMR::PhaseChangeUtilities::call_set_density_callback,
                                                    static_cast<void*>(&set_fluid_properties));

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        RefinementRegion refinement_region(input_db->getDoubleWithDefault("DT_MAX", 1.0));
        if (check_amr || check_source_transfer)
        {
            time_integrator->registerApplyGradientDetectorCallback(&tag_moving_refinement_region, &refinement_region);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        Pointer<RegridCountingIntegrator<AllenCahnHierarchyIntegrator>> counting_integrator = ac_hier_integrator;
        const int initial_reset_count = counting_integrator->getConfigurationResetCount();
        const int initial_mesh_change_count = counting_integrator->getMeshChangeCount();
        const int initial_step = time_integrator->getIntegratorStep();

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

        // Tracking the mass of the PCM.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int rho_idx = var_db->mapVariableAndContextToIndex(rho_cc_var, time_integrator->getCurrentContext());
        const int H_idx = var_db->mapVariableAndContextToIndex(H_var, time_integrator->getCurrentContext());
        const int pcm_mass_idx = var_db->registerClonedPatchDataIndex(H_var, H_idx);

        const int coarsest_ln = 0;
        int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(pcm_mass_idx, loop_time);
        }

        Pointer<HierarchyCellDataOpsReal<NDIM, double>> hier_cc_data_ops =
            new HierarchyCellDataOpsReal<NDIM, double>(patch_hierarchy, coarsest_ln, finest_ln);
        std::ostringstream results;
        if (check_source_transfer)
        {
            check_divergence_source_transfer(
                time_integrator, ac_hier_integrator, patch_hierarchy, refinement_region, results);
            if (IBTK_MPI::sumReduction(counting_integrator->getMeshChangeCount() - initial_mesh_change_count) == 0)
            {
                TBOX_ERROR("Source transfer test did not move the refined region.\n");
            }
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!check_source_transfer && !MathUtilities<double>::equalEps(loop_time, loop_time_end) &&
               time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            if (check_amr)
            {
                check_liquid_fraction_gradient(ac_hier_integrator, patch_hierarchy, lf_var, lf_gradient_var, results);
            }
            loop_time += dt;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            finest_ln = patch_hierarchy->getFinestLevelNumber();
            hier_cc_data_ops->resetLevels(coarsest_ln, finest_ln);
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
                if (!level->checkAllocated(pcm_mass_idx))
                {
                    level->allocatePatchData(pcm_mass_idx, loop_time);
                }
            }
            hier_cc_data_ops->multiply(pcm_mass_idx, rho_idx, H_idx);
            HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy, coarsest_ln, finest_ln);
            const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
            const double mass = hier_cc_data_ops->integral(pcm_mass_idx, wgt_cc_idx);
            if (!check_restart && !check_profiles)
            {
                results << std::setprecision(13) << loop_time << "\t" << mass << "\n";
            }
            else if (check_restart)
            {
                const std::vector<int> cell_indices = {
                    ac_hier_integrator->getVelocityDivergencePatchDataIndex(),
                    var_db->mapVariableAndContextToIndex(lf_gradient_var, ac_hier_integrator->getCurrentContext()),
                    var_db->mapVariableAndContextToIndex(T_var, ac_hier_integrator->getCurrentContext()),
                    var_db->mapVariableAndContextToIndex(lf_var, ac_hier_integrator->getCurrentContext()),
                    var_db->mapVariableAndContextToIndex(H_var, ac_hier_integrator->getCurrentContext()),
                    var_db->mapVariableAndContextToIndex(rho_cc_var, ac_hier_integrator->getCurrentContext()),
                    var_db->mapVariableAndContextToIndex(Cp_var, ac_hier_integrator->getCurrentContext())
                };
                check_restart_fields(
                    patch_hierarchy, cell_indices, {}, iteration_num + 1, loop_time, from_restart, results);
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
        }

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(pcm_mass_idx);
        }

        var_db->removePatchDataIndex(pcm_mass_idx);
        if (check_profiles)
        {
            if (IBTK_MPI::sumReduction(profile_integrator->getProfileCheckCount()) == 0)
            {
                TBOX_ERROR("Interpolation profile was not evaluated.\n");
            }
            results << std::setprecision(13) << "Maximum interpolation derivative error = "
                    << IBTK_MPI::maxReduction(profile_integrator->getProfileMaxError()) << '\n';
        }

        if (check_amr)
        {
            if (patch_hierarchy->getFinestLevelNumber() != 1 ||
                counting_integrator->getConfigurationResetCount() <= initial_reset_count)
            {
                TBOX_ERROR("The phase-change operators were not reset on a refined hierarchy.\n");
            }
            if (IBTK_MPI::sumReduction(counting_integrator->getMeshChangeCount() - initial_mesh_change_count) == 0)
            {
                TBOX_ERROR("The refined region did not move during integration.\n");
            }
        }
        if (check_restart)
        {
            if (time_integrator->getIntegratorStep() <= initial_step)
            {
                TBOX_ERROR("Restart test did not advance the hierarchy.\n");
            }
        }
        PIO::logOnlyNodeZero("output");
        plog << results.str();

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).

    } // cleanup dynamically allocated objects prior to shutdown
} // main
