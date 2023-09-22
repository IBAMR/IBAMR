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
#include <ibamr/CarmanKozenyDragForce.h>
#include <ibamr/EnthalpyHierarchyIntegrator.h>
#include <ibamr/HeavisideForcingFunction.h>
#include <ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/PhaseChangeDivUSourceFunction.h>
#include <ibamr/RelaxationLSMethod.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <ibamr/app_namespaces.h>

// Application
#include "LSLocateInterface.h"
#include "LevelSetInitialCondition.h"
#include "SetFluidProperties.h"
#include "SetLSProperties.h"
#include "TagInterfaceRefinementCells.h"

struct SynchronizeLevelSetCtx
{
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator;
    Pointer<CellVariable<NDIM, double> > ls_var;
    Pointer<CellVariable<NDIM, double> > H_var;
    int num_interface_cells;
};

void
synchronize_levelset_with_heaviside_fcn(int H_current_idx,
                                        Pointer<HierarchyMathOps> hier_math_ops,
                                        int /*integrator_step*/,
                                        double /*time*/,
                                        bool /*initial_time*/,
                                        bool /*regrid_time*/,
                                        void* ctx)
{
    SynchronizeLevelSetCtx* sync_ls_ctx = static_cast<SynchronizeLevelSetCtx*>(ctx);
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int ls_current_idx = var_db->mapVariableAndContextToIndex(
        sync_ls_ctx->ls_var, sync_ls_ctx->adv_diff_hier_integrator->getCurrentContext());

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* patch_dx = patch_geom->getDx();
            double vol_cell = 1.0;
            for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
            const int num_interface_cells = sync_ls_ctx->num_interface_cells;
            const double alpha = num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_current_idx);
            Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(ls_current_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double phi = (*ls_data)(ci);

                (*H_data)(ci) = IBTK::smooth_heaviside(phi, alpha);
            }
        }
    }
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
        Pointer<INSVCStaggeredHierarchyIntegrator> time_integrator = new INSVCStaggeredConservativeHierarchyIntegrator(
            "INSVCStaggeredConservativeHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSVCStaggeredConservativeHierarchyIntegrator"));

        Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator = new EnthalpyHierarchyIntegrator(
            "EnthalpyHierarchyIntegrator", app_initializer->getComponentDatabase("EnthalpyHierarchyIntegrator"));
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

        // register level set
        Pointer<CellVariable<NDIM, double> > ls_var = new CellVariable<NDIM, double>("ls_var");
        adv_diff_integrator->registerTransportedQuantity(ls_var, true);
        adv_diff_integrator->setDiffusionCoefficient(ls_var, 0.0);

        const double initial_liquid_gas_interface_position = input_db->getDouble("INITIAL_INTERFACE_POSITION");
        ;
        Pointer<RelaxationLSMethod> level_set_ops =
            new RelaxationLSMethod("RelaxationLSMethod", app_initializer->getComponentDatabase("RelaxationLSMethod"));
        LSLocateInterface* ptr_LSLocateInterface = new LSLocateInterface(
            "LSLocateInterface", adv_diff_integrator, ls_var, initial_liquid_gas_interface_position);
        level_set_ops->registerInterfaceNeighborhoodLocatingFcn(&callLSLocateInterfaceCallbackFunction,
                                                                static_cast<void*>(ptr_LSLocateInterface));
        SetLSProperties* ptr_SetLSProperties = new SetLSProperties("SetLSProperties", level_set_ops);
        adv_diff_integrator->registerResetFunction(
            ls_var, &callSetLSCallbackFunction, static_cast<void*>(ptr_SetLSProperties));

        // register liquid fraction
        Pointer<CellVariable<NDIM, double> > lf_var = new CellVariable<NDIM, double>("lf_var");
        Pointer<EnthalpyHierarchyIntegrator> enthalpy_hier_integrator = adv_diff_integrator;
        enthalpy_hier_integrator->registerLiquidFractionVariable(lf_var, true);

        // register liquid fraction gradient
        Pointer<CellVariable<NDIM, double> > lf_gradient_var = new CellVariable<NDIM, double>("lf_gradient_var", NDIM);
        enthalpy_hier_integrator->registerLiquidFractionGradientVariable(lf_gradient_var, true);

        // register Heaviside
        Pointer<CellVariable<NDIM, double> > H_var = new CellVariable<NDIM, double>("heaviside_var");
        adv_diff_integrator->registerTransportedQuantity(H_var, true);
        adv_diff_integrator->setDiffusionCoefficient(H_var, 0.0);

        // set Heaviside
        enthalpy_hier_integrator->registerHeavisideVariable(H_var);

        // register temperature
        Pointer<CellVariable<NDIM, double> > T_var = new CellVariable<NDIM, double>("Temperature");
        enthalpy_hier_integrator->registerTemperatureVariable(T_var, true);

        // set Advection velocity.
        adv_diff_integrator->setAdvectionVelocity(ls_var, time_integrator->getAdvectionVelocityVariable());
        adv_diff_integrator->setAdvectionVelocity(H_var, time_integrator->getAdvectionVelocityVariable());
        enthalpy_hier_integrator->setAdvectionVelocity(time_integrator->getAdvectionVelocityVariable());

        const ConvectiveDifferencingType ls_difference_form =
            IBAMR::string_to_enum<ConvectiveDifferencingType>(input_db->getString("LS_CONVECTIVE_FORM"));
        adv_diff_integrator->setConvectiveDifferencingType(ls_var, ls_difference_form);

        const ConvectiveDifferencingType H_difference_form =
            IBAMR::string_to_enum<ConvectiveDifferencingType>(input_db->getString("H_CONVECTIVE_FORM"));
        adv_diff_integrator->setConvectiveDifferencingType(H_var, H_difference_form);

        // set priority.
        adv_diff_integrator->setResetPriority(ls_var, 0);
        adv_diff_integrator->setResetPriority(H_var, 1);

        // set initial conditions for the variables.
        Pointer<CartGridFunction> ls_init =
            new LevelSetInitialCondition("ls_init", initial_liquid_gas_interface_position);
        adv_diff_integrator->setInitialConditions(ls_var, ls_init);

        // Since H is synchronized with ls, the initial conditions for H is not rquired.

        Pointer<CartGridFunction> T_init = new muParserCartGridFunction(
            "T_init", app_initializer->getComponentDatabase("TemperatureInitialConditions"), grid_geometry);
        enthalpy_hier_integrator->setTemperatureInitialCondition(T_var, T_init);

        Pointer<CartGridFunction> lf_init = new muParserCartGridFunction(
            "lf_init", app_initializer->getComponentDatabase("LiquidFractionInitialConditions"), grid_geometry);
        enthalpy_hier_integrator->setLiquidFractionInitialCondition(lf_var, lf_init);

        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            time_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            time_integrator->registerPressureInitialConditions(p_init);
        }

        SynchronizeLevelSetCtx sync_ls_ctx;
        sync_ls_ctx.adv_diff_hier_integrator = adv_diff_integrator;
        sync_ls_ctx.ls_var = ls_var;
        sync_ls_ctx.H_var = H_var;
        sync_ls_ctx.num_interface_cells = input_db->getDouble("NUMBER_OF_INTERFACE_CELLS");

        adv_diff_integrator->registerResetFunction(
            H_var, &synchronize_levelset_with_heaviside_fcn, static_cast<void*>(&sync_ls_ctx));

        // Setup the INS maintained material properties.
        Pointer<SideVariable<NDIM, double> > rho_sc_var = new SideVariable<NDIM, double>("rho_sc_var");
        time_integrator->registerMassDensityVariable(rho_sc_var);

        Pointer<CellVariable<NDIM, double> > mu_var = new CellVariable<NDIM, double>("mu");
        time_integrator->registerViscosityVariable(mu_var);

        Pointer<CellVariable<NDIM, double> > rho_cc_var = new CellVariable<NDIM, double>("rho_cc_var");
        enthalpy_hier_integrator->registerDensityVariable(rho_cc_var, true);

        Pointer<CellVariable<NDIM, double> > Cp_var = new CellVariable<NDIM, double>("Cp");
        enthalpy_hier_integrator->registerSpecificHeatVariable(Cp_var, true);

        // Tag cells for refinement
        const double min_tag_val = input_db->getDouble("MIN_TAG_VAL");
        const double max_tag_val = input_db->getDouble("MAX_TAG_VAL");
        TagInterfaceRefinementCells tagger(enthalpy_hier_integrator, lf_var, lf_gradient_var, min_tag_val, max_tag_val);
        enthalpy_hier_integrator->registerApplyGradientDetectorCallback(
            &callTagInterfaceRefinementCellsCallbackFunction, static_cast<void*>(&tagger));

        // Create Eulerian boundary condition specification objects (when
        // necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();

        RobinBcCoefStrategy<NDIM>* H_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("HeavisideBcCoefs"))
        {
            H_bc_coef = new muParserRobinBcCoefs(
                "H_bc_coef", app_initializer->getComponentDatabase("HeavisideBcCoefs"), grid_geometry);
            adv_diff_integrator->setPhysicalBcCoef(H_var, H_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* T_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("TemperatureBcCoefs"))
        {
            T_bc_coef = new muParserRobinBcCoefs(
                "T_bc_coef", app_initializer->getComponentDatabase("TemperatureBcCoefs"), grid_geometry);
            enthalpy_hier_integrator->setTemperaturePhysicalBcCoef(T_var, T_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* h_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("EnthalpyBcCoefs"))
        {
            h_bc_coef = new muParserRobinBcCoefs(
                "h_bc_coef", app_initializer->getComponentDatabase("EnthalpyBcCoefs"), grid_geometry);
            enthalpy_hier_integrator->setEnthalpyBcCoef(h_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* lf_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("LiquidFractionBcCoefs"))
        {
            lf_bc_coef = new muParserRobinBcCoefs(
                "lf_bc_coef", app_initializer->getComponentDatabase("LiquidFractionBcCoefs"), grid_geometry);
        }

        RobinBcCoefStrategy<NDIM>* k_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("ThermalConductivityBcCoefs"))
        {
            k_bc_coef = new muParserRobinBcCoefs(
                "k_bc_coef", app_initializer->getComponentDatabase("ThermalConductivityBcCoefs"), grid_geometry);
            enthalpy_hier_integrator->registerThermalConductivityBoundaryConditions(k_bc_coef);
        }

        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = NULL;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const std::string bc_coefs_name = "u_bc_coefs_" + std::to_string(d);

                const std::string bc_coefs_db_name = "VelocityBcCoefs_" + std::to_string(d);

                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            time_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        RobinBcCoefStrategy<NDIM>* rho_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("DensityBcCoefs"))
        {
            rho_bc_coef = new muParserRobinBcCoefs(
                "rho_bc_coef", app_initializer->getComponentDatabase("DensityBcCoefs"), grid_geometry);
            time_integrator->registerMassDensityBoundaryConditions(rho_bc_coef);
            enthalpy_hier_integrator->registerMassDensityBoundaryConditions(rho_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* mu_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("ViscosityBcCoefs"))
        {
            mu_bc_coef = new muParserRobinBcCoefs(
                "mu_bc_coef", app_initializer->getComponentDatabase("ViscosityBcCoefs"), grid_geometry);
            time_integrator->registerViscosityBoundaryConditions(mu_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* ls_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("LevelSetBcCoefs"))
        {
            ls_bc_coef = new muParserRobinBcCoefs(
                "ls_bc_coef", app_initializer->getComponentDatabase("LevelSetBcCoefs"), grid_geometry);
            adv_diff_integrator->setPhysicalBcCoef(ls_var, ls_bc_coef);
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
        SetFluidProperties* ptr_SetFluidProperties = new SetFluidProperties("SetFluidProperties",
                                                                            adv_diff_integrator,
                                                                            lf_var,
                                                                            lf_bc_coef,
                                                                            H_var,
                                                                            H_bc_coef,
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

        time_integrator->registerResetFluidDensityFcn(&callSetLiquidSolidGasDensityCallbackFunction,
                                                      static_cast<void*>(ptr_SetFluidProperties));
        time_integrator->registerResetFluidViscosityFcn(&callSetLiquidGasSolidViscosityCallbackFunction,
                                                        static_cast<void*>(ptr_SetFluidProperties));

        enthalpy_hier_integrator->registerResetDiffusionCoefficientFcn(
            &callSetLiquidSolidGasConductivityCallbackFunction, static_cast<void*>(ptr_SetFluidProperties));

        enthalpy_hier_integrator->registerResetSpecificHeatFcn(&callSetLiquidSolidGasSpecificHeatCallbackFunction,
                                                               static_cast<void*>(ptr_SetFluidProperties));

        enthalpy_hier_integrator->registerResetDensityFcn(&callSetLiquidSolidGasDensityCallbackFunction,
                                                          static_cast<void*>(ptr_SetFluidProperties));

        // Register H Div U term in the Heaviside equation.
        Pointer<CellVariable<NDIM, double> > F_var = new CellVariable<NDIM, double>("F");
        adv_diff_integrator->registerSourceTerm(F_var, true);
        Pointer<CartGridFunction> H_forcing_fcn = new HeavisideForcingFunction(
            "H_forcing_fcn", adv_diff_integrator, H_var, time_integrator->getAdvectionVelocityVariable());
        adv_diff_integrator->setSourceTermFunction(F_var, H_forcing_fcn);
        adv_diff_integrator->setSourceTerm(H_var, F_var);

        // Register source term for Div U equation.
        Pointer<CartGridFunction> Div_U_forcing_fcn =
            new PhaseChangeDivUSourceFunction("Div_U_forcing_fcn", enthalpy_hier_integrator);
        time_integrator->registerVelocityDivergenceFunction(Div_U_forcing_fcn);

        // Configure the drag force object to enforce solid velocity to be zero.
        Pointer<CarmanKozenyDragForce> drag_force =
            new CarmanKozenyDragForce("drag_force",
                                      H_var,
                                      lf_var,
                                      adv_diff_integrator,
                                      time_integrator,
                                      app_initializer->getComponentDatabase("CarmanKozenyDragForce"),
                                      /*register_for_restart*/ true);
        time_integrator->registerBrinkmanPenalizationStrategy(drag_force);

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
