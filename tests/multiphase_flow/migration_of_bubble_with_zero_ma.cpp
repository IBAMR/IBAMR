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
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/MarangoniSurfaceTensionForceFunction.h>
#include <ibamr/RelaxationLSMethod.h>
#include <ibamr/SurfaceTensionForceFunction.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <ibamr/app_namespaces.h>

// Application
#include "LSLocateCircularInterface.h"
#include "SetFluidProperties.h"
#include "SetLSProperties.h"

struct MaskSurfaceTensionForceCtx
{
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator;
    Pointer<INSVCStaggeredHierarchyIntegrator> ins_hier_integrator;
    double rho_liquid;
    double rho_gas;
};

void
mask_surface_tension_force(int F_idx,
                           Pointer<HierarchyMathOps> hier_math_ops,
                           int /*integrator_step*/,
                           double /*time*/,
                           double /*current_time*/,
                           double /*new_time*/,
                           void* ctx)
{
    MaskSurfaceTensionForceCtx* mask_surface_tension_force_ctx = static_cast<MaskSurfaceTensionForceCtx*>(ctx);
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    int rho_idx = mask_surface_tension_force_ctx->ins_hier_integrator->getLinearOperatorRhoPatchDataIndex();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();

            Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);
            Pointer<SideData<NDIM, double> > F_data = patch->getPatchData(F_idx);

            for (unsigned int axis = 0; axis < NDIM; axis++)
            {
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                {
                    SideIndex<NDIM> si(it(), axis, SideIndex<NDIM>::Lower);

                    const double multiplier_term =
                        2.0 * (*rho_data)(si) /
                        (mask_surface_tension_force_ctx->rho_liquid + mask_surface_tension_force_ctx->rho_gas);
                    (*F_data)(si) *= multiplier_term;
                }
            }
        }
    }
    return;
}

struct ComputeVariableSurfaceTensionCoefCtx
{
    Pointer<CellVariable<NDIM, double> > T_var;
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator;
    Pointer<INSVCStaggeredHierarchyIntegrator> ins_hier_integrator;
    double sigma0;
    double dsigma_dT0;
    double T_ref;
};

void
compute_surface_tension_coef_function(int F_idx,
                                      Pointer<Patch<NDIM> > patch,
                                      int /*integrator_step*/,
                                      double time,
                                      double /*current_time*/,
                                      double /*new_time*/,
                                      void* ctx)
{
    ComputeVariableSurfaceTensionCoefCtx* compute_variable_surface_tension_coef_ctx =
        static_cast<ComputeVariableSurfaceTensionCoefCtx*>(ctx);

    // parameters
    const double sigma_0 = compute_variable_surface_tension_coef_ctx->sigma0;
    const double dsigma_dT0 = compute_variable_surface_tension_coef_ctx->dsigma_dT0;
    const double T_ref = compute_variable_surface_tension_coef_ctx->T_ref;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int T_scratch_idx = var_db->mapVariableAndContextToIndex(
        compute_variable_surface_tension_coef_ctx->T_var,
        compute_variable_surface_tension_coef_ctx->adv_diff_hier_integrator->getScratchContext());

    const Box<NDIM>& patch_box = patch->getBox();

    Pointer<CellData<NDIM, double> > T_data = patch->getPatchData(T_scratch_idx);
    Pointer<SideData<NDIM, double> > F_data = patch->getPatchData(F_idx);

    for (unsigned int axis = 0; axis < NDIM; axis++)
    {
        for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
        {
            SideIndex<NDIM> si(it(), axis, SideIndex<NDIM>::Lower);

            const double T_sc = 0.5 * ((*T_data)(si.toCell(0)) + (*T_data)(si.toCell(1)));

            const double sigma = sigma_0 + dsigma_dT0 * (T_sc - T_ref);

            (*F_data)(si) *= sigma;
        }
    }
    return;
}

void
compute_marangoni_coef_function(int F_idx,
                                Pointer<Patch<NDIM> > patch,
                                int /*integrator_step*/,
                                double /*time*/,
                                double /*current_time*/,
                                double /*new_time*/,
                                void* ctx)
{
    ComputeVariableSurfaceTensionCoefCtx* compute_variable_surface_tension_coef_ctx =
        static_cast<ComputeVariableSurfaceTensionCoefCtx*>(ctx);

    // parameters
    const double dsigma_dT0 = compute_variable_surface_tension_coef_ctx->dsigma_dT0;

    const Box<NDIM>& patch_box = patch->getBox();

    Pointer<SideData<NDIM, double> > F_data = patch->getPatchData(F_idx);

    for (unsigned int axis = 0; axis < NDIM; axis++)
    {
        for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
        {
            SideIndex<NDIM> si(it(), axis, SideIndex<NDIM>::Lower);

            const double marangoni_coef = dsigma_dT0;
            (*F_data)(si) *= marangoni_coef; // marangoni_coef is constant for this example. So it can be set
                                             // through input file as well instead of
                                             // using this callback function.
        }
    }
    return;
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

        Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
            "AdvDiffSemiImplicitHierarchyIntegrator",
            app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
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
        CircularInterface circle;
        circle.R = input_db->getDouble("R");
        circle.X0[0] = input_db->getDouble("XCOM");
        circle.X0[1] = input_db->getDouble("YCOM");
#if (NDIM == 3)
        circle.X0[2] = input_db->getDouble("ZCOM");
#endif

        // register level set
        Pointer<CellVariable<NDIM, double> > ls_var = new CellVariable<NDIM, double>("ls_var");
        adv_diff_integrator->registerTransportedQuantity(ls_var, true);
        adv_diff_integrator->setDiffusionCoefficient(ls_var, 0.0);

        Pointer<RelaxationLSMethod> level_set_ops =
            new RelaxationLSMethod("RelaxationLSMethod", app_initializer->getComponentDatabase("RelaxationLSMethod"));
        LSLocateCircularInterface* ptr_LSLocateCircularInterface =
            new LSLocateCircularInterface("LSLocateCircularInterface", adv_diff_integrator, ls_var, &circle);
        level_set_ops->registerInterfaceNeighborhoodLocatingFcn(&callLSLocateCircularInterfaceCallbackFunction,
                                                                static_cast<void*>(ptr_LSLocateCircularInterface));
        SetLSProperties* ptr_SetLSProperties = new SetLSProperties("SetLSProperties", NULL, level_set_ops);
        adv_diff_integrator->registerResetFunction(
            ls_var, &callSetGasLSCallbackFunction, static_cast<void*>(ptr_SetLSProperties));

        // register temperature
        Pointer<CellVariable<NDIM, double> > T_var = new CellVariable<NDIM, double>("Temperature");
        adv_diff_integrator->registerTransportedQuantity(T_var, true);

        // set Advection velocity.
        adv_diff_integrator->setAdvectionVelocity(ls_var, time_integrator->getAdvectionVelocityVariable());

        // set priority.
        adv_diff_integrator->setResetPriority(ls_var, 0);
        adv_diff_integrator->setResetPriority(T_var, 1);

        // set initial conditions for the variables.
        if (input_db->keyExists("LevelSetInitialConditions"))
        {
            Pointer<CartGridFunction> ls_init = new muParserCartGridFunction(
                "ls_init", app_initializer->getComponentDatabase("LevelSetInitialConditions"), grid_geometry);
            adv_diff_integrator->setInitialConditions(ls_var, ls_init);
        }

        // set initial conditions for the variables.
        if (input_db->keyExists("TemperatureInitialConditions"))
        {
            Pointer<CartGridFunction> T_init = new muParserCartGridFunction(
                "T_init", app_initializer->getComponentDatabase("TemperatureInitialConditions"), grid_geometry);
            adv_diff_integrator->setInitialConditions(T_var, T_init);
        }

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

        // Setup the INS maintained material properties.
        Pointer<SideVariable<NDIM, double> > rho_sc_var = new SideVariable<NDIM, double>("rho_sc_var");
        time_integrator->registerMassDensityVariable(rho_sc_var);

        Pointer<CellVariable<NDIM, double> > mu_var = new CellVariable<NDIM, double>("mu");
        time_integrator->registerViscosityVariable(mu_var);

        // Create Eulerian boundary condition specification objects (when
        // necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();

        RobinBcCoefStrategy<NDIM>* T_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("TemperatureBcCoefs"))
        {
            T_bc_coef = new muParserRobinBcCoefs(
                "T_bc_coef", app_initializer->getComponentDatabase("TemperatureBcCoefs"), grid_geometry);
            adv_diff_integrator->setPhysicalBcCoef(T_var, T_bc_coef);
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

        // thermophysical properties and parameters.
        const double rho_liquid = input_db->getDouble("RHO_L");
        const double rho_gas = input_db->getDouble("RHO_G");
        const double mu_liquid = input_db->getDouble("MU_L");
        const double mu_gas = input_db->getDouble("MU_G");
        const double sigma_0 = input_db->getDouble("SIGMA_0");
        const double dsigma_dT_0 = input_db->getDouble("DSIGMA_DT_0");
        const double ref_temperature = input_db->getDouble("REFERENCE_TEMPERATURE");
        const double temperature_gradient = input_db->getDouble("TEMPERATURE_GRADIENT");
        const double num_interface_cells = input_db->getDouble("NUMBER_OF_INTERFACE_CELLS");

        // Callback functions can either be registered with the NS integrator, or
        // the advection-diffusion integrator
        SetFluidProperties* ptr_SetFluidProperties = new SetFluidProperties("SetFluidProperties",
                                                                            adv_diff_integrator,
                                                                            ls_var,
                                                                            rho_liquid,
                                                                            rho_gas,
                                                                            mu_liquid,
                                                                            mu_gas,
                                                                            1 /*ls_reinit_interval*/,
                                                                            num_interface_cells);

        time_integrator->registerResetFluidDensityFcn(&callSetFluidDensityCallbackFunction,
                                                      static_cast<void*>(ptr_SetFluidProperties));
        time_integrator->registerResetFluidViscosityFcn(&callSetFluidViscosityCallbackFunction,
                                                        static_cast<void*>(ptr_SetFluidProperties));

        // Register surface tension force.
        Pointer<SurfaceTensionForceFunction> surface_tension_force = new MarangoniSurfaceTensionForceFunction(
            "MarangoniSurfaceTensionForceFunction",
            app_initializer->getComponentDatabase("MarangoniSurfaceTensionForceFunction"),
            adv_diff_integrator,
            ls_var,
            T_var,
            T_bc_coef);

        // Register callback function to multiply the surface tension term with the
        // coefficient.
        MaskSurfaceTensionForceCtx mask_surface_tension_force_ctx;
        mask_surface_tension_force_ctx.ins_hier_integrator = time_integrator;
        mask_surface_tension_force_ctx.rho_liquid = rho_liquid;
        mask_surface_tension_force_ctx.rho_gas = rho_gas;
        mask_surface_tension_force_ctx.adv_diff_hier_integrator = adv_diff_integrator;

        surface_tension_force->registerSurfaceTensionForceMasking(&mask_surface_tension_force,
                                                                  static_cast<void*>(&mask_surface_tension_force_ctx));

        // Register variable coefficient surface tension.
        ComputeVariableSurfaceTensionCoefCtx compute_variable_surface_tension_coef_ctx;
        compute_variable_surface_tension_coef_ctx.T_var = T_var;
        compute_variable_surface_tension_coef_ctx.adv_diff_hier_integrator = adv_diff_integrator;
        compute_variable_surface_tension_coef_ctx.ins_hier_integrator = time_integrator;
        compute_variable_surface_tension_coef_ctx.sigma0 = sigma_0;
        compute_variable_surface_tension_coef_ctx.dsigma_dT0 = dsigma_dT_0;
        compute_variable_surface_tension_coef_ctx.T_ref = ref_temperature;

        surface_tension_force->registerSurfaceTensionCoefficientFunction(
            &compute_surface_tension_coef_function, static_cast<void*>(&compute_variable_surface_tension_coef_ctx));

        // Register variable marangoni coefficient dsigma_dT.
        Pointer<MarangoniSurfaceTensionForceFunction> marangoni_force = surface_tension_force;
        marangoni_force->registerMarangoniCoefficientFunction(
            &compute_marangoni_coef_function, static_cast<void*>(&compute_variable_surface_tension_coef_ctx));

        Pointer<CartGridFunctionSet> eul_forces = new CartGridFunctionSet("eulerian_forces");
        eul_forces->addFunction(surface_tension_force);
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

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int phi_idx = var_db->mapVariableAndContextToIndex(ls_var, adv_diff_integrator->getCurrentContext());
        const int H_cloned_idx = var_db->registerClonedPatchDataIndex(ls_var, phi_idx);

        // Interpolating side-centered velocity to cell-centered
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > U_sc_var =
            time_integrator->getVelocityVariable();
        const int U_sc_idx = var_db->mapVariableAndContextToIndex(U_sc_var, time_integrator->getCurrentContext());
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > U_cc_var;
        U_cc_var = new CellVariable<NDIM, double>("U_cc", NDIM);
        int U_cc_idx = var_db->registerVariableAndContext(U_cc_var, var_db->getContext("U_cc"), 0);

        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > v_var;
        v_var = new CellVariable<NDIM, double>("v_cc");
        int v_idx = var_db->registerVariableAndContext(v_var, var_db->getContext("v_cc"), 0);

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(H_cloned_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(U_cc_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(v_idx, loop_time);
        }

        // File to write to for fluid mass data
        ofstream output_file("output");

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

            HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy, coarsest_ln, finest_ln);

            const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
            HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
            HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);

            hier_math_ops.interp(U_cc_idx, U_cc_var, U_sc_idx, U_sc_var, nullptr, loop_time, true);

            const double U_ref =
                -(2.0 * dsigma_dT_0 * temperature_gradient * circle.R) / (6.0 * mu_liquid + 9.0 * mu_gas);
            const double t_ref = 1.0;

            // Calculate Heaviside function and compute the rise velocity.
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
                    const double alpha = num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

                    Pointer<CellData<NDIM, double> > U_cc_data = patch->getPatchData(U_cc_idx);
                    Pointer<CellData<NDIM, double> > v_data = patch->getPatchData(v_idx);
                    Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(phi_idx);
                    Pointer<CellData<NDIM, double> > H_cloned_data = patch->getPatchData(H_cloned_idx);

                    for (Box<NDIM>::Iterator it(patch_box); it; it++)
                    {
                        CellIndex<NDIM> ci(it());
                        (*v_data)(ci) = (*U_cc_data)(ci, 1) / U_ref; // non_dimensional velocity
                        const double phi = (*phi_data)(ci);
                        (*H_cloned_data)(ci) = IBTK::smooth_heaviside(phi, alpha); // This gives H=1 in the liquid.
                        (*H_cloned_data)(ci) = 1.0 - (*H_cloned_data)(ci);
                    }
                }
            }

            double vol = hier_cc_data_ops.integral(H_cloned_idx, wgt_cc_idx);

            // Get volume weights in the region
            hier_cc_data_ops.multiply(H_cloned_idx, H_cloned_idx, wgt_cc_idx);
            double v_integral = hier_cc_data_ops.integral(v_idx, H_cloned_idx);
            pout << "vol is\t" << vol << "\n";
            pout << "Rise velocity is\t" << v_integral / vol << "\t"
                 << "at time \t" << loop_time << "\t" << loop_time / t_ref << "\n";
            // Write the rise velocity and time to a file
            if (IBTK_MPI::getRank() == 0)
            {
                output_file << std::setprecision(13) << loop_time / t_ref << "\t" << v_integral / vol << std::endl;
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

        output_file.close();

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(H_cloned_idx);
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(U_cc_idx);
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(v_idx);
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

        // Cleanup other dumb pointers
        delete ptr_SetLSProperties;
        delete ptr_SetFluidProperties;
        delete ptr_LSLocateCircularInterface;

    } // cleanup dynamically allocated objects prior to shutdown
} // main
