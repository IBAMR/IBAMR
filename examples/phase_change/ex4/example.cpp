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
#include <ibamr/LevelSetUtilities.h>
#include <ibamr/PhaseChangeDivUSourceFunction.h>
#include <ibamr/PhaseChangeUtilities.h>
#include <ibamr/RelaxationLSMethod.h>
#include <ibamr/SurfaceTensionForceFunction.h>
#include <ibamr/vc_ins_utilities.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <algorithm>

#include "vc_ins_vof_utilities.h"

#include <ibamr/app_namespaces.h>

// Application
#include "LSLocateInterface.h"
#include "LevelSetInitialCondition.h"
#include "LiquidFractionInitialCondition.h"
#include "TemperatureInitialCondition.h"

struct MaskSurfaceTensionForceCtx
{
    Pointer<CellVariable<NDIM, double>> lf_var;
    RobinBcCoefStrategy<NDIM>* lf_bc_coef;
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator;
    Pointer<INSVCStaggeredHierarchyIntegrator> ins_hier_integrator;
    double rho_liquid;
    double rho_gas;
};

void
mask_surface_tension_force(int F_idx,
                           Pointer<HierarchyMathOps> hier_math_ops,
                           int /*integrator_step*/,
                           double time,
                           double /*current_time*/,
                           double /*new_time*/,
                           void* ctx)
{
    MaskSurfaceTensionForceCtx* mask_surface_tension_force_ctx = static_cast<MaskSurfaceTensionForceCtx*>(ctx);
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int lf_new_idx =
        var_db->mapVariableAndContextToIndex(mask_surface_tension_force_ctx->lf_var,
                                             mask_surface_tension_force_ctx->adv_diff_hier_integrator->getNewContext());
    const int lf_scratch_idx = var_db->mapVariableAndContextToIndex(
        mask_surface_tension_force_ctx->lf_var,
        mask_surface_tension_force_ctx->adv_diff_hier_integrator->getScratchContext());

    // ghost cell filling for liquid fraction variable.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> lf_transaction_comps(1);
    lf_transaction_comps[0] = InterpolationTransactionComponent(lf_scratch_idx,
                                                                lf_new_idx,
                                                                "CONSERVATIVE_LINEAR_REFINE",
                                                                false,
                                                                "CONSERVATIVE_COARSEN",
                                                                "LINEAR",
                                                                false,
                                                                mask_surface_tension_force_ctx->lf_bc_coef);

    Pointer<HierarchyGhostCellInterpolation> lf_hier_bdry_fill = new HierarchyGhostCellInterpolation();
    lf_hier_bdry_fill->initializeOperatorState(lf_transaction_comps, patch_hierarchy);
    lf_hier_bdry_fill->fillData(time);

    int rho_idx = mask_surface_tension_force_ctx->ins_hier_integrator->getLinearOperatorRhoPatchDataIndex();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();

            Pointer<CellData<NDIM, double>> lf_data = patch->getPatchData(lf_scratch_idx);
            Pointer<SideData<NDIM, double>> rho_data = patch->getPatchData(rho_idx);
            Pointer<SideData<NDIM, double>> F_data = patch->getPatchData(F_idx);

            for (unsigned int axis = 0; axis < NDIM; axis++)
            {
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                {
                    SideIndex<NDIM> si(it(), axis, SideIndex<NDIM>::Lower);

                    const double lf_sc_value = 0.5 * ((*lf_data)(si.toCell(0)) + (*lf_data)(si.toCell(1)));

                    const double multiplier_term =
                        2.0 * (*rho_data)(si) /
                        (mask_surface_tension_force_ctx->rho_liquid + mask_surface_tension_force_ctx->rho_gas) *
                        lf_sc_value;
                    (*F_data)(si) *= multiplier_term;
                }
            }
        }
    }
}

struct SynchronizePCMVOFWithLSCtx
{
    IBAMR::VCINSVOFUtilities::VOFFromLevelSetInitializer* vof_from_ls = nullptr;
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator;
    Pointer<CellVariable<NDIM, double>> pcm_vof_var;
    Pointer<CellVariable<NDIM, double>> liquid_fraction_var;
};

void
clamp_liquid_fraction_to_pcm_vof(Pointer<PatchHierarchy<NDIM>> hierarchy,
                                 const int pcm_vof_idx,
                                 const int liquid_fraction_idx)
{
    const int finest_ln = hierarchy->getFinestLevelNumber();
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM, double>> pcm_vof_data = patch->getPatchData(pcm_vof_idx);
            Pointer<CellData<NDIM, double>> liquid_fraction_data = patch->getPatchData(liquid_fraction_idx);

#if !defined(NDEBUG)
            TBOX_ASSERT(!pcm_vof_data.isNull());
            TBOX_ASSERT(!liquid_fraction_data.isNull());
#endif

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                const CellIndex<NDIM> ci(it());
                const double C = std::max(0.0, std::min(1.0, (*pcm_vof_data)(ci)));
                const double L = std::max(0.0, std::min(C, (*liquid_fraction_data)(ci)));
                (*pcm_vof_data)(ci) = C;
                (*liquid_fraction_data)(ci) = L;
            }
        }
    }
}

void
synchronize_pcm_vof_with_level_set(int pcm_vof_current_idx,
                                   Pointer<HierarchyMathOps> hier_math_ops,
                                   int /*integrator_step*/,
                                   double time,
                                   bool /*initial_time*/,
                                   bool /*regrid_time*/,
                                   void* ctx)
{
    auto* sync_ctx = static_cast<SynchronizePCMVOFWithLSCtx*>(ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(sync_ctx);
    TBOX_ASSERT(sync_ctx->vof_from_ls);
    TBOX_ASSERT(!sync_ctx->adv_diff_integrator.isNull());
    TBOX_ASSERT(!sync_ctx->pcm_vof_var.isNull());
    TBOX_ASSERT(!sync_ctx->liquid_fraction_var.isNull());
    TBOX_ASSERT(hier_math_ops);
#endif

    // Reset priority guarantees that ls_var has already been reset/reinitialized.
    sync_ctx->vof_from_ls->computeVOFFromLevelSet(time, /*use_new_context=*/false);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int mapped_pcm_vof_current_idx =
        var_db->mapVariableAndContextToIndex(sync_ctx->pcm_vof_var, sync_ctx->adv_diff_integrator->getCurrentContext());
    const int liquid_fraction_current_idx = var_db->mapVariableAndContextToIndex(
        sync_ctx->liquid_fraction_var, sync_ctx->adv_diff_integrator->getCurrentContext());

#if !defined(NDEBUG)
    TBOX_ASSERT(mapped_pcm_vof_current_idx == pcm_vof_current_idx);
#endif

    // liquid_fraction is a whole-cell liquid PCM fraction: enforce 0 <= L <= C.
    clamp_liquid_fraction_to_pcm_vof(
        hier_math_ops->getPatchHierarchy(), pcm_vof_current_idx, liquid_fraction_current_idx);
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

        // register level set
        Pointer<CellVariable<NDIM, double>> ls_var = new CellVariable<NDIM, double>("ls_var");
        adv_diff_integrator->registerTransportedQuantity(ls_var, true);
        adv_diff_integrator->setDiffusionCoefficient(ls_var, 0.0);

        const double initial_gas_pcm_interface_position = input_db->getDouble("INITIAL_GAS_PCM_INTERFACE_POSITION");
        Pointer<RelaxationLSMethod> level_set_ops =
            new RelaxationLSMethod("RelaxationLSMethod", app_initializer->getComponentDatabase("RelaxationLSMethod"));
        LSLocateInterface* ptr_LSLocateInterface =
            new LSLocateInterface("LSLocateInterface", adv_diff_integrator, ls_var, initial_gas_pcm_interface_position);
        level_set_ops->registerInterfaceNeighborhoodLocatingFcn(&callLSLocateInterfaceCallbackFunction,
                                                                static_cast<void*>(ptr_LSLocateInterface));
        IBAMR::LevelSetUtilities::SetLSProperties setSetLSProperties("SetLSProperties", level_set_ops);
        adv_diff_integrator->registerResetFunction(
            ls_var, &IBAMR::LevelSetUtilities::setLSDataPatchHierarchy, static_cast<void*>(&setSetLSProperties));

        // register liquid fraction
        Pointer<CellVariable<NDIM, double>> lf_var = new CellVariable<NDIM, double>("lf_var");
        Pointer<EnthalpyHierarchyIntegrator> enthalpy_hier_integrator = adv_diff_integrator;

        enthalpy_hier_integrator->registerLevelSetVariable(ls_var);

        enthalpy_hier_integrator->registerLiquidFractionVariable(lf_var, true);

        // register specific enthalpy
        Pointer<CellVariable<NDIM, double>> h_var = new CellVariable<NDIM, double>("h_var");
        enthalpy_hier_integrator->registerSpecificEnthalpyVariable(h_var, true);

        Pointer<CellVariable<NDIM, double>> pcm_vof_var = new CellVariable<NDIM, double>("pcm_vof_var");
        adv_diff_integrator->registerTransportedQuantity(pcm_vof_var, true);
        adv_diff_integrator->setDiffusionCoefficient(pcm_vof_var, 0.0);

        enthalpy_hier_integrator->registerHeavisideVariable(pcm_vof_var);

        // register temperature
        Pointer<CellVariable<NDIM, double>> T_var = new CellVariable<NDIM, double>("Temperature");
        enthalpy_hier_integrator->registerTemperatureVariable(T_var, true);

        // set Advection velocity.
        adv_diff_integrator->setAdvectionVelocity(ls_var, time_integrator->getAdvectionVelocityVariable());
        enthalpy_hier_integrator->setAdvectionVelocity(time_integrator->getAdvectionVelocityVariable());

        const ConvectiveDifferencingType ls_difference_form =
            IBAMR::string_to_enum<ConvectiveDifferencingType>(input_db->getString("LS_CONVECTIVE_FORM"));
        adv_diff_integrator->setConvectiveDifferencingType(ls_var, ls_difference_form);

        // adv_diff_integrator->setAdvectionVelocity( pcm_vof_var, time_integrator->getAdvectionVelocityVariable());

        // const ConvectiveDifferencingType pcm_vof_difference_form =
        //     IBAMR::string_to_enum<ConvectiveDifferencingType>(input_db->getString("H_CONVECTIVE_FORM"));
        // adv_diff_integrator->setConvectiveDifferencingType(pcm_vof_var, pcm_vof_difference_form);

        adv_diff_integrator->setResetPriority(ls_var, 0);
        adv_diff_integrator->setResetPriority(pcm_vof_var, 1);

        // set initial conditions for the variables.
        Pointer<CartGridFunction> ls_init = new LevelSetInitialCondition("ls_init", initial_gas_pcm_interface_position);
        adv_diff_integrator->setInitialConditions(ls_var, ls_init);

        Pointer<CartGridFunction> pcm_vof_init =
            new IBAMR::VCINSVOFUtilities::VOFInitialConditionFromLevelSet("pcm_vof_init", ls_init);
        adv_diff_integrator->setInitialConditions(pcm_vof_var, pcm_vof_init);

        IBAMR::VCINSVOFUtilities::VOFFromLevelSetInitializer pcm_vof_from_ls(
            "pcm_vof_from_ls", adv_diff_integrator, ls_var, pcm_vof_var);


        pcm_vof_from_ls.registerIntegrateHierarchyCallback();

        SynchronizePCMVOFWithLSCtx pcm_vof_sync_ctx;
        pcm_vof_sync_ctx.vof_from_ls = &pcm_vof_from_ls;
        pcm_vof_sync_ctx.adv_diff_integrator = adv_diff_integrator;
        pcm_vof_sync_ctx.pcm_vof_var = pcm_vof_var;
        pcm_vof_sync_ctx.liquid_fraction_var = lf_var;

        adv_diff_integrator->registerResetFunction(
            pcm_vof_var, &synchronize_pcm_vof_with_level_set, static_cast<void*>(&pcm_vof_sync_ctx));

        const double initial_liquid_solid_interface_position =
            input_db->getDouble("INITIAL_LIQUID_SOLID_INTERFACE_POSITION");
        const double initial_liquid_temperature = input_db->getDouble("INITIAL_LIQUID_TEMPERATURE");
        const double initial_solid_temperature = input_db->getDouble("INITIAL_SOLID_TEMPERATURE");

        Pointer<CartGridFunction> T_init = new TemperatureInitialCondition(
            "T_init", initial_liquid_solid_interface_position, initial_liquid_temperature, initial_solid_temperature);
        enthalpy_hier_integrator->setTemperatureInitialCondition(T_var, T_init);

        Pointer<CartGridFunction> lf_init =
            new LiquidFractionInitialCondition("lf_init", initial_liquid_solid_interface_position);
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

        // Setup the INS maintained material properties.
        Pointer<SideVariable<NDIM, double>> rho_sc_var = new SideVariable<NDIM, double>("rho_sc_var");
        time_integrator->registerMassDensityVariable(rho_sc_var);

        Pointer<CellVariable<NDIM, double>> mu_var = new CellVariable<NDIM, double>("mu");
        time_integrator->registerViscosityVariable(mu_var);

        Pointer<CellVariable<NDIM, double>> rho_cc_var = new CellVariable<NDIM, double>("rho_cc_var");
        enthalpy_hier_integrator->registerDensityVariable(rho_cc_var, true);

        Pointer<CellVariable<NDIM, double>> Cp_var = new CellVariable<NDIM, double>("Cp");
        enthalpy_hier_integrator->registerSpecificHeatVariable(Cp_var, true);

        // Create Eulerian boundary condition specification objects (when
        // necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();

        std::unique_ptr<RobinBcCoefStrategy<NDIM>> pcm_vof_bc_coef;
        if (!(periodic_shift.min() > 0))
        {
            pcm_vof_bc_coef = std::make_unique<muParserRobinBcCoefs>(
                "pcm_vof_bc_coef", app_initializer->getComponentDatabase("HeavisideBcCoefs"), grid_geometry);
            adv_diff_integrator->setPhysicalBcCoef(pcm_vof_var, pcm_vof_bc_coef.get());
        }

        std::unique_ptr<RobinBcCoefStrategy<NDIM>> T_bc_coef;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("TemperatureBcCoefs"))
        {
            T_bc_coef = std::make_unique<muParserRobinBcCoefs>(
                "T_bc_coef", app_initializer->getComponentDatabase("TemperatureBcCoefs"), grid_geometry);
            enthalpy_hier_integrator->setTemperaturePhysicalBcCoef(T_var, T_bc_coef.get());
        }

        std::unique_ptr<RobinBcCoefStrategy<NDIM>> h_bc_coef;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("EnthalpyBcCoefs"))
        {
            h_bc_coef = std::make_unique<muParserRobinBcCoefs>(
                "h_bc_coef", app_initializer->getComponentDatabase("EnthalpyBcCoefs"), grid_geometry);
            enthalpy_hier_integrator->setEnthalpyBcCoef(h_bc_coef.get());
        }

        std::unique_ptr<RobinBcCoefStrategy<NDIM>> lf_bc_coef;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("LiquidFractionBcCoefs"))
        {
            lf_bc_coef = std::make_unique<muParserRobinBcCoefs>(
                "lf_bc_coef", app_initializer->getComponentDatabase("LiquidFractionBcCoefs"), grid_geometry);
        }

        vector<std::unique_ptr<RobinBcCoefStrategy<NDIM>>> u_bc_coefs(NDIM);
        if (periodic_shift.min() == 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();

                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                u_bc_coefs[d] = std::make_unique<muParserRobinBcCoefs>(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            time_integrator->registerPhysicalBoundaryConditions({
                u_bc_coefs[0].get(), u_bc_coefs[1].get()
#if (NDIM == 3)
                                         ,
                    u_bc_coefs[2].get()
#endif
            });
        }

        std::unique_ptr<RobinBcCoefStrategy<NDIM>> rho_bc_coef;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("DensityBcCoefs"))
        {
            rho_bc_coef = std::make_unique<muParserRobinBcCoefs>(
                "rho_bc_coef", app_initializer->getComponentDatabase("DensityBcCoefs"), grid_geometry);
            time_integrator->registerMassDensityBoundaryConditions(rho_bc_coef.get());
            enthalpy_hier_integrator->registerMassDensityBoundaryConditions(rho_bc_coef.get());
        }

        std::unique_ptr<RobinBcCoefStrategy<NDIM>> mu_bc_coef;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("ViscosityBcCoefs"))
        {
            mu_bc_coef = std::make_unique<muParserRobinBcCoefs>(
                "mu_bc_coef", app_initializer->getComponentDatabase("ViscosityBcCoefs"), grid_geometry);
            time_integrator->registerViscosityBoundaryConditions(mu_bc_coef.get());
        }

        std::unique_ptr<RobinBcCoefStrategy<NDIM>> k_bc_coef;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("ThermalConductivityBcCoefs"))
        {
            k_bc_coef = std::make_unique<muParserRobinBcCoefs>(
                "k_bc_coef", app_initializer->getComponentDatabase("ThermalConductivityBcCoefs"), grid_geometry);
            enthalpy_hier_integrator->registerThermalConductivityBoundaryConditions(k_bc_coef.get());
        }

        std::unique_ptr<RobinBcCoefStrategy<NDIM>> ls_bc_coef;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("LevelSetBcCoefs"))
        {
            ls_bc_coef = std::make_unique<muParserRobinBcCoefs>(
                "ls_bc_coef", app_initializer->getComponentDatabase("LevelSetBcCoefs"), grid_geometry);
            adv_diff_integrator->setPhysicalBcCoef(ls_var, ls_bc_coef.get());
            level_set_ops->registerPhysicalBoundaryCondition(ls_bc_coef.get());
        }

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
        IBAMR::PhaseChangeUtilities::SetFluidProperties setSetFluidProperties("SetFluidProperties",
                                                                              adv_diff_integrator,
                                                                              pcm_vof_var,
                                                                              pcm_vof_bc_coef.get(),
                                                                              lf_var,
                                                                              lf_bc_coef.get(),
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

        time_integrator->registerResetFluidDensityFcn(&IBAMR::PhaseChangeUtilities::callSetDensityCallbackFunction,
                                                      static_cast<void*>(&setSetFluidProperties));
        time_integrator->registerResetFluidViscosityFcn(&IBAMR::PhaseChangeUtilities::callSetViscosityCallbackFunction,
                                                        static_cast<void*>(&setSetFluidProperties));

        enthalpy_hier_integrator->registerResetDiffusionCoefficientFcn(
            &IBAMR::PhaseChangeUtilities::callSetThermalConductivityCallbackFunction,
            static_cast<void*>(&setSetFluidProperties));

        enthalpy_hier_integrator->registerResetSpecificHeatFcn(
            &IBAMR::PhaseChangeUtilities::callSetSpecificHeatCallbackFunction,
            static_cast<void*>(&setSetFluidProperties));

        enthalpy_hier_integrator->registerResetDensityFcn(&IBAMR::PhaseChangeUtilities::callSetDensityCallbackFunction,
                                                          static_cast<void*>(&setSetFluidProperties));

        // Pointer<CellVariable<NDIM, double>> pcm_vof_F_var =
        //     new CellVariable<NDIM, double>(pcm_vof_var->getName() + "_F");
        // adv_diff_integrator->registerSourceTerm(pcm_vof_F_var, true);

        // Pointer<CartGridFunction> pcm_vof_forcing_fcn = new HeavisideForcingFunction(
        //     "pcm_vof_forcing_fcn",
        //     adv_diff_integrator,
        //     pcm_vof_var,
        //     time_integrator->getAdvectionVelocityVariable());
        // adv_diff_integrator->setSourceTermFunction(pcm_vof_F_var, pcm_vof_forcing_fcn);
        // adv_diff_integrator->setSourceTerm(pcm_vof_var, pcm_vof_F_var);

        // Register source term for Div U equation.
        Pointer<CartGridFunction> Div_U_forcing_fcn =
            new PhaseChangeDivUSourceFunction("Div_U_forcing_fcn", enthalpy_hier_integrator);
        time_integrator->registerVelocityDivergenceFunction(Div_U_forcing_fcn);

        // Register surface tension force.
        Pointer<SurfaceTensionForceFunction> surface_tension_force =
            new SurfaceTensionForceFunction("SurfaceTensionForceFunction",
                                            app_initializer->getComponentDatabase("SurfaceTensionForceFunction"),
                                            adv_diff_integrator,
                                            ls_var);

        // Register callback function to multiply the surface tension term with the coefficient.
        MaskSurfaceTensionForceCtx mask_surface_tension_force_ctx;
        mask_surface_tension_force_ctx.ins_hier_integrator = time_integrator;
        mask_surface_tension_force_ctx.rho_liquid = rho_liquid;
        mask_surface_tension_force_ctx.rho_gas = rho_gas;
        mask_surface_tension_force_ctx.adv_diff_hier_integrator = adv_diff_integrator;
        mask_surface_tension_force_ctx.lf_var = lf_var;
        mask_surface_tension_force_ctx.lf_bc_coef = lf_bc_coef.get();

        surface_tension_force->registerSurfaceTensionForceMasking(&mask_surface_tension_force,
                                                                  static_cast<void*>(&mask_surface_tension_force_ctx));

        // Register gravity force.
        std::vector<double> grav_const(NDIM);
        input_db->getDoubleArray("GRAV_CONST", &grav_const[0], NDIM);
        Pointer<CartGridFunction> grav_force =
            new IBAMR::VCINSUtilities::GravityForcing("GravityForcing", time_integrator, grav_const);

        Pointer<CartGridFunctionSet> eul_forces = new CartGridFunctionSet("eulerian_forces");
        eul_forces->addFunction(grav_force);
        eul_forces->addFunction(surface_tension_force);
        time_integrator->registerBodyForceFunction(eul_forces);

        // Configure the drag force object to enforce solid velocity to be zero.
        Pointer<CarmanKozenyDragForce> drag_force =
            new CarmanKozenyDragForce("drag_force",
                                      pcm_vof_var,
                                      lf_var,
                                      adv_diff_integrator,
                                      time_integrator,
                                      app_initializer->getComponentDatabase("CarmanKozenyDragForce"),
                                      /*register_for_restart*/ true);
        time_integrator->registerBrinkmanPenalizationStrategy(drag_force);

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        pcm_vof_from_ls.computeVOFFromLevelSet(time_integrator->getIntegratorTime(), /*use_new_context=*/false);

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int pcm_vof_current_idx =
            var_db->mapVariableAndContextToIndex(pcm_vof_var, adv_diff_integrator->getCurrentContext());
        const int lf_current_idx =
            var_db->mapVariableAndContextToIndex(lf_var, adv_diff_integrator->getCurrentContext());

        clamp_liquid_fraction_to_pcm_vof(patch_hierarchy, pcm_vof_current_idx, lf_current_idx);

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
        // VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int rho_idx = var_db->mapVariableAndContextToIndex(rho_cc_var, adv_diff_integrator->getCurrentContext());
        const int pcm_vof_idx =
            var_db->mapVariableAndContextToIndex(pcm_vof_var, adv_diff_integrator->getCurrentContext());
        const int pcm_mass_idx = var_db->registerClonedPatchDataIndex(pcm_vof_var, pcm_vof_idx);

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(pcm_mass_idx, loop_time);
        }

        Pointer<HierarchyCellDataOpsReal<NDIM, double>> hier_cc_data_ops =
            new HierarchyCellDataOpsReal<NDIM, double>(patch_hierarchy, coarsest_ln, finest_ln);

        std::ofstream pcm_mass_file;
        if (SAMRAI_MPI::getRank() == 0)
        {
            std::string FILE_NAME = input_db->getString("MASS_CONSERVATION_FILE_NAME");
            pcm_mass_file.open(FILE_NAME, ios_base::out | ios_base::app);
            pcm_mass_file.precision(16);
            pcm_mass_file.setf(ios::fixed, ios::floatfield);
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

            hier_cc_data_ops->multiply(pcm_mass_idx, rho_idx, pcm_vof_idx);
            HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy, coarsest_ln, finest_ln);
            const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
            const double mass = hier_cc_data_ops->integral(pcm_mass_idx, wgt_cc_idx);

            if (SAMRAI_MPI::getRank() == 0)
            {
                pcm_mass_file.precision(16);
                pcm_mass_file.setf(ios::fixed, ios::floatfield);
                pcm_mass_file << loop_time << "\t" << mass << "\n";
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

        // Cleanup pointers.
        delete ptr_LSLocateInterface;

        // Close the logging streams.
        if (SAMRAI_MPI::getRank() == 0)
        {
            pcm_mass_file.close();
        }

    } // cleanup dynamically allocated objects prior to shutdown
} // main
