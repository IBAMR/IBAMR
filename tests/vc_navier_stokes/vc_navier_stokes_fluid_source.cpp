// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
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
#include <CellVariable.h>
#include <LoadBalancer.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Applications
#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>
#include <string>

#include "SetFluidProperties.cpp"

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        const std::string discretization_form = app_initializer->getComponentDatabase("Main")->getStringWithDefault(
            "discretization_form", "NON_CONSERVATIVE");
        Pointer<INSVCStaggeredHierarchyIntegrator> time_integrator;
        const bool is_conservative = (discretization_form == "CONSERVATIVE");
        if (discretization_form == "NON_CONSERVATIVE")
        {
            time_integrator = new INSVCStaggeredNonConservativeHierarchyIntegrator(
                "INSVCStaggeredNonConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredNonConservativeHierarchyIntegrator"));
        }
        else if (is_conservative)
        {
            time_integrator = new INSVCStaggeredConservativeHierarchyIntegrator(
                "INSVCStaggeredConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredConservativeHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported discretization_form: " << discretization_form << "\n"
                                                           << "Valid options are: NON_CONSERVATIVE, CONSERVATIVE\n");
        }

        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        if (periodic_shift.min() <= 0)
        {
            TBOX_ERROR("vc_navier_stokes_fluid_source:\n"
                       << "  CartesianGeometry requires periodic boundaries in all directions\n");
        }

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

        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        time_integrator->registerVelocityInitialConditions(u_init);
        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        time_integrator->registerPressureInitialConditions(p_init);

        Pointer<Variable<NDIM>> rho_var;
        if (is_conservative)
        {
            rho_var = new SideVariable<NDIM, double>("rho_var");
        }
        else
        {
            rho_var = new CellVariable<NDIM, double>("rho_var");
        }
        Pointer<CellVariable<NDIM, double>> mu_var = new CellVariable<NDIM, double>("mu_var");
        time_integrator->registerMassDensityVariable(rho_var);
        time_integrator->registerViscosityVariable(mu_var);

        Pointer<CartGridFunction> rho_fcn = new muParserCartGridFunction(
            "rho_fcn", app_initializer->getComponentDatabase("DensityFunction"), grid_geometry);
        Pointer<CartGridFunction> mu_fcn = new muParserCartGridFunction(
            "mu_fcn", app_initializer->getComponentDatabase("ViscosityFunction"), grid_geometry);
        std::unique_ptr<SetFluidProperties> fluid_properties =
            std::make_unique<SetFluidProperties>("SetFluidProperties", rho_fcn, mu_fcn);

        time_integrator->registerMassDensityInitialConditions(rho_fcn);
        time_integrator->registerViscosityInitialConditions(mu_fcn);
        if (!is_conservative)
        {
            time_integrator->registerResetFluidDensityFcn(&callSetFluidDensityCallbackFunction, fluid_properties.get());
        }
        time_integrator->registerResetFluidViscosityFcn(&callSetFluidViscosityCallbackFunction, fluid_properties.get());

        Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
            "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
        time_integrator->registerBodyForceFunction(f_fcn);

        Pointer<CartGridFunction> q_fcn = new muParserCartGridFunction(
            "q_fcn", app_initializer->getComponentDatabase("FluidSourceFunction"), grid_geometry);
        time_integrator->registerFluidSourceFunction(q_fcn);

        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        app_initializer.setNull();

        plog << "Input database:\n";
        input_db->printClassData(plog);

        const double loop_time_end = time_integrator->getEndTime();
        double loop_time = time_integrator->getIntegratorTime();
        double dt = 0.0;
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time = time_integrator->getIntegratorTime();
        }

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const Pointer<Variable<NDIM>> u_var = time_integrator->getVelocityVariable();
        const Pointer<VariableContext> u_ctx = time_integrator->getCurrentContext();
        const int u_idx = var_db->mapVariableAndContextToIndex(u_var, u_ctx);
        const Pointer<SideVariable<NDIM, double>> u_sc_var = u_var;
        if (!u_sc_var)
        {
            TBOX_ERROR("vc_navier_stokes_fluid_source:\n"
                       << "  velocity variable is not side-centered\n");
        }
        const int u_cloned_idx = var_db->registerClonedPatchDataIndex(u_var, u_idx);

        const Pointer<Variable<NDIM>> p_var = time_integrator->getPressureVariable();
        const Pointer<VariableContext> p_ctx = time_integrator->getCurrentContext();
        const int p_idx = var_db->mapVariableAndContextToIndex(p_var, p_ctx);
        const int p_cloned_idx = var_db->registerClonedPatchDataIndex(p_var, p_idx);

        Pointer<VariableContext> check_ctx = var_db->getContext("fluid_source_check");
        Pointer<CellVariable<NDIM, double>> div_u_var = new CellVariable<NDIM, double>("div_u");
        Pointer<CellVariable<NDIM, double>> q_exact_var = new CellVariable<NDIM, double>("q_exact");
        const int div_u_idx = var_db->registerVariableAndContext(div_u_var, check_ctx, IntVector<NDIM>(0));
        const int q_exact_idx = var_db->registerVariableAndContext(q_exact_var, check_ctx, IntVector<NDIM>(0));

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_cloned_idx, loop_time);
            level->allocatePatchData(p_cloned_idx, loop_time);
            level->allocatePatchData(div_u_idx, loop_time);
            level->allocatePatchData(q_exact_idx, loop_time);
        }

        HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(coarsest_ln, finest_ln);
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);

        // Check discrete div u = Q at new_time (loop_time) before subtracting u_init from u_idx.
        hier_math_ops.div(div_u_idx, div_u_var, 1.0, u_idx, u_sc_var, nullptr, loop_time, /*src1_cf_bdry_synch*/ false);
        q_fcn->setDataOnPatchHierarchy(q_exact_idx, q_exact_var, patch_hierarchy, loop_time);
        hier_cc_data_ops.subtract(div_u_idx, div_u_idx, q_exact_idx);
        const double div_l2 = hier_cc_data_ops.L2Norm(div_u_idx, wgt_cc_idx);
        const double div_max = hier_cc_data_ops.maxNorm(div_u_idx, wgt_cc_idx);
        if (!std::isfinite(div_l2) || !std::isfinite(div_max) || div_max > 1.0e-9)
        {
            TBOX_ERROR("vc_navier_stokes_fluid_source:\n"
                       << "  discrete divergence residual exceeds 1.0e-9 at t = " << loop_time << "\n"
                       << "  L2-norm = " << std::setprecision(10) << div_l2 << ", max-norm = " << div_max << "\n");
        }

        u_init->setDataOnPatchHierarchy(u_cloned_idx, u_var, patch_hierarchy, loop_time);
        p_init->setDataOnPatchHierarchy(p_cloned_idx, p_var, patch_hierarchy, loop_time - 0.5 * dt);

        hier_sc_data_ops.subtract(u_idx, u_idx, u_cloned_idx);
        const double u_l1 = hier_sc_data_ops.L1Norm(u_idx, wgt_sc_idx);
        const double u_l2 = hier_sc_data_ops.L2Norm(u_idx, wgt_sc_idx);
        const double u_max = hier_sc_data_ops.maxNorm(u_idx, wgt_sc_idx);

        hier_cc_data_ops.subtract(p_idx, p_idx, p_cloned_idx);
        const double p_l1 = hier_cc_data_ops.L1Norm(p_idx, wgt_cc_idx);
        const double p_l2 = hier_cc_data_ops.L2Norm(p_idx, wgt_cc_idx);
        const double p_max = hier_cc_data_ops.maxNorm(p_idx, wgt_cc_idx);

        if (!std::isfinite(u_l1) || !std::isfinite(u_l2) || !std::isfinite(u_max) || !std::isfinite(p_l1) ||
            !std::isfinite(p_l2) || !std::isfinite(p_max))
        {
            TBOX_ERROR("vc_navier_stokes_fluid_source:\n"
                       << "  computed error norm is nonfinite\n");
        }

        if (IBTK_MPI::getRank() == 0)
        {
            std::ofstream out("output");
            out << "\n"
                << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n"
                << "Computing error norms.\n\n"
                << "Error in u_sc at time " << loop_time << ":\n"
                << "  L1-norm:  " << std::setprecision(10) << u_l1 << "\n"
                << "  L2-norm:  " << u_l2 << "\n"
                << "  max-norm: " << u_max << "\n"
                << "Error in p at time " << loop_time - 0.5 * dt << ":\n"
                << "  L1-norm:  " << p_l1 << "\n"
                << "  L2-norm:  " << p_l2 << "\n"
                << "  max-norm: " << p_max << "\n"
                << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        }

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(u_cloned_idx);
            level->deallocatePatchData(p_cloned_idx);
            level->deallocatePatchData(div_u_idx);
            level->deallocatePatchData(q_exact_idx);
        }
        var_db->removePatchDataIndex(u_cloned_idx);
        var_db->removePatchDataIndex(p_cloned_idx);
        var_db->removePatchDataIndex(div_u_idx);
        var_db->removePatchDataIndex(q_exact_idx);
    }
} // main
