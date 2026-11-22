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

// Headers for basic PETSc objects
#include <petscsys.h>

// Headers for major SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <HierarchySideDataOpsReal.h>
#include <IntVector.h>
#include <LoadBalancer.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RobinBcCoefStrategy.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/INSVCStaggeredConservativeMassMomentumSSPRKIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartSideDoubleCubicCoarsen.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>

// Compare the opt-in constant-density SSPRK3 path against full SSPRK3 integration on a two-level AMR hierarchy.
int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "vc_ssprk3_const_density.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        grid_geometry->addSpatialCoarsenOperator(new CartSideDoubleCubicCoarsen());
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");

        Pointer<SideVariable<NDIM, double>> gcw_guard_var = new SideVariable<NDIM, double>("gcw_guard");
        var_db->registerVariableAndContext(gcw_guard_var, ctx, IntVector<NDIM>(4));

        Pointer<SideVariable<NDIM, double>> v_old_var = new SideVariable<NDIM, double>("v_old");
        Pointer<SideVariable<NDIM, double>> v_cur_var = new SideVariable<NDIM, double>("v_cur");
        Pointer<SideVariable<NDIM, double>> v_new_var = new SideVariable<NDIM, double>("v_new");
        Pointer<SideVariable<NDIM, double>> rho_cur_var = new SideVariable<NDIM, double>("rho_cur");
        Pointer<SideVariable<NDIM, double>> n_fast_var = new SideVariable<NDIM, double>("n_fast");
        Pointer<SideVariable<NDIM, double>> n_ref_var = new SideVariable<NDIM, double>("n_ref");
        Pointer<SideVariable<NDIM, double>> diff_var = new SideVariable<NDIM, double>("diff");

        const int v_old_idx = var_db->registerVariableAndContext(v_old_var, ctx, IntVector<NDIM>(0));
        const int v_cur_idx = var_db->registerVariableAndContext(v_cur_var, ctx, IntVector<NDIM>(0));
        const int v_new_idx = var_db->registerVariableAndContext(v_new_var, ctx, IntVector<NDIM>(0));
        const int rho_cur_idx = var_db->registerVariableAndContext(rho_cur_var, ctx, IntVector<NDIM>(0));
        const int n_fast_idx = var_db->registerVariableAndContext(n_fast_var, ctx, IntVector<NDIM>(0));
        const int n_ref_idx = var_db->registerVariableAndContext(n_ref_var, ctx, IntVector<NDIM>(0));
        const int diff_idx = var_db->registerVariableAndContext(diff_var, ctx, IntVector<NDIM>(0));

        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        int level_number = 0;
        bool done = false;
        while (!done && gridding_algorithm->levelCanBeRefined(level_number))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, 1);
            done = !patch_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }

        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (const int idx : { v_old_idx, v_cur_idx, v_new_idx, rho_cur_idx, n_fast_idx, n_ref_idx, diff_idx })
            {
                level->allocatePatchData(idx, 0.0);
            }
        }

        const double dt = input_db->getDouble("dt");
        muParserCartGridFunction rho_init("rho_init", app_initializer->getComponentDatabase("rho_init"), grid_geometry);
        muParserCartGridFunction u_init("u_init", app_initializer->getComponentDatabase("u_init"), grid_geometry);
        rho_init.setDataOnPatchHierarchy(rho_cur_idx, rho_cur_var, patch_hierarchy, 0.0);
        u_init.setDataOnPatchHierarchy(v_old_idx, v_old_var, patch_hierarchy, -dt);
        u_init.setDataOnPatchHierarchy(v_cur_idx, v_cur_var, patch_hierarchy, 0.0);
        u_init.setDataOnPatchHierarchy(v_new_idx, v_new_var, patch_hierarchy, dt);

        std::vector<std::unique_ptr<RobinBcCoefStrategy<NDIM>>> u_bc_storage;
        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            u_bc_storage.push_back(std::make_unique<muParserRobinBcCoefs>(
                "u_bc_" + std::to_string(d),
                app_initializer->getComponentDatabase("VelocityBcCoefs_" + std::to_string(d)),
                grid_geometry));
            u_bc_coefs[d] = u_bc_storage.back().get();
        }
        std::unique_ptr<RobinBcCoefStrategy<NDIM>> rho_bc = std::make_unique<muParserRobinBcCoefs>(
            "rho_bc", app_initializer->getComponentDatabase("DensityBcCoefs"), grid_geometry);
        std::vector<RobinBcCoefStrategy<NDIM>*> rho_bc_coefs(NDIM, rho_bc.get());

        Pointer<HierarchyMathOps> hier_math_ops = new HierarchyMathOps("hier_math_ops", patch_hierarchy, 0, finest_ln);
        const int h_sc_idx = hier_math_ops->getSideWeightPatchDescriptorIndex();

        auto configure_and_run = [&](const std::string& name,
                                     Pointer<Database> db,
                                     const int n_idx) -> Pointer<INSVCStaggeredConservativeMassMomentumSSPRKIntegrator>
        {
            Pointer<INSVCStaggeredConservativeMassMomentumSSPRKIntegrator> integrator =
                new INSVCStaggeredConservativeMassMomentumSSPRKIntegrator(name, db);
            integrator->setHierarchyMathOps(hier_math_ops);
            integrator->setSideCenteredVelocityBoundaryConditions(u_bc_coefs);
            integrator->setDensityBoundaryConditions(rho_bc_coefs);
            integrator->setTimeInterval(0.0, dt);
            integrator->setSolutionTime(0.5 * dt);
            integrator->initializeSTSIntegrator(patch_hierarchy);
            integrator->setTimeInterval(0.0, dt);
            integrator->setSolutionTime(0.5 * dt);
            integrator->setCycleNumber(1);
            integrator->setPreviousTimeStepSize(dt);
            integrator->setDensityPatchDataIndex(rho_cur_idx);
            integrator->setFluidVelocityPatchDataIndices(v_old_idx, v_cur_idx, v_new_idx);
            integrator->setConvectiveDerivativePatchDataIndex(n_idx);
            integrator->integrate(dt);
            return integrator;
        };

        Pointer<INSVCStaggeredConservativeMassMomentumSSPRKIntegrator> ref_integrator =
            configure_and_run("ref_integrator", input_db->getDatabase("ref_integrator_db"), n_ref_idx);
        Pointer<INSVCStaggeredConservativeMassMomentumSSPRKIntegrator> fast_integrator =
            configure_and_run("fast_integrator", input_db->getDatabase("fast_integrator_db"), n_fast_idx);

        const int rho_ref_idx = ref_integrator->getUpdatedSideCenteredDensityPatchDataIndex();
        const int rho_fast_idx = fast_integrator->getUpdatedSideCenteredDensityPatchDataIndex();

        const long expected_hits = input_db->getInteger("expected_fast_path_hits");
        const long expected_fallbacks = input_db->getInteger("expected_fallback_count");
        const long hits = fast_integrator->getConstantDensityFastPathHits();
        const long fallbacks = fast_integrator->getConstantDensityFallbackCount();
        if (hits != expected_hits || fallbacks != expected_fallbacks)
        {
            TBOX_ERROR("Unexpected fast-path counts: hits=" << hits << ", fallbacks=" << fallbacks << ".\n");
        }
        if (ref_integrator->getConstantDensityFastPathHits() != 0 ||
            ref_integrator->getConstantDensityFallbackCount() != 0)
        {
            TBOX_ERROR("Reference integrator reported nonzero fast-path counts.\n");
        }

        HierarchySideDataOpsReal<NDIM, double> sc_ops(patch_hierarchy, 0, finest_ln);
        sc_ops.subtract(diff_idx, rho_fast_idx, rho_ref_idx);
        const double rho_diff_max = sc_ops.maxNorm(diff_idx);
        sc_ops.subtract(diff_idx, n_fast_idx, n_ref_idx);
        const double n_diff_max = sc_ops.maxNorm(diff_idx);
        const double rho_min = sc_ops.min(rho_fast_idx);
        const double rho_max = sc_ops.max(rho_fast_idx);
        const double rho_l2 = sc_ops.L2Norm(rho_fast_idx, h_sc_idx);
        const double n_l2 = sc_ops.L2Norm(n_fast_idx, h_sc_idx);
        const double n_max = sc_ops.maxNorm(n_fast_idx);

        if (!std::isfinite(rho_diff_max) || !std::isfinite(n_diff_max) || !std::isfinite(rho_min) ||
            !std::isfinite(rho_max) || !std::isfinite(rho_l2) || !std::isfinite(n_l2) || !std::isfinite(n_max))
        {
            TBOX_ERROR("Nonfinite norm observed after SSPRK3 integration.\n");
        }
        if (rho_diff_max > 1.0e-12 || n_diff_max > 1.0e-12)
        {
            TBOX_ERROR("Fast-path and reference integration diverged: rho_diff_max=" << rho_diff_max << ", n_diff_max="
                                                                                     << n_diff_max << ".\n");
        }

        if (IBTK_MPI::getRank() == 0)
        {
            std::ofstream out("output");
            out << std::setprecision(10) << "fast_path_hits = " << hits << '\n'
                << "fallback_count = " << fallbacks << '\n'
                << "rho_min = " << rho_min << ", rho_max = " << rho_max << ", rho_l2 = " << rho_l2 << '\n'
                << "N_l2 = " << n_l2 << ", N_max = " << n_max << '\n';
        }
    }
    return 0;
}
