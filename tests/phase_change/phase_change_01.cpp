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
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/AllenCahnHierarchyIntegrator.h>
#include <ibamr/EnthalpyHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <cmath>
#include <fstream>
#include <iomanip>

#include <ibamr/app_namespaces.h>

namespace
{
// Material property that varies linearly with the liquid fraction.
struct LiquidFractionProperty
{
    Pointer<AdvDiffHierarchyIntegrator> integrator;
    Pointer<CellVariable<NDIM, double>> lf_var;
    double solid_value;
    double liquid_value;
};

void
set_property_from_liquid_fraction(const int property_idx,
                                  Pointer<Variable<NDIM>> /*property_var*/,
                                  Pointer<HierarchyMathOps> hier_math_ops,
                                  const int /*cycle_num*/,
                                  const double time,
                                  const double current_time,
                                  const double new_time,
                                  void* ctx)
{
    LiquidFractionProperty* property = static_cast<LiquidFractionProperty*>(ctx);
    Pointer<VariableContext> lf_ctx;
    if (IBTK::rel_equal_eps(time, current_time))
    {
        lf_ctx = property->integrator->getCurrentContext();
    }
    else if (IBTK::rel_equal_eps(time, new_time))
    {
        lf_ctx = property->integrator->getNewContext();
    }
    else
    {
        TBOX_ERROR("set_property_from_liquid_fraction(): unexpected time " << time << "\n");
    }
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int lf_idx = var_db->mapVariableAndContextToIndex(property->lf_var, lf_ctx);

    Pointer<PatchHierarchy<NDIM>> hierarchy = hier_math_ops->getPatchHierarchy();
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> lf_data = patch->getPatchData(lf_idx);
            Pointer<CellData<NDIM, double>> property_data = patch->getPatchData(property_idx);
            for (Box<NDIM>::Iterator it(patch->getBox()); it; it++)
            {
                const CellIndex<NDIM> ci(it());
                (*property_data)(ci) =
                    property->solid_value + (property->liquid_value - property->solid_value) * (*lf_data)(ci);
            }
        }
    }
}

// Tag cells in the mushy region, where the liquid fraction is neither zero nor one.
struct MushyRegionTagger
{
    Pointer<AdvDiffHierarchyIntegrator> integrator;
    Pointer<CellVariable<NDIM, double>> lf_var;
    double min_lf;
    double max_lf;
};

void
tag_mushy_region(const Pointer<BasePatchHierarchy<NDIM>> hierarchy,
                 const int level_number,
                 const double /*error_data_time*/,
                 const int tag_index,
                 const bool /*initial_time*/,
                 const bool /*uses_richardson_extrapolation_too*/,
                 void* ctx)
{
    MushyRegionTagger* tagger = static_cast<MushyRegionTagger*>(ctx);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int lf_idx = var_db->mapVariableAndContextToIndex(tagger->lf_var, tagger->integrator->getCurrentContext());

    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<CellData<NDIM, int>> tag_data = patch->getPatchData(tag_index);
        Pointer<CellData<NDIM, double>> lf_data = patch->getPatchData(lf_idx);
        for (CellIterator<NDIM> ic(patch->getBox()); ic; ic++)
        {
            const double lf = (*lf_data)(ic());
            if (lf > tagger->min_lf && lf < tagger->max_lf)
            {
                (*tag_data)(ic()) = 1;
            }
        }
    }
}
} // namespace

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

    { // cleanup dynamically allocated objects prior to shutdown

        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "phase_change.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const std::string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        // Create the phase change integrator, either Allen-Cahn or enthalpy.
        const std::string integrator_type = input_db->getString("INTEGRATOR_TYPE");
        Pointer<PhaseChangeHierarchyIntegrator> time_integrator;
        Pointer<AllenCahnHierarchyIntegrator> allen_cahn_integrator;
        if (integrator_type == "ALLEN_CAHN")
        {
            allen_cahn_integrator = new AllenCahnHierarchyIntegrator(
                "AllenCahnHierarchyIntegrator", app_initializer->getComponentDatabase("AllenCahnHierarchyIntegrator"));
            time_integrator = allen_cahn_integrator;
        }
        else if (integrator_type == "ENTHALPY")
        {
            time_integrator = new EnthalpyHierarchyIntegrator(
                "EnthalpyHierarchyIntegrator", app_initializer->getComponentDatabase("EnthalpyHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("unknown INTEGRATOR_TYPE: " << integrator_type << "\n");
        }

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

        // Register the liquid fraction, its gradient, the Heaviside function, and the temperature.
        Pointer<CellVariable<NDIM, double>> lf_var = new CellVariable<NDIM, double>("lf_var");
        time_integrator->registerLiquidFractionVariable(lf_var, true);

        Pointer<CellVariable<NDIM, double>> lf_gradient_var = new CellVariable<NDIM, double>("lf_gradient_var", NDIM);
        time_integrator->registerLiquidFractionGradientVariable(lf_gradient_var, true);

        Pointer<CellVariable<NDIM, double>> H_var = new CellVariable<NDIM, double>("heaviside_var");
        time_integrator->registerTransportedQuantity(H_var, true);
        time_integrator->setDiffusionCoefficient(H_var, 0.0);
        time_integrator->registerHeavisideVariable(H_var);

        Pointer<CellVariable<NDIM, double>> T_var = new CellVariable<NDIM, double>("Temperature");
        time_integrator->registerTemperatureVariable(T_var, true);

        Pointer<CellVariable<NDIM, double>> rho_var = new CellVariable<NDIM, double>("rho_cc_var");
        time_integrator->registerDensityVariable(rho_var, true);

        Pointer<CellVariable<NDIM, double>> Cp_var = new CellVariable<NDIM, double>("Cp");
        time_integrator->registerSpecificHeatVariable(Cp_var, true);

        // Initial conditions.
        Pointer<CartGridFunction> H_init = new muParserCartGridFunction(
            "H_init", app_initializer->getComponentDatabase("HeavisideInitialConditions"), grid_geometry);
        time_integrator->setInitialConditions(H_var, H_init);

        Pointer<CartGridFunction> T_init = new muParserCartGridFunction(
            "T_init", app_initializer->getComponentDatabase("TemperatureInitialConditions"), grid_geometry);
        time_integrator->setTemperatureInitialCondition(T_var, T_init);

        Pointer<CartGridFunction> lf_init = new muParserCartGridFunction(
            "lf_init", app_initializer->getComponentDatabase("LiquidFractionInitialConditions"), grid_geometry);
        time_integrator->setLiquidFractionInitialCondition(lf_var, lf_init);

        // Boundary conditions.
        RobinBcCoefStrategy<NDIM>* H_bc_coef = new muParserRobinBcCoefs(
            "H_bc_coef", app_initializer->getComponentDatabase("HeavisideBcCoefs"), grid_geometry);
        time_integrator->setPhysicalBcCoef(H_var, H_bc_coef);

        RobinBcCoefStrategy<NDIM>* T_bc_coef = new muParserRobinBcCoefs(
            "T_bc_coef", app_initializer->getComponentDatabase("TemperatureBcCoefs"), grid_geometry);
        time_integrator->setTemperaturePhysicalBcCoef(T_var, T_bc_coef);

        RobinBcCoefStrategy<NDIM>* lf_bc_coef = nullptr;
        if (allen_cahn_integrator)
        {
            lf_bc_coef = new muParserRobinBcCoefs(
                "lf_bc_coef", app_initializer->getComponentDatabase("LiquidFractionBcCoefs"), grid_geometry);
            allen_cahn_integrator->setLiquidFractionPhysicalBcCoef(lf_var, lf_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* k_bc_coef = new muParserRobinBcCoefs(
            "k_bc_coef", app_initializer->getComponentDatabase("ThermalConductivityBcCoefs"), grid_geometry);
        time_integrator->registerThermalConductivityBoundaryConditions(k_bc_coef);

        // Density, specific heat, and thermal conductivity vary linearly with the liquid fraction.
        LiquidFractionProperty rho{
            time_integrator, lf_var, input_db->getDouble("RHO_S"), input_db->getDouble("RHO_L")
        };
        LiquidFractionProperty Cp{ time_integrator, lf_var, input_db->getDouble("CP_S"), input_db->getDouble("CP_L") };
        LiquidFractionProperty kappa{
            time_integrator, lf_var, input_db->getDouble("KAPPA_S"), input_db->getDouble("KAPPA_L")
        };
        time_integrator->registerResetDensityFcn(&set_property_from_liquid_fraction, static_cast<void*>(&rho));
        time_integrator->registerResetSpecificHeatFcn(&set_property_from_liquid_fraction, static_cast<void*>(&Cp));
        time_integrator->registerResetDiffusionCoefficientFcn(&set_property_from_liquid_fraction,
                                                              static_cast<void*>(&kappa));

        // Refine where the liquid fraction is between zero and one.
        MushyRegionTagger tagger{
            time_integrator, lf_var, input_db->getDouble("MIN_TAG_VAL"), input_db->getDouble("MAX_TAG_VAL")
        };
        time_integrator->registerApplyGradientDetectorCallback(&tag_mushy_region, static_cast<void*>(&tagger));

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Main time step loop.
        double loop_time = time_integrator->getIntegratorTime();
        const double loop_time_end = time_integrator->getEndTime();
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            const int iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "At beginning of timestep # " << iteration_num << ", simulation time is " << loop_time << "\n";

            const double dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            if (dump_restart_data && (iteration_num + 1) % restart_dump_interval == 0)
            {
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num + 1);
            }
        }

        // Write composite norms of the solution.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);

        std::ofstream out;
        if (IBTK_MPI::getRank() == 0)
        {
            out.open("output");
        }
        out << "Number of levels at time " << loop_time << ": " << finest_ln + 1 << "\n";
        const std::vector<std::pair<std::string, Pointer<CellVariable<NDIM, double>>>> variables = { { "T", T_var },
                                                                                                     { "lf", lf_var } };
        for (const auto& name_and_var : variables)
        {
            const int idx =
                var_db->mapVariableAndContextToIndex(name_and_var.second, time_integrator->getCurrentContext());
            const double L1_norm = hier_cc_data_ops.L1Norm(idx, wgt_cc_idx);
            const double L2_norm = hier_cc_data_ops.L2Norm(idx, wgt_cc_idx);
            const double max_norm = hier_cc_data_ops.maxNorm(idx, wgt_cc_idx);
            if (!std::isfinite(L1_norm) || !std::isfinite(L2_norm) || !std::isfinite(max_norm))
            {
                TBOX_ERROR("nonfinite norm of " << name_and_var.first << "\n");
            }
            out << std::setprecision(10) << name_and_var.first << " norms:\n"
                << "  L1-norm:  " << L1_norm << "\n"
                << "  L2-norm:  " << L2_norm << "\n"
                << "  max-norm: " << max_norm << "\n";
        }

        delete H_bc_coef;
        delete T_bc_coef;
        delete lf_bc_coef;
        delete k_bc_coef;

    } // cleanup dynamically allocated objects prior to shutdown
} // main
