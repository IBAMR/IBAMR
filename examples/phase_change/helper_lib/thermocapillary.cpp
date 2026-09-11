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

#include <ibamr/MarangoniSurfaceTensionForceFunction.h>

#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/HierarchyMathOps.h>

#include <CapillaryForces.h>
#include <VariableDatabase.h>

#include <cmath>

#include "Applications.h"
#include "CoupledApplication.h"
#include "DiagnosticUtilities.h"

#include <ibamr/app_namespaces.h>

namespace PhaseChangeExamples
{
namespace
{
class ThermocapillaryApplication : public CoupledApplication
{
public:
    explicit ThermocapillaryApplication(Pointer<AppInitializer> app_initializer) : CoupledApplication(app_initializer)
    {
        if (d_input_db->keyExists("EnthalpyHierarchyIntegrator"))
        {
            useEnthalpy();
        }
        else
        {
            useAllenCahn();
        }
        d_bubble_radius = d_input_db->getDouble("R");
        registerLevelSet();
        registerLiquidFraction();
        if (d_enthalpy_integrator)
        {
            registerLiquidFractionGradient();
            registerExtrapolatedLiquidFraction();
            registerSpecificEnthalpy();
            registerLevelSetForExtrapolation();
        }
        registerHeaviside();
        registerTransportFields();
        registerMaterialFields();
        if (d_enthalpy_integrator)
        {
            registerLevelSetTagging();
            registerLiquidFractionTagging();
        }
        setHeavisideBoundary();
        setTemperatureBoundary();
        if (d_enthalpy_integrator)
        {
            setEnthalpyBoundary();
        }
        setLiquidFractionBoundary();
        setVelocityBoundary();
        setDensityBoundary();
        setViscosityBoundary();
        if (d_allen_cahn_integrator)
        {
            setSpecificHeatBoundary();
        }
        setConductivityBoundary();
        setLevelSetBoundary();
        if (d_enthalpy_integrator)
        {
            registerMaterialProperties();
        }
        else
        {
            registerMaterialProperties(0.0, 0.0, 0.0, 0.0);
        }
        const double sigma_0 = d_input_db->getDouble("SIGMA_0");
        d_dsigma_dt = d_input_db->getDouble("DSIGMA_DT_0");
        const double ref_temperature = d_input_db->getDouble("REFERENCE_TEMP_SIGMA");
        d_temperature_gradient = d_input_db->getDouble("TEMPERATURE_GRADIENT");
        registerPhaseChangeSources(d_heaviside_var->getName() + "_F");
        // Register surface tension force.
        Pointer<SurfaceTensionForceFunction> surface_tension_force = new MarangoniSurfaceTensionForceFunction(
            "MarangoniSurfaceTensionForceFunction",
            d_app_initializer->getComponentDatabase("MarangoniSurfaceTensionForceFunction"),
            d_phase_change_integrator,
            d_level_set_var,
            d_temperature_var,
            d_temperature_bc_coef.get());

        registerForceMask(surface_tension_force, d_extrapolated_liquid_fraction_var);

        // Register variable coefficient surface tension.
        d_coefficients = std::make_unique<MultiphaseExamples::SurfaceTensionCoefficients>(
            d_temperature_var, d_phase_change_integrator->getScratchContext(), sigma_0, d_dsigma_dt, ref_temperature);

        surface_tension_force->registerSurfaceTensionCoefficientFunction(
            &MultiphaseExamples::SurfaceTensionCoefficients::compute_surface_tension_coef_function,
            static_cast<void*>(d_coefficients.get()));

        // Register variable marangoni coefficient dsigma_dT.
        Pointer<MarangoniSurfaceTensionForceFunction> marangoni_force = surface_tension_force;
        marangoni_force->registerMarangoniCoefficientFunction(
            &MultiphaseExamples::SurfaceTensionCoefficients::compute_marangoni_coef_function,
            static_cast<void*>(d_coefficients.get()));

        Pointer<CartGridFunctionSet> eul_forces = new CartGridFunctionSet("eulerian_forces");
        eul_forces->addFunction(surface_tension_force);
        d_time_integrator->registerBodyForceFunction(eul_forces);
    }

private:
    void setThermalInitialConditions() override
    {
        if (d_input_db->keyExists("TemperatureInitialConditions"))
        {
            setTemperatureInitialCondition();
        }
        if (d_input_db->keyExists("LiquidFractionInitialConditions"))
        {
            setLiquidFractionInitialCondition();
        }
    }
    void initializeDiagnostics(const double loop_time) override
    {
        d_var_db = VariableDatabase<NDIM>::getDatabase();
        d_h_idx =
            d_var_db->mapVariableAndContextToIndex(d_heaviside_var, d_phase_change_integrator->getCurrentContext());
        d_phi_idx =
            d_var_db->mapVariableAndContextToIndex(d_level_set_var, d_phase_change_integrator->getCurrentContext());
        d_h_cloned_idx = d_var_db->registerClonedPatchDataIndex(d_level_set_var, d_phi_idx);

        // Interpolating side-centered velocity to cell-centered
        d_u_sc_var = d_time_integrator->getVelocityVariable();
        d_u_sc_idx = d_var_db->mapVariableAndContextToIndex(d_u_sc_var, d_time_integrator->getCurrentContext());

        d_u_cc_var = new CellVariable<NDIM, double>("U_cc", NDIM);
        d_u_cc_idx = d_var_db->registerVariableAndContext(d_u_cc_var, d_var_db->getContext("U_cc"), 0);

        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> v_var;
        v_var = new CellVariable<NDIM, double>("v_cc");
        d_v_idx = d_var_db->registerVariableAndContext(v_var, d_var_db->getContext("v_cc"), 0);

        d_coarsest_ln = 0;
        d_finest_ln = d_patch_hierarchy->getFinestLevelNumber();
        allocate_diagnostic_data(
            d_patch_hierarchy, { d_h_cloned_idx, d_u_cc_idx, d_v_idx }, d_coarsest_ln, d_finest_ln, loop_time);

        // File to write rise velocity of a bubble.

        if (SAMRAI_MPI::getRank() == 0)
        {
            string output_file_name = d_input_db->getString("OUTPUT_FILE_NAME");
            open_diagnostic_file(d_output_file, output_file_name);
        }
    }
    void postprocessStep(const double loop_time) override
    {
        HierarchyMathOps hier_math_ops("HierarchyMathOps", d_patch_hierarchy, d_coarsest_ln, d_finest_ln);

        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(d_patch_hierarchy, d_coarsest_ln, d_finest_ln);

        allocate_diagnostic_data(
            d_patch_hierarchy, { d_h_cloned_idx, d_u_cc_idx, d_v_idx }, d_coarsest_ln, d_finest_ln, loop_time);

        hier_math_ops.interp(d_u_cc_idx, d_u_cc_var, d_u_sc_idx, d_u_sc_var, nullptr, loop_time, true);

        const double U_ref = std::abs(d_dsigma_dt) * d_temperature_gradient * d_bubble_radius / d_mu_liquid;
        const double t_ref = d_bubble_radius / U_ref;

        // Calculate Heaviside function and compute the rise velocity.
        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = d_patch_hierarchy->getPatchLevel(ln);

            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();

                Pointer<CellData<NDIM, double>> U_cc_data = patch->getPatchData(d_u_cc_idx);
                Pointer<CellData<NDIM, double>> v_data = patch->getPatchData(d_v_idx);
                Pointer<CellData<NDIM, double>> H_cloned_data = patch->getPatchData(d_h_cloned_idx);

                Pointer<CellData<NDIM, double>> H_data = patch->getPatchData(d_h_idx);

                for (Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    CellIndex<NDIM> ci(it());
                    (*v_data)(ci) = (*U_cc_data)(ci, 1) / U_ref; // non_dimensional velocity
                    (*H_cloned_data)(ci) = 1.0 - (*H_data)(ci);
                }
            }
        }

        double vol = hier_cc_data_ops.integral(d_h_cloned_idx, wgt_cc_idx);
        hier_cc_data_ops.multiply(d_h_cloned_idx, d_h_cloned_idx, wgt_cc_idx);
        double v_integral = hier_cc_data_ops.integral(d_v_idx, d_h_cloned_idx);

        if (SAMRAI_MPI::getRank() == 0)
        {
            d_output_file << loop_time / t_ref << "\t" << v_integral / vol << "\n";
        }
    }
    void finalizeDiagnostics() override
    {
        deallocate_diagnostic_data(
            d_patch_hierarchy, { d_h_cloned_idx, d_u_cc_idx, d_v_idx }, d_coarsest_ln, d_finest_ln);
        if (SAMRAI_MPI::getRank() == 0)
        {
            d_output_file.close();
        }
    }
    VariableDatabase<NDIM>* d_var_db = nullptr;
    int d_h_idx = -1;
    int d_phi_idx = -1;
    int d_h_cloned_idx = -1;
    int d_u_sc_idx = -1;
    int d_u_cc_idx = -1;
    int d_v_idx = -1;
    int d_coarsest_ln = -1;
    int d_finest_ln = -1;
    Pointer<SideVariable<NDIM, double>> d_u_sc_var;
    Pointer<CellVariable<NDIM, double>> d_u_cc_var;
    std::ofstream d_output_file;
    double d_bubble_radius = 0.0, d_dsigma_dt = 0.0, d_temperature_gradient = 0.0;
    std::unique_ptr<MultiphaseExamples::SurfaceTensionCoefficients> d_coefficients;
};

std::unique_ptr<CoupledApplication>
create_application(Pointer<AppInitializer> app_initializer)
{
    return std::make_unique<ThermocapillaryApplication>(app_initializer);
}
} // namespace

int
run_thermocapillary(int argc, char* argv[])
{
    return run_coupled(argc, argv, &create_application);
}
} // namespace PhaseChangeExamples
