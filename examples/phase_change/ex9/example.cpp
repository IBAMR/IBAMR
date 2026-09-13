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

#include <ibamr/LaserSourceFunction.h>

#include <ibtk/HierarchyMathOps.h>

#include <CoupledApplication.h>
#include <HierarchyCellDataOpsReal.h>

#include <cmath>

#include <ibamr/app_namespaces.h>

namespace
{
class StefanHeatFlux : public PhaseChangeExamples::CoupledApplication
{
public:
    explicit StefanHeatFlux(Pointer<AppInitializer> app_initializer) : CoupledApplication(app_initializer)
    {
        useEnthalpy();
        registerLevelSet();
        registerLiquidFraction();
        registerLiquidFractionGradient();
        registerSpecificEnthalpy();
        registerHeaviside();
        registerTransportFields();
        registerMaterialFields();
        registerLiquidFractionTagging();
        setHeavisideBoundary();
        setTemperatureBoundary();
        setEnthalpyBoundary();
        setLiquidFractionBoundary();
        setConductivityBoundary();
        setVelocityBoundary();
        setDensityBoundary();
        setViscosityBoundary();
        setLevelSetBoundary();
        registerMaterialProperties();
        registerPhaseChangeSources("F");
        registerSolidDrag();

        Pointer<LaserSourceFunction> laser_source =
            new LaserSourceFunction("LaserSourceFunction",
                                    d_app_initializer->getComponentDatabase("LaserSourceFunction"),
                                    d_enthalpy_integrator,
                                    d_level_set_var);

        d_melt_temperature = d_input_db->getDouble("MELT_TEMPERATURE");
        d_liquid_temperature = d_input_db->getDouble("LIQUID_TEMPERATURE");
        d_alpha_liquid = d_input_db->getDouble("ALPHA_LIQUID");
        d_alpha_solid = d_input_db->getDouble("ALPHA_SOLID");
        d_lambda = d_input_db->getDouble("LAMBDA");

        laser_source->registerHeatFlux(&StefanHeatFlux::compute_heat_flux, this);
        d_enthalpy_integrator->setEnergyEquationSourceTermFunction(laser_source);
    }

private:
    static void compute_heat_flux(const int F_idx,
                                  Pointer<HierarchyMathOps> hier_math_ops,
                                  int /*integrator_step*/,
                                  const double time,
                                  double /*current_time*/,
                                  double /*new_time*/,
                                  void* ctx)
    {
        const auto* application = static_cast<StefanHeatFlux*>(ctx);
        const double kappa_liquid = application->d_kappa_liquid;
        const double T_m = application->d_melt_temperature;
        const double T_0 = application->d_liquid_temperature;
        const double alpha_liquid = application->d_alpha_liquid;
        const double alpha_solid = application->d_alpha_solid;
        const double lambda = application->d_lambda;

        const double flux = -kappa_liquid * (T_m - T_0) /
                            (erf(lambda * sqrt(alpha_solid / alpha_liquid)) * sqrt(M_PI * alpha_liquid * time));

        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = hier_math_ops->getPatchHierarchy();
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);

        hier_cc_data_ops.scale(F_idx, flux, F_idx);
    }

    double d_melt_temperature = 0.0, d_liquid_temperature = 0.0;
    double d_alpha_liquid = 0.0, d_alpha_solid = 0.0, d_lambda = 0.0;
};

std::unique_ptr<PhaseChangeExamples::CoupledApplication>
create_application(Pointer<AppInitializer> app_initializer)
{
    return std::make_unique<StefanHeatFlux>(app_initializer);
}
} // namespace

int
main(int argc, char* argv[])
{
    return PhaseChangeExamples::run_coupled(argc, argv, &create_application);
}
