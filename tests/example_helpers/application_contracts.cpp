// ---------------------------------------------------------------------
//
// Copyright (c) 2026 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTK_MPI.h>

#include <tbox/PIO.h>
#include <tbox/RestartManager.h>

#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIndex.h>
#include <CoupledApplication.h>
#include <DiagnosticUtilities.h>
#include <VariableDatabase.h>

#include <cmath>
#include <iomanip>
#include <sstream>

#include <ibamr/app_namespaces.h>

namespace
{
// Run the shared example lifecycle against a discrete heat-equation solution.
// The x direction is periodic; the y boundaries have zero normal heat flux.
class ApplicationContract : public PhaseChangeExamples::CoupledApplication
{
public:
    explicit ApplicationContract(Pointer<AppInitializer> app_initializer);

private:
    void initializeDiagnostics(double time) override;
    void postprocessStep(double time) override;
    void finalizeDiagnostics() override;

    std::unique_ptr<PhaseChangeExamples::PhaseMassDiagnostic> d_mass;
    std::ostringstream d_results;
    int d_initial_step = -1, d_postprocess_count = 0;
};

ApplicationContract::ApplicationContract(Pointer<AppInitializer> app_initializer) : CoupledApplication(app_initializer)
{
    if (d_input_db->keyExists("EnthalpyHierarchyIntegrator"))
    {
        useEnthalpy();
    }
    else
    {
        useAllenCahn();
    }
    registerLevelSet();
    registerLiquidFraction();
    registerLiquidFractionGradient();
    if (d_enthalpy_integrator)
    {
        registerSpecificEnthalpy();
    }
    registerHeaviside();
    registerTransportFields();
    registerMaterialFields();
    setHeavisideBoundary();
    setTemperatureBoundary();
    if (d_enthalpy_integrator)
    {
        setEnthalpyBoundary();
    }
    else
    {
        setSpecificHeatBoundary();
    }
    setLiquidFractionBoundary();
    setConductivityBoundary();
    setVelocityBoundary();
    setDensityBoundary();
    setViscosityBoundary();
    setLevelSetBoundary();
    registerMaterialProperties();
    registerPhaseChangeSources("F");
}

void
ApplicationContract::initializeDiagnostics(const double time)
{
    d_initial_step = d_time_integrator->getIntegratorStep();
    if (d_initial_step != (RestartManager::getManager()->isFromRestart() ? 2 : 0))
    {
        TBOX_ERROR("Unexpected initial step in shared application\n");
    }
    d_mass = std::make_unique<PhaseChangeExamples::PhaseMassDiagnostic>(
        d_patch_hierarchy, d_phase_change_integrator, d_density_var, d_heaviside_var, time);
}

void
ApplicationContract::postprocessStep(const double time)
{
    ++d_postprocess_count;
    const int step = d_time_integrator->getIntegratorStep();
    if (!d_mass || step != d_initial_step + d_postprocess_count || std::abs(time - step * 0.001) > 1.0e-14)
    {
        TBOX_ERROR("Incorrect diagnostic callback order or timestep\n");
    }
    auto* db = VariableDatabase<NDIM>::getDatabase();
    const auto context = d_phase_change_integrator->getCurrentContext();
    const int T_idx = db->mapVariableAndContextToIndex(d_temperature_var, context);
    const int lf_idx = db->mapVariableAndContextToIndex(d_liquid_fraction_var, context);
    const int rho_idx = db->mapVariableAndContextToIndex(d_density_var, context);
    const int Cp_idx = db->mapVariableAndContextToIndex(d_specific_heat_var, context);
    const double pi = std::acos(-1.0);
    // k/(rho Cp) = 0.06/(2*3), with Crank-Nicolson diffusion on a 16-cell period.
    const double lambda = 4.0 * 0.01 * std::pow(16.0 * std::sin(pi / 16.0), 2.0);
    const double amplitude = 0.1 * std::pow((1.0 - 0.0005 * lambda) / (1.0 + 0.0005 * lambda), step);
    double temperature_error = 0.0, material_error = 0.0;
    for (int ln = 0; ln <= d_patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM>> geom = patch->getPatchGeometry();
            Pointer<CellData<NDIM, double>> T = patch->getPatchData(T_idx);
            Pointer<CellData<NDIM, double>> lf = patch->getPatchData(lf_idx);
            Pointer<CellData<NDIM, double>> H = patch->getPatchData(d_mass->getHeavisideIndex());
            Pointer<CellData<NDIM, double>> rho = patch->getPatchData(rho_idx);
            Pointer<CellData<NDIM, double>> Cp = patch->getPatchData(Cp_idx);
            for (Box<NDIM>::Iterator i(patch->getBox()); i; i++)
            {
                const CellIndex<NDIM> ci(i());
                for (const double value : { (*T)(ci), (*lf)(ci), (*H)(ci), (*rho)(ci), (*Cp)(ci) })
                {
                    if (!std::isfinite(value))
                    {
                        TBOX_ERROR("Nonfinite shared application field\n");
                    }
                }
                const double x = geom->getXLower()[0] + (ci(0) - patch->getBox().lower(0) + 0.5) * geom->getDx()[0];
                temperature_error =
                    std::max(temperature_error, std::abs((*T)(ci) - (2.0 + amplitude * std::cos(2.0 * pi * x))));
                material_error = std::max({ material_error,
                                            std::abs((*lf)(ci)-1.0),
                                            std::abs((*H)(ci)-1.0),
                                            std::abs((*rho)(ci)-2.0),
                                            std::abs((*Cp)(ci)-3.0) });
            }
        }
    }
    HierarchyMathOps math_ops("diagnostic_math", d_patch_hierarchy);
    d_mass->allocateData(time);
    d_mass->multiply();
    const double mass = d_mass->integral(math_ops.getCellWeightPatchDescriptorIndex());
    if (!std::isfinite(mass))
    {
        TBOX_ERROR("Nonfinite phase mass diagnostic\n");
    }
    d_results << std::setprecision(13) << "Step " << step
              << " temperature/material/mass error = " << IBTK_MPI::maxReduction(temperature_error) << ' '
              << IBTK_MPI::maxReduction(material_error) << ' ' << std::abs(mass - 2.0) << '\n';
}

void
ApplicationContract::finalizeDiagnostics()
{
    if (!d_mass || d_postprocess_count != 4 - d_initial_step || d_time_integrator->getIntegratorStep() != 4)
    {
        TBOX_ERROR("Shared application did not complete its diagnostic lifecycle\n");
    }
    d_mass.reset();
    PIO::logOnlyNodeZero("output");
    plog << d_results.str();
}

std::unique_ptr<PhaseChangeExamples::CoupledApplication>
create_application(Pointer<AppInitializer> app_initializer)
{
    return std::make_unique<ApplicationContract>(app_initializer);
}
} // namespace

int
main(int argc, char* argv[])
{
    return PhaseChangeExamples::run_coupled(argc, argv, &create_application);
}
