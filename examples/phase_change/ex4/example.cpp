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

#include <ibtk/HierarchyMathOps.h>

#include <CoupledApplication.h>
#include <DiagnosticUtilities.h>

#include <fstream>

#include <ibamr/app_namespaces.h>

namespace
{
class MeltingCylinder : public PhaseChangeExamples::CoupledApplication
{
public:
    explicit MeltingCylinder(Pointer<AppInitializer> app_initializer) : CoupledApplication(app_initializer)
    {
        useEnthalpy();
        registerLevelSet();
        registerLiquidFraction();
        registerSpecificEnthalpy();
        registerHeaviside();
        registerTransportFields();
        registerMaterialFields();
        setHeavisideBoundary();
        setTemperatureBoundary();
        setEnthalpyBoundary();
        setLiquidFractionBoundary();
        setVelocityBoundary();
        setDensityBoundary();
        setViscosityBoundary();
        setConductivityBoundary();
        setLevelSetBoundary();
        registerMaterialProperties();
        registerPhaseChangeSources(d_heaviside_var->getName() + "_F");
        registerGravityAndSurfaceTension(d_liquid_fraction_var);
        registerSolidDrag();
    }

private:
    void initializeDiagnostics(const double loop_time) override
    {
        d_mass = std::make_unique<PhaseChangeExamples::PhaseMassDiagnostic>(
            d_patch_hierarchy, d_phase_change_integrator, d_density_var, d_heaviside_var, loop_time);
        if (SAMRAI_MPI::getRank() == 0)
        {
            const std::string file_name = d_input_db->getString("MASS_CONSERVATION_FILE_NAME");
            PhaseChangeExamples::open_diagnostic_file(d_pcm_mass_file, file_name);
        }
    }
    void postprocessStep(const double loop_time) override
    {
        d_mass->multiply();
        HierarchyMathOps hier_math_ops("HierarchyMathOps", d_patch_hierarchy, 0, d_mass->getFinestLevel());
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const double mass = d_mass->integral(wgt_cc_idx);

        if (SAMRAI_MPI::getRank() == 0)
        {
            d_pcm_mass_file << loop_time << "\t" << mass << "\n";
        }
    }
    void finalizeDiagnostics() override
    {
        d_mass.reset();
        if (SAMRAI_MPI::getRank() == 0)
        {
            d_pcm_mass_file.close();
        }
    }
    std::unique_ptr<PhaseChangeExamples::PhaseMassDiagnostic> d_mass;
    std::ofstream d_pcm_mass_file;
};

std::unique_ptr<PhaseChangeExamples::CoupledApplication>
create_application(Pointer<AppInitializer> app_initializer)
{
    return std::make_unique<MeltingCylinder>(app_initializer);
}
} // namespace

int
main(int argc, char* argv[])
{
    return PhaseChangeExamples::run_coupled(argc, argv, &create_application);
}
