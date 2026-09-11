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
#include <ibtk/IBTK_MPI.h>

#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIndex.h>
#include <CoupledApplication.h>
#include <DiagnosticUtilities.h>

#include <fstream>

#include <ibamr/app_namespaces.h>

namespace
{
class MeltingBubbles : public PhaseChangeExamples::CoupledApplication
{
public:
    explicit MeltingBubbles(Pointer<AppInitializer> app_initializer) : CoupledApplication(app_initializer)
    {
        useEnthalpy();
        registerLevelSet();
        registerLiquidFraction();
        registerLiquidFractionGradient();
        registerExtrapolatedLiquidFraction();
        registerSpecificEnthalpy();
        registerHeaviside();
        registerLevelSetForExtrapolation();
        registerTransportFields();
        registerMaterialFields();
        registerLevelSetTagging();
        registerLiquidFractionTagging();
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
        registerGravityAndSurfaceTension(d_extrapolated_liquid_fraction_var);
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
        d_mass->allocateData(loop_time);

        double left_bubble_volume = 0.0;
        double right_bubble_volume = 0.0;
        HierarchyMathOps hier_math_ops("HierarchyMathOps", d_patch_hierarchy, 0, d_mass->getFinestLevel());
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        for (int ln = 0; ln <= d_mass->getFinestLevel(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = d_patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
                const double* patch_dx = patch_geom->getDx();
                const double* patch_X_lower = patch_geom->getXLower();
                const hier::Index<NDIM>& patch_lower_idx = patch_box.lower();

                Pointer<CellData<NDIM, double>> H_data = patch->getPatchData(d_mass->getHeavisideIndex());
                Pointer<CellData<NDIM, double>> wgt_data = patch->getPatchData(wgt_cc_idx);
                for (Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    CellIndex<NDIM> ci(it());

                    IBTK::Vector coord = IBTK::Vector::Zero();
                    for (int d = 0; d < NDIM; ++d)
                    {
                        coord[d] =
                            patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
                    }
                    if (coord[0] < 5.0e-3 && coord[1] < 4.0e-3)
                    {
                        left_bubble_volume += (1.0 - (*H_data)(ci)) * (*wgt_data)(ci);
                    }
                    else if (coord[0] > 5.0e-3 && coord[1] < 4.0e-3)
                    {
                        right_bubble_volume += (1.0 - (*H_data)(ci)) * (*wgt_data)(ci);
                    }
                }
            }
        }
        std::vector<double> bubble_integrals{ left_bubble_volume, right_bubble_volume };
        IBTK_MPI::sumReduction(&bubble_integrals[0], bubble_integrals.size());

        d_mass->multiply();
        const double mass = d_mass->integral(wgt_cc_idx);

        if (SAMRAI_MPI::getRank() == 0)
        {
            d_pcm_mass_file << loop_time << "\t" << mass << "\t" << bubble_integrals[0] << "\t" << bubble_integrals[1]
                            << std::endl;
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
    return std::make_unique<MeltingBubbles>(app_initializer);
}
} // namespace

int
main(int argc, char* argv[])
{
    return PhaseChangeExamples::run_coupled(argc, argv, &create_application);
}
