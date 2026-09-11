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

#include <ibamr/LevelSetUtilities.h>

#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTK_MPI.h>

#include <CellData.h>
#include <CellIndex.h>
#include <CoupledApplication.h>
#include <DiagnosticUtilities.h>
#include <HierarchySideDataOpsReal.h>
#include <VariableDatabase.h>

#include <fstream>

#include <ibamr/app_namespaces.h>

namespace
{
class AdvectedBubble : public PhaseChangeExamples::CoupledApplication
{
public:
    explicit AdvectedBubble(Pointer<AppInitializer> app_initializer) : CoupledApplication(app_initializer)
    {
        useEnthalpy();
        registerLevelSet();
        // Lagrange multiplier to conserve mass of the phases.
        std::vector<Pointer<CellVariable<NDIM, double>>> ls_vars{ d_level_set_var };
        d_level_set_fixer = std::make_unique<LevelSetUtilities::LevelSetMassLossFixer>(
            "LevelSetMassLossFixer",
            d_phase_change_integrator,
            ls_vars,
            d_app_initializer->getComponentDatabase("LevelSetMassFixer"),
            /*restart*/ true);
        if (d_input_db->getBoolWithDefault("FIX_MASS_LOSS", false))
        {
            d_time_integrator->registerPostprocessIntegrateHierarchyCallback(
                &LevelSetUtilities::fixMassLoss2PhaseFlows, static_cast<void*>(d_level_set_fixer.get()));
        }
        registerLiquidFraction();
        registerLiquidFractionGradient();
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
        registerLevelSetTagging();
        registerLiquidFractionTagging();
        registerPhaseChangeSources("F");
    }

private:
    void initializeDiagnostics(const double loop_time) override
    {
        // Tracking the total enthalpy and the bubble volume.
        d_var_db = VariableDatabase<NDIM>::getDatabase();
        d_rho_cc_idx =
            d_var_db->mapVariableAndContextToIndex(d_density_var, d_phase_change_integrator->getCurrentContext());
        d_rho_sc_idx = d_time_integrator->getLinearOperatorRhoPatchDataIndex();
        d_h_idx =
            d_var_db->mapVariableAndContextToIndex(d_enthalpy_var, d_phase_change_integrator->getCurrentContext());
        d_total_enthalpy_idx = d_var_db->registerClonedPatchDataIndex(d_enthalpy_var, d_h_idx);
        d_x_mom_idx = d_var_db->registerClonedPatchDataIndex(d_enthalpy_var, d_h_idx);
        d_y_mom_idx = d_var_db->registerClonedPatchDataIndex(d_enthalpy_var, d_h_idx);

        // Interpolating side-centered velocity and density to cell-centered
        d_u_sc_var = d_time_integrator->getVelocityVariable();
        d_u_sc_idx = d_var_db->mapVariableAndContextToIndex(d_u_sc_var, d_time_integrator->getCurrentContext());

        d_u_cc_var = new CellVariable<NDIM, double>("U_cc", NDIM);
        d_u_cc_idx = d_var_db->registerVariableAndContext(d_u_cc_var, d_var_db->getContext("U_cc"), 0);

        d_rho_cc_interp_var = new CellVariable<NDIM, double>("rho_cc_interp", NDIM);
        d_rho_cc_interp_idx =
            d_var_db->registerVariableAndContext(d_rho_cc_interp_var, d_var_db->getContext("rho_cc_interp"), 0);

        d_coarsest_ln = 0;
        d_finest_ln = d_patch_hierarchy->getFinestLevelNumber();
        PhaseChangeExamples::allocate_diagnostic_data(
            d_patch_hierarchy,
            { d_total_enthalpy_idx, d_u_cc_idx, d_rho_cc_interp_idx, d_x_mom_idx, d_y_mom_idx },
            d_coarsest_ln,
            d_finest_ln,
            loop_time);
        d_hier_cc_data_ops = new HierarchyCellDataOpsReal<NDIM, double>(d_patch_hierarchy, d_coarsest_ln, d_finest_ln);
        d_hier_sc_data_ops = new HierarchySideDataOpsReal<NDIM, double>(d_patch_hierarchy, d_coarsest_ln, d_finest_ln);

        if (SAMRAI_MPI::getRank() == 0)
        {
            const std::string file_name = d_input_db->getString("FILE_NAME");
            PhaseChangeExamples::open_diagnostic_file(d_vol_stream, file_name);
        }

        const bool is_from_restart = RestartManager::getManager()->isFromRestart();
        if (!is_from_restart)
        {
            // Target the gas volume in the Newton's iterations
            std::vector<double> H_integrals =
                LevelSetUtilities::computeHeavisideIntegrals2PhaseFlows(d_level_set_fixer->getLevelSetContainer());
            d_level_set_fixer->setInitialVolume(H_integrals[0]);

            // Save the initial volume of gas, liquid, and gas+solid phases, and the Lagrange multiplier in the stream.
            if (SAMRAI_MPI::getRank() == 0)
            {
                d_vol_stream << 0.0 << "\t" << H_integrals[0] << "\t" << H_integrals[1] << "\t" << 0.0 << "\t" << 0.0
                             << std::endl;
            }
        }
    }
    void postprocessStep(const double loop_time) override
    {
        HierarchyMathOps hier_math_ops("HierarchyMathOps", d_patch_hierarchy, d_coarsest_ln, d_finest_ln);
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        std::vector<double> H_integrals =
            LevelSetUtilities::computeHeavisideIntegrals2PhaseFlows(d_level_set_fixer->getLevelSetContainer());

        PhaseChangeExamples::allocate_diagnostic_data(
            d_patch_hierarchy,
            { d_total_enthalpy_idx, d_u_cc_idx, d_rho_cc_interp_idx, d_x_mom_idx, d_y_mom_idx },
            d_coarsest_ln,
            d_finest_ln,
            loop_time);

        hier_math_ops.interp(d_u_cc_idx, d_u_cc_var, d_u_sc_idx, d_u_sc_var, nullptr, loop_time, true);
        hier_math_ops.interp(
            d_rho_cc_interp_idx, d_rho_cc_interp_var, d_rho_sc_idx, d_side_density_var, nullptr, loop_time, true);
        double mass = 0;
        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = d_patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();

                Pointer<CellData<NDIM, double>> U_cc_data = patch->getPatchData(d_u_cc_idx);
                Pointer<CellData<NDIM, double>> rho_cc_interp_data = patch->getPatchData(d_rho_cc_interp_idx);
                Pointer<CellData<NDIM, double>> x_mom_data = patch->getPatchData(d_x_mom_idx);
                Pointer<CellData<NDIM, double>> y_mom_data = patch->getPatchData(d_y_mom_idx);
                Pointer<CellData<NDIM, double>> wgt_cc_data = patch->getPatchData(wgt_cc_idx);
                for (Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    CellIndex<NDIM> ci(it());

                    (*x_mom_data)(ci) = (*rho_cc_interp_data)(ci, 0) * (*U_cc_data)(ci, 0);
                    (*y_mom_data)(ci) = (*rho_cc_interp_data)(ci, 1) * (*U_cc_data)(ci, 1);
                    mass += (*rho_cc_interp_data)(ci, 0) * (*wgt_cc_data)(ci);
                }
            }
        }
        IBTK_MPI::sumReduction(&mass);

        d_hier_cc_data_ops->multiply(d_total_enthalpy_idx, d_rho_cc_idx, d_h_idx);
        const double total_enthalpy = d_hier_cc_data_ops->integral(d_total_enthalpy_idx, wgt_cc_idx);
        const double mass_sc = d_hier_sc_data_ops->integral(d_rho_sc_idx, wgt_sc_idx);
        const double x_mom = d_hier_cc_data_ops->integral(d_x_mom_idx, wgt_cc_idx);
        const double y_mom = d_hier_cc_data_ops->integral(d_y_mom_idx, wgt_cc_idx);
        if (SAMRAI_MPI::getRank() == 0)
        {
            d_vol_stream << loop_time << "\t" << H_integrals[0] << "\t" << H_integrals[1] << "\t"
                         << d_level_set_fixer->getLagrangeMultiplier() << "\t" << mass << "\t" << mass_sc << "\t"
                         << x_mom << "\t" << y_mom << "\t" << total_enthalpy << std::endl;
        }
    }
    void finalizeDiagnostics() override
    {
        PhaseChangeExamples::deallocate_diagnostic_data(
            d_patch_hierarchy,
            { d_total_enthalpy_idx, d_u_cc_idx, d_rho_cc_interp_idx, d_x_mom_idx, d_y_mom_idx },
            d_coarsest_ln,
            d_finest_ln);

        if (SAMRAI_MPI::getRank() == 0)
        {
            d_vol_stream.close();
        }
    }
    std::unique_ptr<LevelSetUtilities::LevelSetMassLossFixer> d_level_set_fixer;
    VariableDatabase<NDIM>* d_var_db = nullptr;
    int d_coarsest_ln = -1;
    int d_finest_ln = -1;
    Pointer<HierarchyCellDataOpsReal<NDIM, double>> d_hier_cc_data_ops;
    int d_rho_cc_idx = -1;
    int d_rho_sc_idx = -1;
    int d_h_idx = -1;
    int d_total_enthalpy_idx = -1;
    int d_x_mom_idx = -1;
    int d_y_mom_idx = -1;
    int d_u_sc_idx = -1;
    int d_u_cc_idx = -1;
    int d_rho_cc_interp_idx = -1;
    Pointer<SideVariable<NDIM, double>> d_u_sc_var;
    Pointer<CellVariable<NDIM, double>> d_u_cc_var;
    Pointer<CellVariable<NDIM, double>> d_rho_cc_interp_var;
    Pointer<HierarchySideDataOpsReal<NDIM, double>> d_hier_sc_data_ops;
    std::ofstream d_vol_stream;
};

std::unique_ptr<PhaseChangeExamples::CoupledApplication>
create_application(Pointer<AppInitializer> app_initializer)
{
    return std::make_unique<AdvectedBubble>(app_initializer);
}
} // namespace

int
main(int argc, char* argv[])
{
    return PhaseChangeExamples::run_coupled(argc, argv, &create_application);
}
