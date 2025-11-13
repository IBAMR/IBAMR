// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2024 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////
#include "ibamr/PhaseChangeUtilities.h"

#include <ibtk/HierarchyMathOps.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{

namespace PhaseChangeUtilities
{
void
callSetDensityCallbackFunction(int rho_idx,
                               Pointer<Variable<NDIM> > rho_var,
                               Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                               const int cycle_num,
                               const double time,
                               const double current_time,
                               const double new_time,
                               void* ctx)
{
    // Set the density from the level set information
    auto ptr_SetFluidProperties = static_cast<SetFluidProperties*>(ctx);
    ptr_SetFluidProperties->setDensityPatchData(
        rho_idx, rho_var, hier_math_ops, cycle_num, time, current_time, new_time);

    return;
} // callSetDensityCallBackFunction

void
callSetThermalConductivityCallbackFunction(int kappa_idx,
                                           Pointer<Variable<NDIM> > kappa_var,
                                           Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                           const int cycle_num,
                                           const double time,
                                           const double current_time,
                                           const double new_time,
                                           void* ctx)
{
    // Set the density from the level set information
    auto ptr_SetFluidProperties = static_cast<SetFluidProperties*>(ctx);
    ptr_SetFluidProperties->setThermalConductivityPatchData(
        kappa_idx, kappa_var, hier_math_ops, cycle_num, time, current_time, new_time);

    return;
} // callSetThermalConductivityCallbackFunction

void
callSetSpecificHeatCallbackFunction(int specific_heat_idx,
                                    Pointer<Variable<NDIM> > specific_heat_var,
                                    Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                    const int cycle_num,
                                    const double time,
                                    const double current_time,
                                    const double new_time,
                                    void* ctx)
{
    // Set the density from the level set information
    auto ptr_SetFluidProperties = static_cast<SetFluidProperties*>(ctx);
    ptr_SetFluidProperties->setSpecificHeatPatchData(
        specific_heat_idx, specific_heat_var, hier_math_ops, cycle_num, time, current_time, new_time);

    return;
} // callSetSpecificHeatCallbackFunction

void
callSetViscosityCallbackFunction(int mu_idx,
                                 Pointer<Variable<NDIM> > mu_var,
                                 Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                 const int cycle_num,
                                 const double time,
                                 const double current_time,
                                 const double new_time,
                                 void* ctx)
{
    // Set the density from the level set information
    auto ptr_SetFluidProperties = static_cast<SetFluidProperties*>(ctx);
    ptr_SetFluidProperties->setViscosityPatchData(
        mu_idx, mu_var, hier_math_ops, cycle_num, time, current_time, new_time);

    return;
} // callSetViscosityCallBackFunction

void
callTagLiquidFractionCellsCallbackFunction(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                           const int level_number,
                                           const double error_data_time,
                                           const int tag_index,
                                           const bool initial_time,
                                           const bool uses_richardson_extrapolation_too,
                                           void* ctx)
{
    // Set the density from the level set information
    auto ptr_tagLiquidFractionCells = static_cast<TagLiquidFractionRefinementCells*>(ctx);
    ptr_tagLiquidFractionCells->tagLiquidFractionCells(
        hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too, ctx);

    return;
} // callTagLiquidFractionCellsCallbackFunction

SetFluidProperties::SetFluidProperties(const std::string& object_name,
                                       Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                       const Pointer<CellVariable<NDIM, double> > H_var,
                                       RobinBcCoefStrategy<NDIM>* H_bc_coef,
                                       const Pointer<CellVariable<NDIM, double> > lf_var,
                                       RobinBcCoefStrategy<NDIM>* lf_bc_coef,
                                       const double rho_liquid,
                                       const double rho_solid,
                                       const double rho_gas,
                                       const double kappa_liquid,
                                       const double kappa_solid,
                                       const double kappa_gas,
                                       const double specific_heat_liquid,
                                       const double specific_heat_solid,
                                       const double specific_heat_gas,
                                       const double mu_liquid,
                                       const double mu_solid,
                                       const double mu_gas)
    : d_object_name(object_name),
      d_adv_diff_solver(adv_diff_solver),
      d_H_var(H_var),
      d_H_bc_coef(H_bc_coef),
      d_lf_var(lf_var),
      d_lf_bc_coef(lf_bc_coef),
      d_rho_liquid(rho_liquid),
      d_rho_solid(rho_solid),
      d_rho_gas(rho_gas),
      d_kappa_liquid(kappa_liquid),
      d_kappa_solid(kappa_solid),
      d_kappa_gas(kappa_gas),
      d_specific_heat_liquid(specific_heat_liquid),
      d_specific_heat_solid(specific_heat_solid),
      d_specific_heat_gas(specific_heat_gas),
      d_mu_liquid(mu_liquid),
      d_mu_solid(mu_solid),
      d_mu_gas(mu_gas)
{
    // Intentionally left blank
    return;
} // SetFluidProperties

SetFluidProperties::SetFluidProperties(const std::string& object_name,
                                       Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                       const Pointer<CellVariable<NDIM, double> > H_var,
                                       RobinBcCoefStrategy<NDIM>* H_bc_coef,
                                       const Pointer<CellVariable<NDIM, double> > lf_var,
                                       RobinBcCoefStrategy<NDIM>* lf_bc_coef,
                                       const double rho_liquid,
                                       const double rho_solid,
                                       const double rho_gas,
                                       const double kappa_liquid,
                                       const double kappa_solid,
                                       const double kappa_gas,
                                       const double specific_heat_liquid,
                                       const double specific_heat_solid,
                                       const double specific_heat_gas)
    : d_object_name(object_name),
      d_adv_diff_solver(adv_diff_solver),
      d_H_var(H_var),
      d_H_bc_coef(H_bc_coef),
      d_lf_var(lf_var),
      d_lf_bc_coef(lf_bc_coef),
      d_rho_liquid(rho_liquid),
      d_rho_solid(rho_solid),
      d_rho_gas(rho_gas),
      d_kappa_liquid(kappa_liquid),
      d_kappa_solid(kappa_solid),
      d_kappa_gas(kappa_gas),
      d_specific_heat_liquid(specific_heat_liquid),
      d_specific_heat_solid(specific_heat_solid),
      d_specific_heat_gas(specific_heat_gas)
{
    // Intentionally left blank
    return;
} // SetFluidProperties

void
SetFluidProperties::setDensityPatchData(int rho_idx,
                                        Pointer<Variable<NDIM> > rho_var,
                                        SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
                                        const int /*cycle_num*/,
                                        const double time,
                                        const double current_time,
                                        const double new_time)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int H_idx = IBTK::invalid_index;
    if (IBTK::rel_equal_eps(time, current_time))
    {
        H_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (IBTK::rel_equal_eps(time, new_time))
    {
        H_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    int lf_idx = IBTK::invalid_index;
    if (IBTK::rel_equal_eps(time, current_time))
    {
        lf_idx = var_db->mapVariableAndContextToIndex(d_lf_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (IBTK::rel_equal_eps(time, new_time))
    {
        lf_idx = var_db->mapVariableAndContextToIndex(d_lf_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    Pointer<CellVariable<NDIM, double> > rho_cc_var = rho_var;
    if (rho_cc_var)
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_idx);
                const Pointer<CellData<NDIM, double> > lf_data = patch->getPatchData(lf_idx);
                Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);

                for (Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    CellIndex<NDIM> ci(it());
                    const double heaviside = (*H_data)(ci);
                    const double liquid_fraction = (*lf_data)(ci);

                    (*rho_data)(ci) = d_rho_gas + (d_rho_solid - d_rho_gas) * heaviside +
                                      (d_rho_liquid - d_rho_solid) * liquid_fraction * heaviside;
                }
            }
        }
    }

    Pointer<SideVariable<NDIM, double> > rho_sc_var = rho_var;
    if (rho_sc_var)
    {
        IntVector<NDIM> cell_ghosts = 1;
        const int H_scratch_idx =
            var_db->registerVariableAndContext(d_H_var, var_db->getContext(d_object_name + "::SCRATCH"), cell_ghosts);
        const int lf_scratch_idx =
            var_db->registerVariableAndContext(d_lf_var, var_db->getContext(d_object_name + "::SCRATCH"), cell_ghosts);

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(H_scratch_idx, time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(lf_scratch_idx, time);
        }
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> transaction_comps(2);
        transaction_comps[0] = InterpolationTransactionComponent(H_scratch_idx,
                                                                 H_idx,
                                                                 "CONSERVATIVE_LINEAR_REFINE",
                                                                 false,
                                                                 "CONSERVATIVE_COARSEN",
                                                                 "LINEAR",
                                                                 false,
                                                                 d_H_bc_coef);

        transaction_comps[1] = InterpolationTransactionComponent(lf_scratch_idx,
                                                                 lf_idx,
                                                                 "CONSERVATIVE_LINEAR_REFINE",
                                                                 false,
                                                                 "CONSERVATIVE_COARSEN",
                                                                 "LINEAR",
                                                                 false,
                                                                 d_lf_bc_coef);

        Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
        hier_bdry_fill->initializeOperatorState(transaction_comps, patch_hierarchy);
        hier_bdry_fill->fillData(time);

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_scratch_idx);
                const Pointer<CellData<NDIM, double> > lf_data = patch->getPatchData(lf_scratch_idx);
                Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                    {
                        SideIndex<NDIM> si(it(), axis, 0);
                        const double heaviside_lower = (*H_data)(si.toCell(0));
                        const double heaviside_upper = (*H_data)(si.toCell(1));
                        const double heaviside = 0.5 * (heaviside_lower + heaviside_upper);

                        const double liquid_fraction_lower = (*lf_data)(si.toCell(0));
                        const double liquid_fraction_upper = (*lf_data)(si.toCell(1));
                        const double liquid_fraction = 0.5 * (liquid_fraction_lower + liquid_fraction_upper);

                        (*rho_data)(si) = d_rho_gas + (d_rho_solid - d_rho_gas) * heaviside +
                                          (d_rho_liquid - d_rho_solid) * liquid_fraction * heaviside;
                    }
                }
            }
        }

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(H_scratch_idx);
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(lf_scratch_idx);
        }
        var_db->removePatchDataIndex(H_scratch_idx);
        var_db->removePatchDataIndex(lf_scratch_idx);
    }

    return;
} // setDensityPatchData

void
SetFluidProperties::setThermalConductivityPatchData(int kappa_idx,
                                                    Pointer<Variable<NDIM> > /*kappa_var*/,
                                                    SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
                                                    const int /*cycle_num*/,
                                                    const double time,
                                                    const double current_time,
                                                    const double new_time)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int H_idx = IBTK::invalid_index;
    if (IBTK::rel_equal_eps(time, current_time))
    {
        H_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (IBTK::rel_equal_eps(time, new_time))
    {
        H_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    int lf_idx = IBTK::invalid_index;
    if (IBTK::rel_equal_eps(time, current_time))
    {
        lf_idx = var_db->mapVariableAndContextToIndex(d_lf_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (IBTK::rel_equal_eps(time, new_time))
    {
        lf_idx = var_db->mapVariableAndContextToIndex(d_lf_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_idx);
            const Pointer<CellData<NDIM, double> > lf_data = patch->getPatchData(lf_idx);
            Pointer<CellData<NDIM, double> > kappa_data = patch->getPatchData(kappa_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                const double heaviside = (*H_data)(ci);
                const double liquid_fraction = (*lf_data)(ci);
                (*kappa_data)(ci) = d_kappa_gas + (d_kappa_solid - d_kappa_gas) * heaviside +
                                    (d_kappa_liquid - d_kappa_solid) * liquid_fraction * heaviside;
            }
        }
    }

    return;
} // setThermalConductivityPatchData

void
SetFluidProperties::setSpecificHeatPatchData(int specific_heat_idx,
                                             Pointer<Variable<NDIM> > /*specific_heat_var*/,
                                             SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
                                             const int /*cycle_num*/,
                                             const double time,
                                             const double current_time,
                                             const double new_time)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int H_idx = IBTK::invalid_index;
    if (IBTK::rel_equal_eps(time, current_time))
    {
        H_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (IBTK::rel_equal_eps(time, new_time))
    {
        H_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    int lf_idx = IBTK::invalid_index;
    if (IBTK::rel_equal_eps(time, current_time))
    {
        lf_idx = var_db->mapVariableAndContextToIndex(d_lf_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (IBTK::rel_equal_eps(time, new_time))
    {
        lf_idx = var_db->mapVariableAndContextToIndex(d_lf_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_idx);
            const Pointer<CellData<NDIM, double> > lf_data = patch->getPatchData(lf_idx);
            Pointer<CellData<NDIM, double> > specific_heat_data = patch->getPatchData(specific_heat_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                const double heaviside = (*H_data)(ci);
                const double liquid_fraction = (*lf_data)(ci);
                (*specific_heat_data)(ci) =
                    d_specific_heat_gas + (d_specific_heat_solid - d_specific_heat_gas) * heaviside +
                    (d_specific_heat_liquid - d_specific_heat_solid) * liquid_fraction * heaviside;
            }
        }
    }

    return;
} // setSpecificHeatPatchData

void
SetFluidProperties::setViscosityPatchData(int mu_idx,
                                          Pointer<Variable<NDIM> > /*mu_var*/,
                                          SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
                                          const int /*cycle_num*/,
                                          const double time,
                                          const double current_time,
                                          const double new_time)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int H_idx = IBTK::invalid_index;
    if (IBTK::rel_equal_eps(time, current_time))
    {
        H_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (IBTK::rel_equal_eps(time, new_time))
    {
        H_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    int lf_idx = IBTK::invalid_index;
    if (IBTK::rel_equal_eps(time, current_time))
    {
        lf_idx = var_db->mapVariableAndContextToIndex(d_lf_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (IBTK::rel_equal_eps(time, new_time))
    {
        lf_idx = var_db->mapVariableAndContextToIndex(d_lf_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_idx);
            const Pointer<CellData<NDIM, double> > lf_data = patch->getPatchData(lf_idx);
            Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                const double heaviside = (*H_data)(ci);
                const double liquid_fraction = (*lf_data)(ci);
                (*mu_data)(ci) = d_mu_gas + (d_mu_solid - d_mu_gas) * heaviside +
                                 (d_mu_liquid - d_mu_solid) * liquid_fraction * heaviside;
            }
        }
    }

    return;
} // setViscosityPatchData

void
TagLiquidFractionRefinementCells::tagLiquidFractionCells(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                                         const int level_number,
                                                         const double /*error_data_time*/,
                                                         const int tag_index,
                                                         const bool initial_time,
                                                         const bool /*uses_richardson_extrapolation_too*/,
                                                         void* /*ctx*/)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT(level_number >= 0);
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
#endif

    // Get the liquid fraction information
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int lf_idx = var_db->mapVariableAndContextToIndex(d_lf_var, d_adv_diff_solver->getCurrentContext());
    const int lf_grad_idx = var_db->mapVariableAndContextToIndex(d_lf_grad_var, d_adv_diff_solver->getCurrentContext());

    // Tag cells based on the value of the level set variable
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, int> > tags_data = patch->getPatchData(tag_index);
        Pointer<CellData<NDIM, double> > lf_data = patch->getPatchData(lf_idx);
        Pointer<CellData<NDIM, double> > lf_grad_data = patch->getPatchData(lf_grad_idx);

        for (CellIterator<NDIM> ic(patch_box); ic; ic++)
        {
            const hier::Index<NDIM>& i = ic();
            const double liquid_fraction = (*lf_data)(i);

            if (initial_time)
            {
                if (liquid_fraction >= d_tag_min_value && liquid_fraction <= d_tag_max_value) (*tags_data)(i) = 1;
            }
            else
            {
                bool non_zero_gradient = false;
                for (int d = 0; d < NDIM; ++d)
                {
                    if (!IBTK::abs_equal_eps((*lf_grad_data)(i, d), 0.0))
                    {
                        non_zero_gradient = true;
                        break;
                    }
                }
                if (non_zero_gradient) (*tags_data)(i) = 1;
            }
        }
    }

    return;
} // tagLiquidFractionCells

//////////////////////////////////////////////////////////////////////////////

} // namespace PhaseChangeUtilities

} // namespace IBAMR
