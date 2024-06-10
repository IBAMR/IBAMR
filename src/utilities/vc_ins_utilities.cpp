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
#include "ibamr/vc_ins_utilities.h"

#include "ibtk/PhysicalBoundaryUtilities.h"
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/ibtk_utilities.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{

namespace VCINSUtilities
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
} // callSetDensityCallBackFunction

SetFluidProperties::SetFluidProperties(const std::string& object_name,
                                       Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                       Pointer<CellVariable<NDIM, double> > ls_var,
                                       const double rho_liquid,
                                       const double rho_gas,
                                       const double mu_liquid,
                                       const double mu_gas,
                                       const double num_interface_cells)
    : d_object_name(object_name),
      d_adv_diff_solver(adv_diff_solver),
      d_ls_gas_var(ls_var),
      d_rho_liquid(rho_liquid),
      d_rho_gas(rho_gas),
      d_mu_liquid(mu_liquid),
      d_mu_gas(mu_gas),
      d_num_gas_interface_cells(num_interface_cells),
      d_num_phases(2)
{
    // Intentionally left blank
    return;
} // SetFluidProperties

SetFluidProperties::SetFluidProperties(const std::string& object_name,
                                       Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                       Pointer<CellVariable<NDIM, double> > ls_gas_var,
                                       Pointer<CellVariable<NDIM, double> > ls_solid_var,
                                       const double rho_liquid,
                                       const double rho_gas,
                                       const double rho_solid,
                                       const double mu_liquid,
                                       const double mu_gas,
                                       const double mu_solid,
                                       const double num_gas_interface_cells,
                                       const double num_solid_interface_cells,
                                       const bool set_mu_solid)
    : d_object_name(object_name),
      d_adv_diff_solver(adv_diff_solver),
      d_ls_gas_var(ls_gas_var),
      d_ls_solid_var(ls_solid_var),
      d_rho_liquid(rho_liquid),
      d_rho_gas(rho_gas),
      d_rho_solid(rho_solid),
      d_mu_liquid(mu_liquid),
      d_mu_gas(mu_gas),
      d_mu_solid(mu_solid),
      d_num_gas_interface_cells(num_gas_interface_cells),
      d_num_solid_interface_cells(num_solid_interface_cells),
      d_set_mu_solid(set_mu_solid),
      d_num_phases(3)
{
    // Intentionally left blank
    return;
} // SetFluidProperties

void
SetFluidProperties::setDensityPatchData2PhaseFlows(int rho_idx,
                                                   Pointer<Variable<NDIM> > rho_var,
                                                   SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
                                                   const int /*cycle_num*/,
                                                   const double time,
                                                   const double current_time,
                                                   const double new_time)
{
    // Get the current level set information
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int ls_idx = IBTK::invalid_index;
    if (IBTK::rel_equal_eps(time, current_time))
    {
        ls_idx = var_db->mapVariableAndContextToIndex(d_ls_gas_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (IBTK::rel_equal_eps(time, new_time))
    {
        ls_idx = var_db->mapVariableAndContextToIndex(d_ls_gas_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    // Set the density based on the level set
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // Normal way to set cell centered density
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
                const double* const patch_dx = patch_geom->getDx();
                double vol_cell = 1.0;
                for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
                const double alpha = d_num_gas_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(ls_idx);
                Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);

                for (Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    CellIndex<NDIM> ci(it());
                    const double phi = (*ls_data)(ci);
                    const double h_phi = IBTK::smooth_heaviside(phi, alpha);

                    (*rho_data)(ci) = d_rho_gas + (d_rho_liquid - d_rho_gas) * h_phi;
                }
            }
        }
    }

    Pointer<SideVariable<NDIM, double> > rho_sc_var = rho_var;
    if (rho_sc_var)
    {
        // Note, this method requires ghost cells to be filled for the level set variable
        RobinBcCoefStrategy<NDIM>* ls_bc_coef = d_adv_diff_solver->getPhysicalBcCoefs(d_ls_gas_var).front();
        IntVector<NDIM> cell_ghosts = 1;
        const int ls_scratch_idx = var_db->registerVariableAndContext(
            d_ls_gas_var, var_db->getContext(d_object_name + "::SCRATCH"), cell_ghosts);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(ls_scratch_idx, time);
        }
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        InterpolationTransactionComponent ls_transaction(ls_scratch_idx,
                                                         ls_idx,
                                                         "CONSERVATIVE_LINEAR_REFINE",
                                                         false,
                                                         "CONSERVATIVE_COARSEN",
                                                         "LINEAR",
                                                         false,
                                                         ls_bc_coef);
        Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
        hier_bdry_fill->initializeOperatorState(ls_transaction, patch_hierarchy);
        hier_bdry_fill->fillData(time);

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const double* const patch_dx = patch_geom->getDx();
                double vol_cell = 1.0;
                for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
                const double alpha = d_num_gas_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(ls_scratch_idx);
                Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);
// Various options to setting side-centered densities
#define SMOOTH_SC_RHO_2_PHASE 1
#define DESJARDINS_SC_RHO_2_PHASE 0
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                    {
                        SideIndex<NDIM> si(it(), axis, 0);
                        const double phi_lower = (*ls_data)(si.toCell(0));
                        const double phi_upper = (*ls_data)(si.toCell(1));
#if (DESJARDINS_SC_RHO_2_PHASE)
                        double h;
                        // Desjardins way to set side-centered density
                        if (phi_lower >= 0.0 && phi_upper >= 0.0)
                        {
                            h = 1.0;
                        }
                        else if (phi_lower < 0.0 && phi_upper < 0.0)
                        {
                            h = 0.0;
                        }
                        else
                        {
                            h = (std::max(phi_lower, 0.0) + std::max(phi_upper, 0.0)) /
                                (std::abs(phi_lower) + std::abs(phi_upper));
                        }
                        (*rho_data)(si) = d_rho_inside + (d_rho_outside - d_rho_inside) * h;
#endif
#if (SMOOTH_SC_RHO_2_PHASE)
                        // Simple average of phi onto side centers and set rho_sc directly
                        const double phi = 0.5 * (phi_lower + phi_upper);
                        const double h_phi = IBTK::smooth_heaviside(phi, alpha);

                        (*rho_data)(si) = d_rho_gas + (d_rho_liquid - d_rho_gas) * h_phi;
#endif
                    }
                }
            }
        }

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(ls_scratch_idx);
        }
        var_db->removePatchDataIndex(ls_scratch_idx);
    }

    return;
} // setDensityPatchData2PhaseFlows

void
SetFluidProperties::setDensityPatchData3PhaseFlows(int rho_idx,
                                                   Pointer<Variable<NDIM> > rho_var,
                                                   SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
                                                   const int /*cycle_num*/,
                                                   const double time,
                                                   const double current_time,
                                                   const double new_time)
{
    // Get the current level set information
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    int ls_solid_idx = IBTK::invalid_index;
    int ls_gas_idx = IBTK::invalid_index;
    if (IBTK::rel_equal_eps(time, current_time))
    {
        ls_solid_idx = var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (IBTK::rel_equal_eps(time, new_time))
    {
        ls_solid_idx = var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    if (IBTK::rel_equal_eps(time, current_time))
    {
        ls_gas_idx = var_db->mapVariableAndContextToIndex(d_ls_gas_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (IBTK::rel_equal_eps(time, new_time))
    {
        ls_gas_idx = var_db->mapVariableAndContextToIndex(d_ls_gas_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // Set the density based on the cell centered level set
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
                const double* const patch_dx = patch_geom->getDx();

                double vol_cell = 1.0;
                for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
                const double alpha = d_num_gas_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));
                const double beta = d_num_solid_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(ls_solid_idx);
                const Pointer<CellData<NDIM, double> > ls_gas_data = patch->getPatchData(ls_gas_idx);
                Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);

                // Li et al, 2015
                for (Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    CellIndex<NDIM> ci(it());
                    const double phi_s = (*ls_solid_data)(ci);
                    const double phi_g = (*ls_gas_data)(ci);

                    const double Hphi_s = IBTK::smooth_heaviside(phi_s, beta);
                    const double Hphi_g = IBTK::smooth_heaviside(phi_g, alpha);

                    // First, compute the density of the "flowing" phases
                    const double rho_flow = (d_rho_liquid - d_rho_gas) * Hphi_g + d_rho_gas;

                    // Next, set the density of the solid phase
                    (*rho_data)(ci) = (rho_flow - d_rho_solid) * Hphi_s + d_rho_solid;
                }
            }
        }
    }

    // Setting side centered density directly
    Pointer<SideVariable<NDIM, double> > rho_sc_var = rho_var;
    if (rho_sc_var)
    {
        // Note, this method requires ghost cells to be filled for the level set variable
        RobinBcCoefStrategy<NDIM>* ls_solid_bc_coef = d_adv_diff_solver->getPhysicalBcCoefs(d_ls_solid_var).front();
        RobinBcCoefStrategy<NDIM>* ls_gas_bc_coef = d_adv_diff_solver->getPhysicalBcCoefs(d_ls_gas_var).front();
        IntVector<NDIM> cell_ghosts = 1;
        const int ls_solid_scratch_idx = var_db->registerVariableAndContext(
            d_ls_solid_var, var_db->getContext(d_object_name + "::SOLID::SCRATCH"), cell_ghosts);
        const int ls_gas_scratch_idx = var_db->registerVariableAndContext(
            d_ls_gas_var, var_db->getContext(d_object_name + "::GAS::SCRATCH"), cell_ghosts);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(ls_solid_scratch_idx, time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(ls_gas_scratch_idx, time);
        }
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> ls_transaction_comps(2);
        ls_transaction_comps[0] = InterpolationTransactionComponent(ls_solid_scratch_idx,
                                                                    ls_solid_idx,
                                                                    "CONSERVATIVE_LINEAR_REFINE",
                                                                    false,
                                                                    "CONSERVATIVE_COARSEN",
                                                                    "LINEAR",
                                                                    false,
                                                                    ls_solid_bc_coef);
        ls_transaction_comps[1] = InterpolationTransactionComponent(ls_gas_scratch_idx,
                                                                    ls_gas_idx,
                                                                    "CONSERVATIVE_LINEAR_REFINE",
                                                                    false,
                                                                    "CONSERVATIVE_COARSEN",
                                                                    "LINEAR",
                                                                    false,
                                                                    ls_gas_bc_coef);
        Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
        hier_bdry_fill->initializeOperatorState(ls_transaction_comps, patch_hierarchy);
        hier_bdry_fill->fillData(time);

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();

                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(ls_solid_scratch_idx);
                const Pointer<CellData<NDIM, double> > ls_gas_data = patch->getPatchData(ls_gas_scratch_idx);
                Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);

                const double* const patch_dx = patch_geom->getDx();
                double vol_cell = 1.0;
                for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
                const double alpha = d_num_gas_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));
                const double beta = d_num_solid_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

// Various options to setting side-centered densities
#define SMOOTH_SC_RHO_3_PHASE 0
#define DESJARDINS_SC_RHO_3_PHASE 0
#define HARMONIC_CC_TO_SC_RHO_3_PHASE 1

                // Compute the indicators for both level sets
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                    {
                        SideIndex<NDIM> si(it(), axis, 0);
                        const double phi_solid_lower = (*ls_solid_data)(si.toCell(0));
                        const double phi_solid_upper = (*ls_solid_data)(si.toCell(1));
                        const double phi_gas_lower = (*ls_gas_data)(si.toCell(0));
                        const double phi_gas_upper = (*ls_gas_data)(si.toCell(1));

#if (DESJARDINS_SC_RHO_3_PHASE)
                        {
                            // SETTING 1: Desjardins way to set side-centered density
                            double h_solid, h_gas;
                            if (phi_solid_lower >= 0.0 && phi_solid_upper >= 0.0)
                            {
                                h_solid = 1.0;
                            }
                            else if (phi_solid_lower <= 0.0 && phi_solid_upper <= 0.0)
                            {
                                h_solid = 0.0;
                            }
                            else
                            {
                                h_solid = (std::max(phi_solid_lower, 0.0) + std::max(phi_solid_upper, 0.0)) /
                                          (std::abs(phi_solid_lower) + std::abs(phi_solid_upper));
                            }
                            if (phi_gas_lower >= 0.0 && phi_gas_upper >= 0.0)
                            {
                                h_gas = 1.0;
                            }
                            else if (phi_gas_lower <= 0.0 && phi_gas_upper <= 0.0)
                            {
                                h_gas = 0.0;
                            }
                            else
                            {
                                h_gas = (std::max(phi_gas_lower, 0.0) + std::max(phi_gas_upper, 0.0)) /
                                        (std::abs(phi_gas_lower) + std::abs(phi_gas_upper));
                            }

                            // First, compute the density of the "flowing" phases
                            const double rho_flow = d_rho_gas + (d_rho_liquid - d_rho_gas) * h_gas;

                            // Next, compute the density taking into account the solid phase
                            (*rho_data)(si) = d_rho_solid + (rho_flow - d_rho_solid) * h_solid;
#endif
#if (HARMONIC_CC_TO_SC_RHO_3_PHASE)
                            // SETTING 2: Set rho on cell centers and harmonic average to side centers
                            const double h_solid_lower = IBTK::smooth_heaviside(phi_solid_lower, beta);
                            const double h_solid_upper = IBTK::smooth_heaviside(phi_solid_upper, beta);

                            const double h_gas_lower = IBTK::smooth_heaviside(phi_gas_lower, alpha);
                            const double h_gas_upper = IBTK::smooth_heaviside(phi_gas_upper, alpha);

                            const double rho_flow_lower = (d_rho_liquid - d_rho_gas) * h_gas_lower + d_rho_gas;
                            const double rho_flow_upper = (d_rho_liquid - d_rho_gas) * h_gas_upper + d_rho_gas;

                            const double rho_full_lower = (rho_flow_lower - d_rho_solid) * h_solid_lower + d_rho_solid;
                            const double rho_full_upper = (rho_flow_upper - d_rho_solid) * h_solid_upper + d_rho_solid;

                            (*rho_data)(si) = 2.0 * rho_full_upper * rho_full_lower / (rho_full_upper + rho_full_lower);
#endif
#if (SMOOTH_SC_RHO_3_PHASE)
                            // SETTING 3: Simple average of phi onto side centers and set rho_sc directly
                            const double phi_solid = 0.5 * (phi_solid_lower + phi_solid_upper);
                            const double phi_gas = 0.5 * (phi_gas_lower + phi_gas_upper);

                            const double h_solid = IBTK::smooth_heaviside(phi_solid, alpha);
                            const double h_gas = IBTK::smooth_heaviside(phi_gas, beta);

                            const double rho_flow = (d_rho_liquid - d_rho_gas) * h_gas + d_rho_gas;
                            const double rho_full = (rho_flow - d_rho_solid) * h_solid + d_rho_solid;

                            (*rho_data)(si) = rho_full;
#endif
                        }
                    }
                }
            }

            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(ls_solid_scratch_idx);
                patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(ls_gas_scratch_idx);
            }
            var_db->removePatchDataIndex(ls_solid_scratch_idx);
            var_db->removePatchDataIndex(ls_gas_scratch_idx);
        }

        return;
    } // setDensityPatchData3PhaseFlows

    void SetFluidProperties::setViscosityPatchData2PhaseFlows(int mu_idx,
                                                              Pointer<Variable<NDIM> > /*mu_var*/,
                                                              SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
                                                              const int /*cycle_num*/,
                                                              const double time,
                                                              const double current_time,
                                                              const double new_time)
    {
        // Get the current level set information
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        int ls_idx = IBTK::invalid_index;
        if (IBTK::rel_equal_eps(time, current_time))
        {
            ls_idx = var_db->mapVariableAndContextToIndex(d_ls_gas_var, d_adv_diff_solver->getCurrentContext());
        }
        else if (IBTK::rel_equal_eps(time, new_time))
        {
            ls_idx = var_db->mapVariableAndContextToIndex(d_ls_gas_var, d_adv_diff_solver->getNewContext());
        }
        else
        {
            TBOX_ERROR("This statement should not be reached");
        }

        // Set the density based on the level set
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
                const double* const patch_dx = patch_geom->getDx();
                double vol_cell = 1.0;
                for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
                const double alpha = d_num_gas_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(ls_idx);
                if (!ls_data) TBOX_ERROR("This statement should not be reached");
                Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_idx);

                for (Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    CellIndex<NDIM> ci(it());
                    const double phi = (*ls_data)(ci);

                    double h_phi = IBTK::smooth_heaviside(phi, alpha);
                    (*mu_data)(ci) = d_mu_gas + (d_mu_liquid - d_mu_gas) * h_phi;
                }
            }
        }

        return;
    } // setViscosityPatchData2PhaseFlows

    void SetFluidProperties::setViscosityPatchData3PhaseFlows(int mu_idx,
                                                              Pointer<Variable<NDIM> > mu_var,
                                                              SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
                                                              const int /*cycle_num*/,
                                                              const double time,
                                                              const double current_time,
                                                              const double new_time)
    {
        // Get the current level set information
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        int ls_solid_idx = IBTK::invalid_index;
        int ls_gas_idx = IBTK::invalid_index;

        if (IBTK::rel_equal_eps(time, current_time))
        {
            ls_solid_idx = var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getCurrentContext());
        }
        else if (IBTK::rel_equal_eps(time, new_time))
        {
            ls_solid_idx = var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getNewContext());
        }
        else
        {
            TBOX_ERROR("This statement should not be reached");
        }

        if (IBTK::rel_equal_eps(time, current_time))
        {
            ls_gas_idx = var_db->mapVariableAndContextToIndex(d_ls_gas_var, d_adv_diff_solver->getCurrentContext());
        }
        else if (IBTK::rel_equal_eps(time, new_time))
        {
            ls_gas_idx = var_db->mapVariableAndContextToIndex(d_ls_gas_var, d_adv_diff_solver->getNewContext());
        }
        else
        {
            TBOX_ERROR("This statement should not be reached");
        }

        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();

        Pointer<CellVariable<NDIM, double> > mu_cc_var = mu_var;
        if (mu_cc_var)
        {
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                    const double* const patch_dx = patch_geom->getDx();
                    double vol_cell = 1.0;
                    for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
                    const double alpha =
                        d_num_gas_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));
                    const double beta =
                        d_num_solid_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

                    const Box<NDIM>& patch_box = patch->getBox();
                    const Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(ls_solid_idx);
                    const Pointer<CellData<NDIM, double> > ls_gas_data = patch->getPatchData(ls_gas_idx);
                    Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_idx);

                    for (Box<NDIM>::Iterator it(patch_box); it; it++)
                    {
                        CellIndex<NDIM> ci(it());
                        const double phi_s = (*ls_solid_data)(ci);
                        const double phi_g = (*ls_gas_data)(ci);
                        const double Hphi_s = IBTK::smooth_heaviside(phi_s, beta);
                        const double Hphi_g = IBTK::smooth_heaviside(phi_g, alpha);

                        // First, compute the viscosity of the "flowing" phases
                        const double mu_flow = (d_mu_liquid - d_mu_gas) * Hphi_g + d_mu_gas;

                        // Next, set the viscosity of the solid phase in the usual way
                        if (d_set_mu_solid)
                            (*mu_data)(ci) = (mu_flow - d_mu_solid) * Hphi_s + d_mu_solid;
                        else
                            (*mu_data)(ci) = mu_flow;
                    }
                }
            }
        }
        else
        {
            // Erroring out if any other centered is used for mu
            TBOX_ERROR("This statement should not have been reached");
        }

        return;
    } // setViscosityPatchData3PhaseFlows

    GravityForcing::GravityForcing(const std::string& object_name,
                                   Pointer<INSVCStaggeredHierarchyIntegrator> ins_hierarchy_integrator,
                                   std::vector<double> grav_const)
        : d_object_name(object_name), d_ins_hierarchy_integrator(ins_hierarchy_integrator), d_grav_const(grav_const)
    {
        d_grav_type = "FULL";
        return;
    } // GravityForcing

    GravityForcing::GravityForcing(const std::string& object_name,
                                   Pointer<AdvDiffHierarchyIntegrator> adv_diff_hierarchy_integrator,
                                   Pointer<CellVariable<NDIM, double> > ls_gas_var,
                                   Pointer<Database> input_db,
                                   std::vector<double> grav_const)
        : d_object_name(object_name),
          d_adv_diff_hierarchy_integrator(adv_diff_hierarchy_integrator),
          d_ls_gas_var(ls_gas_var),
          d_grav_const(grav_const)
    {
        d_grav_type = "FLOW";
        d_rho_neg = input_db->getDouble("rho_neg");
        d_rho_pos = input_db->getDouble("rho_pos");
        d_num_gas_interface_cells = input_db->getDouble("num_interface_cells");

        return;
    } // GravityForcing

    bool GravityForcing::isTimeDependent() const
    {
        return true;
    } // isTimeDependent

    void GravityForcing::setDataOnPatchHierarchy(const int data_idx,
                                                 Pointer<Variable<NDIM> > var,
                                                 Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                 const double data_time,
                                                 const bool initial_time,
                                                 const int coarsest_ln_in,
                                                 const int finest_ln_in)
    {
        const int coarsest_ln = (coarsest_ln_in == IBTK::invalid_level_number ? 0 : coarsest_ln_in);
        const int finest_ln =
            (finest_ln_in == IBTK::invalid_level_number ? hierarchy->getFinestLevelNumber() : finest_ln_in);
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

        if (d_grav_type == "FLOW")
        {
            // Get level set information
            int ls_gas_current_idx = var_db->mapVariableAndContextToIndex(
                d_ls_gas_var, d_adv_diff_hierarchy_integrator->getCurrentContext());
            int ls_gas_new_idx =
                var_db->mapVariableAndContextToIndex(d_ls_gas_var, d_adv_diff_hierarchy_integrator->getNewContext());
            const bool ls_gas_new_is_allocated = d_adv_diff_hierarchy_integrator->isAllocatedPatchData(ls_gas_new_idx);
            int ls_gas_idx = ls_gas_new_is_allocated ? ls_gas_new_idx : ls_gas_current_idx;

            IntVector<NDIM> cell_ghosts = 1;
            d_ls_gas_scratch_idx = var_db->registerVariableAndContext(
                d_ls_gas_var, var_db->getContext(d_object_name + "::LS_GAS_SCRATCH"), cell_ghosts);

#if !defined(NDEBUG)
            TBOX_ASSERT(ls_gas_idx >= 0);
            TBOX_ASSERT(d_ls_gas_scratch_idx >= 0);
#endif

            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                hierarchy->getPatchLevel(ln)->allocatePatchData(d_ls_gas_scratch_idx, data_time);
            }
            using InterpolationTransactionComponent =
                HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
            std::vector<InterpolationTransactionComponent> ls_transaction_comp(1);
            ls_transaction_comp[0] =
                InterpolationTransactionComponent(d_ls_gas_scratch_idx,
                                                  ls_gas_idx,
                                                  "CONSERVATIVE_LINEAR_REFINE",
                                                  false,
                                                  "CONSERVATIVE_COARSEN",
                                                  "LINEAR",
                                                  false,
                                                  d_adv_diff_hierarchy_integrator->getPhysicalBcCoefs(d_ls_gas_var));
            Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
            hier_bdry_fill->initializeOperatorState(ls_transaction_comp, hierarchy);
            hier_bdry_fill->fillData(data_time);
        }

        // Fill data on each patch level
        CartGridFunction::setDataOnPatchHierarchy(
            data_idx, var, hierarchy, data_time, initial_time, coarsest_ln_in, finest_ln_in);

        if (d_grav_type == "FLOW")
        {
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                hierarchy->getPatchLevel(ln)->deallocatePatchData(d_ls_gas_scratch_idx);
            }
            var_db->removePatchDataIndex(d_ls_gas_scratch_idx);
        }

        return;
    } // setDataOnPatchHierarchy

    void GravityForcing::setDataOnPatch(const int data_idx,
                                        Pointer<Variable<NDIM> > /*var*/,
                                        Pointer<Patch<NDIM> > patch,
                                        const double /*data_time*/,
                                        const bool initial_time,
                                        Pointer<PatchLevel<NDIM> > /*patch_level*/)
    {
        Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(data_idx);
        if (initial_time)
        {
            f_data->fillAll(0.0);
            return;
        }

        const Box<NDIM>& patch_box = patch->getBox();
        if (d_grav_type == "FULL")
        {
            // Get interpolated density variable
            const int rho_ins_idx = d_ins_hierarchy_integrator->getLinearOperatorRhoPatchDataIndex();

#if !defined(NDEBUG)
            TBOX_ASSERT(rho_ins_idx >= 0);
#endif

            const Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_ins_idx);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                {
                    SideIndex<NDIM> s_i(it(), axis, SideIndex<NDIM>::Lower);
                    (*f_data)(s_i) = ((*rho_data)(s_i)) * d_grav_const[axis];
                }
            }
        }
        else if (d_grav_type == "FLOW")
        {
            // Set the gravity force. In this version, the gravity force is reconstructed from the flow density field.
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            const Pointer<CellData<NDIM, double> > ls_gas_data = patch->getPatchData(d_ls_gas_scratch_idx);

            double alpha = 1.0;
            for (int d = 0; d < NDIM; ++d) alpha *= patch_dx[d];
            alpha = std::pow(alpha, 1.0 / NDIM) * d_num_gas_interface_cells;

            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                {
                    SideIndex<NDIM> s_i(it(), axis, SideIndex<NDIM>::Lower);

                    // Reconstruct density
                    double phi_gas_lower = (*ls_gas_data)(s_i.toCell(0));
                    double phi_gas_upper = (*ls_gas_data)(s_i.toCell(1));

                    double h_gas_lower, h_gas_upper;
                    h_gas_lower = IBTK::smooth_heaviside(phi_gas_lower, alpha);
                    h_gas_upper = IBTK::smooth_heaviside(phi_gas_upper, alpha);

                    const double rho_flow_lower = (d_rho_pos - d_rho_neg) * h_gas_lower + d_rho_neg;
                    const double rho_flow_upper = (d_rho_pos - d_rho_neg) * h_gas_upper + d_rho_neg;
                    (*f_data)(s_i) = d_grav_const[axis] * 2.0 * (rho_flow_lower * rho_flow_upper) /
                                     (rho_flow_lower + rho_flow_upper);
                }
            }
        }

        return;
    } // setDataOnPatch

    //////////////////////////////////////////////////////////////////////////////

} // namespace VCINSUtilities

} // namespace IBAMR
