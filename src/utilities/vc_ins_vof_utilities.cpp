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
#include <ibamr/vc_ins_vof_utilities.h>

#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/ibtk_utilities.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{

namespace VCINSVOFUtilities
{
static double
vof_alpha(const double& phi, const std::array<double, NDIM>& dPhi)
{
#if (NDIM == 2)
    const double ax = std::abs(dPhi[0]);
    const double ay = std::abs(dPhi[1]);

    const double phi_max = 0.5 * (ax + ay);
    const double phi_mid = 0.5 * std::abs(ax - ay);

    // Degenerate case: If no gradient in the cell -> phi constant in cell
    if (tbox::MathUtilities<double>::equalEps(phi_max, 0.0))
    {
        return IBTK::discontinuous_heaviside(phi);
    }

    // Denominator used in quadratic pieces: phi_max^2 - phi_mid^2 = |Dx||Dy|
    const double denom = std::pow(phi_max, 2) - std::pow(phi_mid, 2);

    // Clamp function to limit values between 0 and 1
    auto clamp = [](double x, double lower, double upper) { return (x < lower) ? lower : (x > upper ? upper : x); };

    // Handle vertical/horizontal cut (denom==0): only linear+saturation remains
    if (tbox::MathUtilities<double>::equalEps(denom, 0.0))
    {
        // Here phi_mid == phi_max and phi_max + phi_mid == gmax = max(|Dx|,|Dy|)
        double gmax = ax;
        if (ay > gmax)
        {
            gmax = ay;
        }

        // Exact linear middle branch with saturation (the only non-empty branch)
        return clamp(0.5 + phi / gmax, 0.0, 1.0);
    }

    // Eq. (28) of Pijl paper: piecewise branches
    double alpha;
    if (phi <= -phi_max)
    {
        alpha = 0.0;
    }
    else if (phi < -phi_mid)
    {
        alpha = 0.5 * (std::pow(phi_max + phi, 2) / denom);
    }
    else if (phi <= phi_mid)
    {
        alpha = 0.5 + phi / (phi_max + phi_mid);
    }
    else if (phi < phi_max)
    {
        alpha = 1.0 - 0.5 * (std::pow(phi_max - phi, 2) / denom);
    }
    else
    {
        alpha = 1.0;
    }

    // Guard against roundoff
    alpha = clamp(alpha, 0.0, 1.0);
    return alpha;

#else
    TBOX_ERROR("vof_alpha() not implemented for NDIM != 2");
    return 0.0;
#endif
}

VOFInitialConditionFromLevelSet::VOFInitialConditionFromLevelSet(const std::string& object_name,
                                                                 SAMRAI::tbox::Pointer<IBTK::CartGridFunction> ls_ic)
    : d_object_name(object_name), d_ls_ic(ls_ic)
{
    if (d_ls_ic.isNull())
    {
        TBOX_ERROR(d_object_name << ": level set initial condition pointer is null.\n");
    }

    // Register scratch storage used to ask the LS initial condition to write \phi on each patch.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_phi_scratch_var = new CellVariable<NDIM, double>(d_object_name + "::phi_scratch");
    const IntVector<NDIM> no_ghosts = 0;
    d_phi_scratch_idx = var_db->registerVariableAndContext(
        d_phi_scratch_var, var_db->getContext(d_object_name + "::PHI_SCRATCH"), no_ghosts);

    return;
}

bool
VOFInitialConditionFromLevelSet::isTimeDependent() const
{
    return d_ls_ic->isTimeDependent();
}

void
VOFInitialConditionFromLevelSet::setDataOnPatch(const int data_idx,
                                                Pointer<Variable<NDIM>> /*var*/,
                                                Pointer<Patch<NDIM>> patch,
                                                const double data_time,
                                                const bool initial_time,
                                                Pointer<PatchLevel<NDIM>> patch_level)
{
    if (!initial_time) return;

    // Ask the LS initial condition to fill a scratch cell-centered field \phi on this patch.
    patch->allocatePatchData(d_phi_scratch_idx, data_time);
    d_ls_ic->setDataOnPatch(d_phi_scratch_idx,
                            d_phi_scratch_var,
                            patch,
                            data_time,
                            /*initial_time=*/true,
                            patch_level);

    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CellData<NDIM, double>> phi_data = patch->getPatchData(d_phi_scratch_idx);
    Pointer<CellData<NDIM, double>> vof_data = patch->getPatchData(data_idx);

    auto get_shift = [](int dir, int shift)
    {
        SAMRAI::hier::Index<NDIM> iv(0);
        iv(dir) = shift;
        return iv;
    };

    // We only have interior \phi values here (no ghost fill). Use one-sided differences near
    // patch boundaries and centered differences in the interior.
    auto in_box = [&patch_box](const CellIndex<NDIM>& idx) -> bool
    {
        const hier::Index<NDIM>& lo = patch_box.lower();
        const hier::Index<NDIM>& hi = patch_box.upper();
        for (int d = 0; d < NDIM; ++d)
        {
            if (idx(d) < lo(d) || idx(d) > hi(d)) return false;
        }
        return true;
    };

    for (Box<NDIM>::Iterator it(patch_box); it; it++)
    {
        const CellIndex<NDIM> ci(it());
        const double phi = (*phi_data)(ci);

        std::array<double, NDIM> dPhi;
        for (int d = 0; d < NDIM; ++d)
        {
            const CellIndex<NDIM> ci_plus(ci + get_shift(d, 1));
            const CellIndex<NDIM> ci_minus(ci + get_shift(d, -1));

            const bool has_plus = in_box(ci_plus);
            const bool has_minus = in_box(ci_minus);

            if (has_plus && has_minus)
            {
                dPhi[d] = 0.5 * ((*phi_data)(ci_plus) - (*phi_data)(ci_minus));
            }
            else if (has_plus)
            {
                dPhi[d] = (*phi_data)(ci_plus)-phi;
            }
            else if (has_minus)
            {
                dPhi[d] = phi - (*phi_data)(ci_minus);
            }
            else
            {
                dPhi[d] = 0.0;
            }
        }

        (*vof_data)(ci) = vof_alpha(phi, dPhi);
    }

    patch->deallocatePatchData(d_phi_scratch_idx);
    return;
}
VOFFromLevelSetInitializer::VOFFromLevelSetInitializer(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_integrator,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> ls_var,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> vof_var)
    : d_object_name(object_name), d_integrator(adv_diff_integrator), d_ls_var(ls_var), d_vof_var(vof_var)
{
    if (d_integrator.isNull())
    {
        TBOX_ERROR(d_object_name << ": adv-diff integrator pointer is null.\n");
    }
    if (d_ls_var.isNull() || d_vof_var.isNull())
    {
        TBOX_ERROR(d_object_name << ": ls_var and/or vof_var pointer is null.\n");
    }
    return;
}

void
VOFFromLevelSetInitializer::registerIntegrateHierarchyCallback()
{
    d_integrator->registerIntegrateHierarchyCallback(&VOFFromLevelSetInitializer::integrateHierarchyCallback,
                                                     static_cast<void*>(this));
    return;
}

void
VOFFromLevelSetInitializer::computeVOFFromLevelSet(const double time, const bool use_new_context)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int ls_idx = var_db->mapVariableAndContextToIndex(
        d_ls_var, use_new_context ? d_integrator->getNewContext() : d_integrator->getCurrentContext());
    const int vof_idx = var_db->mapVariableAndContextToIndex(
        d_vof_var, use_new_context ? d_integrator->getNewContext() : d_integrator->getCurrentContext());
    computeVOFInternal(time, ls_idx, vof_idx);
    return;
}

void
VOFFromLevelSetInitializer::integrateHierarchyCallback(double /*current_time*/,
                                                       double new_time,
                                                       int /*cycle_num*/,
                                                       void* ctx)
{
    auto* self = static_cast<VOFFromLevelSetInitializer*>(ctx);
    self->computeVOFFromLevelSet(new_time, /*use_new_context=*/true);
    return;
}

void
VOFFromLevelSetInitializer::computeVOFInternal(const double time, const int ls_idx, const int vof_idx)
{
    Pointer<PatchHierarchy<NDIM>> patch_hier = d_integrator->getPatchHierarchy();
    const int finest_ln = patch_hier->getFinestLevelNumber();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Fill LS ghost cells into a temporary scratch index so we can take centered differences.
    RobinBcCoefStrategy<NDIM>* ls_bc_coef = d_integrator->getPhysicalBcCoefs(d_ls_var).front();
    const IntVector<NDIM> cell_ghosts = 1;
    const int ls_scratch_idx =
        var_db->registerVariableAndContext(d_ls_var, var_db->getContext("LSToVOF::SCRATCH"), cell_ghosts);
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        patch_hier->getPatchLevel(ln)->allocatePatchData(ls_scratch_idx, time);
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
    hier_bdry_fill->initializeOperatorState(ls_transaction, patch_hier);
    hier_bdry_fill->fillData(time);

    auto get_shift = [](int dir, int shift)
    {
        SAMRAI::hier::Index<NDIM> iv(0);
        iv(dir) = shift;
        return iv;
    };

    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> patch_level = patch_hier->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM, double>> ls_data = patch->getPatchData(ls_scratch_idx);
            Pointer<CellData<NDIM, double>> vof_data = patch->getPatchData(vof_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                const hier::Index<NDIM>& ci = it();

                std::array<double, NDIM> dPhi;
                for (int d = 0; d < NDIM; ++d)
                {
                    const CellIndex<NDIM> ci_plus(ci + get_shift(d, 1));
                    const CellIndex<NDIM> ci_minus(ci + get_shift(d, -1));
                    dPhi[d] = ((*ls_data)(ci_plus) - (*ls_data)(ci_minus)) / 2.0;
                }
                const double phi = (*ls_data)(ci);
                (*vof_data)(ci) = vof_alpha(phi, dPhi);
            }
        }
    }

    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        patch_hier->getPatchLevel(ln)->deallocatePatchData(ls_scratch_idx);
    }
    var_db->removePatchDataIndex(ls_scratch_idx);
    return;
}

void
callSetVOFBasedDensity(int rho_idx,
                       Pointer<Variable<NDIM>> rho_var,
                       Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                       const int cycle_num,
                       const double time,
                       const double current_time,
                       const double new_time,
                       void* ctx)
{
    // Set the density from the level set information
    auto ptr_SetFluidProperties = static_cast<SetVOFBasedFluidProperties*>(ctx);
    ptr_SetFluidProperties->setVOFBasedDensityPatchData(
        rho_idx, rho_var, hier_math_ops, cycle_num, time, current_time, new_time);

    return;
} // callSetVOFBasedDensity

void
callSetVOFBasedViscosity(int mu_idx,
                         Pointer<Variable<NDIM>> mu_var,
                         Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                         const int cycle_num,
                         const double time,
                         const double current_time,
                         const double new_time,
                         void* ctx)
{
    // Set the density from the level set information
    auto ptr_SetFluidProperties = static_cast<SetVOFBasedFluidProperties*>(ctx);
    ptr_SetFluidProperties->setVOFBasedViscosityPatchData(
        mu_idx, mu_var, hier_math_ops, cycle_num, time, current_time, new_time);

    return;
} // callSetVOFBasedViscosity

SetVOFBasedFluidProperties::SetVOFBasedFluidProperties(const std::string& object_name,
                                                       Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                                       Pointer<CellVariable<NDIM, double>> vof_var,
                                                       const double rho_liquid,
                                                       const double rho_gas,
                                                       const double mu_liquid,
                                                       const double mu_gas)
    : d_object_name(object_name),
      d_adv_diff_solver(adv_diff_solver),
      d_vof_var(vof_var),
      d_rho_liquid(rho_liquid),
      d_rho_gas(rho_gas),
      d_mu_liquid(mu_liquid),
      d_mu_gas(mu_gas),
      d_num_phases(2)
{
    // Intentionally left blank
    return;
} // SetFluidProperties

void
SetVOFBasedFluidProperties::setVOFBasedDensityPatchData2PhaseFlows(
    int rho_idx,
    Pointer<Variable<NDIM>> rho_var,
    SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
    const int /*cycle_num*/,
    const double time,
    const double current_time,
    const double new_time)
{
    // Get the current level set information
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int vof_idx = IBTK::invalid_index;
    if (IBTK::rel_equal_eps(time, current_time))
    {
        vof_idx = var_db->mapVariableAndContextToIndex(d_vof_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (IBTK::rel_equal_eps(time, new_time))
    {
        vof_idx = var_db->mapVariableAndContextToIndex(d_vof_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    // Set the density based on the level set
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // Normal way to set cell centered density
    Pointer<CellVariable<NDIM, double>> rho_cc_var = rho_var;
    if (rho_cc_var)
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CellData<NDIM, double>> vof_data = patch->getPatchData(vof_idx);
                Pointer<CellData<NDIM, double>> rho_data = patch->getPatchData(rho_idx);

                for (Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    CellIndex<NDIM> ci(it());
                    const double C = (*vof_data)(ci);

                    (*rho_data)(ci) = d_rho_gas + (d_rho_liquid - d_rho_gas) * C;
                }
            }
        }
    }

    Pointer<SideVariable<NDIM, double>> rho_sc_var = rho_var;
    if (rho_sc_var)
    {
        // Note, this method requires ghost cells to be filled for the vof variable
        RobinBcCoefStrategy<NDIM>* vof_bc_coef = d_adv_diff_solver->getPhysicalBcCoefs(d_vof_var).front();
        int vof_scratch_idx = var_db->mapVariableAndContextToIndex(d_vof_var, d_adv_diff_solver->getScratchContext());
        d_adv_diff_solver->allocatePatchData(vof_scratch_idx, time);

        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        InterpolationTransactionComponent vof_transaction(vof_scratch_idx,
                                                          vof_idx,
                                                          "CONSERVATIVE_LINEAR_REFINE",
                                                          false,
                                                          "CONSERVATIVE_COARSEN",
                                                          "LINEAR",
                                                          false,
                                                          vof_bc_coef);
        Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
        hier_bdry_fill->initializeOperatorState(vof_transaction, patch_hierarchy);
        hier_bdry_fill->fillData(time);

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();

                const Pointer<CellData<NDIM, double>> vof_data = patch->getPatchData(vof_scratch_idx);
                Pointer<SideData<NDIM, double>> rho_data = patch->getPatchData(rho_idx);

                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                    {
                        SideIndex<NDIM> si(it(), axis, 0);
                        const double C_lower = (*vof_data)(si.toCell(0));
                        const double C_upper = (*vof_data)(si.toCell(1));
                        const double C_mid = 0.5 * (C_lower + C_upper);

                        (*rho_data)(si) = d_rho_gas + (d_rho_liquid - d_rho_gas) * C_mid;
                    }
                }
            }
        }
        d_adv_diff_solver->deallocatePatchData(vof_scratch_idx);
    }

    return;
} // setDensityPatchData2PhaseFlows

void
SetVOFBasedFluidProperties::setVOFBasedViscosityPatchData2PhaseFlows(
    int mu_idx,
    Pointer<Variable<NDIM>> /*mu_var*/,
    SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
    const int /*cycle_num*/,
    const double time,
    const double current_time,
    const double new_time)
{
    // Get the VOF information
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int vof_idx = IBTK::invalid_index;
    if (IBTK::rel_equal_eps(time, current_time))
    {
        vof_idx = var_db->mapVariableAndContextToIndex(d_vof_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (IBTK::rel_equal_eps(time, new_time))
    {
        vof_idx = var_db->mapVariableAndContextToIndex(d_vof_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    // Set viscosity based on the VOF function
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CellData<NDIM, double>> vof_data = patch->getPatchData(vof_idx);
            Pointer<CellData<NDIM, double>> mu_data = patch->getPatchData(mu_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                const double C = (*vof_data)(ci);

                (*mu_data)(ci) = d_mu_gas + (d_mu_liquid - d_mu_gas) * C;
            }
        }
    }

    return;
} // setViscosityPatchData2PhaseFlows

VOFBasedGravityForcing::VOFBasedGravityForcing(const std::string& object_name,
                                               Pointer<INSVCStaggeredHierarchyIntegrator> ins_hierarchy_integrator,
                                               std::vector<double> grav_const)
    : d_object_name(object_name), d_ins_hierarchy_integrator(ins_hierarchy_integrator), d_grav_const(grav_const)
{
    d_grav_type = "FULL";
    return;
} // GravityForcing

VOFBasedGravityForcing::VOFBasedGravityForcing(const std::string& object_name,
                                               Pointer<AdvDiffHierarchyIntegrator> adv_diff_hierarchy_integrator,
                                               Pointer<CellVariable<NDIM, double>> vof_var,
                                               Pointer<Database> input_db,
                                               std::vector<double> grav_const)
    : d_object_name(object_name),
      d_adv_diff_hierarchy_integrator(adv_diff_hierarchy_integrator),
      d_vof_var(vof_var),
      d_grav_const(grav_const)
{
    d_grav_type = "FLOW";
    d_rho_neg = input_db->getDouble("rho_neg");
    d_rho_pos = input_db->getDouble("rho_pos");

    return;
} // GravityForcing

bool
VOFBasedGravityForcing::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
VOFBasedGravityForcing::setDataOnPatchHierarchy(const int data_idx,
                                                Pointer<Variable<NDIM>> var,
                                                Pointer<PatchHierarchy<NDIM>> hierarchy,
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
        int vof_current_idx =
            var_db->mapVariableAndContextToIndex(d_vof_var, d_adv_diff_hierarchy_integrator->getCurrentContext());
        int vof_new_idx =
            var_db->mapVariableAndContextToIndex(d_vof_var, d_adv_diff_hierarchy_integrator->getNewContext());
        const bool vof_new_is_allocated = d_adv_diff_hierarchy_integrator->isAllocatedPatchData(vof_new_idx);
        int vof_idx = vof_new_is_allocated ? vof_new_idx : vof_current_idx;

        int d_vof_scratch_idx =
            var_db->mapVariableAndContextToIndex(d_vof_var, d_adv_diff_hierarchy_integrator->getScratchContext());
        d_adv_diff_hierarchy_integrator->allocatePatchData(d_vof_scratch_idx, data_time);
#if !defined(NDEBUG)
        TBOX_ASSERT(vof_idx >= 0);
        TBOX_ASSERT(d_vof_scratch_idx >= 0);
#endif

        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> vof_transaction_comp(1);
        vof_transaction_comp[0] =
            InterpolationTransactionComponent(d_vof_scratch_idx,
                                              vof_idx,
                                              "CONSERVATIVE_LINEAR_REFINE",
                                              false,
                                              "CONSERVATIVE_COARSEN",
                                              "LINEAR",
                                              false,
                                              d_adv_diff_hierarchy_integrator->getPhysicalBcCoefs(d_vof_var));
        Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
        hier_bdry_fill->initializeOperatorState(vof_transaction_comp, hierarchy);
        hier_bdry_fill->fillData(data_time);
    }

    // Fill data on each patch level
    CartGridFunction::setDataOnPatchHierarchy(
        data_idx, var, hierarchy, data_time, initial_time, coarsest_ln_in, finest_ln_in);

    if (d_grav_type == "FLOW")
    {
        d_adv_diff_hierarchy_integrator->deallocatePatchData(d_vof_scratch_idx);
    }

    return;
} // setDataOnPatchHierarchy

void
VOFBasedGravityForcing::setDataOnPatch(const int data_idx,
                                       Pointer<Variable<NDIM>> /*var*/,
                                       Pointer<Patch<NDIM>> patch,
                                       const double /*data_time*/,
                                       const bool initial_time,
                                       Pointer<PatchLevel<NDIM>> /*patch_level*/)
{
    Pointer<SideData<NDIM, double>> f_data = patch->getPatchData(data_idx);
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

        const Pointer<SideData<NDIM, double>> rho_data = patch->getPatchData(rho_ins_idx);
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
        Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
        const Pointer<CellData<NDIM, double>> vof_data = patch->getPatchData(d_vof_scratch_idx);

        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
            {
                SideIndex<NDIM> s_i(it(), axis, SideIndex<NDIM>::Lower);

                // Reconstruct density
                double C_lower = (*vof_data)(s_i.toCell(0));
                double C_upper = (*vof_data)(s_i.toCell(1));

                const double rho_flow_lower = (d_rho_pos - d_rho_neg) * C_lower + d_rho_neg;
                const double rho_flow_upper = (d_rho_pos - d_rho_neg) * C_upper + d_rho_neg;
                (*f_data)(s_i) =
                    d_grav_const[axis] * 2.0 * (rho_flow_lower * rho_flow_upper) / (rho_flow_lower + rho_flow_upper);
            }
        }
    }

    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////

} // namespace VCINSVOFUtilities

} // namespace IBAMR
