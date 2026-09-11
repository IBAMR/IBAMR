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

#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/ibtk_utilities.h>

#include <CapillaryForces.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <VariableDatabase.h>

#include <cmath>

#include "LiquidFractionForceMask.h"

#include <ibamr/app_namespaces.h>

namespace PhaseChangeExamples
{
LiquidFractionForceMask::LiquidFractionForceMask(Pointer<CellVariable<NDIM, double>> liquid_fraction,
                                                 RobinBcCoefStrategy<NDIM>* bc_coef,
                                                 Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator,
                                                 Pointer<INSVCStaggeredHierarchyIntegrator> ins_integrator,
                                                 const double rho_liquid,
                                                 const double rho_gas)
    : d_liquid_fraction(liquid_fraction),
      d_bc_coef(bc_coef),
      d_adv_diff_integrator(adv_diff_integrator),
      d_ins_integrator(ins_integrator),
      d_rho_liquid(rho_liquid),
      d_rho_gas(rho_gas)
{
}

void
LiquidFractionForceMask::mask_surface_tension_force(int F_idx,
                                                    Pointer<HierarchyMathOps> hier_math_ops,
                                                    int /*integrator_step*/,
                                                    double time,
                                                    double /*current_time*/,
                                                    double /*new_time*/,
                                                    void* ctx)
{
    const auto& self = *static_cast<LiquidFractionForceMask*>(ctx);
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = hier_math_ops->getPatchHierarchy();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int lf_new_idx =
        var_db->mapVariableAndContextToIndex(self.d_liquid_fraction, self.d_adv_diff_integrator->getNewContext());
    const int lf_scratch_idx =
        var_db->mapVariableAndContextToIndex(self.d_liquid_fraction, self.d_adv_diff_integrator->getScratchContext());

    // ghost cell filling for liquid fraction variable.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> lf_transaction_comps(1);
    lf_transaction_comps[0] = InterpolationTransactionComponent(lf_scratch_idx,
                                                                lf_new_idx,
                                                                "CONSERVATIVE_LINEAR_REFINE",
                                                                false,
                                                                "CONSERVATIVE_COARSEN",
                                                                "LINEAR",
                                                                false,
                                                                self.d_bc_coef);

    Pointer<HierarchyGhostCellInterpolation> lf_hier_bdry_fill = new HierarchyGhostCellInterpolation();
    lf_hier_bdry_fill->initializeOperatorState(lf_transaction_comps, patch_hierarchy);
    lf_hier_bdry_fill->fillData(time);

    int rho_idx = self.d_ins_integrator->getLinearOperatorRhoPatchDataIndex();

    MultiphaseExamples::mask_force(patch_hierarchy, F_idx, rho_idx, lf_scratch_idx, self.d_rho_liquid, self.d_rho_gas);
}
} // namespace PhaseChangeExamples
