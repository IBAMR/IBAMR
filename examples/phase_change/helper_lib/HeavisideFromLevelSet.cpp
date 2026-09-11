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

#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <VariableDatabase.h>

#include <cmath>

#include "HeavisideFromLevelSet.h"

#include <ibamr/app_namespaces.h>

namespace PhaseChangeExamples
{
HeavisideFromLevelSet::HeavisideFromLevelSet(Pointer<AdvDiffHierarchyIntegrator> integrator,
                                             Pointer<CellVariable<NDIM, double>> ls_var,
                                             const double num_interface_cells)
    : d_integrator(integrator), d_ls_var(ls_var), d_num_interface_cells(num_interface_cells)
{
}

void
HeavisideFromLevelSet::synchronize_levelset_with_heaviside_fcn(int H_current_idx,
                                                               Pointer<HierarchyMathOps> hier_math_ops,
                                                               int /*integrator_step*/,
                                                               double /*time*/,
                                                               bool /*initial_time*/,
                                                               bool /*regrid_time*/,
                                                               void* ctx)
{
    const auto& self = *static_cast<HeavisideFromLevelSet*>(ctx);
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int ls_current_idx =
        var_db->mapVariableAndContextToIndex(self.d_ls_var, self.d_integrator->getCurrentContext());

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const double* patch_dx = patch_geom->getDx();
            double vol_cell = 1.0;
            for (int d = 0; d < NDIM; ++d)
            {
                vol_cell *= patch_dx[d];
            }
            const double num_interface_cells = self.d_num_interface_cells;
            const double alpha = num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

            Pointer<CellData<NDIM, double>> H_data = patch->getPatchData(H_current_idx);
            Pointer<CellData<NDIM, double>> ls_data = patch->getPatchData(ls_current_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double phi = (*ls_data)(ci);

                (*H_data)(ci) = IBTK::smooth_heaviside(phi, alpha);
            }
        }
    }
    return;
}
} // namespace PhaseChangeExamples
