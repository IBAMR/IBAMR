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

#include "CapillaryForces.h"

#include <ibamr/app_namespaces.h>

namespace MultiphaseExamples
{
void
mask_force(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
           const int F_idx,
           const int rho_idx,
           const int lf_scratch_idx,
           const double rho_liquid,
           const double rho_gas)
{
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<SideData<NDIM, double>> rho_data = patch->getPatchData(rho_idx);
            Pointer<SideData<NDIM, double>> F_data = patch->getPatchData(F_idx);
            Pointer<CellData<NDIM, double>> lf_data;
            if (lf_scratch_idx >= 0)
            {
                lf_data = patch->getPatchData(lf_scratch_idx);
            }
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                {
                    const SideIndex<NDIM> si(it(), axis, SideIndex<NDIM>::Lower);
                    double multiplier_term = 2.0 * (*rho_data)(si) / (rho_liquid + rho_gas);
                    if (lf_data)
                    {
                        const double lf_sc_value = 0.5 * ((*lf_data)(si.toCell(0)) + (*lf_data)(si.toCell(1)));
                        multiplier_term *= lf_sc_value;
                    }
                    (*F_data)(si) *= multiplier_term;
                }
            }
        }
    }
}

DensityForceMask::DensityForceMask(Pointer<INSVCStaggeredHierarchyIntegrator> integrator,
                                   const double rho_liquid,
                                   const double rho_gas)
    : d_integrator(integrator), d_rho_liquid(rho_liquid), d_rho_gas(rho_gas)
{
}

void
DensityForceMask::mask_surface_tension_force(const int F_idx,
                                             Pointer<HierarchyMathOps> hier_math_ops,
                                             int /*integrator_step*/,
                                             double /*time*/,
                                             double /*current_time*/,
                                             double /*new_time*/,
                                             void* ctx)
{
    const auto& self = *static_cast<DensityForceMask*>(ctx);
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int rho_idx = self.d_integrator->getLinearOperatorRhoPatchDataIndex();
    mask_force(patch_hierarchy, F_idx, rho_idx, -1, self.d_rho_liquid, self.d_rho_gas);
}

SurfaceTensionCoefficients::SurfaceTensionCoefficients(Pointer<CellVariable<NDIM, double>> temperature,
                                                       Pointer<VariableContext> context,
                                                       const double sigma0,
                                                       const double dsigma_dT0,
                                                       const double T_ref)
    : d_temperature(temperature), d_context(context), d_sigma0(sigma0), d_dsigma_dT0(dsigma_dT0), d_T_ref(T_ref)
{
}

void
SurfaceTensionCoefficients::compute_surface_tension_coef_function(int F_idx,
                                                                  Pointer<Patch<NDIM>> patch,
                                                                  int /*integrator_step*/,
                                                                  double /*time*/,
                                                                  double /*current_time*/,
                                                                  double /*new_time*/,
                                                                  void* ctx)
{
    const auto& self = *static_cast<SurfaceTensionCoefficients*>(ctx);

    // parameters
    const double sigma_0 = self.d_sigma0;
    const double dsigma_dT0 = self.d_dsigma_dT0;
    const double T_ref = self.d_T_ref;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int T_scratch_idx = var_db->mapVariableAndContextToIndex(self.d_temperature, self.d_context);

    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CellData<NDIM, double>> T_data = patch->getPatchData(T_scratch_idx);
    Pointer<SideData<NDIM, double>> F_data = patch->getPatchData(F_idx);

    for (unsigned int axis = 0; axis < NDIM; axis++)
    {
        for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
        {
            SideIndex<NDIM> si(it(), axis, SideIndex<NDIM>::Lower);

            const double T_sc = 0.5 * ((*T_data)(si.toCell(0)) + (*T_data)(si.toCell(1)));
            const double sigma = sigma_0 + dsigma_dT0 * (T_sc - T_ref);
            (*F_data)(si) *= sigma;
        }
    }
    return;
}

void
SurfaceTensionCoefficients::compute_marangoni_coef_function(int F_idx,
                                                            Pointer<Patch<NDIM>> patch,
                                                            int /*integrator_step*/,
                                                            double /*time*/,
                                                            double /*current_time*/,
                                                            double /*new_time*/,
                                                            void* ctx)
{
    const auto& self = *static_cast<SurfaceTensionCoefficients*>(ctx);

    const double dsigma_dT0 = self.d_dsigma_dT0;
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<SideData<NDIM, double>> F_data = patch->getPatchData(F_idx);

    for (unsigned int axis = 0; axis < NDIM; axis++)
    {
        for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
        {
            SideIndex<NDIM> si(it(), axis, SideIndex<NDIM>::Lower);

            const double marangoni_coef = dsigma_dT0;
            (*F_data)(si) *= marangoni_coef; // Since marangoni_coef is constant for this example, it can be set
                                             // through input file as well.
        }
    }
    return;
}
} // namespace MultiphaseExamples
