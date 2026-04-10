// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibamr/INSSmagorinskyTurbulenceModel.h>
#include <ibamr/StokesSpecifications.h>

#include <ibtk/ibtk_utilities.h>

#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIndex.h>
#include <CellIterator.h>
#include <EdgeData.h>
#include <EdgeGeometry.h>
#include <EdgeIterator.h>
#include <NodeData.h>
#include <NodeIndex.h>
#include <NodeIterator.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SideData.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

#include <ibamr/app_namespaces.h> // IWYU pragma: keep

namespace IBAMR
{
INSSmagorinskyTurbulenceModel::INSSmagorinskyTurbulenceModel(std::string object_name, Pointer<Database> input_db)
    : INSTurbulenceModel(std::move(object_name))
{
    if (input_db)
    {
        d_smagorinsky_constant = input_db->getDoubleWithDefault("smagorinsky_constant", d_smagorinsky_constant);
        d_filter_width_scale = input_db->getDoubleWithDefault("filter_width_scale", d_filter_width_scale);
        d_max_turbulent_viscosity =
            input_db->getDoubleWithDefault("max_turbulent_viscosity", d_max_turbulent_viscosity);
    }

    d_kinematics = new INSSGSKinematics(d_object_name + "::Kinematics");
    d_tau_sgs_data = new INSSGSStressData(d_object_name + "::StaggeredStress");
}

void
INSSmagorinskyTurbulenceModel::computeTurbulenceForce(const int F_idx,
                                                      const Pointer<SideVariable<NDIM, double>> /*F_var*/,
                                                      const int U_idx,
                                                      const Pointer<SideVariable<NDIM, double>> /*U_var*/,
                                                      const std::vector<RobinBcCoefStrategy<NDIM>*>& velocity_bc_coefs,
                                                      const Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                      const double data_time,
                                                      const StokesSpecifications& problem_coefs)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> F_data = patch->getPatchData(F_idx);
            F_data->fillAll(0.0);
        }
    }

    const double rho = problem_coefs.getRho();
    if (IBTK::abs_equal_eps(rho, 0.0)) return;

    d_kinematics->fillGhostedVelocity(U_idx, velocity_bc_coefs, hierarchy, data_time);
    d_kinematics->computeCellCenteredStrainRate(hierarchy, data_time);
#if (NDIM == 2)
    d_kinematics->computeNodeCenteredStrainRate(hierarchy, data_time);
#endif
#if (NDIM == 3)
    d_kinematics->computeEdgeCenteredStrainRate(hierarchy, data_time);
#endif
    d_tau_sgs_data->allocatePatchData(data_time, hierarchy);
    d_tau_sgs_data->setToZero(hierarchy);

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> S_data =
                patch->getPatchData(d_kinematics->getCellCenteredStrainRatePatchDataIndex());
            Pointer<CellData<NDIM, double>> tau_diag_data =
                patch->getPatchData(d_tau_sgs_data->getDiagonalStressPatchDataIndex());
            const Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();
            double delta_vol = 1.0;
            for (int d = 0; d < NDIM; ++d) delta_vol *= dx[d];
            const double delta = d_filter_width_scale * std::pow(delta_vol, 1.0 / static_cast<double>(NDIM));
            const double cs_delta_sq = std::pow(d_smagorinsky_constant * delta, 2.0);

            tau_diag_data->fillAll(0.0);
            const Box<NDIM> tau_box = tau_diag_data->getGhostBox();
            for (CellIterator<NDIM> ci(tau_box); ci; ci++)
            {
                const CellIndex<NDIM> idx = *ci;

#if (NDIM == 2)
                const double Sxx = (*S_data)(idx, 0);
                const double Syy = (*S_data)(idx, 1);
                const double Sxy = (*S_data)(idx, 2);
                const double mag_S = std::sqrt(2.0 * (Sxx * Sxx + Syy * Syy + 2.0 * Sxy * Sxy));
                const double mu_t = rho * std::min(cs_delta_sq * mag_S, d_max_turbulent_viscosity);

                (*tau_diag_data)(idx, 0) = 2.0 * mu_t * Sxx;
                (*tau_diag_data)(idx, 1) = 2.0 * mu_t * Syy;
#endif

#if (NDIM == 3)
                const double Sxx = (*S_data)(idx, 0);
                const double Syy = (*S_data)(idx, 1);
                const double Szz = (*S_data)(idx, 2);
                const double Syz = (*S_data)(idx, 3);
                const double Sxz = (*S_data)(idx, 4);
                const double Sxy = (*S_data)(idx, 5);
                const double mag_S =
                    std::sqrt(2.0 * (Sxx * Sxx + Syy * Syy + Szz * Szz + 2.0 * (Syz * Syz + Sxz * Sxz + Sxy * Sxy)));
                const double mu_t = rho * std::min(cs_delta_sq * mag_S, d_max_turbulent_viscosity);

                (*tau_diag_data)(idx, 0) = 2.0 * mu_t * Sxx;
                (*tau_diag_data)(idx, 1) = 2.0 * mu_t * Syy;
                (*tau_diag_data)(idx, 2) = 2.0 * mu_t * Szz;
#endif
            }

#if (NDIM == 2)
            Pointer<NodeData<NDIM, double>> S_node_data =
                patch->getPatchData(d_kinematics->getNodeCenteredStrainRatePatchDataIndex());
            Pointer<NodeData<NDIM, double>> tau_shear_data =
                patch->getPatchData(d_tau_sgs_data->getShearStressPatchDataIndex());
            tau_shear_data->fillAll(0.0);
            const Box<NDIM> tau_node_box = tau_shear_data->getGhostBox();
            for (NodeIterator<NDIM> ni(tau_node_box); ni; ni++)
            {
                const NodeIndex<NDIM> idx = ni();
                const double Sxx = (*S_node_data)(idx, 0);
                const double Syy = (*S_node_data)(idx, 1);
                const double Sxy = (*S_node_data)(idx, 2);
                const double mag_S = std::sqrt(2.0 * (Sxx * Sxx + Syy * Syy + 2.0 * Sxy * Sxy));
                const double mu_t = rho * std::min(cs_delta_sq * mag_S, d_max_turbulent_viscosity);
                (*tau_shear_data)(idx) = 2.0 * mu_t * Sxy;
            }
#endif

#if (NDIM == 3)
            Pointer<EdgeData<NDIM, double>> S_edge_data =
                patch->getPatchData(d_kinematics->getEdgeCenteredStrainRatePatchDataIndex());
            Pointer<EdgeData<NDIM, double>> tau_shear_data =
                patch->getPatchData(d_tau_sgs_data->getShearStressPatchDataIndex());
            tau_shear_data->fillAll(0.0);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                ArrayData<NDIM, double>& tau_array = tau_shear_data->getArrayData(axis);
                const ArrayData<NDIM, double>& S_array = S_edge_data->getArrayData(axis);
                const Box<NDIM> edge_box = EdgeGeometry<NDIM>::toEdgeBox(patch->getBox(), axis);
                const int shear_depth = axis == 0 ? 3 : (axis == 1 ? 4 : 5);
                for (EdgeIterator<NDIM> ei(edge_box, axis); ei; ei++)
                {
                    const hier::Index<NDIM>& i = EdgeIndex<NDIM>(ei());
                    const double Sxx = S_array(i, 0);
                    const double Syy = S_array(i, 1);
                    const double Szz = S_array(i, 2);
                    const double Syz = S_array(i, 3);
                    const double Sxz = S_array(i, 4);
                    const double Sxy = S_array(i, 5);
                    const double mag_S = std::sqrt(
                        2.0 * (Sxx * Sxx + Syy * Syy + Szz * Szz + 2.0 * (Syz * Syz + Sxz * Sxz + Sxy * Sxy)));
                    const double mu_t = rho * std::min(cs_delta_sq * mag_S, d_max_turbulent_viscosity);
                    tau_array(i, 0) = 2.0 * mu_t * S_array(i, shear_depth);
                }
            }
#endif
        }
    }

    d_tau_sgs_data->computeDivergence(F_idx, hierarchy);
    d_tau_sgs_data->deallocatePatchData(hierarchy);
}
} // namespace IBAMR
