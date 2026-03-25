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

#include <ibamr/INSSGSKinematics.h>

#include <ibtk/HierarchyGhostCellInterpolation.h>

#include <CellData.h>
#include <CellIndex.h>
#include <CellIterator.h>
#include <EdgeData.h>
#include <EdgeGeometry.h>
#include <EdgeIndex.h>
#include <EdgeIterator.h>
#include <NodeData.h>
#include <NodeIndex.h>
#include <NodeIterator.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RobinBcCoefStrategy.h>
#include <SideData.h>
#include <SideIndex.h>
#include <VariableDatabase.h>

#include <utility>

#include <ibamr/app_namespaces.h> // IWYU pragma: keep

namespace
{
void
allocate_patch_data(const int idx,
                    const double time,
                    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(idx)) level->allocatePatchData(idx, time);
    }
}

SAMRAI::hier::IntVector<NDIM>
axis_shift(const int axis, const int offset)
{
    SAMRAI::hier::IntVector<NDIM> shift(0);
    shift(axis) = offset;
    return shift;
}
} // namespace

namespace IBAMR
{
INSSGSKinematics::INSSGSKinematics(std::string object_name) : d_object_name(std::move(object_name))
{
    d_U_scratch_var = new SideVariable<NDIM, double>(d_object_name + "::U_scratch");
    d_S_cc_var = new CellVariable<NDIM, double>(d_object_name + "::S_cc", NDIM * (NDIM + 1) / 2);
#if (NDIM == 2)
    d_S_nc_var = new NodeVariable<NDIM, double>(d_object_name + "::S_nc", NDIM * (NDIM + 1) / 2, false);
#endif
#if (NDIM == 3)
    d_S_ec_var = new EdgeVariable<NDIM, double>(d_object_name + "::S_ec", NDIM * (NDIM + 1) / 2);
#endif

    auto* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext(d_object_name + "::CONTEXT");
    d_U_scratch_idx = var_db->registerVariableAndContext(d_U_scratch_var, ctx, IntVector<NDIM>(2));
    d_S_cc_idx = var_db->registerVariableAndContext(d_S_cc_var, ctx, IntVector<NDIM>(1));
#if (NDIM == 2)
    d_S_nc_idx = var_db->registerVariableAndContext(d_S_nc_var, ctx, IntVector<NDIM>(1));
#endif
#if (NDIM == 3)
    d_S_ec_idx = var_db->registerVariableAndContext(d_S_ec_var, ctx, IntVector<NDIM>(1));
#endif
}

void
INSSGSKinematics::fillGhostedVelocity(const int U_idx,
                                      const std::vector<RobinBcCoefStrategy<NDIM>*>& velocity_bc_coefs,
                                      const Pointer<PatchHierarchy<NDIM>> hierarchy,
                                      const double data_time)
{
    allocate_patch_data(d_U_scratch_idx, data_time, hierarchy);

    using ITC = IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    IBTK::HierarchyGhostCellInterpolation ghost_fill;
    ghost_fill.initializeOperatorState(ITC(d_U_scratch_idx,
                                           U_idx,
                                           "CONSERVATIVE_LINEAR_REFINE",
                                           true,
                                           "CONSERVATIVE_COARSEN",
                                           "LINEAR",
                                           false,
                                           velocity_bc_coefs),
                                       hierarchy);
    ghost_fill.fillData(data_time);
}

#if (NDIM == 2)
void
INSSGSKinematics::computeNodeCenteredStrainRate(const Pointer<PatchHierarchy<NDIM>> hierarchy, const double data_time)
{
    allocate_patch_data(d_S_nc_idx, data_time, hierarchy);

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> U_data = patch->getPatchData(d_U_scratch_idx);
            Pointer<NodeData<NDIM, double>> S_data = patch->getPatchData(d_S_nc_idx);
            const Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();

            S_data->fillAll(0.0);
            const Box<NDIM> S_box = S_data->getGhostBox();
            for (NodeIterator<NDIM> ni(S_box); ni; ni++)
            {
                const NodeIndex<NDIM> idx = ni();
                const CellIndex<NDIM> c(idx);
                const CellIndex<NDIM> c_x_minus(idx - axis_shift(0, 1));
                const CellIndex<NDIM> c_y_minus(idx - axis_shift(1, 1));
                const CellIndex<NDIM> c_xy_minus(idx - axis_shift(0, 1) - axis_shift(1, 1));

                const double du0_dx0_upper = ((*U_data)(SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Upper)) -
                                              (*U_data)(SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Lower))) /
                                             dx[0];
                const double du0_dx0_lower = ((*U_data)(SideIndex<NDIM>(c_y_minus, 0, SideIndex<NDIM>::Upper)) -
                                              (*U_data)(SideIndex<NDIM>(c_y_minus, 0, SideIndex<NDIM>::Lower))) /
                                             dx[0];
                const double du0_dx0 = 0.5 * (du0_dx0_upper + du0_dx0_lower);

                const double du1_dx1_right = ((*U_data)(SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Upper)) -
                                              (*U_data)(SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Lower))) /
                                             dx[1];
                const double du1_dx1_left = ((*U_data)(SideIndex<NDIM>(c_x_minus, 1, SideIndex<NDIM>::Upper)) -
                                             (*U_data)(SideIndex<NDIM>(c_x_minus, 1, SideIndex<NDIM>::Lower))) /
                                            dx[1];
                const double du1_dx1 = 0.5 * (du1_dx1_right + du1_dx1_left);

                const double du0_dx1 = ((*U_data)(SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Lower)) -
                                        (*U_data)(SideIndex<NDIM>(c_y_minus, 0, SideIndex<NDIM>::Lower))) /
                                       dx[1];
                const double du1_dx0 = ((*U_data)(SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Lower)) -
                                        (*U_data)(SideIndex<NDIM>(c_x_minus, 1, SideIndex<NDIM>::Lower))) /
                                       dx[0];

                (*S_data)(idx, 0) = du0_dx0;
                (*S_data)(idx, 1) = du1_dx1;
                (*S_data)(idx, 2) = 0.5 * (du0_dx1 + du1_dx0);
            }
        }
    }
}
#endif

#if (NDIM == 3)
void
INSSGSKinematics::computeEdgeCenteredStrainRate(const Pointer<PatchHierarchy<NDIM>> hierarchy, const double data_time)
{
    allocate_patch_data(d_S_ec_idx, data_time, hierarchy);

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> U_data = patch->getPatchData(d_U_scratch_idx);
            Pointer<EdgeData<NDIM, double>> S_data = patch->getPatchData(d_S_ec_idx);
            const Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();

            S_data->fillAll(0.0);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                ArrayData<NDIM, double>& S_array = S_data->getArrayData(axis);
                const Box<NDIM> edge_box = EdgeGeometry<NDIM>::toEdgeBox(patch->getBox(), axis);
                for (EdgeIterator<NDIM> ei(edge_box, axis); ei; ei++)
                {
                    const EdgeIndex<NDIM> e_idx = ei();
                    const hier::Index<NDIM>& i = e_idx;
                    const CellIndex<NDIM> c(i);
                    const CellIndex<NDIM> c_xm(i - axis_shift(0, 1));
                    const CellIndex<NDIM> c_xp(i + axis_shift(0, 1));
                    const CellIndex<NDIM> c_ym(i - axis_shift(1, 1));
                    const CellIndex<NDIM> c_yp(i + axis_shift(1, 1));
                    const CellIndex<NDIM> c_zm(i - axis_shift(2, 1));
                    const CellIndex<NDIM> c_zp(i + axis_shift(2, 1));
                    const CellIndex<NDIM> c_xym(i - axis_shift(0, 1) - axis_shift(1, 1));
                    const CellIndex<NDIM> c_xzm(i - axis_shift(0, 1) - axis_shift(2, 1));
                    const CellIndex<NDIM> c_yzm(i - axis_shift(1, 1) - axis_shift(2, 1));

                    double Sxx = 0.0, Syy = 0.0, Szz = 0.0, Syz = 0.0, Sxz = 0.0, Sxy = 0.0;

                    if (axis == 0)
                    {
                        const double du0_dx0 = ((*U_data)(SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Upper)) +
                                                (*U_data)(SideIndex<NDIM>(c_ym, 0, SideIndex<NDIM>::Upper)) +
                                                (*U_data)(SideIndex<NDIM>(c_zm, 0, SideIndex<NDIM>::Upper)) +
                                                (*U_data)(SideIndex<NDIM>(c_yzm, 0, SideIndex<NDIM>::Upper)) -
                                                (*U_data)(SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_ym, 0, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_zm, 0, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_yzm, 0, SideIndex<NDIM>::Lower))) /
                                               (4.0 * dx[0]);
                        const double du0_dx1 = ((*U_data)(SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Lower)) +
                                                (*U_data)(SideIndex<NDIM>(c_zm, 0, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_ym, 0, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_yzm, 0, SideIndex<NDIM>::Lower))) /
                                               (2.0 * dx[1]);
                        const double du0_dx2 = ((*U_data)(SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Lower)) +
                                                (*U_data)(SideIndex<NDIM>(c_ym, 0, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_zm, 0, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_yzm, 0, SideIndex<NDIM>::Lower))) /
                                               (2.0 * dx[2]);
                        const double du1_dx0 =
                            ((*U_data)(SideIndex<NDIM>(c_xp, 1, SideIndex<NDIM>::Lower)) +
                             (*U_data)(SideIndex<NDIM>(c_xp - axis_shift(2, 1), 1, SideIndex<NDIM>::Lower)) -
                             (*U_data)(SideIndex<NDIM>(c_xm, 1, SideIndex<NDIM>::Lower)) -
                             (*U_data)(SideIndex<NDIM>(c_xzm, 1, SideIndex<NDIM>::Lower))) /
                            (4.0 * dx[0]);
                        const double du1_dx1 = ((*U_data)(SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Upper)) +
                                                (*U_data)(SideIndex<NDIM>(c_zm, 1, SideIndex<NDIM>::Upper)) -
                                                (*U_data)(SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_zm, 1, SideIndex<NDIM>::Lower))) /
                                               (2.0 * dx[1]);
                        const double du1_dx2 = ((*U_data)(SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_zm, 1, SideIndex<NDIM>::Lower))) /
                                               dx[2];
                        const double du2_dx0 =
                            ((*U_data)(SideIndex<NDIM>(c_xp, 2, SideIndex<NDIM>::Lower)) +
                             (*U_data)(SideIndex<NDIM>(c_xp - axis_shift(1, 1), 2, SideIndex<NDIM>::Lower)) -
                             (*U_data)(SideIndex<NDIM>(c_xm, 2, SideIndex<NDIM>::Lower)) -
                             (*U_data)(SideIndex<NDIM>(c_xym, 2, SideIndex<NDIM>::Lower))) /
                            (4.0 * dx[0]);
                        const double du2_dx1 = ((*U_data)(SideIndex<NDIM>(c, 2, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_ym, 2, SideIndex<NDIM>::Lower))) /
                                               dx[1];
                        const double du2_dx2 = ((*U_data)(SideIndex<NDIM>(c, 2, SideIndex<NDIM>::Upper)) +
                                                (*U_data)(SideIndex<NDIM>(c_ym, 2, SideIndex<NDIM>::Upper)) -
                                                (*U_data)(SideIndex<NDIM>(c, 2, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_ym, 2, SideIndex<NDIM>::Lower))) /
                                               (2.0 * dx[2]);

                        Sxx = du0_dx0;
                        Syy = du1_dx1;
                        Szz = du2_dx2;
                        Syz = 0.5 * (du1_dx2 + du2_dx1);
                        Sxz = 0.5 * (du0_dx2 + du2_dx0);
                        Sxy = 0.5 * (du0_dx1 + du1_dx0);
                    }
                    else if (axis == 1)
                    {
                        const double du0_dx0 = ((*U_data)(SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Upper)) +
                                                (*U_data)(SideIndex<NDIM>(c_zm, 0, SideIndex<NDIM>::Upper)) -
                                                (*U_data)(SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_zm, 0, SideIndex<NDIM>::Lower))) /
                                               (2.0 * dx[0]);
                        const double du0_dx1 =
                            ((*U_data)(SideIndex<NDIM>(c_yp, 0, SideIndex<NDIM>::Lower)) +
                             (*U_data)(SideIndex<NDIM>(c_yp - axis_shift(2, 1), 0, SideIndex<NDIM>::Lower)) -
                             (*U_data)(SideIndex<NDIM>(c_ym, 0, SideIndex<NDIM>::Lower)) -
                             (*U_data)(SideIndex<NDIM>(c_yzm, 0, SideIndex<NDIM>::Lower))) /
                            (4.0 * dx[1]);
                        const double du0_dx2 = ((*U_data)(SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_zm, 0, SideIndex<NDIM>::Lower))) /
                                               dx[2];
                        const double du1_dx0 = ((*U_data)(SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Lower)) +
                                                (*U_data)(SideIndex<NDIM>(c_zm, 1, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_xm, 1, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_xzm, 1, SideIndex<NDIM>::Lower))) /
                                               (2.0 * dx[0]);
                        const double du1_dx1 =
                            ((*U_data)(SideIndex<NDIM>(c_xp, 1, SideIndex<NDIM>::Upper)) +
                             (*U_data)(SideIndex<NDIM>(c_xp - axis_shift(2, 1), 1, SideIndex<NDIM>::Upper)) -
                             (*U_data)(SideIndex<NDIM>(c_xm, 1, SideIndex<NDIM>::Lower)) -
                             (*U_data)(SideIndex<NDIM>(c_xzm, 1, SideIndex<NDIM>::Lower))) /
                            (4.0 * dx[1]);
                        const double du1_dx2 = ((*U_data)(SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Lower)) +
                                                (*U_data)(SideIndex<NDIM>(c_xm, 1, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_zm, 1, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_xzm, 1, SideIndex<NDIM>::Lower))) /
                                               (2.0 * dx[2]);
                        const double du2_dx0 = ((*U_data)(SideIndex<NDIM>(c_xp, 2, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_xm, 2, SideIndex<NDIM>::Lower))) /
                                               (2.0 * dx[0]);
                        const double du2_dx1 = ((*U_data)(SideIndex<NDIM>(c, 2, SideIndex<NDIM>::Lower)) +
                                                (*U_data)(SideIndex<NDIM>(c_xm, 2, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_ym, 2, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_xym, 2, SideIndex<NDIM>::Lower))) /
                                               (2.0 * dx[1]);
                        const double du2_dx2 = ((*U_data)(SideIndex<NDIM>(c, 2, SideIndex<NDIM>::Upper)) +
                                                (*U_data)(SideIndex<NDIM>(c_xm, 2, SideIndex<NDIM>::Upper)) -
                                                (*U_data)(SideIndex<NDIM>(c, 2, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_xm, 2, SideIndex<NDIM>::Lower))) /
                                               (2.0 * dx[2]);

                        Sxx = du0_dx0;
                        Syy = du1_dx1;
                        Szz = du2_dx2;
                        Syz = 0.5 * (du1_dx2 + du2_dx1);
                        Sxz = 0.5 * (du0_dx2 + du2_dx0);
                        Sxy = 0.5 * (du0_dx1 + du1_dx0);
                    }
                    else
                    {
                        const double du0_dx0 = ((*U_data)(SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Upper)) +
                                                (*U_data)(SideIndex<NDIM>(c_ym, 0, SideIndex<NDIM>::Upper)) -
                                                (*U_data)(SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_ym, 0, SideIndex<NDIM>::Lower))) /
                                               (2.0 * dx[0]);
                        const double du0_dx1 = ((*U_data)(SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_ym, 0, SideIndex<NDIM>::Lower))) /
                                               dx[1];
                        const double du0_dx2 =
                            ((*U_data)(SideIndex<NDIM>(c_zp, 0, SideIndex<NDIM>::Lower)) +
                             (*U_data)(SideIndex<NDIM>(c_zp - axis_shift(1, 1), 0, SideIndex<NDIM>::Lower)) -
                             (*U_data)(SideIndex<NDIM>(c_zm, 0, SideIndex<NDIM>::Lower)) -
                             (*U_data)(SideIndex<NDIM>(c_yzm, 0, SideIndex<NDIM>::Lower))) /
                            (4.0 * dx[2]);
                        const double du1_dx0 = ((*U_data)(SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_xm, 1, SideIndex<NDIM>::Lower))) /
                                               dx[0];
                        const double du1_dx1 = ((*U_data)(SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Upper)) +
                                                (*U_data)(SideIndex<NDIM>(c_xm, 1, SideIndex<NDIM>::Upper)) -
                                                (*U_data)(SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_xm, 1, SideIndex<NDIM>::Lower))) /
                                               (2.0 * dx[1]);
                        const double du1_dx2 =
                            ((*U_data)(SideIndex<NDIM>(c_zp, 1, SideIndex<NDIM>::Lower)) +
                             (*U_data)(SideIndex<NDIM>(c_zp - axis_shift(0, 1), 1, SideIndex<NDIM>::Lower)) -
                             (*U_data)(SideIndex<NDIM>(c_zm, 1, SideIndex<NDIM>::Lower)) -
                             (*U_data)(SideIndex<NDIM>(c_xzm, 1, SideIndex<NDIM>::Lower))) /
                            (4.0 * dx[2]);
                        const double du2_dx0 = ((*U_data)(SideIndex<NDIM>(c, 2, SideIndex<NDIM>::Lower)) +
                                                (*U_data)(SideIndex<NDIM>(c_ym, 2, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_xm, 2, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_xym, 2, SideIndex<NDIM>::Lower))) /
                                               (2.0 * dx[0]);
                        const double du2_dx1 = ((*U_data)(SideIndex<NDIM>(c, 2, SideIndex<NDIM>::Lower)) +
                                                (*U_data)(SideIndex<NDIM>(c_xm, 2, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_ym, 2, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_xym, 2, SideIndex<NDIM>::Lower))) /
                                               (2.0 * dx[1]);
                        const double du2_dx2 = ((*U_data)(SideIndex<NDIM>(c, 2, SideIndex<NDIM>::Upper)) +
                                                (*U_data)(SideIndex<NDIM>(c_xm, 2, SideIndex<NDIM>::Upper)) +
                                                (*U_data)(SideIndex<NDIM>(c_ym, 2, SideIndex<NDIM>::Upper)) +
                                                (*U_data)(SideIndex<NDIM>(c_xym, 2, SideIndex<NDIM>::Upper)) -
                                                (*U_data)(SideIndex<NDIM>(c, 2, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_xm, 2, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_ym, 2, SideIndex<NDIM>::Lower)) -
                                                (*U_data)(SideIndex<NDIM>(c_xym, 2, SideIndex<NDIM>::Lower))) /
                                               (4.0 * dx[2]);

                        Sxx = du0_dx0;
                        Syy = du1_dx1;
                        Szz = du2_dx2;
                        Syz = 0.5 * (du1_dx2 + du2_dx1);
                        Sxz = 0.5 * (du0_dx2 + du2_dx0);
                        Sxy = 0.5 * (du0_dx1 + du1_dx0);
                    }

                    S_array(i, 0) = Sxx;
                    S_array(i, 1) = Syy;
                    S_array(i, 2) = Szz;
                    S_array(i, 3) = Syz;
                    S_array(i, 4) = Sxz;
                    S_array(i, 5) = Sxy;
                }
            }
        }
    }
}
#endif

void
INSSGSKinematics::computeCellCenteredStrainRate(const Pointer<PatchHierarchy<NDIM>> hierarchy, const double data_time)
{
    allocate_patch_data(d_S_cc_idx, data_time, hierarchy);

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> U_data = patch->getPatchData(d_U_scratch_idx);
            Pointer<CellData<NDIM, double>> S_data = patch->getPatchData(d_S_cc_idx);
            const Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();

            S_data->fillAll(0.0);
            const Box<NDIM> S_box = S_data->getGhostBox();
            for (CellIterator<NDIM> ci(S_box); ci; ci++)
            {
                const CellIndex<NDIM> idx = *ci;

#if (NDIM == 2)
                const CellIndex<NDIM> idx_x_minus(idx - axis_shift(0, 1));
                const CellIndex<NDIM> idx_x_plus(idx + axis_shift(0, 1));
                const CellIndex<NDIM> idx_y_minus(idx - axis_shift(1, 1));
                const CellIndex<NDIM> idx_y_plus(idx + axis_shift(1, 1));

                const double du0_dx0 =
                    ((*U_data)(SideIndex<NDIM>(idx, 0, 1)) - (*U_data)(SideIndex<NDIM>(idx, 0, 0))) / dx[0];
                const double du1_dx1 =
                    ((*U_data)(SideIndex<NDIM>(idx, 1, 1)) - (*U_data)(SideIndex<NDIM>(idx, 1, 0))) / dx[1];
                const double du0_dx1 =
                    ((*U_data)(SideIndex<NDIM>(idx_y_plus, 0, 0)) + (*U_data)(SideIndex<NDIM>(idx_y_plus, 0, 1)) -
                     (*U_data)(SideIndex<NDIM>(idx_y_minus, 0, 0)) - (*U_data)(SideIndex<NDIM>(idx_y_minus, 0, 1))) /
                    (4.0 * dx[1]);
                const double du1_dx0 =
                    ((*U_data)(SideIndex<NDIM>(idx_x_plus, 1, 0)) + (*U_data)(SideIndex<NDIM>(idx_x_plus, 1, 1)) -
                     (*U_data)(SideIndex<NDIM>(idx_x_minus, 1, 0)) - (*U_data)(SideIndex<NDIM>(idx_x_minus, 1, 1))) /
                    (4.0 * dx[0]);

                (*S_data)(idx, 0) = du0_dx0;
                (*S_data)(idx, 1) = du1_dx1;
                (*S_data)(idx, 2) = 0.5 * (du0_dx1 + du1_dx0);
#endif

#if (NDIM == 3)
                const CellIndex<NDIM> idx_x_minus(idx - axis_shift(0, 1));
                const CellIndex<NDIM> idx_x_plus(idx + axis_shift(0, 1));
                const CellIndex<NDIM> idx_y_minus(idx - axis_shift(1, 1));
                const CellIndex<NDIM> idx_y_plus(idx + axis_shift(1, 1));
                const CellIndex<NDIM> idx_z_minus(idx - axis_shift(2, 1));
                const CellIndex<NDIM> idx_z_plus(idx + axis_shift(2, 1));

                const double du0_dx0 =
                    ((*U_data)(SideIndex<NDIM>(idx, 0, 1)) - (*U_data)(SideIndex<NDIM>(idx, 0, 0))) / dx[0];
                const double du1_dx1 =
                    ((*U_data)(SideIndex<NDIM>(idx, 1, 1)) - (*U_data)(SideIndex<NDIM>(idx, 1, 0))) / dx[1];
                const double du2_dx2 =
                    ((*U_data)(SideIndex<NDIM>(idx, 2, 1)) - (*U_data)(SideIndex<NDIM>(idx, 2, 0))) / dx[2];

                const double du1_dx2 =
                    ((*U_data)(SideIndex<NDIM>(idx_z_plus, 1, 0)) + (*U_data)(SideIndex<NDIM>(idx_z_plus, 1, 1)) -
                     (*U_data)(SideIndex<NDIM>(idx_z_minus, 1, 0)) - (*U_data)(SideIndex<NDIM>(idx_z_minus, 1, 1))) /
                    (4.0 * dx[2]);
                const double du2_dx1 =
                    ((*U_data)(SideIndex<NDIM>(idx_y_plus, 2, 0)) + (*U_data)(SideIndex<NDIM>(idx_y_plus, 2, 1)) -
                     (*U_data)(SideIndex<NDIM>(idx_y_minus, 2, 0)) - (*U_data)(SideIndex<NDIM>(idx_y_minus, 2, 1))) /
                    (4.0 * dx[1]);
                const double du0_dx2 =
                    ((*U_data)(SideIndex<NDIM>(idx_z_plus, 0, 0)) + (*U_data)(SideIndex<NDIM>(idx_z_plus, 0, 1)) -
                     (*U_data)(SideIndex<NDIM>(idx_z_minus, 0, 0)) - (*U_data)(SideIndex<NDIM>(idx_z_minus, 0, 1))) /
                    (4.0 * dx[2]);
                const double du2_dx0 =
                    ((*U_data)(SideIndex<NDIM>(idx_x_plus, 2, 0)) + (*U_data)(SideIndex<NDIM>(idx_x_plus, 2, 1)) -
                     (*U_data)(SideIndex<NDIM>(idx_x_minus, 2, 0)) - (*U_data)(SideIndex<NDIM>(idx_x_minus, 2, 1))) /
                    (4.0 * dx[0]);
                const double du0_dx1 =
                    ((*U_data)(SideIndex<NDIM>(idx_y_plus, 0, 0)) + (*U_data)(SideIndex<NDIM>(idx_y_plus, 0, 1)) -
                     (*U_data)(SideIndex<NDIM>(idx_y_minus, 0, 0)) - (*U_data)(SideIndex<NDIM>(idx_y_minus, 0, 1))) /
                    (4.0 * dx[1]);
                const double du1_dx0 =
                    ((*U_data)(SideIndex<NDIM>(idx_x_plus, 1, 0)) + (*U_data)(SideIndex<NDIM>(idx_x_plus, 1, 1)) -
                     (*U_data)(SideIndex<NDIM>(idx_x_minus, 1, 0)) - (*U_data)(SideIndex<NDIM>(idx_x_minus, 1, 1))) /
                    (4.0 * dx[0]);

                (*S_data)(idx, 0) = du0_dx0;
                (*S_data)(idx, 1) = du1_dx1;
                (*S_data)(idx, 2) = du2_dx2;
                (*S_data)(idx, 3) = 0.5 * (du1_dx2 + du2_dx1);
                (*S_data)(idx, 4) = 0.5 * (du0_dx2 + du2_dx0);
                (*S_data)(idx, 5) = 0.5 * (du0_dx1 + du1_dx0);
#endif
            }
        }
    }
}

int
INSSGSKinematics::getGhostedVelocityPatchDataIndex() const
{
    return d_U_scratch_idx;
}

int
INSSGSKinematics::getCellCenteredStrainRatePatchDataIndex() const
{
    return d_S_cc_idx;
}

#if (NDIM == 2)
int
INSSGSKinematics::getNodeCenteredStrainRatePatchDataIndex() const
{
    return d_S_nc_idx;
}
#endif

#if (NDIM == 3)
int
INSSGSKinematics::getEdgeCenteredStrainRatePatchDataIndex() const
{
    return d_S_ec_idx;
}
#endif
} // namespace IBAMR
