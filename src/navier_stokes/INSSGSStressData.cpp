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

#include <ibamr/INSSGSStressData.h>

#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIndex.h>
#include <CellIterator.h>
#include <EdgeData.h>
#include <EdgeIndex.h>
#include <EdgeIterator.h>
#include <NodeData.h>
#include <NodeIndex.h>
#include <NodeIterator.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SideData.h>
#include <SideIndex.h>
#include <VariableDatabase.h>

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

void
deallocate_patch_data(const int idx, const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(idx)) level->deallocatePatchData(idx);
    }
}

SAMRAI::hier::Index<NDIM>
shift_index(const int axis, const int offset)
{
    SAMRAI::hier::Index<NDIM> idx(0);
    idx(axis) = offset;
    return idx;
}
} // namespace

namespace IBAMR
{
INSSGSStressData::INSSGSStressData(std::string object_name) : d_object_name(std::move(object_name))
{
    d_tau_diag_var = new CellVariable<NDIM, double>(d_object_name + "::tau_diag", NDIM);
#if (NDIM == 2)
    d_tau_shear_var = new NodeVariable<NDIM, double>(d_object_name + "::tau_xy", 1, false);
#endif
#if (NDIM == 3)
    d_tau_shear_var = new EdgeVariable<NDIM, double>(d_object_name + "::tau_shear", 1);
#endif

    auto* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext(d_object_name + "::CONTEXT");
    d_tau_diag_idx = var_db->registerVariableAndContext(d_tau_diag_var, ctx, IntVector<NDIM>(1));
    d_tau_shear_idx = var_db->registerVariableAndContext(d_tau_shear_var, ctx, IntVector<NDIM>(1));
}

void
INSSGSStressData::allocatePatchData(const double data_time, const Pointer<PatchHierarchy<NDIM>> hierarchy)
{
    allocate_patch_data(d_tau_diag_idx, data_time, hierarchy);
    allocate_patch_data(d_tau_shear_idx, data_time, hierarchy);
}

void
INSSGSStressData::deallocatePatchData(const Pointer<PatchHierarchy<NDIM>> hierarchy)
{
    deallocate_patch_data(d_tau_diag_idx, hierarchy);
    deallocate_patch_data(d_tau_shear_idx, hierarchy);
}

void
INSSGSStressData::setToZero(const Pointer<PatchHierarchy<NDIM>> hierarchy)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> tau_diag_data = patch->getPatchData(d_tau_diag_idx);
            tau_diag_data->fillAll(0.0);
#if (NDIM == 2)
            Pointer<NodeData<NDIM, double>> tau_shear_data = patch->getPatchData(d_tau_shear_idx);
            tau_shear_data->fillAll(0.0);
#endif
#if (NDIM == 3)
            Pointer<EdgeData<NDIM, double>> tau_shear_data = patch->getPatchData(d_tau_shear_idx);
            tau_shear_data->fillAll(0.0);
#endif
        }
    }
}

void
INSSGSStressData::setFromCellCenteredStressTensor(const int tau_cc_idx, const Pointer<PatchHierarchy<NDIM>> hierarchy)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> tau_cc_data = patch->getPatchData(tau_cc_idx);
            Pointer<CellData<NDIM, double>> tau_diag_data = patch->getPatchData(d_tau_diag_idx);
            const Box<NDIM> cell_box = tau_diag_data->getGhostBox();
            for (CellIterator<NDIM> ci(cell_box); ci; ci++)
            {
                const CellIndex<NDIM> idx = *ci;
#if (NDIM == 2)
                (*tau_diag_data)(idx, 0) = (*tau_cc_data)(idx, 0);
                (*tau_diag_data)(idx, 1) = (*tau_cc_data)(idx, 1);
#endif
#if (NDIM == 3)
                (*tau_diag_data)(idx, 0) = (*tau_cc_data)(idx, 0);
                (*tau_diag_data)(idx, 1) = (*tau_cc_data)(idx, 1);
                (*tau_diag_data)(idx, 2) = (*tau_cc_data)(idx, 2);
#endif
            }

#if (NDIM == 2)
            Pointer<NodeData<NDIM, double>> tau_shear_data = patch->getPatchData(d_tau_shear_idx);
            const Box<NDIM> node_box = tau_shear_data->getGhostBox();
            for (NodeIterator<NDIM> ni(node_box); ni; ni++)
            {
                const NodeIndex<NDIM> idx = ni();
                const CellIndex<NDIM> ll(idx - 1);
                const CellIndex<NDIM> lr(idx - shift_index(1, 1));
                const CellIndex<NDIM> ul(idx - shift_index(0, 1));
                const CellIndex<NDIM> ur(idx);
                (*tau_shear_data)(idx) = 0.25 * ((*tau_cc_data)(ll, 2) + (*tau_cc_data)(lr, 2) + (*tau_cc_data)(ul, 2) +
                                                 (*tau_cc_data)(ur, 2));
            }
#endif

#if (NDIM == 3)
            Pointer<EdgeData<NDIM, double>> tau_shear_data = patch->getPatchData(d_tau_shear_idx);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> edge_box = tau_shear_data->getArrayData(axis).getBox();
                for (EdgeIterator<NDIM> ei(edge_box, axis); ei; ei++)
                {
                    const EdgeIndex<NDIM> e_idx = ei();
                    const hier::Index<NDIM>& i = e_idx;
                    double value = 0.0;
                    if (axis == 0)
                    {
                        value = 0.25 * ((*tau_cc_data)(CellIndex<NDIM>(i - shift_index(1, 1) - shift_index(2, 1)), 3) +
                                        (*tau_cc_data)(CellIndex<NDIM>(i - shift_index(2, 1)), 3) +
                                        (*tau_cc_data)(CellIndex<NDIM>(i - shift_index(1, 1)), 3) +
                                        (*tau_cc_data)(CellIndex<NDIM>(i), 3));
                    }
                    else if (axis == 1)
                    {
                        value = 0.25 * ((*tau_cc_data)(CellIndex<NDIM>(i - shift_index(0, 1) - shift_index(2, 1)), 4) +
                                        (*tau_cc_data)(CellIndex<NDIM>(i - shift_index(2, 1)), 4) +
                                        (*tau_cc_data)(CellIndex<NDIM>(i - shift_index(0, 1)), 4) +
                                        (*tau_cc_data)(CellIndex<NDIM>(i), 4));
                    }
                    else
                    {
                        value = 0.25 * ((*tau_cc_data)(CellIndex<NDIM>(i - shift_index(0, 1) - shift_index(1, 1)), 5) +
                                        (*tau_cc_data)(CellIndex<NDIM>(i - shift_index(1, 1)), 5) +
                                        (*tau_cc_data)(CellIndex<NDIM>(i - shift_index(0, 1)), 5) +
                                        (*tau_cc_data)(CellIndex<NDIM>(i), 5));
                    }
                    tau_shear_data->getArrayData(axis)(i, 0) = value;
                }
            }
#endif
        }
    }
}

void
INSSGSStressData::computeDivergence(const int F_idx, const Pointer<PatchHierarchy<NDIM>> hierarchy) const
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> F_data = patch->getPatchData(F_idx);
            Pointer<CellData<NDIM, double>> tau_diag_data = patch->getPatchData(d_tau_diag_idx);
            const Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();

#if (NDIM == 2)
            Pointer<NodeData<NDIM, double>> tau_shear_data = patch->getPatchData(d_tau_shear_idx);
            for (SideIterator<NDIM> si(patch->getBox(), 0); si; si++)
            {
                const SideIndex<NDIM> s_idx = si();
                const CellIndex<NDIM> c_l = s_idx.toCell(SideIndex<NDIM>::Lower);
                const CellIndex<NDIM> c_u = s_idx.toCell(SideIndex<NDIM>::Upper);
                const NodeIndex<NDIM> n_l(c_u, NodeIndex<NDIM>::LowerLeft);
                const NodeIndex<NDIM> n_u(c_u, NodeIndex<NDIM>::UpperLeft);
                (*F_data)(s_idx) = ((*tau_diag_data)(c_u, 0) - (*tau_diag_data)(c_l, 0)) / dx[0] +
                                   ((*tau_shear_data)(n_u) - (*tau_shear_data)(n_l)) / dx[1];
            }
            for (SideIterator<NDIM> si(patch->getBox(), 1); si; si++)
            {
                const SideIndex<NDIM> s_idx = si();
                const CellIndex<NDIM> c_l = s_idx.toCell(SideIndex<NDIM>::Lower);
                const CellIndex<NDIM> c_u = s_idx.toCell(SideIndex<NDIM>::Upper);
                const NodeIndex<NDIM> n_l(c_u, NodeIndex<NDIM>::LowerLeft);
                const NodeIndex<NDIM> n_u(c_u, NodeIndex<NDIM>::LowerRight);
                (*F_data)(s_idx) = ((*tau_shear_data)(n_u) - (*tau_shear_data)(n_l)) / dx[0] +
                                   ((*tau_diag_data)(c_u, 1) - (*tau_diag_data)(c_l, 1)) / dx[1];
            }
#endif

#if (NDIM == 3)
            Pointer<EdgeData<NDIM, double>> tau_shear_data = patch->getPatchData(d_tau_shear_idx);
            const ArrayData<NDIM, double>& tau_yz_data = tau_shear_data->getArrayData(0);
            const ArrayData<NDIM, double>& tau_xz_data = tau_shear_data->getArrayData(1);
            const ArrayData<NDIM, double>& tau_xy_data = tau_shear_data->getArrayData(2);

            for (SideIterator<NDIM> si(patch->getBox(), 0); si; si++)
            {
                const SideIndex<NDIM> s_idx = si();
                const CellIndex<NDIM> c_l = s_idx.toCell(SideIndex<NDIM>::Lower);
                const CellIndex<NDIM> c_u = s_idx.toCell(SideIndex<NDIM>::Upper);
                const hier::Index<NDIM>& i = c_u;
                (*F_data)(s_idx) = ((*tau_diag_data)(c_u, 0) - (*tau_diag_data)(c_l, 0)) / dx[0] +
                                   (tau_xy_data(i + shift_index(1, 1), 0) - tau_xy_data(i, 0)) / dx[1] +
                                   (tau_xz_data(i + shift_index(2, 1), 0) - tau_xz_data(i, 0)) / dx[2];
            }
            for (SideIterator<NDIM> si(patch->getBox(), 1); si; si++)
            {
                const SideIndex<NDIM> s_idx = si();
                const CellIndex<NDIM> c_l = s_idx.toCell(SideIndex<NDIM>::Lower);
                const CellIndex<NDIM> c_u = s_idx.toCell(SideIndex<NDIM>::Upper);
                const hier::Index<NDIM>& i = c_u;
                (*F_data)(s_idx) = (tau_xy_data(i + shift_index(0, 1), 0) - tau_xy_data(i, 0)) / dx[0] +
                                   ((*tau_diag_data)(c_u, 1) - (*tau_diag_data)(c_l, 1)) / dx[1] +
                                   (tau_yz_data(i + shift_index(2, 1), 0) - tau_yz_data(i, 0)) / dx[2];
            }
            for (SideIterator<NDIM> si(patch->getBox(), 2); si; si++)
            {
                const SideIndex<NDIM> s_idx = si();
                const CellIndex<NDIM> c_l = s_idx.toCell(SideIndex<NDIM>::Lower);
                const CellIndex<NDIM> c_u = s_idx.toCell(SideIndex<NDIM>::Upper);
                const hier::Index<NDIM>& i = c_u;
                (*F_data)(s_idx) = (tau_xz_data(i + shift_index(0, 1), 0) - tau_xz_data(i, 0)) / dx[0] +
                                   (tau_yz_data(i + shift_index(1, 1), 0) - tau_yz_data(i, 0)) / dx[1] +
                                   ((*tau_diag_data)(c_u, 2) - (*tau_diag_data)(c_l, 2)) / dx[2];
            }
#endif
        }
    }
}

int
INSSGSStressData::getDiagonalStressPatchDataIndex() const
{
    return d_tau_diag_idx;
}

int
INSSGSStressData::getShearStressPatchDataIndex() const
{
    return d_tau_shear_idx;
}
} // namespace IBAMR
