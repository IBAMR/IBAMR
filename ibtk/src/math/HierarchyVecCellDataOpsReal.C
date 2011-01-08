// Filename: HierarchyVecCellDataOpsReal.C
// Created on 09 Apr 2010 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include "HierarchyVecCellDataOpsReal.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/VecCellDataFactory.h>
#include <ibtk/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

template<class TYPE>
HierarchyVecCellDataOpsReal<TYPE>::HierarchyVecCellDataOpsReal(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level)
    : HierarchyDataOpsReal<NDIM,TYPE>(),
      d_hierarchy(hierarchy),
      d_coarsest_level(coarsest_level),
      d_finest_level(finest_level),
      d_patch_ops()
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
#endif
    d_hierarchy = hierarchy;
    if ((coarsest_level < 0) || (finest_level < 0))
    {
        if (d_hierarchy->getNumberOfLevels() == 0)
        {
            d_coarsest_level = coarsest_level;
            d_finest_level = finest_level;
        }
        else
        {
            resetLevels(0, d_hierarchy->getFinestLevelNumber());
        }
    }
    else
    {
        resetLevels(coarsest_level, finest_level);
    }
    return;
}// HierarchyVecCellDataOpsReal

template<class TYPE>
HierarchyVecCellDataOpsReal<TYPE>::~HierarchyVecCellDataOpsReal()
{
    // intentionally blank
    return;
}// ~HierarchyVecCellDataOpsReal

template<class TYPE>
void
HierarchyVecCellDataOpsReal<TYPE>::setPatchHierarchy(
    Pointer<PatchHierarchy<NDIM> > hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
#endif
    d_hierarchy = hierarchy;
    return;
}// setPatchHierarchy

template<class TYPE>
void
HierarchyVecCellDataOpsReal<TYPE>::resetLevels(
    const int coarsest_level,
    const int finest_level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
    TBOX_ASSERT((coarsest_level >= 0)
                && (finest_level >= coarsest_level)
                && (finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif
    d_coarsest_level = coarsest_level;
    d_finest_level = finest_level;
    return;
}// resetLevels

template<class TYPE>
const Pointer<PatchHierarchy<NDIM> >
HierarchyVecCellDataOpsReal<TYPE>::getPatchHierarchy() const
{
    return d_hierarchy;
}// getPatchHierarchy

template<class TYPE>
void
HierarchyVecCellDataOpsReal<TYPE>::copyData(
    const int dst_id,
    const int src_id,
    const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
    TBOX_ASSERT((d_coarsest_level >= 0)
                && (d_finest_level >= d_coarsest_level)
                && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif
    for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (typename PatchLevel<NDIM>::Iterator ip(level); ip; ip++)
        {
            Pointer<Patch<NDIM> > p = level->getPatch(ip());

            Pointer<VecCellData<TYPE> > dst = p->getPatchData(dst_id);
            Pointer<VecCellData<TYPE> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!dst.isNull());
#endif
            const Box<NDIM>& box = interior_only ? p->getBox() : dst->getGhostBox();
            d_patch_ops.copyData(dst, src, box);
        }
    }
    return;
}// copyData

template<class TYPE>
void
HierarchyVecCellDataOpsReal<TYPE>::swapData(
    const int data1_id,
    const int data2_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    Pointer<VecCellDataFactory<TYPE> > d1fact = d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data1_id);
    TBOX_ASSERT(!d1fact.isNull());
    Pointer<VecCellDataFactory<TYPE> > d2fact = d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data2_id);
    TBOX_ASSERT(!d2fact.isNull());
    TBOX_ASSERT(d1fact->getDefaultDepth() == d2fact->getDefaultDepth());
    TBOX_ASSERT(d1fact->getGhostCellWidth() == d2fact->getGhostCellWidth());
#endif
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
    TBOX_ASSERT((d_coarsest_level >= 0)
                && (d_finest_level >= d_coarsest_level)
                && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif
    for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (typename PatchLevel<NDIM>::Iterator ip(level); ip; ip++)
        {
            Pointer<Patch<NDIM> > p = level->getPatch(ip());
            d_patch_ops.swapData(p, data1_id, data2_id);
        }
    }
    return;
}// swapData

template<class TYPE>
void
HierarchyVecCellDataOpsReal<TYPE>::printData(
    const int data_id,
    std::ostream& s,
    const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
    TBOX_ASSERT((d_coarsest_level >= 0)
                && (d_finest_level >= d_coarsest_level)
                && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif
    s << "Patch descriptor id = " << data_id << std::endl;
    s << "Factory = " << typeid(*d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data_id)).name() << std::endl;
    for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln)
    {
        s << "Level number = " << ln << std::endl;
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (typename PatchLevel<NDIM>::Iterator ip(level); ip; ip++)
        {
            Pointer<Patch<NDIM> > p = level->getPatch(ip());

            Pointer<VecCellData<TYPE> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!d.isNull());
#endif
            const Box<NDIM>& box = interior_only ? p->getBox() : d->getGhostBox();
            d_patch_ops.printData(d, box, s);
        }
    }
    return;
}// printData

template<class TYPE>
void
HierarchyVecCellDataOpsReal<TYPE>::setToScalar(
    const int data_id,
    const TYPE& alpha,
    const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
    TBOX_ASSERT((d_coarsest_level >= 0)
                && (d_finest_level >= d_coarsest_level)
                && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif
    for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (typename PatchLevel<NDIM>::Iterator ip(level); ip; ip++)
        {
            Pointer<Patch<NDIM> > p = level->getPatch(ip());

            Pointer<VecCellData<TYPE> > d = p->getPatchData(data_id);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!d.isNull());
#endif
            const Box<NDIM>& box = interior_only ? p->getBox() : d->getGhostBox();
            d_patch_ops.setToScalar(d, alpha, box);
        }
    }
    return;
}// setToScalar

template<class TYPE>
void
HierarchyVecCellDataOpsReal<TYPE>::scale(
    const int dst_id,
    const TYPE& alpha,
    const int src_id,
    const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
    TBOX_ASSERT((d_coarsest_level >= 0)
                && (d_finest_level >= d_coarsest_level)
                && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

    for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (typename PatchLevel<NDIM>::Iterator ip(level); ip; ip++)
        {
            Pointer<Patch<NDIM> > p = level->getPatch(ip());

            Pointer<VecCellData<TYPE> > dst = p->getPatchData(dst_id);
            Pointer<VecCellData<TYPE> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!dst.isNull());
#endif
            const Box<NDIM>& box = interior_only ? p->getBox() : dst->getGhostBox();
            d_patch_ops.scale(dst, alpha, src, box);
        }
    }
    return;
}// scale

template<class TYPE>
void
HierarchyVecCellDataOpsReal<TYPE>::addScalar(
    const int dst_id,
    const int src_id,
    const TYPE& alpha,
    const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
    TBOX_ASSERT((d_coarsest_level >= 0)
                && (d_finest_level >= d_coarsest_level)
                && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

    for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (typename PatchLevel<NDIM>::Iterator ip(level); ip; ip++)
        {
            Pointer<Patch<NDIM> > p = level->getPatch(ip());

            Pointer<VecCellData<TYPE> > dst = p->getPatchData(dst_id);
            Pointer<VecCellData<TYPE> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!dst.isNull());
#endif
            const Box<NDIM>& box = interior_only ? p->getBox() : dst->getGhostBox();
            d_patch_ops.addScalar(dst, src, alpha, box);
        }
    }
    return;
}// addScalar

template<class TYPE>
void
HierarchyVecCellDataOpsReal<TYPE>::add(
    const int dst_id,
    const int src1_id,
    const int src2_id,
    const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
    TBOX_ASSERT((d_coarsest_level >= 0)
                && (d_finest_level >= d_coarsest_level)
                && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

    for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (typename PatchLevel<NDIM>::Iterator ip(level); ip; ip++)
        {
            Pointer<Patch<NDIM> > p = level->getPatch(ip());

            Pointer<VecCellData<TYPE> > dst  = p->getPatchData(dst_id);
            Pointer<VecCellData<TYPE> > src1 = p->getPatchData(src1_id);
            Pointer<VecCellData<TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!dst.isNull());
#endif
            const Box<NDIM>& box = interior_only ? p->getBox() : dst->getGhostBox();
            d_patch_ops.add(dst, src1, src2, box);
        }
    }
    return;
}// add

template<class TYPE>
void
HierarchyVecCellDataOpsReal<TYPE>::subtract(
    const int dst_id,
    const int src1_id,
    const int src2_id,
    const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
    TBOX_ASSERT((d_coarsest_level >= 0)
                && (d_finest_level >= d_coarsest_level)
                && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

    for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (typename PatchLevel<NDIM>::Iterator ip(level); ip; ip++)
        {
            Pointer<Patch<NDIM> > p = level->getPatch(ip());

            Pointer<VecCellData<TYPE> > dst  = p->getPatchData(dst_id);
            Pointer<VecCellData<TYPE> > src1 = p->getPatchData(src1_id);
            Pointer<VecCellData<TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!dst.isNull());
#endif
            const Box<NDIM>& box = interior_only ? p->getBox() : dst->getGhostBox();
            d_patch_ops.subtract(dst, src1, src2, box);
        }
    }
    return;
}// subtract

template<class TYPE>
void
HierarchyVecCellDataOpsReal<TYPE>::multiply(
    const int dst_id,
    const int src1_id,
    const int src2_id,
    const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
    TBOX_ASSERT((d_coarsest_level >= 0)
                && (d_finest_level >= d_coarsest_level)
                && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

    for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (typename PatchLevel<NDIM>::Iterator ip(level); ip; ip++)
        {
            Pointer<Patch<NDIM> > p = level->getPatch(ip());

            Pointer<VecCellData<TYPE> > dst  = p->getPatchData(dst_id);
            Pointer<VecCellData<TYPE> > src1 = p->getPatchData(src1_id);
            Pointer<VecCellData<TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!dst.isNull());
#endif
            const Box<NDIM>& box = interior_only ? p->getBox() : dst->getGhostBox();
            d_patch_ops.multiply(dst, src1, src2, box);
        }
    }
    return;
}// multiply

template<class TYPE>
void
HierarchyVecCellDataOpsReal<TYPE>::divide(
    const int dst_id,
    const int src1_id,
    const int src2_id,
    const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
    TBOX_ASSERT((d_coarsest_level >= 0)
                && (d_finest_level >= d_coarsest_level)
                && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

    for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (typename PatchLevel<NDIM>::Iterator ip(level); ip; ip++)
        {
            Pointer<Patch<NDIM> > p = level->getPatch(ip());

            Pointer<VecCellData<TYPE> > dst  = p->getPatchData(dst_id);
            Pointer<VecCellData<TYPE> > src1 = p->getPatchData(src1_id);
            Pointer<VecCellData<TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!dst.isNull());
#endif
            const Box<NDIM>& box = interior_only ? p->getBox() : dst->getGhostBox();
            d_patch_ops.divide(dst, src1, src2, box);
        }
    }
    return;
}// divide

template<class TYPE>
void
HierarchyVecCellDataOpsReal<TYPE>::reciprocal(
    const int dst_id,
    const int src_id,
    const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
    TBOX_ASSERT((d_coarsest_level >= 0)
                && (d_finest_level >= d_coarsest_level)
                && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

    for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (typename PatchLevel<NDIM>::Iterator ip(level); ip; ip++)
        {
            Pointer<Patch<NDIM> > p = level->getPatch(ip());

            Pointer<VecCellData<TYPE> > dst = p->getPatchData(dst_id);
            Pointer<VecCellData<TYPE> > src = p->getPatchData(src_id);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!dst.isNull());
#endif
            const Box<NDIM>& box = interior_only ? p->getBox() : dst->getGhostBox();
            d_patch_ops.reciprocal(dst, src, box);
        }
    }
    return;
}// reciprocal

template<class TYPE>
void
HierarchyVecCellDataOpsReal<TYPE>::linearSum(
    const int dst_id,
    const TYPE& alpha,
    const int src1_id,
    const TYPE& beta,
    const int src2_id,
    const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
    TBOX_ASSERT((d_coarsest_level >= 0)
                && (d_finest_level >= d_coarsest_level)
                && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

    for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (typename PatchLevel<NDIM>::Iterator ip(level); ip; ip++)
        {
            Pointer<Patch<NDIM> > p = level->getPatch(ip());

            Pointer<VecCellData<TYPE> > dst  = p->getPatchData(dst_id);
            Pointer<VecCellData<TYPE> > src1 = p->getPatchData(src1_id);
            Pointer<VecCellData<TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!dst.isNull());
#endif
            const Box<NDIM>& box = interior_only ? p->getBox() : dst->getGhostBox();
            d_patch_ops.linearSum(dst, alpha, src1, beta, src2, box);
        }
    }
    return;
}// linearSum

template<class TYPE>
void
HierarchyVecCellDataOpsReal<TYPE>::axpy(
    const int dst_id,
    const TYPE& alpha,
    const int src1_id,
    const int src2_id,
    const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
    TBOX_ASSERT((d_coarsest_level >= 0)
                && (d_finest_level >= d_coarsest_level)
                && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

    for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (typename PatchLevel<NDIM>::Iterator ip(level); ip; ip++)
        {
            Pointer<Patch<NDIM> > p = level->getPatch(ip());

            Pointer<VecCellData<TYPE> > dst  = p->getPatchData(dst_id);
            Pointer<VecCellData<TYPE> > src1 = p->getPatchData(src1_id);
            Pointer<VecCellData<TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!dst.isNull());
#endif
            const Box<NDIM>& box = interior_only ? p->getBox() : dst->getGhostBox();
            d_patch_ops.axpy(dst, alpha, src1, src2, box);
        }
    }
    return;
}// axpy

template<class TYPE>
void
HierarchyVecCellDataOpsReal<TYPE>::axmy(
    const int dst_id,
    const TYPE& alpha,
    const int src1_id,
    const int src2_id,
    const bool interior_only) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hierarchy.isNull());
    TBOX_ASSERT((d_coarsest_level >= 0)
                && (d_finest_level >= d_coarsest_level)
                && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));
#endif

    for (int ln = d_coarsest_level; ln <= d_finest_level; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (typename PatchLevel<NDIM>::Iterator ip(level); ip; ip++)
        {
            Pointer<Patch<NDIM> > p = level->getPatch(ip());

            Pointer<VecCellData<TYPE> > dst  = p->getPatchData(dst_id);
            Pointer<VecCellData<TYPE> > src1 = p->getPatchData(src1_id);
            Pointer<VecCellData<TYPE> > src2 = p->getPatchData(src2_id);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!dst.isNull());
#endif
            const Box<NDIM>& box = interior_only ? p->getBox() : dst->getGhostBox();
            d_patch_ops.axmy(dst, alpha, src1, src2, box);
        }
    }
    return;
}// axmy

template<class TYPE>
void
HierarchyVecCellDataOpsReal<TYPE>::abs(
    const int dst_id,
    const int src_id,
    const bool interior_only) const
{
    TBOX_ERROR("HierarchyVecCellDataOpsReal<TYPE>::abs() unimplemented" << std::endl);
    return;
}// abs

template<class TYPE>
TYPE
HierarchyVecCellDataOpsReal<TYPE>::min(
    const int data_id,
    const bool interior_only) const
{
    TBOX_ERROR("HierarchyVecCellDataOpsReal<TYPE>::min() unimplemented" << std::endl);
    return 0.0;
}// min

template<class TYPE>
TYPE
HierarchyVecCellDataOpsReal<TYPE>::max(
    const int data_id,
    const bool interior_only) const
{
    TBOX_ERROR("HierarchyVecCellDataOpsReal<TYPE>::max() unimplemented" << std::endl);
    return 0.0;
}// max

template<class TYPE>
void
HierarchyVecCellDataOpsReal<TYPE>::setRandomValues(
    const int data_id,
    const TYPE& width,
    const TYPE& low,
    const bool interior_only) const
{
    TBOX_ERROR("HierarchyVecCellDataOpsReal<TYPE>::setRandomValues() unimplemented" << std::endl);
    return;
}// setRandomValues

template<class TYPE>
int
HierarchyVecCellDataOpsReal<TYPE>::numberOfEntries(
    const int data_id,
    const bool interior_only) const
{
    TBOX_ERROR("HierarchyVecCellDataOpsReal<TYPE>::numberOfEntries() unimplemented" << std::endl);
    return 0;
}// numberOfEntries

template<class TYPE>
double
HierarchyVecCellDataOpsReal<TYPE>::sumControlVolumes(
    const int data_id,
    const int vol_id) const
{
    TBOX_ERROR("HierarchyVecCellDataOpsReal<TYPE>::sumControlVolumes() unimplemented" << std::endl);
    return 0.0;
}// sumControlVolumes

template<class TYPE>
double
HierarchyVecCellDataOpsReal<TYPE>::L1Norm(
    const int data_id,
    const int vol_id,
    bool local_only) const
{
    TBOX_ERROR("HierarchyVecCellDataOpsReal<TYPE>::L1Norm() unimplemented" << std::endl);
    return 0.0;
}// L1Norm

template<class TYPE>
double
HierarchyVecCellDataOpsReal<TYPE>::L2Norm(
    const int data_id,
    const int vol_id,
    bool local_only) const
{
    TBOX_ERROR("HierarchyVecCellDataOpsReal<TYPE>::L2Norm() unimplemented" << std::endl);
    return 0.0;
}// L2Norm

template<class TYPE>
double
HierarchyVecCellDataOpsReal<TYPE>::weightedL2Norm(
    const int data_id,
    const int wgt_id,
    const int vol_id) const
{
    TBOX_ERROR("HierarchyVecCellDataOpsReal<TYPE>::weightedL2Norm() unimplemented" << std::endl);
    return 0.0;
}// weightedL2Norm

template<class TYPE>
double
HierarchyVecCellDataOpsReal<TYPE>::RMSNorm(
    const int data_id,
    const int vol_id) const
{
    TBOX_ERROR("HierarchyVecCellDataOpsReal<TYPE>::RMSNorm() unimplemented" << std::endl);
    return 0.0;
}// RMSNorm

template<class TYPE>
double
HierarchyVecCellDataOpsReal<TYPE>::weightedRMSNorm(
    const int data_id,
    const int wgt_id,
    const int vol_id) const
{
    TBOX_ERROR("HierarchyVecCellDataOpsReal<TYPE>::weightedRMSNorm() unimplemented" << std::endl);
    return 0.0;
}// weightedRMSNorm

template<class TYPE>
double
HierarchyVecCellDataOpsReal<TYPE>::maxNorm(
    const int data_id,
    const int vol_id,
    bool local_only) const
{
    TBOX_ERROR("HierarchyVecCellDataOpsReal<TYPE>::maxNorm() unimplemented" << std::endl);
    return 0.0;
}// maxNorm

template<class TYPE>
TYPE
HierarchyVecCellDataOpsReal<TYPE>::dot(
    const int data1_id,
    const int data2_id,
    const int vol_id,
    bool local_only) const
{
    TBOX_ERROR("HierarchyVecCellDataOpsReal<TYPE>::dot() unimplemented" << std::endl);
    return 0.0;
}// dot

template<class TYPE>
TYPE
HierarchyVecCellDataOpsReal<TYPE>::integral(
    const int data_id,
    const int vol_id) const
{
    TBOX_ERROR("HierarchyVecCellDataOpsReal<TYPE>::integral() unimplemented" << std::endl);
    return 0.0;
}// integral

template<class TYPE>
int
HierarchyVecCellDataOpsReal<TYPE>::computeConstrProdPos(
    const int data1_id,
    const int data2_id,
    const int vol_id) const
{
    TBOX_ERROR("HierarchyVecCellDataOpsReal<TYPE>::computeConstrProdPos() unimplemented" << std::endl);
    return 0;
}// computeConstrProdPos

template<class TYPE>
void
HierarchyVecCellDataOpsReal<TYPE>::compareToScalar(
    const int dst_id,
    const int src_id,
    const TYPE& alpha,
    const int vol_id) const
{
    TBOX_ERROR("HierarchyVecCellDataOpsReal<TYPE>::compareToScalar() unimplemented" << std::endl);
    return;
}// compareToScalar

template<class TYPE>
int
HierarchyVecCellDataOpsReal<TYPE>::testReciprocal(
    const int dst_id,
    const int src_id,
    const int vol_id) const
{
    TBOX_ERROR("HierarchyVecCellDataOpsReal<TYPE>::testReciprocal() unimplemented" << std::endl);
    return 0;
}// testReciprocal

template<class TYPE>
TYPE
HierarchyVecCellDataOpsReal<TYPE>::maxPointwiseDivide(
    const int numer_id,
    const int denom_id,
    bool local_only) const
{
    TBOX_ERROR("HierarchyVecCellDataOpsReal<TYPE>::maxPointwiseDivide() unimplemented" << std::endl);
    return 0.0;
}// maxPointwiseDivide

template<class TYPE>
TYPE
HierarchyVecCellDataOpsReal<TYPE>::minPointwiseDivide(
    const int numer_id,
    const int denom_id,
    bool local_only) const
{
    TBOX_ERROR("HierarchyVecCellDataOpsReal<TYPE>::minPointwiseDivide() unimplemented" << std::endl);
    return 0.0;
}// minPointwiseDivide

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class IBTK::HierarchyVecCellDataOpsReal<double>;
template class Pointer<IBTK::HierarchyVecCellDataOpsReal<double> >;

//////////////////////////////////////////////////////////////////////////////
