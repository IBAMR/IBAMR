// Filename: LNodeIndexData.C
// Last modified: <17.Apr.2007 18:19:51 griffith@box221.cims.nyu.edu>
// Created on 01 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

#include "LNodeIndexData.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/LNodeIndexVariable.h>

// SAMRAI INCLUDES
#include <CellOverlap.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LNodeIndexData::LNodeIndexData(
    const SAMRAI::hier::Box<NDIM>& box,
    const SAMRAI::hier::IntVector<NDIM>& ghosts)
    : SAMRAI::pdat::IndexData<NDIM,LNodeIndexSet>(box,ghosts),
      d_interior_local_indices(),
      d_ghost_local_indices()
{
    // intentionally blank
    return;
}// LNodeIndexData

LNodeIndexData::~LNodeIndexData()
{
    // intentionally blank
    return;
}// ~LNodeIndexData

void
LNodeIndexData::copy(
    const SAMRAI::hier::PatchData<NDIM>& src)
{
    LNodeIndexData* t_src = const_cast<LNodeIndexData*>(dynamic_cast<const LNodeIndexData*>(&src));
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(t_src != NULL);
#endif
    const SAMRAI::hier::Box<NDIM>& ghost_box = getGhostBox();
    const SAMRAI::hier::Box<NDIM>& src_ghost_box = t_src->getGhostBox();
    removeInsideBox(ghost_box*src_ghost_box);
    for (LNodeIndexData::Iterator it(*t_src); it; it++)
    {
        const SAMRAI::hier::Index<NDIM>& i = it.getIndex();
        const LNodeIndexSet& src_node_set = *it;
        if (ghost_box.contains(i)) appendItem(i,src_node_set);
    }
}// copy

void
LNodeIndexData::copy2(
    SAMRAI::hier::PatchData<NDIM>& dst) const
{
    dst.copy(*this);
}// copy

void
LNodeIndexData::copy(
    const SAMRAI::hier::PatchData<NDIM>& src,
    const SAMRAI::hier::BoxOverlap<NDIM>& overlap)
{
    LNodeIndexData* t_src = const_cast<LNodeIndexData*>(dynamic_cast<const LNodeIndexData*>(&src));
    const SAMRAI::pdat::CellOverlap<NDIM>* t_overlap = dynamic_cast<const SAMRAI::pdat::CellOverlap<NDIM>*>(&overlap);
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(t_src != NULL);
    assert(t_overlap != NULL);
#endif
    const SAMRAI::hier::Box<NDIM>& ghost_box = getGhostBox();
    const SAMRAI::hier::IntVector<NDIM>& src_offset = t_overlap->getSourceOffset();
    const SAMRAI::hier::BoxList<NDIM>& dst_box_list = t_overlap->getDestinationBoxList();

    for (SAMRAI::hier::BoxList<NDIM>::Iterator b(dst_box_list); b; b++)
    {
        const SAMRAI::hier::Box<NDIM>& dst_box = b();
        removeInsideBox(ghost_box*dst_box);
    }

    for (SAMRAI::hier::BoxList<NDIM>::Iterator b(dst_box_list); b; b++)
    {
        const SAMRAI::hier::Box<NDIM>& dst_box = b();
        const SAMRAI::hier::Box<NDIM>  src_box(SAMRAI::hier::Box<NDIM>::shift(dst_box, -src_offset));
        for (LNodeIndexData::Iterator it(*t_src); it; it++)
        {
            const SAMRAI::hier::Index<NDIM>& i_src = it.getIndex();
            const SAMRAI::hier::Index<NDIM>& i_dst = i_src+src_offset;
            const LNodeIndexSet& src_node_set = *it;
            if (src_box.contains(i_src) && ghost_box.contains(i_dst))
            {
                LNodeIndexSet new_item;
                new_item.copySourceItem(i_src, src_offset, src_node_set);
                appendItem(i_dst, new_item);
            }
        }
    }
}// copy

void
LNodeIndexData::copy2(
    SAMRAI::hier::PatchData<NDIM>& dst,
    const SAMRAI::hier::BoxOverlap<NDIM>& overlap) const
{
    dst.copy(*this, overlap);
}// copy2

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::LNodeIndexData>;

//////////////////////////////////////////////////////////////////////////////
