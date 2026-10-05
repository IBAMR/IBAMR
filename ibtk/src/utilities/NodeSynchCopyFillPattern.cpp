// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
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

#include <ibtk/NodeSynchCopyFillPattern.h>

#include <tbox/Pointer.h>

#include <Box.h>
#include <BoxGeometry.h>
#include <NodeGeometry.h>
#include <NodeOverlap.h>

#include <algorithm>
#include <string>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const std::string PATTERN_NAME = "NODE_SYNCH_COPY_FILL_PATTERN";
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

NodeSynchCopyFillPattern::NodeSynchCopyFillPattern(const unsigned int axis) : d_axis(axis)
{
    // intentionally blank
    return;
} // NodeSynchCopyFillPattern

Pointer<BoxOverlap<NDIM>>
NodeSynchCopyFillPattern::calculateOverlap(const BoxGeometry<NDIM>& dst_geometry,
                                           const BoxGeometry<NDIM>& src_geometry,
                                           const Box<NDIM>& /*dst_patch_box*/,
                                           const Box<NDIM>& src_mask,
                                           const bool overwrite_interior,
                                           const IntVector<NDIM>& src_offset) const
{
    Pointer<NodeOverlap<NDIM>> box_geom_overlap =
        dst_geometry.calculateOverlap(src_geometry, src_mask, overwrite_interior, src_offset);
#if !defined(NDEBUG)
    TBOX_ASSERT(box_geom_overlap);
#endif
    if (box_geom_overlap->isOverlapEmpty()) return box_geom_overlap;

    auto const t_dst_geometry = dynamic_cast<const NodeGeometry<NDIM>*>(&dst_geometry);
    auto const t_src_geometry = dynamic_cast<const NodeGeometry<NDIM>*>(&src_geometry);
#if !defined(NDEBUG)
    TBOX_ASSERT(t_dst_geometry);
    TBOX_ASSERT(t_src_geometry);
#endif
    const Box<NDIM>& dst_box = t_dst_geometry->getBox();
    const Box<NDIM> src_box = Box<NDIM>::shift(t_src_geometry->getBox(), src_offset);

    // The copy of a shared node that a box holds belongs to the highest cell of
    // the box that touches the node.  In the pass for direction d_axis, a copy
    // whose cell is on the lower side of the node in that direction is replaced
    // by a copy whose cell is on the upper side in that direction and on the
    // same side in each later direction.  The passes for the earlier
    // directions have made all such copies equal, so the result does not
    // depend on which of them is used, and after the last pass every copy has
    // the value of the copy whose cell is the highest, comparing the last
    // coordinate first.
    BoxList<NDIM> dst_boxes;
    if (src_box.upper(d_axis) > dst_box.upper(d_axis))
    {
        // Determine the stencil box.
        Box<NDIM> stencil_box = NodeGeometry<NDIM>::toNodeBox(dst_box);
        stencil_box.lower(d_axis) = stencil_box.upper(d_axis);
        for (unsigned int d = d_axis + 1; d < NDIM; ++d)
        {
            if (src_box.upper(d) != dst_box.upper(d))
            {
                stencil_box.upper(d) = std::min(src_box.upper(d), dst_box.upper(d));
            }
        }

        // Intersect the original overlap boxes with the stencil box.
        const BoxList<NDIM>& box_geom_overlap_boxes = box_geom_overlap->getDestinationBoxList();
        for (BoxList<NDIM>::Iterator it(box_geom_overlap_boxes); it; it++)
        {
            const Box<NDIM> overlap_box = stencil_box * it();
            if (!overlap_box.empty()) dst_boxes.appendItem(overlap_box);
        }
    }
    return new NodeOverlap<NDIM>(dst_boxes, src_offset);
} // calculateOverlap

IntVector<NDIM>&
NodeSynchCopyFillPattern::getStencilWidth()
{
    return d_stencil_width;
} // getStencilWidth

const std::string&
NodeSynchCopyFillPattern::getPatternName() const
{
    return PATTERN_NAME;
} // getPatternName

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
