// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#include "ibtk/NodeSynchCopyFillPattern.h"

#include "Box.h"
#include "BoxGeometry.h"
#include "NodeGeometry.h"
#include "NodeOverlap.h"
#include "tbox/Pointer.h"

#include <string>

#include "ibtk/namespaces.h" // IWYU pragma: keep

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

Pointer<BoxOverlapNd>
NodeSynchCopyFillPattern::calculateOverlap(const BoxGeometryNd& dst_geometry,
                                           const BoxGeometryNd& src_geometry,
                                           const BoxNd& /*dst_patch_box*/,
                                           const BoxNd& src_mask,
                                           const bool overwrite_interior,
                                           const IntVectorNd& src_offset) const
{
    Pointer<NodeOverlapNd> box_geom_overlap =
        dst_geometry.calculateOverlap(src_geometry, src_mask, overwrite_interior, src_offset);
#if !defined(NDEBUG)
    TBOX_ASSERT(box_geom_overlap);
#endif
    if (box_geom_overlap->isOverlapEmpty()) return box_geom_overlap;

    auto const t_dst_geometry = dynamic_cast<const NodeGeometryNd*>(&dst_geometry);
#if !defined(NDEBUG)
    TBOX_ASSERT(t_dst_geometry);
#endif
    BoxListNd dst_boxes;
    bool skip = false;
    for (unsigned int d = 0; d < NDIM && !skip; ++d)
    {
        if (d != d_axis)
        {
            skip = skip || (src_offset(d) != 0);
        }
    }
    if (!skip)
    {
        // Determine the stencil box.
        const BoxNd& dst_box = t_dst_geometry->getBox();
        BoxNd stencil_box = NodeGeometryNd::toNodeBox(dst_box);
        stencil_box.lower(d_axis) = stencil_box.upper(d_axis);

        // Intersect the original overlap boxes with the stencil box.
        const BoxListNd& box_geom_overlap_boxes = box_geom_overlap->getDestinationBoxList();
        for (BoxListNd::Iterator it(box_geom_overlap_boxes); it; it++)
        {
            const BoxNd overlap_box = stencil_box * it();
            if (!overlap_box.empty()) dst_boxes.appendItem(overlap_box);
        }
    }
    return new NodeOverlapNd(dst_boxes, src_offset);
} // calculateOverlap

IntVectorNd&
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
