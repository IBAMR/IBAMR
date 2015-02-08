// Filename: NodeSynchCopyFillPattern.cpp
// Created on 02 Feb 2011 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ostream>
#include <string>

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/BoxList.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#include "SAMRAI/pdat/NodeOverlap.h"
#include "ibtk/NodeSynchCopyFillPattern.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const std::string PATTERN_NAME = "NODE_SYNCH_COPY_FILL_PATTERN";
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

NodeSynchCopyFillPattern::NodeSynchCopyFillPattern(const unsigned int axis) : d_stencil_width(DIM, 1), d_axis(axis)
{
    // intentionally blank
    return;
}

NodeSynchCopyFillPattern::~NodeSynchCopyFillPattern()
{
    // intentionally blank
    return;
}

Pointer<BoxOverlap> NodeSynchCopyFillPattern::calculateOverlap(const BoxGeometry& dst_geometry,
                                                               const BoxGeometry& src_geometry,
                                                               const Box& /*dst_patch_box*/,
                                                               const Box& src_mask,
                                                               const bool overwrite_interior,
                                                               const IntVector& src_offset) const
{
    Pointer<NodeOverlap> box_geom_overlap =
        dst_geometry.calculateOverlap(src_geometry, src_mask, overwrite_interior, src_offset);
    TBOX_ASSERT(box_geom_overlap);
    if (box_geom_overlap->isOverlapEmpty()) return box_geom_overlap;

    const NodeGeometry* const t_dst_geometry = dynamic_cast<const NodeGeometry*>(&dst_geometry);
    TBOX_ASSERT(t_dst_geometry);
    BoxList dst_boxes;
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
        const Box& dst_box = t_dst_geometry->getBox();
        Box stencil_box = NodeGeometry::toNodeBox(dst_box);
        stencil_box.lower()(d_axis) = stencil_box.upper()(d_axis);

        // Intersect the original overlap boxes with the stencil box.
        const BoxList& box_geom_overlap_boxes = box_geom_overlap->getDestinationBoxList();
        for (BoxList::Iterator it(box_geom_overlap_boxes); it; it++)
        {
            const Box overlap_box = stencil_box * it();
            if (!overlap_box.empty()) dst_boxes.appendItem(overlap_box);
        }
    }
    return Pointer<BoxOverlap>(new NodeOverlap(dst_boxes, src_offset));
}

IntVector& NodeSynchCopyFillPattern::getStencilWidth()
{
    return d_stencil_width;
}

const std::string& NodeSynchCopyFillPattern::getPatternName() const
{
    return PATTERN_NAME;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////
}
//////////////////////////////////////////////////////////////////////////////
