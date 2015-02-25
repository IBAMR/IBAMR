// Filename: SideSynchCopyFillPattern.cpp
// Created on 10 Mar 2010 by Boyce Griffith
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
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideOverlap.h"
#include "ibtk/SideSynchCopyFillPattern.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const std::string PATTERN_NAME = "SIDE_SYNCH_COPY_FILL_PATTERN";
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

SideSynchCopyFillPattern::SideSynchCopyFillPattern() : d_stencil_width(DIM, 1)
{
    // intentionally blank
    return;
} // SideSynchCopyFillPattern

SideSynchCopyFillPattern::~SideSynchCopyFillPattern()
{
    // intentionally blank
    return;
} // SideSynchCopyFillPattern

boost::shared_ptr<BoxOverlap> SideSynchCopyFillPattern::calculateOverlap(const BoxGeometry& dst_geometry,
                                                                         const BoxGeometry& src_geometry,
                                                                         const Box& dst_patch_box,
                                                                         const Box& src_mask,
                                                                         const Box& fill_box,
                                                                         const bool overwrite_interior,
                                                                         const Transformation& transformation) const
{
    auto box_geom_overlap = BOOST_CAST<SideOverlap>(
        dst_geometry.calculateOverlap(src_geometry, src_mask, fill_box, overwrite_interior, transformation));
    TBOX_ASSERT(box_geom_overlap);
    if (box_geom_overlap->isOverlapEmpty()) return box_geom_overlap;

    auto t_dst_geometry = CPP_CAST<const SideGeometry*>(&dst_geometry);
    TBOX_ASSERT(t_dst_geometry);

    std::vector<BoxContainer> dst_boxes(NDIM);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        // Determine the stencil box.
        const Box& dst_box = t_dst_geometry->getBox();
        Box stencil_box = SideGeometry::toSideBox(dst_box, axis);
        stencil_box.setLower(axis, stencil_box.upper(axis));

        // Intersect the original overlap boxes with the stencil box.
        const BoxContainer& box_geom_overlap_boxes = box_geom_overlap->getDestinationBoxContainer(axis);
        for (BoxContainer::const_iterator it = box_geom_overlap_boxes.begin(); it != box_geom_overlap_boxes.end(); ++it)
        {
            const Box overlap_box(stencil_box * *it);
            if (!overlap_box.empty()) dst_boxes[axis].push_back(overlap_box);
        }
    }
    return boost::make_shared<SideOverlap>(dst_boxes, src_offset);
} // calculateOverlap

IntVector& SideSynchCopyFillPattern::getStencilWidth()
{
    return d_stencil_width;
} // getStencilWidth

const std::string& SideSynchCopyFillPattern::getPatternName() const
{
    return PATTERN_NAME;
} // getPatternName

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
