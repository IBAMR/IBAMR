// Filename: SideNoCornersFillPattern.cpp
// Created on 09 Mar 2010 by Boyce Griffith
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
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideOverlap.h"
#include "boost/array.hpp"
#include "ibtk/SideNoCornersFillPattern.h"
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
static const std::string PATTERN_NAME = "SIDE_NO_CORNERS_FILL_PATTERN";
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

SideNoCornersFillPattern::SideNoCornersFillPattern(const int stencil_width,
                                                   const bool include_dst_patch_box,
                                                   const bool include_edges_on_dst_level,
                                                   const bool include_edges_on_src_level)
    : d_stencil_width(DIM, stencil_width), d_include_dst_patch_box(include_dst_patch_box),
      d_include_edges_on_dst_level(include_edges_on_dst_level),
      d_include_edges_on_src_level(include_edges_on_src_level), d_target_level_num(-1)
{
    // intentionally blank
    return;
} // SideNoCornersFillPattern

SideNoCornersFillPattern::~SideNoCornersFillPattern()
{
    // intentionally blank
    return;
} // SideNoCornersFillPattern

Pointer<BoxOverlap> SideNoCornersFillPattern::calculateOverlap(const BoxGeometry& dst_geometry,
                                                               const BoxGeometry& src_geometry,
                                                               const Box& /*dst_patch_box*/,
                                                               const Box& src_mask,
                                                               const bool overwrite_interior,
                                                               const IntVector& src_offset) const
{
    Pointer<SideOverlap> box_geom_overlap =
        dst_geometry.calculateOverlap(src_geometry, src_mask, overwrite_interior, src_offset);
    TBOX_ASSERT(box_geom_overlap);
    const SideGeometry* const t_dst_geometry = dynamic_cast<const SideGeometry*>(&dst_geometry);
    TBOX_ASSERT(t_dst_geometry);
    Array<BoxList> dst_boxes(NDIM);
    if (!box_geom_overlap->isOverlapEmpty())
    {
        const Box& dst_box = t_dst_geometry->getBox();
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            const BoxList& box_geom_overlap_boxes = box_geom_overlap->getDestinationBoxList(axis);

            // Determine the stencil boxes with the specified ghost cell width.
            BoxList stencil_boxes;
            if (NDIM == 2 || (!d_include_edges_on_src_level && !d_include_edges_on_dst_level))
            {
                for (unsigned int i = 0; i < NDIM; ++i)
                {
                    Box box = dst_box;
                    box.lower()(i) -= d_stencil_width(i);
                    box.upper()(i) += d_stencil_width(i);
                    stencil_boxes.appendItem(SideGeometry::toSideBox(box, axis));
                }
            }
            else
            {
                for (unsigned int j = 0; j < NDIM; ++j)
                {
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        if (i == j) continue;
                        Box box = dst_box;
                        box.lower()(i) -= d_stencil_width(i);
                        box.upper()(i) += d_stencil_width(i);
                        box.lower()(j) -= d_stencil_width(j);
                        box.upper()(j) += d_stencil_width(j);
                        stencil_boxes.appendItem(SideGeometry::toSideBox(box, axis));
                    }
                }
            }

            // Intersect the overlap boxes with the stencil boxes.
            for (BoxList::Iterator it1(box_geom_overlap_boxes); it1; it1++)
            {
                BoxList overlap_boxes(stencil_boxes);
                overlap_boxes.intersectBoxes(it1());
                for (BoxList::Iterator it2(overlap_boxes); it2; it2++)
                {
                    const Box& overlap_box = it2();
                    if (!overlap_box.empty()) dst_boxes[axis].appendItem(overlap_box);
                }
            }
        }
    }
    return Pointer<BoxOverlap>(new SideOverlap(dst_boxes, src_offset));
} // calculateOverlap

void SideNoCornersFillPattern::setTargetPatchLevelNumber(const int level_num)
{
    d_target_level_num = level_num;
    return;
} // setTargetPatchLevelNumber

IntVector& SideNoCornersFillPattern::getStencilWidth()
{
    return d_stencil_width;
} // getStencilWidth

const std::string& SideNoCornersFillPattern::getPatternName() const
{
    return PATTERN_NAME;
} // getPatternName

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
