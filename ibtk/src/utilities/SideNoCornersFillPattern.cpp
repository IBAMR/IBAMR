// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2025 by the IBAMR developers
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

#include <ibtk/SideNoCornersFillPattern.h>
#include <ibtk/samrai_compatibility_names.h>

#include <SAMRAIBox.h>
#include <SAMRAIBoxGeometry.h>
#include <SAMRAIBoxList.h>
#include <SAMRAIBoxOverlap.h>
#include <SAMRAIIntVector.h>
#include <SAMRAIPointer.h>
#include <SAMRAISideGeometry.h>
#include <SAMRAISideOverlap.h>

#include <array>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const std::string PATTERN_NAME = "SIDE_NO_CORNERS_FILL_PATTERN";
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

SideNoCornersFillPattern::SideNoCornersFillPattern(const int stencil_width, const bool overwrite_interior)
    : d_stencil_width(stencil_width), d_overwrite_interior(overwrite_interior)
{
    // intentionally blank
    return;
} // SideNoCornersFillPattern

SAMRAIPointer<SAMRAIBoxOverlap>
SideNoCornersFillPattern::calculateOverlap(const SAMRAIBoxGeometry& dst_geometry,
                                           const SAMRAIBoxGeometry& src_geometry,
                                           const SAMRAIBox& /*dst_patch_box*/,
                                           const SAMRAIBox& src_mask,
                                           const bool overwrite_interior,
                                           const SAMRAIIntVector& src_offset) const
{
    return dst_geometry.calculateOverlap(src_geometry, src_mask, overwrite_interior, src_offset);
} // calculateOverlap

SAMRAIPointer<SAMRAIBoxOverlap>
SideNoCornersFillPattern::calculateOverlapOnLevel(const SAMRAIBoxGeometry& dst_geometry,
                                                  const SAMRAIBoxGeometry& src_geometry,
                                                  const SAMRAIBox& dst_patch_box,
                                                  const SAMRAIBox& src_mask,
                                                  const bool overwrite_interior,
                                                  const SAMRAIIntVector& src_offset,
                                                  const int dst_level_num,
                                                  const int /*src_level_num*/) const
{
    if (d_target_level_num == dst_level_num)
    {
        SAMRAIPointer<SAMRAISideOverlap> default_overlap = dst_geometry.calculateOverlap(
            dst_geometry, src_geometry, src_mask, overwrite_interior, src_offset, /*retry*/ false);

        // If we are on the target level, we set up the stencil to include overlaps
        // only if this->d_overwrite_interior == true.
        std::array<SAMRAIBoxList, NDIM> dst_boxes;
        for (int axis = 0; axis < NDIM; ++axis)
        {
            dst_boxes[axis] = default_overlap->getDestinationBoxList(axis);
            SAMRAIBoxList stencil_boxes;
            auto dst_side_box = SAMRAISideGeometry::toSideBox(dst_patch_box, axis);
            if (d_overwrite_interior) stencil_boxes.addItem(dst_side_box);
            for (int d = 0; d < NDIM; ++d)
            {
                SAMRAIBox lower_box = dst_side_box;
                lower_box.lower(d) = dst_side_box.lower(d) - d_stencil_width(d);
                lower_box.upper(d) = dst_side_box.lower(d) - 1;
                stencil_boxes.addItem(lower_box);

                SAMRAIBox upper_box = dst_side_box;
                upper_box.lower(d) = dst_side_box.upper(d) + 1;
                upper_box.upper(d) = dst_side_box.upper(d) + d_stencil_width(d);
                stencil_boxes.addItem(upper_box);
            }
            dst_boxes[axis].intersectBoxes(stencil_boxes);
        }
        SAMRAIPointer<SAMRAISideOverlap> overlap = new SAMRAISideOverlap(dst_boxes.data(), src_offset);
        return overlap;
    }
    else
    {
        return calculateOverlap(dst_geometry, src_geometry, dst_patch_box, src_mask, overwrite_interior, src_offset);
    }
} // calculateOverlapOnLevel

void
SideNoCornersFillPattern::setTargetPatchLevelNumber(const int level_num)
{
    d_target_level_num = level_num;
    return;
} // setTargetPatchLevelNumber

SAMRAIIntVector&
SideNoCornersFillPattern::getStencilWidth()
{
    return d_stencil_width;
} // getStencilWidth

const std::string&
SideNoCornersFillPattern::getPatternName() const
{
    return PATTERN_NAME;
} // getPatternName

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
