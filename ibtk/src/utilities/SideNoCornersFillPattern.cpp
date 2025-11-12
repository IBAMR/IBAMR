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

#include "ibtk/SideNoCornersFillPattern.h"

#include "Box.h"
#include "BoxGeometry.h"
#include "BoxList.h"
#include "BoxOverlap.h"
#include "IntVector.h"
#include "SideGeometry.h"
#include "SideOverlap.h"
#include "tbox/Pointer.h"

#include <array>

#include "ibtk/namespaces.h" // IWYU pragma: keep

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

Pointer<BoxOverlap<NDIM> >
SideNoCornersFillPattern::calculateOverlap(const BoxGeometry<NDIM>& dst_geometry,
                                           const BoxGeometry<NDIM>& src_geometry,
                                           const Box<NDIM>& /*dst_patch_box*/,
                                           const Box<NDIM>& src_mask,
                                           const bool overwrite_interior,
                                           const IntVector<NDIM>& src_offset) const
{
    return dst_geometry.calculateOverlap(src_geometry, src_mask, overwrite_interior, src_offset);
} // calculateOverlap

Pointer<BoxOverlap<NDIM> >
SideNoCornersFillPattern::calculateOverlapOnLevel(const BoxGeometry<NDIM>& dst_geometry,
                                                  const BoxGeometry<NDIM>& src_geometry,
                                                  const Box<NDIM>& dst_patch_box,
                                                  const Box<NDIM>& src_mask,
                                                  const bool overwrite_interior,
                                                  const IntVector<NDIM>& src_offset,
                                                  const int dst_level_num,
                                                  const int /*src_level_num*/) const
{
    if (d_target_level_num == dst_level_num)
    {
        Pointer<SideOverlap<NDIM> > default_overlap = dst_geometry.calculateOverlap(
            dst_geometry, src_geometry, src_mask, overwrite_interior, src_offset, /*retry*/ false);

        // If we are on the target level, we set up the stencil to include overlaps
        // only if this->d_overwrite_interior == true.
        std::array<BoxList<NDIM>, NDIM> dst_boxes;
        for (int axis = 0; axis < NDIM; ++axis)
        {
            dst_boxes[axis] = default_overlap->getDestinationBoxList(axis);
            BoxList<NDIM> stencil_boxes;
            auto dst_side_box = SideGeometry<NDIM>::toSideBox(dst_patch_box, axis);
            if (d_overwrite_interior) stencil_boxes.addItem(dst_side_box);
            for (int d = 0; d < NDIM; ++d)
            {
                Box<NDIM> lower_box = dst_side_box;
                lower_box.lower(d) = dst_side_box.lower(d) - d_stencil_width(d);
                lower_box.upper(d) = dst_side_box.lower(d) - 1;
                stencil_boxes.addItem(lower_box);

                Box<NDIM> upper_box = dst_side_box;
                upper_box.lower(d) = dst_side_box.upper(d) + 1;
                upper_box.upper(d) = dst_side_box.upper(d) + d_stencil_width(d);
                stencil_boxes.addItem(upper_box);
            }
            dst_boxes[axis].intersectBoxes(stencil_boxes);
        }
        Pointer<SideOverlap<NDIM> > overlap = new SideOverlap<NDIM>(dst_boxes.data(), src_offset);
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

IntVector<NDIM>&
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
