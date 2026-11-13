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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/RegionFillPattern.h>

#include <CellOverlap.h>
#include <SideOverlap.h>

#include <array>
#include <utility>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

RegionFillPattern::RegionFillPattern(Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> inner,
                                     std::vector<BoxList<NDIM>> regions)
    : d_inner(inner), d_regions(std::move(regions))
{
    // intentionally blank
}

Pointer<BoxOverlap<NDIM>>
RegionFillPattern::calculateOverlap(const BoxGeometry<NDIM>& dst_geometry,
                                    const BoxGeometry<NDIM>& src_geometry,
                                    const Box<NDIM>& dst_patch_box,
                                    const Box<NDIM>& src_mask,
                                    const bool overwrite_interior,
                                    const IntVector<NDIM>& src_offset) const
{
    return restrictOverlap(
        d_inner->calculateOverlap(dst_geometry, src_geometry, dst_patch_box, src_mask, overwrite_interior, src_offset),
        src_offset);
}

Pointer<BoxOverlap<NDIM>>
RegionFillPattern::calculateOverlapOnLevel(const BoxGeometry<NDIM>& dst_geometry,
                                           const BoxGeometry<NDIM>& src_geometry,
                                           const Box<NDIM>& dst_patch_box,
                                           const Box<NDIM>& src_mask,
                                           const bool overwrite_interior,
                                           const IntVector<NDIM>& src_offset,
                                           const int dst_level_num,
                                           const int src_level_num) const
{
    return restrictOverlap(d_inner->calculateOverlapOnLevel(dst_geometry,
                                                            src_geometry,
                                                            dst_patch_box,
                                                            src_mask,
                                                            overwrite_interior,
                                                            src_offset,
                                                            dst_level_num,
                                                            src_level_num),
                           src_offset);
}

void
RegionFillPattern::setTargetPatchLevelNumber(const int level_num)
{
    SAMRAI::xfer::VariableFillPattern<NDIM>::setTargetPatchLevelNumber(level_num);
    d_inner->setTargetPatchLevelNumber(level_num);
}

IntVector<NDIM>&
RegionFillPattern::getStencilWidth()
{
    return d_inner->getStencilWidth();
}

const std::string&
RegionFillPattern::getPatternName() const
{
    static const std::string name = "REGION_FILL_PATTERN";
    return name;
}

/////////////////////////////// PRIVATE //////////////////////////////////////

Pointer<BoxOverlap<NDIM>>
RegionFillPattern::restrictOverlap(Pointer<BoxOverlap<NDIM>> overlap, const IntVector<NDIM>& src_offset) const
{
    Pointer<CellOverlap<NDIM>> cell_overlap = overlap;
    if (cell_overlap)
    {
        BoxList<NDIM> boxes = cell_overlap->getDestinationBoxList();
        boxes.intersectBoxes(d_regions[0]);
        return new CellOverlap<NDIM>(boxes, src_offset);
    }
    Pointer<SideOverlap<NDIM>> side_overlap = overlap;
    std::array<BoxList<NDIM>, NDIM> boxes;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        boxes[axis] = side_overlap->getDestinationBoxList(axis);
        boxes[axis].intersectBoxes(d_regions[axis]);
    }
    return new SideOverlap<NDIM>(boxes.data(), src_offset);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
