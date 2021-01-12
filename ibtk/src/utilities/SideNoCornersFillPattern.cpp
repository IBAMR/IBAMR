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

#include "ibtk/SideNoCornersFillPattern.h"

#include "Box.h"
#include "BoxGeometry.h"
#include "BoxList.h"
#include "BoxOverlap.h"
#include "Index.h"
#include "IntVector.h"
#include "SideGeometry.h"
#include "SideOverlap.h"
#include "tbox/Pointer.h"

#include <array>
#include <string>

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

SideNoCornersFillPattern::SideNoCornersFillPattern(const int stencil_width,
                                                   const bool include_dst_patch_box,
                                                   const bool include_edges_on_dst_level,
                                                   const bool include_edges_on_src_level)
    : d_stencil_width(stencil_width),
      d_include_dst_patch_box(include_dst_patch_box),
      d_include_edges_on_dst_level(include_edges_on_dst_level),
      d_include_edges_on_src_level(include_edges_on_src_level)
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
    Pointer<SideOverlap<NDIM> > box_geom_overlap =
        dst_geometry.calculateOverlap(src_geometry, src_mask, overwrite_interior, src_offset);
#if !defined(NDEBUG)
    TBOX_ASSERT(box_geom_overlap);
#endif
    auto const t_dst_geometry = dynamic_cast<const SideGeometry<NDIM>*>(&dst_geometry);
#if !defined(NDEBUG)
    TBOX_ASSERT(t_dst_geometry);
#endif
    std::array<BoxList<NDIM>, NDIM> dst_boxes;
    if (!box_geom_overlap->isOverlapEmpty())
    {
        const Box<NDIM>& dst_box = t_dst_geometry->getBox();
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            const BoxList<NDIM>& box_geom_overlap_boxes = box_geom_overlap->getDestinationBoxList(axis);

            // Determine the stencil boxes with the specified ghost cell width.
            BoxList<NDIM> stencil_boxes;
            if (NDIM == 2 || (!d_include_edges_on_src_level && !d_include_edges_on_dst_level))
            {
                for (unsigned int i = 0; i < NDIM; ++i)
                {
                    Box<NDIM> box = dst_box;
                    box.lower()(i) -= d_stencil_width(i);
                    box.upper()(i) += d_stencil_width(i);
                    stencil_boxes.appendItem(SideGeometry<NDIM>::toSideBox(box, axis));
                }
            }
            else
            {
                for (unsigned int j = 0; j < NDIM; ++j)
                {
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        if (i == j) continue;
                        Box<NDIM> box = dst_box;
                        box.lower()(i) -= d_stencil_width(i);
                        box.upper()(i) += d_stencil_width(i);
                        box.lower()(j) -= d_stencil_width(j);
                        box.upper()(j) += d_stencil_width(j);
                        stencil_boxes.appendItem(SideGeometry<NDIM>::toSideBox(box, axis));
                    }
                }
            }

            // Intersect the overlap boxes with the stencil boxes.
            for (BoxList<NDIM>::Iterator it1(box_geom_overlap_boxes); it1; it1++)
            {
                BoxList<NDIM> overlap_boxes(stencil_boxes);
                overlap_boxes.intersectBoxes(it1());
                for (BoxList<NDIM>::Iterator it2(overlap_boxes); it2; it2++)
                {
                    const Box<NDIM>& overlap_box = it2();
                    if (!overlap_box.empty()) dst_boxes[axis].appendItem(overlap_box);
                }
            }
        }
    }
    return new SideOverlap<NDIM>(dst_boxes.data(), src_offset);
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
    Pointer<SideOverlap<NDIM> > box_geom_overlap =
        dst_geometry.calculateOverlap(src_geometry, src_mask, overwrite_interior, src_offset);
#if !defined(NDEBUG)
    TBOX_ASSERT(box_geom_overlap);
#endif
    auto const t_dst_geometry = dynamic_cast<const SideGeometry<NDIM>*>(&dst_geometry);
#if !defined(NDEBUG)
    TBOX_ASSERT(t_dst_geometry);
#endif
    std::array<BoxList<NDIM>, NDIM> dst_boxes;
    if (!box_geom_overlap->isOverlapEmpty())
    {
        const Box<NDIM>& dst_box = t_dst_geometry->getBox();
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            const BoxList<NDIM>& box_geom_overlap_boxes = box_geom_overlap->getDestinationBoxList(axis);

            // Determine the stencil boxes with the specified ghost cell width.
            BoxList<NDIM> stencil_boxes;
            if (NDIM == 2 || (!d_include_edges_on_dst_level && dst_level_num == d_target_level_num) ||
                (!d_include_edges_on_src_level && dst_level_num != d_target_level_num))
            {
                for (unsigned int i = 0; i < NDIM; ++i)
                {
                    Box<NDIM> box = dst_box;
                    box.lower()(i) -= d_stencil_width(i);
                    box.upper()(i) += d_stencil_width(i);
                    stencil_boxes.appendItem(SideGeometry<NDIM>::toSideBox(box, axis));
                }
            }
            else
            {
                for (unsigned int j = 0; j < NDIM; ++j)
                {
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        if (i == j) continue;
                        Box<NDIM> box = dst_box;
                        box.lower()(i) -= d_stencil_width(i);
                        box.upper()(i) += d_stencil_width(i);
                        box.lower()(j) -= d_stencil_width(j);
                        box.upper()(j) += d_stencil_width(j);
                        stencil_boxes.appendItem(SideGeometry<NDIM>::toSideBox(box, axis));
                    }
                }
            }

            // Intersect the overlap boxes with the stencil boxes.
            for (BoxList<NDIM>::Iterator it1(box_geom_overlap_boxes); it1; it1++)
            {
                BoxList<NDIM> overlap_boxes(stencil_boxes);
                overlap_boxes.intersectBoxes(it1());
                if (dst_level_num == d_target_level_num && !d_include_dst_patch_box)
                {
                    overlap_boxes.removeIntersections(SideGeometry<NDIM>::toSideBox(dst_patch_box, axis));
                }
                for (BoxList<NDIM>::Iterator it2(overlap_boxes); it2; it2++)
                {
                    const Box<NDIM>& overlap_box = it2();
                    if (!overlap_box.empty()) dst_boxes[axis].appendItem(overlap_box);
                }
            }
        }
    }
    return new SideOverlap<NDIM>(dst_boxes.data(), src_offset);
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
