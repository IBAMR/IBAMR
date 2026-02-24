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

#include "ibtk/FaceSynchCopyFillPattern.h"
#include "ibtk/samrai_compatibility_names.h"

#include "SAMRAIBox.h"
#include "SAMRAIBoxGeometry.h"
#include "SAMRAIBoxList.h"
#include "SAMRAIBoxOverlap.h"
#include "SAMRAIFaceGeometry.h"
#include "SAMRAIFaceOverlap.h"
#include "SAMRAIIntVector.h"
#include "SAMRAIPointer.h"

#include <string>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const std::string PATTERN_NAME = "FACE_SYNCH_COPY_FILL_PATTERN";
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

SAMRAIPointer<SAMRAIBoxOverlap>
FaceSynchCopyFillPattern::calculateOverlap(const SAMRAIBoxGeometry& dst_geometry,
                                           const SAMRAIBoxGeometry& src_geometry,
                                           const SAMRAIBox& /*dst_patch_box*/,
                                           const SAMRAIBox& src_mask,
                                           const bool overwrite_interior,
                                           const SAMRAIIntVector& src_offset) const
{
    SAMRAIPointer<SAMRAIFaceOverlap> box_geom_overlap =
        dst_geometry.calculateOverlap(src_geometry, src_mask, overwrite_interior, src_offset);
#if !defined(NDEBUG)
    TBOX_ASSERT(box_geom_overlap);
#endif
    if (box_geom_overlap->isOverlapEmpty()) return box_geom_overlap;

    auto const t_dst_geometry = dynamic_cast<const SAMRAIFaceGeometry*>(&dst_geometry);
#if !defined(NDEBUG)
    TBOX_ASSERT(t_dst_geometry);
#endif
    SAMRAIBoxList dst_boxes[NDIM];
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        bool skip = false;
        for (unsigned int d = 0; d < NDIM && !skip; ++d)
        {
            if (d != axis)
            {
                skip = skip || (src_offset(d) != 0);
            }
        }
        if (!skip)
        {
            // Determine the stencil box.
            const SAMRAIBox& dst_box = t_dst_geometry->getBox();
            SAMRAIBox stencil_box = SAMRAIFaceGeometry::toFaceBox(dst_box, axis);
            stencil_box.lower(0) = stencil_box.upper(0);

            // Intersect the original overlap boxes with the stencil box.
            const SAMRAIBoxList& box_geom_overlap_boxes = box_geom_overlap->getDestinationBoxList(axis);
            for (SAMRAIBoxList::Iterator it(box_geom_overlap_boxes); it; it++)
            {
                const SAMRAIBox overlap_box = stencil_box * it();
                if (!overlap_box.empty()) dst_boxes[axis].appendItem(overlap_box);
            }
        }
    }
    return new SAMRAIFaceOverlap(dst_boxes, src_offset);
} // calculateOverlap

SAMRAIIntVector&
FaceSynchCopyFillPattern::getStencilWidth()
{
    return d_stencil_width;
} // getStencilWidth

const std::string&
FaceSynchCopyFillPattern::getPatternName() const
{
    return PATTERN_NAME;
} // getPatternName

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
