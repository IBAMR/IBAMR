// Filename: CellNoCornersFillPattern.cpp
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
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/pdat/CellOverlap.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "ibtk/CellNoCornersFillPattern.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const std::string PATTERN_NAME = "CELL_NO_CORNERS_FILL_PATTERN";
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CellNoCornersFillPattern::CellNoCornersFillPattern(const int stencil_width,
                                                   const bool include_dst_patch_box,
                                                   const bool include_edges_on_dst_level,
                                                   const bool include_edges_on_src_level)
    : d_stencil_width(DIM, stencil_width), d_include_dst_patch_box(include_dst_patch_box),
      d_include_edges_on_dst_level(include_edges_on_dst_level), d_include_edges_on_src_level(include_edges_on_src_level)
{
    // intentionally blank
    return;
}

CellNoCornersFillPattern::~CellNoCornersFillPattern()
{
    // intentionally blank
    return;
}

boost::shared_ptr<BoxOverlap> CellNoCornersFillPattern::calculateOverlap(const BoxGeometry& dst_geometry,
                                                                         const BoxGeometry& src_geometry,
                                                                         const Box& dst_patch_box,
                                                                         const Box& src_mask,
                                                                         const Box& fill_box,
                                                                         const bool overwrite_interior,
                                                                         const Transformation& transformation) const
{
    auto box_geom_overlap = BOOST_CAST<CellOverlap>(
        dst_geometry.calculateOverlap(src_geometry, src_mask, fill_box, overwrite_interior, transformation));
    TBOX_ASSERT(box_geom_overlap);
    auto t_dst_geometry = CPP_CAST<const CellGeometry*>(&dst_geometry);
    TBOX_ASSERT(t_dst_geometry);

    BoxContainer dst_boxes;
    if (!box_geom_overlap->isOverlapEmpty())
    {
        const Box& dst_box = t_dst_geometry->getBox();
        const BoxContainer& box_geom_overlap_boxes = box_geom_overlap->getDestinationBoxContainer();

        // Determine the stencil boxes with the specified ghost cell width.
        BoxContainer stencil_boxes;
        if (NDIM == 2 || (!d_include_edges_on_src_level && !d_include_edges_on_dst_level))
        {
            for (unsigned int i = 0; i < NDIM; ++i)
            {
                Box box(dst_box);
                box.setLower(i, box.lower(i) - d_stencil_width(i));
                box.setUpper(i, box.upper(i) + d_stencil_width(i));
                stencil_boxes.push_back(box);
            }
        }
        else
        {
            for (unsigned int j = 0; j < NDIM; ++j)
            {
                for (unsigned int i = 0; i < NDIM; ++i)
                {
                    if (i == j) continue;
                    Box box(dst_box);
                    box.setLower(i, box.lower(i) - d_stencil_width(i));
                    box.setUpper(i, box.upper(i) + d_stencil_width(i));
                    box.setLower(j, box.lower(j) - d_stencil_width(j));
                    box.setUpper(j, box.upper(j) + d_stencil_width(j));
                    stencil_boxes.push_back(box);
                }
            }
        }

        // Intersect the overlap boxes with the stencil boxes.
        for (auto it1 = box_geom_overlap_boxes.begin(); it1 != box_geom_overlap_boxes.end(); ++it1)
        {
            BoxContainer overlap_boxes(stencil_boxes);
            overlap_boxes.intersectBoxes(*it1);
            for (auto it2 = overlap_boxes.begin(); it2 != overlap_boxes.end(); ++it2)
            {
                const Box& overlap_box = *it2;
                if (!overlap_box.empty()) dst_boxes.push_back(overlap_box);
            }
        }
    }
    return boost::make_shared<CellOverlap>(dst_boxes, transformation);
}

IntVector& CellNoCornersFillPattern::getStencilWidth()
{
    return d_stencil_width;
}

const std::string& CellNoCornersFillPattern::getPatternName() const
{
    return PATTERN_NAME;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
