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
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideOverlap.h"
#include "boost/array.hpp"
#include "ibtk/SideNoCornersFillPattern.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "boost/make_shared.hpp"
#include "SAMRAI/tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const std::string PATTERN_NAME = "SIDE_NO_CORNERS_FILL_PATTERN";

void compute_stencil_boxes(std::vector<BoxContainer> stencil_boxes,
                           const IntVector& stencil_width,
                           const Box& dst_cell_box)
{
    const Dimension& dim = dst_cell_box.getDim();
    for (auto axis = 0; axis < dim.getValue(); ++axis)
    {
        const auto dst_box = SideGeometry::toSideBox(dst_cell_box, axis);
        for (auto i = 0; i < dim.getValue(); ++i)
        {
            Box low_box(dst_box);
            low_box.setLower(i, dst_box.lower(i) - stencil_width(i));
            low_box.setUpper(i, dst_box.lower(i) - 1);
            stencil_boxes[axis].push_back(low_box);

            Box high_box(dst_box);
            high_box.setLower(i, dst_box.upper(i) + 1);
            high_box.setUpper(i, dst_box.lower(i) + stencil_width(i));
            stencil_boxes[axis].push_back(high_box);
        }
    }
    return;
}
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

SideNoCornersFillPattern::SideNoCornersFillPattern(const int stencil_width) : d_stencil_width(DIM, stencil_width)
{
    // intentionally blank
    return;
}

SideNoCornersFillPattern::~SideNoCornersFillPattern()
{
    // intentionally blank
    return;
}

boost::shared_ptr<BoxOverlap> SideNoCornersFillPattern::calculateOverlap(const BoxGeometry& dst_geometry,
                                                                         const BoxGeometry& src_geometry,
                                                                         const Box& dst_patch_box,
                                                                         const Box& src_mask,
                                                                         const Box& fill_box,
                                                                         const bool overwrite_interior,
                                                                         const Transformation& transformation) const
{
    const Dimension& dim = dst_patch_box.getDim();
    std::vector<BoxContainer> stencil_boxes(dim.getValue());
    compute_stencil_boxes(stencil_boxes, d_stencil_width, dst_patch_box);
    std::vector<hier::BoxContainer> dst_boxes(dim.getValue());
    auto t_dst = CPP_CAST<const SideGeometry*>(&dst_geometry);
    auto t_src = CPP_CAST<const SideGeometry*>(&src_geometry);
    t_dst->computeDestinationBoxes(dst_boxes, *t_src, src_mask, fill_box, overwrite_interior, transformation);
    for (int d = 0; d < dim.getValue(); ++d)
    {
        dst_boxes[d].intersectBoxes(stencil_boxes[d]);
    }
    return boost::make_shared<SideOverlap>(dst_boxes, transformation);
}

boost::shared_ptr<hier::BoxOverlap>
SideNoCornersFillPattern::computeFillBoxesOverlap(const BoxContainer& fill_boxes,
                                                  const BoxContainer& /*node_fill_boxes*/,
                                                  const Box& patch_box,
                                                  const Box& data_box,
                                                  const PatchDataFactory& /*pdf*/) const
{
    const Dimension& dim = patch_box.getDim();

    // Compute the stencil boxes.
    std::vector<BoxContainer> stencil_boxes(dim.getValue());
    compute_stencil_boxes(stencil_boxes, d_stencil_width, patch_box);

    // Convert overlap_boxes to face-based centerings.
    std::vector<BoxContainer> overlap_boxes(dim.getValue());
    for (auto d = 0; d < dim.getValue(); ++d)
    {
        for (auto b = fill_boxes.begin(), b_end = fill_boxes.end(); b != b_end; ++b)
        {
            overlap_boxes[d].pushBack(SideGeometry::toSideBox(*b, d));
        }
        overlap_boxes[d].intersectBoxes(SideGeometry::toSideBox(data_box, d));
        overlap_boxes[d].intersectBoxes(stencil_boxes[d]);
        overlap_boxes[d].coalesce(); // to prevent redundant faces.
    }

    return boost::make_shared<SideOverlap>(overlap_boxes, Transformation(IntVector::getZero(dim)));
}

IntVector& SideNoCornersFillPattern::getStencilWidth()
{
    return d_stencil_width;
}

const std::string& SideNoCornersFillPattern::getPatternName() const
{
    return PATTERN_NAME;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
