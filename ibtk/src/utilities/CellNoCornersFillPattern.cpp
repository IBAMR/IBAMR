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

void compute_stencil_boxes(BoxContainer& stencil_boxes, const IntVector& stencil_width, const Box& dst_box)
{
    const Dimension& dim = dst_box.getDim();
    for (auto i = 0; i < dim.getValue(); ++i)
    {
        Box low_box(dst_box);
        low_box.setLower(i, dst_box.lower(i) - stencil_width(i));
        low_box.setUpper(i, dst_box.lower(i) - 1);
        stencil_boxes.push_back(low_box);

        Box high_box(dst_box);
        high_box.setLower(i, dst_box.upper(i) + 1);
        high_box.setUpper(i, dst_box.lower(i) + stencil_width(i));
        stencil_boxes.push_back(high_box);
    }
    return;
}
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CellNoCornersFillPattern::CellNoCornersFillPattern(const int stencil_width) : d_stencil_width(DIM, stencil_width)
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
    BoxContainer stencil_boxes;
    compute_stencil_boxes(stencil_boxes, d_stencil_width, dst_patch_box);
    return dst_geometry.calculateOverlap(src_geometry, src_mask, fill_box, overwrite_interior, transformation,
                                         stencil_boxes);
}

boost::shared_ptr<hier::BoxOverlap>
CellNoCornersFillPattern::computeFillBoxesOverlap(const BoxContainer& fill_boxes,
                                                  const BoxContainer& /*node_fill_boxes*/,
                                                  const Box& patch_box,
                                                  const Box& data_box,
                                                  const PatchDataFactory& /*pdf*/) const
{
    BoxContainer stencil_boxes;
    compute_stencil_boxes(stencil_boxes, d_stencil_width, patch_box);
    BoxContainer overlap_boxes(fill_boxes);
    overlap_boxes.intersectBoxes(data_box);
    overlap_boxes.intersectBoxes(stencil_boxes);
    return boost::make_shared<CellOverlap>(overlap_boxes, Transformation(IntVector::getZero(patch_box.getDim())));
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
