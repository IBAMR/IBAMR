// Filename: EdgeSynchCopyFillPattern.cpp
// Created on 02 Feb 2011 by Boyce Griffith
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
#include "SAMRAI/pdat/EdgeGeometry.h"
#include "SAMRAI/pdat/EdgeOverlap.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "ibtk/EdgeSynchCopyFillPattern.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "boost/make_shared.hpp"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const std::string PATTERN_NAME = "EDGE_SYNCH_COPY_FILL_PATTERN";

void compute_stencil_boxes(std::vector<BoxContainer>& stencil_boxes, const unsigned int axis, const Box& dst_box)
{
    Box stencil_box = EdgeGeometry::toEdgeBox(dst_box, axis);
    stencil_box.setLower(axis, stencil_box.upper(axis));
    stencil_boxes[axis].push_back(stencil_box);
    return;
}
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

EdgeSynchCopyFillPattern::EdgeSynchCopyFillPattern(const unsigned int axis)
    : d_stencil_width(IntVector::getOne(DIM)), d_axis(axis)
{
    // intentionally blank
    return;
}

EdgeSynchCopyFillPattern::~EdgeSynchCopyFillPattern()
{
    // intentionally blank
    return;
}

boost::shared_ptr<BoxOverlap> EdgeSynchCopyFillPattern::calculateOverlap(const BoxGeometry& dst_geometry,
                                                                         const BoxGeometry& src_geometry,
                                                                         const Box& dst_patch_box,
                                                                         const Box& src_mask,
                                                                         const Box& fill_box,
                                                                         const bool overwrite_interior,
                                                                         const Transformation& transformation) const
{
    const Dimension& dim = dst_patch_box.getDim();
    std::vector<BoxContainer> stencil_boxes(dim.getValue());
    compute_stencil_boxes(stencil_boxes, d_axis, dst_patch_box);
    std::vector<hier::BoxContainer> dst_boxes(dim.getValue());
    auto t_dst = CPP_CAST<const EdgeGeometry*>(&dst_geometry);
    auto t_src = CPP_CAST<const EdgeGeometry*>(&src_geometry);
    t_dst->computeDestinationBoxes(dst_boxes, *t_src, src_mask, fill_box, overwrite_interior, transformation);
    for (int d = 0; d < dim.getValue(); ++d)
    {
        dst_boxes[d].intersectBoxes(stencil_boxes[d]);
    }
    return boost::make_shared<EdgeOverlap>(dst_boxes, transformation);
}

boost::shared_ptr<BoxOverlap> EdgeSynchCopyFillPattern::computeFillBoxesOverlap(const BoxContainer& fill_boxes,
                                                                                const BoxContainer& /*node_fill_boxes*/,
                                                                                const Box& patch_box,
                                                                                const Box& data_box,
                                                                                const PatchDataFactory& /*pdf*/) const
{
    const Dimension& dim = patch_box.getDim();

    // Compute the stencil boxes.
    std::vector<BoxContainer> stencil_boxes(dim.getValue());
    compute_stencil_boxes(stencil_boxes, d_axis, patch_box);

    // Convert overlap_boxes to edge-based centerings.
    std::vector<BoxContainer> overlap_boxes(dim.getValue());
    for (auto d = 0; d < dim.getValue(); ++d)
    {
        for (auto b = fill_boxes.begin(), b_end = fill_boxes.end(); b != b_end; ++b)
        {
            overlap_boxes[d].pushBack(EdgeGeometry::toEdgeBox(*b, d));
        }
        overlap_boxes[d].intersectBoxes(EdgeGeometry::toEdgeBox(data_box, d));
        overlap_boxes[d].intersectBoxes(stencil_boxes[d]);
        overlap_boxes[d].coalesce(); // to prevent redundant edges.
    }

    return boost::make_shared<EdgeOverlap>(overlap_boxes, Transformation(IntVector::getZero(dim)));
}

IntVector& EdgeSynchCopyFillPattern::getStencilWidth()
{
    return d_stencil_width;
}

const std::string& EdgeSynchCopyFillPattern::getPatternName() const
{
    return PATTERN_NAME;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
