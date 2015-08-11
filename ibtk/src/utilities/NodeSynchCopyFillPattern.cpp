// Filename: NodeSynchCopyFillPattern.cpp
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
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#include "SAMRAI/pdat/NodeOverlap.h"
#include "ibtk/NodeSynchCopyFillPattern.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const std::string PATTERN_NAME = "NODE_SYNCH_COPY_FILL_PATTERN";

void compute_stencil_boxes(BoxContainer& stencil_boxes, const unsigned int axis, const Box& dst_patch_box)
{
    Box stencil_box = NodeGeometry::toNodeBox(dst_patch_box);
    stencil_box.setLower(axis, stencil_box.upper(axis));
    stencil_boxes.push_back(stencil_box);
    return;
}
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

NodeSynchCopyFillPattern::NodeSynchCopyFillPattern(const unsigned int axis) : d_stencil_width(DIM, 1), d_axis(axis)
{
    // intentionally blank
    return;
}

NodeSynchCopyFillPattern::~NodeSynchCopyFillPattern()
{
    // intentionally blank
    return;
}

boost::shared_ptr<BoxOverlap> NodeSynchCopyFillPattern::calculateOverlap(const BoxGeometry& dst_geometry,
                                                                         const BoxGeometry& src_geometry,
                                                                         const Box& dst_patch_box,
                                                                         const Box& src_mask,
                                                                         const Box& fill_box,
                                                                         const bool overwrite_interior,
                                                                         const Transformation& transformation) const
{
    BoxContainer stencil_boxes;
    compute_stencil_boxes(stencil_boxes, d_axis, dst_patch_box);
    BoxContainer dst_boxes;
    auto t_dst = CPP_CAST<const NodeGeometry*>(&dst_geometry);
    auto t_src = CPP_CAST<const NodeGeometry*>(&src_geometry);
    t_dst->computeDestinationBoxes(dst_boxes, *t_src, src_mask, fill_box, overwrite_interior, transformation);
    dst_boxes.intersectBoxes(stencil_boxes);
    return boost::make_shared<NodeOverlap>(dst_boxes, transformation);
}

boost::shared_ptr<hier::BoxOverlap>
NodeSynchCopyFillPattern::computeFillBoxesOverlap(const BoxContainer& fill_boxes,
                                                  const BoxContainer& /*node_fill_boxes*/,
                                                  const Box& patch_box,
                                                  const Box& data_box,
                                                  const PatchDataFactory& /*pdf*/) const
{
    const Dimension& dim = patch_box.getDim();

    // Compute the stencil boxes.
    BoxContainer stencil_boxes;
    compute_stencil_boxes(stencil_boxes, d_axis, patch_box);

    // Convert overlap_boxes to node-based centerings.
    BoxContainer overlap_boxes(fill_boxes);
    for (auto b = overlap_boxes.begin(), b_end = overlap_boxes.end(); b != b_end; ++b)
    {
        b->growUpper(IntVector::getOne(dim));
    }
    overlap_boxes.intersectBoxes(NodeGeometry::toNodeBox(data_box));
    overlap_boxes.intersectBoxes(stencil_boxes);
    overlap_boxes.coalesce(); // to prevent redundant nodes.

    return boost::make_shared<NodeOverlap>(overlap_boxes, Transformation(IntVector::getZero(dim)));
}

IntVector& NodeSynchCopyFillPattern::getStencilWidth()
{
    return d_stencil_width;
}

const std::string& NodeSynchCopyFillPattern::getPatternName() const
{
    return PATTERN_NAME;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////
}
//////////////////////////////////////////////////////////////////////////////
