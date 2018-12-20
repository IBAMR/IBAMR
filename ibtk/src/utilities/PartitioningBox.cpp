// Filename: PartitioningBox.cpp
// Created on 06 Dec 2018 by David Wells
//
// Copyright (c) 2018, Boyce Griffith, David Wells
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
#include <ibtk/ibtk_utilities.h>
#include <ibtk/PartitioningBox.h>

#include <CartesianPatchGeometry.h>

#include <algorithm>
#include <type_traits>
#include <utility>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////
namespace IBTK
{
PartitioningBox::PartitioningBox()
{
    for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
    {
        d_bounding_points.first[dim_n] = std::numeric_limits<double>::infinity();
        d_bounding_points.second[dim_n] = std::numeric_limits<double>::infinity();
    }
} // PartitioningBox

PartitioningBox::PartitioningBox(const Point &bottom_point,
                                 const Point &top_point)
    : d_bounding_points(bottom_point, top_point)
{
    for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
    {
        // Do not permit negative volume boxes
        TBOX_ASSERT(bottom()[dim_n] <= top()[dim_n]);
    }
} // PartitioningBox

PartitioningBox::PartitioningBox(const SAMRAI::geom::CartesianPatchGeometry<NDIM> &patch)
{
    std::copy(patch.getXLower(), patch.getXLower() + NDIM,
              d_bounding_points.first.data());
    std::copy(patch.getXUpper(), patch.getXUpper() + NDIM,
              d_bounding_points.second.data());
    for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
    {
        // Do not permit negative volume boxes
        TBOX_ASSERT(bottom()[dim_n] <= top()[dim_n]);
    }
} // PartitioningBox

double PartitioningBox::volume() const
{
    // The bounds may be infinite: check for point equality first to avoid
    // evaluating oo - oo
    for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
    {
        if (d_bounding_points.first[dim_n] == d_bounding_points.second[dim_n])
            return 0.0;
    }
    double vol = 1.0;
    for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
    {
        vol *= d_bounding_points.second[dim_n] - d_bounding_points.first[dim_n];
    }
    return vol;
} // volume

namespace
{
    std::vector<IBTK::PartitioningBox>
    get_partitioning_boxes(const SAMRAI::hier::PatchHierarchy<NDIM> &patch_hierarchy)
    {
        std::vector<IBTK::PartitioningBox> boxes;
        const int finest_level = patch_hierarchy.getFinestLevelNumber();
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level
            = patch_hierarchy.getPatchLevel(finest_level);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const SAMRAI::hier::Patch<NDIM>& patch = *level->getPatch(p());
            SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geometry
                = patch.getPatchGeometry();
            boxes.emplace_back(*patch_geometry);
        }
        return boxes;
    }
}

PartitioningBoxes::PartitioningBoxes(const SAMRAI::hier::PatchHierarchy<NDIM> &hierarchy)
{
    const auto boxes = get_partitioning_boxes(hierarchy);
    *this = PartitioningBoxes(boxes.begin(), boxes.end());
}

}
//////////////////////////////////////////////////////////////////////////////
