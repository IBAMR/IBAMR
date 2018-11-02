// Filename: BoundingBox.h
// Created on 30 Oct 2018 by David Wells
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

#ifndef included_IBTK_ibtk_boundingbox
#define included_IBTK_ibtk_boundingbox

/////////////////////////////// INCLUDES /////////////////////////////////////
#include <ibtk/ibtk_utilities.h>

#include <CartesianPatchGeometry.h>

#include <algorithm>
#include <utility>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
// TODO: move all code to a proper source file

/*!
 * \brief Class BoundingBox implements an NDIM-dimensional bounding box
 * defined by two points.
 *
 *
*/
class BoundingBox
{
public:
    // Constructor.
    BoundingBox(const Point &bottom_point,
                const Point &top_point)
        : d_bounding_points(bottom_point, top_point)
    {
        for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
        {
            // Do not permit negative or zero volume boxes
            TBOX_ASSERT(bottom()[dim_n] < top()[dim_n]);
        }
    }

    // Constructor, starting from a SAMRAI data type.
    BoundingBox(const SAMRAI::geom::CartesianPatchGeometry<NDIM> &patch)
    {
        std::copy(patch.getXLower(), patch.getXLower() + NDIM,
                  d_bounding_points.first.data());
        std::copy(patch.getXUpper(), patch.getXUpper() + NDIM,
                  d_bounding_points.second.data());

        for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
        {
            // Do not permit negative or zero volume boxes
            TBOX_ASSERT(bottom()[dim_n] < top()[dim_n]);
        }
    }

    const Point &bottom() const
    {
        return d_bounding_points.first;
    }

    const Point &top() const
    {
        return d_bounding_points.second;
    }

    bool point_inside(const Point &point) const
    {
        for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
            if (!(bottom()[dim_n] <= point[dim_n] && point[dim_n] < top()[dim_n]))
                return false;
        return true;
    }

protected:
    std::pair<Point, Point> d_bounding_points;
};



/*!
 * \brief Class BoundingBoxes implements an NDIM-dimensional bounding box
 * defined by two points.
 *
 *
*/
class BoundingBoxes
{
public:
    BoundingBoxes(const std::vector<BoundingBox> &boxes)
        : d_boxes(boxes)
    {
        TBOX_ASSERT(0 < d_boxes.size());
        for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
        {
            d_bounding_points.first[dim_n] = std::min_element(
                d_boxes.begin(), d_boxes.end(),
                [=](const BoundingBox &a,
                    const BoundingBox &b) -> bool
                {
                    return a.bottom()[dim_n] < b.bottom()[dim_n];
                })->bottom()[dim_n];

            d_bounding_points.second[dim_n] = std::max_element(
                d_boxes.begin(), d_boxes.end(),
                [=](const BoundingBox &a,
                    const BoundingBox &b) -> bool
                {
                    return a.top()[dim_n] < b.top()[dim_n];
                })->top()[dim_n];
        }
    }

    const Point &bottom() const
    {
        return d_bounding_points.first;
    }

    const Point &top() const
    {
        return d_bounding_points.second;
    }

    bool point_inside(const Point &point) const
    {
        // start with the cheaper check:
        for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
            if (!(bottom()[dim_n] <= point[dim_n] && point[dim_n] < top()[dim_n]))
                return false;

        // TODO it would be vastly more efficient to store the collection of
        // boxes in a tree structure so that we can search in logarithmic time
        // instead of linear time.
        for (const BoundingBox &box : d_boxes)
            if (!box.point_inside(point))
                return false;
        return true;
    }

protected:
    std::pair<Point, Point> d_bounding_points;
    std::vector<BoundingBox> d_boxes;
};
}

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_ibtk_boundingbox
