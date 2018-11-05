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
#include <type_traits>
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
 * It is possible for the bounding box to be 'empty': in that case the bounds
 * are [oo, oo).
 */
class BoundingBox
{
public:
    // Default Constructor: an 'empty' bounding box.
    BoundingBox()
    {
        for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
        {
            d_bounding_points.first[dim_n] = std::numeric_limits<double>::infinity();
            d_bounding_points.second[dim_n] = std::numeric_limits<double>::infinity();
        }
    }

    // Constructor.
    BoundingBox(const Point &bottom_point,
                const Point &top_point)
        : d_bounding_points(bottom_point, top_point)
    {
        for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
        {
            // Do not permit negative volume boxes, but zero volume is okay
            TBOX_ASSERT(bottom()[dim_n] <= top()[dim_n]);
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
            TBOX_ASSERT(bottom()[dim_n] <= top()[dim_n]);
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

    double volume() const
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
    }

protected:
    std::pair<Point, Point> d_bounding_points;
};



/*!
 * \brief Class BoundingBoxes stores a set of bounding boxes and can check if
 * a point is in the set of bounding boxes or not in a more optimized way than
 * just looping over a std::vector<BoundingBox>.
 *
 * It is possible for the bounding box to be 'empty': in that case the bounds
 * are [oo, oo).
 */
class BoundingBoxes
{
public:
    // Default Constructor: an 'empty' bounding box.
    BoundingBoxes()
    {
        for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
        {
            d_bounding_points.first[dim_n] = std::numeric_limits<double>::infinity();
            d_bounding_points.second[dim_n] = std::numeric_limits<double>::infinity();
        }
    }

    template <typename ForwardIterator>
    BoundingBoxes(const ForwardIterator begin, const ForwardIterator end)
        : d_boxes(begin, end)
    {
        static_assert(std::is_same<decltype(*begin), BoundingBox &>::value ||
                      std::is_same<decltype(*begin), const BoundingBox &>::value,
                      "The iterators should point to BoundingBoxes");
        if (begin != end)
        {
            for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
            {
                d_bounding_points.first[dim_n] = std::min_element(
                    begin, end,
                    [=](const BoundingBox &a,
                        const BoundingBox &b) -> bool
                    {
                        return a.bottom()[dim_n] < b.bottom()[dim_n];
                    })->bottom()[dim_n];

                d_bounding_points.second[dim_n] = std::max_element(
                    begin, end,
                    [=](const BoundingBox &a,
                        const BoundingBox &b) -> bool
                    {
                        return a.top()[dim_n] < b.top()[dim_n];
                    })->top()[dim_n];
            }
        }
        else
        {
            for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
            {
                d_bounding_points.first[dim_n] = std::numeric_limits<double>::infinity();
                d_bounding_points.second[dim_n] = std::numeric_limits<double>::infinity();
            }
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
