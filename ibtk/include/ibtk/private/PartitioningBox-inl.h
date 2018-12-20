// Filename: PartitioningBox-inl.h
// Created on 06 Dec 2018 by David Wells
//
// Copyright (c) 2018, Boyce Griffith
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

#ifndef included_IBTK_PartitioningBox_inl_h
#define included_IBTK_PartitioningBox_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/PartitioningBox.h"

#include <algorithm>
#include <limits>
#include <type_traits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////
inline const Point& PartitioningBox::bottom() const
{
    return d_bounding_points.first;
} // bottom

inline const Point& PartitioningBox::top() const
{
    return d_bounding_points.second;
} // top

inline bool PartitioningBox::contains(const Point &point) const
{
    for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
        if (!(bottom()[dim_n] <= point[dim_n] && point[dim_n] < top()[dim_n]))
            return false;
    return true;
} // contains

template <typename ForwardIterator>
inline PartitioningBoxes::PartitioningBoxes(const ForwardIterator begin, const ForwardIterator end)
    : d_boxes(begin, end)
{
    static_assert(std::is_same<decltype(*begin), PartitioningBox &>::value ||
                  std::is_same<decltype(*begin), const PartitioningBox &>::value,
                  "The iterators should point to PartitioningBoxes");
    if (begin != end)
    {
        IBTK::Point bottom;
        IBTK::Point top;
        for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
        {
            bottom[dim_n] = std::min_element(
                begin, end,
                [=](const PartitioningBox &a,
                    const PartitioningBox &b) -> bool
                {
                    return a.bottom()[dim_n] < b.bottom()[dim_n];
                })->bottom()[dim_n];
            top[dim_n] = std::max_element(
                begin, end,
                [=](const PartitioningBox &a,
                    const PartitioningBox &b) -> bool
                {
                    return a.top()[dim_n] < b.top()[dim_n];
                })->top()[dim_n];
        }
        d_bounding_partitioning_box = PartitioningBox(bottom, top);
    }
} // PartitioningBoxes

inline const Point& PartitioningBoxes::bottom() const
{
    return d_bounding_partitioning_box.bottom();
} // bottom

inline const Point& PartitioningBoxes::top() const
{
    return d_bounding_partitioning_box.top();
} // top

inline bool PartitioningBoxes::contains(const Point &point) const
{
    if (d_bounding_partitioning_box.contains(point))
    {
        // TODO it would be vastly more efficient to store the collection of
        // boxes in a tree structure so that we can search in logarithmic time
        // instead of linear time.
        for (const PartitioningBox& box : d_boxes)
        {
            if (box.contains(point))
            {
                return true;
            }
        }
    }
    return false;
} // contains

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_PartitioningBox_inl_h
