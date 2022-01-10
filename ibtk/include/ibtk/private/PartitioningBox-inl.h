// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_PartitioningBox_inl_h
#define included_IBTK_PartitioningBox_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/PartitioningBox.h"

#include <algorithm>
#include <limits>
#include <type_traits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////
inline const Point&
PartitioningBox::bottom() const
{
    return d_bounding_points.first;
} // bottom

inline const Point&
PartitioningBox::top() const
{
    return d_bounding_points.second;
} // top

inline bool
PartitioningBox::contains(const Point& point) const
{
    for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
        if (!(bottom()[dim_n] <= point[dim_n] && point[dim_n] < top()[dim_n])) return false;
    return true;
} // contains

template <typename ForwardIterator>
inline PartitioningBoxes::PartitioningBoxes(const ForwardIterator begin, const ForwardIterator end)
    : d_boxes(begin, end)
{
    static_assert(std::is_same<decltype(*begin), PartitioningBox&>::value ||
                      std::is_same<decltype(*begin), const PartitioningBox&>::value,
                  "The iterators should point to PartitioningBoxes");
    if (begin != end)
    {
        IBTK::Point bottom;
        IBTK::Point top;
        for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
        {
            bottom[dim_n] = std::min_element(begin,
                                             end,
                                             [=](const PartitioningBox& a, const PartitioningBox& b) -> bool
                                             { return a.bottom()[dim_n] < b.bottom()[dim_n]; })
                                ->bottom()[dim_n];
            top[dim_n] = std::max_element(begin,
                                          end,
                                          [=](const PartitioningBox& a, const PartitioningBox& b) -> bool
                                          { return a.top()[dim_n] < b.top()[dim_n]; })
                             ->top()[dim_n];
        }
        d_bounding_partitioning_box = PartitioningBox(bottom, top);
    }
} // PartitioningBoxes

inline const Point&
PartitioningBoxes::bottom() const
{
    return d_bounding_partitioning_box.bottom();
} // bottom

inline const Point&
PartitioningBoxes::top() const
{
    return d_bounding_partitioning_box.top();
} // top

inline const PartitioningBox*
PartitioningBoxes::begin() const
{
    return d_boxes.data();
}

inline const PartitioningBox*
PartitioningBoxes::end() const
{
    return d_boxes.data() + d_boxes.size();
}

inline bool
PartitioningBoxes::contains(const Point& point) const
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
