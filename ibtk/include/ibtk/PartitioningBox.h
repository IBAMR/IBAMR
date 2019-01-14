// Filename: PartitioningBox.h
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
#ifndef included_IBTK_ibtk_partitioningbox
#define included_IBTK_ibtk_partitioningbox

/////////////////////////////// INCLUDES /////////////////////////////////////
#include <ibtk/ibtk_utilities.h>

#include <CartesianPatchGeometry.h>
#include <PatchHierarchy.h>

#include <utility>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////
namespace IBTK
{
// TODO: move all code to a proper source file

/*!
 * \brief Class PartitioningBox implements an NDIM-dimensional bounding box
 * defined by two points. Unlike a standard bounding box, a PartitioningBox is
 * an <code>NDIM</code>-dimensional tensor product of half-open intervals:
 * i.e., it is a half-open box. This property allows one to partition a domain
 * into a set of boxes.
 *
 * It is possible for a partitioning box to be 'empty': in that case the
 * bounds are [oo, oo) in each dimension.
 */
class PartitioningBox
{
public:
    /// Default Constructor: an 'empty' partitioning box.
    PartitioningBox();

    /// Constructor.
    PartitioningBox(const Point &bottom_point,
                    const Point &top_point);

    /// Constructor, starting from a SAMRAI data type.
    PartitioningBox(const SAMRAI::geom::CartesianPatchGeometry<NDIM> &patch);

    /// Get the bottom left corner of the box.
    const Point &bottom() const;

    /// Get the top right corner of the box. Note that this point is not in
    /// the box.
    const Point &top() const;

    /// Return true if the point is inside the box and false otherwise.
    bool contains(const Point &point) const;

    /// Return the volume of the box. If the box is empty then this is zero.
    double volume() const;

protected:
    /// bottom left and top right corners of the bounding box.
    std::pair<Point, Point> d_bounding_points;
};

/*!
 * \brief Class PartitioningBoxes stores a set of bounding boxes and can check
 * if a point is in the set of partitioning boxes or not in a more optimized
 * way than just looping over a std::vector<PartitioningBox>.
 *
 * It is possible for the set of bounding boxes to be empty: in that case the
 * bounds are [oo, oo) in each dimension.
 */
class PartitioningBoxes
{
public:
    /// Default constructor: an empty collection of partitioning boxes with
    /// bounds at oo.
    PartitioningBoxes() = default;

    /// Constructor. Creates an object from a sequence of partitioning boxes.
    template <typename ForwardIterator>
    PartitioningBoxes(const ForwardIterator begin, const ForwardIterator end);

    /// Constructor. Uses the finest level boxes on the provided
    /// PatchHierarchy to create a set of bounding boxes.
    PartitioningBoxes(const SAMRAI::hier::PatchHierarchy<NDIM> &hierarchy);

    /// Get the bottom left corner of the partitioning box bounding all other
    /// partitioning boxes.
    const Point &bottom() const;

    /// Get the top right corner of the partitioning box bounding all other
    /// partitioning boxes. Note that this point is not in the collection of
    /// boxes.
    const Point &top() const;

    /// Return true if the point is inside the set of boxes and false
    /// otherwise.
    bool contains(const Point &point) const;

    /// Return a pointer to the first partitioning box.
    const PartitioningBox * begin() const;

    /// Return a pointer to one past the end of the end of the partitioning
    /// box array.
    const PartitioningBox * end() const;

protected:
    /// The partitioning box that bounds all other partitioning boxes.
    PartitioningBox d_bounding_partitioning_box;

    /// The set of bounding boxes.
    std::vector<PartitioningBox> d_boxes;
};
} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/PartitioningBox-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_ibtk_partitioningbox
