// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
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

#ifndef included_IBTK_BoxArrayIndex
#define included_IBTK_BoxArrayIndex

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <Box.h>
#include <BoxArray.h>
#include <IntVector.h>
#include <PatchLevel.h>

#include <array>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
namespace detail
{
/*!
 * \brief Class BoxArrayIndex finds the boxes of a BoxArray, and of periodic or other shifted images of them, that
 * intersect a query box.
 *
 * findOverlapIndices() returns the indices in ascending order; findIntersections() orders the intersections by offset
 * and then by box. The boxes of the array may overlap and may be empty. The object stores the boxes, so the array that
 * it was built from may change afterward.
 */
class BoxArrayIndex
{
public:
    /*!
     * \brief Default constructor. The object indexes no boxes.
     */
    BoxArrayIndex() = default;

    /*!
     * \brief Constructor. The object indexes each box of boxes shifted by each of the offsets. An offset of zero
     * indexes the box itself.
     */
    BoxArrayIndex(const SAMRAI::hier::BoxArray<NDIM>& boxes, const std::vector<SAMRAI::hier::IntVector<NDIM>>& offsets);

    /*!
     * \brief Return the boxes of the array.
     */
    const SAMRAI::hier::BoxArray<NDIM>& getBoxes() const;

    /*!
     * \brief Return the offsets.
     */
    const std::vector<SAMRAI::hier::IntVector<NDIM>>& getOffsets() const;

    /*!
     * \brief Set indices to the indices, in ascending order and without repeats, of the boxes of the array such that
     * the box shifted by at least one offset intersects query.
     */
    void findOverlapIndices(std::vector<int>& indices, const SAMRAI::hier::Box<NDIM>& query) const;

    /*!
     * \brief Set intersections to the nonempty intersections of query with each box of the array shifted by each
     * offset, ordered by offset and then by index of the box.
     */
    void findIntersections(std::vector<SAMRAI::hier::Box<NDIM>>& intersections,
                           const SAMRAI::hier::Box<NDIM>& query) const;

private:
    /*!
     * \brief Set entries to the numbers, in ascending order, of the entries that intersect query.
     */
    void collect(std::vector<int>& entries, const SAMRAI::hier::Box<NDIM>& query) const;

    /*!
     * \name The boxes and offsets, and the entries: entry k * d_boxes.getNumberOfBoxes() + i is box i shifted by
     * offset k.
     */
    //\{
    SAMRAI::hier::BoxArray<NDIM> d_boxes;
    std::vector<SAMRAI::hier::IntVector<NDIM>> d_offsets;
    std::vector<SAMRAI::hier::Box<NDIM>> d_entries;
    //\}

    /*!
     * \name A uniform grid of buckets over the bounding box of the nonempty entries. The entries of bucket b are
     * d_bucket_entries[d_bucket_start[b]] up to d_bucket_entries[d_bucket_start[b + 1] - 1]. An entry is in every
     * bucket that it overlaps.
     */
    //\{
    SAMRAI::hier::Box<NDIM> d_bounds;
    std::array<int, NDIM> d_bucket_size{};
    std::array<int, NDIM> d_bucket_count{};
    std::array<int, NDIM> d_bucket_stride{};
    std::vector<int> d_bucket_start, d_bucket_entries;
    //\}
};

/*!
 * \brief Return an index of the region of the next coarser level that a finer level refines: the boxes of the finer
 * level coarsened to the next coarser level, and their periodic images shifted by plus or minus the periodic shift of
 * the coarser level along one axis at a time.
 */
BoxArrayIndex refined_region_index(const SAMRAI::hier::PatchLevel<NDIM>& finer_level,
                                   const SAMRAI::hier::IntVector<NDIM>& periodic_shift);

/*!
 * \brief Return an index of the boxes of a level, level.getBoxes(), and all of their periodic images, including images
 * along several axes. The zero offset comes first. The index of a box is the number of its patch.
 */
BoxArrayIndex periodic_level_box_index(const SAMRAI::hier::PatchLevel<NDIM>& level);
} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_BoxArrayIndex
