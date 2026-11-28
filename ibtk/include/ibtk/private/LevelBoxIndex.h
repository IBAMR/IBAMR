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

#ifndef included_IBTK_LevelBoxIndex
#define included_IBTK_LevelBoxIndex

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <tbox/Array.h>

#include <Box.h>
#include <BoxArray.h>
#include <BoxList.h>
#include <PatchLevel.h>

#include <array>
#include <vector>

namespace IBTK
{
namespace detail
{
/*!
 * \brief Class LevelBoxIndex answers the queries of the BoxTree of a level in time proportional to the number of boxes
 * near the query box, instead of the number of boxes of the level.
 *
 * The tree holds the boxes of the level and their periodic images, each image with the index of its box. The class
 * indexes the same entries. Its results are identical to those of the tree, box for box and in the same order: it ranks
 * the entries in the order in which the tree lists them and sorts the entries that a query finds by that rank. The
 * entries must be disjoint. The object keeps work buffers that its queries reuse, so a query is not thread safe.
 */
class LevelBoxIndex
{
public:
    /*!
     * \brief Constructor. The object indexes the box tree of the level.
     */
    explicit LevelBoxIndex(SAMRAI::hier::PatchLevel<NDIM>& level);

    /*!
     * \brief Remove the boxes of the level and their periodic images from the boxes of fragments, as
     * BoxTree::removeIntersections() does.
     */
    void removeIntersections(SAMRAI::hier::BoxList<NDIM>& fragments);

    /*!
     * \brief Set indices to the sorted indices of the boxes of the level whose box or periodic image intersects query,
     * as BoxTree::findOverlapIndices() does.
     */
    void findOverlapIndices(SAMRAI::tbox::Array<int>& indices, const SAMRAI::hier::Box<NDIM>& query);

private:
    //! Set d_candidates to the entries that intersect query, in no particular order.
    void collect(const SAMRAI::hier::Box<NDIM>& query);

    //! The entries of the tree: the boxes of the level and their periodic images, each with the index of its box.
    SAMRAI::hier::BoxArray<NDIM> d_entries;
    SAMRAI::tbox::Array<int> d_entry_index;

    //! The position of each entry in the order in which the tree lists the entries.
    std::vector<int> d_rank;

    //! A uniform grid of buckets over the bounding box of the entries. The entries of bucket b are
    //! d_bucket_entries[d_bucket_start[b]] up to d_bucket_entries[d_bucket_start[b + 1] - 1].
    SAMRAI::hier::Box<NDIM> d_bounds;
    std::array<int, NDIM> d_bucket_size{};
    std::array<int, NDIM> d_bucket_count{};
    std::array<int, NDIM> d_bucket_stride{};
    std::vector<int> d_bucket_start, d_bucket_entries;

    //! Work buffers: the entries that the last query found, and for each entry the number of the last query that
    //! visited it.
    std::vector<int> d_candidates, d_visited;
    int d_query = 0;
    std::vector<SAMRAI::hier::Box<NDIM>> d_fragments;
};
} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_LevelBoxIndex
