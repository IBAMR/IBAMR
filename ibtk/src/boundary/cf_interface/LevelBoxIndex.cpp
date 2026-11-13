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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/private/LevelBoxIndex.h>

#include <tbox/Utilities.h>

#include <BoxGraphUtilities.h>
#include <BoxTree.h>

#include <algorithm>
#include <limits>
#include <map>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
namespace detail
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

LevelBoxIndex::LevelBoxIndex(PatchLevel<NDIM>& level)
{
    // The entries are those of the tree: each box followed by its periodic images.
    BoxGraphUtilities<NDIM>::makeBoxesPlusPeriodicBoxes(
        d_entries, d_entry_index, level.getBoxes(), level.getShiftsForLevel());
    const int num_boxes = d_entries.getNumberOfBoxes();
    d_rank.assign(num_boxes, -1);
    d_visited.assign(num_boxes, 0);
    d_bucket_size.fill(1);
    d_bucket_count.fill(1);
    d_bucket_stride.fill(1);
    d_bucket_start.assign(2, 0);
    if (num_boxes == 0)
    {
        return;
    }

    // The bounding box, and the largest box extent along each axis.
    d_bounds = d_entries[0];
    std::array<int, NDIM> max_extent{};
    max_extent.fill(1);
    for (int i = 0; i < num_boxes; ++i)
    {
        d_bounds += d_entries[i];
        for (int d = 0; d < NDIM; ++d)
        {
            max_extent[d] = std::max(max_extent[d], d_entries[i].numberCells(d));
        }
    }

    // The entries in the order of the tree. Every query lists its entries in this order. Each entry intersects the
    // bounding box, and the entries are disjoint, so the lower corner identifies the entry.
    BoxList<NDIM> listed;
    level.getBoxTree()->findOverlapBoxes(listed, d_bounds);
    if (listed.size() != num_boxes)
    {
        TBOX_ERROR("LevelBoxIndex::LevelBoxIndex():\n"
                   << "  the tree lists " << listed.size() << " boxes but the level has " << num_boxes
                   << " boxes and periodic images\n");
    }
    std::map<std::array<int, NDIM>, int> index_of_lower;
    for (int i = 0; i < num_boxes; ++i)
    {
        std::array<int, NDIM> lower{};
        for (int d = 0; d < NDIM; ++d)
        {
            lower[d] = d_entries[i].lower(d);
        }
        if (!index_of_lower.emplace(lower, i).second)
        {
            TBOX_ERROR("LevelBoxIndex::LevelBoxIndex():\n"
                       << "  two boxes of the level or their periodic images have the lower corner of " << d_entries[i]
                       << "\n");
        }
    }
    int rank = 0;
    for (typename BoxList<NDIM>::Iterator it(listed); it; it++, ++rank)
    {
        std::array<int, NDIM> lower{};
        for (int d = 0; d < NDIM; ++d)
        {
            lower[d] = (*it).lower(d);
        }
        const auto found = index_of_lower.find(lower);
        if (found == index_of_lower.end() || !(d_entries[found->second] == *it))
        {
            TBOX_ERROR("LevelBoxIndex::LevelBoxIndex():\n"
                       << "  the tree lists the box " << *it
                       << " that is not a box of the level or a periodic image\n");
        }
        if (d_rank[found->second] != -1)
        {
            TBOX_ERROR("LevelBoxIndex::LevelBoxIndex():\n"
                       << "  the tree lists the box " << *it << " twice\n");
        }
        d_rank[found->second] = rank;
    }

    // The buckets are as large as the largest boxes, so a box overlaps at most two buckets along each axis. Doubling
    // the buckets keeps their number within a few times the number of boxes when the boxes are far apart.
    d_bucket_size = max_extent;
    const long max_buckets = std::max(64L, 4L * num_boxes);
    while (true)
    {
        long buckets = 1;
        for (int d = 0; d < NDIM; ++d)
        {
            buckets *= (d_bounds.numberCells(d) + d_bucket_size[d] - 1) / d_bucket_size[d];
        }
        if (buckets <= max_buckets)
        {
            break;
        }
        for (int d = 0; d < NDIM; ++d)
        {
            d_bucket_size[d] *= 2;
        }
    }
    for (int d = 0; d < NDIM; ++d)
    {
        d_bucket_count[d] = (d_bounds.numberCells(d) + d_bucket_size[d] - 1) / d_bucket_size[d];
        d_bucket_stride[d] = d == 0 ? 1 : d_bucket_stride[d - 1] * d_bucket_count[d - 1];
    }
    const int num_buckets = d_bucket_stride[NDIM - 1] * d_bucket_count[NDIM - 1];

    // Register each box in every bucket that it overlaps, in two passes: count, then fill.
    d_bucket_start.assign(num_buckets + 1, 0);
    for (int pass = 0; pass < 2; ++pass)
    {
        std::vector<int> next(d_bucket_start.begin(), d_bucket_start.end() - 1);
        for (int i = 0; i < num_boxes; ++i)
        {
            std::array<int, NDIM> first{};
            std::array<int, NDIM> last{};
            std::array<int, NDIM> c{};
            for (int d = 0; d < NDIM; ++d)
            {
                first[d] = (d_entries[i].lower(d) - d_bounds.lower(d)) / d_bucket_size[d];
                last[d] = (d_entries[i].upper(d) - d_bounds.lower(d)) / d_bucket_size[d];
            }
            c = first;
            while (true)
            {
                int bucket = 0;
                for (int d = 0; d < NDIM; ++d)
                {
                    bucket += c[d] * d_bucket_stride[d];
                }
                if (pass == 0)
                {
                    ++d_bucket_start[bucket + 1];
                }
                else
                {
                    d_bucket_entries[next[bucket]++] = i;
                }
                int d = 0;
                for (; d < NDIM; ++d)
                {
                    if (++c[d] <= last[d])
                    {
                        break;
                    }
                    c[d] = first[d];
                }
                if (d == NDIM)
                {
                    break;
                }
            }
        }
        if (pass == 0)
        {
            for (int b = 0; b < num_buckets; ++b)
            {
                d_bucket_start[b + 1] += d_bucket_start[b];
            }
            d_bucket_entries.assign(d_bucket_start[num_buckets], 0);
        }
    }
} // LevelBoxIndex

void
LevelBoxIndex::removeIntersections(BoxList<NDIM>& fragments)
{
    d_fragments.clear();
    for (typename BoxList<NDIM>::Iterator it(fragments); it; it++)
    {
        d_fragments.push_back(*it);
    }
    fragments.clearItems();
    for (const Box<NDIM>& fragment : d_fragments)
    {
        collect(fragment);
        if (d_candidates.empty())
        {
            fragments.appendItem(fragment);
        }
        else
        {
            // Removing a list of boxes removes them one at a time in order, which is what the tree does.
            std::sort(d_candidates.begin(),
                      d_candidates.end(),
                      [this](const int a, const int b) { return d_rank[a] < d_rank[b]; });
            BoxList<NDIM> burst(fragment);
            for (const int index : d_candidates)
            {
                burst.removeIntersections(d_entries[index]);
            }
            fragments.catenateItems(burst);
        }
    }
} // removeIntersections

void
LevelBoxIndex::findOverlapIndices(tbox::Array<int>& indices, const Box<NDIM>& query)
{
    collect(query);
    for (int& entry : d_candidates)
    {
        entry = d_entry_index[entry];
    }
    std::sort(d_candidates.begin(), d_candidates.end());
    d_candidates.erase(std::unique(d_candidates.begin(), d_candidates.end()), d_candidates.end());
    indices.resizeArray(static_cast<int>(d_candidates.size()));
    for (std::size_t j = 0; j < d_candidates.size(); ++j)
    {
        indices[static_cast<int>(j)] = d_candidates[j];
    }
} // findOverlapIndices

/////////////////////////////// PRIVATE //////////////////////////////////////

void
LevelBoxIndex::collect(const Box<NDIM>& query)
{
    d_candidates.clear();
    if (d_entries.getNumberOfBoxes() == 0 || query.empty())
    {
        return;
    }

    // The buckets that the query overlaps. The query can lie partly or fully outside the grid.
    std::array<int, NDIM> first{};
    std::array<int, NDIM> last{};
    std::array<int, NDIM> c{};
    for (int d = 0; d < NDIM; ++d)
    {
        const int lower = std::max(query.lower(d), d_bounds.lower(d));
        const int upper = std::min(query.upper(d), d_bounds.upper(d));
        if (lower > upper)
        {
            return;
        }
        first[d] = (lower - d_bounds.lower(d)) / d_bucket_size[d];
        last[d] = (upper - d_bounds.lower(d)) / d_bucket_size[d];
    }

    if (d_query == std::numeric_limits<int>::max())
    {
        std::fill(d_visited.begin(), d_visited.end(), 0);
        d_query = 0;
    }
    ++d_query;

    c = first;
    while (true)
    {
        int bucket = 0;
        for (int d = 0; d < NDIM; ++d)
        {
            bucket += c[d] * d_bucket_stride[d];
        }
        for (int j = d_bucket_start[bucket]; j < d_bucket_start[bucket + 1]; ++j)
        {
            const int index = d_bucket_entries[j];
            if (d_visited[index] != d_query)
            {
                d_visited[index] = d_query;
                if (query.intersects(d_entries[index]))
                {
                    d_candidates.push_back(index);
                }
            }
        }
        int d = 0;
        for (; d < NDIM; ++d)
        {
            if (++c[d] <= last[d])
            {
                break;
            }
            c[d] = first[d];
        }
        if (d == NDIM)
        {
            break;
        }
    }
} // collect

} // namespace detail
} // namespace IBTK
