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

#include <ibtk/private/BoxArrayIndex.h>

#include <GridGeometry.h>

#include <algorithm>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
namespace detail
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Call action(bucket, c) for each bucket with grid coordinates c in the range from first to last along each axis. The
// buckets are numbered with the first axis varying fastest.
template <class Action>
void
for_each_bucket(const std::array<int, NDIM>& first,
                const std::array<int, NDIM>& last,
                const std::array<int, NDIM>& stride,
                Action action)
{
    std::array<int, NDIM> c = first;
    while (true)
    {
        int bucket = 0;
        for (int d = 0; d < NDIM; ++d)
        {
            bucket += c[d] * stride[d];
        }
        action(bucket, c);
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
            return;
        }
    }
}
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

BoxArrayIndex::BoxArrayIndex(const BoxArray<NDIM>& boxes, const std::vector<IntVector<NDIM>>& offsets)
    : d_boxes(boxes), d_offsets(offsets)
{
    const int num_boxes = d_boxes.getNumberOfBoxes();
    d_entries.reserve(d_offsets.size() * num_boxes);
    for (const IntVector<NDIM>& offset : d_offsets)
    {
        for (int i = 0; i < num_boxes; ++i)
        {
            d_entries.push_back(Box<NDIM>::shift(d_boxes[i], offset));
        }
    }

    // The bounding box of the nonempty entries and their largest extent along each axis. An empty box has no cells and
    // does not change the bounding box.
    std::array<int, NDIM> max_extent{};
    max_extent.fill(1);
    for (const Box<NDIM>& entry : d_entries)
    {
        d_bounds += entry;
        for (int d = 0; d < NDIM; ++d)
        {
            max_extent[d] = std::max(max_extent[d], entry.numberCells(d));
        }
    }

    // The buckets are as large as the largest entries, so an entry overlaps at most two buckets along each axis.
    // Doubling the buckets keeps their number within a few times the number of entries when the entries are far apart.
    d_bucket_size = max_extent;
    const long max_buckets = std::max(64L, 4L * static_cast<long>(d_entries.size()));
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

    // Register each nonempty entry in every bucket that it overlaps, in two passes: count, then fill.
    d_bucket_start.assign(num_buckets + 1, 0);
    for (int pass = 0; pass < 2; ++pass)
    {
        std::vector<int> next(d_bucket_start.begin(), d_bucket_start.end() - 1);
        for (std::size_t e = 0; e < d_entries.size(); ++e)
        {
            const Box<NDIM>& entry = d_entries[e];
            if (entry.empty())
            {
                continue;
            }
            std::array<int, NDIM> first{};
            std::array<int, NDIM> last{};
            for (int d = 0; d < NDIM; ++d)
            {
                first[d] = (entry.lower(d) - d_bounds.lower(d)) / d_bucket_size[d];
                last[d] = (entry.upper(d) - d_bounds.lower(d)) / d_bucket_size[d];
            }
            for_each_bucket(first,
                            last,
                            d_bucket_stride,
                            [&](const int bucket, const std::array<int, NDIM>&)
                            {
                                if (pass == 0)
                                {
                                    ++d_bucket_start[bucket + 1];
                                }
                                else
                                {
                                    d_bucket_entries[next[bucket]++] = static_cast<int>(e);
                                }
                            });
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
} // BoxArrayIndex

const BoxArray<NDIM>&
BoxArrayIndex::getBoxes() const
{
    return d_boxes;
} // getBoxes

const std::vector<IntVector<NDIM>>&
BoxArrayIndex::getOffsets() const
{
    return d_offsets;
} // getOffsets

void
BoxArrayIndex::findOverlapIndices(std::vector<int>& indices, const Box<NDIM>& query) const
{
    collect(indices, query);
    const int num_boxes = d_boxes.getNumberOfBoxes();
    for (int& entry : indices)
    {
        entry %= num_boxes;
    }
    std::sort(indices.begin(), indices.end());
    indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
} // findOverlapIndices

void
BoxArrayIndex::findIntersections(std::vector<Box<NDIM>>& intersections, const Box<NDIM>& query) const
{
    std::vector<int> entries;
    collect(entries, query);
    intersections.clear();
    for (const int entry : entries)
    {
        intersections.push_back(query * d_entries[entry]);
    }
} // findIntersections

BoxArrayIndex
refined_region_index(const PatchLevel<NDIM>& finer_level, const IntVector<NDIM>& periodic_shift)
{
    BoxArray<NDIM> boxes = finer_level.getBoxes();
    boxes.coarsen(finer_level.getRatioToCoarserLevel());
    std::vector<IntVector<NDIM>> offsets(1, IntVector<NDIM>(0));
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        if (periodic_shift(axis) != 0)
        {
            for (int sgn = -1; sgn <= 1; sgn += 2)
            {
                IntVector<NDIM> offset(0);
                offset(axis) = sgn * periodic_shift(axis);
                offsets.push_back(offset);
            }
        }
    }
    return BoxArrayIndex(boxes, offsets);
} // refined_region_index

BoxArrayIndex
periodic_level_box_index(const PatchLevel<NDIM>& level)
{
    const IntVector<NDIM> periodic_shift = level.getGridGeometry()->getPeriodicShift(level.getRatio());
    std::vector<IntVector<NDIM>> offsets(1, IntVector<NDIM>(0));
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (periodic_shift(d) == 0)
        {
            continue;
        }
        const auto n = offsets.size();
        for (int image = -1; image <= 1; image += 2)
        {
            for (std::size_t k = 0; k < n; ++k)
            {
                IntVector<NDIM> offset = offsets[k];
                offset(d) = image * periodic_shift(d);
                offsets.push_back(offset);
            }
        }
    }
    return BoxArrayIndex(level.getBoxes(), offsets);
} // periodic_level_box_index

/////////////////////////////// PRIVATE //////////////////////////////////////

void
BoxArrayIndex::collect(std::vector<int>& entries, const Box<NDIM>& query) const
{
    entries.clear();
    if (d_bounds.empty() || query.empty())
    {
        return;
    }

    // The buckets that the query overlaps. The query can lie partly or fully outside the grid.
    std::array<int, NDIM> first{};
    std::array<int, NDIM> last{};
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

    // An entry that overlaps several of these buckets is reported once, from the bucket that holds the lower corner of
    // its intersection with the query.
    for_each_bucket(first,
                    last,
                    d_bucket_stride,
                    [&](const int bucket, const std::array<int, NDIM>& c)
                    {
                        for (int j = d_bucket_start[bucket]; j < d_bucket_start[bucket + 1]; ++j)
                        {
                            const int entry = d_bucket_entries[j];
                            const Box<NDIM>& box = d_entries[entry];
                            if (!query.intersects(box))
                            {
                                continue;
                            }
                            bool owner = true;
                            for (int d = 0; d < NDIM && owner; ++d)
                            {
                                const int corner = std::max(query.lower(d), box.lower(d));
                                owner = (corner - d_bounds.lower(d)) / d_bucket_size[d] == c[d];
                            }
                            if (owner)
                            {
                                entries.push_back(entry);
                            }
                        }
                    });
    std::sort(entries.begin(), entries.end());
} // collect

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
