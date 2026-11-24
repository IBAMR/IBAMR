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

#include <ibtk/IBTKInit.h>
#include <ibtk/private/BoxArrayIndex.h>

#include <tbox/Utilities.h>

#include <BoxArray.h>

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <vector>

using namespace SAMRAI;

// Compare IBTK::detail::BoxArrayIndex with a scan of all boxes and offsets for random boxes, which may be empty and may
// overlap, and for query boxes inside, outside, and straddling the boxes.

namespace
{
// A small generator with the same sequence on every platform.
class Random
{
public:
    int next(const int lower, const int upper)
    {
        d_state = d_state * 6364136223846793005ULL + 1442695040888963407ULL;
        return lower + static_cast<int>((d_state >> 33) % static_cast<std::uint64_t>(upper - lower + 1));
    }

private:
    std::uint64_t d_state = 12345;
};

hier::Box<NDIM>
random_box(Random& random, const int lower_min, const int lower_max, const int extent_max)
{
    hier::Index<NDIM> lower, upper;
    for (int d = 0; d < NDIM; ++d)
    {
        lower(d) = random.next(lower_min, lower_max);
        upper(d) = lower(d) + random.next(-1, extent_max - 1);
    }
    return hier::Box<NDIM>(lower, upper);
}
} // namespace

int
main(int argc, char** argv)
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    std::ofstream out("output");

    Random random;
    const int period = 30;
    long num_queries = 0, num_indices = 0, num_intersections = 0;
    for (int trial = 0; trial < 60; ++trial)
    {
        // Boxes in [0, period) and their images along each axis, along one axis only, or no images.
        const int num_boxes = random.next(0, 40);
        hier::BoxArray<NDIM> boxes(num_boxes);
        for (int i = 0; i < num_boxes; ++i)
        {
            boxes[i] = random_box(random, 0, period - 1, 9);
        }
        std::vector<hier::IntVector<NDIM>> offsets(1, hier::IntVector<NDIM>(0));
        const int images = trial % 3;
        for (int d = 0; d < NDIM && images > 0; ++d)
        {
            if (images == 1 && d > 0)
            {
                break;
            }
            for (int sgn = -1; sgn <= 1; sgn += 2)
            {
                hier::IntVector<NDIM> offset(0);
                offset(d) = sgn * period;
                offsets.push_back(offset);
            }
        }
        const IBTK::detail::BoxArrayIndex index(boxes, offsets);

        for (int q = 0; q < 50; ++q)
        {
            hier::Box<NDIM> query;
            switch (q % 5)
            {
            case 0: // inside
                query = random_box(random, 0, period - 1, 4);
                break;
            case 1: // outside the boxes and their images
                query = random_box(random, 3 * period, 4 * period, 6);
                break;
            case 2: // straddling the edge of the region of the boxes
                query = random_box(random, -6, 2, 12);
                break;
            case 3: // large, over many boxes and images
                query = random_box(random, -period, period, 3 * period);
                break;
            default: // possibly empty
                query = random_box(random, -period, 2 * period, 3);
                break;
            }

            std::vector<int> expected_indices;
            std::vector<hier::Box<NDIM>> expected_intersections;
            for (const hier::IntVector<NDIM>& offset : offsets)
            {
                for (int i = 0; i < num_boxes; ++i)
                {
                    const hier::Box<NDIM> image = hier::Box<NDIM>::shift(boxes[i], offset);
                    if (query.intersects(image))
                    {
                        expected_indices.push_back(i);
                        expected_intersections.push_back(query * image);
                    }
                }
            }
            std::sort(expected_indices.begin(), expected_indices.end());
            expected_indices.erase(std::unique(expected_indices.begin(), expected_indices.end()),
                                   expected_indices.end());

            std::vector<int> indices;
            index.findOverlapIndices(indices, query);
            std::vector<hier::Box<NDIM>> intersections;
            index.findIntersections(intersections, query);
            bool same = indices == expected_indices && intersections.size() == expected_intersections.size();
            for (std::size_t j = 0; same && j < intersections.size(); ++j)
            {
                same = intersections[j] == expected_intersections[j];
            }
            if (!same)
            {
                TBOX_ERROR("BoxArrayIndex differs from a scan of the boxes in trial " << trial << ", query " << q
                                                                                      << ": " << query << "\n");
            }

            ++num_queries;
            num_indices += static_cast<long>(indices.size());
            num_intersections += static_cast<long>(intersections.size());
        }
    }
    out << "queries " << num_queries << "\n"
        << "indices " << num_indices << "\n"
        << "intersections " << num_intersections << "\n";
    return 0;
}
