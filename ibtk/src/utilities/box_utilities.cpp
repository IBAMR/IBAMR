// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_box_utilities
#define included_IBTK_box_utilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "tbox/ReferenceCounter.h"
#include "tbox/Utilities.h"

#include <Box.h>

#include <algorithm>
#include <iterator>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace IBTK
{
using namespace SAMRAI;
// Index::operator< does not use lexical comparison so its not useful for
// sorting: correctly implement lexical sorting here for use in, e.g.,
// std::set.
struct IndexPairLexical
{
    bool operator()(const std::pair<hier::Index<NDIM>, hier::Index<NDIM> >& a,
                    const std::pair<hier::Index<NDIM>, hier::Index<NDIM> >& b) const
    {
        int a_indices[2 * NDIM];
        int b_indices[2 * NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            a_indices[d] = a.first(d);
            a_indices[d + NDIM] = a.second(d);
            b_indices[d] = b.first(d);
            b_indices[d + NDIM] = b.second(d);
        }
        return std::lexicographical_compare(
            std::begin(a_indices), std::end(a_indices), std::begin(b_indices), std::end(b_indices));
    }
};

// Compare the number of cells between two different faces and return true
// if @p a has more cells or false if @p b has more cells.
struct IndexPairMagnitudeGreater
{
    bool operator()(const std::pair<hier::Index<NDIM>, hier::Index<NDIM> >& a,
                    const std::pair<hier::Index<NDIM>, hier::Index<NDIM> >& b) const
    {
        long a_magnitude = 1;
        long b_magnitude = 1;
        for (int d = 0; d < NDIM; ++d)
        {
            // only one entry in each product will be nonzero
            a_magnitude *= std::max<long>(1, a.second(d) - a.first(d));
            b_magnitude *= std::max<long>(1, b.second(d) - b.first(d));
        }
        return a_magnitude > b_magnitude;
    }
};

// Box::operator< does not use lexical comparison so its not useful for
// sorting: implement lexical comparison for use in, e.g., std::set.
struct BoxLexical
{
    bool operator()(const hier::Box<NDIM>& a, const hier::Box<NDIM>& b) const
    {
        int a_indices[2 * NDIM];
        int b_indices[2 * NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            a_indices[d] = a.lower(d);
            a_indices[d + NDIM] = a.upper(d);
            b_indices[d] = b.lower(d);
            b_indices[d + NDIM] = b.upper(d);
        }
        return std::lexicographical_compare(
            std::begin(a_indices), std::end(a_indices), std::begin(b_indices), std::end(b_indices));
    }
};

std::vector<SAMRAI::hier::Box<NDIM> >
merge_boxes_by_longest_edge(const std::vector<SAMRAI::hier::Box<NDIM> >& input_boxes)
{
    // Naught to do with no boxes
    if (input_boxes.size() == 0) return input_boxes;

    std::set<hier::Box<NDIM>, BoxLexical> boxes;
    // In this function we interpret the boxes as partitioning boxes
    // rather than bounding boxes: i.e., the upper bounds in each
    // coordinate direction are one past the end. This is helpful because
    // it makes these boxes with common faces have equal indices on that
    // face.
    for (const auto& box : input_boxes)
    {
        hier::Box<NDIM> new_box(box);
        const hier::Index<NDIM> ones(1);
        new_box.growUpper(ones);
        boxes.insert(new_box);
    }

    int round_n = 0;
    while (true)
    {
        std::map<std::pair<hier::Index<NDIM>, hier::Index<NDIM> >,
                 std::set<hier::Box<NDIM>, BoxLexical>,
                 IndexPairLexical>
            face_to_boxes;
        for (const hier::Box<NDIM>& box : boxes)
        {
            // 1. collect boxes by common face:
            const hier::Index<NDIM> lower = box.lower();
            const hier::Index<NDIM> upper = box.upper();
            if (NDIM == 2)
            {
                const std::pair<hier::Index<NDIM>, hier::Index<NDIM> > left = { { lower(0), lower(1) },
                                                                                { lower(0), upper(1) } };
                const std::pair<hier::Index<NDIM>, hier::Index<NDIM> > right = { { upper(0), lower(1) },
                                                                                 { upper(0), upper(1) } };
                const std::pair<hier::Index<NDIM>, hier::Index<NDIM> > bottom = { { lower(0), lower(1) },
                                                                                  { upper(0), lower(1) } };
                const std::pair<hier::Index<NDIM>, hier::Index<NDIM> > top = { { lower(0), upper(1) },
                                                                               { upper(0), upper(1) } };
                face_to_boxes[left].insert(box);
                face_to_boxes[right].insert(box);
                face_to_boxes[bottom].insert(box);
                face_to_boxes[top].insert(box);
            }
            else if (NDIM == 3)
            {
                const std::pair<hier::Index<NDIM>, hier::Index<NDIM> > left = { { lower(0), lower(1), lower(2) },
                                                                                { lower(0), upper(1), upper(2) } };
                const std::pair<hier::Index<NDIM>, hier::Index<NDIM> > right = { { upper(0), lower(1), lower(2) },
                                                                                 { upper(0), upper(1), upper(2) } };
                const std::pair<hier::Index<NDIM>, hier::Index<NDIM> > front = { { lower(0), lower(1), lower(2) },
                                                                                 { upper(0), lower(1), upper(2) } };
                const std::pair<hier::Index<NDIM>, hier::Index<NDIM> > back = { { lower(0), upper(1), lower(2) },
                                                                                { upper(0), upper(1), upper(2) } };
                const std::pair<hier::Index<NDIM>, hier::Index<NDIM> > bottom = { { lower(0), lower(1), lower(2) },
                                                                                  { upper(0), upper(1), lower(2) } };
                const std::pair<hier::Index<NDIM>, hier::Index<NDIM> > top = { { lower(0), lower(1), upper(2) },
                                                                               { upper(0), upper(1), upper(2) } };
                face_to_boxes[left].insert(box);
                face_to_boxes[right].insert(box);
                face_to_boxes[front].insert(box);
                face_to_boxes[back].insert(box);
                face_to_boxes[bottom].insert(box);
                face_to_boxes[top].insert(box);
            }
            else
            {
                TBOX_ERROR("not implemented");
            }
        }

        // 1.5. Remove any entries that correspond to just a single face
        {
            std::vector<std::pair<hier::Index<NDIM>, hier::Index<NDIM> > > faces_to_remove;
            for (const auto& pair : face_to_boxes)
                if (pair.second.size() < 2) faces_to_remove.push_back(pair.first);
            for (const auto& face : faces_to_remove) face_to_boxes.erase(face);
        }

        if (face_to_boxes.size() == 0) break;

        // 2. sort faces by size in *descending* order:
        std::vector<std::pair<hier::Index<NDIM>, hier::Index<NDIM> > > faces;
        for (const auto& pair : face_to_boxes) faces.push_back(pair.first);
        // Since some faces have equal magnitudes but don't refer to the same
        // slice of index space, use a stable sort to get platform-independent
        // iteration order in step 3:
        std::stable_sort(faces.begin(), faces.end(), IndexPairMagnitudeGreater());

        // 3. merge boxes.
        std::set<hier::Box<NDIM>, BoxLexical> removed_boxes;
        for (auto& face : faces)
        {
            const std::set<hier::Box<NDIM>, BoxLexical>& face_boxes = face_to_boxes[face];
            TBOX_ASSERT(face_boxes.size() == 2);
            bool skip_current_face = false;
            for (const hier::Box<NDIM>& box : face_boxes)
            {
                // make sure we don't merge the same box more than
                // once per round. If we have already merged one
                // of these boxes then just run the outer loop
                // again.
                if (removed_boxes.find(box) != removed_boxes.end())
                {
                    skip_current_face = true;
                    break;
                }
                else
                    removed_boxes.insert(box);
            }
            if (skip_current_face) continue;

            // remove the boxes with the current face from the list of
            // boxes on the current processor
            for (const hier::Box<NDIM>& box : face_boxes)
            {
                const auto it = boxes.find(box);
                TBOX_ASSERT(it != boxes.end());
                boxes.erase(it);
            }
            // add the merged box to the list of boxes
            const hier::Box<NDIM> box_1 = *face_boxes.begin();
            const hier::Box<NDIM> box_2 = *(++face_boxes.begin());
            const hier::Box<NDIM> union_box = box_1 + box_2;
            boxes.insert(union_box);
            const hier::Box<NDIM> intersect_box = box_1 * box_2;
            // double check that the boxes really overlap over exactly
            // one face
            TBOX_ASSERT(intersect_box.size() != 0);
            TBOX_ASSERT(box_1.size() + box_2.size() - intersect_box.size() == union_box.size());
        }
        ++round_n;
    }

    // convert back to the usual box format by undoing the box expansion
    // above.
    std::vector<hier::Box<NDIM> > result;
    const hier::Index<NDIM> negative_ones(-1);
    for (const hier::Box<NDIM>& box : boxes)
    {
        hier::Box<NDIM> new_box = box;
        new_box.growUpper(negative_ones);
        result.push_back(new_box);
    }

#ifndef NDEBUG
    // Unless we have a fancy data structure there is no linear time
    // algorithm for checking that we have exactly the same index
    // space. Do a quick sanity check that at least the volumes match:
    {
        long input_volume = 0;
        for (const hier::Box<NDIM>& box : input_boxes) input_volume += box.size();
        long output_volume = 0;
        for (const hier::Box<NDIM>& box : result) output_volume += box.size();
        TBOX_ASSERT(input_volume == output_volume);
    }
#endif

    return result;
}
} // namespace IBTK

#endif
