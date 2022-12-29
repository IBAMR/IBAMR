// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2022 by the IBAMR developers
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

#include "ibtk/IBTK_MPI.h"
#include "ibtk/ParallelEdgeMap.h"

#include <tbox/Utilities.h>

#include <algorithm>
#include <map>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

int
ParallelEdgeMap::addEdge(const std::pair<int, int>& link, int mastr_idx)
{
    if (mastr_idx == invalid_index)
    {
        mastr_idx = std::min(link.first, link.second);
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(mastr_idx == link.first || mastr_idx == link.second);
#endif
    d_pending_additions.insert(std::make_pair(mastr_idx, link));
    return mastr_idx;
} // addEdge

void
ParallelEdgeMap::removeEdge(const std::pair<int, int>& link, int mastr_idx)
{
    if (mastr_idx == invalid_index)
    {
        mastr_idx = std::min(link.first, link.second);
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(mastr_idx == link.first || mastr_idx == link.second);
#endif
    d_pending_removals.insert(std::make_pair(mastr_idx, link));
    return;
} // removeEdge

void
ParallelEdgeMap::communicateData()
{
    const int size = IBTK_MPI::getNodes();
    const int rank = IBTK_MPI::getRank();

    std::vector<int> num_additions_and_removals(2 * size, 0);
    num_additions_and_removals[2 * rank] = static_cast<int>(d_pending_additions.size());
    num_additions_and_removals[2 * rank + 1] = static_cast<int>(d_pending_removals.size());
    IBTK_MPI::sumReduction(&num_additions_and_removals[0], 2 * size);

    int num_transactions = 0, offset = 0;
    for (int k = 0; k < size; ++k)
    {
        int size_k = num_additions_and_removals[2 * k] + num_additions_and_removals[2 * k + 1];
        num_transactions += size_k;
        if (k < rank)
        {
            offset += size_k;
        }
    }

    if (num_transactions == 0) return;

    static const int SIZE = 3;
    std::vector<int> transactions(SIZE * num_transactions, 0);
    for (auto cit = d_pending_additions.begin(); cit != d_pending_additions.end(); ++cit, ++offset)
    {
        transactions[SIZE * offset] = cit->first;
        transactions[SIZE * offset + 1] = cit->second.first;
        transactions[SIZE * offset + 2] = cit->second.second;
    }
    for (auto cit = d_pending_removals.begin(); cit != d_pending_removals.end(); ++cit, ++offset)
    {
        transactions[SIZE * offset] = cit->first;
        transactions[SIZE * offset + 1] = cit->second.first;
        transactions[SIZE * offset + 2] = cit->second.second;
    }
    IBTK_MPI::sumReduction(&transactions[0], SIZE * num_transactions);

    offset = 0;
    for (int k = 0; k < size; ++k)
    {
        for (int t = 0; t < num_additions_and_removals[2 * k]; ++t, ++offset)
        {
            const int mastr_idx = transactions[SIZE * offset];
            const int idx1 = transactions[SIZE * offset + 1];
            const int idx2 = transactions[SIZE * offset + 2];
            const std::pair<int, int> link = std::make_pair(idx1, idx2);
            d_pending_additions.insert(std::make_pair(mastr_idx, link));
        }
        for (int t = 0; t < num_additions_and_removals[2 * k + 1]; ++t, ++offset)
        {
            const int mastr_idx = transactions[SIZE * offset];
            const int idx1 = transactions[SIZE * offset + 1];
            const int idx2 = transactions[SIZE * offset + 2];
            const std::pair<int, int> link = std::make_pair(idx1, idx2);
            d_pending_removals.insert(std::make_pair(mastr_idx, link));
        }
    }

    using multimap_iterator = std::multimap<int, std::pair<int, int> >::iterator;
    for (const auto& pending_addition : d_pending_additions)
    {
        d_edge_map.insert(std::make_pair(pending_addition.first, pending_addition.second));
    }

    for (const auto& pending_removal : d_pending_removals)
    {
        int mastr_idx = pending_removal.first;
        const std::pair<int, int>& link = pending_removal.second;

        bool found_link = false;

        std::pair<multimap_iterator, multimap_iterator> range = d_edge_map.equal_range(mastr_idx);
        for (auto it = range.first; it != range.second && !found_link; ++it)
        {
            if (it->second == link)
            {
                found_link = true;
                d_edge_map.erase(it);
            }
        }

        if (!found_link)
        {
            const int idx1 = link.first;
            const int idx2 = link.second;
            if (mastr_idx == idx1)
            {
                mastr_idx = idx2;
            }
            else
            {
                mastr_idx = idx1;
            }

            std::pair<multimap_iterator, multimap_iterator> range = d_edge_map.equal_range(mastr_idx);
            for (auto it = range.first; it != range.second && !found_link; ++it)
            {
                if (it->second == link)
                {
                    found_link = true;
                    d_edge_map.erase(it);
                }
            }
        }
    }

    d_pending_additions.clear();
    d_pending_removals.clear();
    return;
} // communicateData

const std::multimap<int, std::pair<int, int> >&
ParallelEdgeMap::getEdgeMap() const
{
    return d_edge_map;
} // getEdgeMap

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
