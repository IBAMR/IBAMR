// Filename: ParallelEdgeMap.cpp
// Created on 28 Jun 2010 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <algorithm>
#include <map>
#include <ostream>
#include <utility>
#include <vector>

#include "ibtk/ParallelEdgeMap.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

ParallelEdgeMap::ParallelEdgeMap() : d_edge_map(), d_pending_additions(), d_pending_removals()
{
    // intentionally blank
    return;
} // ParallelEdgeMap

ParallelEdgeMap::~ParallelEdgeMap()
{
    // intentionally blank
    return;
} // ~ParallelEdgeMap

int
ParallelEdgeMap::addEdge(const std::pair<int, int>& link, int mastr_idx)
{
    if (mastr_idx == -1)
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
    if (mastr_idx == -1)
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
    const int size = SAMRAI_MPI::getNodes();
    const int rank = SAMRAI_MPI::getRank();

    std::vector<int> num_additions_and_removals(2 * size, 0);
    num_additions_and_removals[2 * rank] = static_cast<int>(d_pending_additions.size());
    num_additions_and_removals[2 * rank + 1] = static_cast<int>(d_pending_removals.size());
    SAMRAI_MPI::sumReduction(&num_additions_and_removals[0], 2 * size);

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
    for (std::multimap<int, std::pair<int, int> >::const_iterator cit = d_pending_additions.begin();
         cit != d_pending_additions.end();
         ++cit, ++offset)
    {
        transactions[SIZE * offset] = cit->first;
        transactions[SIZE * offset + 1] = cit->second.first;
        transactions[SIZE * offset + 2] = cit->second.second;
    }
    for (std::multimap<int, std::pair<int, int> >::const_iterator cit = d_pending_removals.begin();
         cit != d_pending_removals.end();
         ++cit, ++offset)
    {
        transactions[SIZE * offset] = cit->first;
        transactions[SIZE * offset + 1] = cit->second.first;
        transactions[SIZE * offset + 2] = cit->second.second;
    }
    SAMRAI_MPI::sumReduction(&transactions[0], SIZE * num_transactions);

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

    typedef std::multimap<int, std::pair<int, int> >::iterator multimap_iterator;
    typedef std::multimap<int, std::pair<int, int> >::const_iterator multimap_const_iterator;
    for (multimap_const_iterator cit = d_pending_additions.begin(); cit != d_pending_additions.end(); ++cit)
    {
        d_edge_map.insert(std::make_pair(cit->first, cit->second));
    }

    typedef std::multimap<int, std::pair<int, int> >::const_iterator multimap_const_iterator;
    for (multimap_const_iterator cit = d_pending_removals.begin(); cit != d_pending_removals.end(); ++cit)
    {
        int mastr_idx = cit->first;
        const std::pair<int, int>& link = cit->second;

        bool found_link = false;

        std::pair<multimap_iterator, multimap_iterator> range = d_edge_map.equal_range(mastr_idx);
        for (multimap_iterator it = range.first; it != range.second && !found_link; ++it)
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
            for (multimap_iterator it = range.first; it != range.second && !found_link; ++it)
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
