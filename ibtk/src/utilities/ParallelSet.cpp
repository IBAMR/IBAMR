// Filename: ParallelSet.cpp
// Created on 03 Mar 2011 by Boyce Griffith
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

#include <ostream>
#include <set>
#include <vector>

#include "ibtk/ParallelSet.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

ParallelSet::ParallelSet() : d_set(), d_pending_additions(), d_pending_removals()
{
    // intentionally blank
    return;
} // ParallelSet

ParallelSet::ParallelSet(const ParallelSet& from)
    : d_set(from.d_set), d_pending_additions(from.d_pending_additions), d_pending_removals(from.d_pending_removals)
{
    // intentionally blank
    return;
} // ParallelSet

ParallelSet::~ParallelSet()
{
    // intentionally blank
    return;
} // ~ParallelSet

ParallelSet& ParallelSet::operator=(const ParallelSet& that)
{
    if (this != &that)
    {
        d_set = that.d_set;
        d_pending_additions = that.d_pending_additions;
        d_pending_removals = that.d_pending_removals;
    }
    return *this;
} // operator=

void
ParallelSet::addItem(const int key)
{
    d_pending_additions.push_back(key);
    return;
} // addItem

void
ParallelSet::removeItem(const int key)
{
    d_pending_removals.push_back(key);
    return;
} // removeItem

void
ParallelSet::communicateData()
{
    const int size = SAMRAI_MPI::getNodes();
    const int rank = SAMRAI_MPI::getRank();

    // Add items from the set.
    if (SAMRAI_MPI::maxReduction(static_cast<int>(d_pending_additions.size())) > 0)
    {
        // Determine how many keys have been registered for addition on each
        // process.
        std::vector<int> num_additions(size, 0);
        num_additions[rank] = static_cast<int>(d_pending_additions.size());
        SAMRAI_MPI::sumReduction(&num_additions[0], size);

        // Broadcast data from each process.
        for (int sending_proc = 0; sending_proc < size; ++sending_proc)
        {
            int num_keys = num_additions[sending_proc];
            if (num_keys == 0) continue;
            if (sending_proc == rank)
            {
// Pack and broadcast data on process sending_proc.
#if !defined(NDEBUG)
                TBOX_ASSERT(static_cast<int>(d_pending_additions.size()) == num_keys);
#endif
                SAMRAI_MPI::bcast(&d_pending_additions[0], num_keys, sending_proc);
                for (int k = 0; k < num_keys; ++k)
                {
                    d_set.insert(d_pending_additions[k]);
                }
            }
            else
            {
                // Receive and unpack data broadcast from process sending_proc.
                std::vector<int> keys_received(num_additions[sending_proc]);
                SAMRAI_MPI::bcast(&keys_received[0], num_keys, sending_proc);
                for (int k = 0; k < num_keys; ++k)
                {
                    d_set.insert(keys_received[k]);
                }
            }
        }

        // Clear the set of pending additions.
        d_pending_additions.clear();
    }

    // Remove items from the set.
    if (SAMRAI_MPI::maxReduction(static_cast<int>(d_pending_removals.size())) > 0)
    {
        // Determine how many keys have been registered for removal on each
        // process.
        std::vector<int> num_removals(size, 0);
        num_removals[rank] = static_cast<int>(d_pending_removals.size());
        SAMRAI_MPI::sumReduction(&num_removals[0], size);

        // Broadcast data from each process.
        for (int sending_proc = 0; sending_proc < size; ++sending_proc)
        {
            int num_keys = num_removals[sending_proc];
            if (num_keys == 0) continue;
            if (sending_proc == rank)
            {
// Pack and broadcast data on process sending_proc.
#if !defined(NDEBUG)
                TBOX_ASSERT(static_cast<int>(d_pending_removals.size()) == num_keys);
#endif
                SAMRAI_MPI::bcast(&d_pending_removals[0], num_keys, sending_proc);
                for (int k = 0; k < num_keys; ++k)
                {
                    d_set.erase(d_pending_removals[k]);
                }
            }
            else
            {
                // Receive and unpack data broadcast from process sending_proc.
                std::vector<int> keys_received(num_removals[sending_proc]);
                SAMRAI_MPI::bcast(&keys_received[0], num_keys, sending_proc);
                for (int k = 0; k < num_keys; ++k)
                {
                    d_set.erase(keys_received[k]);
                }
            }
        }

        // Clear the set of pending removals.
        d_pending_removals.clear();
    }
    return;
} // communicateData

const std::set<int>&
ParallelSet::getSet() const
{
    return d_set;
} // getSet

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
