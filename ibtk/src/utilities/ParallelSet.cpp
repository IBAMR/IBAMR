// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2021 by the IBAMR developers
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
#include "ibtk/ParallelSet.h"

#include <tbox/Utilities.h>

#include <set>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

ParallelSet&
ParallelSet::operator=(const ParallelSet& that)
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
    const int size = IBTK_MPI::getNodes();
    const int rank = IBTK_MPI::getRank();

    // Add items from the set.
    if (IBTK_MPI::maxReduction(static_cast<int>(d_pending_additions.size())) > 0)
    {
        // Determine how many keys have been registered for addition on each
        // process.
        std::vector<int> num_additions(size, 0);
        num_additions[rank] = static_cast<int>(d_pending_additions.size());
        IBTK_MPI::sumReduction(&num_additions[0], size);

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
                IBTK_MPI::bcast(&d_pending_additions[0], num_keys, sending_proc);
                for (int k = 0; k < num_keys; ++k)
                {
                    d_set.insert(d_pending_additions[k]);
                }
            }
            else
            {
                // Receive and unpack data broadcast from process sending_proc.
                std::vector<int> keys_received(num_additions[sending_proc]);
                IBTK_MPI::bcast(&keys_received[0], num_keys, sending_proc);
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
    if (IBTK_MPI::maxReduction(static_cast<int>(d_pending_removals.size())) > 0)
    {
        // Determine how many keys have been registered for removal on each
        // process.
        std::vector<int> num_removals(size, 0);
        num_removals[rank] = static_cast<int>(d_pending_removals.size());
        IBTK_MPI::sumReduction(&num_removals[0], size);

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
                IBTK_MPI::bcast(&d_pending_removals[0], num_keys, sending_proc);
                for (int k = 0; k < num_keys; ++k)
                {
                    d_set.erase(d_pending_removals[k]);
                }
            }
            else
            {
                // Receive and unpack data broadcast from process sending_proc.
                std::vector<int> keys_received(num_removals[sending_proc]);
                IBTK_MPI::bcast(&keys_received[0], num_keys, sending_proc);
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
