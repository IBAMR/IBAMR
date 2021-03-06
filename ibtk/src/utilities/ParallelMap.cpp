// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#include "ibtk/FixedSizedStream.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/ParallelMap.h"
#include "ibtk/Streamable.h"
#include "ibtk/StreamableManager.h"

#include "IntVector.h"
#include "tbox/Pointer.h"

#include <map>
#include <utility>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

ParallelMap&
ParallelMap::operator=(const ParallelMap& that)
{
    if (this != &that)
    {
        d_map = that.d_map;
        d_pending_additions = that.d_pending_additions;
        d_pending_removals = that.d_pending_removals;
    }
    return *this;
} // operator=

void
ParallelMap::addItem(const int key, const tbox::Pointer<Streamable> item)
{
    d_pending_additions.insert(std::make_pair(key, item));
    return;
} // addItem

void
ParallelMap::removeItem(const int key)
{
    d_pending_removals.push_back(key);
    return;
} // removeItem

void
ParallelMap::communicateData()
{
    const int size = IBTK_MPI::getNodes();
    const int rank = IBTK_MPI::getRank();

    // Add items to the map.
    if (IBTK_MPI::maxReduction(static_cast<int>(d_pending_additions.size())) > 0)
    {
        StreamableManager* streamable_manager = StreamableManager::getManager();

        // Determine how many keys have been registered for addition on each
        // process.
        std::vector<int> num_additions(size, 0);
        num_additions[rank] = static_cast<int>(d_pending_additions.size());
        IBTK_MPI::sumReduction(&num_additions[0], size);

        // Get the local values to send and determine the amount of data to be
        // broadcast by each process.
        std::vector<int> keys_to_send;
        std::vector<tbox::Pointer<Streamable> > data_items_to_send;
        for (const auto& pending_addition : d_pending_additions)
        {
            keys_to_send.push_back(pending_addition.first);
            data_items_to_send.push_back(pending_addition.second);
        }
        std::vector<int> data_sz(size, 0);
        data_sz[rank] = static_cast<int>(tbox::AbstractStream::sizeofInt() * keys_to_send.size() +
                                         streamable_manager->getDataStreamSize(data_items_to_send));
        IBTK_MPI::sumReduction(&data_sz[0], size);

        // Broadcast data from each process.
        for (int sending_proc = 0; sending_proc < size; ++sending_proc)
        {
            int num_keys = num_additions[sending_proc];
            if (num_keys == 0) continue;
            if (sending_proc == rank)
            {
                // Pack and broadcast data on process sending_proc.
                FixedSizedStream stream(data_sz[sending_proc]);
                stream.pack(&keys_to_send[0], static_cast<int>(keys_to_send.size()));
                streamable_manager->packStream(stream, data_items_to_send);
                int data_size = stream.getCurrentSize();
#if !defined(NDEBUG)
                TBOX_ASSERT(static_cast<int>(d_pending_additions.size()) == num_keys);
                TBOX_ASSERT(data_size == data_sz[sending_proc]);
#endif
                IBTK_MPI::bcast(static_cast<char*>(stream.getBufferStart()), data_size, sending_proc);
                for (int k = 0; k < num_keys; ++k)
                {
                    d_map[keys_to_send[k]] = data_items_to_send[k];
                }
            }
            else
            {
                // Receive and unpack data broadcast from process sending_proc.
                std::vector<char> buffer(data_sz[sending_proc]);
                int data_size = data_sz[sending_proc];
                IBTK_MPI::bcast(&buffer[0], data_size, sending_proc);
#if !defined(NDEBUG)
                TBOX_ASSERT(data_size == data_sz[sending_proc]);
#endif
                FixedSizedStream stream(&buffer[0], data_size);
                std::vector<int> keys_received(num_keys);
                stream.unpack(&keys_received[0], num_keys);
                std::vector<tbox::Pointer<Streamable> > data_items_received;
                hier::IntVector<NDIM> offset = 0;
                streamable_manager->unpackStream(stream, offset, data_items_received);
#if !defined(NDEBUG)
                TBOX_ASSERT(keys_received.size() == data_items_received.size());
#endif
                for (int k = 0; k < num_keys; ++k)
                {
                    d_map[keys_received[k]] = data_items_received[k];
                }
            }
        }

        // Clear the set of pending additions.
        d_pending_additions.clear();
    }

    // Remove items from the map.
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
                    d_map.erase(d_pending_removals[k]);
                }
            }
            else
            {
                // Receive and unpack data broadcast from process sending_proc.
                std::vector<int> keys_received(num_removals[sending_proc]);
                IBTK_MPI::bcast(&keys_received[0], num_keys, sending_proc);
                for (int k = 0; k < num_keys; ++k)
                {
                    d_map.erase(keys_received[k]);
                }
            }
        }

        // Clear the set of pending removals.
        d_pending_removals.clear();
    }
    return;
} // communicateData

const std::map<int, SAMRAI::tbox::Pointer<Streamable> >&
ParallelMap::getMap() const
{
    return d_map;
} // getMap

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
