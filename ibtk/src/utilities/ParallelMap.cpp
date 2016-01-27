// Filename: ParallelMap.cpp
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

#include <map>
#include <ostream>
#include <utility>
#include <vector>

#include "IntVector.h"
#include "ibtk/FixedSizedStream.h"
#include "ibtk/ParallelMap.h"
#include "ibtk/Streamable.h"
#include "ibtk/StreamableManager.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/AbstractStream.h"
#include "tbox/Pointer.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

ParallelMap::ParallelMap() : d_map(), d_pending_additions(), d_pending_removals()
{
    // intentionally blank
    return;
} // ParallelMap

ParallelMap::ParallelMap(const ParallelMap& from)
    : d_map(from.d_map), d_pending_additions(from.d_pending_additions), d_pending_removals(from.d_pending_removals)
{
    // intentionally blank
    return;
} // ParallelMap

ParallelMap::~ParallelMap()
{
    // intentionally blank
    return;
} // ~ParallelMap

ParallelMap& ParallelMap::operator=(const ParallelMap& that)
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
    const int size = SAMRAI_MPI::getNodes();
    const int rank = SAMRAI_MPI::getRank();

    // Add items to the map.
    if (SAMRAI_MPI::maxReduction(static_cast<int>(d_pending_additions.size())) > 0)
    {
        StreamableManager* streamable_manager = StreamableManager::getManager();

        // Determine how many keys have been registered for addition on each
        // process.
        std::vector<int> num_additions(size, 0);
        num_additions[rank] = static_cast<int>(d_pending_additions.size());
        SAMRAI_MPI::sumReduction(&num_additions[0], size);

        // Get the local values to send and determine the amount of data to be
        // broadcast by each process.
        std::vector<int> keys_to_send;
        std::vector<tbox::Pointer<Streamable> > data_items_to_send;
        for (std::map<int, tbox::Pointer<Streamable> >::const_iterator cit = d_pending_additions.begin();
             cit != d_pending_additions.end();
             ++cit)
        {
            keys_to_send.push_back(cit->first);
            data_items_to_send.push_back(cit->second);
        }
        std::vector<int> data_sz(size, 0);
        data_sz[rank] = static_cast<int>(tbox::AbstractStream::sizeofInt() * keys_to_send.size() +
                                         streamable_manager->getDataStreamSize(data_items_to_send));
        SAMRAI_MPI::sumReduction(&data_sz[0], size);

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
                SAMRAI_MPI::bcast(static_cast<char*>(stream.getBufferStart()), data_size, sending_proc);
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
                SAMRAI_MPI::bcast(&buffer[0], data_size, sending_proc);
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
                    d_map.erase(d_pending_removals[k]);
                }
            }
            else
            {
                // Receive and unpack data broadcast from process sending_proc.
                std::vector<int> keys_received(num_removals[sending_proc]);
                SAMRAI_MPI::bcast(&keys_received[0], num_keys, sending_proc);
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
