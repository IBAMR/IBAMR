// Filename: StreamableManager.C
// Created on 14 Jun 2004 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#include "StreamableManager.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/efficient_add_or_update.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <tbox/SAMRAI_MPI.h>
#include <tbox/ShutdownRegistry.h>

// C++ STDLIB INCLUDES
#include <algorithm>
#include <functional>
#include <numeric>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
struct StreamableGetDataStreamSizeSum
    : std::binary_function<size_t,Pointer<Streamable>,size_t>
{
    inline size_t
    operator()(
        size_t size_so_far,
        const Pointer<Streamable>& data) const
        {
            return size_so_far+StreamableManager::getManager()->getDataStreamSize(data);
        }
};

class StreamablePackStream
    : public std::unary_function<Pointer<Streamable>,void>
{
public:
    inline
    StreamablePackStream(
        AbstractStream* const stream)
        : d_stream(stream)
        {
            return;
        }

    inline void
    operator()(
        Pointer<Streamable>& data) const
        {
            StreamableManager::getManager()->packStream(*d_stream,data);
            return;
        }

private:
    AbstractStream* const d_stream;
};

class StreamableUnpackStream
    : public std::unary_function<void,Pointer<Streamable> >
{
public:
    inline
    StreamableUnpackStream(
        AbstractStream* const stream,
        const IntVector<NDIM>& offset)
        : d_stream(stream),
          d_offset(offset)
        {
            return;
        }

    inline Pointer<Streamable>
    operator()() const
        {
            Pointer<Streamable> data_out;
            StreamableManager::getManager()->unpackStream(
                *d_stream,d_offset,data_out);
            return data_out;
        }

private:
    AbstractStream* const d_stream;
    const IntVector<NDIM>& d_offset;
};
}

StreamableManager* StreamableManager::s_data_manager_instance = NULL;
bool StreamableManager::s_registered_callback = false;
int StreamableManager::s_current_id_number = 0;
const int StreamableManager::s_unregistered_id_number = -1;
unsigned char StreamableManager::s_shutdown_priority = 200;

StreamableManager*
StreamableManager::getManager()
{
    if (s_data_manager_instance == NULL)
    {
        s_data_manager_instance = new StreamableManager();
    }
    if (!s_registered_callback)
    {
        ShutdownRegistry::registerShutdownRoutine(
            freeManager,s_shutdown_priority);
        s_registered_callback = true;
    }
    return s_data_manager_instance;
}// getManager

void
StreamableManager::freeManager()
{
    delete s_data_manager_instance;
    s_data_manager_instance = NULL;
    return;
}// freeManager

int
StreamableManager::getUnregisteredID()
{
    return s_unregistered_id_number;
}// getUnregisteredID

/////////////////////////////// PUBLIC ///////////////////////////////////////

bool
StreamableManager::checkFactoryRegistration(
    Pointer<StreamableFactory> factory)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!factory.isNull());
#endif
    return d_factory_map.count(factory->getStreamableClassID()) == 1;
}// checkFactoryRegistration

int
StreamableManager::registerFactory(
    Pointer<StreamableFactory> factory)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!factory.isNull());
    TBOX_ASSERT(factory->getStreamableClassID() == getUnregisteredID());
#endif
    // These barriers ensure that each factory is assigned the same class ID
    // number on each MPI process.
    SAMRAI_MPI::barrier();
    const int factory_id = createUniqueID();
    SAMRAI_MPI::barrier();
    factory->setStreamableClassID(factory_id);
    efficient_add_or_update(d_factory_map, factory_id, factory);
    return factory_id;
}// registerFactory

size_t
StreamableManager::getDataStreamSize(
    const std::vector<Pointer<Streamable> >& data_items) const
{
    return std::accumulate(data_items.begin(), data_items.end(),
                           AbstractStream::sizeofInt(),
                           StreamableGetDataStreamSizeSum());
}//getDataStreamSize

void
StreamableManager::packStream(
    AbstractStream& stream,
    std::vector<Pointer<Streamable> >& data_items)
{
    const int num_data = data_items.size();;
    stream.pack(&num_data,1);
    for_each(&data_items[0],&data_items[0]+data_items.size(),StreamablePackStream(&stream));
    return;
}// packStream

void
StreamableManager::unpackStream(
    AbstractStream& stream,
    const IntVector<NDIM>& offset,
    std::vector<Pointer<Streamable> >& data_items)
{
    int num_data;
    stream.unpack(&num_data,1);
    data_items.resize(num_data);
    generate(data_items.begin(),data_items.end(),StreamableUnpackStream(&stream,offset));
    std::vector<Pointer<Streamable> >(data_items).swap(data_items); // trim-to-fit
    return;
}// unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

StreamableManager::StreamableManager()
    : d_factory_map()
{
    // intentionally blank
    return;
}// StreamableManager

StreamableManager::~StreamableManager()
{
    d_factory_map.clear();
    return;
}// ~StreamableManager

int
StreamableManager::createUniqueID()
{
    return s_current_id_number++;
}// createUniqueID

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
