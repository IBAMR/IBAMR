// Filename: StashableManager.C
// Created on 14 Jun 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)
// Last modified: <16.Apr.2007 05:33:25 boyce@bigboy.nyconnect.com>

#include "StashableManager.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// STOOLS INCLUDES
#include <stools/STOOLS_Utilities.h>

// SAMRAI INCLUDES
#include <tbox/MPI.h>
#include <tbox/ShutdownRegistry.h>

// C++ STDLIB INCLUDES
#include <algorithm>
#include <functional>
#include <numeric>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
struct StashableGetDataStreamSizeSum
    : std::binary_function<size_t,SAMRAI::tbox::Pointer<Stashable>,size_t>
{
    inline size_t
    operator()(
        size_t size_so_far,
        const SAMRAI::tbox::Pointer<Stashable>& data) const
        {
            return size_so_far+StashableManager::getManager()->getDataStreamSize(data);
        }
};

class StashablePackStream
    : public std::unary_function<SAMRAI::tbox::Pointer<Stashable>,void>
{
public:
    inline
    StashablePackStream(
        SAMRAI::tbox::AbstractStream* const stream)
        : d_stream(stream)
        {
            return;
        }

    inline void
    operator()(
        SAMRAI::tbox::Pointer<Stashable>& data) const
        {
            StashableManager::getManager()->packStream(*d_stream,data);
            return;
        }

private:
    SAMRAI::tbox::AbstractStream* const d_stream;
};

class StashableUnpackStream
    : public std::unary_function<void,SAMRAI::tbox::Pointer<Stashable> >
{
public:
    inline
    StashableUnpackStream(
        SAMRAI::tbox::AbstractStream* const stream,
        const SAMRAI::hier::IntVector<NDIM>& offset)
        : d_stream(stream),
          d_offset(offset)
        {
            return;
        }

    inline SAMRAI::tbox::Pointer<Stashable>
    operator()() const
        {
            SAMRAI::tbox::Pointer<Stashable> data_out;
            StashableManager::getManager()->unpackStream(
                *d_stream,d_offset,data_out);
            return data_out;
        }

private:
    SAMRAI::tbox::AbstractStream* const d_stream;
    const SAMRAI::hier::IntVector<NDIM>& d_offset;
};
}

StashableManager* StashableManager::s_data_manager_instance = NULL;
bool StashableManager::s_registered_callback = false;
int StashableManager::s_current_id_number = 0;
const int StashableManager::s_unregistered_number = -1;
unsigned char StashableManager::s_shutdown_priority = 200;

StashableManager*
StashableManager::getManager()
{
    if (s_data_manager_instance == NULL)
    {
        s_data_manager_instance = new StashableManager();
    }
    if (!s_registered_callback)
    {
        SAMRAI::tbox::ShutdownRegistry::registerShutdownRoutine(
            freeManager,s_shutdown_priority);
        s_registered_callback = true;
    }
    return s_data_manager_instance;
}// getManager

void
StashableManager::freeManager()
{
    delete s_data_manager_instance;
    s_data_manager_instance = NULL;
    return;
}// freeManager

int
StashableManager::getUnregisteredID()
{
    return s_unregistered_number;
}// getUnregisteredID

/////////////////////////////// PUBLIC ///////////////////////////////////////

bool
StashableManager::checkFactoryRegistration(
    SAMRAI::tbox::Pointer<StashableFactory> factory)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!factory.isNull());
#endif
    return d_factory_map.count(factory->getStashableID()) == 1;
}// checkFactoryRegistration

int
StashableManager::registerFactory(
    SAMRAI::tbox::Pointer<StashableFactory> factory)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!factory.isNull());
    assert(factory->getStashableID() == getUnregisteredID());
#endif
    SAMRAI::tbox::MPI::barrier();  // ensure all processes use the same stashable ID
                                   // for the registered class
    const int factory_id = getUniqueID();
    factory->setStashableID(factory_id);
    STOOLS::efficient_add_or_update(d_factory_map, factory_id, factory);
    return factory_id;
}// registerFactory

size_t
StashableManager::getDataStreamSize(
    const std::vector<SAMRAI::tbox::Pointer<Stashable> >& stash_data) const
{
    return std::accumulate(stash_data.begin(), stash_data.end(),
                           SAMRAI::tbox::AbstractStream::sizeofInt(),
                           StashableGetDataStreamSizeSum());
}//getDataStreamSize

void
StashableManager::packStream(
    SAMRAI::tbox::AbstractStream& stream,
    std::vector<SAMRAI::tbox::Pointer<Stashable> >& stash_data)
{
    const int num_data = stash_data.size();
    stream.pack(&num_data,1);
    for_each(stash_data.begin(),stash_data.end(),
             StashablePackStream(&stream));
    return;
}// packStream

void
StashableManager::unpackStream(
    SAMRAI::tbox::AbstractStream& stream,
    const SAMRAI::hier::IntVector<NDIM>& offset,
    std::vector<SAMRAI::tbox::Pointer<Stashable> >& stash_data)
{
    int num_data;
    stream.unpack(&num_data,1);
    stash_data.resize(num_data);
    generate(stash_data.begin(),stash_data.end(),
             StashableUnpackStream(&stream,offset));
    std::vector<SAMRAI::tbox::Pointer<Stashable> >(stash_data).swap(stash_data); // trim-to-fit
    return;
}// unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

StashableManager::StashableManager()
    : d_factory_map()
{
    // intentionally blank
    return;
}// StashableManager

StashableManager::~StashableManager()
{
    d_factory_map.clear();
    return;
}// ~StashableManager

int
StashableManager::getUniqueID()
{
    return s_current_id_number++;
}// getUniqueID

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
