// Filename: StashableManager.I
// Created on 14 Jun 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)
// Last modified: <02.Oct.2006 10:37:45 boyce@boyce-griffiths-powerbook-g4-15.local>

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#endif

// C++ STDLIB INCLUDES
#include <algorithm>
#include <cassert>
#include <functional>
#include <numeric>
#include <sstream>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
size_t
StashableManager::getDataStreamSize(
    const SAMRAI::tbox::Pointer<Stashable>& stash_data) const
{
    return SAMRAI::tbox::AbstractStream::sizeofInt() + stash_data->getDataStreamSize();
}//getDataStreamSize

inline
void
StashableManager::packStream(
    SAMRAI::tbox::AbstractStream& stream,
    SAMRAI::tbox::Pointer<Stashable>& stash_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!stash_data.isNull());
#endif
    const int stashable_id = stash_data->getStashableID();
    stream.pack(&stashable_id,1);
    stash_data->packStream(stream);
    return;
}// packStream

inline
void
StashableManager::unpackStream(
    SAMRAI::tbox::AbstractStream& stream,
    const SAMRAI::hier::IntVector<NDIM>& offset,
    SAMRAI::tbox::Pointer<Stashable>& stash_data)
{
    int stashable_id;
    stream.unpack(&stashable_id,1);
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_factory_map.count(stashable_id) == 1);
#endif
    stash_data = d_factory_map[stashable_id]->unpackStream(stream,offset);
    return;
}// unpackStream

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
