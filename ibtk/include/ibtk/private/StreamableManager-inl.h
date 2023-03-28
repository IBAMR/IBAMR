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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_StreamableManager_inl_h
#define included_IBTK_StreamableManager_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/Streamable.h"
#include "ibtk/StreamableManager.h"

#include "tbox/AbstractStream.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline size_t
StreamableManager::getDataStreamSize(const SAMRAI::tbox::Pointer<Streamable> data_item) const
{
    return SAMRAI::tbox::AbstractStream::sizeofInt() + data_item->getDataStreamSize();
} // getDataStreamSize

inline size_t
StreamableManager::getDataStreamSize(const std::vector<SAMRAI::tbox::Pointer<Streamable> >& data_items) const
{
    size_t size = SAMRAI::tbox::AbstractStream::sizeofInt();
    for (const auto& data_item : data_items)
    {
        size += getDataStreamSize(data_item);
    }
    return size;
} // getDataStreamSize

inline void
StreamableManager::packStream(SAMRAI::tbox::AbstractStream& stream, SAMRAI::tbox::Pointer<Streamable> data_item)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(data_item);
#endif
    const int streamable_id = data_item->getStreamableClassID();
    stream.pack(&streamable_id, 1);
    data_item->packStream(stream);
    return;
} // packStream

inline void
StreamableManager::packStream(SAMRAI::tbox::AbstractStream& stream,
                              std::vector<SAMRAI::tbox::Pointer<Streamable> >& data_items)
{
    const int num_data = static_cast<int>(data_items.size());
    stream.pack(&num_data, 1);
    for (auto& data_item : data_items)
    {
        packStream(stream, data_item);
    }
    return;
} // packStream

inline SAMRAI::tbox::Pointer<Streamable>
StreamableManager::unpackStream(SAMRAI::tbox::AbstractStream& stream, const SAMRAI::hier::IntVector<NDIM>& offset)
{
    int streamable_id;
    stream.unpack(&streamable_id, 1);
#if !defined(NDEBUG)
    TBOX_ASSERT(d_factory_map.count(streamable_id) == 1);
#endif
    return d_factory_map[streamable_id]->unpackStream(stream, offset);
} // unpackStream

inline void
StreamableManager::unpackStream(SAMRAI::tbox::AbstractStream& stream,
                                const SAMRAI::hier::IntVector<NDIM>& offset,
                                std::vector<SAMRAI::tbox::Pointer<Streamable> >& data_items)
{
    int num_data;
    stream.unpack(&num_data, 1);
    data_items.resize(num_data);
    for (auto& data_item : data_items)
    {
        data_item = unpackStream(stream, offset);
    }
    std::vector<SAMRAI::tbox::Pointer<Streamable> >(data_items).swap(data_items); // trim-to-fit
    return;
} // unpackStream

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_StreamableManager_inl_h
