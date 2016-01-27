// Filename: StreamableManager-inl.h
// Created on 14 Jun 2004 by Boyce Griffith
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

#ifndef included_StreamableManager_inl_h
#define included_StreamableManager_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

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
    for (unsigned int k = 0; k < data_items.size(); ++k)
    {
        size += getDataStreamSize(data_items[k]);
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
    for (unsigned int k = 0; k < data_items.size(); ++k)
    {
        packStream(stream, data_items[k]);
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
    for (unsigned int k = 0; k < data_items.size(); ++k)
    {
        data_items[k] = unpackStream(stream, offset);
    }
    std::vector<SAMRAI::tbox::Pointer<Streamable> >(data_items).swap(data_items); // trim-to-fit
    return;
} // unpackStream

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_StreamableManager_inl_h
