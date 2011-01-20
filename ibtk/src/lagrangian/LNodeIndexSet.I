// Filename: LNodeIndexSet.I
// Created on 29 Feb 2004 by Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/FixedSizedStream.h>

// C++ STDLIB INCLUDES
#include <algorithm>
#include <functional>
#include <numeric>
#include <sstream>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
LNodeIndexSet::LNodeIndexSet()
    : d_set(),
      d_offset(0)
{
    // intentionally blank
    return;
}// LNodeIndexSet

inline
LNodeIndexSet::LNodeIndexSet(
    const LNodeIndexSet& from)
    : d_set(from.d_set),
      d_offset(from.d_offset)
{
    // intentionally blank
    return;
}// LNodeIndexSet

inline
LNodeIndexSet::~LNodeIndexSet()
{
    d_set.clear();
    return;
}// ~LNodeIndexSet

inline LNodeIndexSet&
LNodeIndexSet::operator=(
    const LNodeIndexSet& that)
{
    if (this == &that) return *this;  // check for self-assignment

    d_set    = that.d_set;
    d_offset = that.d_offset;
    return *this;
}// operator=

inline LNodeIndexSet::const_iterator
LNodeIndexSet::begin() const
{
    return d_set.begin();
}// begin

inline LNodeIndexSet::iterator
LNodeIndexSet::begin()
{
    return d_set.begin();
}// begin

inline LNodeIndexSet::const_iterator
LNodeIndexSet::end() const
{
    return d_set.end();
}// end

inline LNodeIndexSet::iterator
LNodeIndexSet::end()
{
    return d_set.end();
}// end

inline LNodeIndexSet::size_type
LNodeIndexSet::size() const
{
    return d_set.size();
}// size

inline bool
LNodeIndexSet::empty() const
{
    return d_set.empty();
}// empty

inline void
LNodeIndexSet::push_back(
    const LNodeIndexSet::value_type& value)
{
    d_set.push_back(value);
    return;
}// push_back

inline const SAMRAI::hier::IntVector<NDIM>&
LNodeIndexSet::getPeriodicOffset() const
{
    return d_offset;
}// getPeriodicOffset

inline void
LNodeIndexSet::setPeriodicOffset(
    const SAMRAI::hier::IntVector<NDIM>& offset)
{
    d_offset = offset;
    return;
}// setPeriodicOffset

inline void
LNodeIndexSet::copySourceItem(
    const SAMRAI::hier::Index<NDIM>& src_index,
    const SAMRAI::hier::IntVector<NDIM>& src_offset,
    const LNodeIndexSet& src_item)
{
    d_set    = src_item.d_set;
    d_offset = src_offset;
    return;
}// copySourceItem

inline size_t
LNodeIndexSet::getDataStreamSize() const
{
    return std::accumulate(d_set.begin(),d_set.end(),
                           SAMRAI::tbox::AbstractStream::sizeofInt(),
                           LNodeIndexGetDataStreamSizeSum());
}// getDataStreamSize

inline void
LNodeIndexSet::packStream(
    SAMRAI::tbox::AbstractStream& stream) const
{
    const int num_idx = d_set.size();
    stream.pack(&num_idx,1);
    for_each(d_set.begin(),d_set.end(), LNodeIndexPackStream(&stream));
    return;
}// packStream

inline void
LNodeIndexSet::unpackStream(
    SAMRAI::tbox::AbstractStream& stream,
    const SAMRAI::hier::IntVector<NDIM>& offset)
{
    int num_idx;
    stream.unpack(&num_idx,1);
    d_set.resize(num_idx);
    generate(d_set.begin(),d_set.end(),LNodeIndexUnpackStream(&stream,offset));

    d_offset = offset;
    return;
}// unpackStream

inline void
LNodeIndexSet::putToDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& database) const
{
    const size_t data_sz = getDataStreamSize();
    FixedSizedStream stream(data_sz);
    packStream(stream);
    database->putInteger("data_sz", data_sz);
    database->putCharArray("data", static_cast<char*>(stream.getBufferStart()), data_sz);
    database->putIntegerArray("d_offset", d_offset, NDIM);
    return;
}// putToDatabase

inline void
LNodeIndexSet::getFromDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& database)
{
    database->getIntegerArray("d_offset", d_offset, NDIM);

    const size_t data_sz = database->getInteger("data_sz");
    std::vector<char> data(data_sz);
    database->getCharArray("data", &data[0], data_sz);
    FixedSizedStream stream(&data[0], data_sz);
    unpackStream(stream, d_offset);
    return;
}// getFromDatabase

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
