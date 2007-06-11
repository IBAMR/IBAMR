// Filename: LNodeIndexSet.I
// Last modified: <10.Jun.2007 22:18:25 boyce@bigboy.nyconnect.com>
// Created on 29 Feb 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/StashableStream.h>

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
    StashableStream stream(data_sz, StashableStream::Write);
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
    StashableStream stream(&data[0], data_sz, StashableStream::Read);
    unpackStream(stream, d_offset);
    return;
}// getFromDatabase

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
