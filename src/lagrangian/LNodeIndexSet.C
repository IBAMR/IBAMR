// Filename: LNodeIndexSet.C
// Created on 29 Feb 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)
// Last modified: <25.Oct.2006 18:27:52 boyce@bigboy.nyconnect.com>

#include "LNodeIndexSet.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/StashableStream.h>

// C++ STDLIB INCLUDES
#include <algorithm>
#include <functional>
#include <numeric>
#include <sstream>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

namespace
{
struct LNodeIndexGetDataStreamSizeSum
    : binary_function<size_t,SAMRAI::tbox::Pointer<LNodeIndex>,size_t>
{
    size_t operator()(
        size_t size_so_far,
        const SAMRAI::tbox::Pointer<LNodeIndex>& index) const
        {
            return size_so_far+index->getDataStreamSize();
        }
};
}

size_t
LNodeIndexSet::getDataStreamSize() const
{
    return accumulate(d_set.begin(),d_set.end(),
                      SAMRAI::tbox::AbstractStream::sizeofInt(),
                      LNodeIndexGetDataStreamSizeSum());
}// getDataStreamSize

namespace
{
class LNodeIndexPackStream
    : public unary_function<SAMRAI::tbox::Pointer<LNodeIndex>,void>
{
public:
    LNodeIndexPackStream(
        SAMRAI::tbox::AbstractStream* const stream)
        : d_stream(stream)
        {
            return;
        }

    void operator()(
        const SAMRAI::tbox::Pointer<LNodeIndex>& index) const
        {
            index->packStream(*d_stream);
            return;
        }
private:
    SAMRAI::tbox::AbstractStream* const d_stream;
};
}

void
LNodeIndexSet::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
    const int num_idx = d_set.size();
    stream.pack(&num_idx,1);
    for_each(d_set.begin(),d_set.end(), LNodeIndexPackStream(&stream));
    return;
}// packStream

namespace
{
class LNodeIndexUnpackStream
    : public unary_function<void,SAMRAI::tbox::Pointer<LNodeIndex> >
{
public:
    LNodeIndexUnpackStream(
        SAMRAI::tbox::AbstractStream* const stream,
        const SAMRAI::hier::IntVector<NDIM>& offset)
        : d_stream(stream),
          d_offset(offset)
        {
            return;
        }

    SAMRAI::tbox::Pointer<LNodeIndex> operator()() const
        {
            SAMRAI::tbox::Pointer<LNodeIndex> index_out = new LNodeIndex();
            index_out->unpackStream(*d_stream,d_offset);
            return index_out;
        }
private:
    SAMRAI::tbox::AbstractStream* const d_stream;
    const SAMRAI::hier::IntVector<NDIM>& d_offset;
};
}

void
LNodeIndexSet::unpackStream(
    SAMRAI::tbox::AbstractStream& stream,
    const SAMRAI::hier::IntVector<NDIM>& offset)
{
    int num_idx;
    stream.unpack(&num_idx,1);
    d_set.resize(num_idx);
    generate(d_set.begin(),d_set.end(),
             LNodeIndexUnpackStream(&stream,offset));
    trimToFit();

    d_offset = offset;
    return;
}// unpackStream

void
LNodeIndexSet::putToDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& database)
{
    const size_t data_sz = getDataStreamSize();
    static const bool use_xdr = false;
    StashableStream stream(data_sz, StashableStream::Write, use_xdr);
    packStream(stream);
    database->putInteger("data_sz", data_sz);
    database->putCharArray("data", static_cast<char*>(stream.getBufferStart()), data_sz);
    database->putIntegerArray("d_offset", d_offset, NDIM);
    return;
}// putToDatabase

void
LNodeIndexSet::getFromDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& database)
{
    database->getIntegerArray("d_offset", d_offset, NDIM);

    const size_t data_sz = database->getInteger("data_sz");
    vector<char> data(data_sz);
    database->getCharArray("data", &data[0], data_sz);
    static const bool use_xdr = false;
    StashableStream stream(&data[0], data_sz, StashableStream::Read, use_xdr);
    unpackStream(stream, d_offset);
    return;
}// getFromDatabase

/////////////////////////////// PRIVATE //////////////////////////////////////

namespace
{
struct LNodeIndexLessThan
    : binary_function<SAMRAI::tbox::Pointer<LNodeIndex>,SAMRAI::tbox::Pointer<LNodeIndex>,bool>
{
    bool operator()(
        const SAMRAI::tbox::Pointer<LNodeIndex>& lhs,
        const SAMRAI::tbox::Pointer<LNodeIndex>& rhs) const
        {
            return *lhs < *rhs;
        }
};
}

void
LNodeIndexSet::reorderCollection()
{
    sort(d_set.begin(),d_set.end(),LNodeIndexLessThan());
    return;
}// sortSet

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <IndexData.C>
#include <IndexDataFactory.C>
#include <IndexVariable.C>
#include <tbox/Array.C>
#include <tbox/List.C>
#include <tbox/Pointer.C>

template class SAMRAI::tbox::Pointer<IBAMR::LNodeIndexSet>;

#if (NDIM == 1)
template class SAMRAI::pdat::IndexData<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexDataFactory<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexDataNode<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexIterator<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexVariable<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::tbox::Array<IBAMR::LNodeIndexSet>;
template class SAMRAI::tbox::Array<SAMRAI::pdat::IndexDataNode<1,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::List<SAMRAI::pdat::IndexDataNode<1,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::ListIterator<SAMRAI::pdat::IndexDataNode<1,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::ListNode<SAMRAI::pdat::IndexDataNode<1,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<1,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexVariable<1,IBAMR::LNodeIndexSet> >;
#endif

#if (NDIM == 2)
template class SAMRAI::pdat::IndexData<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexDataFactory<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexDataNode<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexIterator<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexVariable<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::tbox::Array<IBAMR::LNodeIndexSet>;
template class SAMRAI::tbox::Array<SAMRAI::pdat::IndexDataNode<2,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::List<SAMRAI::pdat::IndexDataNode<2,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::ListIterator<SAMRAI::pdat::IndexDataNode<2,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::ListNode<SAMRAI::pdat::IndexDataNode<2,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<2,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexVariable<2,IBAMR::LNodeIndexSet> >;
#endif

#if (NDIM == 3)
template class SAMRAI::pdat::IndexData<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexDataFactory<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexDataNode<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexIterator<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexVariable<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::tbox::Array<IBAMR::LNodeIndexSet>;
template class SAMRAI::tbox::Array<SAMRAI::pdat::IndexDataNode<3,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::List<SAMRAI::pdat::IndexDataNode<3,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::ListIterator<SAMRAI::pdat::IndexDataNode<3,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::ListNode<SAMRAI::pdat::IndexDataNode<3,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<3,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexVariable<3,IBAMR::LNodeIndexSet> >;
#endif

//////////////////////////////////////////////////////////////////////////////
