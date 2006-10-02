//
// LNodeIndexSet.C
//
// Created on 29 Feb 2004
//         by Boyce Griffith (boyce@bigboy.speakeasy.net).
//
// Last modified: <29.Jun.2005 17:11:51 boyce@mstu1.cims.nyu.edu>
//

#include "LNodeIndexSet.h"

// STL INCLUDES
//
#include <algorithm>
#include <functional>
#include <numeric>
#include <sstream>

// SAMRAI-tools INCLUDES
//
#include "SerializableStream.h"

/////////////////////////////// INLINE ///////////////////////////////////////

#ifdef DEBUG_NO_INLINE
#include "LNodeIndexSet.I"
#endif

/////////////////////////////// PUBLIC ///////////////////////////////////////

namespace
{
    struct LNodeIndexGetDataStreamSizeSum
        : binary_function<size_t,tbox::Pointer<LNodeIndex>,size_t>
    {
        size_t operator()(
            size_t size_so_far,
            const tbox::Pointer<LNodeIndex>& index) const
            {
                return(size_so_far+index->getDataStreamSize());
            }
    };        
}

size_t LNodeIndexSet::getDataStreamSize() const
{
    return(accumulate(d_set.begin(),d_set.end(),
                      tbox::AbstractStream::sizeofInt(),
                      LNodeIndexGetDataStreamSizeSum()));
}// getDataStreamSize

namespace
{
    class LNodeIndexPackStream
        : public unary_function<tbox::Pointer<LNodeIndex>,void>
    {
    public:
        LNodeIndexPackStream(
            tbox::AbstractStream* const stream)
            : d_stream(stream)
            {
                return;
            }

        void operator()(
            const tbox::Pointer<LNodeIndex>& index) const
            {
                index->packStream(*d_stream);
                return;
            }
    private:
        tbox::AbstractStream* const d_stream;
    };
}

void LNodeIndexSet::packStream(
    tbox::AbstractStream& stream)
{
    const int num_idx = d_set.size();
    stream.pack(&num_idx,1);
    for_each(d_set.begin(),d_set.end(), LNodeIndexPackStream(&stream));
    return;
}// packStream

namespace
{
    class LNodeIndexUnpackStream
        : public unary_function<void,tbox::Pointer<LNodeIndex> >
    {
    public:
        LNodeIndexUnpackStream(
            tbox::AbstractStream* const stream,
            const hier::IntVector<NDIM>& offset)
            : d_stream(stream),
              d_offset(offset)
            {
                return;
            }

        tbox::Pointer<LNodeIndex> operator()() const
            {
                tbox::Pointer<LNodeIndex> index_out = new LNodeIndex();
                index_out->unpackStream(*d_stream,d_offset);
                return(index_out);
            }
    private:
        tbox::AbstractStream* const d_stream;
        const hier::IntVector<NDIM>& d_offset;
    };
}

void LNodeIndexSet::unpackStream(
    tbox::AbstractStream& stream,
    const hier::IntVector<NDIM>& offset)
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

void LNodeIndexSet::putToDatabase(
    tbox::Pointer<tbox::Database>& database)
{
    const size_t data_sz = getDataStreamSize();
    static const bool use_xdr = false;
    SerializableStream stream(data_sz, SerializableStream::Write, use_xdr);
    packStream(stream);
    database->putInteger("data_sz", data_sz);
    database->putCharArray("data", (char*)stream.getBufferStart(), data_sz);
    database->putIntegerArray("d_offset", d_offset, NDIM);
    return;
}// putToDatabase

void LNodeIndexSet::getFromDatabase(
    tbox::Pointer<tbox::Database>& database)
{
    database->getIntegerArray("d_offset", d_offset, NDIM);
    
    const size_t data_sz = database->getInteger("data_sz");
    vector<char> data(data_sz);
    database->getCharArray("data", &data[0], data_sz);
    static const bool use_xdr = false;
    SerializableStream stream(&data[0], data_sz,
                              SerializableStream::Read, use_xdr);
    unpackStream(stream, d_offset);
    return;
}// getFromDatabase

/////////////////////////////// PRIVATE //////////////////////////////////////

namespace
{
    struct LNodeIndexLessThan
        : binary_function<tbox::Pointer<LNodeIndex>,tbox::Pointer<LNodeIndex>,bool>
    {
        bool operator()(
            const tbox::Pointer<LNodeIndex>& lhs,
            const tbox::Pointer<LNodeIndex>& rhs) const
            {
                return(*lhs < *rhs);
            }
    };        
}

void LNodeIndexSet::reorderCollection()
{
    sort(d_set.begin(),d_set.end(),LNodeIndexLessThan());
    return;
}// sortSet

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#ifndef LACKS_EXPLICIT_TEMPLATE_INSTANTIATION

#include "IndexData.C"
#include "IndexDataFactory.C"
#include "IndexVariable.C"
#include "tbox/Array.C"
#include "tbox/List.C"
#include "tbox/Pointer.C"

//////////////////////////////////////////////////////////////////////
///
/// These declarations are required to use the LNodeIndexSet class.
///
//////////////////////////////////////////////////////////////////////

template class tbox::Pointer<LNodeIndexSet>;

#if (NDIM == 1)
template class pdat::IndexData<1,LNodeIndexSet>;
template class pdat::IndexDataFactory<1,LNodeIndexSet>;
template class pdat::IndexDataNode<1,LNodeIndexSet>;
template class pdat::IndexIterator<1,LNodeIndexSet>;
template class pdat::IndexVariable<1,LNodeIndexSet>;
template class tbox::Array<LNodeIndexSet>;
template class tbox::Array<pdat::IndexDataNode<1,LNodeIndexSet> >;
template class tbox::List<pdat::IndexDataNode<1,LNodeIndexSet> >;
template class tbox::ListIterator<pdat::IndexDataNode<1,LNodeIndexSet> >;
template class tbox::ListNode<pdat::IndexDataNode<1,LNodeIndexSet> >;
template class tbox::Pointer<pdat::IndexData<1,LNodeIndexSet> >;
template class tbox::Pointer<pdat::IndexVariable<1,LNodeIndexSet> >;
#endif

#if (NDIM == 2)
template class pdat::IndexData<2,LNodeIndexSet>;
template class pdat::IndexDataFactory<2,LNodeIndexSet>;
template class pdat::IndexDataNode<2,LNodeIndexSet>;
template class pdat::IndexIterator<2,LNodeIndexSet>;
template class pdat::IndexVariable<2,LNodeIndexSet>;
template class tbox::Array<LNodeIndexSet>;
template class tbox::Array<pdat::IndexDataNode<2,LNodeIndexSet> >;
template class tbox::List<pdat::IndexDataNode<2,LNodeIndexSet> >;
template class tbox::ListIterator<pdat::IndexDataNode<2,LNodeIndexSet> >;
template class tbox::ListNode<pdat::IndexDataNode<2,LNodeIndexSet> >;
template class tbox::Pointer<pdat::IndexData<2,LNodeIndexSet> >;
template class tbox::Pointer<pdat::IndexVariable<2,LNodeIndexSet> >;
#endif

#if (NDIM == 3)
template class pdat::IndexData<3,LNodeIndexSet>;
template class pdat::IndexDataFactory<3,LNodeIndexSet>;
template class pdat::IndexDataNode<3,LNodeIndexSet>;
template class pdat::IndexIterator<3,LNodeIndexSet>;
template class pdat::IndexVariable<3,LNodeIndexSet>;
template class tbox::Array<LNodeIndexSet>;
template class tbox::Array<pdat::IndexDataNode<3,LNodeIndexSet> >;
template class tbox::List<pdat::IndexDataNode<3,LNodeIndexSet> >;
template class tbox::ListIterator<pdat::IndexDataNode<3,LNodeIndexSet> >;
template class tbox::ListNode<pdat::IndexDataNode<3,LNodeIndexSet> >;
template class tbox::Pointer<pdat::IndexData<3,LNodeIndexSet> >;
template class tbox::Pointer<pdat::IndexVariable<3,LNodeIndexSet> >;
#endif

#endif

//////////////////////////////////////////////////////////////////////////////
