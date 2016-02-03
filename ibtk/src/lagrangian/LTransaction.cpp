// Filename: LTransaction.cpp
// Created on 03 Mar 2010 by Boyce Griffith
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

#include <ostream>
#include <vector>

#include "IntVector.h"
#include "ibtk/LSet.h"
#include "ibtk/LTransaction.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/AbstractStream.h"

namespace IBTK
{
class LMarker;
class LNode;
class LNodeIndex;
} // namespace IBTK

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

template <class T>
LTransaction<T>::LTransaction(const int src_proc, const int dst_proc)
    : d_src_item_set(), d_src_proc(src_proc), d_outgoing_bytes(0), d_dst_item_set(), d_dst_proc(dst_proc)
{
    // intentionally blank
    return;
} // LTransaction

template <class T>
LTransaction<T>::LTransaction(const int src_proc,
                              const int dst_proc,
                              const std::vector<LTransactionComponent>& src_item_set)
    : d_src_item_set(src_item_set), d_src_proc(src_proc), d_outgoing_bytes(0), d_dst_item_set(), d_dst_proc(dst_proc)
{
    d_outgoing_bytes = AbstractStream::sizeofInt();
    for (typename std::vector<LTransactionComponent>::const_iterator cit = d_src_item_set.begin();
         cit != d_src_item_set.end();
         ++cit)
    {
        d_outgoing_bytes += cit->item->getDataStreamSize() + NDIM * AbstractStream::sizeofDouble();
    }
    return;
} // LTransaction

template <class T>
LTransaction<T>::~LTransaction()
{
    // intentionally blank
    return;
} // ~LTransaction

template <class T>
bool
LTransaction<T>::canEstimateIncomingMessageSize()
{
    return false;
} // canEstimateIncomingMessageSize

template <class T>
int
LTransaction<T>::computeIncomingMessageSize()
{
    return 0;
} // computeIncomingMessageSize

template <class T>
int
LTransaction<T>::computeOutgoingMessageSize()
{
    return d_outgoing_bytes;
} // computeOutgoingMessageSize

template <class T>
int
LTransaction<T>::getSourceProcessor()
{
    return d_src_proc;
} // getSourceProcessor

template <class T>
int
LTransaction<T>::getDestinationProcessor()
{
    return d_dst_proc;
} // getDestinationProcessor

template <class T>
void
LTransaction<T>::packStream(AbstractStream& stream)
{
    stream << static_cast<int>(d_src_item_set.size());
    for (typename std::vector<LTransactionComponent>::iterator it = d_src_item_set.begin(); it != d_src_item_set.end();
         ++it)
    {
        typename LSet<T>::value_type& item = it->item;
        item->packStream(stream);
        const Point& posn = it->posn;
        stream.pack(posn.data(), NDIM);
    }
    return;
} // packStream

template <class T>
void
LTransaction<T>::unpackStream(AbstractStream& stream)
{
    static const IntVector<NDIM> periodic_offset = 0;
    int num_items;
    stream >> num_items;
    d_dst_item_set.resize(num_items);
    for (typename std::vector<LTransactionComponent>::iterator it = d_dst_item_set.begin(); it != d_dst_item_set.end();
         ++it)
    {
        it->item->unpackStream(stream, periodic_offset);
        Point& posn = it->posn;
        stream.unpack(posn.data(), NDIM);
    }
    return;
} // unpackStream

template <class T>
void
LTransaction<T>::copyLocalData()
{
    d_dst_item_set = d_src_item_set;
    return;
} // copyLocalData

template <class T>
void
LTransaction<T>::printClassData(std::ostream& stream) const
{
    stream << "LNodeIndex Transaction" << std::endl;
    stream << "   number of outgoing indices: " << d_src_item_set.size() << std::endl;
    stream << "   outgoing processor rank:    " << d_src_proc << std::endl;
    stream << "   outgoing bytes:             " << d_outgoing_bytes << std::endl;
    stream << "   number of incoming indices: " << d_dst_item_set.size() << std::endl;
    stream << "   incoming processor rank:    " << d_dst_proc << std::endl;
    return;
} // printClassData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

template class IBTK::LTransaction<IBTK::LMarker>;
template class IBTK::LTransaction<IBTK::LNode>;
template class IBTK::LTransaction<IBTK::LNodeIndex>;

//////////////////////////////////////////////////////////////////////////////
