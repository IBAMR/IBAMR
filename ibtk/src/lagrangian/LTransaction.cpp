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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/LNode.h"
#include "ibtk/LNodeIndex.h"
#include "ibtk/LSet.h"
#include "ibtk/LTransaction.h"
#include "ibtk/ibtk_utilities.h"

#include "IntVector.h"
#include "tbox/AbstractStream.h"

#include <algorithm>
#include <ostream>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

template <class T>
LTransaction<T>::LTransaction(const int src_proc, const int dst_proc) : d_src_proc(src_proc), d_dst_proc(dst_proc)
{
    // intentionally blank
    return;
} // LTransaction

template <class T>
LTransaction<T>::LTransaction(const int src_proc, const int dst_proc, std::vector<LTransactionComponent> src_item_set)
    : d_src_item_set(std::move(src_item_set)), d_src_proc(src_proc), d_dst_proc(dst_proc)
{
    d_outgoing_bytes = AbstractStream::sizeofInt();
    for (auto cit = d_src_item_set.begin(); cit != d_src_item_set.end(); ++cit)
    {
        d_outgoing_bytes += cit->item->getDataStreamSize() + NDIM * AbstractStream::sizeofDouble();
    }
    return;
} // LTransaction

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
    for (auto it = d_src_item_set.begin(); it != d_src_item_set.end(); ++it)
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
    for (auto it = d_dst_item_set.begin(); it != d_dst_item_set.end(); ++it)
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

template class IBTK::LTransaction<IBTK::LNode>;
template class IBTK::LTransaction<IBTK::LNodeIndex>;

//////////////////////////////////////////////////////////////////////////////
