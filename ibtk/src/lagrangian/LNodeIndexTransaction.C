// Filename: LNodeIndexTransaction.C
// Created on 03 Mar 2010 by Boyce Griffith
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

#include "LNodeIndexTransaction.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/namespaces.h>

// C++ STDLIB INCLUDES
#include <algorithm>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LNodeIndexTransaction::LNodeIndexTransaction(
    const int src_proc,
    const int dst_proc)
    : d_src_index_set(),
      d_src_proc(src_proc),
      d_outgoing_bytes(0),
      d_dst_index_set(),
      d_dst_proc(dst_proc)
{
    // intentionally blank
    return;
}// LNodeIndexTransaction

LNodeIndexTransaction::LNodeIndexTransaction(
    const int src_proc,
    const int dst_proc,
    const std::vector<std::pair<Pointer<LNodeIndex>,std::vector<double> > >& src_index_set)
    : d_src_index_set(src_index_set),
      d_src_proc(src_proc),
      d_outgoing_bytes(0),
      d_dst_index_set(),
      d_dst_proc(dst_proc)
{
    d_outgoing_bytes = AbstractStream::sizeofInt();
    for (std::vector<std::pair<Pointer<LNodeIndex>,std::vector<double> > >::const_iterator cit = d_src_index_set.begin();
         cit != d_src_index_set.end(); ++cit)
    {
        d_outgoing_bytes += (*cit).first->getDataStreamSize() + NDIM*AbstractStream::sizeofDouble();
    }
    return;
}// LNodeIndexTransaction

LNodeIndexTransaction::~LNodeIndexTransaction()
{
    // intentionally blank
    return;
}// ~LNodeIndexTransaction

bool
LNodeIndexTransaction::canEstimateIncomingMessageSize()
{
    return false;
}// canEstimateIncomingMessageSize

int
LNodeIndexTransaction::computeIncomingMessageSize()
{
    return 0;
}// computeIncomingMessageSize

int
LNodeIndexTransaction::computeOutgoingMessageSize()
{
    return d_outgoing_bytes;
}// computeOutgoingMessageSize

int
LNodeIndexTransaction::getSourceProcessor()
{
    return d_src_proc;
}// getSourceProcessor

int
LNodeIndexTransaction::getDestinationProcessor()
{
    return d_dst_proc;
}// getDestinationProcessor

void
LNodeIndexTransaction::packStream(
    AbstractStream& stream)
{
    stream << int(d_src_index_set.size());
    for (std::vector<std::pair<Pointer<LNodeIndex>,std::vector<double> > >::const_iterator cit = d_src_index_set.begin();
         cit != d_src_index_set.end(); ++cit)
    {
        const Pointer<LNodeIndex>& idx = (*cit).first;
        idx->packStream(stream);
        const std::vector<double>& posn = (*cit).second;
        stream.pack(&posn[0],NDIM);
    }
    return;
}// packStream

void
LNodeIndexTransaction::unpackStream(
    AbstractStream& stream)
{
    static const IntVector<NDIM> periodic_offset = 0;
    int num_idxs;
    stream >> num_idxs;
    d_dst_index_set.resize(num_idxs,std::make_pair(Pointer<LNodeIndex>(NULL),std::vector<double>(NDIM)));
    for (std::vector<std::pair<Pointer<LNodeIndex>,std::vector<double> > >::iterator it = d_dst_index_set.begin();
         it != d_dst_index_set.end(); ++it)
    {
        (*it).first = new LNodeIndex();
        Pointer<LNodeIndex> idx = (*it).first;
        idx->unpackStream(stream, periodic_offset);
        std::vector<double>& posn = (*it).second;
        stream.unpack(&posn[0],NDIM);
    }
    return;
}// unpackStream

void
LNodeIndexTransaction::copyLocalData()
{
    d_dst_index_set = d_src_index_set;
    return;
}// copyLocalData

void
LNodeIndexTransaction::printClassData(
    std::ostream& stream) const
{
    stream << "LNodeIndex Transaction"                                    << std::endl;
    stream << "   number of outgoing indices: " << d_src_index_set.size() << std::endl;
    stream << "   outgoing processor rank:    " << d_src_proc             << std::endl;
    stream << "   outgoing bytes:             " << d_outgoing_bytes       << std::endl;
    stream << "   number of incoming indices: " << d_dst_index_set.size() << std::endl;
    stream << "   incoming processor rank:    " << d_dst_proc             << std::endl;
    return;
}// printClassData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
