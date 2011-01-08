// Filename: SideSynchCopyTransaction.C
// Created on 17 Dec 2009 by Boyce Griffith
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

#include "SideSynchCopyTransaction.h"

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

// SAMRAI INCLUDES
#include <tbox/SAMRAI_MPI.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

const RefineClasses<NDIM>::Data** SideSynchCopyTransaction::s_refine_items =
    (const RefineClasses<NDIM>::Data**)NULL;
int SideSynchCopyTransaction::s_num_refine_items = 0;

void
SideSynchCopyTransaction::setRefineItems(
    const RefineClasses<NDIM>::Data** refine_items,
    int num_refine_items)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(refine_items != (const RefineClasses<NDIM>::Data**)NULL);
    TBOX_ASSERT(num_refine_items >= 0);
#endif
    s_refine_items = refine_items;
    s_num_refine_items = num_refine_items;
    return;
}// setRefineItems

void
SideSynchCopyTransaction::unsetRefineItems()
{
    s_refine_items = (const RefineClasses<NDIM>::Data**)NULL;
    s_num_refine_items = 0;
    return;
}// unsetRefineItems

/////////////////////////////// PUBLIC ///////////////////////////////////////

SideSynchCopyTransaction::SideSynchCopyTransaction(
    Pointer<PatchLevel<NDIM> > dst_level,
    Pointer<PatchLevel<NDIM> > src_level,
    Pointer<BoxOverlap<NDIM> > overlap,
    int dst_patch,
    int src_patch,
    int refine_item_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!dst_level.isNull());
    TBOX_ASSERT(!src_level.isNull());
    TBOX_ASSERT(!overlap.isNull());
    TBOX_ASSERT(dst_patch >= 0 && dst_patch < dst_level->getNumberOfPatches());
    TBOX_ASSERT(src_patch >= 0 && src_patch < src_level->getNumberOfPatches());
    TBOX_ASSERT(refine_item_id >= 0);
    // Note: s_num_refine_items cannot be used at this point!
#endif
    d_dst_level        = dst_level;
    d_src_level        = src_level;
    d_overlap          = overlap;
    d_dst_patch        = dst_patch;
    d_src_patch        = src_patch;
    d_refine_item_id   = refine_item_id;
    d_incoming_bytes   = 0;
    d_outgoing_bytes   = 0;
    return;
}// SideSynchCopyTransaction

SideSynchCopyTransaction::~SideSynchCopyTransaction()
{
    // intentionally blank
    return;
}// ~SideSynchCopyTransaction

bool
SideSynchCopyTransaction::canEstimateIncomingMessageSize()
{
    bool can_estimate = false;
    if (getSourceProcessor() == SAMRAI_MPI::getRank())
    {
        can_estimate =
            d_src_level->getPatch(d_src_patch)
            ->getPatchData(s_refine_items[d_refine_item_id]->
                           d_src)
            ->canEstimateStreamSizeFromBox();
    }
    else
    {
        can_estimate =
            d_dst_level->getPatch(d_dst_patch)
            ->getPatchData(s_refine_items[d_refine_item_id]->
                           d_scratch)
            ->canEstimateStreamSizeFromBox();
    }
    return can_estimate;
}// canEstimateIncomingMessageSize

int
SideSynchCopyTransaction::computeIncomingMessageSize()
{
    d_incoming_bytes =
        d_dst_level->getPatch(d_dst_patch)
        ->getPatchData(s_refine_items[d_refine_item_id]->
                       d_scratch)
        ->getDataStreamSize(*d_overlap);
    return d_incoming_bytes;
}// computeIncomingMessageSize

int
SideSynchCopyTransaction::computeOutgoingMessageSize()
{
    d_outgoing_bytes =
        d_src_level->getPatch(d_src_patch)
        ->getPatchData(s_refine_items[d_refine_item_id]->
                       d_src)
        ->getDataStreamSize(*d_overlap);
    return d_outgoing_bytes;
}// computeOutgoingMessageSize

int
SideSynchCopyTransaction::getSourceProcessor()
{
    return d_src_level->getMappingForPatch(d_src_patch);
}// getSourceProcessor

int
SideSynchCopyTransaction::getDestinationProcessor()
{
    return d_dst_level->getMappingForPatch(d_dst_patch);
}// getDestinationProcessor

void
SideSynchCopyTransaction::packStream(
    AbstractStream& stream)
{
    d_src_level->getPatch(d_src_patch)
        ->getPatchData(s_refine_items[d_refine_item_id]->
                       d_src)
        ->packStream(stream, *d_overlap);
    return;
}// packStream

void
SideSynchCopyTransaction::unpackStream(
    AbstractStream& stream)
{
    Pointer<SideData<NDIM,double> > dst_data =
        d_dst_level->getPatch(d_dst_patch)->getPatchData(s_refine_items[d_refine_item_id]->d_scratch);
    SideData<NDIM,double> src_data(
        dst_data->getBox(), dst_data->getDepth(), dst_data->getGhostCellWidth());
    src_data.unpackStream(stream, *d_overlap);
    synchronizeData(dst_data, Pointer<SideData<NDIM,double> >(&src_data,false));
    return;
}// unpackStream

void
SideSynchCopyTransaction::copyLocalData()
{
    Pointer<SideData<NDIM,double> > dst_data =
        d_dst_level->getPatch(d_dst_patch)->getPatchData(s_refine_items[d_refine_item_id]->d_scratch);
    Pointer<SideData<NDIM,double> > src_data =
        d_src_level->getPatch(d_src_patch)->getPatchData(s_refine_items[d_refine_item_id]->d_src);
    synchronizeData(dst_data, src_data);
    return;
}// copyLocalData

void
SideSynchCopyTransaction::printClassData(
    std::ostream& stream) const
{
    stream << "Side Synch Copy Transaction"                       << std::endl;
    stream << "   refine item array:      "
           << (RefineClasses<NDIM>::Data**)s_refine_items
           << std::endl;
    stream << "   num refine items:       " << s_num_refine_items << std::endl;
    stream << "   destination patch:      " << d_dst_patch        << std::endl;
    stream << "   source patch:           " << d_src_patch        << std::endl;
    stream << "   refine item id:         " << d_refine_item_id   << std::endl;
    stream << "   destination patch data: "
           << s_refine_items[d_refine_item_id]->d_scratch         << std::endl;
    stream << "   source patch data:      "
           << s_refine_items[d_refine_item_id]->d_src             << std::endl;
    stream << "   incoming bytes:         " << d_incoming_bytes   << std::endl;
    stream << "   outgoing bytes:         " << d_outgoing_bytes   << std::endl;
    stream << "   destination level:      "
           << (PatchLevel<NDIM>*)d_src_level        << std::endl;
    stream << "   source level:           "
           << (PatchLevel<NDIM>*)d_src_level        << std::endl;
    stream << "   overlap:                "                       << std::endl;
    d_overlap->print(stream);
    return;
}// printClassData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
SideSynchCopyTransaction::synchronizeData(
    Pointer<SideData<NDIM,double> > dst_data,
    Pointer<SideData<NDIM,double> > src_data)
{
    Pointer<SideOverlap<NDIM> > overlap = d_overlap;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!dst_data.isNull());
    TBOX_ASSERT(!src_data.isNull());
    TBOX_ASSERT(!overlap.isNull());
#endif
    const IntVector<NDIM>& src_offset = overlap->getSourceOffset();
    for (int d = 0; d < NDIM; ++d)
    {
        const BoxList<NDIM>& box_list = overlap->getDestinationBoxList(d);
        dst_data->getArrayData(d).copy(src_data->getArrayData(d), box_list, src_offset);
    }
    return;
}// synchronizeData

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
