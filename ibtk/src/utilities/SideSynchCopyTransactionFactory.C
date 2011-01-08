// Filename: SideSynchCopyTransactionFactory.C
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

#include "SideSynchCopyTransactionFactory.h"

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
#include <ibtk/SideSynchCopyTransaction.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <tbox/ArenaManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

SideSynchCopyTransactionFactory::SideSynchCopyTransactionFactory()
    : d_refine_items(NULL),
      d_number_refine_items(0)
{
    // intentionally blank
    return;
}// SideSynchCopyTransactionFactory

SideSynchCopyTransactionFactory::~SideSynchCopyTransactionFactory()
{
    // intentionally blank
}// ~SideSynchCopyTransactionFactory

void
SideSynchCopyTransactionFactory::setRefineItems(
    const RefineClasses<NDIM>::Data** refine_items,
    int num_refine_items)
{
    SideSynchCopyTransaction::setRefineItems(refine_items, num_refine_items);
    d_refine_items = refine_items;
    d_number_refine_items = num_refine_items;
    return;
}// setRefineItems

void
SideSynchCopyTransactionFactory::unsetRefineItems()
{
    SideSynchCopyTransaction::unsetRefineItems();
    d_refine_items = (const RefineClasses<NDIM>::Data**)NULL;
    d_number_refine_items = 0;
    return;
}// unsetRefineItems

Pointer<Transaction>
SideSynchCopyTransactionFactory::allocate(
    Pointer<PatchLevel<NDIM> > dst_level,
    Pointer<PatchLevel<NDIM> > src_level,
    Pointer<BoxOverlap<NDIM> > overlap,
    int dst_patch_id,
    int src_patch_id,
    int ritem_id,
    const Box<NDIM>& box,
    bool use_time_interpolation,
    Pointer<Arena> pool) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(dst_level == src_level);
#endif
    if (pool.isNull())
    {
        pool = ArenaManager::getManager()->getStandardAllocator();
    }
    SideSynchCopyTransaction* transaction = new (pool) SideSynchCopyTransaction(
        dst_level, src_level, overlap, dst_patch_id, src_patch_id, ritem_id);
    return Pointer<Transaction>(transaction, pool);
}// allocate

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
