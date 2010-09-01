// Filename: IBInstrumentationSpecFactory.C
// Created on 11 Jun 2007 by Boyce Griffith
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

#include "IBInstrumentationSpecFactory.h"

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
#include <ibamr/IBInstrumentationSpec.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/StashableManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

int IBInstrumentationSpecFactory::s_stashable_id = -1;

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBInstrumentationSpecFactory::IBInstrumentationSpecFactory()
{
    setStashableID(StashableManager::getUnregisteredID());
    return;
}// IBInstrumentationSpecFactory

IBInstrumentationSpecFactory::~IBInstrumentationSpecFactory()
{
    // intentionally blank
    return;
}// ~IBInstrumentationSpecFactory

int
IBInstrumentationSpecFactory::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

void
IBInstrumentationSpecFactory::setStashableID(
    const int stashable_id)
{
    s_stashable_id = stashable_id;
    return;
}// setStashableID

Pointer<Stashable>
IBInstrumentationSpecFactory::unpackStream(
    AbstractStream& stream,
    const IntVector<NDIM>& offset)
{
    int master_idx;
    stream.unpack(&master_idx,1);
    int meter_idx;
    stream.unpack(&meter_idx,1);
    int node_idx;
    stream.unpack(&node_idx,1);
    return new IBInstrumentationSpec(master_idx,meter_idx,node_idx);
}// unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBInstrumentationSpecFactory>;

//////////////////////////////////////////////////////////////////////////////
