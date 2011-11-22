// Filename: IBSourceSpecFactory.C
// Created on 28 Apr 2011 by Boyce Griffith
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

#include "IBSourceSpecFactory.h"

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
#include <ibamr/IBSourceSpec.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/StreamableManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

int IBSourceSpecFactory::s_class_id = StreamableManager::getUnregisteredID();

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBSourceSpecFactory::IBSourceSpecFactory()
{
    setStreamableClassID(StreamableManager::getUnregisteredID());
    return;
}// IBSourceSpecFactory

IBSourceSpecFactory::~IBSourceSpecFactory()
{
    // intentionally blank
    return;
}// ~IBSourceSpecFactory

int
IBSourceSpecFactory::getStreamableClassID() const
{
    return s_class_id;
}// getStreamableClassID

void
IBSourceSpecFactory::setStreamableClassID(
    const int class_id)
{
    s_class_id = class_id;
    return;
}// setStreamableClassID

Pointer<Streamable>
IBSourceSpecFactory::unpackStream(
    AbstractStream& stream,
    const IntVector<NDIM>& /*offset*/)
{
    Pointer<IBSourceSpec> ret_val = new IBSourceSpec();
    stream.unpack(&ret_val->d_master_idx,1);
    stream.unpack(&ret_val->d_source_idx,1);
    return ret_val;
}// unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
