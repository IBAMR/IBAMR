// Filename: IBInstrumentationSpecFactory.cpp
// Created on 11 Jun 2007 by Boyce Griffith
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

#include "ibamr/IBInstrumentationSpec.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/Streamable.h"
#include "ibtk/StreamableManager.h"
#include "tbox/AbstractStream.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class IntVector;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

IBInstrumentationSpec::Factory::Factory()
{
    setStreamableClassID(StreamableManager::getUnregisteredID());
    return;
} // Factory

IBInstrumentationSpec::Factory::~Factory()
{
    // intentionally blank
    return;
} // ~Factory

int
IBInstrumentationSpec::Factory::getStreamableClassID() const
{
    return STREAMABLE_CLASS_ID;
} // getStreamableClassID

void
IBInstrumentationSpec::Factory::setStreamableClassID(const int class_id)
{
    STREAMABLE_CLASS_ID = class_id;
    return;
} // setStreamableClassID

Pointer<Streamable>
IBInstrumentationSpec::Factory::unpackStream(AbstractStream& stream, const IntVector<NDIM>& /*offset*/)
{
    Pointer<IBInstrumentationSpec> ret_val = new IBInstrumentationSpec();
    stream.unpack(&ret_val->d_master_idx, 1);
    stream.unpack(&ret_val->d_meter_idx, 1);
    stream.unpack(&ret_val->d_node_idx, 1);
    return ret_val;
} // unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
