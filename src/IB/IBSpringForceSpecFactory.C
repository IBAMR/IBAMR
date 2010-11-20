// Filename: IBSpringForceSpecFactory.C
// Created on 14 Jul 2004 by Boyce Griffith
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

#include "IBSpringForceSpecFactory.h"

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
#include <ibamr/IBSpringForceSpec.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/StreamableManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

int IBSpringForceSpecFactory::s_class_id = -1;

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBSpringForceSpecFactory::IBSpringForceSpecFactory()
{
    setStreamableClassID(StreamableManager::getUnregisteredID());
    return;
}// IBSpringForceSpecFactory

IBSpringForceSpecFactory::~IBSpringForceSpecFactory()
{
    // intentionally blank
    return;
}// ~IBSpringForceSpecFactory

int
IBSpringForceSpecFactory::getStreamableClassID() const
{
    return s_class_id;
}// getStreamableClassID

void
IBSpringForceSpecFactory::setStreamableClassID(
    const int class_id)
{
    s_class_id = class_id;
    return;
}// setStreamableClassID

Pointer<Streamable>
IBSpringForceSpecFactory::unpackStream(
    AbstractStream& stream,
    const IntVector<NDIM>& offset)
{
    int num_springs;
    stream.unpack(&num_springs,1);
    int master_idx;
    stream.unpack(&master_idx,1);
    std::vector<int> slave_idxs(num_springs);
    stream.unpack(&slave_idxs[0],num_springs);
    std::vector<int> force_fcn_idxs(num_springs);
    stream.unpack(&force_fcn_idxs[0],num_springs);
    std::vector<double> stiffnesses(num_springs);
    stream.unpack(&stiffnesses[0],num_springs);
    std::vector<double> rest_lengths(num_springs);
    stream.unpack(&rest_lengths[0],num_springs);
    return new IBSpringForceSpec(master_idx,slave_idxs,force_fcn_idxs,stiffnesses,rest_lengths);
}// unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBSpringForceSpecFactory>;

//////////////////////////////////////////////////////////////////////////////
