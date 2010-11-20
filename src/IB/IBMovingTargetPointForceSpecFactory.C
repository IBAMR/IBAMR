// Filename: IBMovingTargetPointForceSpecFactory.C
// Created on 14 Aug 2008 by Boyce Griffith
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

#include "IBMovingTargetPointForceSpecFactory.h"

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
#include <ibamr/IBMovingTargetPointForceSpec.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/StreamableManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

int IBMovingTargetPointForceSpecFactory::s_class_id = -1;

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBMovingTargetPointForceSpecFactory::IBMovingTargetPointForceSpecFactory()
{
    setStreamableClassID(StreamableManager::getUnregisteredID());
    return;
}// IBMovingTargetPointForceSpecFactory

IBMovingTargetPointForceSpecFactory::~IBMovingTargetPointForceSpecFactory()
{
    // intentionally blank
    return;
}// ~IBMovingTargetPointForceSpecFactory

int
IBMovingTargetPointForceSpecFactory::getStreamableClassID() const
{
    return s_class_id;
}// getStreamableClassID

void
IBMovingTargetPointForceSpecFactory::setStreamableClassID(
    const int class_id)
{
    s_class_id = class_id;
    return;
}// setStreamableClassID

Pointer<Streamable>
IBMovingTargetPointForceSpecFactory::unpackStream(
    AbstractStream& stream,
    const IntVector<NDIM>& offset)
{
    int mastr_idx;
    stream.unpack(&mastr_idx,1);
    double kappa_target;
    stream.unpack(&kappa_target,1);
    double eta_target;
    stream.unpack(&eta_target,1);
    int spec_fcn_idx;
    stream.unpack(&spec_fcn_idx,1);
    std::vector<double> periodic_shift(NDIM);
    stream.unpack(&periodic_shift[0],NDIM);
    return new IBMovingTargetPointForceSpec(mastr_idx,kappa_target,eta_target,spec_fcn_idx,periodic_shift);
}// unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBMovingTargetPointForceSpecFactory>;

//////////////////////////////////////////////////////////////////////////////
