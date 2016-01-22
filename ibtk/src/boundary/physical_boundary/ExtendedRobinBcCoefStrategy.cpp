// Filename: ExtendedRobinBcCoefStrategy.cpp
// Created on 16 May 2007 by Boyce Griffith
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

#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

ExtendedRobinBcCoefStrategy::ExtendedRobinBcCoefStrategy() : d_target_data_idx(-1), d_homogeneous_bc(false)
{
    // intentionally blank
    return;
} // ExtendedRobinBcCoefStrategy

ExtendedRobinBcCoefStrategy::~ExtendedRobinBcCoefStrategy()
{
    // intentionally blank
    return;
} // ~ExtendedRobinBcCoefStrategy

void
ExtendedRobinBcCoefStrategy::setTargetPatchDataIndex(int target_data_idx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(target_data_idx >= 0);
#endif
    d_target_data_idx = target_data_idx;
    return;
} // setTargetPatchDataIndex

void
ExtendedRobinBcCoefStrategy::clearTargetPatchDataIndex()
{
    d_target_data_idx = -1;
    return;
} // clearTargetPatchDataIndex

void
ExtendedRobinBcCoefStrategy::setHomogeneousBc(bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
} // setHomogeneousBc

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
