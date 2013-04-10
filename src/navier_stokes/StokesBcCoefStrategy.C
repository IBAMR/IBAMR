// Filename: StokesBcCoefStrategy.C
// Created on 04 Sep 2012 by Boyce Griffith
//
// Copyright (c) 2002-2013, Boyce Griffith
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

#include "StokesBcCoefStrategy.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

StokesBcCoefStrategy::StokesBcCoefStrategy()
    : d_problem_coefs(NULL),
      d_u_target_data_idx(-1),
      d_p_target_data_idx(-1)
{
    // intentionally blank
    return;
}// StokesBcCoefStrategy

StokesBcCoefStrategy::~StokesBcCoefStrategy()
{
    // intentionally blank
    return;
}// ~StokesBcCoefStrategy

void
StokesBcCoefStrategy::setStokesSpecifications(
    const StokesSpecifications* const problem_coefs)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(problem_coefs);
#endif
    d_problem_coefs = problem_coefs;
    return;
}// setStokesSpecifications

void
StokesBcCoefStrategy::setTargetVelocityPatchDataIndex(
    int u_target_data_idx)
{
    d_u_target_data_idx = u_target_data_idx;
    return;
}// setTargetVelocityPatchDataIndex

void
StokesBcCoefStrategy::clearTargetVelocityPatchDataIndex()
{
    d_u_target_data_idx = -1;
    return;
}// clearTargetVelocityPatchDataIndex

void
StokesBcCoefStrategy::setTargetPressurePatchDataIndex(
    int p_target_data_idx)
{
    d_p_target_data_idx = p_target_data_idx;
    return;
}// setPressurePatchDataIndex

void
StokesBcCoefStrategy::clearTargetPressurePatchDataIndex()
{
    d_p_target_data_idx = -1;
    return;
}// clearPressurePatchDataIndex

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
