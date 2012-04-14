// Filename: LinearOperator.C
// Created on 14 Sep 2003 by Boyce Griffith
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

#include "LinearOperator.h"

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

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LinearOperator::LinearOperator(
    bool homogeneous_bc)
    : d_homogeneous_bc(homogeneous_bc),
      d_correcting_rhs(false)
{
    // intentionally blank
    return;
}// LinearOperator()

LinearOperator::~LinearOperator()
{
    deallocateOperatorState();
    return;
}// ~LinearOperator()

void
LinearOperator::setHomogeneousBc(
    const bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

bool
LinearOperator::getHomogeneousBc()
{
    return d_homogeneous_bc;
}// getHomogeneousBc

void
LinearOperator::modifyRhsForInhomogeneousBc(
    SAMRAIVectorReal<NDIM,double>& /*y*/)
{
    TBOX_WARNING("LinearOperator::modifyRhsForInhomogeneousBc() not implemented for this operator" << std::endl);
    return;
}// modifyRhsForInhomogeneousBc

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
