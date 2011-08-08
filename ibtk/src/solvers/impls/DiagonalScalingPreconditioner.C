// Filename: DiagonalScalingPreconditioner.C
// Created on 18 Apr 2009 by Boyce Griffith
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

#include "DiagonalScalingPreconditioner.h"

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

// C++ STDLIB INCLUDES
#include <ostream>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

DiagonalScalingPreconditioner::DiagonalScalingPreconditioner(
    const double alpha)
    : d_alpha(alpha)
{
    // intentionally blank
    return;
}// DiagonalScalingPreconditioner()

DiagonalScalingPreconditioner::~DiagonalScalingPreconditioner()
{
    // intentionally blank
    return;
}// ~DiagonalScalingPreconditioner()

void
DiagonalScalingPreconditioner::setScalingFactor(
    const double alpha)
{
    d_alpha = alpha;
    return;
}// setScalingFactor

bool
DiagonalScalingPreconditioner::solveSystem(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& b)
{
    x.scale(d_alpha,Pointer<SAMRAIVectorReal<NDIM,double> >(&b,false));
    return true;
}// solveSystem

void
DiagonalScalingPreconditioner::initializeSolverState(
    const SAMRAIVectorReal<NDIM,double>& /*x*/,
    const SAMRAIVectorReal<NDIM,double>& /*b*/)
{
    // intentionally blank
    return;
}// initializeSolverState

void
DiagonalScalingPreconditioner::deallocateSolverState()
{
    // intentionally blank
    return;
}// deallocateSolverState

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBTK::DiagonalScalingPreconditioner>;

//////////////////////////////////////////////////////////////////////////////
