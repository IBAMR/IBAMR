// Filename: GeneralOperator.C
// Created on 18 Nov 2003 by Boyce Griffith
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

#include "GeneralOperator.h"

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

GeneralOperator::GeneralOperator()
    : d_hier_math_ops(NULL),
      d_hier_math_ops_external(false)
{
    // intentionally blank
    return;
}// GeneralOperator()

GeneralOperator::~GeneralOperator()
{
    deallocateOperatorState();
    return;
}// ~GeneralOperator()

void
GeneralOperator::setHierarchyMathOps(
    Pointer<HierarchyMathOps> hier_math_ops)
{
    d_hier_math_ops = hier_math_ops;
    d_hier_math_ops_external = !d_hier_math_ops.isNull();
    return;
}// setHierarchyMathOps

void
GeneralOperator::applyAdd(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& y,
    SAMRAIVectorReal<NDIM,double>& z)
{
    // Guard against the case that y == z.
    Pointer<SAMRAIVectorReal<NDIM,double> > zz = z.cloneVector(z.getName());
    zz->allocateVectorData();
    zz->copyVector(Pointer<SAMRAIVectorReal<NDIM,double> >(&z,false));
    apply(x,*zz);
    z.add(Pointer<SAMRAIVectorReal<NDIM,double> >(&y,false),zz);
    zz->freeVectorComponents();
    zz.setNull();
    return;
}// applyAdd

void
GeneralOperator::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& /*in*/,
    const SAMRAIVectorReal<NDIM,double>& /*out*/)
{
    // intentionally blank
    return;
}// initializeOperatorState

void
GeneralOperator::deallocateOperatorState()
{
    // intentionally blank
    return;
}// deallocateOperatorState

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
