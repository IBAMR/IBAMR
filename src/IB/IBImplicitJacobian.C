// Filename: IBImplicitJacobian.C
// Created on 30 Aug 2010 by Boyce Griffith
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

#include "IBImplicitJacobian.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScSAMRAIVectorReal.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBImplicitJacobian::IBImplicitJacobian(
    Pointer<INSStaggeredStokesOperator> stokes_op,
    Pointer<JacobianOperator> ib_SJSstar_op)
    : JacobianOperator(false),
      d_is_initialized(false),
      d_stokes_op(stokes_op),
      d_ib_SJSstar_op(ib_SJSstar_op)
{
    // intentionally blank
    return;
}// IBImplicitJacobian()

IBImplicitJacobian::~IBImplicitJacobian()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
}// ~IBImplicitJacobian()

void
IBImplicitJacobian::formJacobian(
    SAMRAIVectorReal<NDIM,double>& x)
{
    d_ib_SJSstar_op->formJacobian(x);
    return;
}// formJacobiam

Pointer<SAMRAIVectorReal<NDIM,double> >
IBImplicitJacobian::getBaseVector() const
{
    return d_ib_SJSstar_op->getBaseVector();
}// getBaseVector

void
IBImplicitJacobian::apply(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& y)
{
    if (!d_is_initialized) initializeOperatorState(x, y);

    // Apply the linear part of the operator with homogeneous boundary
    // conditions.
    static const bool homogeneous_bc = true;
    d_stokes_op->apply(homogeneous_bc, x, y);

    // Apply the nonlinear part of the operator.
    SAMRAIVectorReal<NDIM,double> x_u(x.getName(), x.getPatchHierarchy(), x.getCoarsestLevelNumber(), x.getFinestLevelNumber());
    x_u.addComponent(x.getComponentVariable(0), x.getComponentDescriptorIndex(0), x.getControlVolumeIndex(0));
    SAMRAIVectorReal<NDIM,double> y_u(y.getName(), y.getPatchHierarchy(), y.getCoarsestLevelNumber(), y.getFinestLevelNumber());
    y_u.addComponent(y.getComponentVariable(0), y.getComponentDescriptorIndex(0), y.getControlVolumeIndex(0));
    d_ib_SJSstar_op->applyAdd(x_u, y_u, y_u);
    return;
}// apply

void
IBImplicitJacobian::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& in,
    const SAMRAIVectorReal<NDIM,double>& out)
{
    if (d_is_initialized) deallocateOperatorState();

    d_stokes_op->initializeOperatorState(in, out);

    SAMRAIVectorReal<NDIM,double> in_u(in.getName(), in.getPatchHierarchy(), in.getCoarsestLevelNumber(), in.getFinestLevelNumber());
    in_u.addComponent(in.getComponentVariable(0), in.getComponentDescriptorIndex(0), in.getControlVolumeIndex(0));
    SAMRAIVectorReal<NDIM,double> out_u(out.getName(), out.getPatchHierarchy(), out.getCoarsestLevelNumber(), out.getFinestLevelNumber());
    out_u.addComponent(out.getComponentVariable(0), out.getComponentDescriptorIndex(0), out.getControlVolumeIndex(0));
    d_ib_SJSstar_op->initializeOperatorState(in_u, out_u);

    d_is_initialized = true;
    return;
}// initializeOperatorState

void
IBImplicitJacobian::deallocateOperatorState()
{
    if (!d_is_initialized) return;
    d_stokes_op->deallocateOperatorState();
    d_ib_SJSstar_op->deallocateOperatorState();
    d_is_initialized = false;
    return;
}// deallocateOperatorState

void
IBImplicitJacobian::enableLogging(
    bool enabled)
{
    // intentionally blank
    return;
}// enableLogging

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBImplicitJacobian>;

//////////////////////////////////////////////////////////////////////////////
