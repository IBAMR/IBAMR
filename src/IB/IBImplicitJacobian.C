// Filename: IBImplicitJacobian.C
// Created on 30 Aug 2010 by Boyce Griffith
//
// Copyright (c) 2002-2010 Boyce Griffith
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

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
