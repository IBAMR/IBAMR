// Filename: IBImplicitModHelmholtzOperator.C
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

#include "IBImplicitModHelmholtzOperator.h"

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
#include <ibamr/namespaces.h>

// SAMRAI INCLUDES
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Pointer<Timer> t_apply;
static Pointer<Timer> t_initialize_operator_state;
static Pointer<Timer> t_deallocate_operator_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBImplicitModHelmholtzOperator::IBImplicitModHelmholtzOperator(
    Pointer<SCLaplaceOperator> helmholtz_op,
    Pointer<LinearOperator> ib_SJSstar_op)
    : d_is_initialized(false),
      d_helmholtz_op(helmholtz_op),
      d_ib_SJSstar_op(ib_SJSstar_op)
{
    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_apply                     = TimerManager::getManager()->getTimer("IBAMR::IBImplicitModHelmholtzOperator::apply()");
        t_initialize_operator_state = TimerManager::getManager()->getTimer("IBAMR::IBImplicitModHelmholtzOperator::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer("IBAMR::IBImplicitModHelmholtzOperator::deallocateOperatorState()");
        timers_need_init = false;
    }
    return;
}// IBImplicitModHelmholtzOperator

IBImplicitModHelmholtzOperator::~IBImplicitModHelmholtzOperator()
{
    deallocateOperatorState();
    return;
}// ~IBImplicitModHelmholtzOperator

void
IBImplicitModHelmholtzOperator::apply(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& y)
{
    t_apply->start();

    SAMRAIVectorReal<NDIM,double> x_u(x.getName(), x.getPatchHierarchy(), x.getCoarsestLevelNumber(), x.getFinestLevelNumber());
    x_u.addComponent(x.getComponentVariable(0), x.getComponentDescriptorIndex(0), x.getControlVolumeIndex(0));
    SAMRAIVectorReal<NDIM,double> y_u(y.getName(), y.getPatchHierarchy(), y.getCoarsestLevelNumber(), y.getFinestLevelNumber());
    y_u.addComponent(y.getComponentVariable(0), y.getComponentDescriptorIndex(0), y.getControlVolumeIndex(0));

    // Apply the linear part of the operator with homogeneous boundary
    // conditions.
    d_helmholtz_op->apply(x_u, y_u);

    // Apply the nonlinear part of the operator.
    d_ib_SJSstar_op->applyAdd(x_u, y_u, y_u);

    t_apply->stop();
    return;
}// apply

void
IBImplicitModHelmholtzOperator::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& in,
    const SAMRAIVectorReal<NDIM,double>& out)
{
    t_initialize_operator_state->start();

    if (d_is_initialized) deallocateOperatorState();

    SAMRAIVectorReal<NDIM,double> in_u(in.getName(), in.getPatchHierarchy(), in.getCoarsestLevelNumber(), in.getFinestLevelNumber());
    in_u.addComponent(in.getComponentVariable(0), in.getComponentDescriptorIndex(0), in.getControlVolumeIndex(0));
    SAMRAIVectorReal<NDIM,double> out_u(out.getName(), out.getPatchHierarchy(), out.getCoarsestLevelNumber(), out.getFinestLevelNumber());
    out_u.addComponent(out.getComponentVariable(0), out.getComponentDescriptorIndex(0), out.getControlVolumeIndex(0));

    d_helmholtz_op->initializeOperatorState(in_u, out_u);
//  d_ib_SJSstar_op->initializeOperatorState(in_u, out_u);

    d_is_initialized = true;

    t_initialize_operator_state->stop();
    return;
}// initializeOperatorState

void
IBImplicitModHelmholtzOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    t_deallocate_operator_state->start();

    d_helmholtz_op->deallocateOperatorState();
//  d_ib_SJSstar_op->deallocateOperatorState();

    d_is_initialized = false;

    t_deallocate_operator_state->stop();
    return;
}// deallocateOperatorState

void
IBImplicitModHelmholtzOperator::enableLogging(
    bool enabled)
{
    // intentionally blank
    return;
}// enableLogging

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBImplicitModHelmholtzOperator>;

//////////////////////////////////////////////////////////////////////////////
