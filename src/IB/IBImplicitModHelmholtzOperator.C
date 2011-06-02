// Filename: IBImplicitModHelmholtzOperator.C
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
#include <ibamr/ibamr_utilities.h>
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
static Timer* t_apply;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
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
    IBAMR_DO_ONCE(
        t_apply                     = TimerManager::getManager()->getTimer("IBAMR::IBImplicitModHelmholtzOperator::apply()");
        t_initialize_operator_state = TimerManager::getManager()->getTimer("IBAMR::IBImplicitModHelmholtzOperator::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer("IBAMR::IBImplicitModHelmholtzOperator::deallocateOperatorState()");
                  );
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
    IBAMR_TIMER_START(t_apply);

    SAMRAIVectorReal<NDIM,double> x_u(x.getName(), x.getPatchHierarchy(), x.getCoarsestLevelNumber(), x.getFinestLevelNumber());
    x_u.addComponent(x.getComponentVariable(0), x.getComponentDescriptorIndex(0), x.getControlVolumeIndex(0));
    SAMRAIVectorReal<NDIM,double> y_u(y.getName(), y.getPatchHierarchy(), y.getCoarsestLevelNumber(), y.getFinestLevelNumber());
    y_u.addComponent(y.getComponentVariable(0), y.getComponentDescriptorIndex(0), y.getControlVolumeIndex(0));

    // Apply the linear part of the operator with homogeneous boundary
    // conditions.
    d_helmholtz_op->apply(x_u, y_u);

    // Apply the nonlinear part of the operator.
    d_ib_SJSstar_op->applyAdd(x_u, y_u, y_u);

    IBAMR_TIMER_STOP(t_apply);
    return;
}// apply

void
IBImplicitModHelmholtzOperator::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& in,
    const SAMRAIVectorReal<NDIM,double>& out)
{
    IBAMR_TIMER_START(t_initialize_operator_state);

    if (d_is_initialized) deallocateOperatorState();

    SAMRAIVectorReal<NDIM,double> in_u(in.getName(), in.getPatchHierarchy(), in.getCoarsestLevelNumber(), in.getFinestLevelNumber());
    in_u.addComponent(in.getComponentVariable(0), in.getComponentDescriptorIndex(0), in.getControlVolumeIndex(0));
    SAMRAIVectorReal<NDIM,double> out_u(out.getName(), out.getPatchHierarchy(), out.getCoarsestLevelNumber(), out.getFinestLevelNumber());
    out_u.addComponent(out.getComponentVariable(0), out.getComponentDescriptorIndex(0), out.getControlVolumeIndex(0));

    d_helmholtz_op->initializeOperatorState(in_u, out_u);
//  d_ib_SJSstar_op->initializeOperatorState(in_u, out_u);

    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_operator_state);
    return;
}// initializeOperatorState

void
IBImplicitModHelmholtzOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_operator_state);

    d_helmholtz_op->deallocateOperatorState();
//  d_ib_SJSstar_op->deallocateOperatorState();

    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_operator_state);
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
