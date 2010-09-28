// Filename: IBImplicitOperator.C
// Created on 26 Aug 2010 by Boyce Griffith
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

#include "IBImplicitOperator.h"

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

// IBTK INCLUDES
#include <ibtk/SideDataSynchronization.h> // XXXX

// SAMRAI INCLUDES
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

// C++ STDLIB INCLUDES
#include <limits>

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

IBImplicitOperator::IBImplicitOperator(
    Pointer<INSStaggeredStokesOperator> stokes_op,
    Pointer<IBImplicitSFROperator> ib_SFR_op)
    : d_is_initialized(false),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_dt(std::numeric_limits<double>::quiet_NaN()),
      d_stokes_op(stokes_op),
      d_ib_SFR_op(ib_SFR_op)
{
    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_apply                     = TimerManager::getManager()->getTimer("IBAMR::IBImplicitOperator::apply()");
        t_initialize_operator_state = TimerManager::getManager()->getTimer("IBAMR::IBImplicitOperator::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer("IBAMR::IBImplicitOperator::deallocateOperatorState()");
        timers_need_init = false;
    }
    return;
}// IBImplicitOperator

IBImplicitOperator::~IBImplicitOperator()
{
    deallocateOperatorState();
    return;
}// ~IBImplicitOperator

void
IBImplicitOperator::setTimeInterval(
    const double current_time,
    const double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    d_dt = d_new_time-d_current_time;
    d_stokes_op->setTimeInterval(current_time, new_time);
    d_ib_SFR_op->setTimeInterval(current_time, new_time);
    return;
}// setTimeInterval

void
IBImplicitOperator::apply(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& y)
{
    t_apply->start();

    // Apply the linear part of the operator with inhomogeneous boundary
    // conditions.
    static const bool homogeneous_bc = false;
    d_stokes_op->apply(homogeneous_bc, x, y);

    // Apply the nonlinear part of the operator.
    static const bool zero_y_before_spread = false;
    d_ib_SFR_op->apply(zero_y_before_spread, x, y);

    typedef SideDataSynchronization::SynchronizationTransactionComponent SynchronizationTransactionComponent; // XXXX
    const int y_idx = y.getComponentDescriptorIndex(0);
    SynchronizationTransactionComponent y_synch_transaction = SynchronizationTransactionComponent(y_idx, "CONSERVATIVE_COARSEN");
    Pointer<SideDataSynchronization> side_synch_op = new SideDataSynchronization();
    side_synch_op->initializeOperatorState(y_synch_transaction, y.getPatchHierarchy());
    side_synch_op->synchronizeData(0.0);

    t_apply->stop();
    return;
}// apply

void
IBImplicitOperator::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& in,
    const SAMRAIVectorReal<NDIM,double>& out)
{
    t_initialize_operator_state->start();

    if (d_is_initialized) deallocateOperatorState();

    d_stokes_op->initializeOperatorState(in, out);
    d_ib_SFR_op->initializeOperatorState(in, out);

    d_is_initialized = true;

    t_initialize_operator_state->stop();
    return;
}// initializeOperatorState

void
IBImplicitOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    t_deallocate_operator_state->start();

    d_stokes_op->deallocateOperatorState();
    d_ib_SFR_op->deallocateOperatorState();

    d_is_initialized = false;

    t_deallocate_operator_state->stop();
    return;
}// deallocateOperatorState

void
IBImplicitOperator::enableLogging(
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
template class Pointer<IBAMR::IBImplicitOperator>;

//////////////////////////////////////////////////////////////////////////////
