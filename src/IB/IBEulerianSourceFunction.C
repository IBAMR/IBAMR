// Filename: IBEulerianSourceFunction.C
// Created on 18 Jun 2005 by Boyce Griffith
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

#include "IBEulerianSourceFunction.h"

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
#include <CellData.h>
#include <tbox/MathUtilities.h>

// C++ STDLIB INCLUDES
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

////////////////////////////// PUBLIC ///////////////////////////////////////

IBEulerianSourceFunction::IBEulerianSourceFunction(
    const std::string& object_name,
    const int Q_current_idx,
    const int Q_new_idx,
    const int Q_half_idx)
    : CartGridFunction(object_name),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_Q_current_idx(Q_current_idx),
      d_Q_new_idx(Q_new_idx),
      d_Q_half_idx(Q_half_idx)
{
    // intentionally blank
    return;
}// IBEulerianSourceFunction

IBEulerianSourceFunction::~IBEulerianSourceFunction()
{
    // intentionally blank
    return;
}// ~IBEulerianSourceFunction

void
IBEulerianSourceFunction::setTimeInterval(
    const double current_time,
    const double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    return;
}// setTimeInterval

bool
IBEulerianSourceFunction::isTimeDependent() const
{
    return true;
}// isTimeDependent

void
IBEulerianSourceFunction::setDataOnPatch(
    const int data_idx,
    Pointer<Variable<NDIM> > /*var*/,
    Pointer<Patch<NDIM> > patch,
    const double data_time,
    const bool initial_time,
    Pointer<PatchLevel<NDIM> > /*level*/)
{
    Pointer<CellData<NDIM,double> > q_data = patch->getPatchData(data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!q_data.isNull());
#endif
    if (initial_time)
    {
        q_data->fillAll(0.0);
    }
    else if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        if (d_Q_current_idx != -1)
        {
            Pointer<CellData<NDIM,double> > q_current_data = patch->getPatchData(d_Q_current_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!q_current_data.isNull());
#endif
            q_data->copy(*q_current_data);
        }
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        if (d_Q_new_idx != -1)
        {
            Pointer<CellData<NDIM,double> > q_new_data = patch->getPatchData(d_Q_new_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!q_new_data.isNull());
#endif
            q_data->copy(*q_new_data);
        }
    }
    else if (MathUtilities<double>::equalEps(data_time, 0.5*(d_current_time+d_new_time)))
    {
        if (d_Q_half_idx != -1)
        {
            Pointer<CellData<NDIM,double> > q_half_data = patch->getPatchData(d_Q_half_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(!q_half_data.isNull());
#endif
            q_data->copy(*q_half_data);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name << "::setDataOnPatch():\n"
                   << "  data time " << data_time << " is not the current, new, or half time." << std::endl);
    }
    return;
}// setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
