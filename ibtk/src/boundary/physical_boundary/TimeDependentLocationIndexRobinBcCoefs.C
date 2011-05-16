// Filename: TimeDependentLocationIndexRobinBcCoefs.C
// Created on 30 Sep 2006 by Boyce Griffith
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

#include "TimeDependentLocationIndexRobinBcCoefs.h"

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

// SAMRAI INCLUDES
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

TimeDependentLocationIndexRobinBcCoefs::TimeDependentLocationIndexRobinBcCoefs(
    const std::string& object_name,
    Pointer<Database> database)
    : d_object_name(object_name),
      d_bc_coef(new LocationIndexRobinBcCoefs<NDIM>(
                    object_name+"_bc_coef",
                    Pointer<Database>(NULL))),
      d_tau(2*NDIM,std::numeric_limits<double>::quiet_NaN()),
      d_a_coef(2*NDIM,std::numeric_limits<double>::quiet_NaN()),
      d_b_coef(2*NDIM,std::numeric_limits<double>::quiet_NaN()),
      d_initial_g_coef(2*NDIM,std::numeric_limits<double>::quiet_NaN()),
      d_final_g_coef(2*NDIM,std::numeric_limits<double>::quiet_NaN())
{
    if (!database.isNull())
    {
        for (int i = 0; i < 2*NDIM; ++i)
        {
            char buf[256];
            sprintf(buf, "boundary_%d", i);
            std::string name(buf);
            if (database->isString(name))
            {
                d_tau[i] = 1.0;
                d_a_coef[i] = 1.0;
                d_b_coef[i] = 0.0;
                d_initial_g_coef[i] = 0.0;
                d_final_g_coef[i] = 0.0;

                Array<std::string> specs = database->getStringArray(name);
                for (int k = 0; k < specs.size(); k += 2)
                {
                    if (specs[k] == "time_constant")
                    {
                        const double tau = atof(specs[k+1].c_str());
                        setBoundaryTimeConstant(i, tau);
                    }
                    else if (specs[k] == "initial_value")
                    {
                        const double initial_value = atof(specs[k+1].c_str());
                        setBoundaryInitialValue(i, initial_value);
                    }
                    else if (specs[k] == "final_value")
                    {
                        const double final_value = atof(specs[k+1].c_str());
                        setBoundaryFinalValue(i, final_value);
                    }
                    else if (specs[k] == "initial_slope")
                    {
                        const double initial_slope = atof(specs[k+1].c_str());
                        setBoundaryInitialSlope(i, initial_slope);
                    }
                    else if (specs[k] == "final_slope")
                    {
                        const double final_slope = atof(specs[k+1].c_str());
                        setBoundaryFinalSlope(i, final_slope);
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name + "::TimeDependentLocationIndexRobinBcCoefs():\n"
                                   << "  unrecognized key: " << specs[k] << "\n"
                                   << "  valid keys are: time_constant, initial_value, final_value, initial_slope, final_slope" << std::endl);
                    }
                }
            }
        }
    }
    return;
}// TimeDependentLocationIndexRobinBcCoefs

TimeDependentLocationIndexRobinBcCoefs::~TimeDependentLocationIndexRobinBcCoefs()
{
    delete d_bc_coef;
    return;
}// ~TimeDependentLocationIndexRobinBcCoefs

void
TimeDependentLocationIndexRobinBcCoefs::setBoundaryTimeConstant(
    const unsigned int location_index,
    const double tau)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(location_index < 2*NDIM);
    TBOX_ASSERT(tau > 0.0);
#endif
    d_tau[location_index] = tau;
    return;
}// setBoundaryTimeConstant

void
TimeDependentLocationIndexRobinBcCoefs::setBoundaryInitialValue(
    const unsigned int location_index,
    const double value)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(location_index < 2*NDIM);
#endif
    d_a_coef[location_index] = 1.0;
    d_b_coef[location_index] = 0.0;
    d_initial_g_coef[location_index] = value;
    return;
}// setBoundaryInitialValue

void
TimeDependentLocationIndexRobinBcCoefs::setBoundaryFinalValue(
    const unsigned int location_index,
    const double value)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(location_index < 2*NDIM);
#endif
    d_a_coef[location_index] = 1.0;
    d_b_coef[location_index] = 0.0;
    d_final_g_coef[location_index] = value;
    return;
}// setBoundaryFinalValue

void
TimeDependentLocationIndexRobinBcCoefs::setBoundaryInitialSlope(
    const unsigned int location_index,
    const double slope)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(location_index < 2*NDIM);
#endif
    d_a_coef[location_index] = 0.0;
    d_b_coef[location_index] = 1.0;
    d_initial_g_coef[location_index] = slope;
    return;
}// setBoundaryInitialSlope

void
TimeDependentLocationIndexRobinBcCoefs::setBoundaryFinalSlope(
    const unsigned int location_index,
    const double slope)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(location_index < 2*NDIM);
#endif
    d_a_coef[location_index] = 0.0;
    d_b_coef[location_index] = 1.0;
    d_final_g_coef[location_index] = slope;
    return;
}// setBoundaryFinalSlope

void
TimeDependentLocationIndexRobinBcCoefs::setBcCoefs(
    Pointer<ArrayData<NDIM,double> >& acoef_data,
    Pointer<ArrayData<NDIM,double> >& bcoef_data,
    Pointer<ArrayData<NDIM,double> >& gcoef_data,
    const Pointer<Variable<NDIM> >& variable,
    const Patch<NDIM>& patch,
    const BoundaryBox<NDIM>& bdry_box,
    double fill_time) const
{
    for (unsigned int location_index = 0; location_index < 2*NDIM; ++location_index)
    {
        const double tau = d_tau[location_index];
        const double a = d_a_coef[location_index];
        const double b = d_b_coef[location_index];
        const double g_init = d_initial_g_coef[location_index];
        const double g_final = d_final_g_coef[location_index];
        const double g_fill_time = g_init + (g_final-g_init)*(1.0/(1.0+tanh(2.0)))*(tanh(4.0*fill_time/tau-2.0)+tanh(2.0));
        d_bc_coef->setRawCoefficients(location_index, a, b, g_fill_time);
    }
    d_bc_coef->setBcCoefs(acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
    return;
}// setBcCoefs

IntVector<NDIM>
TimeDependentLocationIndexRobinBcCoefs::numberOfExtensionsFillable() const
{
    return d_bc_coef->numberOfExtensionsFillable();
}// numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
