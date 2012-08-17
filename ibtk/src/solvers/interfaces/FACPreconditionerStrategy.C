// Filename: FACPreconditionerStrategy.C
// Created on 10 Sep 2010 by Boyce Griffith
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

#include "FACPreconditionerStrategy.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

FACPreconditionerStrategy::FACPreconditionerStrategy(
    const std::string& object_name,
    bool homogeneous_bc)
    : d_object_name(object_name),
      d_homogeneous_bc(homogeneous_bc),
      d_solution_time(std::numeric_limits<double>::quiet_NaN()),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_is_initialized(false)
{
    // intentionally blank
    return;
}// FACPreconditionerStrategy

FACPreconditionerStrategy::~FACPreconditionerStrategy()
{
    // intentionally blank
    return;
}// ~FACPreconditionerStrategy

void
FACPreconditionerStrategy::setFACPreconditioner(
    ConstPointer<FACPreconditioner> preconditioner)
{
    d_preconditioner = preconditioner;
    return;
}// setPreconditioner

const std::string&
FACPreconditionerStrategy::getName() const
{
    return d_object_name;
}// getName

void
FACPreconditionerStrategy::setHomogeneousBc(
    bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

bool
FACPreconditionerStrategy::getHomogeneousBc() const
{
    return d_homogeneous_bc;
}// getHomogeneousBc

void
FACPreconditionerStrategy::setSolutionTime(
    double solution_time)
{
    d_solution_time = solution_time;
    return;
}// setSolutionTime

double
FACPreconditionerStrategy::getSolutionTime() const
{
    return d_solution_time;
}// getSolutionTime

void
FACPreconditionerStrategy::setTimeInterval(
    double current_time,
    double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    return;
}// setTimeInterval

std::pair<double,double>
FACPreconditionerStrategy::getTimeInterval() const
{
    return std::make_pair(d_current_time,d_new_time);
}// getTimeInterval

void
FACPreconditionerStrategy::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& /*solution*/,
    const SAMRAIVectorReal<NDIM,double>& /*rhs*/)
{
    // intentionally blank
    return;
}// initializeOperatorState

void
FACPreconditionerStrategy::deallocateOperatorState()
{
    // intentionally blank
    return;
}// deallocateOperatorState

/////////////////////////////// PROTECTED ////////////////////////////////////

Pointer<SAMRAIVectorReal<NDIM,double> >
FACPreconditionerStrategy::getLevelSAMRAIVectorReal(
    const SAMRAIVectorReal<NDIM,double>& vec,
    int level_num) const
{
    std::ostringstream name_str;
    name_str << vec.getName() << "::level_" << level_num;
    Pointer<SAMRAIVectorReal<NDIM,double> > level_vec = new SAMRAIVectorReal<NDIM,double>(name_str.str(), vec.getPatchHierarchy(), level_num, level_num);
    for (int comp = 0; comp < vec.getNumberOfComponents(); ++comp)
    {
        level_vec->addComponent(vec.getComponentVariable(comp), vec.getComponentDescriptorIndex(comp), vec.getControlVolumeIndex(comp));
    }
    return level_vec;
}// getLevelSAMRAIVectorReal

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
