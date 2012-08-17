// Filename: GeneralSolver.C
// Created on 07 Apr 2012 by Boyce Griffith
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

#include "GeneralSolver.h"

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
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

GeneralSolver::GeneralSolver(
    const std::string& object_name,
    bool homogeneous_bc)
    : d_object_name(object_name),
      d_is_initialized(false),
      d_homogeneous_bc(homogeneous_bc),
      d_solution_time(std::numeric_limits<double>::quiet_NaN()),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_hier_math_ops(NULL),
      d_hier_math_ops_external(false),
      d_enable_logging(false)
{
    // intentionally blank
    return;
}// GeneralSolver()

GeneralSolver::~GeneralSolver()
{
    // intentionally blank
    return;
}// ~GeneralSolver()

const std::string&
GeneralSolver::getName() const
{
    return d_object_name;
}// getName

void
GeneralSolver::setHomogeneousBc(
    bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

bool
GeneralSolver::getHomogeneousBc() const
{
    return d_homogeneous_bc;
}// getHomogeneousBc

void
GeneralSolver::setSolutionTime(
    double solution_time)
{
    d_solution_time = solution_time;
    return;
}// setSolutionTime

double
GeneralSolver::getSolutionTime() const
{
    return d_solution_time;
}// getSolutionTime

void
GeneralSolver::setTimeInterval(
    double current_time,
    double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    return;
}// setTimeInterval

std::pair<double,double>
GeneralSolver::getTimeInterval() const
{
    return std::make_pair(d_current_time,d_new_time);
}// getTimeInterval

void
GeneralSolver::setHierarchyMathOps(
    Pointer<HierarchyMathOps> hier_math_ops)
{
    d_hier_math_ops = hier_math_ops;
    d_hier_math_ops_external = !d_hier_math_ops.isNull();
    return;
}// setHierarchyMathOps

Pointer<HierarchyMathOps>
GeneralSolver::getHierarchyMathOps() const
{
    return d_hier_math_ops;
}// getHierarchyMathOps

void
GeneralSolver::initializeSolverState(
    const SAMRAIVectorReal<NDIM,double>& /*u*/,
    const SAMRAIVectorReal<NDIM,double>& /*r*/)
{
    // intentionally blank
    return;
}// initializeSolverState

void
GeneralSolver::deallocateSolverState()
{
    // intentionally blank
    return;
}// deallocateSolverState

void
GeneralSolver::enableLogging(
    bool enabled)
{
    d_enable_logging = enabled;
    return;
}// enableLogging

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
