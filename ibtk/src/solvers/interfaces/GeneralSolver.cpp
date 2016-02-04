// Filename: GeneralSolver.cpp
// Created on 07 Apr 2012 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <limits>
#include <ostream>
#include <string>
#include <utility>

#include "ibtk/GeneralSolver.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

GeneralSolver::GeneralSolver()
    : d_object_name("unitialized"),
      d_is_initialized(false),
      d_homogeneous_bc(false),
      d_solution_time(std::numeric_limits<double>::quiet_NaN()),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_rel_residual_tol(0.0),
      d_abs_residual_tol(0.0),
      d_max_iterations(100),
      d_current_iterations(0),
      d_current_residual_norm(std::numeric_limits<double>::quiet_NaN()),
      d_hier_math_ops(NULL),
      d_hier_math_ops_external(false),
      d_enable_logging(false)
{
    // intentionally blank
    return;
} // GeneralSolver()

GeneralSolver::~GeneralSolver()
{
    // intentionally blank
    return;
} // ~GeneralSolver()

const std::string&
GeneralSolver::getName() const
{
    return d_object_name;
} // getName

bool
GeneralSolver::getIsInitialized() const
{
    return d_is_initialized;
} // getIsInitialized

void
GeneralSolver::setHomogeneousBc(bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
} // setHomogeneousBc

bool
GeneralSolver::getHomogeneousBc() const
{
    return d_homogeneous_bc;
} // getHomogeneousBc

void
GeneralSolver::setSolutionTime(double solution_time)
{
    d_solution_time = solution_time;
    return;
} // setSolutionTime

double
GeneralSolver::getSolutionTime() const
{
    return d_solution_time;
} // getSolutionTime

void
GeneralSolver::setTimeInterval(double current_time, double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    return;
} // setTimeInterval

std::pair<double, double>
GeneralSolver::getTimeInterval() const
{
    return std::make_pair(d_current_time, d_new_time);
} // getTimeInterval

double
GeneralSolver::getDt() const
{
    return d_new_time - d_current_time;
} // getDt

void
GeneralSolver::setHierarchyMathOps(Pointer<HierarchyMathOps> hier_math_ops)
{
    d_hier_math_ops = hier_math_ops;
    d_hier_math_ops_external = d_hier_math_ops;
    return;
} // setHierarchyMathOps

Pointer<HierarchyMathOps>
GeneralSolver::getHierarchyMathOps() const
{
    return d_hier_math_ops;
} // getHierarchyMathOps

void
GeneralSolver::initializeSolverState(const SAMRAIVectorReal<NDIM, double>& /*u*/,
                                     const SAMRAIVectorReal<NDIM, double>& /*r*/)
{
    d_is_initialized = true;
    return;
} // initializeSolverState

void
GeneralSolver::deallocateSolverState()
{
    d_is_initialized = false;
    return;
} // deallocateSolverState

void
GeneralSolver::setMaxIterations(int max_iterations)
{
    d_max_iterations = max_iterations;
    return;
} // setMaxIterations

int
GeneralSolver::getMaxIterations() const
{
    return d_max_iterations;
} // getMaxIterations

void
GeneralSolver::setAbsoluteTolerance(double abs_residual_tol)
{
    d_abs_residual_tol = abs_residual_tol;
    return;
} // setAbsoluteTolerance

double
GeneralSolver::getAbsoluteTolerance() const
{
    return d_abs_residual_tol;
} // getAbsoluteTolerance

void
GeneralSolver::setRelativeTolerance(double rel_residual_tol)
{
    d_rel_residual_tol = rel_residual_tol;
    return;
} // setRelativeTolerance

double
GeneralSolver::getRelativeTolerance() const
{
    return d_rel_residual_tol;
} // getRelativeTolerance

int
GeneralSolver::getNumIterations() const
{
    return d_current_iterations;
} // getNumIterations

double
GeneralSolver::getResidualNorm() const
{
    return d_current_residual_norm;
} // getResidualNorm

void
GeneralSolver::setLoggingEnabled(bool enable_logging)
{
    d_enable_logging = enable_logging;
    return;
} // setLoggingEnabled

bool
GeneralSolver::getLoggingEnabled() const
{
    return d_enable_logging;
} // getLoggingEnabled

void
GeneralSolver::printClassData(std::ostream& stream)
{
    stream << "\n"
           << "object_name = " << d_object_name << "\n"
           << "is_initialized = " << d_is_initialized << "\n"
           << "homogeneous_bc = " << d_homogeneous_bc << "\n"
           << "solution_time = " << d_solution_time << "\n"
           << "current_time = " << d_current_time << "\n"
           << "new_time = " << d_new_time << "\n"
           << "rel_residual_tol = " << d_rel_residual_tol << "\n"
           << "abs_residual_tol = " << d_abs_residual_tol << "\n"
           << "max_iterations = " << d_max_iterations << "\n"
           << "current_iterations = " << d_current_iterations << "\n"
           << "current_residual_norm = " << d_current_residual_norm << "\n"
           << "hier_math_ops = " << d_hier_math_ops.getPointer() << "\n"
           << "hier_math_ops_external = " << d_hier_math_ops_external << "\n"
           << "enable_logging = " << d_enable_logging << "\n";
    return;
} // printClassData

/////////////////////////////// PROTECTED ////////////////////////////////////

void
GeneralSolver::init(const std::string& object_name, const bool homogeneous_bc)
{
    d_object_name = object_name;
    d_homogeneous_bc = homogeneous_bc;
    initSpecialized(object_name, homogeneous_bc);
    return;
} // init

void
GeneralSolver::initSpecialized(const std::string& /*object_name*/, const bool /*homogeneous_bc*/)
{
    // intentionally blank
    return;
} // init

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
