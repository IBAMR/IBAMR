// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/GeneralSolver.h"
#include "ibtk/HierarchyMathOps.h"

#include "tbox/Pointer.h"

#include <ostream>
#include <string>
#include <utility>

#include "ibtk/namespaces.h" // IWYU pragma: keep

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
