// Filename: NewtonKrylovSolver.C
// Created on 18 Nov 2003 by Boyce Griffith
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

#include "NewtonKrylovSolver.h"

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

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

NewtonKrylovSolver::NewtonKrylovSolver(
    const std::string& object_name,
    bool homogeneous_bc)
    : GeneralSolver(object_name, homogeneous_bc),
      d_F(NULL),
      d_J(NULL),
      d_krylov_solver(NULL),
      d_x(NULL),
      d_b(NULL),
      d_r(NULL),
      d_rel_residual_tol(1.0e-8),
      d_abs_residual_tol(1.0e-50),
      d_solution_tol(1.0e-8),
      d_max_iterations(50),
      d_max_evaluations(10000),
      d_current_iterations(0),
      d_current_linear_iterations(0),
      d_current_residual_norm(std::numeric_limits<double>::quiet_NaN())
{
    // intentionally blank
    return;
}// NewtonKrylovSolver()

NewtonKrylovSolver::~NewtonKrylovSolver()
{
    // intentionally blank
    return;
}// ~NewtonKrylovSolver()

void
NewtonKrylovSolver::setHierarchyMathOps(
    Pointer<HierarchyMathOps> hier_math_ops)
{
    NewtonKrylovSolver::setHierarchyMathOps(hier_math_ops);
    if (d_F) d_F->setHierarchyMathOps(d_hier_math_ops);
    if (d_J) d_J->setHierarchyMathOps(d_hier_math_ops);
    if (d_krylov_solver) d_krylov_solver->setHierarchyMathOps(d_hier_math_ops);
    return;
}// setHierarchyMathOps

void
NewtonKrylovSolver::setOperator(
    Pointer<GeneralOperator> F)
{
    Pointer<GeneralOperator> F_old = d_F;
    d_F = F;
    d_F->setHomogeneousBc(d_homogeneous_bc);
    d_F->setSolutionTime(d_solution_time);
    d_F->setTimeInterval(d_current_time, d_new_time);
    if (d_is_initialized && (d_F != F_old) && d_F) d_F->initializeOperatorState(*d_x, *d_b);
    return;
}// setOperator

Pointer<GeneralOperator>
NewtonKrylovSolver::getOperator() const
{
    return d_F;
}// getOperator

void
NewtonKrylovSolver::setJacobian(
    Pointer<JacobianOperator> J)
{
    Pointer<JacobianOperator> J_old = d_J;
    d_J = J;
    d_J->setHomogeneousBc(true);
    d_J->setSolutionTime(d_solution_time);
    d_J->setTimeInterval(d_current_time, d_new_time);
    if (d_is_initialized && (d_J != J_old) && d_J) d_J->initializeOperatorState(*d_x, *d_b);
    return;
}// setJacobian

Pointer<JacobianOperator>
NewtonKrylovSolver::getJacobian() const
{
    return d_J;
}// getJacobian

Pointer<KrylovLinearSolver>
NewtonKrylovSolver::getLinearSolver() const
{
    return d_krylov_solver;
}// getLinearSolver

void
NewtonKrylovSolver::setMaxIterations(
    int max_iterations)
{
    d_max_iterations = max_iterations;
    return;
}// setMaxIterations

int
NewtonKrylovSolver::getMaxIterations() const
{
    return d_max_iterations;
}// getMaxIterations

void
NewtonKrylovSolver::setMaxEvaluations(
    int max_evaluations)
{
    d_max_evaluations = max_evaluations;
    return;
}// setMaxEvaluations

int
NewtonKrylovSolver::getMaxEvaluations() const
{
    return d_max_evaluations;
}// getMaxEvaluations

void
NewtonKrylovSolver::setAbsoluteTolerance(
    double abs_residual_tol)
{
    d_abs_residual_tol = abs_residual_tol;
    return;
}// setAbsoluteTolerance

double
NewtonKrylovSolver::getAbsoluteTolerance() const
{
    return d_abs_residual_tol;
}// getAbsoluteTolerance

void
NewtonKrylovSolver::setRelativeTolerance(
    double rel_residual_tol)
{
    d_rel_residual_tol = rel_residual_tol;
    return;
}// setRelativeTolerance

double
NewtonKrylovSolver::getRelativeTolerance() const
{
    return d_rel_residual_tol;
}// getRelativeTolerance

void
NewtonKrylovSolver::setSolutionTolerance(
    double solution_tol)
{
    d_solution_tol = solution_tol;
    return;
}// setSolutionTolerance

double
NewtonKrylovSolver::getSolutionTolerance() const
{
    return d_solution_tol;
}// getSolutionTolerance

int
NewtonKrylovSolver::getNumIterations() const
{
    return d_current_iterations;
}// getNumIterations

int
NewtonKrylovSolver::getNumLinearIterations() const
{
    return d_current_linear_iterations;
}// getNumLinearIterations

double
NewtonKrylovSolver::getResidualNorm() const
{
    return d_current_residual_norm;
}// getResidualNorm

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
