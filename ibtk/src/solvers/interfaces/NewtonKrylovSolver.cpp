// Filename: NewtonKrylovSolver.cpp
// Created on 18 Nov 2003 by Boyce Griffith
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

#include "IntVector.h"
#include "SAMRAIVectorReal.h"
#include "ibtk/GeneralOperator.h"
#include "ibtk/GeneralSolver.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/JacobianOperator.h"
#include "ibtk/KrylovLinearSolver.h"
#include "ibtk/NewtonKrylovSolver.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Pointer.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

NewtonKrylovSolver::NewtonKrylovSolver()
    : d_F(NULL),
      d_J(NULL),
      d_krylov_solver(NULL),
      d_x(NULL),
      d_b(NULL),
      d_r(NULL),
      d_max_evaluations(10000),
      d_solution_tol(1.0e-8),
      d_current_linear_iterations(0)
{
    d_max_iterations = 50;
    d_rel_residual_tol = 1.0e-8;
    d_abs_residual_tol = 1.0e-50;
    return;
} // NewtonKrylovSolver()

NewtonKrylovSolver::~NewtonKrylovSolver()
{
    // intentionally blank
    return;
} // ~NewtonKrylovSolver()

void
NewtonKrylovSolver::setHierarchyMathOps(Pointer<HierarchyMathOps> hier_math_ops)
{
    NewtonKrylovSolver::setHierarchyMathOps(hier_math_ops);
    if (d_F) d_F->setHierarchyMathOps(d_hier_math_ops);
    if (d_J) d_J->setHierarchyMathOps(d_hier_math_ops);
    if (d_krylov_solver) d_krylov_solver->setHierarchyMathOps(d_hier_math_ops);
    return;
} // setHierarchyMathOps

void
NewtonKrylovSolver::setHomogeneousBc(const bool homogeneous_bc)
{
    GeneralSolver::setHomogeneousBc(homogeneous_bc);
    if (d_F) d_F->setHomogeneousBc(homogeneous_bc);
    if (d_krylov_solver) d_krylov_solver->setHomogeneousBc(homogeneous_bc);
    return;
} // setHomogeneousBc

void
NewtonKrylovSolver::setSolutionTime(const double solution_time)
{
    GeneralSolver::setSolutionTime(solution_time);
    if (d_F) d_F->setSolutionTime(solution_time);
    if (d_J) d_J->setSolutionTime(solution_time);
    if (d_krylov_solver) d_krylov_solver->setSolutionTime(solution_time);
    return;
} // setSolutionTime

void
NewtonKrylovSolver::setTimeInterval(const double current_time, const double new_time)
{
    GeneralSolver::setTimeInterval(current_time, new_time);
    if (d_F) d_F->setTimeInterval(current_time, new_time);
    if (d_J) d_J->setTimeInterval(current_time, new_time);
    if (d_krylov_solver) d_krylov_solver->setTimeInterval(current_time, new_time);
    return;
} // setTimeInterval

void
NewtonKrylovSolver::setOperator(Pointer<GeneralOperator> F)
{
    Pointer<GeneralOperator> F_old = d_F;
    d_F = F;
    d_F->setHomogeneousBc(d_homogeneous_bc);
    d_F->setSolutionTime(d_solution_time);
    d_F->setTimeInterval(d_current_time, d_new_time);
    if (d_is_initialized && (d_F != F_old) && d_F) d_F->initializeOperatorState(*d_x, *d_b);
    return;
} // setOperator

Pointer<GeneralOperator>
NewtonKrylovSolver::getOperator() const
{
    return d_F;
} // getOperator

void
NewtonKrylovSolver::setJacobian(Pointer<JacobianOperator> J)
{
    Pointer<JacobianOperator> J_old = d_J;
    d_J = J;
    d_J->setHomogeneousBc(true);
    d_J->setSolutionTime(d_solution_time);
    d_J->setTimeInterval(d_current_time, d_new_time);
    if (d_is_initialized && (d_J != J_old) && d_J) d_J->initializeOperatorState(*d_x, *d_b);
    return;
} // setJacobian

Pointer<JacobianOperator>
NewtonKrylovSolver::getJacobian() const
{
    return d_J;
} // getJacobian

Pointer<KrylovLinearSolver>
NewtonKrylovSolver::getLinearSolver() const
{
    return d_krylov_solver;
} // getLinearSolver

void
NewtonKrylovSolver::setMaxEvaluations(int max_evaluations)
{
    d_max_evaluations = max_evaluations;
    return;
} // setMaxEvaluations

int
NewtonKrylovSolver::getMaxEvaluations() const
{
    return d_max_evaluations;
} // getMaxEvaluations

void
NewtonKrylovSolver::setSolutionTolerance(double solution_tol)
{
    d_solution_tol = solution_tol;
    return;
} // setSolutionTolerance

double
NewtonKrylovSolver::getSolutionTolerance() const
{
    return d_solution_tol;
} // getSolutionTolerance

int
NewtonKrylovSolver::getNumLinearIterations() const
{
    return d_current_linear_iterations;
} // getNumLinearIterations

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
