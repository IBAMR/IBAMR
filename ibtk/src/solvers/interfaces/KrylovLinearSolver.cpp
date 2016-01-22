// Filename: KrylovLinearSolver.cpp
// Created on 09 Sep 2003 by Boyce Griffith
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

#include "ibtk/KrylovLinearSolver.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

KrylovLinearSolver::KrylovLinearSolver() : d_A(NULL), d_pc_solver(NULL), d_x(NULL), d_b(NULL)
{
    // intentionally blank
    return;
} // KrylovLinearSolver()

KrylovLinearSolver::~KrylovLinearSolver()
{
    // intentionally blank
    return;
} // ~KrylovLinearSolver()

void
KrylovLinearSolver::setHierarchyMathOps(Pointer<HierarchyMathOps> hier_math_ops)
{
    LinearSolver::setHierarchyMathOps(hier_math_ops);
    if (d_A) d_A->setHierarchyMathOps(d_hier_math_ops);
    if (d_pc_solver) d_pc_solver->setHierarchyMathOps(d_hier_math_ops);
    return;
} // setHierarchyMathOps

void
KrylovLinearSolver::setHomogeneousBc(bool homogeneous_bc)
{
    LinearSolver::setHomogeneousBc(homogeneous_bc);
    if (d_A) d_A->setHomogeneousBc(homogeneous_bc);
    return;
} // setHomogeneousBc

void
KrylovLinearSolver::setSolutionTime(const double solution_time)
{
    LinearSolver::setSolutionTime(solution_time);
    if (d_A) d_A->setSolutionTime(solution_time);
    if (d_pc_solver) d_pc_solver->setSolutionTime(solution_time);
    return;
} // setSolutionTime

void
KrylovLinearSolver::setTimeInterval(const double current_time, const double new_time)
{
    LinearSolver::setTimeInterval(current_time, new_time);
    if (d_A) d_A->setTimeInterval(current_time, new_time);
    if (d_pc_solver) d_pc_solver->setTimeInterval(current_time, new_time);
    return;
} // setTimeInterval

void
KrylovLinearSolver::setOperator(Pointer<LinearOperator> A)
{
    Pointer<LinearOperator> A_old = d_A;
    d_A = A;
    if (d_A)
    {
        d_A->setHomogeneousBc(d_homogeneous_bc);
        d_A->setSolutionTime(d_solution_time);
        d_A->setTimeInterval(d_current_time, d_new_time);
        if (d_is_initialized && (d_A != A_old) && d_A)
        {
            d_A->initializeOperatorState(*d_x, *d_b);
        }
    }
    return;
} // setOperator

Pointer<LinearOperator>
KrylovLinearSolver::getOperator() const
{
    return d_A;
} // getOperator

void
KrylovLinearSolver::setPreconditioner(Pointer<LinearSolver> pc_solver)
{
    Pointer<LinearSolver> pc_solver_old = d_pc_solver;
    d_pc_solver = pc_solver;
    if (d_pc_solver)
    {
        d_pc_solver->setHomogeneousBc(true);
        d_pc_solver->setSolutionTime(d_solution_time);
        d_pc_solver->setTimeInterval(d_current_time, d_new_time);
        if (d_is_initialized && (d_pc_solver != pc_solver_old) && d_pc_solver)
        {
            d_pc_solver->initializeSolverState(*d_b, *d_b);
        }
    }
    return;
} // setPreconditioner

Pointer<LinearSolver>
KrylovLinearSolver::getPreconditioner() const
{
    return d_pc_solver;
} // getPreconditioner

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
