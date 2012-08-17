// Filename: KrylovLinearSolverWrapper.C
// Created on 17 Aug 2012 by Boyce Griffith
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

#include "KrylovLinearSolverWrapper.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

KrylovLinearSolverWrapper::KrylovLinearSolverWrapper(
    Pointer<KrylovLinearSolver> krylov_solver)
    : LinearSolver(krylov_solver->getName(), krylov_solver->getHomogeneousBc()),
      KrylovLinearSolver(krylov_solver->getName(), krylov_solver->getHomogeneousBc()),
      d_krylov_solver(krylov_solver)
{
    // Configure the wrapper object to have the same state as the wrapped
    // object.

    // First, we use all of the setter methods.
    KrylovLinearSolver::setHomogeneousBc(krylov_solver->getHomogeneousBc());
    KrylovLinearSolver::setSolutionTime(krylov_solver->getSolutionTime());
    KrylovLinearSolver::setTimeInterval(krylov_solver->getTimeInterval().first, krylov_solver->getTimeInterval().second);
    KrylovLinearSolver::setHierarchyMathOps(krylov_solver->getHierarchyMathOps());
    KrylovLinearSolver::setNullspace(krylov_solver->getNullspaceContainsConstantVector(), krylov_solver->getNullspaceBasisVectors());
    KrylovLinearSolver::setOperator(krylov_solver->getOperator());
    KrylovLinearSolver::setPreconditioner(krylov_solver->getPreconditioner());
    KrylovLinearSolver::setInitialGuessNonzero(krylov_solver->getInitialGuessNonzero());
    KrylovLinearSolver::setMaxIterations(krylov_solver->getMaxIterations());
    KrylovLinearSolver::setAbsoluteTolerance(krylov_solver->getAbsoluteTolerance());
    KrylovLinearSolver::setRelativeTolerance(krylov_solver->getRelativeTolerance());
    KrylovLinearSolver::setLoggingEnabled(krylov_solver->getLoggingEnabled());

    // Second, we set those fields for which there are no setter methods.
    d_is_initialized = krylov_solver->getIsInitialized();
    d_current_its = krylov_solver->getNumIterations();
    d_current_residual_norm = krylov_solver->getResidualNorm();
    return;
}// KrylovLinearSolverWrapper

KrylovLinearSolverWrapper::~KrylovLinearSolverWrapper()
{
    // intentionally blank
    return;
}// ~KrylovLinearSolverWrapper

bool
KrylovLinearSolverWrapper::getIsInitialized() const
{
    return d_krylov_solver->getIsInitialized();
}// getIsInitialized

void
KrylovLinearSolverWrapper::setHomogeneousBc(
    bool homogeneous_bc)
{
    KrylovLinearSolver::setHomogeneousBc(homogeneous_bc);
    d_krylov_solver->setHomogeneousBc(d_homogeneous_bc);
    return;
}// setHomogeneousBc

bool
KrylovLinearSolverWrapper::getHomogeneousBc() const
{
    return d_krylov_solver->getHomogeneousBc();
}// getHomogeneousBc

void
KrylovLinearSolverWrapper::setSolutionTime(
    double solution_time)
{
    KrylovLinearSolver::setSolutionTime(solution_time);
    d_krylov_solver->setSolutionTime(d_solution_time);
    return;
}// setSolutionTime

double
KrylovLinearSolverWrapper::getSolutionTime() const
{
    return d_krylov_solver->getSolutionTime();
}// getSolutionTime

void
KrylovLinearSolverWrapper::setTimeInterval(
    double current_time,
    double new_time)
{
    KrylovLinearSolver::setTimeInterval(current_time, new_time);
    d_krylov_solver->setTimeInterval(d_current_time, d_new_time);
    return;
}// setTimeInterval

std::pair<double,double>
KrylovLinearSolverWrapper::getTimeInterval() const
{
    return d_krylov_solver->getTimeInterval();
}// getTimeInterval

void
KrylovLinearSolverWrapper::setHierarchyMathOps(
    Pointer<HierarchyMathOps> hier_math_ops)
{
    KrylovLinearSolver::setHierarchyMathOps(hier_math_ops);
    d_krylov_solver->setHierarchyMathOps(d_hier_math_ops);
    return;
}// setHierarchyMathOps

Pointer<HierarchyMathOps>
KrylovLinearSolverWrapper::getHierarchyMathOps() const
{
    return d_krylov_solver->getHierarchyMathOps();
}// getHierarchyMathOps

bool
KrylovLinearSolverWrapper::solveSystem(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& b)
{
    return d_krylov_solver->solveSystem(x,b);
}// solveSystem

void
KrylovLinearSolverWrapper::initializeSolverState(
    const SAMRAIVectorReal<NDIM,double>& x,
    const SAMRAIVectorReal<NDIM,double>& b)
{
    KrylovLinearSolver::initializeSolverState(x,b);
    d_krylov_solver->initializeSolverState(x,b);
    return;
}// initializeSolverState

void
KrylovLinearSolverWrapper::deallocateSolverState()
{
    KrylovLinearSolver::deallocateSolverState();
    d_krylov_solver->deallocateSolverState();
    return;
}// deallocateSolverState

void
KrylovLinearSolverWrapper::setNullspace(
    bool contains_constant_vector,
    const std::vector<Pointer<SAMRAIVectorReal<NDIM,double> > >& nullspace_basis_vecs)
{
    KrylovLinearSolver::setNullspace(contains_constant_vector, nullspace_basis_vecs);
    d_krylov_solver->setNullspace(d_nullspace_contains_constant_vector, d_nullspace_basis_vecs);
    return;
}// setNullspace

bool
KrylovLinearSolverWrapper::getNullspaceContainsConstantVector() const
{
    return d_krylov_solver->getNullspaceContainsConstantVector();
}// getNullspaceContainsConstantVector

const std::vector<Pointer<SAMRAIVectorReal<NDIM,double> > >&
KrylovLinearSolverWrapper::getNullspaceBasisVectors() const
{
    return d_krylov_solver->getNullspaceBasisVectors();
}// getNullspaceBasisVectors

void
KrylovLinearSolverWrapper::setOperator(
    Pointer<LinearOperator> A)
{
    // We would prefer to use KrylovLinearSolver::setOperator(A) here, but
    // KrylovLinearSolver::setOperator(A) will call
    // A::initializeOperatorState(), as will d_krylov_solver->setOperator(A).
    d_A = A;
    d_krylov_solver->setOperator(d_A);
    return;
}// setOperator

Pointer<LinearOperator>
KrylovLinearSolverWrapper::getOperator() const
{
    return d_krylov_solver->getOperator();
}// getOperator

void
KrylovLinearSolverWrapper::setPreconditioner(
    Pointer<LinearSolver> pc_solver)
{
    // We would prefer to use KrylovLinearSolver::setPreconditioner(pc_solver)
    // here, but KrylovLinearSolver::setPreconditioner(pc_solver) will call
    // pc_solver::initializeSolverState(), as will
    // d_krylov_solver->setPreconditioner(pc_solver).
    d_pc_solver = pc_solver;
    d_krylov_solver->setPreconditioner(d_pc_solver);
    return;
}// setPreconditioner

Pointer<LinearSolver>
KrylovLinearSolverWrapper::getPreconditioner() const
{
    return d_krylov_solver->getPreconditioner();
}// getPreconditioner

void
KrylovLinearSolverWrapper::setInitialGuessNonzero(
    bool initial_guess_nonzero)
{
    KrylovLinearSolver::setInitialGuessNonzero(initial_guess_nonzero);
    d_krylov_solver->setInitialGuessNonzero(d_initial_guess_nonzero);
    return;
}// setInitialGuessNonzero

bool
KrylovLinearSolverWrapper::getInitialGuessNonzero() const
{
    return d_krylov_solver->getInitialGuessNonzero();
}// getInitialGuessNonzero

void
KrylovLinearSolverWrapper::setMaxIterations(
    int max_iterations)
{
    KrylovLinearSolver::setMaxIterations(max_iterations);
    d_krylov_solver->setMaxIterations(d_max_iterations);
    return;
}// setMaxIterations

int
KrylovLinearSolverWrapper::getMaxIterations() const
{
    return d_krylov_solver->getMaxIterations();
}// getMaxIterations

void
KrylovLinearSolverWrapper::setAbsoluteTolerance(
    double abs_residual_tol)
{
    KrylovLinearSolver::setAbsoluteTolerance(abs_residual_tol);
    d_krylov_solver->setAbsoluteTolerance(d_abs_residual_tol);
    return;
}// setAbsoluteTolerance

double
KrylovLinearSolverWrapper::getAbsoluteTolerance() const
{
    return d_krylov_solver->getAbsoluteTolerance();
}// getAbsoluteTolerance

void
KrylovLinearSolverWrapper::setRelativeTolerance(
    double rel_residual_tol)
{
    KrylovLinearSolver::setRelativeTolerance(rel_residual_tol);
    d_krylov_solver->setRelativeTolerance(d_rel_residual_tol);
    return;
}// setRelativeTolerance

double
KrylovLinearSolverWrapper::getRelativeTolerance() const
{
    return d_krylov_solver->getRelativeTolerance();
}// getRelativeTolerance

int
KrylovLinearSolverWrapper::getNumIterations() const
{
    return d_krylov_solver->getNumIterations();
}// getNumIterations

double
KrylovLinearSolverWrapper::getResidualNorm() const
{
    return d_krylov_solver->getResidualNorm();
}// getResidualNorm

void
KrylovLinearSolverWrapper::setLoggingEnabled(
    bool enable_logging)
{
    KrylovLinearSolver::setLoggingEnabled(enable_logging);
    d_krylov_solver->setLoggingEnabled(d_enable_logging);
    return;
}// setLoggingEnabled

bool
KrylovLinearSolverWrapper::getLoggingEnabled() const
{
    return d_krylov_solver->getLoggingEnabled();
}// getLoggingEnabled

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
