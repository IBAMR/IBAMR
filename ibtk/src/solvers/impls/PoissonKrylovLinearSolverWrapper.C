// Filename: PoissonKrylovLinearSolverWrapper.C
// Created on 13 Aug 2012 by Boyce Griffith
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

#include "PoissonKrylovLinearSolverWrapper.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/LaplaceOperator.h>
#include <ibtk/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

PoissonKrylovLinearSolverWrapper::PoissonKrylovLinearSolverWrapper(
    Pointer<KrylovLinearSolver> krylov_solver)
    : LinearSolver(krylov_solver->getName(), krylov_solver->getHomogeneousBc()),
      KrylovLinearSolver(krylov_solver->getName(), krylov_solver->getHomogeneousBc()),
      PoissonSolver(krylov_solver->getName(), krylov_solver->getHomogeneousBc()),
      d_krylov_solver(krylov_solver)
{
    setSolutionTime(krylov_solver->getSolutionTime());
    setTimeInterval(krylov_solver->getTimeInterval().first, krylov_solver->getTimeInterval().second);
    return;
}// PoissonKrylovLinearSolverWrapper

PoissonKrylovLinearSolverWrapper::~PoissonKrylovLinearSolverWrapper()
{
    // intentionally blank
    return;
}// ~PoissonKrylovLinearSolverWrapper

void
PoissonKrylovLinearSolverWrapper::setPoissonSpecifications(
    const PoissonSpecifications& poisson_spec)
{
    PoissonSolver::setPoissonSpecifications(poisson_spec);
    Pointer<LaplaceOperator> p_operator = getOperator();
    if (!p_operator.isNull()) p_operator->setPoissonSpecifications(d_poisson_spec);
    Pointer<PoissonSolver> p_preconditioner = getPreconditioner();
    if (!p_preconditioner.isNull()) p_preconditioner->setPoissonSpecifications(d_poisson_spec);
    return;
}// setPoissonSpecifications

void
PoissonKrylovLinearSolverWrapper::setPhysicalBcCoef(
    RobinBcCoefStrategy<NDIM>* bc_coef)
{
    PoissonSolver::setPhysicalBcCoef(bc_coef);
    Pointer<LaplaceOperator> p_operator = getOperator();
    if (!p_operator.isNull()) p_operator->setPhysicalBcCoefs(d_bc_coefs);
    Pointer<PoissonSolver> p_preconditioner = getPreconditioner();
    if (!p_preconditioner.isNull()) p_preconditioner->setPhysicalBcCoefs(d_bc_coefs);
    return;
}// setPhysicalBcCoef

void
PoissonKrylovLinearSolverWrapper::setPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
    PoissonSolver::setPhysicalBcCoefs(bc_coefs);
    Pointer<LaplaceOperator> p_operator = getOperator();
    if (!p_operator.isNull()) p_operator->setPhysicalBcCoefs(d_bc_coefs);
    Pointer<PoissonSolver> p_preconditioner = getPreconditioner();
    if (!p_preconditioner.isNull()) p_preconditioner->setPhysicalBcCoefs(d_bc_coefs);
    return;
}// setPhysicalBcCoefs

void
PoissonKrylovLinearSolverWrapper::setHomogeneousBc(
    bool homogeneous_bc)
{
    KrylovLinearSolver::setHomogeneousBc(homogeneous_bc);
    d_krylov_solver->setHomogeneousBc(d_homogeneous_bc);
    return;
}// setHomogeneousBc

bool
PoissonKrylovLinearSolverWrapper::getHomogeneousBc() const
{
    return d_krylov_solver->getHomogeneousBc();
}// getHomogeneousBc

void
PoissonKrylovLinearSolverWrapper::setSolutionTime(
    double solution_time)
{
    KrylovLinearSolver::setSolutionTime(solution_time);
    d_krylov_solver->setSolutionTime(d_solution_time);
    return;
}// setSolutionTime

double
PoissonKrylovLinearSolverWrapper::getSolutionTime() const
{
    return d_krylov_solver->getSolutionTime();
}// getSolutionTime

void
PoissonKrylovLinearSolverWrapper::setTimeInterval(
    double current_time,
    double new_time)
{
    KrylovLinearSolver::setTimeInterval(current_time, new_time);
    d_krylov_solver->setTimeInterval(d_current_time, d_new_time);
    return;
}// setTimeInterval

std::pair<double,double>
PoissonKrylovLinearSolverWrapper::getTimeInterval() const
{
    return d_krylov_solver->getTimeInterval();
}// getTimeInterval

void
PoissonKrylovLinearSolverWrapper::setHierarchyMathOps(
    Pointer<HierarchyMathOps> hier_math_ops)
{
    KrylovLinearSolver::setHierarchyMathOps(hier_math_ops);
    d_krylov_solver->setHierarchyMathOps(d_hier_math_ops);
    return;
}// setHierarchyMathOps

Pointer<HierarchyMathOps>
PoissonKrylovLinearSolverWrapper::getHierarchyMathOps() const
{
    return d_krylov_solver->getHierarchyMathOps();
}// getHierarchyMathOps

bool
PoissonKrylovLinearSolverWrapper::solveSystem(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& b)
{
    return d_krylov_solver->solveSystem(x,b);
}// solveSystem

void
PoissonKrylovLinearSolverWrapper::initializeSolverState(
    const SAMRAIVectorReal<NDIM,double>& x,
    const SAMRAIVectorReal<NDIM,double>& b)
{
    KrylovLinearSolver::initializeSolverState(x,b);
    d_krylov_solver->initializeSolverState(x,b);
    return;
}// initializeSolverState

void
PoissonKrylovLinearSolverWrapper::deallocateSolverState()
{
    KrylovLinearSolver::deallocateSolverState();
    d_krylov_solver->deallocateSolverState();
    return;
}// deallocateSolverState

void
PoissonKrylovLinearSolverWrapper::setNullspace(
    bool contains_constant_vector,
    const std::vector<Pointer<SAMRAIVectorReal<NDIM,double> > >& nullspace_basis_vecs)
{
    KrylovLinearSolver::setNullspace(contains_constant_vector, nullspace_basis_vecs);
    d_krylov_solver->setNullspace(d_nullspace_contains_constant_vector, d_nullspace_basis_vecs);
    return;
}// setNullspace

void
PoissonKrylovLinearSolverWrapper::setOperator(
    Pointer<LinearOperator> A)
{
    d_A = A;
    d_krylov_solver->setOperator(d_A);
    return;
}// setOperator

Pointer<LinearOperator>
PoissonKrylovLinearSolverWrapper::getOperator() const
{
    return d_krylov_solver->getOperator();
}// getOperator

void
PoissonKrylovLinearSolverWrapper::setPreconditioner(
    Pointer<LinearSolver> pc_solver)
{
    d_pc_solver = pc_solver;
    d_krylov_solver->setPreconditioner(d_pc_solver);
    return;
}// setPreconditioner

Pointer<LinearSolver>
PoissonKrylovLinearSolverWrapper::getPreconditioner() const
{
    return d_krylov_solver->getPreconditioner();
}// getPreconditioner

void
PoissonKrylovLinearSolverWrapper::setInitialGuessNonzero(
    bool initial_guess_nonzero)
{
    KrylovLinearSolver::setInitialGuessNonzero(initial_guess_nonzero);
    d_krylov_solver->setInitialGuessNonzero(d_initial_guess_nonzero);
    return;
}// setInitialGuessNonzero

bool
PoissonKrylovLinearSolverWrapper::getInitialGuessNonzero() const
{
    return d_krylov_solver->getInitialGuessNonzero();
}// getInitialGuessNonzero

void
PoissonKrylovLinearSolverWrapper::setMaxIterations(
    int max_iterations)
{
    KrylovLinearSolver::setMaxIterations(max_iterations);
    d_krylov_solver->setMaxIterations(d_max_iterations);
    return;
}// setMaxIterations

int
PoissonKrylovLinearSolverWrapper::getMaxIterations() const
{
    return d_krylov_solver->getMaxIterations();
}// getMaxIterations

void
PoissonKrylovLinearSolverWrapper::setAbsoluteTolerance(
    double abs_residual_tol)
{
    KrylovLinearSolver::setAbsoluteTolerance(abs_residual_tol);
    d_krylov_solver->setAbsoluteTolerance(d_abs_residual_tol);
    return;
}// setAbsoluteTolerance

double
PoissonKrylovLinearSolverWrapper::getAbsoluteTolerance() const
{
    return d_krylov_solver->getAbsoluteTolerance();
}// getAbsoluteTolerance

void
PoissonKrylovLinearSolverWrapper::setRelativeTolerance(
    double rel_residual_tol)
{
    KrylovLinearSolver::setRelativeTolerance(rel_residual_tol);
    d_krylov_solver->setRelativeTolerance(d_rel_residual_tol);
    return;
}// setRelativeTolerance

double
PoissonKrylovLinearSolverWrapper::getRelativeTolerance() const
{
    return d_krylov_solver->getRelativeTolerance();
}// getRelativeTolerance

int
PoissonKrylovLinearSolverWrapper::getNumIterations() const
{
    return d_krylov_solver->getNumIterations();
}// getNumIterations

double
PoissonKrylovLinearSolverWrapper::getResidualNorm() const
{
    return d_krylov_solver->getResidualNorm();
}// getResidualNorm

void
PoissonKrylovLinearSolverWrapper::enableLogging(
    bool enabled)
{
    KrylovLinearSolver::enableLogging(enabled);
    d_krylov_solver->enableLogging(d_enable_logging);
    return;
}// enableLogging

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
