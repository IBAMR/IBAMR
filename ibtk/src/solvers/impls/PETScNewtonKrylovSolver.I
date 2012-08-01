// Filename: PETScNewtonKrylovSolver.I
// Created on 26 Nov 2003 by Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScSAMRAIVectorReal.h>

// SAMRAI INCLUDES
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline const SNES&
PETScNewtonKrylovSolver::getPETScSNES() const
{
    return d_petsc_snes;
}// getPETScSNES

inline void
PETScNewtonKrylovSolver::setHomogeneousBc(
    const bool homogeneous_bc)
{
    NewtonKrylovSolver::setHomogeneousBc(homogeneous_bc);
    if (!d_F.isNull()) d_F->setHomogeneousBc(homogeneous_bc);
    if (!d_krylov_solver.isNull()) d_krylov_solver->setHomogeneousBc(homogeneous_bc);
    return;
}// setHomogeneousBc

inline void
PETScNewtonKrylovSolver::setSolutionTime(
    const double solution_time)
{
    NewtonKrylovSolver::setSolutionTime(solution_time);
    if (!d_F.isNull()) d_F->setSolutionTime(solution_time);
    if (!d_J.isNull()) d_J->setSolutionTime(solution_time);
    if (!d_krylov_solver.isNull()) d_krylov_solver->setSolutionTime(solution_time);
    return;
}// setSolutionTime

inline void
PETScNewtonKrylovSolver::setTimeInterval(
    const double current_time,
    const double new_time)
{
    NewtonKrylovSolver::setTimeInterval(current_time, new_time);
    if (!d_F.isNull()) d_F->setTimeInterval(current_time, new_time);
    if (!d_J.isNull()) d_J->setTimeInterval(current_time, new_time);
    if (!d_krylov_solver.isNull()) d_krylov_solver->setTimeInterval(current_time, new_time);
    return;
}// setTimeInterval

inline void
PETScNewtonKrylovSolver::setOperator(
    SAMRAI::tbox::Pointer<GeneralOperator> F)
{
    d_user_provided_function = true;
    SAMRAI::tbox::Pointer<GeneralOperator> F_old = d_F;
    d_F = F;
    d_F->setHomogeneousBc(d_homogeneous_bc);
    d_F->setSolutionTime(d_solution_time);
    d_F->setTimeInterval(d_current_time, d_new_time);
    if (d_is_initialized && (d_F != F_old) && !d_F.isNull()) d_F->initializeOperatorState(*d_solver_x, *d_solver_b);
    resetSNESFunction();
    return;
}// setOperator

inline SAMRAI::tbox::Pointer<GeneralOperator>
PETScNewtonKrylovSolver::getOperator() const
{
    return d_F;
}// getOperator

inline SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> >
PETScNewtonKrylovSolver::getSolutionVector() const
{
    Vec petsc_x;
    int ierr = SNESGetSolution(d_petsc_snes, &petsc_x); IBTK_CHKERRQ(ierr);
    return PETScSAMRAIVectorReal::getSAMRAIVector(petsc_x);
}// getSolutionVector

inline SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> >
PETScNewtonKrylovSolver::getFunctionVector() const
{
    Vec petsc_f;
    int ierr = SNESGetFunction(d_petsc_snes, &petsc_f, PETSC_NULL, PETSC_NULL); IBTK_CHKERRQ(ierr);
    return PETScSAMRAIVectorReal::getSAMRAIVector(petsc_f);
}// getFunctionVector

inline void
PETScNewtonKrylovSolver::setJacobian(
    SAMRAI::tbox::Pointer<JacobianOperator> J)
{
    d_user_provided_jacobian = true;
    SAMRAI::tbox::Pointer<JacobianOperator> J_old = d_J;
    d_J = J;
    d_J->setHomogeneousBc(true);
    d_J->setSolutionTime(d_solution_time);
    d_J->setTimeInterval(d_current_time, d_new_time);
    if (d_is_initialized && (d_J != J_old) && !d_J.isNull()) d_J->initializeOperatorState(*d_solver_x, *d_solver_b);
    resetSNESJacobian();
    return;
}// setJacobian

inline SAMRAI::tbox::Pointer<JacobianOperator>
PETScNewtonKrylovSolver::getJacobian() const
{
    return d_J;
}// getJacobian

inline SAMRAI::tbox::Pointer<KrylovLinearSolver>
PETScNewtonKrylovSolver::getLinearSolver() const
{
    return d_krylov_solver;
}// getLinearSolver

inline void
PETScNewtonKrylovSolver::setMaxIterations(
    int max_iterations)
{
    d_max_iterations = max_iterations;
    resetSNESOptions();
    return;
}// setMaxIterations

inline int
PETScNewtonKrylovSolver::getMaxIterations() const
{
    return d_max_iterations;
}// getMaxIterations

inline void
PETScNewtonKrylovSolver::setMaxEvaluations(
    int max_evaluations)
{
    d_max_evaluations = max_evaluations;
    resetSNESOptions();
    return;
}// setMaxEvaluations

inline int
PETScNewtonKrylovSolver::getMaxEvaluations() const
{
    return d_max_evaluations;
}// getMaxEvaluations

inline void
PETScNewtonKrylovSolver::setAbsoluteTolerance(
    double abs_residual_tol)
{
    d_abs_residual_tol = abs_residual_tol;
    resetSNESOptions();
    return;
}// setAbsoluteTolerance

inline double
PETScNewtonKrylovSolver::getAbsoluteTolerance() const
{
    return d_abs_residual_tol;
}// getAbsoluteTolerance

inline void
PETScNewtonKrylovSolver::setRelativeTolerance(
    double rel_residual_tol)
{
    d_rel_residual_tol = rel_residual_tol;
    resetSNESOptions();
    return;
}// setRelativeTolerance

inline double
PETScNewtonKrylovSolver::getRelativeTolerance() const
{
    return d_rel_residual_tol;
}// getRelativeTolerance

inline void
PETScNewtonKrylovSolver::setSolutionTolerance(
    double solution_tol)
{
    d_solution_tol = solution_tol;
    resetSNESOptions();
    return;
}// setSolutionTolerance

inline double
PETScNewtonKrylovSolver::getSolutionTolerance() const
{
    return d_solution_tol;
}// getSolutionTolerance

inline int
PETScNewtonKrylovSolver::getNumIterations() const
{
    return d_current_its;
}// getNumIterations

inline int
PETScNewtonKrylovSolver::getNumLinearIterations() const
{
    return d_current_lits;
}// getNumLinearIterations

inline double
PETScNewtonKrylovSolver::getResidualNorm() const
{
    return d_current_residual_norm;
}// getResidualNorm

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
