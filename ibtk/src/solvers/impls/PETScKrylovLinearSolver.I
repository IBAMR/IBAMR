// Filename: PETScKrylovLinearSolver.I
// Created on 16 Sep 2003 by Boyce Griffith
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

// SAMRAI INCLUDES
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline void
PETScKrylovLinearSolver::setKSPType(
    const std::string& ksp_type)
{
    d_ksp_type = ksp_type;
    resetKSPOptions();
    return;
}// setKSPType

inline KSP
PETScKrylovLinearSolver::getPETScKSP() const
{
    return d_petsc_ksp;
}// getPETScKSP

inline void
PETScKrylovLinearSolver::setTimeInterval(
    const double current_time,
    const double new_time)
{
    if (!d_pc_solver.isNull()) d_pc_solver->setTimeInterval(current_time, new_time);
    return;
}// setTimeInterval

inline void
PETScKrylovLinearSolver::setOperator(
    SAMRAI::tbox::Pointer<LinearOperator> A)
{
    SAMRAI::tbox::Pointer<LinearOperator> A_old = d_A;
    d_A = A;
    d_user_provided_mat = true;
    if (d_is_initialized && (d_A != A_old) && !d_A.isNull())
    {
        d_A->initializeOperatorState(*d_solver_x, *d_solver_b);
    }
    resetKSPOperators();
    return;
}// setOperator

inline SAMRAI::tbox::Pointer<LinearOperator>
PETScKrylovLinearSolver::getOperator() const
{
    return d_A;
}// getOperator

inline void
PETScKrylovLinearSolver::setPreconditioner(
    SAMRAI::tbox::Pointer<LinearSolver> pc_solver)
{
    SAMRAI::tbox::Pointer<LinearSolver> pc_solver_old = d_pc_solver;
    d_pc_solver = pc_solver;
    d_user_provided_pc = true;
    if (d_is_initialized && (d_pc_solver != pc_solver_old) && !d_pc_solver.isNull())
    {
        d_pc_solver->initializeSolverState(*d_solver_b, *d_solver_b);
    }
    resetKSPPC();
    return;
}// setPreconditioner

inline SAMRAI::tbox::Pointer<LinearSolver>
PETScKrylovLinearSolver::getPreconditioner() const
{
    return d_pc_solver;
}// getPreconditioner

inline void
PETScKrylovLinearSolver::setNullspace(
    const bool contains_constant_vector,
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > nullspace_basis_vec)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(contains_constant_vector || !nullspace_basis_vec.isNull());
#endif
    if (!nullspace_basis_vec.isNull())
    {
        setNullspace(contains_constant_vector, std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > >(1, nullspace_basis_vec));
    }
    else
    {
        setNullspace(contains_constant_vector, std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > >());
    }
    return;
}// setNullspace

inline void
PETScKrylovLinearSolver::setNullspace(
    const bool contains_constant_vector,
    const std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > >& nullspace_basis_vecs)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(contains_constant_vector || !nullspace_basis_vecs.empty());
    for (unsigned int k = 0; k < nullspace_basis_vecs.size(); ++k)
    {
        TBOX_ASSERT(!nullspace_basis_vecs[k].isNull());
    }
#endif
    deallocateNullspaceData();
    d_nullsp_contains_constant_vector = contains_constant_vector;
    d_solver_nullsp_basis = nullspace_basis_vecs;
    resetKSPNullspace();
    return;
}// setNullspace

inline void
PETScKrylovLinearSolver::setInitialGuessNonzero(
    bool initial_guess_nonzero)
{
    d_initial_guess_nonzero = initial_guess_nonzero;
    resetKSPOptions();
    return;
}// setInitialGuessNonzero

inline bool
PETScKrylovLinearSolver::getInitialGuessNonzero() const
{
    return d_initial_guess_nonzero;
}// getInitialGuessNonzero

inline void
PETScKrylovLinearSolver::setMaxIterations(
    int max_iterations)
{
    d_max_iterations = max_iterations;
    resetKSPOptions();
    return;
}// setMaxIterations

inline int
PETScKrylovLinearSolver::getMaxIterations() const
{
    return d_max_iterations;
}// getMaxIterations

inline void
PETScKrylovLinearSolver::setAbsoluteTolerance(
    double abs_residual_tol)
{
    d_abs_residual_tol = abs_residual_tol;
    resetKSPOptions();
    return;
}// setAbsoluteTolerance

inline double
PETScKrylovLinearSolver::getAbsoluteTolerance() const
{
    return d_abs_residual_tol;
}// getAbsoluteTolerance

inline void
PETScKrylovLinearSolver::setRelativeTolerance(
    double rel_residual_tol)
{
    d_rel_residual_tol = rel_residual_tol;
    resetKSPOptions();
    return;
}// setRelativeTolerance

inline double
PETScKrylovLinearSolver::getRelativeTolerance() const
{
    return d_rel_residual_tol;
}// getRelativeTolerance

inline int
PETScKrylovLinearSolver::getNumIterations() const
{
    return d_current_its;
}// getNumIterations

inline double
PETScKrylovLinearSolver::getResidualNorm() const
{
    return d_current_residual_norm;
}// getResidualNorm

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
