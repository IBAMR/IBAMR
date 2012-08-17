// Filename: StokesKrylovLinearSolverWrapper.C
// Created on 16 Aug 2012 by Boyce Griffith
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

#include "StokesKrylovLinearSolverWrapper.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/StokesOperator.h>
#include <ibamr/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

StokesKrylovLinearSolverWrapper::StokesKrylovLinearSolverWrapper(
    Pointer<KrylovLinearSolver> krylov_solver)
    : LinearSolver(krylov_solver->getName(), krylov_solver->getHomogeneousBc()),
      KrylovLinearSolverWrapper(krylov_solver),
      StokesSolver(krylov_solver->getName(), krylov_solver->getHomogeneousBc())
{
    // intentionally blank
    return;
}// StokesKrylovLinearSolverWrapper

StokesKrylovLinearSolverWrapper::~StokesKrylovLinearSolverWrapper()
{
    // intentionally blank
    return;
}// ~StokesKrylovLinearSolverWrapper

void
StokesKrylovLinearSolverWrapper::setVelocityPoissonSpecifications(
    const PoissonSpecifications& U_problem_coefs)
{
    StokesSolver::setVelocityPoissonSpecifications(U_problem_coefs);
    Pointer<StokesOperator> p_operator = getOperator();
    if (!p_operator.isNull()) p_operator->setVelocityPoissonSpecifications(d_U_problem_coefs);
    Pointer<StokesSolver> p_preconditioner = getPreconditioner();
    if (!p_preconditioner.isNull()) p_preconditioner->setVelocityPoissonSpecifications(d_U_problem_coefs);
    return;
}// setVelocityPoissonSpecifications

void
StokesKrylovLinearSolverWrapper::setPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
    RobinBcCoefStrategy<NDIM>* P_bc_coef)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(U_bc_coefs.size() == NDIM);
#endif
    StokesSolver::setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
    Pointer<StokesOperator> p_operator = getOperator();
    if (!p_operator.isNull()) p_operator->setPhysicalBcCoefs(d_U_bc_coefs, d_P_bc_coef);
    Pointer<StokesSolver> p_preconditioner = getPreconditioner();
    if (!p_preconditioner.isNull()) p_preconditioner->setPhysicalBcCoefs(d_U_bc_coefs, d_P_bc_coef);
    return;
}// setPhysicalBcCoefs

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
