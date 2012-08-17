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
      KrylovLinearSolverWrapper(krylov_solver),
      PoissonSolver(krylov_solver->getName(), krylov_solver->getHomogeneousBc())
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

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
