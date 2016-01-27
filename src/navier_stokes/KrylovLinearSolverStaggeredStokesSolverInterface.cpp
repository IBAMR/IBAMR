// Filename: KrylovLinearSolverStaggeredStokesSolverInterface.cpp
// Created on 16 Aug 2012 by Boyce Griffith
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

#include <ostream>
#include <vector>

#include "PoissonSpecifications.h"
#include "ibamr/KrylovLinearSolverStaggeredStokesSolverInterface.h"
#include "ibamr/StaggeredStokesOperator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/KrylovLinearSolver.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

KrylovLinearSolverStaggeredStokesSolverInterface::KrylovLinearSolverStaggeredStokesSolverInterface()
{
    // intentionally blank
    return;
} // KrylovLinearSolverStaggeredStokesSolverInterface

KrylovLinearSolverStaggeredStokesSolverInterface::~KrylovLinearSolverStaggeredStokesSolverInterface()
{
    // intentionally blank
    return;
} // ~KrylovLinearSolverStaggeredStokesSolverInterface

void
KrylovLinearSolverStaggeredStokesSolverInterface::setVelocityPoissonSpecifications(
    const PoissonSpecifications& U_problem_coefs)
{
    KrylovLinearSolver* p_this = dynamic_cast<KrylovLinearSolver*>(this);
#if !defined(NDEBUG)
    TBOX_ASSERT(p_this);
#endif
    StaggeredStokesSolver::setVelocityPoissonSpecifications(U_problem_coefs);
    Pointer<StaggeredStokesOperator> p_operator = p_this->getOperator();
    if (p_operator) p_operator->setVelocityPoissonSpecifications(d_U_problem_coefs);
    Pointer<StaggeredStokesSolver> p_preconditioner = p_this->getPreconditioner();
    if (p_preconditioner) p_preconditioner->setVelocityPoissonSpecifications(d_U_problem_coefs);
    return;
} // setVelocityPoissonSpecifications

void
KrylovLinearSolverStaggeredStokesSolverInterface::setComponentsHaveNullspace(const bool has_velocity_nullspace,
                                                                             const bool has_pressure_nullspace)
{
    KrylovLinearSolver* p_this = dynamic_cast<KrylovLinearSolver*>(this);
#if !defined(NDEBUG)
    TBOX_ASSERT(p_this);
#endif
    StaggeredStokesSolver::setComponentsHaveNullspace(has_velocity_nullspace, has_pressure_nullspace);
    Pointer<StaggeredStokesSolver> p_preconditioner = p_this->getPreconditioner();
    if (p_preconditioner) p_preconditioner->setComponentsHaveNullspace(has_velocity_nullspace, has_pressure_nullspace);

    return;
} // setComponentsHaveNullspace

void
KrylovLinearSolverStaggeredStokesSolverInterface::setPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
    RobinBcCoefStrategy<NDIM>* P_bc_coef)
{
    KrylovLinearSolver* p_this = dynamic_cast<KrylovLinearSolver*>(this);
#if !defined(NDEBUG)
    TBOX_ASSERT(p_this);
#endif
    StaggeredStokesSolver::setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
    Pointer<StaggeredStokesOperator> p_operator = p_this->getOperator();
    if (p_operator) p_operator->setPhysicalBcCoefs(d_U_bc_coefs, d_P_bc_coef);
    Pointer<StaggeredStokesSolver> p_preconditioner = p_this->getPreconditioner();
    if (p_preconditioner) p_preconditioner->setPhysicalBcCoefs(d_U_bc_coefs, d_P_bc_coef);
    return;
} // setPhysicalBcCoefs

void
KrylovLinearSolverStaggeredStokesSolverInterface::setPhysicalBoundaryHelper(
    Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper)
{
    KrylovLinearSolver* p_this = dynamic_cast<KrylovLinearSolver*>(this);
#if !defined(NDEBUG)
    TBOX_ASSERT(p_this);
#endif
    StaggeredStokesSolver::setPhysicalBoundaryHelper(bc_helper);
    Pointer<StaggeredStokesOperator> p_operator = p_this->getOperator();
    if (p_operator) p_operator->setPhysicalBoundaryHelper(d_bc_helper);
    Pointer<StaggeredStokesSolver> p_preconditioner = p_this->getPreconditioner();
    if (p_preconditioner) p_preconditioner->setPhysicalBoundaryHelper(d_bc_helper);
    return;
} // setPhysicalBoundaryHelper

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
