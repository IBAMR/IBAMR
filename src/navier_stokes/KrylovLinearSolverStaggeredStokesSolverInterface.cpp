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

#include "SAMRAI/solv/PoissonSpecifications.h"
#include "ibamr/KrylovLinearSolverStaggeredStokesSolverInterface.h"
#include "ibamr/StaggeredStokesOperator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/KrylovLinearSolver.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI
{
namespace solv
{

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
}

KrylovLinearSolverStaggeredStokesSolverInterface::~KrylovLinearSolverStaggeredStokesSolverInterface()
{
    // intentionally blank
    return;
}

void KrylovLinearSolverStaggeredStokesSolverInterface::setVelocityPoissonSpecifications(
    const PoissonSpecifications& U_problem_coefs)
{
    auto p_this = CPP_CAST<KrylovLinearSolver*>(this);
    StaggeredStokesSolver::setVelocityPoissonSpecifications(U_problem_coefs);
    auto p_operator = boost::dynamic_pointer_cast<StaggeredStokesOperator>(p_this->getOperator());
    if (p_operator) p_operator->setVelocityPoissonSpecifications(d_U_problem_coefs);
    auto p_preconditioner = boost::dynamic_pointer_cast<StaggeredStokesSolver>(p_this->getPreconditioner());
    if (p_preconditioner) p_preconditioner->setVelocityPoissonSpecifications(d_U_problem_coefs);
    return;
}

void KrylovLinearSolverStaggeredStokesSolverInterface::setPhysicalBcCoefs(
    const std::vector<boost::shared_ptr<RobinBcCoefStrategy>>& U_bc_coefs,
    const boost::shared_ptr<RobinBcCoefStrategy>& P_bc_coef)
{
    auto p_this = CPP_CAST<KrylovLinearSolver*>(this);
    StaggeredStokesSolver::setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
    auto p_operator = boost::dynamic_pointer_cast<StaggeredStokesOperator>(p_this->getOperator());
    if (p_operator) p_operator->setPhysicalBcCoefs(d_U_bc_coefs, d_P_bc_coef);
    auto p_preconditioner = boost::dynamic_pointer_cast<StaggeredStokesSolver>(p_this->getPreconditioner());
    if (p_preconditioner) p_preconditioner->setPhysicalBcCoefs(d_U_bc_coefs, d_P_bc_coef);
    return;
}

void KrylovLinearSolverStaggeredStokesSolverInterface::setPhysicalBoundaryHelper(
    const boost::shared_ptr<StaggeredStokesPhysicalBoundaryHelper>& bc_helper)
{
    auto p_this = CPP_CAST<KrylovLinearSolver*>(this);
    StaggeredStokesSolver::setPhysicalBoundaryHelper(bc_helper);
    auto p_operator = boost::dynamic_pointer_cast<StaggeredStokesOperator>(p_this->getOperator());
    if (p_operator) p_operator->setPhysicalBoundaryHelper(d_bc_helper);
    auto p_preconditioner = boost::dynamic_pointer_cast<StaggeredStokesSolver>(p_this->getPreconditioner());
    if (p_preconditioner) p_preconditioner->setPhysicalBoundaryHelper(d_bc_helper);
    return;
}

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
