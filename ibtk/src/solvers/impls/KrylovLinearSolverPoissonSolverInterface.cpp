// Filename: KrylovLinearSolverPoissonSolverInterface.cpp
// Created on 13 Aug 2012 by Boyce Griffith
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
#include "ibtk/KrylovLinearSolver.h"
#include "ibtk/KrylovLinearSolverPoissonSolverInterface.h"
#include "ibtk/LaplaceOperator.h"
#include "ibtk/PoissonSolver.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
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

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

KrylovLinearSolverPoissonSolverInterface::KrylovLinearSolverPoissonSolverInterface()
{
    // intentionally blank
    return;
} // KrylovLinearSolverPoissonSolverInterface

KrylovLinearSolverPoissonSolverInterface::~KrylovLinearSolverPoissonSolverInterface()
{
    // intentionally blank
    return;
} // ~KrylovLinearSolverPoissonSolverInterface

void
KrylovLinearSolverPoissonSolverInterface::setPoissonSpecifications(const PoissonSpecifications& poisson_spec)
{
    KrylovLinearSolver* p_this = dynamic_cast<KrylovLinearSolver*>(this);
#if !defined(NDEBUG)
    TBOX_ASSERT(p_this);
#endif
    PoissonSolver::setPoissonSpecifications(poisson_spec);
    Pointer<LaplaceOperator> p_operator = p_this->getOperator();
    if (p_operator) p_operator->setPoissonSpecifications(d_poisson_spec);
    Pointer<PoissonSolver> p_preconditioner = p_this->getPreconditioner();
    if (p_preconditioner) p_preconditioner->setPoissonSpecifications(d_poisson_spec);
    return;
} // setPoissonSpecifications

void
KrylovLinearSolverPoissonSolverInterface::setPhysicalBcCoef(RobinBcCoefStrategy<NDIM>* bc_coef)
{
    KrylovLinearSolver* p_this = dynamic_cast<KrylovLinearSolver*>(this);
#if !defined(NDEBUG)
    TBOX_ASSERT(p_this);
#endif
    PoissonSolver::setPhysicalBcCoef(bc_coef);
    Pointer<LaplaceOperator> p_operator = p_this->getOperator();
    if (p_operator) p_operator->setPhysicalBcCoefs(d_bc_coefs);
    Pointer<PoissonSolver> p_preconditioner = p_this->getPreconditioner();
    if (p_preconditioner) p_preconditioner->setPhysicalBcCoefs(d_bc_coefs);
    return;
} // setPhysicalBcCoef

void
KrylovLinearSolverPoissonSolverInterface::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
    KrylovLinearSolver* p_this = dynamic_cast<KrylovLinearSolver*>(this);
#if !defined(NDEBUG)
    TBOX_ASSERT(p_this);
#endif
    PoissonSolver::setPhysicalBcCoefs(bc_coefs);
    Pointer<LaplaceOperator> p_operator = p_this->getOperator();
    if (p_operator) p_operator->setPhysicalBcCoefs(d_bc_coefs);
    Pointer<PoissonSolver> p_preconditioner = p_this->getPreconditioner();
    if (p_preconditioner) p_preconditioner->setPhysicalBcCoefs(d_bc_coefs);
    return;
} // setPhysicalBcCoefs

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
