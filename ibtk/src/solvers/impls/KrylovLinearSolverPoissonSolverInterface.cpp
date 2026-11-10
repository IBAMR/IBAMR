// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/KrylovLinearSolver.h>
#include <ibtk/KrylovLinearSolverPoissonSolverInterface.h>
#include <ibtk/LaplaceOperator.h>
#include <ibtk/PoissonSolver.h>

#include <PoissonSpecifications.h>

#include <vector>

#include <ibtk/namespaces.h> // IWYU pragma: keep

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

void
KrylovLinearSolverPoissonSolverInterface::setPoissonSpecifications(const PoissonSpecifications& poisson_spec)
{
    auto p_this = dynamic_cast<KrylovLinearSolver*>(this);
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
    auto p_this = dynamic_cast<KrylovLinearSolver*>(this);
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
    auto p_this = dynamic_cast<KrylovLinearSolver*>(this);
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

bool
KrylovLinearSolverPoissonSolverInterface::supportsBdryInterpType(const std::string& bdry_interp_type) const
{
    auto p_this = dynamic_cast<const KrylovLinearSolver*>(this);
#if !defined(NDEBUG)
    TBOX_ASSERT(p_this);
#endif
    Pointer<LaplaceOperator> p_operator = p_this->getOperator();
    return p_operator ? p_operator->supportsBdryInterpType(bdry_interp_type) :
                        PoissonSolver::supportsBdryInterpType(bdry_interp_type);
} // supportsBdryInterpType

void
KrylovLinearSolverPoissonSolverInterface::setBdryInterpType(const std::string& bdry_interp_type)
{
    auto p_this = dynamic_cast<KrylovLinearSolver*>(this);
#if !defined(NDEBUG)
    TBOX_ASSERT(p_this);
#endif
    PoissonSolver::setBdryInterpType(bdry_interp_type);
    Pointer<LaplaceOperator> p_operator = p_this->getOperator();
    if (p_operator) p_operator->setBdryInterpType(d_bdry_interp_type);
    return;
} // setBdryInterpType

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
