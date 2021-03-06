// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#include "ibamr/KrylovLinearSolverStaggeredStokesSolverInterface.h"
#include "ibamr/StaggeredStokesOperator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/StaggeredStokesSolver.h"

#include "ibtk/KrylovLinearSolver.h"

#include "PoissonSpecifications.h"
#include "tbox/Pointer.h"

#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

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

void
KrylovLinearSolverStaggeredStokesSolverInterface::setVelocityPoissonSpecifications(
    const PoissonSpecifications& U_problem_coefs)
{
    auto p_this = dynamic_cast<KrylovLinearSolver*>(this);
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
    auto p_this = dynamic_cast<KrylovLinearSolver*>(this);
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
    auto p_this = dynamic_cast<KrylovLinearSolver*>(this);
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
    auto p_this = dynamic_cast<KrylovLinearSolver*>(this);
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
