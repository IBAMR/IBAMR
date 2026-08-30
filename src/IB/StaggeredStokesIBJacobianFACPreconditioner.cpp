// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
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

#include <ibamr/IBImplicitStrategy.h>
#include <ibamr/StaggeredStokesIBJacobianFACPreconditioner.h>
#include <ibamr/StaggeredStokesIBLevelRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/FACPreconditionerStrategy.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesIBJacobianFACPreconditioner::StaggeredStokesIBJacobianFACPreconditioner(
    const std::string& object_name,
    Pointer<IBTK::FACPreconditionerStrategy> fac_strategy,
    Pointer<Database> input_db,
    const std::string& default_options_prefix)
    : FACPreconditioner(object_name, fac_strategy, input_db, default_options_prefix)
{
    // intentionally blank
    return;
} // StaggeredStokesIBJacobianFACPreconditioner

void
StaggeredStokesIBJacobianFACPreconditioner::setVelocityPoissonSpecifications(
    const PoissonSpecifications& U_problem_coefs)
{
    StaggeredStokesSolver::setVelocityPoissonSpecifications(U_problem_coefs);
    Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_strategy = getIBFACPreconditionerStrategy();
    if (fac_strategy) fac_strategy->setVelocityPoissonSpecifications(U_problem_coefs);
    return;
} // setVelocityPoissonSpecifications

void
StaggeredStokesIBJacobianFACPreconditioner::setComponentsHaveNullSpace(const bool has_velocity_nullspace,
                                                                       const bool has_pressure_nullspace)
{
    StaggeredStokesSolver::setComponentsHaveNullSpace(has_velocity_nullspace, has_pressure_nullspace);
    Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_strategy = getIBFACPreconditionerStrategy();
    if (fac_strategy)
    {
        fac_strategy->setComponentsHaveNullSpace(d_has_velocity_nullspace, d_has_pressure_nullspace);
    }
    return;
} // setComponentsHaveNullSpace

void
StaggeredStokesIBJacobianFACPreconditioner::setPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
    RobinBcCoefStrategy<NDIM>* P_bc_coef)
{
    StaggeredStokesSolver::setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
    Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_strategy = getIBFACPreconditionerStrategy();
    if (fac_strategy) fac_strategy->setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
    return;
} // setPhysicalBcCoefs

void
StaggeredStokesIBJacobianFACPreconditioner::setPhysicalBoundaryHelper(
    Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper)
{
    StaggeredStokesSolver::setPhysicalBoundaryHelper(bc_helper);
    Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_strategy = getIBFACPreconditionerStrategy();
    if (fac_strategy) fac_strategy->setPhysicalBoundaryHelper(d_bc_helper);
    return;
} // setPhysicalBoundaryHelper

void
StaggeredStokesIBJacobianFACPreconditioner::setIBTimeSteppingType(const TimeSteppingType time_stepping_type)
{
    Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_strategy = getIBFACPreconditionerStrategy();
    if (fac_strategy) fac_strategy->setIBTimeSteppingType(time_stepping_type);
    return;
} // setIBTimeSteppingType

void
StaggeredStokesIBJacobianFACPreconditioner::setIBForceJacobian(Mat& A_mat)
{
    Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_strategy = getIBFACPreconditionerStrategy();
    if (fac_strategy) fac_strategy->setIBForceJacobian(A_mat);
    return;
} // setIBForceJacobian

void
StaggeredStokesIBJacobianFACPreconditioner::setIBInterpOp(Mat& J_mat)
{
    Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_strategy = getIBFACPreconditionerStrategy();
    if (fac_strategy) fac_strategy->setIBInterpOp(J_mat);
    return;
} // setIBInterpOp

void
StaggeredStokesIBJacobianFACPreconditioner::setIBImplicitStrategy(Pointer<IBImplicitStrategy> ib_implicit_ops)
{
    d_ib_implicit_ops = ib_implicit_ops;
    return;
} // setIBImplicitStrategy

void
StaggeredStokesIBJacobianFACPreconditioner::initializeSolverState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                                                                  const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b)
{
    if (d_ib_implicit_ops)
    {
        d_ib_implicit_ops->setUseFixedLEOperators(true);
        d_ib_implicit_ops->updateFixedLEOperators();
    }
    FACPreconditioner::initializeSolverState(x, b);
    return;
} // initializeSolverState

Pointer<StaggeredStokesIBLevelRelaxationFACOperator>
StaggeredStokesIBJacobianFACPreconditioner::getIBFACPreconditionerStrategy() const
{
    return d_fac_strategy;
} // getIBFACPreconditionerStrategy

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
