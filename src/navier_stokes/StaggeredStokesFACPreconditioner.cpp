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

#include "ibamr/StaggeredStokesFACPreconditioner.h"
#include "ibamr/StaggeredStokesFACPreconditionerStrategy.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"

#include "ibtk/FACPreconditionerStrategy.h"

#include "tbox/Database.h"

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

StaggeredStokesFACPreconditioner::StaggeredStokesFACPreconditioner(const std::string& object_name,
                                                                   Pointer<FACPreconditionerStrategy> fac_strategy,
                                                                   Pointer<Database> input_db,
                                                                   const std::string& default_options_prefix)
    : FACPreconditioner(object_name, fac_strategy, input_db, default_options_prefix)
{
    // intentionally blank
    return;
} // StaggeredStokesFACPreconditioner

void
StaggeredStokesFACPreconditioner::setVelocityPoissonSpecifications(const PoissonSpecifications& U_problem_coefs)
{
    StaggeredStokesSolver::setVelocityPoissonSpecifications(U_problem_coefs);
    Pointer<StaggeredStokesFACPreconditionerStrategy> p_fac_strategy = d_fac_strategy;
    if (p_fac_strategy) p_fac_strategy->setVelocityPoissonSpecifications(U_problem_coefs);
    return;
} // setVelocityPoissonSpecifications

void
StaggeredStokesFACPreconditioner::setComponentsHaveNullspace(const bool has_velocity_nullspace,
                                                             const bool has_pressure_nullspace)
{
    StaggeredStokesSolver::setComponentsHaveNullspace(has_velocity_nullspace, has_pressure_nullspace);
    Pointer<StaggeredStokesFACPreconditionerStrategy> p_fac_strategy = d_fac_strategy;
    if (p_fac_strategy) p_fac_strategy->setComponentsHaveNullspace(d_has_velocity_nullspace, d_has_pressure_nullspace);

    return;
} // setComponentsHaveNullspace

void
StaggeredStokesFACPreconditioner::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
                                                     RobinBcCoefStrategy<NDIM>* P_bc_coef)
{
    StaggeredStokesSolver::setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
    Pointer<StaggeredStokesFACPreconditionerStrategy> p_fac_strategy = d_fac_strategy;
    if (p_fac_strategy) p_fac_strategy->setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
    return;
} // setPhysicalBcCoefs

void
StaggeredStokesFACPreconditioner::setPhysicalBoundaryHelper(Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper)
{
    StaggeredStokesSolver::setPhysicalBoundaryHelper(bc_helper);
    Pointer<StaggeredStokesFACPreconditionerStrategy> p_fac_strategy = d_fac_strategy;
    if (p_fac_strategy) p_fac_strategy->setPhysicalBoundaryHelper(d_bc_helper);
    return;
} // setPhysicalBoundaryHelper

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
