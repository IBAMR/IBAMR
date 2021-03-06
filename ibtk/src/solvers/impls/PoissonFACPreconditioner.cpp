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

#include "ibtk/FACPreconditioner.h"
#include "ibtk/FACPreconditionerStrategy.h"
#include "ibtk/GeneralSolver.h"
#include "ibtk/PoissonFACPreconditioner.h"
#include "ibtk/PoissonFACPreconditionerStrategy.h"
#include "ibtk/PoissonSolver.h"

#include "PoissonSpecifications.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <string>
#include <utility>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

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

PoissonFACPreconditioner::PoissonFACPreconditioner(const std::string& object_name,
                                                   Pointer<PoissonFACPreconditionerStrategy> fac_strategy,
                                                   Pointer<Database> input_db,
                                                   std::string default_options_prefix)
    : FACPreconditioner(object_name, fac_strategy, input_db, std::move(default_options_prefix))
{
    GeneralSolver::init(object_name, /*homogeneous_bc*/ true);
    return;
} // PoissonFACPreconditioner

void
PoissonFACPreconditioner::setPoissonSpecifications(const PoissonSpecifications& poisson_spec)
{
    PoissonSolver::setPoissonSpecifications(poisson_spec);
    Pointer<PoissonFACPreconditionerStrategy> p_fac_strategy = d_fac_strategy;
    if (p_fac_strategy) p_fac_strategy->setPoissonSpecifications(d_poisson_spec);
    return;
} // setPoissonSpecifications

void
PoissonFACPreconditioner::setPhysicalBcCoef(RobinBcCoefStrategy<NDIM>* bc_coef)
{
    PoissonSolver::setPhysicalBcCoef(bc_coef);
    Pointer<PoissonFACPreconditionerStrategy> p_fac_strategy = d_fac_strategy;
    if (p_fac_strategy) p_fac_strategy->setPhysicalBcCoefs(d_bc_coefs);
    return;
} // setPhysicalBcCoef

void
PoissonFACPreconditioner::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
    PoissonSolver::setPhysicalBcCoefs(bc_coefs);
    Pointer<PoissonFACPreconditionerStrategy> p_fac_strategy = d_fac_strategy;
    if (p_fac_strategy) p_fac_strategy->setPhysicalBcCoefs(d_bc_coefs);
    return;
} // setPhysicalBcCoefs

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
