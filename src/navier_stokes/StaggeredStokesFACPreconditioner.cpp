// Filename: StaggeredStokesFACPreconditioner.cpp
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

#include "ibamr/StaggeredStokesFACPreconditioner.h"
#include "ibamr/StaggeredStokesFACPreconditionerStrategy.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/FACPreconditionerStrategy.h"
#include "tbox/Database.h"

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

StaggeredStokesFACPreconditioner::~StaggeredStokesFACPreconditioner()
{
    // intentionally blank
    return;
} // ~StaggeredStokesFACPreconditioner

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
