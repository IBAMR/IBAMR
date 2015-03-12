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
#include "SAMRAI/tbox/Database.h"

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

StaggeredStokesFACPreconditioner::StaggeredStokesFACPreconditioner(
    const std::string& object_name,
    boost::shared_ptr<FACPreconditionerStrategy> fac_strategy,
    boost::shared_ptr<Database> input_db,
    const std::string& default_options_prefix)
    : FACPreconditioner(object_name, fac_strategy, input_db, default_options_prefix)
{
    // intentionally blank
    return;
}

StaggeredStokesFACPreconditioner::~StaggeredStokesFACPreconditioner()
{
    // intentionally blank
    return;
}

void StaggeredStokesFACPreconditioner::setVelocityPoissonSpecifications(const PoissonSpecifications& U_problem_coefs)
{
    StaggeredStokesSolver::setVelocityPoissonSpecifications(U_problem_coefs);
    auto p_fac_strategy = boost::dynamic_pointer_cast<StaggeredStokesFACPreconditionerStrategy>(d_fac_strategy);
    if (p_fac_strategy) p_fac_strategy->setVelocityPoissonSpecifications(U_problem_coefs);
    return;
}

void StaggeredStokesFACPreconditioner::setPhysicalBcCoefs(const std::vector<boost::shared_ptr<RobinBcCoefStrategy>>& U_bc_coefs,
                                                          boost::shared_ptr<RobinBcCoefStrategy> P_bc_coef)
{
    StaggeredStokesSolver::setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
    auto p_fac_strategy = boost::dynamic_pointer_cast<StaggeredStokesFACPreconditionerStrategy>(d_fac_strategy);
    if (p_fac_strategy) p_fac_strategy->setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
    return;
}

void StaggeredStokesFACPreconditioner::setPhysicalBoundaryHelper(
    boost::shared_ptr<StaggeredStokesPhysicalBoundaryHelper> bc_helper)
{
    StaggeredStokesSolver::setPhysicalBoundaryHelper(bc_helper);
    auto p_fac_strategy = boost::dynamic_pointer_cast<StaggeredStokesFACPreconditionerStrategy>(d_fac_strategy);
    if (p_fac_strategy) p_fac_strategy->setPhysicalBoundaryHelper(d_bc_helper);
    return;
}

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
