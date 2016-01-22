// Filename: PoissonFACPreconditioner.cpp
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

#include "ibtk/PoissonFACPreconditioner.h"
#include "ibtk/FACPreconditionerStrategy.h"
#include "ibtk/GeneralSolver.h"
#include "ibtk/PoissonFACPreconditionerStrategy.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
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

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

PoissonFACPreconditioner::PoissonFACPreconditioner(const std::string& object_name,
                                                   Pointer<PoissonFACPreconditionerStrategy> fac_strategy,
                                                   Pointer<Database> input_db,
                                                   const std::string& default_options_prefix)
    : FACPreconditioner(object_name, fac_strategy, input_db, default_options_prefix)
{
    GeneralSolver::init(object_name, /*homogeneous_bc*/ true);
    return;
} // PoissonFACPreconditioner

PoissonFACPreconditioner::~PoissonFACPreconditioner()
{
    // intentionally blank
    return;
} // ~PoissonFACPreconditioner

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
