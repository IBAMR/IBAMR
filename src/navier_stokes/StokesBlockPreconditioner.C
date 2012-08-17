// Filename: StokesBlockPreconditioner.C
// Created on 16 Aug 2012 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#include "StokesBlockPreconditioner.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// SAMRAI INCLUDES
#include <LocationIndexRobinBcCoefs.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

StokesBlockPreconditioner::StokesBlockPreconditioner(
    const std::string& object_name,
    bool needs_velocity_solver,
    bool needs_pressure_solver,
    bool homogeneous_bc)
    : StokesSolver(object_name, homogeneous_bc),
      d_needs_velocity_solver(needs_velocity_solver),
      d_velocity_solver(),
      d_P_problem_coefs(object_name+"::P_problem_coefs"),
      d_needs_pressure_solver(needs_pressure_solver),
      d_pressure_solver()
{
    // intentionally blank
    return;
}// StokesBlockPreconditioner()

StokesBlockPreconditioner::~StokesBlockPreconditioner()
{
    // intentionally blank
    return;
}// ~StokesBlockPreconditioner()

bool
StokesBlockPreconditioner::needsVelocitySubdomainSolver() const
{
    return d_needs_velocity_solver;
}// needsVelocitySubdomainSolver

void
StokesBlockPreconditioner::setVelocitySubdomainSolver(
    Pointer<PoissonSolver> velocity_solver)
{
    IBAMR_DO_ONCE(
        if (!needsVelocitySubdomainSolver())
        {
            pout << d_object_name << "::setVelocitySubdomainSolver():\n"
                 << "WARNING: implementation does not require velocity subdomain solver\n";
        }
                  );
    d_velocity_solver = velocity_solver;
    return;
}// setVelocitySubdomainSolver

void
StokesBlockPreconditioner::setVelocityPoissonSpecifications(
    const PoissonSpecifications& U_problem_coefs)
{
    StokesSolver::setVelocityPoissonSpecifications(U_problem_coefs);
    if (!d_velocity_solver.isNull()) d_velocity_solver->setPoissonSpecifications(U_problem_coefs);
    return;
}// setVelocityPoissonSpecifications

bool
StokesBlockPreconditioner::needsPressureSubdomainSolver() const
{
    return d_needs_pressure_solver;
}// needsPressureSubdomainSolver

void
StokesBlockPreconditioner::setPressureSubdomainSolver(
    Pointer<PoissonSolver> pressure_solver)
{
    IBAMR_DO_ONCE(
        if (!needsPressureSubdomainSolver())
        {
            pout << d_object_name << "::setPressureSubdomainSolver():\n"
                 << "WARNING: implementation does not require pressure subdomain solver\n";
        }
                  );
    d_pressure_solver = pressure_solver;
    return;
}// setPressureSubdomainSolver

void
StokesBlockPreconditioner::setPressurePoissonSpecifications(
    const PoissonSpecifications& P_problem_coefs)
{
    d_P_problem_coefs = P_problem_coefs;
    if (!d_pressure_solver.isNull()) d_pressure_solver->setPoissonSpecifications(P_problem_coefs);
    return;
}// setPressurePoissonSpecifications

void
StokesBlockPreconditioner::setPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
    RobinBcCoefStrategy<NDIM>* P_bc_coef)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(U_bc_coefs.size() == NDIM);
#endif
    StokesSolver::setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
    if (!d_velocity_solver.isNull()) d_velocity_solver->setPhysicalBcCoefs(d_U_bc_coefs);
    if (!d_pressure_solver.isNull()) d_pressure_solver->setPhysicalBcCoef(d_P_bc_coef);
    return;
}// setPhysicalBcCoefs

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
