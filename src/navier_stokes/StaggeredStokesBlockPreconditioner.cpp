// Filename: StaggeredStokesBlockPreconditioner.cpp
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

#include <stddef.h>
#include <ostream>
#include <string>
#include <vector>

#include "SAMRAI/math/HierarchyDataOpsManager.h"
#include "SAMRAI/math/HierarchyDataOpsReal.h"
#include "SAMRAI/hier/IntVector.h"

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/solv/PoissonSpecifications.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "ibamr/StaggeredStokesBlockPreconditioner.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/PoissonSolver.h"
#include "SAMRAI/tbox/PIO.h"

#include "SAMRAI/tbox/Utilities.h"

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

StaggeredStokesBlockPreconditioner::StaggeredStokesBlockPreconditioner(bool needs_velocity_solver,
                                                                       bool needs_pressure_solver)
    : d_needs_velocity_solver(needs_velocity_solver), d_velocity_solver(), d_P_problem_coefs("P_problem_coefs"),
      d_needs_pressure_solver(needs_pressure_solver), d_pressure_solver(), d_hierarchy(NULL), d_coarsest_ln(-1),
      d_finest_ln(-1), d_velocity_data_ops(NULL), d_pressure_data_ops(NULL), d_velocity_wgt_idx(-1),
      d_pressure_wgt_idx(-1), d_hier_math_ops(NULL)
{
    // intentionally blank
    return;
}

StaggeredStokesBlockPreconditioner::~StaggeredStokesBlockPreconditioner()
{
    // intentionally blank
    return;
}

bool StaggeredStokesBlockPreconditioner::needsVelocitySubdomainSolver() const
{
    return d_needs_velocity_solver;
}

void StaggeredStokesBlockPreconditioner::setVelocitySubdomainSolver(
    const boost::shared_ptr<PoissonSolver>& velocity_solver)
{
    IBAMR_DO_ONCE(if (!needsVelocitySubdomainSolver())
                  {
                      pout << d_object_name << "::setVelocitySubdomainSolver():\n"
                           << "WARNING: implementation does not require velocity subdomain solver\n";
                  });
    d_velocity_solver = velocity_solver;
    return;
}

void StaggeredStokesBlockPreconditioner::setVelocityPoissonSpecifications(const PoissonSpecifications& U_problem_coefs)
{
    StaggeredStokesSolver::setVelocityPoissonSpecifications(U_problem_coefs);
    if (d_velocity_solver) d_velocity_solver->setPoissonSpecifications(U_problem_coefs);
    return;
}

bool StaggeredStokesBlockPreconditioner::needsPressureSubdomainSolver() const
{
    return d_needs_pressure_solver;
}

void StaggeredStokesBlockPreconditioner::setPressureSubdomainSolver(
    const boost::shared_ptr<PoissonSolver>& pressure_solver)
{
    IBAMR_DO_ONCE(if (!needsPressureSubdomainSolver())
                  {
                      pout << d_object_name << "::setPressureSubdomainSolver():\n"
                           << "WARNING: implementation does not require pressure subdomain solver\n";
                  });
    d_pressure_solver = pressure_solver;
    return;
}

void StaggeredStokesBlockPreconditioner::setPressurePoissonSpecifications(const PoissonSpecifications& P_problem_coefs)
{
    d_P_problem_coefs = P_problem_coefs;
    if (d_pressure_solver) d_pressure_solver->setPoissonSpecifications(P_problem_coefs);
    return;
}

void StaggeredStokesBlockPreconditioner::setPhysicalBcCoefs(
    const std::vector<boost::shared_ptr<RobinBcCoefStrategy>>& U_bc_coefs,
    const boost::shared_ptr<RobinBcCoefStrategy>& P_bc_coef)
{
    TBOX_ASSERT(U_bc_coefs.size() == NDIM);
    StaggeredStokesSolver::setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
    if (d_velocity_solver) d_velocity_solver->setPhysicalBcCoefs(d_U_bc_coefs);
    if (d_pressure_solver) d_pressure_solver->setPhysicalBcCoef(d_P_bc_coef);
    return;
}

void StaggeredStokesBlockPreconditioner::initializeSolverState(const SAMRAIVectorReal<double>& x,
                                                               const SAMRAIVectorReal<double>& b)
{
    // Get the hierarchy configuration.
    d_hierarchy = x.getPatchHierarchy();
    d_coarsest_ln = x.getCoarsestLevelNumber();
    d_finest_ln = x.getFinestLevelNumber();
    TBOX_ASSERT(d_hierarchy == b.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == b.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == b.getFinestLevelNumber());

    // Setup hierarchy operators.
    auto hier_ops_manager = HierarchyDataOpsManager::getManager();

    d_velocity_data_ops = hier_ops_manager->getOperationsDouble(x.getComponentVariable(0), d_hierarchy, true);
    d_velocity_data_ops->setPatchHierarchy(d_hierarchy);
    d_velocity_data_ops->resetLevels(d_coarsest_ln, d_finest_ln);
    d_velocity_wgt_idx = x.getControlVolumeIndex(0);

    d_pressure_data_ops = hier_ops_manager->getOperationsDouble(x.getComponentVariable(1), d_hierarchy, true);
    d_pressure_data_ops->setPatchHierarchy(d_hierarchy);
    d_pressure_data_ops->resetLevels(d_coarsest_ln, d_finest_ln);
    d_pressure_wgt_idx = x.getControlVolumeIndex(1);

    d_hier_math_ops = boost::make_shared<HierarchyMathOps>(d_object_name + "::HierarchyMathOps", d_hierarchy,
                                                           d_coarsest_ln, d_finest_ln);
    return;
}

void StaggeredStokesBlockPreconditioner::deallocateSolverState()
{
    d_velocity_data_ops.reset();
    d_pressure_data_ops.reset();
    d_velocity_wgt_idx = -1;
    d_pressure_wgt_idx = -1;
    d_hier_math_ops.reset();
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

void StaggeredStokesBlockPreconditioner::correctNullspace(const boost::shared_ptr<SAMRAIVectorReal<double>>& U_vec,
                                                          const boost::shared_ptr<SAMRAIVectorReal<double>>& P_vec)
{
    auto p_velocity_solver = boost::dynamic_pointer_cast<LinearSolver>(d_velocity_solver);
    if (p_velocity_solver)
    {
        const auto& U_nul_vecs = p_velocity_solver->getNullspaceBasisVectors();
        if (!U_nul_vecs.empty())
        {
            for (unsigned int k = 0; k < U_nul_vecs.size(); ++k)
            {
                const double alpha = U_vec->dot(U_nul_vecs[k]) / U_nul_vecs[k]->dot(U_nul_vecs[k]);
                U_vec->axpy(-alpha, U_nul_vecs[k], U_vec);
            }
        }
        TBOX_ASSERT(!p_velocity_solver->getNullspaceContainsConstantVector());
    }

    auto p_pressure_solver = boost::dynamic_pointer_cast<LinearSolver>(d_pressure_solver);
    if (p_pressure_solver)
    {
        if (p_pressure_solver->getNullspaceContainsConstantVector())
        {
            const int P_idx = P_vec->getComponentDescriptorIndex(0);
            const double volume = d_hier_math_ops->getVolumeOfPhysicalDomain();
            auto p_pressure_data_ops =
                boost::dynamic_pointer_cast<HierarchyCellDataOpsReal<double>>(d_pressure_data_ops);
            const double P_mean = (1.0 / volume) * p_pressure_data_ops->integral(P_idx, d_pressure_wgt_idx);
            d_pressure_data_ops->addScalar(P_idx, P_idx, -P_mean);
        }
        TBOX_ASSERT(p_pressure_solver->getNullspaceBasisVectors().empty());
    }
    return;
}

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
