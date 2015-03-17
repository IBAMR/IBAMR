// Filename: IBStrategy.cpp
// Created on 21 Sep 2011 by Boyce Griffith
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
#include <algorithm>
#include <ostream>
#include <string>
#include <vector>

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/math/HierarchyDataOpsReal.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/IBStrategy.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyMathOps.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Database.h"

#include "SAMRAI/tbox/Utilities.h"

namespace IBAMR
{
class INSHierarchyIntegrator;
} // namespace IBAMR

namespace SAMRAI
{
namespace mesh
{

class ChopAndPackLoadBalancer;
} // namespace mesh
namespace xfer
{

class CoarsenPatchStrategy;

class CoarsenSchedule;

class RefinePatchStrategy;

class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBStrategy::IBStrategy() : d_ib_solver(NULL), d_use_fixed_coupling_ops(false)
{
    // intentionally blank
    return;
}

IBStrategy::~IBStrategy()
{
    // intentionally blank
    return;
}

void IBStrategy::registerIBHierarchyIntegrator(IBHierarchyIntegrator* ib_solver)
{
    d_ib_solver = ib_solver;
    return;
}

void IBStrategy::registerEulerianVariables()
{
    // intentionally blank
    return;
}

void IBStrategy::registerEulerianCommunicationAlgorithms()
{
    // intentionally blank
    return;
}

void IBStrategy::setupTagBuffer(std::vector<int>& tag_buffer, const boost::shared_ptr<PatchHierarchy>& hierarchy) const
{
    const int finest_hier_ln = hierarchy->getMaxNumberOfLevels() - 1;
    const auto tsize = tag_buffer.size();
    tag_buffer.resize(finest_hier_ln);
    for (auto i = tsize; i < finest_hier_ln; ++i) tag_buffer[i] = 0;
    const int gcw = getMinimumGhostCellWidth().max();
    for (int i = 0; i < tag_buffer.size(); ++i)
    {
        tag_buffer[i] = std::max(tag_buffer[i], gcw);
    }
    return;
}

void IBStrategy::preprocessIntegrateData(double /*current_time*/, double /*new_time*/, int /*num_cycles*/)
{
    // intentionally blank
    return;
}

void IBStrategy::postprocessIntegrateData(double /*current_time*/, double /*new_time*/, int /*num_cycles*/)
{
    // intentionally blank
    return;
}

void IBStrategy::setUseFixedLEOperators(bool use_fixed_coupling_ops)
{
    d_use_fixed_coupling_ops = use_fixed_coupling_ops;
    return;
}

void IBStrategy::updateFixedLEOperators()
{
    TBOX_ERROR("IBStrategy::updateFixedLEOperators(): unimplemented\n");
    return;
}

bool IBStrategy::hasFluidSources() const
{
    return false;
}

void IBStrategy::computeLagrangianFluidSource(double /*data_time*/)
{
    // intentionally blank
    return;
}

void IBStrategy::spreadFluidSource(int /*q_data_idx*/,
                                   const std::vector<boost::shared_ptr<RefineSchedule> >& /*q_prolongation_scheds*/,
                                   double /*data_time*/)
{
    // intentionally blank
    return;
}

void IBStrategy::interpolatePressure(int /*p_data_idx*/,
                                     const std::vector<boost::shared_ptr<CoarsenSchedule> >& /*p_synch_scheds*/,
                                     const std::vector<boost::shared_ptr<RefineSchedule> >& /*p_ghost_fill_scheds*/,
                                     double /*data_time*/)
{
    // intentionally blank
    return;
}

void IBStrategy::preprocessSolveFluidEquations(double /*current_time*/, double /*new_time*/, int /*cycle_num*/)
{
    // intentionally blank
    return;
}

void IBStrategy::postprocessSolveFluidEquations(double /*current_time*/, double /*new_time*/, int /*cycle_num*/)
{
    // intentionally blank
    return;
}

void IBStrategy::postprocessData()
{
    // intentionally blank
    return;
}

void IBStrategy::initializePatchHierarchy(
    const boost::shared_ptr<PatchHierarchy>& /*hierarchy*/,
    const boost::shared_ptr<GriddingAlgorithm>& /*gridding_alg*/,
    int /*u_data_idx*/,
    const std::vector<boost::shared_ptr<CoarsenSchedule> >& /*u_synch_scheds*/,
    const std::vector<boost::shared_ptr<RefineSchedule> >& /*u_ghost_fill_scheds*/,
    int /*integrator_step*/,
    double /*init_data_time*/,
    bool /*initial_time*/)
{
    // intentionally blank
    return;
}

void IBStrategy::registerLoadBalancer(const boost::shared_ptr<ChopAndPackLoadBalancer>& /*load_balancer*/,
                                      int /*workload_data_idx*/)
{
    // intentionally blank
    return;
}

void IBStrategy::updateWorkloadEstimates(const boost::shared_ptr<PatchHierarchy>& /*hierarchy*/,
                                         int /*workload_data_idx*/)
{
    // intentionally blank
    return;
}

void IBStrategy::beginDataRedistribution(const boost::shared_ptr<PatchHierarchy>& /*hierarchy*/,
                                         const boost::shared_ptr<GriddingAlgorithm>& /*gridding_alg*/)
{
    // intentionally blank
    return;
}

void IBStrategy::endDataRedistribution(const boost::shared_ptr<PatchHierarchy>& /*hierarchy*/,
                                       const boost::shared_ptr<GriddingAlgorithm>& /*gridding_alg*/)
{
    // intentionally blank
    return;
}

void IBStrategy::initializeLevelData(const boost::shared_ptr<PatchHierarchy>& /*hierarchy*/,
                                     int /*level_number*/,
                                     double /*init_data_time*/,
                                     bool /*can_be_refined*/,
                                     bool /*initial_time*/,
                                     const boost::shared_ptr<PatchLevel>& /*old_level*/,
                                     bool /*allocate_data*/)
{
    // intentionally blank
    return;
}

void IBStrategy::resetHierarchyConfiguration(const boost::shared_ptr<PatchHierarchy>& /*hierarchy*/,
                                             int /*coarsest_level*/,
                                             int /*finest_level*/)
{
    // intentionally blank
    return;
}

void IBStrategy::applyGradientDetector(const boost::shared_ptr<PatchHierarchy>& /*hierarchy*/,
                                       int /*level_number*/,
                                       double /*error_data_time*/,
                                       int /*tag_index*/,
                                       bool /*initial_time*/,
                                       bool /*uses_richardson_extrapolation_too*/)
{
    // intentionally blank
    return;
}

void IBStrategy::putToRestart(const boost::shared_ptr<Database>& /*db*/) const
{
    // intentionally blank
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

INSHierarchyIntegrator* IBStrategy::getINSHierarchyIntegrator() const
{
    TBOX_ASSERT(d_ib_solver);
    return d_ib_solver->d_ins_hier_integrator.get();
}

boost::shared_ptr<HierarchyDataOpsReal<double> > IBStrategy::getVelocityHierarchyDataOps() const
{
    TBOX_ASSERT(d_ib_solver);
    return d_ib_solver->d_hier_velocity_data_ops;
}

boost::shared_ptr<HierarchyDataOpsReal<double> > IBStrategy::getPressureHierarchyDataOps() const
{
    TBOX_ASSERT(d_ib_solver);
    return d_ib_solver->d_hier_pressure_data_ops;
}

boost::shared_ptr<HierarchyMathOps> IBStrategy::getHierarchyMathOps() const
{
    TBOX_ASSERT(d_ib_solver);
    return d_ib_solver->d_hier_math_ops;
}

void IBStrategy::registerVariable(int& current_idx,
                                  int& new_idx,
                                  int& scratch_idx,
                                  const boost::shared_ptr<Variable>& variable,
                                  const IntVector& scratch_ghosts,
                                  const std::string& coarsen_name,
                                  const std::string& refine_name,
                                  const boost::shared_ptr<CartGridFunction>& init_fcn)
{
    TBOX_ASSERT(d_ib_solver);
    d_ib_solver->registerVariable(current_idx, new_idx, scratch_idx, variable, scratch_ghosts, coarsen_name,
                                  refine_name, init_fcn);
    return;
}

void IBStrategy::registerVariable(int& idx,
                                  const boost::shared_ptr<Variable>& variable,
                                  const IntVector& ghosts,
                                  const boost::shared_ptr<VariableContext>& ctx)
{
    TBOX_ASSERT(d_ib_solver);
    d_ib_solver->registerVariable(idx, variable, ghosts, ctx);
    return;
}

void IBStrategy::registerGhostfillRefineAlgorithm(const std::string& name,
                                                  const boost::shared_ptr<RefineAlgorithm>& ghostfill_alg,
                                                  RefinePatchStrategy* ghostfill_patch_strategy)
{
    d_ib_solver->registerGhostfillRefineAlgorithm(name, ghostfill_alg, ghostfill_patch_strategy);
}

void IBStrategy::registerProlongRefineAlgorithm(const std::string& name,
                                                const boost::shared_ptr<RefineAlgorithm>& prolong_alg,
                                                RefinePatchStrategy* prolong_patch_strategy)
{
    d_ib_solver->registerProlongRefineAlgorithm(name, prolong_alg, prolong_patch_strategy);
}

void IBStrategy::registerCoarsenAlgorithm(const std::string& name,
                                          const boost::shared_ptr<CoarsenAlgorithm>& coarsen_alg,
                                          CoarsenPatchStrategy* coarsen_patch_strategy)
{
    d_ib_solver->registerCoarsenAlgorithm(name, coarsen_alg, coarsen_patch_strategy);
}

boost::shared_ptr<RefineAlgorithm> IBStrategy::getGhostfillRefineAlgorithm(const std::string& name) const
{
    return d_ib_solver->getGhostfillRefineAlgorithm(name);
}

boost::shared_ptr<RefineAlgorithm> IBStrategy::getProlongRefineAlgorithm(const std::string& name) const
{
    return d_ib_solver->getProlongRefineAlgorithm(name);
}

boost::shared_ptr<CoarsenAlgorithm> IBStrategy::getCoarsenAlgorithm(const std::string& name) const
{
    return d_ib_solver->getCoarsenAlgorithm(name);
}

const std::vector<boost::shared_ptr<RefineSchedule> >&
IBStrategy::getGhostfillRefineSchedules(const std::string& name) const
{
    return d_ib_solver->getGhostfillRefineSchedules(name);
}

const std::vector<boost::shared_ptr<RefineSchedule> >&
IBStrategy::getProlongRefineSchedules(const std::string& name) const
{
    return d_ib_solver->getProlongRefineSchedules(name);
}

const std::vector<boost::shared_ptr<CoarsenSchedule> >& IBStrategy::getCoarsenSchedules(const std::string& name) const
{
    return d_ib_solver->getCoarsenSchedules(name);
}

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
