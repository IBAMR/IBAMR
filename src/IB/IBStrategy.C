// Filename: IBStrategy.C
// Created on 21 Sep 2011 by Boyce Griffith
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

#include "IBStrategy.h"

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
#include <ibamr/IBHierarchyIntegrator.h>
#include <ibamr/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBStrategy::IBStrategy()
    : d_ib_solver(NULL),
      d_use_fixed_coupling_ops(false)
{
    // intentionally blank
    return;
}// IBStrategy

IBStrategy::~IBStrategy()
{
    // intentionally blank
    return;
}// ~IBStrategy

void
IBStrategy::registerIBHierarchyIntegrator(
    IBHierarchyIntegrator* ib_solver)
{
    d_ib_solver = ib_solver;
    return;
}// registerIBHierarchyIntegrator

void
IBStrategy::registerEulerianVariables()
{
    // intentionally blank
    return;
}// registerEulerianVariables

void
IBStrategy::registerEulerianCommunicationAlgorithms()
{
    // intentionally blank
    return;
}// registerEulerianCommunicationAlgorithms

void
IBStrategy::preprocessIntegrateData(
    double /*current_time*/,
    double /*new_time*/,
    int /*num_cycles*/)
{
    // intentionally blank
    return;
}// preprocessIntegrateData

void
IBStrategy::postprocessIntegrateData(
    double /*current_time*/,
    double /*new_time*/,
    int /*num_cycles*/)
{
    // intentionally blank
    return;
}// postprocessIntegrateData

void
IBStrategy::setUseFixedLEOperators(
    bool use_fixed_coupling_ops)
{
    d_use_fixed_coupling_ops = use_fixed_coupling_ops;
    return;
}// setUseFixedLEOperators

void
IBStrategy::updateFixedLEOperators()
{
    TBOX_ERROR("IBStrategy::updateFixedLEOperators(): unimplemented\n");
    return;
}// updateFixedLEOperators

void
IBStrategy::computeLagrangianForceJacobianNonzeroStructure(
    std::vector<int>& /*d_nnz*/,
    std::vector<int>& /*o_nnz*/)
{
    TBOX_ERROR("IBStrategy::computeLagrangianForceJacobianNonzeroStructure(): unimplemented\n");
    return;
}// computeLagrangianForceJacobianNonzeroStructure

void
IBStrategy::computeLagrangianForceJacobian(
    Mat& /*J_mat*/,
    MatAssemblyType /*assembly_type*/,
    double /*X_coef*/,
    double /*U_coef*/,
    double /*data_time*/)
{
    TBOX_ERROR("IBStrategy::computeLagrangianForceJacobian(): unimplemented\n");
    return;
}// computeLagrangianForceJacobian

void
IBStrategy::applyLagrangianForceJacobian(
    int /*f_data_idx*/,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& /*f_prolongation_scheds*/,
    int /*u_data_idx*/,
    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*u_synch_scheds*/,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& /*u_ghost_fill_scheds*/,
    double /*data_time*/,
    Mat& /*J_mat*/)
{
    TBOX_ERROR("IBStrategy::applyLagrangianForceJacobian(): unimplemented\n");
    return;
}// applyLagrangianForceJacobian

bool
IBStrategy::hasFluidSources() const
{
    return false;
}// hasFluidSources

void
IBStrategy::computeLagrangianFluidSource(
    double /*data_time*/)
{
    // intentionally blank
    return;
}// computeLagrangianFluidSource

void
IBStrategy::spreadFluidSource(
    int /*q_data_idx*/,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& /*q_prolongation_scheds*/,
    double /*data_time*/)
{
    // intentionally blank
    return;
}// spreadFluidSource

void
IBStrategy::interpolatePressure(
    int /*p_data_idx*/,
    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*p_synch_scheds*/,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& /*p_ghost_fill_scheds*/,
    double /*data_time*/)
{
    // intentionally blank
    return;
}// interpolatePressure

void
IBStrategy::preprocessSolveFluidEquations(
    double /*current_time*/,
    double /*new_time*/,
    int /*cycle_num*/)
{
    // intentionally blank
    return;
}// preprocessSolveFluidEquations

void
IBStrategy::postprocessSolveFluidEquations(
    double /*current_time*/,
    double /*new_time*/,
    int /*cycle_num*/)
{
    // intentionally blank
    return;
}// postprocessSolveFluidEquations

void
IBStrategy::postprocessData()
{
    // intentionally blank
    return;
}// postprocessData

void
IBStrategy::initializePatchHierarchy(
    Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/,
    int /*u_data_idx*/,
    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*u_synch_scheds*/,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& /*u_ghost_fill_scheds*/,
    int /*integrator_step*/,
    double /*init_data_time*/,
    bool /*initial_time*/)
{
    // intentionally blank
    return;
}// initializePatchHierarchy

void
IBStrategy::registerLoadBalancer(
    Pointer<LoadBalancer<NDIM> > /*load_balancer*/,
    int /*workload_data_idx*/)
{
    // intentionally blank
    return;
}// registerLoadBalancer

void
IBStrategy::updateWorkloadEstimates(
    Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    int /*workload_data_idx*/)
{
    // intentionally blank
    return;
}// updateWorkloadEstimates

void
IBStrategy::beginDataRedistribution(
    Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    // intentionally blank
    return;
}// beginDataRedistribution

void
IBStrategy::endDataRedistribution(
    Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    // intentionally blank
    return;
}// endDataRedistribution

void
IBStrategy::initializeLevelData(
    Pointer<BasePatchHierarchy<NDIM> > /*hierarchy*/,
    int /*level_number*/,
    double /*init_data_time*/,
    bool /*can_be_refined*/,
    bool /*initial_time*/,
    Pointer<BasePatchLevel<NDIM> > /*old_level*/,
    bool /*allocate_data*/)
{
    // intentionally blank
    return;
}// initializeLevelData

void
IBStrategy::resetHierarchyConfiguration(
    Pointer<BasePatchHierarchy<NDIM> > /*hierarchy*/,
    int /*coarsest_level*/,
    int /*finest_level*/)
{
    // intentionally blank
    return;
}// resetHierarchyConfiguration

void
IBStrategy::applyGradientDetector(
    Pointer<BasePatchHierarchy<NDIM> > /*hierarchy*/,
    int /*level_number*/,
    double /*error_data_time*/,
    int /*tag_index*/,
    bool /*initial_time*/,
    bool /*uses_richardson_extrapolation_too*/)
{
    // intentionally blank
    return;
}// applyGradientDetector

void
IBStrategy::putToDatabase(
    Pointer<Database> /*db*/)
{
    // intentionally blank
    return;
}// putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

Pointer<HierarchyDataOpsReal<NDIM,double> >
IBStrategy::getVelocityHierarchyDataOps() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_ib_solver != NULL);
#endif
    return d_ib_solver->d_hier_velocity_data_ops;
}// getVelocityHierarchyDataOps

Pointer<HierarchyDataOpsReal<NDIM,double> >
IBStrategy::getPressureHierarchyDataOps() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_ib_solver != NULL);
#endif
    return d_ib_solver->d_hier_pressure_data_ops;
}// getPressureHierarchyDataOps

Pointer<HierarchyMathOps>
IBStrategy::getHierarchyMathOps() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_ib_solver != NULL);
#endif
    return d_ib_solver->d_hier_math_ops;
}// getHierarchyMathOps

void
IBStrategy::registerVariable(
    int& current_idx,
    int& new_idx,
    int& scratch_idx,
    Pointer<Variable<NDIM> > variable,
    const IntVector<NDIM>& scratch_ghosts,
    const std::string& coarsen_name,
    const std::string& refine_name,
    Pointer<CartGridFunction> init_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_ib_solver != NULL);
#endif
    d_ib_solver->registerVariable(current_idx, new_idx, scratch_idx, variable, scratch_ghosts, coarsen_name, refine_name, init_fcn);
    return;
}// registerVariable

void
IBStrategy::registerVariable(
    int& idx,
    Pointer<Variable<NDIM> > variable,
    const IntVector<NDIM>& ghosts,
    Pointer<VariableContext> ctx)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_ib_solver != NULL);
#endif
    d_ib_solver->registerVariable(idx, variable, ghosts, ctx);
    return;
}// registerVariable

void
IBStrategy::registerGhostfillRefineAlgorithm(
    const std::string& name,
    Pointer<RefineAlgorithm<NDIM> > ghostfill_alg,
    RefinePatchStrategy<NDIM>* ghostfill_patch_strategy)
{
    d_ib_solver->registerGhostfillRefineAlgorithm(name, ghostfill_alg, ghostfill_patch_strategy);
}// registerGhostfillRefineAlgorithm

void
IBStrategy::registerProlongRefineAlgorithm(
    const std::string& name,
    Pointer<RefineAlgorithm<NDIM> > prolong_alg,
    RefinePatchStrategy<NDIM>* prolong_patch_strategy)
{
    d_ib_solver->registerProlongRefineAlgorithm(name, prolong_alg, prolong_patch_strategy);
}// registerProlongRefineAlgorithm

void
IBStrategy::registerCoarsenAlgorithm(
    const std::string& name,
    Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg,
    CoarsenPatchStrategy<NDIM>* coarsen_patch_strategy)
{
    d_ib_solver->registerCoarsenAlgorithm(name, coarsen_alg, coarsen_patch_strategy);
}// registerCoarsenAlgorithm

Pointer<RefineAlgorithm<NDIM> >
IBStrategy::getGhostfillRefineAlgorithm(
    const std::string& name) const
{
    return d_ib_solver->getGhostfillRefineAlgorithm(name);
}// getGhostfillRefineAlgorithm

Pointer<RefineAlgorithm<NDIM> >
IBStrategy::getProlongRefineAlgorithm(
    const std::string& name) const
{
    return d_ib_solver->getProlongRefineAlgorithm(name);
}// getProlongRefineAlgorithm

Pointer<CoarsenAlgorithm<NDIM> >
IBStrategy::getCoarsenAlgorithm(
    const std::string& name) const
{
    return d_ib_solver->getCoarsenAlgorithm(name);
}// getCoarsenAlgorithm

const std::vector<Pointer<RefineSchedule<NDIM> > >&
IBStrategy::getGhostfillRefineSchedules(
    const std::string& name) const
{
    return d_ib_solver->getGhostfillRefineSchedules(name);
}// getGhostfillRefineSchedules

const std::vector<Pointer<RefineSchedule<NDIM> > >&
IBStrategy::getProlongRefineSchedules(
    const std::string& name) const
{
    return d_ib_solver->getProlongRefineSchedules(name);
}// getProlongRefineSchedules

const std::vector<Pointer<CoarsenSchedule<NDIM> > >&
IBStrategy::getCoarsenSchedules(
    const std::string& name) const
{
    return d_ib_solver->getCoarsenSchedules(name);
}// getCoarsenSchedules

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
