// Filename: IBLevelSetMethod.cpp
// Created on 12 Oct 2018 by Amneet Bhalla
//
// Copyright (c) 2002-2018, Amneet Bhalla
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

#include <string>
#include <vector>

#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "PatchHierarchy.h"
#include "ibamr/IBLevelSetMethod.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "tbox/Database.h"
#include "tbox/Pointer.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
IBLevelSetMethod::deactivateIBFEMethod()
{
    d_ibfe_method_ops.setNull();
    d_ibfe_method_ops = NULL;
    return;
} // deactivateIBFEMethod

IBLevelSetMethod::~IBLevelSetMethod()
{
    // intentionally blank
    return;
} // ~IBLevelSetMethod

void
IBLevelSetMethod::registerIBHierarchyIntegrator(IBHierarchyIntegrator* ib_solver)
{
    // Register IBHierarchyIntegrator with the two IB ops.
    if (d_ib_method_ops) d_ib_method_ops->registerIBHierarchyIntegrator(ib_solver);
    if (d_ibfe_method_ops) d_ibfe_method_ops->registerIBHierarchyIntegrator(ib_solver);

    return;
} // registerIBHierarchyIntegrator

void
IBLevelSetMethod::registerEulerianVariables()
{
    if (d_ib_method_ops) d_ib_method_ops->registerEulerianVariables();
    if (d_ibfe_method_ops) d_ibfe_method_ops->registerEulerianVariables();

    return;
} // registerEulerianVariables

void
IBLevelSetMethod::registerEulerianCommunicationAlgorithms()
{
    if (d_ib_method_ops) d_ib_method_ops->registerEulerianCommunicationAlgorithms();
    if (d_ibfe_method_ops) d_ibfe_method_ops->registerEulerianCommunicationAlgorithms();

    return;
} // registerEulerianCommunicationAlgorithms

const IntVector<NDIM>&
IBLevelSetMethod::getMinimumGhostCellWidth() const
{
    IntVector<NDIM> gcw = 0;
    if (d_ib_method_ops) gcw = IntVector<NDIM>::max(gcw, d_ib_method_ops->getMinimumGhostCellWidth());
    if (d_ibfe_method_ops) gcw = IntVector<NDIM>::max(gcw, d_ibfe_method_ops->getMinimumGhostCellWidth());
    return gcw;
} // getMinimumGhostCellWidth

void
IBLevelSetMethod::setupTagBuffer(Array<int>& tag_buffer, Pointer<GriddingAlgorithm<NDIM> > gridding_alg) const
{
    if (d_ib_method_ops) d_ib_method_ops->setupTagBuffer(tag_buffer, gridding_alg);
    if (d_ibfe_method_ops) d_ibfe_method_ops->setupTagBuffer(tag_buffer, gridding_alg);

    return;
} // setupTagBuffer

void
IBLevelSetMethod::preprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    if (d_ib_method_ops) d_ib_method_ops->preprocessIntegrateData(current_time, new_time, num_cycles);
    if (d_ibfe_method_ops) d_ibfe_method_ops->preprocessIntegrateData(current_time, new_time, num_cycles);

    return;
} // preprocessIntegrateData

void
IBLevelSetMethod::postprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    if (d_ib_method_ops) d_ib_method_ops->postprocessIntegrateData(current_time, new_time, num_cycles);
    if (d_ibfe_method_ops) d_ibfe_method_ops->postprocessIntegrateData(current_time, new_time, num_cycles);

    return;
} // postprocessIntegrateData

void
IBLevelSetMethod::updateFixedLEOperators()
{
    if (d_ib_method_ops) d_ib_method_ops->updateFixedLEOperators();
    if (d_ibfe_method_ops) d_ibfe_method_ops->updateFixedLEOperators();

    return;
} // updateFixedLEOperators

void
IBLevelSetMethod::interpolateVelocity(int u_data_idx,
                                      const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
                                      const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                      double data_time)
{
    if (d_ib_method_ops)
        d_ib_method_ops->interpolateVelocity(u_data_idx, u_synch_scheds, u_ghost_fill_scheds, data_time);
    if (d_ibfe_method_ops)
        d_ibfe_method_ops->interpolateVelocity(u_data_idx, u_synch_scheds, u_ghost_fill_scheds, data_time);
    return;
} // interpolateVelocity

void
IBLevelSetMethod::IBLevelSetMethod::forwardEulerStep(double current_time, double new_time)
{
    if (d_ib_method_ops) d_ib_method_ops->forwardEulerStep(current_time, new_time);
    if (d_ibfe_method_ops) d_ibfe_method_ops->forwardEulerStep(current_time, new_time);
    return;
} // eulerStep

void
IBLevelSetMethod::midpointStep(double current_time, double new_time)
{
    if (d_ib_method_ops) d_ib_method_ops->midpointStep(current_time, new_time);
    if (d_ibfe_method_ops) d_ibfe_method_ops->midpointStep(current_time, new_time);
    return;
} // midpointStep

void
IBLevelSetMethod::trapezoidalStep(double current_time, double new_time)
{
    if (d_ib_method_ops) d_ib_method_ops->trapezoidalStep(current_time, new_time);
    if (d_ibfe_method_ops) d_ibfe_method_ops->trapezoidalStep(current_time, new_time);
    return;
} // trapezoidalStep

void
IBLevelSetMethod::computeLagrangianForce(double data_time)
{
    if (d_ib_method_ops) d_ib_method_ops->computeLagrangianForce(data_time);
    if (d_ibfe_method_ops) d_ibfe_method_ops->computeLagrangianForce(data_time);
    return;
} // computeLagrangianForce

void
IBLevelSetMethod::spreadForce(int f_data_idx,
                              RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                              const std::vector<Pointer<RefineSchedule<NDIM> > >& f_prolongation_scheds,
                              double data_time)
{
    if (d_ib_method_ops) d_ib_method_ops->spreadForce(f_data_idx, f_phys_bdry_op, f_prolongation_scheds, data_time);
    if (d_ibfe_method_ops) d_ibfe_method_ops->spreadForce(f_data_idx, f_phys_bdry_op, f_prolongation_scheds, data_time);
    return;
} // spreadForce

bool
IBLevelSetMethod::hasFluidSources() const
{
    bool has_fluid_sources = false;
    if (d_ib_method_ops) has_fluid_sources = has_fluid_sources || d_ib_method_ops->hasFluidSources();
    if (d_ibfe_method_ops) has_fluid_sources = has_fluid_sources || d_ibfe_method_ops->hasFluidSources();
    return has_fluid_sources;
} // hasFluidSources

void
IBLevelSetMethod::computeLagrangianFluidSource(double data_time)
{
    if (d_ib_method_ops) d_ib_method_ops->computeLagrangianFluidSource(data_time);
    if (d_ibfe_method_ops) d_ibfe_method_ops->computeLagrangianFluidSource(data_time);
    return;
} // computeLagrangianFluidSource

void
IBLevelSetMethod::spreadFluidSource(int q_data_idx,
                                    RobinPhysBdryPatchStrategy* q_phys_bdry_op,
                                    const std::vector<Pointer<RefineSchedule<NDIM> > >& q_prolongation_scheds,
                                    double data_time)
{
    if (d_ib_method_ops)
        d_ib_method_ops->spreadFluidSource(q_data_idx, q_phys_bdry_op, q_prolongation_scheds, data_time);
    if (d_ibfe_method_ops)
        d_ibfe_method_ops->spreadFluidSource(q_data_idx, q_phys_bdry_op, q_prolongation_scheds, data_time);
    return;
} // spreadFluidSource

void
IBLevelSetMethod::interpolatePressure(int p_data_idx,
                                      const std::vector<Pointer<CoarsenSchedule<NDIM> > >& p_synch_scheds,
                                      const std::vector<Pointer<RefineSchedule<NDIM> > >& p_ghost_fill_scheds,
                                      double data_time)
{
    if (d_ib_method_ops)
        d_ib_method_ops->interpolatePressure(p_data_idx, p_synch_scheds, p_ghost_fill_scheds, data_time);
    if (d_ibfe_method_ops)
        d_ibfe_method_ops->interpolatePressure(p_data_idx, p_synch_scheds, p_ghost_fill_scheds, data_time);
    return;
} // interpolatePressure

void
IBLevelSetMethod::preprocessSolveFluidEquations(double current_time, double new_time, int cycle_num)
{
    if (d_ib_method_ops) d_ib_method_ops->preprocessSolveFluidEquations(current_time, new_time, cycle_num);
    if (d_ibfe_method_ops) d_ibfe_method_ops->preprocessSolveFluidEquations(current_time, new_time, cycle_num);
    return;
} // preprocessSolveFluidEquations

void
IBLevelSetMethod::postprocessSolveFluidEquations(double current_time, double new_time, int cycle_num)
{
    if (d_ib_method_ops) d_ib_method_ops->postprocessSolveFluidEquations(current_time, new_time, cycle_num);
    if (d_ibfe_method_ops) d_ibfe_method_ops->postprocessSolveFluidEquations(current_time, new_time, cycle_num);
    return;
} // postprocessSolveFluidEquations

void
IBLevelSetMethod::postprocessData()
{
    if (d_ib_method_ops) d_ib_method_ops->postprocessData();
    if (d_ibfe_method_ops) d_ibfe_method_ops->postprocessData();
    return;
} // postprocessData

void
IBLevelSetMethod::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                           Pointer<GriddingAlgorithm<NDIM> > gridding_alg,
                                           int u_data_idx,
                                           const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
                                           const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                           int integrator_step,
                                           double init_data_time,
                                           bool initial_time)
{
    if (d_ib_method_ops)
        d_ib_method_ops->initializePatchHierarchy(hierarchy,
                                                  gridding_alg,
                                                  u_data_idx,
                                                  u_synch_scheds,
                                                  u_ghost_fill_scheds,
                                                  integrator_step,
                                                  init_data_time,
                                                  initial_time);

    if (d_ibfe_method_ops)
        d_ibfe_method_ops->initializePatchHierarchy(hierarchy,
                                                    gridding_alg,
                                                    u_data_idx,
                                                    u_synch_scheds,
                                                    u_ghost_fill_scheds,
                                                    integrator_step,
                                                    init_data_time,
                                                    initial_time);

    return;
} // initializePatchHierarchy

void
IBLevelSetMethod::registerLoadBalancer(Pointer<LoadBalancer<NDIM> > load_balancer, int workload_data_idx)
{
    if (d_ib_method_ops) d_ib_method_ops->registerLoadBalancer(load_balancer, workload_data_idx);
    if (d_ibfe_method_ops) d_ibfe_method_ops->registerLoadBalancer(load_balancer, workload_data_idx);
    return;
} // registerLoadBalancer

void
IBLevelSetMethod::updateWorkloadEstimates(Pointer<PatchHierarchy<NDIM> > hierarchy, int workload_data_idx)
{
    if (d_ib_method_ops) d_ib_method_ops->updateWorkloadEstimates(hierarchy, workload_data_idx);
    if (d_ibfe_method_ops) d_ibfe_method_ops->updateWorkloadEstimates(hierarchy, workload_data_idx);
    return;
} // updateWorkloadEstimates

void
IBLevelSetMethod::beginDataRedistribution(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                          Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_ib_method_ops) d_ib_method_ops->beginDataRedistribution(hierarchy, gridding_alg);
    if (d_ibfe_method_ops) d_ibfe_method_ops->beginDataRedistribution(hierarchy, gridding_alg);
    return;
} // beginDataRedistribution

void
IBLevelSetMethod::endDataRedistribution(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                        Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_ib_method_ops) d_ib_method_ops->endDataRedistribution(hierarchy, gridding_alg);
    if (d_ibfe_method_ops) d_ibfe_method_ops->endDataRedistribution(hierarchy, gridding_alg);
    return;
} // endDataRedistribution

void
IBLevelSetMethod::initializeLevelData(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                      int level_number,
                                      double init_data_time,
                                      bool can_be_refined,
                                      bool initial_time,
                                      Pointer<BasePatchLevel<NDIM> > old_level,
                                      bool allocate_data)
{
    if (d_ib_method_ops)
        d_ib_method_ops->initializeLevelData(
            hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);

    if (d_ibfe_method_ops)
        d_ibfe_method_ops->initializeLevelData(
            hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);

    return;
} // initializeLevelData

void
IBLevelSetMethod::resetHierarchyConfiguration(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                              int coarsest_level,
                                              int finest_level)
{
    if (d_ib_method_ops) d_ib_method_ops->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);
    if (d_ibfe_method_ops) d_ibfe_method_ops->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);
    return;
} // resetHierarchyConfiguration

void
IBLevelSetMethod::applyGradientDetector(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                        int level_number,
                                        double error_data_time,
                                        int tag_index,
                                        bool initial_time,
                                        bool uses_richardson_extrapolation_too)
{
    if (d_ib_method_ops)
        d_ib_method_ops->applyGradientDetector(
            hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    if (d_ibfe_method_ops)
        d_ibfe_method_ops->applyGradientDetector(
            hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);

    return;
} // applyGradientDetector

void
IBLevelSetMethod::putToDatabase(Pointer<Database> db)
{
    if (d_ib_method_ops) d_ib_method_ops->putToDatabase(db);
    if (d_ibfe_method_ops) d_ibfe_method_ops->putToDatabase(db);
    return;
} // putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
