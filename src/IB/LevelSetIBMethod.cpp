// Filename: LevelSetIBMethod.cpp
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
#include "ibamr/LevelSetIBMethod.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "tbox/Database.h"
#include "tbox/Pointer.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LevelSetIBMethod::~LevelSetIBMethod()
{
    // intentionally blank
    return;
} // ~LevelSetIBMethod

void
LevelSetIBMethod::registerIBHierarchyIntegrator(IBHierarchyIntegrator* ib_solver)
{
    // Register IBHierarchyIntegrator with the two IB ops.
    d_ib_method_ops->registerIBHierarchyIntegrator(ib_solver);
    d_ibfe_method_ops->registerIBHierarchyIntegrator(ib_solver);

    return;
} // registerIBHierarchyIntegrator

void
LevelSetIBMethod::registerEulerianVariables()
{
    d_ib_method_ops->registerEulerianVariables();
    d_ibfe_method_ops->registerEulerianVariables();

    return;
} // registerEulerianVariables

void
LevelSetIBMethod::registerEulerianCommunicationAlgorithms()
{
    d_ib_method_ops->registerEulerianCommunicationAlgorithms();
    d_ibfe_method_ops->registerEulerianCommunicationAlgorithms();

    return;
} // registerEulerianCommunicationAlgorithms

const IntVector<NDIM>&
LevelSetIBMethod::getMinimumGhostCellWidth() const
{
    return IntVector<NDIM>::max(d_ibfe_method_ops->getMinimumGhostCellWidth(),
                                d_ib_method_ops->getMinimumGhostCellWidth());
} // getMinimumGhostCellWidth

void
LevelSetIBMethod::setupTagBuffer(Array<int>& tag_buffer, Pointer<GriddingAlgorithm<NDIM> > gridding_alg) const
{
    d_ib_method_ops->setupTagBuffer(tag_buffer, gridding_alg);
    d_ibfe_method_ops->setupTagBuffer(tag_buffer, gridding_alg);

    return;
} // setupTagBuffer

void
LevelSetIBMethod::preprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    d_ib_method_ops->preprocessIntegrateData(current_time, new_time, num_cycles);
    d_ibfe_method_ops->preprocessIntegrateData(current_time, new_time, num_cycles);

    return;
} // preprocessIntegrateData

void
LevelSetIBMethod::postprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    d_ib_method_ops->postprocessIntegrateData(current_time, new_time, num_cycles);
    d_ibfe_method_ops->postprocessIntegrateData(current_time, new_time, num_cycles);

    return;
} // postprocessIntegrateData

void
LevelSetIBMethod::updateFixedLEOperators()
{
    d_ib_method_ops->updateFixedLEOperators();
    d_ibfe_method_ops->updateFixedLEOperators();

    return;
} // updateFixedLEOperators

void
LevelSetIBMethod::interpolateVelocity(int u_data_idx,
                                      const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
                                      const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                      double data_time)
{
    d_ib_method_ops->interpolateVelocity(u_data_idx, u_synch_scheds, u_ghost_fill_scheds, data_time);
    d_ibfe_method_ops->interpolateVelocity(u_data_idx, u_synch_scheds, u_ghost_fill_scheds, data_time);
    return;
} // interpolateVelocity

void
LevelSetIBMethod::LevelSetIBMethod::forwardEulerStep(double current_time, double new_time)
{
    d_ib_method_ops->forwardEulerStep(current_time, new_time);
    d_ibfe_method_ops->forwardEulerStep(current_time, new_time);
    return;
} // eulerStep

void
LevelSetIBMethod::midpointStep(double current_time, double new_time)
{
    d_ib_method_ops->midpointStep(current_time, new_time);
    d_ibfe_method_ops->midpointStep(current_time, new_time);
    return;
} // midpointStep

void
LevelSetIBMethod::trapezoidalStep(double current_time, double new_time)
{
    d_ib_method_ops->trapezoidalStep(current_time, new_time);
    d_ibfe_method_ops->trapezoidalStep(current_time, new_time);
    return;
} // trapezoidalStep

void
LevelSetIBMethod::computeLagrangianForce(double data_time)
{
    d_ib_method_ops->computeLagrangianForce(data_time);
    d_ibfe_method_ops->computeLagrangianForce(data_time);
    return;
} // computeLagrangianForce

void
LevelSetIBMethod::spreadForce(int f_data_idx,
                              RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                              const std::vector<Pointer<RefineSchedule<NDIM> > >& f_prolongation_scheds,
                              double data_time)
{
    d_ib_method_ops->spreadForce(f_data_idx, f_phys_bdry_op, f_prolongation_scheds, data_time);
    d_ibfe_method_ops->spreadForce(f_data_idx, f_phys_bdry_op, f_prolongation_scheds, data_time);

    return;
} // spreadForce

bool
LevelSetIBMethod::hasFluidSources() const
{
    return (d_ib_method_ops->hasFluidSources() || d_ibfe_method_ops->hasFluidSources());

} // hasFluidSources

void
LevelSetIBMethod::computeLagrangianFluidSource(double data_time)
{
    d_ib_method_ops->computeLagrangianFluidSource(data_time);
    d_ibfe_method_ops->computeLagrangianFluidSource(data_time);

    return;
} // computeLagrangianFluidSource

void
LevelSetIBMethod::spreadFluidSource(int q_data_idx,
                                    RobinPhysBdryPatchStrategy* q_phys_bdry_op,
                                    const std::vector<Pointer<RefineSchedule<NDIM> > >& q_prolongation_scheds,
                                    double data_time)
{
    d_ib_method_ops->spreadFluidSource(q_data_idx, q_phys_bdry_op, q_prolongation_scheds, data_time);
    d_ibfe_method_ops->spreadFluidSource(q_data_idx, q_phys_bdry_op, q_prolongation_scheds, data_time);
    return;
} // spreadFluidSource

void
LevelSetIBMethod::interpolatePressure(int p_data_idx,
                                      const std::vector<Pointer<CoarsenSchedule<NDIM> > >& p_synch_scheds,
                                      const std::vector<Pointer<RefineSchedule<NDIM> > >& p_ghost_fill_scheds,
                                      double data_time)
{
    d_ib_method_ops->interpolatePressure(p_data_idx, p_synch_scheds, p_ghost_fill_scheds, data_time);
    d_ibfe_method_ops->interpolatePressure(p_data_idx, p_synch_scheds, p_ghost_fill_scheds, data_time);
    return;
} // interpolatePressure

void
LevelSetIBMethod::preprocessSolveFluidEquations(double current_time, double new_time, int cycle_num)
{
    d_ib_method_ops->preprocessSolveFluidEquations(current_time, new_time, cycle_num);
    d_ibfe_method_ops->preprocessSolveFluidEquations(current_time, new_time, cycle_num);
    return;
} // preprocessSolveFluidEquations

void
LevelSetIBMethod::postprocessSolveFluidEquations(double current_time, double new_time, int cycle_num)
{
    d_ib_method_ops->postprocessSolveFluidEquations(current_time, new_time, cycle_num);
    d_ibfe_method_ops->postprocessSolveFluidEquations(current_time, new_time, cycle_num);
    return;
} // postprocessSolveFluidEquations

void
LevelSetIBMethod::postprocessData()
{
    d_ib_method_ops->postprocessData();
    d_ibfe_method_ops->postprocessData();
    return;
} // postprocessData

void
LevelSetIBMethod::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                           Pointer<GriddingAlgorithm<NDIM> > gridding_alg,
                                           int u_data_idx,
                                           const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
                                           const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                           int integrator_step,
                                           double init_data_time,
                                           bool initial_time)
{
    d_ib_method_ops->initializePatchHierarchy(hierarchy,
                                              gridding_alg,
                                              u_data_idx,
                                              u_synch_scheds,
                                              u_ghost_fill_scheds,
                                              integrator_step,
                                              init_data_time,
                                              initial_time);

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
LevelSetIBMethod::registerLoadBalancer(Pointer<LoadBalancer<NDIM> > load_balancer, int workload_data_idx)
{
    d_ib_method_ops->registerLoadBalancer(load_balancer, workload_data_idx);
    d_ibfe_method_ops->registerLoadBalancer(load_balancer, workload_data_idx);
    return;
} // registerLoadBalancer

void
LevelSetIBMethod::updateWorkloadEstimates(Pointer<PatchHierarchy<NDIM> > hierarchy, int workload_data_idx)
{
    d_ib_method_ops->updateWorkloadEstimates(hierarchy, workload_data_idx);
    d_ibfe_method_ops->updateWorkloadEstimates(hierarchy, workload_data_idx);
    return;
} // updateWorkloadEstimates

void
LevelSetIBMethod::beginDataRedistribution(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                          Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    d_ib_method_ops->beginDataRedistribution(hierarchy, gridding_alg);
    d_ibfe_method_ops->beginDataRedistribution(hierarchy, gridding_alg);
    return;
} // beginDataRedistribution

void
LevelSetIBMethod::endDataRedistribution(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                        Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    d_ib_method_ops->endDataRedistribution(hierarchy, gridding_alg);
    d_ibfe_method_ops->endDataRedistribution(hierarchy, gridding_alg);
    return;
} // endDataRedistribution

void
LevelSetIBMethod::initializeLevelData(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                      int level_number,
                                      double init_data_time,
                                      bool can_be_refined,
                                      bool initial_time,
                                      Pointer<BasePatchLevel<NDIM> > old_level,
                                      bool allocate_data)
{
    d_ib_method_ops->initializeLevelData(
        hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);

    d_ibfe_method_ops->initializeLevelData(
        hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);

    return;
} // initializeLevelData

void
LevelSetIBMethod::resetHierarchyConfiguration(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                              int coarsest_level,
                                              int finest_level)
{
    d_ib_method_ops->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);
    d_ibfe_method_ops->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);
    return;
} // resetHierarchyConfiguration

void
LevelSetIBMethod::applyGradientDetector(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                        int level_number,
                                        double error_data_time,
                                        int tag_index,
                                        bool initial_time,
                                        bool uses_richardson_extrapolation_too)
{
    d_ib_method_ops->applyGradientDetector(
        hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    d_ibfe_method_ops->applyGradientDetector(
        hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);

    return;
} // applyGradientDetector

void
LevelSetIBMethod::putToDatabase(Pointer<Database> db)
{
    d_ib_method_ops->putToDatabase(db);
    d_ibFE_method_ops->putToDatabase(db);
    return;
} // putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
