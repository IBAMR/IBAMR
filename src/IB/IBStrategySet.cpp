// Filename: IBStrategySet.cpp
// Created on 08 Mar 2012 by Boyce Griffith
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

#include <string>
#include <vector>

#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "PatchHierarchy.h"
#include "ibamr/IBStrategy.h"
#include "ibamr/IBStrategySet.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "tbox/Database.h"
#include "tbox/Pointer.h"

namespace IBTK
{
class RobinPhysBdryPatchStrategy;
} // namespace IBTK

namespace IBAMR
{
class IBHierarchyIntegrator;
} // namespace IBAMR
namespace SAMRAI
{
namespace tbox
{
template <class TYPE>
class Array;
} // namespace tbox
namespace xfer
{
template <int DIM>
class CoarsenSchedule;
template <int DIM>
class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBStrategySet::~IBStrategySet()
{
    // intentionally blank
    return;
} // ~IBStrategySet

void
IBStrategySet::registerIBHierarchyIntegrator(IBHierarchyIntegrator* ib_solver)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->registerIBHierarchyIntegrator(ib_solver);
    }
    return;
} // registerIBHierarchyIntegrator

void
IBStrategySet::registerEulerianVariables()
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->registerEulerianVariables();
    }
    return;
} // registerEulerianVariables

void
IBStrategySet::registerEulerianCommunicationAlgorithms()
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->registerEulerianCommunicationAlgorithms();
    }
    return;
} // registerEulerianCommunicationAlgorithms

const IntVector<NDIM>&
IBStrategySet::getMinimumGhostCellWidth() const
{
    static IntVector<NDIM> ghost_cell_width = 0;
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        ghost_cell_width = IntVector<NDIM>::max(ghost_cell_width, (*cit)->getMinimumGhostCellWidth());
    }
    return ghost_cell_width;
} // getMinimumGhostCellWidth

void
IBStrategySet::setupTagBuffer(Array<int>& tag_buffer, Pointer<GriddingAlgorithm<NDIM> > gridding_alg) const
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->setupTagBuffer(tag_buffer, gridding_alg);
    }
    return;
} // setupTagBuffer

void
IBStrategySet::preprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->preprocessIntegrateData(current_time, new_time, num_cycles);
    }
    return;
} // preprocessIntegrateData

void
IBStrategySet::postprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->postprocessIntegrateData(current_time, new_time, num_cycles);
    }
    return;
} // postprocessIntegrateData

void
IBStrategySet::updateFixedLEOperators()
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->updateFixedLEOperators();
    }
    return;
} // updateFixedLEOperators

void
IBStrategySet::interpolateVelocity(int u_data_idx,
                                   const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
                                   const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                   double data_time)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->interpolateVelocity(u_data_idx, u_synch_scheds, u_ghost_fill_scheds, data_time);
    }
    return;
} // interpolateVelocity

void
IBStrategySet::IBStrategySet::eulerStep(double current_time, double new_time)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->eulerStep(current_time, new_time);
    }
    return;
} // eulerStep

void
IBStrategySet::midpointStep(double current_time, double new_time)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->midpointStep(current_time, new_time);
    }
    return;
} // midpointStep

void
IBStrategySet::trapezoidalStep(double current_time, double new_time)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->trapezoidalStep(current_time, new_time);
    }
    return;
} // trapezoidalStep

void
IBStrategySet::computeLagrangianForce(double data_time)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->computeLagrangianForce(data_time);
    }
    return;
} // computeLagrangianForce

void
IBStrategySet::spreadForce(int f_data_idx,
                           RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                           const std::vector<Pointer<RefineSchedule<NDIM> > >& f_prolongation_scheds,
                           double data_time)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->spreadForce(f_data_idx, f_phys_bdry_op, f_prolongation_scheds, data_time);
    }
    return;
} // spreadForce

bool
IBStrategySet::hasFluidSources() const
{
    bool has_fluid_sources = false;
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        has_fluid_sources = has_fluid_sources || (*cit)->hasFluidSources();
    }
    return has_fluid_sources;
} // hasFluidSources

void
IBStrategySet::computeLagrangianFluidSource(double data_time)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->computeLagrangianFluidSource(data_time);
    }
    return;
} // computeLagrangianFluidSource

void
IBStrategySet::spreadFluidSource(int q_data_idx,
                                 const std::vector<Pointer<RefineSchedule<NDIM> > >& q_prolongation_scheds,
                                 double data_time)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->spreadFluidSource(q_data_idx, q_prolongation_scheds, data_time);
    }
    return;
} // spreadFluidSource

void
IBStrategySet::interpolatePressure(int p_data_idx,
                                   const std::vector<Pointer<CoarsenSchedule<NDIM> > >& p_synch_scheds,
                                   const std::vector<Pointer<RefineSchedule<NDIM> > >& p_ghost_fill_scheds,
                                   double data_time)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->interpolatePressure(p_data_idx, p_synch_scheds, p_ghost_fill_scheds, data_time);
    }
    return;
} // interpolatePressure

void
IBStrategySet::preprocessSolveFluidEquations(double current_time, double new_time, int cycle_num)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->preprocessSolveFluidEquations(current_time, new_time, cycle_num);
    }
    return;
} // preprocessSolveFluidEquations

void
IBStrategySet::postprocessSolveFluidEquations(double current_time, double new_time, int cycle_num)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->postprocessSolveFluidEquations(current_time, new_time, cycle_num);
    }
    return;
} // postprocessSolveFluidEquations

void
IBStrategySet::postprocessData()
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->postprocessData();
    }
    return;
} // postprocessData

void
IBStrategySet::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                        Pointer<GriddingAlgorithm<NDIM> > gridding_alg,
                                        int u_data_idx,
                                        const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
                                        const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                        int integrator_step,
                                        double init_data_time,
                                        bool initial_time)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->initializePatchHierarchy(hierarchy,
                                         gridding_alg,
                                         u_data_idx,
                                         u_synch_scheds,
                                         u_ghost_fill_scheds,
                                         integrator_step,
                                         init_data_time,
                                         initial_time);
    }
    return;
} // initializePatchHierarchy

void
IBStrategySet::registerLoadBalancer(Pointer<LoadBalancer<NDIM> > load_balancer, int workload_data_idx)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->registerLoadBalancer(load_balancer, workload_data_idx);
    }
    return;
} // registerLoadBalancer

void
IBStrategySet::updateWorkloadEstimates(Pointer<PatchHierarchy<NDIM> > hierarchy, int workload_data_idx)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->updateWorkloadEstimates(hierarchy, workload_data_idx);
    }
    return;
} // updateWorkloadEstimates

void
IBStrategySet::beginDataRedistribution(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                       Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->beginDataRedistribution(hierarchy, gridding_alg);
    }
    return;
} // beginDataRedistribution

void
IBStrategySet::endDataRedistribution(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                     Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->endDataRedistribution(hierarchy, gridding_alg);
    }
    return;
} // endDataRedistribution

void
IBStrategySet::initializeLevelData(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                   int level_number,
                                   double init_data_time,
                                   bool can_be_refined,
                                   bool initial_time,
                                   Pointer<BasePatchLevel<NDIM> > old_level,
                                   bool allocate_data)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->initializeLevelData(
            hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);
    }
    return;
} // initializeLevelData

void
IBStrategySet::resetHierarchyConfiguration(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                           int coarsest_level,
                                           int finest_level)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);
    }
    return;
} // resetHierarchyConfiguration

void
IBStrategySet::applyGradientDetector(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                     int level_number,
                                     double error_data_time,
                                     int tag_index,
                                     bool initial_time,
                                     bool uses_richardson_extrapolation_too)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->applyGradientDetector(
            hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    }
    return;
} // applyGradientDetector

void
IBStrategySet::putToDatabase(Pointer<Database> db)
{
    for (std::vector<Pointer<IBStrategy> >::const_iterator cit = d_strategy_set.begin(); cit != d_strategy_set.end();
         ++cit)
    {
        (*cit)->putToDatabase(db);
    }
    return;
} // putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
