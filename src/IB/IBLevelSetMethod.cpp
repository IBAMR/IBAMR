// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/IBLevelSetMethod.h"

#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <string>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBLevelSetMethod::IBLevelSetMethod(Pointer<IBStrategy> ib_method_ops, Pointer<IBStrategy> ibfe_method_ops)
    : d_ib_method_ops(ib_method_ops), d_ibfe_method_ops(ibfe_method_ops)
{
    // intentionally blank
    return;
} // IBLevelSetMethod

void
IBLevelSetMethod::deactivateIBFEMethod()
{
    d_ibfe_method_ops.setNull();
    d_ibfe_method_ops = nullptr;
    d_ibfe_method_deactivated = true;
    return;
} // deactivateIBFEMethod

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
    static IntVector<NDIM> gcw = 0;
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
IBLevelSetMethod::addWorkloadEstimate(Pointer<PatchHierarchy<NDIM> > hierarchy, const int workload_data_idx)
{
    if (d_ib_method_ops) d_ib_method_ops->addWorkloadEstimate(hierarchy, workload_data_idx);
    if (d_ibfe_method_ops) d_ibfe_method_ops->addWorkloadEstimate(hierarchy, workload_data_idx);
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
