// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
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

#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/IBStrategy.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep

#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyMathOps.h"

#include "BasePatchLevel.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenPatchStrategy.h"
#include "GriddingAlgorithm.h"
#include "HierarchyDataOpsReal.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "RefineAlgorithm.h"
#include "RefinePatchStrategy.h"
#include "Variable.h"
#include "VariableContext.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace IBTK
{
class RobinPhysBdryPatchStrategy;
} // namespace IBTK

namespace IBAMR
{
class INSHierarchyIntegrator;
} // namespace IBAMR

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BasePatchHierarchy;
} // namespace hier

namespace mesh
{
template <int DIM>
class LoadBalancer;
} // namespace mesh
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

void
IBStrategy::registerIBHierarchyIntegrator(IBHierarchyIntegrator* ib_solver)
{
    d_ib_solver = ib_solver;
    return;
} // registerIBHierarchyIntegrator

void
IBStrategy::registerEulerianVariables()
{
    // intentionally blank
    return;
} // registerEulerianVariables

void
IBStrategy::registerEulerianCommunicationAlgorithms()
{
    // intentionally blank
    return;
} // registerEulerianCommunicationAlgorithms

void
IBStrategy::setupTagBuffer(Array<int>& tag_buffer, Pointer<GriddingAlgorithm<NDIM> > gridding_alg) const
{
    const int finest_hier_ln = gridding_alg->getMaxLevels() - 1;
    const int tsize = tag_buffer.size();
    tag_buffer.resizeArray(finest_hier_ln);
    for (int i = tsize; i < finest_hier_ln; ++i) tag_buffer[i] = 0;
    const int gcw = getMinimumGhostCellWidth().max();
    for (int i = 0; i < tag_buffer.size(); ++i)
    {
        tag_buffer[i] = std::max(tag_buffer[i], gcw);
    }
    return;
} // setupTagBuffer

void
IBStrategy::inactivateLagrangianStructure(int /*structure_number*/, int /*level_number*/)
{
    TBOX_ERROR("IBStrategy::inactivateLagrangianStructure(): unimplemented\n");
    return;
} // inactivateLagrangianStructure

void
IBStrategy::activateLagrangianStructure(int /*structure_number*/, int /*level_number*/)
{
    TBOX_ERROR("IBStrategy::activateLagrangianStructure(): unimplemented\n");
    return;
} // activateLagrangianStructure

bool
IBStrategy::getLagrangianStructureIsActivated(int /*structure_number*/, int /*level_number*/) const
{
    TBOX_ERROR("IBStrategy::getLagrangianStructureIsActivated(): unimplemented\n");
    return true;
} // activateLagrangianStructure

void
IBStrategy::preprocessIntegrateData(double /*current_time*/, double /*new_time*/, int /*num_cycles*/)
{
    // intentionally blank
    return;
} // preprocessIntegrateData

void
IBStrategy::postprocessIntegrateData(double /*current_time*/, double /*new_time*/, int /*num_cycles*/)
{
    // intentionally blank
    return;
} // postprocessIntegrateData

void
IBStrategy::setUseFixedLEOperators(bool use_fixed_coupling_ops)
{
    d_use_fixed_coupling_ops = use_fixed_coupling_ops;
    return;
} // setUseFixedLEOperators

void
IBStrategy::updateFixedLEOperators()
{
    TBOX_ERROR("IBStrategy::updateFixedLEOperators(): unimplemented\n");
    return;
} // updateFixedLEOperators

void
IBStrategy::backwardEulerStep(double /*current_time*/, double /*new_time*/)
{
    TBOX_ERROR("IBStrategy::backwardEulerStep(): unimplemented\n");
    return;
} // backwardEulerStep

bool
IBStrategy::hasFluidSources() const
{
    return false;
} // hasFluidSources

void
IBStrategy::computeLagrangianFluidSource(double /*data_time*/)
{
    // intentionally blank
    return;
} // computeLagrangianFluidSource

void
IBStrategy::spreadFluidSource(int /*q_data_idx*/,
                              RobinPhysBdryPatchStrategy* /*q_phys_bdry_op*/,
                              const std::vector<Pointer<RefineSchedule<NDIM> > >& /*q_prolongation_scheds*/,
                              double /*data_time*/)
{
    // intentionally blank
    return;
} // spreadFluidSource

void
IBStrategy::interpolatePressure(int /*p_data_idx*/,
                                const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*p_synch_scheds*/,
                                const std::vector<Pointer<RefineSchedule<NDIM> > >& /*p_ghost_fill_scheds*/,
                                double /*data_time*/)
{
    // intentionally blank
    return;
} // interpolatePressure

void
IBStrategy::preprocessSolveFluidEquations(double /*current_time*/, double /*new_time*/, int /*cycle_num*/)
{
    // intentionally blank
    return;
} // preprocessSolveFluidEquations

void
IBStrategy::postprocessSolveFluidEquations(double /*current_time*/, double /*new_time*/, int /*cycle_num*/)
{
    // intentionally blank
    return;
} // postprocessSolveFluidEquations

void
IBStrategy::postprocessData()
{
    // intentionally blank
    return;
} // postprocessData

void
IBStrategy::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
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
} // initializePatchHierarchy

void
IBStrategy::registerLoadBalancer(Pointer<LoadBalancer<NDIM> > /*load_balancer*/, int /*workload_data_idx*/)
{
    IBAMR_DEPRECATED_MEMBER_FUNCTION1("IBStrategy", "registerLoadBalancer");
    // intentionally blank
    return;
} // registerLoadBalancer

void
IBStrategy::addWorkloadEstimate(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/, const int /*workload_data_idx*/)
{
    // intentionally blank
    return;
} // addWorkloadEstimate

void IBStrategy::beginDataRedistribution(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                         Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    // intentionally blank
    return;
} // beginDataRedistribution

void IBStrategy::endDataRedistribution(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                       Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    // intentionally blank
    return;
} // endDataRedistribution

void
IBStrategy::initializeLevelData(Pointer<BasePatchHierarchy<NDIM> > /*hierarchy*/,
                                int /*level_number*/,
                                double /*init_data_time*/,
                                bool /*can_be_refined*/,
                                bool /*initial_time*/,
                                Pointer<BasePatchLevel<NDIM> > /*old_level*/,
                                bool /*allocate_data*/)
{
    // intentionally blank
    return;
} // initializeLevelData

void
IBStrategy::resetHierarchyConfiguration(Pointer<BasePatchHierarchy<NDIM> > /*hierarchy*/,
                                        int /*coarsest_level*/,
                                        int /*finest_level*/)
{
    // intentionally blank
    return;
} // resetHierarchyConfiguration

void
IBStrategy::applyGradientDetector(Pointer<BasePatchHierarchy<NDIM> > /*hierarchy*/,
                                  int /*level_number*/,
                                  double /*error_data_time*/,
                                  int /*tag_index*/,
                                  bool /*initial_time*/,
                                  bool /*uses_richardson_extrapolation_too*/)
{
    // intentionally blank
    return;
} // applyGradientDetector

void IBStrategy::putToDatabase(Pointer<Database> /*db*/)
{
    // intentionally blank
    return;
} // putToDatabase

void
IBStrategy::registerMaskingPatchDataIndex(int mask_data_idx)
{
    d_mask_data_idx = mask_data_idx;
    return;
} // registerMaskingPatchDataIndex

int
IBStrategy::getMaskingPatchDataIndex()
{
    return d_mask_data_idx;
} // getMaskingPatchDataIndex

/////////////////////////////// PROTECTED ////////////////////////////////////

INSHierarchyIntegrator*
IBStrategy::getINSHierarchyIntegrator() const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_ib_solver);
#endif
    return d_ib_solver->d_ins_hier_integrator;
} // getINSHierarchyIntegrator

Pointer<HierarchyDataOpsReal<NDIM, double> >
IBStrategy::getVelocityHierarchyDataOps() const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_ib_solver);
#endif
    return d_ib_solver->d_hier_velocity_data_ops;
} // getVelocityHierarchyDataOps

Pointer<HierarchyDataOpsReal<NDIM, double> >
IBStrategy::getPressureHierarchyDataOps() const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_ib_solver);
#endif
    return d_ib_solver->d_hier_pressure_data_ops;
} // getPressureHierarchyDataOps

Pointer<HierarchyMathOps>
IBStrategy::getHierarchyMathOps() const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_ib_solver);
#endif
    return d_ib_solver->d_hier_math_ops;
} // getHierarchyMathOps

void
IBStrategy::registerVariable(int& current_idx,
                             int& new_idx,
                             int& scratch_idx,
                             Pointer<Variable<NDIM> > variable,
                             const IntVector<NDIM>& scratch_ghosts,
                             const std::string& coarsen_name,
                             const std::string& refine_name,
                             Pointer<CartGridFunction> init_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_ib_solver);
#endif
    d_ib_solver->registerVariable(
        current_idx, new_idx, scratch_idx, variable, scratch_ghosts, coarsen_name, refine_name, init_fcn);
    return;
} // registerVariable

void
IBStrategy::registerVariable(int& idx,
                             Pointer<Variable<NDIM> > variable,
                             const IntVector<NDIM>& ghosts,
                             Pointer<VariableContext> ctx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_ib_solver);
#endif
    d_ib_solver->registerVariable(idx, variable, ghosts, ctx);
    return;
} // registerVariable

void
IBStrategy::registerGhostfillRefineAlgorithm(const std::string& name,
                                             Pointer<RefineAlgorithm<NDIM> > ghostfill_alg,
                                             std::unique_ptr<RefinePatchStrategy<NDIM> > ghostfill_patch_strategy)
{
    d_ib_solver->registerGhostfillRefineAlgorithm(name, ghostfill_alg, std::move(ghostfill_patch_strategy));
} // registerGhostfillRefineAlgorithm

void
IBStrategy::registerProlongRefineAlgorithm(const std::string& name,
                                           Pointer<RefineAlgorithm<NDIM> > prolong_alg,
                                           std::unique_ptr<RefinePatchStrategy<NDIM> > prolong_patch_strategy)
{
    d_ib_solver->registerProlongRefineAlgorithm(name, prolong_alg, std::move(prolong_patch_strategy));
} // registerProlongRefineAlgorithm

void
IBStrategy::registerCoarsenAlgorithm(const std::string& name,
                                     Pointer<CoarsenAlgorithm<NDIM> > coarsen_alg,
                                     std::unique_ptr<CoarsenPatchStrategy<NDIM> > coarsen_patch_strategy)
{
    d_ib_solver->registerCoarsenAlgorithm(name, coarsen_alg, std::move(coarsen_patch_strategy));
} // registerCoarsenAlgorithm

Pointer<RefineAlgorithm<NDIM> >
IBStrategy::getGhostfillRefineAlgorithm(const std::string& name) const
{
    return d_ib_solver->getGhostfillRefineAlgorithm(name);
} // getGhostfillRefineAlgorithm

Pointer<RefineAlgorithm<NDIM> >
IBStrategy::getProlongRefineAlgorithm(const std::string& name) const
{
    return d_ib_solver->getProlongRefineAlgorithm(name);
} // getProlongRefineAlgorithm

Pointer<CoarsenAlgorithm<NDIM> >
IBStrategy::getCoarsenAlgorithm(const std::string& name) const
{
    return d_ib_solver->getCoarsenAlgorithm(name);
} // getCoarsenAlgorithm

const std::vector<Pointer<RefineSchedule<NDIM> > >&
IBStrategy::getGhostfillRefineSchedules(const std::string& name) const
{
    return d_ib_solver->getGhostfillRefineSchedules(name);
} // getGhostfillRefineSchedules

const std::vector<Pointer<RefineSchedule<NDIM> > >&
IBStrategy::getProlongRefineSchedules(const std::string& name) const
{
    return d_ib_solver->getProlongRefineSchedules(name);
} // getProlongRefineSchedules

const std::vector<Pointer<CoarsenSchedule<NDIM> > >&
IBStrategy::getCoarsenSchedules(const std::string& name) const
{
    return d_ib_solver->getCoarsenSchedules(name);
} // getCoarsenSchedules

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
