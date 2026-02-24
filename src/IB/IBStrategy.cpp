// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2025 by the IBAMR developers
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

#include <ibamr/IBHierarchyIntegrator.h>
#include <ibamr/IBStrategy.h>
#include <ibamr/ibamr_utilities.h>

#include <ibtk/CartGridFunction.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/samrai_compatibility_names.h>

#include <SAMRAIArray.h>
#include <SAMRAIBasePatchHierarchy.h>
#include <SAMRAIBasePatchLevel.h>
#include <SAMRAICoarsenAlgorithm.h>
#include <SAMRAICoarsenPatchStrategy.h>
#include <SAMRAICoarsenSchedule.h>
#include <SAMRAIDatabase.h>
#include <SAMRAIGriddingAlgorithm.h>
#include <SAMRAIHierarchyDataOpsReal.h>
#include <SAMRAIIntVector.h>
#include <SAMRAILoadBalancer.h>
#include <SAMRAIPatchHierarchy.h>
#include <SAMRAIPointer.h>
#include <SAMRAIRefineAlgorithm.h>
#include <SAMRAIRefinePatchStrategy.h>
#include <SAMRAIRefineSchedule.h>
#include <SAMRAIUtilities.h>
#include <SAMRAIVariable.h>
#include <SAMRAIVariableContext.h>

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <ibamr/namespaces.h> // IWYU pragma: keep

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
IBStrategy::setupTagBuffer(SAMRAIArray<int>& tag_buffer, SAMRAIPointer<SAMRAIGriddingAlgorithm> gridding_alg) const
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

double
IBStrategy::getMaxPointDisplacement() const
{
    TBOX_ERROR(
        "IBStrategy::getMaxPointDisplacement(): This method, which permits regridding based on the maximum "
        "displacement of the structure, is not yet implemented by the inheriting class. You must use regridding "
        "based on the fluid CFL estimate instead. this is typically done by setting regrid_fluid_cfl_interval "
        "and *not* setting regrid_structure_cfl_interval in the input database for IBHierarchyIntegrator.\n");
    return std::numeric_limits<double>::max();
} // getMaxPointDisplacement

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
IBStrategy::setUseMultistepTimeStepping(const unsigned int /*n_previous_steps*/)
{
    TBOX_ERROR("IBStrategy::setUseMultistepTimeStepping(): unimplemented\n");
    return;
} // setUseMultistepTimeStepping

void
IBStrategy::backwardEulerStep(double /*current_time*/, double /*new_time*/)
{
    TBOX_ERROR("IBStrategy::backwardEulerStep(): unimplemented\n");
    return;
} // backwardEulerStep

void
IBStrategy::AB2Step(double /*current_time*/, double /*new_time*/)
{
    TBOX_ERROR("IBStrategy::AB2Step(): unimplemented\n");
    return;
} // AB2Step

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
                              const std::vector<SAMRAIPointer<SAMRAIRefineSchedule>>& /*q_prolongation_scheds*/,
                              double /*data_time*/)
{
    // intentionally blank
    return;
} // spreadFluidSource

void
IBStrategy::interpolatePressure(int /*p_data_idx*/,
                                const std::vector<SAMRAIPointer<SAMRAICoarsenSchedule>>& /*p_synch_scheds*/,
                                const std::vector<SAMRAIPointer<SAMRAIRefineSchedule>>& /*p_ghost_fill_scheds*/,
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
IBStrategy::initializePatchHierarchy(SAMRAIPointer<SAMRAIPatchHierarchy> /*hierarchy*/,
                                     SAMRAIPointer<SAMRAIGriddingAlgorithm> /*gridding_alg*/,
                                     int /*u_data_idx*/,
                                     const std::vector<SAMRAIPointer<SAMRAICoarsenSchedule>>& /*u_synch_scheds*/,
                                     const std::vector<SAMRAIPointer<SAMRAIRefineSchedule>>& /*u_ghost_fill_scheds*/,
                                     int /*integrator_step*/,
                                     double /*init_data_time*/,
                                     bool /*initial_time*/)
{
    // intentionally blank
    return;
} // initializePatchHierarchy

void
IBStrategy::registerLoadBalancer(SAMRAIPointer<SAMRAILoadBalancer> /*load_balancer*/, int /*workload_data_idx*/)
{
    // intentionally blank
    return;
} // registerLoadBalancer

void
IBStrategy::addWorkloadEstimate(SAMRAIPointer<SAMRAIPatchHierarchy> /*hierarchy*/, const int /*workload_data_idx*/)
{
    // intentionally blank
    return;
} // addWorkloadEstimate

void
IBStrategy::beginDataRedistribution(SAMRAIPointer<SAMRAIPatchHierarchy> /*hierarchy*/,
                                    SAMRAIPointer<SAMRAIGriddingAlgorithm> /*gridding_alg*/)
{
    // intentionally blank
    return;
} // beginDataRedistribution

void
IBStrategy::endDataRedistribution(SAMRAIPointer<SAMRAIPatchHierarchy> /*hierarchy*/,
                                  SAMRAIPointer<SAMRAIGriddingAlgorithm> /*gridding_alg*/)
{
    // intentionally blank
    return;
} // endDataRedistribution

void
IBStrategy::initializeLevelData(SAMRAIPointer<SAMRAIBasePatchHierarchy> /*hierarchy*/,
                                int /*level_number*/,
                                double /*init_data_time*/,
                                bool /*can_be_refined*/,
                                bool /*initial_time*/,
                                SAMRAIPointer<SAMRAIBasePatchLevel> /*old_level*/,
                                bool /*allocate_data*/)
{
    // intentionally blank
    return;
} // initializeLevelData

void
IBStrategy::resetHierarchyConfiguration(SAMRAIPointer<SAMRAIBasePatchHierarchy> /*hierarchy*/,
                                        int /*coarsest_level*/,
                                        int /*finest_level*/)
{
    // intentionally blank
    return;
} // resetHierarchyConfiguration

void
IBStrategy::applyGradientDetector(SAMRAIPointer<SAMRAIBasePatchHierarchy> /*hierarchy*/,
                                  int /*level_number*/,
                                  double /*error_data_time*/,
                                  int /*tag_index*/,
                                  bool /*initial_time*/,
                                  bool /*uses_richardson_extrapolation_too*/)
{
    // intentionally blank
    return;
} // applyGradientDetector

void
IBStrategy::putToDatabase(SAMRAIPointer<SAMRAIDatabase> /*db*/)
{
    // intentionally blank
    return;
} // putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

INSHierarchyIntegrator*
IBStrategy::getINSHierarchyIntegrator() const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_ib_solver);
#endif
    return d_ib_solver->d_ins_hier_integrator;
} // getINSHierarchyIntegrator

SAMRAIPointer<SAMRAIHierarchyDataOpsReal<double>>
IBStrategy::getVelocityHierarchyDataOps() const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_ib_solver);
#endif
    return d_ib_solver->d_hier_velocity_data_ops;
} // getVelocityHierarchyDataOps

SAMRAIPointer<SAMRAIHierarchyDataOpsReal<double>>
IBStrategy::getPressureHierarchyDataOps() const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_ib_solver);
#endif
    return d_ib_solver->d_hier_pressure_data_ops;
} // getPressureHierarchyDataOps

SAMRAIPointer<HierarchyMathOps>
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
                             SAMRAIPointer<SAMRAIVariable> variable,
                             const SAMRAIIntVector& scratch_ghosts,
                             const std::string& coarsen_name,
                             const std::string& refine_name,
                             SAMRAIPointer<CartGridFunction> init_fcn,
                             const bool register_for_restart)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_ib_solver);
#endif
    d_ib_solver->registerVariable(current_idx,
                                  new_idx,
                                  scratch_idx,
                                  variable,
                                  scratch_ghosts,
                                  coarsen_name,
                                  refine_name,
                                  init_fcn,
                                  register_for_restart);
    return;
} // registerVariable

void
IBStrategy::registerVariable(int& idx,
                             SAMRAIPointer<SAMRAIVariable> variable,
                             const SAMRAIIntVector& ghosts,
                             SAMRAIPointer<SAMRAIVariableContext> ctx,
                             const bool register_for_restart)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_ib_solver);
#endif
    d_ib_solver->registerVariable(idx, variable, ghosts, ctx, register_for_restart);
    return;
} // registerVariable

void
IBStrategy::registerGhostfillRefineAlgorithm(const std::string& name,
                                             SAMRAIPointer<SAMRAIRefineAlgorithm> ghostfill_alg,
                                             std::unique_ptr<SAMRAIRefinePatchStrategy> ghostfill_patch_strategy)
{
    d_ib_solver->registerGhostfillRefineAlgorithm(name, ghostfill_alg, std::move(ghostfill_patch_strategy));
} // registerGhostfillRefineAlgorithm

void
IBStrategy::registerProlongRefineAlgorithm(const std::string& name,
                                           SAMRAIPointer<SAMRAIRefineAlgorithm> prolong_alg,
                                           std::unique_ptr<SAMRAIRefinePatchStrategy> prolong_patch_strategy)
{
    d_ib_solver->registerProlongRefineAlgorithm(name, prolong_alg, std::move(prolong_patch_strategy));
} // registerProlongRefineAlgorithm

void
IBStrategy::registerCoarsenAlgorithm(const std::string& name,
                                     SAMRAIPointer<SAMRAICoarsenAlgorithm> coarsen_alg,
                                     std::unique_ptr<SAMRAICoarsenPatchStrategy> coarsen_patch_strategy)
{
    d_ib_solver->registerCoarsenAlgorithm(name, coarsen_alg, std::move(coarsen_patch_strategy));
} // registerCoarsenAlgorithm

SAMRAIPointer<SAMRAIRefineAlgorithm>
IBStrategy::getGhostfillRefineAlgorithm(const std::string& name) const
{
    return d_ib_solver->getGhostfillRefineAlgorithm(name);
} // getGhostfillRefineAlgorithm

SAMRAIPointer<SAMRAIRefineAlgorithm>
IBStrategy::getProlongRefineAlgorithm(const std::string& name) const
{
    return d_ib_solver->getProlongRefineAlgorithm(name);
} // getProlongRefineAlgorithm

SAMRAIPointer<SAMRAICoarsenAlgorithm>
IBStrategy::getCoarsenAlgorithm(const std::string& name) const
{
    return d_ib_solver->getCoarsenAlgorithm(name);
} // getCoarsenAlgorithm

const std::vector<SAMRAIPointer<SAMRAIRefineSchedule>>&
IBStrategy::getGhostfillRefineSchedules(const std::string& name) const
{
    return d_ib_solver->getGhostfillRefineSchedules(name);
} // getGhostfillRefineSchedules

const std::vector<SAMRAIPointer<SAMRAIRefineSchedule>>&
IBStrategy::getProlongRefineSchedules(const std::string& name) const
{
    return d_ib_solver->getProlongRefineSchedules(name);
} // getProlongRefineSchedules

const std::vector<SAMRAIPointer<SAMRAICoarsenSchedule>>&
IBStrategy::getCoarsenSchedules(const std::string& name) const
{
    return d_ib_solver->getCoarsenSchedules(name);
} // getCoarsenSchedules

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
