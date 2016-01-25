// Filename: PoissonFACPreconditionerStrategy.cpp
// Created on 27 Apr 2012 by Boyce Griffith
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

#include "CartesianGridGeometry.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "CoarsenSchedule.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchyDataOpsReal.h"
#include "IntVector.h"
#include "LocationIndexRobinBcCoefs.h"
#include "MultiblockDataTranslator.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PoissonSpecifications.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefinePatchStrategy.h"
#include "RefineSchedule.h"
#include "RobinBcCoefStrategy.h"
#include "SAMRAIVectorReal.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "VariableFillPattern.h"
#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"
#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/FACPreconditionerStrategy.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/PoissonFACPreconditionerStrategy.h"
#include "ibtk/RefinePatchStrategySet.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_restrict_residual;
static Timer* t_prolong_error;
static Timer* t_prolong_error_and_correct;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

PoissonFACPreconditionerStrategy::PoissonFACPreconditionerStrategy(const std::string& object_name,
                                                                   Pointer<Variable<NDIM> > scratch_var,
                                                                   const int ghost_cell_width,
                                                                   const Pointer<Database> input_db,
                                                                   const std::string& default_options_prefix)
    : FACPreconditionerStrategy(object_name),
      d_poisson_spec(object_name + "::poisson_spec"),
      d_default_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_bc_coef", Pointer<Database>(NULL))),
      d_bc_coefs(1, d_default_bc_coef),
      d_gcw(ghost_cell_width),
      d_solution(NULL),
      d_rhs(NULL),
      d_hierarchy(),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_level_data_ops(),
      d_level_bdry_fill_ops(),
      d_level_math_ops(),
      d_in_initialize_operator_state(false),
      d_coarsest_reset_ln(-1),
      d_finest_reset_ln(-1),
      d_smoother_type("DEFAULT"),
      d_prolongation_method("DEFAULT"),
      d_restriction_method("DEFAULT"),
      d_coarse_solver_type("DEFAULT"),
      d_coarse_solver_default_options_prefix(default_options_prefix + "_coarse"),
      d_coarse_solver_rel_residual_tol(1.0e-5),
      d_coarse_solver_abs_residual_tol(1.0e-50),
      d_coarse_solver_max_iterations(10),
      d_context(NULL),
      d_bc_op(NULL),
      d_cf_bdry_op(),
      d_op_stencil_fill_pattern(),
      d_prolongation_refine_operator(),
      d_prolongation_refine_patch_strategy(),
      d_prolongation_refine_algorithm(),
      d_prolongation_refine_schedules(),
      d_restriction_coarsen_operator(),
      d_restriction_coarsen_algorithm(),
      d_restriction_coarsen_schedules(),
      d_ghostfill_nocoarse_refine_algorithm(),
      d_ghostfill_nocoarse_refine_schedules(),
      d_synch_refine_algorithm(),
      d_synch_refine_schedules()
{
    // Initialize the Poisson specifications.
    d_poisson_spec.setCZero();
    d_poisson_spec.setDConstant(-1.0);

    // Setup a default boundary condition object that specifies homogeneous
    // Dirichlet boundary conditions.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        LocationIndexRobinBcCoefs<NDIM>* p_default_bc_coef =
            dynamic_cast<LocationIndexRobinBcCoefs<NDIM>*>(d_default_bc_coef);
        p_default_bc_coef->setBoundaryValue(2 * d, 0.0);
        p_default_bc_coef->setBoundaryValue(2 * d + 1, 0.0);
    }

    // Get values from the input database.
    if (input_db)
    {
        if (input_db->keyExists("smoother_type")) d_smoother_type = input_db->getString("smoother_type");
        if (input_db->keyExists("prolongation_method"))
            d_prolongation_method = input_db->getString("prolongation_method");
        if (input_db->keyExists("restriction_method")) d_restriction_method = input_db->getString("restriction_method");
        if (input_db->keyExists("coarse_solver_type")) d_coarse_solver_type = input_db->getString("coarse_solver_type");
        if (input_db->keyExists("coarse_solver_rel_residual_tol"))
            d_coarse_solver_rel_residual_tol = input_db->getDouble("coarse_solver_rel_residual_tol");
        if (input_db->keyExists("coarse_solver_abs_residual_tol"))
            d_coarse_solver_abs_residual_tol = input_db->getDouble("coarse_solver_abs_residual_tol");
        if (input_db->keyExists("coarse_solver_max_iterations"))
            d_coarse_solver_max_iterations = input_db->getInteger("coarse_solver_max_iterations");
    }

    // Setup scratch variables.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(d_object_name + "::CONTEXT");
    const IntVector<NDIM> ghosts = d_gcw;
    if (var_db->checkVariableExists(scratch_var->getName()))
    {
        scratch_var = var_db->getVariable(scratch_var->getName());
        d_scratch_idx = var_db->mapVariableAndContextToIndex(scratch_var, d_context);
        var_db->removePatchDataIndex(d_scratch_idx);
    }
    d_scratch_idx = var_db->registerVariableAndContext(scratch_var, d_context, ghosts);

    // Setup Timers.
    IBTK_DO_ONCE(
        t_restrict_residual =
            TimerManager::getManager()->getTimer("IBTK::PoissonFACPreconditionerStrategy::restrictResidual()");
        t_prolong_error =
            TimerManager::getManager()->getTimer("IBTK::PoissonFACPreconditionerStrategy::prolongError()");
        t_prolong_error_and_correct =
            TimerManager::getManager()->getTimer("IBTK::PoissonFACPreconditionerStrategy::prolongErrorAndCorrect()");
        t_initialize_operator_state =
            TimerManager::getManager()->getTimer("IBTK::PoissonFACPreconditionerStrategy::initializeOperatorState()");
        t_deallocate_operator_state =
            TimerManager::getManager()->getTimer("IBTK::PoissonFACPreconditionerStrategy::deallocateOperatorState()"););
    return;
} // PoissonFACPreconditionerStrategy

PoissonFACPreconditionerStrategy::~PoissonFACPreconditionerStrategy()
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::~PoissonFACPreconditionerStrategy()\n"
                                 << "  subclass must call deallocateOperatorState in subclass destructor"
                                 << std::endl);
    }
    delete d_default_bc_coef;
    d_default_bc_coef = NULL;
    return;
} // ~PoissonFACPreconditionerStrategy

void
PoissonFACPreconditionerStrategy::setPoissonSpecifications(const PoissonSpecifications& poisson_spec)
{
    d_poisson_spec = poisson_spec;
    return;
} // setPoissonSpecifications

void
PoissonFACPreconditionerStrategy::setPhysicalBcCoef(RobinBcCoefStrategy<NDIM>* const bc_coef)
{
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(1, bc_coef));
    return;
} // setPhysicalBcCoef

void
PoissonFACPreconditionerStrategy::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
    d_bc_coefs.resize(bc_coefs.size());
    for (unsigned int l = 0; l < bc_coefs.size(); ++l)
    {
        if (bc_coefs[l])
        {
            d_bc_coefs[l] = bc_coefs[l];
        }
        else
        {
            d_bc_coefs[l] = d_default_bc_coef;
        }
    }
    return;
} // setPhysicalBcCoefs

void
PoissonFACPreconditionerStrategy::setResetLevels(const int coarsest_ln, const int finest_ln)
{
#if !defined(NDEBUG)
    TBOX_ASSERT((coarsest_ln == -1 && finest_ln == -1) || (coarsest_ln >= 0 && finest_ln >= coarsest_ln));
#endif
    if (d_is_initialized)
    {
        d_coarsest_reset_ln = coarsest_ln;
        d_finest_reset_ln = finest_ln;
    }
    return;
} // setResetLevels

void
PoissonFACPreconditionerStrategy::setCoarseSolverMaxIterations(int coarse_solver_max_iterations)
{
    d_coarse_solver_max_iterations = coarse_solver_max_iterations;
    return;
} // setCoarseSolverMaxIterations

void
PoissonFACPreconditionerStrategy::setCoarseSolverAbsoluteTolerance(double coarse_solver_abs_residual_tol)
{
    d_coarse_solver_abs_residual_tol = coarse_solver_abs_residual_tol;
    return;
} // setCoarseSolverAbsoluteTolerance

void
PoissonFACPreconditionerStrategy::setCoarseSolverRelativeTolerance(double coarse_solver_rel_residual_tol)
{
    d_coarse_solver_rel_residual_tol = coarse_solver_rel_residual_tol;
    return;
} // setCoarseSolverRelativeTolerance

void
PoissonFACPreconditionerStrategy::setProlongationMethod(const std::string& prolongation_method)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setProlongationMethod()\n"
                                 << "  cannot be called while operator state is initialized"
                                 << std::endl);
    }
    d_prolongation_method = prolongation_method;
    return;
} // setProlongationMethod

void
PoissonFACPreconditionerStrategy::setRestrictionMethod(const std::string& restriction_method)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setRestrictionMethod()\n"
                                 << "  cannot be called while operator state is initialized"
                                 << std::endl);
    }
    d_restriction_method = restriction_method;
    return;
} // setRestrictionMethod

void
PoissonFACPreconditionerStrategy::restrictResidual(const SAMRAIVectorReal<NDIM, double>& src,
                                                   SAMRAIVectorReal<NDIM, double>& dst,
                                                   int dst_ln)
{
    IBTK_TIMER_START(t_restrict_residual);

    const int src_idx = src.getComponentDescriptorIndex(0);
    const int dst_idx = dst.getComponentDescriptorIndex(0);

    if (src_idx != dst_idx)
    {
        d_level_data_ops[dst_ln]->copyData(dst_idx, src_idx, /*interior_only*/ false);
    }
    xeqScheduleRestriction(dst_idx, src_idx, dst_ln);

    IBTK_TIMER_STOP(t_restrict_residual);
    return;
} // restrictResidual

void
PoissonFACPreconditionerStrategy::prolongError(const SAMRAIVectorReal<NDIM, double>& src,
                                               SAMRAIVectorReal<NDIM, double>& dst,
                                               int dst_ln)
{
    IBTK_TIMER_START(t_prolong_error);

    const int dst_idx = dst.getComponentDescriptorIndex(0);
    const int src_idx = src.getComponentDescriptorIndex(0);

    // Prolong the correction from the coarse src level data into the fine level
    // dst data.
    xeqScheduleProlongation(dst_idx, src_idx, dst_ln);

    IBTK_TIMER_STOP(t_prolong_error);
    return;
} // prolongError

void
PoissonFACPreconditionerStrategy::prolongErrorAndCorrect(const SAMRAIVectorReal<NDIM, double>& src,
                                                         SAMRAIVectorReal<NDIM, double>& dst,
                                                         int dst_ln)
{
    IBTK_TIMER_START(t_prolong_error_and_correct);

    const int dst_idx = dst.getComponentDescriptorIndex(0);
    const int src_idx = src.getComponentDescriptorIndex(0);

    // Prolong the correction from the coarse level src data into the fine level
    // scratch data and then correct the fine level dst data.
    if (src_idx != dst_idx)
    {
        d_level_data_ops[dst_ln - 1]->add(dst_idx, dst_idx, src_idx, /*interior_only*/ false);
    }
    xeqScheduleProlongation(d_scratch_idx, src_idx, dst_ln);
    d_level_data_ops[dst_ln]->add(dst_idx, dst_idx, d_scratch_idx, /*interior_only*/ false);

    IBTK_TIMER_STOP(t_prolong_error_and_correct);
    return;
} // prolongErrorAndCorrect

void
PoissonFACPreconditionerStrategy::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& solution,
                                                          const SAMRAIVectorReal<NDIM, double>& rhs)
{
    IBTK_TIMER_START(t_initialize_operator_state);

    d_in_initialize_operator_state = true;

    // Cache the level range to be reset.
    //
    // NOTE: We cannot use d_coarsest_reset_ln and d_finest_reset_ln since those
    // values are reset by deallocateOperatorState().
    const int coarsest_reset_ln =
        (d_coarsest_reset_ln != -1 && d_finest_reset_ln != -1 ? d_coarsest_reset_ln :
                                                                solution.getCoarsestLevelNumber());
    const int finest_reset_ln =
        (d_coarsest_reset_ln != -1 && d_finest_reset_ln != -1 ? d_finest_reset_ln : solution.getFinestLevelNumber());

    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Setup solution and rhs vectors.
    d_solution = solution.cloneVector(solution.getName());
    d_rhs = rhs.cloneVector(rhs.getName());

    Pointer<Variable<NDIM> > sol_var = d_solution->getComponentVariable(0);
    const int sol_idx = d_solution->getComponentDescriptorIndex(0);
    const int rhs_idx = d_rhs->getComponentDescriptorIndex(0);

    // Reset the hierarchy configuration.
    d_hierarchy = solution.getPatchHierarchy();
    d_coarsest_ln = solution.getCoarsestLevelNumber();
    d_finest_ln = solution.getFinestLevelNumber();

    // Perform implementation-specific initialization.
    initializeOperatorStateSpecialized(solution, rhs, coarsest_reset_ln, finest_reset_ln);
#if !defined(NDEBUG)
    TBOX_ASSERT(d_bc_op);
    TBOX_ASSERT(d_cf_bdry_op);
#endif

    // Setup level operators.
    d_level_data_ops.resize(d_finest_ln + 1);
    d_level_bdry_fill_ops.resize(d_finest_ln + 1, NULL);
    d_level_math_ops.resize(d_finest_ln + 1, NULL);
    HierarchyDataOpsManager<NDIM>* hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln); ln <= finest_reset_ln; ++ln)
    {
        d_level_data_ops[ln] = hier_data_ops_manager->getOperationsDouble(sol_var,
                                                                          d_hierarchy,
                                                                          /*get_unique*/ true);
        d_level_data_ops[ln]->resetLevels(ln, ln);
        d_level_bdry_fill_ops[ln].setNull();
        d_level_math_ops[ln].setNull();
    }

    // Allocate scratch data.
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln); ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_scratch_idx)) level->allocatePatchData(d_scratch_idx);
    }

    // Get the transfer operators.
    Pointer<CartesianGridGeometry<NDIM> > geometry = d_hierarchy->getGridGeometry();
    d_prolongation_refine_operator = geometry->lookupRefineOperator(sol_var, d_prolongation_method);
    d_restriction_coarsen_operator = geometry->lookupCoarsenOperator(sol_var, d_restriction_method);
    d_cf_bdry_op->setConsistentInterpolationScheme(false);
    d_cf_bdry_op->setPatchDataIndex(d_scratch_idx);
    d_cf_bdry_op->setPatchHierarchy(d_hierarchy);

    // Make space for saving communication schedules.  There is no need to
    // delete the old schedules first because we have deallocated the solver
    // state above.
    std::vector<RefinePatchStrategy<NDIM>*> prolongation_refine_patch_strategies;
    prolongation_refine_patch_strategies.push_back(d_cf_bdry_op);
    prolongation_refine_patch_strategies.push_back(d_bc_op);
    d_prolongation_refine_patch_strategy = new RefinePatchStrategySet(
        prolongation_refine_patch_strategies.begin(), prolongation_refine_patch_strategies.end(), false);

    d_prolongation_refine_schedules.resize(d_finest_ln + 1);
    d_restriction_coarsen_schedules.resize(d_finest_ln);
    d_ghostfill_nocoarse_refine_schedules.resize(d_finest_ln + 1);
    d_synch_refine_schedules.resize(d_finest_ln + 1);

    d_prolongation_refine_algorithm = new RefineAlgorithm<NDIM>();
    d_restriction_coarsen_algorithm = new CoarsenAlgorithm<NDIM>();
    d_ghostfill_nocoarse_refine_algorithm = new RefineAlgorithm<NDIM>();
    d_synch_refine_algorithm = new RefineAlgorithm<NDIM>();

    d_prolongation_refine_algorithm->registerRefine(
        d_scratch_idx, sol_idx, d_scratch_idx, d_prolongation_refine_operator, d_op_stencil_fill_pattern);
    d_restriction_coarsen_algorithm->registerCoarsen(d_scratch_idx, rhs_idx, d_restriction_coarsen_operator);
    d_ghostfill_nocoarse_refine_algorithm->registerRefine(
        sol_idx, sol_idx, sol_idx, Pointer<RefineOperator<NDIM> >(), d_op_stencil_fill_pattern);
    d_synch_refine_algorithm->registerRefine(
        sol_idx, sol_idx, sol_idx, Pointer<RefineOperator<NDIM> >(), d_synch_fill_pattern);

    // TODO: Here we take a pessimistic approach and are recreating refine schedule for
    // (coarsest_reset_ln - 1) level as well.
    for (int dst_ln = std::max(d_coarsest_ln + 1, coarsest_reset_ln - 1); dst_ln <= finest_reset_ln; ++dst_ln)
    {
        d_prolongation_refine_schedules[dst_ln] =
            d_prolongation_refine_algorithm->createSchedule(d_hierarchy->getPatchLevel(dst_ln),
                                                            Pointer<PatchLevel<NDIM> >(),
                                                            dst_ln - 1,
                                                            d_hierarchy,
                                                            d_prolongation_refine_patch_strategy.getPointer());
    }

    for (int dst_ln = coarsest_reset_ln; dst_ln < std::min(finest_reset_ln + 1, d_finest_ln); ++dst_ln)
    {
        d_restriction_coarsen_schedules[dst_ln] = d_restriction_coarsen_algorithm->createSchedule(
            d_hierarchy->getPatchLevel(dst_ln), d_hierarchy->getPatchLevel(dst_ln + 1));
    }

    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        d_ghostfill_nocoarse_refine_schedules[ln] =
            d_ghostfill_nocoarse_refine_algorithm->createSchedule(d_hierarchy->getPatchLevel(ln), d_bc_op.getPointer());
        d_synch_refine_schedules[ln] = d_synch_refine_algorithm->createSchedule(d_hierarchy->getPatchLevel(ln));
    }

    // Indicate that the operator is initialized.
    d_is_initialized = true;
    d_in_initialize_operator_state = false;

    IBTK_TIMER_STOP(t_initialize_operator_state);
    return;
} // initializeOperatorState

void
PoissonFACPreconditionerStrategy::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_deallocate_operator_state);

    const int coarsest_reset_ln =
        (d_in_initialize_operator_state && (d_coarsest_reset_ln != -1) && (d_finest_reset_ln != -1)) ?
            d_coarsest_reset_ln :
            d_coarsest_ln;
    const int finest_reset_ln =
        (d_in_initialize_operator_state && (d_coarsest_reset_ln != -1) && (d_finest_reset_ln != -1)) ?
            d_finest_reset_ln :
            d_finest_ln;
    deallocateOperatorStateSpecialized(coarsest_reset_ln, finest_reset_ln);

    // Deallocate scratch data.
    for (int ln = coarsest_reset_ln; ln <= std::min(d_finest_ln, finest_reset_ln); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_scratch_idx)) level->deallocatePatchData(d_scratch_idx);
    }

    // Delete the solution and rhs vectors.
    d_solution->freeVectorComponents();
    d_solution.setNull();

    d_rhs->freeVectorComponents();
    d_rhs.setNull();

    // Only fully deallocate operator data when we are not reinitializing the
    // operator.
    if (!d_in_initialize_operator_state)
    {
        d_hierarchy.setNull();
        d_coarsest_ln = -1;
        d_finest_ln = -1;

        d_level_bdry_fill_ops.clear();
        d_level_math_ops.clear();

        d_prolongation_refine_operator.setNull();
        d_prolongation_refine_patch_strategy.setNull();
        d_prolongation_refine_algorithm.setNull();
        d_prolongation_refine_schedules.resize(0);

        d_restriction_coarsen_operator.setNull();
        d_restriction_coarsen_algorithm.setNull();
        d_restriction_coarsen_schedules.resize(0);

        d_ghostfill_nocoarse_refine_algorithm.setNull();
        d_ghostfill_nocoarse_refine_schedules.resize(0);

        d_synch_refine_algorithm.setNull();
        d_synch_refine_schedules.resize(0);
    }

    // Clear the "reset level" range.
    d_coarsest_reset_ln = -1;
    d_finest_reset_ln = -1;

    // Indicate that the operator is not initialized.
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_operator_state);
    return;
} // deallocateOperatorState

/////////////////////////////// PROTECTED ////////////////////////////////////

void
PoissonFACPreconditionerStrategy::xeqScheduleProlongation(const int dst_idx, const int src_idx, const int dst_ln)
{
    d_cf_bdry_op->setPatchDataIndex(dst_idx);
    d_bc_op->setPatchDataIndex(dst_idx);
    d_bc_op->setPhysicalBcCoefs(d_bc_coefs);
    d_bc_op->setHomogeneousBc(true);
    for (unsigned int k = 0; k < d_bc_coefs.size(); ++k)
    {
        ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[k]);
        if (extended_bc_coef)
        {
            extended_bc_coef->setTargetPatchDataIndex(dst_idx);
            extended_bc_coef->setHomogeneousBc(true);
        }
    }
    RefineAlgorithm<NDIM> refiner;
    refiner.registerRefine(dst_idx, src_idx, dst_idx, d_prolongation_refine_operator, d_op_stencil_fill_pattern);
    refiner.resetSchedule(d_prolongation_refine_schedules[dst_ln]);
    d_prolongation_refine_schedules[dst_ln]->fillData(d_solution_time);
    d_prolongation_refine_algorithm->resetSchedule(d_prolongation_refine_schedules[dst_ln]);
    for (unsigned int k = 0; k < d_bc_coefs.size(); ++k)
    {
        ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[k]);
        if (extended_bc_coef) extended_bc_coef->clearTargetPatchDataIndex();
    }
    return;
} // xeqScheduleProlongation

void
PoissonFACPreconditionerStrategy::xeqScheduleRestriction(const int dst_idx, const int src_idx, const int dst_ln)
{
    CoarsenAlgorithm<NDIM> coarsener;
    coarsener.registerCoarsen(dst_idx, src_idx, d_restriction_coarsen_operator);
    coarsener.resetSchedule(d_restriction_coarsen_schedules[dst_ln]);
    d_restriction_coarsen_schedules[dst_ln]->coarsenData();
    d_restriction_coarsen_algorithm->resetSchedule(d_restriction_coarsen_schedules[dst_ln]);
    return;
} // xeqScheduleRestriction

void
PoissonFACPreconditionerStrategy::xeqScheduleGhostFillNoCoarse(const int dst_idx, const int dst_ln)
{
    d_bc_op->setPatchDataIndex(dst_idx);
    d_bc_op->setPhysicalBcCoefs(d_bc_coefs);
    d_bc_op->setHomogeneousBc(true);
    for (unsigned int k = 0; k < d_bc_coefs.size(); ++k)
    {
        ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[k]);
        if (extended_bc_coef)
        {
            extended_bc_coef->setTargetPatchDataIndex(dst_idx);
            extended_bc_coef->setHomogeneousBc(true);
        }
    }
    RefineAlgorithm<NDIM> refiner;
    refiner.registerRefine(dst_idx, dst_idx, dst_idx, Pointer<RefineOperator<NDIM> >(), d_op_stencil_fill_pattern);
    refiner.resetSchedule(d_ghostfill_nocoarse_refine_schedules[dst_ln]);
    d_ghostfill_nocoarse_refine_schedules[dst_ln]->fillData(d_solution_time);
    d_ghostfill_nocoarse_refine_algorithm->resetSchedule(d_ghostfill_nocoarse_refine_schedules[dst_ln]);
    for (unsigned int k = 0; k < d_bc_coefs.size(); ++k)
    {
        ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[k]);
        if (extended_bc_coef) extended_bc_coef->clearTargetPatchDataIndex();
    }
    return;
} // xeqScheduleGhostFillNoCoarse

void
PoissonFACPreconditionerStrategy::xeqScheduleDataSynch(const int dst_idx, const int dst_ln)
{
    RefineAlgorithm<NDIM> refiner;
    refiner.registerRefine(dst_idx, dst_idx, dst_idx, Pointer<RefineOperator<NDIM> >(), d_synch_fill_pattern);
    refiner.resetSchedule(d_synch_refine_schedules[dst_ln]);
    d_synch_refine_schedules[dst_ln]->fillData(d_solution_time);
    d_synch_refine_algorithm->resetSchedule(d_synch_refine_schedules[dst_ln]);
    return;
} // xeqScheduleDataSynch

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
