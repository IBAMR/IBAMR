// Filename: INSStaggeredFACPreconditionerStrategy.C
// Created on 18 Apr 2012 by Boyce Griffith
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

#include "INSStaggeredFACPreconditionerStrategy.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/CartCellDoubleCubicCoarsen.h>
#include <ibtk/CartSideDoubleCubicCoarsen.h>
#include <ibtk/RefinePatchStrategySet.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
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

INSStaggeredFACPreconditionerStrategy::INSStaggeredFACPreconditionerStrategy(
    const std::string& object_name,
    const int ghost_cell_width,
    const Pointer<Database> input_db)
    : d_object_name(object_name),
      d_is_initialized(false),
      d_current_time(0.0),
      d_new_time(0.0),
      d_gcw(ghost_cell_width),
      d_solution(NULL),
      d_rhs(NULL),
      d_hierarchy(),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_hier_bdry_fill_ops(),
      d_hier_math_ops(),
      d_in_initialize_operator_state(false),
      d_coarsest_reset_ln(-1),
      d_finest_reset_ln(-1),
      d_smoother_choice("additive"),
      d_U_prolongation_method("CONSTANT_REFINE"),
      d_P_prolongation_method("LINEAR_REFINE"),
      d_U_restriction_method("CONSERVATIVE_COARSEN"),
      d_P_restriction_method("CONSERVATIVE_COARSEN"),
      d_coarse_solver_choice("block_jacobi"),
      d_coarse_solver_tol(1.0e-6),
      d_coarse_solver_max_its(10),
      d_context(NULL),
      d_side_scratch_idx(-1),
      d_cell_scratch_idx(-1),
      d_U_cf_bdry_op(),
      d_P_cf_bdry_op(),
      d_U_op_stencil_fill_pattern(),
      d_P_op_stencil_fill_pattern(),
      d_U_synch_fill_pattern(),
      d_U_prolongation_refine_operator(),
      d_P_prolongation_refine_operator(),
      d_prolongation_refine_patch_strategy(),
      d_prolongation_refine_algorithm(),
      d_prolongation_refine_schedules(),
      d_U_restriction_coarsen_operator(),
      d_P_restriction_coarsen_operator(),
      d_restriction_coarsen_algorithm(),
      d_restriction_coarsen_schedules(),
      d_ghostfill_nocoarse_refine_algorithm(),
      d_ghostfill_nocoarse_refine_schedules(),
      d_synch_refine_algorithm(),
      d_synch_refine_schedules()
{
    // Get values from the input database.
    if (!input_db.isNull())
    {
        d_smoother_choice = input_db->getStringWithDefault("smoother_choice", d_smoother_choice);
        d_U_prolongation_method = input_db->getStringWithDefault("U_prolongation_method", d_U_prolongation_method);
        d_P_prolongation_method = input_db->getStringWithDefault("P_prolongation_method", d_P_prolongation_method);
        d_U_restriction_method = input_db->getStringWithDefault("U_restriction_method", d_U_restriction_method);
        d_P_restriction_method = input_db->getStringWithDefault("P_restriction_method", d_P_restriction_method);
        d_coarse_solver_choice = input_db->getStringWithDefault("coarse_solver_choice", d_coarse_solver_choice);
        d_coarse_solver_tol = input_db->getDoubleWithDefault("coarse_solver_tolerance", d_coarse_solver_tol);
        d_coarse_solver_max_its = input_db->getIntegerWithDefault("coarse_solver_max_iterations", d_coarse_solver_max_its);
    }

    // Setup scratch variables.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(d_object_name+"::CONTEXT");
    const IntVector<NDIM> side_ghosts = d_gcw;
    Pointer<SideVariable<NDIM,double> > side_scratch_var = new SideVariable<NDIM,double>(d_object_name+"::side_scratch");
    if (var_db->checkVariableExists(side_scratch_var->getName()))
    {
        side_scratch_var = var_db->getVariable(side_scratch_var->getName());
        d_side_scratch_idx = var_db->mapVariableAndContextToIndex(side_scratch_var, d_context);
        var_db->removePatchDataIndex(d_side_scratch_idx);
    }
    d_side_scratch_idx = var_db->registerVariableAndContext(side_scratch_var, d_context, side_ghosts);
    const IntVector<NDIM> cell_ghosts = d_gcw;
    Pointer<CellVariable<NDIM,double> > cell_scratch_var = new CellVariable<NDIM,double>(d_object_name+"::cell_scratch");
    if (var_db->checkVariableExists(cell_scratch_var->getName()))
    {
        cell_scratch_var = var_db->getVariable(cell_scratch_var->getName());
        d_cell_scratch_idx = var_db->mapVariableAndContextToIndex(cell_scratch_var, d_context);
        var_db->removePatchDataIndex(d_cell_scratch_idx);
    }
    d_cell_scratch_idx = var_db->registerVariableAndContext(cell_scratch_var, d_context, cell_ghosts);

    // Setup Timers.
    IBAMR_DO_ONCE(
        t_restrict_residual         = TimerManager::getManager()->getTimer("INSStaggeredFACPreconditionerStrategy::restrictResidual()");
        t_prolong_error             = TimerManager::getManager()->getTimer("INSStaggeredFACPreconditionerStrategy::prolongError()");
        t_prolong_error_and_correct = TimerManager::getManager()->getTimer("INSStaggeredFACPreconditionerStrategy::prolongErrorAndCorrect()");
        t_initialize_operator_state = TimerManager::getManager()->getTimer("INSStaggeredFACPreconditionerStrategy::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer("INSStaggeredFACPreconditionerStrategy::deallocateOperatorState()");
                  );
    return;
}// INSStaggeredFACPreconditionerStrategy

INSStaggeredFACPreconditionerStrategy::~INSStaggeredFACPreconditionerStrategy()
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::~INSStaggeredFACPreconditionerStrategy()\n"
                   << "  subclass must call deallocateOperatorState in subclass destructor" << std::endl);
    }
    return;
}// ~INSStaggeredFACPreconditionerStrategy

void
INSStaggeredFACPreconditionerStrategy::setTimeInterval(
    const double current_time,
    const double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    return;
}// setTimeInterval

void
INSStaggeredFACPreconditionerStrategy::setResetLevels(
    const int coarsest_ln,
    const int finest_ln)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((coarsest_ln == -1 && finest_ln == -1) ||
                (coarsest_ln >=  0 && finest_ln >= coarsest_ln));
#endif
    if (d_is_initialized)
    {
        d_coarsest_reset_ln = coarsest_ln;
        d_finest_reset_ln = finest_ln;
    }
    return;
}// setResetLevels

void
INSStaggeredFACPreconditionerStrategy::setSmootherChoice(
    const std::string& smoother_choice)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setSmootherChoice()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_smoother_choice = smoother_choice;
    return;
}// setSmootherChoice

void
INSStaggeredFACPreconditionerStrategy::setCoarsestLevelSolverChoice(
    const std::string& coarse_solver_choice)
{
    d_coarse_solver_choice = coarse_solver_choice;
    return;
}// setCoarsestLevelSolverChoice

void
INSStaggeredFACPreconditionerStrategy::setCoarsestLevelSolverTolerance(
    double coarse_solver_tol)
{
    d_coarse_solver_tol = coarse_solver_tol;
    return;
}// setCoarsestLevelSolverTolerance

void
INSStaggeredFACPreconditionerStrategy::setCoarsestLevelSolverMaxIterations(
    int coarse_solver_max_its)
{
    d_coarse_solver_max_its = coarse_solver_max_its;
    return;
}// setCoarsestLevelSolverMaxIterations

void
INSStaggeredFACPreconditionerStrategy::setProlongationMethods(
    const std::string& U_prolongation_method,
    const std::string& P_prolongation_method)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setProlongationMethods()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_U_prolongation_method = U_prolongation_method;
    d_P_prolongation_method = P_prolongation_method;
    return;
}// setProlongationMethods

void
INSStaggeredFACPreconditionerStrategy::setRestrictionMethods(
    const std::string& U_restriction_method,
    const std::string& P_restriction_method)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setRestrictionMethods()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_U_restriction_method = U_restriction_method;
    d_P_restriction_method = P_restriction_method;
    return;
}// setRestrictionMethods

void
INSStaggeredFACPreconditionerStrategy::restrictResidual(
    const SAMRAIVectorReal<NDIM,double>& src,
    SAMRAIVectorReal<NDIM,double>& dst,
    int dst_ln)
{
    IBAMR_TIMER_START(t_restrict_residual);

    const int U_src_idx = src.getComponentDescriptorIndex(0);
    const int P_src_idx = src.getComponentDescriptorIndex(1);
    const std::pair<int,int> src_idxs = std::make_pair(U_src_idx,P_src_idx);

    const int U_dst_idx = dst.getComponentDescriptorIndex(0);
    const int P_dst_idx = dst.getComponentDescriptorIndex(1);
    const std::pair<int,int> dst_idxs = std::make_pair(U_dst_idx,P_dst_idx);

    if (U_src_idx != U_dst_idx)
    {
        HierarchySideDataOpsReal<NDIM,double> hier_sc_data_ops(d_hierarchy, dst_ln, dst_ln);
        static const bool interior_only = false;
        hier_sc_data_ops.copyData(U_dst_idx, U_src_idx, interior_only);
    }

    if (P_src_idx != P_dst_idx)
    {
        HierarchyCellDataOpsReal<NDIM,double> hier_sc_data_ops(d_hierarchy, dst_ln, dst_ln);
        static const bool interior_only = false;
        hier_sc_data_ops.copyData(P_dst_idx, P_src_idx, interior_only);
    }

    xeqScheduleRestriction(dst_idxs, src_idxs, dst_ln);

    IBAMR_TIMER_STOP(t_restrict_residual);
    return;
}// restrictResidual

void
INSStaggeredFACPreconditionerStrategy::prolongError(
    const SAMRAIVectorReal<NDIM,double>& src,
    SAMRAIVectorReal<NDIM,double>& dst,
    int dst_ln)
{
    IBAMR_TIMER_START(t_prolong_error);

    const int U_src_idx = src.getComponentDescriptorIndex(0);
    const int P_src_idx = src.getComponentDescriptorIndex(1);
    const std::pair<int,int> src_idxs = std::make_pair(U_src_idx,P_src_idx);

    const int U_dst_idx = dst.getComponentDescriptorIndex(0);
    const int P_dst_idx = dst.getComponentDescriptorIndex(1);
    const std::pair<int,int> dst_idxs = std::make_pair(U_dst_idx,P_dst_idx);

    // Refine the correction from the coarse level src data directly into the
    // fine level error.
    xeqScheduleProlongation(dst_idxs, src_idxs, dst_ln);

    IBAMR_TIMER_STOP(t_prolong_error);
    return;
}// prolongError

void
INSStaggeredFACPreconditionerStrategy::prolongErrorAndCorrect(
    const SAMRAIVectorReal<NDIM,double>& src,
    SAMRAIVectorReal<NDIM,double>& dst,
    int dst_ln)
{
    IBAMR_TIMER_START(t_prolong_error_and_correct);

    const int U_src_idx = src.getComponentDescriptorIndex(0);
    const int P_src_idx = src.getComponentDescriptorIndex(1);
    const std::pair<int,int> src_idxs = std::make_pair(U_src_idx,P_src_idx);

    const int U_dst_idx = dst.getComponentDescriptorIndex(0);
    const int P_dst_idx = dst.getComponentDescriptorIndex(1);

    const std::pair<int,int> scratch_idxs = std::make_pair(d_side_scratch_idx,d_cell_scratch_idx);

    // Prolong the correction from the coarse level src data into the fine level
    // scratch data and then correct the fine level dst data.
    static const bool interior_only = false;
    if (U_src_idx != U_dst_idx)
    {
        HierarchySideDataOpsReal<NDIM,double> hier_sc_data_ops_coarse(d_hierarchy, dst_ln-1, dst_ln-1);
        hier_sc_data_ops_coarse.add(U_dst_idx, U_dst_idx, U_src_idx, interior_only);
    }
    if (P_src_idx != P_dst_idx)
    {
        HierarchyCellDataOpsReal<NDIM,double> hier_sc_data_ops_coarse(d_hierarchy, dst_ln-1, dst_ln-1);
        hier_sc_data_ops_coarse.add(P_dst_idx, P_dst_idx, P_src_idx, interior_only);
    }
    xeqScheduleProlongation(scratch_idxs, src_idxs, dst_ln);
    HierarchySideDataOpsReal<NDIM,double> hier_sc_data_ops_fine(d_hierarchy, dst_ln, dst_ln);
    hier_sc_data_ops_fine.add(U_dst_idx, U_dst_idx, d_side_scratch_idx, interior_only);
    HierarchyCellDataOpsReal<NDIM,double> hier_cc_data_ops_fine(d_hierarchy, dst_ln, dst_ln);
    hier_cc_data_ops_fine.add(P_dst_idx, P_dst_idx, d_cell_scratch_idx, interior_only);

    IBAMR_TIMER_STOP(t_prolong_error_and_correct);
    return;
}// prolongErrorAndCorrect

void
INSStaggeredFACPreconditionerStrategy::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& solution,
    const SAMRAIVectorReal<NDIM,double>& rhs)
{
    IBAMR_TIMER_START(t_initialize_operator_state);

    d_in_initialize_operator_state = true;

    // Cache the level range to be reset.
    //
    // NOTE: We cannot use d_coarsest_reset_ln and d_finest_reset_ln since those
    // values are reset by deallocateOperatorState().
    const int coarsest_reset_ln =
        (d_coarsest_reset_ln != -1 && d_finest_reset_ln != -1
         ? d_coarsest_reset_ln
         : solution.getCoarsestLevelNumber());
    const int finest_reset_ln =
        (d_coarsest_reset_ln != -1 && d_finest_reset_ln != -1
         ? d_finest_reset_ln
         : solution.getFinestLevelNumber());

    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Setup solution and rhs vectors.
    d_solution = solution.cloneVector(solution.getName());
    d_solution->allocateVectorData();

    d_rhs = rhs.cloneVector(rhs.getName());
    d_rhs->allocateVectorData();

    // Reset the hierarchy configuration.
    d_hierarchy   = solution.getPatchHierarchy();
    d_coarsest_ln = solution.getCoarsestLevelNumber();
    d_finest_ln   = solution.getFinestLevelNumber();

    // Perform implementation-specific initialization.
    initializeOperatorStateSpecialized(solution, rhs, coarsest_reset_ln, finest_reset_ln);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_U_bc_op.isNull());
    TBOX_ASSERT(!d_P_bc_op.isNull());
    TBOX_ASSERT(!d_U_cf_bdry_op.isNull());
    TBOX_ASSERT(!d_P_cf_bdry_op.isNull());
    TBOX_ASSERT(!d_U_op_stencil_fill_pattern.isNull());
    TBOX_ASSERT(!d_P_op_stencil_fill_pattern.isNull());
    TBOX_ASSERT(!d_U_synch_fill_pattern.isNull());
#endif

    // Setup level operators.
    d_hier_bdry_fill_ops.resize(d_finest_ln+1, NULL);
    d_hier_math_ops.resize(d_finest_ln+1, NULL);
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln); ln <= finest_reset_ln; ++ln)
    {
        d_hier_bdry_fill_ops[ln].setNull();
        d_hier_math_ops[ln].setNull();
    }

    // Allocate scratch data.
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln); ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_side_scratch_idx)) level->allocatePatchData(d_side_scratch_idx);
        if (!level->checkAllocated(d_cell_scratch_idx)) level->allocatePatchData(d_cell_scratch_idx);
    }

    // Get the transfer operators.
    Pointer<CartesianGridGeometry<NDIM> > geometry = d_hierarchy->getGridGeometry();
    IBAMR_DO_ONCE(
        geometry->addSpatialCoarsenOperator(new CartSideDoubleCubicCoarsen());
        geometry->addSpatialCoarsenOperator(new CartCellDoubleCubicCoarsen());
                  );
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > var;

    var_db->mapIndexToVariable(d_side_scratch_idx, var);
    d_U_prolongation_refine_operator = geometry->lookupRefineOperator(var, d_U_prolongation_method);

    d_U_cf_bdry_op->setConsistentInterpolationScheme(false);
    d_U_cf_bdry_op->setPatchDataIndex(d_side_scratch_idx);
    d_U_cf_bdry_op->setPatchHierarchy(d_hierarchy);

    var_db->mapIndexToVariable(d_cell_scratch_idx, var);
    d_P_prolongation_refine_operator = geometry->lookupRefineOperator(var, d_P_prolongation_method);

    d_P_cf_bdry_op->setConsistentInterpolationScheme(false);
    d_P_cf_bdry_op->setPatchDataIndex(d_cell_scratch_idx);
    d_P_cf_bdry_op->setPatchHierarchy(d_hierarchy);

    var_db->mapIndexToVariable(d_side_scratch_idx, var);
    d_U_restriction_coarsen_operator = geometry->lookupCoarsenOperator(var, d_U_restriction_method);

    var_db->mapIndexToVariable(d_cell_scratch_idx, var);
    d_P_restriction_coarsen_operator = geometry->lookupCoarsenOperator(var, d_P_restriction_method);

    // Make space for saving communication schedules.  There is no need to
    // delete the old schedules first because we have deallocated the solver
    // state above.
    std::vector<RefinePatchStrategy<NDIM>*> prolongation_refine_patch_strategies;
    prolongation_refine_patch_strategies.push_back(d_U_cf_bdry_op);
    prolongation_refine_patch_strategies.push_back(d_P_cf_bdry_op);
    prolongation_refine_patch_strategies.push_back(d_U_bc_op);
    prolongation_refine_patch_strategies.push_back(d_P_bc_op);
    d_prolongation_refine_patch_strategy = new RefinePatchStrategySet(
        prolongation_refine_patch_strategies.begin(), prolongation_refine_patch_strategies.end(), false);

    d_prolongation_refine_schedules.resize(d_finest_ln+1);
    d_restriction_coarsen_schedules.resize(d_finest_ln+1);
    d_ghostfill_nocoarse_refine_schedules.resize(d_finest_ln+1);
    d_synch_refine_schedules.resize(d_finest_ln+1);

    d_prolongation_refine_algorithm = new RefineAlgorithm<NDIM>();
    d_restriction_coarsen_algorithm = new CoarsenAlgorithm<NDIM>();
    d_ghostfill_nocoarse_refine_algorithm = new RefineAlgorithm<NDIM>();
    d_synch_refine_algorithm = new RefineAlgorithm<NDIM>();

    d_prolongation_refine_algorithm->registerRefine(
        d_side_scratch_idx,
        solution.getComponentDescriptorIndex(0),
        d_side_scratch_idx,
        d_U_prolongation_refine_operator,
        d_U_op_stencil_fill_pattern);
    d_prolongation_refine_algorithm->registerRefine(
        d_cell_scratch_idx,
        solution.getComponentDescriptorIndex(1),
        d_cell_scratch_idx,
        d_P_prolongation_refine_operator,
        d_P_op_stencil_fill_pattern);

    d_restriction_coarsen_algorithm->registerCoarsen(
        d_side_scratch_idx,
        rhs.getComponentDescriptorIndex(0),
        d_U_restriction_coarsen_operator);
    d_restriction_coarsen_algorithm->registerCoarsen(
        d_cell_scratch_idx,
        rhs.getComponentDescriptorIndex(1),
        d_P_restriction_coarsen_operator);

    d_ghostfill_nocoarse_refine_algorithm->registerRefine(
        solution.getComponentDescriptorIndex(0),
        solution.getComponentDescriptorIndex(0),
        solution.getComponentDescriptorIndex(0),
        Pointer<RefineOperator<NDIM> >(),
        d_U_op_stencil_fill_pattern);
    d_ghostfill_nocoarse_refine_algorithm->registerRefine(
        solution.getComponentDescriptorIndex(1),
        solution.getComponentDescriptorIndex(1),
        solution.getComponentDescriptorIndex(1),
        Pointer<RefineOperator<NDIM> >(),
        d_P_op_stencil_fill_pattern);

    d_synch_refine_algorithm->registerRefine(
        solution.getComponentDescriptorIndex(0),
        solution.getComponentDescriptorIndex(0),
        solution.getComponentDescriptorIndex(0),
        Pointer<RefineOperator<NDIM> >(),
        d_U_synch_fill_pattern);

    std::vector<RefinePatchStrategy<NDIM>*> bc_op_ptrs(2);
    bc_op_ptrs[0] = d_U_bc_op;
    bc_op_ptrs[1] = d_P_bc_op;
    d_U_P_bc_op = new RefinePatchStrategySet(bc_op_ptrs.begin(), bc_op_ptrs.end(), false);

    for (int dst_ln = d_coarsest_ln+1; dst_ln <= d_finest_ln; ++dst_ln)
    {
        d_prolongation_refine_schedules[dst_ln] =
            d_prolongation_refine_algorithm->createSchedule(
                d_hierarchy->getPatchLevel(dst_ln),
                Pointer<PatchLevel<NDIM> >(),
                dst_ln-1, d_hierarchy, d_prolongation_refine_patch_strategy.getPointer());

        d_ghostfill_nocoarse_refine_schedules[dst_ln] =
            d_ghostfill_nocoarse_refine_algorithm->createSchedule(
                d_hierarchy->getPatchLevel(dst_ln), d_U_P_bc_op);

        d_synch_refine_schedules[dst_ln] =
            d_synch_refine_algorithm->createSchedule(
                d_hierarchy->getPatchLevel(dst_ln));
    }

    d_ghostfill_nocoarse_refine_schedules[d_coarsest_ln] =
        d_ghostfill_nocoarse_refine_algorithm->createSchedule(
            d_hierarchy->getPatchLevel(d_coarsest_ln), d_U_P_bc_op);

    d_synch_refine_schedules[d_coarsest_ln] =
        d_synch_refine_algorithm->createSchedule(
            d_hierarchy->getPatchLevel(d_coarsest_ln));

    for (int dst_ln = d_coarsest_ln; dst_ln < d_finest_ln; ++dst_ln)
    {
        d_restriction_coarsen_schedules[dst_ln] =
            d_restriction_coarsen_algorithm->createSchedule(
                d_hierarchy->getPatchLevel(dst_ln  ),
                d_hierarchy->getPatchLevel(dst_ln+1));
    }

    // Indicate that the operator is initialized.
    d_is_initialized = true;
    d_in_initialize_operator_state = false;

    IBAMR_TIMER_STOP(t_initialize_operator_state);
    return;
}// initializeOperatorState

void
INSStaggeredFACPreconditionerStrategy::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_operator_state);

    const int coarsest_reset_ln =
        (d_in_initialize_operator_state &&
         (d_coarsest_reset_ln != -1) && (d_finest_reset_ln != -1))
        ? d_coarsest_reset_ln : d_coarsest_ln;
    const int finest_reset_ln =
        (d_in_initialize_operator_state &&
         (d_coarsest_reset_ln != -1) && (d_finest_reset_ln != -1))
        ? d_finest_reset_ln : d_finest_ln;
    deallocateOperatorStateSpecialized(coarsest_reset_ln, finest_reset_ln);

    // Deallocate scratch data.
    for (int ln = coarsest_reset_ln; ln <= std::min(d_finest_ln,finest_reset_ln); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_side_scratch_idx)) level->deallocatePatchData(d_side_scratch_idx);
    }

    // Delete the solution and rhs vectors.
    d_solution->resetLevels(d_solution->getCoarsestLevelNumber(), std::min(d_solution->getFinestLevelNumber(),d_hierarchy->getFinestLevelNumber()));
    d_solution->freeVectorComponents();
    d_solution.setNull();

    d_rhs->resetLevels(d_rhs->getCoarsestLevelNumber(), std::min(d_rhs->getFinestLevelNumber(),d_hierarchy->getFinestLevelNumber()));
    d_rhs->freeVectorComponents();
    d_rhs.setNull();

    // Only fully deallocate operator data when we are not reinitializing the
    // operator.
    if (!d_in_initialize_operator_state)
    {
        d_hierarchy.setNull();
        d_coarsest_ln = -1;
        d_finest_ln   = -1;

        d_hier_bdry_fill_ops.clear();
        d_hier_math_ops.clear();

        d_U_prolongation_refine_operator    .setNull();
        d_P_prolongation_refine_operator    .setNull();
        d_prolongation_refine_patch_strategy.setNull();
        d_prolongation_refine_algorithm     .setNull();
        d_prolongation_refine_schedules     .resize(0);

        d_U_restriction_coarsen_operator.setNull();
        d_P_restriction_coarsen_operator.setNull();
        d_restriction_coarsen_algorithm .setNull();
        d_restriction_coarsen_schedules .resize(0);

        d_ghostfill_nocoarse_refine_algorithm.setNull();
        d_ghostfill_nocoarse_refine_schedules.resize(0);

        d_synch_refine_algorithm.setNull();
        d_synch_refine_schedules.resize(0);
    }

    delete d_U_P_bc_op;

    // Clear the "reset level" range.
    d_coarsest_reset_ln = -1;
    d_finest_reset_ln   = -1;

    // Indicate that the operator is not initialized.
    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_operator_state);
    return;
}// deallocateOperatorState

/////////////////////////////// PROTECTED ////////////////////////////////////

void
INSStaggeredFACPreconditionerStrategy::xeqScheduleProlongation(
    const std::pair<int,int>& dst_idxs,
    const std::pair<int,int>& src_idxs,
    const int dst_ln)
{
    const int U_dst_idx = dst_idxs.first;
    const int U_src_idx = src_idxs.first;
    d_U_bc_op->setPatchDataIndex(U_dst_idx);
    d_U_bc_op->setHomogeneousBc(true);
    d_U_cf_bdry_op->setPatchDataIndex(U_dst_idx);

    const int P_dst_idx = dst_idxs.second;
    const int P_src_idx = src_idxs.second;
    d_P_bc_op->setPatchDataIndex(P_dst_idx);
    d_P_bc_op->setHomogeneousBc(true);
    d_P_cf_bdry_op->setPatchDataIndex(P_dst_idx);

    RefineAlgorithm<NDIM> refine_alg;
    refine_alg.registerRefine(U_dst_idx, U_src_idx, U_dst_idx, d_U_prolongation_refine_operator, d_U_op_stencil_fill_pattern);
    refine_alg.registerRefine(P_dst_idx, P_src_idx, P_dst_idx, d_P_prolongation_refine_operator, d_P_op_stencil_fill_pattern);
    refine_alg.resetSchedule(d_prolongation_refine_schedules[dst_ln]);
    d_prolongation_refine_schedules[dst_ln]->fillData(d_new_time);
    d_prolongation_refine_algorithm->resetSchedule(d_prolongation_refine_schedules[dst_ln]);
    return;
}// xeqScheduleProlongation

void
INSStaggeredFACPreconditionerStrategy::xeqScheduleRestriction(
    const std::pair<int,int>& dst_idxs,
    const std::pair<int,int>& src_idxs,
    const int dst_ln)
{
    const int U_dst_idx = dst_idxs.first;
    const int U_src_idx = src_idxs.first;

    const int P_dst_idx = dst_idxs.second;
    const int P_src_idx = src_idxs.second;

    CoarsenAlgorithm<NDIM> coarsen_alg;
    coarsen_alg.registerCoarsen(U_dst_idx, U_src_idx, d_U_restriction_coarsen_operator);
    coarsen_alg.registerCoarsen(P_dst_idx, P_src_idx, d_P_restriction_coarsen_operator);
    coarsen_alg.resetSchedule(d_restriction_coarsen_schedules[dst_ln]);
    d_restriction_coarsen_schedules[dst_ln]->coarsenData();
    d_restriction_coarsen_algorithm->resetSchedule(d_restriction_coarsen_schedules[dst_ln]);
    return;
}// xeqScheduleRestriction

void
INSStaggeredFACPreconditionerStrategy::xeqScheduleGhostFillNoCoarse(
    const std::pair<int,int>& dst_idxs,
    const int dst_ln)
{
    const int U_dst_idx = dst_idxs.first;
    d_U_bc_op->setPatchDataIndex(U_dst_idx);
    d_U_bc_op->setHomogeneousBc(true);

    const int P_dst_idx = dst_idxs.second;
    d_P_bc_op->setPatchDataIndex(P_dst_idx);
    d_P_bc_op->setHomogeneousBc(true);

    RefineAlgorithm<NDIM> refine_alg;
    refine_alg.registerRefine(U_dst_idx, U_dst_idx, U_dst_idx, Pointer<RefineOperator<NDIM> >(), d_U_op_stencil_fill_pattern);
    refine_alg.registerRefine(P_dst_idx, P_dst_idx, P_dst_idx, Pointer<RefineOperator<NDIM> >(), d_P_op_stencil_fill_pattern);
    refine_alg.resetSchedule(d_ghostfill_nocoarse_refine_schedules[dst_ln]);
    d_ghostfill_nocoarse_refine_schedules[dst_ln]->fillData(d_new_time);
    d_ghostfill_nocoarse_refine_algorithm->resetSchedule(d_ghostfill_nocoarse_refine_schedules[dst_ln]);
    return;
}// xeqScheduleGhostFillNoCoarse

void
INSStaggeredFACPreconditionerStrategy::xeqScheduleDataSynch(
    const int U_dst_idx,
    const int dst_ln)
{
    RefineAlgorithm<NDIM> refine_alg;
    refine_alg.registerRefine(U_dst_idx, U_dst_idx, U_dst_idx, Pointer<RefineOperator<NDIM> >(), d_U_synch_fill_pattern);
    refine_alg.resetSchedule(d_synch_refine_schedules[dst_ln]);
    d_synch_refine_schedules[dst_ln]->fillData(d_new_time);
    d_synch_refine_algorithm->resetSchedule(d_synch_refine_schedules[dst_ln]);
    return;
}// xeqScheduleDataSynch

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
