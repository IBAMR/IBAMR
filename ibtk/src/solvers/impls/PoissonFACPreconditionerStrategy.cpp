// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2023 by the IBAMR developers
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

#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"
#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/FACPreconditionerStrategy.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/PoissonFACPreconditionerStrategy.h"
#include "ibtk/RefinePatchStrategySet.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/samrai_compatibility_names.h"

#include "MultiblockDataTranslator.h"
#include "SAMRAIBox.h"
#include "SAMRAICartesianGridGeometry.h"
#include "SAMRAICoarsenAlgorithm.h"
#include "SAMRAICoarsenOperator.h"
#include "SAMRAICoarsenSchedule.h"
#include "SAMRAIDatabase.h"
#include "SAMRAIHierarchyDataOpsManager.h"
#include "SAMRAIHierarchyDataOpsReal.h"
#include "SAMRAIIntVector.h"
#include "SAMRAILocationIndexRobinBcCoefs.h"
#include "SAMRAIPatchDescriptor.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAIPoissonSpecifications.h"
#include "SAMRAIRefineAlgorithm.h"
#include "SAMRAIRefineOperator.h"
#include "SAMRAIRefinePatchStrategy.h"
#include "SAMRAIRefineSchedule.h"
#include "SAMRAIRobinBcCoefStrategy.h"
#include "SAMRAISAMRAIVectorReal.h"
#include "SAMRAITimer.h"
#include "SAMRAITimerManager.h"
#include "SAMRAIUtilities.h"
#include "SAMRAIVariable.h"
#include "SAMRAIVariableContext.h"
#include "SAMRAIVariableDatabase.h"
#include "SAMRAIVariableFillPattern.h"

#include <algorithm>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static SAMRAITimer* t_restrict_residual;
static SAMRAITimer* t_prolong_error;
static SAMRAITimer* t_prolong_error_and_correct;
static SAMRAITimer* t_initialize_operator_state;
static SAMRAITimer* t_deallocate_operator_state;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

PoissonFACPreconditionerStrategy::PoissonFACPreconditionerStrategy(std::string object_name,
                                                                   SAMRAIPointer<SAMRAIVariable> scratch_var,
                                                                   const int ghost_cell_width,
                                                                   const SAMRAIPointer<SAMRAIDatabase> input_db,
                                                                   const std::string& default_options_prefix)
    : FACPreconditionerStrategy(std::move(object_name)),
      d_poisson_spec(d_object_name + "::poisson_spec"),
      d_default_bc_coef(new SAMRAILocationIndexRobinBcCoefs(d_object_name + "::default_bc_coef",
                                                            SAMRAIPointer<SAMRAIDatabase>(nullptr))),
      d_bc_coefs(1, d_default_bc_coef.get()),
      d_gcw(ghost_cell_width),
      d_coarse_solver_default_options_prefix(default_options_prefix + "_coarse")
{
    // Initialize the Poisson specifications.
    d_poisson_spec.setCZero();
    d_poisson_spec.setDConstant(-1.0);

    // Setup a default boundary condition object that specifies homogeneous
    // Dirichlet boundary conditions.
    auto p_default_bc_coef = dynamic_cast<SAMRAILocationIndexRobinBcCoefs*>(d_default_bc_coef.get());
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        p_default_bc_coef->setBoundaryValue(2 * d, 0.0);
        p_default_bc_coef->setBoundaryValue(2 * d + 1, 0.0);
    }

    // Get values from the input database.
    if (input_db)
    {
        if (input_db->keyExists("ghost_cell_width"))
        {
            int gcw = input_db->getInteger("ghost_cell_width");
            for (int d = 0; d < NDIM; ++d) d_gcw(d) = gcw;
        }
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
    SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
    d_context = var_db->getContext(d_object_name + "::CONTEXT");
    const SAMRAIIntVector ghosts = d_gcw;
    if (var_db->checkVariableExists(scratch_var->getName()))
    {
        scratch_var = var_db->getVariable(scratch_var->getName());
        d_scratch_idx = var_db->mapVariableAndContextToIndex(scratch_var, d_context);
        var_db->removePatchDataIndex(d_scratch_idx);
    }
    d_scratch_idx = var_db->registerVariableAndContext(scratch_var, d_context, ghosts);

    // Setup Timers.
    IBTK_DO_ONCE(t_restrict_residual = SAMRAITimerManager::getManager()->getTimer(
                     "IBTK::PoissonFACPreconditionerStrategy::restrictResidual()");
                 t_prolong_error = SAMRAITimerManager::getManager()->getTimer(
                     "IBTK::PoissonFACPreconditionerStrategy::prolongError()");
                 t_prolong_error_and_correct = SAMRAITimerManager::getManager()->getTimer(
                     "IBTK::PoissonFACPreconditionerStrategy::prolongErrorAndCorrect()");
                 t_initialize_operator_state = SAMRAITimerManager::getManager()->getTimer(
                     "IBTK::PoissonFACPreconditionerStrategy::initializeOperatorState()");
                 t_deallocate_operator_state = SAMRAITimerManager::getManager()->getTimer(
                     "IBTK::PoissonFACPreconditionerStrategy::deallocateOperatorState()"););
    return;
} // PoissonFACPreconditionerStrategy

PoissonFACPreconditionerStrategy::~PoissonFACPreconditionerStrategy()
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::~PoissonFACPreconditionerStrategy()\n"
                                 << "  subclass must call deallocateOperatorState in subclass destructor" << std::endl);
    }
    return;
} // ~PoissonFACPreconditionerStrategy

void
PoissonFACPreconditionerStrategy::setPoissonSpecifications(const SAMRAIPoissonSpecifications& poisson_spec)
{
    d_poisson_spec = poisson_spec;
    return;
} // setPoissonSpecifications

void
PoissonFACPreconditionerStrategy::setPhysicalBcCoef(SAMRAIRobinBcCoefStrategy* const bc_coef)
{
    setPhysicalBcCoefs(std::vector<SAMRAIRobinBcCoefStrategy*>(1, bc_coef));
    return;
} // setPhysicalBcCoef

void
PoissonFACPreconditionerStrategy::setPhysicalBcCoefs(const std::vector<SAMRAIRobinBcCoefStrategy*>& bc_coefs)
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
            d_bc_coefs[l] = d_default_bc_coef.get();
        }
    }
    return;
} // setPhysicalBcCoefs

void
PoissonFACPreconditionerStrategy::setResetLevels(const int coarsest_ln, const int finest_ln)
{
#if !defined(NDEBUG)
    TBOX_ASSERT((coarsest_ln == invalid_level_number && finest_ln == invalid_level_number) ||
                (coarsest_ln >= 0 && finest_ln >= coarsest_ln));
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
                                 << "  cannot be called while operator state is initialized" << std::endl);
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
                                 << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_restriction_method = restriction_method;
    return;
} // setRestrictionMethod

void
PoissonFACPreconditionerStrategy::setToZero(SAMRAISAMRAIVectorReal<double>& vec, int level_num)
{
    const int data_idx = vec.getComponentDescriptorIndex(0);
    d_level_data_ops[level_num]->setToScalar(data_idx, 0.0, /*interior_only*/ false);
    return;
} // setToZero

void
PoissonFACPreconditionerStrategy::restrictResidual(const SAMRAISAMRAIVectorReal<double>& src,
                                                   SAMRAISAMRAIVectorReal<double>& dst,
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
PoissonFACPreconditionerStrategy::prolongError(const SAMRAISAMRAIVectorReal<double>& src,
                                               SAMRAISAMRAIVectorReal<double>& dst,
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
PoissonFACPreconditionerStrategy::prolongErrorAndCorrect(const SAMRAISAMRAIVectorReal<double>& src,
                                                         SAMRAISAMRAIVectorReal<double>& dst,
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
PoissonFACPreconditionerStrategy::initializeOperatorState(const SAMRAISAMRAIVectorReal<double>& solution,
                                                          const SAMRAISAMRAIVectorReal<double>& rhs)
{
    IBTK_TIMER_START(t_initialize_operator_state);

    d_in_initialize_operator_state = true;

    // Cache the level range to be reset.
    //
    // NOTE: We cannot use d_coarsest_reset_ln and d_finest_reset_ln since those
    // values are reset by deallocateOperatorState().
    const int coarsest_reset_ln =
        (d_coarsest_reset_ln != invalid_level_number && d_finest_reset_ln != invalid_level_number ?
             d_coarsest_reset_ln :
             solution.getCoarsestLevelNumber());
    const int finest_reset_ln =
        (d_coarsest_reset_ln != invalid_level_number && d_finest_reset_ln != invalid_level_number ?
             d_finest_reset_ln :
             solution.getFinestLevelNumber());

    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Setup solution and rhs vectors.
    d_solution = solution.cloneVector(solution.getName());
    d_rhs = rhs.cloneVector(rhs.getName());

    SAMRAIPointer<SAMRAIVariable> sol_var = d_solution->getComponentVariable(0);
    const int sol_idx = d_solution->getComponentDescriptorIndex(0);
    const int rhs_idx = d_rhs->getComponentDescriptorIndex(0);

    // Reset the hierarchy configuration.
    d_hierarchy = solution.getPatchHierarchy();
    d_coarsest_ln = solution.getCoarsestLevelNumber();
    d_finest_ln = solution.getFinestLevelNumber();

#if !defined(NDEBUG)
    // To prevent very obtuse errors later on, we check that the rhs index has the correct ghost cell width.
    SAMRAIPointer<SAMRAIPatchDescriptor> pd = SAMRAIVariableDatabase::getDatabase()->getPatchDescriptor();
    const SAMRAIIntVector& gcw = pd->getPatchDataFactory(rhs_idx)->getGhostCellWidth();
    if (gcw != d_gcw)
    {
        TBOX_ERROR(
            d_object_name +
                "::initializeOperatorState(): RHS index does not have the correct ghost width. RHS has ghost width of "
            << gcw << ". Expecting a ghost cell width of " << d_gcw << ".\n");
    }
#endif

    // Perform implementation-specific initialization.
    initializeOperatorStateSpecialized(solution, rhs, coarsest_reset_ln, finest_reset_ln);
#if !defined(NDEBUG)
    TBOX_ASSERT(d_bc_op);
    TBOX_ASSERT(d_cf_bdry_op);
#endif

    // Setup level operators.
    d_level_data_ops.resize(d_finest_ln + 1);
    d_level_bdry_fill_ops.resize(d_finest_ln + 1, nullptr);
    d_level_math_ops.resize(d_finest_ln + 1, nullptr);
    SAMRAIHierarchyDataOpsManager* hier_data_ops_manager = SAMRAIHierarchyDataOpsManager::getManager();
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
        SAMRAIPointer<SAMRAIPatchLevel> level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_scratch_idx)) level->allocatePatchData(d_scratch_idx);
    }

    // Get the transfer operators.
    SAMRAIPointer<SAMRAICartesianGridGeometry> geometry = d_hierarchy->getGridGeometry();
    d_prolongation_refine_operator = geometry->lookupRefineOperator(sol_var, d_prolongation_method);
    d_restriction_coarsen_operator = geometry->lookupCoarsenOperator(sol_var, d_restriction_method);
    d_cf_bdry_op->setConsistentInterpolationScheme(false);
    d_cf_bdry_op->setPatchDataIndex(d_scratch_idx);
    d_cf_bdry_op->setPatchHierarchy(d_hierarchy);

    // Make space for saving communication schedules.  There is no need to
    // delete the old schedules first because we have deallocated the solver
    // state above.
    std::vector<SAMRAIRefinePatchStrategy*> prolongation_refine_patch_strategies;
    prolongation_refine_patch_strategies.push_back(d_cf_bdry_op);
    prolongation_refine_patch_strategies.push_back(d_bc_op);
    d_prolongation_refine_patch_strategy = new RefinePatchStrategySet(
        prolongation_refine_patch_strategies.begin(), prolongation_refine_patch_strategies.end(), false);

    d_prolongation_refine_schedules.resize(d_finest_ln + 1);
    d_restriction_coarsen_schedules.resize(d_finest_ln);
    d_ghostfill_nocoarse_refine_schedules.resize(d_finest_ln + 1);
    d_synch_refine_schedules.resize(d_finest_ln + 1);

    d_prolongation_refine_algorithm = new SAMRAIRefineAlgorithm();
    d_restriction_coarsen_algorithm = new SAMRAICoarsenAlgorithm();
    d_ghostfill_nocoarse_refine_algorithm = new SAMRAIRefineAlgorithm();
    d_synch_refine_algorithm = new SAMRAIRefineAlgorithm();

    d_prolongation_refine_algorithm->registerRefine(
        d_scratch_idx, sol_idx, d_scratch_idx, d_prolongation_refine_operator, d_op_stencil_fill_pattern);
    d_restriction_coarsen_algorithm->registerCoarsen(d_scratch_idx, rhs_idx, d_restriction_coarsen_operator);
    d_ghostfill_nocoarse_refine_algorithm->registerRefine(
        sol_idx, sol_idx, sol_idx, SAMRAIPointer<SAMRAIRefineOperator>(), d_op_stencil_fill_pattern);
    d_synch_refine_algorithm->registerRefine(
        sol_idx, sol_idx, sol_idx, SAMRAIPointer<SAMRAIRefineOperator>(), d_synch_fill_pattern);

    // TODO: Here we take a pessimistic approach and are recreating refine schedule for
    // (coarsest_reset_ln - 1) level as well.
    for (int dst_ln = std::max(d_coarsest_ln + 1, coarsest_reset_ln - 1); dst_ln <= finest_reset_ln; ++dst_ln)
    {
        d_prolongation_refine_schedules[dst_ln] =
            d_prolongation_refine_algorithm->createSchedule(d_hierarchy->getPatchLevel(dst_ln),
                                                            SAMRAIPointer<SAMRAIPatchLevel>(),
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

    const int coarsest_reset_ln = (d_in_initialize_operator_state && (d_coarsest_reset_ln != invalid_level_number) &&
                                   (d_finest_reset_ln != invalid_level_number)) ?
                                      d_coarsest_reset_ln :
                                      d_coarsest_ln;
    const int finest_reset_ln = (d_in_initialize_operator_state && (d_coarsest_reset_ln != invalid_level_number) &&
                                 (d_finest_reset_ln != invalid_level_number)) ?
                                    d_finest_reset_ln :
                                    d_finest_ln;
    deallocateOperatorStateSpecialized(coarsest_reset_ln, finest_reset_ln);

    // Deallocate scratch data.
    for (int ln = coarsest_reset_ln;
         ln <= std::min({ d_finest_ln, finest_reset_ln, d_hierarchy->getFinestLevelNumber() });
         ++ln)
    {
        SAMRAIPointer<SAMRAIPatchLevel> level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_scratch_idx)) level->deallocatePatchData(d_scratch_idx);
    }

    // Delete the solution and rhs vectors.
    free_vector_components(*d_solution);
    d_solution.setNull();

    free_vector_components(*d_rhs);
    d_rhs.setNull();

    // Only fully deallocate operator data when we are not reinitializing the
    // operator.
    if (!d_in_initialize_operator_state)
    {
        d_hierarchy.setNull();
        d_coarsest_ln = invalid_level_number;
        d_finest_ln = invalid_level_number;

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
    d_coarsest_reset_ln = invalid_level_number;
    d_finest_reset_ln = invalid_level_number;

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
    for (const auto& bc_coef : d_bc_coefs)
    {
        auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coef);
        if (extended_bc_coef)
        {
            extended_bc_coef->setTargetPatchDataIndex(dst_idx);
            extended_bc_coef->setHomogeneousBc(true);
        }
    }
    SAMRAIRefineAlgorithm refiner;
    refiner.registerRefine(dst_idx, src_idx, dst_idx, d_prolongation_refine_operator, d_op_stencil_fill_pattern);
    refiner.resetSchedule(d_prolongation_refine_schedules[dst_ln]);
    d_prolongation_refine_schedules[dst_ln]->fillData(d_solution_time);
    d_prolongation_refine_algorithm->resetSchedule(d_prolongation_refine_schedules[dst_ln]);
    for (const auto& bc_coef : d_bc_coefs)
    {
        auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coef);
        if (extended_bc_coef) extended_bc_coef->clearTargetPatchDataIndex();
    }
    return;
} // xeqScheduleProlongation

void
PoissonFACPreconditionerStrategy::xeqScheduleRestriction(const int dst_idx, const int src_idx, const int dst_ln)
{
    SAMRAICoarsenAlgorithm coarsener;
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
    for (const auto& bc_coef : d_bc_coefs)
    {
        auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coef);
        if (extended_bc_coef)
        {
            extended_bc_coef->setTargetPatchDataIndex(dst_idx);
            extended_bc_coef->setHomogeneousBc(true);
        }
    }
    SAMRAIRefineAlgorithm refiner;
    refiner.registerRefine(dst_idx, dst_idx, dst_idx, SAMRAIPointer<SAMRAIRefineOperator>(), d_op_stencil_fill_pattern);
    refiner.resetSchedule(d_ghostfill_nocoarse_refine_schedules[dst_ln]);
    d_ghostfill_nocoarse_refine_schedules[dst_ln]->fillData(d_solution_time);
    d_ghostfill_nocoarse_refine_algorithm->resetSchedule(d_ghostfill_nocoarse_refine_schedules[dst_ln]);
    for (const auto& bc_coef : d_bc_coefs)
    {
        auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coef);
        if (extended_bc_coef) extended_bc_coef->clearTargetPatchDataIndex();
    }
    return;
} // xeqScheduleGhostFillNoCoarse

void
PoissonFACPreconditionerStrategy::xeqScheduleDataSynch(const int dst_idx, const int dst_ln)
{
    SAMRAIRefineAlgorithm refiner;
    refiner.registerRefine(dst_idx, dst_idx, dst_idx, SAMRAIPointer<SAMRAIRefineOperator>(), d_synch_fill_pattern);
    refiner.resetSchedule(d_synch_refine_schedules[dst_ln]);
    d_synch_refine_schedules[dst_ln]->fillData(d_solution_time);
    d_synch_refine_algorithm->resetSchedule(d_synch_refine_schedules[dst_ln]);
    return;
} // xeqScheduleDataSynch

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
