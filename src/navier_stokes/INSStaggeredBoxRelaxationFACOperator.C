// Filename: INSStaggeredBoxRelaxationFACOperator.C
// Last modified: <11.Jun.2010 18:33:16 griffith@boyce-griffiths-mac-pro.local>
// Created on 11 Jun 2010 by Boyce Griffith (griffith@boyce-griffiths-mac-pro.local)

#include "INSStaggeredBoxRelaxationFACOperator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/CartCellDoubleQuadraticCFInterpolation.h>
#include <ibtk/CartSideDoubleCubicCoarsen.h>
#include <ibtk/CartSideDoubleQuadraticCFInterpolation.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/NormOps.h>
#include <ibtk/RefinePatchStrategySet.h>
#include <ibtk/SideNoCornersFillPattern.h>
#include <ibtk/SideSynchCopyFillPattern.h>
#include <ibtk/SideSynchCopyTransactionFactory.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>

// C++ STDLIB INCLUDES
#include <algorithm>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_restrict_solution;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_restrict_residual;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_prolong_error_and_correct;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_smooth_error;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_solve_coarsest_level;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_compute_composite_residual_on_level;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_compute_residual_norm;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_operator_state;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_deallocate_operator_state;

// Number of ghosts cells used for each variable quantity.
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);
static const int SIDEG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredBoxRelaxationFACOperator::INSStaggeredBoxRelaxationFACOperator(
    const std::string& object_name,
    const INSCoefs& problem_coefs,
    const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db)
    : d_object_name(object_name),
      d_is_initialized(false),
      d_solution(NULL),
      d_rhs(NULL),
      d_gcw(SIDEG),
      d_patch_vec_e(),
      d_patch_vec_f(),
      d_patch_mat(),
      d_patch_bc_box_overlap(),
      d_hierarchy(),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_in_initialize_operator_state(false),
      d_coarsest_reset_ln(-1),
      d_finest_reset_ln(-1),
      d_problem_coefs(problem_coefs),
      d_smoother_choice("additive"),
      d_U_prolongation_method("CONSTANT_REFINE"),
      d_P_prolongation_method("LINEAR_REFINE"),
      d_U_restriction_method("CUBIC_COARSEN"),
      d_P_restriction_method("CONSERVATIVE_COARSEN"),
      d_preconditioner(NULL),
      d_fac_max_cycles(1),
      d_fac_uses_presmoothing(false),
      d_fac_initial_guess_nonzero(false),
      d_skip_restrict_sol(true),
      d_skip_restrict_residual(false),
      d_coarse_solver_choice("block_jacobi"),
      d_coarse_solver_tol(1.0e-6),
      d_coarse_solver_max_its(10),
      d_context(NULL),
      d_side_scratch_idx(-1),
      d_cell_scratch_idx(-1),
      d_U_bc_op(NULL),
      d_default_U_bc_coef(new SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>(d_object_name+"::default_U_bc_coef", SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>(NULL))),
      d_U_bc_coefs(std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(NDIM,d_default_U_bc_coef)),
      d_P_bc_op(NULL),
      d_default_P_bc_coef(new SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>(d_object_name+"::default_P_bc_coef", SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>(NULL))),
      d_P_bc_coef(d_default_P_bc_coef),
      d_homogeneous_bc(false),
      d_current_time(0.0),
      d_new_time(0.0),
      d_dt(0.0),
      d_U_op_stencil_fill_pattern(),
      d_P_op_stencil_fill_pattern(),
      d_U_synch_fill_pattern(),
      d_U_prolongation_refine_operator(),
      d_P_prolongation_refine_operator(),
      d_U_cf_bdry_op(),
      d_P_cf_bdry_op(),
      d_prolongation_refine_patch_strategy(),
      d_prolongation_refine_algorithm(),
      d_prolongation_refine_schedules(),
      d_U_urestriction_coarsen_operator(),
      d_P_urestriction_coarsen_operator(),
      d_urestriction_coarsen_algorithm(),
      d_urestriction_coarsen_schedules(),
      d_U_rrestriction_coarsen_operator(),
      d_P_rrestriction_coarsen_operator(),
      d_rrestriction_coarsen_algorithm(),
      d_rrestriction_coarsen_schedules(),
      d_U_ghostfill_refine_operator(),
      d_P_ghostfill_refine_operator(),
      d_ghostfill_refine_algorithm(),
      d_ghostfill_refine_schedules(),
      d_U_ghostfill_nocoarse_refine_operator(),
      d_P_ghostfill_nocoarse_refine_operator(),
      d_ghostfill_nocoarse_refine_algorithm(),
      d_ghostfill_nocoarse_refine_schedules(),
      d_U_synch_refine_operator(),
      d_U_synch_refine_algorithm(),
      d_U_synch_refine_schedules()
{
    // Get values from the input database.
    if (!input_db.isNull())
    {
        d_smoother_choice = input_db->getStringWithDefault("smoother_choice", d_smoother_choice);
        d_U_prolongation_method = input_db->getStringWithDefault("U_prolongation_method", d_U_prolongation_method);
        d_P_prolongation_method = input_db->getStringWithDefault("P_prolongation_method", d_P_prolongation_method);
        d_U_restriction_method = input_db->getStringWithDefault("U_restriction_method", d_U_restriction_method);
        d_P_restriction_method = input_db->getStringWithDefault("P_restriction_method", d_P_restriction_method);
        d_fac_max_cycles = input_db->getIntegerWithDefault("fac_max_cycles", d_fac_max_cycles);
        d_fac_uses_presmoothing = input_db->getBoolWithDefault("fac_uses_presmoothing", d_fac_uses_presmoothing);
        d_fac_initial_guess_nonzero = input_db->getBoolWithDefault("fac_initial_guess_nonzero", d_fac_initial_guess_nonzero);
        d_skip_restrict_sol = input_db->getBoolWithDefault("skip_restrict_sol", d_skip_restrict_sol);
        d_skip_restrict_residual = input_db->getBoolWithDefault("skip_restrict_residual", d_skip_restrict_residual);
        d_coarse_solver_choice = input_db->getStringWithDefault("coarse_solver_choice", d_coarse_solver_choice);
        d_coarse_solver_tol = input_db->getDoubleWithDefault("coarse_solver_tolerance", d_coarse_solver_tol);
        d_coarse_solver_max_its = input_db->getIntegerWithDefault("coarse_solver_max_iterations", d_coarse_solver_max_its);
    }
    sanityCheck();

    // Create the hypre solver, if needed.
    setCoarsestLevelSolverChoice(d_coarse_solver_choice);

    // Setup a default boundary condition object that specifies homogeneous
    // Dirichlet boundary conditions for the velocity and homogeneous Neumann
    // boundary conditions for the pressure.
    for (int d = 0; d < NDIM; ++d)
    {
        d_default_U_bc_coef->setBoundaryValue(2*d  ,0.0);
        d_default_U_bc_coef->setBoundaryValue(2*d+1,0.0);
        d_default_P_bc_coef->setBoundarySlope(2*d  ,0.0);
        d_default_P_bc_coef->setBoundarySlope(2*d+1,0.0);
    }

    // Initialize the boundary conditions objects.
    setHomogeneousBc(d_homogeneous_bc);
    setPhysicalBcCoefs(std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(NDIM,d_default_U_bc_coef),d_default_P_bc_coef);

    // Setup scratch variables.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(d_object_name+"::CONTEXT");

    const SAMRAI::hier::IntVector<NDIM> side_ghosts = d_gcw;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > side_scratch_var =
        new SAMRAI::pdat::SideVariable<NDIM,double>(d_object_name+"::side_scratch");
    if (var_db->checkVariableExists(side_scratch_var->getName()))
    {
        side_scratch_var = var_db->getVariable(side_scratch_var->getName());
        d_side_scratch_idx = var_db->mapVariableAndContextToIndex(side_scratch_var, d_context);
        var_db->removePatchDataIndex(d_side_scratch_idx);
    }
    d_side_scratch_idx = var_db->registerVariableAndContext(side_scratch_var, d_context, side_ghosts);

    const SAMRAI::hier::IntVector<NDIM> cell_ghosts = d_gcw;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > cell_scratch_var =
        new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::cell_scratch");
    if (var_db->checkVariableExists(cell_scratch_var->getName()))
    {
        cell_scratch_var = var_db->getVariable(cell_scratch_var->getName());
        d_cell_scratch_idx = var_db->mapVariableAndContextToIndex(cell_scratch_var, d_context);
        var_db->removePatchDataIndex(d_cell_scratch_idx);
    }
    d_cell_scratch_idx = var_db->registerVariableAndContext(cell_scratch_var, d_context, cell_ghosts);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_restrict_solution                   = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBTK::INSStaggeredBoxRelaxationFACOperator::restrictSolution()");
        t_restrict_residual                   = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBTK::INSStaggeredBoxRelaxationFACOperator::restrictResidual()");
        t_prolong_error_and_correct           = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBTK::INSStaggeredBoxRelaxationFACOperator::prolongErrorAndCorrect()");
        t_smooth_error                        = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBTK::INSStaggeredBoxRelaxationFACOperator::smoothError()");
        t_solve_coarsest_level                = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBTK::INSStaggeredBoxRelaxationFACOperator::solveCoarsestLevel()");
        t_compute_composite_residual_on_level = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBTK::INSStaggeredBoxRelaxationFACOperator::computeCompositeResidualOnLevel()");
        t_compute_residual_norm               = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBTK::INSStaggeredBoxRelaxationFACOperator::computeResidualNorm()");
        t_initialize_operator_state           = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBTK::INSStaggeredBoxRelaxationFACOperator::initializeOperatorState()");
        t_deallocate_operator_state           = SAMRAI::tbox::TimerManager::getManager()->getTimer("IBTK::INSStaggeredBoxRelaxationFACOperator::deallocateOperatorState()");
        timers_need_init = false;
    }
    return;
}// INSStaggeredBoxRelaxationFACOperator

INSStaggeredBoxRelaxationFACOperator::~INSStaggeredBoxRelaxationFACOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    delete d_default_U_bc_coef;
    delete d_default_P_bc_coef;
    return;
}// ~INSStaggeredBoxRelaxationFACOperator

void
INSStaggeredBoxRelaxationFACOperator::setPhysicalBcCoefs(
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(U_bc_coefs.size() == NDIM);
#endif
    d_U_bc_coefs.resize(U_bc_coefs.size());
    for (unsigned l = 0; l < U_bc_coefs.size(); ++l)
    {
        if (U_bc_coefs[l] != NULL)
        {
            d_U_bc_coefs[l] = U_bc_coefs[l];
        }
        else
        {
            d_U_bc_coefs[l] = d_default_U_bc_coef;
        }
    }

    if (P_bc_coef != NULL)
    {
        d_P_bc_coef = P_bc_coef;
    }
    else
    {
        d_P_bc_coef = d_default_P_bc_coef;
    }
    return;
}// setPhysicalBcCoefs

void
INSStaggeredBoxRelaxationFACOperator::setHomogeneousBc(
    const bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

void
INSStaggeredBoxRelaxationFACOperator::setTimeInterval(
    const double current_time,
    const double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    d_dt = d_new_time-d_current_time;
    return;
}// setTimeInterval

void
INSStaggeredBoxRelaxationFACOperator::setPreconditioner(
    const SAMRAI::solv::FACPreconditioner<NDIM>* preconditioner)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setPreconditioner()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_preconditioner = preconditioner;
    sanityCheck();
    return;
}// setPreconditioner

void
INSStaggeredBoxRelaxationFACOperator::setResetLevels(
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
INSStaggeredBoxRelaxationFACOperator::setGhostCellWidth(
    const SAMRAI::hier::IntVector<NDIM>& ghost_cell_width)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setGhostCellWidth()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    if (ghost_cell_width.min() == 0)
    {
        TBOX_ERROR(d_object_name << "::setGhostCellWidth()\n"
                   << "  ghost_cell_width.min() must be greater than zero" << std::endl);
    }
    d_gcw = ghost_cell_width;
    sanityCheck();

    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var;

    var_db->mapIndexToVariable(d_side_scratch_idx, var);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > side_scratch_var = var;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!side_scratch_var.isNull());
#endif
    var_db->removePatchDataIndex(d_side_scratch_idx);
    const SAMRAI::hier::IntVector<NDIM> side_ghosts = d_gcw;
    d_side_scratch_idx = var_db->registerVariableAndContext(side_scratch_var, d_context, side_ghosts);
#ifdef DEBUG_CHECK_ASSERTIONS
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideDataFactory<NDIM,double> > side_scratch_pdat_fac = var_db->getPatchDescriptor()->getPatchDataFactory(d_side_scratch_idx);
    TBOX_ASSERT(!side_scratch_pdat_fac.isNull());
    TBOX_ASSERT(side_scratch_pdat_fac->getGhostCellWidth() == d_gcw);
#endif

    var_db->mapIndexToVariable(d_cell_scratch_idx, var);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > cell_scratch_var = var;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!cell_scratch_var.isNull());
#endif
    var_db->removePatchDataIndex(d_cell_scratch_idx);
    const SAMRAI::hier::IntVector<NDIM> cell_ghosts = d_gcw;
    d_cell_scratch_idx = var_db->registerVariableAndContext(cell_scratch_var, d_context, cell_ghosts);
#ifdef DEBUG_CHECK_ASSERTIONS
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataFactory<NDIM,double> > cell_scratch_pdat_fac = var_db->getPatchDescriptor()->getPatchDataFactory(d_cell_scratch_idx);
    TBOX_ASSERT(!cell_scratch_pdat_fac.isNull());
    TBOX_ASSERT(cell_scratch_pdat_fac->getGhostCellWidth() == d_gcw);
#endif
    return;
}// setGhostCellWidth

void
INSStaggeredBoxRelaxationFACOperator::setSmootherChoice(
    const std::string& smoother_choice)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setSmootherChoice()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_smoother_choice = smoother_choice;
    sanityCheck();
    return;
}// setSmootherChoice

void
INSStaggeredBoxRelaxationFACOperator::setCoarsestLevelSolverChoice(
    const std::string& coarse_solver_choice)
{
    d_coarse_solver_choice = coarse_solver_choice;
    sanityCheck();
    return;
}// setCoarsestLevelSolverChoice

void
INSStaggeredBoxRelaxationFACOperator::setCoarsestLevelSolverTolerance(
    double coarse_solver_tol)
{
    d_coarse_solver_tol = coarse_solver_tol;
    sanityCheck();
    return;
}// setCoarsestLevelSolverTolerance

void
INSStaggeredBoxRelaxationFACOperator::setCoarsestLevelSolverMaxIterations(
    int coarse_solver_max_its)
{
    d_coarse_solver_max_its = coarse_solver_max_its;
    sanityCheck();
    return;
}// setCoarsestLevelSolverMaxIterations

void
INSStaggeredBoxRelaxationFACOperator::setProlongationMethods(
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
    sanityCheck();
    return;
}// setProlongationMethods

void
INSStaggeredBoxRelaxationFACOperator::setRestrictionMethods(
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
    sanityCheck();
    return;
}// setRestrictionMethods

void
INSStaggeredBoxRelaxationFACOperator::setFACPreconditionerMaxCycles(
    int fac_max_cycles)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setFACPreconditionerMaxCycles()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_fac_max_cycles = fac_max_cycles;
    sanityCheck();
    return;
}// setFACPreconditionerMaxCycles

void
INSStaggeredBoxRelaxationFACOperator::setFACPreconditionerUsesPresmoothing(
    bool fac_uses_presmoothing)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setFACPreconditionerUsesPresmoothing()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_fac_uses_presmoothing = fac_uses_presmoothing;
    sanityCheck();
    return;
}// setFACPreconditionerUsesPresmoothing

void
INSStaggeredBoxRelaxationFACOperator::setFACPreconditionerInitialGuessNonzero(
    bool fac_initial_guess_nonzero)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setFACPreconditionerInitialGuessNonzero()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_fac_initial_guess_nonzero = fac_initial_guess_nonzero;
    sanityCheck();
    return;
}// setFACPreconditionerInitialGuessNonzero

void
INSStaggeredBoxRelaxationFACOperator::setSkipRestrictSolution(
    bool skip_restrict_sol)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setSkipRestrictSolution()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_skip_restrict_sol = skip_restrict_sol;
    sanityCheck();
    return;
}// setSkipRestrictSolution

void
INSStaggeredBoxRelaxationFACOperator::setSkipRestrictResidual(
    bool skip_restrict_residual)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setSkipRestrictResidual()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_skip_restrict_residual = skip_restrict_residual;
    sanityCheck();
    return;
}// setSkipRestrictResidual

///
///  The following routines:
///
///      restrictSolution(),
///      restrictResidual(),
///      prolongErrorAndCorrect(),
///      smoothError(),
///      solveCoarsestLevel(),
///      computeCompositeResidualOnLevel(),
///      computeResidualNorm(),
///      initializeOperatorState(),
///      deallocateOperatorState()
///
///  are concrete implementations of functions declared in the
///  SAMRAI::solv::FACOperatorStrategy abstract base class.
///

void
INSStaggeredBoxRelaxationFACOperator::restrictSolution(
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& src,
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& dst,
    int dst_ln)
{
    if (d_skip_restrict_sol) return;

    t_restrict_solution->start();

    const int src_idx = src.getComponentDescriptorIndex(0);
    const int dst_idx = dst.getComponentDescriptorIndex(0);

    xeqScheduleURestriction(dst_idx, src_idx, dst_ln);

    static const bool homogeneous_bc = true;
    if (dst_ln == d_coarsest_ln)
    {
        xeqScheduleGhostFillNoCoarse(dst_idx, dst_ln, homogeneous_bc);
    }
    else
    {
        xeqScheduleGhostFill(dst_idx, dst_ln, homogeneous_bc);
    }

    t_restrict_solution->stop();
    return;
}// restrictSolution

void
INSStaggeredBoxRelaxationFACOperator::restrictResidual(
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& src,
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& dst,
    int dst_ln)
{
    if (d_skip_restrict_residual) return;

    t_restrict_residual->start();

    const int src_idx = src.getComponentDescriptorIndex(0);
    const int dst_idx = dst.getComponentDescriptorIndex(0);

    xeqScheduleRRestriction(dst_idx, src_idx, dst_ln);

    t_restrict_residual->stop();
    return;
}// restrictResidual

void
INSStaggeredBoxRelaxationFACOperator::prolongErrorAndCorrect(
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& src,
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& dst,
    int dst_ln)
{
    t_prolong_error_and_correct->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    const bool dst_data_is_prolonged_src_data = !d_fac_uses_presmoothing && !d_fac_initial_guess_nonzero && d_fac_max_cycles == 1;
    TBOX_ASSERT(dst_data_is_prolonged_src_data);
#endif

    // Refine the correction from the coarse level src data directly into the
    // fine level error.
    const int dst_idx = dst.getComponentDescriptorIndex(0);
    const int src_idx = src.getComponentDescriptorIndex(0);
    static const bool homogeneous_bc = true;
    xeqScheduleProlongation(dst_idx, src_idx, dst_ln, homogeneous_bc);

    t_prolong_error_and_correct->stop();
    return;
}// prolongErrorAndCorrect

void
INSStaggeredBoxRelaxationFACOperator::smoothError(
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& error,
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& residual,
    int level_num,
    int num_sweeps)
{
    if (num_sweeps == 0) return;

    t_smooth_error->start();

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    const int error_idx = error.getComponentDescriptorIndex(0);
    const int scratch_idx = d_side_scratch_idx;

    // Cache coarse-fine interface ghost cell values in the "scratch" data.
    if (level_num > d_coarsest_ln)
    {
        SAMRAI::math::HierarchySideDataOpsReal<NDIM,double> hierarchy_sc_data_ops(d_hierarchy, level_num, level_num);
        hierarchy_sc_data_ops.copyData(scratch_idx, error_idx, false);
    }

    // Smooth the error by the specified number of sweeps.
    for (int isweep = 0; isweep < num_sweeps; ++isweep)
    {
        // Re-fill ghost cell data as needed.
        if (level_num > d_coarsest_ln)
        {
            if (isweep > 0)
            {
                // Copy the coarse-fine interface ghost cell values which are
                // cached in the scratch data into the error data.
                for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                    const int patch_num = patch->getPatchNumber();
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> >   error_data = error.getComponentPatchData(0, *patch);
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > scratch_data = patch->getPatchData(scratch_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
                    const SAMRAI::hier::Box<NDIM>& ghost_box = error_data->getGhostBox();
                    TBOX_ASSERT(ghost_box == scratch_data->getGhostBox());
                    TBOX_ASSERT(  error_data->getGhostCellWidth() == d_gcw);
                    TBOX_ASSERT(scratch_data->getGhostCellWidth() == d_gcw);
#endif
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        error_data->getArrayData(axis).copy(
                            scratch_data->getArrayData(axis),
                            d_patch_bc_box_overlap[level_num][patch_num][axis],
                            SAMRAI::hier::IntVector<NDIM>(0));
                    }
                }

                // Fill the non-coarse-fine interface ghost cell values.
                static const bool homogeneous_bc = true;
                xeqScheduleGhostFillNoCoarse(error_idx, level_num, homogeneous_bc);
            }

            // Complete the coarse-fine interface interpolation by computing the
            // normal extension.
            d_U_cf_bdry_op->setPatchDataIndex(error_idx);
            const SAMRAI::hier::IntVector<NDIM>& ratio = level->getRatioToCoarserLevel();
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill = patch->getPatchData(error_idx)->getGhostCellWidth();
                d_U_cf_bdry_op->computeNormalExtension(*patch, ratio, ghost_width_to_fill);
            }
        }
        else if (isweep > 0)
        {
            static const bool homogeneous_bc = true;
            xeqScheduleGhostFillNoCoarse(error_idx, level_num, homogeneous_bc);
        }

        // Smooth the error on the patches.
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const int patch_num = patch->getPatchNumber();
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> >    error_data = error   .getComponentPatchData(0, *patch);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > residual_data = residual.getComponentPatchData(0, *patch);
#ifdef DEBUG_CHECK_ASSERTIONS
            const SAMRAI::hier::Box<NDIM>& ghost_box = error_data->getGhostBox();
            TBOX_ASSERT(ghost_box == residual_data->getGhostBox());
            TBOX_ASSERT(   error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(residual_data->getGhostCellWidth() == d_gcw);
#endif
            for (int axis = 0; axis < NDIM; ++axis)
            {
                // Reset ghost cell values in the copy of the residual data so
                // that patch boundary conditions are properly handled.
                residual_data->getArrayData(axis).copy(
                    error_data->getArrayData(axis),
                    d_patch_bc_box_overlap[level_num][patch_num][axis],
                    SAMRAI::hier::IntVector<NDIM>(0));

                // Setup the PETSc Vec wrappers for the given patch data and axis.
                int ierr;

                Vec& e = d_patch_vec_e[level_num][patch_num][axis];
                Vec& f = d_patch_vec_f[level_num][patch_num][axis];

                ierr = VecPlaceArray(e,    error_data->getPointer(axis,0));  IBTK_CHKERRQ(ierr);
                ierr = VecPlaceArray(f, residual_data->getPointer(axis,0));  IBTK_CHKERRQ(ierr);

                // Smooth the error on the patch using PETSc.  Here, we are
                // approximately solving
                //
                //     Ae = f
                //
                // using an iteration of the form
                //
                //     e <- e + PC(f - Ae) = e + PC(r) = e + x.
                //
                // Presently, we simply employ symmetric Gauss-Seidel as the
                // patch smoother.
                static const double omega = 1.0;
                static const double shift = 0.0;
                static const int its = 1;
                Mat& A = d_patch_mat[level_num][patch_num][axis];
                ierr = MatRelax(A, f, omega, SOR_SYMMETRIC_SWEEP, shift, its, its, e);  IBTK_CHKERRQ(ierr);

                // Reset the PETSc Vec wrappers.
                ierr = VecResetArray(e);  IBTK_CHKERRQ(ierr);
                ierr = VecResetArray(f);  IBTK_CHKERRQ(ierr);
            }
        }
    }

    // Synchronize data along patch boundaries.
    xeqScheduleSideDataSynch(error_idx, level_num);

    t_smooth_error->stop();
    return;
}// smoothError

int
INSStaggeredBoxRelaxationFACOperator::solveCoarsestLevel(
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& error,
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& residual,
    int coarsest_ln)
{
    t_solve_coarsest_level->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(coarsest_ln == d_coarsest_ln);
#endif

    smoothError(error, residual, coarsest_ln, d_coarse_solver_max_its);

    // Re-fill ghost values.
    const int error_idx = error.getComponentDescriptorIndex(0);
    static const bool homogeneous_bc = true;
    xeqScheduleGhostFillNoCoarse(error_idx, coarsest_ln, homogeneous_bc);

    t_solve_coarsest_level->stop();
    static const int ret_val = 1;
    return ret_val;
}// solveCoarsestLevel

void
INSStaggeredBoxRelaxationFACOperator::computeCompositeResidualOnLevel(
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& residual,
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& solution,
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& rhs,
    int level_num,
    bool error_equation_indicator)
{
    t_compute_composite_residual_on_level->start();

    if (!d_fac_uses_presmoothing && !d_fac_initial_guess_nonzero && d_fac_max_cycles == 1)
    {
        // The residual needs to be computed in two different cases:
        //
        // - before each FAC sweep commences, and
        // - after performing any presmoothing sweeps
        //
        // If the FAC preconditioner (a) does not use presmoothing, (b) uses a
        // zero initial guess, and (c) only employs one FAC sweep (only one FAC
        // sweep is needed, for instance, when the preconditioner is being used
        // in conjunction with a Krylov subspace method), then we can simply set
        // the residual equal to the right hand side.
        //
        // IMPORTANT NOTE: We assume that the FAC algorithm being used does not
        // use the composite residual in postsweeps.
        static const bool interior_only = false;
        SAMRAI::math::HierarchySideDataOpsReal<NDIM,double> hierarchy_sc_data_ops(d_hierarchy, level_num, level_num);
        const int dst_idx = residual.getComponentDescriptorIndex(0);
        const int src_idx = rhs.getComponentDescriptorIndex(0);
        if (dst_idx != src_idx)
        {
            hierarchy_sc_data_ops.copyData(dst_idx, src_idx, interior_only);
        }
    }
    else
    {
        TBOX_ERROR("INSStaggeredBoxRelaxationFACOperator::computeResidualOnLevel()\n"
                   << "  current implementation cannot compute residuals,\n"
                   << "  consequently, we require that the FAC algorithm:\n"
                   << "     (1) does not use pre-smoothing,\n"
                   << "     (2) does not use a non-zero initial guess, and\n"
                   << "     (3) only uses one cycle\n"
                   << "  thus the implemented solver is really only suitable for use as a preconditioner." << std::endl);
    }

    t_compute_composite_residual_on_level->stop();
    return;
}// computeCompositeResidualOnLevel

double
INSStaggeredBoxRelaxationFACOperator::computeResidualNorm(
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& residual,
    int fine_ln,
    int coarse_ln)
{
    t_compute_residual_norm->start();

    if (coarse_ln != residual.getCoarsestLevelNumber() ||
        fine_ln   != residual.getFinestLevelNumber())
    {
        TBOX_ERROR("INSStaggeredBoxRelaxationFACOperator::computeResidualNorm()\n"
                   << "  residual can only be computed over the range of levels\n"
                   << "  defined in the SAMRAIVectorReal residual vector" << std::endl);
    }

    t_compute_residual_norm->stop();
    return IBTK::NormOps::L2Norm(&residual);
}// computeResidualNorm

void
INSStaggeredBoxRelaxationFACOperator::initializeOperatorState(
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& solution,
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& rhs)
{
    t_initialize_operator_state->start();
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

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > solution_var = solution.getComponentVariable(0);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >      rhs_var =      rhs.getComponentVariable(0);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideDataFactory<NDIM,double> > solution_pdat_fac = solution_var->getPatchDataFactory();
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideDataFactory<NDIM,double> >      rhs_pdat_fac =      rhs_var->getPatchDataFactory();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!solution_var.isNull());
    TBOX_ASSERT(!rhs_var.isNull());
    TBOX_ASSERT(!solution_pdat_fac.isNull());
    TBOX_ASSERT(!rhs_pdat_fac.isNull());
#endif

    if (solution_pdat_fac->getDefaultDepth() != rhs_pdat_fac->getDefaultDepth())
    {
        TBOX_ERROR("INSStaggeredBoxRelaxationFACOperator::initializeOperatorState()\n"
                   << "  solution and rhs vectors must have the same data depths\n"
                   << "  solution data depth = " << solution_pdat_fac->getDefaultDepth() << "\n"
                   << "  rhs      data depth = " << rhs_pdat_fac     ->getDefaultDepth() << std::endl);
    }

    if (solution_pdat_fac->getDefaultDepth() != 1)
    {
        TBOX_ERROR("INSStaggeredBoxRelaxationFACOperator::initializeOperatorState()\n"
                   << "  velocity data must have depth of 1" << std::endl);
    }

    // Reset the hierarchy configuration.
    d_hierarchy   = solution.getPatchHierarchy();
    d_coarsest_ln = solution.getCoarsestLevelNumber();
    d_finest_ln   = solution.getFinestLevelNumber();

    // Allocate scratch data.
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln); ln <= finest_reset_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_side_scratch_idx)) level->allocatePatchData(d_side_scratch_idx);
    }

    // Get the transfer operators.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > geometry = d_hierarchy->getGridGeometry();
    static bool need_to_add_cubic_coarsen = true;
    if (need_to_add_cubic_coarsen)
    {
        geometry->addSpatialCoarsenOperator(new IBTK::CartSideDoubleCubicCoarsen());
        need_to_add_cubic_coarsen = false;
    }

    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var;

    var_db->mapIndexToVariable(d_side_scratch_idx, var);
    d_U_prolongation_refine_operator = geometry->lookupRefineOperator(var, d_U_prolongation_method);

    d_U_cf_bdry_op = new IBTK::CartSideDoubleQuadraticCFInterpolation();
    d_U_cf_bdry_op->setConsistentInterpolationScheme(false);
    d_U_cf_bdry_op->setPatchDataIndex(d_side_scratch_idx);
    d_U_cf_bdry_op->setPatchHierarchy(d_hierarchy);

    var_db->mapIndexToVariable(d_cell_scratch_idx, var);
    d_P_prolongation_refine_operator = geometry->lookupRefineOperator(var, d_P_prolongation_method);

    d_P_cf_bdry_op = new IBTK::CartCellDoubleQuadraticCFInterpolation();
    d_P_cf_bdry_op->setConsistentInterpolationScheme(false);
    d_P_cf_bdry_op->setPatchDataIndex(d_cell_scratch_idx);
    d_P_cf_bdry_op->setPatchHierarchy(d_hierarchy);

    var_db->mapIndexToVariable(d_side_scratch_idx, var);
    d_U_urestriction_coarsen_operator = geometry->lookupCoarsenOperator(var, d_U_restriction_method);
    d_U_rrestriction_coarsen_operator = d_U_urestriction_coarsen_operator;

    var_db->mapIndexToVariable(d_cell_scratch_idx, var);
    d_P_urestriction_coarsen_operator = geometry->lookupCoarsenOperator(var, d_P_restriction_method);
    d_P_rrestriction_coarsen_operator = d_P_urestriction_coarsen_operator;

    var_db->mapIndexToVariable(d_side_scratch_idx, var);
    d_U_ghostfill_refine_operator = geometry->lookupRefineOperator(var, d_U_prolongation_method);
    d_U_ghostfill_nocoarse_refine_operator = NULL;

    var_db->mapIndexToVariable(d_cell_scratch_idx, var);
    d_P_ghostfill_refine_operator = geometry->lookupRefineOperator(var, d_P_prolongation_method);
    d_P_ghostfill_nocoarse_refine_operator = NULL;

    d_U_synch_refine_operator = NULL;

    // Make space for saving communication schedules.  There is no need to
    // delete the old schedules first because we have deallocated the solver
    // state above.
    d_U_bc_op = new IBTK::CartSideRobinPhysBdryOp(d_side_scratch_idx, d_U_bc_coefs, false);

#if (NDIM == 2)
    d_U_op_stencil_fill_pattern = new IBTK::SideNoCornersFillPattern(SIDEG);
#endif
#if (NDIM == 3)
    d_U_op_stencil_fill_pattern = new IBTK::SideNoCornersFillPattern(SIDEG,false);
#endif
    d_U_synch_fill_pattern = new IBTK::SideSynchCopyFillPattern();

    std::vector<SAMRAI::xfer::RefinePatchStrategy<NDIM>*> prolongation_refine_patch_strategies;
    prolongation_refine_patch_strategies.push_back(d_U_cf_bdry_op);
    prolongation_refine_patch_strategies.push_back(d_U_bc_op);
    d_prolongation_refine_patch_strategy = new IBTK::RefinePatchStrategySet(
        prolongation_refine_patch_strategies.begin(), prolongation_refine_patch_strategies.end(), false);

    d_prolongation_refine_schedules.resize(d_finest_ln+1);
    d_urestriction_coarsen_schedules.resize(d_finest_ln+1);
    d_rrestriction_coarsen_schedules.resize(d_finest_ln+1);
    d_ghostfill_refine_schedules.resize(d_finest_ln+1);
    d_ghostfill_nocoarse_refine_schedules.resize(d_finest_ln+1);
    d_U_synch_refine_schedules.resize(d_finest_ln+1);

    d_prolongation_refine_algorithm = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    d_urestriction_coarsen_algorithm = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
    d_rrestriction_coarsen_algorithm = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
    d_ghostfill_refine_algorithm = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    d_ghostfill_nocoarse_refine_algorithm = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
    d_U_synch_refine_algorithm = new SAMRAI::xfer::RefineAlgorithm<NDIM>();

    d_prolongation_refine_algorithm->registerRefine(
        d_side_scratch_idx,
        solution.getComponentDescriptorIndex(0),
        d_side_scratch_idx,
        d_U_prolongation_refine_operator,
        d_U_op_stencil_fill_pattern);
    d_urestriction_coarsen_algorithm->registerCoarsen(
        d_side_scratch_idx,
        d_side_scratch_idx,
        d_U_urestriction_coarsen_operator);
    d_rrestriction_coarsen_algorithm->registerCoarsen(
        d_side_scratch_idx,
        d_side_scratch_idx,
        d_U_rrestriction_coarsen_operator);
    d_ghostfill_refine_algorithm->registerRefine(
        d_side_scratch_idx,
        d_side_scratch_idx,
        d_side_scratch_idx,
        d_U_ghostfill_refine_operator,
        d_U_op_stencil_fill_pattern);
    d_ghostfill_nocoarse_refine_algorithm->registerRefine(
        d_side_scratch_idx,
        d_side_scratch_idx,
        d_side_scratch_idx,
        d_U_ghostfill_nocoarse_refine_operator,
        d_U_op_stencil_fill_pattern);
    d_U_synch_refine_algorithm->registerRefine(
        solution.getComponentDescriptorIndex(0),
        solution.getComponentDescriptorIndex(0),
        solution.getComponentDescriptorIndex(0),
        d_U_synch_refine_operator,
        d_U_synch_fill_pattern);

    TBOX_ASSERT(false);  // need unified bc op!

    for (int dst_ln = d_coarsest_ln+1; dst_ln <= d_finest_ln; ++dst_ln)
    {
        d_prolongation_refine_schedules[dst_ln] =
            d_prolongation_refine_algorithm->createSchedule(
                d_hierarchy->getPatchLevel(dst_ln),
                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(),
                dst_ln-1, d_hierarchy, d_prolongation_refine_patch_strategy.getPointer());

        d_ghostfill_refine_schedules[dst_ln] =
            d_ghostfill_refine_algorithm->createSchedule(
                d_hierarchy->getPatchLevel(dst_ln),
                dst_ln-1, d_hierarchy, d_U_bc_op.getPointer());

        d_ghostfill_nocoarse_refine_schedules[dst_ln] =
            d_ghostfill_nocoarse_refine_algorithm->createSchedule(
                d_hierarchy->getPatchLevel(dst_ln), d_U_bc_op.getPointer());

        d_U_synch_refine_schedules[dst_ln] =
            d_U_synch_refine_algorithm->createSchedule(
                d_hierarchy->getPatchLevel(dst_ln), NULL, new IBTK::SideSynchCopyTransactionFactory());
    }

    d_ghostfill_nocoarse_refine_schedules[d_coarsest_ln] =
        d_ghostfill_nocoarse_refine_algorithm->createSchedule(
            d_hierarchy->getPatchLevel(d_coarsest_ln), d_U_bc_op.getPointer());

    d_U_synch_refine_schedules[d_coarsest_ln] =
        d_U_synch_refine_algorithm->createSchedule(
            d_hierarchy->getPatchLevel(d_coarsest_ln), NULL, new IBTK::SideSynchCopyTransactionFactory());

    for (int dst_ln = d_coarsest_ln; dst_ln < d_finest_ln; ++dst_ln)
    {
        if (!d_skip_restrict_sol)
        {
            d_urestriction_coarsen_schedules[dst_ln] =
                d_urestriction_coarsen_algorithm->createSchedule(
                    d_hierarchy->getPatchLevel(dst_ln  ),
                    d_hierarchy->getPatchLevel(dst_ln+1));
        }
        if (!d_skip_restrict_residual)
        {
            d_rrestriction_coarsen_schedules[dst_ln] =
                d_rrestriction_coarsen_algorithm->createSchedule(
                    d_hierarchy->getPatchLevel(dst_ln  ),
                    d_hierarchy->getPatchLevel(dst_ln+1));
        }
    }

    // Initialize all cached PETSc data.
    d_patch_vec_e.resize(d_finest_ln+1);
    d_patch_vec_f.resize(d_finest_ln+1);
    d_patch_mat.resize(d_finest_ln+1);
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const int patch_num = patch->getPatchNumber();

            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            const SAMRAI::hier::Box<NDIM>  ghost_box = SAMRAI::hier::Box<NDIM>::grow(patch_box, d_gcw);

            std::vector<Vec>& e = d_patch_vec_e[ln][patch_num];
            std::vector<Vec>& f = d_patch_vec_f[ln][patch_num];
            std::vector<Mat>& A = d_patch_mat  [ln][patch_num];
            e.resize(NDIM);
            f.resize(NDIM);
            A.resize(NDIM);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                int ierr;

                const SAMRAI::hier::Box<NDIM> axis_ghost_box = SAMRAI::pdat::SideGeometry<NDIM>::toSideBox(ghost_box,axis);
                const int size = axis_ghost_box.size();

                ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, size, PETSC_NULL, &e[axis]);  IBTK_CHKERRQ(ierr);
                ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, size, PETSC_NULL, &f[axis]);  IBTK_CHKERRQ(ierr);

                buildPatchOperator(A[axis], d_problem_coefs, patch, axis, d_gcw);
            }
        }
    }

    // Get overlap information for setting patch boundary conditions.
    d_patch_bc_box_overlap.resize(d_finest_ln+1);
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const SAMRAI::hier::BoxArray<NDIM>& box_array = level->getBoxes();
        const SAMRAI::tbox::Array<int>& local_ids = level->getProcessorMapping().getLocalIndices();
        for (int k = 0; k < local_ids.size(); ++k)
        {
            const int local_id = local_ids[k];
            d_patch_bc_box_overlap[ln][local_id].resize(NDIM);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                const SAMRAI::hier::Box<NDIM> side_box = SAMRAI::pdat::SideGeometry<NDIM>::toSideBox(box_array[local_id],axis);
                const SAMRAI::hier::Box<NDIM> side_ghost_box = SAMRAI::hier::Box<NDIM>::grow(side_box, 1);
                d_patch_bc_box_overlap[ln][local_id][axis] = SAMRAI::hier::BoxList<NDIM>(side_ghost_box);
                d_patch_bc_box_overlap[ln][local_id][axis].removeIntersections(side_box);
            }
        }
    }

    // Indicate that the operator is initialized.
    d_is_initialized = true;

    d_in_initialize_operator_state = false;
    t_initialize_operator_state->stop();
    return;
}// initializeOperatorState

void
INSStaggeredBoxRelaxationFACOperator::deallocateOperatorState()
{
    if (d_is_initialized && !d_in_initialize_operator_state &&
        (d_coarsest_reset_ln != -1) && (d_finest_reset_ln != -1))
    {
        return;
    }

    t_deallocate_operator_state->start();

    if (d_is_initialized)
    {
        const int coarsest_reset_ln =
            (d_in_initialize_operator_state &&
             (d_coarsest_reset_ln != -1) && (d_finest_reset_ln != -1))
            ? d_coarsest_reset_ln : d_coarsest_ln;
        const int finest_reset_ln =
            (d_in_initialize_operator_state &&
             (d_coarsest_reset_ln != -1) && (d_finest_reset_ln != -1))
            ? d_finest_reset_ln : d_finest_ln;

        for (int ln = coarsest_reset_ln; ln <= std::min(d_finest_ln, finest_reset_ln); ++ln)
        {
            if (ln <= d_hierarchy->getFinestLevelNumber())
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                if (level->checkAllocated(d_side_scratch_idx)) level->deallocatePatchData(d_side_scratch_idx);
            }

            int ierr;

            for (std::map<int,std::vector<Vec> >::iterator it = d_patch_vec_e[ln].begin();
                 it != d_patch_vec_e[ln].end(); ++it)
            {
                std::vector<Vec>& e = (*it).second;
                for (unsigned k = 0; k < e.size(); ++k)
                {
                    ierr = VecDestroy(e[k]);  IBTK_CHKERRQ(ierr);
                }
            }
            d_patch_vec_e[ln].clear();
            for (std::map<int,std::vector<Vec> >::iterator it = d_patch_vec_f[ln].begin();
                 it != d_patch_vec_f[ln].end(); ++it)
            {
                std::vector<Vec>& f = (*it).second;
                for (unsigned k = 0; k < f.size(); ++k)
                {
                    ierr = VecDestroy(f[k]);  IBTK_CHKERRQ(ierr);
                }
            }
            d_patch_vec_f[ln].clear();
            for (std::map<int,std::vector<Mat> >::iterator it = d_patch_mat[ln].begin();
                 it != d_patch_mat[ln].end(); ++it)
            {
                std::vector<Mat>& A = (*it).second;
                for (unsigned k = 0; k < A.size(); ++k)
                {
                    ierr = MatDestroy(A[k]);  IBTK_CHKERRQ(ierr);
                }
            }
            d_patch_mat[ln].clear();
         }

        // Delete the solution and rhs vectors.
        d_solution->resetLevels(d_solution->getCoarsestLevelNumber(), std::min(d_solution->getFinestLevelNumber(),d_hierarchy->getFinestLevelNumber()));
        d_solution->freeVectorComponents();
        d_solution.setNull();

        d_rhs->resetLevels(d_rhs->getCoarsestLevelNumber(), std::min(d_rhs->getFinestLevelNumber(),d_hierarchy->getFinestLevelNumber()));
        d_rhs->freeVectorComponents();
        d_rhs.setNull();

        if (!d_in_initialize_operator_state ||
            (d_coarsest_reset_ln == -1) || (d_finest_reset_ln == -1))
        {
            d_patch_vec_e.resize(0);
            d_patch_vec_f.resize(0);
            d_patch_mat.resize(0);
            d_patch_bc_box_overlap.resize(0);

            d_hierarchy.setNull();
            d_coarsest_ln = -1;
            d_finest_ln   = -1;

            d_U_prolongation_refine_operator    .setNull();
            d_U_cf_bdry_op                      .setNull();
            d_prolongation_refine_patch_strategy.setNull();
            d_prolongation_refine_algorithm     .setNull();
            d_prolongation_refine_schedules     .resize(0);

            d_U_urestriction_coarsen_operator.setNull();
            d_urestriction_coarsen_algorithm .setNull();
            d_urestriction_coarsen_schedules .resize(0);

            d_U_rrestriction_coarsen_operator.setNull();
            d_rrestriction_coarsen_algorithm .setNull();
            d_rrestriction_coarsen_schedules .resize(0);

            d_U_ghostfill_refine_operator.setNull();
            d_ghostfill_refine_algorithm .setNull();
            d_ghostfill_refine_schedules .resize(0);

            d_U_ghostfill_nocoarse_refine_operator.setNull();
            d_ghostfill_nocoarse_refine_algorithm .setNull();
            d_ghostfill_nocoarse_refine_schedules .resize(0);

            d_U_synch_refine_operator .setNull();
            d_U_synch_refine_algorithm.setNull();
            d_U_synch_refine_schedules.resize(0);
        }
    }

    // Clear the "reset level" range.
    d_coarsest_reset_ln = -1;
    d_finest_reset_ln   = -1;

    // Indicate that the operator is not initialized.
    d_is_initialized = false;

    t_deallocate_operator_state->stop();
    return;
}// deallocateOperatorState

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
INSStaggeredBoxRelaxationFACOperator::xeqScheduleProlongation(
    const int dst_idx,
    const int src_idx,
    const int dst_ln,
    const bool homogeneous_bc)
{
    d_U_bc_op->setPatchDataIndex(dst_idx);
    d_U_bc_op->setPhysicalBcCoefs(d_U_bc_coefs);
    d_U_bc_op->setHomogeneousBc(homogeneous_bc);
    d_U_cf_bdry_op->setPatchDataIndex(dst_idx);

    SAMRAI::xfer::RefineAlgorithm<NDIM> refiner;
    refiner.registerRefine(dst_idx, src_idx, dst_idx, d_U_prolongation_refine_operator, d_U_op_stencil_fill_pattern);
    refiner.resetSchedule(d_prolongation_refine_schedules[dst_ln]);
    d_prolongation_refine_schedules[dst_ln]->fillData(d_new_time);
    d_prolongation_refine_algorithm->resetSchedule(d_prolongation_refine_schedules[dst_ln]);
    return;
}// xeqScheduleProlongation

void
INSStaggeredBoxRelaxationFACOperator::xeqScheduleURestriction(
    const int dst_idx,
    const int src_idx,
    const int dst_ln)
{
    SAMRAI::xfer::CoarsenAlgorithm<NDIM> coarsener;
    coarsener.registerCoarsen(dst_idx, src_idx, d_U_urestriction_coarsen_operator);
    coarsener.resetSchedule(d_urestriction_coarsen_schedules[dst_ln]);
    d_urestriction_coarsen_schedules[dst_ln]->coarsenData();
    d_urestriction_coarsen_algorithm->resetSchedule(d_urestriction_coarsen_schedules[dst_ln]);
    return;
}// xeqScheduleURestriction

void
INSStaggeredBoxRelaxationFACOperator::xeqScheduleRRestriction(
    const int dst_idx,
    const int src_idx,
    const int dst_ln)
{
    SAMRAI::xfer::CoarsenAlgorithm<NDIM> coarsener;
    coarsener.registerCoarsen(dst_idx, src_idx, d_U_rrestriction_coarsen_operator);
    coarsener.resetSchedule(d_rrestriction_coarsen_schedules[dst_ln]);
    d_rrestriction_coarsen_schedules[dst_ln]->coarsenData();
    d_rrestriction_coarsen_algorithm->resetSchedule(d_rrestriction_coarsen_schedules[dst_ln]);
    return;
}// xeqScheduleRRestriction

void
INSStaggeredBoxRelaxationFACOperator::xeqScheduleGhostFill(
    const int dst_idx,
    const int dst_ln,
    const bool homogeneous_bc)
{
    d_U_bc_op->setPatchDataIndex(dst_idx);
    d_U_bc_op->setPhysicalBcCoefs(d_U_bc_coefs);
    d_U_bc_op->setHomogeneousBc(homogeneous_bc);

    SAMRAI::xfer::RefineAlgorithm<NDIM> refiner;
    refiner.registerRefine(dst_idx, dst_idx, dst_idx, d_U_ghostfill_refine_operator, d_U_op_stencil_fill_pattern);
    refiner.resetSchedule(d_ghostfill_refine_schedules[dst_ln]);
    d_ghostfill_refine_schedules[dst_ln]->fillData(d_new_time);
    d_ghostfill_refine_algorithm->resetSchedule(d_ghostfill_refine_schedules[dst_ln]);
    return;
}// xeqScheduleGhostFill

void
INSStaggeredBoxRelaxationFACOperator::xeqScheduleGhostFillNoCoarse(
    const int dst_idx,
    const int dst_ln,
    const bool homogeneous_bc)
{
    d_U_bc_op->setPatchDataIndex(dst_idx);
    d_U_bc_op->setPhysicalBcCoefs(d_U_bc_coefs);
    d_U_bc_op->setHomogeneousBc(homogeneous_bc);

    SAMRAI::xfer::RefineAlgorithm<NDIM> refiner;
    refiner.registerRefine(dst_idx, dst_idx, dst_idx, d_U_ghostfill_nocoarse_refine_operator, d_U_op_stencil_fill_pattern);
    refiner.resetSchedule(d_ghostfill_nocoarse_refine_schedules[dst_ln]);
    d_ghostfill_nocoarse_refine_schedules[dst_ln]->fillData(d_new_time);
    d_ghostfill_nocoarse_refine_algorithm->resetSchedule(d_ghostfill_nocoarse_refine_schedules[dst_ln]);
    return;
}// xeqScheduleGhostFillNoCoarse

void
INSStaggeredBoxRelaxationFACOperator::xeqScheduleSideDataSynch(
    const int dst_idx,
    const int dst_ln)
{
    SAMRAI::xfer::RefineAlgorithm<NDIM> refiner;
    refiner.registerRefine(dst_idx, dst_idx, dst_idx, d_U_synch_refine_operator, d_U_synch_fill_pattern);
    refiner.resetSchedule(d_U_synch_refine_schedules[dst_ln]);
    d_U_synch_refine_schedules[dst_ln]->fillData(d_new_time);
    d_U_synch_refine_algorithm->resetSchedule(d_U_synch_refine_schedules[dst_ln]);
    return;
}// xeqScheduleSideDataSynch

void
INSStaggeredBoxRelaxationFACOperator::buildPatchOperator(
    Mat& A,
    const INSCoefs& problem_coefs,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
    const int component_axis,
    const SAMRAI::hier::IntVector<NDIM>& ghost_cell_width)
{
    const double C = 0.0; // XXXX
    const double D = 0.0;

    int ierr;

    // Allocate a PETSc matrix for the patch operator.
    const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
    const SAMRAI::hier::Box<NDIM> side_box = SAMRAI::pdat::SideGeometry<NDIM>::toSideBox(patch_box, component_axis);
    const SAMRAI::hier::Box<NDIM> ghost_box = SAMRAI::hier::Box<NDIM>::grow(side_box, ghost_cell_width);
    const int size = ghost_box.size();
    SAMRAI::hier::IntVector<NDIM> component_axis_shift = 0;
    component_axis_shift(component_axis) = -1;

    static const int stencil_sz = 2*NDIM+1;

    SAMRAI::hier::BoxList<NDIM> ghost_boxes(ghost_box);
    ghost_boxes.removeIntersections(side_box);
    std::vector<int> nnz(size, stencil_sz);
    for (SAMRAI::hier::BoxList<NDIM>::Iterator bl(ghost_boxes); bl; bl++)
    {
        for (SAMRAI::hier::Box<NDIM>::Iterator b(bl()); b; b++)
        {
            nnz[ghost_box.offset(b())] = 1;
        }
    }
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, size, size, PETSC_DEFAULT, &nnz[0], &A);  IBTK_CHKERRQ(ierr);

    // Set some general matrix options.
#ifdef DEBUG_CHECK_ASSERTIONS
    ierr = MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);    IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);  IBTK_CHKERRQ(ierr);
#endif

    // Setup the finite difference stencil.  The stencil order is chosen to
    // optimize performance when setting the matrix coefficients.
    static const int x_axis = 0; (void) x_axis;
    static const int y_axis = 1; (void) y_axis;
    static const int z_axis = 2; (void) z_axis;
    std::vector<int> num_cells(NDIM);
    for (int d = 0; d < NDIM; ++d)
    {
        num_cells[d] = ghost_box.numberCells(d);
    }
    std::vector<int> mat_stencil(stencil_sz);
#if (NDIM == 2)
    mat_stencil[0] = -num_cells[x_axis]; // ylower
    mat_stencil[1] = -1;                 // xlower
    mat_stencil[2] = 0;
    mat_stencil[3] = +1;                 // xupper
    mat_stencil[4] = +num_cells[x_axis]; // yupper
#endif
#if (NDIM == 3)
    mat_stencil[0] = -num_cells[x_axis]*num_cells[y_axis]; // zlower
    mat_stencil[1] = -num_cells[x_axis];                   // ylower
    mat_stencil[2] = -1;                                   // xlower
    mat_stencil[3] = 0;
    mat_stencil[4] = +1;                                   // xupper
    mat_stencil[5] = +num_cells[x_axis];                   // yupper
    mat_stencil[6] = +num_cells[x_axis]*num_cells[y_axis]; // zupper
#endif

    // Set the matrix coefficients to correspond to the standard finite
    // difference stencil for the Laplace operator.
    //
    // Note that boundary conditions at both physical boundaries and at
    // coarse-fine interfaces are implicitly treated by setting ghost cell
    // values appropriately.  Thus the matrix coefficients are independent of
    // any boundary conditions.
    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    for (SAMRAI::hier::Box<NDIM>::Iterator b(side_box); b; b++)
    {
        const SAMRAI::hier::Index<NDIM>& i = b();

        std::vector<double> mat_vals(stencil_sz,0.0);
        mat_vals[NDIM] = C;
        for (int axis = 0; axis < NDIM; ++axis)
        {
            const double& h = dx[axis];
            {
                mat_vals[NDIM-axis-1] += D/(h*h);
                mat_vals[NDIM       ] -= D/(h*h);
            }
            {
                mat_vals[NDIM+axis+1] += D/(h*h);
                mat_vals[NDIM       ] -= D/(h*h);
            }
        }

        static const int m = 1;
        static const int n = stencil_sz;
        std::vector<int> idxn(stencil_sz);
        const int idxm = ghost_box.offset(i);

        std::transform(mat_stencil.begin(), mat_stencil.end(), idxn.begin(), std::bind2nd(std::plus<int>(), idxm));
        ierr = MatSetValues(A, m, &idxm, n, &idxn[0], &mat_vals[0], INSERT_VALUES);  IBTK_CHKERRQ(ierr);
    }

    // Set the entries in the ghost cell region so that ghost cell values are
    // not modified by the smoother.
    for (SAMRAI::hier::BoxList<NDIM>::Iterator bl(ghost_boxes); bl; bl++)
    {
        for (SAMRAI::hier::Box<NDIM>::Iterator b(bl()); b; b++)
        {
            const int i = ghost_box.offset(b());
            ierr = MatSetValue(A, i, i, 1.0, INSERT_VALUES);  IBTK_CHKERRQ(ierr);
        }
    }

    // Assemble the matrices.
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    return;
}// buildPatchOperator

void
INSStaggeredBoxRelaxationFACOperator::sanityCheck()
{
    if (d_gcw.min() <= 0)
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  ghost_cell_width.min() must be greater than zero" << std::endl);
    }

    if (d_smoother_choice != "additive")
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  unknown smoother type: " << d_smoother_choice << "\n"
                   << "  valid choices are: additive" << std::endl);
    }

    if (d_coarse_solver_choice != "block_jacobi" &&
        d_coarse_solver_choice != "hypre")
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  unknown coarse solver type: " << d_coarse_solver_choice << "\n"
                   << "  valid choices are: block_jacobi, hypre" << std::endl);
    }

    if (d_coarse_solver_tol < 0.0)
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  invalid coarse solver tolerance: " << d_coarse_solver_tol << std::endl);
    }

    if (d_coarse_solver_max_its <= 0)
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  invalid coarse solver maximum iterations: " << d_coarse_solver_max_its << std::endl);
    }

    for (unsigned l = 0; l < d_U_bc_coefs.size(); ++l)
    {
        if (d_U_bc_coefs[l] == NULL)
        {
            TBOX_ERROR(d_object_name << ":\n"
                       << "  invalid velocity physical bc object at depth = " << l << std::endl);
        }
    }

    if (d_P_bc_coef == NULL)
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  invalid pressure physical bc object" << std::endl);
    }

    if (d_fac_max_cycles <= 0)
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  invalid value for fac_max_cycles: " << d_fac_max_cycles << std::endl);
    }
    return;
}// sanityCheck

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::INSStaggeredBoxRelaxationFACOperator>;

//////////////////////////////////////////////////////////////////////////////
