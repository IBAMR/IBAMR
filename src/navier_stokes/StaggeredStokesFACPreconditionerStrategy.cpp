// Filename: StaggeredStokesFACPreconditionerStrategy.cpp
// Created on 18 Apr 2012 by Boyce Griffith
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
#include <utility>
#include <vector>

#include "CartesianGridGeometry.h"
#include "CellVariable.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "CoarsenSchedule.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
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
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "VariableFillPattern.h"
#include "ibamr/StaggeredStokesFACPreconditionerStrategy.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibamr/StaggeredStokesSolverManager.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CartCellDoubleCubicCoarsen.h"
#include "ibtk/CartCellDoubleQuadraticCFInterpolation.h"
#include "ibtk/CartCellRobinPhysBdryOp.h"
#include "ibtk/CartSideDoubleCubicCoarsen.h"
#include "ibtk/CartSideDoubleQuadraticCFInterpolation.h"
#include "ibtk/CartSideRobinPhysBdryOp.h"
#include "ibtk/CellNoCornersFillPattern.h"
#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"
#include "ibtk/FACPreconditionerStrategy.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/PETScKrylovLinearSolver.h"
#include "ibtk/PETScLevelSolver.h"
#include "ibtk/RefinePatchStrategySet.h"
#include "ibtk/SideNoCornersFillPattern.h"
#include "ibtk/SideSynchCopyFillPattern.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE = "NONE";
static const bool USE_CF_INTERPOLATION = true;
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries; used only to evaluate
// composite grid residuals.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells; used only to evaluate composite grid residuals.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Timers.
static Timer* t_restrict_residual;
static Timer* t_prolong_error;
static Timer* t_prolong_error_and_correct;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesFACPreconditionerStrategy::StaggeredStokesFACPreconditionerStrategy(
    const std::string& object_name,
    const int ghost_cell_width,
    const Pointer<Database> input_db,
    const std::string& default_options_prefix)
    : FACPreconditionerStrategy(object_name),
      d_U_problem_coefs(object_name + "::U_problem_coefs"),
      d_default_U_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_U_bc_coef", Pointer<Database>(NULL))),
      d_U_bc_coefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_U_bc_coef)),
      d_default_P_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_P_bc_coef", Pointer<Database>(NULL))),
      d_P_bc_coef(d_default_P_bc_coef),
      d_bc_helper(NULL),
      d_gcw(ghost_cell_width),
      d_solution(NULL),
      d_rhs(NULL),
      d_hierarchy(),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_level_bdry_fill_ops(),
      d_level_math_ops(),
      d_in_initialize_operator_state(false),
      d_coarsest_reset_ln(-1),
      d_finest_reset_ln(-1),
      d_smoother_type("ADDITIVE"),
      d_U_prolongation_method("CONSTANT_REFINE"),
      d_P_prolongation_method("LINEAR_REFINE"),
      d_U_restriction_method("CONSERVATIVE_COARSEN"),
      d_P_restriction_method("CONSERVATIVE_COARSEN"),
      d_coarse_solver_init_subclass(false),
      d_coarse_solver_type("LEVEL_SMOOTHER"),
      d_coarse_solver_default_options_prefix(default_options_prefix + "level_0_"),
      d_coarse_solver_rel_residual_tol(1.0e-5),
      d_coarse_solver_abs_residual_tol(1.0e-50),
      d_coarse_solver_max_iterations(10),
      d_coarse_solver(),
      d_coarse_solver_db(),
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
    // set some default values.
    d_has_velocity_nullspace = false;
    d_has_pressure_nullspace = false;

    // Get values from the input database.
    if (input_db)
    {
        if (input_db->keyExists("smoother_type")) d_smoother_type = input_db->getString("smoother_type");
        if (input_db->keyExists("U_prolongation_method"))
            d_U_prolongation_method = input_db->getString("U_prolongation_method");
        if (input_db->keyExists("P_prolongation_method"))
            d_P_prolongation_method = input_db->getString("P_prolongation_method");
        if (input_db->keyExists("U_restriction_method"))
            d_U_restriction_method = input_db->getString("U_restriction_method");
        if (input_db->keyExists("P_restriction_method"))
            d_P_restriction_method = input_db->getString("P_restriction_method");
        if (input_db->keyExists("coarse_solver_type")) d_coarse_solver_type = input_db->getString("coarse_solver_type");
        if (input_db->keyExists("coarse_solver_rel_residual_tol"))
            d_coarse_solver_rel_residual_tol = input_db->getDouble("coarse_solver_rel_residual_tol");
        if (input_db->keyExists("coarse_solver_abs_residual_tol"))
            d_coarse_solver_abs_residual_tol = input_db->getDouble("coarse_solver_abs_residual_tol");
        if (input_db->keyExists("coarse_solver_max_iterations"))
            d_coarse_solver_max_iterations = input_db->getInteger("coarse_solver_max_iterations");
        if (input_db->isDatabase("coarse_solver_db")) d_coarse_solver_db = input_db->getDatabase("coarse_solver_db");
    }

    // Configure the coarse level solver.
    setCoarseSolverType(d_coarse_solver_type);

    // Setup scratch variables.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(d_object_name + "::CONTEXT");
    const IntVector<NDIM> side_ghosts = d_gcw;
    Pointer<SideVariable<NDIM, double> > side_scratch_var =
        new SideVariable<NDIM, double>(d_object_name + "::side_scratch");
    if (var_db->checkVariableExists(side_scratch_var->getName()))
    {
        side_scratch_var = var_db->getVariable(side_scratch_var->getName());
        d_side_scratch_idx = var_db->mapVariableAndContextToIndex(side_scratch_var, d_context);
        var_db->removePatchDataIndex(d_side_scratch_idx);
    }
    d_side_scratch_idx = var_db->registerVariableAndContext(side_scratch_var, d_context, side_ghosts);
    const IntVector<NDIM> cell_ghosts = d_gcw;
    Pointer<CellVariable<NDIM, double> > cell_scratch_var =
        new CellVariable<NDIM, double>(d_object_name + "::cell_scratch");
    if (var_db->checkVariableExists(cell_scratch_var->getName()))
    {
        cell_scratch_var = var_db->getVariable(cell_scratch_var->getName());
        d_cell_scratch_idx = var_db->mapVariableAndContextToIndex(cell_scratch_var, d_context);
        var_db->removePatchDataIndex(d_cell_scratch_idx);
    }
    d_cell_scratch_idx = var_db->registerVariableAndContext(cell_scratch_var, d_context, cell_ghosts);

    // Setup Timers.
    IBAMR_DO_ONCE(
        t_restrict_residual =
            TimerManager::getManager()->getTimer("StaggeredStokesFACPreconditionerStrategy::restrictResidual()");
        t_prolong_error =
            TimerManager::getManager()->getTimer("StaggeredStokesFACPreconditionerStrategy::prolongError()");
        t_prolong_error_and_correct =
            TimerManager::getManager()->getTimer("StaggeredStokesFACPreconditionerStrategy::prolongErrorAndCorrect()");
        t_initialize_operator_state =
            TimerManager::getManager()->getTimer("StaggeredStokesFACPreconditionerStrategy::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer(
            "StaggeredStokesFACPreconditionerStrategy::deallocateOperatorState()"););
    return;
} // StaggeredStokesFACPreconditionerStrategy

StaggeredStokesFACPreconditionerStrategy::~StaggeredStokesFACPreconditionerStrategy()
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::~StaggeredStokesFACPreconditionerStrategy()\n"
                                 << "  subclass must call deallocateOperatorState in subclass destructor"
                                 << std::endl);
    }
    delete d_default_U_bc_coef;
    d_default_U_bc_coef = NULL;
    delete d_default_P_bc_coef;
    d_default_P_bc_coef = NULL;
    return;
} // ~StaggeredStokesFACPreconditionerStrategy

void
StaggeredStokesFACPreconditionerStrategy::setVelocityPoissonSpecifications(const PoissonSpecifications& U_problem_coefs)
{
    d_U_problem_coefs = U_problem_coefs;
    return;
} // setVelocityPoissonSpecifications

void
StaggeredStokesFACPreconditionerStrategy::setComponentsHaveNullspace(const bool has_velocity_nullspace,
                                                                     const bool has_pressure_nullspace)
{
    d_has_velocity_nullspace = has_velocity_nullspace;
    d_has_pressure_nullspace = has_pressure_nullspace;

    return;
} // setComponentsHaveNullspace

void
StaggeredStokesFACPreconditionerStrategy::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
                                                             RobinBcCoefStrategy<NDIM>* P_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(U_bc_coefs.size() == NDIM);
#endif
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (U_bc_coefs[d])
        {
            d_U_bc_coefs[d] = U_bc_coefs[d];
        }
        else
        {
            d_U_bc_coefs[d] = d_default_U_bc_coef;
        }
    }

    if (P_bc_coef)
    {
        d_P_bc_coef = P_bc_coef;
    }
    else
    {
        d_P_bc_coef = d_default_P_bc_coef;
    }

    return;
} // setPhysicalBcCoefs

void
StaggeredStokesFACPreconditionerStrategy::setPhysicalBoundaryHelper(
    Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_helper);
#endif
    d_bc_helper = bc_helper;
    return;
} // setPhysicalBoundaryHelper

void
StaggeredStokesFACPreconditionerStrategy::setResetLevels(const int coarsest_ln, const int finest_ln)
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
StaggeredStokesFACPreconditionerStrategy::setSmootherType(const std::string& smoother_type)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setSmootherType()\n"
                                 << "  cannot be called while operator state is initialized"
                                 << std::endl);
    }
    d_smoother_type = smoother_type;
    return;
} // setSmootherType

void
StaggeredStokesFACPreconditionerStrategy::setCoarseSolverType(const std::string& coarse_solver_type)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setCoarseSolverType():\n"
                                 << "  cannot be called while operator state is initialized"
                                 << std::endl);
    }
    if (d_coarse_solver_type != coarse_solver_type) d_coarse_solver.setNull();
    d_coarse_solver_type = coarse_solver_type;
    return;
} // setCoarseSolverType

void
StaggeredStokesFACPreconditionerStrategy::setCoarseSolverMaxIterations(int coarse_solver_max_iterations)
{
    d_coarse_solver_max_iterations = coarse_solver_max_iterations;
    return;
} // setCoarseSolverMaxIterations

void
StaggeredStokesFACPreconditionerStrategy::setCoarseSolverAbsoluteTolerance(double coarse_solver_abs_residual_tol)
{
    d_coarse_solver_abs_residual_tol = coarse_solver_abs_residual_tol;
    return;
} // setCoarseSolverAbsoluteTolerance

void
StaggeredStokesFACPreconditionerStrategy::setCoarseSolverRelativeTolerance(double coarse_solver_rel_residual_tol)
{
    d_coarse_solver_rel_residual_tol = coarse_solver_rel_residual_tol;
    return;
} // setCoarseSolverRelativeTolerance

void
StaggeredStokesFACPreconditionerStrategy::setProlongationMethods(const std::string& U_prolongation_method,
                                                                 const std::string& P_prolongation_method)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setProlongationMethods()\n"
                                 << "  cannot be called while operator state is initialized"
                                 << std::endl);
    }
    d_U_prolongation_method = U_prolongation_method;
    d_P_prolongation_method = P_prolongation_method;
    return;
} // setProlongationMethods

void
StaggeredStokesFACPreconditionerStrategy::setRestrictionMethods(const std::string& U_restriction_method,
                                                                const std::string& P_restriction_method)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setRestrictionMethods()\n"
                                 << "  cannot be called while operator state is initialized"
                                 << std::endl);
    }
    d_U_restriction_method = U_restriction_method;
    d_P_restriction_method = P_restriction_method;
    return;
} // setRestrictionMethods

void
StaggeredStokesFACPreconditionerStrategy::restrictResidual(const SAMRAIVectorReal<NDIM, double>& src,
                                                           SAMRAIVectorReal<NDIM, double>& dst,
                                                           int dst_ln)
{
    IBAMR_TIMER_START(t_restrict_residual);

    const int U_src_idx = src.getComponentDescriptorIndex(0);
    const int P_src_idx = src.getComponentDescriptorIndex(1);
    const std::pair<int, int> src_idxs = std::make_pair(U_src_idx, P_src_idx);

    const int U_dst_idx = dst.getComponentDescriptorIndex(0);
    const int P_dst_idx = dst.getComponentDescriptorIndex(1);
    const std::pair<int, int> dst_idxs = std::make_pair(U_dst_idx, P_dst_idx);

    if (U_src_idx != U_dst_idx)
    {
        HierarchySideDataOpsReal<NDIM, double> level_sc_data_ops(d_hierarchy, dst_ln, dst_ln);
        static const bool interior_only = false;
        level_sc_data_ops.copyData(U_dst_idx, U_src_idx, interior_only);
    }
    if (P_src_idx != P_dst_idx)
    {
        HierarchyCellDataOpsReal<NDIM, double> level_cc_data_ops(d_hierarchy, dst_ln, dst_ln);
        static const bool interior_only = false;
        level_cc_data_ops.copyData(P_dst_idx, P_src_idx, interior_only);
    }
    xeqScheduleRestriction(dst_idxs, src_idxs, dst_ln);

    IBAMR_TIMER_STOP(t_restrict_residual);
    return;
} // restrictResidual

void
StaggeredStokesFACPreconditionerStrategy::prolongError(const SAMRAIVectorReal<NDIM, double>& src,
                                                       SAMRAIVectorReal<NDIM, double>& dst,
                                                       int dst_ln)
{
    IBAMR_TIMER_START(t_prolong_error);

    const int U_src_idx = src.getComponentDescriptorIndex(0);
    const int P_src_idx = src.getComponentDescriptorIndex(1);
    const std::pair<int, int> src_idxs = std::make_pair(U_src_idx, P_src_idx);

    const int U_dst_idx = dst.getComponentDescriptorIndex(0);
    const int P_dst_idx = dst.getComponentDescriptorIndex(1);
    const std::pair<int, int> dst_idxs = std::make_pair(U_dst_idx, P_dst_idx);

    // Refine the correction from the coarse level src data directly into the
    // fine level error.
    xeqScheduleProlongation(dst_idxs, src_idxs, dst_ln);

    IBAMR_TIMER_STOP(t_prolong_error);
    return;
} // prolongError

void
StaggeredStokesFACPreconditionerStrategy::prolongErrorAndCorrect(const SAMRAIVectorReal<NDIM, double>& src,
                                                                 SAMRAIVectorReal<NDIM, double>& dst,
                                                                 int dst_ln)
{
    IBAMR_TIMER_START(t_prolong_error_and_correct);

    const int U_src_idx = src.getComponentDescriptorIndex(0);
    const int P_src_idx = src.getComponentDescriptorIndex(1);
    const std::pair<int, int> src_idxs = std::make_pair(U_src_idx, P_src_idx);

    const int U_dst_idx = dst.getComponentDescriptorIndex(0);
    const int P_dst_idx = dst.getComponentDescriptorIndex(1);

    const std::pair<int, int> scratch_idxs = std::make_pair(d_side_scratch_idx, d_cell_scratch_idx);

    // Prolong the correction from the coarse level src data into the fine level
    // scratch data and then correct the fine level dst data.
    static const bool interior_only = false;
    if (U_src_idx != U_dst_idx)
    {
        HierarchySideDataOpsReal<NDIM, double> level_sc_data_ops_coarse(d_hierarchy, dst_ln - 1, dst_ln - 1);
        level_sc_data_ops_coarse.add(U_dst_idx, U_dst_idx, U_src_idx, interior_only);
    }
    if (P_src_idx != P_dst_idx)
    {
        HierarchyCellDataOpsReal<NDIM, double> level_cc_data_ops_coarse(d_hierarchy, dst_ln - 1, dst_ln - 1);
        level_cc_data_ops_coarse.add(P_dst_idx, P_dst_idx, P_src_idx, interior_only);
    }
    xeqScheduleProlongation(scratch_idxs, src_idxs, dst_ln);
    HierarchySideDataOpsReal<NDIM, double> level_sc_data_ops_fine(d_hierarchy, dst_ln, dst_ln);
    level_sc_data_ops_fine.add(U_dst_idx, U_dst_idx, d_side_scratch_idx, interior_only);
    HierarchyCellDataOpsReal<NDIM, double> level_cc_data_ops_fine(d_hierarchy, dst_ln, dst_ln);
    level_cc_data_ops_fine.add(P_dst_idx, P_dst_idx, d_cell_scratch_idx, interior_only);

    IBAMR_TIMER_STOP(t_prolong_error_and_correct);
    return;
} // prolongErrorAndCorrect

bool
StaggeredStokesFACPreconditionerStrategy::solveCoarsestLevel(SAMRAIVectorReal<NDIM, double>& error,
                                                             const SAMRAIVectorReal<NDIM, double>& residual,
                                                             int coarsest_ln)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(coarsest_ln == d_coarsest_ln);
#endif
    if (!d_coarse_solver)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(d_coarse_solver_type == "LEVEL_SMOOTHER");
#endif
        smoothError(error, residual, coarsest_ln, d_coarse_solver_max_iterations, false, false);
    }
    else
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(d_coarse_solver);
#endif
        d_coarse_solver->setSolutionTime(d_solution_time);
        d_coarse_solver->setTimeInterval(d_current_time, d_new_time);
        d_coarse_solver->setMaxIterations(d_coarse_solver_max_iterations);
        d_coarse_solver->setAbsoluteTolerance(d_coarse_solver_abs_residual_tol);
        d_coarse_solver->setRelativeTolerance(d_coarse_solver_rel_residual_tol);
        d_coarse_solver->setComponentsHaveNullspace(d_has_velocity_nullspace, d_has_pressure_nullspace);
        LinearSolver* p_coarse_solver = dynamic_cast<LinearSolver*>(d_coarse_solver.getPointer());

        if (p_coarse_solver)
        {
            bool initial_guess_nonzero = true;
            PETScKrylovLinearSolver* p_petsc_solver = dynamic_cast<PETScKrylovLinearSolver*>(p_coarse_solver);
            PETScLevelSolver* p_petsc_level_solver = dynamic_cast<PETScLevelSolver*>(p_coarse_solver);
            if (p_petsc_solver || p_petsc_level_solver)
            {
                const KSP& petsc_ksp =
                    p_petsc_solver ? p_petsc_solver->getPETScKSP() : p_petsc_level_solver->getPETScKSP();
                KSPType ksp_type;
                KSPGetType(petsc_ksp, &ksp_type);

                if (!strcmp(ksp_type, "preonly")) initial_guess_nonzero = false;
            }
            p_coarse_solver->setInitialGuessNonzero(initial_guess_nonzero);
        }

        d_coarse_solver->solveSystem(*getLevelSAMRAIVectorReal(error, d_coarsest_ln),
                                     *getLevelSAMRAIVectorReal(residual, d_coarsest_ln));
    }
    return true;
} // solveCoarsestLevel

void
StaggeredStokesFACPreconditionerStrategy::computeResidual(SAMRAIVectorReal<NDIM, double>& residual,
                                                          const SAMRAIVectorReal<NDIM, double>& solution,
                                                          const SAMRAIVectorReal<NDIM, double>& rhs,
                                                          int coarsest_level_num,
                                                          int finest_level_num)
{
    const int U_res_idx = residual.getComponentDescriptorIndex(0);
    const int U_sol_idx = solution.getComponentDescriptorIndex(0);
    const int U_rhs_idx = rhs.getComponentDescriptorIndex(0);

    const Pointer<SideVariable<NDIM, double> > U_res_sc_var = residual.getComponentVariable(0);
    const Pointer<SideVariable<NDIM, double> > U_sol_sc_var = solution.getComponentVariable(0);
    const Pointer<SideVariable<NDIM, double> > U_rhs_sc_var = rhs.getComponentVariable(0);

    const int P_res_idx = residual.getComponentDescriptorIndex(1);
    const int P_sol_idx = solution.getComponentDescriptorIndex(1);
    const int P_rhs_idx = rhs.getComponentDescriptorIndex(1);

    const Pointer<CellVariable<NDIM, double> > P_res_cc_var = residual.getComponentVariable(1);
    const Pointer<CellVariable<NDIM, double> > P_sol_cc_var = solution.getComponentVariable(1);
    const Pointer<CellVariable<NDIM, double> > P_rhs_cc_var = rhs.getComponentVariable(1);

    // Fill ghost-cell values.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    Pointer<VariableFillPattern<NDIM> > sc_fill_pattern = new SideNoCornersFillPattern(d_gcw, false, false, true);
    Pointer<VariableFillPattern<NDIM> > cc_fill_pattern = new CellNoCornersFillPattern(d_gcw, false, false, true);
    InterpolationTransactionComponent U_scratch_component(U_sol_idx,
                                                          DATA_REFINE_TYPE,
                                                          USE_CF_INTERPOLATION,
                                                          DATA_COARSEN_TYPE,
                                                          BDRY_EXTRAP_TYPE,
                                                          CONSISTENT_TYPE_2_BDRY,
                                                          d_U_bc_coefs,
                                                          sc_fill_pattern);
    InterpolationTransactionComponent P_scratch_component(P_sol_idx,
                                                          DATA_REFINE_TYPE,
                                                          USE_CF_INTERPOLATION,
                                                          DATA_COARSEN_TYPE,
                                                          BDRY_EXTRAP_TYPE,
                                                          CONSISTENT_TYPE_2_BDRY,
                                                          d_P_bc_coef,
                                                          cc_fill_pattern);
    std::vector<InterpolationTransactionComponent> U_P_components(2);
    U_P_components[0] = U_scratch_component;
    U_P_components[1] = P_scratch_component;
    if (!d_level_bdry_fill_ops[finest_level_num])
    {
        d_level_bdry_fill_ops[finest_level_num] = new HierarchyGhostCellInterpolation();
        d_level_bdry_fill_ops[finest_level_num]->initializeOperatorState(
            U_P_components, d_hierarchy, coarsest_level_num, finest_level_num);
    }
    else
    {
        d_level_bdry_fill_ops[finest_level_num]->resetTransactionComponents(U_P_components);
    }
    d_level_bdry_fill_ops[finest_level_num]->setHomogeneousBc(true);
    d_level_bdry_fill_ops[finest_level_num]->fillData(d_new_time);

    InterpolationTransactionComponent default_U_scratch_component(d_solution->getComponentDescriptorIndex(0),
                                                                  DATA_REFINE_TYPE,
                                                                  USE_CF_INTERPOLATION,
                                                                  DATA_COARSEN_TYPE,
                                                                  BDRY_EXTRAP_TYPE,
                                                                  CONSISTENT_TYPE_2_BDRY,
                                                                  d_U_bc_coefs,
                                                                  sc_fill_pattern);

    InterpolationTransactionComponent default_P_scratch_component(d_solution->getComponentDescriptorIndex(1),
                                                                  DATA_REFINE_TYPE,
                                                                  USE_CF_INTERPOLATION,
                                                                  DATA_COARSEN_TYPE,
                                                                  BDRY_EXTRAP_TYPE,
                                                                  CONSISTENT_TYPE_2_BDRY,
                                                                  d_P_bc_coef,
                                                                  cc_fill_pattern);
    std::vector<InterpolationTransactionComponent> default_U_P_components(2);
    default_U_P_components[0] = default_U_scratch_component;
    default_U_P_components[1] = default_P_scratch_component;
    d_level_bdry_fill_ops[finest_level_num]->resetTransactionComponents(default_U_P_components);

    // Compute the residual, r = f - A*u.
    if (!d_level_math_ops[finest_level_num])
    {
        std::ostringstream stream;
        stream << d_object_name << "::level_math_ops_" << finest_level_num;
        d_level_math_ops[finest_level_num] =
            new HierarchyMathOps(stream.str(), d_hierarchy, coarsest_level_num, finest_level_num);
    }
    d_level_math_ops[finest_level_num]->grad(U_res_idx,
                                             U_res_sc_var,
                                             /*cf_bdry_synch*/ true,
                                             1.0,
                                             P_sol_idx,
                                             P_sol_cc_var,
                                             NULL,
                                             d_new_time);
    d_level_math_ops[finest_level_num]->laplace(U_res_idx,
                                                U_res_sc_var,
                                                d_U_problem_coefs,
                                                U_sol_idx,
                                                U_sol_sc_var,
                                                NULL,
                                                d_new_time,
                                                1.0,
                                                U_res_idx,
                                                U_res_sc_var);
    HierarchySideDataOpsReal<NDIM, double> level_sc_data_ops(d_hierarchy, coarsest_level_num, finest_level_num);
    level_sc_data_ops.axpy(U_res_idx, -1.0, U_res_idx, U_rhs_idx, false);
    d_level_math_ops[finest_level_num]->div(P_res_idx,
                                            P_res_cc_var,
                                            -1.0,
                                            U_sol_idx,
                                            U_sol_sc_var,
                                            NULL,
                                            d_new_time,
                                            /*cf_bdry_synch*/ true);
    HierarchyCellDataOpsReal<NDIM, double> level_cc_data_ops(d_hierarchy, coarsest_level_num, finest_level_num);
    level_cc_data_ops.axpy(P_res_idx, -1.0, P_res_idx, P_rhs_idx, false);
    return;
} // computeResidual

void
StaggeredStokesFACPreconditionerStrategy::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& solution,
                                                                  const SAMRAIVectorReal<NDIM, double>& rhs)
{
    IBAMR_TIMER_START(t_initialize_operator_state);

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

    // Reset the hierarchy configuration.
    d_hierarchy = solution.getPatchHierarchy();
    d_coarsest_ln = solution.getCoarsestLevelNumber();
    d_finest_ln = solution.getFinestLevelNumber();

    // Setup boundary condition handling objects.
    d_U_bc_op = new CartSideRobinPhysBdryOp(d_side_scratch_idx, d_U_bc_coefs, false);
    d_P_bc_op = new CartCellRobinPhysBdryOp(d_cell_scratch_idx, d_P_bc_coef, false);
    d_U_cf_bdry_op = new CartSideDoubleQuadraticCFInterpolation();
    d_P_cf_bdry_op = new CartCellDoubleQuadraticCFInterpolation();
    d_U_op_stencil_fill_pattern = new SideNoCornersFillPattern(d_gcw, false, false, false);
    d_P_op_stencil_fill_pattern = new CellNoCornersFillPattern(d_gcw, false, false, false);
    d_U_synch_fill_pattern = new SideSynchCopyFillPattern();

    // Initialize the coarse level solvers when needed.
    if (!d_coarse_solver_init_subclass && coarsest_reset_ln == d_coarsest_ln &&
        d_coarse_solver_type != "LEVEL_SMOOTHER")
    {
        if (!d_coarse_solver)
        {
            d_coarse_solver =
                StaggeredStokesSolverManager::getManager()->allocateSolver(d_coarse_solver_type,
                                                                           d_object_name + "::coarse_solver",
                                                                           d_coarse_solver_db,
                                                                           d_coarse_solver_default_options_prefix);
        }
#if !defined(NDEBUG)
        TBOX_ASSERT(d_coarse_solver);
#endif
        // Note that since the coarse level solver is solving for the error, it
        // must always employ homogeneous boundary conditions.
        d_coarse_solver->setSolutionTime(d_solution_time);
        d_coarse_solver->setTimeInterval(d_current_time, d_new_time);
        d_coarse_solver->setVelocityPoissonSpecifications(d_U_problem_coefs);
        d_coarse_solver->setPhysicalBcCoefs(d_U_bc_coefs, d_P_bc_coef);
        d_coarse_solver->setPhysicalBoundaryHelper(d_bc_helper);
        d_coarse_solver->setHomogeneousBc(true);
        d_coarse_solver->setComponentsHaveNullspace(d_has_velocity_nullspace, d_has_pressure_nullspace);
        d_coarse_solver->initializeSolverState(*getLevelSAMRAIVectorReal(*d_solution, d_coarsest_ln),
                                               *getLevelSAMRAIVectorReal(*d_rhs, d_coarsest_ln));
    }

    // Perform implementation-specific initialization.
    initializeOperatorStateSpecialized(solution, rhs, coarsest_reset_ln, finest_reset_ln);

    // Setup level operators.
    d_level_bdry_fill_ops.resize(d_finest_ln + 1, NULL);
    d_level_math_ops.resize(d_finest_ln + 1, NULL);
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln); ln <= finest_reset_ln; ++ln)
    {
        d_level_bdry_fill_ops[ln].setNull();
        d_level_math_ops[ln].setNull();
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
    IBAMR_DO_ONCE(geometry->addSpatialCoarsenOperator(new CartSideDoubleCubicCoarsen());
                  geometry->addSpatialCoarsenOperator(new CartCellDoubleCubicCoarsen()););
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

    d_prolongation_refine_schedules.resize(d_finest_ln + 1);
    d_restriction_coarsen_schedules.resize(d_finest_ln + 1);
    d_ghostfill_nocoarse_refine_schedules.resize(d_finest_ln + 1);
    d_synch_refine_schedules.resize(d_finest_ln + 1);

    d_prolongation_refine_algorithm = new RefineAlgorithm<NDIM>();
    d_restriction_coarsen_algorithm = new CoarsenAlgorithm<NDIM>();
    d_ghostfill_nocoarse_refine_algorithm = new RefineAlgorithm<NDIM>();
    d_synch_refine_algorithm = new RefineAlgorithm<NDIM>();

    d_prolongation_refine_algorithm->registerRefine(d_side_scratch_idx,
                                                    solution.getComponentDescriptorIndex(0),
                                                    d_side_scratch_idx,
                                                    d_U_prolongation_refine_operator,
                                                    d_U_op_stencil_fill_pattern);
    d_prolongation_refine_algorithm->registerRefine(d_cell_scratch_idx,
                                                    solution.getComponentDescriptorIndex(1),
                                                    d_cell_scratch_idx,
                                                    d_P_prolongation_refine_operator,
                                                    d_P_op_stencil_fill_pattern);

    d_restriction_coarsen_algorithm->registerCoarsen(
        d_side_scratch_idx, rhs.getComponentDescriptorIndex(0), d_U_restriction_coarsen_operator);
    d_restriction_coarsen_algorithm->registerCoarsen(
        d_cell_scratch_idx, rhs.getComponentDescriptorIndex(1), d_P_restriction_coarsen_operator);

    d_ghostfill_nocoarse_refine_algorithm->registerRefine(solution.getComponentDescriptorIndex(0),
                                                          solution.getComponentDescriptorIndex(0),
                                                          solution.getComponentDescriptorIndex(0),
                                                          Pointer<RefineOperator<NDIM> >(),
                                                          d_U_op_stencil_fill_pattern);
    d_ghostfill_nocoarse_refine_algorithm->registerRefine(solution.getComponentDescriptorIndex(1),
                                                          solution.getComponentDescriptorIndex(1),
                                                          solution.getComponentDescriptorIndex(1),
                                                          Pointer<RefineOperator<NDIM> >(),
                                                          d_P_op_stencil_fill_pattern);

    d_synch_refine_algorithm->registerRefine(solution.getComponentDescriptorIndex(0),
                                             solution.getComponentDescriptorIndex(0),
                                             solution.getComponentDescriptorIndex(0),
                                             Pointer<RefineOperator<NDIM> >(),
                                             d_U_synch_fill_pattern);

    std::vector<RefinePatchStrategy<NDIM>*> bc_op_ptrs(2);
    bc_op_ptrs[0] = d_U_bc_op;
    bc_op_ptrs[1] = d_P_bc_op;
    d_U_P_bc_op = new RefinePatchStrategySet(bc_op_ptrs.begin(), bc_op_ptrs.end(), false);

    for (int dst_ln = d_coarsest_ln + 1; dst_ln <= d_finest_ln; ++dst_ln)
    {
        d_prolongation_refine_schedules[dst_ln] =
            d_prolongation_refine_algorithm->createSchedule(d_hierarchy->getPatchLevel(dst_ln),
                                                            Pointer<PatchLevel<NDIM> >(),
                                                            dst_ln - 1,
                                                            d_hierarchy,
                                                            d_prolongation_refine_patch_strategy.getPointer());

        d_ghostfill_nocoarse_refine_schedules[dst_ln] =
            d_ghostfill_nocoarse_refine_algorithm->createSchedule(d_hierarchy->getPatchLevel(dst_ln), d_U_P_bc_op);

        d_synch_refine_schedules[dst_ln] = d_synch_refine_algorithm->createSchedule(d_hierarchy->getPatchLevel(dst_ln));
    }

    d_ghostfill_nocoarse_refine_schedules[d_coarsest_ln] =
        d_ghostfill_nocoarse_refine_algorithm->createSchedule(d_hierarchy->getPatchLevel(d_coarsest_ln), d_U_P_bc_op);

    d_synch_refine_schedules[d_coarsest_ln] =
        d_synch_refine_algorithm->createSchedule(d_hierarchy->getPatchLevel(d_coarsest_ln));

    for (int dst_ln = d_coarsest_ln; dst_ln < d_finest_ln; ++dst_ln)
    {
        d_restriction_coarsen_schedules[dst_ln] = d_restriction_coarsen_algorithm->createSchedule(
            d_hierarchy->getPatchLevel(dst_ln), d_hierarchy->getPatchLevel(dst_ln + 1));
    }

    // Indicate that the operator is initialized.
    d_is_initialized = true;
    d_in_initialize_operator_state = false;

    IBAMR_TIMER_STOP(t_initialize_operator_state);
    return;
} // initializeOperatorState

void
StaggeredStokesFACPreconditionerStrategy::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_operator_state);

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
        if (level->checkAllocated(d_side_scratch_idx)) level->deallocatePatchData(d_side_scratch_idx);
        if (level->checkAllocated(d_cell_scratch_idx)) level->deallocatePatchData(d_cell_scratch_idx);
    }

    // Delete the solution and rhs vectors.
    d_solution->resetLevels(d_solution->getCoarsestLevelNumber(),
                            std::min(d_solution->getFinestLevelNumber(), d_hierarchy->getFinestLevelNumber()));
    d_solution->freeVectorComponents();
    d_solution.setNull();

    d_rhs->resetLevels(d_rhs->getCoarsestLevelNumber(),
                       std::min(d_rhs->getFinestLevelNumber(), d_hierarchy->getFinestLevelNumber()));
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

        if (d_coarse_solver) d_coarse_solver->deallocateSolverState();

        d_U_prolongation_refine_operator.setNull();
        d_P_prolongation_refine_operator.setNull();
        d_prolongation_refine_patch_strategy.setNull();
        d_prolongation_refine_algorithm.setNull();
        d_prolongation_refine_schedules.resize(0);

        d_U_restriction_coarsen_operator.setNull();
        d_P_restriction_coarsen_operator.setNull();
        d_restriction_coarsen_algorithm.setNull();
        d_restriction_coarsen_schedules.resize(0);

        d_ghostfill_nocoarse_refine_algorithm.setNull();
        d_ghostfill_nocoarse_refine_schedules.resize(0);

        d_synch_refine_algorithm.setNull();
        d_synch_refine_schedules.resize(0);
    }

    delete d_U_P_bc_op;

    // Clear the "reset level" range.
    d_coarsest_reset_ln = -1;
    d_finest_reset_ln = -1;

    // Indicate that the operator is not initialized.
    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_operator_state);
    return;
} // deallocateOperatorState

/////////////////////////////// PROTECTED ////////////////////////////////////

void
StaggeredStokesFACPreconditionerStrategy::xeqScheduleProlongation(const std::pair<int, int>& dst_idxs,
                                                                  const std::pair<int, int>& src_idxs,
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
    refine_alg.registerRefine(
        U_dst_idx, U_src_idx, U_dst_idx, d_U_prolongation_refine_operator, d_U_op_stencil_fill_pattern);
    refine_alg.registerRefine(
        P_dst_idx, P_src_idx, P_dst_idx, d_P_prolongation_refine_operator, d_P_op_stencil_fill_pattern);
    refine_alg.resetSchedule(d_prolongation_refine_schedules[dst_ln]);
    d_prolongation_refine_schedules[dst_ln]->fillData(d_new_time);
    d_prolongation_refine_algorithm->resetSchedule(d_prolongation_refine_schedules[dst_ln]);
    return;
} // xeqScheduleProlongation

void
StaggeredStokesFACPreconditionerStrategy::xeqScheduleRestriction(const std::pair<int, int>& dst_idxs,
                                                                 const std::pair<int, int>& src_idxs,
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
} // xeqScheduleRestriction

void
StaggeredStokesFACPreconditionerStrategy::xeqScheduleGhostFillNoCoarse(const std::pair<int, int>& dst_idxs,
                                                                       const int dst_ln)
{
    const int U_dst_idx = dst_idxs.first;
    d_U_bc_op->setPatchDataIndex(U_dst_idx);
    d_U_bc_op->setHomogeneousBc(true);

    const int P_dst_idx = dst_idxs.second;
    d_P_bc_op->setPatchDataIndex(P_dst_idx);
    d_P_bc_op->setHomogeneousBc(true);

    RefineAlgorithm<NDIM> refine_alg;
    refine_alg.registerRefine(
        U_dst_idx, U_dst_idx, U_dst_idx, Pointer<RefineOperator<NDIM> >(), d_U_op_stencil_fill_pattern);
    refine_alg.registerRefine(
        P_dst_idx, P_dst_idx, P_dst_idx, Pointer<RefineOperator<NDIM> >(), d_P_op_stencil_fill_pattern);
    refine_alg.resetSchedule(d_ghostfill_nocoarse_refine_schedules[dst_ln]);
    d_ghostfill_nocoarse_refine_schedules[dst_ln]->fillData(d_new_time);
    d_ghostfill_nocoarse_refine_algorithm->resetSchedule(d_ghostfill_nocoarse_refine_schedules[dst_ln]);
    return;
} // xeqScheduleGhostFillNoCoarse

void
StaggeredStokesFACPreconditionerStrategy::xeqScheduleDataSynch(const int U_dst_idx, const int dst_ln)
{
    RefineAlgorithm<NDIM> refine_alg;
    refine_alg.registerRefine(
        U_dst_idx, U_dst_idx, U_dst_idx, Pointer<RefineOperator<NDIM> >(), d_U_synch_fill_pattern);
    refine_alg.resetSchedule(d_synch_refine_schedules[dst_ln]);
    d_synch_refine_schedules[dst_ln]->fillData(d_new_time);
    d_synch_refine_algorithm->resetSchedule(d_synch_refine_schedules[dst_ln]);
    return;
} // xeqScheduleDataSynch

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
