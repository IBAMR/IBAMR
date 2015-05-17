// Filename: StaggeredStokesIBLevelRelaxationFACOperator.cpp
// Created on 22 Mar 2015 by Amneet Bhalla
//
// Copyright (c) 2002-2015, Amneet Bhalla and Boyce Griffith
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
#include "ibamr/StaggeredStokesIBLevelRelaxationFACOperator.h"
#include "ibamr/StaggeredStokesPETScLevelSolver.h"
#include "ibamr/StaggeredStokesPETScVecUtilities.h"
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
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/PETScMatUtilities.h"
#include "ibtk/PETScVecUtilities.h"
#include "ibtk/RefinePatchStrategySet.h"
#include "ibtk/SideNoCornersFillPattern.h"
#include "ibtk/SideSynchCopyFillPattern.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#include "petscksp.h"
#include "boost/array.hpp"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{

// Ghost cell width for pressure
static const int SIDEG = 1;
static const int CELLG = 1;
static const int NOGHOST = 0;

// Timers.
static Timer* t_compute_residual;
static Timer* t_restrict_residual;
static Timer* t_prolong_error;
static Timer* t_prolong_error_and_correct;
static Timer* t_smooth_error;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesIBLevelRelaxationFACOperator::StaggeredStokesIBLevelRelaxationFACOperator(
    const std::string& object_name,
    const Pointer<Database> input_db,
    const std::string& default_options_prefix)
    : FACPreconditionerStrategy(object_name), d_U_problem_coefs(object_name + "::U_problem_coefs")
{
    // Set some default values.
    d_default_U_bc_coef =
        new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_U_bc_coef", Pointer<Database>(NULL)),
    d_U_bc_coefs.resize(NDIM, d_default_U_bc_coef);
    d_default_P_bc_coef =
        new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_P_bc_coef", Pointer<Database>(NULL));
    d_P_bc_coef = d_default_P_bc_coef;
    d_coarsest_ln = -1;
    d_finest_ln = -1;
    d_in_initialize_operator_state = false;
    d_coarsest_reset_ln = -1;
    d_finest_reset_ln = -1;
    d_smoother_type = "ADDITIVE";
    d_U_prolongation_method = "CONSTANT_REFINE";
    d_P_prolongation_method = "LINEAR_REFINE";
    d_U_restriction_method = "CONSERVATIVE_COARSEN";
    d_P_restriction_method = "CONSERVATIVE_COARSEN";
    d_coarse_solver_type = "PETSC_LEVEL_SOLVER";
    d_coarse_solver_default_options_prefix = default_options_prefix + "level_0_";
    d_coarse_solver_rel_residual_tol = 1.0e-5;
    d_coarse_solver_abs_residual_tol = 1.0e-50;
    d_coarse_solver_max_iterations = 10;
    d_level_solver_type = "PETSC_LEVEL_SOLVER";
    d_level_solver_default_options_prefix = default_options_prefix + "level_";
    d_level_solver_rel_residual_tol = 1.0e-5;
    d_level_solver_abs_residual_tol = 1.0e-50;
    d_level_solver_max_iterations = 10;
    d_has_velocity_nullspace = false;
    d_has_pressure_nullspace = false;
    d_side_scratch_idx = -1;
    d_cell_scratch_idx = -1;
    d_u_dof_index_idx = -1;
    d_p_dof_index_idx = -1;

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
        if (input_db->keyExists("level_solver_type")) d_level_solver_type = input_db->getString("level_solver_type");
        if (input_db->keyExists("level_solver_rel_residual_tol"))
            d_level_solver_rel_residual_tol = input_db->getDouble("level_solver_rel_residual_tol");
        if (input_db->keyExists("level_solver_abs_residual_tol"))
            d_level_solver_abs_residual_tol = input_db->getDouble("level_solver_abs_residual_tol");
        if (input_db->keyExists("level_solver_max_iterations"))
            d_level_solver_max_iterations = input_db->getInteger("level_solver_max_iterations");
        if (input_db->isDatabase("level_solver_db"))
        {
            d_level_solver_db = input_db->getDatabase("level_solver_db");
        }
    }

    // Setup scratch variables.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(d_object_name + "::CONTEXT");
    const IntVector<NDIM> side_ghosts = SIDEG;
    Pointer<SideVariable<NDIM, double> > side_scratch_var =
        new SideVariable<NDIM, double>(d_object_name + "::side_scratch");
    if (var_db->checkVariableExists(side_scratch_var->getName()))
    {
        side_scratch_var = var_db->getVariable(side_scratch_var->getName());
        d_side_scratch_idx = var_db->mapVariableAndContextToIndex(side_scratch_var, d_context);
        var_db->removePatchDataIndex(d_side_scratch_idx);
    }
    d_side_scratch_idx = var_db->registerVariableAndContext(side_scratch_var, d_context, side_ghosts);
    const IntVector<NDIM> cell_ghosts = CELLG;
    Pointer<CellVariable<NDIM, double> > cell_scratch_var =
        new CellVariable<NDIM, double>(d_object_name + "::cell_scratch");
    if (var_db->checkVariableExists(cell_scratch_var->getName()))
    {
        cell_scratch_var = var_db->getVariable(cell_scratch_var->getName());
        d_cell_scratch_idx = var_db->mapVariableAndContextToIndex(cell_scratch_var, d_context);
        var_db->removePatchDataIndex(d_cell_scratch_idx);
    }
    d_cell_scratch_idx = var_db->registerVariableAndContext(cell_scratch_var, d_context, cell_ghosts);

    // Construct the DOF index variable/context.
    d_u_dof_index_var = new SideVariable<NDIM, int>(object_name + "::u_dof_index");
    if (var_db->checkVariableExists(d_u_dof_index_var->getName()))
    {
        d_u_dof_index_var = var_db->getVariable(d_u_dof_index_var->getName());
        d_u_dof_index_idx = var_db->mapVariableAndContextToIndex(d_u_dof_index_var, d_context);
        var_db->removePatchDataIndex(d_u_dof_index_idx);
    }
    d_u_dof_index_idx = var_db->registerVariableAndContext(d_u_dof_index_var, d_context, NOGHOST);
    d_p_dof_index_var = new CellVariable<NDIM, int>(object_name + "::p_dof_index");
    if (var_db->checkVariableExists(d_p_dof_index_var->getName()))
    {
        d_p_dof_index_var = var_db->getVariable(d_p_dof_index_var->getName());
        d_p_dof_index_idx = var_db->mapVariableAndContextToIndex(d_p_dof_index_var, d_context);
        var_db->removePatchDataIndex(d_p_dof_index_idx);
    }
    d_p_dof_index_idx = var_db->registerVariableAndContext(d_p_dof_index_var, d_context, NOGHOST);
    // Setup Timers.
    IBAMR_DO_ONCE(
        t_compute_residual =
            TimerManager::getManager()->getTimer("StaggeredStokesIBLevelRelaxationFACOperator::computeResidual()");
        t_restrict_residual =
            TimerManager::getManager()->getTimer("StaggeredStokesIBLevelRelaxationFACOperator::restrictResidual()");
        t_prolong_error =
            TimerManager::getManager()->getTimer("StaggeredStokesIBLevelRelaxationFACOperator::prolongError()");
        t_prolong_error_and_correct = TimerManager::getManager()->getTimer(
            "StaggeredStokesIBLevelRelaxationFACOperator::prolongErrorAndCorrect()");
        t_smooth_error =
            TimerManager::getManager()->getTimer("IBAMR::StaggeredStokesLevelRelaxationFACOperator::smoothError()");
        t_initialize_operator_state = TimerManager::getManager()->getTimer(
            "StaggeredStokesIBLevelRelaxationFACOperator::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer(
            "StaggeredStokesIBLevelRelaxationFACOperator::deallocateOperatorState()"););
    return;
} // StaggeredStokesIBLevelRelaxationFACOperator

StaggeredStokesIBLevelRelaxationFACOperator::~StaggeredStokesIBLevelRelaxationFACOperator()
{
    if (d_is_initialized)
    {
        deallocateOperatorState();
    }
    delete d_default_U_bc_coef;
    d_default_U_bc_coef = NULL;
    delete d_default_P_bc_coef;
    d_default_P_bc_coef = NULL;
    return;
} // ~StaggeredStokesIBLevelRelaxationFACOperator

void StaggeredStokesIBLevelRelaxationFACOperator::setVelocityPoissonSpecifications(
    const PoissonSpecifications& U_problem_coefs)
{
    d_U_problem_coefs = U_problem_coefs;
    return;
} // setVelocityPoissonSpecifications

void StaggeredStokesIBLevelRelaxationFACOperator::setComponentsHaveNullspace(const bool has_velocity_nullspace,
                                                                             const bool has_pressure_nullspace)
{
    d_has_velocity_nullspace = has_velocity_nullspace;
    d_has_pressure_nullspace = has_pressure_nullspace;

    return;
} // setComponentsHaveNullspace

void StaggeredStokesIBLevelRelaxationFACOperator::setPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
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

void StaggeredStokesIBLevelRelaxationFACOperator::setPhysicalBoundaryHelper(
    Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_helper);
#endif
    d_bc_helper = bc_helper;
    return;
} // setPhysicalBoundaryHelper

void StaggeredStokesIBLevelRelaxationFACOperator::setResetLevels(const int coarsest_ln, const int finest_ln)
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

void StaggeredStokesIBLevelRelaxationFACOperator::setSmootherType(const std::string& smoother_type)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setSmootherType()\n"
                                 << "  cannot be called while operator state is initialized" << std::endl);
    }
    if (d_level_solver_type != smoother_type)
    {
        d_level_solvers.clear();
    }
    d_level_solver_type = smoother_type;
    d_smoother_type = smoother_type;

    return;
} // setSmootherType

void StaggeredStokesIBLevelRelaxationFACOperator::setCoarseSolverType(const std::string& coarse_solver_type)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setCoarseSolverType():\n"
                                 << "  cannot be called while operator state is initialized" << std::endl);
    }
    if (d_coarse_solver_type != coarse_solver_type) d_coarse_solver.setNull();
    d_coarse_solver_type = coarse_solver_type;
    if (!d_coarse_solver)
    {
        d_coarse_solver = StaggeredStokesSolverManager::getManager()->allocateSolver(
            d_coarse_solver_type, d_object_name + "::coarse_solver", d_coarse_solver_db,
            d_coarse_solver_default_options_prefix);
    }
    return;
} // setCoarseSolverType

void StaggeredStokesIBLevelRelaxationFACOperator::setCoarseSolverMaxIterations(int coarse_solver_max_iterations)
{
    d_coarse_solver_max_iterations = coarse_solver_max_iterations;
    return;
} // setCoarseSolverMaxIterations

void StaggeredStokesIBLevelRelaxationFACOperator::setCoarseSolverAbsoluteTolerance(
    double coarse_solver_abs_residual_tol)
{
    d_coarse_solver_abs_residual_tol = coarse_solver_abs_residual_tol;
    return;
} // setCoarseSolverAbsoluteTolerance

void StaggeredStokesIBLevelRelaxationFACOperator::setCoarseSolverRelativeTolerance(
    double coarse_solver_rel_residual_tol)
{
    d_coarse_solver_rel_residual_tol = coarse_solver_rel_residual_tol;
    return;
} // setCoarseSolverRelativeTolerance

void StaggeredStokesIBLevelRelaxationFACOperator::setProlongationMethods(const std::string& U_prolongation_method,
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
} // setProlongationMethods

void StaggeredStokesIBLevelRelaxationFACOperator::setRestrictionMethods(const std::string& U_restriction_method,
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
} // setRestrictionMethods

void StaggeredStokesIBLevelRelaxationFACOperator::setIBForceJacobian(Mat& A)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setIBForceJacobian()\n"
                                 << "  cannot be called while operator state is initialized" << std::endl);
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(A);
#endif

    d_A_mat = A;
    return;
} // setIBForceJacobian

void StaggeredStokesIBLevelRelaxationFACOperator::setIBInterpOp(Mat& J)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setIBInterpOp()\n"
                                 << "  cannot be called while operator state is initialized" << std::endl);
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(J);
#endif

    d_J_mat = J;
    return;
} // setIBInterpOp

Pointer<StaggeredStokesPETScLevelSolver>
StaggeredStokesIBLevelRelaxationFACOperator::getStaggeredStokesIBLevelSolver(const int ln) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(ln >= 0 && ln <= d_finest_ln);
#endif
    return d_level_solvers[ln];
} // getStaggeredStokesIBLevelSolver

Mat StaggeredStokesIBLevelRelaxationFACOperator::getGalerkinElasticityLevelOp(const int ln) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(ln >= 0 && ln <= d_finest_ln);
#endif
    return d_SAJ_mat[ln];
} // getGalerkinElasticityLevelOp

Mat StaggeredStokesIBLevelRelaxationFACOperator::getProlongationOp(const int ln) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(ln >= 0 && ln < d_finest_ln);
#endif
    return d_prolongation_mat[ln];
} // getProlongationOp

Vec StaggeredStokesIBLevelRelaxationFACOperator::getRestrictionScalingOp(const int ln) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(ln >= 0 && ln < d_finest_ln);
#endif
    return d_scale_restriction_mat[ln];
} // getRestrictionScalingOp

void StaggeredStokesIBLevelRelaxationFACOperator::restrictResidual(const SAMRAIVectorReal<NDIM, double>& src,
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

void StaggeredStokesIBLevelRelaxationFACOperator::prolongError(const SAMRAIVectorReal<NDIM, double>& src,
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

void StaggeredStokesIBLevelRelaxationFACOperator::prolongErrorAndCorrect(const SAMRAIVectorReal<NDIM, double>& src,
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

bool StaggeredStokesIBLevelRelaxationFACOperator::solveCoarsestLevel(SAMRAIVectorReal<NDIM, double>& error,
                                                                     const SAMRAIVectorReal<NDIM, double>& residual,
                                                                     int coarsest_ln)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(coarsest_ln == d_coarsest_ln);
#endif

    d_coarse_solver->setSolutionTime(d_solution_time);
    d_coarse_solver->setTimeInterval(d_current_time, d_new_time);
    d_coarse_solver->setMaxIterations(d_coarse_solver_max_iterations);
    d_coarse_solver->setAbsoluteTolerance(d_coarse_solver_abs_residual_tol);
    d_coarse_solver->setRelativeTolerance(d_coarse_solver_rel_residual_tol);
    d_coarse_solver->setComponentsHaveNullspace(d_has_velocity_nullspace, d_has_pressure_nullspace);

    bool initial_guess_nonzero = true;
    const KSP& petsc_ksp = d_coarse_solver->getPETScKSP();
    KSPType ksp_type;
    KSPGetType(petsc_ksp, &ksp_type);
    if (!strcmp(ksp_type, "preonly")) initial_guess_nonzero = false;
    d_coarse_solver->setInitialGuessNonzero(initial_guess_nonzero);

    d_coarse_solver->solveSystem(*getLevelSAMRAIVectorReal(error, d_coarsest_ln),
                                 *getLevelSAMRAIVectorReal(residual, d_coarsest_ln));

    return true;
} // solveCoarsestLevel

void StaggeredStokesIBLevelRelaxationFACOperator::computeResidual(SAMRAIVectorReal<NDIM, double>& residual,
                                                                  const SAMRAIVectorReal<NDIM, double>& solution,
                                                                  const SAMRAIVectorReal<NDIM, double>& rhs,
                                                                  int coarsest_level_num,
                                                                  int finest_level_num)
{
    IBAMR_TIMER_START(t_compute_residual);

    const int U_res_idx = residual.getComponentDescriptorIndex(0);
    const int U_sol_idx = solution.getComponentDescriptorIndex(0);
    const int U_rhs_idx = rhs.getComponentDescriptorIndex(0);

    const int P_res_idx = residual.getComponentDescriptorIndex(1);
    const int P_sol_idx = solution.getComponentDescriptorIndex(1);
    const int P_rhs_idx = rhs.getComponentDescriptorIndex(1);

    // Compute the residual, r = f - A*u.
    int ierr;
    int rank = SAMRAI_MPI::getRank();
    for (int ln = coarsest_level_num; ln <= finest_level_num; ++ln)
    {
        Vec solution_vec;
        Vec residual_vec;
        Vec rhs_vec;
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[ln][rank], PETSC_DETERMINE, &solution_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[ln][rank], PETSC_DETERMINE, &residual_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[ln][rank], PETSC_DETERMINE, &rhs_vec);
        IBTK_CHKERRQ(ierr);

        StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(solution_vec, U_sol_idx, d_u_dof_index_idx, P_sol_idx,
                                                              d_p_dof_index_idx, level);
        StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(rhs_vec, U_rhs_idx, d_u_dof_index_idx, P_rhs_idx,
                                                              d_p_dof_index_idx, level);

        const KSP& level_ksp =
            (ln == d_coarsest_ln) ? d_coarse_solver->getPETScKSP() : d_level_solvers[ln]->getPETScKSP();
        Mat A;
        ierr = KSPGetOperators(level_ksp, &A, PETSC_NULL);
        IBTK_CHKERRQ(ierr);
        ierr = MatMult(A, solution_vec, residual_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecAYPX(residual_vec, -1.0, rhs_vec);
        IBTK_CHKERRQ(ierr);

        StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(residual_vec, U_res_idx, d_u_dof_index_idx, P_res_idx,
                                                                d_p_dof_index_idx, level, d_synch_refine_schedules[ln],
                                                                d_ghostfill_nocoarse_refine_schedules[ln]);

        ierr = VecDestroy(&solution_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&rhs_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&residual_vec);
        IBTK_CHKERRQ(ierr);
    }

    IBAMR_TIMER_STOP(t_compute_residual);

    return;
} // computeResidual

void StaggeredStokesIBLevelRelaxationFACOperator::smoothError(SAMRAIVectorReal<NDIM, double>& error,
                                                              const SAMRAIVectorReal<NDIM, double>& residual,
                                                              int level_num,
                                                              int num_sweeps,
                                                              bool /*performing_pre_sweeps*/,
                                                              bool /*performing_post_sweeps*/)
{
    if (num_sweeps == 0) return;

    IBAMR_TIMER_START(t_smooth_error);

    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    const int U_error_idx = error.getComponentDescriptorIndex(0);
    const int P_error_idx = error.getComponentDescriptorIndex(1);
    const int U_scratch_idx = d_side_scratch_idx;
    const int P_scratch_idx = d_cell_scratch_idx;

    Pointer<SAMRAIVectorReal<NDIM, double> > e_level = getLevelSAMRAIVectorReal(error, level_num);
    Pointer<SAMRAIVectorReal<NDIM, double> > r_level = getLevelSAMRAIVectorReal(residual, level_num);

    // Cache coarse-fine interface ghost cell values in the "scratch" data.
    if (level_num > d_coarsest_ln && num_sweeps > 1)
    {
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM, double> > U_error_data = error.getComponentPatchData(0, *patch);
            Pointer<SideData<NDIM, double> > U_scratch_data = patch->getPatchData(U_scratch_idx);
#if !defined(NDEBUG)
            const Box<NDIM>& U_ghost_box = U_error_data->getGhostBox();
            TBOX_ASSERT(U_ghost_box == U_scratch_data->getGhostBox());
            TBOX_ASSERT(U_error_data->getGhostCellWidth() == SIDEG);
            TBOX_ASSERT(U_scratch_data->getGhostCellWidth() == SIDEG);
#endif
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                U_scratch_data->getArrayData(axis).copy(U_error_data->getArrayData(axis),
                                                        d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                                                        IntVector<NDIM>(0));
            }

            Pointer<CellData<NDIM, double> > P_error_data = error.getComponentPatchData(1, *patch);
            Pointer<CellData<NDIM, double> > P_scratch_data = patch->getPatchData(P_scratch_idx);
#if !defined(NDEBUG)
            const Box<NDIM>& P_ghost_box = P_error_data->getGhostBox();
            TBOX_ASSERT(P_ghost_box == P_scratch_data->getGhostBox());
            TBOX_ASSERT(P_error_data->getGhostCellWidth() == CELLG);
            TBOX_ASSERT(P_scratch_data->getGhostCellWidth() == CELLG);
#endif
            P_scratch_data->getArrayData().copy(P_error_data->getArrayData(),
                                                d_patch_cell_bc_box_overlap[level_num][patch_counter],
                                                IntVector<NDIM>(0));
        }
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
                int patch_counter = 0;
                for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());

                    Pointer<SideData<NDIM, double> > U_error_data = error.getComponentPatchData(0, *patch);
                    Pointer<SideData<NDIM, double> > U_scratch_data = patch->getPatchData(U_scratch_idx);
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        U_error_data->getArrayData(axis)
                            .copy(U_scratch_data->getArrayData(axis),
                                  d_patch_side_bc_box_overlap[level_num][patch_counter][axis], IntVector<NDIM>(0));
                    }

                    Pointer<CellData<NDIM, double> > P_error_data = error.getComponentPatchData(1, *patch);
                    Pointer<CellData<NDIM, double> > P_scratch_data = patch->getPatchData(P_scratch_idx);
                    P_error_data->getArrayData().copy(P_scratch_data->getArrayData(),
                                                      d_patch_cell_bc_box_overlap[level_num][patch_counter],
                                                      IntVector<NDIM>(0));
                }
            }

            // Complete the coarse-fine interface interpolation by computing the
            // normal extension. This uses latest interior values on fine patch
            // but cached tangential values at the fine locations of the
            // coarse-fine interface
            d_U_cf_bdry_op->setPatchDataIndex(U_error_idx);
            d_P_cf_bdry_op->setPatchDataIndex(P_error_idx);
            const IntVector<NDIM>& ratio = level->getRatioToCoarserLevel();
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                d_U_cf_bdry_op->computeNormalExtension(*patch, ratio, SIDEG);
                d_P_cf_bdry_op->computeNormalExtension(*patch, ratio, CELLG);
            }
        }

        // Smooth the error on the level.
        Pointer<StaggeredStokesPETScLevelSolver> level_solver = d_level_solvers[level_num];
        level_solver->setSolutionTime(d_solution_time);
        level_solver->setTimeInterval(d_current_time, d_new_time);
        level_solver->setMaxIterations(d_level_solver_max_iterations);
        level_solver->setAbsoluteTolerance(d_level_solver_abs_residual_tol);
        level_solver->setRelativeTolerance(d_level_solver_rel_residual_tol);
        level_solver->setComponentsHaveNullspace(d_has_velocity_nullspace, d_has_pressure_nullspace);

        bool initial_guess_nonzero = true;
        const KSP& petsc_ksp = level_solver->getPETScKSP();
        KSPType ksp_type;
        KSPGetType(petsc_ksp, &ksp_type);
        if (!strcmp(ksp_type, "preonly")) initial_guess_nonzero = false;
        level_solver->setInitialGuessNonzero(initial_guess_nonzero);

        level_solver->solveSystem(*e_level, *r_level);
    }

    IBAMR_TIMER_STOP(t_smooth_error);

    return;
} // smoothError

void StaggeredStokesIBLevelRelaxationFACOperator::initializeOperatorState(
    const SAMRAIVectorReal<NDIM, double>& solution,
    const SAMRAIVectorReal<NDIM, double>& rhs)
{
    IBAMR_TIMER_START(t_initialize_operator_state);

    d_in_initialize_operator_state = true;
    int ierr;

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
    d_U_op_stencil_fill_pattern = NULL;
    d_P_op_stencil_fill_pattern = new CellNoCornersFillPattern(CELLG, false, false, false);
    d_U_synch_fill_pattern = new SideSynchCopyFillPattern();

    // Construct patch level DOFs.
    d_num_dofs_per_proc.resize(d_finest_ln + 1);
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln - 1); ln <= std::min(d_finest_ln, finest_reset_ln); ++ln)
    {
        // Allocate DOF index data.
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_u_dof_index_idx)) level->allocatePatchData(d_u_dof_index_idx);
        if (!level->checkAllocated(d_p_dof_index_idx)) level->allocatePatchData(d_p_dof_index_idx);

        // Construct DOF indices and SAMRAI to PETSc ordering.
        StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(d_num_dofs_per_proc[ln], d_u_dof_index_idx,
                                                                        d_p_dof_index_idx, level);
    }

    // Setup application ordering for the velocity DOFs.
    d_u_app_ordering.resize(d_finest_ln + 1);
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln - 1); ln <= std::min(d_finest_ln - 1, finest_reset_ln);
         ++ln)
    {

        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        PETScVecUtilities::constructPatchLevelAO(d_u_app_ordering[ln], d_num_dofs_per_proc[ln], d_u_dof_index_idx,
                                                 level, /*ao offset*/ 0);
    }

    // Construct prolongation matrix and scaling for restriction matrix for various levels.
    d_prolongation_mat.resize(d_finest_ln + 1);
    d_scale_restriction_mat.resize(d_finest_ln + 1);
    for (int ln = std::min(d_finest_ln - 1, finest_reset_ln); ln >= std::max(d_coarsest_ln, coarsest_reset_ln - 1);
         --ln)
    {

        Pointer<PatchLevel<NDIM> > fine_level = d_hierarchy->getPatchLevel(ln + 1);
        Pointer<PatchLevel<NDIM> > coarse_level = d_hierarchy->getPatchLevel(ln);

        PETScMatUtilities::constructProlongationOp(d_prolongation_mat[ln], d_u_dof_index_idx,
                                                   d_num_dofs_per_proc[ln + 1], d_num_dofs_per_proc[ln], fine_level,
                                                   coarse_level, d_u_app_ordering[ln]);

        // Get the scaling for restriction operator R = P^T.
        PETScMatUtilities::constructRestrictionScalingOp(d_prolongation_mat[ln], d_scale_restriction_mat[ln]);
    }

    // Compute SAJ operator for various patch levels.
    d_SAJ_mat.resize(d_finest_ln + 1);
    for (int ln = std::min(d_finest_ln, finest_reset_ln); ln >= std::max(d_coarsest_ln, coarsest_reset_ln - 1); --ln)
    {
        if (ln == d_finest_ln)
        {
            ierr = MatPtAP(d_A_mat, d_J_mat, MAT_INITIAL_MATRIX, 1.0, &d_SAJ_mat[ln]);
            IBTK_CHKERRQ(ierr);

            // Compute the scale for the spreading operator.
            Pointer<PatchLevel<NDIM> > finest_level = d_hierarchy->getPatchLevel(d_finest_ln);
            Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
            const double* const dx0 = grid_geom->getDx();
            IntVector<NDIM> ratio = finest_level->getRatio();
            double spread_scale = -0.25 * (d_new_time - d_current_time);
            for (unsigned d = 0; d < NDIM; ++d) spread_scale *= ratio(d) / dx0[d];

            ierr = MatScale(d_SAJ_mat[ln], spread_scale);
            IBTK_CHKERRQ(ierr);
        }
        else
        {
            ierr = MatPtAP(d_SAJ_mat[ln + 1], d_prolongation_mat[ln], MAT_INITIAL_MATRIX, 1.0, &d_SAJ_mat[ln]);
            IBTK_CHKERRQ(ierr);
            ierr = MatDiagonalScale(d_SAJ_mat[ln], d_scale_restriction_mat[ln], PETSC_NULL);
            IBTK_CHKERRQ(ierr);
        }
    }

    // Initialize the coarse level solver when needed.
    if (coarsest_reset_ln == d_coarsest_ln)
    {
        if (!d_coarse_solver)
        {
            setCoarseSolverType(d_coarse_solver_type);
        }

        // Note that since the coarse level solver is solving for the error, it
        // must always employ homogeneous boundary conditions.
        d_coarse_solver->setSolutionTime(d_solution_time);
        d_coarse_solver->setTimeInterval(d_current_time, d_new_time);
        d_coarse_solver->setVelocityPoissonSpecifications(d_U_problem_coefs);
        d_coarse_solver->setPhysicalBcCoefs(d_U_bc_coefs, d_P_bc_coef);
        d_coarse_solver->setPhysicalBoundaryHelper(d_bc_helper);
        d_coarse_solver->setHomogeneousBc(true);
        d_coarse_solver->setComponentsHaveNullspace(d_has_velocity_nullspace, d_has_pressure_nullspace);
        d_coarse_solver->addLinearOperator(d_SAJ_mat[d_coarsest_ln]);
        d_coarse_solver->initializeSolverState(*getLevelSAMRAIVectorReal(*d_solution, d_coarsest_ln),
                                               *getLevelSAMRAIVectorReal(*d_rhs, d_coarsest_ln));
    }

    // Initialize the fine level solvers when needed.
    d_level_solvers.resize(d_finest_ln + 1);
    for (int ln = std::max(1, coarsest_reset_ln); ln <= finest_reset_ln; ++ln)
    {
        Pointer<StaggeredStokesPETScLevelSolver>& level_solver = d_level_solvers[ln];
        if (!level_solver)
        {
            std::ostringstream level_solver_stream;
            level_solver_stream << d_level_solver_default_options_prefix;
            level_solver = StaggeredStokesSolverManager::getManager()->allocateSolver(
                d_level_solver_type, d_object_name + "::level_solver", d_level_solver_db, level_solver_stream.str());
        }

        level_solver->setSolutionTime(d_solution_time);
        level_solver->setTimeInterval(d_current_time, d_new_time);
        level_solver->setVelocityPoissonSpecifications(d_U_problem_coefs);
        level_solver->setPhysicalBcCoefs(d_U_bc_coefs, d_P_bc_coef);
        level_solver->setPhysicalBoundaryHelper(d_bc_helper);
        level_solver->setHomogeneousBc(true);
        level_solver->setComponentsHaveNullspace(d_has_velocity_nullspace, d_has_pressure_nullspace);
        level_solver->addLinearOperator(d_SAJ_mat[ln]);
        level_solver->initializeSolverState(*getLevelSAMRAIVectorReal(*d_solution, ln),
                                            *getLevelSAMRAIVectorReal(*d_rhs, ln));
    }

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

    d_prolongation_refine_algorithm->registerRefine(d_side_scratch_idx, solution.getComponentDescriptorIndex(0),
                                                    d_side_scratch_idx, d_U_prolongation_refine_operator,
                                                    d_U_op_stencil_fill_pattern);
    d_prolongation_refine_algorithm->registerRefine(d_cell_scratch_idx, solution.getComponentDescriptorIndex(1),
                                                    d_cell_scratch_idx, d_P_prolongation_refine_operator,
                                                    d_P_op_stencil_fill_pattern);

    d_restriction_coarsen_algorithm->registerCoarsen(d_side_scratch_idx, rhs.getComponentDescriptorIndex(0),
                                                     d_U_restriction_coarsen_operator);
    d_restriction_coarsen_algorithm->registerCoarsen(d_cell_scratch_idx, rhs.getComponentDescriptorIndex(1),
                                                     d_P_restriction_coarsen_operator);

    d_ghostfill_nocoarse_refine_algorithm->registerRefine(
        solution.getComponentDescriptorIndex(0), solution.getComponentDescriptorIndex(0),
        solution.getComponentDescriptorIndex(0), Pointer<RefineOperator<NDIM> >(), d_U_op_stencil_fill_pattern);
    d_ghostfill_nocoarse_refine_algorithm->registerRefine(
        solution.getComponentDescriptorIndex(1), solution.getComponentDescriptorIndex(1),
        solution.getComponentDescriptorIndex(1), Pointer<RefineOperator<NDIM> >(), d_P_op_stencil_fill_pattern);

    d_synch_refine_algorithm->registerRefine(
        solution.getComponentDescriptorIndex(0), solution.getComponentDescriptorIndex(0),
        solution.getComponentDescriptorIndex(0), Pointer<RefineOperator<NDIM> >(), d_U_synch_fill_pattern);

    std::vector<RefinePatchStrategy<NDIM>*> bc_op_ptrs(2);
    bc_op_ptrs[0] = d_U_bc_op;
    bc_op_ptrs[1] = d_P_bc_op;
    d_U_P_bc_op = new RefinePatchStrategySet(bc_op_ptrs.begin(), bc_op_ptrs.end(), false);

    for (int dst_ln = d_coarsest_ln + 1; dst_ln <= d_finest_ln; ++dst_ln)
    {
        d_prolongation_refine_schedules[dst_ln] = d_prolongation_refine_algorithm->createSchedule(
            d_hierarchy->getPatchLevel(dst_ln), Pointer<PatchLevel<NDIM> >(), dst_ln - 1, d_hierarchy,
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

    // Get overlap information for setting patch boundary conditions.
    d_patch_side_bc_box_overlap.resize(d_finest_ln + 1);
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
        d_patch_side_bc_box_overlap[ln].resize(num_local_patches);
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                const Box<NDIM> side_ghost_box = Box<NDIM>::grow(side_box, 1);
                d_patch_side_bc_box_overlap[ln][patch_counter][axis] = BoxList<NDIM>(side_ghost_box);
                d_patch_side_bc_box_overlap[ln][patch_counter][axis].removeIntersections(side_box);
            }
        }
    }

    d_patch_cell_bc_box_overlap.resize(d_finest_ln + 1);
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
        d_patch_cell_bc_box_overlap[ln].resize(num_local_patches);

        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Box<NDIM>& ghost_box = Box<NDIM>::grow(patch_box, 1);

            d_patch_cell_bc_box_overlap[ln][patch_counter] = BoxList<NDIM>(ghost_box);
            d_patch_cell_bc_box_overlap[ln][patch_counter].removeIntersections(patch_box);
        }
    }

    // Indicate that the operator is initialized.
    d_is_initialized = true;
    d_in_initialize_operator_state = false;

    IBAMR_TIMER_STOP(t_initialize_operator_state);
    return;
} // initializeOperatorState

void StaggeredStokesIBLevelRelaxationFACOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_operator_state);

    int ierr;
    const int coarsest_reset_ln =
        (d_in_initialize_operator_state && (d_coarsest_reset_ln != -1) && (d_finest_reset_ln != -1)) ?
            d_coarsest_reset_ln :
            d_coarsest_ln;
    const int finest_reset_ln =
        (d_in_initialize_operator_state && (d_coarsest_reset_ln != -1) && (d_finest_reset_ln != -1)) ?
            d_finest_reset_ln :
            d_finest_ln;

    // Deallocate level solvers and overlap boxes for side and cell data.
    for (int ln = coarsest_reset_ln; ln <= std::min(d_finest_ln, finest_reset_ln); ++ln)
    {
        if (d_level_solvers[ln]) d_level_solvers[ln]->deallocateSolverState();
        d_patch_side_bc_box_overlap[ln].resize(0);
        d_patch_cell_bc_box_overlap[ln].resize(0);
    }

    // Deallocate application ordering
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln - 1); ln <= std::min(d_finest_ln - 1, finest_reset_ln);
         ++ln)
    {
        ierr = AODestroy(&d_u_app_ordering[ln]);
        IBTK_CHKERRQ(ierr);
        d_u_app_ordering[ln] = PETSC_NULL;
    }

    // Deallocate prolongation mat and restriction scaling mat
    for (int ln = std::min(d_finest_ln - 1, finest_reset_ln); ln >= std::max(d_coarsest_ln, coarsest_reset_ln - 1);
         --ln)
    {
        ierr = MatDestroy(&d_prolongation_mat[ln]);
        IBTK_CHKERRQ(ierr);
        d_prolongation_mat[ln] = PETSC_NULL;

        ierr = VecDestroy(&d_scale_restriction_mat[ln]);
        IBTK_CHKERRQ(ierr);
        d_scale_restriction_mat[ln] = PETSC_NULL;
    }

    // Deallocate SAJ mat.
    for (int ln = std::min(d_finest_ln, finest_reset_ln); ln >= std::max(d_coarsest_ln, coarsest_reset_ln - 1); --ln)
    {
        ierr = MatDestroy(&d_SAJ_mat[ln]);
        IBTK_CHKERRQ(ierr);
        d_SAJ_mat[ln] = PETSC_NULL;
    }

    // Deallocate DOF index data.
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln - 1); ln <= std::min(d_finest_ln, finest_reset_ln); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_u_dof_index_idx)) level->deallocatePatchData(d_u_dof_index_idx);
        if (level->checkAllocated(d_p_dof_index_idx)) level->deallocatePatchData(d_p_dof_index_idx);
    }

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

void StaggeredStokesIBLevelRelaxationFACOperator::allocateScratchData()
{
    if (d_solution) d_solution->allocateVectorData();
    if (d_rhs) d_rhs->allocateVectorData();
    return;
}

void StaggeredStokesIBLevelRelaxationFACOperator::deallocateScratchData()
{
    if (d_solution) d_solution->deallocateVectorData();
    if (d_rhs) d_rhs->deallocateVectorData();
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void StaggeredStokesIBLevelRelaxationFACOperator::xeqScheduleProlongation(const std::pair<int, int>& dst_idxs,
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
    refine_alg.registerRefine(U_dst_idx, U_src_idx, U_dst_idx, d_U_prolongation_refine_operator,
                              d_U_op_stencil_fill_pattern);
    refine_alg.registerRefine(P_dst_idx, P_src_idx, P_dst_idx, d_P_prolongation_refine_operator,
                              d_P_op_stencil_fill_pattern);
    refine_alg.resetSchedule(d_prolongation_refine_schedules[dst_ln]);
    d_prolongation_refine_schedules[dst_ln]->fillData(d_new_time);
    d_prolongation_refine_algorithm->resetSchedule(d_prolongation_refine_schedules[dst_ln]);
    return;
} // xeqScheduleProlongation

void StaggeredStokesIBLevelRelaxationFACOperator::xeqScheduleRestriction(const std::pair<int, int>& dst_idxs,
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

void StaggeredStokesIBLevelRelaxationFACOperator::xeqScheduleGhostFillNoCoarse(const std::pair<int, int>& dst_idxs,
                                                                               const int dst_ln)
{
    const int U_dst_idx = dst_idxs.first;
    d_U_bc_op->setPatchDataIndex(U_dst_idx);
    d_U_bc_op->setHomogeneousBc(true);

    const int P_dst_idx = dst_idxs.second;
    d_P_bc_op->setPatchDataIndex(P_dst_idx);
    d_P_bc_op->setHomogeneousBc(true);

    RefineAlgorithm<NDIM> refine_alg;
    refine_alg.registerRefine(U_dst_idx, U_dst_idx, U_dst_idx, Pointer<RefineOperator<NDIM> >(),
                              d_U_op_stencil_fill_pattern);
    refine_alg.registerRefine(P_dst_idx, P_dst_idx, P_dst_idx, Pointer<RefineOperator<NDIM> >(),
                              d_P_op_stencil_fill_pattern);
    refine_alg.resetSchedule(d_ghostfill_nocoarse_refine_schedules[dst_ln]);
    d_ghostfill_nocoarse_refine_schedules[dst_ln]->fillData(d_new_time);
    d_ghostfill_nocoarse_refine_algorithm->resetSchedule(d_ghostfill_nocoarse_refine_schedules[dst_ln]);
    return;
} // xeqScheduleGhostFillNoCoarse

void StaggeredStokesIBLevelRelaxationFACOperator::xeqScheduleDataSynch(const int U_dst_idx, const int dst_ln)
{
    RefineAlgorithm<NDIM> refine_alg;
    refine_alg.registerRefine(U_dst_idx, U_dst_idx, U_dst_idx, Pointer<RefineOperator<NDIM> >(),
                              d_U_synch_fill_pattern);
    refine_alg.resetSchedule(d_synch_refine_schedules[dst_ln]);
    d_synch_refine_schedules[dst_ln]->fillData(d_new_time);
    d_synch_refine_algorithm->resetSchedule(d_synch_refine_schedules[dst_ln]);
    return;
} // xeqScheduleDataSynch

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
