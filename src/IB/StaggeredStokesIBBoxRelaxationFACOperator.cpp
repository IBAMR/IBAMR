// Filename: StaggeredStokesIBBoxRelaxationFACOperator.cpp
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
#include <limits>
#include <numeric>
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
#include "ibamr/StaggeredStokesIBBoxRelaxationFACOperator.h"
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
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesIBBoxRelaxationFACOperator::StaggeredStokesIBBoxRelaxationFACOperator(
    const std::string& object_name,
    const Pointer<Database> input_db,
    const std::string& default_options_prefix)
    : StaggeredStokesFACPreconditionerStrategy(object_name, std::max(SIDEG, CELLG), input_db, default_options_prefix),
      d_level_solver_type("PETSC_LEVEL_SOLVER"),
      d_level_solver_default_options_prefix(default_options_prefix + "level_"),
      d_level_solver_abs_residual_tol(1.0e-50),
      d_level_solver_rel_residual_tol(1.0e-5),
      d_level_solver_max_iterations(10)
{
    // Indicate that this subclass handles initializaing the coarse-grid solver.
    d_coarse_solver_init_subclass = true;

    // Get values from the input database.
    if (input_db)
    {
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

    // Construct the DOF index variable/context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
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
    return;
} // StaggeredStokesIBBoxRelaxationFACOperator

StaggeredStokesIBBoxRelaxationFACOperator::~StaggeredStokesIBBoxRelaxationFACOperator()
{
    if (d_is_initialized)
    {
        deallocateOperatorState();
    }
    return;
} // ~StaggeredStokesIBBoxRelaxationFACOperator

void
StaggeredStokesIBBoxRelaxationFACOperator::setIBForceJacobian(Mat& A)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setIBForceJacobian()\n"
                                 << "  cannot be called while operator state is initialized"
                                 << std::endl);
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(A);
#endif
    d_A_mat = A;
    return;
} // setIBForceJacobian

void
StaggeredStokesIBBoxRelaxationFACOperator::setIBInterpOp(Mat& J)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setIBInterpOp()\n"
                                 << "  cannot be called while operator state is initialized"
                                 << std::endl);
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(J);
#endif
    d_J_mat = J;
    return;
} // setIBInterpOp

Pointer<StaggeredStokesPETScLevelSolver>
StaggeredStokesIBBoxRelaxationFACOperator::getStaggeredStokesPETScLevelSolver(const int ln) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(ln >= 0 && ln <= d_finest_ln);
#endif
    return d_level_solvers[ln];
} // getStaggeredStokesLevelSolver

Mat
StaggeredStokesIBBoxRelaxationFACOperator::getEulerianElasticityLevelOp(const int ln) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(ln >= 0 && ln <= d_finest_ln);
#endif
    return d_SAJ_mat[ln];
} // getEulerianElasticityLevelOp

Mat
StaggeredStokesIBBoxRelaxationFACOperator::getProlongationOp(const int ln) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(ln >= 0 && ln < d_finest_ln);
#endif
    return d_prolongation_mat[ln];
} // getProlongationOp

Vec
StaggeredStokesIBBoxRelaxationFACOperator::getRestrictionScalingOp(const int ln) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(ln >= 0 && ln < d_finest_ln);
#endif
    return d_scale_restriction_mat[ln];
} // getRestrictionScalingOp

void
StaggeredStokesIBBoxRelaxationFACOperator::computeResidual(SAMRAIVectorReal<NDIM, double>& residual,
                                                           const SAMRAIVectorReal<NDIM, double>& solution,
                                                           const SAMRAIVectorReal<NDIM, double>& rhs,
                                                           int coarsest_level_num,
                                                           int finest_level_num)
{
    StaggeredStokesFACPreconditionerStrategy::computeResidual(
        residual, solution, rhs, coarsest_level_num, finest_level_num);

    const int U_res_idx = residual.getComponentDescriptorIndex(0);
    const int U_sol_idx = solution.getComponentDescriptorIndex(0);

    const int P_res_idx = residual.getComponentDescriptorIndex(1);
    const int P_sol_idx = solution.getComponentDescriptorIndex(1);

    // Update the residual, r = f - A*u, to include the IB part of the operator.
    int ierr;
    int rank = SAMRAI_MPI::getRank();
    for (int ln = coarsest_level_num; ln <= finest_level_num; ++ln)
    {
        Vec solution_vec, residual_vec;
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[ln][rank], PETSC_DETERMINE, &solution_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[ln][rank], PETSC_DETERMINE, &residual_vec);
        IBTK_CHKERRQ(ierr);

        StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            solution_vec, U_sol_idx, d_u_dof_index_idx, P_sol_idx, d_p_dof_index_idx, level);
        StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            residual_vec, U_res_idx, d_u_dof_index_idx, P_res_idx, d_p_dof_index_idx, level);

        ierr = VecScale(residual_vec, -1.0);
        IBTK_CHKERRQ(ierr);
        ierr = MatMultAdd(d_SAJ_mat[ln], solution_vec, residual_vec, residual_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecScale(residual_vec, -1.0);
        IBTK_CHKERRQ(ierr);

        StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(
            residual_vec, U_res_idx, d_u_dof_index_idx, P_res_idx, d_p_dof_index_idx, level, NULL, NULL);
        xeqScheduleDataSynch(U_res_idx, ln);
        xeqScheduleGhostFillNoCoarse(std::make_pair(U_res_idx, P_res_idx), ln);

        ierr = VecDestroy(&solution_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&residual_vec);
        IBTK_CHKERRQ(ierr);
    }
    return;
} // computeResidual

void
StaggeredStokesIBBoxRelaxationFACOperator::smoothError(SAMRAIVectorReal<NDIM, double>& error,
                                                       const SAMRAIVectorReal<NDIM, double>& residual,
                                                       int level_num,
                                                       int num_sweeps,
                                                       bool performing_pre_sweeps,
                                                       bool performing_post_sweeps)
{
    if (num_sweeps == 0) return;

    if (d_smoother_type == "MULTIPLICATIVE")
    {
        smoothErrorMultiplicative(
            error, residual, level_num, num_sweeps, performing_pre_sweeps, performing_post_sweeps);
    }
    else if (d_smoother_type == "REDBLACK")
    {
        smoothErrorRedBlackMultiplicative(
            error, residual, level_num, num_sweeps, performing_pre_sweeps, performing_post_sweeps);
    }
    else if (d_smoother_type == "ADDITIVE")
    {
        smoothErrorAdditive(error, residual, level_num, num_sweeps, performing_pre_sweeps, performing_post_sweeps);
    }
    else
    {
        TBOX_ERROR(d_object_name << "::smoothError()\n"
                                 << "Unknown smooth error option specified. Available options are ADDITIVE, "
                                    "MULTIPLICATIVE, REDBLACK"
                                 << std::endl);
    }

    return;
} // smoothError

/////////////////////////////// PROTECTED ////////////////////////////////////

void
StaggeredStokesIBBoxRelaxationFACOperator::initializeOperatorStateSpecialized(
    const SAMRAIVectorReal<NDIM, double>& /*solution*/,
    const SAMRAIVectorReal<NDIM, double>& /*rhs*/,
    const int coarsest_reset_ln,
    const int finest_reset_ln)
{
    int ierr;

    // Construct patch level DOFs.
    d_num_dofs_per_proc.resize(d_finest_ln + 1);
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln - 1); ln <= std::min(d_finest_ln, finest_reset_ln); ++ln)
    {
        // Allocate DOF index data.
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_u_dof_index_idx)) level->allocatePatchData(d_u_dof_index_idx);
        if (!level->checkAllocated(d_p_dof_index_idx)) level->allocatePatchData(d_p_dof_index_idx);

        // Construct DOF indices and SAMRAI to PETSc ordering.
        StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
            d_num_dofs_per_proc[ln], d_u_dof_index_idx, d_p_dof_index_idx, level);
    }

    // Setup application ordering for the velocity DOFs.
    d_u_app_ordering.resize(d_finest_ln + 1);
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln - 1); ln <= std::min(d_finest_ln - 1, finest_reset_ln);
         ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        PETScVecUtilities::constructPatchLevelAO(
            d_u_app_ordering[ln], d_num_dofs_per_proc[ln], d_u_dof_index_idx, level, /*ao offset*/ 0);
    }

    // Construct prolongation matrix and scaling for restriction matrix for various levels.
    d_prolongation_mat.resize(d_finest_ln + 1);
    d_scale_restriction_mat.resize(d_finest_ln + 1);
    for (int ln = std::min(d_finest_ln - 1, finest_reset_ln); ln >= std::max(d_coarsest_ln, coarsest_reset_ln - 1);
         --ln)
    {
        Pointer<PatchLevel<NDIM> > fine_level = d_hierarchy->getPatchLevel(ln + 1);
        Pointer<PatchLevel<NDIM> > coarse_level = d_hierarchy->getPatchLevel(ln);
        PETScMatUtilities::constructProlongationOp(d_prolongation_mat[ln],
                                                   d_u_dof_index_idx,
                                                   d_num_dofs_per_proc[ln + 1],
                                                   d_num_dofs_per_proc[ln],
                                                   fine_level,
                                                   coarse_level,
                                                   d_u_app_ordering[ln]);
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
            ierr = MatDiagonalScale(d_SAJ_mat[ln], d_scale_restriction_mat[ln], NULL);
            IBTK_CHKERRQ(ierr);
        }
    }

    // Initialize the level solvers when needed.
    if (d_coarse_solver_init_subclass && coarsest_reset_ln == d_coarsest_ln && d_coarse_solver_type != "LEVEL_SMOOTHER")
    {
        if (!d_coarse_solver)
        {
            d_coarse_solver =
                StaggeredStokesSolverManager::getManager()->allocateSolver(d_coarse_solver_type,
                                                                           d_object_name + "::coarse_solver",
                                                                           d_coarse_solver_db,
                                                                           d_coarse_solver_default_options_prefix);
        }
        d_coarse_solver->setSolutionTime(d_solution_time);
        d_coarse_solver->setTimeInterval(d_current_time, d_new_time);
        d_coarse_solver->setVelocityPoissonSpecifications(d_U_problem_coefs);
        d_coarse_solver->setPhysicalBcCoefs(d_U_bc_coefs, d_P_bc_coef);
        d_coarse_solver->setPhysicalBoundaryHelper(d_bc_helper);
        d_coarse_solver->setHomogeneousBc(true);
        d_coarse_solver->setComponentsHaveNullspace(d_has_velocity_nullspace, d_has_pressure_nullspace);
        Pointer<StaggeredStokesPETScLevelSolver> p_coarse_solver = d_coarse_solver;
        if (p_coarse_solver)
            p_coarse_solver->addLinearOperator(d_SAJ_mat[d_coarsest_ln]);
        else
            TBOX_ERROR("no mechanism for specifying SAJ part of operator!");
        d_coarse_solver->initializeSolverState(*getLevelSAMRAIVectorReal(*d_solution, d_coarsest_ln),
                                               *getLevelSAMRAIVectorReal(*d_rhs, d_coarsest_ln));
    }

    d_level_solvers.resize(d_finest_ln + 1);
    for (int ln = std::max(d_coarsest_ln + 1, coarsest_reset_ln); ln <= finest_reset_ln; ++ln)
    {
        Pointer<StaggeredStokesPETScLevelSolver>& level_solver = d_level_solvers[ln];
        if (!level_solver)
        {
            std::ostringstream level_solver_stream;
            level_solver_stream << d_level_solver_default_options_prefix << ln << "_";
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
        Pointer<StaggeredStokesPETScLevelSolver> p_level_solver = level_solver;
        if (p_level_solver)
            p_level_solver->addLinearOperator(d_SAJ_mat[ln]);
        else
            TBOX_ERROR("no mechanism for specifying SAJ part of operator!");
        level_solver->initializeSolverState(*getLevelSAMRAIVectorReal(*d_solution, ln),
                                            *getLevelSAMRAIVectorReal(*d_rhs, ln));
    }
    d_level_solvers[d_coarsest_ln] = d_coarse_solver;
    if (!d_level_solvers[d_coarsest_ln])
    {
        TBOX_ERROR("no mechanism for extracting coarse-grid solver Mat!");
    }

    // Get PETSc Mat from level solvers.
    d_level_mat.resize(d_finest_ln + 1);
    d_diagonal_level_mat.resize(d_finest_ln + 1);
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln - 1); ln <= std::min(d_finest_ln, finest_reset_ln); ++ln)
    {
        const KSP& level_ksp = d_level_solvers[ln]->getPETScKSP();
        ierr = KSPGetOperators(level_ksp, &d_level_mat[ln], NULL);
        IBTK_CHKERRQ(ierr);
        ierr = MatGetDiagonalBlock(d_level_mat[ln], &d_diagonal_level_mat[ln]);
        IBTK_CHKERRQ(ierr);
    }

    // Create subdomain matrices and KSP.
    d_subdomain_row_is.resize(d_finest_ln + 1);
    d_subdomain_col_is.resize(d_finest_ln + 1);
    d_no_subdomains.resize(d_finest_ln + 1);
    d_subdomain_mat.resize(d_finest_ln + 1);
    d_subdomain_bc_mat.resize(d_finest_ln + 1);
    d_subdomain_ksp.resize(d_finest_ln + 1);
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln - 1); ln <= std::min(d_finest_ln, finest_reset_ln); ++ln)
    {
        Pointer<StaggeredStokesPETScLevelSolver>& level_solver = d_level_solvers[ln];
        level_solver->getMSMSubdomains(&d_subdomain_row_is[ln], &d_subdomain_col_is[ln]);

        d_no_subdomains[ln] = static_cast<int>(d_subdomain_row_is[ln]->size());

        ierr = MatGetSubMatrices(d_level_mat[ln],
                                 d_no_subdomains[ln],
                                 &((*d_subdomain_row_is[ln])[0]),
                                 &((*d_subdomain_row_is[ln])[0]),
                                 MAT_INITIAL_MATRIX,
                                 &d_subdomain_mat[ln]);
        IBTK_CHKERRQ(ierr);

        ierr = MatGetSubMatrices(d_level_mat[ln],
                                 d_no_subdomains[ln],
                                 &((*d_subdomain_row_is[ln])[0]),
                                 &((*d_subdomain_col_is[ln])[0]),
                                 MAT_INITIAL_MATRIX,
                                 &d_subdomain_bc_mat[ln]);
        IBTK_CHKERRQ(ierr);

        // Set up subdomain KSPs
        d_subdomain_ksp[ln].resize(d_no_subdomains[ln]);
        for (int subdomain = 0; subdomain < d_no_subdomains[ln]; ++subdomain)
        {
            KSP& sub_ksp = d_subdomain_ksp[ln][subdomain];
            Mat& sub_mat = d_subdomain_mat[ln][subdomain];
            ierr = KSPCreate(PETSC_COMM_SELF, &sub_ksp);
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetOperators(sub_ksp, sub_mat, sub_mat);
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetReusePreconditioner(sub_ksp, PETSC_TRUE);
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetType(sub_ksp, KSPPREONLY);
            IBTK_CHKERRQ(ierr);
            PC sub_pc;
            ierr = KSPGetPC(sub_ksp, &sub_pc);
            IBTK_CHKERRQ(ierr);
            ierr = PCSetType(sub_pc, PCLU);
            IBTK_CHKERRQ(ierr);
            ierr = PCFactorReorderForNonzeroDiagonal(sub_pc, std::numeric_limits<double>::epsilon());
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetUp(sub_ksp);
            IBTK_CHKERRQ(ierr);
        }
    }

    // Create red and black subdomain matrices and KSP.
    d_red_subdomain_row_is.resize(d_finest_ln + 1);
    d_red_subdomain_col_is.resize(d_finest_ln + 1);
    d_no_red_subdomains.resize(d_finest_ln + 1);
    d_red_subdomain_mat.resize(d_finest_ln + 1);
    d_red_subdomain_bc_mat.resize(d_finest_ln + 1);
    d_red_subdomain_ksp.resize(d_finest_ln + 1);

    d_black_subdomain_row_is.resize(d_finest_ln + 1);
    d_black_subdomain_col_is.resize(d_finest_ln + 1);
    d_no_black_subdomains.resize(d_finest_ln + 1);
    d_black_subdomain_mat.resize(d_finest_ln + 1);
    d_black_subdomain_bc_mat.resize(d_finest_ln + 1);
    d_black_subdomain_ksp.resize(d_finest_ln + 1);
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln - 1); ln <= std::min(d_finest_ln, finest_reset_ln); ++ln)
    {
        Pointer<StaggeredStokesPETScLevelSolver>& level_solver = d_level_solvers[ln];
        level_solver->getMSMSubdomains(&d_red_subdomain_row_is[ln],
                                       &d_red_subdomain_col_is[ln],
                                       &d_black_subdomain_row_is[ln],
                                       &d_black_subdomain_col_is[ln]);
        d_no_red_subdomains[ln] = static_cast<int>(d_red_subdomain_row_is[ln]->size());
        d_no_black_subdomains[ln] = static_cast<int>(d_black_subdomain_row_is[ln]->size());

        // Extract red submatrices
        ierr = MatGetSubMatrices(d_level_mat[ln],
                                 d_no_red_subdomains[ln],
                                 &((*d_red_subdomain_row_is[ln])[0]),
                                 &((*d_red_subdomain_row_is[ln])[0]),
                                 MAT_INITIAL_MATRIX,
                                 &d_red_subdomain_mat[ln]);
        IBTK_CHKERRQ(ierr);
        ierr = MatGetSubMatrices(d_level_mat[ln],
                                 d_no_red_subdomains[ln],
                                 &((*d_red_subdomain_row_is[ln])[0]),
                                 &((*d_red_subdomain_col_is[ln])[0]),
                                 MAT_INITIAL_MATRIX,
                                 &d_red_subdomain_bc_mat[ln]);
        IBTK_CHKERRQ(ierr);

        // Extract black submatrices
        ierr = MatGetSubMatrices(d_level_mat[ln],
                                 d_no_black_subdomains[ln],
                                 &((*d_black_subdomain_row_is[ln])[0]),
                                 &((*d_black_subdomain_row_is[ln])[0]),
                                 MAT_INITIAL_MATRIX,
                                 &d_black_subdomain_mat[ln]);
        IBTK_CHKERRQ(ierr);
        ierr = MatGetSubMatrices(d_level_mat[ln],
                                 d_no_black_subdomains[ln],
                                 &((*d_black_subdomain_row_is[ln])[0]),
                                 &((*d_black_subdomain_col_is[ln])[0]),
                                 MAT_INITIAL_MATRIX,
                                 &d_black_subdomain_bc_mat[ln]);
        IBTK_CHKERRQ(ierr);

        // Set up red and black subdomain KSPs
        d_red_subdomain_ksp[ln].resize(d_no_red_subdomains[ln]);
        for (int red_subdomain = 0; red_subdomain < d_no_red_subdomains[ln]; ++red_subdomain)
        {
            KSP& sub_ksp = d_red_subdomain_ksp[ln][red_subdomain];
            Mat& sub_mat = d_red_subdomain_mat[ln][red_subdomain];
            ierr = KSPCreate(PETSC_COMM_SELF, &sub_ksp);
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetOperators(sub_ksp, sub_mat, sub_mat);
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetReusePreconditioner(sub_ksp, PETSC_TRUE);
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetType(sub_ksp, KSPPREONLY);
            IBTK_CHKERRQ(ierr);
            PC sub_pc;
            ierr = KSPGetPC(sub_ksp, &sub_pc);
            IBTK_CHKERRQ(ierr);
            ierr = PCSetType(sub_pc, PCLU);
            IBTK_CHKERRQ(ierr);
            ierr = PCFactorReorderForNonzeroDiagonal(sub_pc, std::numeric_limits<double>::epsilon());
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetUp(sub_ksp);
            IBTK_CHKERRQ(ierr);
        }

        d_black_subdomain_ksp[ln].resize(d_no_black_subdomains[ln]);
        for (int black_subdomain = 0; black_subdomain < d_no_black_subdomains[ln]; ++black_subdomain)
        {
            KSP& sub_ksp = d_black_subdomain_ksp[ln][black_subdomain];
            Mat& sub_mat = d_black_subdomain_mat[ln][black_subdomain];
            ierr = KSPCreate(PETSC_COMM_SELF, &sub_ksp);
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetOperators(sub_ksp, sub_mat, sub_mat);
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetReusePreconditioner(sub_ksp, PETSC_TRUE);
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetType(sub_ksp, KSPPREONLY);
            IBTK_CHKERRQ(ierr);
            PC sub_pc;
            ierr = KSPGetPC(sub_ksp, &sub_pc);
            IBTK_CHKERRQ(ierr);
            ierr = PCSetType(sub_pc, PCLU);
            IBTK_CHKERRQ(ierr);
            ierr = PCFactorReorderForNonzeroDiagonal(sub_pc, std::numeric_limits<double>::epsilon());
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetUp(sub_ksp);
            IBTK_CHKERRQ(ierr);
        }
    }

    return;
} // initializeOperatorStateSpecialized

void
StaggeredStokesIBBoxRelaxationFACOperator::deallocateOperatorStateSpecialized(const int coarsest_reset_ln,
                                                                              const int finest_reset_ln)
{
    int ierr;

    // Deallocate level solvers.
    if (d_coarse_solver == d_level_solvers[d_coarsest_ln]) d_coarse_solver.setNull();
    for (int ln = coarsest_reset_ln; ln <= std::min(d_finest_ln, finest_reset_ln); ++ln)
    {
        if (d_level_solvers[ln]) d_level_solvers[ln]->deallocateSolverState();
    }

    // Deallocate application ordering
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln - 1); ln <= std::min(d_finest_ln - 1, finest_reset_ln);
         ++ln)
    {
        ierr = AODestroy(&d_u_app_ordering[ln]);
        IBTK_CHKERRQ(ierr);
        d_u_app_ordering[ln] = NULL;
    }

    // Deallocate prolongation mat and restriction scaling mat
    for (int ln = std::min(d_finest_ln - 1, finest_reset_ln); ln >= std::max(d_coarsest_ln, coarsest_reset_ln - 1);
         --ln)
    {
        ierr = MatDestroy(&d_prolongation_mat[ln]);
        IBTK_CHKERRQ(ierr);
        d_prolongation_mat[ln] = NULL;

        ierr = VecDestroy(&d_scale_restriction_mat[ln]);
        IBTK_CHKERRQ(ierr);
        d_scale_restriction_mat[ln] = NULL;
    }

    // Deallocate SAJ mat.
    for (int ln = std::min(d_finest_ln, finest_reset_ln); ln >= std::max(d_coarsest_ln, coarsest_reset_ln - 1); --ln)
    {
        ierr = MatDestroy(&d_SAJ_mat[ln]);
        IBTK_CHKERRQ(ierr);
        d_SAJ_mat[ln] = NULL;
    }

    // Deallocate DOF index data.
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln - 1); ln <= std::min(d_finest_ln, finest_reset_ln); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_u_dof_index_idx)) level->deallocatePatchData(d_u_dof_index_idx);
        if (level->checkAllocated(d_p_dof_index_idx)) level->deallocatePatchData(d_p_dof_index_idx);
    }

    // Deallocate subdomain matrices and KSP.
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln - 1); ln <= std::min(d_finest_ln, finest_reset_ln); ++ln)
    {
        ierr = MatDestroyMatrices(d_no_subdomains[ln], &d_subdomain_mat[ln]);
        IBTK_CHKERRQ(ierr);
        d_subdomain_mat[ln] = NULL;
        ierr = MatDestroyMatrices(d_no_subdomains[ln], &d_subdomain_bc_mat[ln]);
        IBTK_CHKERRQ(ierr);
        d_subdomain_bc_mat[ln] = NULL;

        ierr = MatDestroyMatrices(d_no_red_subdomains[ln], &d_red_subdomain_mat[ln]);
        IBTK_CHKERRQ(ierr);
        d_red_subdomain_mat[ln] = NULL;
        ierr = MatDestroyMatrices(d_no_red_subdomains[ln], &d_red_subdomain_bc_mat[ln]);
        IBTK_CHKERRQ(ierr);
        d_red_subdomain_bc_mat[ln] = NULL;

        ierr = MatDestroyMatrices(d_no_black_subdomains[ln], &d_black_subdomain_mat[ln]);
        IBTK_CHKERRQ(ierr);
        d_black_subdomain_mat[ln] = NULL;
        ierr = MatDestroyMatrices(d_no_black_subdomains[ln], &d_black_subdomain_bc_mat[ln]);
        IBTK_CHKERRQ(ierr);
        d_black_subdomain_bc_mat[ln] = NULL;

        for (int subdomain = 0; subdomain < d_no_subdomains[ln]; ++subdomain)
        {
            KSP& sub_ksp = d_subdomain_ksp[ln][subdomain];
            ierr = KSPDestroy(&sub_ksp);
            IBTK_CHKERRQ(ierr);
            d_subdomain_ksp[ln][subdomain] = NULL;
        }

        for (int red_subdomain = 0; red_subdomain < d_no_red_subdomains[ln]; ++red_subdomain)
        {
            KSP& sub_ksp = d_red_subdomain_ksp[ln][red_subdomain];
            ierr = KSPDestroy(&sub_ksp);
            IBTK_CHKERRQ(ierr);
            d_red_subdomain_ksp[ln][red_subdomain] = NULL;
        }

        for (int black_subdomain = 0; black_subdomain < d_no_black_subdomains[ln]; ++black_subdomain)
        {
            KSP& sub_ksp = d_black_subdomain_ksp[ln][black_subdomain];
            ierr = KSPDestroy(&sub_ksp);
            IBTK_CHKERRQ(ierr);
            d_black_subdomain_ksp[ln][black_subdomain] = NULL;
        }
    }
    return;
} // deallocateOperatorStateSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
StaggeredStokesIBBoxRelaxationFACOperator::smoothErrorMultiplicative(SAMRAIVectorReal<NDIM, double>& error,
                                                                     const SAMRAIVectorReal<NDIM, double>& residual,
                                                                     int level_num,
                                                                     int num_sweeps,
                                                                     bool /*performing_pre_sweeps*/,
                                                                     bool /*performing_post_sweeps*/)
{
    if (num_sweeps == 0) return;

    int ierr;
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    const int U_error_idx = error.getComponentDescriptorIndex(0);
    const int P_error_idx = error.getComponentDescriptorIndex(1);
    const int U_residual_idx = residual.getComponentDescriptorIndex(0);
    const int P_residual_idx = residual.getComponentDescriptorIndex(1);

    // Create some Vecs.
    Mat& level_mat = d_level_mat[level_num];
    Mat& diagonal_level_mat = d_diagonal_level_mat[level_num];
    Vec err_vec, res_vec, new_res_vec;
    ierr = MatCreateVecs(level_mat, &err_vec, &res_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(res_vec, &new_res_vec);
    IBTK_CHKERRQ(ierr);
    Vec err_vec_local, new_res_vec_local;
    ierr = MatCreateVecs(diagonal_level_mat, &err_vec_local, &new_res_vec_local);
    IBTK_CHKERRQ(ierr);

    // Copy patch level SAMRAI data to PETSc Vec
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        err_vec, U_error_idx, d_u_dof_index_idx, P_error_idx, d_p_dof_index_idx, level);
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        res_vec, U_residual_idx, d_u_dof_index_idx, P_residual_idx, d_p_dof_index_idx, level);

    // DOFs on this processor.
    const int mpi_rank = SAMRAI_MPI::getRank();
    const std::vector<int>& dofs_level = d_num_dofs_per_proc[level_num];
    const int first_local_dof = std::accumulate(dofs_level.begin(), dofs_level.begin() + mpi_rank, 0);

    // Smooth the error by the specified number of sweeps.
    for (int isweep = 0; isweep < num_sweeps; ++isweep)
    {
        // Extract the local vectors.
        PetscScalar *err_array, *new_res_array;
        ierr = VecGetArray(err_vec, &err_array);
        IBTK_CHKERRQ(ierr);
        ierr = VecGetArray(new_res_vec, &new_res_array);
        IBTK_CHKERRQ(ierr);

        // Update residual
        // Do globally: r = f - A*u
        ierr = VecScale(err_vec, -1.0);
        IBTK_CHKERRQ(ierr);
        ierr = MatMultAdd(level_mat, err_vec, res_vec, new_res_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecScale(err_vec, -1.0);
        IBTK_CHKERRQ(ierr);

        // Do locally: r = r + A*u
        ierr = VecGetLocalVector(err_vec, err_vec_local);
        IBTK_CHKERRQ(ierr);
        ierr = VecGetLocalVector(new_res_vec, new_res_vec_local);
        IBTK_CHKERRQ(ierr);

        ierr = MatMultAdd(diagonal_level_mat, err_vec_local, new_res_vec_local, new_res_vec_local);
        IBTK_CHKERRQ(ierr);

        ierr = VecRestoreLocalVector(err_vec, err_vec_local);
        IBTK_CHKERRQ(ierr);
        ierr = VecRestoreLocalVector(new_res_vec, new_res_vec_local);
        IBTK_CHKERRQ(ierr);

        // Smooth error on subdomains.
        const int n_subdomains = d_no_subdomains[level_num];
        for (int subdomain = 0; subdomain < n_subdomains; ++subdomain)
        {
            Mat& bc_mat = d_subdomain_bc_mat[level_num][subdomain];
            KSP& sub_ksp = d_subdomain_ksp[level_num][subdomain];

            // Get local and non-local DOFs.
            // Non-local DOFs are non-subdomains DOFs that are on this processor.
            IS& local_dofs = (*d_subdomain_row_is[level_num])[subdomain];
            IS& nonlocal_dofs = (*d_subdomain_col_is[level_num])[subdomain];

            PetscInt is_size_local, is_size_nonlocal;
            ierr = ISGetLocalSize(local_dofs, &is_size_local);
            IBTK_CHKERRQ(ierr);
            ierr = ISGetLocalSize(nonlocal_dofs, &is_size_nonlocal);
            IBTK_CHKERRQ(ierr);
            std::vector<int> local_indices(is_size_local);
            for (int i = 0; i < is_size_local; ++i)
            {
                local_indices[i] = i;
            }
            std::vector<int> nonlocal_indices(is_size_nonlocal);
            for (int i = 0; i < is_size_nonlocal; ++i)
            {
                nonlocal_indices[i] = i;
            }

            const PetscInt *is_local_array, *is_nonlocal_array;
            ierr = ISGetIndices(local_dofs, &is_local_array);
            IBTK_CHKERRQ(ierr);
            ierr = ISGetIndices(nonlocal_dofs, &is_nonlocal_array);
            IBTK_CHKERRQ(ierr);

            // Create little block Vecs for solving local sytem and modifying RHS.
            Vec u, r, e;
            ierr = MatCreateVecs(bc_mat, &u, &r);
            IBTK_CHKERRQ(ierr);
            ierr = VecDuplicate(r, &e);
            IBTK_CHKERRQ(ierr);

#if !defined(NDEBUG)
            PetscInt vec_size;

            ierr = VecGetSize(r, &vec_size);
            IBTK_CHKERRQ(ierr);
            TBOX_ASSERT(is_size_local == vec_size);

            ierr = VecGetSize(u, &vec_size);
            IBTK_CHKERRQ(ierr);
            TBOX_ASSERT(is_size_nonlocal == vec_size);
#endif

            // Copy to/from various Vecs
            std::vector<double> local_values(is_size_local);
            std::vector<double> nonlocal_values(is_size_nonlocal);

            for (int i = 0; i < is_size_local; ++i)
            {
                local_values[i] = new_res_array[is_local_array[i] - first_local_dof];
            }
            ierr = VecSetValues(r, is_size_local, &local_indices[0], &local_values[0], INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyBegin(r);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyEnd(r);
            IBTK_CHKERRQ(ierr);

            for (int i = 0; i < is_size_nonlocal; ++i)
            {
                nonlocal_values[i] = err_array[is_nonlocal_array[i] - first_local_dof];
            }
            ierr = VecSetValues(u, is_size_nonlocal, &nonlocal_indices[0], &nonlocal_values[0], INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyBegin(u);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyEnd(u);
            IBTK_CHKERRQ(ierr);

            // Modify RHS for "boundary" conditions.
            ierr = VecScale(u, -1.0);
            IBTK_CHKERRQ(ierr);
            ierr = MatMultAdd(bc_mat, u, r, r);
            IBTK_CHKERRQ(ierr);

            // Do the local solve.
            ierr = KSPSolve(sub_ksp, r, e);
            IBTK_CHKERRQ(ierr);

            // Update error vector.
            ierr = VecGetValues(e, is_size_local, &local_indices[0], &local_values[0]);
            IBTK_CHKERRQ(ierr);
            for (int i = 0; i < is_size_local; ++i)
            {
                err_array[is_local_array[i] - first_local_dof] = local_values[i];
            }

            // Restore IS indices.
            ierr = ISRestoreIndices(local_dofs, &is_local_array);
            IBTK_CHKERRQ(ierr);
            ierr = ISRestoreIndices(nonlocal_dofs, &is_nonlocal_array);
            IBTK_CHKERRQ(ierr);

            // Destroy temporary PETSc objects
            ierr = VecDestroy(&u);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&r);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&e);
            IBTK_CHKERRQ(ierr);
        }

        ierr = VecRestoreArray(err_vec, &err_array);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyBegin(err_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyEnd(err_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecRestoreArray(new_res_vec, &new_res_array);
        IBTK_CHKERRQ(ierr);
    }

    StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(
        err_vec, U_error_idx, d_u_dof_index_idx, P_error_idx, d_p_dof_index_idx, level, NULL, NULL);
    ierr = VecDestroy(&err_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&res_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&new_res_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&err_vec_local);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&new_res_vec_local);
    IBTK_CHKERRQ(ierr);

    // Synchronize data along patch boundaries.
    xeqScheduleDataSynch(U_error_idx, level_num);

    // Fill ghost cells.
    const std::pair<int, int> error_idxs = std::make_pair(U_error_idx, P_error_idx);
    xeqScheduleGhostFillNoCoarse(error_idxs, level_num);

    return;
} // smoothErrorMultiplicative

void
StaggeredStokesIBBoxRelaxationFACOperator::smoothErrorRedBlackMultiplicative(
    SAMRAIVectorReal<NDIM, double>& error,
    const SAMRAIVectorReal<NDIM, double>& residual,
    int level_num,
    int num_sweeps,
    bool /*performing_pre_sweeps*/,
    bool /*performing_post_sweeps*/)
{
    if (num_sweeps == 0) return;

    int ierr;
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    const int U_error_idx = error.getComponentDescriptorIndex(0);
    const int P_error_idx = error.getComponentDescriptorIndex(1);
    const int U_residual_idx = residual.getComponentDescriptorIndex(0);
    const int P_residual_idx = residual.getComponentDescriptorIndex(1);

    // Create some Vecs.
    Mat& level_mat = d_level_mat[level_num];
    Mat& diagonal_level_mat = d_diagonal_level_mat[level_num];

    Vec err_vec, res_vec, new_res_vec;
    PetscScalar *err_array, *new_res_array;
    ierr = MatCreateVecs(level_mat, &err_vec, &res_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(res_vec, &new_res_vec);
    IBTK_CHKERRQ(ierr);

    Vec err_vec_local, new_res_vec_local;
    ierr = MatCreateVecs(diagonal_level_mat, &err_vec_local, &new_res_vec_local);
    IBTK_CHKERRQ(ierr);

    // Copy patch level SAMRAI data to PETSc Vec
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        err_vec, U_error_idx, d_u_dof_index_idx, P_error_idx, d_p_dof_index_idx, level);
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        res_vec, U_residual_idx, d_u_dof_index_idx, P_residual_idx, d_p_dof_index_idx, level);

    // DOFs on this processor.
    const int mpi_rank = SAMRAI_MPI::getRank();
    const std::vector<int>& dofs_level = d_num_dofs_per_proc[level_num];
    const int first_local_dof = std::accumulate(dofs_level.begin(), dofs_level.begin() + mpi_rank, 0);

    // Smooth the error by the specified number of sweeps.
    for (int isweep = 0; isweep < num_sweeps; ++isweep)
    {
        // Update residual for off-processor DOFs
        // a) Do globally: r = f - A*u
        ierr = VecScale(err_vec, -1.0);
        IBTK_CHKERRQ(ierr);
        ierr = MatMultAdd(level_mat, err_vec, res_vec, new_res_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecScale(err_vec, -1.0);
        IBTK_CHKERRQ(ierr);

        // b) Do locally: r = r + A*u
        ierr = VecGetLocalVector(err_vec, err_vec_local);
        IBTK_CHKERRQ(ierr);
        ierr = VecGetLocalVector(new_res_vec, new_res_vec_local);
        IBTK_CHKERRQ(ierr);
        ierr = MatMultAdd(diagonal_level_mat, err_vec_local, new_res_vec_local, new_res_vec_local);
        IBTK_CHKERRQ(ierr);
        ierr = VecRestoreLocalVector(err_vec, err_vec_local);
        IBTK_CHKERRQ(ierr);
        ierr = VecRestoreLocalVector(new_res_vec, new_res_vec_local);
        IBTK_CHKERRQ(ierr);

        // Extract the local vectors.
        ierr = VecGetArray(err_vec, &err_array);
        IBTK_CHKERRQ(ierr);
        ierr = VecGetArray(new_res_vec, &new_res_array);
        IBTK_CHKERRQ(ierr);

        // Smooth error on red subdomains.
        const int n_red_subdomains = d_no_red_subdomains[level_num];
        for (int red_subdomain = 0; red_subdomain < n_red_subdomains; ++red_subdomain)
        {
            Mat& bc_mat = d_red_subdomain_bc_mat[level_num][red_subdomain];
            KSP& sub_ksp = d_red_subdomain_ksp[level_num][red_subdomain];

            // Get local (to box) and non-local (non-box) DOFs.
            // Non-local DOFs are non-subdomains or the non-box DOFs that
            // are still on this processor.
            IS& local_dofs = (*d_red_subdomain_row_is[level_num])[red_subdomain];
            IS& nonlocal_dofs = (*d_red_subdomain_col_is[level_num])[red_subdomain];

            PetscInt is_size_local, is_size_nonlocal;
            ierr = ISGetLocalSize(local_dofs, &is_size_local);
            IBTK_CHKERRQ(ierr);
            ierr = ISGetLocalSize(nonlocal_dofs, &is_size_nonlocal);
            IBTK_CHKERRQ(ierr);
            std::vector<int> local_indices(is_size_local);
            for (int i = 0; i < is_size_local; ++i)
            {
                local_indices[i] = i;
            }
            std::vector<int> nonlocal_indices(is_size_nonlocal);
            for (int i = 0; i < is_size_nonlocal; ++i)
            {
                nonlocal_indices[i] = i;
            }

            const PetscInt *is_local_array, *is_nonlocal_array;
            ierr = ISGetIndices(local_dofs, &is_local_array);
            IBTK_CHKERRQ(ierr);
            ierr = ISGetIndices(nonlocal_dofs, &is_nonlocal_array);
            IBTK_CHKERRQ(ierr);

            // Create little block Vecs for solving local sytem and modifying RHS.
            Vec u, r, e;
            ierr = MatCreateVecs(bc_mat, &u, &r);
            IBTK_CHKERRQ(ierr);
            ierr = VecDuplicate(r, &e);
            IBTK_CHKERRQ(ierr);

#if !defined(NDEBUG)
            PetscInt vec_size;

            ierr = VecGetSize(r, &vec_size);
            IBTK_CHKERRQ(ierr);
            TBOX_ASSERT(is_size_local == vec_size);

            ierr = VecGetSize(u, &vec_size);
            IBTK_CHKERRQ(ierr);
            TBOX_ASSERT(is_size_nonlocal == vec_size);
#endif

            // Copy to/from various Vecs
            std::vector<double> local_values(is_size_local);
            std::vector<double> nonlocal_values(is_size_nonlocal);

            for (int i = 0; i < is_size_local; ++i)
            {
                local_values[i] = new_res_array[is_local_array[i] - first_local_dof];
            }
            ierr = VecSetValues(r, is_size_local, &local_indices[0], &local_values[0], INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyBegin(r);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyEnd(r);
            IBTK_CHKERRQ(ierr);

            for (int i = 0; i < is_size_nonlocal; ++i)
            {
                nonlocal_values[i] = err_array[is_nonlocal_array[i] - first_local_dof];
            }
            ierr = VecSetValues(u, is_size_nonlocal, &nonlocal_indices[0], &nonlocal_values[0], INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyBegin(u);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyEnd(u);
            IBTK_CHKERRQ(ierr);

            // Modify RHS for "boundary" conditions.
            ierr = VecScale(u, -1.0);
            IBTK_CHKERRQ(ierr);
            ierr = MatMultAdd(bc_mat, u, r, r);
            IBTK_CHKERRQ(ierr);

            // Do the local solve.
            ierr = KSPSolve(sub_ksp, r, e);
            IBTK_CHKERRQ(ierr);

            // Update error vector.
            ierr = VecGetValues(e, is_size_local, &local_indices[0], &local_values[0]);
            IBTK_CHKERRQ(ierr);
            for (int i = 0; i < is_size_local; ++i)
            {
                err_array[is_local_array[i] - first_local_dof] = local_values[i];
            }

            // Restore IS indices.
            ierr = ISRestoreIndices(local_dofs, &is_local_array);
            IBTK_CHKERRQ(ierr);
            ierr = ISRestoreIndices(nonlocal_dofs, &is_nonlocal_array);
            IBTK_CHKERRQ(ierr);

            // Destroy temporary PETSc objects
            ierr = VecDestroy(&u);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&r);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&e);
            IBTK_CHKERRQ(ierr);
        }

        // Assemble error vector for contribution from red subdomains.
        ierr = VecRestoreArray(err_vec, &err_array);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyBegin(err_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyEnd(err_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecRestoreArray(new_res_vec, &new_res_array);
        IBTK_CHKERRQ(ierr);

        // Update residual for black subdomains

        // Do globally: r = f - A*u
        ierr = VecScale(err_vec, -1.0);
        IBTK_CHKERRQ(ierr);
        ierr = MatMultAdd(level_mat, err_vec, res_vec, new_res_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecScale(err_vec, -1.0);
        IBTK_CHKERRQ(ierr);

        // Do locally: r = r + A*u
        ierr = VecGetLocalVector(err_vec, err_vec_local);
        IBTK_CHKERRQ(ierr);
        ierr = VecGetLocalVector(new_res_vec, new_res_vec_local);
        IBTK_CHKERRQ(ierr);
        ierr = MatMultAdd(diagonal_level_mat, err_vec_local, new_res_vec_local, new_res_vec_local);
        IBTK_CHKERRQ(ierr);
        ierr = VecRestoreLocalVector(err_vec, err_vec_local);
        IBTK_CHKERRQ(ierr);
        ierr = VecRestoreLocalVector(new_res_vec, new_res_vec_local);
        IBTK_CHKERRQ(ierr);

        // Extract the local vectors.
        ierr = VecGetArray(err_vec, &err_array);
        IBTK_CHKERRQ(ierr);
        ierr = VecGetArray(new_res_vec, &new_res_array);
        IBTK_CHKERRQ(ierr);

        // Smooth error on black subdomains.
        const int n_black_subdomains = d_no_black_subdomains[level_num];
        for (int black_subdomain = 0; black_subdomain < n_black_subdomains; ++black_subdomain)
        {
            Mat& bc_mat = d_black_subdomain_bc_mat[level_num][black_subdomain];
            KSP& sub_ksp = d_black_subdomain_ksp[level_num][black_subdomain];

            // Get local (to box) and non-local (non-box) DOFs.
            // Non-local DOFs are non-subdomains or the non-box DOFs that
            // are still on this processor.
            IS& local_dofs = (*d_black_subdomain_row_is[level_num])[black_subdomain];
            IS& nonlocal_dofs = (*d_black_subdomain_col_is[level_num])[black_subdomain];

            PetscInt is_size_local, is_size_nonlocal;
            ierr = ISGetLocalSize(local_dofs, &is_size_local);
            IBTK_CHKERRQ(ierr);
            ierr = ISGetLocalSize(nonlocal_dofs, &is_size_nonlocal);
            IBTK_CHKERRQ(ierr);
            std::vector<int> local_indices(is_size_local);
            for (int i = 0; i < is_size_local; ++i)
            {
                local_indices[i] = i;
            }
            std::vector<int> nonlocal_indices(is_size_nonlocal);
            for (int i = 0; i < is_size_nonlocal; ++i)
            {
                nonlocal_indices[i] = i;
            }

            const PetscInt *is_local_array, *is_nonlocal_array;
            ierr = ISGetIndices(local_dofs, &is_local_array);
            IBTK_CHKERRQ(ierr);
            ierr = ISGetIndices(nonlocal_dofs, &is_nonlocal_array);
            IBTK_CHKERRQ(ierr);

            // Create little block Vecs for solving local sytem and modifying RHS.
            Vec u, r, e;
            ierr = MatCreateVecs(bc_mat, &u, &r);
            IBTK_CHKERRQ(ierr);
            ierr = VecDuplicate(r, &e);
            IBTK_CHKERRQ(ierr);

#if !defined(NDEBUG)
            PetscInt vec_size;

            ierr = VecGetSize(r, &vec_size);
            IBTK_CHKERRQ(ierr);
            TBOX_ASSERT(is_size_local == vec_size);

            ierr = VecGetSize(u, &vec_size);
            IBTK_CHKERRQ(ierr);
            TBOX_ASSERT(is_size_nonlocal == vec_size);
#endif

            // Copy to/from various Vecs
            std::vector<double> local_values(is_size_local);
            std::vector<double> nonlocal_values(is_size_nonlocal);

            for (int i = 0; i < is_size_local; ++i)
            {
                local_values[i] = new_res_array[is_local_array[i] - first_local_dof];
            }
            ierr = VecSetValues(r, is_size_local, &local_indices[0], &local_values[0], INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyBegin(r);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyEnd(r);
            IBTK_CHKERRQ(ierr);

            for (int i = 0; i < is_size_nonlocal; ++i)
            {
                nonlocal_values[i] = err_array[is_nonlocal_array[i] - first_local_dof];
            }
            ierr = VecSetValues(u, is_size_nonlocal, &nonlocal_indices[0], &nonlocal_values[0], INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyBegin(u);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyEnd(u);
            IBTK_CHKERRQ(ierr);

            // Modify RHS for "boundary" conditions.
            ierr = VecScale(u, -1.0);
            IBTK_CHKERRQ(ierr);
            ierr = MatMultAdd(bc_mat, u, r, r);
            IBTK_CHKERRQ(ierr);

            // Do the local solve.
            ierr = KSPSolve(sub_ksp, r, e);
            IBTK_CHKERRQ(ierr);

            // Update error vector.
            ierr = VecGetValues(e, is_size_local, &local_indices[0], &local_values[0]);
            IBTK_CHKERRQ(ierr);
            for (int i = 0; i < is_size_local; ++i)
            {
                err_array[is_local_array[i] - first_local_dof] = local_values[i];
            }

            // Restore IS indices.
            ierr = ISRestoreIndices(local_dofs, &is_local_array);
            IBTK_CHKERRQ(ierr);
            ierr = ISRestoreIndices(nonlocal_dofs, &is_nonlocal_array);
            IBTK_CHKERRQ(ierr);

            // Destroy temporary PETSc objects
            ierr = VecDestroy(&u);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&r);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&e);
            IBTK_CHKERRQ(ierr);
        }

        // Assemble error vector for contribution from black subdomains.
        ierr = VecRestoreArray(err_vec, &err_array);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyBegin(err_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyEnd(err_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecRestoreArray(new_res_vec, &new_res_array);
        IBTK_CHKERRQ(ierr);
    }

    StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(
        err_vec, U_error_idx, d_u_dof_index_idx, P_error_idx, d_p_dof_index_idx, level, NULL, NULL);
    ierr = VecDestroy(&err_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&res_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&new_res_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&err_vec_local);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&new_res_vec_local);
    IBTK_CHKERRQ(ierr);

    // Synchronize data along patch boundaries.
    xeqScheduleDataSynch(U_error_idx, level_num);

    // Fill ghost cells.
    const std::pair<int, int> error_idxs = std::make_pair(U_error_idx, P_error_idx);
    xeqScheduleGhostFillNoCoarse(error_idxs, level_num);

    return;
} // smoothErrorRedBlackMultiplicative

void
StaggeredStokesIBBoxRelaxationFACOperator::smoothErrorAdditive(SAMRAIVectorReal<NDIM, double>& error,
                                                               const SAMRAIVectorReal<NDIM, double>& residual,
                                                               int level_num,
                                                               int num_sweeps,
                                                               bool /*performing_pre_sweeps*/,
                                                               bool /*performing_post_sweeps*/)
{
    if (num_sweeps == 0) return;

    int ierr;
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    const int U_error_idx = error.getComponentDescriptorIndex(0);
    const int P_error_idx = error.getComponentDescriptorIndex(1);
    const int U_residual_idx = residual.getComponentDescriptorIndex(0);
    const int P_residual_idx = residual.getComponentDescriptorIndex(1);

    // Create some Vecs.
    Mat& level_mat = d_level_mat[level_num];
    Vec err_vec, old_err_vec, res_vec;
    ierr = MatCreateVecs(level_mat, &err_vec, &res_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(err_vec, &old_err_vec);
    IBTK_CHKERRQ(ierr);

    // Copy patch level SAMRAI data to PETSc Vec
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        err_vec, U_error_idx, d_u_dof_index_idx, P_error_idx, d_p_dof_index_idx, level);
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        res_vec, U_residual_idx, d_u_dof_index_idx, P_residual_idx, d_p_dof_index_idx, level);

    // DOFs on this processor.
    const int mpi_rank = SAMRAI_MPI::getRank();
    const std::vector<int>& dofs_level = d_num_dofs_per_proc[level_num];
    const int first_local_dof = std::accumulate(dofs_level.begin(), dofs_level.begin() + mpi_rank, 0);

    // Smooth the error by the specified number of sweeps.
    for (int isweep = 0; isweep < num_sweeps; ++isweep)
    {
        ierr = VecCopy(err_vec, old_err_vec);
        IBTK_CHKERRQ(ierr);

        // Extract the local vectors.
        PetscScalar *err_array, *old_err_array, *res_array;
        ierr = VecGetArray(err_vec, &err_array);
        IBTK_CHKERRQ(ierr);
        ierr = VecGetArray(old_err_vec, &old_err_array);
        IBTK_CHKERRQ(ierr);
        ierr = VecGetArray(res_vec, &res_array);
        IBTK_CHKERRQ(ierr);

        // Smooth error on subdomains.
        const int n_subdomains = d_no_subdomains[level_num];
        for (int subdomain = 0; subdomain < n_subdomains; ++subdomain)
        {
            Mat& bc_mat = d_subdomain_bc_mat[level_num][subdomain];
            KSP& sub_ksp = d_subdomain_ksp[level_num][subdomain];

            // Get local and non-local DOFs.
            // Non-local DOFs are non-subdomains DOFs that are on this processor.
            IS& local_dofs = (*d_subdomain_row_is[level_num])[subdomain];
            IS& nonlocal_dofs = (*d_subdomain_col_is[level_num])[subdomain];

            PetscInt is_size_local, is_size_nonlocal;
            ierr = ISGetLocalSize(local_dofs, &is_size_local);
            IBTK_CHKERRQ(ierr);
            ierr = ISGetLocalSize(nonlocal_dofs, &is_size_nonlocal);
            IBTK_CHKERRQ(ierr);
            std::vector<int> local_indices(is_size_local);
            for (int i = 0; i < is_size_local; ++i)
            {
                local_indices[i] = i;
            }
            std::vector<int> nonlocal_indices(is_size_nonlocal);
            for (int i = 0; i < is_size_nonlocal; ++i)
            {
                nonlocal_indices[i] = i;
            }

            const PetscInt *is_local_array, *is_nonlocal_array;
            ierr = ISGetIndices(local_dofs, &is_local_array);
            IBTK_CHKERRQ(ierr);
            ierr = ISGetIndices(nonlocal_dofs, &is_nonlocal_array);
            IBTK_CHKERRQ(ierr);

            // Create little block Vecs for solving local sytem and modifying RHS.
            Vec u, r, e;
            ierr = MatCreateVecs(bc_mat, &u, &r);
            IBTK_CHKERRQ(ierr);
            ierr = VecDuplicate(r, &e);
            IBTK_CHKERRQ(ierr);

#if !defined(NDEBUG)
            PetscInt vec_size;

            ierr = VecGetSize(r, &vec_size);
            IBTK_CHKERRQ(ierr);
            TBOX_ASSERT(is_size_local == vec_size);

            ierr = VecGetSize(u, &vec_size);
            IBTK_CHKERRQ(ierr);
            TBOX_ASSERT(is_size_nonlocal == vec_size);
#endif

            // Copy to/from various Vecs
            std::vector<double> local_values(is_size_local);
            std::vector<double> nonlocal_values(is_size_nonlocal);

            for (int i = 0; i < is_size_local; ++i)
            {
                local_values[i] = res_array[is_local_array[i] - first_local_dof];
            }
            ierr = VecSetValues(r, is_size_local, &local_indices[0], &local_values[0], INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyBegin(r);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyEnd(r);
            IBTK_CHKERRQ(ierr);

            for (int i = 0; i < is_size_nonlocal; ++i)
            {
                nonlocal_values[i] = old_err_array[is_nonlocal_array[i] - first_local_dof];
            }
            ierr = VecSetValues(u, is_size_nonlocal, &nonlocal_indices[0], &nonlocal_values[0], INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyBegin(u);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyEnd(u);
            IBTK_CHKERRQ(ierr);

            // Modify RHS for "boundary" conditions.
            ierr = VecScale(u, -1.0);
            IBTK_CHKERRQ(ierr);
            ierr = MatMultAdd(bc_mat, u, r, r);
            IBTK_CHKERRQ(ierr);

            // Do the local solve.
            ierr = KSPSolve(sub_ksp, r, e);
            IBTK_CHKERRQ(ierr);

            // Update error vector.
            ierr = VecGetValues(e, is_size_local, &local_indices[0], &local_values[0]);
            IBTK_CHKERRQ(ierr);
            for (int i = 0; i < is_size_local; ++i)
            {
                err_array[is_local_array[i] - first_local_dof] = local_values[i];
            }

            // Restore IS indices.
            ierr = ISRestoreIndices(local_dofs, &is_local_array);
            IBTK_CHKERRQ(ierr);
            ierr = ISRestoreIndices(nonlocal_dofs, &is_nonlocal_array);
            IBTK_CHKERRQ(ierr);

            // Destroy temporary PETSc objects
            ierr = VecDestroy(&u);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&r);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&e);
            IBTK_CHKERRQ(ierr);
        }

        ierr = VecRestoreArray(err_vec, &err_array);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyBegin(err_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyEnd(err_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecRestoreArray(old_err_vec, &old_err_array);
    }

    StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(
        err_vec, U_error_idx, d_u_dof_index_idx, P_error_idx, d_p_dof_index_idx, level, NULL, NULL);
    ierr = VecDestroy(&err_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&res_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&old_err_vec);
    IBTK_CHKERRQ(ierr);

    // Synchronize data along patch boundaries.
    xeqScheduleDataSynch(U_error_idx, level_num);

    // Fill ghost cells.
    const std::pair<int, int> error_idxs = std::make_pair(U_error_idx, P_error_idx);
    xeqScheduleGhostFillNoCoarse(error_idxs, level_num);

    return;
} // smoothErrorAdditive

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
