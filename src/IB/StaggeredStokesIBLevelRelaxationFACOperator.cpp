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
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesIBLevelRelaxationFACOperator::StaggeredStokesIBLevelRelaxationFACOperator(
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
} // StaggeredStokesIBLevelRelaxationFACOperator

StaggeredStokesIBLevelRelaxationFACOperator::~StaggeredStokesIBLevelRelaxationFACOperator()
{
    if (d_is_initialized)
    {
        deallocateOperatorState();
    }
    return;
} // ~StaggeredStokesIBLevelRelaxationFACOperator

void
StaggeredStokesIBLevelRelaxationFACOperator::setIBForceJacobian(Mat& A)
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
StaggeredStokesIBLevelRelaxationFACOperator::setIBInterpOp(Mat& J)
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
StaggeredStokesIBLevelRelaxationFACOperator::getStaggeredStokesPETScLevelSolver(const int ln) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(ln >= 0 && ln <= d_finest_ln);
#endif
    return d_level_solvers[ln];
} // getStaggeredStokesIBLevelSolver

Mat
StaggeredStokesIBLevelRelaxationFACOperator::getEulerianElasticityLevelOp(const int ln) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(ln >= 0 && ln <= d_finest_ln);
#endif
    return d_SAJ_mat[ln];
} // getEulerianElasticityLevelOp

Mat
StaggeredStokesIBLevelRelaxationFACOperator::getProlongationOp(const int ln) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(ln >= 0 && ln < d_finest_ln);
#endif
    return d_prolongation_mat[ln];
} // getProlongationOp

Vec
StaggeredStokesIBLevelRelaxationFACOperator::getRestrictionScalingOp(const int ln) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(ln >= 0 && ln < d_finest_ln);
#endif
    return d_scale_restriction_mat[ln];
} // getRestrictionScalingOp

void
StaggeredStokesIBLevelRelaxationFACOperator::computeResidual(SAMRAIVectorReal<NDIM, double>& residual,
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
StaggeredStokesIBLevelRelaxationFACOperator::smoothError(SAMRAIVectorReal<NDIM, double>& error,
                                                         const SAMRAIVectorReal<NDIM, double>& residual,
                                                         int level_num,
                                                         int num_sweeps,
                                                         bool /*performing_pre_sweeps*/,
                                                         bool /*performing_post_sweeps*/)
{
    if (num_sweeps == 0) return;

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
                                  d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                                  IntVector<NDIM>(0));
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
    return;
} // smoothError

/////////////////////////////// PROTECTED ////////////////////////////////////

void
StaggeredStokesIBLevelRelaxationFACOperator::initializeOperatorStateSpecialized(
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

    return;
} // initializeOperatorStateSpecialized

void
StaggeredStokesIBLevelRelaxationFACOperator::deallocateOperatorStateSpecialized(const int coarsest_reset_ln,
                                                                                const int finest_reset_ln)
{
    int ierr;

    // Deallocate level solvers and overlap boxes for side and cell data.
    if (d_coarse_solver == d_level_solvers[d_coarsest_ln]) d_coarse_solver.setNull();
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
    return;
} // deallocateOperatorStateSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
