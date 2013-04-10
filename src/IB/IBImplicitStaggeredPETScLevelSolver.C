// Filename: IBImplicitStaggeredPETScLevelSolver.C
// Created on 16 Apr 2012 by Boyce Griffith
//
// Copyright (c) 2002-2013, Boyce Griffith
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

#include "IBImplicitStaggeredPETScLevelSolver.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/StaggeredStokesPETScMatUtilities.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>
#include <ibamr/namespaces.h>
#include <ibamr/ibamr_utilities.h>

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScMatUtilities.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <HierarchyDataOpsManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);
static const int SIDEG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBImplicitStaggeredPETScLevelSolver::IBImplicitStaggeredPETScLevelSolver(
    const std::string& object_name,
    Pointer<Database> input_db,
    const std::string& default_options_prefix)
    : d_problem_coefs(),
      d_J_mat(NULL),
      d_interp_fcn(NULL),
      d_interp_stencil(0),
      d_X_vec(NULL),
      d_default_u_bc_coef(new LocationIndexRobinBcCoefs<NDIM>(d_object_name+"::default_u_bc_coef", Pointer<Database>(NULL))),
      d_u_bc_coefs(),
      d_context(NULL),
      d_u_dof_index_idx(-1),
      d_p_dof_index_idx(-1),
      d_u_dof_index_var(NULL),
      d_p_dof_index_var(NULL),
      d_data_synch_sched(NULL),
      d_ghost_fill_sched(NULL)
{
    GeneralSolver::init(object_name, /*homogeneous_bc*/ false);
    PETScLevelSolver::init(input_db, default_options_prefix);

    // Setup a default boundary condition object that specifies homogeneous
    // Dirichlet boundary conditions.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_default_u_bc_coef->setBoundaryValue(2*d  ,0.0);
        d_default_u_bc_coef->setBoundaryValue(2*d+1,0.0);
    }

#if 0  // XXXX
    // Initialize the boundary conditions objects.
    setHomogeneousBc(d_homogeneous_bc);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (u_bc_coefs[d])
        {
            d_u_bc_coefs[d] = u_bc_coefs[d];
        }
        else
        {
            d_u_bc_coefs[d] = d_default_u_bc_coef;
        }
    }
#endif

    // Construct the DOF index variable/context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(d_object_name + "::CONTEXT");
    d_u_dof_index_var = new SideVariable<NDIM,int>(d_object_name + "::u_dof_index");
    if (var_db->checkVariableExists(d_u_dof_index_var->getName()))
    {
        d_u_dof_index_var = var_db->getVariable(d_u_dof_index_var->getName());
        d_u_dof_index_idx = var_db->mapVariableAndContextToIndex(d_u_dof_index_var, d_context);
        var_db->removePatchDataIndex(d_u_dof_index_idx);
    }
    d_u_dof_index_idx = var_db->registerVariableAndContext(d_u_dof_index_var, d_context, SIDEG);
    d_p_dof_index_var = new CellVariable<NDIM,int>(d_object_name + "::p_dof_index");
    if (var_db->checkVariableExists(d_p_dof_index_var->getName()))
    {
        d_p_dof_index_var = var_db->getVariable(d_p_dof_index_var->getName());
        d_p_dof_index_idx = var_db->mapVariableAndContextToIndex(d_p_dof_index_var, d_context);
        var_db->removePatchDataIndex(d_p_dof_index_idx);
    }
    d_p_dof_index_idx = var_db->registerVariableAndContext(d_p_dof_index_var, d_context, CELLG);
    return;
}// IBImplicitStaggeredPETScLevelSolver

IBImplicitStaggeredPETScLevelSolver::~IBImplicitStaggeredPETScLevelSolver()
{
    if (d_is_initialized) deallocateSolverState();
    delete d_default_u_bc_coef;
    return;
}// ~IBImplicitStaggeredPETScLevelSolver

void
IBImplicitStaggeredPETScLevelSolver::initializeOperator()
{
    int ierr;

    // Setup PETSc objects.
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_num);
    PETScMatUtilities::constructPatchLevelSCInterpOp(d_R_mat, d_interp_fcn, d_interp_stencil, *d_X_vec, d_num_dofs_per_proc, d_u_dof_index_idx, level);
    Pointer<CartesianGridGeometry<NDIM> > geometry = d_hierarchy->getGridGeometry();
    const IntVector<NDIM>& ratio = d_hierarchy->getPatchLevel(d_level_num)->getRatio();
    const double* const dx_coarsest = geometry->getDx();
    double vol = 1.0;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        vol *= dx_coarsest[d]/static_cast<double>(ratio(d));
    }
    ierr = MatPtAP(*d_J_mat, d_R_mat, MAT_INITIAL_MATRIX, 5.0, &d_RtJR_mat); IBTK_CHKERRQ(ierr);
    ierr = MatDuplicate(d_stokes_mat, MAT_COPY_VALUES, &d_petsc_mat); IBTK_CHKERRQ(ierr);
    ierr = MatAXPY(d_petsc_mat, 1.0/vol, d_RtJR_mat, DIFFERENT_NONZERO_PATTERN); IBTK_CHKERRQ(ierr);
    ierr = MatDuplicate(d_petsc_mat, MAT_COPY_VALUES, &d_petsc_pc); IBTK_CHKERRQ(ierr);
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<HierarchyDataOpsInteger<NDIM> > hier_p_dof_index_ops = hier_ops_manager->getOperationsInteger(d_p_dof_index_var, d_hierarchy, true);
    hier_p_dof_index_ops->resetLevels(d_level_num, d_level_num);
    const int min_p_idx = hier_p_dof_index_ops->min(d_p_dof_index_idx);  // NOTE: HierarchyDataOpsInteger::max() is broken
    ierr = MatZeroRowsColumns(d_petsc_pc, 1, &min_p_idx, 1.0, NULL, NULL); IBTK_CHKERRQ(ierr);
    d_petsc_ksp_ops_flag = DIFFERENT_NONZERO_PATTERN;
    ierr = KSPSetOperators(d_petsc_ksp, d_petsc_mat, d_petsc_pc, d_petsc_ksp_ops_flag); IBTK_CHKERRQ(ierr);
    return;
}// initializeOperator

void
IBImplicitStaggeredPETScLevelSolver::updateOperator()
{
    int ierr;

    // Setup PETSc objects.
    Pointer<CartesianGridGeometry<NDIM> > geometry = d_hierarchy->getGridGeometry();
    const IntVector<NDIM>& ratio = d_hierarchy->getPatchLevel(d_level_num)->getRatio();
    const double* const dx_coarsest = geometry->getDx();
    double vol = 1.0;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        vol *= dx_coarsest[d]/static_cast<double>(ratio(d));
    }
    ierr = MatPtAP(*d_J_mat, d_R_mat, MAT_REUSE_MATRIX, 5.0, &d_RtJR_mat); IBTK_CHKERRQ(ierr);
    ierr = MatZeroEntries(d_petsc_mat); IBTK_CHKERRQ(ierr);
    ierr = MatAXPY(d_petsc_mat, 1.0, d_stokes_mat, SAME_NONZERO_PATTERN); IBTK_CHKERRQ(ierr);
    ierr = MatAXPY(d_petsc_mat, 1.0/vol, d_RtJR_mat, SAME_NONZERO_PATTERN); IBTK_CHKERRQ(ierr);
    ierr = MatCopy(d_petsc_mat, d_petsc_pc, SAME_NONZERO_PATTERN); IBTK_CHKERRQ(ierr);
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<HierarchyDataOpsInteger<NDIM> > hier_p_dof_index_ops = hier_ops_manager->getOperationsInteger(d_p_dof_index_var, d_hierarchy, true);
    hier_p_dof_index_ops->resetLevels(d_level_num, d_level_num);
    const int min_p_idx = hier_p_dof_index_ops->min(d_p_dof_index_idx);  // NOTE: HierarchyDataOpsInteger::max() is broken
    ierr = MatZeroRowsColumns(d_petsc_pc, 1, &min_p_idx, 1.0, NULL, NULL); IBTK_CHKERRQ(ierr);
    d_petsc_ksp_ops_flag = DIFFERENT_NONZERO_PATTERN;
    ierr = KSPSetOperators(d_petsc_ksp, d_petsc_mat, d_petsc_pc, d_petsc_ksp_ops_flag); IBTK_CHKERRQ(ierr);
    return;
}// updateOperator

/////////////////////////////// PROTECTED ////////////////////////////////////

void
IBImplicitStaggeredPETScLevelSolver::initializeSolverStateSpecialized(
    const SAMRAIVectorReal<NDIM,double>& x,
    const SAMRAIVectorReal<NDIM,double>& /*b*/)
{
    // Allocate DOF index data.
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_num);
    if (!level->checkAllocated(d_u_dof_index_idx)) level->allocatePatchData(d_u_dof_index_idx);
    if (!level->checkAllocated(d_p_dof_index_idx)) level->allocatePatchData(d_p_dof_index_idx);

    // Setup PETSc objects.
    int ierr;
    StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(d_num_dofs_per_proc, d_u_dof_index_idx, d_p_dof_index_idx, level);
    const int mpi_rank = SAMRAI_MPI::getRank();
    ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[mpi_rank], PETSC_DETERMINE, &d_petsc_x); IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[mpi_rank], PETSC_DETERMINE, &d_petsc_b); IBTK_CHKERRQ(ierr);
    TBOX_ASSERT(false);  // XXXX
//  StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(d_stokes_mat, &d_problem_coefs, d_u_bc_coefs, d_new_time, d_num_dofs_per_proc, d_u_dof_index_idx, d_p_dof_index_idx, level);
    ierr = MatDuplicate(d_stokes_mat, MAT_COPY_VALUES, &d_petsc_mat); IBTK_CHKERRQ(ierr);
    ierr = MatDuplicate(d_petsc_mat, MAT_COPY_VALUES, &d_petsc_pc); IBTK_CHKERRQ(ierr);
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<HierarchyDataOpsInteger<NDIM> > hier_p_dof_index_ops = hier_ops_manager->getOperationsInteger(d_p_dof_index_var, d_hierarchy, true);
    hier_p_dof_index_ops->resetLevels(d_level_num, d_level_num);
    const int min_p_idx = hier_p_dof_index_ops->min(d_p_dof_index_idx);  // NOTE: HierarchyDataOpsInteger::max() is broken
    ierr = MatZeroRowsColumns(d_petsc_pc, 1, &min_p_idx, 1.0, NULL, NULL); IBTK_CHKERRQ(ierr);
    d_petsc_ksp_ops_flag = DIFFERENT_NONZERO_PATTERN;
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    d_data_synch_sched = StaggeredStokesPETScVecUtilities::constructDataSynchSchedule(u_idx, p_idx, level);
    d_ghost_fill_sched = StaggeredStokesPETScVecUtilities::constructGhostFillSchedule(u_idx, p_idx, level);
    return;
}// initializeSolverStateSpecialized

void
IBImplicitStaggeredPETScLevelSolver::deallocateSolverStateSpecialized()
{
    // Deallocate DOF index data.
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_num);
    if (level->checkAllocated(d_u_dof_index_idx)) level->deallocatePatchData(d_u_dof_index_idx);
    if (level->checkAllocated(d_p_dof_index_idx)) level->deallocatePatchData(d_p_dof_index_idx);
    return;
}// deallocateSolverStateSpecialized

void
IBImplicitStaggeredPETScLevelSolver::copyToPETScVec(
    Vec& petsc_x,
    SAMRAIVectorReal<NDIM,double>& x,
    Pointer<PatchLevel<NDIM> > patch_level)
{
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(petsc_x, u_idx, d_u_dof_index_idx, p_idx, d_p_dof_index_idx, patch_level);
    return;
}// copyToPETScVec

void
IBImplicitStaggeredPETScLevelSolver::copyFromPETScVec(
    Vec& petsc_x,
    SAMRAIVectorReal<NDIM,double>& x,
    Pointer<PatchLevel<NDIM> > patch_level)
{
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(petsc_x, u_idx, d_u_dof_index_idx, p_idx, d_p_dof_index_idx, patch_level, d_data_synch_sched, d_ghost_fill_sched);
    return;
}// copyFromPETScVec

void
IBImplicitStaggeredPETScLevelSolver::setupKSPVecs(
    Vec& petsc_x,
    Vec& petsc_b,
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& b,
    Pointer<PatchLevel<NDIM> > patch_level)
{
    if (!d_initial_guess_nonzero) copyToPETScVec(petsc_x, x, patch_level);
    copyToPETScVec(petsc_b, b, patch_level);
    return;
}// setupKSPVecs

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
