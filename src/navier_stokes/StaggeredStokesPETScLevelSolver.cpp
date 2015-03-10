// Filename: StaggeredStokesPETScLevelSolver.cpp
// Created on 08 Sep 2010 by Boyce Griffith
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
#include <ostream>
#include <string>
#include <vector>

#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/math/HierarchyDataOpsInteger.h"
#include "SAMRAI/math/HierarchyDataOpsManager.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "ibamr/StaggeredStokesPETScLevelSolver.h"
#include "ibamr/StaggeredStokesPETScMatUtilities.h"
#include "ibamr/StaggeredStokesPETScVecUtilities.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/GeneralSolver.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/PETScLevelSolver.h"
#include "petscmat.h"
#include "petscsys.h"
#include "petscvec.h"
#include "SAMRAI/tbox/Database.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;
static const int SIDEG = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesPETScLevelSolver::StaggeredStokesPETScLevelSolver(const std::string& object_name,
                                                                 boost::shared_ptr<Database> input_db,
                                                                 const std::string& default_options_prefix)
    : d_context(NULL), d_u_dof_index_idx(-1), d_p_dof_index_idx(-1), d_u_dof_index_var(NULL), d_p_dof_index_var(NULL),
      d_data_synch_sched(NULL), d_ghost_fill_sched(NULL)
{
    GeneralSolver::init(object_name, /*homogeneous_bc*/ false);
    PETScLevelSolver::init(input_db, default_options_prefix);

    // Construct the DOF index variable/context.
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    d_context = var_db->getContext(object_name + "::CONTEXT");
    d_u_dof_index_var = boost::make_shared<SideVariable<int> >(DIM, object_name + "::u_dof_index");
    if (var_db->checkVariableExists(d_u_dof_index_var->getName()))
    {
        d_u_dof_index_var = BOOST_CAST<SideVariable<int> >(var_db->getVariable(d_u_dof_index_var->getName()));
        TBOX_ASSERT(d_u_dof_index_var);
        d_u_dof_index_idx = var_db->mapVariableAndContextToIndex(d_u_dof_index_var, d_context);
        var_db->removePatchDataIndex(d_u_dof_index_idx);
    }
    d_u_dof_index_idx = var_db->registerVariableAndContext(d_u_dof_index_var, d_context, IntVector(DIM, SIDEG));
    d_p_dof_index_var = boost::make_shared<CellVariable<int> >(DIM, object_name + "::p_dof_index");
    if (var_db->checkVariableExists(d_p_dof_index_var->getName()))
    {
        d_p_dof_index_var = BOOST_CAST<CellVariable<int> >(var_db->getVariable(d_p_dof_index_var->getName()));
        TBOX_ASSERT(d_p_dof_index_var);
        d_p_dof_index_idx = var_db->mapVariableAndContextToIndex(d_p_dof_index_var, d_context);
        var_db->removePatchDataIndex(d_p_dof_index_idx);
    }
    d_p_dof_index_idx = var_db->registerVariableAndContext(d_p_dof_index_var, d_context, IntVector(DIM, CELLG));
    return;
}

StaggeredStokesPETScLevelSolver::~StaggeredStokesPETScLevelSolver()
{
    if (d_is_initialized) deallocateSolverState();
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

void StaggeredStokesPETScLevelSolver::initializeSolverStateSpecialized(const SAMRAIVectorReal<double>& x,
                                                                       const SAMRAIVectorReal<double>& /*b*/)
{
    // Allocate DOF index data.
    auto level =d_hierarchy->getPatchLevel(d_level_num);
    if (!level->checkAllocated(d_u_dof_index_idx)) level->allocatePatchData(d_u_dof_index_idx);
    if (!level->checkAllocated(d_p_dof_index_idx)) level->allocatePatchData(d_p_dof_index_idx);

    // Setup PETSc objects.
    int ierr;
    StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
        d_num_dofs_per_proc, d_u_dof_index_idx, d_p_dof_index_idx, level);
    tbox::SAMRAI_MPI comm(PETSC_COMM_WORLD);
    const int mpi_rank = comm.getRank();
    ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[mpi_rank], PETSC_DETERMINE, &d_petsc_x);
    IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[mpi_rank], PETSC_DETERMINE, &d_petsc_b);
    IBTK_CHKERRQ(ierr);
    StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(d_petsc_mat,
                                                                     d_U_problem_coefs,
                                                                     d_U_bc_coefs,
                                                                     d_new_time,
                                                                     d_num_dofs_per_proc,
                                                                     d_u_dof_index_idx,
                                                                     d_p_dof_index_idx,
                                                                     level);
    ierr = MatDuplicate(d_petsc_mat, MAT_COPY_VALUES, &d_petsc_pc);
    IBTK_CHKERRQ(ierr);
    HierarchyDataOpsManager* hier_ops_manager = HierarchyDataOpsManager::getManager();
    boost::shared_ptr<HierarchyDataOpsInteger> hier_p_dof_index_ops =
        hier_ops_manager->getOperationsInteger(d_p_dof_index_var, d_hierarchy, true);
    hier_p_dof_index_ops->resetLevels(d_level_num, d_level_num);
    const int min_p_idx =
        hier_p_dof_index_ops->min(d_p_dof_index_idx); // NOTE: HierarchyDataOpsInteger::max() is broken
    ierr = MatZeroRowsColumns(d_petsc_pc, 1, &min_p_idx, 1.0, NULL, NULL);
    IBTK_CHKERRQ(ierr);
    d_petsc_ksp_ops_flag = SAME_PRECONDITIONER;
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    d_data_synch_sched = StaggeredStokesPETScVecUtilities::constructDataSynchSchedule(u_idx, p_idx, level);
    d_ghost_fill_sched = StaggeredStokesPETScVecUtilities::constructGhostFillSchedule(u_idx, p_idx, level);
    return;
}

void StaggeredStokesPETScLevelSolver::deallocateSolverStateSpecialized()
{
    // Deallocate DOF index data.
    auto level =d_hierarchy->getPatchLevel(d_level_num);
    if (level->checkAllocated(d_u_dof_index_idx)) level->deallocatePatchData(d_u_dof_index_idx);
    if (level->checkAllocated(d_p_dof_index_idx)) level->deallocatePatchData(d_p_dof_index_idx);
    return;
}

void StaggeredStokesPETScLevelSolver::copyToPETScVec(Vec& petsc_x,
                                                     SAMRAIVectorReal<double>& x,
                                                     boost::shared_ptr<PatchLevel> patch_level)
{
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        petsc_x, u_idx, d_u_dof_index_idx, p_idx, d_p_dof_index_idx, patch_level);
    return;
}

void StaggeredStokesPETScLevelSolver::copyFromPETScVec(Vec& petsc_x,
                                                       SAMRAIVectorReal<double>& x,
                                                       boost::shared_ptr<PatchLevel> patch_level)
{
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(petsc_x,
                                                            u_idx,
                                                            d_u_dof_index_idx,
                                                            p_idx,
                                                            d_p_dof_index_idx,
                                                            patch_level,
                                                            d_data_synch_sched,
                                                            d_ghost_fill_sched);
    return;
}

void StaggeredStokesPETScLevelSolver::setupKSPVecs(Vec& petsc_x,
                                                   Vec& petsc_b,
                                                   SAMRAIVectorReal<double>& x,
                                                   SAMRAIVectorReal<double>& b,
                                                   boost::shared_ptr<PatchLevel> patch_level)
{
    if (!d_initial_guess_nonzero) copyToPETScVec(petsc_x, x, patch_level);
    copyToPETScVec(petsc_b, b, patch_level);
    return;
}

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
