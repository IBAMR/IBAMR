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

#include "CellVariable.h"
#include "HierarchyDataOpsInteger.h"
#include "HierarchyDataOpsManager.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RefineSchedule.h"
#include "SAMRAIVectorReal.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "ibamr/StaggeredStokesPETScLevelSolver.h"
#include "ibamr/StaggeredStokesPETScMatUtilities.h"
#include "ibamr/StaggeredStokesPETScVecUtilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/GeneralSolver.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/PETScLevelSolver.h"
#include "ibtk/PoissonUtilities.h"
#include "petscmat.h"
#include "petscsys.h"
#include "petscvec.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/SAMRAI_MPI.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;
static const int SIDEG = 1;
static const int NOGHOST = 0;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesPETScLevelSolver::StaggeredStokesPETScLevelSolver(const std::string& object_name,
                                                                 Pointer<Database> input_db,
                                                                 const std::string& default_options_prefix)
    : d_context(NULL),
      d_u_dof_index_idx(-1),
      d_p_dof_index_idx(-1),
      d_u_nullspace_idx(-1),
      d_p_nullspace_idx(-1),
      d_u_dof_index_var(NULL),
      d_u_nullspace_var(NULL),
      d_p_dof_index_var(NULL),
      d_p_nullspace_var(NULL),
      d_data_synch_sched(NULL),
      d_ghost_fill_sched(NULL)
{
    GeneralSolver::init(object_name, /*homogeneous_bc*/ false);
    PETScLevelSolver::init(input_db, default_options_prefix);

    // Construct the DOF index variable/context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(object_name + "::CONTEXT");
    d_u_dof_index_var = new SideVariable<NDIM, int>(object_name + "::u_dof_index");
    if (var_db->checkVariableExists(d_u_dof_index_var->getName()))
    {
        d_u_dof_index_var = var_db->getVariable(d_u_dof_index_var->getName());
        d_u_dof_index_idx = var_db->mapVariableAndContextToIndex(d_u_dof_index_var, d_context);
        var_db->removePatchDataIndex(d_u_dof_index_idx);
    }
    const int u_gcw = std::max(d_overlap_size.max(), SIDEG);
    d_u_dof_index_idx = var_db->registerVariableAndContext(d_u_dof_index_var, d_context, u_gcw);
    d_p_dof_index_var = new CellVariable<NDIM, int>(object_name + "::p_dof_index");
    if (var_db->checkVariableExists(d_p_dof_index_var->getName()))
    {
        d_p_dof_index_var = var_db->getVariable(d_p_dof_index_var->getName());
        d_p_dof_index_idx = var_db->mapVariableAndContextToIndex(d_p_dof_index_var, d_context);
        var_db->removePatchDataIndex(d_p_dof_index_idx);
    }
    const int p_gcw = std::max(d_overlap_size.max(), CELLG);
    d_p_dof_index_idx = var_db->registerVariableAndContext(d_p_dof_index_var, d_context, p_gcw);

    // Construct the nullspace variable/index.
    d_u_nullspace_var = new SideVariable<NDIM, double>(object_name + "::u_nullspace_var");
    if (var_db->checkVariableExists(d_u_nullspace_var->getName()))
    {
        d_u_nullspace_var = var_db->getVariable(d_u_nullspace_var->getName());
        d_u_nullspace_idx = var_db->mapVariableAndContextToIndex(d_u_nullspace_var, d_context);
        var_db->removePatchDataIndex(d_u_nullspace_idx);
    }
    d_u_nullspace_idx = var_db->registerVariableAndContext(d_u_nullspace_var, d_context, NOGHOST);
    d_p_nullspace_var = new CellVariable<NDIM, double>(object_name + "::p_nullspace_var");
    if (var_db->checkVariableExists(d_p_nullspace_var->getName()))
    {
        d_p_nullspace_var = var_db->getVariable(d_p_nullspace_var->getName());
        d_p_nullspace_idx = var_db->mapVariableAndContextToIndex(d_p_nullspace_var, d_context);
        var_db->removePatchDataIndex(d_p_nullspace_idx);
    }
    d_p_nullspace_idx = var_db->registerVariableAndContext(d_p_nullspace_var, d_context, NOGHOST);

    return;
} // StaggeredStokesPETScLevelSolver

StaggeredStokesPETScLevelSolver::~StaggeredStokesPETScLevelSolver()
{
    if (d_is_initialized) deallocateSolverState();
    return;
} // ~StaggeredStokesPETScLevelSolver

/////////////////////////////// PROTECTED ////////////////////////////////////

void
StaggeredStokesPETScLevelSolver::initializeSolverStateSpecialized(const SAMRAIVectorReal<NDIM, double>& x,
                                                                  const SAMRAIVectorReal<NDIM, double>& /*b*/)
{
    // Allocate DOF index data.
    if (!d_level->checkAllocated(d_u_dof_index_idx)) d_level->allocatePatchData(d_u_dof_index_idx);
    if (!d_level->checkAllocated(d_p_dof_index_idx)) d_level->allocatePatchData(d_p_dof_index_idx);
    StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
        d_num_dofs_per_proc, d_u_dof_index_idx, d_p_dof_index_idx, d_level);

    // Setup PETSc objects.
    int ierr;
    const int mpi_rank = SAMRAI_MPI::getRank();
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
                                                                     d_level);
    if (d_petsc_extern_mat)
    {
        ierr = MatAXPY(d_petsc_mat, 1.0, d_petsc_extern_mat, DIFFERENT_NONZERO_PATTERN);
        IBTK_CHKERRQ(ierr);
    }
    d_petsc_pc = d_petsc_mat;

    // Construct subdomains for ASM and MSM preconditioner.
    StaggeredStokesPETScMatUtilities::constructPatchLevelASMSubdomains(d_overlap_is,
                                                                       d_nonoverlap_is,
                                                                       d_box_size,
                                                                       d_overlap_size,
                                                                       d_num_dofs_per_proc,
                                                                       d_u_dof_index_idx,
                                                                       d_p_dof_index_idx,
                                                                       d_level,
                                                                       d_cf_boundary);
    StaggeredStokesPETScMatUtilities::constructPatchLevelMSMSubdomains(d_subdomain_row_is,
                                                                       d_subdomain_col_is,
                                                                       d_box_size,
                                                                       d_overlap_size,
                                                                       d_num_dofs_per_proc,
                                                                       d_u_dof_index_idx,
                                                                       d_p_dof_index_idx,
                                                                       d_level);
    StaggeredStokesPETScMatUtilities::constructPatchLevelMSMSubdomains(d_red_subdomain_row_is,
                                                                       d_red_subdomain_col_is,
                                                                       d_black_subdomain_row_is,
                                                                       d_black_subdomain_col_is,
                                                                       d_box_size,
                                                                       d_overlap_size,
                                                                       d_num_dofs_per_proc,
                                                                       d_u_dof_index_idx,
                                                                       d_p_dof_index_idx,
                                                                       d_level);

    // Set IS'es for field split preconditioner.
    StaggeredStokesPETScMatUtilities::constructPatchLevelFields(
        d_field_is, d_field_name, d_num_dofs_per_proc, d_u_dof_index_idx, d_p_dof_index_idx, d_level);

    // Set pressure nullspace if the level covers the entire domain.
    if (d_has_pressure_nullspace)
    {
        bool level_covers_entire_domain = d_level_num == 0;
        if (d_level_num > 0)
        {
            int local_cf_bdry_box_size = 0;
            for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = d_level->getPatch(p());
                const Array<BoundaryBox<NDIM> >& type_1_cf_bdry =
                    d_cf_boundary->getBoundaries(patch->getPatchNumber(), /* boundary type */ 1);
                local_cf_bdry_box_size += type_1_cf_bdry.size();
            }
            level_covers_entire_domain = SAMRAI_MPI::sumReduction(local_cf_bdry_box_size) == 0;
        }

        if (level_covers_entire_domain)
        {
            // Allocate pressure nullspace data.
            if (!d_level->checkAllocated(d_u_nullspace_idx)) d_level->allocatePatchData(d_u_nullspace_idx);
            if (!d_level->checkAllocated(d_p_nullspace_idx)) d_level->allocatePatchData(d_p_nullspace_idx);

            Pointer<SAMRAIVectorReal<NDIM, double> > nullspace_vec = new SAMRAIVectorReal<NDIM, double>(
                d_object_name + "nullspace_vec", d_hierarchy, d_level_num, d_level_num);
            nullspace_vec->addComponent(d_u_nullspace_var, d_u_nullspace_idx);
            nullspace_vec->addComponent(d_p_nullspace_var, d_p_nullspace_idx);
            for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = d_level->getPatch(p());
                Pointer<SideData<NDIM, double> > u_patch_data = nullspace_vec->getComponentPatchData(0, *patch);
                u_patch_data->fill(0.0);
                Pointer<CellData<NDIM, double> > p_patch_data = nullspace_vec->getComponentPatchData(1, *patch);
                p_patch_data->fill(1.0);
            }

            LinearSolver::setNullspace(/*const vec*/ false,
                                       std::vector<Pointer<SAMRAIVectorReal<NDIM, double> > >(1, nullspace_vec));
        }
    }

    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    d_data_synch_sched = StaggeredStokesPETScVecUtilities::constructDataSynchSchedule(u_idx, p_idx, d_level);
    d_ghost_fill_sched = StaggeredStokesPETScVecUtilities::constructGhostFillSchedule(u_idx, p_idx, d_level);
    return;
} // initializeSolverStateSpecialized

void
StaggeredStokesPETScLevelSolver::deallocateSolverStateSpecialized()
{
    // Deallocate DOF index data.
    if (d_level->checkAllocated(d_u_dof_index_idx)) d_level->deallocatePatchData(d_u_dof_index_idx);
    if (d_level->checkAllocated(d_p_dof_index_idx)) d_level->deallocatePatchData(d_p_dof_index_idx);
    return;
} // deallocateSolverStateSpecialized

void
StaggeredStokesPETScLevelSolver::copyToPETScVec(Vec& petsc_x, SAMRAIVectorReal<NDIM, double>& x)
{
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        petsc_x, u_idx, d_u_dof_index_idx, p_idx, d_p_dof_index_idx, d_level);
    return;
} // copyToPETScVec

void
StaggeredStokesPETScLevelSolver::copyFromPETScVec(Vec& petsc_x, SAMRAIVectorReal<NDIM, double>& x)
{
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(
        petsc_x, u_idx, d_u_dof_index_idx, p_idx, d_p_dof_index_idx, d_level, d_data_synch_sched, d_ghost_fill_sched);
    return;
} // copyFromPETScVec

void
StaggeredStokesPETScLevelSolver::setupKSPVecs(Vec& petsc_x,
                                              Vec& petsc_b,
                                              SAMRAIVectorReal<NDIM, double>& x,
                                              SAMRAIVectorReal<NDIM, double>& b)
{
    if (d_initial_guess_nonzero) copyToPETScVec(petsc_x, x);
    const bool level_zero = (d_level_num == 0);
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    const int f_idx = b.getComponentDescriptorIndex(0);
    const int h_idx = b.getComponentDescriptorIndex(1);
    Pointer<SideVariable<NDIM, double> > f_var = b.getComponentVariable(0);
    Pointer<CellVariable<NDIM, double> > h_var = b.getComponentVariable(1);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int f_adj_idx = var_db->registerClonedPatchDataIndex(f_var, f_idx);
    int h_adj_idx = var_db->registerClonedPatchDataIndex(h_var, h_idx);
    d_level->allocatePatchData(f_adj_idx);
    d_level->allocatePatchData(h_adj_idx);
    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = d_level->getPatch(p());
        Pointer<PatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_idx);
        Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(f_idx);
        Pointer<CellData<NDIM, double> > h_data = patch->getPatchData(h_idx);
        Pointer<SideData<NDIM, double> > f_adj_data = patch->getPatchData(f_adj_idx);
        Pointer<CellData<NDIM, double> > h_adj_data = patch->getPatchData(h_adj_idx);
        f_adj_data->copy(*f_data);
        h_adj_data->copy(*h_data);
        const bool at_physical_bdry = pgeom->intersectsPhysicalBoundary();
        // TODO: should we be using target data idx's here?
        StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
            d_U_bc_coefs, d_P_bc_coef, u_idx, p_idx, d_homogeneous_bc);
        if (at_physical_bdry)
        {
            PoissonUtilities::adjustRHSAtPhysicalBoundary(
                *f_adj_data, patch, d_U_problem_coefs, d_U_bc_coefs, d_solution_time, d_homogeneous_bc);
            d_bc_helper->enforceNormalVelocityBoundaryConditions(
                f_adj_idx, h_adj_idx, d_U_bc_coefs, d_solution_time, d_homogeneous_bc, d_level_num, d_level_num);
        }
        const Array<BoundaryBox<NDIM> >& type_1_cf_bdry =
            level_zero ? Array<BoundaryBox<NDIM> >() :
                         d_cf_boundary->getBoundaries(patch->getPatchNumber(), /* boundary type */ 1);
        const bool at_cf_bdry = type_1_cf_bdry.size() > 0;
        if (at_cf_bdry)
        {
            PoissonUtilities::adjustRHSAtCoarseFineBoundary(
                *f_adj_data, *u_data, patch, d_U_problem_coefs, type_1_cf_bdry);
        }
    }

    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        petsc_b, f_adj_idx, d_u_dof_index_idx, h_adj_idx, d_p_dof_index_idx, d_level);

    d_level->deallocatePatchData(f_adj_idx);
    d_level->deallocatePatchData(h_adj_idx);
    var_db->removePatchDataIndex(f_adj_idx);
    var_db->removePatchDataIndex(h_adj_idx);

    copyToPETScVec(petsc_b, b);
    return;
} // setupKSPVecs

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
