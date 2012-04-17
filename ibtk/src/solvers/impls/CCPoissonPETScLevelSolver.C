// Filename: CCPoissonPETScLevelSolver.C
// Created on 30 Aug 2010 by Boyce Griffith
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

#include "CCPoissonPETScLevelSolver.h"

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
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScVecUtilities.h>
#include <ibtk/PoissonUtilities.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CCPoissonPETScLevelSolver::CCPoissonPETScLevelSolver(
    const std::string& object_name,
    Pointer<Database> input_db)
    : PETScLevelSolver(object_name, input_db),
      d_poisson_spec(d_object_name+"::Poisson specs"),
      d_default_bc_coef(new LocationIndexRobinBcCoefs<NDIM>(d_object_name+"::default_bc_coef", Pointer<Database>(NULL))),
      d_bc_coefs(),
      d_homogeneous_bc(true),
      d_apply_time(0.0),
      d_context(NULL),
      d_dof_index_idx(-1),
      d_dof_index_var(NULL),
      d_data_synch_sched(NULL),
      d_ghost_fill_sched(NULL)
{
    // Initialize the Poisson specifications.
    d_poisson_spec.setCZero();
    d_poisson_spec.setDConstant(-1.0);

    // Setup a default boundary condition object that specifies homogeneous
    // Dirichlet boundary conditions.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_default_bc_coef->setBoundaryValue(2*d  ,0.0);
        d_default_bc_coef->setBoundaryValue(2*d+1,0.0);
    }

    // Initialize the boundary conditions objects.
    setHomogeneousBc(d_homogeneous_bc);
    setPhysicalBcCoef(d_default_bc_coef);

    // Construct the DOF index variable/context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(d_object_name + "::CONTEXT");
    d_dof_index_var = new CellVariable<NDIM,int>(d_object_name + "::dof_index");
    if (var_db->checkVariableExists(d_dof_index_var->getName()))
    {
        d_dof_index_var = var_db->getVariable(d_dof_index_var->getName());
        d_dof_index_idx = var_db->mapVariableAndContextToIndex(d_dof_index_var, d_context);
        var_db->removePatchDataIndex(d_dof_index_idx);
    }
    d_dof_index_idx = var_db->registerVariableAndContext(d_dof_index_var, d_context, CELLG);
    return;
}// CCPoissonPETScLevelSolver

CCPoissonPETScLevelSolver::~CCPoissonPETScLevelSolver()
{
    delete d_default_bc_coef;
    return;
}// ~CCPoissonPETScLevelSolver

void
CCPoissonPETScLevelSolver::setPoissonSpecifications(
    const SAMRAI::solv::PoissonSpecifications& poisson_spec)
{
    d_poisson_spec = poisson_spec;
    return;
}// setPoissonSpecifications

void
CCPoissonPETScLevelSolver::setPhysicalBcCoef(
    RobinBcCoefStrategy<NDIM>* const bc_coef)
{
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(1,bc_coef));
    return;
}// setPhysicalBcCoef

void
CCPoissonPETScLevelSolver::setPhysicalBcCoefs(
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& bc_coefs)
{
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(&bc_coefs[0],&bc_coefs[0]+NDIM));
    return;
}// setPhysicalBcCoefs

void
CCPoissonPETScLevelSolver::setPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
    d_bc_coefs.resize(bc_coefs.size());
    for (unsigned int l = 0; l < bc_coefs.size(); ++l)
    {
        if (bc_coefs[l] != NULL)
        {
            d_bc_coefs[l] = bc_coefs[l];
        }
        else
        {
            d_bc_coefs[l] = d_default_bc_coef;
        }
    }
    return;
}// setPhysicalBcCoefs

void
CCPoissonPETScLevelSolver::setHomogeneousBc(
    bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

void
CCPoissonPETScLevelSolver::setTime(
    const double time)
{
    d_apply_time = time;
    return;
}// setTime

/////////////////////////////// PROTECTED ////////////////////////////////////

void
CCPoissonPETScLevelSolver::initializeSolverStateSpecialized(
    const SAMRAIVectorReal<NDIM,double>& x,
    const SAMRAIVectorReal<NDIM,double>& /*b*/)
{
    // Allocate DOF index data.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int x_idx = x.getComponentDescriptorIndex(0);
    Pointer<CellDataFactory<NDIM,double> > x_fac = var_db->getPatchDescriptor()->getPatchDataFactory(x_idx);
    const int depth = x_fac->getDefaultDepth();
    Pointer<CellDataFactory<NDIM,int> > dof_index_fac = var_db->getPatchDescriptor()->getPatchDataFactory(d_dof_index_idx);
    dof_index_fac->setDefaultDepth(depth);
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_num);
    if (!level->checkAllocated(d_dof_index_idx)) level->allocatePatchData(d_dof_index_idx);

    // Setup PETSc objects.
    int ierr;
    PETScVecUtilities::constructPatchLevelDOFIndices(d_num_dofs_per_proc, d_dof_index_idx, level);
    const int mpi_rank = SAMRAI_MPI::getRank();
    ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[mpi_rank], PETSC_DETERMINE, &d_petsc_x); IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[mpi_rank], PETSC_DETERMINE, &d_petsc_b); IBTK_CHKERRQ(ierr);
    PETScMatUtilities::constructPatchLevelCCLaplaceOp(d_petsc_mat, d_poisson_spec, d_bc_coefs, d_apply_time, d_num_dofs_per_proc, d_dof_index_idx, level);
    d_data_synch_sched = PETScVecUtilities::constructDataSynchSchedule(x_idx, level);
    d_ghost_fill_sched = PETScVecUtilities::constructGhostFillSchedule(x_idx, level);
    return;
}// initializeSolverStateSpecialized

void
CCPoissonPETScLevelSolver::deallocateSolverStateSpecialized()
{
    // Deallocate DOF index data.
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_num);
    if (level->checkAllocated(d_dof_index_idx)) level->deallocatePatchData(d_dof_index_idx);
    return;
}// deallocateSolverStateSpecialized

void
CCPoissonPETScLevelSolver::copyToPETScVecs(
    Vec& petsc_x,
    Vec& petsc_b,
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& b,
    Pointer<PatchLevel<NDIM> > patch_level)
{
    if (!d_initial_guess_nonzero)
    {
        const int x_idx = x.getComponentDescriptorIndex(0);
        PETScVecUtilities::copyToPatchLevelVec(petsc_x, x_idx, d_dof_index_idx, patch_level);
    }
    const int b_idx = b.getComponentDescriptorIndex(0);
    Pointer<CellVariable<NDIM,double> > b_var = b.getComponentVariable(0);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int b_adj_idx = var_db->registerClonedPatchDataIndex(b_var, b_idx);
    patch_level->allocatePatchData(b_adj_idx);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<CellData<NDIM,double> > b_data = patch->getPatchData(b_idx);
        Pointer<CellData<NDIM,double> > b_adj_data = patch->getPatchData(b_adj_idx);
        b_adj_data->copy(*b_data);
        if (!patch->getPatchGeometry()->intersectsPhysicalBoundary()) continue;
        PoissonUtilities::adjustCCBoundaryRhsEntries(patch, *b_adj_data, d_poisson_spec, d_bc_coefs, d_apply_time, d_homogeneous_bc);
    }
    PETScVecUtilities::copyToPatchLevelVec(petsc_b, b_adj_idx, d_dof_index_idx, patch_level);
    patch_level->deallocatePatchData(b_adj_idx);
    var_db->removePatchDataIndex(b_adj_idx);
    return;
}// copyToPETScVecs

void
CCPoissonPETScLevelSolver::copyFromPETScVec(
    Vec& petsc_x,
    SAMRAIVectorReal<NDIM,double>& x,
    Pointer<PatchLevel<NDIM> > patch_level)
{
    const int x_idx = x.getComponentDescriptorIndex(0);
    PETScVecUtilities::copyFromPatchLevelVec(petsc_x, x_idx, d_dof_index_idx, patch_level, d_data_synch_sched, d_ghost_fill_sched);
    return;
}// copyFromPETScVec

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
