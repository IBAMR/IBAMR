// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2021 by the IBAMR developers
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

#include "ibtk/GeneralSolver.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/PETScMatUtilities.h"
#include "ibtk/PETScVecUtilities.h"
#include "ibtk/PoissonUtilities.h"
#include "ibtk/SCPoissonPETScLevelSolver.h"

#include "BoundaryBox.h"
#include "Box.h"
#include "CoarseFineBoundary.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchDescriptor.h"
#include "PatchGeometry.h"
#include "PatchLevel.h"
#include "RefineSchedule.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "SideDataFactory.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Array.h"

#include "petscvec.h"
#include <petsclog.h>

#include <algorithm>
#include <set>
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
// Number of ghosts cells used for each variable quantity.
static const int SIDEG = 1;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

SCPoissonPETScLevelSolver::SCPoissonPETScLevelSolver(const std::string& object_name,
                                                     Pointer<Database> input_db,
                                                     std::string default_options_prefix)
{
    // Configure solver.
    GeneralSolver::init(object_name, /*homogeneous_bc*/ false);
    PETScLevelSolver::init(input_db, std::move(default_options_prefix));

    // Construct the DOF index variable/context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(object_name + "::CONTEXT");
    d_dof_index_var = new SideVariable<NDIM, int>(object_name + "::dof_index");
    if (var_db->checkVariableExists(d_dof_index_var->getName()))
    {
        d_dof_index_var = var_db->getVariable(d_dof_index_var->getName());
        d_dof_index_idx = var_db->mapVariableAndContextToIndex(d_dof_index_var, d_context);
        var_db->removePatchDataIndex(d_dof_index_idx);
    }
    const int gcw = std::max(d_overlap_size.max(), SIDEG);
    d_dof_index_idx = var_db->registerVariableAndContext(d_dof_index_var, d_context, gcw);
    return;
} // SCPoissonPETScLevelSolver

SCPoissonPETScLevelSolver::~SCPoissonPETScLevelSolver()
{
    if (d_is_initialized) deallocateSolverState();
    return;
} // ~SCPoissonPETScLevelSolver

/////////////////////////////// PROTECTED ////////////////////////////////////

void
SCPoissonPETScLevelSolver::generateASMSubdomains(std::vector<std::set<int> >& /*overlap_is*/,
                                                 std::vector<std::set<int> >& /*nonoverlap_is*/)
{
    // Construct subdomains for ASM and MSM preconditioner, indexed directly by PETSc IS.
    PETScMatUtilities::constructPatchLevelASMSubdomains(d_overlap_is,
                                                        d_nonoverlap_is,
                                                        d_box_size,
                                                        d_overlap_size,
                                                        d_num_dofs_per_proc,
                                                        d_dof_index_idx,
                                                        d_level,
                                                        d_cf_boundary);
    return;
} // generateASMSubdomains

void
SCPoissonPETScLevelSolver::initializeSolverStateSpecialized(const SAMRAIVectorReal<NDIM, double>& x,
                                                            const SAMRAIVectorReal<NDIM, double>& /*b*/)
{
    // Allocate DOF index data.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int x_idx = x.getComponentDescriptorIndex(0);
    Pointer<SideDataFactory<NDIM, double> > x_fac = var_db->getPatchDescriptor()->getPatchDataFactory(x_idx);
    const int depth = x_fac->getDefaultDepth();
    Pointer<SideDataFactory<NDIM, int> > dof_index_fac =
        var_db->getPatchDescriptor()->getPatchDataFactory(d_dof_index_idx);
    dof_index_fac->setDefaultDepth(depth);
    if (!d_level->checkAllocated(d_dof_index_idx)) d_level->allocatePatchData(d_dof_index_idx);
    PETScVecUtilities::constructPatchLevelDOFIndices(d_num_dofs_per_proc, d_dof_index_idx, d_level);

    // Setup PETSc objects.
    int ierr;
    const int mpi_rank = IBTK_MPI::getRank();
    ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[mpi_rank], PETSC_DETERMINE, &d_petsc_x);
    IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[mpi_rank], PETSC_DETERMINE, &d_petsc_b);
    IBTK_CHKERRQ(ierr);
    PETScMatUtilities::constructPatchLevelSCLaplaceOp(
        d_petsc_mat, d_poisson_spec, d_bc_coefs, d_solution_time, d_num_dofs_per_proc, d_dof_index_idx, d_level);
    d_petsc_pc = d_petsc_mat;

    // Setup SAMRAI communication objects.
    d_data_synch_sched = PETScVecUtilities::constructDataSynchSchedule(x_idx, d_level);
    d_ghost_fill_sched = PETScVecUtilities::constructGhostFillSchedule(x_idx, d_level);
    return;
} // initializeSolverStateSpecialized

void
SCPoissonPETScLevelSolver::deallocateSolverStateSpecialized()
{
    // Deallocate DOF index data.
    if (d_level->checkAllocated(d_dof_index_idx)) d_level->deallocatePatchData(d_dof_index_idx);
    return;
} // deallocateSolverStateSpecialized

void
SCPoissonPETScLevelSolver::copyToPETScVec(Vec& petsc_x, SAMRAIVectorReal<NDIM, double>& x)
{
    const int x_idx = x.getComponentDescriptorIndex(0);
    PETScVecUtilities::copyToPatchLevelVec(petsc_x, x_idx, d_dof_index_idx, d_level);
    return;
} // copyToPETScVec

void
SCPoissonPETScLevelSolver::copyFromPETScVec(Vec& petsc_x, SAMRAIVectorReal<NDIM, double>& x)
{
    const int x_idx = x.getComponentDescriptorIndex(0);
    PETScVecUtilities::copyFromPatchLevelVec(
        petsc_x, x_idx, d_dof_index_idx, d_level, d_data_synch_sched, d_ghost_fill_sched);
    return;
} // copyFromPETScVec

void
SCPoissonPETScLevelSolver::setupKSPVecs(Vec& petsc_x,
                                        Vec& petsc_b,
                                        SAMRAIVectorReal<NDIM, double>& x,
                                        SAMRAIVectorReal<NDIM, double>& b)
{
    if (d_initial_guess_nonzero) copyToPETScVec(petsc_x, x);
    const bool level_zero = (d_level_num == 0);
    const int x_idx = x.getComponentDescriptorIndex(0);
    const int b_idx = b.getComponentDescriptorIndex(0);
    const auto b_adj_idx = d_cached_eulerian_data.getCachedPatchDataIndex(b_idx);
    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = d_level->getPatch(p());
        Pointer<PatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        Pointer<SideData<NDIM, double> > x_data = patch->getPatchData(x_idx);
        Pointer<SideData<NDIM, double> > b_data = patch->getPatchData(b_idx);
        Pointer<SideData<NDIM, double> > b_adj_data = patch->getPatchData(b_adj_idx);
        b_adj_data->copy(*b_data);
        const bool at_physical_bdry = pgeom->intersectsPhysicalBoundary();
        if (at_physical_bdry)
        {
            PoissonUtilities::adjustRHSAtPhysicalBoundary(
                *b_adj_data, patch, d_poisson_spec, d_bc_coefs, d_solution_time, d_homogeneous_bc);
        }
        const Array<BoundaryBox<NDIM> >& type_1_cf_bdry =
            level_zero ? Array<BoundaryBox<NDIM> >() :
                         d_cf_boundary->getBoundaries(patch->getPatchNumber(), /* boundary type */ 1);
        const bool at_cf_bdry = type_1_cf_bdry.size() > 0;
        if (at_cf_bdry)
        {
            PoissonUtilities::adjustRHSAtCoarseFineBoundary(
                *b_adj_data, *x_data, patch, d_poisson_spec, type_1_cf_bdry);
        }
    }
    PETScVecUtilities::copyToPatchLevelVec(petsc_b, b_adj_idx, d_dof_index_idx, d_level);
    return;
} // setupKSPVecs

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
