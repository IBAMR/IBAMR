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

#include "ibamr/AcousticStreamingPETScMatUtilities.h"
#include "ibamr/AcousticStreamingPETScVecUtilities.h"
#include "ibamr/FOAcousticStreamingPETScLevelSolver.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"

#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/GeneralSolver.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/PETScLevelSolver.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "ibtk/PoissonUtilities.h"

#include "BoundaryBox.h"
#include "CellData.h"
#include "CellVariable.h"
#include "CoarseFineBoundary.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RefineSchedule.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include "petscvec.h"
#include <petsclog.h>

#include <algorithm>
#include <string>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int GCW = 1;

static const int REAL = 0;
static const int IMAG = 1;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

FOAcousticStreamingPETScLevelSolver::~FOAcousticStreamingPETScLevelSolver()
{
    if (d_is_initialized) deallocateSolverState();
    return;
} // ~FOAcousticStreamingPETScLevelSolver

/////////////////////////////// PROTECTED ////////////////////////////////////

FOAcousticStreamingPETScLevelSolver::FOAcousticStreamingPETScLevelSolver(const std::string& object_name,
                                                                         Pointer<Database> input_db,
                                                                         Pointer<Variable<NDIM> > u_dof_idx_var,
                                                                         Pointer<Variable<NDIM> > p_dof_idx_var,
                                                                         const std::string& default_options_prefix)
    : d_u_dof_index_var(u_dof_idx_var), d_p_dof_index_var(p_dof_idx_var)
{
    GeneralSolver::init(object_name, /*homogeneous_bc*/ false);
    PETScLevelSolver::init(input_db, default_options_prefix);

    if (input_db)
    {
        if (input_db->keyExists("acoustic_frequency")) d_omega = 2.0 * M_PI * input_db->getDouble("acoustic_frequency");
        if (input_db->keyExists("frequency")) d_omega = 2.0 * M_PI * input_db->getDouble("frequency");
        if (input_db->keyExists("acoustic_angular_frequency"))
            d_omega = input_db->getDouble("acoustic_angular_frequency");
        if (input_db->keyExists("angular_frequency")) d_omega = input_db->getDouble("angular_frequency");
        if (input_db->keyExists("sound_speed")) d_sound_speed = input_db->getDouble("sound_speed");
    }

    // Construct the DOF index variable/context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(object_name + "::CONTEXT");

    if (var_db->checkVariableExists(d_u_dof_index_var->getName()))
    {
        d_u_dof_index_var = var_db->getVariable(d_u_dof_index_var->getName());
        d_u_dof_index_idx = var_db->mapVariableAndContextToIndex(d_u_dof_index_var, d_context);
        var_db->removePatchDataIndex(d_u_dof_index_idx);
    }
    const int u_gcw = std::max(d_overlap_size.max(), GCW);
    d_u_dof_index_idx = var_db->registerVariableAndContext(d_u_dof_index_var, d_context, u_gcw);

    if (var_db->checkVariableExists(d_p_dof_index_var->getName()))
    {
        d_p_dof_index_var = var_db->getVariable(d_p_dof_index_var->getName());
        d_p_dof_index_idx = var_db->mapVariableAndContextToIndex(d_p_dof_index_var, d_context);
        var_db->removePatchDataIndex(d_p_dof_index_idx);
    }
    const int p_gcw = std::max(d_overlap_size.max(), GCW);
    d_p_dof_index_idx = var_db->registerVariableAndContext(d_p_dof_index_var, d_context, p_gcw);

    return;
} // FOAcousticStreamingPETScLevelSolver

void
FOAcousticStreamingPETScLevelSolver::initializeSolverStateSpecialized(const SAMRAIVectorReal<NDIM, double>& x,
                                                                      const SAMRAIVectorReal<NDIM, double>& /*b*/)
{
    // Allocate DOF index data.
    if (!d_level->checkAllocated(d_u_dof_index_idx)) d_level->allocatePatchData(d_u_dof_index_idx);
    if (!d_level->checkAllocated(d_p_dof_index_idx)) d_level->allocatePatchData(d_p_dof_index_idx);
    AcousticStreamingPETScVecUtilities::constructPatchLevelDOFIndices(
        d_num_dofs_per_proc, d_u_dof_index_idx, d_p_dof_index_idx, d_level);

    // Setup PETSc objects.
    int ierr;
    const int mpi_rank = IBTK_MPI::getRank();
    ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[mpi_rank], PETSC_DETERMINE, &d_petsc_x);
    IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[mpi_rank], PETSC_DETERMINE, &d_petsc_b);
    IBTK_CHKERRQ(ierr);

    AcousticStreamingPETScMatUtilities::constructPatchLevelFOAcousticStreamingOp(d_petsc_mat,
                                                                                 d_omega,
                                                                                 d_sound_speed,
                                                                                 d_rho_idx,
                                                                                 d_mu_idx,
                                                                                 d_lambda_idx,
                                                                                 d_chi_idx,
                                                                                 d_U_bc_coefs,
                                                                                 d_P_bc_coefs,
                                                                                 d_new_time,
                                                                                 d_num_dofs_per_proc,
                                                                                 d_u_dof_index_idx,
                                                                                 d_p_dof_index_idx,
                                                                                 d_level,
                                                                                 d_mu_interp_type);

    d_petsc_pc = d_petsc_mat;

    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    d_data_synch_sched = AcousticStreamingPETScVecUtilities::constructDataSynchSchedule(u_idx, p_idx, d_level);
    d_ghost_fill_sched = AcousticStreamingPETScVecUtilities::constructGhostFillSchedule(u_idx, p_idx, d_level);
    return;
} // initializeSolverStateSpecialized

void
FOAcousticStreamingPETScLevelSolver::deallocateSolverStateSpecialized()
{
    // Deallocate DOF index data.
    if (d_level->checkAllocated(d_u_dof_index_idx)) d_level->deallocatePatchData(d_u_dof_index_idx);
    if (d_level->checkAllocated(d_p_dof_index_idx)) d_level->deallocatePatchData(d_p_dof_index_idx);
    return;
} // deallocateSolverStateSpecialized

void
FOAcousticStreamingPETScLevelSolver::copyToPETScVec(Vec& petsc_x, SAMRAIVectorReal<NDIM, double>& x)
{
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    AcousticStreamingPETScVecUtilities::copyToPatchLevelVec(
        petsc_x, u_idx, d_u_dof_index_idx, p_idx, d_p_dof_index_idx, d_level);
    return;
} // copyToPETScVec

void
FOAcousticStreamingPETScLevelSolver::copyFromPETScVec(Vec& petsc_x, SAMRAIVectorReal<NDIM, double>& x)
{
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    AcousticStreamingPETScVecUtilities::copyFromPatchLevelVec(
        petsc_x, u_idx, d_u_dof_index_idx, p_idx, d_p_dof_index_idx, d_level, d_data_synch_sched, d_ghost_fill_sched);
    return;
} // copyFromPETScVec

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
