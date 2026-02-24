// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2021 by the IBAMR developers
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

#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/PETScMatUtilities.h"
#include "ibtk/PETScVecUtilities.h"
#include "ibtk/PoissonUtilities.h"
#include "ibtk/SCPoissonPETScLevelSolver.h"
#include "ibtk/VCSCViscousPETScLevelSolver.h"
#include "ibtk/ibtk_enums.h"
#include "ibtk/samrai_compatibility_names.h"

#include "MultiblockDataTranslator.h"
#include "SAMRAIArray.h"
#include "SAMRAIBoundaryBox.h"
#include "SAMRAIBox.h"
#include "SAMRAICoarseFineBoundary.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchDescriptor.h"
#include "SAMRAIPatchGeometry.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAIRefineSchedule.h"
#include "SAMRAISAMRAIVectorReal.h"
#include "SAMRAISideData.h"
#include "SAMRAISideDataFactory.h"
#include "SAMRAIVariableDatabase.h"

#include "petscvec.h"
#include <petsclog.h>

#include <string>
#include <utility>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

VCSCViscousPETScLevelSolver::VCSCViscousPETScLevelSolver(std::string object_name,
                                                         SAMRAIPointer<Database> input_db,
                                                         std::string default_options_prefix)
    : SCPoissonPETScLevelSolver(std::move(object_name), input_db, std::move(default_options_prefix))
{
    // Set a default interpolation type.
    d_mu_interp_type = VC_HARMONIC_INTERP;
    return;
} // VCSCViscousPETScLevelSolver

VCSCViscousPETScLevelSolver::~VCSCViscousPETScLevelSolver()
{
    if (d_is_initialized) deallocateSolverState();
    return;
} // ~VCSCViscousPETScLevelSolver

/////////////////////////////// PROTECTED ////////////////////////////////////

void
VCSCViscousPETScLevelSolver::initializeSolverStateSpecialized(const SAMRAISAMRAIVectorReal<double>& x,
                                                              const SAMRAISAMRAIVectorReal<double>& /*b*/)
{
    // Allocate DOF index data.
    SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
    const int x_idx = x.getComponentDescriptorIndex(0);
    SAMRAIPointer<SAMRAISideDataFactory<double>> x_fac = var_db->getPatchDescriptor()->getPatchDataFactory(x_idx);
    const int depth = x_fac->getDefaultDepth();
    SAMRAIPointer<SAMRAISideDataFactory<int>> dof_index_fac =
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
    const double alpha = 1.0;
    const double beta = 1.0;
    PETScMatUtilities::constructPatchLevelVCSCViscousOp(d_petsc_mat,
                                                        d_poisson_spec,
                                                        alpha,
                                                        beta,
                                                        d_bc_coefs,
                                                        d_solution_time,
                                                        d_num_dofs_per_proc,
                                                        d_dof_index_idx,
                                                        d_level,
                                                        d_mu_interp_type);

    d_petsc_pc = d_petsc_mat;

    // Setup SAMRAI communication objects.
    d_data_synch_sched = PETScVecUtilities::constructDataSynchSchedule(x_idx, d_level);
    d_ghost_fill_sched = PETScVecUtilities::constructGhostFillSchedule(x_idx, d_level);
    return;
} // initializeSolverStateSpecialized

void
VCSCViscousPETScLevelSolver::setupKSPVecs(Vec& petsc_x,
                                          Vec& petsc_b,
                                          SAMRAISAMRAIVectorReal<double>& x,
                                          SAMRAISAMRAIVectorReal<double>& b)
{
    if (d_initial_guess_nonzero) copyToPETScVec(petsc_x, x);
    const bool level_zero = (d_level_num == 0);
    const int x_idx = x.getComponentDescriptorIndex(0);
    const int b_idx = b.getComponentDescriptorIndex(0);
    const auto b_adj_idx = d_cached_eulerian_data.getCachedPatchDataIndex(b_idx);
    for (SAMRAIPatchLevel::Iterator p(d_level); p; p++)
    {
        SAMRAIPointer<SAMRAIPatch> patch = d_level->getPatch(p());
        SAMRAIPointer<SAMRAIPatchGeometry> pgeom = patch->getPatchGeometry();
        SAMRAIPointer<SAMRAISideData<double>> x_data = patch->getPatchData(x_idx);
        SAMRAIPointer<SAMRAISideData<double>> b_data = patch->getPatchData(b_idx);
        SAMRAIPointer<SAMRAISideData<double>> b_adj_data = patch->getPatchData(b_adj_idx);
        b_adj_data->copy(*b_data);
        const bool at_physical_bdry = pgeom->intersectsPhysicalBoundary();
        if (at_physical_bdry)
        {
            PoissonUtilities::adjustVCSCViscousOpRHSAtPhysicalBoundary(*b_adj_data,
                                                                       patch,
                                                                       d_poisson_spec,
                                                                       1.0,
                                                                       d_bc_coefs,
                                                                       d_solution_time,
                                                                       d_homogeneous_bc,
                                                                       d_mu_interp_type);
        }
        const SAMRAIArray<SAMRAIBoundaryBox>& type_1_cf_bdry =
            level_zero ? SAMRAIArray<SAMRAIBoundaryBox>() :
                         d_cf_boundary->getBoundaries(patch->getPatchNumber(), /* boundary type */ 1, d_mu_interp_type);
        const bool at_cf_bdry = type_1_cf_bdry.size() > 0;
        if (at_cf_bdry)
        {
            PoissonUtilities::adjustVCSCViscousOpRHSAtCoarseFineBoundary(
                *b_adj_data, *x_data, patch, d_poisson_spec, 1.0, type_1_cf_bdry);
        }
    }
    PETScVecUtilities::copyToPatchLevelVec(petsc_b, b_adj_idx, d_dof_index_idx, d_level);
    return;
} // setupKSPVecs

void
VCSCViscousPETScLevelSolver::setViscosityInterpolationType(const IBTK::VCInterpType mu_interp_type)
{
    d_mu_interp_type = mu_interp_type;
    return;
} // setViscosityInterpolationType

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
