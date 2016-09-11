// Filename: PETScLevelSolver.cpp
// Created on 16 Apr 2012 by Boyce Griffith
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

#include <math.h>
#include <stddef.h>
#include <ostream>
#include <string>
#include <vector>

#include "IntVector.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "SAMRAIVectorReal.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/PETScLevelSolver.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "petsc/private/petscimpl.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscsys.h"
#include "petscvec.h"
#include "petscviewerhdf5.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_solve_system;
static Timer* t_initialize_solver_state;
static Timer* t_deallocate_solver_state;

void
generate_petsc_is_from_std_is(std::vector<std::set<int> >& overlap_std,
                              std::vector<std::set<int> >& nonoverlap_std,
                              std::vector<IS>& overlap_petsc,
                              std::vector<IS>& nonoverlap_petsc)
{
    // Destroy old IS'es and generate new ones.
    int ierr;
    for (unsigned int k = 0; k < overlap_petsc.size(); ++k)
    {
        ierr = ISDestroy(&overlap_petsc[k]);
        IBTK_CHKERRQ(ierr);
    }
    overlap_petsc.clear();
    for (unsigned int k = 0; k < nonoverlap_petsc.size(); ++k)
    {
        ierr = ISDestroy(&nonoverlap_petsc[k]);
        IBTK_CHKERRQ(ierr);
    }
    nonoverlap_petsc.clear();

    const int n_overlap_subdomains = static_cast<int>(overlap_std.size());
    overlap_petsc.resize(n_overlap_subdomains);
    for (int k = 0; k < n_overlap_subdomains; ++k)
    {
        PetscInt* overlap_dof_arr;
        const int n_overlap_dofs = static_cast<int>(overlap_std[k].size());
        ierr = PetscMalloc1(n_overlap_dofs, &overlap_dof_arr);
        IBTK_CHKERRQ(ierr);
        std::copy(overlap_std[k].begin(), overlap_std[k].end(), overlap_dof_arr);
        ierr = ISCreateGeneral(PETSC_COMM_SELF, n_overlap_dofs, overlap_dof_arr, PETSC_OWN_POINTER, &overlap_petsc[k]);
        IBTK_CHKERRQ(ierr);
    }

    const int n_nonoverlap_subdomains = static_cast<int>(nonoverlap_std.size());
    nonoverlap_petsc.resize(n_nonoverlap_subdomains);
    for (int k = 0; k < n_nonoverlap_subdomains; ++k)
    {
        PetscInt* nonoverlap_dof_arr;
        const int n_nonoverlap_dofs = static_cast<int>(nonoverlap_std[k].size());
        ierr = PetscMalloc1(n_nonoverlap_dofs, &nonoverlap_dof_arr);
        IBTK_CHKERRQ(ierr);
        std::copy(nonoverlap_std[k].begin(), nonoverlap_std[k].end(), nonoverlap_dof_arr);
        ierr = ISCreateGeneral(
            PETSC_COMM_SELF, n_nonoverlap_dofs, nonoverlap_dof_arr, PETSC_OWN_POINTER, &nonoverlap_petsc[k]);
        IBTK_CHKERRQ(ierr);
    }

    return;
} // generate_petsc_is_from_std_is
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

PETScLevelSolver::PETScLevelSolver()
    : d_hierarchy(),
      d_level_num(-1),
      d_ksp_type(KSPGMRES),
      d_shell_pc_type(""),
      d_options_prefix(""),
      d_petsc_ksp(NULL),
      d_petsc_mat(NULL),
      d_petsc_pc(NULL),
      d_petsc_extern_mat(NULL),
      d_petsc_x(NULL),
      d_petsc_b(NULL)
{
    // Setup default options.
    d_max_iterations = 10000;
    d_abs_residual_tol = 1.0e-50;
    d_rel_residual_tol = 1.0e-5;
    d_ksp_type = KSPGMRES;
    d_pc_type = PCILU;
    d_initial_guess_nonzero = true;
    d_enable_logging = false;
    d_box_size = 2;
    d_overlap_size = 1;

    // Setup Timers.
    IBTK_DO_ONCE(t_solve_system = TimerManager::getManager()->getTimer("IBTK::PETScLevelSolver::solveSystem()");
                 t_initialize_solver_state =
                     TimerManager::getManager()->getTimer("IBTK::PETScLevelSolver::initializeSolverState()");
                 t_deallocate_solver_state =
                     TimerManager::getManager()->getTimer("IBTK::PETScLevelSolver::deallocateSolverState()"););
    return;
} // PETScLevelSolver

PETScLevelSolver::~PETScLevelSolver()
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::~PETScLevelSolver()\n"
                                 << "  subclass must call deallocateSolverState in subclass destructor"
                                 << std::endl);
    }

    int ierr;
    for (size_t i = 0; i < d_nonoverlap_is.size(); ++i)
    {
        ierr = ISDestroy(&d_nonoverlap_is[i]);
        IBTK_CHKERRQ(ierr);
    }
    for (size_t i = 0; i < d_overlap_is.size(); ++i)
    {
        ierr = ISDestroy(&d_overlap_is[i]);
        IBTK_CHKERRQ(ierr);
    }
    return;
} // ~PETScLevelSolver

void
PETScLevelSolver::setKSPType(const std::string& ksp_type)
{
    d_ksp_type = ksp_type;
    return;
} // setKSPType

void
PETScLevelSolver::setOptionsPrefix(const std::string& options_prefix)
{
    d_options_prefix = options_prefix;
    return;
} // setOptionsPrefix

const KSP&
PETScLevelSolver::getPETScKSP() const
{
    return d_petsc_ksp;
} // getPETScKSP

void
PETScLevelSolver::getASMSubdomains(std::vector<IS>** nonoverlapping_subdomains,
                                   std::vector<IS>** overlapping_subdomains)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_is_initialized);
#endif
    *nonoverlapping_subdomains = &d_nonoverlap_is;
    *overlapping_subdomains = &d_overlap_is;
    return;
} // getASMSubdomains

void
PETScLevelSolver::setNullspace(bool contains_constant_vec,
                               const std::vector<Pointer<SAMRAIVectorReal<NDIM, double> > >& nullspace_basis_vecs)
{
    LinearSolver::setNullspace(contains_constant_vec, nullspace_basis_vecs);
    if (d_is_initialized) setupNullspace();
    return;
} // setNullspace

bool
PETScLevelSolver::solveSystem(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& b)
{
    IBTK_TIMER_START(t_solve_system);

    int ierr;

    if (d_enable_logging) plog << d_object_name << "::solveSystem():" << std::endl;

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x, b);

    // Configure solver.
    ierr = KSPSetTolerances(d_petsc_ksp, d_rel_residual_tol, d_abs_residual_tol, PETSC_DEFAULT, d_max_iterations);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(d_petsc_ksp, d_initial_guess_nonzero ? PETSC_TRUE : PETSC_FALSE);
    IBTK_CHKERRQ(ierr);

    // Solve the system.
    setupKSPVecs(d_petsc_x, d_petsc_b, x, b);
    ierr = KSPSolve(d_petsc_ksp, d_petsc_b, d_petsc_x);
    IBTK_CHKERRQ(ierr);
    copyFromPETScVec(d_petsc_x, x);

    // Log solver info.
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(d_petsc_ksp, &reason);
    IBTK_CHKERRQ(ierr);
    const bool converged = reason > 0;
    if (d_enable_logging)
    {
        plog << d_object_name << "::solveSystem(): solver " << (converged ? "converged" : "diverged") << "\n"
             << "iterations = " << d_current_iterations << "\n"
             << "residual norm = " << d_current_residual_norm << std::endl;
    }

    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();

    IBTK_TIMER_STOP(t_solve_system);
    return converged;
} // solveSystem

void
PETScLevelSolver::initializeSolverState(const SAMRAIVectorReal<NDIM, double>& x,
                                        const SAMRAIVectorReal<NDIM, double>& b)
{
    IBTK_TIMER_START(t_initialize_solver_state);

// Rudimentary error checking.
#if !defined(NDEBUG)
    if (x.getNumberOfComponents() != b.getNumberOfComponents())
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  vectors must have the same number of components"
                                 << std::endl);
    }

    const Pointer<PatchHierarchy<NDIM> >& patch_hierarchy = x.getPatchHierarchy();
    if (patch_hierarchy != b.getPatchHierarchy())
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  vectors must have the same hierarchy"
                                 << std::endl);
    }

    const int coarsest_ln = x.getCoarsestLevelNumber();
    if (coarsest_ln < 0)
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  coarsest level number must not be negative"
                                 << std::endl);
    }
    if (coarsest_ln != b.getCoarsestLevelNumber())
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  vectors must have same coarsest level number"
                                 << std::endl);
    }

    const int finest_ln = x.getFinestLevelNumber();
    if (finest_ln < coarsest_ln)
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  finest level number must be >= coarsest level number"
                                 << std::endl);
    }
    if (finest_ln != b.getFinestLevelNumber())
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  vectors must have same finest level number"
                                 << std::endl);
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!patch_hierarchy->getPatchLevel(ln))
        {
            TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                     << "  hierarchy level "
                                     << ln
                                     << " does not exist"
                                     << std::endl);
        }
    }

    if (coarsest_ln != finest_ln)
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  coarsest_ln != finest_ln in PETScLevelSolver"
                                 << std::endl);
    }
#endif
    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized) deallocateSolverState();

    // Get the hierarchy information.
    d_hierarchy = x.getPatchHierarchy();
    d_level_num = x.getCoarsestLevelNumber();
    TBOX_ASSERT(d_level_num == x.getFinestLevelNumber());
    d_level = d_hierarchy->getPatchLevel(d_level_num);
    if (d_level_num > 0)
    {
        d_cf_boundary = new CoarseFineBoundary<NDIM>(*d_hierarchy, d_level_num, IntVector<NDIM>(1));
    }

    // Perform specialized operations to initialize solver state();
    initializeSolverStateSpecialized(x, b);

    // Setup PETSc objects.
    int ierr;
    ierr = KSPCreate(PETSC_COMM_WORLD, &d_petsc_ksp);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetOperators(d_petsc_ksp, d_petsc_mat, d_petsc_pc);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetReusePreconditioner(d_petsc_ksp, PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetType(d_petsc_ksp, d_ksp_type.c_str());
    IBTK_CHKERRQ(ierr);
    PetscBool initial_guess_nonzero = d_initial_guess_nonzero ? PETSC_TRUE : PETSC_FALSE;
    ierr = KSPSetInitialGuessNonzero(d_petsc_ksp, initial_guess_nonzero);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetTolerances(d_petsc_ksp, d_rel_residual_tol, d_abs_residual_tol, PETSC_DEFAULT, d_max_iterations);
    IBTK_CHKERRQ(ierr);

    // Setup KSP PC.
    PC ksp_pc;
    ierr = KSPGetPC(d_petsc_ksp, &ksp_pc);
    IBTK_CHKERRQ(ierr);
    PCType pc_type = d_pc_type.c_str();
    ierr = PCSetType(ksp_pc, pc_type);
    IBTK_CHKERRQ(ierr);
    if (d_options_prefix != "")
    {
        ierr = KSPSetOptionsPrefix(d_petsc_ksp, d_options_prefix.c_str());
        IBTK_CHKERRQ(ierr);
    }
    ierr = KSPSetFromOptions(d_petsc_ksp);
    IBTK_CHKERRQ(ierr);

    // Reset class data structure to correspond to command-line options.
    ierr = KSPGetTolerances(d_petsc_ksp, &d_rel_residual_tol, &d_abs_residual_tol, NULL, &d_max_iterations);
    IBTK_CHKERRQ(ierr);
    ierr = PCGetType(ksp_pc, &pc_type);
    IBTK_CHKERRQ(ierr);
    d_pc_type = pc_type;

    // Set the nullspace.
    if (d_nullspace_contains_constant_vec || !d_nullspace_basis_vecs.empty()) setupNullspace();

    // Setup the preconditioner.
    if (d_pc_type == "asm")
    {
        // Generate user-defined subdomains.
        std::vector<std::set<int> > overlap_is, nonoverlap_is;
        generateASMSubdomains(overlap_is, nonoverlap_is);
        generate_petsc_is_from_std_is(overlap_is, nonoverlap_is, d_overlap_is, d_nonoverlap_is);

        int num_subdomains = static_cast<int>(overlap_is.size());
        if (num_subdomains == 0)
        {
            IS is;
            ierr = ISCreateGeneral(PETSC_COMM_SELF, 0, NULL, PETSC_OWN_POINTER, &is);
            IBTK_CHKERRQ(ierr);
            ierr = PCASMSetLocalSubdomains(ksp_pc, 1, &is, &is);
            IBTK_CHKERRQ(ierr);
            ierr = ISDestroy(&is);
            IBTK_CHKERRQ(ierr);
        }
        else
        {
            ierr = PCASMSetLocalSubdomains(ksp_pc, num_subdomains, &d_overlap_is[0], &d_nonoverlap_is[0]);
            IBTK_CHKERRQ(ierr);
        }
    }

    if (d_pc_type == "fieldsplit")
    {
        std::vector<std::set<int> > field_is;
        std::vector<std::string> field_name;
        generateFieldSplitSubdomains(field_name, field_is);
        d_field_name = field_name;
        const int n_fields = static_cast<int>(field_is.size());

        // Destroy old IS'es and generate new ones.
        for (unsigned int k = 0; k < d_field_is.size(); ++k)
        {
            ierr = ISDestroy(&d_field_is[k]);
            IBTK_CHKERRQ(ierr);
        }
        d_field_is.clear();

        d_field_is.resize(n_fields);
        for (int k = 0; k < n_fields; ++k)
        {
            PetscInt* field_dof_arr;
            const int n_field_dofs = static_cast<int>(field_is[k].size());
            ierr = PetscMalloc1(n_field_dofs, &field_dof_arr);
            IBTK_CHKERRQ(ierr);
            std::copy(field_is[k].begin(), field_is[k].end(), field_dof_arr);
            ierr = ISCreateGeneral(PETSC_COMM_WORLD, n_field_dofs, field_dof_arr, PETSC_OWN_POINTER, &d_field_is[k]);
            IBTK_CHKERRQ(ierr);
            ierr = PCFieldSplitSetIS(ksp_pc, d_field_name[k].c_str(), d_field_is[k]);
            IBTK_CHKERRQ(ierr);
        }
    }

    if (d_pc_type == "shell")
    {
        Mat diagonal_mat_block;
        ierr = MatGetDiagonalBlock(d_petsc_mat, &diagonal_mat_block);
        IBTK_CHKERRQ(ierr);
        ierr = MatCreateVecs(diagonal_mat_block, &d_local_x, &d_local_y);
        IBTK_CHKERRQ(ierr);

        // Generate user-defined subdomains.
        std::vector<std::set<int> > overlap_is, nonoverlap_is;
        generateASMSubdomains(overlap_is, nonoverlap_is);
        d_n_local_subdomains = static_cast<int>(overlap_is.size());
        d_n_subdomains_max = SAMRAI_MPI::maxReduction(d_n_local_subdomains);
        generate_petsc_is_from_std_is(overlap_is, nonoverlap_is, d_overlap_is, d_nonoverlap_is);

        // Get the local submatrices.
        ierr = MatGetSubMatrices(
            d_petsc_mat, d_n_local_subdomains, &d_overlap_is[0], &d_overlap_is[0], MAT_INITIAL_MATRIX, &d_sub_mat);
        IBTK_CHKERRQ(ierr);

        // Setup data for communicating values between local and global representations.
        d_local_overlap_is.resize(d_n_subdomains_max);
        d_local_nonoverlap_is.resize(d_n_subdomains_max);
        d_restriction.resize(d_n_subdomains_max);
        d_prolongation.resize(d_n_subdomains_max);
        d_sub_x.resize(d_n_subdomains_max);
        d_sub_y.resize(d_n_subdomains_max);
#if !defined(NDEBUG)
        std::set<int> idxs;
#endif
        for (int i = 0; i < d_n_subdomains_max; ++i)
        {
            int overlap_is_size = 0, nonoverlap_is_size = 0;
            PetscInt *overlap_indices = NULL, *nonoverlap_indices = NULL;
            if (i < d_n_local_subdomains)
            {
                ierr = ISGetLocalSize(d_overlap_is[i], &overlap_is_size);
                IBTK_CHKERRQ(ierr);
                const int* overlap_is_arr;
                ierr = ISGetIndices(d_overlap_is[i], &overlap_is_arr);

                ierr = ISGetLocalSize(d_nonoverlap_is[i], &nonoverlap_is_size);
                IBTK_CHKERRQ(ierr);
                const int* nonoverlap_is_arr;
                ierr = ISGetIndices(d_nonoverlap_is[i], &nonoverlap_is_arr);

                PetscMalloc(overlap_is_size * sizeof(PetscInt), &overlap_indices);
                PetscMalloc(nonoverlap_is_size * sizeof(PetscInt), &nonoverlap_indices);
                int ii = 0, jj = 0;
                for (; ii < overlap_is_size; ++ii)
                {
                    overlap_indices[ii] = ii;

                    // Keep the local indices of nonoverlap DOFs in an array.
                    // Since we have sorted IS'es, it is easier to locate contigous nonoverlap DOFs.
                    if (jj < nonoverlap_is_size && overlap_is_arr[ii] == nonoverlap_is_arr[jj])
                    {
#if !defined(NDEBUG)
                        TBOX_ASSERT(idxs.find(overlap_is_arr[ii]) == idxs.end());
                        idxs.insert(overlap_is_arr[ii]);
#endif
                        nonoverlap_indices[jj] = ii;
                        ++jj;
                    }
                }
                TBOX_ASSERT(ii == overlap_is_size);
                TBOX_ASSERT(jj == nonoverlap_is_size);

                ierr = ISRestoreIndices(d_overlap_is[i], &overlap_is_arr);
                IBTK_CHKERRQ(ierr);
                ierr = ISRestoreIndices(d_nonoverlap_is[i], &nonoverlap_is_arr);
                IBTK_CHKERRQ(ierr);

                ierr = MatCreateVecs(d_sub_mat[i], &d_sub_x[i], &d_sub_y[i]);
            }
            else
            {
                ierr = VecCreateSeq(PETSC_COMM_SELF, 0, &d_sub_x[i]);
                IBTK_CHKERRQ(ierr);
                ierr = VecCreateSeq(PETSC_COMM_SELF, 0, &d_sub_y[i]);
                IBTK_CHKERRQ(ierr);
            }

            ierr = ISCreateGeneral(
                PETSC_COMM_WORLD, overlap_is_size, overlap_indices, PETSC_OWN_POINTER, &d_local_overlap_is[i]);
            IBTK_CHKERRQ(ierr);
            ierr = ISCreateGeneral(
                PETSC_COMM_WORLD, nonoverlap_is_size, nonoverlap_indices, PETSC_OWN_POINTER, &d_local_nonoverlap_is[i]);
            IBTK_CHKERRQ(ierr);

            IS& overlap_is = (i < d_n_local_subdomains ? d_overlap_is[i] : d_local_overlap_is[i]);
            IS& nonoverlap_is = (i < d_n_local_subdomains ? d_nonoverlap_is[i] : d_local_nonoverlap_is[i]);
            ierr = VecScatterCreate(d_petsc_x, overlap_is, d_sub_x[i], d_local_overlap_is[i], &d_restriction[i]);
            IBTK_CHKERRQ(ierr);
            ierr = VecScatterCreate(d_sub_y[i], d_local_nonoverlap_is[i], d_petsc_b, nonoverlap_is, &d_prolongation[i]);
            IBTK_CHKERRQ(ierr);
        }
#if !defined(NDEBUG)
        int n_local_dofs;
        VecGetSize(d_local_x, &n_local_dofs);
        TBOX_ASSERT(n_local_dofs == static_cast<int>(idxs.size()));
#endif
        if (d_shell_pc_type == "multiplicative")
        {
            PetscInt n_lo, n_hi;
            ierr = VecGetOwnershipRange(d_petsc_x, &n_lo, &n_hi);
            IBTK_CHKERRQ(ierr);
            IS local_idx;
            ierr = ISCreateStride(PETSC_COMM_WORLD, n_hi - n_lo, n_lo, 1, &local_idx);
            IBTK_CHKERRQ(ierr);
            std::vector<IS> local_idxs(d_n_local_subdomains, local_idx);
            ierr = MatGetSubMatrices(d_petsc_mat,
                                     d_n_local_subdomains,
                                     d_n_local_subdomains ? &d_overlap_is[0] : NULL,
                                     d_n_local_subdomains ? &local_idxs[0] : NULL,
                                     MAT_INITIAL_MATRIX,
                                     &d_sub_bc_mat);
            IBTK_CHKERRQ(ierr);
            for (int i = 0; i < d_n_local_subdomains; ++i)
            {
                ierr = MatScale(d_sub_bc_mat[i], -1.0);
                IBTK_CHKERRQ(ierr);
            }
            ierr = ISDestroy(&local_idx);
            IBTK_CHKERRQ(ierr);
        }

        // Set up subdomain KSPs
        d_sub_ksp.resize(d_n_local_subdomains);
        for (int i = 0; i < d_n_local_subdomains; ++i)
        {
            KSP& sub_ksp = d_sub_ksp[i];
            Mat& sub_mat = d_sub_mat[i];
            ierr = KSPCreate(PETSC_COMM_SELF, &sub_ksp);
            IBTK_CHKERRQ(ierr);
            std::string sub_prefix = d_options_prefix + "_sub";
            ierr = KSPSetOptionsPrefix(sub_ksp, sub_prefix.c_str());
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetOperators(sub_ksp, sub_mat, sub_mat);
            IBTK_CHKERRQ(ierr);

            // Set default configuraiton.
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

            // Set from options.
            ierr = KSPSetFromOptions(sub_ksp);
            IBTK_CHKERRQ(ierr);

            // Always use a zero initial guess.
            ierr = KSPSetInitialGuessNonzero(sub_ksp, PETSC_FALSE);
            IBTK_CHKERRQ(ierr);
        }
        ierr = PCSetType(ksp_pc, PCSHELL);
        IBTK_CHKERRQ(ierr);
        ierr = PCShellSetContext(ksp_pc, static_cast<void*>(this));
        IBTK_CHKERRQ(ierr);
        if (d_shell_pc_type == "additive")
        {
            ierr = PCShellSetApply(ksp_pc, PETScLevelSolver::PCApply_Additive);
            IBTK_CHKERRQ(ierr);
            std::string pc_name = d_options_prefix + "PC_Additive";
            ierr = PCShellSetName(ksp_pc, pc_name.c_str());
            IBTK_CHKERRQ(ierr);
        }
        else if (d_shell_pc_type == "multiplicative")
        {
            ierr = PCShellSetApply(ksp_pc, PETScLevelSolver::PCApply_Multiplicative);
            IBTK_CHKERRQ(ierr);
            std::string pc_name = d_options_prefix + "PC_Multiplicative";
            ierr = PCShellSetName(ksp_pc, pc_name.c_str());
            IBTK_CHKERRQ(ierr);
        }
        else
        {
            TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::initializeSolverState()\n"
                                     << "Unknown PCSHELL specified. Supported PCSHELL types are additive and "
                                        "multiplicative."
                                     << std::endl);
        }
    }

    // Indicate that the solver is initialized.
    d_is_initialized = true;

    IBTK_TIMER_STOP(t_initialize_solver_state);
    return;
} // initializeSolverState

void
PETScLevelSolver::deallocateSolverState()
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_deallocate_solver_state);

    // Perform specialized operations to deallocate solver state.
    deallocateSolverStateSpecialized();

    // Deallocate PETSc objects.
    int ierr;
    ierr = KSPDestroy(&d_petsc_ksp);
    IBTK_CHKERRQ(ierr);
    if (d_petsc_mat != d_petsc_pc)
    {
        ierr = MatDestroy(&d_petsc_pc);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatDestroy(&d_petsc_mat);
    IBTK_CHKERRQ(ierr);
    if (d_nullspace_contains_constant_vec || !d_nullspace_basis_vecs.empty())
    {
        ierr = MatNullSpaceDestroy(&d_petsc_nullsp);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecDestroy(&d_petsc_x);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&d_petsc_b);
    IBTK_CHKERRQ(ierr);

    // Deallocate PETSc objects for shell preconditioner.
    if (d_pc_type == "shell")
    {
        for (int i = 0; i < d_n_local_subdomains; ++i)
        {
            ierr = KSPDestroy(&d_sub_ksp[i]);
            IBTK_CHKERRQ(ierr);
        }
        for (int i = 0; i < d_n_subdomains_max; ++i)
        {
            ierr = ISDestroy(&d_local_overlap_is[i]);
            IBTK_CHKERRQ(ierr);
            ierr = ISDestroy(&d_local_nonoverlap_is[i]);
            IBTK_CHKERRQ(ierr);
            ierr = VecScatterDestroy(&d_prolongation[i]);
            IBTK_CHKERRQ(ierr);
            ierr = VecScatterDestroy(&d_restriction[i]);
            IBTK_CHKERRQ(ierr);
        }
        ierr = MatDestroyMatrices(d_n_local_subdomains, &d_sub_mat);
        IBTK_CHKERRQ(ierr);
        if (d_shell_pc_type == "multiplicative" && d_n_local_subdomains > 0)
        {
            ierr = MatDestroyMatrices(d_n_local_subdomains, &d_sub_bc_mat);
            IBTK_CHKERRQ(ierr);
        }
        d_sub_mat = NULL;
        ierr = VecDestroy(&d_local_x);
        IBTK_CHKERRQ(ierr);
        d_local_x = NULL;
        ierr = VecDestroy(&d_local_y);
        IBTK_CHKERRQ(ierr);
        d_local_y = NULL;
        d_n_local_subdomains = 0;
        d_n_subdomains_max = 0;

        d_local_overlap_is.clear();
        d_local_nonoverlap_is.clear();
        d_restriction.clear();
        d_prolongation.clear();
        d_sub_ksp.clear();
        d_sub_x.clear();
        d_sub_x.clear();
    }

    d_petsc_ksp = NULL;
    d_petsc_mat = NULL;
    d_petsc_x = NULL;
    d_petsc_b = NULL;

    // Indicate that the solver is NOT initialized.
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_solver_state);
    return;
} // deallocateSolverState

void
PETScLevelSolver::addLinearOperator(Mat& op)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_is_initialized);
    TBOX_ASSERT(op);
#endif
    d_petsc_extern_mat = op;
    return;
} // addLinearOperator

/////////////////////////////// PROTECTED ////////////////////////////////////

void
PETScLevelSolver::init(Pointer<Database> input_db, const std::string& default_options_prefix)
{
    d_options_prefix = default_options_prefix;
    if (input_db)
    {
        if (input_db->keyExists("options_prefix")) d_options_prefix = input_db->getString("options_prefix");
        if (input_db->keyExists("enable_logging")) d_enable_logging = input_db->getBool("enable_logging");
        if (input_db->keyExists("max_iterations")) d_max_iterations = input_db->getInteger("max_iterations");
        if (input_db->keyExists("abs_residual_tol")) d_abs_residual_tol = input_db->getDouble("abs_residual_tol");
        if (input_db->keyExists("rel_residual_tol")) d_rel_residual_tol = input_db->getDouble("rel_residual_tol");
        if (input_db->keyExists("ksp_type")) d_ksp_type = input_db->getString("ksp_type");
        if (input_db->keyExists("pc_type")) d_pc_type = input_db->getString("pc_type");
        if (input_db->keyExists("shell_pc_type")) d_shell_pc_type = input_db->getString("shell_pc_type");
        if (input_db->keyExists("initial_guess_nonzero"))
            d_initial_guess_nonzero = input_db->getBool("initial_guess_nonzero");
        if (input_db->keyExists("subdomain_box_size"))
            input_db->getIntegerArray("subdomain_box_size", d_box_size, NDIM);
        if (input_db->keyExists("subdomain_overlap_size"))
            input_db->getIntegerArray("subdomain_overlap_size", d_overlap_size, NDIM);
    }
    return;
} // init

void
PETScLevelSolver::generateASMSubdomains(std::vector<std::set<int> >& overlap_is,
                                        std::vector<std::set<int> >& nonoverlap_is)
{
    TBOX_ERROR("PETScLevelSolver::generateASMSubdomains(): Subclasses need to generate ASM subdomains. \n");

    return;
} // generateASMSubdomains

void
PETScLevelSolver::generateFieldSplitSubdomains(std::vector<std::string>& field_names,
                                               std::vector<std::set<int> >& field_is)
{
    TBOX_ERROR(
        "PETScLevelSolver::generateFieldSplitSubdomains(): Subclasses need to generate FieldSplit subdomains. \n");

    return;
} // generateFieldSplitSubdomains

void
PETScLevelSolver::setupNullspace()
{
    int ierr;
    std::vector<Vec> petsc_nullspace_basis_vecs(d_nullspace_basis_vecs.size());
    for (unsigned k = 0; k < d_nullspace_basis_vecs.size(); ++k)
    {
        Vec& petsc_nullspace_vec = petsc_nullspace_basis_vecs[k];
        ierr = MatCreateVecs(d_petsc_mat, NULL, &petsc_nullspace_vec);
        IBTK_CHKERRQ(ierr);
        copyToPETScVec(petsc_nullspace_vec, *d_nullspace_basis_vecs[k]);
        double norm;
        ierr = VecNorm(petsc_nullspace_vec, NORM_2, &norm);
        IBTK_CHKERRQ(ierr);
        ierr = VecScale(petsc_nullspace_vec, 1.0 / norm);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,
                              d_nullspace_contains_constant_vec ? PETSC_TRUE : PETSC_FALSE,
                              static_cast<int>(petsc_nullspace_basis_vecs.size()),
                              (petsc_nullspace_basis_vecs.empty() ? NULL : &petsc_nullspace_basis_vecs[0]),
                              &d_petsc_nullsp);
    IBTK_CHKERRQ(ierr);
    ierr = MatSetNullSpace(d_petsc_mat, d_petsc_nullsp);
    IBTK_CHKERRQ(ierr);
    for (unsigned k = 0; k < d_nullspace_basis_vecs.size(); ++k)
    {
        ierr = VecDestroy(&petsc_nullspace_basis_vecs[k]);
        IBTK_CHKERRQ(ierr);
    }
    return;
} // setupNullspace

/////////////////////////////// PRIVATE //////////////////////////////////////

PetscErrorCode
PETScLevelSolver::PCApply_Additive(PC pc, Vec x, Vec y)
{
    int ierr;
    void* ctx;
    ierr = PCShellGetContext(pc, &ctx);
    IBTK_CHKERRQ(ierr);
    PETScLevelSolver* solver = static_cast<PETScLevelSolver*>(ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(solver);
#endif
    const int n_local_subdomains = solver->d_n_local_subdomains;
    const int n_subdomains_max = solver->d_n_subdomains_max;
    std::vector<VecScatter>& restriction = solver->d_restriction;
    std::vector<VecScatter>& prolongation = solver->d_prolongation;
    std::vector<KSP>& sub_ksp = solver->d_sub_ksp;
    std::vector<Vec>& sub_x = solver->d_sub_x;
    std::vector<Vec>& sub_y = solver->d_sub_y;

    // Restrict the global vector to the local vectors, solve the local systems, and
    // prolong the data back into the global vector.
    for (int i = 0; i < n_subdomains_max; ++i)
    {
        ierr = VecScatterBegin(restriction[i], x, sub_x[i], INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
    }
    for (int i = 0; i < n_subdomains_max; ++i)
    {
        ierr = VecScatterEnd(restriction[i], x, sub_x[i], INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
        if (i < n_local_subdomains)
        {
            ierr = KSPSolve(sub_ksp[i], sub_x[i], sub_y[i]);
            IBTK_CHKERRQ(ierr);
        }
        ierr = VecScatterBegin(prolongation[i], sub_y[i], y, INSERT_VALUES, SCATTER_FORWARD_LOCAL);
        IBTK_CHKERRQ(ierr);
    }
    for (int i = 0; i < n_subdomains_max; ++i)
    {
        ierr = VecScatterEnd(prolongation[i], sub_y[i], y, INSERT_VALUES, SCATTER_FORWARD_LOCAL);
        IBTK_CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} // PCApply_Additive

PetscErrorCode
PETScLevelSolver::PCApply_Multiplicative(PC pc, Vec x, Vec y)
{
    int ierr;
    void* ctx;
    ierr = PCShellGetContext(pc, &ctx);
    IBTK_CHKERRQ(ierr);
    PETScLevelSolver* solver = static_cast<PETScLevelSolver*>(ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(solver);
#endif
    ierr = VecZeroEntries(y);
    IBTK_CHKERRQ(ierr);
    Vec local_y = solver->d_local_y;
    const int n_local_subdomains = solver->d_n_local_subdomains;
    const int n_subdomains_max = solver->d_n_subdomains_max;
    std::vector<VecScatter>& restriction = solver->d_restriction;
    std::vector<VecScatter>& prolongation = solver->d_prolongation;
    std::vector<KSP>& sub_ksp = solver->d_sub_ksp;
    Mat* sub_bc_mat = solver->d_sub_bc_mat;
    std::vector<Vec>& sub_x = solver->d_sub_x;
    std::vector<Vec>& sub_y = solver->d_sub_y;

    // Restrict the global vector to the local vectors, solve the local systems, and
    // prolong the data back into the global vector.
    for (int i = 0; i < n_subdomains_max; ++i)
    {
        ierr = VecScatterBegin(restriction[i], x, sub_x[i], INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
    }
    for (int i = 0; i < n_subdomains_max; ++i)
    {
        ierr = VecScatterEnd(restriction[i], x, sub_x[i], INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
        if (i < n_local_subdomains)
        {
            if (i > 0)
            {
                ierr = VecGetLocalVectorRead(y, local_y);
                IBTK_CHKERRQ(ierr);
                ierr = MatMultAdd(sub_bc_mat[i], local_y, sub_x[i], sub_x[i]);
                IBTK_CHKERRQ(ierr);
                ierr = VecRestoreLocalVectorRead(y, local_y);
                IBTK_CHKERRQ(ierr);
            }
            ierr = KSPSolve(sub_ksp[i], sub_x[i], sub_y[i]);
            IBTK_CHKERRQ(ierr);
        }
        ierr = VecScatterBegin(prolongation[i], sub_y[i], y, INSERT_VALUES, SCATTER_FORWARD_LOCAL);
        IBTK_CHKERRQ(ierr);
        ierr = VecScatterEnd(prolongation[i], sub_y[i], y, INSERT_VALUES, SCATTER_FORWARD_LOCAL);
        IBTK_CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} // PCApply_Multiplicative

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
