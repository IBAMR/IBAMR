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
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

PETScLevelSolver::PETScLevelSolver()
    : d_hierarchy(),
      d_level_num(-1),
      d_use_ksp_as_smoother(false),
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
PETScLevelSolver::getMSMSubdomains(std::vector<IS>** rows_subdomains, std::vector<IS>** cols_subdomains)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_is_initialized);
#endif
    *rows_subdomains = &d_subdomain_row_is;
    *cols_subdomains = &d_subdomain_col_is;

    return;
} // getMSMSubdomains

void
PETScLevelSolver::getMSMSubdomains(std::vector<IS>** red_rows_subdomains,
                                   std::vector<IS>** red_cols_subdomains,
                                   std::vector<IS>** black_rows_subdomains,
                                   std::vector<IS>** black_cols_subdomains)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_is_initialized);
#endif
    *red_rows_subdomains = &d_red_subdomain_row_is;
    *red_cols_subdomains = &d_red_subdomain_col_is;
    *black_rows_subdomains = &d_black_subdomain_row_is;
    *black_cols_subdomains = &d_black_subdomain_col_is;

    return;
} // getMSMSubdomains

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

    // For level smoothers skip computing norms.
    if (d_use_ksp_as_smoother)
    {
        ierr = KSPSetNormType(d_petsc_ksp, KSP_NORM_NONE);
        IBTK_CHKERRQ(ierr);
    }

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
        int num_subdomains = static_cast<int>(d_overlap_is.size());
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
        const int n_fields = static_cast<int>(d_field_is.size());
        for (int k = 0; k < n_fields; ++k)
        {
            ierr = PCFieldSplitSetIS(ksp_pc, d_field_name[k].c_str(), d_field_is[k]);
            IBTK_CHKERRQ(ierr);
        }
    }

    if (d_pc_type == "shell")
    {
        ierr = MatGetDiagonalBlock(d_petsc_mat, &d_diagonal_mat);
        IBTK_CHKERRQ(ierr);

        if (d_shell_pc_type == "additive")
        {
            d_no_subdomains = static_cast<int>(d_subdomain_row_is.size());
            ierr = MatGetSubMatrices(d_petsc_mat,
                                     d_no_subdomains,
                                     &d_subdomain_row_is[0],
                                     &d_subdomain_row_is[0],
                                     MAT_INITIAL_MATRIX,
                                     &d_subdomain_mat);
            IBTK_CHKERRQ(ierr);
            ierr = MatGetSubMatrices(d_petsc_mat,
                                     d_no_subdomains,
                                     &d_subdomain_row_is[0],
                                     &d_subdomain_col_is[0],
                                     MAT_INITIAL_MATRIX,
                                     &d_subdomain_bc_mat);
            IBTK_CHKERRQ(ierr);

            // Set up subdomain KSPs
            d_subdomain_ksp.resize(d_no_subdomains);
            for (int subdomain = 0; subdomain < d_no_subdomains; ++subdomain)
            {
                KSP& sub_ksp = d_subdomain_ksp[subdomain];
                Mat& sub_mat = d_subdomain_mat[subdomain];
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

            ierr = PCSetType(ksp_pc, PCSHELL);
            IBTK_CHKERRQ(ierr);
            ierr = PCShellSetContext(ksp_pc, static_cast<void*>(this));
            IBTK_CHKERRQ(ierr);
            ierr = PCShellSetApply(ksp_pc, PETScLevelSolver::PCApply_Additive);
            IBTK_CHKERRQ(ierr);
            std::string pc_name = d_options_prefix + "PCApply_Additive";
            ierr = PCShellSetName(ksp_pc, pc_name.c_str());
            IBTK_CHKERRQ(ierr);
        }
        else if (d_shell_pc_type == "multiplicative")
        {
            d_no_red_subdomains = static_cast<int>(d_red_subdomain_row_is.size());
            d_no_black_subdomains = static_cast<int>(d_black_subdomain_row_is.size());

            ierr = MatGetSubMatrices(d_petsc_mat,
                                     d_no_red_subdomains,
                                     &d_red_subdomain_row_is[0],
                                     &d_red_subdomain_row_is[0],
                                     MAT_INITIAL_MATRIX,
                                     &d_red_subdomain_mat);
            IBTK_CHKERRQ(ierr);
            ierr = MatGetSubMatrices(d_petsc_mat,
                                     d_no_red_subdomains,
                                     &d_red_subdomain_row_is[0],
                                     &d_red_subdomain_col_is[0],
                                     MAT_INITIAL_MATRIX,
                                     &d_red_subdomain_bc_mat);
            IBTK_CHKERRQ(ierr);
            ierr = MatGetSubMatrices(d_petsc_mat,
                                     d_no_black_subdomains,
                                     &d_black_subdomain_row_is[0],
                                     &d_black_subdomain_row_is[0],
                                     MAT_INITIAL_MATRIX,
                                     &d_black_subdomain_mat);
            IBTK_CHKERRQ(ierr);
            ierr = MatGetSubMatrices(d_petsc_mat,
                                     d_no_black_subdomains,
                                     &d_black_subdomain_row_is[0],
                                     &d_black_subdomain_col_is[0],
                                     MAT_INITIAL_MATRIX,
                                     &d_black_subdomain_bc_mat);
            IBTK_CHKERRQ(ierr);

            // Set up red subdomain KSPs
            d_red_subdomain_ksp.resize(d_no_red_subdomains);
            for (int red_subdomain = 0; red_subdomain < d_no_red_subdomains; ++red_subdomain)
            {
                KSP& sub_ksp = d_red_subdomain_ksp[red_subdomain];
                Mat& sub_mat = d_red_subdomain_mat[red_subdomain];
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

            // Set up black subdomain KSPs
            d_black_subdomain_ksp.resize(d_no_black_subdomains);
            for (int black_subdomain = 0; black_subdomain < d_no_black_subdomains; ++black_subdomain)
            {
                KSP& sub_ksp = d_black_subdomain_ksp[black_subdomain];
                Mat& sub_mat = d_black_subdomain_mat[black_subdomain];
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

            ierr = PCSetType(ksp_pc, PCSHELL);
            IBTK_CHKERRQ(ierr);
            ierr = PCShellSetContext(ksp_pc, static_cast<void*>(this));
            IBTK_CHKERRQ(ierr);
            ierr = PCShellSetApply(ksp_pc, PETScLevelSolver::PCApply_Multiplicative);
            IBTK_CHKERRQ(ierr);
            std::string pc_name = d_options_prefix + "PCApply_Multiplicative";
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
    if (d_pc_type == "shell" && d_shell_pc_type == "additive")
    {
        ierr = MatDestroyMatrices(d_no_subdomains, &d_subdomain_mat);
        IBTK_CHKERRQ(ierr);
        d_subdomain_mat = NULL;
        ierr = MatDestroyMatrices(d_no_subdomains, &d_subdomain_bc_mat);
        IBTK_CHKERRQ(ierr);
        d_subdomain_bc_mat = NULL;

        for (int subdomain = 0; subdomain < d_no_subdomains; ++subdomain)
        {
            KSP& sub_ksp = d_subdomain_ksp[subdomain];
            ierr = KSPDestroy(&sub_ksp);
            IBTK_CHKERRQ(ierr);
            d_subdomain_ksp[subdomain] = NULL;
        }
        d_no_subdomains = std::numeric_limits<int>::min();
    }

    if (d_pc_type == "shell" && d_shell_pc_type == "multiplicative")
    {
        ierr = MatDestroyMatrices(d_no_red_subdomains, &d_red_subdomain_mat);
        IBTK_CHKERRQ(ierr);
        d_red_subdomain_mat = NULL;
        ierr = MatDestroyMatrices(d_no_red_subdomains, &d_red_subdomain_bc_mat);
        IBTK_CHKERRQ(ierr);
        d_red_subdomain_bc_mat = NULL;

        ierr = MatDestroyMatrices(d_no_black_subdomains, &d_black_subdomain_mat);
        IBTK_CHKERRQ(ierr);
        d_black_subdomain_mat = NULL;
        ierr = MatDestroyMatrices(d_no_black_subdomains, &d_black_subdomain_bc_mat);
        IBTK_CHKERRQ(ierr);
        d_black_subdomain_bc_mat = NULL;

        for (int red_subdomain = 0; red_subdomain < d_no_red_subdomains; ++red_subdomain)
        {
            KSP& sub_ksp = d_red_subdomain_ksp[red_subdomain];
            ierr = KSPDestroy(&sub_ksp);
            IBTK_CHKERRQ(ierr);
            d_red_subdomain_ksp[red_subdomain] = NULL;
        }

        for (int black_subdomain = 0; black_subdomain < d_no_black_subdomains; ++black_subdomain)
        {
            KSP& sub_ksp = d_black_subdomain_ksp[black_subdomain];
            ierr = KSPDestroy(&sub_ksp);
            IBTK_CHKERRQ(ierr);
            d_black_subdomain_ksp[black_subdomain] = NULL;
        }
        d_no_red_subdomains = std::numeric_limits<int>::min();
        d_no_black_subdomains = std::numeric_limits<int>::min();
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
        if (input_db->keyExists("use_ksp_as_smoother"))
            d_use_ksp_as_smoother = input_db->getBool("use_ksp_as_smoother");
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
        double dot;
        ierr = VecDot(petsc_nullspace_vec, petsc_nullspace_vec, &dot);
        IBTK_CHKERRQ(ierr);
        ierr = VecScale(petsc_nullspace_vec, 1.0 / sqrt(dot));
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

    // Make sure initial value of y is zero upon entry to preconditioner.
    ierr = VecSet(y, 0.0);
    IBTK_CHKERRQ(ierr);

    // Apply the preconditioner.
    Vec res_vec;
    ierr = VecDuplicate(x, &res_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecCopy(x, res_vec);
    IBTK_CHKERRQ(ierr);

    PetscInt first_local_dof;
    ierr = VecGetOwnershipRange(y, &first_local_dof, NULL);
    IBTK_CHKERRQ(ierr);

    PetscScalar *err_array, *res_array;
    ierr = VecGetArray(y, &err_array);
    IBTK_CHKERRQ(ierr);
    ierr = VecGetArray(res_vec, &res_array);
    IBTK_CHKERRQ(ierr);

    const int& n_subdomains = solver->d_no_subdomains;
    for (int subdomain = 0; subdomain < n_subdomains; ++subdomain)
    {
        Mat& sub_mat = solver->d_subdomain_mat[subdomain];
        KSP& sub_ksp = solver->d_subdomain_ksp[subdomain];

        // Get local (box) DOFs.
        IS& local_dofs = solver->d_subdomain_row_is[subdomain];
        PetscInt is_size_local;
        ierr = ISGetLocalSize(local_dofs, &is_size_local);
        IBTK_CHKERRQ(ierr);
        std::vector<int> local_indices(is_size_local);
        for (int i = 0; i < is_size_local; ++i)
        {
            local_indices[i] = i;
        }
        const PetscInt* is_local_array;
        ierr = ISGetIndices(local_dofs, &is_local_array);
        IBTK_CHKERRQ(ierr);

        // Create little block Vecs for solving local sytem.
        Vec r, e;
        ierr = MatCreateVecs(sub_mat, &e, &r);
        IBTK_CHKERRQ(ierr);

        // Copy to/from various Vecs
        std::vector<double> local_values(is_size_local);
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

        // Destroy temporary PETSc objects
        ierr = VecDestroy(&r);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&e);
        IBTK_CHKERRQ(ierr);
    }

    ierr = VecRestoreArray(y, &err_array);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyBegin(y);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(y);
    IBTK_CHKERRQ(ierr);

    // Destroy temporary objects.
    ierr = VecDestroy(&res_vec);
    IBTK_CHKERRQ(ierr);

    // Reflect change in PETSc state.
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    IBTK_CHKERRQ(ierr);

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

    // Make sure initial value of y is zero upon entry to preconditioner.
    ierr = VecSet(y, 0.0);
    IBTK_CHKERRQ(ierr);

    // Get some temporary Vecs.
    Vec new_res_vec;
    ierr = VecDuplicate(x, &new_res_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecCopy(x, new_res_vec);
    IBTK_CHKERRQ(ierr);

    Vec err_vec_local, new_res_vec_local;
    ierr = MatCreateVecs(solver->d_diagonal_mat, &err_vec_local, &new_res_vec_local);
    IBTK_CHKERRQ(ierr);

    // Extract the local vectors.
    PetscScalar *err_array, *new_res_array;
    ierr = VecGetArray(y, &err_array);
    IBTK_CHKERRQ(ierr);
    ierr = VecGetArray(new_res_vec, &new_res_array);
    IBTK_CHKERRQ(ierr);

    PetscInt first_local_dof;
    ierr = VecGetOwnershipRange(y, &first_local_dof, NULL);
    IBTK_CHKERRQ(ierr);

    // Apply the preconditioner
    // Smooth error on red subdomains.
    const int& n_red_subdomains = solver->d_no_red_subdomains;
    for (int red_subdomain = 0; red_subdomain < n_red_subdomains; ++red_subdomain)
    {
        Mat& bc_mat = solver->d_red_subdomain_bc_mat[red_subdomain];
        KSP& sub_ksp = solver->d_red_subdomain_ksp[red_subdomain];

        // Get local (to box) and non-local (non-box) DOFs.
        // Non-local DOFs are non-subdomains or the non-box DOFs that
        // are still on this processor.
        IS& local_dofs = solver->d_red_subdomain_row_is[red_subdomain];
        IS& nonlocal_dofs = solver->d_red_subdomain_col_is[red_subdomain];

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

    // Assemble error vector for contribution from red subdomains
    // and restore arrays.
    ierr = VecRestoreArray(y, &err_array);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyBegin(y);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(y);
    IBTK_CHKERRQ(ierr);
    ierr = VecRestoreArray(new_res_vec, &new_res_array);
    IBTK_CHKERRQ(ierr);

    // Update residual for black subdomains
    // Do globally: r = f - A*u
    ierr = VecScale(y, -1.0);
    IBTK_CHKERRQ(ierr);
    ierr = MatMultAdd(solver->d_petsc_mat, y, x, new_res_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecScale(y, -1.0);
    IBTK_CHKERRQ(ierr);

    // Do locally: r = r + A*u
    ierr = VecGetLocalVector(y, err_vec_local);
    IBTK_CHKERRQ(ierr);
    ierr = VecGetLocalVector(new_res_vec, new_res_vec_local);
    IBTK_CHKERRQ(ierr);
    ierr = MatMultAdd(solver->d_diagonal_mat, err_vec_local, new_res_vec_local, new_res_vec_local);
    IBTK_CHKERRQ(ierr);
    ierr = VecRestoreLocalVector(y, err_vec_local);
    IBTK_CHKERRQ(ierr);
    ierr = VecRestoreLocalVector(new_res_vec, new_res_vec_local);
    IBTK_CHKERRQ(ierr);

    // Extract the local vectors.
    ierr = VecGetArray(y, &err_array);
    IBTK_CHKERRQ(ierr);
    ierr = VecGetArray(new_res_vec, &new_res_array);
    IBTK_CHKERRQ(ierr);

    // Smooth error on black subdomains.
    const int n_black_subdomains = solver->d_no_black_subdomains;
    for (int black_subdomain = 0; black_subdomain < n_black_subdomains; ++black_subdomain)
    {
        Mat& bc_mat = solver->d_black_subdomain_bc_mat[black_subdomain];
        KSP& sub_ksp = solver->d_black_subdomain_ksp[black_subdomain];

        // Get local (to box) and non-local (non-box) DOFs.
        // Non-local DOFs are non-subdomains or the non-box DOFs that
        // are still on this processor.
        IS& local_dofs = solver->d_black_subdomain_row_is[black_subdomain];
        IS& nonlocal_dofs = solver->d_black_subdomain_col_is[black_subdomain];

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
    ierr = VecRestoreArray(y, &err_array);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyBegin(y);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(y);
    IBTK_CHKERRQ(ierr);
    ierr = VecRestoreArray(new_res_vec, &new_res_array);
    IBTK_CHKERRQ(ierr);

    // Destroy temporary vectors
    ierr = VecDestroy(&new_res_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&err_vec_local);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&new_res_vec_local);
    IBTK_CHKERRQ(ierr);

    // Reflect change in PETSc state.
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    IBTK_CHKERRQ(ierr);

    PetscFunctionReturn(0);
} // PCApply_Multiplicative

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
