// Filename: PETScLevelSolver.C
// Created on 16 Apr 2012 by Boyce Griffith
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

#include "PETScLevelSolver.h"

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

// SAMRAI INCLUDES
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

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

// Number of ghosts cells used for each variable quantity.
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

PETScLevelSolver::PETScLevelSolver(
    const std::string& object_name,
    Pointer<Database> input_db)
    : d_object_name(object_name),
      d_is_initialized(false),
      d_hierarchy(),
      d_level_num(-1),
      d_options_prefix(""),
      d_petsc_ksp(PETSC_NULL),
      d_petsc_mat(PETSC_NULL),
      d_petsc_x(PETSC_NULL),
      d_petsc_b(PETSC_NULL),
      d_max_iterations(10),
      d_abs_residual_tol(0.0),
      d_rel_residual_tol(1.0e-6),
      d_initial_guess_nonzero(false),
      d_current_its(-1),
      d_current_residual_norm(-1.0),
      d_enable_logging(false)
{
    // Get values from the input database.
    if (!input_db.isNull())
    {
        d_options_prefix = input_db->getStringWithDefault("options_prefix", d_options_prefix);
        d_max_iterations = input_db->getIntegerWithDefault("max_iterations", d_max_iterations);
        d_abs_residual_tol = input_db->getDoubleWithDefault("absolute_residual_tol", d_abs_residual_tol);
        d_rel_residual_tol = input_db->getDoubleWithDefault("relative_residual_tol", d_rel_residual_tol);
        d_initial_guess_nonzero = input_db->getBoolWithDefault("initial_guess_nonzero", d_initial_guess_nonzero);
        d_enable_logging = input_db->getBoolWithDefault("enable_logging", d_enable_logging);
    }

    // Setup Timers.
    IBTK_DO_ONCE(
        t_solve_system            = TimerManager::getManager()->getTimer("IBTK::PETScLevelSolver::solveSystem()");
        t_initialize_solver_state = TimerManager::getManager()->getTimer("IBTK::PETScLevelSolver::initializeSolverState()");
        t_deallocate_solver_state = TimerManager::getManager()->getTimer("IBTK::PETScLevelSolver::deallocateSolverState()");
                 );
    return;
}// PETScLevelSolver

PETScLevelSolver::~PETScLevelSolver()
{
    if (d_is_initialized) deallocateSolverState();
    return;
}// ~PETScLevelSolver

bool
PETScLevelSolver::solveSystem(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& b)
{
    IBTK_TIMER_START(t_solve_system);

    int ierr;

    if (d_enable_logging) plog << d_object_name << "::solveSystem():" << std::endl;

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x,b);

    // Configure solver.
    ierr = KSPSetTolerances(d_petsc_ksp, d_rel_residual_tol, d_abs_residual_tol, PETSC_DEFAULT, d_max_iterations); IBTK_CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(d_petsc_ksp, d_initial_guess_nonzero ? PETSC_TRUE : PETSC_FALSE); IBTK_CHKERRQ(ierr);

    // Solve the system.
    Pointer<PatchLevel<NDIM> > patch_level = d_hierarchy->getPatchLevel(d_level_num);
    copyToPETScVecs(d_petsc_x, d_petsc_b, x, b, patch_level);
    ierr = KSPSolve(d_petsc_ksp, d_petsc_b, d_petsc_x); IBTK_CHKERRQ(ierr);
    copyFromPETScVec(d_petsc_x, x, patch_level);

    // Log solver info.
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(d_petsc_ksp, &reason); IBTK_CHKERRQ(ierr);
    const bool converged = reason > 0;
    if (d_enable_logging)
    {
        plog << d_object_name << "::solveSystem(): solver " << (converged ? "converged" : "diverged") << "\n"
             << "iterations = " << d_current_its << "\n"
             << "residual norm = " << d_current_residual_norm << std::endl;
    }

    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();

    IBTK_TIMER_STOP(t_solve_system);
    return converged;
}// solveSystem

void
PETScLevelSolver::initializeSolverState(
    const SAMRAIVectorReal<NDIM,double>& x,
    const SAMRAIVectorReal<NDIM,double>& b)
{
    IBTK_TIMER_START(t_initialize_solver_state);

    // Rudimentary error checking.
#ifdef DEBUG_CHECK_ASSERTIONS
    if (x.getNumberOfComponents() != b.getNumberOfComponents())
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                   << "  vectors must have the same number of components" << std::endl);
    }

    const Pointer<PatchHierarchy<NDIM> >& patch_hierarchy = x.getPatchHierarchy();
    if (patch_hierarchy != b.getPatchHierarchy())
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                   << "  vectors must have the same hierarchy" << std::endl);
    }

    const int coarsest_ln = x.getCoarsestLevelNumber();
    if (coarsest_ln < 0)
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                   << "  coarsest level number must not be negative" << std::endl);
    }
    if (coarsest_ln != b.getCoarsestLevelNumber())
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                   << "  vectors must have same coarsest level number" << std::endl);
    }

    const int finest_ln = x.getFinestLevelNumber();
    if (finest_ln < coarsest_ln)
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                   << "  finest level number must be >= coarsest level number" << std::endl);
    }
    if (finest_ln != b.getFinestLevelNumber())
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                   << "  vectors must have same finest level number" << std::endl);
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (patch_hierarchy->getPatchLevel(ln).isNull())
        {
            TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                       << "  hierarchy level " << ln << " does not exist" << std::endl);
        }
    }

    if (coarsest_ln != finest_ln)
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                   << "  coarsest_ln != finest_ln in PETScLevelSolver" << std::endl);
    }
#endif
    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized) deallocateSolverState();

    // Get the hierarchy information.
    d_hierarchy = x.getPatchHierarchy();
    d_level_num = x.getCoarsestLevelNumber();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_level_num == x.getFinestLevelNumber());
#endif

    // Perform specialized operations to initialize solver state();
    initializeSolverStateSpecialized(x,b);

    // Setup PETSc objects.
    int ierr;
    ierr = KSPCreate(PETSC_COMM_WORLD, &d_petsc_ksp); IBTK_CHKERRQ(ierr);
    ierr = KSPSetOperators(d_petsc_ksp, d_petsc_mat, d_petsc_mat, SAME_PRECONDITIONER); IBTK_CHKERRQ(ierr);
    if (!d_options_prefix.empty())
    {
        ierr = KSPSetOptionsPrefix(d_petsc_ksp, d_options_prefix.c_str()); IBTK_CHKERRQ(ierr);
    }
    ierr = KSPSetFromOptions(d_petsc_ksp); IBTK_CHKERRQ(ierr);

    // Indicate that the solver is initialized.
    d_is_initialized = true;

    IBTK_TIMER_STOP(t_initialize_solver_state);
    return;
}// initializeSolverState

void
PETScLevelSolver::deallocateSolverState()
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_deallocate_solver_state);

    // Perform specialized operations to deallocate solver state.
    deallocateSolverStateSpecialized();

    // Deallocate PETSc objects.
    int ierr;
    ierr = KSPDestroy(&d_petsc_ksp); IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&d_petsc_mat); IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&d_petsc_x); IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&d_petsc_b); IBTK_CHKERRQ(ierr);

    d_petsc_ksp = PETSC_NULL;
    d_petsc_mat = PETSC_NULL;
    d_petsc_x = PETSC_NULL;
    d_petsc_b = PETSC_NULL;

    // Indicate that the solver is NOT initialized.
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_solver_state);
    return;
}// deallocateSolverState

void
PETScLevelSolver::enableLogging(
    bool enabled)
{
    d_enable_logging = enabled;
    return;
}// enableLogging

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
