// Filename: PETScNewtonKrylovSolver.C
// Created on 26 Nov 2003 by Boyce Griffith
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

#include "PETScNewtonKrylovSolver.h"

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
#include <ibtk/PETScKrylovLinearSolver.h>
#include <ibtk/PETScSNESFunctionGOWrapper.h>
#include <ibtk/PETScSNESJacobianJOWrapper.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <PatchHierarchy.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

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

PETScNewtonKrylovSolver::PETScNewtonKrylovSolver(
    const std::string& name,
    const std::string& options_prefix,
    MPI_Comm petsc_comm)
    : d_object_name(name),
      d_is_initialized(false),
      d_reinitializing_solver(false),
      d_do_log(false),
      d_solver_x(NULL),
      d_solver_b(NULL),
      d_solver_r(NULL),
      d_petsc_x(PETSC_NULL),
      d_petsc_b(PETSC_NULL),
      d_petsc_r(PETSC_NULL),
      d_options_prefix(options_prefix),
      d_petsc_comm(petsc_comm),
      d_petsc_snes(PETSC_NULL),
      d_petsc_jac (PETSC_NULL),
      d_managing_petsc_snes(true),
      d_user_provided_function(false),
      d_user_provided_jacobian(false),
      d_F(NULL),
      d_J(NULL),
      d_krylov_solver(NULL),
      d_abs_residual_tol(PETSC_DEFAULT),
      d_rel_residual_tol(PETSC_DEFAULT),
      d_solution_tol(PETSC_DEFAULT),
      d_max_iterations(PETSC_DEFAULT),
      d_max_evaluations(PETSC_DEFAULT),
      d_current_its(0),
      d_current_lits(0),
      d_current_residual_norm(0.0)
{
    // Common constructor functionality.
    common_ctor();
    return;
}// PETScNewtonKrylovSolver()

PETScNewtonKrylovSolver::PETScNewtonKrylovSolver(
    const std::string& name,
    const SNES& petsc_snes,
    const std::string& options_prefix)
    : d_object_name(name),
      d_is_initialized(false),
      d_reinitializing_solver(false),
      d_do_log(false),
      d_solver_x(NULL),
      d_solver_b(NULL),
      d_solver_r(NULL),
      d_petsc_x(PETSC_NULL),
      d_petsc_b(PETSC_NULL),
      d_petsc_r(PETSC_NULL),
      d_options_prefix(options_prefix),
      d_petsc_comm(PETSC_COMM_WORLD),
      d_petsc_snes(petsc_snes),
      d_petsc_jac (PETSC_NULL),
      d_managing_petsc_snes(false),
      d_user_provided_function(false),
      d_user_provided_jacobian(false),
      d_F(NULL),
      d_J(NULL),
      d_krylov_solver(NULL),
      d_abs_residual_tol(PETSC_DEFAULT),
      d_rel_residual_tol(PETSC_DEFAULT),
      d_solution_tol(PETSC_DEFAULT),
      d_max_iterations(PETSC_DEFAULT),
      d_max_evaluations(PETSC_DEFAULT),
      d_current_its(0),
      d_current_lits(0),
      d_current_residual_norm(0.0)
{
    int ierr;

    // Set d_petsc_comm to be the MPI communicator used by the supplied SNES.
    ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(d_petsc_snes), &d_petsc_comm); IBTK_CHKERRQ(ierr);

    // Create an GeneralOperator wrapper to correspond to the SNES function.
    PetscErrorCode (*petsc_snes_form_func)(SNES,Vec,Vec,void*);
    void* petsc_snes_func_ctx;
    ierr = SNESGetFunction(d_petsc_snes, PETSC_NULL, &petsc_snes_form_func, &petsc_snes_func_ctx); IBTK_CHKERRQ(ierr);
    d_F = new PETScSNESFunctionGOWrapper(
        d_object_name+"::SNESFunction Wrapper",
        d_petsc_snes, petsc_snes_form_func, petsc_snes_func_ctx);

    // Create a JacobianOperator wrapper to correspond to the SNES Jacobian.
    PetscErrorCode (*petsc_snes_form_jac)(SNES,Vec,Mat*,Mat*,MatStructure*,void*);
    void* petsc_snes_jac_ctx;
    ierr = SNESGetJacobian(d_petsc_snes, PETSC_NULL, PETSC_NULL, &petsc_snes_form_jac, &petsc_snes_jac_ctx); IBTK_CHKERRQ(ierr);
    d_J = new PETScSNESJacobianJOWrapper(
        d_object_name+"::SNESJacobian Wrapper",
        d_petsc_snes, petsc_snes_form_jac, petsc_snes_jac_ctx);

    // Create a KrylovLinearSolver wrapper to correspond to the KSP.
    KSP petsc_ksp;
    ierr = SNESGetKSP(d_petsc_snes, &petsc_ksp); IBTK_CHKERRQ(ierr);
    d_krylov_solver = new PETScKrylovLinearSolver(
        d_object_name+"::KSP Wrapper", petsc_ksp, d_options_prefix);

    // Reset the member state variables to correspond to the values used by the
    // SNES object.
    ierr = SNESGetTolerances(d_petsc_snes,
                             &d_abs_residual_tol, // absolute residual tol
                             &d_rel_residual_tol, // relative residual tol
                             &d_solution_tol,     // change between solutions
                             &d_max_iterations,   // max Newton iterations
                             &d_max_evaluations); // max function evaluations
    IBTK_CHKERRQ(ierr);

    // Common constructor functionality.
    common_ctor();
    return;
}// PETScNewtonKrylovSolver()

PETScNewtonKrylovSolver::~PETScNewtonKrylovSolver()
{
    if (d_is_initialized) deallocateSolverState();

    // Delete allocated PETSc solver components.
    int ierr;
    if (d_petsc_jac != PETSC_NULL)
    {
        ierr = MatDestroy(d_petsc_jac); IBTK_CHKERRQ(ierr);
        d_petsc_jac = PETSC_NULL;
    }
    if (d_managing_petsc_snes && d_petsc_snes != PETSC_NULL)
    {
        ierr = SNESDestroy(d_petsc_snes); IBTK_CHKERRQ(ierr);
        d_petsc_snes = PETSC_NULL;
    }
    return;
}// ~PETScNewtonKrylovSolver()

bool
PETScNewtonKrylovSolver::solveSystem(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& b)
{
    t_solve_system->start();

    int ierr;

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x,b);

    // Solve the system using a PETSc SNES object.
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_x, Pointer<SAMRAIVectorReal<NDIM,double> >(&x,false));
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_b, Pointer<SAMRAIVectorReal<NDIM,double> >(&b,false));

    ierr = SNESSolve(d_petsc_snes, d_petsc_b, d_petsc_x); IBTK_CHKERRQ(ierr);
    ierr = SNESGetIterationNumber(d_petsc_snes, &d_current_its); IBTK_CHKERRQ(ierr);
    ierr = SNESGetLinearSolveIterations(d_petsc_snes, &d_current_lits); IBTK_CHKERRQ(ierr);
    ierr = SNESGetFunctionNorm(d_petsc_snes, &d_current_residual_norm); IBTK_CHKERRQ(ierr);

    // Determine the convergence reason.
    SNESConvergedReason reason;
    ierr = SNESGetConvergedReason(d_petsc_snes, &reason); IBTK_CHKERRQ(ierr);
    const bool converged = (static_cast<int>(reason) > 0);
    if (d_do_log) reportSNESConvergedReason(reason, plog);

    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();

    t_solve_system->stop();
    return converged;
}// solveSystem

void
PETScNewtonKrylovSolver::initializeSolverState(
    const SAMRAIVectorReal<NDIM,double>& x,
    const SAMRAIVectorReal<NDIM,double>& b)
{
    t_initialize_solver_state->start();

    int ierr;

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
                       <<"  hierarchy level " << ln << " does not exist" << std::endl);
        }
    }
#endif
    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized)
    {
        d_reinitializing_solver = true;
        deallocateSolverState();
    }

    // Create the SNES solver.
    if (d_managing_petsc_snes)
    {
        ierr = SNESCreate(d_petsc_comm, &d_petsc_snes); IBTK_CHKERRQ(ierr);
        resetSNESOptions();

        // Create a KrylovLinearSolver wrapper to correspond to the KSP employed
        // by the SNES solver.
        KSP petsc_ksp;
        ierr = SNESGetKSP(d_petsc_snes, &petsc_ksp); IBTK_CHKERRQ(ierr);
        ierr = KSPSetType(petsc_ksp, KSPFGMRES); IBTK_CHKERRQ(ierr);
        PC petsc_pc;
        ierr = KSPGetPC(petsc_ksp, &petsc_pc); IBTK_CHKERRQ(ierr);
        ierr = PCSetType(petsc_pc, PCNONE); IBTK_CHKERRQ(ierr);

        d_krylov_solver = new PETScKrylovLinearSolver(
            d_object_name+"::KSP Wrapper", petsc_ksp, d_options_prefix);
    }

    // Setup solution and rhs vectors.
    //
    // NOTE: Vector data are allocated only for the residual vector r.
    d_solver_x = x.cloneVector(x.getName());
    d_petsc_x = PETScSAMRAIVectorReal::createPETScVector(d_solver_x, d_petsc_comm);

    d_solver_b = b.cloneVector(b.getName());
    d_petsc_b = PETScSAMRAIVectorReal::createPETScVector(d_solver_b, d_petsc_comm);

    d_solver_r = b.cloneVector(b.getName());
    d_solver_r->allocateVectorData();
    d_petsc_r = PETScSAMRAIVectorReal::createPETScVector(d_solver_r, d_petsc_comm);

    // Setup the nonlinear operator.
    if (!d_F.isNull()) d_F->initializeOperatorState(*d_solver_x, *d_solver_b);
    if (d_managing_petsc_snes || d_user_provided_function) resetSNESFunction();

    // Setup the Jacobian.
    if (!d_J.isNull()) d_J->initializeOperatorState(*d_solver_x, *d_solver_b);
    if (d_managing_petsc_snes || d_user_provided_jacobian) resetSNESJacobian();

    // Setup the Krylov solver.
    if (!d_krylov_solver.isNull()) d_krylov_solver->initializeSolverState(*d_solver_x, *d_solver_b);

    // Set the SNES options from the PETSc options database.
    if (!d_options_prefix.empty())
    {
        ierr = SNESSetOptionsPrefix(d_petsc_snes, d_options_prefix.c_str()); IBTK_CHKERRQ(ierr);
    }
    ierr = SNESSetFromOptions(d_petsc_snes); IBTK_CHKERRQ(ierr);

    // Reset the member state variables to correspond to the values used by the
    // SNES object.
    ierr = SNESGetTolerances(d_petsc_snes,
                             &d_abs_residual_tol, // absolute residual tol
                             &d_rel_residual_tol, // relative residual tol
                             &d_solution_tol,     // change between solutions
                             &d_max_iterations,   // max Newton iterations
                             &d_max_evaluations); // max function evaluations
    IBTK_CHKERRQ(ierr);

    // Indicate that the solver is initialized.
    d_reinitializing_solver = false;
    d_is_initialized = true;

    t_initialize_solver_state->stop();
    return;
}// initializeSolverState

void
PETScNewtonKrylovSolver::deallocateSolverState()
{
    if (!d_is_initialized) return;

    t_deallocate_solver_state->start();

    int ierr;

    // Deallocate the linear solver and operator states only if we are not
    // re-initializing the Newton solver.
    if (!d_reinitializing_solver)
    {
        if (!d_krylov_solver.isNull()) d_krylov_solver->deallocateSolverState();
        if (!d_J.isNull()) d_J->deallocateOperatorState();
        if (!d_F.isNull()) d_F->deallocateOperatorState();
    }

    // Delete the solution and rhs vectors.
    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_x);
    d_petsc_x = PETSC_NULL;
    d_solver_x->resetLevels(d_solver_x->getCoarsestLevelNumber(), std::min(d_solver_x->getFinestLevelNumber(),d_solver_x->getPatchHierarchy()->getFinestLevelNumber()));
    d_solver_x->freeVectorComponents();
    d_solver_x.setNull();

    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_b);
    d_petsc_b = PETSC_NULL;
    d_solver_b->resetLevels(d_solver_b->getCoarsestLevelNumber(), std::min(d_solver_b->getFinestLevelNumber(),d_solver_b->getPatchHierarchy()->getFinestLevelNumber()));
    d_solver_b->freeVectorComponents();
    d_solver_b.setNull();

    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_r);
    d_petsc_r = PETSC_NULL;
    d_solver_r->resetLevels(d_solver_r->getCoarsestLevelNumber(), std::min(d_solver_r->getFinestLevelNumber(),d_solver_r->getPatchHierarchy()->getFinestLevelNumber()));
    d_solver_r->freeVectorComponents();
    d_solver_r.setNull();

    // Destroy the SNES solver.
    if (d_managing_petsc_snes)
    {
        ierr = SNESDestroy(d_petsc_snes);  IBTK_CHKERRQ(ierr);
        d_petsc_snes = PETSC_NULL;
    }

    // Indicate that the solver is NOT initialized.
    d_is_initialized = false;

    t_deallocate_solver_state->stop();
    return;
}// deallocateSolverState

void
PETScNewtonKrylovSolver::enableLogging(
    bool enabled)
{
    d_do_log = enabled;
    return;
}// enableLogging

/////////////////////////////// PRIVATE //////////////////////////////////////

void
PETScNewtonKrylovSolver::common_ctor()
{
    // Setup Timers.
    IBTK_DO_ONCE(
        t_solve_system            = TimerManager::getManager()->getTimer("IBTK::PETScNewtonKrylovSolver::solveSystem()");
        t_initialize_solver_state = TimerManager::getManager()->getTimer("IBTK::PETScNewtonKrylovSolver::initializeOperatorState()");
        t_deallocate_solver_state = TimerManager::getManager()->getTimer("IBTK::PETScNewtonKrylovSolver::deallocateOperatorState()");
                 );
    return;
}// common_ctor

void
PETScNewtonKrylovSolver::reportSNESConvergedReason(
    const SNESConvergedReason& reason,
    std::ostream& os) const
{
    switch (static_cast<int>(reason))
    {
        case SNES_CONVERGED_FNORM_ABS:
            os << d_object_name << ": converged: |F| less than specified absolute tolerance.\n";
            break;
        case SNES_CONVERGED_FNORM_RELATIVE:
            os << d_object_name << ": converged: |F| less than specified relative tolerance.\n";
            break;
        case SNES_CONVERGED_PNORM_RELATIVE:
            os << d_object_name << ": converged: |P| less than specified relative tolerance.\n";
            break;
        case SNES_CONVERGED_TR_DELTA:
            os << d_object_name << ": converged: trust-region delta.\n";
            break;
        case SNES_DIVERGED_FUNCTION_COUNT:
            os << d_object_name << ": diverged: exceeded maximum number of function evaluations.\n";
            break;
        case SNES_DIVERGED_FNORM_NAN:
            os << d_object_name << ": diverged: |F| is NaN.\n";
            break;
        case SNES_DIVERGED_MAX_IT:
            os << d_object_name << ": diverged: exceeded maximum number of iterations.\n";
            break;
        case SNES_DIVERGED_LS_FAILURE:
            os << d_object_name << ": diverged: line-search failure.\n";
            break;
        case SNES_DIVERGED_LOCAL_MIN:
            os << d_object_name << ": diverged: attained non-zero local minimum.\n";
            break;
        case SNES_CONVERGED_ITERATING:
            os << d_object_name << ": iterating.\n";
            break;
        default:
            os << d_object_name << ": unknown completion code " << static_cast<int>(reason) << " reported.\n";
            break;
    }
    return;
}// reportSNESConvergedReason

void
PETScNewtonKrylovSolver::resetSNESOptions()
{
    if (d_petsc_snes != PETSC_NULL)
    {
        int ierr = SNESSetTolerances(d_petsc_snes,
                                     d_abs_residual_tol, // absolute residual tol
                                     d_rel_residual_tol, // relative residual tol
                                     d_solution_tol,     // change between solutions
                                     d_max_iterations,   // max Newton iterations
                                     d_max_evaluations); // max function evaluations
        IBTK_CHKERRQ(ierr);
    }
    return;
}// resetSNESSNESOptions

void
PETScNewtonKrylovSolver::resetSNESFunction()
{
    if (d_petsc_snes != PETSC_NULL)
    {
        int ierr;
        ierr = SNESSetFunction(d_petsc_snes, d_petsc_r, PETScNewtonKrylovSolver::FormFunction_SAMRAI, static_cast<void*>(this)); IBTK_CHKERRQ(ierr);
    }
    return;
}// resetSNESFunction

void
PETScNewtonKrylovSolver::resetSNESJacobian()
{
    if (d_petsc_snes != PETSC_NULL)
    {
        // Create and configure the Jacobian matrix.
        int ierr;
        if (d_petsc_jac != PETSC_NULL)
        {
            ierr = MatDestroy(d_petsc_jac); IBTK_CHKERRQ(ierr);
            d_petsc_jac = PETSC_NULL;
        }
        if (d_J.isNull() || !d_user_provided_jacobian)
        {
            ierr = MatCreateMFFD(d_petsc_comm, 0, 0, 0, 0, &d_petsc_jac); IBTK_CHKERRQ(ierr);
            ierr = MatMFFDSetFunction(d_petsc_jac,
                                      reinterpret_cast<PetscErrorCode(*)(void*, Vec, Vec)>(SNESComputeFunction),
                                      d_petsc_snes); IBTK_CHKERRQ(ierr);
            if (!d_options_prefix.empty())
            {
                ierr = MatSetOptionsPrefix(d_petsc_jac, d_options_prefix.c_str()); IBTK_CHKERRQ(ierr);
            }
            ierr = MatSetFromOptions(d_petsc_jac); IBTK_CHKERRQ(ierr);
        }
        else
        {
            ierr = MatCreateShell(d_petsc_comm,
                                  0, 0,
                                  PETSC_DETERMINE, PETSC_DETERMINE,
                                  static_cast<void*>(this),
                                  &d_petsc_jac); IBTK_CHKERRQ(ierr);
            ierr = MatShellSetOperation(
                d_petsc_jac, MATOP_MULT              , reinterpret_cast<void(*)(void)>(PETScNewtonKrylovSolver::MatVecMult_SAMRAI            )); IBTK_CHKERRQ(ierr);
            ierr = MatShellSetOperation(
                d_petsc_jac, MATOP_MULT_ADD          , reinterpret_cast<void(*)(void)>(PETScNewtonKrylovSolver::MatVecMultAdd_SAMRAI         )); IBTK_CHKERRQ(ierr);
#if 0
            ierr = MatShellSetOperation(
                d_petsc_jac, MATOP_MULT_TRANSPOSE    , reinterpret_cast<void(*)(void)>(PETScNewtonKrylovSolver::MatVecMultTranspose_SAMRAI   )); IBTK_CHKERRQ(ierr);
            ierr = MatShellSetOperation(
                d_petsc_jac, MATOP_MULT_TRANSPOSE_ADD, reinterpret_cast<void(*)(void)>(PETScNewtonKrylovSolver::MatVecMultTransposeAdd_SAMRAI)); IBTK_CHKERRQ(ierr);
#endif
            ierr = MatShellSetOperation(
                d_petsc_jac, MATOP_GET_VECS          , reinterpret_cast<void(*)(void)>(PETScNewtonKrylovSolver::MatGetVecs_SAMRAI            )); IBTK_CHKERRQ(ierr);
        }

        // Reset the configuration of the PETSc SNES object.
        ierr = SNESSetJacobian(d_petsc_snes, d_petsc_jac, d_petsc_jac,
                               PETScNewtonKrylovSolver::FormJacobian_SAMRAI,
                               static_cast<void*>(this)); IBTK_CHKERRQ(ierr);
    }
    return;
}// resetSNESJacobian

PetscErrorCode
PETScNewtonKrylovSolver::FormFunction_SAMRAI(
    SNES snes,
    Vec x,
    Vec f,
    void* p_ctx)
{
    PETScNewtonKrylovSolver* newton_solver = static_cast<PETScNewtonKrylovSolver*>(p_ctx);
#if (DEBUG_CHECK_ASSERTIONS)
    TBOX_ASSERT(newton_solver != NULL);
    TBOX_ASSERT(!newton_solver->d_F.isNull());
#endif
    int ierr;
    newton_solver->d_F->apply(*PETScSAMRAIVectorReal::getSAMRAIVector(x),
                              *PETScSAMRAIVectorReal::getSAMRAIVector(f));
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(f)); IBTK_CHKERRQ(ierr);
    PetscFunctionReturn(0);
}// FormFunction_SAMRAI

PetscErrorCode
PETScNewtonKrylovSolver::FormJacobian_SAMRAI(
    SNES snes,
    Vec x,
    Mat* A,
    Mat* B,
    MatStructure* mat_structure,
    void* p_ctx)
{
    (void) snes;
    (void) B;
    PETScNewtonKrylovSolver* newton_solver = static_cast<PETScNewtonKrylovSolver*>(p_ctx);
#if (DEBUG_CHECK_ASSERTIONS)
    TBOX_ASSERT(newton_solver != NULL);
#endif
    if (newton_solver->d_J.isNull())
    {
        int ierr;
        Vec u, f;
        ierr = SNESGetSolution(snes, &u); IBTK_CHKERRQ(ierr);
        ierr = SNESGetFunction(snes, &f, PETSC_NULL, PETSC_NULL); IBTK_CHKERRQ(ierr);
        ierr = MatMFFDSetBase(*A, u, f); IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY); IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY); IBTK_CHKERRQ(ierr);
    }
    else
    {
        newton_solver->d_J->formJacobian(*PETScSAMRAIVectorReal::getSAMRAIVector(x));
    }
    PetscFunctionReturn(0);
}// FormJacobian_SAMRAI

PetscErrorCode
PETScNewtonKrylovSolver::MatVecMult_SAMRAI(
    Mat A,
    Vec x,
    Vec y)
{
    int ierr;
    void* p_ctx;
    ierr = MatShellGetContext(A, &p_ctx); IBTK_CHKERRQ(ierr);
    PETScNewtonKrylovSolver* newton_solver = static_cast<PETScNewtonKrylovSolver*>(p_ctx);
#if (DEBUG_CHECK_ASSERTIONS)
    TBOX_ASSERT(newton_solver != NULL);
    TBOX_ASSERT(!newton_solver->d_J.isNull());
#endif
    newton_solver->d_J->apply(*PETScSAMRAIVectorReal::getSAMRAIVector(x),
                              *PETScSAMRAIVectorReal::getSAMRAIVector(y));
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y)); IBTK_CHKERRQ(ierr);
    PetscFunctionReturn(0);
}// MatVecMult_SAMRAI

PetscErrorCode
PETScNewtonKrylovSolver::MatVecMultAdd_SAMRAI(
    Mat A,
    Vec x,
    Vec y,
    Vec z)
{
    int ierr;
    void* p_ctx;
    ierr = MatShellGetContext(A, &p_ctx); IBTK_CHKERRQ(ierr);
    PETScNewtonKrylovSolver* newton_solver = static_cast<PETScNewtonKrylovSolver*>(p_ctx);
#if (DEBUG_CHECK_ASSERTIONS)
    TBOX_ASSERT(newton_solver != NULL);
    TBOX_ASSERT(!newton_solver->d_J.isNull());
#endif
    newton_solver->d_J->applyAdd(*PETScSAMRAIVectorReal::getSAMRAIVector(x),
                                 *PETScSAMRAIVectorReal::getSAMRAIVector(y),
                                 *PETScSAMRAIVectorReal::getSAMRAIVector(z));
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(z)); IBTK_CHKERRQ(ierr);
    PetscFunctionReturn(0);
}// MatVecMultAdd_SAMRAI

PetscErrorCode
PETScNewtonKrylovSolver::MatVecMultTranspose_SAMRAI(
    Mat A,
    Vec x,
    Vec y)
{
    int ierr;
    void* p_ctx;
    ierr = MatShellGetContext(A, &p_ctx); IBTK_CHKERRQ(ierr);
    PETScNewtonKrylovSolver* newton_solver = static_cast<PETScNewtonKrylovSolver*>(p_ctx);
#if (DEBUG_CHECK_ASSERTIONS)
    TBOX_ASSERT(newton_solver != NULL);
    TBOX_ASSERT(!newton_solver->d_J.isNull());
#endif
    if (newton_solver->d_J->isSymmetric())
    {
        PetscFunctionReturn(PETScNewtonKrylovSolver::MatVecMult_SAMRAI(A,x,y));
    }
    newton_solver->d_J->applyAdjoint(*PETScSAMRAIVectorReal::getSAMRAIVector(x),
                                     *PETScSAMRAIVectorReal::getSAMRAIVector(y));
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y)); IBTK_CHKERRQ(ierr);
    PetscFunctionReturn(0);
}// MatVecMultTranspose_SAMRAI

PetscErrorCode
PETScNewtonKrylovSolver::MatVecMultTransposeAdd_SAMRAI(
    Mat A,
    Vec x,
    Vec y,
    Vec z)
{
    int ierr;
    void* p_ctx;
    ierr = MatShellGetContext(A, &p_ctx); IBTK_CHKERRQ(ierr);
    PETScNewtonKrylovSolver* newton_solver = static_cast<PETScNewtonKrylovSolver*>(p_ctx);
#if (DEBUG_CHECK_ASSERTIONS)
    TBOX_ASSERT(newton_solver != NULL);
    TBOX_ASSERT(!newton_solver->d_J.isNull());
#endif
    if (newton_solver->d_J->isSymmetric())
    {
        PetscFunctionReturn(PETScNewtonKrylovSolver::MatVecMultAdd_SAMRAI(A,x,y,z));
    }
    newton_solver->d_J->applyAdjointAdd(*PETScSAMRAIVectorReal::getSAMRAIVector(x),
                                        *PETScSAMRAIVectorReal::getSAMRAIVector(y),
                                        *PETScSAMRAIVectorReal::getSAMRAIVector(z));
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(z)); IBTK_CHKERRQ(ierr);
    PetscFunctionReturn(0);
}// MatVecMultTransposeAdd_SAMRAI

PetscErrorCode
PETScNewtonKrylovSolver::MatGetVecs_SAMRAI(
    Mat A,
    Vec* right,
    Vec* left)
{
    int ierr;
    void* p_ctx;
    ierr = MatShellGetContext(A, &p_ctx); IBTK_CHKERRQ(ierr);
    PETScNewtonKrylovSolver* newton_solver = static_cast<PETScNewtonKrylovSolver*>(p_ctx);
#if (DEBUG_CHECK_ASSERTIONS)
    TBOX_ASSERT(newton_solver != NULL);
#endif
    if (right != PETSC_NULL)
    {
        // vector that the matrix can be multiplied against
        ierr = VecDuplicate(newton_solver->d_petsc_x, right); IBTK_CHKERRQ(ierr);
        ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(*right)); IBTK_CHKERRQ(ierr);
    }
    if (left != PETSC_NULL)
    {
        // vector that the matrix vector product can be stored in
        ierr = VecDuplicate(newton_solver->d_petsc_b, left); IBTK_CHKERRQ(ierr);
        ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(*left)); IBTK_CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}// MatGetVecs_SAMRAI

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBTK::PETScNewtonKrylovSolver>;

//////////////////////////////////////////////////////////////////////////////
