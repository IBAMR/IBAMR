// Filename: PETScKrylovLinearSolver.C
// Created on 16 Sep 2003 by Boyce Griffith
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

#include "PETScKrylovLinearSolver.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// PETSc INCLUDES
#include <private/pcimpl.h>

// IBTK INCLUDES
#include <ibtk/PETScSAMRAIVectorReal.h>
#include <ibtk/PETScMatLOWrapper.h>
#include <ibtk/PETScPCLSWrapper.h>
#include <ibtk/IBTK_CHKERRQ.h>
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

PETScKrylovLinearSolver::PETScKrylovLinearSolver(
    const std::string& object_name,
    const std::string& options_prefix,
    MPI_Comm petsc_comm)
    : d_object_name(object_name),
      d_ksp_type("gmres"),
      d_pc_shell_types(),
      d_pc_shell_type("none"),
      d_is_initialized(false),
      d_reinitializing_solver(false),
      d_do_log(false),
      d_solver_x(NULL),
      d_solver_b(NULL),
      d_petsc_x(PETSC_NULL),
      d_petsc_b(PETSC_NULL),
      d_options_prefix(options_prefix),
      d_petsc_comm  (petsc_comm),
      d_petsc_ksp   (PETSC_NULL),
      d_petsc_mat   (PETSC_NULL),
      d_petsc_nullsp(PETSC_NULL),
      d_managing_petsc_ksp(true),
      d_user_provided_mat(false),
      d_user_provided_pc (false),
      d_A(NULL),
      d_pc_solver(NULL),
      d_nullsp_contains_constant_vector(false),
      d_solver_nullsp_constant(NULL),
      d_solver_nullsp_basis(),
      d_petsc_nullsp_constant(PETSC_NULL),
      d_petsc_nullsp_basis(),
      d_solver_has_attached_nullsp(false),
      d_initial_guess_nonzero(false),
      d_rel_residual_tol(PETSC_DEFAULT),
      d_abs_residual_tol(PETSC_DEFAULT),
      d_divergence_tol(PETSC_DEFAULT),
      d_max_iterations(PETSC_DEFAULT),
      d_current_its(0),
      d_current_residual_norm(0.0)
{
    // Common constructor functionality.
    common_ctor();
    return;
}// PETScKrylovLinearSolver()

PETScKrylovLinearSolver::PETScKrylovLinearSolver(
    const std::string& object_name,
    const KSP& petsc_ksp,
    const std::string& options_prefix)
    : d_object_name(object_name),
      d_ksp_type("none"),
      d_is_initialized(false),
      d_reinitializing_solver(false),
      d_do_log(false),
      d_solver_x(NULL),
      d_solver_b(NULL),
      d_petsc_x(PETSC_NULL),
      d_petsc_b(PETSC_NULL),
      d_options_prefix(options_prefix),
      d_petsc_comm  (PETSC_COMM_WORLD),
      d_petsc_ksp   (petsc_ksp),
      d_petsc_mat   (PETSC_NULL),
      d_petsc_nullsp(PETSC_NULL),
      d_managing_petsc_ksp(false),
      d_user_provided_mat(false),
      d_user_provided_pc (false),
      d_A(NULL),
      d_pc_solver(NULL),
      d_nullsp_contains_constant_vector(false),
      d_solver_nullsp_constant(NULL),
      d_solver_nullsp_basis(),
      d_petsc_nullsp_constant(PETSC_NULL),
      d_petsc_nullsp_basis(),
      d_solver_has_attached_nullsp(false),
      d_initial_guess_nonzero(false),
      d_rel_residual_tol(PETSC_DEFAULT),
      d_abs_residual_tol(PETSC_DEFAULT),
      d_divergence_tol(PETSC_DEFAULT),
      d_max_iterations(PETSC_DEFAULT),
      d_current_its(0),
      d_current_residual_norm(0.0)
{
    if (d_petsc_ksp != NULL) resetWrappedKSP(d_petsc_ksp);

    // Common constructor functionality.
    common_ctor();
    return;
}// PETScKrylovLinearSolver()

PETScKrylovLinearSolver::~PETScKrylovLinearSolver()
{
    if (d_is_initialized) deallocateSolverState();

    // Delete allocated PETSc solver components.
    int ierr;
    if (d_petsc_mat != PETSC_NULL)
    {
        ierr = MatDestroy(&d_petsc_mat); IBTK_CHKERRQ(ierr);
        d_petsc_mat = PETSC_NULL;
    }
    if (d_managing_petsc_ksp && d_petsc_ksp != PETSC_NULL)
    {
        ierr = KSPDestroy(&d_petsc_ksp); IBTK_CHKERRQ(ierr);
        d_petsc_ksp = PETSC_NULL;
    }
    return;
}// ~PETScKrylovLinearSolver()

void
PETScKrylovLinearSolver::setValidPCShellTypes(
    const std::vector<std::string>& pc_shell_types)
{
    d_pc_shell_types = pc_shell_types;
    return;
}// setValidPCShellTypes

const std::string&
PETScKrylovLinearSolver::getPCShellType() const
{
    return d_pc_shell_type;
}// getPCShellType

bool
PETScKrylovLinearSolver::solveSystem(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& b)
{
    IBTK_TIMER_START(t_solve_system);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_A.isNull());
#endif
    int ierr;

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x,b);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_petsc_ksp != PETSC_NULL);
#endif

    // Solve the system using a PETSc KSP object.
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_x, Pointer<SAMRAIVectorReal<NDIM,double> >(&x,false));
    d_solver_b->copyVector(Pointer<SAMRAIVectorReal<NDIM,double> >(&b,false));
    d_A->modifyRhsForInhomogeneousBc(*d_solver_b);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(d_petsc_b)); IBTK_CHKERRQ(ierr);

    ierr = KSPSolve(d_petsc_ksp, d_petsc_b, d_petsc_x); IBTK_CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(d_petsc_ksp, &d_current_its); IBTK_CHKERRQ(ierr);
    ierr = KSPGetResidualNorm(d_petsc_ksp, &d_current_residual_norm); IBTK_CHKERRQ(ierr);

    // Determine the convergence reason.
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(d_petsc_ksp, &reason); IBTK_CHKERRQ(ierr);
    const bool converged = (static_cast<int>(reason) > 0);
    if (d_do_log) reportKSPConvergedReason(reason, plog);

    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();

    IBTK_TIMER_STOP(t_solve_system);
    return converged;
}// solveSystem

void
PETScKrylovLinearSolver::initializeSolverState(
    const SAMRAIVectorReal<NDIM,double>& x,
    const SAMRAIVectorReal<NDIM,double>& b)
{
    IBTK_TIMER_START(t_initialize_solver_state);

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
                       << "  hierarchy level " << ln << " does not exist" << std::endl);
        }
    }
#endif
    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized)
    {
        d_reinitializing_solver = true;
        deallocateSolverState();
    }

    // Create the KSP solver.
    if (d_managing_petsc_ksp)
    {
        ierr = KSPCreate(d_petsc_comm, &d_petsc_ksp); IBTK_CHKERRQ(ierr);
        resetKSPOptions();
    }
    else if (d_petsc_ksp == PETSC_NULL)
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                   << "  cannot initialize solver state for wrapped PETSc KSP object if the wrapped object is NULL" << std::endl);
    }

    // Setup solution and rhs vectors.
    d_solver_x = x.cloneVector(x.getName());
    d_petsc_x = PETScSAMRAIVectorReal::createPETScVector(d_solver_x, d_petsc_comm);

    d_solver_b = b.cloneVector(b.getName());
    d_solver_b->allocateVectorData();
    d_petsc_b = PETScSAMRAIVectorReal::createPETScVector(d_solver_b, d_petsc_comm);

    // Initialize the linear operator and preconditioner objects.
    if (!d_A.isNull()) d_A->initializeOperatorState(*d_solver_x, *d_solver_b);
    if (d_managing_petsc_ksp || d_user_provided_mat) resetKSPOperators();

    if (!d_pc_solver.isNull()) d_pc_solver->initializeSolverState(*d_solver_x, *d_solver_b);
    if (d_managing_petsc_ksp || d_user_provided_pc) resetKSPPC();

    // Set the KSP options from the PETSc options database.
    if (!d_options_prefix.empty())
    {
        ierr = KSPSetOptionsPrefix(d_petsc_ksp, d_options_prefix.c_str()); IBTK_CHKERRQ(ierr);
    }
    ierr = KSPSetFromOptions(d_petsc_ksp); IBTK_CHKERRQ(ierr);

    // Reset the member state variables to correspond to the values used by the
    // KSP object.  (Command-line options always take precedence.)
    const KSPType ksp_type;
    ierr = KSPGetType(d_petsc_ksp, &ksp_type); IBTK_CHKERRQ(ierr);
    d_ksp_type = ksp_type;
    PetscBool initial_guess_nonzero;
    ierr = KSPGetInitialGuessNonzero(d_petsc_ksp, &initial_guess_nonzero); IBTK_CHKERRQ(ierr);
    d_initial_guess_nonzero = (initial_guess_nonzero == PETSC_TRUE);
    ierr = KSPGetTolerances(d_petsc_ksp, &d_rel_residual_tol, &d_abs_residual_tol, &d_divergence_tol, &d_max_iterations); IBTK_CHKERRQ(ierr);

    // Configure the nullspace object.
    resetKSPNullspace();

    // Indicate that the solver is initialized.
    d_reinitializing_solver = false;
    d_is_initialized = true;

    IBTK_TIMER_STOP(t_initialize_solver_state);
    return;
}// initializeSolverState

void
PETScKrylovLinearSolver::deallocateSolverState()
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_deallocate_solver_state);

    int ierr;

    // Deallocate the operator and preconditioner states only if we are not
    // re-initializing the solver.
    if (!d_reinitializing_solver)
    {
        if (!d_pc_solver.isNull()) d_pc_solver->deallocateSolverState();
        if (!d_A.isNull()) d_A->deallocateOperatorState();
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

    // Deallocate the nullspace object.
    deallocateNullspaceData();

    // Destroy the KSP solver.
    if (d_managing_petsc_ksp)
    {
        ierr = KSPDestroy(&d_petsc_ksp);  IBTK_CHKERRQ(ierr);
        d_petsc_ksp = PETSC_NULL;
        d_solver_has_attached_nullsp = false;
    }

    // Indicate that the solver is NOT initialized.
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_solver_state);
    return;
}// deallocateSolverState

void
PETScKrylovLinearSolver::enableLogging(
    bool enabled)
{
    d_do_log = enabled;
    return;
}// enableLogging

/////////////////////////////// PRIVATE //////////////////////////////////////

void
PETScKrylovLinearSolver::common_ctor()
{
    // Setup Timers.
    IBTK_DO_ONCE(
        t_solve_system            = TimerManager::getManager()->getTimer("IBTK::PETScKrylovLinearSolver::solveSystem()");
        t_initialize_solver_state = TimerManager::getManager()->getTimer("IBTK::PETScKrylovLinearSolver::initializeSolverState()");
        t_deallocate_solver_state = TimerManager::getManager()->getTimer("IBTK::PETScKrylovLinearSolver::deallocateSolverState()");
                 );
    return;
}// common_ctor

void
PETScKrylovLinearSolver::reportKSPConvergedReason(
    const KSPConvergedReason& reason,
    std::ostream& os) const
{
    switch (static_cast<int>(reason))
    {
        case KSP_CONVERGED_RTOL:
            os << d_object_name << ": converged: |Ax-b| <= rtol*|b| --- residual norm is less than specified relative tolerance.\n";
            break;
        case KSP_CONVERGED_ATOL:
            os << d_object_name << ": converged: |Ax-b| <= atol --- residual norm is less than specified absolute tolerance.\n";
            break;
        case KSP_CONVERGED_ITS:
            os << d_object_name << ": converged: single iteration of KSPPREONLY.\n";
            break;
        case KSP_CONVERGED_STEP_LENGTH:
            os << d_object_name << ": converged: step size less than specified tolerance.\n";
            break;
        case KSP_DIVERGED_NULL:
            os << d_object_name << ": diverged: null.\n";
            break;
        case KSP_DIVERGED_ITS:
            os << d_object_name << ": diverged: reached maximum number of iterations before any convergence criteria were satisfied.\n";
            break;
        case KSP_DIVERGED_DTOL:
            os << d_object_name << ": diverged: |Ax-b| >= dtol*|b| --- residual is greater than specified divergence tolerance.\n";
            break;
        case KSP_DIVERGED_BREAKDOWN:
            os << d_object_name << ": diverged: breakdown in the Krylov method.\n";
            break;
        case KSP_DIVERGED_BREAKDOWN_BICG:
            os << d_object_name << ": diverged: breakdown in the bi-congugate gradient method.\n";
            break;
        case KSP_DIVERGED_NONSYMMETRIC:
            os << d_object_name << ": diverged: it appears the operator or preconditioner is not symmetric, but this Krylov method (KSPCG, KSPMINRES, KSPCR) requires symmetry\n";
            break;
        case KSP_DIVERGED_INDEFINITE_PC:
            os << d_object_name << ": diverged: it appears the preconditioner is indefinite (has both positive and negative eigenvalues), but this Krylov method (KSPCG) requires it to be positive definite.\n";
            break;
        case KSP_CONVERGED_ITERATING:
            os << d_object_name << ": iterating: KSPSolve() is still running.\n";
            break;
        default:
            os << d_object_name << ": unknown completion code " << static_cast<int>(reason) << " reported.\n";
            break;
    }
    return;
}// reportKSPConvergedReason

void
PETScKrylovLinearSolver::resetWrappedKSP(
    KSP& petsc_ksp)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_managing_petsc_ksp);
#endif
    d_petsc_ksp = petsc_ksp;
    if (d_petsc_ksp == PETSC_NULL) return;
    int ierr;

    // Set d_petsc_comm to be the MPI communicator used by the supplied KSP.
    ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(d_petsc_ksp), &d_petsc_comm); IBTK_CHKERRQ(ierr);

    // Set d_ksp_type to correspond to the KSP type used by the supplied KSP.
    const KSPType ksp_type;
    ierr = KSPGetType(d_petsc_ksp, &ksp_type); IBTK_CHKERRQ(ierr);
    d_ksp_type = std::string(ksp_type);

    // Setup operators and preconditioners.
    if (d_user_provided_mat) resetKSPOperators();
    else
    {
        // Create a LinearOperator wrapper to correspond to the PETSc Mat used
        // by the KSP.
        Mat petsc_mat;
        ierr = KSPGetOperators(d_petsc_ksp, &petsc_mat, PETSC_NULL, PETSC_NULL); IBTK_CHKERRQ(ierr);
        d_A = new PETScMatLOWrapper(d_object_name+"::Mat Wrapper", petsc_mat);
        d_A->setHomogeneousBc(d_homogeneous_bc);
        d_A->setSolutionTime(d_solution_time);
        d_A->setTimeInterval(d_current_time, d_new_time);
    }

    if (d_user_provided_pc) resetKSPPC();
    else
    {
        // Create a LinearSolver wrapper to correspond to the PETSc PC used by
        // the KSP.
        PC petsc_pc;
        ierr = KSPGetPC(d_petsc_ksp, &petsc_pc); IBTK_CHKERRQ(ierr);
        d_pc_solver = new PETScPCLSWrapper(d_object_name+"::PC Wrapper", petsc_pc);
        d_pc_solver->setHomogeneousBc(true);
        d_pc_solver->setSolutionTime(d_solution_time);
        d_pc_solver->setTimeInterval(d_current_time, d_new_time);
    }

    // Reset the member state variables to correspond to the values used by the
    // KSP object.
    PetscBool initial_guess_nonzero;
    ierr = KSPGetInitialGuessNonzero(d_petsc_ksp, &initial_guess_nonzero); IBTK_CHKERRQ(ierr);
    d_initial_guess_nonzero = (initial_guess_nonzero == PETSC_TRUE);
    ierr = KSPGetTolerances(d_petsc_ksp, &d_rel_residual_tol, &d_abs_residual_tol, &d_divergence_tol, &d_max_iterations); IBTK_CHKERRQ(ierr);
    return;
}// resetWrappedKSP

void
PETScKrylovLinearSolver::resetKSPOptions()
{
    if (d_petsc_ksp == PETSC_NULL) return;
    int ierr;
    const KSPType ksp_type = d_ksp_type.c_str();
    ierr = KSPSetType(d_petsc_ksp, ksp_type); IBTK_CHKERRQ(ierr);
    std::string ksp_type_name(ksp_type);
    if (ksp_type_name.find("gmres") != std::string::npos)
    {
        ierr = KSPGMRESSetCGSRefinementType(d_petsc_ksp, KSP_GMRES_CGS_REFINE_IFNEEDED); IBTK_CHKERRQ(ierr);
    }
    PetscBool initial_guess_nonzero = (d_initial_guess_nonzero ? PETSC_TRUE : PETSC_FALSE);
    ierr = KSPSetInitialGuessNonzero(d_petsc_ksp, initial_guess_nonzero); IBTK_CHKERRQ(ierr);
    ierr = KSPSetTolerances(d_petsc_ksp, d_rel_residual_tol, d_abs_residual_tol, d_divergence_tol, d_max_iterations); IBTK_CHKERRQ(ierr);
    return;
}// resetKSPOptions

void
PETScKrylovLinearSolver::resetKSPOperators()
{
    int ierr;

    // Create and configure the MatShell object.
    if (d_petsc_mat != PETSC_NULL)
    {
        const MatType mat_type;
        ierr = MatGetType(d_petsc_mat, &mat_type); IBTK_CHKERRQ(ierr);
        if (strcmp(mat_type,MATSHELL))
        {
            ierr = MatDestroy(&d_petsc_mat); IBTK_CHKERRQ(ierr);
            d_petsc_mat = PETSC_NULL;
        }
    }
    if (d_petsc_mat == PETSC_NULL)
    {
        ierr = MatCreateShell(d_petsc_comm, 0, 0, 0, 0, static_cast<void*>(this), &d_petsc_mat); IBTK_CHKERRQ(ierr);
    }
    ierr = MatShellSetOperation(d_petsc_mat, MATOP_MULT    , reinterpret_cast<void(*)(void)>(PETScKrylovLinearSolver::MatVecMult_SAMRAI   )); IBTK_CHKERRQ(ierr);
    ierr = MatShellSetOperation(d_petsc_mat, MATOP_MULT_ADD, reinterpret_cast<void(*)(void)>(PETScKrylovLinearSolver::MatVecMultAdd_SAMRAI)); IBTK_CHKERRQ(ierr);
    ierr = MatShellSetOperation(d_petsc_mat, MATOP_GET_VECS, reinterpret_cast<void(*)(void)>(PETScKrylovLinearSolver::MatGetVecs_SAMRAI   )); IBTK_CHKERRQ(ierr);

    // Reset the configuration of the PETSc KSP object.
    if (d_petsc_ksp != PETSC_NULL)
    {
        ierr = KSPSetOperators(d_petsc_ksp, d_petsc_mat, d_petsc_mat, SAME_PRECONDITIONER); IBTK_CHKERRQ(ierr);
    }
    return;
}// resetKSPOperators

void
PETScKrylovLinearSolver::resetKSPPC()
{
    if (d_petsc_ksp == PETSC_NULL) return;
    int ierr;

    // Determine the preconditioner type to use.
    static const size_t len = 255;
    char pc_type_str[len];
    PetscBool flg;
    ierr = PetscOptionsGetString(d_options_prefix.c_str(), "-pc_type", pc_type_str, len, &flg);  IBTK_CHKERRQ(ierr);
    std::string pc_type = "shell";
    if (flg)
    {
        pc_type = std::string(pc_type_str);
    }

    if (!(pc_type == "none" || pc_type == "shell"))
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                   << "  valid values for -" << d_options_prefix << "pc_type are: none, shell" << std::endl);
    }

    PC petsc_pc;
    ierr = KSPGetPC(d_petsc_ksp, &petsc_pc); IBTK_CHKERRQ(ierr);
    if (pc_type == "none" || d_pc_solver.isNull())
    {
        ierr = PCSetType(petsc_pc, PCNONE); IBTK_CHKERRQ(ierr);
    }
    else if (pc_type == "shell" && !d_pc_solver.isNull())
    {
        ierr = PCSetType(petsc_pc, PCSHELL); IBTK_CHKERRQ(ierr);
        ierr = PCShellSetContext(petsc_pc, static_cast<void*>(this)); IBTK_CHKERRQ(ierr);
        ierr = PCShellSetApply(petsc_pc, PETScKrylovLinearSolver::PCApply_SAMRAI); IBTK_CHKERRQ(ierr);
        ierr = PCShellSetName(petsc_pc, "PETScKrylovLinearSolver PC"); IBTK_CHKERRQ(ierr);
        petsc_pc->ops->setfromoptions = PCShellSetFromOptions_SAMRAI;
    }
    else
    {
        TBOX_ERROR("this statement should not be reached!\n");
    }
    return;
}// resetKSPPC

void
PETScKrylovLinearSolver::resetKSPNullspace()
{
    if (d_petsc_ksp == PETSC_NULL) return;
    int ierr;
    PetscBool flg;
    ierr = PetscOptionsHasName(d_options_prefix.c_str(), "-ksp_constant_null_space", &flg); IBTK_CHKERRQ(ierr);
    if (flg == PETSC_TRUE) d_nullsp_contains_constant_vector = true;
    if (d_nullsp_contains_constant_vector || !d_solver_nullsp_basis.empty())
    {
        std::vector<Vec> nullspace_vecs;
        nullspace_vecs.reserve(1+d_solver_nullsp_basis.size());
        if (d_nullsp_contains_constant_vector)
        {
            d_solver_nullsp_constant = d_solver_x->cloneVector(d_solver_x->getName());
            d_solver_nullsp_constant->allocateVectorData();
            d_petsc_nullsp_constant = PETScSAMRAIVectorReal::createPETScVector(d_solver_nullsp_constant, d_petsc_comm);
            ierr = VecSet(d_petsc_nullsp_constant, 1.0); IBTK_CHKERRQ(ierr);
            nullspace_vecs.push_back(d_petsc_nullsp_constant);
        }

        d_petsc_nullsp_basis.resize(d_solver_nullsp_basis.size());
        for (unsigned int k = 0; k < d_solver_nullsp_basis.size(); ++k)
        {
            d_petsc_nullsp_basis[k] = PETScSAMRAIVectorReal::createPETScVector(d_solver_nullsp_basis[k], d_petsc_comm);
            nullspace_vecs.push_back(d_petsc_nullsp_basis[k]);
        }

        for (unsigned int k = 0; k < nullspace_vecs.size(); ++k)
        {
            Vec petsc_nvec = nullspace_vecs[k];
            double dot;
            ierr = VecDot(petsc_nvec, petsc_nvec, &dot); IBTK_CHKERRQ(ierr);
            ierr = VecScale(petsc_nvec, 1.0/sqrt(dot)); IBTK_CHKERRQ(ierr);
        }

        static const PetscBool has_cnst = PETSC_FALSE;
        ierr = MatNullSpaceCreate(d_petsc_comm, has_cnst, nullspace_vecs.size(), &nullspace_vecs[0], &d_petsc_nullsp); IBTK_CHKERRQ(ierr);
        ierr = KSPSetNullSpace(d_petsc_ksp, d_petsc_nullsp); IBTK_CHKERRQ(ierr);
        d_solver_has_attached_nullsp = true;
    }
    else if (d_solver_has_attached_nullsp)
    {
        TBOX_ERROR(d_object_name << "::resetKSPNullspace():\n"
                   << "  it is not possible to remove the nullspace from a PETSc KSP object\n");
    }
    return;
}// resetKSPNullspace

void
PETScKrylovLinearSolver::deallocateNullspaceData()
{
    int ierr;

    if (d_petsc_nullsp != PETSC_NULL)
    {
        ierr = MatNullSpaceDestroy(&d_petsc_nullsp); IBTK_CHKERRQ(ierr);
        d_petsc_nullsp = PETSC_NULL;
    }

    if (d_petsc_nullsp_constant != PETSC_NULL)
    {
        PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_nullsp_constant);
        d_petsc_nullsp_constant = PETSC_NULL;
        d_solver_nullsp_constant->resetLevels(d_solver_nullsp_constant->getCoarsestLevelNumber(), std::min(d_solver_nullsp_constant->getFinestLevelNumber(),d_solver_nullsp_constant->getPatchHierarchy()->getFinestLevelNumber()));
        d_solver_nullsp_constant->freeVectorComponents();
        d_solver_nullsp_constant.setNull();
    }

    for (unsigned int k = 0; k < d_petsc_nullsp_basis.size(); ++k)
    {
        PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_nullsp_basis[k]);
    }
    d_petsc_nullsp_basis.clear();
    return;
}// deallocateNullspaceData

PetscErrorCode
PETScKrylovLinearSolver::PCShellSetFromOptions_SAMRAI(
    PC pc)
{
    PetscErrorCode ierr;
    void* p_ctx;
    ierr = PCShellGetContext(pc, &p_ctx); IBTK_CHKERRQ(ierr);
    PETScKrylovLinearSolver* krylov_solver = static_cast<PETScKrylovLinearSolver*>(p_ctx);
    if (!krylov_solver->d_pc_shell_types.empty())
    {
        std::vector<const char*> pc_shell_types(krylov_solver->d_pc_shell_types.size());
        for (unsigned int k = 0; k < krylov_solver->d_pc_shell_types.size(); ++k)
        {
            pc_shell_types[k] = krylov_solver->d_pc_shell_types[k].c_str();
        }
        int pc_shell_type;
        PetscBool pc_shell_type_set;
        ierr = PetscOptionsEList("-pc_shell_type","Type of shell preconditioner","PCSetFromOptions",&pc_shell_types[0],pc_shell_types.size(),pc_shell_types[0],&pc_shell_type,&pc_shell_type_set); IBTK_CHKERRQ(ierr);
        if (pc_shell_type_set)
        {
            krylov_solver->d_pc_shell_type = krylov_solver->d_pc_shell_types[pc_shell_type];
        }
        else
        {
            krylov_solver->d_pc_shell_type = krylov_solver->d_pc_shell_types[0];
        }
    }
    PetscFunctionReturn(0);
}// PCShellSetFromOptions_SAMRAI

PetscErrorCode
PETScKrylovLinearSolver::MatVecMult_SAMRAI(
    Mat A,
    Vec x,
    Vec y)
{
    int ierr;
    void* p_ctx;
    ierr = MatShellGetContext(A, &p_ctx); IBTK_CHKERRQ(ierr);
    PETScKrylovLinearSolver* krylov_solver = static_cast<PETScKrylovLinearSolver*>(p_ctx);
#if (DEBUG_CHECK_ASSERTIONS)
    TBOX_ASSERT(krylov_solver != NULL);
    TBOX_ASSERT(!krylov_solver->d_A.isNull());
#endif
    krylov_solver->d_A->apply(*PETScSAMRAIVectorReal::getSAMRAIVector(x), *PETScSAMRAIVectorReal::getSAMRAIVector(y));
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y)); IBTK_CHKERRQ(ierr);
    PetscFunctionReturn(0);
}// MatVecMult_SAMRAI

PetscErrorCode
PETScKrylovLinearSolver::MatVecMultAdd_SAMRAI(
    Mat A,
    Vec x,
    Vec y,
    Vec z)
{
    int ierr;
    void* p_ctx;
    ierr = MatShellGetContext(A, &p_ctx); IBTK_CHKERRQ(ierr);
    PETScKrylovLinearSolver* krylov_solver = static_cast<PETScKrylovLinearSolver*>(p_ctx);
#if (DEBUG_CHECK_ASSERTIONS)
    TBOX_ASSERT(krylov_solver != NULL);
    TBOX_ASSERT(!krylov_solver->d_A.isNull());
#endif
    krylov_solver->d_A->applyAdd(*PETScSAMRAIVectorReal::getSAMRAIVector(x), *PETScSAMRAIVectorReal::getSAMRAIVector(y), *PETScSAMRAIVectorReal::getSAMRAIVector(z));
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(z)); IBTK_CHKERRQ(ierr);
    PetscFunctionReturn(0);
}// MatVecMultAdd_SAMRAI

PetscErrorCode
PETScKrylovLinearSolver::MatGetVecs_SAMRAI(
    Mat A,
    Vec* right,
    Vec* left)
{
    int ierr;
    void* p_ctx;
    ierr = MatShellGetContext(A, &p_ctx); IBTK_CHKERRQ(ierr);
    PETScKrylovLinearSolver* krylov_solver = static_cast<PETScKrylovLinearSolver*>(p_ctx);
#if (DEBUG_CHECK_ASSERTIONS)
    TBOX_ASSERT(krylov_solver != NULL);
#endif
    if (right != PETSC_NULL)
    {
        // vector that the matrix can be multiplied against
        ierr = VecDuplicate(krylov_solver->d_petsc_x, right); IBTK_CHKERRQ(ierr);
        ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(*right)); IBTK_CHKERRQ(ierr);
    }
    if (left != PETSC_NULL)
    {
        // vector that the matrix vector product can be stored in
        ierr = VecDuplicate(krylov_solver->d_petsc_b, left); IBTK_CHKERRQ(ierr);
        ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(*left)); IBTK_CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}// MatGetVecs_SAMRAI

PetscErrorCode
PETScKrylovLinearSolver::PCApply_SAMRAI(
    PC pc,
    Vec x,
    Vec y)
{
    int ierr;
    void* ctx;
    ierr = PCShellGetContext(pc, &ctx); IBTK_CHKERRQ(ierr);
    PETScKrylovLinearSolver* krylov_solver = static_cast<PETScKrylovLinearSolver*>(ctx);
#if (DEBUG_CHECK_ASSERTIONS)
    TBOX_ASSERT(krylov_solver != NULL);
    TBOX_ASSERT(!krylov_solver->d_pc_solver.isNull());
#endif

    // Indicate that the initial guess should be zero.
    const bool pc_initial_guess_nonzero = krylov_solver->d_pc_solver->getInitialGuessNonzero();
    krylov_solver->d_pc_solver->setInitialGuessNonzero(false);

    // Apply the preconditioner.
    krylov_solver->d_pc_solver->solveSystem(*PETScSAMRAIVectorReal::getSAMRAIVector(y), *PETScSAMRAIVectorReal::getSAMRAIVector(x));
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y)); IBTK_CHKERRQ(ierr);

    // Reset the configuration of the preconditioner object.
    krylov_solver->d_pc_solver->setInitialGuessNonzero(pc_initial_guess_nonzero);
    PetscFunctionReturn(0);
}// PCApply_SAMRAI

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
