// Filename: PETScNewtonKrylovSolver.cpp
// Created on 26 Nov 2003 by Boyce Griffith
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

#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "PatchHierarchy.h"
#include "SAMRAIVectorReal.h"
#include "ibtk/GeneralOperator.h"
#include "ibtk/GeneralSolver.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/JacobianOperator.h"
#include "ibtk/KrylovLinearSolver.h"
#include "ibtk/LinearOperator.h"
#include "ibtk/NewtonKrylovSolver.h"
#include "ibtk/PETScKrylovLinearSolver.h"
#include "ibtk/PETScNewtonKrylovSolver.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "ibtk/PETScSNESFunctionGOWrapper.h"
#include "ibtk/PETScSNESJacobianJOWrapper.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "mpi.h"
#include "petscerror.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscsnes.h"
#include "petscsys.h"
#include "petscvec.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
// IWYU pragma: no_include "petsc-private/petscimpl.h"

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

PETScNewtonKrylovSolver::PETScNewtonKrylovSolver(const std::string& object_name,
                                                 Pointer<Database> input_db,
                                                 const std::string& default_options_prefix,
                                                 MPI_Comm petsc_comm)
    : d_reinitializing_solver(false),
      d_petsc_x(NULL),
      d_petsc_b(NULL),
      d_petsc_r(NULL),
      d_options_prefix(default_options_prefix),
      d_petsc_comm(petsc_comm),
      d_petsc_snes(NULL),
      d_petsc_jac(NULL),
      d_managing_petsc_snes(true),
      d_user_provided_function(false),
      d_user_provided_jacobian(false)
{
    // Setup default values.
    GeneralSolver::init(object_name, /*homogeneous_bc*/ false);
    d_options_prefix = default_options_prefix;
    d_max_iterations = 50;
    d_abs_residual_tol = 1.0e-50;
    d_rel_residual_tol = 1.0e-8;
    d_solution_tol = 1.0e-8;
    d_enable_logging = false;

    // Get values from the input database.
    if (input_db)
    {
        if (input_db->keyExists("options_prefix")) d_options_prefix = input_db->getString("options_prefix");
        if (input_db->keyExists("max_iterations")) d_max_iterations = input_db->getInteger("max_iterations");
        if (input_db->keyExists("abs_residual_tol")) d_abs_residual_tol = input_db->getDouble("abs_residual_tol");
        if (input_db->keyExists("rel_residual_tol")) d_rel_residual_tol = input_db->getDouble("rel_residual_tol");
        if (input_db->keyExists("solution_tol")) d_solution_tol = input_db->getDouble("solution_tol");
        if (input_db->keyExists("enable_logging")) d_enable_logging = input_db->getBool("enable_logging");
    }

    // Common constructor functionality.
    common_ctor();
    return;
} // PETScNewtonKrylovSolver()

PETScNewtonKrylovSolver::PETScNewtonKrylovSolver(const std::string& object_name, const SNES& petsc_snes)
    : d_reinitializing_solver(false),
      d_petsc_x(NULL),
      d_petsc_b(NULL),
      d_petsc_r(NULL),
      d_options_prefix(""),
      d_petsc_comm(PETSC_COMM_WORLD),
      d_petsc_snes(petsc_snes),
      d_petsc_jac(NULL),
      d_managing_petsc_snes(false),
      d_user_provided_function(false),
      d_user_provided_jacobian(false)
{
    GeneralSolver::init(object_name, /*homogeneous_bc*/ false);
    if (d_petsc_snes) resetWrappedSNES(d_petsc_snes);
    common_ctor();
    return;
} // PETScNewtonKrylovSolver()

PETScNewtonKrylovSolver::~PETScNewtonKrylovSolver()
{
    if (d_is_initialized) deallocateSolverState();

    // Delete allocated PETSc solver components.
    int ierr;
    if (d_petsc_jac)
    {
        ierr = MatDestroy(&d_petsc_jac);
        IBTK_CHKERRQ(ierr);
        d_petsc_jac = NULL;
    }
    if (d_managing_petsc_snes && d_petsc_snes)
    {
        ierr = SNESDestroy(&d_petsc_snes);
        IBTK_CHKERRQ(ierr);
        d_petsc_snes = NULL;
    }
    return;
} // ~PETScNewtonKrylovSolver()

void
PETScNewtonKrylovSolver::setOptionsPrefix(const std::string& options_prefix)
{
    d_options_prefix = options_prefix;
    return;
} // setOptionsPrefix

const SNES&
PETScNewtonKrylovSolver::getPETScSNES() const
{
    return d_petsc_snes;
} // getPETScSNES

void
PETScNewtonKrylovSolver::setOperator(Pointer<GeneralOperator> F)
{
    NewtonKrylovSolver::setOperator(F);
    d_user_provided_function = true;
    resetSNESFunction();
    return;
} // setOperator

Pointer<SAMRAIVectorReal<NDIM, double> >
PETScNewtonKrylovSolver::getSolutionVector() const
{
    Vec petsc_x;
    int ierr = SNESGetSolution(d_petsc_snes, &petsc_x);
    IBTK_CHKERRQ(ierr);
    return PETScSAMRAIVectorReal::getSAMRAIVector(petsc_x);
} // getSolutionVector

Pointer<SAMRAIVectorReal<NDIM, double> >
PETScNewtonKrylovSolver::getFunctionVector() const
{
    Vec petsc_f;
    int ierr = SNESGetFunction(d_petsc_snes, &petsc_f, NULL, NULL);
    IBTK_CHKERRQ(ierr);
    return PETScSAMRAIVectorReal::getSAMRAIVector(petsc_f);
} // getFunctionVector

void
PETScNewtonKrylovSolver::setJacobian(Pointer<JacobianOperator> J)
{
    NewtonKrylovSolver::setJacobian(J);
    d_user_provided_jacobian = true;
    resetSNESJacobian();
    return;
} // setJacobian

bool
PETScNewtonKrylovSolver::solveSystem(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& b)
{
    IBTK_TIMER_START(t_solve_system);

    int ierr;

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x, b);
#if !defined(NDEBUG)
    TBOX_ASSERT(d_petsc_snes);
#endif
    resetSNESOptions();
    Pointer<PETScKrylovLinearSolver> p_krylov_solver = d_krylov_solver;
    if (p_krylov_solver) p_krylov_solver->resetKSPOptions();

    // Allocate scratch data.
    d_b->allocateVectorData();
    d_r->allocateVectorData();

    // Solve the system using a PETSc SNES object.
    d_b->copyVector(Pointer<SAMRAIVectorReal<NDIM, double> >(&b, false));
    d_F->setHomogeneousBc(d_homogeneous_bc);
    d_F->modifyRhsForBcs(*d_b);
    Pointer<LinearOperator> A = d_F;
    if (A) A->setHomogeneousBc(true);
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_x, Pointer<SAMRAIVectorReal<NDIM, double> >(&x, false));
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_b, d_b);
    ierr = SNESSolve(d_petsc_snes, d_petsc_b, d_petsc_x);
    if (A) A->setHomogeneousBc(d_homogeneous_bc);
    d_F->imposeSolBcs(x);

    // Get iterations counts and residual norm.
    IBTK_CHKERRQ(ierr);
    ierr = SNESGetIterationNumber(d_petsc_snes, &d_current_iterations);
    IBTK_CHKERRQ(ierr);
    ierr = SNESGetLinearSolveIterations(d_petsc_snes, &d_current_linear_iterations);
    IBTK_CHKERRQ(ierr);
    Vec residual;
    ierr = SNESGetFunction(d_petsc_snes, &residual, NULL, NULL);
    IBTK_CHKERRQ(ierr);
    ierr = VecNorm(residual, NORM_2, &d_current_residual_norm);
    IBTK_CHKERRQ(ierr);

    // Determine the convergence reason.
    SNESConvergedReason reason;
    ierr = SNESGetConvergedReason(d_petsc_snes, &reason);
    IBTK_CHKERRQ(ierr);
    const bool converged = (static_cast<int>(reason) > 0);
    if (d_enable_logging) reportSNESConvergedReason(reason, plog);

    // Deallocate scratch data.
    d_b->deallocateVectorData();
    d_r->deallocateVectorData();

    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();

    IBTK_TIMER_STOP(t_solve_system);
    return converged;
} // solveSystem

void
PETScNewtonKrylovSolver::initializeSolverState(const SAMRAIVectorReal<NDIM, double>& x,
                                               const SAMRAIVectorReal<NDIM, double>& b)
{
    IBTK_TIMER_START(t_initialize_solver_state);

    int ierr;

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
        ierr = SNESCreate(d_petsc_comm, &d_petsc_snes);
        IBTK_CHKERRQ(ierr);
        resetSNESOptions();
    }
    else if (!d_petsc_snes)
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  cannot initialize solver state for wrapped PETSc SNES "
                                    "object if the wrapped object is NULL"
                                 << std::endl);
    }

    // Setup solution and rhs vectors.
    d_x = x.cloneVector(x.getName());
    d_petsc_x = PETScSAMRAIVectorReal::createPETScVector(d_x, d_petsc_comm);

    d_b = b.cloneVector(b.getName());
    d_petsc_b = PETScSAMRAIVectorReal::createPETScVector(d_b, d_petsc_comm);

    d_r = b.cloneVector(b.getName());
    d_petsc_r = PETScSAMRAIVectorReal::createPETScVector(d_r, d_petsc_comm);

    // Setup the nonlinear operator.
    if (d_F) d_F->initializeOperatorState(*d_x, *d_b);
    if (d_managing_petsc_snes || d_user_provided_function) resetSNESFunction();

    // Setup the Jacobian.
    if (d_J) d_J->initializeOperatorState(*d_x, *d_b);
    if (d_managing_petsc_snes || d_user_provided_jacobian) resetSNESJacobian();

    // Set the SNES options from the PETSc options database.
    if (d_options_prefix != "")
    {
        ierr = SNESSetOptionsPrefix(d_petsc_snes, d_options_prefix.c_str());
        IBTK_CHKERRQ(ierr);
    }
    ierr = SNESSetFromOptions(d_petsc_snes);
    IBTK_CHKERRQ(ierr);

    // Reset the member state variables to correspond to the values used by the
    // SNES object.  (Command-line options always take precedence.)
    ierr = SNESGetTolerances(
        d_petsc_snes, &d_abs_residual_tol, &d_rel_residual_tol, &d_solution_tol, &d_max_iterations, &d_max_evaluations);
    IBTK_CHKERRQ(ierr);

    // Setup the KrylovLinearSolver wrapper to correspond to the KSP employed by
    // the SNES solver.
    KSP petsc_ksp;
    ierr = SNESGetKSP(d_petsc_snes, &petsc_ksp);
    IBTK_CHKERRQ(ierr);
    Pointer<PETScKrylovLinearSolver> p_krylov_solver = d_krylov_solver;
    if (p_krylov_solver) p_krylov_solver->resetWrappedKSP(petsc_ksp);

    // Setup the Krylov solver.
    if (d_krylov_solver) d_krylov_solver->initializeSolverState(*d_x, *d_b);

    // Indicate that the solver is initialized.
    d_reinitializing_solver = false;
    d_is_initialized = true;

    IBTK_TIMER_STOP(t_initialize_solver_state);
    return;
} // initializeSolverState

void
PETScNewtonKrylovSolver::deallocateSolverState()
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_deallocate_solver_state);

    int ierr;

    // Deallocate the linear solver and operator states only if we are not
    // re-initializing the Newton solver.
    if (!d_reinitializing_solver)
    {
        if (d_krylov_solver) d_krylov_solver->deallocateSolverState();
        if (d_J) d_J->deallocateOperatorState();
        if (d_F) d_F->deallocateOperatorState();
    }

    // Delete the solution and rhs vectors.
    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_x);
    d_petsc_x = NULL;
    d_x->freeVectorComponents();
    d_x.setNull();

    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_b);
    d_petsc_b = NULL;
    d_b->freeVectorComponents();
    d_b.setNull();

    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_r);
    d_petsc_r = NULL;
    d_r->freeVectorComponents();
    d_r.setNull();

    // Destroy the SNES solver.
    if (d_managing_petsc_snes)
    {
        ierr = SNESDestroy(&d_petsc_snes);
        IBTK_CHKERRQ(ierr);
        d_petsc_snes = NULL;
    }

    // Indicate that the solver is NOT initialized.
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_solver_state);
    return;
} // deallocateSolverState

/////////////////////////////// PRIVATE //////////////////////////////////////

void
PETScNewtonKrylovSolver::common_ctor()
{
    // Setup linear solver wrapper.
    KSP petsc_ksp = NULL;
    d_krylov_solver = new PETScKrylovLinearSolver(d_object_name + "::KSP Wrapper", petsc_ksp);
    d_krylov_solver->setHomogeneousBc(d_homogeneous_bc);
    d_krylov_solver->setSolutionTime(d_solution_time);
    d_krylov_solver->setTimeInterval(d_current_time, d_new_time);

    // Setup Timers.
    IBTK_DO_ONCE(t_solve_system = TimerManager::getManager()->getTimer("IBTK::PETScNewtonKrylovSolver::solveSystem()");
                 t_initialize_solver_state =
                     TimerManager::getManager()->getTimer("IBTK::PETScNewtonKrylovSolver::initializeOperatorState()");
                 t_deallocate_solver_state =
                     TimerManager::getManager()->getTimer("IBTK::PETScNewtonKrylovSolver::deallocateOperatorState()"););
    return;
} // common_ctor

void
PETScNewtonKrylovSolver::reportSNESConvergedReason(const SNESConvergedReason& reason, std::ostream& os) const
{
    switch (static_cast<int>(reason))
    {
    case SNES_CONVERGED_FNORM_ABS:
        os << d_object_name << ": converged: |F| less than specified absolute tolerance.\n";
        break;
    case SNES_CONVERGED_FNORM_RELATIVE:
        os << d_object_name << ": converged: |F| less than specified relative tolerance.\n";
        break;
    case SNES_CONVERGED_SNORM_RELATIVE:
        os << d_object_name << ": converged: step size less than specified relative tolerance.\n";
        break;
    case SNES_CONVERGED_ITS:
        os << d_object_name << ": converged: maximum number of iterations reached.\n";
        break;
    case SNES_CONVERGED_TR_DELTA:
        os << d_object_name << ": converged: trust-region delta.\n";
        break;
    case SNES_DIVERGED_FUNCTION_DOMAIN:
        os << d_object_name << ": diverged: new x location passed to the function is not in "
                               "the function domain.\n";
        break;
    case SNES_DIVERGED_FUNCTION_COUNT:
        os << d_object_name << ": diverged: exceeded maximum number of function evaluations.\n";
        break;
    case SNES_DIVERGED_LINEAR_SOLVE:
        os << d_object_name << ": diverged: the linear solve failed.\n";
        break;
    case SNES_DIVERGED_FNORM_NAN:
        os << d_object_name << ": diverged: |F| is NaN.\n";
        break;
    case SNES_DIVERGED_MAX_IT:
        os << d_object_name << ": diverged: exceeded maximum number of iterations.\n";
        break;
    case SNES_DIVERGED_LINE_SEARCH:
        os << d_object_name << ": diverged: line-search failure.\n";
        break;
    case SNES_DIVERGED_INNER:
        os << d_object_name << ": diverged: inner solve failed.\n";
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
} // reportSNESConvergedReason

void
PETScNewtonKrylovSolver::resetWrappedSNES(SNES& petsc_snes)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_managing_petsc_snes);
#endif
    d_petsc_snes = petsc_snes;
    if (!d_petsc_snes) return;

    int ierr;

    // Set d_petsc_comm to be the MPI communicator used by the supplied SNES.
    ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(d_petsc_snes), &d_petsc_comm);
    IBTK_CHKERRQ(ierr);

    // Set d_options_prefix to correspond to that used by the supplied SNES.
    const char* options_prefix;
    ierr = SNESGetOptionsPrefix(d_petsc_snes, &options_prefix);
    IBTK_CHKERRQ(ierr);
    d_options_prefix = options_prefix;

    // Setup the function and Jacobian.
    if (d_user_provided_function)
        resetSNESFunction();
    else
    {
        // Create an GeneralOperator wrapper to correspond to the SNES function.
        PetscErrorCode (*petsc_snes_form_func)(SNES, Vec, Vec, void*);
        void* petsc_snes_func_ctx;
        ierr = SNESGetFunction(d_petsc_snes, NULL, &petsc_snes_form_func, &petsc_snes_func_ctx);
        IBTK_CHKERRQ(ierr);
        d_F = new PETScSNESFunctionGOWrapper(
            d_object_name + "::SNESFunction Wrapper", d_petsc_snes, petsc_snes_form_func, petsc_snes_func_ctx);
        d_F->setHomogeneousBc(d_homogeneous_bc);
        d_F->setSolutionTime(d_solution_time);
        d_F->setTimeInterval(d_current_time, d_new_time);
    }

    if (d_user_provided_jacobian)
        resetSNESJacobian();
    else
    {
        // Create a JacobianOperator wrapper to correspond to the SNES Jacobian.
        PetscErrorCode (*petsc_snes_form_jac)(SNES, Vec, Mat, Mat, void*);
        void* petsc_snes_jac_ctx;
        ierr = SNESGetJacobian(d_petsc_snes, NULL, NULL, &petsc_snes_form_jac, &petsc_snes_jac_ctx);
        IBTK_CHKERRQ(ierr);
        d_J = new PETScSNESJacobianJOWrapper(
            d_object_name + "::SNESJacobian Wrapper", d_petsc_snes, petsc_snes_form_jac, petsc_snes_jac_ctx);
        d_J->setHomogeneousBc(true);
        d_J->setSolutionTime(d_solution_time);
        d_J->setTimeInterval(d_current_time, d_new_time);
    }

    // Setup the KrylovLinearSolver wrapper to use the KSP associated with the
    // wrapped SNES object.
    KSP petsc_ksp;
    ierr = SNESGetKSP(d_petsc_snes, &petsc_ksp);
    IBTK_CHKERRQ(ierr);
    Pointer<PETScKrylovLinearSolver> p_krylov_solver = d_krylov_solver;
    if (p_krylov_solver) p_krylov_solver->resetWrappedKSP(petsc_ksp);

    // Reset the member state variables to correspond to the values used by the
    // SNES object.
    ierr = SNESGetTolerances(
        d_petsc_snes, &d_abs_residual_tol, &d_rel_residual_tol, &d_solution_tol, &d_max_iterations, &d_max_evaluations);
    IBTK_CHKERRQ(ierr);
    return;
} // resetWrappedSNES

void
PETScNewtonKrylovSolver::resetSNESOptions()
{
    if (!d_petsc_snes) return;
    int ierr = SNESSetTolerances(
        d_petsc_snes, d_abs_residual_tol, d_rel_residual_tol, d_solution_tol, d_max_iterations, d_max_evaluations);
    IBTK_CHKERRQ(ierr);
    return;
} // resetSNESSNESOptions

void
PETScNewtonKrylovSolver::resetSNESFunction()
{
    if (!d_petsc_snes) return;
    int ierr = SNESSetFunction(
        d_petsc_snes, d_petsc_r, PETScNewtonKrylovSolver::FormFunction_SAMRAI, static_cast<void*>(this));
    IBTK_CHKERRQ(ierr);
    return;
} // resetSNESFunction

void
PETScNewtonKrylovSolver::resetSNESJacobian()
{
    if (!d_petsc_snes) return;
    int ierr;

    // Create and configure the Jacobian matrix.
    if (d_petsc_jac)
    {
        ierr = MatDestroy(&d_petsc_jac);
        IBTK_CHKERRQ(ierr);
        d_petsc_jac = NULL;
    }
    if (d_J && d_user_provided_jacobian)
    {
        ierr = MatCreateShell(
            d_petsc_comm, 1, 1, PETSC_DETERMINE, PETSC_DETERMINE, static_cast<void*>(this), &d_petsc_jac);
        IBTK_CHKERRQ(ierr);
        ierr = MatShellSetOperation(
            d_petsc_jac, MATOP_MULT, reinterpret_cast<void (*)(void)>(PETScNewtonKrylovSolver::MatVecMult_SAMRAI));
        IBTK_CHKERRQ(ierr);
    }
    else
    {
        ierr = MatCreateMFFD(d_petsc_comm, 1, 1, PETSC_DETERMINE, PETSC_DETERMINE, &d_petsc_jac);
        IBTK_CHKERRQ(ierr);
        ierr = MatMFFDSetFunction(
            d_petsc_jac, reinterpret_cast<PetscErrorCode (*)(void*, Vec, Vec)>(SNESComputeFunction), d_petsc_snes);
        IBTK_CHKERRQ(ierr);
        if (!d_options_prefix.empty())
        {
            ierr = MatSetOptionsPrefix(d_petsc_jac, d_options_prefix.c_str());
            IBTK_CHKERRQ(ierr);
        }
        ierr = MatSetFromOptions(d_petsc_jac);
        IBTK_CHKERRQ(ierr);
    }

    // Reset the configuration of the PETSc SNES object.
    ierr = SNESSetJacobian(
        d_petsc_snes, d_petsc_jac, d_petsc_jac, PETScNewtonKrylovSolver::FormJacobian_SAMRAI, static_cast<void*>(this));
    IBTK_CHKERRQ(ierr);
    return;
} // resetSNESJacobian

PetscErrorCode
PETScNewtonKrylovSolver::FormFunction_SAMRAI(SNES /*snes*/, Vec x, Vec f, void* p_ctx)
{
    int ierr;
    PETScNewtonKrylovSolver* newton_solver = static_cast<PETScNewtonKrylovSolver*>(p_ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(newton_solver);
    TBOX_ASSERT(newton_solver->d_F);
#endif
    newton_solver->d_F->apply(*PETScSAMRAIVectorReal::getSAMRAIVector(x), *PETScSAMRAIVectorReal::getSAMRAIVector(f));
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(f));
    IBTK_CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // FormFunction_SAMRAI

PetscErrorCode
PETScNewtonKrylovSolver::FormJacobian_SAMRAI(SNES snes, Vec x, Mat A, Mat /*B*/, void* p_ctx)
{
    int ierr;
    PETScNewtonKrylovSolver* newton_solver = static_cast<PETScNewtonKrylovSolver*>(p_ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(newton_solver);
#endif
    if (newton_solver->d_J)
    {
        newton_solver->d_J->formJacobian(*PETScSAMRAIVectorReal::getSAMRAIVector(x));
    }
    else
    {
        Vec u, f;
        ierr = SNESGetSolution(snes, &u);
        IBTK_CHKERRQ(ierr);
        ierr = SNESGetFunction(snes, &f, NULL, NULL);
        IBTK_CHKERRQ(ierr);
        ierr = MatMFFDSetBase(A, u, f);
        IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} // FormJacobian_SAMRAI

PetscErrorCode
PETScNewtonKrylovSolver::MatVecMult_SAMRAI(Mat A, Vec x, Vec y)
{
    int ierr;
    void* p_ctx;
    ierr = MatShellGetContext(A, &p_ctx);
    IBTK_CHKERRQ(ierr);
    PETScNewtonKrylovSolver* newton_solver = static_cast<PETScNewtonKrylovSolver*>(p_ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(newton_solver);
    TBOX_ASSERT(newton_solver->d_J);
#endif
    newton_solver->d_J->apply(*PETScSAMRAIVectorReal::getSAMRAIVector(x), *PETScSAMRAIVectorReal::getSAMRAIVector(y));
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    IBTK_CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // MatVecMult_SAMRAI

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
