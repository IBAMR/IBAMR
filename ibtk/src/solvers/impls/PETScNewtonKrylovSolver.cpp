// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#include "ibtk/GeneralOperator.h"
#include "ibtk/GeneralSolver.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/JacobianOperator.h"
#include "ibtk/KrylovLinearSolver.h"
#include "ibtk/LinearOperator.h"
#include "ibtk/NewtonKrylovSolver.h"
#include "ibtk/PETScKrylovLinearSolver.h"
#include "ibtk/PETScNewtonKrylovSolver.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "ibtk/PETScSNESFunctionGOWrapper.h"
#include "ibtk/PETScSNESJacobianJOWrapper.h"
#include "ibtk/solver_utilities.h"

#include "Box.h"
#include "MultiblockDataTranslator.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"

#include "petscksp.h"
#include "petscmat.h"
#include "petscsnes.h"
#include "petscsys.h"
#include "petscvec.h"
#include "petscversion.h"

#include <mpi.h>

#include <ostream>
#include <string>
#include <utility>

#include "ibtk/namespaces.h" // IWYU pragma: keep

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
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

PETScNewtonKrylovSolver::PETScNewtonKrylovSolver(std::string object_name,
                                                 Pointer<Database> input_db,
                                                 std::string default_options_prefix,
                                                 MPI_Comm petsc_comm)
    : d_options_prefix(std::move(default_options_prefix)), d_petsc_comm(petsc_comm)
{
    // Setup default values.
    GeneralSolver::init(std::move(object_name), /*homogeneous_bc*/ false);
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
}

PETScNewtonKrylovSolver::PETScNewtonKrylovSolver(std::string object_name, SNES petsc_snes)
    : d_petsc_snes(std::move(petsc_snes))
{
    GeneralSolver::init(std::move(object_name), /*homogeneous_bc*/ false);
    if (d_petsc_snes) resetWrappedSNES(d_petsc_snes);
    common_ctor();
}

PETScNewtonKrylovSolver::~PETScNewtonKrylovSolver()
{
    if (d_is_initialized) deallocateSolverState();

    // Delete allocated PETSc solver components.
    int ierr;
    if (d_petsc_jac)
    {
        ierr = MatDestroy(&d_petsc_jac);
        IBTK_CHKERRQ(ierr);
        d_petsc_jac = nullptr;
    }
    if (d_managing_petsc_snes && d_petsc_snes)
    {
        ierr = SNESDestroy(&d_petsc_snes);
        IBTK_CHKERRQ(ierr);
        d_petsc_snes = nullptr;
    }
}

void
PETScNewtonKrylovSolver::setOptionsPrefix(const std::string& options_prefix)
{
    d_options_prefix = options_prefix;
}

const SNES&
PETScNewtonKrylovSolver::getPETScSNES() const
{
    return d_petsc_snes;
}

void
PETScNewtonKrylovSolver::setOperator(Pointer<GeneralOperator> F)
{
    NewtonKrylovSolver::setOperator(F);
    d_user_provided_function = true;
    resetSNESFunction();
}

Pointer<SAMRAIVectorReal<NDIM, double> >
PETScNewtonKrylovSolver::getSolutionVector() const
{
    Vec petsc_x;
    int ierr = SNESGetSolution(d_petsc_snes, &petsc_x);
    IBTK_CHKERRQ(ierr);
    Pointer<SAMRAIVectorReal<NDIM, double> > samrai_x;
    PETScSAMRAIVectorReal::getSAMRAIVectorRead(petsc_x, &samrai_x);
    Pointer<SAMRAIVectorReal<NDIM, double> > samrai_x_ptr(samrai_x, false);
    PETScSAMRAIVectorReal::restoreSAMRAIVectorRead(petsc_x, &samrai_x);
    return samrai_x_ptr;
}

Pointer<SAMRAIVectorReal<NDIM, double> >
PETScNewtonKrylovSolver::getFunctionVector() const
{
    Vec petsc_f;
    int ierr = SNESGetFunction(d_petsc_snes, &petsc_f, nullptr, nullptr);
    IBTK_CHKERRQ(ierr);
    Pointer<SAMRAIVectorReal<NDIM, double> > samrai_f;
    PETScSAMRAIVectorReal::getSAMRAIVectorRead(petsc_f, &samrai_f);
    Pointer<SAMRAIVectorReal<NDIM, double> > samrai_f_ptr(samrai_f, false);
    PETScSAMRAIVectorReal::restoreSAMRAIVectorRead(petsc_f, &samrai_f);
    return samrai_f_ptr;
}

void
PETScNewtonKrylovSolver::setJacobian(Pointer<JacobianOperator> J)
{
    NewtonKrylovSolver::setJacobian(J);
    d_user_provided_jacobian = true;
    resetSNESJacobian();
}

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
    ierr = SNESGetFunction(d_petsc_snes, &residual, nullptr, nullptr);
    IBTK_CHKERRQ(ierr);
    ierr = VecNorm(residual, NORM_2, &d_current_residual_norm);
    IBTK_CHKERRQ(ierr);

    // Determine the convergence reason.
    SNESConvergedReason reason;
    ierr = SNESGetConvergedReason(d_petsc_snes, &reason);
    IBTK_CHKERRQ(ierr);
    const bool converged = (static_cast<int>(reason) > 0);
    if (d_enable_logging) reportPETScSNESConvergedReason(d_object_name, reason, plog);

    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();

    IBTK_TIMER_STOP(t_solve_system);
    return converged;
}

void
PETScNewtonKrylovSolver::initializeSolverState(const SAMRAIVectorReal<NDIM, double>& x,
                                               const SAMRAIVectorReal<NDIM, double>& b)
{
    IBTK_TIMER_START(t_initialize_solver_state);

    int ierr;

#if !defined(NDEBUG)
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
        if (!patch_hierarchy->getPatchLevel(ln))
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

    // Allocate scratch data.
    d_b->allocateVectorData();
    d_r->allocateVectorData();

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
}

void
PETScNewtonKrylovSolver::deallocateSolverState()
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_deallocate_solver_state);

    // Deallocate the linear solver and operator states only if we are not
    // re-initializing the Newton solver.
    if (!d_reinitializing_solver)
    {
        if (d_krylov_solver) d_krylov_solver->deallocateSolverState();
        if (d_J) d_J->deallocateOperatorState();
        if (d_F) d_F->deallocateOperatorState();
    }

    // Deallocate scratch data.
    d_b->deallocateVectorData();
    d_r->deallocateVectorData();

    // Delete the solution and rhs vectors.
    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_x);
    d_petsc_x = nullptr;
    d_x->freeVectorComponents();
    d_x.setNull();

    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_b);
    d_petsc_b = nullptr;
    d_b->freeVectorComponents();
    d_b.setNull();

    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_r);
    d_petsc_r = nullptr;
    d_r->freeVectorComponents();
    d_r.setNull();

    // Destroy the SNES solver.
    if (d_managing_petsc_snes)
    {
        int ierr = SNESDestroy(&d_petsc_snes);
        IBTK_CHKERRQ(ierr);
        d_petsc_snes = nullptr;
    }

    // Indicate that the solver is NOT initialized.
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_solver_state);
}

/////////////////////////////// PRIVATE //////////////////////////////////////

void
PETScNewtonKrylovSolver::common_ctor()
{
    // Setup linear solver wrapper.
    KSP petsc_ksp = nullptr;
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
}

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
        ierr = SNESGetFunction(d_petsc_snes, nullptr, &petsc_snes_form_func, &petsc_snes_func_ctx);
        IBTK_CHKERRQ(ierr);
        d_F = new PETScSNESFunctionGOWrapper(
            d_object_name + "::SNESFunction Wrapper", d_petsc_snes, petsc_snes_form_func, petsc_snes_func_ctx);
        d_F->setHomogeneousBc(d_homogeneous_bc);
        d_F->setSolutionTime(d_solution_time);
        d_F->setTimeInterval(d_current_time, d_new_time);
    }

    if (d_user_provided_jacobian)
    {
        resetSNESJacobian();
    }
    else
    {
        // Create a JacobianOperator wrapper to correspond to the SNES Jacobian.
        PetscErrorCode (*petsc_snes_form_jac)(SNES, Vec, Mat, Mat, void*);
        void* petsc_snes_jac_ctx;
        ierr = SNESGetJacobian(d_petsc_snes, nullptr, nullptr, &petsc_snes_form_jac, &petsc_snes_jac_ctx);
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
}

void
PETScNewtonKrylovSolver::resetSNESOptions()
{
    if (!d_petsc_snes) return;
    int ierr = SNESSetTolerances(
        d_petsc_snes, d_abs_residual_tol, d_rel_residual_tol, d_solution_tol, d_max_iterations, d_max_evaluations);
    IBTK_CHKERRQ(ierr);
}

void
PETScNewtonKrylovSolver::resetSNESFunction()
{
    if (!d_petsc_snes) return;
    int ierr = SNESSetFunction(
        d_petsc_snes, d_petsc_r, PETScNewtonKrylovSolver::FormFunction_SAMRAI, static_cast<void*>(this));
    IBTK_CHKERRQ(ierr);
}

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
        d_petsc_jac = nullptr;
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
}

PetscErrorCode
PETScNewtonKrylovSolver::FormFunction_SAMRAI(SNES /*snes*/, Vec x, Vec f, void* p_ctx)
{
    auto newton_solver = static_cast<PETScNewtonKrylovSolver*>(p_ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(newton_solver);
    TBOX_ASSERT(newton_solver->d_F);
#endif
    Pointer<SAMRAIVectorReal<NDIM, double> > samrai_x, samrai_f;
    PETScSAMRAIVectorReal::getSAMRAIVectorRead(x, &samrai_x);
    PETScSAMRAIVectorReal::getSAMRAIVector(f, &samrai_f);
    newton_solver->d_F->apply(*samrai_x, *samrai_f);
    PETScSAMRAIVectorReal::restoreSAMRAIVectorRead(x, &samrai_x);
    PETScSAMRAIVectorReal::restoreSAMRAIVector(f, &samrai_f);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScNewtonKrylovSolver::FormJacobian_SAMRAI(SNES snes, Vec x, Mat A, Mat /*B*/, void* p_ctx)
{
    auto newton_solver = static_cast<PETScNewtonKrylovSolver*>(p_ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(newton_solver);
#endif
    if (newton_solver->d_J)
    {
        Pointer<SAMRAIVectorReal<NDIM, double> > samrai_x;
        PETScSAMRAIVectorReal::getSAMRAIVectorRead(x, &samrai_x);
        newton_solver->d_J->formJacobian(*samrai_x);
        PETScSAMRAIVectorReal::restoreSAMRAIVectorRead(x, &samrai_x);
    }
    else
    {
        Vec u, f;
        int ierr = SNESGetSolution(snes, &u);
        CHKERRQ(ierr);
        ierr = SNESGetFunction(snes, &f, nullptr, nullptr);
        CHKERRQ(ierr);
        ierr = MatMFFDSetBase(A, u, f);
        CHKERRQ(ierr);
        ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
        CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScNewtonKrylovSolver::MatVecMult_SAMRAI(Mat A, Vec x, Vec y)
{
    int ierr;
    void* p_ctx;
    ierr = MatShellGetContext(A, &p_ctx);
    CHKERRQ(ierr);
    auto newton_solver = static_cast<PETScNewtonKrylovSolver*>(p_ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(newton_solver);
    TBOX_ASSERT(newton_solver->d_J);
#endif
    Pointer<SAMRAIVectorReal<NDIM, double> > samrai_x, samrai_y;
    PETScSAMRAIVectorReal::getSAMRAIVectorRead(x, &samrai_x);
    PETScSAMRAIVectorReal::getSAMRAIVector(y, &samrai_y);
    newton_solver->d_J->apply(*samrai_x, *samrai_y);
    PETScSAMRAIVectorReal::restoreSAMRAIVectorRead(x, &samrai_x);
    PETScSAMRAIVectorReal::restoreSAMRAIVector(y, &samrai_y);
    PetscFunctionReturn(0);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
