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

#include "ibtk/GeneralSolver.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/KrylovLinearSolver.h"
#include "ibtk/LinearOperator.h"
#include "ibtk/PETScAugmentedKrylovLinearSolver.h"
#include "ibtk/PETScMatLOWrapper.h"
#include "ibtk/PETScPCLSWrapper.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "ibtk/solver_utilities.h"

#include "Box.h"
#include "MultiblockDataTranslator.h"
#include "PatchHierarchy.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/Utilities.h"

#include "petscksp.h"
#include "petscmat.h"
#include "petscpc.h"
#include "petscpctypes.h"
#include "petscvec.h"
#include <petsclog.h>

#include <mpi.h>

#include <algorithm>
#include <cstring>
#include <ostream>
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
// Timers.
static Timer* t_solve_system;
static Timer* t_initialize_solver_state;
static Timer* t_deallocate_solver_state;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

PETScAugmentedKrylovLinearSolver::PETScAugmentedKrylovLinearSolver(std::string object_name,
                                                                   Pointer<Database> input_db,
                                                                   std::string default_options_prefix,
                                                                   MPI_Comm petsc_comm)
    : d_ksp_type(KSPGMRES), d_options_prefix(std::move(default_options_prefix)), d_petsc_comm(petsc_comm)
{
    // Setup default values.
    GeneralSolver::init(std::move(object_name), /*homogeneous_bc*/ false);
    d_max_iterations = 10000;
    d_abs_residual_tol = 1.0e-50;
    d_rel_residual_tol = 1.0e-5;
    d_ksp_type = KSPGMRES;
    d_initial_guess_nonzero = true;
    d_enable_logging = false;

    // Get values from the input database.
    if (input_db)
    {
        if (input_db->keyExists("options_prefix")) d_options_prefix = input_db->getString("options_prefix");
        if (input_db->keyExists("max_iterations")) d_max_iterations = input_db->getInteger("max_iterations");
        if (input_db->keyExists("abs_residual_tol")) d_abs_residual_tol = input_db->getDouble("abs_residual_tol");
        if (input_db->keyExists("rel_residual_tol")) d_rel_residual_tol = input_db->getDouble("rel_residual_tol");
        if (input_db->keyExists("ksp_type")) d_ksp_type = input_db->getString("ksp_type");
        if (input_db->keyExists("initial_guess_nonzero"))
            d_initial_guess_nonzero = input_db->getBool("initial_guess_nonzero");
        if (input_db->keyExists("enable_logging")) d_enable_logging = input_db->getBool("enable_logging");
    }

    // Common constructor functionality.
    common_ctor();
    return;
} // PETScAugmentedKrylovLinearSolver()

PETScAugmentedKrylovLinearSolver::~PETScAugmentedKrylovLinearSolver()
{
    if (d_is_initialized) deallocateSolverState();

    // Delete allocated PETSc solver components.
    int ierr;
    if (d_petsc_mat)
    {
        ierr = MatDestroy(&d_petsc_mat);
        IBTK_CHKERRQ(ierr);
        d_petsc_mat = nullptr;
    }
    return;
} // ~PETScAugmentedKrylovLinearSolver()

void
PETScAugmentedKrylovLinearSolver::setKSPType(const std::string& ksp_type)
{
    d_ksp_type = ksp_type;
    return;
} // setKSPType

void
PETScAugmentedKrylovLinearSolver::setOptionsPrefix(const std::string& options_prefix)
{
    d_options_prefix = options_prefix;
    return;
} // setOptionsPrefix

void
PETScAugmentedKrylovLinearSolver::setAugmentedRHS(const Vec& vec)
{
    d_aug_b = vec;
}

void
PETScAugmentedKrylovLinearSolver::setInitialGuess(const Vec& vec)
{
    d_aug_vec = vec;
}

const KSP&
PETScAugmentedKrylovLinearSolver::getPETScKSP() const
{
    return d_petsc_ksp;
} // getPETScKSP

const Vec&
PETScAugmentedKrylovLinearSolver::getAugmentedVec() const
{
    return d_aug_x;
}

void
PETScAugmentedKrylovLinearSolver::setOperator(Pointer<LinearOperator> A)
{
#ifndef NDEBUG
    Pointer<PETScLinearAugmentedOperator> A_op = A;
    if (!A_op)
        TBOX_ERROR(d_object_name + "::setOperator: The linear operator must be a PETScLinearAugmentedOperator\n");
#endif
    KrylovLinearSolver::setOperator(A);
    return;
} // setOperator

void
PETScAugmentedKrylovLinearSolver::setPreconditioner(Pointer<LinearSolver> pc_solver)
{
    TBOX_WARNING(d_object_name + "::setPreconditioner: preconditioner not set up for augmented systems.\n");
    return;
} // setPreconditioner

void
PETScAugmentedKrylovLinearSolver::setNullspace(
    const bool /*contains_constant_vec*/,
    const std::vector<Pointer<SAMRAIVectorReal<NDIM, double> > >& /*nullspace_basis_vecs*/)
{
    TBOX_WARNING(d_object_name + "::setNullspace: nullspace not set up for augmented systems.\n");
    return;
} // setNullspace

bool
PETScAugmentedKrylovLinearSolver::solveSystem(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& b)
{
    IBTK_TIMER_START(t_solve_system);

    Pointer<PETScLinearAugmentedOperator> A_op = d_A;
#if !defined(NDEBUG)
    TBOX_ASSERT(d_A);
    TBOX_ASSERT(A_op);
#endif

    int ierr;

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x, b);
#if !defined(NDEBUG)
    TBOX_ASSERT(d_petsc_ksp);
#endif
    resetKSPOptions();

    // Solve the system using a PETSc KSP object.
    d_b->copyVector(Pointer<SAMRAIVectorReal<NDIM, double> >(&b, false));
    A_op->setHomogeneousBc(d_homogeneous_bc);
    A_op->setAugmentedRhsForBcs(d_aug_b);
    A_op->modifyRhsForBcs(*d_b);
    d_A->setHomogeneousBc(true);
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_eul_x, Pointer<SAMRAIVectorReal<NDIM, double> >(&x, false));
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_eul_b, d_b);
    // Create nested vec
    ierr = VecCopy(d_aug_vec, d_aug_x);
    IBTK_CHKERRQ(ierr);
    Vec petsc_b, petsc_x;
    std::array<Vec, 2> b_blocks = { d_eul_b, d_aug_b }, x_blocks = { d_eul_x, d_aug_x };
    ierr = VecCreateNest(d_petsc_comm, /*num nested blocks*/ 2, nullptr, b_blocks.data(), &petsc_b);
    IBTK_CHKERRQ(ierr);
    ierr = VecCreateNest(d_petsc_comm, /*num nested blocks*/ 2, nullptr, x_blocks.data(), &petsc_x);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSolve(d_petsc_ksp, petsc_b, petsc_x);
    IBTK_CHKERRQ(ierr);
    A_op->setHomogeneousBc(d_homogeneous_bc);
    A_op->imposeSolBcs(x);

    // Get iterations count and residual norm.
    ierr = KSPGetIterationNumber(d_petsc_ksp, &d_current_iterations);
    IBTK_CHKERRQ(ierr);
    ierr = KSPGetResidualNorm(d_petsc_ksp, &d_current_residual_norm);
    IBTK_CHKERRQ(ierr);
    A_op->setHomogeneousBc(d_homogeneous_bc);

    // Determine the convergence reason.
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(d_petsc_ksp, &reason);
    IBTK_CHKERRQ(ierr);
    const bool converged = (static_cast<int>(reason) > 0);
    reportPETScKSPConvergedReason(d_object_name, reason, plog);

    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();

    IBTK_TIMER_STOP(t_solve_system);
    return converged;
} // solveSystem

void
PETScAugmentedKrylovLinearSolver::initializeSolverState(const SAMRAIVectorReal<NDIM, double>& x,
                                                        const SAMRAIVectorReal<NDIM, double>& b)
{
    IBTK_TIMER_START(t_initialize_solver_state);

    int ierr;

// Rudimentary error checking.
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

    ierr = KSPCreate(d_petsc_comm, &d_petsc_ksp);
    IBTK_CHKERRQ(ierr);
    resetKSPOptions();

    // Set PC type to be none for KSP
    PC petsc_pc;
    ierr = KSPGetPC(d_petsc_ksp, &petsc_pc);
    IBTK_CHKERRQ(ierr);
    ierr = PCSetType(petsc_pc, PCNONE);
    IBTK_CHKERRQ(ierr);

    // Setup solution and rhs vectors.
    d_x = x.cloneVector(x.getName());
    d_eul_x = PETScSAMRAIVectorReal::createPETScVector(d_x, d_petsc_comm);

    d_b = b.cloneVector(b.getName());
    d_eul_b = PETScSAMRAIVectorReal::createPETScVector(d_b, d_petsc_comm);

    // Setup augmented solution
    if (d_aug_vec)
    {
        ierr = VecDuplicate(d_aug_vec, &d_aug_x);
        IBTK_CHKERRQ(ierr);
        ierr = VecCopy(d_aug_vec, d_aug_x);
        IBTK_CHKERRQ(ierr);
    }
    else
    {
        ierr = VecDuplicate(d_aug_b, &d_aug_x);
        IBTK_CHKERRQ(ierr);
        // Setup zero initial guess
        ierr = VecZeroEntries(d_aug_x);
        IBTK_CHKERRQ(ierr);
    }

    // Allocate scratch data.
    d_b->allocateVectorData();

    // Initialize the linear operator and preconditioner objects.
    // Before we initialize, we need to set up augmented vector
    Pointer<PETScLinearAugmentedOperator> A_aug = d_A;
    A_aug->setAugmentedVec(d_aug_x);
    d_A->initializeOperatorState(*d_x, *d_b);
    resetKSPOperators();

    if (d_pc_solver) d_pc_solver->initializeSolverState(*d_x, *d_b);
    resetKSPPC();

    // Set the KSP options from the PETSc options database.
    if (d_options_prefix != "")
    {
        ierr = KSPSetOptionsPrefix(d_petsc_ksp, d_options_prefix.c_str());
        IBTK_CHKERRQ(ierr);
    }
    ierr = KSPSetFromOptions(d_petsc_ksp);
    IBTK_CHKERRQ(ierr);

    // Reset the member state variables to correspond to the values used by the
    // KSP object.  (Command-line options always take precedence.)
    const char* ksp_type;
    ierr = KSPGetType(d_petsc_ksp, &ksp_type);
    IBTK_CHKERRQ(ierr);
    d_ksp_type = ksp_type;
    PetscBool initial_guess_nonzero;
    ierr = KSPGetInitialGuessNonzero(d_petsc_ksp, &initial_guess_nonzero);
    IBTK_CHKERRQ(ierr);
    d_initial_guess_nonzero = (initial_guess_nonzero == PETSC_TRUE);
    ierr = KSPGetTolerances(d_petsc_ksp, &d_rel_residual_tol, &d_abs_residual_tol, nullptr, &d_max_iterations);
    IBTK_CHKERRQ(ierr);

    // Configure the nullspace object.
    resetMatNullspace();

    // Indicate that the solver is initialized.
    d_reinitializing_solver = false;
    d_is_initialized = true;

    IBTK_TIMER_STOP(t_initialize_solver_state);
    return;
} // initializeSolverState

void
PETScAugmentedKrylovLinearSolver::deallocateSolverState()
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_deallocate_solver_state);

    // Deallocate the operator and preconditioner states only if we are not
    // re-initializing the solver.
    if (!d_reinitializing_solver)
    {
        if (d_pc_solver) d_pc_solver->deallocateSolverState();
        if (d_A) d_A->deallocateOperatorState();
    }

    // Dealocate scratch data.
    d_b->deallocateVectorData();

    // Delete the solution and rhs vectors.
    PETScSAMRAIVectorReal::destroyPETScVector(d_eul_x);
    d_eul_x = nullptr;
    d_x->freeVectorComponents();
    d_x.setNull();

    PETScSAMRAIVectorReal::destroyPETScVector(d_eul_b);
    d_eul_b = nullptr;
    d_b->freeVectorComponents();
    d_b.setNull();

    // Deallocate the nullspace object.
    deallocateNullspaceData();

    // Destroy the KSP solver.
    int ierr = KSPDestroy(&d_petsc_ksp);
    IBTK_CHKERRQ(ierr);
    d_petsc_ksp = nullptr;
    d_solver_has_attached_nullspace = false;

    // Indicate that the solver is NOT initialized.
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_solver_state);
    return;
} // deallocateSolverState

/////////////////////////////// PRIVATE //////////////////////////////////////

void
PETScAugmentedKrylovLinearSolver::common_ctor()
{
    // Setup Timers.
    IBTK_DO_ONCE(
        t_solve_system = TimerManager::getManager()->getTimer("IBTK::PETScAugmentedKrylovLinearSolver::solveSystem()");
        t_initialize_solver_state =
            TimerManager::getManager()->getTimer("IBTK::PETScAugmentedKrylovLinearSolver::initializeSolverState()");
        t_deallocate_solver_state =
            TimerManager::getManager()->getTimer("IBTK::PETScAugmentedKrylovLinearSolver::deallocateSolverState()"););
    return;
} // common_ctor

void
PETScAugmentedKrylovLinearSolver::resetKSPOptions()
{
    if (!d_petsc_ksp) return;
    int ierr;
    const KSPType ksp_type = d_ksp_type.c_str();
    ierr = KSPSetType(d_petsc_ksp, ksp_type);
    IBTK_CHKERRQ(ierr);
    std::string ksp_type_name(ksp_type);
    if (ksp_type_name.find("gmres") != std::string::npos)
    {
        ierr = KSPGMRESSetCGSRefinementType(d_petsc_ksp, KSP_GMRES_CGS_REFINE_IFNEEDED);
        IBTK_CHKERRQ(ierr);
    }
    PetscBool initial_guess_nonzero = (d_initial_guess_nonzero ? PETSC_TRUE : PETSC_FALSE);
    ierr = KSPSetInitialGuessNonzero(d_petsc_ksp, initial_guess_nonzero);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetTolerances(d_petsc_ksp, d_rel_residual_tol, d_abs_residual_tol, PETSC_DEFAULT, d_max_iterations);
    IBTK_CHKERRQ(ierr);
    return;
} // resetKSPOptions

void
PETScAugmentedKrylovLinearSolver::resetKSPOperators()
{
    int ierr;

    // Create and configure the MatShell object.
    if (d_petsc_mat)
    {
        const char* mat_type;
        ierr = MatGetType(d_petsc_mat, &mat_type);
        IBTK_CHKERRQ(ierr);
        if (std::strcmp(mat_type, MATSHELL))
        {
            ierr = MatDestroy(&d_petsc_mat);
            IBTK_CHKERRQ(ierr);
            d_petsc_mat = nullptr;
        }
    }
    if (!d_petsc_mat)
    {
        // Note we need the augmented vectors set at this point
        int aug_size;
        ierr = VecGetLocalSize(d_aug_b, &aug_size);
        IBTK_CHKERRQ(ierr);
        ierr = MatCreateShell(d_petsc_comm,
                              aug_size + 1,
                              aug_size + 1,
                              PETSC_DETERMINE,
                              PETSC_DETERMINE,
                              static_cast<void*>(this),
                              &d_petsc_mat);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatShellSetOperation(
        d_petsc_mat, MATOP_MULT, reinterpret_cast<void (*)(void)>(PETScAugmentedKrylovLinearSolver::MatVecMult_SAMRAI));
    IBTK_CHKERRQ(ierr);

    // Reset the configuration of the PETSc KSP object.
    if (d_petsc_ksp)
    {
        ierr = KSPSetOperators(d_petsc_ksp, d_petsc_mat, d_petsc_mat);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSetReusePreconditioner(d_petsc_ksp, PETSC_TRUE);
        IBTK_CHKERRQ(ierr);
    }
    return;
} // resetKSPOperators

void
PETScAugmentedKrylovLinearSolver::resetKSPPC()
{
    // intentionally blank
    return;
} // resetKSPPC

void
PETScAugmentedKrylovLinearSolver::resetMatNullspace()
{
    // intentionally blank
    return;
} // resetMatNullspace

void
PETScAugmentedKrylovLinearSolver::deallocateNullspaceData()
{
    // intentionally blank
    return;
} // deallocateNullspaceData

PetscErrorCode
PETScAugmentedKrylovLinearSolver::MatVecMult_SAMRAI(Mat A, Vec x, Vec y)
{
    void* p_ctx;
    int ierr = MatShellGetContext(A, &p_ctx);
    CHKERRQ(ierr);
    auto krylov_solver = static_cast<PETScAugmentedKrylovLinearSolver*>(p_ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(krylov_solver);
    TBOX_ASSERT(krylov_solver->d_A);
#endif
    Pointer<PETScLinearAugmentedOperator> A_op = krylov_solver->d_A;
#if !defined(NDEBUG)
    TBOX_ASSERT(A_op);
#endif
    // Pull out the two blocks
    // First element : Eulerian Vec
    // Second element: Augmented Vec

    std::array<Vec, 2> x_vecs, y_vecs;
    for (unsigned int i = 0; i < 2; ++i)
    {
        ierr = VecNestGetSubVec(x, i, &x_vecs[i]);
        CHKERRQ(ierr);
        ierr = VecNestGetSubVec(y, i, &y_vecs[i]);
        CHKERRQ(ierr);
    }

    Pointer<SAMRAIVectorReal<NDIM, double> > samrai_x, samrai_y;
    PETScSAMRAIVectorReal::getSAMRAIVectorRead(x_vecs[0], &samrai_x);
    PETScSAMRAIVectorReal::getSAMRAIVector(y_vecs[0], &samrai_y);
    A_op->setAugmentedVec(x_vecs[1]);
    A_op->apply(*samrai_x, *samrai_y);
    const Vec& temp_y_vec = A_op->getAugmentedVec();
    ierr = VecCopy(temp_y_vec, y_vecs[1]);
    IBTK_CHKERRQ(ierr);
    PETScSAMRAIVectorReal::restoreSAMRAIVectorRead(x_vecs[0], &samrai_x);
    PETScSAMRAIVectorReal::restoreSAMRAIVector(y_vecs[0], &samrai_y);

    PetscFunctionReturn(0);
} // MatVecMult_SAMRAI

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
