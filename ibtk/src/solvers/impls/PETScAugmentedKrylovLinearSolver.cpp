// Filename: PETScAugmentedKrylovLinearSolver.cpp
// Created on 04 Apr 2016 by Boyce Griffith
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

#include <algorithm>
#include <math.h>
#include <ostream>
#include <stddef.h>
#include <string.h>
#include <string>
#include <vector>

#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "PatchHierarchy.h"
#include "SAMRAIVectorReal.h"
#include "ibtk/GeneralSolver.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/KrylovLinearSolver.h"
#include "ibtk/LinearOperator.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/PETScAugmentedKrylovLinearSolver.h"
#include "ibtk/PETScMatLOWrapper.h"
#include "ibtk/PETScMultiVec.h"
#include "ibtk/PETScPCLSWrapper.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "mpi.h"
#include "petscerror.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscoptions.h"
#include "petscpc.h"
#include "petscsys.h"
#include "petscvec.h"
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

PETScAugmentedKrylovLinearSolver::PETScAugmentedKrylovLinearSolver(const std::string& object_name,
                                                                   Pointer<Database> input_db,
                                                                   const std::string& default_options_prefix,
                                                                   MPI_Comm petsc_comm)
    : d_ksp_type(KSPGMRES),
      d_reinitializing_solver(false),
      d_petsc_x(NULL),
      d_petsc_b(NULL),
      d_options_prefix(default_options_prefix),
      d_petsc_comm(petsc_comm),
      d_petsc_ksp(NULL),
      d_petsc_mat(NULL)
{
    // Setup default values.
    GeneralSolver::init(object_name, /*homogeneous_bc*/ false);
    d_options_prefix = default_options_prefix;
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
        d_petsc_mat = NULL;
    }
    if (d_petsc_ksp)
    {
        ierr = KSPDestroy(&d_petsc_ksp);
        IBTK_CHKERRQ(ierr);
        d_petsc_ksp = NULL;
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

const KSP&
PETScAugmentedKrylovLinearSolver::getPETScKSP() const
{
    return d_petsc_ksp;
} // getPETScKSP

void
PETScAugmentedKrylovLinearSolver::setOperator(Pointer<LinearOperator> A)
{
    KrylovLinearSolver::setOperator(A);
    resetKSPOperators();
    return;
} // setOperator

void
PETScAugmentedKrylovLinearSolver::setPreconditioner(Pointer<LinearSolver> pc_solver)
{
    KrylovLinearSolver::setPreconditioner(pc_solver);
    resetKSPPC();
    return;
} // setPreconditioner

#if 0
void
PETScAugmentedKrylovLinearSolver::setNullspace(
    const bool contains_constant_vec,
    const std::vector<Pointer<SAMRAIVectorReal<NDIM, double> > >& nullspace_basis_vecs)
{
    deallocateNullspaceData();
    KrylovLinearSolver::setNullspace(contains_constant_vec, nullspace_basis_vecs);
    resetMatNullspace();
    return;
} // setNullspace
#endif

bool
PETScAugmentedKrylovLinearSolver::solveSystem(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& b)
{
    IBTK_TIMER_START(t_solve_system);

#if !defined(NDEBUG)
    TBOX_ASSERT(d_A);
#endif
    int ierr;

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x, b);
#if !defined(NDEBUG)
    TBOX_ASSERT(d_petsc_ksp);
#endif
    resetKSPOptions();

    // Allocate scratch data.
    d_b->allocateVectorData();

    // Solve the system using a PETSc KSP object.
    d_b->copyVector(Pointer<SAMRAIVectorReal<NDIM, double> >(&b, false));
    d_A->setHomogeneousBc(d_homogeneous_bc);
    d_A->modifyRhsForBcs(*d_b);
    d_A->setHomogeneousBc(true);
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_x, Pointer<SAMRAIVectorReal<NDIM, double> >(&x, false));
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_b, d_b);
    d_A->setHomogeneousBc(d_homogeneous_bc);
    d_A->imposeSolBcs(x);

    int n_vecs = 2;
    Vec x_vecs[] = { d_petsc_x, extra_x_vec };
    Vec multi_x;
    ierr = VecCreateMultiVec(d_petsc_comm, n_vecs, x_vecs, &multi_x);
    IBTK_CHKERRQ(ierr);
    Vec b_vecs[] = { d_petsc_b, extra_b_vec };
    Vec multi_b;
    ierr = VecCreateMultiVec(d_petsc_comm, n_vecs, b_vecs, &multi_b);
    IBTK_CHKERRQ(ierr);

    ierr = KSPSolve(d_petsc_ksp, multi_b, multi_x);
    IBTK_CHKERRQ(ierr);

    // Get iterations count and residual norm.
    ierr = KSPGetIterationNumber(d_petsc_ksp, &d_current_iterations);
    IBTK_CHKERRQ(ierr);
    ierr = KSPGetResidualNorm(d_petsc_ksp, &d_current_residual_norm);
    IBTK_CHKERRQ(ierr);
    d_A->setHomogeneousBc(d_homogeneous_bc);

    // Determine the convergence reason.
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(d_petsc_ksp, &reason);
    IBTK_CHKERRQ(ierr);
    const bool converged = (static_cast<int>(reason) > 0);
    if (d_enable_logging) reportKSPConvergedReason(reason, plog);

    ierr = VecDestroy(&multi_x);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&multi_b);
    IBTK_CHKERRQ(ierr);

    // Dealocate scratch data.
    d_b->deallocateVectorData();

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

    // Create the KSP solver.
    ierr = KSPCreate(d_petsc_comm, &d_petsc_ksp);
    IBTK_CHKERRQ(ierr);
    resetKSPOptions();

    // Setup solution and rhs vectors.
    d_x = x.cloneVector(x.getName());
    d_petsc_x = PETScSAMRAIVectorReal::createPETScVector(d_x, d_petsc_comm);

    d_b = b.cloneVector(b.getName());
    d_petsc_b = PETScSAMRAIVectorReal::createPETScVector(d_b, d_petsc_comm);

    // Initialize the linear operator and preconditioner objects.
    if (d_A) d_A->initializeOperatorState(*d_x, *d_b);
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
    ierr = KSPGetTolerances(d_petsc_ksp, &d_rel_residual_tol, &d_abs_residual_tol, NULL, &d_max_iterations);
    IBTK_CHKERRQ(ierr);

#if 0
    // Configure the nullspace object.
    resetMatNullspace();
#endif

    // Initialize augmented system data.
    initializeSolverStateAugmented(x, b);

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

    int ierr;

    // Deallocate the operator and preconditioner states only if we are not
    // re-initializing the solver.
    if (!d_reinitializing_solver)
    {
        if (d_pc_solver) d_pc_solver->deallocateSolverState();
        if (d_A) d_A->deallocateOperatorState();
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

#if 0
    // Deallocate the nullspace object.
    deallocateNullspaceData();
#endif

    // Destroy the KSP solver.
    ierr = KSPDestroy(&d_petsc_ksp);
    IBTK_CHKERRQ(ierr);
    d_petsc_ksp = NULL;

    // Indicate that the solver is NOT initialized.
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_solver_state);
    return;
} // deallocateSolverState

/////////////////////////////// PROTECTED ////////////////////////////////////

void
PETScAugmentedKrylovLinearSolver::initializeSolverStateAugmented(const SAMRAIVectorReal<NDIM, double>& x,
                                                                 const SAMRAIVectorReal<NDIM, double>& b)
{
    // intentionally blank
    return;
} // initializeSolverStateAugmented

void
PETScAugmentedKrylovLinearSolver::PCApplyAugmented(Vec x, Vec y)
{
    PetscErrorCode ierr;
    ierr = VecCopy(x, y);
    IBTK_CHKERRQ(ierr);
    return;
} // PCApplyAugmented

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
PETScAugmentedKrylovLinearSolver::reportKSPConvergedReason(const KSPConvergedReason& reason, std::ostream& os) const
{
    switch (static_cast<int>(reason))
    {
    case KSP_CONVERGED_RTOL:
        os << d_object_name << ": converged: |Ax-b| <= rtol*|b| --- residual norm is less than "
                               "specified relative tolerance.\n";
        break;
    case KSP_CONVERGED_ATOL:
        os << d_object_name << ": converged: |Ax-b| <= atol --- residual norm is less than "
                               "specified absolute tolerance.\n";
        break;
    case KSP_CONVERGED_ITS:
        os << d_object_name << ": converged: maximum number of iterations reached.\n";
        break;
    case KSP_CONVERGED_STEP_LENGTH:
        os << d_object_name << ": converged: step size less than specified tolerance.\n";
        break;
    case KSP_DIVERGED_NULL:
        os << d_object_name << ": diverged: null.\n";
        break;
    case KSP_DIVERGED_ITS:
        os << d_object_name << ": diverged: reached maximum number of iterations before any "
                               "convergence criteria were satisfied.\n";
        break;
    case KSP_DIVERGED_DTOL:
        os << d_object_name << ": diverged: |Ax-b| >= dtol*|b| --- residual is greater than "
                               "specified divergence tolerance.\n";
        break;
    case KSP_DIVERGED_BREAKDOWN:
        os << d_object_name << ": diverged: breakdown in the Krylov method.\n";
        break;
    case KSP_DIVERGED_BREAKDOWN_BICG:
        os << d_object_name << ": diverged: breakdown in the bi-congugate gradient method.\n";
        break;
    case KSP_DIVERGED_NONSYMMETRIC:
        os << d_object_name << ": diverged: it appears the operator or preconditioner is not "
                               "symmetric, but this Krylov method (KSPCG, KSPMINRES, KSPCR) "
                               "requires symmetry\n";
        break;
    case KSP_DIVERGED_INDEFINITE_PC:
        os << d_object_name << ": diverged: it appears the preconditioner is indefinite (has both "
                               "positive and negative eigenvalues), but this Krylov method (KSPCG) "
                               "requires it to be positive definite.\n";
        break;
    case KSP_CONVERGED_ITERATING:
        os << d_object_name << ": iterating: KSPSolve() is still running.\n";
        break;
    default:
        os << d_object_name << ": unknown completion code " << static_cast<int>(reason) << " reported.\n";
        break;
    }
    return;
} // reportKSPConvergedReason

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
        if (strcmp(mat_type, MATSHELL))
        {
            ierr = MatDestroy(&d_petsc_mat);
            IBTK_CHKERRQ(ierr);
            d_petsc_mat = NULL;
        }
    }
    if (!d_petsc_mat)
    {
        ierr = MatCreateShell(d_petsc_comm,
                              1 + (SAMRAI_MPI::getRank() == 0 ? 3 : 0),
                              1 + (SAMRAI_MPI::getRank() == 0 ? 3 : 0),
                              PETSC_DETERMINE,
                              PETSC_DETERMINE,
                              static_cast<void*>(this),
                              &d_petsc_mat); // XXXX dirty hack
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
    if (!d_petsc_ksp) return;
    int ierr;

    // Determine the preconditioner type to use.
    static const size_t len = 255;
    char pc_type_str[len];
    PetscBool flg;
#if (!PETSC_VERSION_RELEASE)
    ierr = PetscOptionsGetString(NULL, d_options_prefix.c_str(), "-pc_type", pc_type_str, len, &flg);
#else
    ierr = PetscOptionsGetString(d_options_prefix.c_str(), "-pc_type", pc_type_str, len, &flg);
#endif
    IBTK_CHKERRQ(ierr);
    std::string pc_type = "shell";
    if (flg)
    {
        pc_type = std::string(pc_type_str);
    }

    if (!(pc_type == "none" || pc_type == "shell"))
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  valid values for -"
                                 << d_options_prefix
                                 << "pc_type are: none, shell"
                                 << std::endl);
    }

    PC petsc_pc;
    ierr = KSPGetPC(d_petsc_ksp, &petsc_pc);
    IBTK_CHKERRQ(ierr);
    if (pc_type == "none" || !d_pc_solver)
    {
        ierr = PCSetType(petsc_pc, PCNONE);
        IBTK_CHKERRQ(ierr);
    }
    else if (pc_type == "shell" && d_pc_solver)
    {
        ierr = PCSetType(petsc_pc, PCSHELL);
        IBTK_CHKERRQ(ierr);
        ierr = PCShellSetContext(petsc_pc, static_cast<void*>(this));
        IBTK_CHKERRQ(ierr);
        ierr = PCShellSetApply(petsc_pc, PETScAugmentedKrylovLinearSolver::PCApply_SAMRAI);
        IBTK_CHKERRQ(ierr);
        ierr = PCShellSetName(petsc_pc, d_pc_solver->getName().c_str());
        IBTK_CHKERRQ(ierr);
    }
    else
    {
        TBOX_ERROR("this statement should not be reached!\n");
    }
    return;
} // resetKSPPC

PetscErrorCode
PETScAugmentedKrylovLinearSolver::MatVecMult_SAMRAI(Mat A, Vec x, Vec y)
{
    int ierr;
    void* p_ctx;
    ierr = MatShellGetContext(A, &p_ctx);
    IBTK_CHKERRQ(ierr);
    PETScAugmentedKrylovLinearSolver* krylov_solver = static_cast<PETScAugmentedKrylovLinearSolver*>(p_ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(krylov_solver);
    TBOX_ASSERT(krylov_solver->d_A);
#endif

    // Apply the augmented system.
    Vec* x_vecs;
    ierr = VecMultiVecGetSubVecs(x, &x_vecs);
    IBTK_CHKERRQ(ierr);
    Vec* y_vecs;
    ierr = VecMultiVecGetSubVecs(y, &y_vecs);
    IBTK_CHKERRQ(ierr);
    krylov_solver->d_A->apply(*PETScSAMRAIVectorReal::getSAMRAIVector(x_vecs[0]),
                              *PETScSAMRAIVectorReal::getSAMRAIVector(y_vecs[0]));
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y_vecs[0]));
    IBTK_CHKERRQ(ierr);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    IBTK_CHKERRQ(ierr);
    krylov_solver->MatVecMultAugmented(x, y);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    IBTK_CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // MatVecMult_SAMRAI

PetscErrorCode
PETScAugmentedKrylovLinearSolver::PCApply_SAMRAI(PC pc, Vec x, Vec y)
{
    int ierr;
    void* ctx;
    ierr = PCShellGetContext(pc, &ctx);
    IBTK_CHKERRQ(ierr);
    PETScAugmentedKrylovLinearSolver* krylov_solver = static_cast<PETScAugmentedKrylovLinearSolver*>(ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(krylov_solver);
    TBOX_ASSERT(krylov_solver->d_pc_solver);
#endif

    // Indicate that the initial guess should be zero.
    const bool pc_initial_guess_nonzero = krylov_solver->d_pc_solver->getInitialGuessNonzero();
    krylov_solver->d_pc_solver->setInitialGuessNonzero(false);

    // Apply the preconditioner.
    Vec* x_vecs;
    ierr = VecMultiVecGetSubVecs(x, &x_vecs);
    IBTK_CHKERRQ(ierr);
    Vec* y_vecs;
    ierr = VecMultiVecGetSubVecs(y, &y_vecs);
    IBTK_CHKERRQ(ierr);
    krylov_solver->d_pc_solver->solveSystem(*PETScSAMRAIVectorReal::getSAMRAIVector(y_vecs[0]),
                                            *PETScSAMRAIVectorReal::getSAMRAIVector(x_vecs[0]));
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y_vecs[0]));
    IBTK_CHKERRQ(ierr);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    IBTK_CHKERRQ(ierr);
    krylov_solver->PCApplyAugmented(x, y);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    IBTK_CHKERRQ(ierr);

    // Reset the configuration of the preconditioner object.
    krylov_solver->d_pc_solver->setInitialGuessNonzero(pc_initial_guess_nonzero);
    PetscFunctionReturn(0);
} // PCApply_SAMRAI

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
