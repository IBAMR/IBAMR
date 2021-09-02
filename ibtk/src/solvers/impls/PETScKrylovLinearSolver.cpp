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
#include "ibtk/PETScKrylovLinearSolver.h"
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

PETScKrylovLinearSolver::PETScKrylovLinearSolver(std::string object_name,
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
} // PETScKrylovLinearSolver()

PETScKrylovLinearSolver::PETScKrylovLinearSolver(std::string object_name, const KSP& petsc_ksp)
    : d_petsc_comm(PETSC_COMM_WORLD), d_petsc_ksp(petsc_ksp), d_managing_petsc_ksp(false)
{
    GeneralSolver::init(std::move(object_name), /*homogeneous_bc*/ false);
    if (d_petsc_ksp) resetWrappedKSP(d_petsc_ksp);
    common_ctor();
    return;
} // PETScKrylovLinearSolver()

PETScKrylovLinearSolver::~PETScKrylovLinearSolver()
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
    if (d_managing_petsc_ksp && d_petsc_ksp)
    {
        ierr = KSPDestroy(&d_petsc_ksp);
        IBTK_CHKERRQ(ierr);
        d_petsc_ksp = nullptr;
    }
    return;
} // ~PETScKrylovLinearSolver()

void
PETScKrylovLinearSolver::setKSPType(const std::string& ksp_type)
{
    d_ksp_type = ksp_type;
    return;
} // setKSPType

void
PETScKrylovLinearSolver::setOptionsPrefix(const std::string& options_prefix)
{
    d_options_prefix = options_prefix;
    return;
} // setOptionsPrefix

const KSP&
PETScKrylovLinearSolver::getPETScKSP() const
{
    return d_petsc_ksp;
} // getPETScKSP

void
PETScKrylovLinearSolver::setOperator(Pointer<LinearOperator> A)
{
    KrylovLinearSolver::setOperator(A);
    d_user_provided_mat = true;
    resetKSPOperators();
    return;
} // setOperator

void
PETScKrylovLinearSolver::setPreconditioner(Pointer<LinearSolver> pc_solver)
{
    KrylovLinearSolver::setPreconditioner(pc_solver);
    d_user_provided_pc = true;
    resetKSPPC();
    return;
} // setPreconditioner

void
PETScKrylovLinearSolver::setNullspace(
    const bool contains_constant_vec,
    const std::vector<Pointer<SAMRAIVectorReal<NDIM, double> > >& nullspace_basis_vecs)
{
    deallocateNullspaceData();
    KrylovLinearSolver::setNullspace(contains_constant_vec, nullspace_basis_vecs);
    resetMatNullspace();
    return;
} // setNullspace

bool
PETScKrylovLinearSolver::solveSystem(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& b)
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

    // Solve the system using a PETSc KSP object.
    d_b->copyVector(Pointer<SAMRAIVectorReal<NDIM, double> >(&b, false));
    d_A->setHomogeneousBc(d_homogeneous_bc);
    d_A->modifyRhsForBcs(*d_b);
    d_A->setHomogeneousBc(true);
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_x, Pointer<SAMRAIVectorReal<NDIM, double> >(&x, false));
    PETScSAMRAIVectorReal::replaceSAMRAIVector(d_petsc_b, d_b);
    ierr = KSPSolve(d_petsc_ksp, d_petsc_b, d_petsc_x);
    IBTK_CHKERRQ(ierr);
    d_A->setHomogeneousBc(d_homogeneous_bc);
    d_A->imposeSolBcs(x);

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
    if (d_enable_logging) reportPETScKSPConvergedReason(d_object_name, reason, plog);

    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();

    IBTK_TIMER_STOP(t_solve_system);
    return converged;
} // solveSystem

void
PETScKrylovLinearSolver::initializeSolverState(const SAMRAIVectorReal<NDIM, double>& x,
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

    // Create the KSP solver.
    if (d_managing_petsc_ksp)
    {
        ierr = KSPCreate(d_petsc_comm, &d_petsc_ksp);
        IBTK_CHKERRQ(ierr);
        resetKSPOptions();
    }
    else if (!d_petsc_ksp)
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  cannot initialize solver state for wrapped PETSc KSP object "
                                    "if the wrapped object is NULL"
                                 << std::endl);
    }

    // Setup solution and rhs vectors.
    d_x = x.cloneVector(x.getName());
    d_petsc_x = PETScSAMRAIVectorReal::createPETScVector(d_x, d_petsc_comm);

    d_b = b.cloneVector(b.getName());
    d_petsc_b = PETScSAMRAIVectorReal::createPETScVector(d_b, d_petsc_comm);

    // Allocate scratch data.
    d_b->allocateVectorData();

    // Initialize the linear operator and preconditioner objects.
    if (d_A) d_A->initializeOperatorState(*d_x, *d_b);
    if (d_managing_petsc_ksp || d_user_provided_mat) resetKSPOperators();

    if (d_pc_solver) d_pc_solver->initializeSolverState(*d_x, *d_b);
    if (d_managing_petsc_ksp || d_user_provided_pc) resetKSPPC();

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
PETScKrylovLinearSolver::deallocateSolverState()
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
    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_x);
    d_petsc_x = nullptr;
    d_x->freeVectorComponents();
    d_x.setNull();

    PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_b);
    d_petsc_b = nullptr;
    d_b->freeVectorComponents();
    d_b.setNull();

    // Deallocate the nullspace object.
    deallocateNullspaceData();

    // Destroy the KSP solver.
    if (d_managing_petsc_ksp)
    {
        int ierr = KSPDestroy(&d_petsc_ksp);
        IBTK_CHKERRQ(ierr);
        d_petsc_ksp = nullptr;
        d_solver_has_attached_nullspace = false;
    }

    // Indicate that the solver is NOT initialized.
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_solver_state);
    return;
} // deallocateSolverState

/////////////////////////////// PRIVATE //////////////////////////////////////

void
PETScKrylovLinearSolver::common_ctor()
{
    // Setup Timers.
    IBTK_DO_ONCE(t_solve_system = TimerManager::getManager()->getTimer("IBTK::PETScKrylovLinearSolver::solveSystem()");
                 t_initialize_solver_state =
                     TimerManager::getManager()->getTimer("IBTK::PETScKrylovLinearSolver::initializeSolverState()");
                 t_deallocate_solver_state =
                     TimerManager::getManager()->getTimer("IBTK::PETScKrylovLinearSolver::deallocateSolverState()"););
    return;
} // common_ctor

void
PETScKrylovLinearSolver::resetWrappedKSP(KSP& petsc_ksp)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_managing_petsc_ksp);
#endif
    d_petsc_ksp = petsc_ksp;
    if (!d_petsc_ksp) return;
    int ierr;

    // Set d_petsc_comm to be the MPI communicator used by the supplied KSP.
    ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(d_petsc_ksp), &d_petsc_comm);
    IBTK_CHKERRQ(ierr);

    // Set d_ksp_type to correspond to the KSP type used by the supplied KSP.
    const char* ksp_type;
    ierr = KSPGetType(d_petsc_ksp, &ksp_type);
    IBTK_CHKERRQ(ierr);
    d_ksp_type = std::string(ksp_type);

    // Set d_options_prefix to correspond to that used by the supplied KSP.
    const char* options_prefix;
    ierr = KSPGetOptionsPrefix(d_petsc_ksp, &options_prefix);
    IBTK_CHKERRQ(ierr);
    d_options_prefix = options_prefix;

    // Setup operators and preconditioners.
    if (d_user_provided_mat)
    {
        resetKSPOperators();
    }
    else
    {
        // Create a LinearOperator wrapper to correspond to the PETSc Mat used
        // by the KSP.
        Mat petsc_mat;
        ierr = KSPGetOperators(d_petsc_ksp, &petsc_mat, nullptr);
        IBTK_CHKERRQ(ierr);
        d_A = new PETScMatLOWrapper(d_object_name + "::Mat Wrapper", petsc_mat);
        d_A->setHomogeneousBc(d_homogeneous_bc);
        d_A->setSolutionTime(d_solution_time);
        d_A->setTimeInterval(d_current_time, d_new_time);
    }

    if (d_user_provided_pc)
    {
        resetKSPPC();
    }
    else
    {
        // Create a LinearSolver wrapper to correspond to the PETSc PC used by
        // the KSP.
        PC petsc_pc;
        ierr = KSPGetPC(d_petsc_ksp, &petsc_pc);
        IBTK_CHKERRQ(ierr);
        d_pc_solver = new PETScPCLSWrapper(d_object_name + "::PC Wrapper", petsc_pc);
        d_pc_solver->setHomogeneousBc(true);
        d_pc_solver->setSolutionTime(d_solution_time);
        d_pc_solver->setTimeInterval(d_current_time, d_new_time);
    }

    // Reset the member state variables to correspond to the values used by the
    // KSP object.
    PetscBool initial_guess_nonzero;
    ierr = KSPGetInitialGuessNonzero(d_petsc_ksp, &initial_guess_nonzero);
    IBTK_CHKERRQ(ierr);
    d_initial_guess_nonzero = (initial_guess_nonzero == PETSC_TRUE);
    ierr = KSPGetTolerances(d_petsc_ksp, &d_rel_residual_tol, &d_abs_residual_tol, nullptr, &d_max_iterations);
    IBTK_CHKERRQ(ierr);
    return;
} // resetWrappedKSP

void
PETScKrylovLinearSolver::resetKSPOptions()
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
PETScKrylovLinearSolver::resetKSPOperators()
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
        ierr = MatCreateShell(
            d_petsc_comm, 1, 1, PETSC_DETERMINE, PETSC_DETERMINE, static_cast<void*>(this), &d_petsc_mat);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatShellSetOperation(
        d_petsc_mat, MATOP_MULT, reinterpret_cast<void (*)(void)>(PETScKrylovLinearSolver::MatVecMult_SAMRAI));
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
PETScKrylovLinearSolver::resetKSPPC()
{
    if (!d_petsc_ksp) return;
    int ierr;

    // Determine the preconditioner type to use.
    static const size_t len = 255;
    char pc_type_str[len];
    PetscBool flg;
    ierr = PetscOptionsGetString(nullptr, d_options_prefix.c_str(), "-pc_type", pc_type_str, len, &flg);
    IBTK_CHKERRQ(ierr);
    std::string pc_type = "shell";
    if (flg)
    {
        pc_type = std::string(pc_type_str);
    }

    if (!(pc_type == "none" || pc_type == "shell"))
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  valid values for -" << d_options_prefix << "pc_type are: none, shell"
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
        ierr = PCShellSetApply(petsc_pc, PETScKrylovLinearSolver::PCApply_SAMRAI);
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

void
PETScKrylovLinearSolver::resetMatNullspace()
{
    if (!d_petsc_ksp) return;
    int ierr;
    PetscBool flg;
    ierr = PetscOptionsHasName(nullptr, d_options_prefix.c_str(), "-ksp_constant_null_space", &flg);
    IBTK_CHKERRQ(ierr);
    if (flg == PETSC_TRUE) d_nullspace_contains_constant_vec = true;
    if (d_nullspace_contains_constant_vec || !d_nullspace_basis_vecs.empty())
    {
        std::vector<Vec> nullspace_vecs;
        nullspace_vecs.reserve(1 + d_nullspace_basis_vecs.size());
        if (d_nullspace_contains_constant_vec)
        {
            d_nullspace_constant_vec = d_x->cloneVector(d_x->getName());
            d_nullspace_constant_vec->allocateVectorData();
            d_petsc_nullspace_constant_vec =
                PETScSAMRAIVectorReal::createPETScVector(d_nullspace_constant_vec, d_petsc_comm);
            ierr = VecSet(d_petsc_nullspace_constant_vec, 1.0);
            IBTK_CHKERRQ(ierr);
            nullspace_vecs.push_back(d_petsc_nullspace_constant_vec);
        }

        d_petsc_nullspace_basis_vecs.resize(d_nullspace_basis_vecs.size());
        for (unsigned int k = 0; k < d_nullspace_basis_vecs.size(); ++k)
        {
            d_petsc_nullspace_basis_vecs[k] =
                PETScSAMRAIVectorReal::createPETScVector(d_nullspace_basis_vecs[k], d_petsc_comm);
            nullspace_vecs.push_back(d_petsc_nullspace_basis_vecs[k]);
        }

        for (const auto& petsc_nvec : nullspace_vecs)
        {
            double norm;
            ierr = VecNorm(petsc_nvec, NORM_2, &norm);
            IBTK_CHKERRQ(ierr);
            ierr = VecScale(petsc_nvec, 1.0 / norm);
            IBTK_CHKERRQ(ierr);
        }

        static const PetscBool has_cnst = PETSC_FALSE;
        ierr = MatNullSpaceCreate(
            d_petsc_comm, has_cnst, static_cast<int>(nullspace_vecs.size()), &nullspace_vecs[0], &d_petsc_nullsp);
        IBTK_CHKERRQ(ierr);
        ierr = MatSetNullSpace(d_petsc_mat, d_petsc_nullsp);
        IBTK_CHKERRQ(ierr);
        d_solver_has_attached_nullspace = true;
    }
    else if (d_solver_has_attached_nullspace)
    {
        static const PetscBool has_cnst = PETSC_FALSE;
        ierr = MatNullSpaceCreate(d_petsc_comm, has_cnst, 0, nullptr, &d_petsc_nullsp);
        IBTK_CHKERRQ(ierr);
        ierr = MatSetNullSpace(d_petsc_mat, d_petsc_nullsp);
        IBTK_CHKERRQ(ierr);
    }
    return;
} // resetMatNullspace

void
PETScKrylovLinearSolver::deallocateNullspaceData()
{
    // Deallocate PETSc data structures.
    if (d_petsc_nullsp)
    {
        int ierr = MatNullSpaceDestroy(&d_petsc_nullsp);
        IBTK_CHKERRQ(ierr);
        d_petsc_nullsp = nullptr;
        if (d_petsc_mat)
        {
            ierr = MatSetNullSpace(d_petsc_mat, nullptr);
            IBTK_CHKERRQ(ierr);
        }
    }
    if (d_petsc_nullspace_constant_vec)
    {
        PETScSAMRAIVectorReal::destroyPETScVector(d_petsc_nullspace_constant_vec);
    }
    for (const auto& petsc_nullspace_basis_vec : d_petsc_nullspace_basis_vecs)
    {
        PETScSAMRAIVectorReal::destroyPETScVector(petsc_nullspace_basis_vec);
    }
    d_petsc_nullspace_constant_vec = nullptr;
    d_petsc_nullspace_basis_vecs.clear();

    // Clear SAMRAI data structures.
    if (d_nullspace_constant_vec)
    {
        d_nullspace_constant_vec->resetLevels(
            0,
            std::min(d_nullspace_constant_vec->getFinestLevelNumber(),
                     d_nullspace_constant_vec->getPatchHierarchy()->getFinestLevelNumber()));
        d_nullspace_constant_vec->deallocateVectorData();
        d_nullspace_constant_vec->freeVectorComponents();
    }
    d_nullspace_constant_vec.setNull();
    return;
} // deallocateNullspaceData

PetscErrorCode
PETScKrylovLinearSolver::MatVecMult_SAMRAI(Mat A, Vec x, Vec y)
{
    void* p_ctx;
    int ierr = MatShellGetContext(A, &p_ctx);
    CHKERRQ(ierr);
    auto krylov_solver = static_cast<PETScKrylovLinearSolver*>(p_ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(krylov_solver);
    TBOX_ASSERT(krylov_solver->d_A);
#endif
    Pointer<SAMRAIVectorReal<NDIM, double> > samrai_x, samrai_y;
    PETScSAMRAIVectorReal::getSAMRAIVectorRead(x, &samrai_x);
    PETScSAMRAIVectorReal::getSAMRAIVector(y, &samrai_y);
    krylov_solver->d_A->apply(*samrai_x, *samrai_y);
    PETScSAMRAIVectorReal::restoreSAMRAIVectorRead(x, &samrai_x);
    PETScSAMRAIVectorReal::restoreSAMRAIVector(y, &samrai_y);
    PetscFunctionReturn(0);
} // MatVecMult_SAMRAI

PetscErrorCode
PETScKrylovLinearSolver::PCApply_SAMRAI(PC pc, Vec x, Vec y)
{
    int ierr;
    void* ctx;
    ierr = PCShellGetContext(pc, &ctx);
    CHKERRQ(ierr);
    auto krylov_solver = static_cast<PETScKrylovLinearSolver*>(ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(krylov_solver);
    TBOX_ASSERT(krylov_solver->d_pc_solver);
#endif

    // Indicate that the initial guess should be zero.
    const bool pc_initial_guess_nonzero = krylov_solver->d_pc_solver->getInitialGuessNonzero();
    krylov_solver->d_pc_solver->setInitialGuessNonzero(false);

    // Apply the preconditioner.
    Pointer<SAMRAIVectorReal<NDIM, double> > samrai_x, samrai_y;
    PETScSAMRAIVectorReal::getSAMRAIVectorRead(x, &samrai_x);
    PETScSAMRAIVectorReal::getSAMRAIVector(y, &samrai_y);
    krylov_solver->d_pc_solver->solveSystem(*samrai_y, *samrai_x);
    PETScSAMRAIVectorReal::restoreSAMRAIVectorRead(x, &samrai_x);
    PETScSAMRAIVectorReal::restoreSAMRAIVector(y, &samrai_y);

    // Reset the configuration of the preconditioner object.
    krylov_solver->d_pc_solver->setInitialGuessNonzero(pc_initial_guess_nonzero);
    PetscFunctionReturn(0);
} // PCApply_SAMRAI

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
