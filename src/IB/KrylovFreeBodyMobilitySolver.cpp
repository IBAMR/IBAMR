// Filename: KrylovFreeBodyMobilitySolver.cpp
// Created on 16 Aug 2015 by Amneet Bhalla and Bakytzhan Kallemov.
//
// Copyright (c) 2002-2017, Amneet Bhalla, Bakytzhan Kallemov, and
// Boyce Griffith.
//
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
//    * Neither the name of The University of North Carolina nor the names of its
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <limits>

#include "ibamr/CIBMobilitySolver.h"
#include "ibamr/CIBStrategy.h"
#include "ibamr/DirectMobilitySolver.h"
#include "ibamr/KrylovFreeBodyMobilitySolver.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h"
#include "ibtk/ibtk_utilities.h"
#include "petsc/private/petscimpl.h"
#include "tbox/TimerManager.h"

namespace IBAMR
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

KrylovFreeBodyMobilitySolver::KrylovFreeBodyMobilitySolver(std::string object_name,
                                                           Pointer<Database> input_db,
                                                           std::string default_options_prefix,
                                                           Pointer<CIBStrategy> cib_strategy,
                                                           MPI_Comm petsc_comm)
    : d_object_name(std::move(object_name)),
      d_options_prefix(std::move(default_options_prefix)),
      d_petsc_comm(petsc_comm),
      d_cib_strategy(cib_strategy)
{
    // Get values from the input database.
    if (input_db)
    {
        if (input_db->keyExists("options_prefix")) d_options_prefix = input_db->getString("options_prefix");
        if (input_db->keyExists("max_iterations")) d_max_iterations = input_db->getInteger("max_iterations");
        if (input_db->keyExists("abs_residual_tol")) d_abs_residual_tol = input_db->getDouble("abs_residual_tol");
        if (input_db->keyExists("rel_residual_tol")) d_rel_residual_tol = input_db->getDouble("rel_residual_tol");
        if (input_db->keyExists("ksp_type")) d_ksp_type = input_db->getString("ksp_type");

        if (input_db->keyExists("pc_type")) d_pc_type = input_db->getString("pc_type");
        if (input_db->keyExists("initial_guess_nonzero"))
            d_initial_guess_nonzero = input_db->getBool("initial_guess_nonzero");
        if (input_db->keyExists("enable_logging")) d_enable_logging = input_db->getBool("enable_logging");
    }

    IBAMR_DO_ONCE(
        t_solve_system = TimerManager::getManager()->getTimer("IBAMR::KrylovFreeBodyMobilitySolver::solveSystem()");
        t_initialize_solver_state =
            TimerManager::getManager()->getTimer("IBAMR::KrylovFreeBodyMobilitySolver::initializeSolverState()");
        t_deallocate_solver_state =
            TimerManager::getManager()->getTimer("IBAMR::KrylovFreeBodyMobilitySolver::deallocateSolverState()"););
} // KrylovFreeBodyMobilitySolver

KrylovFreeBodyMobilitySolver::~KrylovFreeBodyMobilitySolver()
{
    if (d_is_initialized) deallocateSolverState();

    // Delete allocated PETSc solver components.
    if (d_petsc_mat)
    {
        MatDestroy(&d_petsc_mat);
        d_petsc_mat = nullptr;
    }
    if (d_petsc_ksp)
    {
        KSPDestroy(&d_petsc_ksp);
        d_petsc_ksp = nullptr;
    }
} // ~KrylovFreeBodyMobilitySolver

void
KrylovFreeBodyMobilitySolver::setMobilitySolver(Pointer<CIBMobilitySolver> mobility_solver)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(mobility_solver);
#endif
    d_mobility_solver = mobility_solver;
} // setMobilitySolver

void
KrylovFreeBodyMobilitySolver::setInterpScale(const double interp_scale)
{
    d_interp_scale = interp_scale;
} // setInterpScale

void
KrylovFreeBodyMobilitySolver::setSpreadScale(const double spread_scale)
{
    d_spread_scale = spread_scale;
} // setSpreadScale

void
KrylovFreeBodyMobilitySolver::setStokesSpecifications(const StokesSpecifications& stokes_spec)
{
    d_rho = stokes_spec.getRho();
    d_mu = stokes_spec.getMu();
} // setStokesSpecifications

const KSP&
KrylovFreeBodyMobilitySolver::getPETScKSP() const
{
    return d_petsc_ksp;
} // getPETScKSP

void
KrylovFreeBodyMobilitySolver::setKSPType(const std::string& ksp_type)
{
    d_ksp_type = ksp_type;
} // setKSPType

void
KrylovFreeBodyMobilitySolver::setOptionsPrefix(const std::string& options_prefix)
{
    d_options_prefix = options_prefix;
} // setOptionsPrefix

void
KrylovFreeBodyMobilitySolver::setSolutionTime(double solution_time)
{
    d_solution_time = solution_time;
} // setSolutionTime

void
KrylovFreeBodyMobilitySolver::setTimeInterval(double current_time, double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    d_dt = new_time - current_time;
} // setTimeInterval

bool
KrylovFreeBodyMobilitySolver::solveSystem(Vec x, Vec b)
{
    IBAMR_TIMER_START(t_solve_system);

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x, b);

    VecCopy(b, d_petsc_b);
    KSPSolve(d_petsc_ksp, d_petsc_b, x);
    KSPGetIterationNumber(d_petsc_ksp, &d_current_iterations);
    KSPGetResidualNorm(d_petsc_ksp, &d_current_residual_norm);

    // Determine the convergence reason.
    KSPConvergedReason reason;
    KSPGetConvergedReason(d_petsc_ksp, &reason);
    const bool converged = (static_cast<int>(reason) > 0);
    if (d_enable_logging) reportKSPConvergedReason(reason, plog);

    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();

    IBAMR_TIMER_STOP(t_solve_system);
    return converged;
} // solveSystem

void
KrylovFreeBodyMobilitySolver::initializeSolverState(Vec /*x*/, Vec b)
{
    IBAMR_TIMER_START(t_initialize_solver_state);

    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized)
    {
        d_reinitializing_solver = true;
        deallocateSolverState();
    }

    // Generate RHS and temporary vectors for storing Lagrange multiplier
    // and rigid body velocity.
    Vec* vb;
    VecNestGetSubVecs(b, nullptr, &vb);
    VecDuplicate(vb[2], &d_petsc_b);
    VecDuplicate(vb[1], &d_petsc_temp_f);
    VecDuplicate(vb[1], &d_petsc_temp_v);

    // Initialize PETSc KSP
    initializeKSP();

    // Indicate that the solver is initialized.
    d_reinitializing_solver = false;
    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_solver_state);
} // initializeSolverState

void
KrylovFreeBodyMobilitySolver::deallocateSolverState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_solver_state);

    // Destroy RHS and temp vecs.
    VecDestroy(&d_petsc_b);
    VecDestroy(&d_petsc_temp_f);
    VecDestroy(&d_petsc_temp_v);
    d_petsc_temp_v = nullptr;
    d_petsc_temp_f = nullptr;
    d_petsc_b = nullptr;

    // Destroy the KSP solver.
    destroyKSP();

    // Indicate that the solver is NOT initialized
    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_solver_state);
} // deallocateSolverState

/////////////////////////////// PRIVATE //////////////////////////////////////

void
KrylovFreeBodyMobilitySolver::reportKSPConvergedReason(const KSPConvergedReason& reason, std::ostream& os) const
{
    switch (static_cast<int>(reason))
    {
    case KSP_CONVERGED_RTOL:
        os << d_object_name
           << ": converged: |Ax-b| <= rtol*|b| --- residual norm is less than specified relative tolerance.\n";
        break;
    case KSP_CONVERGED_ATOL:
        os << d_object_name
           << ": converged: |Ax-b| <= atol --- residual norm is less than specified absolute tolerance.\n";
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
        os << d_object_name
           << ": diverged: reached maximum number of iterations before any convergence criteria were satisfied.\n";
        break;
    case KSP_DIVERGED_DTOL:
        os << d_object_name
           << ": diverged: |Ax-b| >= dtol*|b| --- residual is greater than specified divergence tolerance.\n";
        break;
    case KSP_DIVERGED_BREAKDOWN:
        os << d_object_name << ": diverged: breakdown in the Krylov method.\n";
        break;
    case KSP_DIVERGED_BREAKDOWN_BICG:
        os << d_object_name << ": diverged: breakdown in the bi-congugate gradient method.\n";
        break;
    case KSP_DIVERGED_NONSYMMETRIC:
        os << d_object_name
           << ": diverged: it appears the operator or preconditioner is not symmetric, but this "
              "Krylov method (KSPCG, KSPMINRES, KSPCR) requires symmetry\n";
        break;
    case KSP_DIVERGED_INDEFINITE_PC:
        os << d_object_name
           << ": diverged: it appears the preconditioner is indefinite (has both positive and "
              "negative eigenvalues), but this Krylov method (KSPCG) requires it to be positive "
              "definite.\n";
        break;
    case KSP_CONVERGED_ITERATING:
        os << d_object_name << ": iterating: KSPSolve() is still running.\n";
        break;
    default:
        os << d_object_name << ": unknown completion code " << static_cast<int>(reason) << " reported.\n";
        break;
    }
} // reportKSPConvergedReason

void
KrylovFreeBodyMobilitySolver::initializeKSP()
{
    // Create the KSP solver.
    KSPCreate(d_petsc_comm, &d_petsc_ksp);
    resetKSPOptions();
    resetKSPOperators();
    resetKSPPC();

    // Set the KSP options from the PETSc options database.
    if (d_options_prefix != "")
    {
        KSPSetOptionsPrefix(d_petsc_ksp, d_options_prefix.c_str());
    }
    KSPSetFromOptions(d_petsc_ksp);

    // Reset the member state variables to correspond to the values used by the
    // KSP object.  (Command-line options always take precedence.)
    KSPType ksp_type;
    KSPGetType(d_petsc_ksp, (const char**)&ksp_type);
    d_ksp_type = ksp_type;
    PetscBool initial_guess_nonzero;
    KSPGetInitialGuessNonzero(d_petsc_ksp, &initial_guess_nonzero);
    d_initial_guess_nonzero = (initial_guess_nonzero == PETSC_TRUE);
    KSPGetTolerances(d_petsc_ksp, &d_rel_residual_tol, &d_abs_residual_tol, nullptr, &d_max_iterations);
} // initializeKSP

void
KrylovFreeBodyMobilitySolver::destroyKSP()
{
    if (!d_petsc_ksp) return;
    KSPDestroy(&d_petsc_ksp);
    d_petsc_ksp = nullptr;
} // destroyKSP

void
KrylovFreeBodyMobilitySolver::resetKSPOptions()
{
    if (!d_petsc_ksp) return;
    const KSPType ksp_type = d_ksp_type.c_str();
    KSPSetType(d_petsc_ksp, ksp_type);
    std::string ksp_type_name(ksp_type);
    if (ksp_type_name.find("gmres") != std::string::npos)
    {
        KSPGMRESSetCGSRefinementType(d_petsc_ksp, KSP_GMRES_CGS_REFINE_IFNEEDED);
    }
    PetscBool initial_guess_nonzero = (d_initial_guess_nonzero ? PETSC_TRUE : PETSC_FALSE);
    KSPSetInitialGuessNonzero(d_petsc_ksp, initial_guess_nonzero);
    KSPSetTolerances(d_petsc_ksp, d_rel_residual_tol, d_abs_residual_tol, PETSC_DEFAULT, d_max_iterations);

    // Set KSP monitor routine.
    if (d_enable_logging)
    {
        KSPMonitorCancel(d_petsc_ksp);
        KSPMonitorSet(d_petsc_ksp,
                      reinterpret_cast<PetscErrorCode (*)(KSP, PetscInt, PetscReal, void*)>(
                          KrylovFreeBodyMobilitySolver::monitorKSP),
                      nullptr,
                      nullptr);
    }
} // resetKSPOptions

void
KrylovFreeBodyMobilitySolver::resetKSPOperators()
{
    // Create and configure the MatShell object.
    if (d_petsc_mat)
    {
        MatDestroy(&d_petsc_mat);
        d_petsc_mat = nullptr;
    }
    if (!d_petsc_mat)
    {
        int n;
        VecGetLocalSize(d_petsc_b, &n);
        MatCreateShell(d_petsc_comm, n, n, PETSC_DETERMINE, PETSC_DETERMINE, static_cast<void*>(this), &d_petsc_mat);
    }
    MatShellSetOperation(
        d_petsc_mat, MATOP_MULT, reinterpret_cast<void (*)(void)>(KrylovFreeBodyMobilitySolver::MatVecMult_KFBMSolver));

    // Reset the configuration of the PETSc KSP object.
    if (d_petsc_ksp)
    {
        KSPSetOperators(d_petsc_ksp, d_petsc_mat, d_petsc_mat);
        KSPSetReusePreconditioner(d_petsc_ksp, PETSC_TRUE);
    }
} // resetKSPOperators

void
KrylovFreeBodyMobilitySolver::resetKSPPC()
{
    if (!d_petsc_ksp) return;

    // Determine the preconditioner type to use.
    static const size_t len = 255;
    char pc_type_str[len];
    PetscBool flg;
    PetscOptionsGetString(nullptr, d_options_prefix.c_str(), "-pc_type", pc_type_str, len, &flg);
    std::string pc_type = d_pc_type;
    if (flg)
    {
        pc_type = std::string(pc_type_str);
    }

    if (!(pc_type == "none" || pc_type == "shell"))
    {
        TBOX_ERROR(d_object_name << "::resetKSPPC()\n"
                                 << "  valid values for -" << d_options_prefix << "pc_type are: none, shell"
                                 << std::endl);
    }

    PC petsc_pc;
    KSPGetPC(d_petsc_ksp, &petsc_pc);
    if (pc_type == "none")
    {
        PCSetType(petsc_pc, PCNONE);
    }
    else if (pc_type == "shell")
    {
        PCSetType(petsc_pc, PCSHELL);
        PCShellSetContext(petsc_pc, static_cast<void*>(this));
        PCShellSetApply(petsc_pc, KrylovFreeBodyMobilitySolver::PCApply_KFBMSolver);
    }
    else
    {
        TBOX_ERROR("This statement should not be reached!\n");
    }
} // resetKSPPC

PetscErrorCode
KrylovFreeBodyMobilitySolver::MatVecMult_KFBMSolver(Mat A, Vec x, Vec y)
{
    void* p_ctx;
    MatShellGetContext(A, &p_ctx);
    auto solver = static_cast<KrylovFreeBodyMobilitySolver*>(p_ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(solver);
#endif
    // a) Set rigid body velocity
    VecSet(solver->d_petsc_temp_v, 0.0); // so that lambda for specified-kinematics part is zero.
    solver->d_cib_strategy->setRigidBodyVelocity(x,
                                                 solver->d_petsc_temp_v,
                                                 /*only_free_dofs*/ true,
                                                 /*only_imposed_dofs*/ false);

    // b) Solve the mobility problem.
    solver->d_mobility_solver->solveMobilitySystem(solver->d_petsc_temp_f, solver->d_petsc_temp_v);

    // c) Compute the generalized force and torque.
    solver->d_cib_strategy->computeNetRigidGeneralizedForce(solver->d_petsc_temp_f,
                                                            y,
                                                            /*only_free_dofs*/ true,
                                                            /*only_imposed_dofs*/ false);
    PetscFunctionReturn(0);
} // MatVecMult_KFBMSolver

// Routine to apply FreeBodyMobility preconditioner
PetscErrorCode
KrylovFreeBodyMobilitySolver::PCApply_KFBMSolver(PC pc, Vec x, Vec y)
{
    // Here we are trying to the solve the problem of the type: Py = x for y.
    void* ctx;
    PCShellGetContext(pc, &ctx);
    auto solver = static_cast<KrylovFreeBodyMobilitySolver*>(ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(solver);
#endif
    DirectMobilitySolver* direct_solver;
    solver->d_mobility_solver->getMobilitySolvers(nullptr, &direct_solver, nullptr);
    direct_solver->solveBodySystem(y, x);
    VecScale(y, 1.0 / (solver->d_interp_scale * solver->d_spread_scale));
    PetscFunctionReturn(0);
} // PCApply_KFBMSolver

// Routine to log output of KrylovFreeBodyMobilitySolver
PetscErrorCode
KrylovFreeBodyMobilitySolver::monitorKSP(KSP ksp, int it, PetscReal rnorm, void* /*mctx*/)
{
    Vec resid, rhs;
    PetscReal truenorm, bnorm;
    char print_normtype[256];
    KSPNormType ksp_normtype;

    KSPBuildResidual(ksp, nullptr, nullptr, &resid);
    VecNorm(resid, NORM_2, &truenorm);
    VecDestroy(&resid);
    KSPGetRhs(ksp, &rhs);
    KSPGetNormType(ksp, &ksp_normtype);
    VecNorm(rhs, NORM_2, &bnorm);
    PetscStrncpy(print_normtype, KSPNormTypes[ksp_normtype], sizeof(print_normtype));
    PetscStrtolower(print_normtype);

    if (it == 0)
    {
        tbox::plog << "\n\n         Residual norms for -KFBMInv_ksp" << std::endl;
    }

    std::streamsize old_precision = tbox::plog.precision(16);
    tbox::plog << std::scientific << it << " KFBMInv_KSP " << print_normtype << " resid norm " << rnorm
               << " true resid norm " << truenorm << " ||r(i)||/||b|| " << truenorm / bnorm << std::endl;
    tbox::plog.precision(old_precision);
    PetscFunctionReturn(0);
} // monitorKSP

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
