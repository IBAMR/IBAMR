// Filename: FreeBodyMobilitySolver.C


/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h"
#include "ibtk/PETScMultiVec.h"
#include "tbox/TimerManager.h"
#include "ibamr/CIBStrategy.h"
#include "ibamr/CIBMobilitySolver.h"
#include "ibamr/FreeBodyMobilitySolver.h"
#include "petsc-private/petscimpl.h"
#include <limits>

namespace IBAMR
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

FreeBodyMobilitySolver::FreeBodyMobilitySolver(
    const std::string& object_name,
    Pointer<Database> input_db,
    const std::string& default_options_prefix,
    Pointer<CIBStrategy> cib_strategy,
    MPI_Comm petsc_comm)
    : d_object_name(object_name),
      d_ksp_type(KSPGMRES),
      d_pc_type("none"),
      d_is_initialized(false),
      d_reinitializing_solver(false),
      d_petsc_b(PETSC_NULL),
      d_petsc_temp_v(PETSC_NULL),
      d_petsc_temp_f(PETSC_NULL),
      d_options_prefix(default_options_prefix),
      d_petsc_comm  (petsc_comm),
      d_petsc_ksp   (PETSC_NULL),
      d_petsc_mat   (PETSC_NULL),
      d_cib_strategy(cib_strategy),
      d_MInv(NULL),
      d_num_rigid_parts(cib_strategy->getNumberOfRigidStructures()),
      d_current_time(std::numeric_limits<double>::signaling_NaN()),
      d_new_time(std::numeric_limits<double>::signaling_NaN()),
      d_dt(std::numeric_limits<double>::signaling_NaN())
{
    // Some default values of the Krylov solver.
    d_max_iterations = 10000;
    d_abs_residual_tol = 1.0e-50;
    d_rel_residual_tol = 1.0e-5;
    d_initial_guess_nonzero = false;
    d_enable_logging = false;

    // Get values from the input database.
    if (input_db)
    {
        if (input_db->keyExists("options_prefix"))        d_options_prefix        = input_db->getString("options_prefix");
        if (input_db->keyExists("max_iterations"))        d_max_iterations        = input_db->getInteger("max_iterations");
        if (input_db->keyExists("abs_residual_tol"))      d_abs_residual_tol      = input_db->getDouble("abs_residual_tol");
        if (input_db->keyExists("rel_residual_tol"))      d_rel_residual_tol      = input_db->getDouble("rel_residual_tol");
        if (input_db->keyExists("ksp_type"))              d_ksp_type              = input_db->getString("ksp_type");

	if (input_db->keyExists("pc_type"))               d_pc_type               = input_db->getString("pc_type");
	if (d_pc_type=="shell") 
	{
	    d_body_effect_hrad = input_db->getDouble("effective_body_hydrad");
	    d_mu=1.0;
	}
        if (input_db->keyExists("initial_guess_nonzero")) d_initial_guess_nonzero = input_db->getBool("initial_guess_nonzero");
        if (input_db->keyExists("enable_logging"))        d_enable_logging        = input_db->getBool("enable_logging");
    }
    
    IBTK_DO_ONCE(
        t_solve_system            = TimerManager::getManager()->getTimer("IBTK::FreeBodyMobilitySolver::solveSystem()");
        t_initialize_solver_state = TimerManager::getManager()->getTimer
                                                              ("IBTK::FreeBodyMobilitySolver::initializeSolverState()");
        t_deallocate_solver_state = TimerManager::getManager()->getTimer
                                                              ("IBTK::FreeBodyMobilitySolver::deallocateSolverState()");
                 );    
    return;
}// FreeBodyMobilitySolver


FreeBodyMobilitySolver::~FreeBodyMobilitySolver()
{
    if (d_is_initialized) deallocateSolverState();

    // Delete allocated PETSc solver components.
    int ierr;
    if (d_petsc_mat)
    {
        ierr = MatDestroy(&d_petsc_mat); IBTK_CHKERRQ(ierr);
        d_petsc_mat = PETSC_NULL;
    }
    if (d_petsc_ksp)
    {
        ierr = KSPDestroy(&d_petsc_ksp); IBTK_CHKERRQ(ierr);
        d_petsc_ksp = PETSC_NULL;
    }
    return;
}// ~FreeBodyMobilitySolver

void
FreeBodyMobilitySolver::setMobilitySolver(
    Pointer<CIBMobilitySolver> minv)
{

#if !defined(NDEBUG) 
    TBOX_ASSERT(dminv);
#endif

    d_MInv = minv;

    return;
}// setDirectMobilitySolver

void
FreeBodyMobilitySolver::setSolutionTime(
    double /*solution_time*/)
{
    // intentionally left blank.
    return;
}// setSolutionTime

void
FreeBodyMobilitySolver::setTimeInterval(
    double current_time,
    double new_time)
{
    d_current_time          = current_time;
    d_new_time              = new_time;
    d_dt                    = new_time-current_time;
     
    return;
}// setTimeInterval

bool
FreeBodyMobilitySolver::solveSystem(
    Vec x,
    Vec b)
{
    IBTK_TIMER_START(t_solve_system);
    int ierr;
    
    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;

    if (deallocate_after_solve) initializeSolverState(x,b);
    Vec *vb;
    IBTK::VecMultiVecGetSubVecs(b,&vb); 
    d_petsc_temp_v=vb[0];
    d_petsc_temp_f=vb[1];
    VecCopy (vb[2],d_petsc_b);

    ierr = KSPSolve(d_petsc_ksp, d_petsc_b, x);                          IBTK_CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(d_petsc_ksp, &d_current_iterations);    IBTK_CHKERRQ(ierr);
    ierr = KSPGetResidualNorm   (d_petsc_ksp, &d_current_residual_norm); IBTK_CHKERRQ(ierr);

    // Determine the convergence reason.
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(d_petsc_ksp, &reason); IBTK_CHKERRQ(ierr);
    const bool converged = (static_cast<int>(reason) > 0);
    if (d_enable_logging) reportKSPConvergedReason(reason, plog);
    
    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();

    IBTK_TIMER_STOP(t_solve_system);
    return converged;
}// solveSystem

void
FreeBodyMobilitySolver::initializeSolverState(
    Vec /*x*/,
    Vec b)
{
    IBTK_TIMER_START(t_initialize_solver_state);

    // Rudimentary error checking.
    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized)
    {
        d_reinitializing_solver = true;
        deallocateSolverState();
    }
    
    // Generate RHS and temporary vector for storing Lagrange multiplier.
    Vec *vb;
    IBTK::VecMultiVecGetSubVecs(b,&vb); 
    VecDuplicate(vb[2],&d_petsc_b); 
    
    // Initialize PETSc KSP
    initializeKSP();
  
    // Indicate that the solver is initialized.
    d_reinitializing_solver = false;
    d_is_initialized = true;

    IBTK_TIMER_STOP(t_initialize_solver_state);
    return;
}// initializeSolverState

void
FreeBodyMobilitySolver::deallocateSolverState()
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_deallocate_solver_state);
      
    // Destroy RHS and temp vecs.
    VecDestroy(&d_petsc_b);
    d_petsc_temp_v    = PETSC_NULL;
    d_petsc_temp_f    = PETSC_NULL;
    d_petsc_b    = PETSC_NULL;

    // Destroy the KSP solver.
    destroyKSP();

    // Indicate that the solver is NOT initialized
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_solver_state);

    return;
}// deallocateSolverState

/////////////////////////////// PRIVATE //////////////////////////////////////

void
FreeBodyMobilitySolver::setKSPType(
    const std::string& ksp_type)
{
    d_ksp_type = ksp_type;
    return;
}// setKSPType

void
FreeBodyMobilitySolver::setOptionsPrefix(
    const std::string& options_prefix)
{
    d_options_prefix = options_prefix;
    return;
}// setOptionsPrefix

const KSP&
FreeBodyMobilitySolver::getPETScKSP() const
{
    return d_petsc_ksp;
}// getPETScKSP

void
FreeBodyMobilitySolver::reportKSPConvergedReason(
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
            os << d_object_name << ": converged: maximum number of iterations reached.\n";
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
FreeBodyMobilitySolver::initializeKSP()
{
    // Create the KSP solver.
    int ierr;
    ierr = KSPCreate(d_petsc_comm, &d_petsc_ksp); IBTK_CHKERRQ(ierr);
    resetKSPOptions();
    resetKSPOperators();
    resetKSPPC();
    
    // Set the KSP options from the PETSc options database.
    if (d_options_prefix != "")
    {
        ierr = KSPSetOptionsPrefix(d_petsc_ksp, d_options_prefix.c_str()); IBTK_CHKERRQ(ierr);
    }
    ierr = KSPSetFromOptions(d_petsc_ksp); IBTK_CHKERRQ(ierr);    
    
    // Reset the member state variables to correspond to the values used by the
    // KSP object.  (Command-line options always take precedence.)
    KSPType ksp_type;
    ierr = KSPGetType(d_petsc_ksp, (const char**) &ksp_type); IBTK_CHKERRQ(ierr);
    d_ksp_type = ksp_type;
    PetscBool initial_guess_nonzero;
    ierr = KSPGetInitialGuessNonzero(d_petsc_ksp, &initial_guess_nonzero); IBTK_CHKERRQ(ierr);
    d_initial_guess_nonzero = (initial_guess_nonzero == PETSC_TRUE);
    ierr = KSPGetTolerances(d_petsc_ksp, &d_rel_residual_tol, &d_abs_residual_tol, PETSC_NULL, &d_max_iterations); IBTK_CHKERRQ(ierr);
   
    return;
}// initializeKSP

void
FreeBodyMobilitySolver::destroyKSP()
{
    int ierr = KSPDestroy(&d_petsc_ksp);  IBTK_CHKERRQ(ierr);
    d_petsc_ksp = PETSC_NULL;

    return;   
}// destroyKSP

void
FreeBodyMobilitySolver::resetKSPOptions()
{
    if (!d_petsc_ksp) return;
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
    ierr = KSPSetTolerances(d_petsc_ksp, d_rel_residual_tol, d_abs_residual_tol, PETSC_DEFAULT, d_max_iterations); IBTK_CHKERRQ(ierr);

    //Set KSP monitor routine.
    if (d_enable_logging)
    {
        ierr = KSPMonitorCancel(d_petsc_ksp); IBTK_CHKERRQ(ierr);
        ierr = KSPMonitorSet(d_petsc_ksp,reinterpret_cast<PetscErrorCode(*)(KSP,PetscInt,PetscReal,void*)>
               (FreeBodyMobilitySolver::monitorKSP),PETSC_NULL,PETSC_NULL) ; IBTK_CHKERRQ(ierr);
    }
    return;
}// resetKSPOptions

void
FreeBodyMobilitySolver::resetKSPOperators()
{
    int ierr;

    // Create and configure the MatShell object.
    if (d_petsc_mat)
    {
        ierr = MatDestroy(&d_petsc_mat); IBTK_CHKERRQ(ierr);
        d_petsc_mat = PETSC_NULL;
    }
    if (!d_petsc_mat)
    {
        int n;
	ierr = VecGetLocalSize(d_petsc_b,&n); IBTK_CHKERRQ(ierr);
        ierr = MatCreateShell(d_petsc_comm, n, n, PETSC_DETERMINE, PETSC_DETERMINE, static_cast<void*>(this), &d_petsc_mat); IBTK_CHKERRQ(ierr);
    }
    ierr = MatShellSetOperation(d_petsc_mat, MATOP_MULT, reinterpret_cast<void(*)(void)>
            (FreeBodyMobilitySolver::MatVecMult_KBMInv)); IBTK_CHKERRQ(ierr);

    // Reset the configuration of the PETSc KSP object.
    if (d_petsc_ksp)
    {
      KSPSetOperators(d_petsc_ksp, d_petsc_mat, d_petsc_mat);	
      KSPSetReusePreconditioner(d_petsc_ksp, PETSC_TRUE);

    }
    return;
}// resetKSPOperators

void
FreeBodyMobilitySolver::resetKSPPC()
{
    if (!d_petsc_ksp) return;
    int ierr;

    // Determine the preconditioner type to use.
    static const size_t len = 255;
    char pc_type_str[len];
    PetscBool flg;
    ierr = PetscOptionsGetString(d_options_prefix.c_str(), "-pc_type", pc_type_str, len, &flg);  IBTK_CHKERRQ(ierr);
    std::string pc_type = d_pc_type;
    if (flg)
    {
        pc_type = std::string(pc_type_str);
    }

    if (!(pc_type == "none" || pc_type == "shell"))
    {
        TBOX_ERROR(d_object_name << "::resetKSPPC()\n"
                   << "  valid values for -" << d_options_prefix << "pc_type are: none, shell" << std::endl);
    }

    PC petsc_pc;
    ierr = KSPGetPC(d_petsc_ksp, &petsc_pc); IBTK_CHKERRQ(ierr);
    if (pc_type == "none")
    {
        ierr = PCSetType(petsc_pc, PCNONE); IBTK_CHKERRQ(ierr);
    }
    else if (pc_type == "shell")
    {
        const std::string pc_name = d_object_name + pc_type;
        ierr = PCSetType(petsc_pc, PCSHELL); IBTK_CHKERRQ(ierr);
        ierr = PCShellSetContext(petsc_pc, static_cast<void*>(this)); IBTK_CHKERRQ(ierr);
        ierr = PCShellSetApply(petsc_pc, FreeBodyMobilitySolver::PCApply_KBMInv); IBTK_CHKERRQ(ierr);
    }
    else
    {
        TBOX_ERROR("This statement should not be reached!\n");
    }
    return;
}// resetKSPPC

PetscErrorCode
FreeBodyMobilitySolver::MatVecMult_KBMInv(
    Mat A,
    Vec x,
    Vec y)
{
    int ierr;
    void* p_ctx;
    ierr = MatShellGetContext(A, &p_ctx); IBTK_CHKERRQ(ierr);
    FreeBodyMobilitySolver* solver = static_cast<FreeBodyMobilitySolver*>(p_ctx);

    Vec* mv_U;
    IBTK::VecMultiVecGetSubVecs(x, &mv_U);

    for (unsigned part = 0, k = 0; part < solver->d_num_rigid_parts; ++part)
    {
	if (solver->d_cib_strategy->getSolveRigidBodyVelocity(part)) 
	{
	    solver->d_cib_strategy->setRigidBodyVelocity(part, mv_U[k], solver->d_petsc_temp_v);
	    ++k;
	}
    }

    solver->d_MInv->solveMobilitySystem(solver->d_petsc_temp_f,solver->d_petsc_temp_v);


    solver->d_cib_strategy->computeNetRigidGeneralizedForce(solver->d_petsc_temp_f, y, /*only_free_parts*/ true,
						    /*only_imposed_parts*/ false);

    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y)); IBTK_CHKERRQ(ierr);   
    PetscFunctionReturn(0);

}// MatVecMult_KBMInv

// Routine to apply DirectMobility preconditioner 
PetscErrorCode
FreeBodyMobilitySolver::PCApply_KBMInv(
    PC pc,
    Vec x,
    Vec y)
{
    int ierr;
    void* ctx;
    ierr = PCShellGetContext(pc, &ctx); IBTK_CHKERRQ(ierr);
    FreeBodyMobilitySolver* solver = static_cast<FreeBodyMobilitySolver*>(ctx);

    int free_comps;
    Vec* vx, *vy;
    IBTK::VecMultiVecGetSubVecs(x, &vx);
    IBTK::VecMultiVecGetSubVecs(y, &vy);
    IBTK::VecMultiVecGetNumberOfSubVecs(x, &free_comps);
  

    const double Mob_tt=1.0/(6.0*M_PI*solver->d_body_effect_hrad*solver->d_mu); //translational self-mobility
    const double Mob_rr=1.0/(8.0*M_PI*solver->d_mu*pow(solver->d_body_effect_hrad,3));//rotational self-mobility


    for (int part = 0; part < free_comps; ++part)
    {
	RigidDOFVector Ur, Fr;

	solver->d_cib_strategy->vecToRDV(vx[part],Fr); 

	for (int d = 0; d<NDIM; ++d)
	{
	    Ur[d] = Mob_tt*Fr[d];
#if (NDIM==3)
	    Ur[d+3] = Mob_rr*Fr[d+3];
	}
#elif (NDIM==2)
        }    
            Ur[2] = Mob_rr*Fr[2];
#endif
       solver->d_cib_strategy->rdvToVec(Ur,vy[part]); 
    }

    PetscFunctionReturn(0);
}// PCApply_KBMInv

// Routine to log output of FreeBodyMobilitySolver
PetscErrorCode
FreeBodyMobilitySolver::monitorKSP(
    KSP ksp, 
    int it, 
    PetscReal rnorm, 
    void* /*mctx*/)
{
    Vec            resid, rhs;
    PetscReal      truenorm,bnorm;
    char           print_normtype[256];
    KSPNormType    ksp_normtype;
   
    KSPBuildResidual(ksp,NULL,NULL,&resid);
    VecNorm(resid,NORM_2,&truenorm);
    VecDestroy(&resid);
    KSPGetRhs(ksp,&rhs);
    KSPGetNormType(ksp,&ksp_normtype);
    VecNorm(rhs,NORM_2,&bnorm);
    PetscStrncpy(print_normtype,KSPNormTypes[ksp_normtype],sizeof(print_normtype));
    PetscStrtolower(print_normtype);

    if (it == 0)
    {
        tbox::plog << "\n\n         Residual norms for -FBMInv_ksp" << std::endl;
    } 
    
    std::streamsize old_precision = tbox::plog.precision(16);
    tbox::plog << std::scientific << it << " FBMInv_KSP " << print_normtype  << " resid norm " 
               << (double)rnorm << " true resid norm "   << (double)truenorm << " ||r(i)||/||b|| " 
               << (double)(truenorm/bnorm)<< std::endl;
    tbox::plog.precision(old_precision);
    
   return(0);
  
}//monitorKSP


//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
