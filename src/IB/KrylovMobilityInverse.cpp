// Filename: KrylovMobilityInverse.C
// Created on 28 Oct 2013 by Amneet Bhalla 

#include "KrylovMobilityInverse.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h"
#include "ibtk/PETScMultiVec.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "tbox/TimerManager.h"
#include "ibtk/PoissonSolver.h"
#include "ibamr/INSStaggeredHierarchyIntegrator.h"
#include "limits"
#include "InterpOperator.h"
#include "SpreadOperator.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibamr/StaggeredStokesSolverManager.h"
#include "ibtk/SCPoissonSolverManager.h"
#include "ibtk/CCPoissonSolverManager.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/NewtonKrylovSolver.h"
#include "ibamr/StaggeredStokesBlockPreconditioner.h"
#include "DirectMobilityInverse.h"
#include "cIBMethod.h"

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

KrylovMobilityInverse::KrylovMobilityInverse(
    const std::string& object_name,
    Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator,
    Pointer<cIBMethod> cib_method,
    Pointer<Database> input_db,
    const std::string& default_options_prefix,
    MPI_Comm petsc_comm)
    : d_object_name(object_name),
      d_ksp_type(KSPGMRES),
      d_pc_type("none"),
      d_is_initialized(false),
      d_reinitializing_solver(false),
      d_petsc_x(PETSC_NULL),
      d_petsc_b(PETSC_NULL),
      d_options_prefix(default_options_prefix),
      d_petsc_comm  (petsc_comm),
      d_petsc_ksp   (PETSC_NULL),
      d_petsc_mat   (PETSC_NULL),
      d_samrai_temp(2,Pointer<SAMRAIVectorReal<NDIM,PetscScalar> >(NULL)),
      d_ins_integrator(navier_stokes_integrator),
      d_cib_method(cib_method),
      d_S(NULL),
      d_J(NULL),
      d_LInv(NULL),
      d_nul_vecs(),
      d_U_nul_vecs(),
      d_normalize_pressure(false),
      d_normalize_velocity(false),
      d_current_time(std::numeric_limits<double>::signaling_NaN()),
      d_new_time(std::numeric_limits<double>::signaling_NaN()),
      d_dt(std::numeric_limits<double>::signaling_NaN()),
      d_scale_interp(1.0),
      d_scale_spread(1.0),
      d_reg_mob_factor(0.0)
{
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
        if (input_db->keyExists("initial_guess_nonzero")) d_initial_guess_nonzero = input_db->getBool("initial_guess_nonzero");
        if (input_db->keyExists("normalize_pressure"))    d_normalize_pressure    = input_db->getBool("normalize_pressure");
        if (input_db->keyExists("normalize_velocity"))    d_normalize_velocity    = input_db->getBool("normalize_velocity");
        if (input_db->keyExists("enable_logging"))        d_enable_logging        = input_db->getBool("enable_logging");
    }
    
    // Create the Stokes solver (LInv) for the linear operator.
    // Create databases for setting up LInv solver.
    std::string stokes_solver_type     = StaggeredStokesSolverManager::PETSC_KRYLOV_SOLVER;
    Pointer<Database> stokes_solver_db = NULL;
    if (input_db->keyExists("stokes_solver_type"))
    {
	stokes_solver_type = input_db->getString("stokes_solver_type");
	if (input_db->keyExists("stokes_solver_db")) stokes_solver_db = input_db->getDatabase("stokes_solver_db");
    }
    if (!stokes_solver_db) 
    {
	stokes_solver_db = new MemoryDatabase("stokes_solver_db");
	stokes_solver_db->putString("ksp_type", "fgmres");
    }

    std::string stokes_precond_type     = StaggeredStokesSolverManager::DEFAULT_BLOCK_PRECONDITIONER;
    Pointer<Database> stokes_precond_db = NULL;
    if (input_db->keyExists("stokes_precond_type"))
    {
	stokes_precond_type = input_db->getString("stokes_precond_type");
	if (input_db->keyExists("stokes_precond_db")) stokes_precond_db = input_db->getDatabase("stokes_precond_db");
    }
    if (!stokes_precond_db)
    {
	stokes_precond_db = new MemoryDatabase("stokes_precond_db");
	stokes_precond_db->putInteger("max_iterations", 1);
    }

    std::string velocity_solver_type     = IBTK::SCPoissonSolverManager::PETSC_KRYLOV_SOLVER;
    Pointer<Database> velocity_solver_db = NULL;
    if (input_db->keyExists("velocity_solver_type") ) 
    {
	velocity_solver_type  = input_db->getString("velocity_solver_type");
	if (input_db->keyExists("velocity_solver_db")) velocity_solver_db = input_db->getDatabase("velocity_solver_db");
    }
    if (!velocity_solver_db)
    {
	velocity_solver_db = new MemoryDatabase("velocity_solver_db");
	velocity_solver_db->putString("ksp_type", "richardson");
	velocity_solver_db->putInteger("max_iterations", 10);
	velocity_solver_db->putDouble("rel_residual_tol", 1.0e-1);
    }

    std::string velocity_precond_type     = IBTK::SCPoissonSolverManager::DEFAULT_FAC_PRECONDITIONER;
    Pointer<Database> velocity_precond_db = NULL;
    if (input_db->keyExists("velocity_precond_type")) 
    {
	velocity_precond_type = input_db->getString("velocity_precond_type");
	if (input_db->keyExists("velocity_precond_db")) velocity_precond_db = input_db->getDatabase("velocity_precond_db");
    }
    if (!velocity_precond_db)
    {
	velocity_precond_db = new MemoryDatabase("velocity_precond_db");
	velocity_precond_db->putInteger("max_iterations", 1);
    }

    std::string pressure_solver_type     = IBTK::CCPoissonSolverManager::PETSC_KRYLOV_SOLVER;;
    Pointer<Database> pressure_solver_db = NULL;
    if (input_db->keyExists("pressure_solver_type")) 
    {
	pressure_solver_type  = input_db->getString("pressure_solver_type");
	if (input_db->keyExists("pressure_solver_db")) pressure_solver_db = input_db->getDatabase("pressure_solver_db");
    }
    if (!pressure_solver_db)
    {
	pressure_solver_db = new MemoryDatabase("pressure_solver_db");
	pressure_solver_db->putString("ksp_type", "richardson");
	pressure_solver_db->putInteger("max_iterations", 10);
	pressure_solver_db->putDouble("rel_residual_tol", 1.0e-1);
    }

    std::string pressure_precond_type     = IBTK::CCPoissonSolverManager::DEFAULT_FAC_PRECONDITIONER;
    Pointer<Database> pressure_precond_db = NULL;
    if (input_db->keyExists("pressure_precond_type")) 
    {
	pressure_precond_type = input_db->getString("pressure_precond_type");
	if (input_db->keyExists("pressure_precond_db")) pressure_precond_db = input_db->getDatabase("pressure_precond_db");
    }
    if (!pressure_precond_db)
    {
	pressure_precond_db = new MemoryDatabase("pressure_precond_db");
	pressure_precond_db->putInteger("max_iterations", 1);
    }

    // Create LInv.
    d_LInv = StaggeredStokesSolverManager::getManager()->allocateSolver(
	stokes_solver_type , d_object_name+"::pc_stokes_solver" , stokes_solver_db , "KM_LInv_"   ,
	stokes_precond_type, d_object_name+"::pc_stokes_precond", stokes_precond_db, "KM_LInv_pc_");

    // Create velocity solver.
    d_velocity_solver = IBTK::SCPoissonSolverManager::getManager()->allocateSolver(
	velocity_solver_type , d_object_name+"::velocity_solver",  velocity_solver_db , "KM_LInv_velocity_"   ,
	velocity_precond_type, d_object_name+"::velocity_precond", velocity_precond_db, "KM_LInv_velocity_pc_");

    // Create pressure solver.
    d_pressure_solver = IBTK::CCPoissonSolverManager::getManager()->allocateSolver(
	pressure_solver_type , d_object_name+"::pressure_solver" , pressure_solver_db , "KM_LInv_pressure_"   ,
	pressure_precond_type, d_object_name+"::pressure_precond", pressure_precond_db, "KM_LInv_pressure_pc_");
    
    // Register Poisson specification 
    const StokesSpecifications& stokes_spec = *d_ins_integrator->getStokesSpecifications();
    const double rho                        = stokes_spec.getRho(); 
    PoissonSpecifications P_problem_coefs(d_object_name+"::P_problem_coefs");
    P_problem_coefs.setCZero();
    P_problem_coefs.setDConstant(rho == 0.0 ? -1.0 : -1.0/rho);
    d_pressure_solver->setPoissonSpecifications(P_problem_coefs);

    // Register velocity and pressure solvers with LInv.
    Pointer<IBTK::LinearSolver> p_stokes_linear_solver = d_LInv;
    if (!p_stokes_linear_solver)
    {
	Pointer<IBTK::NewtonKrylovSolver> p_stokes_newton_solver = d_LInv;
	if (p_stokes_newton_solver) p_stokes_linear_solver = p_stokes_newton_solver->getLinearSolver();
    }
    if (p_stokes_linear_solver)
    {
	Pointer<StaggeredStokesBlockPreconditioner> p_stokes_block_pc = p_stokes_linear_solver;
	if (!p_stokes_block_pc)
	{
	    Pointer<IBTK::KrylovLinearSolver> p_stokes_krylov_solver = p_stokes_linear_solver;
	    if (p_stokes_krylov_solver) p_stokes_block_pc      = p_stokes_krylov_solver->getPreconditioner();
	}
	if (p_stokes_block_pc)
	{
	    if (p_stokes_block_pc->needsVelocitySubdomainSolver())
	    {
		p_stokes_block_pc->setVelocitySubdomainSolver(d_velocity_solver);
	    }
	    if (p_stokes_block_pc->needsPressureSubdomainSolver())
	    {
		p_stokes_block_pc->setPressureSubdomainSolver(d_pressure_solver);
                p_stokes_block_pc->setPressurePoissonSpecifications(P_problem_coefs);
	    }
	}
    }
    
    IBTK_DO_ONCE(
        t_solve_system            = TimerManager::getManager()->getTimer("IBTK::KrylovMobilityInverse::solveSystem()");
        t_initialize_solver_state = TimerManager::getManager()->getTimer("IBTK::KrylovMobilityInverse::initializeSolverState()");
        t_deallocate_solver_state = TimerManager::getManager()->getTimer("IBTK::KrylovMobilityInverse::deallocateSolverState()");
                 );    
    return;
}// KrylovMobilityInverse


KrylovMobilityInverse::~KrylovMobilityInverse()
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
}// ~KrylovMobilityInverse

void
KrylovMobilityInverse::setInterpScaleFactor(
    const double beta)
{
    d_scale_interp = beta;
    return;
}// setInterpScaleFactor

void
KrylovMobilityInverse::setSpreadScaleFactor(
    const double gamma)
{
    d_scale_spread = gamma;
    return;
}// setSpreadScaleFactor  

void
KrylovMobilityInverse::setRegularizeMobilityFactor(
    const double delta)
{
    d_reg_mob_factor = delta;
    return;
}// setRegularizeMobilityFactor

void
KrylovMobilityInverse::setKSPType(
    const std::string& ksp_type)
{
    d_ksp_type = ksp_type;
    return;
}// setKSPType

void
KrylovMobilityInverse::setOptionsPrefix(
    const std::string& options_prefix)
{
    d_options_prefix = options_prefix;
    return;
}// setOptionsPrefix

const KSP&
KrylovMobilityInverse::getPETScKSP() const
{
    return d_petsc_ksp;
}// getPETScKSP

void
KrylovMobilityInverse::setLinearOperators(
    Pointer<InterpOperator> J,
    Pointer<SpreadOperator> S)
{
#if !defined(NDEBUG) 
    TBOX_ASSERT(J);
    TBOX_ASSERT(S);
#endif
    
    d_J = J;
    d_S = S;

    return;
}// setLinearOperators

void
KrylovMobilityInverse::setPreconditioner(
    Pointer<DirectMobilityInverse> pc)
{
#if !defined(NDEBUG) 
    TBOX_ASSERT(pc);
#endif

    d_DMInv = pc;
    return;
}// setPreconditioner

void
KrylovMobilityInverse::setVelocityPoissonSpecifications(
    const PoissonSpecifications& u_problem_coefs)
{
    d_LInv->setVelocityPoissonSpecifications(u_problem_coefs); 
    d_velocity_solver->setPoissonSpecifications(u_problem_coefs);

    return;
}// setVelocityPoissonSpecifications

void
KrylovMobilityInverse::setSolutionTime(
    double solution_time)
{
    d_LInv->setSolutionTime(solution_time); 

    return;
}// setSolutionTime

void
KrylovMobilityInverse::setTimeInterval(
    double current_time,
    double new_time)
{
    d_current_time          = current_time;
    d_new_time              = new_time;
    d_dt                    = new_time-current_time;
    const double half_time  = current_time + 0.5*d_dt;

    d_LInv->setTimeInterval(current_time, new_time); 
    d_velocity_solver->setSolutionTime(new_time);
    d_pressure_solver->setSolutionTime(half_time);  
    d_velocity_solver->setTimeInterval(current_time, new_time);
    d_pressure_solver->setTimeInterval(current_time, new_time);
     
    return;
}// setTimeInterval

void
KrylovMobilityInverse::setPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    RobinBcCoefStrategy<NDIM>* p_bc_coef)
{
    d_LInv->setPhysicalBcCoefs(u_bc_coefs, p_bc_coef); 
    d_velocity_solver->setPhysicalBcCoefs(d_ins_integrator->getIntermediateVelocityBoundaryConditions());
    d_pressure_solver->setPhysicalBcCoef(d_ins_integrator->getProjectionBoundaryConditions());

    return;
}// setPhysicalBcCoefs

void
KrylovMobilityInverse::setPhysicalBoundaryHelper(
    Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper)
{
    d_LInv->setPhysicalBoundaryHelper(bc_helper); 
    
    return;
}// setPhysicalBoundaryHelper

bool
KrylovMobilityInverse::solveSystem(
    Vec x,
    Vec b)
{
    IBTK_TIMER_START(t_solve_system);
    int ierr;
    
    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;

    if (deallocate_after_solve) initializeSolverState(x,b);
#if !defined(NDEBUG)
    TBOX_ASSERT(d_petsc_ksp);
#endif
    
    // Create d_petsc_x Vec with memory pointing to input x
    int comps;
    Vec *vx;
    ierr = IBTK::VecMultiVecGetSubVecs(x,&vx); IBTK_CHKERRQ(ierr);
    ierr = IBTK::VecMultiVecGetNumberOfSubVecs(x,&comps); IBTK_CHKERRQ(ierr);
    ierr = IBTK::VecCreateMultiVec(d_petsc_comm,comps,vx,&d_petsc_x); IBTK_CHKERRQ(ierr);       
    
    // Copy the RHS Vec into d_petsc_b   
    ierr = VecCopy(b,d_petsc_b); IBTK_CHKERRQ(ierr);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(d_petsc_b)); IBTK_CHKERRQ(ierr);

    // Solve the system using a PETSc KSP object.
    ierr = KSPSolve(d_petsc_ksp, d_petsc_b, d_petsc_x); IBTK_CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(d_petsc_ksp, &d_current_iterations); IBTK_CHKERRQ(ierr);
    ierr = KSPGetResidualNorm   (d_petsc_ksp, &d_current_residual_norm); IBTK_CHKERRQ(ierr);

    // Determine the convergence reason.
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(d_petsc_ksp, &reason); IBTK_CHKERRQ(ierr);
    const bool converged = (static_cast<int>(reason) > 0);
    if (d_enable_logging) reportKSPConvergedReason(reason, plog);

    // Destroy d_petsc_x MultiVec
    ierr = VecDestroy(&d_petsc_x); IBTK_CHKERRQ(ierr);
    d_petsc_x = PETSC_NULL;   
    
    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();

    IBTK_TIMER_STOP(t_solve_system);
    return converged;
}// solveSystem

void
KrylovMobilityInverse::initializeSolverState(
    Vec x,
    Vec b)
{
    IBTK_TIMER_START(t_initialize_solver_state);

    int ierr;

    // Rudimentary error checking.
    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized)
    {
        d_reinitializing_solver = true;
        deallocateSolverState();
    }

    int comps;
    Vec *vx, *vb;
    ierr = IBTK::VecMultiVecGetSubVecs(x,&vx);            IBTK_CHKERRQ(ierr);
    ierr = IBTK::VecMultiVecGetSubVecs(b,&vb);            IBTK_CHKERRQ(ierr);
    ierr = IBTK::VecMultiVecGetNumberOfSubVecs(x,&comps); IBTK_CHKERRQ(ierr);
    
    // Create the temporary storage for spreading and Stokes solve operation.
    for (int i = 0; i < 2; ++i)
    {
        d_samrai_temp[i]  = IBTK::PETScSAMRAIVectorReal::getSAMRAIVector(vx[0])->cloneVector("");         
        d_samrai_temp[i]->allocateVectorData();
    }
    
    // Create Vecs to be used in the KSP object.
    Vec* r;
    ierr = PetscMalloc(sizeof(Vec),&r);                          IBTK_CHKERRQ(ierr); 
    ierr = VecDuplicate(vb[comps-1],&r[0]);                      IBTK_CHKERRQ(ierr);
    ierr = IBTK::VecCreateMultiVec(d_petsc_comm,1,r,&d_petsc_b); IBTK_CHKERRQ(ierr);
     
    // Initialize PETSc KSP
    initializeKSP();
  
    // Initialize LInv (Stokes solver) required in the mobility matrix. 
    initializeStokesSolver(*IBTK::PETScSAMRAIVectorReal::getSAMRAIVector(vx[0]),*IBTK::PETScSAMRAIVectorReal::getSAMRAIVector(vb[0]));

    // Indicate that the solver is initialized.
    d_reinitializing_solver = false;
    d_is_initialized = true;

    IBTK_TIMER_STOP(t_initialize_solver_state);
    return;
}// initializeSolverState

void
KrylovMobilityInverse::deallocateSolverState()
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_deallocate_solver_state);
    int ierr;

    // Deallocate the operator and preconditioner states only if we are not
    // re-initializing the solver.
    if (!d_reinitializing_solver)
    {
        d_LInv->deallocateSolverState(); 
    }
    
    // Delete the temporary storage for spreading and Stokes solve operation.
    for (int i = 0; i < 2; ++i)
    {
        d_samrai_temp[i]->resetLevels(d_samrai_temp[i]->getCoarsestLevelNumber(), std::min(d_samrai_temp[i]->getFinestLevelNumber(),d_samrai_temp[i]->getPatchHierarchy()->getFinestLevelNumber()));
        d_samrai_temp[i]->freeVectorComponents();
        d_samrai_temp[i].setNull();
    }
    
    int comps;
    Vec* vb;
    ierr = IBTK::VecMultiVecGetSubVecs(d_petsc_b,&vb); IBTK_CHKERRQ(ierr);
    ierr = IBTK::VecMultiVecGetNumberOfSubVecs(d_petsc_b,&comps); IBTK_CHKERRQ(ierr);
    for (int i = 0; i < comps; ++i) 
    {  
        ierr = VecDestroy(&vb[i]); IBTK_CHKERRQ(ierr);
    }
    ierr = PetscFree(vb); IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&d_petsc_b); IBTK_CHKERRQ(ierr);
    d_petsc_x = PETSC_NULL;
    d_petsc_b = PETSC_NULL;

    // Destroy the KSP solver.
    destroyKSP();

    // Indicate that the solver is NOT initialized
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_solver_state);
    return;
}// deallocateSolverState

/////////////////////////////// PRIVATE //////////////////////////////////////

void
KrylovMobilityInverse::initializeStokesSolver(
    const SAMRAIVectorReal<NDIM,double>& sol_vec,
    const SAMRAIVectorReal<NDIM,double>& rhs_vec)
{
    Pointer<PatchHierarchy<NDIM> > patch_hier = sol_vec.getPatchHierarchy();
    const int coarsest_ln                     = sol_vec.getCoarsestLevelNumber();
    const int finest_ln                       = sol_vec.getFinestLevelNumber();

    //Set the nullspace of the LInv and subdomain solvers
    const double rho = d_ins_integrator->getStokesSpecifications()->getRho();
    const bool has_velocity_nullspace = d_normalize_velocity && MathUtilities<double>::equalEps(rho, 0.0);
    const bool has_pressure_nullspace = d_normalize_pressure;

    for (unsigned int k = 0; k < d_nul_vecs.size(); ++k)
    {
        if (d_nul_vecs[k]) d_nul_vecs[k]->freeVectorComponents();
    }
    const int n_nul_vecs = (has_pressure_nullspace ? 1 : 0) + (has_velocity_nullspace ? NDIM : 0);
    d_nul_vecs.resize(n_nul_vecs);

    for (unsigned int k = 0; k < d_U_nul_vecs.size(); ++k)
    {
        if (d_U_nul_vecs[k]) d_U_nul_vecs[k]->freeVectorComponents();
    }
    const int n_U_nul_vecs = (has_velocity_nullspace ? NDIM : 0);
    d_U_nul_vecs.resize(n_U_nul_vecs);

    if (has_velocity_nullspace)
    {
        for (unsigned int k = 0; k < NDIM; ++k)
        {
            std::ostringstream stream;
            stream << k;
            d_nul_vecs[k] = sol_vec.cloneVector(d_object_name+"::nul_vec_U_"+stream.str());
            d_nul_vecs[k]->allocateVectorData(d_current_time);
            d_nul_vecs[k]->setToScalar(0.0);

            SAMRAIVectorReal<NDIM,double> svr_u(d_object_name+"::U_nul_vec_U_"+stream.str(), patch_hier, coarsest_ln, finest_ln); 
            svr_u.addComponent(sol_vec.getComponentVariable(0),sol_vec.getComponentDescriptorIndex(0),sol_vec.getControlVolumeIndex(0));
            
            d_U_nul_vecs[k] = svr_u.cloneVector(svr_u.getName());
            d_U_nul_vecs[k]->allocateVectorData(d_current_time);
            d_U_nul_vecs[k]->setToScalar(0.0);
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = patch_hier->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    Pointer<SideData<NDIM,double> > nul_data = patch->getPatchData(d_nul_vecs[k]->getComponentDescriptorIndex(0));
                    nul_data->getArrayData(k).fillAll(1.0);
                    Pointer<SideData<NDIM,double> > U_nul_data = patch->getPatchData(d_U_nul_vecs[k]->getComponentDescriptorIndex(0));
                    U_nul_data->getArrayData(k).fillAll(1.0);
                }
            }
        }
    }

    if (has_pressure_nullspace)
    {
        d_nul_vecs.back() = sol_vec.cloneVector(d_object_name+"::nul_vec_p");
        d_nul_vecs.back()->allocateVectorData(d_current_time);
        
        HierarchySideDataOpsReal<NDIM,double> side_ops(patch_hier);
        HierarchyCellDataOpsReal<NDIM,double> cell_ops(patch_hier);
        side_ops.setToScalar(d_nul_vecs.back()->getComponentDescriptorIndex(0), 0.0);
        cell_ops.setToScalar(d_nul_vecs.back()->getComponentDescriptorIndex(1), 1.0);
    }

    // Initialize the velocity and pressure sub-domain solvers
    const int x_u_idx = sol_vec.getComponentDescriptorIndex(0);
    const int x_p_idx = sol_vec.getComponentDescriptorIndex(1);
    const int b_u_idx = rhs_vec.getComponentDescriptorIndex(0);
    const int b_p_idx = rhs_vec.getComponentDescriptorIndex(1);
    
    Pointer<SideVariable<NDIM,double> > x_u_sc_var = sol_vec.getComponentVariable(0);
    Pointer<CellVariable<NDIM,double> > x_p_cc_var = sol_vec.getComponentVariable(1);
    Pointer<SideVariable<NDIM,double> > b_u_sc_var = rhs_vec.getComponentVariable(0);
    Pointer<CellVariable<NDIM,double> > b_p_cc_var = rhs_vec.getComponentVariable(1);
    
    SAMRAIVectorReal<NDIM,double> x_u_vec(d_object_name+"::x_u_vec",patch_hier,coarsest_ln,finest_ln);
    SAMRAIVectorReal<NDIM,double> b_u_vec(d_object_name+"::b_u_vec",patch_hier,coarsest_ln,finest_ln);
    SAMRAIVectorReal<NDIM,double> x_p_vec(d_object_name+"::x_p_vec",patch_hier,coarsest_ln,finest_ln);
    SAMRAIVectorReal<NDIM,double> b_p_vec(d_object_name+"::b_p_vec",patch_hier,coarsest_ln,finest_ln);

    x_u_vec.addComponent(x_u_sc_var,x_u_idx,sol_vec.getControlVolumeIndex(0));
    b_u_vec.addComponent(b_u_sc_var,b_u_idx,rhs_vec.getControlVolumeIndex(0));
    x_p_vec.addComponent(x_p_cc_var,x_p_idx,sol_vec.getControlVolumeIndex(1)); 
    b_p_vec.addComponent(b_p_cc_var,b_p_idx,rhs_vec.getControlVolumeIndex(1));     

    IBTK::LinearSolver* p_velocity_solver = dynamic_cast<IBTK::LinearSolver*>(d_velocity_solver.getPointer());
    if (p_velocity_solver) 
    {
        p_velocity_solver->setInitialGuessNonzero(false);
        if (has_velocity_nullspace) p_velocity_solver->setNullspace(false, d_U_nul_vecs);
    }

    IBTK::LinearSolver* p_pressure_solver = dynamic_cast<IBTK::LinearSolver*>(d_pressure_solver.getPointer());
    if (p_pressure_solver)  
    {
        p_pressure_solver->setInitialGuessNonzero(false);
        if (has_pressure_nullspace) p_pressure_solver->setNullspace(true);
    }
    
    d_velocity_solver->initializeSolverState(x_u_vec,b_u_vec);
    d_pressure_solver->initializeSolverState(x_p_vec,b_p_vec);

    // Initialize LInv (Stokes solver for the mobility matrix).
    IBTK::LinearSolver* p_stokes_linear_solver = dynamic_cast<IBTK::LinearSolver*>(d_LInv.getPointer());
    if (!p_stokes_linear_solver)
    {
        IBTK::NewtonKrylovSolver* p_stokes_newton_solver = dynamic_cast<IBTK::NewtonKrylovSolver*>(d_LInv.getPointer());
        if (p_stokes_newton_solver) p_stokes_linear_solver = p_stokes_newton_solver->getLinearSolver().getPointer();
    }
    if (p_stokes_linear_solver)
    {
        StaggeredStokesBlockPreconditioner* p_stokes_block_pc = dynamic_cast<StaggeredStokesBlockPreconditioner*>(p_stokes_linear_solver);
        if (!p_stokes_block_pc)
        {
            IBTK::KrylovLinearSolver* p_stokes_krylov_solver = dynamic_cast<IBTK::KrylovLinearSolver*>(p_stokes_linear_solver);
            if (p_stokes_krylov_solver) p_stokes_block_pc = dynamic_cast<StaggeredStokesBlockPreconditioner*>(p_stokes_krylov_solver->getPreconditioner().getPointer());
        }
        if (p_stokes_block_pc)
        {
            p_stokes_block_pc->setPhysicalBcCoefs(d_ins_integrator->getIntermediateVelocityBoundaryConditions(),
                d_ins_integrator->getProjectionBoundaryConditions());
        }
   
        p_stokes_linear_solver->setInitialGuessNonzero(false); // In preconditioner initial guess has to be zero.
        if (has_velocity_nullspace || has_pressure_nullspace) p_stokes_linear_solver->setNullspace(false, d_nul_vecs);
        p_stokes_linear_solver->initializeSolverState(sol_vec, rhs_vec);
    }

    return;
}// initializeStokesSolver

void
KrylovMobilityInverse::reportKSPConvergedReason(
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
KrylovMobilityInverse::initializeKSP()
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
KrylovMobilityInverse::destroyKSP()
{
    int ierr = KSPDestroy(&d_petsc_ksp);  IBTK_CHKERRQ(ierr);
    d_petsc_ksp = PETSC_NULL;
    return;   
}// destroyKSP

void
KrylovMobilityInverse::resetKSPOptions()
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
               (KrylovMobilityInverse::monitorKSP),PETSC_NULL,PETSC_NULL) ; IBTK_CHKERRQ(ierr);
    }
    return;
}// resetKSPOptions

void
KrylovMobilityInverse::resetKSPOperators()
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
    ierr = MatShellSetOperation(d_petsc_mat, MATOP_MULT    , reinterpret_cast<void(*)(void)>(KrylovMobilityInverse::MatVecMult_KMInv   )); IBTK_CHKERRQ(ierr);

    // Reset the configuration of the PETSc KSP object.
    if (d_petsc_ksp)
    {
        ierr = KSPSetOperators(d_petsc_ksp, d_petsc_mat, d_petsc_mat, SAME_PRECONDITIONER); IBTK_CHKERRQ(ierr);
    }
    return;
}// resetKSPOperators

void
KrylovMobilityInverse::resetKSPPC()
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
        ierr = PCShellSetApply(petsc_pc, KrylovMobilityInverse::PCApply_KMInv); IBTK_CHKERRQ(ierr);
    }
    else
    {
        TBOX_ERROR("This statement should not be reached!\n");
    }
    return;
}// resetKSPPC

PetscErrorCode
KrylovMobilityInverse::MatVecMult_KMInv(
    Mat A,
    Vec x,
    Vec y)
{
    int ierr;
    void* p_ctx;
    ierr = MatShellGetContext(A, &p_ctx); IBTK_CHKERRQ(ierr);
    KrylovMobilityInverse* solver = static_cast<KrylovMobilityInverse*>(p_ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(solver);
    TBOX_ASSERT(solver->d_petsc_mat);
#endif

    // Get the components of x and y
    Vec *vx, *vy;
    ierr = IBTK::VecMultiVecGetSubVecs(x,&vx); IBTK_CHKERRQ(ierr);
    ierr = IBTK::VecMultiVecGetSubVecs(y,&vy); IBTK_CHKERRQ(ierr);

    // Temporary vector for regularization.
    Vec W;
    VecDuplicate(vx[0],&W);

    // Get the regulator vector for individual Lagrangian markers.
    const int finest_ln            = solver->d_cib_method->getStructuresLevelNumber();
    Pointer<IBTK::LData> regulator = solver->d_cib_method->getLDataManager()->getLData("regulator", finest_ln);
    Vec regulator_petsc_vec        = regulator->getVec();
    
    // Use homogeneous BCs with Stokes solver in the preconditioner.
    dynamic_cast<IBTK::LinearSolver*>(solver->d_LInv.getPointer())->setHomogeneousBc(true);

    solver->d_samrai_temp[0]->setToScalar(0.0);
    solver->d_S->apply(vx[0],*solver->d_samrai_temp[0]);
    solver->d_LInv->solveSystem(*solver->d_samrai_temp[1],*solver->d_samrai_temp[0]);
    solver->d_J->apply(*solver->d_samrai_temp[1],vy[0]);
    ierr = VecPointwiseMult(W,regulator_petsc_vec,vx[0]); IBTK_CHKERRQ(ierr);
    ierr = VecAXPY(vy[0],solver->d_scale_interp * solver->d_reg_mob_factor,W); IBTK_CHKERRQ(ierr);  
    ierr = VecDestroy(&W); IBTK_CHKERRQ(ierr);  
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y)); IBTK_CHKERRQ(ierr);   
    PetscFunctionReturn(0);

}// MatVecMult_KMInv

// Routine to apply DirectMobility preconditioner 
PetscErrorCode
KrylovMobilityInverse::PCApply_KMInv(
    PC pc,
    Vec x,
    Vec y)
{
    // Here we are solving the equation of the type : Py = x; where P is the preconditioner
    int ierr;
    void* ctx;
    ierr = PCShellGetContext(pc, &ctx); IBTK_CHKERRQ(ierr);
    KrylovMobilityInverse* solver = static_cast<KrylovMobilityInverse*>(ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(solver);
    TBOX_ASSERT(solver->d_DMInv);
#endif      
 
    // Get some constants
    static const double gamma = solver->d_scale_spread;
    static const double beta  = solver->d_scale_interp;

    Vec* vx, *vy;
    ierr = IBTK::VecMultiVecGetSubVecs(x,&vx);            IBTK_CHKERRQ(ierr);
    ierr = IBTK::VecMultiVecGetSubVecs(y,&vy);            IBTK_CHKERRQ(ierr);
    
    // Solve for y using DirectMobilityInverse
    solver->d_DMInv->solveSystem(y,x);
    ierr = VecScale(y,1.0/(gamma*beta));                  IBTK_CHKERRQ(ierr);

    PetscFunctionReturn(0);

}// PCApply_KMInv

// Routine to log output of KrylovMobilityInverse
PetscErrorCode
KrylovMobilityInverse::monitorKSP(
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

    if(it == 0)
    {
        tbox::plog << "\n\n         Residual norms for -KMInv_ksp" << std::endl;
    } 
    
    std::streamsize old_precision = tbox::plog.precision(16);
    tbox::plog << std::scientific << it << " KMInv_KSP " << print_normtype << " resid norm " << (double)rnorm << " true resid norm " << (double)truenorm << 
    " ||r(i)||/||b|| " << (double)(truenorm/bnorm)<< std::endl;
    tbox::plog.precision(old_precision);
    
   return(0);
  
}//monitorKSP


//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
