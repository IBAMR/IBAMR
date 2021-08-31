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

#include "ibamr/CIBMobilitySolver.h"
#include "ibamr/CIBSaddlePointSolver.h"
#include "ibamr/CIBStaggeredStokesOperator.h"
#include "ibamr/CIBStrategy.h"
#include "ibamr/IBStrategy.h"
#include "ibamr/INSStaggeredHierarchyIntegrator.h"
#include "ibamr/StaggeredStokesBlockPreconditioner.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibamr/StaggeredStokesSolverManager.h"
#include "ibamr/StokesSpecifications.h"

#include "ibtk/CCPoissonSolverManager.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/KrylovLinearSolver.h"
#include "ibtk/LinearOperator.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/NewtonKrylovSolver.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "ibtk/PoissonSolver.h"
#include "ibtk/SCPoissonSolverManager.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/solver_utilities.h"

#include "ArrayData.h"
#include "CellVariable.h"
#include "CoarsenSchedule.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PoissonSpecifications.h"
#include "RefineSchedule.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableFillPattern.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/MemoryDatabase.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

#include "petscksp.h"
#include "petscmat.h"
#include "petscpc.h"
#include "petscpctypes.h"
#include "petscvec.h"
#include <petsclog.h>
#include <petscsys.h>

#include <algorithm>
#include <ostream>
#include <string>
#include <vector>

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE = "NONE";
static const bool USE_CF_INTERPOLATION = true;
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Timers.
static Timer* t_solve_system;
static Timer* t_initialize_solver_state;
static Timer* t_deallocate_solver_state;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

CIBSaddlePointSolver::CIBSaddlePointSolver(std::string object_name,
                                           Pointer<Database> input_db,
                                           Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator,
                                           Pointer<CIBStrategy> cib_strategy,
                                           std::string default_options_prefix,
                                           MPI_Comm petsc_comm)
    : d_object_name(std::move(object_name)),
      d_options_prefix(std::move(default_options_prefix)),
      d_petsc_comm(petsc_comm),
      d_ins_integrator(navier_stokes_integrator),
      d_cib_strategy(cib_strategy),
      d_num_rigid_parts(cib_strategy->getNumberOfRigidStructures())
{
    // Get values from the input database.
    if (input_db) getFromInput(input_db);

    // Create the linear operator for the extended Stokes (Krylov) solver.
    d_A = new CIBStaggeredStokesOperator(d_object_name + "CIBStaggeredStokesOperator",
                                         d_cib_strategy,
                                         /*homogeneous_bc*/ false);
    d_A->setInterpScaleFactor(d_scale_interp);
    d_A->setSpreadScaleFactor(d_scale_spread);
    d_A->setRegularizeMobilityFactor(d_reg_mob_factor);
    d_A->setNormalizeSpreadForce(d_normalize_spread_force);

    // Create the mobility solver and pass-in some parameters.
    d_mob_solver =
        new CIBMobilitySolver(d_object_name + "CIBMobilitySolver", input_db, d_ins_integrator, d_cib_strategy);
    d_mob_solver->setInterpScale(d_scale_interp);
    d_mob_solver->setSpreadScale(d_scale_spread);
    d_mob_solver->setRegularizeMobilityScale(d_reg_mob_factor);
    d_mob_solver->setNormalizeSpreadForce(d_normalize_spread_force);

    // Create the Stokes solver (LInv) for the preconditioner.
    // Create databases for setting up LInv solver.
    Pointer<Database> LInv_db = input_db->getDatabase("PCStokesSolver");
    if (LInv_db->keyExists("normalize_pressure")) d_normalize_pressure = LInv_db->getBool("normalize_pressure");
    if (LInv_db->keyExists("normalize_velocity")) d_normalize_velocity = LInv_db->getBool("normalize_velocity");

    std::string stokes_solver_type = StaggeredStokesSolverManager::PETSC_KRYLOV_SOLVER;
    Pointer<Database> stokes_solver_db = nullptr;
    if (LInv_db->keyExists("stokes_solver_type"))
    {
        stokes_solver_type = LInv_db->getString("stokes_solver_type");
        if (LInv_db->keyExists("stokes_solver_db")) stokes_solver_db = LInv_db->getDatabase("stokes_solver_db");
    }
    if (!stokes_solver_db)
    {
        stokes_solver_db = new MemoryDatabase("stokes_solver_db");
        stokes_solver_db->putString("ksp_type", "fgmres");
    }

    std::string stokes_precond_type = StaggeredStokesSolverManager::DEFAULT_BLOCK_PRECONDITIONER;
    Pointer<Database> stokes_precond_db = nullptr;
    if (LInv_db->keyExists("stokes_precond_type"))
    {
        stokes_precond_type = LInv_db->getString("stokes_precond_type");
        if (LInv_db->keyExists("stokes_precond_db")) stokes_precond_db = LInv_db->getDatabase("stokes_precond_db");
    }
    if (!stokes_precond_db)
    {
        stokes_precond_db = new MemoryDatabase("stokes_precond_db");
        stokes_precond_db->putInteger("max_iterations", 1);
    }

    std::string velocity_solver_type = IBTK::SCPoissonSolverManager::PETSC_KRYLOV_SOLVER;
    Pointer<Database> velocity_solver_db = nullptr;
    if (LInv_db->keyExists("velocity_solver_type"))
    {
        velocity_solver_type = LInv_db->getString("velocity_solver_type");
        if (LInv_db->keyExists("velocity_solver_db")) velocity_solver_db = LInv_db->getDatabase("velocity_solver_db");
    }
    if (!velocity_solver_db)
    {
        velocity_solver_db = new MemoryDatabase("velocity_solver_db");
        velocity_solver_db->putString("ksp_type", "richardson");
        velocity_solver_db->putInteger("max_iterations", 10);
        velocity_solver_db->putDouble("rel_residual_tol", 1.0e-1);
    }

    std::string velocity_precond_type = IBTK::SCPoissonSolverManager::DEFAULT_FAC_PRECONDITIONER;
    Pointer<Database> velocity_precond_db = nullptr;
    if (LInv_db->keyExists("velocity_precond_type"))
    {
        velocity_precond_type = LInv_db->getString("velocity_precond_type");
        if (LInv_db->keyExists("velocity_precond_db"))
            velocity_precond_db = LInv_db->getDatabase("velocity_precond_db");
    }
    if (!velocity_precond_db)
    {
        velocity_precond_db = new MemoryDatabase("velocity_precond_db");
        velocity_precond_db->putInteger("max_iterations", 1);
    }

    std::string pressure_solver_type = IBTK::CCPoissonSolverManager::PETSC_KRYLOV_SOLVER;
    Pointer<Database> pressure_solver_db = nullptr;
    if (LInv_db->keyExists("pressure_solver_type"))
    {
        pressure_solver_type = LInv_db->getString("pressure_solver_type");
        if (LInv_db->keyExists("pressure_solver_db")) pressure_solver_db = LInv_db->getDatabase("pressure_solver_db");
    }
    if (!pressure_solver_db)
    {
        pressure_solver_db = new MemoryDatabase("pressure_solver_db");
        pressure_solver_db->putString("ksp_type", "richardson");
        pressure_solver_db->putInteger("max_iterations", 10);
        pressure_solver_db->putDouble("rel_residual_tol", 1.0e-1);
    }

    std::string pressure_precond_type = IBTK::CCPoissonSolverManager::DEFAULT_FAC_PRECONDITIONER;
    Pointer<Database> pressure_precond_db = nullptr;
    if (LInv_db->keyExists("pressure_precond_type"))
    {
        pressure_precond_type = LInv_db->getString("pressure_precond_type");
        if (LInv_db->keyExists("pressure_precond_db"))
            pressure_precond_db = LInv_db->getDatabase("pressure_precond_db");
    }
    if (!pressure_precond_db)
    {
        pressure_precond_db = new MemoryDatabase("pressure_precond_db");
        pressure_precond_db->putInteger("max_iterations", 1);
    }

    // Create LInv.
    d_LInv = StaggeredStokesSolverManager::getManager()->allocateSolver(stokes_solver_type,
                                                                        d_object_name + "::pc_stokes_solver",
                                                                        stokes_solver_db,
                                                                        "LInv_",
                                                                        stokes_precond_type,
                                                                        d_object_name + "::pc_stokes_precond",
                                                                        stokes_precond_db,
                                                                        "LInv_pc_");

    // Create velocity solver.
    d_velocity_solver = IBTK::SCPoissonSolverManager::getManager()->allocateSolver(velocity_solver_type,
                                                                                   d_object_name + "::velocity_solver",
                                                                                   velocity_solver_db,
                                                                                   "LInv_velocity_",
                                                                                   velocity_precond_type,
                                                                                   d_object_name + "::velocity_precond",
                                                                                   velocity_precond_db,
                                                                                   "LInv_velocity_pc_");

    // Create pressure solver.
    d_pressure_solver = IBTK::CCPoissonSolverManager::getManager()->allocateSolver(pressure_solver_type,
                                                                                   d_object_name + "::pressure_solver",
                                                                                   pressure_solver_db,
                                                                                   "LInv_pressure_",
                                                                                   pressure_precond_type,
                                                                                   d_object_name + "::pressure_precond",
                                                                                   pressure_precond_db,
                                                                                   "LInv_pressure_pc_");

    // Register Poisson specification with pressure sub-domain solver.
    const StokesSpecifications& stokes_spec = *d_ins_integrator->getStokesSpecifications();
    const double rho = stokes_spec.getRho();
    PoissonSpecifications P_problem_coefs(d_object_name + "::P_problem_coefs");
    P_problem_coefs.setCZero();
    P_problem_coefs.setDConstant(rho == 0.0 ? -1.0 : -1.0 / rho);
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
            if (p_stokes_krylov_solver) p_stokes_block_pc = p_stokes_krylov_solver->getPreconditioner();
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

    IBTK_DO_ONCE(t_solve_system = TimerManager::getManager()->getTimer("IBAMR::CIBSaddlePointSolver::solveSystem()");
                 t_initialize_solver_state =
                     TimerManager::getManager()->getTimer("IBAMR::CIBSaddlePointSolver::initializeSolverState()");
                 t_deallocate_solver_state =
                     TimerManager::getManager()->getTimer("IBAMR::CIBSaddlePointSolver::deallocateSolverState()"););
    return;
} // CIBSaddlePointSolver

CIBSaddlePointSolver::~CIBSaddlePointSolver()
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
    return;
} // ~CIBSaddlePointSolver

void
CIBSaddlePointSolver::setKSPType(const std::string& ksp_type)
{
    d_ksp_type = ksp_type;
    return;
} // setKSPType

void
CIBSaddlePointSolver::setOptionsPrefix(const std::string& options_prefix)
{
    d_options_prefix = options_prefix;
    return;
} // setOptionsPrefix

void
CIBSaddlePointSolver::setVelocityPoissonSpecifications(const PoissonSpecifications& u_problem_coefs)
{
    d_A->setVelocityPoissonSpecifications(u_problem_coefs);
    d_LInv->setVelocityPoissonSpecifications(u_problem_coefs);
    d_velocity_solver->setPoissonSpecifications(u_problem_coefs);
    d_mob_solver->setVelocityPoissonSpecifications(u_problem_coefs);

    return;
} // setVelocityPoissonSpecifications

void
CIBSaddlePointSolver::setSolutionTime(double solution_time)
{
    d_A->setSolutionTime(solution_time);
    d_LInv->setSolutionTime(solution_time);
    d_mob_solver->setSolutionTime(solution_time);

    return;
} // setSolutionTime

void
CIBSaddlePointSolver::setTimeInterval(double current_time, double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    const double half_time = 0.5 * (new_time + current_time);

    d_A->setTimeInterval(current_time, new_time);
    d_LInv->setTimeInterval(current_time, new_time);
    d_mob_solver->setTimeInterval(current_time, new_time);

    d_velocity_solver->setSolutionTime(new_time);
    d_pressure_solver->setSolutionTime(half_time);

    d_velocity_solver->setTimeInterval(current_time, new_time);
    d_pressure_solver->setTimeInterval(current_time, new_time);

    return;
} // setTimeInterval

void
CIBSaddlePointSolver::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                         RobinBcCoefStrategy<NDIM>* p_bc_coef)
{
    d_u_bc_coefs = u_bc_coefs;
    d_A->setPhysicalBcCoefs(d_u_bc_coefs, p_bc_coef);
    d_LInv->setPhysicalBcCoefs(d_u_bc_coefs, p_bc_coef);
    d_velocity_solver->setPhysicalBcCoefs(d_ins_integrator->getIntermediateVelocityBoundaryConditions());
    d_pressure_solver->setPhysicalBcCoef(d_ins_integrator->getProjectionBoundaryConditions());
    d_mob_solver->setPhysicalBcCoefs(d_u_bc_coefs, p_bc_coef);

    return;
} // setPhysicalBcCoefs

void
CIBSaddlePointSolver::setPhysicalBoundaryHelper(Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper)
{
    d_A->setPhysicalBoundaryHelper(bc_helper);
    d_LInv->setPhysicalBoundaryHelper(bc_helper);
    d_mob_solver->setPhysicalBoundaryHelper(bc_helper);

    return;
} // setPhysicalBoundaryHelper

Pointer<IBTK::LinearOperator>
CIBSaddlePointSolver::getA() const
{
    return d_A;
} // getA

Pointer<StaggeredStokesSolver>
CIBSaddlePointSolver::getStokesSolver() const
{
    return d_LInv;
} // getStokesSolver

Pointer<CIBMobilitySolver>
CIBSaddlePointSolver::getCIBMobilitySolver() const
{
    return d_mob_solver;

} // getCIBMobilitySolver

bool
CIBSaddlePointSolver::solveSystem(Vec x, Vec b)
{
    IBTK_TIMER_START(t_solve_system);

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x, b);

    d_petsc_x = x;
    VecCopy(b, d_petsc_b);

    // Modify RHS for inhomogeneous Bcs.
    d_A->setHomogeneousBc(false);
    d_A->modifyRhsForBcs(d_petsc_b);
    d_A->setHomogeneousBc(true);

    // Solve the system.
    KSPSolve(d_petsc_ksp, d_petsc_b, d_petsc_x);
    KSPGetIterationNumber(d_petsc_ksp, &d_current_iterations);
    KSPGetResidualNorm(d_petsc_ksp, &d_current_residual_norm);

    // Impose Solution Bcs.
    d_A->setHomogeneousBc(false);
    d_A->imposeSolBcs(d_petsc_x);

    // Determine the convergence reason.
    KSPConvergedReason reason;
    KSPGetConvergedReason(d_petsc_ksp, &reason);
    const bool converged = (static_cast<int>(reason) > 0);
    if (d_enable_logging) reportPETScKSPConvergedReason(d_object_name, reason, plog);

    // Invalidate d_petsc_x Vec.
    d_petsc_x = nullptr;

    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();

    IBTK_TIMER_STOP(t_solve_system);
    return converged;
} // solveSystem

void
CIBSaddlePointSolver::initializeSolverState(Vec x, Vec b)
{
    IBTK_TIMER_START(t_initialize_solver_state);

    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized)
    {
        d_reinitializing_solver = true;
        deallocateSolverState();
    }

    // Create RHS Vec to be used in KSP object
    VecDuplicate(b, &d_petsc_b);

    // Get components of x and b.
    Vec *vx, *vb;
    VecNestGetSubVecs(x, nullptr, &vx);
    VecNestGetSubVecs(b, nullptr, &vb);

    Pointer<SAMRAIVectorReal<NDIM, double> > vx0, vb0;
    IBTK::PETScSAMRAIVectorReal::getSAMRAIVectorRead(vx[0], &vx0);
    IBTK::PETScSAMRAIVectorReal::getSAMRAIVectorRead(vb[0], &vb0);
    d_hierarchy = vx0->getPatchHierarchy();
    const int u_data_idx = vx0->getComponentDescriptorIndex(0);
    const int coarsest_ln = vx0->getCoarsestLevelNumber();
    const int finest_ln = vx0->getFinestLevelNumber();

    // Initialize various operators and solvers.
    d_A->initializeOperatorState(*vx0, *vb0);
    initializeStokesSolver(*vx0, *vb0);

    // Restore the SAMRAI vectors.
    IBTK::PETScSAMRAIVectorReal::restoreSAMRAIVectorRead(vx[0], &vx0);
    IBTK::PETScSAMRAIVectorReal::restoreSAMRAIVectorRead(vb[0], &vb0);

    // Initialize the mobility solver.
    d_mob_solver->initializeSolverState(x, b);

    // Initialize the KSP.
    initializeKSP();

    // Setup the interpolation transaction information.
    d_fill_pattern = nullptr;
    using InterpolationTransactionComponent = IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    InterpolationTransactionComponent component(u_data_idx,
                                                DATA_REFINE_TYPE,
                                                USE_CF_INTERPOLATION,
                                                DATA_COARSEN_TYPE,
                                                BDRY_EXTRAP_TYPE,
                                                CONSISTENT_TYPE_2_BDRY,
                                                d_u_bc_coefs,
                                                d_fill_pattern);
    d_transaction_comps.push_back(component);

    // Initialize the interpolation operators.
    d_hier_bdry_fill = new IBTK::HierarchyGhostCellInterpolation();
    d_hier_bdry_fill->initializeOperatorState(d_transaction_comps, d_hierarchy, coarsest_ln, finest_ln);

    // Indicate that the solver is initialized.
    d_reinitializing_solver = false;
    d_is_initialized = true;

    IBTK_TIMER_STOP(t_initialize_solver_state);

    return;
} // initializeSolverState

void
CIBSaddlePointSolver::deallocateSolverState()
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_deallocate_solver_state);

    // Deallocate the operator and preconditioner states only if we are not
    // re-initializing the solver.
    if (!d_reinitializing_solver)
    {
        d_A->deallocateOperatorState();
        d_LInv->deallocateSolverState();
        d_mob_solver->deallocateSolverState();
    }

    // Delete the solution and RHS vectors.
    VecDestroy(&d_petsc_b);
    d_petsc_x = nullptr;
    d_petsc_b = nullptr;

    // Destroy the KSP solver.
    destroyKSP();

    // Deallocate the interpolation operators.
    d_hier_bdry_fill->deallocateOperatorState();
    d_hier_bdry_fill.setNull();
    d_transaction_comps.clear();
    d_fill_pattern.setNull();

    // Indicate that the solver is NOT initialized.
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_solver_state);
    return;
} // deallocateSolverState

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CIBSaddlePointSolver::getFromInput(Pointer<Database> input_db)
{
    // Get some options.
    if (input_db->keyExists("options_prefix")) d_options_prefix = input_db->getString("options_prefix");
    if (input_db->keyExists("max_iterations")) d_max_iterations = input_db->getInteger("max_iterations");
    if (input_db->keyExists("abs_residual_tol")) d_abs_residual_tol = input_db->getDouble("abs_residual_tol");
    if (input_db->keyExists("rel_residual_tol")) d_rel_residual_tol = input_db->getDouble("rel_residual_tol");
    if (input_db->keyExists("ksp_type")) d_ksp_type = input_db->getString("ksp_type");
    if (input_db->keyExists("pc_type")) d_pc_type = input_db->getString("pc_type");
    if (input_db->keyExists("initial_guess_nonzero"))
        d_initial_guess_nonzero = input_db->getBool("initial_guess_nonzero");
    if (input_db->keyExists("enable_logging")) d_enable_logging = input_db->getBool("enable_logging");
    if (input_db->keyExists("scale_interp_operator")) d_scale_interp = input_db->getDouble("scale_interp_operator");
    if (input_db->keyExists("scale_spread_operator")) d_scale_spread = input_db->getDouble("scale_spread_operator");
    if (input_db->keyExists("regularize_mob_factor")) d_reg_mob_factor = input_db->getDouble("regularize_mob_factor");
    if (input_db->keyExists("normalize_spread_force"))
        d_normalize_spread_force = input_db->getBool("normalize_spread_force");

    return;
} // getFromInput

void
CIBSaddlePointSolver::initializeKSP()
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

    return;
} // initializeKSP

void
CIBSaddlePointSolver::initializeStokesSolver(const SAMRAIVectorReal<NDIM, double>& sol_vec,
                                             const SAMRAIVectorReal<NDIM, double>& rhs_vec)
{
    Pointer<PatchHierarchy<NDIM> > patch_hier = sol_vec.getPatchHierarchy();
    const int coarsest_ln = sol_vec.getCoarsestLevelNumber();
    const int finest_ln = sol_vec.getFinestLevelNumber();

    // Set the nullspace of the LInv and subdomain solvers
    const double rho = d_ins_integrator->getStokesSpecifications()->getRho();
    const bool has_velocity_nullspace = d_normalize_velocity && MathUtilities<double>::equalEps(rho, 0.0);
    const bool has_pressure_nullspace = d_normalize_pressure;

    for (const auto& nul_vec : d_nul_vecs)
    {
        if (nul_vec) nul_vec->freeVectorComponents();
    }
    const int n_nul_vecs = (has_pressure_nullspace ? 1 : 0) + (has_velocity_nullspace ? NDIM : 0);
    d_nul_vecs.resize(n_nul_vecs);

    for (const auto& U_nul_vec : d_U_nul_vecs)
    {
        if (U_nul_vec) U_nul_vec->freeVectorComponents();
    }
    const int n_U_nul_vecs = (has_velocity_nullspace ? NDIM : 0);
    d_U_nul_vecs.resize(n_U_nul_vecs);

    if (has_velocity_nullspace)
    {
        for (unsigned int k = 0; k < NDIM; ++k)
        {
            d_nul_vecs[k] = sol_vec.cloneVector(d_object_name + "::nul_vec_U_" + std::to_string(k));
            d_nul_vecs[k]->allocateVectorData(d_current_time);
            d_nul_vecs[k]->setToScalar(0.0);

            SAMRAIVectorReal<NDIM, double> svr_u(
                d_object_name + "::U_nul_vec_U_" + std::to_string(k), patch_hier, coarsest_ln, finest_ln);
            svr_u.addComponent(sol_vec.getComponentVariable(0),
                               sol_vec.getComponentDescriptorIndex(0),
                               sol_vec.getControlVolumeIndex(0));

            d_U_nul_vecs[k] = svr_u.cloneVector(svr_u.getName());
            d_U_nul_vecs[k]->allocateVectorData(d_current_time);
            d_U_nul_vecs[k]->setToScalar(0.0);
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = patch_hier->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    Pointer<SideData<NDIM, double> > nul_data =
                        patch->getPatchData(d_nul_vecs[k]->getComponentDescriptorIndex(0));
                    nul_data->getArrayData(k).fillAll(1.0);
                    Pointer<SideData<NDIM, double> > U_nul_data =
                        patch->getPatchData(d_U_nul_vecs[k]->getComponentDescriptorIndex(0));
                    U_nul_data->getArrayData(k).fillAll(1.0);
                }
            }
        }
    }

    if (has_pressure_nullspace)
    {
        d_nul_vecs.back() = sol_vec.cloneVector(d_object_name + "::nul_vec_p");
        d_nul_vecs.back()->allocateVectorData(d_current_time);

        HierarchySideDataOpsReal<NDIM, double> side_ops(patch_hier);
        HierarchyCellDataOpsReal<NDIM, double> cell_ops(patch_hier);
        side_ops.setToScalar(d_nul_vecs.back()->getComponentDescriptorIndex(0), 0.0);
        cell_ops.setToScalar(d_nul_vecs.back()->getComponentDescriptorIndex(1), 1.0);
    }

    // Initialize the velocity and pressure sub-domain solvers
    const int x_u_idx = sol_vec.getComponentDescriptorIndex(0);
    const int x_p_idx = sol_vec.getComponentDescriptorIndex(1);
    const int b_u_idx = rhs_vec.getComponentDescriptorIndex(0);
    const int b_p_idx = rhs_vec.getComponentDescriptorIndex(1);

    Pointer<SideVariable<NDIM, double> > x_u_sc_var = sol_vec.getComponentVariable(0);
    Pointer<CellVariable<NDIM, double> > x_p_cc_var = sol_vec.getComponentVariable(1);
    Pointer<SideVariable<NDIM, double> > b_u_sc_var = rhs_vec.getComponentVariable(0);
    Pointer<CellVariable<NDIM, double> > b_p_cc_var = rhs_vec.getComponentVariable(1);

    SAMRAIVectorReal<NDIM, double> x_u_vec(d_object_name + "::x_u_vec", patch_hier, coarsest_ln, finest_ln);
    SAMRAIVectorReal<NDIM, double> b_u_vec(d_object_name + "::b_u_vec", patch_hier, coarsest_ln, finest_ln);
    SAMRAIVectorReal<NDIM, double> x_p_vec(d_object_name + "::x_p_vec", patch_hier, coarsest_ln, finest_ln);
    SAMRAIVectorReal<NDIM, double> b_p_vec(d_object_name + "::b_p_vec", patch_hier, coarsest_ln, finest_ln);

    x_u_vec.addComponent(x_u_sc_var, x_u_idx, sol_vec.getControlVolumeIndex(0));
    b_u_vec.addComponent(b_u_sc_var, b_u_idx, rhs_vec.getControlVolumeIndex(0));
    x_p_vec.addComponent(x_p_cc_var, x_p_idx, sol_vec.getControlVolumeIndex(1));
    b_p_vec.addComponent(b_p_cc_var, b_p_idx, rhs_vec.getControlVolumeIndex(1));

    auto p_velocity_solver = dynamic_cast<IBTK::LinearSolver*>(d_velocity_solver.getPointer());
    if (p_velocity_solver)
    {
        p_velocity_solver->setInitialGuessNonzero(false);
        if (has_velocity_nullspace) p_velocity_solver->setNullspace(false, d_U_nul_vecs);
    }

    auto p_pressure_solver = dynamic_cast<IBTK::LinearSolver*>(d_pressure_solver.getPointer());
    if (p_pressure_solver)
    {
        p_pressure_solver->setInitialGuessNonzero(false);
        if (has_pressure_nullspace) p_pressure_solver->setNullspace(true);
    }

    d_velocity_solver->initializeSolverState(x_u_vec, b_u_vec);
    d_pressure_solver->initializeSolverState(x_p_vec, b_p_vec);

    // Initialize LInv (Stokes solver for the preconditioner).
    auto p_stokes_linear_solver = dynamic_cast<IBTK::LinearSolver*>(d_LInv.getPointer());
    if (!p_stokes_linear_solver)
    {
        auto p_stokes_newton_solver = dynamic_cast<IBTK::NewtonKrylovSolver*>(d_LInv.getPointer());
        if (p_stokes_newton_solver)
        {
            p_stokes_linear_solver = p_stokes_newton_solver->getLinearSolver().getPointer();
        }
    }
    if (p_stokes_linear_solver)
    {
        auto p_stokes_block_pc = dynamic_cast<StaggeredStokesBlockPreconditioner*>(p_stokes_linear_solver);
        if (!p_stokes_block_pc)
        {
            auto p_stokes_krylov_solver = dynamic_cast<IBTK::KrylovLinearSolver*>(p_stokes_linear_solver);
            if (p_stokes_krylov_solver)
            {
                p_stokes_block_pc = dynamic_cast<StaggeredStokesBlockPreconditioner*>(
                    p_stokes_krylov_solver->getPreconditioner().getPointer());
            }
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
} // initializeStokesSolver

void
CIBSaddlePointSolver::destroyKSP()
{
    int ierr = KSPDestroy(&d_petsc_ksp);
    IBTK_CHKERRQ(ierr);
    d_petsc_ksp = nullptr;
    return;
} // destroyKSP

void
CIBSaddlePointSolver::resetKSPOptions()
{
    if (!d_petsc_ksp) return;
    const KSPType ksp_type = d_ksp_type.c_str();
    KSPSetType(d_petsc_ksp, ksp_type);
    std::string ksp_type_name(ksp_type);
    if (ksp_type_name.find("gmres") != std::string::npos)
    {
        KSPGMRESSetOrthogonalization(d_petsc_ksp, KSPGMRESModifiedGramSchmidtOrthogonalization);
        KSPGMRESSetCGSRefinementType(d_petsc_ksp, KSP_GMRES_CGS_REFINE_ALWAYS);
    }

    PetscBool initial_guess_nonzero = (d_initial_guess_nonzero ? PETSC_TRUE : PETSC_FALSE);
    KSPSetInitialGuessNonzero(d_petsc_ksp, initial_guess_nonzero);
    KSPSetTolerances(d_petsc_ksp, d_rel_residual_tol, d_abs_residual_tol, PETSC_DEFAULT, d_max_iterations);

    // Set KSP monitor routine.
    if (d_enable_logging)
    {
        KSPMonitorCancel(d_petsc_ksp);
        KSPMonitorSet(
            d_petsc_ksp,
            reinterpret_cast<PetscErrorCode (*)(KSP, PetscInt, PetscReal, void*)>(CIBSaddlePointSolver::monitorKSP),
            nullptr,
            nullptr);
    }
    return;
} // resetKSPOptions

void
CIBSaddlePointSolver::resetKSPOperators()
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
        d_petsc_mat, MATOP_MULT, reinterpret_cast<void (*)(void)>(CIBSaddlePointSolver::MatVecMult_SaddlePoint));

    // Reset the configuration of the PETSc KSP object.
    if (d_petsc_ksp)
    {
        KSPSetOperators(d_petsc_ksp, d_petsc_mat, d_petsc_mat);
        KSPSetReusePreconditioner(d_petsc_ksp, PETSC_TRUE);
    }
    return;
} // resetKSPOperators

void
CIBSaddlePointSolver::resetKSPPC()
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
        PCShellSetApply(petsc_pc, CIBSaddlePointSolver::PCApply_SaddlePoint);
    }
    else
    {
        TBOX_ERROR("CIBSaddlePointSolver::resetKSPPC() This statement should not be reached!\n");
    }
    return;
} // resetKSPPC

PetscErrorCode
CIBSaddlePointSolver::MatVecMult_SaddlePoint(Mat A, Vec x, Vec y)
{
    void* p_ctx;
    MatShellGetContext(A, &p_ctx);
    auto solver = static_cast<CIBSaddlePointSolver*>(p_ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(solver);
    TBOX_ASSERT(solver->d_A);
#endif
    solver->d_A->apply(x, y);
    PetscFunctionReturn(0);
} // MatVecMult_SaddlePoint

// Exact-Schur Complement PC
PetscErrorCode
CIBSaddlePointSolver::PCApply_SaddlePoint(PC pc, Vec x, Vec y)
{
    // Here we are solving the equation of the type : Py = x
    // in which P is the preconditioner.
    void* ctx;
    PCShellGetContext(pc, &ctx);
    auto solver = static_cast<CIBSaddlePointSolver*>(ctx);
    Pointer<IBStrategy> ib_method_ops = solver->d_cib_strategy;

#if !defined(NDEBUG)
    TBOX_ASSERT(solver);
    TBOX_ASSERT(solver->d_LInv);
    TBOX_ASSERT(solver->d_mob_solver);
    TBOX_ASSERT(ib_method_ops);
#endif

    // Get some constants
    static const double gamma = solver->d_scale_spread;
    static const double beta = solver->d_scale_interp;
    const double half_time = 0.5 * (solver->d_new_time + solver->d_current_time);

    int total_comps, free_comps = 0;
    Vec *vx, *vy;
    VecNestGetSubVecs(x, nullptr, &vx);
    VecNestGetSubVecs(y, &total_comps, &vy);
    VecGetSize(vx[2], &free_comps);

    Pointer<SAMRAIVectorReal<NDIM, double> > vx0, vy0;
    IBTK::PETScSAMRAIVectorReal::getSAMRAIVectorRead(vx[0], &vx0);
    IBTK::PETScSAMRAIVectorReal::getSAMRAIVector(vy[0], &vy0);

    // Get the individual components.
    Pointer<SAMRAIVectorReal<NDIM, double> > g_h = vx0->cloneVector("");
    g_h->allocateVectorData();
    g_h->copyVector(vx0);
    Pointer<SAMRAIVectorReal<NDIM, double> > u_p = vy0;

    Vec W, Lambda, F_tilde;
    W = vx[1];
    Lambda = vy[1];

    // Create temporary vectors for storage.
    // U is the interpolated velocity and delU is the slip velocity.
    Vec U, delU;
    VecDuplicate(W, &U);

    // 1) (u,p) = L^-1 (g,h)
    dynamic_cast<IBTK::LinearSolver*>(solver->d_LInv.getPointer())->setInitialGuessNonzero(false);
    dynamic_cast<IBTK::LinearSolver*>(solver->d_LInv.getPointer())->setHomogeneousBc(true);
    solver->d_LInv->solveSystem(*u_p, *g_h);

    // 2a) Fill ghost cells of u
    int u_data_idx = u_p->getComponentDescriptorIndex(0);
    using InterpolationTransactionComponent = IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps;
    InterpolationTransactionComponent u_component(u_data_idx,
                                                  DATA_REFINE_TYPE,
                                                  USE_CF_INTERPOLATION,
                                                  DATA_COARSEN_TYPE,
                                                  BDRY_EXTRAP_TYPE,
                                                  CONSISTENT_TYPE_2_BDRY,
                                                  solver->d_u_bc_coefs,
                                                  solver->d_fill_pattern);
    transaction_comps.push_back(u_component);
    solver->d_hier_bdry_fill->resetTransactionComponents(transaction_comps);
    static const bool homogeneous_bc = true;
    solver->d_hier_bdry_fill->setHomogeneousBc(homogeneous_bc);
    solver->d_hier_bdry_fill->fillData(half_time);
    solver->d_hier_bdry_fill->resetTransactionComponents(solver->d_transaction_comps);

    // 2b) U = J u + W.
    solver->d_cib_strategy->setInterpolatedVelocityVector(U, half_time);
    ib_method_ops->interpolateVelocity(u_data_idx,
                                       std::vector<Pointer<CoarsenSchedule<NDIM> > >(),
                                       std::vector<Pointer<RefineSchedule<NDIM> > >(),
                                       half_time);

    solver->d_cib_strategy->getInterpolatedVelocity(U, half_time, beta);
    VecAXPY(U, 1.0, W);

    // 3) Calculate the slip velocity
    if (free_comps)
    {
        VecDuplicate(vx[2], &F_tilde);
        VecDuplicate(W, &delU);

        // 3a) lambda = M^-1(U)
        solver->d_mob_solver->solveMobilitySystem(Lambda, U);

        // 3b) F_tilde = F + T(lambda)
        solver->d_cib_strategy->computeNetRigidGeneralizedForce(Lambda,
                                                                F_tilde,
                                                                /*only_free_dofs*/ true,
                                                                /*only_imposed_dofs*/ false);
        VecAXPY(F_tilde, 1.0, vx[2]);

        // 3c) U_rigid = N^-1(F_tilde)
        solver->d_mob_solver->solveBodyMobilitySystem(vy[2], F_tilde);

        // 3d) delU = T*(U_rigid) - U
        VecSet(delU, 0.0);
        solver->d_cib_strategy->setRigidBodyVelocity(vy[2], delU, /*only_free_dofs*/ true, /*only_imposed_dofs*/ false);
        VecAXPY(delU, -1.0, U);
    }
    else
    {
        delU = U;
        VecScale(delU, -1.0);
    }

    // 4) lambda  = M^-1(delta_U)
    solver->d_mob_solver->solveMobilitySystem(Lambda, delU);

    // 5) (u,p)   = L^-1(S[lambda]+g, h)
    const int g_data_idx = g_h->getComponentDescriptorIndex(0);
    solver->d_cib_strategy->setConstraintForce(Lambda, half_time, gamma);
    ib_method_ops->spreadForce(g_data_idx, nullptr, std::vector<Pointer<RefineSchedule<NDIM> > >(), half_time);
    if (solver->d_normalize_spread_force)
    {
        solver->d_cib_strategy->subtractMeanConstraintForce(Lambda, g_data_idx, gamma);
    }

    // Solve using previous u_p as a guess for Krylov solvers.
    if (solver->d_ksp_type == "preonly")
    {
        dynamic_cast<IBTK::LinearSolver*>(solver->d_LInv.getPointer())->setInitialGuessNonzero(false);
    }
    else
    {
        dynamic_cast<IBTK::LinearSolver*>(solver->d_LInv.getPointer())->setInitialGuessNonzero(true);
    }
    dynamic_cast<IBTK::LinearSolver*>(solver->d_LInv.getPointer())->setHomogeneousBc(true);
    solver->d_LInv->solveSystem(*u_p, *g_h);

    // Destroy temporary vectors
    g_h->resetLevels(g_h->getCoarsestLevelNumber(),
                     std::min(g_h->getFinestLevelNumber(), g_h->getPatchHierarchy()->getFinestLevelNumber()));
    g_h->freeVectorComponents();
    g_h.setNull();
    VecDestroy(&U);
    if (free_comps)
    {
        VecDestroy(&delU);
        VecDestroy(&F_tilde);
    }

    IBTK::PETScSAMRAIVectorReal::restoreSAMRAIVectorRead(vx[0], &vx0);
    IBTK::PETScSAMRAIVectorReal::restoreSAMRAIVector(vy[0], &vy0);

    PetscFunctionReturn(0);

} // PCApply_SaddlePoint

// Routine to log output of CIBSaddlePointSolver
PetscErrorCode
CIBSaddlePointSolver::monitorKSP(KSP ksp, int it, PetscReal rnorm, void* /*mctx*/)
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
        tbox::plog << "\n\n         Residual norms for -SP_ksp" << std::endl;
    }

    std::streamsize old_precision = tbox::plog.precision(16);
    tbox::plog << std::scientific << it << " SP_KSP " << print_normtype << " resid norm " << rnorm
               << " true resid norm " << truenorm << " ||r(i)||/||b|| " << truenorm / bnorm << std::endl;
    tbox::plog.precision(old_precision);

    PetscFunctionReturn(0);

} // monitorKSP
//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
