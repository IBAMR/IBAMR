// ---------------------------------------------------------------------
//
// Copyright (c) 2015 - 2020 by the IBAMR developers
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

#include "ibamr/CIBStrategy.h"
#include "ibamr/IBStrategy.h"
#include "ibamr/INSStaggeredHierarchyIntegrator.h"
#include "ibamr/KrylovMobilitySolver.h"
#include "ibamr/StaggeredStokesBlockPreconditioner.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibamr/StaggeredStokesSolverManager.h"
#include "ibamr/StokesSpecifications.h"

#include "ibtk/CCPoissonSolverManager.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/KrylovLinearSolver.h"
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
#include <memory>
#include <ostream>
#include <string>
#include <utility>
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

KrylovMobilitySolver::KrylovMobilitySolver(std::string object_name,
                                           Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator,
                                           Pointer<CIBStrategy> cib_strategy,
                                           Pointer<Database> input_db,
                                           std::string default_options_prefix,
                                           MPI_Comm petsc_comm)
    : d_object_name(std::move(object_name)),
      d_options_prefix(std::move(default_options_prefix)),
      d_petsc_comm(petsc_comm),
      d_samrai_temp(2, Pointer<SAMRAIVectorReal<NDIM, PetscScalar> >(nullptr)),
      d_ins_integrator(navier_stokes_integrator),
      d_cib_strategy(cib_strategy)
{
    // Get values from the input database.
    if (input_db) getFromInput(input_db);

    // Create the Stokes solver (LInv) for the linear operator.
    // Create databases for setting up LInv solver.
    std::string stokes_solver_type = StaggeredStokesSolverManager::PETSC_KRYLOV_SOLVER;
    Pointer<Database> stokes_solver_db = nullptr;
    if (input_db->keyExists("stokes_solver_type"))
    {
        stokes_solver_type = input_db->getString("stokes_solver_type");
        if (input_db->keyExists("stokes_solver_db"))
        {
            stokes_solver_db = input_db->getDatabase("stokes_solver_db");
        }
    }
    if (!stokes_solver_db)
    {
        stokes_solver_db = new MemoryDatabase("stokes_solver_db");
        stokes_solver_db->putString("ksp_type", "fgmres");
    }

    std::string stokes_precond_type = StaggeredStokesSolverManager::DEFAULT_BLOCK_PRECONDITIONER;
    Pointer<Database> stokes_precond_db = nullptr;
    if (input_db->keyExists("stokes_precond_type"))
    {
        stokes_precond_type = input_db->getString("stokes_precond_type");
        if (input_db->keyExists("stokes_precond_db"))
        {
            stokes_precond_db = input_db->getDatabase("stokes_precond_db");
        }
    }
    if (!stokes_precond_db)
    {
        stokes_precond_db = new MemoryDatabase("stokes_precond_db");
        stokes_precond_db->putInteger("max_iterations", 1);
    }

    std::string velocity_solver_type = IBTK::SCPoissonSolverManager::PETSC_KRYLOV_SOLVER;
    Pointer<Database> velocity_solver_db = nullptr;
    if (input_db->keyExists("velocity_solver_type"))
    {
        velocity_solver_type = input_db->getString("velocity_solver_type");
        if (input_db->keyExists("velocity_solver_db"))
        {
            velocity_solver_db = input_db->getDatabase("velocity_solver_db");
        }
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
    if (input_db->keyExists("velocity_precond_type"))
    {
        velocity_precond_type = input_db->getString("velocity_precond_type");
        if (input_db->keyExists("velocity_precond_db"))
        {
            velocity_precond_db = input_db->getDatabase("velocity_precond_db");
        }
    }
    if (!velocity_precond_db)
    {
        velocity_precond_db = new MemoryDatabase("velocity_precond_db");
        velocity_precond_db->putInteger("max_iterations", 1);
    }

    std::string pressure_solver_type = IBTK::CCPoissonSolverManager::PETSC_KRYLOV_SOLVER;
    ;
    Pointer<Database> pressure_solver_db = nullptr;
    if (input_db->keyExists("pressure_solver_type"))
    {
        pressure_solver_type = input_db->getString("pressure_solver_type");
        if (input_db->keyExists("pressure_solver_db"))
        {
            pressure_solver_db = input_db->getDatabase("pressure_solver_db");
        }
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
    if (input_db->keyExists("pressure_precond_type"))
    {
        pressure_precond_type = input_db->getString("pressure_precond_type");
        if (input_db->keyExists("pressure_precond_db"))
        {
            pressure_precond_db = input_db->getDatabase("pressure_precond_db");
        }
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
                                                                        "KM_LInv_",
                                                                        stokes_precond_type,
                                                                        d_object_name + "::pc_stokes_precond",
                                                                        stokes_precond_db,
                                                                        "KM_LInv_pc_");

    // Create velocity solver.
    d_velocity_solver = IBTK::SCPoissonSolverManager::getManager()->allocateSolver(velocity_solver_type,
                                                                                   d_object_name + "::velocity_solver",
                                                                                   velocity_solver_db,
                                                                                   "KM_LInv_velocity_",
                                                                                   velocity_precond_type,
                                                                                   d_object_name + "::velocity_precond",
                                                                                   velocity_precond_db,
                                                                                   "KM_LInv_velocity_pc_");

    // Create pressure solver.
    d_pressure_solver = IBTK::CCPoissonSolverManager::getManager()->allocateSolver(pressure_solver_type,
                                                                                   d_object_name + "::pressure_solver",
                                                                                   pressure_solver_db,
                                                                                   "KM_LInv_pressure_",
                                                                                   pressure_precond_type,
                                                                                   d_object_name + "::pressure_precond",
                                                                                   pressure_precond_db,
                                                                                   "KM_LInv_pressure_pc_");

    // Register Poisson specification
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
        if (p_stokes_newton_solver)
        {
            p_stokes_linear_solver = p_stokes_newton_solver->getLinearSolver();
        }
    }
    if (p_stokes_linear_solver)
    {
        Pointer<StaggeredStokesBlockPreconditioner> p_stokes_block_pc = p_stokes_linear_solver;
        if (!p_stokes_block_pc)
        {
            Pointer<IBTK::KrylovLinearSolver> p_stokes_krylov_solver = p_stokes_linear_solver;
            if (p_stokes_krylov_solver)
            {
                p_stokes_block_pc = p_stokes_krylov_solver->getPreconditioner();
            }
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

    IBTK_DO_ONCE(t_solve_system = TimerManager::getManager()->getTimer("IBTK::KrylovMobilitySolver::solveSystem()");
                 t_initialize_solver_state =
                     TimerManager::getManager()->getTimer("IBTK::KrylovMobilitySolver::initializeSolverState()");
                 t_deallocate_solver_state =
                     TimerManager::getManager()->getTimer("IBTK::KrylovMobilitySolver::deallocateSolverState()"););
} // KrylovMobilitySolver

KrylovMobilitySolver::~KrylovMobilitySolver()
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
} // ~KrylovMobilitySolver

void
KrylovMobilitySolver::setInterpScale(const double scale_interp)
{
    d_scale_interp = scale_interp;
} // setInterpScaleFactor

void
KrylovMobilitySolver::setSpreadScale(const double scale_spread)
{
    d_scale_spread = scale_spread;
} // setSpreadScaleFactor

void
KrylovMobilitySolver::setRegularizeMobilityScale(const double scale_reg_mob)
{
    d_reg_mob_factor = scale_reg_mob;
} // setRegularizeMobilityFactor

void
KrylovMobilitySolver::setNormalizeSpreadForce(const bool normalize_force)
{
    d_normalize_spread_force = normalize_force;
} // setNormalizeSpreadForce

void
KrylovMobilitySolver::setKSPType(const std::string& ksp_type)
{
    d_ksp_type = ksp_type;
} // setKSPType

void
KrylovMobilitySolver::setOptionsPrefix(const std::string& options_prefix)
{
    d_options_prefix = options_prefix;
} // setOptionsPrefix

const KSP&
KrylovMobilitySolver::getPETScKSP() const
{
    return d_petsc_ksp;
} // getPETScKSP

Pointer<StaggeredStokesSolver>
KrylovMobilitySolver::getStokesSolver() const
{
    return d_LInv;
} // getStokesSolver

void
KrylovMobilitySolver::setVelocityPoissonSpecifications(const PoissonSpecifications& u_problem_coefs)
{
    d_LInv->setVelocityPoissonSpecifications(u_problem_coefs);
    d_velocity_solver->setPoissonSpecifications(u_problem_coefs);
} // setVelocityPoissonSpecifications

void
KrylovMobilitySolver::setSolutionTime(double solution_time)
{
    d_LInv->setSolutionTime(solution_time);
} // setSolutionTime

void
KrylovMobilitySolver::setTimeInterval(double current_time, double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    double dt = new_time - current_time;
    const double half_time = current_time + 0.5 * dt;
    d_LInv->setTimeInterval(current_time, new_time);
    d_velocity_solver->setSolutionTime(new_time);
    d_pressure_solver->setSolutionTime(half_time);
    d_velocity_solver->setTimeInterval(current_time, new_time);
    d_pressure_solver->setTimeInterval(current_time, new_time);
} // setTimeInterval

void
KrylovMobilitySolver::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                         RobinBcCoefStrategy<NDIM>* p_bc_coef)
{
    d_u_bc_coefs = u_bc_coefs;
    d_LInv->setPhysicalBcCoefs(d_u_bc_coefs, p_bc_coef);
    d_velocity_solver->setPhysicalBcCoefs(d_ins_integrator->getIntermediateVelocityBoundaryConditions());
    d_pressure_solver->setPhysicalBcCoef(d_ins_integrator->getProjectionBoundaryConditions());
} // setPhysicalBcCoefs

void
KrylovMobilitySolver::setPhysicalBoundaryHelper(Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper)
{
    d_LInv->setPhysicalBoundaryHelper(bc_helper);
} // setPhysicalBoundaryHelper

bool
KrylovMobilitySolver::solveSystem(Vec x, Vec b)
{
    IBTK_TIMER_START(t_solve_system);

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x, b);

#if !defined(NDEBUG)
    TBOX_ASSERT(d_petsc_ksp);
#endif

    d_petsc_x = x;
    VecCopy(b, d_petsc_b);

    // Solve the system using a PETSc KSP object.
    KSPSolve(d_petsc_ksp, d_petsc_b, d_petsc_x);
    KSPGetIterationNumber(d_petsc_ksp, &d_current_iterations);
    KSPGetResidualNorm(d_petsc_ksp, &d_current_residual_norm);

    // Determine the convergence reason.
    KSPConvergedReason reason;
    KSPGetConvergedReason(d_petsc_ksp, &reason);
    const bool converged = (static_cast<int>(reason) > 0);
    if (d_enable_logging) reportPETScKSPConvergedReason(d_object_name, reason, plog);

    // Deallocate the solver, when necessary.
    d_petsc_x = nullptr;
    if (deallocate_after_solve) deallocateSolverState();

    IBTK_TIMER_STOP(t_solve_system);
    return converged;
} // solveSystem

void
KrylovMobilitySolver::initializeSolverState(Vec x, Vec b)
{
    IBTK_TIMER_START(t_initialize_solver_state);

    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized)
    {
        d_reinitializing_solver = true;
        deallocateSolverState();
    }

    // Get the Eulerian and Lagrangian components.
    Vec *vx, *vb;
    VecNestGetSubVecs(x, nullptr, &vx);
    VecNestGetSubVecs(b, nullptr, &vb);
    Pointer<SAMRAIVectorReal<NDIM, double> > vx0, vb0;

    // Create the RHS Vec to be used in the KSP object.
    VecDuplicate(vb[1], &d_petsc_b);

    // Create the temporary storage for spreading and Stokes solve operation.
    IBTK::PETScSAMRAIVectorReal::getSAMRAIVectorRead(vx[0], &vx0);
    for (int i = 0; i < 2; ++i)
    {
        d_samrai_temp[i] = vx0->cloneVector("");
        d_samrai_temp[i]->allocateVectorData();
    }
    IBTK::PETScSAMRAIVectorReal::restoreSAMRAIVectorRead(vx[0], &vx0);

    // Initialize PETSc KSP
    initializeKSP();

    // Initialize LInv (Stokes solver) required in the mobility matrix.
    IBTK::PETScSAMRAIVectorReal::getSAMRAIVector(vx[0], &vx0);
    IBTK::PETScSAMRAIVectorReal::getSAMRAIVector(vb[0], &vb0);
    initializeStokesSolver(*vx0, *vb0);

    // Get hierarchy information.
    d_hierarchy = vx0->getPatchHierarchy();
    const int u_data_idx = vx0->getComponentDescriptorIndex(0);
    const int coarsest_ln = vx0->getCoarsestLevelNumber();
    const int finest_ln = vx0->getFinestLevelNumber();
    IBTK::PETScSAMRAIVectorReal::restoreSAMRAIVector(vx[0], &vx0);
    IBTK::PETScSAMRAIVectorReal::restoreSAMRAIVector(vb[0], &vb0);

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
} // initializeSolverState

void
KrylovMobilitySolver::deallocateSolverState()
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_deallocate_solver_state);

    // Deallocate the operator and preconditioner states only if we are not
    // re-initializing the solver.
    if (!d_reinitializing_solver)
    {
        d_LInv->deallocateSolverState();
    }

    // Delete the temporary storage for spreading and Stokes solve operation.
    for (int i = 0; i < 2; ++i)
    {
        d_samrai_temp[i]->resetLevels(d_samrai_temp[i]->getCoarsestLevelNumber(),
                                      std::min(d_samrai_temp[i]->getFinestLevelNumber(),
                                               d_samrai_temp[i]->getPatchHierarchy()->getFinestLevelNumber()));
        d_samrai_temp[i]->freeVectorComponents();
        d_samrai_temp[i].setNull();
    }

    VecDestroy(&d_petsc_b);
    d_petsc_x = nullptr;
    d_petsc_b = nullptr;

    // Destroy the KSP solver.
    destroyKSP();

    // Deallocate the interpolation operators.
    d_hier_bdry_fill->deallocateOperatorState();
    d_hier_bdry_fill.setNull();
    d_transaction_comps.clear();

    // Indicate that the solver is NOT initialized
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_solver_state);
} // deallocateSolverState

/////////////////////////////// PRIVATE //////////////////////////////////////

void
KrylovMobilitySolver::getFromInput(Pointer<Database> input_db)
{
    if (input_db->keyExists("options_prefix")) d_options_prefix = input_db->getString("options_prefix");
    if (input_db->keyExists("max_iterations")) d_max_iterations = input_db->getInteger("max_iterations");
    if (input_db->keyExists("abs_residual_tol")) d_abs_residual_tol = input_db->getDouble("abs_residual_tol");
    if (input_db->keyExists("rel_residual_tol")) d_rel_residual_tol = input_db->getDouble("rel_residual_tol");
    if (input_db->keyExists("ksp_type")) d_ksp_type = input_db->getString("ksp_type");
    if (input_db->keyExists("pc_type")) d_pc_type = input_db->getString("pc_type");
    if (input_db->keyExists("initial_guess_nonzero"))
        d_initial_guess_nonzero = input_db->getBool("initial_guess_nonzero");
    if (input_db->keyExists("normalize_pressure")) d_normalize_pressure = input_db->getBool("normalize_pressure");
    if (input_db->keyExists("normalize_velocity")) d_normalize_velocity = input_db->getBool("normalize_velocity");
    if (input_db->keyExists("enable_logging")) d_enable_logging = input_db->getBool("enable_logging");
} // getFromInput

void
KrylovMobilitySolver::initializeStokesSolver(const SAMRAIVectorReal<NDIM, double>& sol_vec,
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

    // Initialize LInv (Stokes solver for the mobility matrix).
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

        // In preconditioner initial guess has to be zero.
        p_stokes_linear_solver->setInitialGuessNonzero(false);
        if (has_velocity_nullspace || has_pressure_nullspace)
        {
            p_stokes_linear_solver->setNullspace(false, d_nul_vecs);
        }
        p_stokes_linear_solver->initializeSolverState(sol_vec, rhs_vec);
    }
} // initializeStokesSolver

void
KrylovMobilitySolver::initializeKSP()
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
KrylovMobilitySolver::destroyKSP()
{
    KSPDestroy(&d_petsc_ksp);
    d_petsc_ksp = nullptr;
} // destroyKSP

void
KrylovMobilitySolver::resetKSPOptions()
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
        KSPMonitorSet(
            d_petsc_ksp,
            reinterpret_cast<PetscErrorCode (*)(KSP, PetscInt, PetscReal, void*)>(KrylovMobilitySolver::monitorKSP),
            nullptr,
            nullptr);
    }
} // resetKSPOptions

void
KrylovMobilitySolver::resetKSPOperators()
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
        d_petsc_mat, MATOP_MULT, reinterpret_cast<void (*)(void)>(KrylovMobilitySolver::MatVecMult_KMInv));

    // Reset the configuration of the PETSc KSP object.
    if (d_petsc_ksp)
    {
        KSPSetOperators(d_petsc_ksp, d_petsc_mat, d_petsc_mat);
        KSPSetReusePreconditioner(d_petsc_ksp, PETSC_TRUE);
    }
} // resetKSPOperators

void
KrylovMobilitySolver::resetKSPPC()
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
        PCShellSetApply(petsc_pc, KrylovMobilitySolver::PCApply_KMInv);
    }
    else
    {
        TBOX_ERROR("This statement should not be reached!\n");
    }
} // resetKSPPC

PetscErrorCode
KrylovMobilitySolver::MatVecMult_KMInv(Mat A, Vec x, Vec y)
{
    void* p_ctx;
    MatShellGetContext(A, &p_ctx);
    auto solver = static_cast<KrylovMobilitySolver*>(p_ctx);
    Pointer<IBStrategy> ib_method_ops = solver->d_cib_strategy;

#if !defined(NDEBUG)
    TBOX_ASSERT(solver);
    TBOX_ASSERT(solver->d_petsc_mat);
    TBOX_ASSERT(ib_method_ops);
#endif

    // Some constants
    static const double gamma = solver->d_scale_spread;
    static const double beta = solver->d_scale_interp;
    static const double delta = solver->d_reg_mob_factor;
    const double half_time = 0.5 * (solver->d_new_time + solver->d_current_time);

    // Use homogeneous BCs with Stokes solver in the preconditioner.
    dynamic_cast<IBTK::LinearSolver*>(solver->d_LInv.getPointer())->setHomogeneousBc(true);

    // Set y:= [J L^-1 S + \delta]x
    // 1) Spread force.
    solver->d_samrai_temp[0]->setToScalar(0.0);
    solver->d_cib_strategy->setConstraintForce(x, half_time, gamma);
    ib_method_ops->spreadForce(solver->d_samrai_temp[0]->getComponentDescriptorIndex(0),
                               nullptr,
                               std::vector<Pointer<RefineSchedule<NDIM> > >(),
                               half_time);
    if (solver->d_normalize_spread_force)
    {
        solver->d_cib_strategy->subtractMeanConstraintForce(
            x, solver->d_samrai_temp[0]->getComponentDescriptorIndex(0), gamma);
    }
    // 2) Solve Stokes system.
    solver->d_LInv->solveSystem(*solver->d_samrai_temp[1], *solver->d_samrai_temp[0]);

    // 3a) Fill velocity ghost cells.
    int u_data_idx = solver->d_samrai_temp[1]->getComponentDescriptorIndex(0);
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

    // 3b) Interpolate velocity
    solver->d_cib_strategy->setInterpolatedVelocityVector(y, half_time);
    ib_method_ops->interpolateVelocity(u_data_idx,
                                       std::vector<Pointer<CoarsenSchedule<NDIM> > >(),
                                       std::vector<Pointer<RefineSchedule<NDIM> > >(),
                                       half_time);
    solver->d_cib_strategy->getInterpolatedVelocity(y, half_time, beta);

    // 4) Regularize mobility.
    if (!MathUtilities<double>::equalEps(delta, 0.0))
    {
        Vec D;
        VecDuplicate(x, &D);
        solver->d_cib_strategy->computeMobilityRegularization(D, x, delta);
        VecAXPY(y, beta, D);
        VecDestroy(&D);
    }

    PetscFunctionReturn(0);
} // MatVecMult_KMInv

// Routine to apply DirectMobility preconditioner
PetscErrorCode KrylovMobilitySolver::PCApply_KMInv(PC /*pc*/, Vec /*x*/, Vec /*y*/)
{
    TBOX_ERROR(
        "KrylovMobilitySolver::PCApply_KMInv(). Shell Preconditioner for KrylovMobilitySolver not implemented.\n");
    PetscFunctionReturn(0);
} // PCApply_KMInv

// Routine to log output of KrylovMobilitySolver
PetscErrorCode
KrylovMobilitySolver::monitorKSP(KSP ksp, int it, PetscReal rnorm, void* /*mctx*/)
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
        tbox::plog << "\n\n         Residual norms for -KMInv_ksp" << std::endl;
    }

    std::streamsize old_precision = tbox::plog.precision(16);
    tbox::plog << std::scientific << it << " KMInv_KSP " << print_normtype << " resid norm " << rnorm
               << " true resid norm " << truenorm << " ||r(i)||/||b|| " << truenorm / bnorm << std::endl;

    tbox::plog.precision(old_precision);

    PetscFunctionReturn(0);
} // monitorKSP

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
