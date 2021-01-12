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

#include "ibamr/CIBMobilitySolver.h"
#include "ibamr/CIBStrategy.h"
#include "ibamr/DirectMobilitySolver.h"
#include "ibamr/INSStaggeredHierarchyIntegrator.h"
#include "ibamr/KrylovFreeBodyMobilitySolver.h"
#include "ibamr/KrylovMobilitySolver.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/ibamr_utilities.h"

#include "PoissonSpecifications.h"
#include "tbox/Database.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_solve_mobility_system;
static Timer* t_solve_body_mobility_system;
static Timer* t_initialize_solver_state;
static Timer* t_deallocate_solver_state;
} // namespace

////////////////////////////// PUBLIC ////////////////////////////////////////

CIBMobilitySolver::CIBMobilitySolver(std::string object_name,
                                     Pointer<Database> input_db,
                                     Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator,
                                     Pointer<CIBStrategy> cib_strategy)
    : d_object_name(std::move(object_name)),
      d_num_rigid_parts(cib_strategy->getNumberOfRigidStructures()),
      d_cib_strategy(cib_strategy)
{
    // Get from input.
    if (input_db) getFromInput(input_db);

    // Create the mobility solvers.
    if (d_mobility_solver_type == KRYLOV)
    {
        d_krylov_mob_solver = new KrylovMobilitySolver(d_object_name + "KrylovMobilitySolver",
                                                       navier_stokes_integrator,
                                                       d_cib_strategy,
                                                       input_db->getDatabase("KrylovMobilitySolver"),
                                                       "KMInv_");
    }
    if (d_mobility_solver_type == DIRECT)
    {
        d_direct_mob_solver = new DirectMobilitySolver(
            d_object_name + "DirectMobilitySolver", input_db->getDatabase("DirectMobilitySolver"), cib_strategy);
        d_direct_mob_solver->setStokesSpecifications(*navier_stokes_integrator->getStokesSpecifications());
    }

    // Create solver for free parts moving under external forces.
    {
        d_krylov_freebody_mob_solver =
            new KrylovFreeBodyMobilitySolver(d_object_name + "FreeBodyMobilitySolver",
                                             input_db->getDatabase("KrylovFreeBodyMobilitySolver"),
                                             "KFBMInv_",
                                             cib_strategy);
        d_krylov_freebody_mob_solver->setMobilitySolver(this);
        d_krylov_freebody_mob_solver->setStokesSpecifications(*navier_stokes_integrator->getStokesSpecifications());
    }

    IBAMR_DO_ONCE(t_solve_mobility_system =
                      TimerManager::getManager()->getTimer("IBAMR::CIBMobilitySolver::solveMobilitySystem()");
                  t_solve_body_mobility_system =
                      TimerManager::getManager()->getTimer("IBAMR::CIBMobilitySolver::solveBodyMobilitySystem()");
                  t_initialize_solver_state =
                      TimerManager::getManager()->getTimer("IBAMR::CIBMobilitySolver::initializeSolverState()");
                  t_deallocate_solver_state =
                      TimerManager::getManager()->getTimer("IBAMR::CIBMobilitySolver::deallocateSolverState()"););

    return;
} // CIBMobilitySolver

CIBMobilitySolver::~CIBMobilitySolver()
{
    return;
} // ~CIBMobilitySolver

void
CIBMobilitySolver::setInterpScale(const double interp_scale)
{
    d_interp_scale = interp_scale;
    if (d_mobility_solver_type == KRYLOV)
    {
        d_krylov_mob_solver->setInterpScale(interp_scale);
    }
    d_krylov_freebody_mob_solver->setInterpScale(interp_scale);

    return;
} // setInterpScale

void
CIBMobilitySolver::setSpreadScale(const double spread_scale)
{
    d_spread_scale = spread_scale;
    if (d_mobility_solver_type == KRYLOV)
    {
        d_krylov_mob_solver->setSpreadScale(spread_scale);
    }
    d_krylov_freebody_mob_solver->setSpreadScale(spread_scale);

    return;
} // setSpreadScale

void
CIBMobilitySolver::setRegularizeMobilityScale(const double reg_mob_scale)
{
    d_reg_mob_scale = reg_mob_scale;
    if (d_mobility_solver_type == KRYLOV)
    {
        d_krylov_mob_solver->setRegularizeMobilityScale(reg_mob_scale);
    }
    return;
} // setRegularizeMobilityScale

void
CIBMobilitySolver::setNormalizeSpreadForce(const bool normalize_spread_force)
{
    if (d_mobility_solver_type == KRYLOV)
    {
        d_krylov_mob_solver->setNormalizeSpreadForce(normalize_spread_force);
    }
    return;
} // setNormalizeSpreadForce

void
CIBMobilitySolver::setSolutionTime(const double solution_time)
{
    d_solution_time = solution_time;

    if (d_mobility_solver_type == KRYLOV)
    {
        d_krylov_mob_solver->setSolutionTime(d_solution_time);
    }
    else if (d_mobility_solver_type == DIRECT)
    {
        d_direct_mob_solver->setSolutionTime(d_solution_time);
    }
    else
    {
        TBOX_ERROR("This statement should not be reached\n");
    }
    d_krylov_freebody_mob_solver->setSolutionTime(solution_time);

    return;
} // setSolutionTime

void
CIBMobilitySolver::setTimeInterval(double current_time, double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;

    if (d_mobility_solver_type == KRYLOV)
    {
        d_krylov_mob_solver->setTimeInterval(d_current_time, d_new_time);
    }
    else if (d_mobility_solver_type == DIRECT)
    {
        d_direct_mob_solver->setTimeInterval(d_current_time, d_new_time);
    }
    else
    {
        TBOX_ERROR("This statement should not be reached\n");
    }
    d_krylov_freebody_mob_solver->setTimeInterval(current_time, new_time);

    return;
} // setTimeInterval

void
CIBMobilitySolver::setVelocityPoissonSpecifications(const PoissonSpecifications& u_problem_coefs)
{
    if (d_mobility_solver_type == KRYLOV)
    {
        d_krylov_mob_solver->setVelocityPoissonSpecifications(u_problem_coefs);
    }
    return;
} // setVelocityPoissonSpecifications

void
CIBMobilitySolver::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                      RobinBcCoefStrategy<NDIM>* p_bc_coef)
{
    if (d_mobility_solver_type == KRYLOV)
    {
        d_krylov_mob_solver->setPhysicalBcCoefs(u_bc_coefs, p_bc_coef);
    }
    return;
} // setPhysicalBcCoefs

void
CIBMobilitySolver::setPhysicalBoundaryHelper(Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper)
{
    if (d_mobility_solver_type == KRYLOV)
    {
        d_krylov_mob_solver->setPhysicalBoundaryHelper(bc_helper);
    }

    return;
} // setPhysicalBoundaryHelper

void
CIBMobilitySolver::getMobilitySolvers(KrylovMobilitySolver** km_solver,
                                      DirectMobilitySolver** dm_solver,
                                      KrylovFreeBodyMobilitySolver** fbm_solver)
{
    if (km_solver)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(d_mobility_solver_type == KRYLOV);
#endif
        *km_solver = d_krylov_mob_solver.getPointer();
    }

    if (dm_solver)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(d_mobility_solver_type == DIRECT);
#endif
        *dm_solver = d_direct_mob_solver.getPointer();
    }

    if (fbm_solver)
    {
        *fbm_solver = d_krylov_freebody_mob_solver.getPointer();
    }

    return;

} // getMobilitySolvers

void
CIBMobilitySolver::initializeSolverState(Vec x, Vec b)
{
    IBAMR_TIMER_START(t_initialize_solver_state);

    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized)
    {
        d_reinitializing_solver = true;
        deallocateSolverState();
    }

    if (d_mobility_solver_type == KRYLOV)
    {
        d_krylov_mob_solver->initializeSolverState(x, b);
    }
    else if (d_mobility_solver_type == DIRECT)
    {
        d_direct_mob_solver->initializeSolverState(x, b);
    }
    else
    {
        TBOX_ERROR("CIBMobilitySolver::initializeSolverState() Unknown mobility solver type" << std::endl);
    }

    d_has_free_parts = false;
    for (unsigned part = 0; part < d_num_rigid_parts && !d_has_free_parts; ++part)
    {
        int num_free_dofs;
        d_cib_strategy->getSolveRigidBodyVelocity(part, num_free_dofs);
        if (num_free_dofs)
        {
            d_has_free_parts = true;
        }
    }
    if (d_has_free_parts)
    {
        d_krylov_freebody_mob_solver->initializeSolverState(x, b);
    }

    // Indicate that the solver is initialized.
    d_reinitializing_solver = false;
    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_solver_state);

    return;
} // initializeSolverState

void
CIBMobilitySolver::deallocateSolverState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_solver_state);

    // Deallocate the operator and preconditioner states only if we are not
    // re-initializing the solver.
    if (!d_reinitializing_solver)
    {
        if (d_mobility_solver_type == KRYLOV)
        {
            d_krylov_mob_solver->deallocateSolverState();
        }
        else if (d_mobility_solver_type == DIRECT)
        {
            d_direct_mob_solver->deallocateSolverState();
        }
        else
        {
            TBOX_ERROR("CIBMobilitySolver::deallocateSolverState() Unknown mobility "
                       << " solver type encountered." << std::endl);
        }

        if (d_has_free_parts)
        {
            d_krylov_freebody_mob_solver->deallocateSolverState();
        }
    }

    // Indicate that the solver is NOT initialized.
    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_solver_state);

    return;
} // deallocateSolverState

bool
CIBMobilitySolver::solveMobilitySystem(Vec x, Vec b)
{
    IBAMR_TIMER_START(t_solve_mobility_system);

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x, b);

    // Solve for x.
    bool converged = false;
    if (d_mobility_solver_type == KRYLOV)
    {
        converged = d_krylov_mob_solver->solveSystem(x, b);
    }
    else if (d_mobility_solver_type == DIRECT)
    {
        converged = d_direct_mob_solver->solveSystem(x, b);
        const double scale = 1.0 / (d_interp_scale * d_spread_scale);
        VecScale(x, scale);
    }
    else
    {
        TBOX_ERROR("This statment should not be reached\n");
    }

    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();

    IBAMR_TIMER_STOP(t_solve_mobility_system);
    return converged;
} // solveMobilitySystem

bool
CIBMobilitySolver::solveBodyMobilitySystem(Vec x, Vec b)
{
    IBAMR_TIMER_START(t_solve_body_mobility_system);

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x, b);

    // Solve for x.
    bool converged = d_krylov_freebody_mob_solver->solveSystem(x, b);

    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();

    IBAMR_TIMER_STOP(t_solve_body_mobility_system);

    return converged;
} // solveBodyMobilitySystem

////////////////////////////// PRIVATE ///////////////////////////////////////

void
CIBMobilitySolver::getFromInput(Pointer<Database> input_db)
{
    // Get the mobility solver type.
    const std::string solver_type = input_db->getString("mobility_solver_type");
    if (solver_type == "DIRECT")
    {
        d_mobility_solver_type = DIRECT;
    }
    else if (solver_type == "KRYLOV")
    {
        d_mobility_solver_type = KRYLOV;
    }
    else
    {
        TBOX_ERROR("CIBMobilitySolver::getFromInput() Unknown mobility solver type = " << solver_type << " provided."
                                                                                       << std::endl);
    }

    return;
} // getFromInput

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR
