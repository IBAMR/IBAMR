// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2022 by the IBAMR developers
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

#include "ibamr/PETScKrylovStaggeredStokesSolver.h"
#include "ibamr/StaggeredStokesBlockFactorizationPreconditioner.h"
#include "ibamr/StaggeredStokesLevelRelaxationFACOperator.h"
#include "ibamr/StaggeredStokesOperator.h"
#include "ibamr/StaggeredStokesPETScLevelSolver.h"
#include "ibamr/StaggeredStokesProjectionPreconditioner.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibamr/StaggeredStokesSolverManager.h"

#include "ibtk/KrylovLinearSolver.h"
#include "ibtk/LinearOperator.h"
#include "ibtk/LinearSolver.h"

#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/Utilities.h"

#include <map>
#include <ostream>
#include <string>
#include <utility>

#include "ibamr/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

const std::string StaggeredStokesSolverManager::UNDEFINED = "UNDEFINED";
const std::string StaggeredStokesSolverManager::DEFAULT_KRYLOV_SOLVER = "DEFAULT_KRYLOV_SOLVER";
const std::string StaggeredStokesSolverManager::PETSC_KRYLOV_SOLVER = "PETSC_KRYLOV_SOLVER";
const std::string StaggeredStokesSolverManager::DEFAULT_BLOCK_PRECONDITIONER = "DEFAULT_BLOCK_PRECONDITIONER";
const std::string StaggeredStokesSolverManager::BLOCK_FACTORIZATION_PRECONDITIONER =
    "BLOCK_FACTORIZATION_PRECONDITIONER";
const std::string StaggeredStokesSolverManager::PROJECTION_PRECONDITIONER = "PROJECTION_PRECONDITIONER";
const std::string StaggeredStokesSolverManager::DEFAULT_FAC_PRECONDITIONER = "DEFAULT_FAC_PRECONDITIONER";
const std::string StaggeredStokesSolverManager::BOX_RELAXATION_FAC_PRECONDITIONER = "BOX_RELAXATION_FAC_PRECONDITIONER";
const std::string StaggeredStokesSolverManager::LEVEL_RELAXATION_FAC_PRECONDITIONER =
    "LEVEL_RELAXATION_FAC_PRECONDITIONER";
const std::string StaggeredStokesSolverManager::DEFAULT_LEVEL_SOLVER = "DEFAULT_LEVEL_SOLVER";
const std::string StaggeredStokesSolverManager::PETSC_LEVEL_SOLVER = "PETSC_LEVEL_SOLVER";

StaggeredStokesSolverManager* StaggeredStokesSolverManager::s_solver_manager_instance = nullptr;
bool StaggeredStokesSolverManager::s_registered_callback = false;
unsigned char StaggeredStokesSolverManager::s_shutdown_priority = 200;

StaggeredStokesSolverManager*
StaggeredStokesSolverManager::getManager()
{
    if (!s_solver_manager_instance)
    {
        s_solver_manager_instance = new StaggeredStokesSolverManager();
    }
    if (!s_registered_callback)
    {
        ShutdownRegistry::registerShutdownRoutine(freeManager, s_shutdown_priority);
        s_registered_callback = true;
    }
    return s_solver_manager_instance;
} // getManager

void
StaggeredStokesSolverManager::freeManager()
{
    delete s_solver_manager_instance;
    s_solver_manager_instance = nullptr;
    return;
} // freeManager

namespace
{
Pointer<StaggeredStokesSolver>
allocate_petsc_krylov_solver(const std::string& object_name,
                             Pointer<Database> input_db,
                             const std::string& default_options_prefix)
{
    Pointer<PETScKrylovStaggeredStokesSolver> krylov_solver =
        new PETScKrylovStaggeredStokesSolver(object_name, input_db, default_options_prefix);
    krylov_solver->setOperator(new StaggeredStokesOperator(object_name + "::StokesOperator", true, input_db));
    return krylov_solver;
} // allocate_petsc_krylov_solver
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

Pointer<StaggeredStokesSolver>
StaggeredStokesSolverManager::allocateSolver(const std::string& solver_type,
                                             const std::string& solver_object_name,
                                             Pointer<Database> solver_input_db,
                                             const std::string& solver_default_options_prefix) const
{
    auto it = d_solver_maker_map.find(solver_type);
    if (it == d_solver_maker_map.end())
    {
        TBOX_ERROR("CCPoissonSolverManager::allocateSolver():\n"
                   << "  unrecognized solver type: " << solver_type << "\n");
    }
    return (it->second)(solver_object_name, solver_input_db, solver_default_options_prefix);
} // allocateSolver

Pointer<StaggeredStokesSolver>
StaggeredStokesSolverManager::allocateSolver(const std::string& solver_type,
                                             const std::string& solver_object_name,
                                             Pointer<Database> solver_input_db,
                                             const std::string& solver_default_options_prefix,
                                             const std::string& precond_type,
                                             const std::string& precond_object_name,
                                             Pointer<Database> precond_input_db,
                                             const std::string& precond_default_options_prefix,
                                             const std::string& sub_precond_type,
                                             const std::string& sub_precond_object_name,
                                             Pointer<Database> sub_precond_input_db,
                                             const std::string& sub_precond_default_options_prefix) const
{
    Pointer<StaggeredStokesSolver> solver =
        allocateSolver(solver_type, solver_object_name, solver_input_db, solver_default_options_prefix);
    Pointer<KrylovLinearSolver> p_solver = solver;
    if (p_solver && !precond_type.empty())
    {
        Pointer<StaggeredStokesSolver> precond = allocateSolver(precond_type,
                                                                precond_object_name,
                                                                precond_input_db,
                                                                precond_default_options_prefix,
                                                                sub_precond_type,
                                                                sub_precond_object_name,
                                                                sub_precond_input_db,
                                                                sub_precond_default_options_prefix);
        if (precond) p_solver->setPreconditioner(precond);
    }
    return solver;
} // allocateSolver

void
StaggeredStokesSolverManager::registerSolverFactoryFunction(const std::string& solver_type, SolverMaker solver_maker)
{
    if (d_solver_maker_map.find(solver_type) != d_solver_maker_map.end())
    {
        pout << "StaggeredStokesSolverManager::registerSolverFactoryFunction():\n"
             << "  NOTICE: overriding initialization function for solver_type = " << solver_type << "\n";
    }
    d_solver_maker_map[solver_type] = solver_maker;
    return;
} // registerSolverFactoryFunction

/////////////////////////////// PROTECTED ////////////////////////////////////

StaggeredStokesSolverManager::StaggeredStokesSolverManager()
{
    registerSolverFactoryFunction(DEFAULT_KRYLOV_SOLVER, allocate_petsc_krylov_solver);
    registerSolverFactoryFunction(PETSC_KRYLOV_SOLVER, allocate_petsc_krylov_solver);
    registerSolverFactoryFunction(DEFAULT_BLOCK_PRECONDITIONER,
                                  StaggeredStokesProjectionPreconditioner::allocate_solver);
    registerSolverFactoryFunction(BLOCK_FACTORIZATION_PRECONDITIONER,
                                  StaggeredStokesBlockFactorizationPreconditioner::allocate_solver);
    registerSolverFactoryFunction(PROJECTION_PRECONDITIONER, StaggeredStokesProjectionPreconditioner::allocate_solver);
    registerSolverFactoryFunction(DEFAULT_FAC_PRECONDITIONER,
                                  StaggeredStokesLevelRelaxationFACOperator::allocate_solver);
    registerSolverFactoryFunction(LEVEL_RELAXATION_FAC_PRECONDITIONER,
                                  StaggeredStokesLevelRelaxationFACOperator::allocate_solver);
    registerSolverFactoryFunction(DEFAULT_LEVEL_SOLVER, StaggeredStokesPETScLevelSolver::allocate_solver);
    registerSolverFactoryFunction(PETSC_LEVEL_SOLVER, StaggeredStokesPETScLevelSolver::allocate_solver);
    return;
} // StaggeredStokesSolverManager

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
