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

#include "ibtk/NewtonKrylovSolver.h"
#include "ibtk/NewtonKrylovSolverManager.h"
#include "ibtk/PETScNewtonKrylovSolver.h"

#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/Utilities.h"

#include <map>
#include <ostream>
#include <string>
#include <utility>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

const std::string NewtonKrylovSolverManager::UNDEFINED = "UNDEFINED";
const std::string NewtonKrylovSolverManager::DEFAULT = "DEFAULT";
const std::string NewtonKrylovSolverManager::PETSC = "PETSC";

NewtonKrylovSolverManager* NewtonKrylovSolverManager::s_solver_manager_instance = nullptr;
bool NewtonKrylovSolverManager::s_registered_callback = false;
unsigned char NewtonKrylovSolverManager::s_shutdown_priority = 200;

NewtonKrylovSolverManager*
NewtonKrylovSolverManager::getManager()
{
    if (!s_solver_manager_instance)
    {
        s_solver_manager_instance = new NewtonKrylovSolverManager();
    }
    if (!s_registered_callback)
    {
        ShutdownRegistry::registerShutdownRoutine(freeManager, s_shutdown_priority);
        s_registered_callback = true;
    }
    return s_solver_manager_instance;
} // getManager

void
NewtonKrylovSolverManager::freeManager()
{
    delete s_solver_manager_instance;
    s_solver_manager_instance = nullptr;
    return;
} // freeManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

Pointer<NewtonKrylovSolver>
NewtonKrylovSolverManager::allocateSolver(const std::string& solver_type,
                                          const std::string& solver_object_name,
                                          Pointer<Database> solver_input_db,
                                          const std::string& solver_default_options_prefix) const
{
    auto it = d_solver_maker_map.find(solver_type);
    if (it == d_solver_maker_map.end())
    {
        TBOX_ERROR("NewtonKrylovSolverManager::allocateSolver():\n"
                   << "  unrecognized solver type: " << solver_type << "\n");
    }
    return (it->second)(solver_object_name, solver_input_db, solver_default_options_prefix);
} // allocateSolver

void
NewtonKrylovSolverManager::registerSolverFactoryFunction(const std::string& solver_type, SolverMaker solver_maker)
{
    if (d_solver_maker_map.find(solver_type) != d_solver_maker_map.end())
    {
        pout << "NewtonKrylovSolverManager::registerSolverFactoryFunction():\n"
             << "  NOTICE: overriding initialization function for solver_type = " << solver_type << "\n";
    }
    d_solver_maker_map[solver_type] = solver_maker;
    return;
} // registerSolverFactoryFunction

/////////////////////////////// PROTECTED ////////////////////////////////////

NewtonKrylovSolverManager::NewtonKrylovSolverManager() : d_solver_maker_map()
{
    registerSolverFactoryFunction(DEFAULT, PETScNewtonKrylovSolver::allocate_solver);
    registerSolverFactoryFunction(PETSC, PETScNewtonKrylovSolver::allocate_solver);
    return;
} // NewtonKrylovSolverManager

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
