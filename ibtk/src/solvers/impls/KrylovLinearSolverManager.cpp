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

#include "ibtk/KrylovLinearSolver.h"
#include "ibtk/KrylovLinearSolverManager.h"
#include "ibtk/PETScKrylovLinearSolver.h"

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

const std::string KrylovLinearSolverManager::UNDEFINED = "UNDEFINED";
const std::string KrylovLinearSolverManager::DEFAULT = "DEFAULT";
const std::string KrylovLinearSolverManager::PETSC = "PETSC";

KrylovLinearSolverManager* KrylovLinearSolverManager::s_solver_manager_instance = nullptr;
bool KrylovLinearSolverManager::s_registered_callback = false;
unsigned char KrylovLinearSolverManager::s_shutdown_priority = 200;

KrylovLinearSolverManager*
KrylovLinearSolverManager::getManager()
{
    if (!s_solver_manager_instance)
    {
        s_solver_manager_instance = new KrylovLinearSolverManager();
    }
    if (!s_registered_callback)
    {
        ShutdownRegistry::registerShutdownRoutine(freeManager, s_shutdown_priority);
        s_registered_callback = true;
    }
    return s_solver_manager_instance;
} // getManager

void
KrylovLinearSolverManager::freeManager()
{
    delete s_solver_manager_instance;
    s_solver_manager_instance = nullptr;
    return;
} // freeManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

Pointer<KrylovLinearSolver>
KrylovLinearSolverManager::allocateSolver(const std::string& solver_type,
                                          const std::string& solver_object_name,
                                          Pointer<Database> solver_input_db,
                                          const std::string& solver_default_options_prefix) const
{
    auto it = d_solver_maker_map.find(solver_type);
    if (it == d_solver_maker_map.end())
    {
        TBOX_ERROR("KrylovLinearSolverManager::allocateSolver():\n"
                   << "  unrecognized solver type: " << solver_type << "\n");
    }
    return (it->second)(solver_object_name, solver_input_db, solver_default_options_prefix);
} // allocateSolver

void
KrylovLinearSolverManager::registerSolverFactoryFunction(const std::string& solver_type, SolverMaker solver_maker)
{
    if (d_solver_maker_map.find(solver_type) != d_solver_maker_map.end())
    {
        pout << "KrylovLinearSolverManager::registerSolverFactoryFunction():\n"
             << "  NOTICE: overriding initialization function for solver_type = " << solver_type << "\n";
    }
    d_solver_maker_map[solver_type] = solver_maker;
    return;
} // registerSolverFactoryFunction

/////////////////////////////// PROTECTED ////////////////////////////////////

KrylovLinearSolverManager::KrylovLinearSolverManager() : d_solver_maker_map()
{
    registerSolverFactoryFunction(DEFAULT, PETScKrylovLinearSolver::allocate_solver);
    registerSolverFactoryFunction(PETSC, PETScKrylovLinearSolver::allocate_solver);
    return;
} // KrylovLinearSolverManager

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
