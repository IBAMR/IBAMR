// Filename: NewtonKrylovSolverManager.cpp
// Created on 13 Aug 2012 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

#include <stddef.h>
#include <map>
#include <ostream>
#include <string>
#include <utility>

#include "ibtk/NewtonKrylovSolver.h"
#include "ibtk/NewtonKrylovSolverManager.h"
#include "ibtk/PETScNewtonKrylovSolver.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

const std::string NewtonKrylovSolverManager::UNDEFINED = "UNDEFINED";
const std::string NewtonKrylovSolverManager::DEFAULT = "DEFAULT";
const std::string NewtonKrylovSolverManager::PETSC = "PETSC";

NewtonKrylovSolverManager* NewtonKrylovSolverManager::s_solver_manager_instance = NULL;
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
    s_solver_manager_instance = NULL;
    return;
} // freeManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

Pointer<NewtonKrylovSolver>
NewtonKrylovSolverManager::allocateSolver(const std::string& solver_type,
                                          const std::string& solver_object_name,
                                          Pointer<Database> solver_input_db,
                                          const std::string& solver_default_options_prefix) const
{
    std::map<std::string, SolverMaker>::const_iterator it = d_solver_maker_map.find(solver_type);
    if (it == d_solver_maker_map.end())
    {
        TBOX_ERROR("NewtonKrylovSolverManager::allocateSolver():\n"
                   << "  unrecognized solver type: "
                   << solver_type
                   << "\n");
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

NewtonKrylovSolverManager::~NewtonKrylovSolverManager()
{
    // intentionally blank
    return;
} // ~NewtonKrylovSolverManager

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
