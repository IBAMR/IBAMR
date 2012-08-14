// Filename: CCPoissonSolverManager.C
// Created on 13 Aug 2012 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
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

#include "CCPoissonSolverManager.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

CCPoissonSolverManager* CCPoissonSolverManager::s_solver_manager_instance = NULL;
bool CCPoissonSolverManager::s_registered_callback = false;
unsigned char CCPoissonSolverManager::s_shutdown_priority = 200;

CCPoissonSolverManager*
CCPoissonSolverManager::getManager()
{
    if (s_solver_manager_instance == NULL)
    {
        s_solver_manager_instance = new CCPoissonSolverManager();
    }
    if (!s_registered_callback)
    {
        ShutdownRegistry::registerShutdownRoutine(freeManager, s_shutdown_priority);
        s_registered_callback = true;
    }
    return s_solver_manager_instance;
}// getManager

void
CCPoissonSolverManager::freeManager()
{
    if (s_solver_manager_instance) delete s_solver_manager_instance;
    s_solver_manager_instance = NULL;
    return;
}// freeManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

Pointer<PoissonSolver>
CCPoissonSolverManager::allocateSolver(
    const std::string& solver_type,
    const std::string& solver_object_name,
    Pointer<Database> solver_input_db) const
{
    std::map<std::string,SolverMaker>::const_iterator it = d_solver_maker_map.find(solver_type);
    if (it == d_solver_maker_map.end())
    {
        TBOX_ERROR("CCPoissonSolverManager::allocateSolver():\n"
                   << "  unrecognized solver type: " << solver_type << "\n");
    }
    return (it->second)(solver_object_name, solver_input_db);
}// allocateSolver

void
CCPoissonSolverManager::registerSolverFactoryFunction(
    const std::string& solver_type,
    SolverMaker solver_maker)
{
    d_solver_maker_map[solver_type] = solver_maker;
    return;
}// registerSolverFactoryFunction

/////////////////////////////// PROTECTED ////////////////////////////////////

CCPoissonSolverManager::CCPoissonSolverManager()
    : d_solver_maker_map()
{
    return;
}// CCPoissonSolverManager

CCPoissonSolverManager::~CCPoissonSolverManager()
{
    // intentionally blank
    return;
}// ~CCPoissonSolverManager

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
