// Filename: INSStaggeredStokesSolverManager.C
// Created on 16 Aug 2012 by Boyce Griffith
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

#include "INSStaggeredStokesSolverManager.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// To avoid compiler warnings related to redefinition of MPICH_SKIP_MPICXX.
#ifdef MPICH_SKIP_MPICXX
#undef MPICH_SKIP_MPICXX
#endif

// IBAMR INCLUDES
#include <ibamr/INSStaggeredBoxRelaxationFACOperator.h>
#include <ibamr/INSStaggeredPETScLevelSolver.h>
#include <ibamr/INSStaggeredStokesOperator.h>
#include <ibamr/StokesKrylovLinearSolverWrapper.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/KrylovLinearSolverManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

INSStaggeredStokesSolverManager* INSStaggeredStokesSolverManager::s_solver_manager_instance = NULL;
bool INSStaggeredStokesSolverManager::s_registered_callback = false;
unsigned char INSStaggeredStokesSolverManager::s_shutdown_priority = 200;

INSStaggeredStokesSolverManager*
INSStaggeredStokesSolverManager::getManager()
{
    if (s_solver_manager_instance == NULL)
    {
        s_solver_manager_instance = new INSStaggeredStokesSolverManager();
    }
    if (!s_registered_callback)
    {
        ShutdownRegistry::registerShutdownRoutine(freeManager, s_shutdown_priority);
        s_registered_callback = true;
    }
    return s_solver_manager_instance;
}// getManager

void
INSStaggeredStokesSolverManager::freeManager()
{
    if (s_solver_manager_instance) delete s_solver_manager_instance;
    s_solver_manager_instance = NULL;
    return;
}// freeManager

namespace
{
Pointer<StokesSolver>
allocate_default_krylov_solver(
    const std::string& solver_object_name,
    Pointer<Database> solver_input_db)
{
    Pointer<KrylovLinearSolver> krylov_solver = KrylovLinearSolverManager::getManager()->allocateSolver("DEFAULT", solver_object_name, solver_input_db);
    krylov_solver->setOperator(new INSStaggeredStokesOperator(solver_object_name+"::StokesOperator"));
    return new StokesKrylovLinearSolverWrapper(krylov_solver);
}// allocate_default_krylov_solver

Pointer<StokesSolver>
allocate_petsc_krylov_solver(
    const std::string& solver_object_name,
    Pointer<Database> solver_input_db)
{
    Pointer<KrylovLinearSolver> krylov_solver = KrylovLinearSolverManager::getManager()->allocateSolver("PETSC_KRYLOV_LINEAR_SOLVER", solver_object_name, solver_input_db);
    krylov_solver->setOperator(new INSStaggeredStokesOperator(solver_object_name+"::StokesOperator"));
    return new StokesKrylovLinearSolverWrapper(krylov_solver);
}// allocate_petsc_krylov_solver

Pointer<StokesSolver>
allocate_box_relaxation_fac_preconditioner(
    const std::string& solver_object_name,
    Pointer<Database> solver_input_db)
{
    Pointer<StokesFACPreconditionerStrategy> fac_operator = new INSStaggeredBoxRelaxationFACOperator(solver_object_name+"::FACOperator", solver_input_db);
    return new StokesFACPreconditioner(solver_object_name, fac_operator);
}// allocate_box_relaxation_fac_preconditioner
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

Pointer<StokesSolver>
INSStaggeredStokesSolverManager::allocateSolver(
    const std::string& solver_type,
    const std::string& solver_object_name,
    Pointer<Database> solver_input_db) const
{
    std::map<std::string,SolverMaker>::const_iterator it = d_solver_maker_map.find(solver_type);
    if (it == d_solver_maker_map.end())
    {
        TBOX_ERROR("INSStaggeredStokesSolverManager::allocateSolver():\n"
                   << "  unrecognized solver type: " << solver_type << "\n");
    }
    return (it->second)(solver_object_name, solver_input_db);
}// allocateSolver

Pointer<StokesSolver>
INSStaggeredStokesSolverManager::allocateSolver(
    const std::string& solver_type,
    const std::string& solver_object_name,
    Pointer<Database> solver_input_db,
    const std::string& precond_type,
    const std::string& precond_object_name,
    Pointer<Database> precond_input_db) const
{
    Pointer<StokesSolver> solver = allocateSolver(solver_type, solver_object_name, solver_input_db);
    Pointer<KrylovLinearSolver> p_solver = solver;
    if (!p_solver.isNull())
    {
        p_solver->setPreconditioner(allocateSolver(precond_type, precond_object_name, precond_input_db));
    }
    return solver;
}// allocateSolver

void
INSStaggeredStokesSolverManager::registerSolverFactoryFunction(
    const std::string& solver_type,
    SolverMaker solver_maker)
{
    d_solver_maker_map[solver_type] = solver_maker;
    return;
}// registerSolverFactoryFunction

/////////////////////////////// PROTECTED ////////////////////////////////////

INSStaggeredStokesSolverManager::INSStaggeredStokesSolverManager()
    : d_solver_maker_map()
{
    d_solver_maker_map[INSStaggeredPETScLevelSolver::SOLVER_TYPE_NAME] = INSStaggeredPETScLevelSolver::allocate_solver;
    return;
}// INSStaggeredStokesSolverManager

INSStaggeredStokesSolverManager::~INSStaggeredStokesSolverManager()
{
    // intentionally blank
    return;
}// ~INSStaggeredStokesSolverManager

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
