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

// To avoid compiler warnings related to redefinition of MPICH_SKIP_MPICXX.
#ifdef MPICH_SKIP_MPICXX
#undef MPICH_SKIP_MPICXX
#endif

// IBTK INCLUDES
#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/CCPoissonHypreLevelSolver.h>
#include <ibtk/CCPoissonPETScLevelSolver.h>
#include <ibtk/CCPoissonPointRelaxationFACOperator.h>
#include <ibtk/KrylovLinearSolverManager.h>
#include <ibtk/PoissonKrylovLinearSolverWrapper.h>
#include <ibtk/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

const std::string CCPoissonSolverManager::DEFAULT_KRYLOV_SOLVER               = "DEFAULT_KRYLOV_SOLVER";
const std::string CCPoissonSolverManager::PETSC_KRYLOV_SOLVER                 = "PETSC_KRYLOV_SOLVER";
const std::string CCPoissonSolverManager::DEFAULT_FAC_PRECONDITIONER          = "DEFAULT_FAC_PRECONDITIONER";
const std::string CCPoissonSolverManager::POINT_RELAXATION_FAC_PRECONDITIONER = "POINT_RELAXATION_FAC_PRECONDITIONER";
const std::string CCPoissonSolverManager::DEFAULT_LEVEL_SOLVER                = "DEFAULT_LEVEL_SOLVER";
const std::string CCPoissonSolverManager::HYPRE_LEVEL_SOLVER                  = "HYPRE_LEVEL_SOLVER";
const std::string CCPoissonSolverManager::PETSC_LEVEL_SOLVER                  = "PETSC_LEVEL_SOLVER";

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

namespace
{
Pointer<PoissonSolver>
allocate_default_krylov_solver(
    const std::string& solver_object_name,
    Pointer<Database> solver_input_db)
{
    Pointer<KrylovLinearSolver> krylov_solver =
        KrylovLinearSolverManager::getManager()->allocateSolver(
            KrylovLinearSolverManager::DEFAULT, solver_object_name, solver_input_db);
    Pointer<PoissonKrylovLinearSolverWrapper> krylov_wrapper = new PoissonKrylovLinearSolverWrapper(krylov_solver);
    krylov_wrapper->setOperator(new CCLaplaceOperator(solver_object_name+"::laplace_operator"));
    return krylov_wrapper;
}// allocate_default_krylov_solver

Pointer<PoissonSolver>
allocate_petsc_krylov_solver(
    const std::string& solver_object_name,
    Pointer<Database> solver_input_db)
{
    Pointer<KrylovLinearSolver> krylov_solver =
        KrylovLinearSolverManager::getManager()->allocateSolver(
            KrylovLinearSolverManager::PETSC, solver_object_name, solver_input_db);
    Pointer<PoissonKrylovLinearSolverWrapper> krylov_wrapper = new PoissonKrylovLinearSolverWrapper(krylov_solver);
    krylov_wrapper->setOperator(new CCLaplaceOperator(solver_object_name+"::laplace_operator"));
    return krylov_wrapper;
}// allocate_petsc_krylov_solver
}

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

Pointer<PoissonSolver>
CCPoissonSolverManager::allocateSolver(
    const std::string& solver_type,
    const std::string& solver_object_name,
    Pointer<Database> solver_input_db,
    const std::string& precond_type,
    const std::string& precond_object_name,
    Pointer<Database> precond_input_db) const
{
    Pointer<PoissonSolver> solver = allocateSolver(solver_type, solver_object_name, solver_input_db);
    Pointer<KrylovLinearSolver> p_solver = solver;
    if (!p_solver.isNull())
    {
        p_solver->setPreconditioner(allocateSolver(precond_type, precond_object_name, precond_input_db));
    }
    return solver;
}// allocateSolver

void
CCPoissonSolverManager::registerSolverFactoryFunction(
    const std::string& solver_type,
    SolverMaker solver_maker)
{
    if (d_solver_maker_map.find(solver_type) != d_solver_maker_map.end())
    {
        pout << "CCPoissonSolverManager::registerSolverFactoryFunction():\n"
             << "  NOTICE: overriding initialization function for solver_type = " << solver_type << "\n";
    }
    d_solver_maker_map[solver_type] = solver_maker;
    return;
}// registerSolverFactoryFunction

/////////////////////////////// PROTECTED ////////////////////////////////////

CCPoissonSolverManager::CCPoissonSolverManager()
    : d_solver_maker_map()
{
    registerSolverFactoryFunction(DEFAULT_KRYLOV_SOLVER              , allocate_default_krylov_solver);
    registerSolverFactoryFunction(PETSC_KRYLOV_SOLVER                , allocate_petsc_krylov_solver);
    registerSolverFactoryFunction(DEFAULT_FAC_PRECONDITIONER         , CCPoissonPointRelaxationFACOperator::allocate_solver);
    registerSolverFactoryFunction(POINT_RELAXATION_FAC_PRECONDITIONER, CCPoissonPointRelaxationFACOperator::allocate_solver);
    registerSolverFactoryFunction(DEFAULT_LEVEL_SOLVER               , CCPoissonHypreLevelSolver::allocate_solver);
    registerSolverFactoryFunction(HYPRE_LEVEL_SOLVER                 , CCPoissonHypreLevelSolver::allocate_solver);
    registerSolverFactoryFunction(PETSC_LEVEL_SOLVER                 , CCPoissonPETScLevelSolver::allocate_solver);
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
