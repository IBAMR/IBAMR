// Filename: SCPoissonSolverManager.cpp
// Created on 14 Aug 2012 by Boyce Griffith
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

#include "ibtk/KrylovLinearSolver.h"
#include "ibtk/LinearOperator.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/PETScKrylovPoissonSolver.h"
#include "ibtk/PoissonSolver.h"
#include "ibtk/SCLaplaceOperator.h"
#include "ibtk/SCPoissonHypreLevelSolver.h"
#include "ibtk/SCPoissonPETScLevelSolver.h"
#include "ibtk/SCPoissonPointRelaxationFACOperator.h"
#include "ibtk/SCPoissonSolverManager.h"
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

const std::string SCPoissonSolverManager::UNDEFINED = "UNDEFINED";
const std::string SCPoissonSolverManager::DEFAULT_KRYLOV_SOLVER = "DEFAULT_KRYLOV_SOLVER";
const std::string SCPoissonSolverManager::PETSC_KRYLOV_SOLVER = "PETSC_KRYLOV_SOLVER";
const std::string SCPoissonSolverManager::DEFAULT_FAC_PRECONDITIONER = "DEFAULT_FAC_PRECONDITIONER";
const std::string SCPoissonSolverManager::POINT_RELAXATION_FAC_PRECONDITIONER = "POINT_RELAXATION_FAC_PRECONDITIONER";
const std::string SCPoissonSolverManager::DEFAULT_LEVEL_SOLVER = "DEFAULT_LEVEL_SOLVER";
const std::string SCPoissonSolverManager::HYPRE_LEVEL_SOLVER = "HYPRE_LEVEL_SOLVER";
const std::string SCPoissonSolverManager::PETSC_LEVEL_SOLVER = "PETSC_LEVEL_SOLVER";

SCPoissonSolverManager* SCPoissonSolverManager::s_solver_manager_instance = NULL;
bool SCPoissonSolverManager::s_registered_callback = false;
unsigned char SCPoissonSolverManager::s_shutdown_priority = 200;

SCPoissonSolverManager*
SCPoissonSolverManager::getManager()
{
    if (!s_solver_manager_instance)
    {
        s_solver_manager_instance = new SCPoissonSolverManager();
    }
    if (!s_registered_callback)
    {
        ShutdownRegistry::registerShutdownRoutine(freeManager, s_shutdown_priority);
        s_registered_callback = true;
    }
    return s_solver_manager_instance;
} // getManager

void
SCPoissonSolverManager::freeManager()
{
    delete s_solver_manager_instance;
    s_solver_manager_instance = NULL;
    return;
} // freeManager

namespace
{
Pointer<PoissonSolver>
allocate_petsc_krylov_solver(const std::string& object_name,
                             Pointer<Database> input_db,
                             const std::string& default_options_prefix)
{
    Pointer<PETScKrylovPoissonSolver> krylov_solver =
        new PETScKrylovPoissonSolver(object_name, input_db, default_options_prefix);
    krylov_solver->setOperator(new SCLaplaceOperator(object_name + "::laplace_operator"));
    return krylov_solver;
} // allocate_petsc_krylov_solver
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

Pointer<PoissonSolver>
SCPoissonSolverManager::allocateSolver(const std::string& solver_type,
                                       const std::string& solver_object_name,
                                       Pointer<Database> solver_input_db,
                                       const std::string& solver_default_options_prefix) const
{
    std::map<std::string, SolverMaker>::const_iterator it = d_solver_maker_map.find(solver_type);
    if (it == d_solver_maker_map.end())
    {
        TBOX_ERROR("SCPoissonSolverManager::allocateSolver():\n"
                   << "  unrecognized solver type: "
                   << solver_type
                   << "\n");
    }
    return (it->second)(solver_object_name, solver_input_db, solver_default_options_prefix);
} // allocateSolver

Pointer<PoissonSolver>
SCPoissonSolverManager::allocateSolver(const std::string& solver_type,
                                       const std::string& solver_object_name,
                                       Pointer<Database> solver_input_db,
                                       const std::string& solver_default_options_prefix,
                                       const std::string& precond_type,
                                       const std::string& precond_object_name,
                                       Pointer<Database> precond_input_db,
                                       const std::string& precond_default_options_prefix) const
{
    Pointer<PoissonSolver> solver =
        allocateSolver(solver_type, solver_object_name, solver_input_db, solver_default_options_prefix);
    Pointer<KrylovLinearSolver> p_solver = solver;
    if (p_solver)
    {
        p_solver->setPreconditioner(
            allocateSolver(precond_type, precond_object_name, precond_input_db, precond_default_options_prefix));
    }
    return solver;
} // allocateSolver

void
SCPoissonSolverManager::registerSolverFactoryFunction(const std::string& solver_type, SolverMaker solver_maker)
{
    if (d_solver_maker_map.find(solver_type) != d_solver_maker_map.end())
    {
        pout << "SCPoissonSolverManager::registerSolverFactoryFunction():\n"
             << "  NOTICE: overriding initialization function for solver_type = " << solver_type << "\n";
    }
    d_solver_maker_map[solver_type] = solver_maker;
    return;
} // registerSolverFactoryFunction

/////////////////////////////// PROTECTED ////////////////////////////////////

SCPoissonSolverManager::SCPoissonSolverManager() : d_solver_maker_map()
{
    registerSolverFactoryFunction(DEFAULT_KRYLOV_SOLVER, allocate_petsc_krylov_solver);
    registerSolverFactoryFunction(PETSC_KRYLOV_SOLVER, allocate_petsc_krylov_solver);
    registerSolverFactoryFunction(DEFAULT_FAC_PRECONDITIONER, SCPoissonPointRelaxationFACOperator::allocate_solver);
    registerSolverFactoryFunction(POINT_RELAXATION_FAC_PRECONDITIONER,
                                  SCPoissonPointRelaxationFACOperator::allocate_solver);
    registerSolverFactoryFunction(DEFAULT_LEVEL_SOLVER, SCPoissonHypreLevelSolver::allocate_solver);
    registerSolverFactoryFunction(HYPRE_LEVEL_SOLVER, SCPoissonHypreLevelSolver::allocate_solver);
    registerSolverFactoryFunction(PETSC_LEVEL_SOLVER, SCPoissonPETScLevelSolver::allocate_solver);
    return;
} // SCPoissonSolverManager

SCPoissonSolverManager::~SCPoissonSolverManager()
{
    // intentionally blank
    return;
} // ~SCPoissonSolverManager

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
