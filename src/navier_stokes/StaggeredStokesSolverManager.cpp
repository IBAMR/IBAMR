// Filename: StaggeredStokesSolverManager.cpp
// Created on 16 Aug 2012 by Boyce Griffith
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

#include "ibamr/PETScKrylovStaggeredStokesSolver.h"
#include "ibamr/StaggeredStokesBlockFactorizationPreconditioner.h"
#include "ibamr/StaggeredStokesBoxRelaxationFACOperator.h"
#include "ibamr/StaggeredStokesFACPreconditioner.h"
#include "ibamr/StaggeredStokesFACPreconditionerStrategy.h"
#include "ibamr/StaggeredStokesLevelRelaxationFACOperator.h"
#include "ibamr/StaggeredStokesOperator.h"
#include "ibamr/StaggeredStokesPETScLevelSolver.h"
#include "ibamr/StaggeredStokesProjectionPreconditioner.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibamr/StaggeredStokesSolverManager.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/FACPreconditionerStrategy.h"
#include "ibtk/KrylovLinearSolver.h"
#include "ibtk/LinearOperator.h"
#include "ibtk/LinearSolver.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/Utilities.h"

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

StaggeredStokesSolverManager* StaggeredStokesSolverManager::s_solver_manager_instance = NULL;
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
    s_solver_manager_instance = NULL;
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
    krylov_solver->setOperator(new StaggeredStokesOperator(object_name + "::StokesOperator"));
    return krylov_solver;
} // allocate_petsc_krylov_solver

Pointer<StaggeredStokesSolver>
allocate_box_relaxation_fac_preconditioner(const std::string& object_name,
                                           Pointer<Database> input_db,
                                           const std::string& default_options_prefix)
{
    Pointer<StaggeredStokesFACPreconditionerStrategy> fac_operator =
        new StaggeredStokesBoxRelaxationFACOperator(object_name + "::FACOperator", input_db, default_options_prefix);
    return new StaggeredStokesFACPreconditioner(object_name, fac_operator, input_db, default_options_prefix);
} // allocate_box_relaxation_fac_preconditioner
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

Pointer<StaggeredStokesSolver>
StaggeredStokesSolverManager::allocateSolver(const std::string& solver_type,
                                             const std::string& solver_object_name,
                                             Pointer<Database> solver_input_db,
                                             const std::string& solver_default_options_prefix) const
{
    std::map<std::string, SolverMaker>::const_iterator it = d_solver_maker_map.find(solver_type);
    if (it == d_solver_maker_map.end())
    {
        TBOX_ERROR("StaggeredStokesSolverManager::allocateSolver():\n"
                   << "  unrecognized solver type: "
                   << solver_type
                   << "\n");
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
                                             const std::string& precond_default_options_prefix) const
{
    Pointer<StaggeredStokesSolver> solver =
        allocateSolver(solver_type, solver_object_name, solver_input_db, solver_default_options_prefix);
    Pointer<KrylovLinearSolver> p_solver = solver;
    if (p_solver)
        p_solver->setPreconditioner(
            allocateSolver(precond_type, precond_object_name, precond_input_db, precond_default_options_prefix));
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

StaggeredStokesSolverManager::StaggeredStokesSolverManager() : d_solver_maker_map()
{
    registerSolverFactoryFunction(DEFAULT_KRYLOV_SOLVER, allocate_petsc_krylov_solver);
    registerSolverFactoryFunction(PETSC_KRYLOV_SOLVER, allocate_petsc_krylov_solver);
    registerSolverFactoryFunction(DEFAULT_BLOCK_PRECONDITIONER,
                                  StaggeredStokesProjectionPreconditioner::allocate_solver);
    registerSolverFactoryFunction(BLOCK_FACTORIZATION_PRECONDITIONER,
                                  StaggeredStokesBlockFactorizationPreconditioner::allocate_solver);
    registerSolverFactoryFunction(PROJECTION_PRECONDITIONER, StaggeredStokesProjectionPreconditioner::allocate_solver);
    registerSolverFactoryFunction(DEFAULT_FAC_PRECONDITIONER, allocate_box_relaxation_fac_preconditioner);
    registerSolverFactoryFunction(BOX_RELAXATION_FAC_PRECONDITIONER, allocate_box_relaxation_fac_preconditioner);
    registerSolverFactoryFunction(LEVEL_RELAXATION_FAC_PRECONDITIONER,
                                  StaggeredStokesLevelRelaxationFACOperator::allocate_solver);
    registerSolverFactoryFunction(DEFAULT_LEVEL_SOLVER, StaggeredStokesPETScLevelSolver::allocate_solver);
    registerSolverFactoryFunction(PETSC_LEVEL_SOLVER, StaggeredStokesPETScLevelSolver::allocate_solver);
    return;
} // StaggeredStokesSolverManager

StaggeredStokesSolverManager::~StaggeredStokesSolverManager()
{
    // intentionally blank
    return;
} // ~StaggeredStokesSolverManager

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
