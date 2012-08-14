// Filename: SCPoissonSolverManager.C
// Created on 14 Aug 2012 by Boyce Griffith
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

#include "SCPoissonSolverManager.h"

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
#include <ibtk/KrylovLinearSolverManager.h>
#include <ibtk/PoissonFACPreconditioner.h>
#include <ibtk/PoissonKrylovLinearSolverWrapper.h>
#include <ibtk/SCLaplaceOperator.h>
#include <ibtk/SCPoissonHypreLevelSolver.h>
#include <ibtk/SCPoissonPETScLevelSolver.h>
#include <ibtk/SCPoissonPointRelaxationFACOperator.h>
#include <ibtk/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

SCPoissonSolverManager* SCPoissonSolverManager::s_solver_manager_instance = NULL;
bool SCPoissonSolverManager::s_registered_callback = false;
unsigned char SCPoissonSolverManager::s_shutdown_priority = 200;

SCPoissonSolverManager*
SCPoissonSolverManager::getManager()
{
    if (s_solver_manager_instance == NULL)
    {
        s_solver_manager_instance = new SCPoissonSolverManager();
    }
    if (!s_registered_callback)
    {
        ShutdownRegistry::registerShutdownRoutine(freeManager, s_shutdown_priority);
        s_registered_callback = true;
    }
    return s_solver_manager_instance;
}// getManager

void
SCPoissonSolverManager::freeManager()
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
    PoissonSpecifications poisson_spec(solver_object_name+"::PoissonSpecifications");
    poisson_spec.setCConstant(0.0);
    poisson_spec.setDConstant(-1.0);
    blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM> bc_coefs;  for (unsigned int d = 0; d < NDIM; ++d) bc_coefs[d] = NULL;
    Pointer<KrylovLinearSolver> krylov_solver = KrylovLinearSolverManager::getManager()->allocateSolver("DEFAULT", solver_object_name, solver_input_db);
    krylov_solver->setOperator(new SCLaplaceOperator(solver_object_name+"::LaplaceOperator", poisson_spec, bc_coefs));
    return new PoissonKrylovLinearSolverWrapper(krylov_solver);
}// allocate_default_krylov_solver

Pointer<PoissonSolver>
allocate_petsc_krylov_solver(
    const std::string& solver_object_name,
    Pointer<Database> solver_input_db)
{
    PoissonSpecifications poisson_spec(solver_object_name+"::PoissonSpecifications");
    poisson_spec.setCConstant(0.0);
    poisson_spec.setDConstant(-1.0);
    blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM> bc_coefs;  for (unsigned int d = 0; d < NDIM; ++d) bc_coefs[d] = NULL;
    Pointer<KrylovLinearSolver> krylov_solver = KrylovLinearSolverManager::getManager()->allocateSolver("PETSC_KRYLOV_LINEAR_SOLVER", solver_object_name, solver_input_db);
    krylov_solver->setOperator(new SCLaplaceOperator(solver_object_name+"::LaplaceOperator", poisson_spec, bc_coefs));
    return new PoissonKrylovLinearSolverWrapper(krylov_solver);
}// allocate_petsc_krylov_solver

Pointer<PoissonSolver>
allocate_point_relaxation_fac_preconditioner(
    const std::string& solver_object_name,
    Pointer<Database> solver_input_db)
{
    PoissonSpecifications poisson_spec(solver_object_name+"::PoissonSpecifications");
    poisson_spec.setCConstant(0.0);
    poisson_spec.setDConstant(-1.0);
    Pointer<PoissonFACPreconditionerStrategy> fac_operator = new SCPoissonPointRelaxationFACOperator(solver_object_name+"::FACOperator", solver_input_db);
    fac_operator->setPoissonSpecifications(poisson_spec);
    return new PoissonFACPreconditioner(solver_object_name, fac_operator);
}// allocate_point_relaxation_fac_preconditioner
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

Pointer<PoissonSolver>
SCPoissonSolverManager::allocateSolver(
    const std::string& solver_type,
    const std::string& solver_object_name,
    Pointer<Database> solver_input_db) const
{
    std::map<std::string,SolverMaker>::const_iterator it = d_solver_maker_map.find(solver_type);
    if (it == d_solver_maker_map.end())
    {
        TBOX_ERROR("SCPoissonSolverManager::allocateSolver():\n"
                   << "  unrecognized solver type: " << solver_type << "\n");
    }
    return (it->second)(solver_object_name, solver_input_db);
}// allocateSolver

void
SCPoissonSolverManager::registerSolverFactoryFunction(
    const std::string& solver_type,
    SolverMaker solver_maker)
{
    d_solver_maker_map[solver_type] = solver_maker;
    return;
}// registerSolverFactoryFunction

/////////////////////////////// PROTECTED ////////////////////////////////////

SCPoissonSolverManager::SCPoissonSolverManager()
    : d_solver_maker_map()
{
    d_solver_maker_map["DEFAULT_KRYLOV_LINEAR_SOLVER"       ] = allocate_default_krylov_solver;
    d_solver_maker_map["PETSC_KRYLOV_LINEAR_SOLVER"         ] = allocate_petsc_krylov_solver;
    d_solver_maker_map["DEFAULT_FAC_PRECONDITIONER"         ] = allocate_point_relaxation_fac_preconditioner;
    d_solver_maker_map["POINT_RELAXATION_FAC_PRECONDITIONER"] = allocate_point_relaxation_fac_preconditioner;
    d_solver_maker_map["DEFAULT_LEVEL_SOLVER"               ] = SCPoissonHypreLevelSolver::allocate_solver;
    d_solver_maker_map["HYPRE_LEVEL_SOLVER"                 ] = SCPoissonHypreLevelSolver::allocate_solver;
    d_solver_maker_map["PETSC_LEVEL_SOLVER"                 ] = SCPoissonPETScLevelSolver::allocate_solver;
    return;
}// SCPoissonSolverManager

SCPoissonSolverManager::~SCPoissonSolverManager()
{
    // intentionally blank
    return;
}// ~SCPoissonSolverManager

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
