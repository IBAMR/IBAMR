// Filename: CCPoissonPETScLevelSolver.C
// Created on 30 Aug 2010 by Boyce Griffith
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

#include "CCPoissonPETScLevelSolver.h"

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
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScVecUtilities.h>
#include <ibtk/PoissonUtilities.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_solve_system;
static Timer* t_initialize_solver_state;
static Timer* t_deallocate_solver_state;

// Number of ghosts cells used for each variable quantity.
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CCPoissonPETScLevelSolver::CCPoissonPETScLevelSolver(
    const std::string& object_name,
    Pointer<Database> input_db)
    : d_object_name(object_name),
      d_is_initialized(false),
      d_hierarchy(),
      d_level_num(-1),
      d_poisson_spec(d_object_name+"::Poisson specs"),
      d_default_bc_coef(new LocationIndexRobinBcCoefs<NDIM>(
                            d_object_name+"::default_bc_coef", Pointer<Database>(NULL))),
      d_bc_coefs(),
      d_homogeneous_bc(true),
      d_apply_time(0.0),
      d_options_prefix(""),
      d_petsc_ksp(PETSC_NULL),
      d_petsc_mat(PETSC_NULL),
      d_petsc_x(PETSC_NULL),
      d_petsc_b(PETSC_NULL),
      d_max_iterations(10),
      d_abs_residual_tol(0.0),
      d_rel_residual_tol(1.0e-6),
      d_initial_guess_nonzero(false),
      d_current_its(-1),
      d_current_residual_norm(-1.0),
      d_context(NULL),
      d_dof_index_idx(-1),
      d_dof_index_var(NULL),
      d_ghost_fill_sched(NULL),
      d_enable_logging(false)
{
    // Get values from the input database.
    if (!input_db.isNull())
    {
        d_options_prefix = input_db->getStringWithDefault("options_prefix", d_options_prefix);
        d_max_iterations = input_db->getIntegerWithDefault("max_iterations", d_max_iterations);
        d_abs_residual_tol = input_db->getDoubleWithDefault("absolute_residual_tol", d_abs_residual_tol);
        d_rel_residual_tol = input_db->getDoubleWithDefault("relative_residual_tol", d_rel_residual_tol);
        d_initial_guess_nonzero = input_db->getBoolWithDefault("initial_guess_nonzero", d_initial_guess_nonzero);
        d_enable_logging = input_db->getBoolWithDefault("enable_logging", d_enable_logging);
    }

    // Initialize the Poisson specifications.
    d_poisson_spec.setCZero();
    d_poisson_spec.setDConstant(-1.0);

    // Setup a default boundary condition object that specifies homogeneous
    // Dirichlet boundary conditions.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_default_bc_coef->setBoundaryValue(2*d  ,0.0);
        d_default_bc_coef->setBoundaryValue(2*d+1,0.0);
    }

    // Initialize the boundary conditions objects.
    setHomogeneousBc(d_homogeneous_bc);
    setPhysicalBcCoef(d_default_bc_coef);

    // Construct the DOF index variable/context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(d_object_name + "::CONTEXT");
    d_dof_index_var = new CellVariable<NDIM,int>(d_object_name + "::dof_index");
    if (var_db->checkVariableExists(d_dof_index_var->getName()))
    {
        d_dof_index_var = var_db->getVariable(d_dof_index_var->getName());
        d_dof_index_idx = var_db->mapVariableAndContextToIndex(d_dof_index_var, d_context);
        var_db->removePatchDataIndex(d_dof_index_idx);
    }
    d_dof_index_idx = var_db->registerVariableAndContext(d_dof_index_var, d_context, CELLG);

    // Setup Timers.
    IBTK_DO_ONCE(
        t_solve_system            = TimerManager::getManager()->getTimer("IBTK::CCPoissonPETScLevelSolver::solveSystem()");
        t_initialize_solver_state = TimerManager::getManager()->getTimer("IBTK::CCPoissonPETScLevelSolver::initializeSolverState()");
        t_deallocate_solver_state = TimerManager::getManager()->getTimer("IBTK::CCPoissonPETScLevelSolver::deallocateSolverState()");
                 );
    return;
}// CCPoissonPETScLevelSolver

CCPoissonPETScLevelSolver::~CCPoissonPETScLevelSolver()
{
    if (d_is_initialized) deallocateSolverState();
    delete d_default_bc_coef;
    return;
}// ~CCPoissonPETScLevelSolver

void
CCPoissonPETScLevelSolver::setPoissonSpecifications(
    const SAMRAI::solv::PoissonSpecifications& poisson_spec)
{
    d_poisson_spec = poisson_spec;
    return;
}// setPoissonSpecifications

void
CCPoissonPETScLevelSolver::setPhysicalBcCoef(
    RobinBcCoefStrategy<NDIM>* const bc_coef)
{
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(1,bc_coef));
    return;
}// setPhysicalBcCoef

void
CCPoissonPETScLevelSolver::setPhysicalBcCoefs(
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& bc_coefs)
{
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(&bc_coefs[0],&bc_coefs[0]+NDIM));
    return;
}// setPhysicalBcCoefs

void
CCPoissonPETScLevelSolver::setPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
    d_bc_coefs.resize(bc_coefs.size());
    for (unsigned int l = 0; l < bc_coefs.size(); ++l)
    {
        if (bc_coefs[l] != NULL)
        {
            d_bc_coefs[l] = bc_coefs[l];
        }
        else
        {
            d_bc_coefs[l] = d_default_bc_coef;
        }
    }
    return;
}// setPhysicalBcCoefs

void
CCPoissonPETScLevelSolver::setHomogeneousBc(
    bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

void
CCPoissonPETScLevelSolver::setTime(
    const double time)
{
    d_apply_time = time;
    return;
}// setTime

bool
CCPoissonPETScLevelSolver::solveSystem(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& b)
{
    IBTK_TIMER_START(t_solve_system);

    int ierr;

    if (d_enable_logging) plog << d_object_name << "::solveSystem():" << std::endl;

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x,b);

    // Configure solver.
    ierr = KSPSetTolerances(d_petsc_ksp, d_rel_residual_tol, d_abs_residual_tol, PETSC_DEFAULT, d_max_iterations); IBTK_CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(d_petsc_ksp, d_initial_guess_nonzero ? PETSC_TRUE : PETSC_FALSE); IBTK_CHKERRQ(ierr);

    // Solve the system.
    Pointer<PatchLevel<NDIM> > patch_level = d_hierarchy->getPatchLevel(d_level_num);
    const int x_idx = x.getComponentDescriptorIndex(0);
    Pointer<CellVariable<NDIM,double> > x_var = x.getComponentVariable(0);
    const int b_idx = b.getComponentDescriptorIndex(0);
    Pointer<CellVariable<NDIM,double> > b_var = b.getComponentVariable(0);
    if (!d_initial_guess_nonzero) PETScVecUtilities::copyToPatchLevelVec(d_petsc_x, x_idx, d_dof_index_idx, patch_level);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int b_adj_idx = var_db->registerClonedPatchDataIndex(b_var, b_idx);
    patch_level->allocatePatchData(b_adj_idx);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<CellData<NDIM,double> > b_data = patch->getPatchData(b_idx);
        Pointer<CellData<NDIM,double> > b_adj_data = patch->getPatchData(b_adj_idx);
        b_adj_data->copy(*b_data);
        if (!patch->getPatchGeometry()->intersectsPhysicalBoundary()) continue;
        PoissonUtilities::adjustCCBoundaryRhsEntries(patch, *b_adj_data, d_poisson_spec, d_bc_coefs, d_apply_time, d_homogeneous_bc);
    }
    PETScVecUtilities::copyToPatchLevelVec(d_petsc_b, b_adj_idx, d_dof_index_idx, patch_level);
    patch_level->deallocatePatchData(b_adj_idx);
    var_db->removePatchDataIndex(b_adj_idx);
    ierr = KSPSolve(d_petsc_ksp, d_petsc_b, d_petsc_x); IBTK_CHKERRQ(ierr);
    PETScVecUtilities::copyFromPatchLevelVec(d_petsc_x, x_idx, d_dof_index_idx, patch_level, d_ghost_fill_sched);

    // Log solver info.
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(d_petsc_ksp, &reason); IBTK_CHKERRQ(ierr);
    const bool converged = reason > 0;
    if (d_enable_logging)
    {
        plog << d_object_name << "::solveSystem(): solver " << (converged ? "converged" : "diverged") << "\n"
             << "iterations = " << d_current_its << "\n"
             << "residual norm = " << d_current_residual_norm << std::endl;
    }

    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();

    IBTK_TIMER_STOP(t_solve_system);
    return converged;
}// solveSystem

void
CCPoissonPETScLevelSolver::initializeSolverState(
    const SAMRAIVectorReal<NDIM,double>& x,
    const SAMRAIVectorReal<NDIM,double>& b)
{
    IBTK_TIMER_START(t_initialize_solver_state);

    // Rudimentary error checking.
#ifdef DEBUG_CHECK_ASSERTIONS
    if (x.getNumberOfComponents() != b.getNumberOfComponents())
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                   << "  vectors must have the same number of components" << std::endl);
    }

    const Pointer<PatchHierarchy<NDIM> >& patch_hierarchy = x.getPatchHierarchy();
    if (patch_hierarchy != b.getPatchHierarchy())
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                   << "  vectors must have the same hierarchy" << std::endl);
    }

    const int coarsest_ln = x.getCoarsestLevelNumber();
    if (coarsest_ln < 0)
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                   << "  coarsest level number must not be negative" << std::endl);
    }
    if (coarsest_ln != b.getCoarsestLevelNumber())
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                   << "  vectors must have same coarsest level number" << std::endl);
    }

    const int finest_ln = x.getFinestLevelNumber();
    if (finest_ln < coarsest_ln)
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                   << "  finest level number must be >= coarsest level number" << std::endl);
    }
    if (finest_ln != b.getFinestLevelNumber())
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                   << "  vectors must have same finest level number" << std::endl);
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (patch_hierarchy->getPatchLevel(ln).isNull())
        {
            TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                       << "  hierarchy level " << ln << " does not exist" << std::endl);
        }
    }

    if (coarsest_ln != finest_ln)
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                   << "  coarsest_ln != finest_ln in CCPoissonPETScLevelSolver" << std::endl);
    }
#endif
    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized) deallocateSolverState();

    // Get the hierarchy information.
    d_hierarchy = x.getPatchHierarchy();
    d_level_num = x.getCoarsestLevelNumber();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_level_num == x.getFinestLevelNumber());
#endif

    // Allocate DOF index data.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int x_idx = x.getComponentDescriptorIndex(0);
    Pointer<CellDataFactory<NDIM,double> > x_fac = var_db->getPatchDescriptor()->getPatchDataFactory(x_idx);
    const int depth = x_fac->getDefaultDepth();
    Pointer<CellDataFactory<NDIM,int> > dof_index_fac = var_db->getPatchDescriptor()->getPatchDataFactory(d_dof_index_idx);
    dof_index_fac->setDefaultDepth(depth);
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_num);
    if (!level->checkAllocated(d_dof_index_idx)) level->allocatePatchData(d_dof_index_idx);

    // Setup PETSc objects.
    int ierr;
    PETScVecUtilities::constructPatchLevelDOFIndices(d_num_dofs_per_proc, d_dof_index_idx, level);
    const int mpi_rank = SAMRAI_MPI::getRank();
    ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[mpi_rank], PETSC_DETERMINE, &d_petsc_x); IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[mpi_rank], PETSC_DETERMINE, &d_petsc_b); IBTK_CHKERRQ(ierr);
    PETScMatUtilities::constructPatchLevelCCLaplaceOp(d_petsc_mat, d_poisson_spec, d_bc_coefs, d_apply_time, d_num_dofs_per_proc, d_dof_index_idx, level);
    ierr = KSPCreate(PETSC_COMM_WORLD, &d_petsc_ksp); IBTK_CHKERRQ(ierr);
    ierr = KSPSetOperators(d_petsc_ksp, d_petsc_mat, d_petsc_mat, SAME_PRECONDITIONER); IBTK_CHKERRQ(ierr);
    if (!d_options_prefix.empty())
    {
        ierr = KSPSetOptionsPrefix(d_petsc_ksp, d_options_prefix.c_str()); IBTK_CHKERRQ(ierr);
    }
    ierr = KSPSetFromOptions(d_petsc_ksp); IBTK_CHKERRQ(ierr);
    d_ghost_fill_sched = PETScVecUtilities::constructGhostFillSchedule(x_idx, level);

    // Indicate that the solver is initialized.
    d_is_initialized = true;

    IBTK_TIMER_STOP(t_initialize_solver_state);
    return;
}// initializeSolverState

void
CCPoissonPETScLevelSolver::deallocateSolverState()
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_deallocate_solver_state);

    // Deallocate PETSc objects.
    int ierr;
    ierr = KSPDestroy(&d_petsc_ksp); IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&d_petsc_mat); IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&d_petsc_x); IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&d_petsc_b); IBTK_CHKERRQ(ierr);
    d_ghost_fill_sched.setNull();

    d_petsc_ksp = PETSC_NULL;
    d_petsc_mat = PETSC_NULL;
    d_petsc_x = PETSC_NULL;
    d_petsc_b = PETSC_NULL;

    // Deallocate DOF index data.
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_num);
    if (level->checkAllocated(d_dof_index_idx)) level->deallocatePatchData(d_dof_index_idx);

    // Indicate that the solver is NOT initialized.
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_solver_state);
    return;
}// deallocateSolverState

void
CCPoissonPETScLevelSolver::enableLogging(
    bool enabled)
{
    d_enable_logging = enabled;
    return;
}// enableLogging

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
