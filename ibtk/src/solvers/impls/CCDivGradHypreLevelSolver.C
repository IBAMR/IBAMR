// Filename: CCDivGradHypreLevelSolver.C
// Created on 20 Feb 2010 by Boyce Griffith
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

#include "CCDivGradHypreLevelSolver.h"

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
#include <ibtk/ibtk_utilities.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_solve_system;
static Timer* t_solve_system_hypre;
static Timer* t_initialize_solver_state;
static Timer* t_deallocate_solver_state;

// hypre solver options.
static const int RAP_TYPE_GALERKIN             = 0;
static const int RAP_TYPE_NON_GALERKIN_PARFLOW = 1;
static const int RAP_TYPE_GALERKIN_GENERAL     = 2;

static const int RELAX_TYPE_JACOBI                       = 0;
static const int RELAX_TYPE_WEIGHTED_JACOBI              = 1;
static const int RELAX_TYPE_RB_GAUSS_SEIDEL              = 2;
static const int RELAX_TYPE_RB_GAUSS_SEIDEL_NONSYMMETRIC = 3;

struct IndexComp
    : std::binary_function<Index<NDIM>,Index<NDIM>,bool>
{
    bool operator()(
        const Index<NDIM>& lhs,
        const Index<NDIM>& rhs) const
        {
            return (lhs(0) < rhs(0)
#if (NDIM > 1)
                    || (lhs(0) == rhs(0) && lhs(1) < rhs(1))
#if (NDIM > 2)
                    || (lhs(0) == rhs(0) && lhs(1) == rhs(1) && lhs(2) < rhs(2))
#endif
#endif
                    );
        }// operator()
};
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CCDivGradHypreLevelSolver::CCDivGradHypreLevelSolver(
    const std::string& object_name,
    Pointer<Database> input_db)
    : d_object_name(object_name),
      d_is_initialized(false),
      d_hierarchy(),
      d_level_num(-1),
      d_grid(NULL),
      d_stencil(NULL),
      d_matrix(NULL),
      d_rhs_vec(NULL),
      d_sol_vec(NULL),
      d_solver(NULL),
      d_precond(NULL),
      d_solver_type("PFMG"),
      d_precond_type("none"),
      d_max_iterations(10),
      d_abs_residual_tol(0.0),
      d_rel_residual_tol(1.0e-6),
      d_initial_guess_nonzero(false),
      d_rel_change(0),
      d_num_pre_relax_steps(1),
      d_num_post_relax_steps(1),
      d_memory_use(0),
      d_rap_type(RAP_TYPE_GALERKIN),
      d_relax_type(RELAX_TYPE_WEIGHTED_JACOBI),
      d_skip_relax(1),
      d_two_norm(1),
      d_current_its(-1),
      d_current_residual_norm(-1.0),
      d_enable_logging(false)
{
    if (NDIM == 1 || NDIM > 3)
    {
        TBOX_ERROR(d_object_name << "::CCDivGradHypreLevelSolver()"
                   << "  hypre solvers are only provided for 2D and 3D problems" << std::endl);
    }

    // Get values from the input database.
    if (!input_db.isNull())
    {
        d_enable_logging = input_db->getBoolWithDefault("enable_logging", d_enable_logging);
        d_solver_type = input_db->getStringWithDefault("solver_type", d_solver_type);
        d_precond_type = input_db->getStringWithDefault("precond_type", d_precond_type);
        d_max_iterations = input_db->getIntegerWithDefault("max_iterations", d_max_iterations);
        d_abs_residual_tol = input_db->getDoubleWithDefault("absolute_residual_tol", d_abs_residual_tol);
        d_rel_residual_tol = input_db->getDoubleWithDefault("relative_residual_tol", d_rel_residual_tol);
        d_initial_guess_nonzero = input_db->getBoolWithDefault("initial_guess_nonzero", d_initial_guess_nonzero);
        d_rel_change = input_db->getIntegerWithDefault("rel_change", d_rel_change);

        if (d_solver_type == "SMG" || d_precond_type == "SMG" || d_solver_type == "PFMG" || d_precond_type == "PFMG")
        {
            d_num_pre_relax_steps = input_db->getIntegerWithDefault("num_pre_relax_steps", d_num_pre_relax_steps);
            d_num_post_relax_steps = input_db->getIntegerWithDefault("num_post_relax_steps", d_num_post_relax_steps);
        }

        if (d_solver_type == "SMG" || d_precond_type == "SMG")
        {
            d_memory_use = input_db->getIntegerWithDefault("memory_use", d_memory_use);
        }

        if (d_solver_type == "PFMG" || d_precond_type == "PFMG")
        {
            d_rap_type = input_db->getIntegerWithDefault("rap_type", d_rap_type);
            d_relax_type = input_db->getIntegerWithDefault("relax_type", d_relax_type);
            d_skip_relax = input_db->getIntegerWithDefault("skip_relax", d_skip_relax);
        }

        if (d_solver_type == "PCG")
        {
            d_two_norm = input_db->getIntegerWithDefault("two_norm", d_two_norm);
        }
    }

    // Setup Timers.
    IBTK_DO_ONCE(
        t_solve_system            = TimerManager::getManager()->getTimer("IBTK::CCDivGradHypreLevelSolver::solveSystem()");
        t_solve_system_hypre      = TimerManager::getManager()->getTimer("IBTK::CCDivGradHypreLevelSolver::solveSystem()[hypre]");
        t_initialize_solver_state = TimerManager::getManager()->getTimer("IBTK::CCDivGradHypreLevelSolver::initializeSolverState()");
        t_deallocate_solver_state = TimerManager::getManager()->getTimer("IBTK::CCDivGradHypreLevelSolver::deallocateSolverState()");
                 );
    return;
}// CCDivGradHypreLevelSolver

CCDivGradHypreLevelSolver::~CCDivGradHypreLevelSolver()
{
    if (d_is_initialized) deallocateSolverState();
    return;
}// ~CCDivGradHypreLevelSolver

bool
CCDivGradHypreLevelSolver::solveSystem(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& b)
{
    SAMRAI_MPI::barrier();
    t_solve_system->start();

    if (d_enable_logging) plog << d_object_name << "::solveSystem():" << std::endl;

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x,b);

    // Solve the system using the hypre solver.
    static const int comp = 0;
    const int x_idx = x.getComponentDescriptorIndex(comp);
    const int b_idx = b.getComponentDescriptorIndex(comp);

    bool converged = true;
    IntVector<NDIM> chkbrd_mode_id;
#if (NDIM > 2)
    for (chkbrd_mode_id(2) = 0; chkbrd_mode_id(2) < 2; ++chkbrd_mode_id(2))
    {
#endif
        for (chkbrd_mode_id(1) = 0; chkbrd_mode_id(1) < 2; ++chkbrd_mode_id(1))
        {
            for (chkbrd_mode_id(0) = 0; chkbrd_mode_id(0) < 2; ++chkbrd_mode_id(0))
            {
                bool converged_mode = solveSystem(x_idx, b_idx, chkbrd_mode_id);
                if (d_enable_logging)
                {
                    plog << d_object_name << "::solveSystem(): solver " << (converged_mode ? "converged" : "diverged") << "\n"
                         << "chkbrd_mode_id = " << chkbrd_mode_id << "\n"
                         << "iterations = " << d_current_its << "\n"
                         << "residual norm = " << d_current_residual_norm << std::endl;
                }
                converged = converged && converged_mode;
            }
        }
#if (NDIM > 2)
    }
#endif

    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();

    t_solve_system->stop();
    return converged;
}// solveSystem

void
CCDivGradHypreLevelSolver::initializeSolverState(
    const SAMRAIVectorReal<NDIM,double>& x,
    const SAMRAIVectorReal<NDIM,double>& b)
{
    SAMRAI_MPI::barrier();
    t_initialize_solver_state->start();

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
                   << "  coarsest_ln != finest_ln in CCDivGradHypreLevelSolver" << std::endl);
    }
#endif
    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized) deallocateSolverState();

    // Get the hierarchy information.
    d_hierarchy = x.getPatchHierarchy();
    d_level_num = x.getCoarsestLevelNumber();

    // Allocate and initialize the hypre data structures.
    allocateHypreData();
    setMatrixCoefficients();
    setupHypreSolver();

    // Indicate that the solver is initialized.
    d_is_initialized = true;

    t_initialize_solver_state->stop();
    return;
}// initializeSolverState

void
CCDivGradHypreLevelSolver::deallocateSolverState()
{
    if (!d_is_initialized) return;

    SAMRAI_MPI::barrier();
    t_deallocate_solver_state->start();

    // Deallocate the hypre data structures.
    destroyHypreSolver();
    deallocateHypreData();

    // Indicate that the solver is NOT initialized.
    d_is_initialized = false;

    t_deallocate_solver_state->stop();
    return;
}// deallocateSolverState

void
CCDivGradHypreLevelSolver::enableLogging(
    bool enabled)
{
    d_enable_logging = enabled;
    return;
}// enableLogging

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CCDivGradHypreLevelSolver::allocateHypreData()
{
    // Get the MPI communicator.
#ifdef HAVE_MPI
    MPI_Comm communicator = SAMRAI_MPI::getCommunicator();
#else
    MPI_Comm communicator;
#endif

    // Setup the hypre grid.
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_num);
    Pointer<CartesianGridGeometry<NDIM> > grid_geometry = d_hierarchy->getGridGeometry();
    const IntVector<NDIM>& ratio = level->getRatio();
    const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift(ratio);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(periodic_shift.min() > 0);
#endif

    HYPRE_StructGridCreate(communicator, NDIM, &d_grid);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        const Box<NDIM> patch_box = Box<NDIM>::coarsen(level->getPatch(p())->getBox(),2);
        Index<NDIM> lower = patch_box.lower();
        Index<NDIM> upper = patch_box.upper();
        HYPRE_StructGridSetExtents(d_grid, lower, upper);
    }

    int hypre_periodic_shift[3];
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        hypre_periodic_shift[d] = periodic_shift(d)/2;
    }
    for (int d = NDIM; d < 3; ++d)
    {
        hypre_periodic_shift[d] = 0;
    }
    HYPRE_StructGridSetPeriodic(d_grid, hypre_periodic_shift);
    HYPRE_StructGridAssemble(d_grid);

    // Allocate stencil data and set stencil offsets.
    static const int stencil_sz = 2*NDIM+1;
#if (NDIM == 2)
    int stencil_offsets[stencil_sz][2] = {
        { -1, 0 }, { 0, -1}, { +1, 0}, { 0, +1 }, { 0, 0 }
    };
#endif
#if (NDIM == 3)
    int stencil_offsets[stencil_sz][3] = {
        { -1,  0,  0}, { 0,  -1,  0}, { 0,  0,  -1},
        { +1,  0,  0}, { 0,  +1,  0}, { 0,  0,  +1},
        { 0,  0,  0}
    };
#endif
    HYPRE_StructStencilCreate(NDIM, stencil_sz, &d_stencil);
    for (int s = 0; s < stencil_sz; ++s)
    {
        HYPRE_StructStencilSetElement(d_stencil, s, stencil_offsets[s]);
    }

    // Allocate the hypre matrix.
#if (NDIM == 2)
    int full_ghosts[2*3] = { 1, 1, 1, 1, 0, 0 };
#endif
#if (NDIM == 3)
    int full_ghosts[2*3] = { 1, 1, 1, 1, 1, 1 };
#endif
    int   no_ghosts[2*3] = { 0, 0, 0, 0, 0, 0 };

    HYPRE_StructMatrixCreate(communicator, d_grid, d_stencil, &d_matrix);
    HYPRE_StructMatrixSetNumGhost(d_matrix, full_ghosts);
    HYPRE_StructMatrixSetSymmetric(d_matrix, 0);
    HYPRE_StructMatrixInitialize(d_matrix);

    // Allocate the hypre vectors.
    HYPRE_StructVectorCreate(communicator, d_grid, &d_sol_vec);
    HYPRE_StructVectorSetNumGhost(d_sol_vec, full_ghosts);
    HYPRE_StructVectorInitialize(d_sol_vec);

    HYPRE_StructVectorCreate(communicator, d_grid, &d_rhs_vec);
    HYPRE_StructVectorSetNumGhost(d_rhs_vec, no_ghosts);
    HYPRE_StructVectorInitialize(d_rhs_vec);
    return;
}// allocateHypreData

void
CCDivGradHypreLevelSolver::setMatrixCoefficients()
{
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_num);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM> patch_box = Box<NDIM>::coarsen(patch->getBox(),2);
        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();

        const int stencil_sz = 2*NDIM+1;
        int stencil_indices[stencil_sz];
        for (int i = 0; i < stencil_sz; ++i)
        {
            stencil_indices[i] = i;
        }

        std::vector<double> mat_vals(stencil_sz,0.0);
        for (Box<NDIM>::Iterator b(patch_box); b; b++)
        {
            Index<NDIM> i = b();

            const SideIndex<NDIM> ixlower(
                i, SideIndex<NDIM>::X, SideIndex<NDIM>::Lower);
            mat_vals[0] = -1.0/(2.0*dx[0]*2.0*dx[0]);

            const SideIndex<NDIM> iylower(
                i, SideIndex<NDIM>::Y, SideIndex<NDIM>::Lower);
            mat_vals[1] = -1.0/(2.0*dx[1]*2.0*dx[1]);
#if (NDIM == 3)
            SideIndex<NDIM> izlower(
                i, SideIndex<NDIM>::Z, SideIndex<NDIM>::Lower);
            mat_vals[2] = -1.0/(2.0*dx[2]*2.0*dx[2]);
#endif
            const SideIndex<NDIM> ixupper(
                i, SideIndex<NDIM>::X, SideIndex<NDIM>::Upper);
            mat_vals[NDIM+0] = -1.0/(2.0*dx[0]*2.0*dx[0]);

            const SideIndex<NDIM> iyupper(
                i, SideIndex<NDIM>::Y, SideIndex<NDIM>::Upper);
            mat_vals[NDIM+1] = -1.0/(2.0*dx[1]*2.0*dx[1]);
#if (NDIM == 3)
            SideIndex<NDIM> izupper(
                i, SideIndex<NDIM>::Z, SideIndex<NDIM>::Upper);
            mat_vals[NDIM+2] = -1.0/(2.0*dx[2]*2.0*dx[2]);
#endif
            mat_vals[2*NDIM] = 0.0;
            for (int k = 0; k < 2*NDIM; ++k)
            {
                mat_vals[2*NDIM] -= mat_vals[k];
            }

            HYPRE_StructMatrixSetValues(d_matrix, i, stencil_sz, stencil_indices, &mat_vals[0]);
        }
    }

    // Assemble the hypre matrix.
    HYPRE_StructMatrixAssemble(d_matrix);
    return;
}// setMatrixCoefficients

void
CCDivGradHypreLevelSolver::setupHypreSolver()
{
    // Get the MPI communicator.
#ifdef HAVE_MPI
    MPI_Comm communicator = SAMRAI_MPI::getCommunicator();
#else
    MPI_Comm communicator;
#endif

    // When using a Krylov method, setup the preconditioner.
    if (d_solver_type == "PCG" || d_solver_type == "GMRES" || d_solver_type == "FlexGMRES" || d_solver_type == "LGMRES" || d_solver_type == "BiCGSTAB")
    {
        if (d_precond_type == "PFMG")
        {
            HYPRE_StructPFMGCreate(communicator, &d_precond);
            HYPRE_StructPFMGSetMaxIter(d_precond, 1);
            HYPRE_StructPFMGSetTol(d_precond, 0.0);
            HYPRE_StructPFMGSetZeroGuess(d_precond);
            HYPRE_StructPFMGSetRAPType(d_precond, d_rap_type);
            HYPRE_StructPFMGSetRelaxType(d_precond, d_relax_type);
            HYPRE_StructPFMGSetNumPreRelax(d_precond, d_num_pre_relax_steps);
            HYPRE_StructPFMGSetNumPostRelax(d_precond, d_num_post_relax_steps);
            HYPRE_StructPFMGSetSkipRelax(d_precond, d_skip_relax);
        }
        else if (d_precond_type == "SMG")
        {
            HYPRE_StructSMGCreate(communicator, &d_precond);
            HYPRE_StructSMGSetMaxIter(d_precond, 1);
            HYPRE_StructSMGSetTol(d_precond, 0.0);
            HYPRE_StructSMGSetZeroGuess(d_precond);
            HYPRE_StructSMGSetMemoryUse(d_precond, d_memory_use);
            HYPRE_StructSMGSetNumPreRelax(d_precond, d_num_pre_relax_steps);
            HYPRE_StructSMGSetNumPostRelax(d_precond, d_num_post_relax_steps);
        }
        else if (d_precond_type == "Jacobi")
        {
            HYPRE_StructJacobiCreate(communicator, &d_precond);
            HYPRE_StructJacobiSetMaxIter(d_precond, 2);
            HYPRE_StructJacobiSetTol(d_precond, 0.0);
            HYPRE_StructJacobiSetZeroGuess(d_precond);
        }
    }

    // Setup the solver.
    if (d_solver_type == "PFMG")
    {
        HYPRE_StructPFMGCreate(communicator, &d_solver);
        HYPRE_StructPFMGSetMaxIter(d_solver, d_max_iterations);
        HYPRE_StructPFMGSetTol(d_solver, d_rel_residual_tol);
        HYPRE_StructPFMGSetRelChange(d_solver, d_rel_change);
        HYPRE_StructPFMGSetRAPType(d_solver, d_rap_type);
        HYPRE_StructPFMGSetRelaxType(d_solver, d_relax_type);
        HYPRE_StructPFMGSetNumPreRelax(d_solver, d_num_pre_relax_steps);
        HYPRE_StructPFMGSetNumPostRelax(d_solver, d_num_post_relax_steps);
        HYPRE_StructPFMGSetSkipRelax(d_solver, d_skip_relax);
        if (d_initial_guess_nonzero)
        {
            HYPRE_StructPFMGSetNonZeroGuess(d_solver);
        }
        else
        {
            HYPRE_StructPFMGSetZeroGuess(d_solver);
        }
        HYPRE_StructPFMGSetup(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
    }
    else if (d_solver_type == "SMG")
    {
        HYPRE_StructSMGCreate(communicator, &d_solver);
        HYPRE_StructSMGSetMaxIter(d_solver, d_max_iterations);
        HYPRE_StructSMGSetTol(d_solver, d_rel_residual_tol);
        HYPRE_StructSMGSetRelChange(d_solver, d_rel_change);
        HYPRE_StructSMGSetMemoryUse(d_solver, d_memory_use);
        HYPRE_StructSMGSetNumPreRelax(d_solver, d_num_pre_relax_steps);
        HYPRE_StructSMGSetNumPostRelax(d_solver, d_num_post_relax_steps);
        if (d_initial_guess_nonzero)
        {
            HYPRE_StructSMGSetNonZeroGuess(d_solver);
        }
        else
        {
            HYPRE_StructSMGSetZeroGuess(d_solver);
        }
        HYPRE_StructSMGSetup(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
    }
    else if (d_solver_type == "PCG")
    {
        HYPRE_StructPCGCreate(communicator, &d_solver);
        HYPRE_StructPCGSetMaxIter(d_solver, d_max_iterations);
        HYPRE_StructPCGSetTol(d_solver, d_rel_residual_tol);
        HYPRE_StructPCGSetAbsoluteTol(d_solver, d_abs_residual_tol);
        HYPRE_StructPCGSetTwoNorm(d_solver, d_two_norm);
        HYPRE_StructPCGSetRelChange(d_solver, d_rel_change);
        if (d_precond_type == "PFMG")
        {
            HYPRE_StructPCGSetPrecond(d_solver,
                                      HYPRE_StructPFMGSolve,
                                      HYPRE_StructPFMGSetup,
                                      d_precond);
        }
        else if (d_precond_type == "SMG")
        {
            HYPRE_StructPCGSetPrecond(d_solver,
                                      HYPRE_StructSMGSolve,
                                      HYPRE_StructSMGSetup,
                                      d_precond);
        }
        else if (d_precond_type == "Jacobi")
        {
            HYPRE_StructPCGSetPrecond(d_solver,
                                      HYPRE_StructJacobiSolve,
                                      HYPRE_StructJacobiSetup,
                                      d_precond);
        }
        else if (d_precond_type == "diagonal_scaling")
        {
            d_precond = NULL;
            HYPRE_StructPCGSetPrecond(d_solver,
                                      HYPRE_StructDiagScale,
                                      HYPRE_StructDiagScaleSetup,
                                      d_precond);
        }
        else if (d_precond_type == "none")
        {
            d_precond = NULL;
        }
        else
        {
            TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                       << "  unknown preconditioner type: " << d_precond_type << std::endl);
        }
        HYPRE_StructPCGSetup(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
    }
    else if (d_solver_type == "GMRES")
    {
        HYPRE_StructGMRESCreate(communicator, &d_solver);
        HYPRE_StructGMRESSetMaxIter(d_solver, d_max_iterations);
        HYPRE_StructGMRESSetTol(d_solver, d_rel_residual_tol);
        HYPRE_StructGMRESSetAbsoluteTol(d_solver, d_abs_residual_tol);
        if (d_precond_type == "PFMG")
        {
            HYPRE_StructGMRESSetPrecond(d_solver,
                                        HYPRE_StructPFMGSolve,
                                        HYPRE_StructPFMGSetup,
                                        d_precond);
        }
        else if (d_precond_type == "SMG")
        {
            HYPRE_StructGMRESSetPrecond(d_solver,
                                        HYPRE_StructSMGSolve,
                                        HYPRE_StructSMGSetup,
                                        d_precond);
        }
        else if (d_precond_type == "Jacobi")
        {
            HYPRE_StructGMRESSetPrecond(d_solver,
                                        HYPRE_StructJacobiSolve,
                                        HYPRE_StructJacobiSetup,
                                        d_precond);
        }
        else if (d_precond_type == "diagonal_scaling")
        {
            d_precond = NULL;
            HYPRE_StructGMRESSetPrecond(d_solver,
                                        HYPRE_StructDiagScale,
                                        HYPRE_StructDiagScaleSetup,
                                        d_precond);
        }
        else if (d_precond_type == "none")
        {
            d_precond = NULL;
        }
        else
        {
            TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                       << "  unknown preconditioner type: " << d_precond_type << std::endl);
        }
        HYPRE_StructGMRESSetup(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
    }
    else if (d_solver_type == "FlexGMRES")
    {
        HYPRE_StructFlexGMRESCreate(communicator, &d_solver);
        HYPRE_StructFlexGMRESSetMaxIter(d_solver, d_max_iterations);
        HYPRE_StructFlexGMRESSetTol(d_solver, d_rel_residual_tol);
        HYPRE_StructFlexGMRESSetAbsoluteTol(d_solver, d_abs_residual_tol);
        if (d_precond_type == "PFMG")
        {
            HYPRE_StructFlexGMRESSetPrecond(d_solver,
                                            HYPRE_StructPFMGSolve,
                                            HYPRE_StructPFMGSetup,
                                            d_precond);
        }
        else if (d_precond_type == "SMG")
        {
            HYPRE_StructFlexGMRESSetPrecond(d_solver,
                                            HYPRE_StructSMGSolve,
                                            HYPRE_StructSMGSetup,
                                            d_precond);
        }
        else if (d_precond_type == "Jacobi")
        {
            HYPRE_StructFlexGMRESSetPrecond(d_solver,
                                            HYPRE_StructJacobiSolve,
                                            HYPRE_StructJacobiSetup,
                                            d_precond);
        }
        else if (d_precond_type == "diagonal_scaling")
        {
            d_precond = NULL;
            HYPRE_StructFlexGMRESSetPrecond(d_solver,
                                            HYPRE_StructDiagScale,
                                            HYPRE_StructDiagScaleSetup,
                                            d_precond);
        }
        else if (d_precond_type == "none")
        {
            d_precond = NULL;
        }
        else
        {
            TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                       << "  unknown preconditioner type: " << d_precond_type << std::endl);
        }
        HYPRE_StructFlexGMRESSetup(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
    }
    else if (d_solver_type == "LGMRES")
    {
        HYPRE_StructLGMRESCreate(communicator, &d_solver);
        HYPRE_StructLGMRESSetMaxIter(d_solver, d_max_iterations);
        HYPRE_StructLGMRESSetTol(d_solver, d_rel_residual_tol);
        HYPRE_StructLGMRESSetAbsoluteTol(d_solver, d_abs_residual_tol);
        if (d_precond_type == "PFMG")
        {
            HYPRE_StructLGMRESSetPrecond(d_solver,
                                         HYPRE_StructPFMGSolve,
                                         HYPRE_StructPFMGSetup,
                                         d_precond);
        }
        else if (d_precond_type == "SMG")
        {
            HYPRE_StructLGMRESSetPrecond(d_solver,
                                         HYPRE_StructSMGSolve,
                                         HYPRE_StructSMGSetup,
                                         d_precond);
        }
        else if (d_precond_type == "Jacobi")
        {
            HYPRE_StructLGMRESSetPrecond(d_solver,
                                         HYPRE_StructJacobiSolve,
                                         HYPRE_StructJacobiSetup,
                                         d_precond);
        }
        else if (d_precond_type == "diagonal_scaling")
        {
            d_precond = NULL;
            HYPRE_StructLGMRESSetPrecond(d_solver,
                                         HYPRE_StructDiagScale,
                                         HYPRE_StructDiagScaleSetup,
                                         d_precond);
        }
        else if (d_precond_type == "none")
        {
            d_precond = NULL;
        }
        else
        {
            TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                       << "  unknown preconditioner type: " << d_precond_type << std::endl);
        }
        HYPRE_StructLGMRESSetup(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
    }
    else if (d_solver_type == "BiCGSTAB")
    {
        HYPRE_StructBiCGSTABCreate(communicator, &d_solver);
        HYPRE_StructBiCGSTABSetMaxIter(d_solver, d_max_iterations);
        HYPRE_StructBiCGSTABSetTol(d_solver, d_rel_residual_tol);
        HYPRE_StructBiCGSTABSetAbsoluteTol(d_solver, d_abs_residual_tol);
        if (d_precond_type == "PFMG")
        {
            HYPRE_StructBiCGSTABSetPrecond(d_solver,
                                           HYPRE_StructPFMGSolve,
                                           HYPRE_StructPFMGSetup,
                                           d_precond);
        }
        else if (d_precond_type == "SMG")
        {
            HYPRE_StructBiCGSTABSetPrecond(d_solver,
                                           HYPRE_StructSMGSolve,
                                           HYPRE_StructSMGSetup,
                                           d_precond);
        }
        else if (d_precond_type == "Jacobi")
        {
            HYPRE_StructBiCGSTABSetPrecond(d_solver,
                                           HYPRE_StructJacobiSolve,
                                           HYPRE_StructJacobiSetup,
                                           d_precond);
        }
        else if (d_precond_type == "diagonal_scaling")
        {
            d_precond = NULL;
            HYPRE_StructBiCGSTABSetPrecond(d_solver,
                                           HYPRE_StructDiagScale,
                                           HYPRE_StructDiagScaleSetup,
                                           d_precond);
        }
        else if (d_precond_type == "none")
        {
            d_precond = NULL;
        }
        else
        {
            TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                       << "  unknown preconditioner type: " << d_precond_type << std::endl);
        }
        HYPRE_StructBiCGSTABSetup(d_solver,d_matrix, d_rhs_vec, d_sol_vec);
    }
    else
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                   << "  unknown solver type: " << d_solver_type << std::endl);
    }
    return;
}// setupHypreSolver

bool
CCDivGradHypreLevelSolver::solveSystem(
    const int x_idx,
    const int b_idx,
    const IntVector<NDIM>& chkbrd_mode_id)
{
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_num);

    // Copy solution and right-hand-side data to hypre structures.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM,double> > x_data = patch->getPatchData(x_idx);
        copyToHypre(d_sol_vec, x_data, patch_box, chkbrd_mode_id);
        Pointer<CellData<NDIM,double> > b_data = patch->getPatchData(b_idx);
        copyToHypre(d_rhs_vec, b_data, patch_box, chkbrd_mode_id);
    }

    // Assemble the hypre vectors.
    HYPRE_StructVectorAssemble(d_sol_vec);
    HYPRE_StructVectorAssemble(d_rhs_vec);

    // Solve the system.
    SAMRAI_MPI::barrier();
    t_solve_system_hypre->start();

    if (d_solver_type == "PFMG")
    {
        HYPRE_StructPFMGSetMaxIter(d_solver, d_max_iterations);
        HYPRE_StructPFMGSetTol(d_solver, d_rel_residual_tol);
        HYPRE_StructPFMGSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_StructPFMGGetNumIterations(d_solver, &d_current_its);
        HYPRE_StructPFMGGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "SMG")
    {
        HYPRE_StructSMGSetMaxIter(d_solver, d_max_iterations);
        HYPRE_StructSMGSetTol(d_solver, d_rel_residual_tol);
        HYPRE_StructSMGSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_StructSMGGetNumIterations(d_solver, &d_current_its);
        HYPRE_StructSMGGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "PCG")
    {
        HYPRE_StructPCGSetMaxIter(d_solver, d_max_iterations);
        HYPRE_StructPCGSetTol(d_solver, d_rel_residual_tol);
        HYPRE_StructPCGSetAbsoluteTol(d_solver, d_abs_residual_tol);
        HYPRE_StructPCGSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_StructPCGGetNumIterations(d_solver, &d_current_its);
        HYPRE_StructPCGGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "GMRES")
    {
        HYPRE_StructGMRESSetMaxIter(d_solver, d_max_iterations);
        HYPRE_StructGMRESSetTol(d_solver, d_rel_residual_tol);
        HYPRE_StructGMRESSetAbsoluteTol(d_solver, d_abs_residual_tol);
        HYPRE_StructGMRESSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_StructGMRESGetNumIterations(d_solver, &d_current_its);
        HYPRE_StructGMRESGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "FlexGMRES")
    {
        HYPRE_StructFlexGMRESSetMaxIter(d_solver, d_max_iterations);
        HYPRE_StructFlexGMRESSetTol(d_solver, d_rel_residual_tol);
        HYPRE_StructFlexGMRESSetAbsoluteTol(d_solver, d_abs_residual_tol);
        HYPRE_StructFlexGMRESSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_StructFlexGMRESGetNumIterations(d_solver, &d_current_its);
        HYPRE_StructFlexGMRESGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "LGMRES")
    {
        HYPRE_StructLGMRESSetMaxIter(d_solver, d_max_iterations);
        HYPRE_StructLGMRESSetTol(d_solver, d_rel_residual_tol);
        HYPRE_StructLGMRESSetAbsoluteTol(d_solver, d_abs_residual_tol);
        HYPRE_StructLGMRESSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_StructLGMRESGetNumIterations(d_solver, &d_current_its);
        HYPRE_StructLGMRESGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "BiCGSTAB")
    {
        HYPRE_StructBiCGSTABSetMaxIter(d_solver, d_max_iterations);
        HYPRE_StructBiCGSTABSetTol(d_solver, d_rel_residual_tol);
        HYPRE_StructBiCGSTABSetAbsoluteTol(d_solver, d_abs_residual_tol);
        HYPRE_StructBiCGSTABSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_StructBiCGSTABGetNumIterations(d_solver, &d_current_its);
        HYPRE_StructBiCGSTABGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }

    t_solve_system_hypre->stop();

    // Pull the solution vector out of the hypre structures.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM,double> > x_data = patch->getPatchData(x_idx);
        copyFromHypre(x_data, d_sol_vec, patch_box, chkbrd_mode_id);
    }
    return (d_current_residual_norm <= d_rel_residual_tol || d_current_residual_norm <= d_abs_residual_tol);
}// solveSystem

void
CCDivGradHypreLevelSolver::copyToHypre(
    HYPRE_StructVector vector,
    const Pointer<CellData<NDIM,double> >& src_data,
    const Box<NDIM>& box,
    const IntVector<NDIM>& chkbrd_mode_id)
{
    static const int ratio = 2;
    Box<NDIM> coarsened_box = Box<NDIM>::coarsen(box,ratio);
    CellData<NDIM,double> hypre_data(coarsened_box,1,0);
    for (Box<NDIM>::Iterator b(coarsened_box); b; b++)
    {
        const Index<NDIM>& coarse_i = b();
        Index<NDIM> fine_i = coarse_i*ratio;
        hypre_data(coarse_i) = (*src_data)(fine_i+chkbrd_mode_id);
    }
    HYPRE_StructVectorSetBoxValues(vector,coarsened_box.lower(),coarsened_box.upper(),hypre_data.getPointer());
    return;
}// copyToHypre

void
CCDivGradHypreLevelSolver::copyFromHypre(
    Pointer<CellData<NDIM,double> >& dst_data,
    HYPRE_StructVector vector,
    const Box<NDIM>& box,
    const IntVector<NDIM>& chkbrd_mode_id)
{
    static const int ratio = 2;
    Box<NDIM> coarsened_box = Box<NDIM>::coarsen(box,ratio);
    CellData<NDIM,double> hypre_data(coarsened_box,1,0);
    HYPRE_StructVectorGetBoxValues(vector,coarsened_box.lower(),coarsened_box.upper(),hypre_data.getPointer());
    for (Box<NDIM>::Iterator b(coarsened_box); b; b++)
    {
        const Index<NDIM>& coarse_i = b();
        const Index<NDIM> fine_i = coarse_i*ratio;
        (*dst_data)(fine_i+chkbrd_mode_id) = hypre_data(coarse_i);
    }
    return;
}// copyFromHypre

void
CCDivGradHypreLevelSolver::destroyHypreSolver()
{
    // Destroy the solver.
    if (d_solver_type == "PFMG")
    {
        HYPRE_StructPFMGDestroy(d_solver);
    }
    else if (d_solver_type == "SMG")
    {
        HYPRE_StructSMGDestroy(d_solver);
    }
    else if (d_solver_type == "PCG")
    {
        HYPRE_StructPCGDestroy(d_solver);
    }
    else if (d_solver_type == "GMRES")
    {
        HYPRE_StructGMRESDestroy(d_solver);
    }
    else if (d_solver_type == "FlexGMRES")
    {
        HYPRE_StructFlexGMRESDestroy(d_solver);
    }
    else if (d_solver_type == "LGMRES")
    {
        HYPRE_StructLGMRESDestroy(d_solver);
    }
    else if (d_solver_type == "BiCGSTAB")
    {
        HYPRE_StructBiCGSTABDestroy(d_solver);
    }

    // When using a Krylov method, destroy the preconditioner.
    if (d_solver_type == "PCG" || d_solver_type == "GMRES" || d_solver_type == "FlexGMRES" || d_solver_type == "LGMRES" || d_solver_type == "BiCGSTAB")
    {
        if (d_precond_type == "PFMG")
        {
            HYPRE_StructPFMGDestroy(d_precond);
        }
        else if (d_precond_type == "SMG")
        {
            HYPRE_StructSMGDestroy(d_precond);
        }
        else if (d_precond_type == "Jacobi")
        {
            HYPRE_StructJacobiDestroy(d_precond);
        }
    }

    // Set the solver and preconditioner pointers to NULL.
    d_solver  = NULL;
    d_precond = NULL;
    return;
}// destroyHypreSolver

void
CCDivGradHypreLevelSolver::deallocateHypreData()
{
    if (d_stencil) HYPRE_StructStencilDestroy(d_stencil);
    if (d_grid   ) HYPRE_StructGridDestroy(d_grid);
    if (d_matrix ) HYPRE_StructMatrixDestroy(d_matrix);
    if (d_sol_vec) HYPRE_StructVectorDestroy(d_sol_vec);
    if (d_rhs_vec) HYPRE_StructVectorDestroy(d_rhs_vec);
    d_grid    = NULL;
    d_stencil = NULL;
    d_matrix  = NULL;
    d_sol_vec = NULL;
    d_rhs_vec = NULL;
    return;
}// deallocateHypreData

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBTK::CCDivGradHypreLevelSolver>;

//////////////////////////////////////////////////////////////////////////////
