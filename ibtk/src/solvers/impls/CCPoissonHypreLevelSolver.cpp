// Filename: CCPoissonHypreLevelSolver.cpp
// Created on 30 May 2005 by Boyce Griffith
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
#include <algorithm>
#include <functional>
#include <map>
#include <ostream>
#include <string>
#include <vector>

#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellDataFactory.h"
#include "CellIndex.h"
#include "HYPRE_struct_ls.h"
#include "HYPRE_struct_mv.h"
#include "Index.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchDescriptor.h"
#include "PatchGeometry.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PoissonSpecifications.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "SideDataFactory.h"
#include "SideIndex.h"
#include "VariableDatabase.h"
#include "ibtk/CCPoissonHypreLevelSolver.h"
#include "ibtk/GeneralSolver.h"
#include "ibtk/PoissonUtilities.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "mpi.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

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
enum HypreStructRAPType
{
    RAP_TYPE_GALERKIN = 0,
    RAP_TYPE_NON_GALERKIN_PARFLOW = 1,
    RAP_TYPE_GALERKIN_GENERAL = 2
};

enum HypreStructRelaxType
{
    RELAX_TYPE_JACOBI = 0,
    RELAX_TYPE_WEIGHTED_JACOBI = 1,
    RELAX_TYPE_RB_GAUSS_SEIDEL = 2,
    RELAX_TYPE_RB_GAUSS_SEIDEL_NONSYMMETRIC = 3
};

struct IndexComp : std::binary_function<Index<NDIM>, Index<NDIM>, bool>
{
    bool operator()(const Index<NDIM>& lhs, const Index<NDIM>& rhs) const
    {
        return (lhs(0) < rhs(0)
#if (NDIM > 1)
                ||
                (lhs(0) == rhs(0) && lhs(1) < rhs(1))
#if (NDIM > 2)
                ||
                (lhs(0) == rhs(0) && lhs(1) == rhs(1) && lhs(2) < rhs(2))
#endif
#endif
                    );
    } // operator()
};
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CCPoissonHypreLevelSolver::CCPoissonHypreLevelSolver(const std::string& object_name,
                                                     Pointer<Database> input_db,
                                                     const std::string& /*default_options_prefix*/)
    : d_hierarchy(),
      d_level_num(-1),
      d_grid_aligned_anisotropy(true),
      d_depth(0),
      d_grid(NULL),
      d_stencil(NULL),
      d_matrices(),
      d_rhs_vecs(),
      d_sol_vecs(),
      d_solvers(),
      d_preconds(),
      d_solver_type("PFMG"),
      d_precond_type("none"),
      d_rel_change(0),
      d_num_pre_relax_steps(1),
      d_num_post_relax_steps(1),
      d_memory_use(0),
      d_rap_type(RAP_TYPE_GALERKIN),
      d_relax_type(RELAX_TYPE_WEIGHTED_JACOBI),
      d_skip_relax(1),
      d_two_norm(1)
{
    if (NDIM == 1 || NDIM > 3)
    {
        TBOX_ERROR(d_object_name << "::CCPoissonHypreLevelSolver()"
                                 << "  hypre solvers are only provided for 2D and 3D problems"
                                 << std::endl);
    }

    // Setup default options.
    GeneralSolver::init(object_name, /*homogeneous_bc*/ false);
    d_initial_guess_nonzero = false;
    d_rel_residual_tol = 1.0e-5;
    d_abs_residual_tol = 1.0e-50;
    d_max_iterations = 25;

    // Get values from the input database.
    if (input_db)
    {
        if (input_db->keyExists("enable_logging")) d_enable_logging = input_db->getBool("enable_logging");
        if (input_db->keyExists("solver_type")) d_solver_type = input_db->getString("solver_type");
        if (input_db->keyExists("precond_type")) d_precond_type = input_db->getString("precond_type");
        if (input_db->keyExists("max_iterations")) d_max_iterations = input_db->getInteger("max_iterations");
        if (input_db->keyExists("abs_residual_tol")) d_abs_residual_tol = input_db->getDouble("abs_residual_tol");
        if (input_db->keyExists("rel_residual_tol")) d_rel_residual_tol = input_db->getDouble("rel_residual_tol");
        if (input_db->keyExists("initial_guess_nonzero"))
            d_initial_guess_nonzero = input_db->getBool("initial_guess_nonzero");
        if (input_db->keyExists("rel_change")) d_rel_change = input_db->getInteger("rel_change");

        if (d_solver_type == "SMG" || d_precond_type == "SMG" || d_solver_type == "PFMG" || d_precond_type == "PFMG")
        {
            if (input_db->keyExists("num_pre_relax_steps"))
                d_num_pre_relax_steps = input_db->getInteger("num_pre_relax_steps");
            if (input_db->keyExists("num_post_relax_steps"))
                d_num_post_relax_steps = input_db->getInteger("num_post_relax_steps");
        }

        if (d_solver_type == "SMG" || d_precond_type == "SMG")
        {
            if (input_db->keyExists("memory_use")) d_memory_use = input_db->getInteger("memory_use");
        }

        if (d_solver_type == "PFMG" || d_precond_type == "PFMG")
        {
            if (input_db->keyExists("rap_type")) d_rap_type = input_db->getInteger("rap_type");
            if (input_db->keyExists("relax_type")) d_relax_type = input_db->getInteger("relax_type");
            if (input_db->keyExists("skip_relax")) d_skip_relax = input_db->getInteger("skip_relax");
        }

        if (d_solver_type == "PCG")
        {
            if (input_db->keyExists("two_norm")) d_two_norm = input_db->getInteger("two_norm");
        }
    }

    // Setup Timers.
    IBTK_DO_ONCE(t_solve_system =
                     TimerManager::getManager()->getTimer("IBTK::CCPoissonHypreLevelSolver::solveSystem()");
                 t_solve_system_hypre =
                     TimerManager::getManager()->getTimer("IBTK::CCPoissonHypreLevelSolver::solveSystem()[hypre]");
                 t_initialize_solver_state =
                     TimerManager::getManager()->getTimer("IBTK::CCPoissonHypreLevelSolver::initializeSolverState()");
                 t_deallocate_solver_state =
                     TimerManager::getManager()->getTimer("IBTK::CCPoissonHypreLevelSolver::deallocateSolverState()"););
    return;
} // CCPoissonHypreLevelSolver

CCPoissonHypreLevelSolver::~CCPoissonHypreLevelSolver()
{
    if (d_is_initialized) deallocateSolverState();
    return;
} // ~CCPoissonHypreLevelSolver

bool
CCPoissonHypreLevelSolver::solveSystem(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& b)
{
    IBTK_TIMER_START(t_solve_system);

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x, b);

    // Ensure the initial guess is zero when appropriate.  (hypre does not
    // reliably honor the SetZeroGuess settings.)
    if (!d_initial_guess_nonzero) x.setToScalar(0.0, /*interior_only*/ false);

    // Solve the system using the hypre solver.
    static const int comp = 0;
    const int x_idx = x.getComponentDescriptorIndex(comp);
    const int b_idx = b.getComponentDescriptorIndex(comp);
    const bool converged = solveSystem(x_idx, b_idx);

    // Log solver info.
    if (d_enable_logging)
    {
        plog << d_object_name << "::solveSystem(): solver " << (converged ? "converged" : "diverged") << "\n"
             << "iterations = " << d_current_iterations << "\n"
             << "residual norm = " << d_current_residual_norm << std::endl;
    }

    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();

    IBTK_TIMER_STOP(t_solve_system);
    return converged;
} // solveSystem

void
CCPoissonHypreLevelSolver::initializeSolverState(const SAMRAIVectorReal<NDIM, double>& x,
                                                 const SAMRAIVectorReal<NDIM, double>& b)
{
    IBTK_TIMER_START(t_initialize_solver_state);

#if !defined(NDEBUG)
    // Rudimentary error checking.
    if (x.getNumberOfComponents() != b.getNumberOfComponents())
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  vectors must have the same number of components"
                                 << std::endl);
    }

    const Pointer<PatchHierarchy<NDIM> >& patch_hierarchy = x.getPatchHierarchy();
    if (patch_hierarchy != b.getPatchHierarchy())
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  vectors must have the same hierarchy"
                                 << std::endl);
    }

    const int coarsest_ln = x.getCoarsestLevelNumber();
    if (coarsest_ln < 0)
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  coarsest level number must not be negative"
                                 << std::endl);
    }
    if (coarsest_ln != b.getCoarsestLevelNumber())
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  vectors must have same coarsest level number"
                                 << std::endl);
    }

    const int finest_ln = x.getFinestLevelNumber();
    if (finest_ln < coarsest_ln)
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  finest level number must be >= coarsest level number"
                                 << std::endl);
    }
    if (finest_ln != b.getFinestLevelNumber())
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  vectors must have same finest level number"
                                 << std::endl);
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!patch_hierarchy->getPatchLevel(ln))
        {
            TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                     << "  hierarchy level "
                                     << ln
                                     << " does not exist"
                                     << std::endl);
        }
    }

    if (coarsest_ln != finest_ln)
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  coarsest_ln != finest_ln in CCPoissonHypreLevelSolver"
                                 << std::endl);
    }
#else
    NULL_USE(b);
#endif
    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized) deallocateSolverState();

    // Get the hierarchy information.
    d_hierarchy = x.getPatchHierarchy();
    d_level_num = x.getCoarsestLevelNumber();
    TBOX_ASSERT(d_level_num == x.getFinestLevelNumber());
    d_level = d_hierarchy->getPatchLevel(d_level_num);
    if (d_level_num > 0)
    {
        d_cf_boundary = new CoarseFineBoundary<NDIM>(*d_hierarchy, d_level_num, IntVector<NDIM>(1));
    }

    // Allocate and initialize the hypre data structures.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int x_idx = x.getComponentDescriptorIndex(0);
    Pointer<CellDataFactory<NDIM, double> > x_fac = var_db->getPatchDescriptor()->getPatchDataFactory(x_idx);
    d_depth = x_fac->getDefaultDepth();
    if (d_poisson_spec.dIsConstant())
    {
        d_grid_aligned_anisotropy = true;
    }
    else
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<SideDataFactory<NDIM, double> > pdat_factory =
            var_db->getPatchDescriptor()->getPatchDataFactory(d_poisson_spec.getDPatchDataId());
#if !defined(NDEBUG)
        TBOX_ASSERT(pdat_factory);
#endif
        d_grid_aligned_anisotropy = pdat_factory->getDefaultDepth() == 1;
    }
    allocateHypreData();
    if (d_grid_aligned_anisotropy)
    {
        setMatrixCoefficients_aligned();
    }
    else
    {
        setMatrixCoefficients_nonaligned();
    }
    setupHypreSolver();

    // Indicate that the solver is initialized.
    d_is_initialized = true;

    IBTK_TIMER_STOP(t_initialize_solver_state);
    return;
} // initializeSolverState

void
CCPoissonHypreLevelSolver::deallocateSolverState()
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_deallocate_solver_state);

    // Deallocate the hypre data structures.
    destroyHypreSolver();
    deallocateHypreData();

    // Indicate that the solver is NOT initialized.
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_solver_state);
    return;
} // deallocateSolverState

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CCPoissonHypreLevelSolver::allocateHypreData()
{
    // Get the MPI communicator.
    MPI_Comm communicator = SAMRAI_MPI::getCommunicator();

    // Setup the hypre grid.
    Pointer<CartesianGridGeometry<NDIM> > grid_geometry = d_hierarchy->getGridGeometry();
    const IntVector<NDIM>& ratio = d_level->getRatio();
    const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift(ratio);

    HYPRE_StructGridCreate(communicator, NDIM, &d_grid);
    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        const Box<NDIM>& patch_box = d_level->getPatch(p())->getBox();
        Index<NDIM> lower = patch_box.lower();
        Index<NDIM> upper = patch_box.upper();
        HYPRE_StructGridSetExtents(d_grid, lower, upper);
    }

    int hypre_periodic_shift[3];
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        hypre_periodic_shift[d] = periodic_shift(d);
    }
    for (int d = NDIM; d < 3; ++d)
    {
        hypre_periodic_shift[d] = 0;
    }
    HYPRE_StructGridSetPeriodic(d_grid, hypre_periodic_shift);
    HYPRE_StructGridAssemble(d_grid);

    // Allocate stencil data and set stencil offsets.
    if (d_grid_aligned_anisotropy)
    {
        static const int stencil_sz = 2 * NDIM + 1;
        d_stencil_offsets.resize(stencil_sz);
        std::fill(d_stencil_offsets.begin(), d_stencil_offsets.end(), Index<NDIM>(0));
        for (unsigned int axis = 0, stencil_index = 1; axis < NDIM; ++axis)
        {
            for (int side = 0; side <= 1; ++side, ++stencil_index)
            {
                d_stencil_offsets[stencil_index](axis) = (side == 0 ? -1 : +1);
            }
        }
        HYPRE_StructStencilCreate(NDIM, stencil_sz, &d_stencil);
        for (int s = 0; s < stencil_sz; ++s)
        {
            HYPRE_StructStencilSetElement(d_stencil, s, d_stencil_offsets[s]);
        }
    }
    else
    {
        static const int stencil_sz = (NDIM == 2) ? 9 : 19;
        d_stencil_offsets.resize(stencil_sz);
        std::fill(d_stencil_offsets.begin(), d_stencil_offsets.end(), Index<NDIM>(0));
        int stencil_index = 0;
#if (NDIM == 3)
        for (int z_offset = -1; z_offset <= 1; ++z_offset)
#endif
        {
            for (int y_offset = -1; y_offset <= 1; ++y_offset)
            {
                for (int x_offset = -1; x_offset <= 1; ++x_offset)
                {
#if (NDIM == 3)
                    // No full-corner coupling in 3D.
                    if (x_offset == 0 || y_offset == 0 || z_offset == 0)
                    {
#endif
                        d_stencil_offsets[stencil_index](0) = x_offset;
                        d_stencil_offsets[stencil_index](1) = y_offset;
#if (NDIM == 3)
                        d_stencil_offsets[stencil_index](2) = z_offset;
#endif
                        ++stencil_index;
#if (NDIM == 3)
                    }
#endif
                }
            }
        }

        HYPRE_StructStencilCreate(NDIM, stencil_sz, &d_stencil);
        for (int s = 0; s < stencil_sz; ++s)
        {
            HYPRE_StructStencilSetElement(d_stencil, s, d_stencil_offsets[s]);
        }
    }

// Allocate the hypre matrices.
#if (NDIM == 2)
    int full_ghosts[2 * 3] = { 1, 1, 1, 1, 0, 0 };
#endif
#if (NDIM == 3)
    int full_ghosts[2 * 3] = { 1, 1, 1, 1, 1, 1 };
#endif
    int no_ghosts[2 * 3] = { 0, 0, 0, 0, 0, 0 };
    d_matrices.resize(d_depth);
    for (unsigned int k = 0; k < d_depth; ++k)
    {
        HYPRE_StructMatrixCreate(communicator, d_grid, d_stencil, &d_matrices[k]);
        HYPRE_StructMatrixSetNumGhost(d_matrices[k], full_ghosts);
        HYPRE_StructMatrixSetSymmetric(d_matrices[k], 0);
        HYPRE_StructMatrixInitialize(d_matrices[k]);
    }

    // Allocate the hypre vectors.
    d_sol_vecs.resize(d_depth);
    d_rhs_vecs.resize(d_depth);
    for (unsigned int k = 0; k < d_depth; ++k)
    {
        HYPRE_StructVectorCreate(communicator, d_grid, &d_sol_vecs[k]);
        HYPRE_StructVectorSetNumGhost(d_sol_vecs[k], full_ghosts);
        HYPRE_StructVectorInitialize(d_sol_vecs[k]);

        HYPRE_StructVectorCreate(communicator, d_grid, &d_rhs_vecs[k]);
        HYPRE_StructVectorSetNumGhost(d_rhs_vecs[k], no_ghosts);
        HYPRE_StructVectorInitialize(d_rhs_vecs[k]);
    }
    return;
} // allocateHypreData

void
CCPoissonHypreLevelSolver::setMatrixCoefficients_aligned()
{
    // Set matrix entries and copy them to the hypre matrix structures.
    const int stencil_sz = static_cast<int>(d_stencil_offsets.size());
    std::vector<int> stencil_indices(stencil_sz);
    for (int i = 0; i < stencil_sz; ++i)
    {
        stencil_indices[i] = i;
    }
    std::vector<double> mat_vals(stencil_sz, 0.0);
    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = d_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        CellData<NDIM, double> matrix_coefs(patch_box, stencil_sz, IntVector<NDIM>(0));
        for (unsigned int k = 0; k < d_depth; ++k)
        {
            PoissonUtilities::computeMatrixCoefficients(
                matrix_coefs, patch, d_stencil_offsets, d_poisson_spec, d_bc_coefs[k], d_solution_time);
            for (Box<NDIM>::Iterator b(patch_box); b; b++)
            {
                Index<NDIM> i = b();
                for (int j = 0; j < stencil_sz; ++j)
                {
                    mat_vals[j] = matrix_coefs(i, j);
                }
                HYPRE_StructMatrixSetValues(d_matrices[k], i, stencil_sz, &stencil_indices[0], &mat_vals[0]);
            }
        }
    }

    // Assemble the hypre matrices.
    for (unsigned int k = 0; k < d_depth; ++k)
    {
        HYPRE_StructMatrixAssemble(d_matrices[k]);
    }
    return;
} // setMatrixCoefficients_aligned

void
CCPoissonHypreLevelSolver::setMatrixCoefficients_nonaligned()
{
    static const IntVector<NDIM> no_ghosts = 0;
    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = d_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();

        // Compute all matrix coefficients.
        //
        // NOTE: Here we assume that no flux boundary conditions are imposed at
        // the physical domain.
        Pointer<CellData<NDIM, double> > C_data;
        if (!d_poisson_spec.cIsZero() && !d_poisson_spec.cIsConstant())
        {
            C_data = patch->getPatchData(d_poisson_spec.getCPatchDataId());
            if (!C_data)
            {
                TBOX_ERROR(d_object_name << "::setMatrixCoefficients_nonaligned()\n"
                                         << "  to solve (C u + div D grad u) = f with non-constant C,\n"
                                         << "  C must be cell-centered double precision data"
                                         << std::endl);
            }
        }
        else
        {
            C_data = new CellData<NDIM, double>(patch_box, 1, no_ghosts);
            if (d_poisson_spec.cIsZero())
                C_data->fill(0.0);
            else
                C_data->fill(d_poisson_spec.getCConstant());
        }

        Pointer<SideData<NDIM, double> > D_data;
        if (!d_poisson_spec.dIsConstant())
        {
            D_data = patch->getPatchData(d_poisson_spec.getDPatchDataId());
        }

        if (!D_data)
        {
            TBOX_ERROR(d_object_name << "::setMatrixCoefficients_nonaligned()\n"
                                     << "  to solve C u + div D grad u = f with non-constant D,\n"
                                     << "  D must be side-centered double precision data"
                                     << std::endl);
        }

        // Setup the finite difference stencil.
        static const int stencil_sz = (NDIM == 2 ? 9 : 19);
        int stencil_indices[stencil_sz];
        for (int i = 0; i < stencil_sz; ++i)
        {
            stencil_indices[i] = i;
        }

        std::map<Index<NDIM>, int, IndexComp> stencil_index_map;
        int stencil_index = 0;
#if (NDIM == 3)
        for (int z_offset = -1; z_offset <= 1; ++z_offset)
#endif
        {
            for (int y_offset = -1; y_offset <= 1; ++y_offset)
            {
                for (int x_offset = -1; x_offset <= 1; ++x_offset)
                {
#if (NDIM == 3)
                    // No full-corner coupling in 3D.
                    if (x_offset == 0 || y_offset == 0 || z_offset == 0)
                    {
#endif
#if (NDIM == 2)
                        const Index<NDIM> i(x_offset, y_offset);
#endif
#if (NDIM == 3)
                        const Index<NDIM> i(x_offset, y_offset, z_offset);
#endif
                        stencil_index_map[i] = stencil_index++;
#if (NDIM == 3)
                    }
#endif
                }
            }
        }

        // Set the matrix coefficients to correspond to a second-order accurate
        // finite difference stencil for the Laplace operator.
        for (Box<NDIM>::Iterator b(patch_box); b; b++)
        {
            Index<NDIM> i = b();
            static const Index<NDIM> i_stencil_center(0);
            const int stencil_center = stencil_index_map[i_stencil_center];

            std::vector<double> mat_vals(stencil_sz, 0.0);
            mat_vals[stencil_center] = (*C_data)(i);

            // The grid aligned part of the stencil (normal derivatives).
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const double& h = dx[axis];
                {
                    Index<NDIM> i_stencil_lower(0);
                    --i_stencil_lower[axis];
                    const int stencil_lower = stencil_index_map[i_stencil_lower];

                    const SideIndex<NDIM> ilower(i, axis, SideIndex<NDIM>::Lower);
                    const double& D_lower = (*D_data)(ilower, axis);
                    mat_vals[stencil_lower] += D_lower / (h * h);
                    mat_vals[stencil_center] -= D_lower / (h * h);
                }
                {
                    Index<NDIM> i_stencil_upper(0);
                    ++i_stencil_upper[axis];
                    const int stencil_upper = stencil_index_map[i_stencil_upper];

                    const SideIndex<NDIM> iupper(i, axis, SideIndex<NDIM>::Upper);
                    const double& D_upper = (*D_data)(iupper, axis);
                    mat_vals[stencil_upper] += D_upper / (h * h);
                    mat_vals[stencil_center] -= D_upper / (h * h);
                }
            }

            // The non-grid aligned part of the stencil (transverse derivatives).
            for (unsigned int norm_axis = 0; norm_axis < NDIM; ++norm_axis)
            {
                const double& norm_h = dx[norm_axis];
                for (unsigned int trans_axis = 0; trans_axis < NDIM; ++trans_axis)
                {
                    if (norm_axis == trans_axis) break;
                    const double& trans_h = dx[trans_axis];

                    static const int lower = 0;
                    const bool skip_lower_side_fluxes =
                        patch->getPatchGeometry()->getTouchesRegularBoundary(norm_axis, lower) &&
                        i(norm_axis) == patch_box.lower()(norm_axis);

                    static const int upper = 1;
                    const bool skip_upper_side_fluxes =
                        patch->getPatchGeometry()->getTouchesRegularBoundary(norm_axis, upper) &&
                        i(norm_axis) == patch_box.upper()(norm_axis);

                    // Lower side fluxes.
                    if (!skip_lower_side_fluxes)
                    {
                        const SideIndex<NDIM> ilower(i, norm_axis, SideIndex<NDIM>::Lower);
                        for (int norm_shift = -1; norm_shift <= 0; ++norm_shift)
                        {
                            for (int trans_shift = -1; trans_shift <= 1; trans_shift += 2)
                            {
                                Index<NDIM> i_stencil(0);
                                i_stencil[norm_axis] += norm_shift;
                                i_stencil[trans_axis] += trans_shift;
                                const int stencil_index = stencil_index_map[i_stencil];
                                if (trans_shift == 1)
                                {
                                    mat_vals[stencil_index] -= 0.25 * (*D_data)(ilower, norm_axis) / (norm_h * trans_h);
                                }
                                else
                                {
                                    mat_vals[stencil_index] += 0.25 * (*D_data)(ilower, norm_axis) / (norm_h * trans_h);
                                }
                            }
                        }
                    }

                    // Upper side fluxes.
                    if (!skip_upper_side_fluxes)
                    {
                        const SideIndex<NDIM> iupper(i, norm_axis, SideIndex<NDIM>::Upper);
                        for (int norm_shift = 0; norm_shift <= 1; ++norm_shift)
                        {
                            for (int trans_shift = -1; trans_shift <= 1; trans_shift += 2)
                            {
                                Index<NDIM> i_stencil(0);
                                i_stencil[norm_axis] += norm_shift;
                                i_stencil[trans_axis] += trans_shift;
                                const int stencil_index = stencil_index_map[i_stencil];
                                if (trans_shift == 1)
                                {
                                    mat_vals[stencil_index] += 0.25 * (*D_data)(iupper, norm_axis) / (norm_h * trans_h);
                                }
                                else
                                {
                                    mat_vals[stencil_index] -= 0.25 * (*D_data)(iupper, norm_axis) / (norm_h * trans_h);
                                }
                            }
                        }
                    }
                }
            }

            for (unsigned int k = 0; k < d_depth; ++k)
            {
                HYPRE_StructMatrixSetValues(d_matrices[k], i, stencil_sz, stencil_indices, &mat_vals[0]);
            }
        }
    }

    // Assemble the hypre matrices.
    for (unsigned int k = 0; k < d_depth; ++k)
    {
        HYPRE_StructMatrixAssemble(d_matrices[k]);
    }
    return;
} // setMatrixCoefficients_nonaligned

void
CCPoissonHypreLevelSolver::setupHypreSolver()
{
    // Get the MPI communicator.
    MPI_Comm communicator = SAMRAI_MPI::getCommunicator();

    d_solvers.resize(d_depth);
    d_preconds.resize(d_depth);
    for (unsigned int k = 0; k < d_depth; ++k)
    {
        // When using a Krylov method, setup the preconditioner.
        if (d_solver_type == "PCG" || d_solver_type == "GMRES" || d_solver_type == "FlexGMRES" ||
            d_solver_type == "LGMRES" || d_solver_type == "BiCGSTAB")
        {
            if (d_precond_type == "PFMG")
            {
                HYPRE_StructPFMGCreate(communicator, &d_preconds[k]);
                HYPRE_StructPFMGSetMaxIter(d_preconds[k], 1);
                HYPRE_StructPFMGSetTol(d_preconds[k], 0.0);
                HYPRE_StructPFMGSetZeroGuess(d_preconds[k]);
                HYPRE_StructPFMGSetRAPType(d_preconds[k], d_rap_type);
                HYPRE_StructPFMGSetRelaxType(d_preconds[k], d_relax_type);
                HYPRE_StructPFMGSetNumPreRelax(d_preconds[k], d_num_pre_relax_steps);
                HYPRE_StructPFMGSetNumPostRelax(d_preconds[k], d_num_post_relax_steps);
                HYPRE_StructPFMGSetSkipRelax(d_preconds[k], d_skip_relax);
            }
            else if (d_precond_type == "SMG")
            {
                HYPRE_StructSMGCreate(communicator, &d_preconds[k]);
                HYPRE_StructSMGSetMaxIter(d_preconds[k], 1);
                HYPRE_StructSMGSetTol(d_preconds[k], 0.0);
                HYPRE_StructSMGSetZeroGuess(d_preconds[k]);
                HYPRE_StructSMGSetMemoryUse(d_preconds[k], d_memory_use);
                HYPRE_StructSMGSetNumPreRelax(d_preconds[k], d_num_pre_relax_steps);
                HYPRE_StructSMGSetNumPostRelax(d_preconds[k], d_num_post_relax_steps);
            }
            else if (d_precond_type == "Jacobi")
            {
                HYPRE_StructJacobiCreate(communicator, &d_preconds[k]);
                HYPRE_StructJacobiSetMaxIter(d_preconds[k], 2);
                HYPRE_StructJacobiSetTol(d_preconds[k], 0.0);
                HYPRE_StructJacobiSetZeroGuess(d_preconds[k]);
            }
        }

        // Setup the solver.
        if (d_solver_type == "PFMG")
        {
            HYPRE_StructPFMGCreate(communicator, &d_solvers[k]);
            HYPRE_StructPFMGSetMaxIter(d_solvers[k], d_max_iterations);
            HYPRE_StructPFMGSetTol(d_solvers[k], d_rel_residual_tol);
            HYPRE_StructPFMGSetRelChange(d_solvers[k], d_rel_change);
            HYPRE_StructPFMGSetRAPType(d_solvers[k], d_rap_type);
            HYPRE_StructPFMGSetRelaxType(d_solvers[k], d_relax_type);
            HYPRE_StructPFMGSetNumPreRelax(d_solvers[k], d_num_pre_relax_steps);
            HYPRE_StructPFMGSetNumPostRelax(d_solvers[k], d_num_post_relax_steps);
            HYPRE_StructPFMGSetSkipRelax(d_solvers[k], d_skip_relax);
            if (d_initial_guess_nonzero)
            {
                HYPRE_StructPFMGSetNonZeroGuess(d_solvers[k]);
            }
            else
            {
                HYPRE_StructPFMGSetZeroGuess(d_solvers[k]);
            }
            HYPRE_StructPFMGSetup(d_solvers[k], d_matrices[k], d_rhs_vecs[k], d_sol_vecs[k]);
        }
        else if (d_solver_type == "SMG")
        {
            HYPRE_StructSMGCreate(communicator, &d_solvers[k]);
            HYPRE_StructSMGSetMaxIter(d_solvers[k], d_max_iterations);
            HYPRE_StructSMGSetTol(d_solvers[k], d_rel_residual_tol);
            HYPRE_StructSMGSetRelChange(d_solvers[k], d_rel_change);
            HYPRE_StructSMGSetMemoryUse(d_solvers[k], d_memory_use);
            HYPRE_StructSMGSetNumPreRelax(d_solvers[k], d_num_pre_relax_steps);
            HYPRE_StructSMGSetNumPostRelax(d_solvers[k], d_num_post_relax_steps);
            if (d_initial_guess_nonzero)
            {
                HYPRE_StructSMGSetNonZeroGuess(d_solvers[k]);
            }
            else
            {
                HYPRE_StructSMGSetZeroGuess(d_solvers[k]);
            }
            HYPRE_StructSMGSetup(d_solvers[k], d_matrices[k], d_rhs_vecs[k], d_sol_vecs[k]);
        }
        else if (d_solver_type == "PCG")
        {
            HYPRE_StructPCGCreate(communicator, &d_solvers[k]);
            HYPRE_StructPCGSetMaxIter(d_solvers[k], d_max_iterations);
            HYPRE_StructPCGSetTol(d_solvers[k], d_rel_residual_tol);
            HYPRE_StructPCGSetAbsoluteTol(d_solvers[k], d_abs_residual_tol);
            HYPRE_StructPCGSetTwoNorm(d_solvers[k], d_two_norm);
            HYPRE_StructPCGSetRelChange(d_solvers[k], d_rel_change);
            if (d_precond_type == "PFMG")
            {
                HYPRE_StructPCGSetPrecond(d_solvers[k], HYPRE_StructPFMGSolve, HYPRE_StructPFMGSetup, d_preconds[k]);
            }
            else if (d_precond_type == "SMG")
            {
                HYPRE_StructPCGSetPrecond(d_solvers[k], HYPRE_StructSMGSolve, HYPRE_StructSMGSetup, d_preconds[k]);
            }
            else if (d_precond_type == "Jacobi")
            {
                HYPRE_StructPCGSetPrecond(
                    d_solvers[k], HYPRE_StructJacobiSolve, HYPRE_StructJacobiSetup, d_preconds[k]);
            }
            else if (d_precond_type == "diagonal_scaling")
            {
                d_preconds[k] = NULL;
                HYPRE_StructPCGSetPrecond(
                    d_solvers[k], HYPRE_StructDiagScale, HYPRE_StructDiagScaleSetup, d_preconds[k]);
            }
            else if (d_precond_type == "none")
            {
                d_preconds[k] = NULL;
            }
            else
            {
                TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                         << "  unknown preconditioner type: "
                                         << d_precond_type
                                         << std::endl);
            }
            HYPRE_StructPCGSetup(d_solvers[k], d_matrices[k], d_rhs_vecs[k], d_sol_vecs[k]);
        }
        else if (d_solver_type == "GMRES")
        {
            HYPRE_StructGMRESCreate(communicator, &d_solvers[k]);
            HYPRE_StructGMRESSetMaxIter(d_solvers[k], d_max_iterations);
            HYPRE_StructGMRESSetTol(d_solvers[k], d_rel_residual_tol);
            HYPRE_StructGMRESSetAbsoluteTol(d_solvers[k], d_abs_residual_tol);
            if (d_precond_type == "PFMG")
            {
                HYPRE_StructGMRESSetPrecond(d_solvers[k], HYPRE_StructPFMGSolve, HYPRE_StructPFMGSetup, d_preconds[k]);
            }
            else if (d_precond_type == "SMG")
            {
                HYPRE_StructGMRESSetPrecond(d_solvers[k], HYPRE_StructSMGSolve, HYPRE_StructSMGSetup, d_preconds[k]);
            }
            else if (d_precond_type == "Jacobi")
            {
                HYPRE_StructGMRESSetPrecond(
                    d_solvers[k], HYPRE_StructJacobiSolve, HYPRE_StructJacobiSetup, d_preconds[k]);
            }
            else if (d_precond_type == "diagonal_scaling")
            {
                d_preconds[k] = NULL;
                HYPRE_StructGMRESSetPrecond(
                    d_solvers[k], HYPRE_StructDiagScale, HYPRE_StructDiagScaleSetup, d_preconds[k]);
            }
            else if (d_precond_type == "none")
            {
                d_preconds[k] = NULL;
            }
            else
            {
                TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                         << "  unknown preconditioner type: "
                                         << d_precond_type
                                         << std::endl);
            }
            HYPRE_StructGMRESSetup(d_solvers[k], d_matrices[k], d_rhs_vecs[k], d_sol_vecs[k]);
        }
        else if (d_solver_type == "FlexGMRES")
        {
            HYPRE_StructFlexGMRESCreate(communicator, &d_solvers[k]);
            HYPRE_StructFlexGMRESSetMaxIter(d_solvers[k], d_max_iterations);
            HYPRE_StructFlexGMRESSetTol(d_solvers[k], d_rel_residual_tol);
            HYPRE_StructFlexGMRESSetAbsoluteTol(d_solvers[k], d_abs_residual_tol);
            if (d_precond_type == "PFMG")
            {
                HYPRE_StructFlexGMRESSetPrecond(
                    d_solvers[k], HYPRE_StructPFMGSolve, HYPRE_StructPFMGSetup, d_preconds[k]);
            }
            else if (d_precond_type == "SMG")
            {
                HYPRE_StructFlexGMRESSetPrecond(
                    d_solvers[k], HYPRE_StructSMGSolve, HYPRE_StructSMGSetup, d_preconds[k]);
            }
            else if (d_precond_type == "Jacobi")
            {
                HYPRE_StructFlexGMRESSetPrecond(
                    d_solvers[k], HYPRE_StructJacobiSolve, HYPRE_StructJacobiSetup, d_preconds[k]);
            }
            else if (d_precond_type == "diagonal_scaling")
            {
                d_preconds[k] = NULL;
                HYPRE_StructFlexGMRESSetPrecond(
                    d_solvers[k], HYPRE_StructDiagScale, HYPRE_StructDiagScaleSetup, d_preconds[k]);
            }
            else if (d_precond_type == "none")
            {
                d_preconds[k] = NULL;
            }
            else
            {
                TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                         << "  unknown preconditioner type: "
                                         << d_precond_type
                                         << std::endl);
            }
            HYPRE_StructFlexGMRESSetup(d_solvers[k], d_matrices[k], d_rhs_vecs[k], d_sol_vecs[k]);
        }
        else if (d_solver_type == "LGMRES")
        {
            HYPRE_StructLGMRESCreate(communicator, &d_solvers[k]);
            HYPRE_StructLGMRESSetMaxIter(d_solvers[k], d_max_iterations);
            HYPRE_StructLGMRESSetTol(d_solvers[k], d_rel_residual_tol);
            HYPRE_StructLGMRESSetAbsoluteTol(d_solvers[k], d_abs_residual_tol);
            if (d_precond_type == "PFMG")
            {
                HYPRE_StructLGMRESSetPrecond(d_solvers[k], HYPRE_StructPFMGSolve, HYPRE_StructPFMGSetup, d_preconds[k]);
            }
            else if (d_precond_type == "SMG")
            {
                HYPRE_StructLGMRESSetPrecond(d_solvers[k], HYPRE_StructSMGSolve, HYPRE_StructSMGSetup, d_preconds[k]);
            }
            else if (d_precond_type == "Jacobi")
            {
                HYPRE_StructLGMRESSetPrecond(
                    d_solvers[k], HYPRE_StructJacobiSolve, HYPRE_StructJacobiSetup, d_preconds[k]);
            }
            else if (d_precond_type == "diagonal_scaling")
            {
                d_preconds[k] = NULL;
                HYPRE_StructLGMRESSetPrecond(
                    d_solvers[k], HYPRE_StructDiagScale, HYPRE_StructDiagScaleSetup, d_preconds[k]);
            }
            else if (d_precond_type == "none")
            {
                d_preconds[k] = NULL;
            }
            else
            {
                TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                         << "  unknown preconditioner type: "
                                         << d_precond_type
                                         << std::endl);
            }
            HYPRE_StructLGMRESSetup(d_solvers[k], d_matrices[k], d_rhs_vecs[k], d_sol_vecs[k]);
        }
        else if (d_solver_type == "BiCGSTAB")
        {
            HYPRE_StructBiCGSTABCreate(communicator, &d_solvers[k]);
            HYPRE_StructBiCGSTABSetMaxIter(d_solvers[k], d_max_iterations);
            HYPRE_StructBiCGSTABSetTol(d_solvers[k], d_rel_residual_tol);
            HYPRE_StructBiCGSTABSetAbsoluteTol(d_solvers[k], d_abs_residual_tol);
            if (d_precond_type == "PFMG")
            {
                HYPRE_StructBiCGSTABSetPrecond(
                    d_solvers[k], HYPRE_StructPFMGSolve, HYPRE_StructPFMGSetup, d_preconds[k]);
            }
            else if (d_precond_type == "SMG")
            {
                HYPRE_StructBiCGSTABSetPrecond(d_solvers[k], HYPRE_StructSMGSolve, HYPRE_StructSMGSetup, d_preconds[k]);
            }
            else if (d_precond_type == "Jacobi")
            {
                HYPRE_StructBiCGSTABSetPrecond(
                    d_solvers[k], HYPRE_StructJacobiSolve, HYPRE_StructJacobiSetup, d_preconds[k]);
            }
            else if (d_precond_type == "diagonal_scaling")
            {
                d_preconds[k] = NULL;
                HYPRE_StructBiCGSTABSetPrecond(
                    d_solvers[k], HYPRE_StructDiagScale, HYPRE_StructDiagScaleSetup, d_preconds[k]);
            }
            else if (d_precond_type == "none")
            {
                d_preconds[k] = NULL;
            }
            else
            {
                TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                         << "  unknown preconditioner type: "
                                         << d_precond_type
                                         << std::endl);
            }
            HYPRE_StructBiCGSTABSetup(d_solvers[k], d_matrices[k], d_rhs_vecs[k], d_sol_vecs[k]);
        }
        else
        {
            TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                     << "  unknown solver type: "
                                     << d_solver_type
                                     << std::endl);
        }
    }
    return;
} // setupHypreSolver

bool
CCPoissonHypreLevelSolver::solveSystem(const int x_idx, const int b_idx)
{
    const bool level_zero = (d_level_num == 0);

    // Modify right-hand-side data to account for boundary conditions and copy
    // solution and right-hand-side data to hypre structures.
    const IntVector<NDIM> ghosts = 1;
    const IntVector<NDIM> no_ghosts = 0;
    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = d_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

        // Copy the solution data into the hypre vector, including ghost cell
        // values
        const Box<NDIM> x_ghost_box = Box<NDIM>::grow(patch_box, 1);
        Pointer<CellData<NDIM, double> > x_data = patch->getPatchData(x_idx);
        copyToHypre(d_sol_vecs, *x_data, x_ghost_box);

        // Modify the right-hand-side data to account for any inhomogeneous
        // boundary conditions and copy the right-hand-side into the hypre
        // vector.
        Pointer<CellData<NDIM, double> > b_data = patch->getPatchData(b_idx);
        const Array<BoundaryBox<NDIM> >& type_1_cf_bdry =
            level_zero ? Array<BoundaryBox<NDIM> >() :
                         d_cf_boundary->getBoundaries(patch->getPatchNumber(), /* boundary type */ 1);
        const bool at_physical_bdry = pgeom->intersectsPhysicalBoundary();
        const bool at_cf_bdry = type_1_cf_bdry.size() > 0;
        if (at_physical_bdry || at_cf_bdry)
        {
            CellData<NDIM, double> b_adj_data(b_data->getBox(), b_data->getDepth(), b_data->getGhostCellWidth());
            b_adj_data.copy(*b_data);
            if (at_physical_bdry)
            {
                PoissonUtilities::adjustRHSAtPhysicalBoundary(
                    b_adj_data, patch, d_poisson_spec, d_bc_coefs, d_solution_time, d_homogeneous_bc);
            }
            if (at_cf_bdry)
            {
                PoissonUtilities::adjustRHSAtCoarseFineBoundary(
                    b_adj_data, *x_data, patch, d_poisson_spec, type_1_cf_bdry);
            }
            copyToHypre(d_rhs_vecs, b_adj_data, patch_box);
        }
        else
        {
            copyToHypre(d_rhs_vecs, *b_data, patch_box);
        }
    }

    for (unsigned int k = 0; k < d_depth; ++k)
    {
        // Assemble the hypre vectors.
        HYPRE_StructVectorAssemble(d_sol_vecs[k]);
        HYPRE_StructVectorAssemble(d_rhs_vecs[k]);

        // Solve the system.
        IBTK_TIMER_START(t_solve_system_hypre);

        if (d_solver_type == "PFMG")
        {
            HYPRE_StructPFMGSetMaxIter(d_solvers[k], d_max_iterations);
            HYPRE_StructPFMGSetTol(d_solvers[k], d_rel_residual_tol);
            if (d_initial_guess_nonzero)
            {
                HYPRE_StructPFMGSetNonZeroGuess(d_solvers[k]);
            }
            else
            {
                HYPRE_StructPFMGSetZeroGuess(d_solvers[k]);
            }
            HYPRE_StructPFMGSolve(d_solvers[k], d_matrices[k], d_rhs_vecs[k], d_sol_vecs[k]);
            HYPRE_StructPFMGGetNumIterations(d_solvers[k], &d_current_iterations);
            HYPRE_StructPFMGGetFinalRelativeResidualNorm(d_solvers[k], &d_current_residual_norm);
        }
        else if (d_solver_type == "SMG")
        {
            HYPRE_StructSMGSetMaxIter(d_solvers[k], d_max_iterations);
            HYPRE_StructSMGSetTol(d_solvers[k], d_rel_residual_tol);
            if (d_initial_guess_nonzero)
            {
                HYPRE_StructSMGSetNonZeroGuess(d_solvers[k]);
            }
            else
            {
                HYPRE_StructSMGSetZeroGuess(d_solvers[k]);
            }
            HYPRE_StructSMGSolve(d_solvers[k], d_matrices[k], d_rhs_vecs[k], d_sol_vecs[k]);
            HYPRE_StructSMGGetNumIterations(d_solvers[k], &d_current_iterations);
            HYPRE_StructSMGGetFinalRelativeResidualNorm(d_solvers[k], &d_current_residual_norm);
        }
        else if (d_solver_type == "PCG")
        {
            HYPRE_StructPCGSetMaxIter(d_solvers[k], d_max_iterations);
            HYPRE_StructPCGSetTol(d_solvers[k], d_rel_residual_tol);
            HYPRE_StructPCGSetAbsoluteTol(d_solvers[k], d_abs_residual_tol);
            HYPRE_StructPCGSolve(d_solvers[k], d_matrices[k], d_rhs_vecs[k], d_sol_vecs[k]);
            HYPRE_StructPCGGetNumIterations(d_solvers[k], &d_current_iterations);
            HYPRE_StructPCGGetFinalRelativeResidualNorm(d_solvers[k], &d_current_residual_norm);
        }
        else if (d_solver_type == "GMRES")
        {
            HYPRE_StructGMRESSetMaxIter(d_solvers[k], d_max_iterations);
            HYPRE_StructGMRESSetTol(d_solvers[k], d_rel_residual_tol);
            HYPRE_StructGMRESSetAbsoluteTol(d_solvers[k], d_abs_residual_tol);
            HYPRE_StructGMRESSolve(d_solvers[k], d_matrices[k], d_rhs_vecs[k], d_sol_vecs[k]);
            HYPRE_StructGMRESGetNumIterations(d_solvers[k], &d_current_iterations);
            HYPRE_StructGMRESGetFinalRelativeResidualNorm(d_solvers[k], &d_current_residual_norm);
        }
        else if (d_solver_type == "FlexGMRES")
        {
            HYPRE_StructFlexGMRESSetMaxIter(d_solvers[k], d_max_iterations);
            HYPRE_StructFlexGMRESSetTol(d_solvers[k], d_rel_residual_tol);
            HYPRE_StructFlexGMRESSetAbsoluteTol(d_solvers[k], d_abs_residual_tol);
            HYPRE_StructFlexGMRESSolve(d_solvers[k], d_matrices[k], d_rhs_vecs[k], d_sol_vecs[k]);
            HYPRE_StructFlexGMRESGetNumIterations(d_solvers[k], &d_current_iterations);
            HYPRE_StructFlexGMRESGetFinalRelativeResidualNorm(d_solvers[k], &d_current_residual_norm);
        }
        else if (d_solver_type == "LGMRES")
        {
            HYPRE_StructLGMRESSetMaxIter(d_solvers[k], d_max_iterations);
            HYPRE_StructLGMRESSetTol(d_solvers[k], d_rel_residual_tol);
            HYPRE_StructLGMRESSetAbsoluteTol(d_solvers[k], d_abs_residual_tol);
            HYPRE_StructLGMRESSolve(d_solvers[k], d_matrices[k], d_rhs_vecs[k], d_sol_vecs[k]);
            HYPRE_StructLGMRESGetNumIterations(d_solvers[k], &d_current_iterations);
            HYPRE_StructLGMRESGetFinalRelativeResidualNorm(d_solvers[k], &d_current_residual_norm);
        }
        else if (d_solver_type == "BiCGSTAB")
        {
            HYPRE_StructBiCGSTABSetMaxIter(d_solvers[k], d_max_iterations);
            HYPRE_StructBiCGSTABSetTol(d_solvers[k], d_rel_residual_tol);
            HYPRE_StructBiCGSTABSetAbsoluteTol(d_solvers[k], d_abs_residual_tol);
            HYPRE_StructBiCGSTABSolve(d_solvers[k], d_matrices[k], d_rhs_vecs[k], d_sol_vecs[k]);
            HYPRE_StructBiCGSTABGetNumIterations(d_solvers[k], &d_current_iterations);
            HYPRE_StructBiCGSTABGetFinalRelativeResidualNorm(d_solvers[k], &d_current_residual_norm);
        }
    }
    IBTK_TIMER_STOP(t_solve_system_hypre);

    // Pull the solution vector out of the hypre structures.
    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = d_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, double> > x_data = patch->getPatchData(x_idx);
        copyFromHypre(*x_data, d_sol_vecs, patch_box);
    }
    return (d_current_residual_norm <= d_rel_residual_tol || d_current_residual_norm <= d_abs_residual_tol);
} // solveSystem

void
CCPoissonHypreLevelSolver::copyToHypre(const std::vector<HYPRE_StructVector>& vectors,
                                       const CellData<NDIM, double>& src_data,
                                       const Box<NDIM>& box)
{
    Index<NDIM> lower = box.lower();
    Index<NDIM> upper = box.upper();
    if (src_data.getGhostBox() == box)
    {
        for (unsigned int k = 0; k < d_depth; ++k)
        {
            HYPRE_StructVectorSetBoxValues(vectors[k], lower, upper, const_cast<double*>(src_data.getPointer(k)));
        }
    }
    else
    {
        CellData<NDIM, double> hypre_data(box, 1, 0);
        for (unsigned int k = 0; k < d_depth; ++k)
        {
            hypre_data.copyDepth(0, src_data, k);
            HYPRE_StructVectorSetBoxValues(vectors[k], lower, upper, hypre_data.getPointer());
        }
    }
    return;
} // copyToHypre

void
CCPoissonHypreLevelSolver::copyFromHypre(CellData<NDIM, double>& dst_data,
                                         const std::vector<HYPRE_StructVector>& vectors,
                                         const Box<NDIM>& box)
{
    Index<NDIM> lower = box.lower();
    Index<NDIM> upper = box.upper();
    if (dst_data.getGhostBox() == box)
    {
        for (unsigned int k = 0; k < d_depth; ++k)
        {
            HYPRE_StructVectorGetBoxValues(vectors[k], lower, upper, dst_data.getPointer(k));
        }
    }
    else
    {
        CellData<NDIM, double> hypre_data(box, 1, 0);
        for (unsigned int k = 0; k < d_depth; ++k)
        {
            HYPRE_StructVectorGetBoxValues(vectors[k], lower, upper, hypre_data.getPointer());
            dst_data.copyDepth(k, hypre_data, 0);
        }
    }
    return;
} // copyFromHypre

void
CCPoissonHypreLevelSolver::destroyHypreSolver()
{
    for (unsigned int k = 0; k < d_depth; ++k)
    {
        // Destroy the solver.
        if (d_solver_type == "PFMG")
        {
            HYPRE_StructPFMGDestroy(d_solvers[k]);
        }
        else if (d_solver_type == "SMG")
        {
            HYPRE_StructSMGDestroy(d_solvers[k]);
        }
        else if (d_solver_type == "PCG")
        {
            HYPRE_StructPCGDestroy(d_solvers[k]);
        }
        else if (d_solver_type == "GMRES")
        {
            HYPRE_StructGMRESDestroy(d_solvers[k]);
        }
        else if (d_solver_type == "FlexGMRES")
        {
            HYPRE_StructFlexGMRESDestroy(d_solvers[k]);
        }
        else if (d_solver_type == "LGMRES")
        {
            HYPRE_StructLGMRESDestroy(d_solvers[k]);
        }
        else if (d_solver_type == "BiCGSTAB")
        {
            HYPRE_StructBiCGSTABDestroy(d_solvers[k]);
        }

        // When using a Krylov method, destroy the preconditioner.
        if (d_solver_type == "PCG" || d_solver_type == "GMRES" || d_solver_type == "FlexGMRES" ||
            d_solver_type == "LGMRES" || d_solver_type == "BiCGSTAB")
        {
            if (d_precond_type == "PFMG")
            {
                HYPRE_StructPFMGDestroy(d_preconds[k]);
            }
            else if (d_precond_type == "SMG")
            {
                HYPRE_StructSMGDestroy(d_preconds[k]);
            }
            else if (d_precond_type == "Jacobi")
            {
                HYPRE_StructJacobiDestroy(d_preconds[k]);
            }
        }

        // Set the solver and preconditioner pointers to NULL.
        d_solvers[k] = NULL;
        d_preconds[k] = NULL;
    }
    return;
} // destroyHypreSolver

void
CCPoissonHypreLevelSolver::deallocateHypreData()
{
    if (d_grid) HYPRE_StructGridDestroy(d_grid);
    if (d_stencil) HYPRE_StructStencilDestroy(d_stencil);
    d_grid = NULL;
    d_stencil = NULL;
    for (unsigned int k = 0; k < d_depth; ++k)
    {
        if (d_matrices[k]) HYPRE_StructMatrixDestroy(d_matrices[k]);
        if (d_sol_vecs[k]) HYPRE_StructVectorDestroy(d_sol_vecs[k]);
        if (d_rhs_vecs[k]) HYPRE_StructVectorDestroy(d_rhs_vecs[k]);
        d_matrices[k] = NULL;
        d_sol_vecs[k] = NULL;
        d_rhs_vecs[k] = NULL;
    }
    return;
} // deallocateHypreData

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
