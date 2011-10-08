// Filename: CCPoissonHypreLevelSolver.C
// Created on 30 May 2005 by Boyce Griffith
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

#include "CCPoissonHypreLevelSolver.h"

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
#include <ibtk/ExtendedRobinBcCoefStrategy.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <ArrayDataBasicOps.h>
#include <Index.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <OutersideData.h>
#include <SideDataFactory.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

// C++ STDLIB INCLUDES
#include <map>
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

CCPoissonHypreLevelSolver::CCPoissonHypreLevelSolver(
    const std::string& object_name,
    Pointer<Database> input_db)
    : d_object_name(object_name),
      d_is_initialized(false),
      d_hierarchy(),
      d_level_num(-1),
      d_poisson_spec(d_object_name+"::Poisson specs"),
      d_grid_aligned_anisotropy(true),
      d_default_bc_coef(new LocationIndexRobinBcCoefs<NDIM>(
                            d_object_name+"::default_bc_coef", Pointer<Database>(NULL))),
      d_bc_coef(NULL),
      d_homogeneous_bc(true),
      d_apply_time(0.0),
      d_depth(0),
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
        TBOX_ERROR(d_object_name << "::CCPoissonHypreLevelSolver()"
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

    // Initialize the Poisson specifications.
    d_poisson_spec.setCZero();
    d_poisson_spec.setDConstant(-1.0);
    d_grid_aligned_anisotropy = true;

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

    // Setup Timers.
    IBTK_DO_ONCE(
        t_solve_system            = TimerManager::getManager()->getTimer("IBTK::CCPoissonHypreLevelSolver::solveSystem()");
        t_solve_system_hypre      = TimerManager::getManager()->getTimer("IBTK::CCPoissonHypreLevelSolver::solveSystem()[hypre]");
        t_initialize_solver_state = TimerManager::getManager()->getTimer("IBTK::CCPoissonHypreLevelSolver::initializeSolverState()");
        t_deallocate_solver_state = TimerManager::getManager()->getTimer("IBTK::CCPoissonHypreLevelSolver::deallocateSolverState()");
                 );
    return;
}// CCPoissonHypreLevelSolver

CCPoissonHypreLevelSolver::~CCPoissonHypreLevelSolver()
{
    if (d_is_initialized) deallocateSolverState();
    delete d_default_bc_coef;
    return;
}// ~CCPoissonHypreLevelSolver

void
CCPoissonHypreLevelSolver::setPoissonSpecifications(
    const PoissonSpecifications& poisson_spec)
{
    d_poisson_spec = poisson_spec;
    if (d_poisson_spec.dIsVariable())
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<SideDataFactory<NDIM,double> > pdat_factory =
            var_db->getPatchDescriptor()->getPatchDataFactory(
                d_poisson_spec.getDPatchDataId());
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(!pdat_factory.isNull());
#endif
        d_grid_aligned_anisotropy = pdat_factory->getDefaultDepth() == 1;
    }
    else
    {
        d_grid_aligned_anisotropy = true;
    }
    return;
}// setPoissonSpecifications

void
CCPoissonHypreLevelSolver::setPhysicalBcCoef(
    RobinBcCoefStrategy<NDIM>* bc_coef)
{
    if (bc_coef != NULL)
    {
        d_bc_coef = bc_coef;
    }
    else
    {
        d_bc_coef = d_default_bc_coef;
        setHomogeneousBc(true);
    }
    return;
}// setPhysicalBcCoef

void
CCPoissonHypreLevelSolver::setHomogeneousBc(
    bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

void
CCPoissonHypreLevelSolver::setTime(
    const double time)
{
    d_apply_time = time;
    return;
}// setTime

void
CCPoissonHypreLevelSolver::setDataDepth(
    const int depth)
{
    d_depth = depth;
    return;
}// setDataDepth

bool
CCPoissonHypreLevelSolver::solveSystem(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& b)
{
    IBTK_TIMER_START(t_solve_system);

    if (d_enable_logging) plog << d_object_name << "::solveSystem():" << std::endl;

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x,b);

    // Solve the system using the hypre solver.
    static const int comp = 0;
    const int x_idx = x.getComponentDescriptorIndex(comp);
    const int b_idx = b.getComponentDescriptorIndex(comp);
    const bool converged = solveSystem(x_idx, b_idx);

    // Log solver info.
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
CCPoissonHypreLevelSolver::initializeSolverState(
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
                   << "  coarsest_ln != finest_ln in CCPoissonHypreLevelSolver" << std::endl);
    }
#else
    NULL_USE(b);
#endif
    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized) deallocateSolverState();

    // Get the hierarchy information.
    d_hierarchy = x.getPatchHierarchy();
    d_level_num = x.getCoarsestLevelNumber();

    // Allocate and initialize the hypre data structures.
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
}// initializeSolverState

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
}// deallocateSolverState

void
CCPoissonHypreLevelSolver::enableLogging(
    bool enabled)
{
    d_enable_logging = enabled;
    return;
}// enableLogging

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CCPoissonHypreLevelSolver::allocateHypreData()
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

    HYPRE_StructGridCreate(communicator, NDIM, &d_grid);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        const Box<NDIM>& patch_box = level->getPatch(p())->getBox();
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
    }
    else
    {
        static const int stencil_sz = (NDIM == 2) ? 9 : 19;
        int stencil_offsets[stencil_sz][NDIM];
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
                        stencil_offsets[stencil_index][0] = x_offset;
                        stencil_offsets[stencil_index][1] = y_offset;
#if (NDIM == 3)
                        stencil_offsets[stencil_index][2] = z_offset;
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
            HYPRE_StructStencilSetElement(d_stencil, s, stencil_offsets[s]);
        }
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
CCPoissonHypreLevelSolver::setMatrixCoefficients_aligned()
{
    ArrayDataBasicOps<NDIM,double> array_ops;
    static const IntVector<NDIM> no_ghosts = 0;
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_num);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();

        // Compute all off-diagonal matrix coefficients for all cell sides,
        // including those that touch the physical boundary, but temporarily
        // ignoring physical boundary conditions.
        Pointer<SideData<NDIM,double> > D_data;
        if (!d_poisson_spec.dIsConstant())
        {
            D_data = patch->getPatchData(d_poisson_spec.getDPatchDataId());
            if (D_data.isNull())
            {
                TBOX_ERROR(d_object_name << "::setMatrixCoefficients_aligned()\n"
                           << "  to solve C u + div D grad u = f with non-constant D,\n"
                           << "  D must be side-centered double precision data" << std::endl);
            }
        }
        else
        {
            D_data = new SideData<NDIM,double>(patch_box, 1, no_ghosts);
            D_data->fill(d_poisson_spec.getDConstant());
        }

        SideData<NDIM,double> off_diagonal(patch_box, 1, no_ghosts);
        off_diagonal.fill(0.0);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
            array_ops.scale(off_diagonal.getArrayData(axis),
                            1.0/(dx[axis]*dx[axis]),
                            D_data->getArrayData(axis),
                            side_box);
        }

        // Compute the diagonal matrix coefficients.
        Pointer<CellData<NDIM,double> > C_data;
        if (!d_poisson_spec.cIsZero() && !d_poisson_spec.cIsConstant())
        {
            C_data = patch->getPatchData(d_poisson_spec.getCPatchDataId());
            if (C_data.isNull())
            {
                TBOX_ERROR(d_object_name << "::setMatrixCoefficients_aligned()\n"
                           << "  to solve (C u + div D grad u) = f with non-constant C,\n"
                           << "  C must be cell-centered double precision data" << std::endl);
            }
        }
        else
        {
            C_data = new CellData<NDIM,double>(patch_box, 1, no_ghosts);
            if (d_poisson_spec.cIsZero()) C_data->fill(0.0);
            else C_data->fill(d_poisson_spec.getCConstant());
        }

        CellData<NDIM,double> diagonal(patch_box, 1, no_ghosts);
        diagonal.copy(*C_data);

        for (Box<NDIM>::Iterator b(patch_box); b; b++)
        {
            Index<NDIM> i = b();
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const SideIndex<NDIM> ilower(i, axis, SideIndex<NDIM>::Lower);
                diagonal(i) -= off_diagonal(ilower);
                const SideIndex<NDIM> iupper(i, axis, SideIndex<NDIM>::Upper);
                diagonal(i) -= off_diagonal(iupper);
            }
        }

        // Modify the diagonal and off-diagonal entries to account for
        // homogeneous boundary conditions.
        //
        // Here, we follow the same linear extrapolation approach implemented in
        // class CartesianRobinBcHelper.  Namely, with u_i
        // denoting the interior cell, u_o denoting the ghost cell, and u_b and
        // u_n denoting the value and normal derivative of u at the boundary,
        //
        //     u_b = (u_i + u_o)/2   and   u_n = (u_o - u_i)/h
        //
        // Now, if
        //
        //     a*u_b + b*u_n = 0
        //
        // then
        //
        //     u_o = -((a*h - 2*b)/(a*h + 2*b))*u_i
        const Array<BoundaryBox<NDIM> > physical_codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
        const int n_physical_codim1_boxes = physical_codim1_boxes.size();
        for (int n = 0; n < n_physical_codim1_boxes; ++n)
        {
            const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
            const BoundaryBox<NDIM> trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
            const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            ArrayData<NDIM,double> acoef_data(bc_coef_box, 1);
            ArrayData<NDIM,double> bcoef_data(bc_coef_box, 1);

            Pointer<ArrayData<NDIM,double> > acoef_data_ptr(&acoef_data, false);
            Pointer<ArrayData<NDIM,double> > bcoef_data_ptr(&bcoef_data, false);
            Pointer<ArrayData<NDIM,double> > gcoef_data_ptr(NULL);

            ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coef);
            if (extended_bc_coef != NULL) extended_bc_coef->setHomogeneousBc(true);
            d_bc_coef->setBcCoefs(
                acoef_data_ptr, bcoef_data_ptr, gcoef_data_ptr, NULL,
                *patch, trimmed_bdry_box, d_apply_time);
            if (extended_bc_coef != NULL) extended_bc_coef->setHomogeneousBc(d_homogeneous_bc);

            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis =  location_index / 2;
            const bool bdry_lower_side = (location_index % 2) == 0;
            const bool bdry_upper_side = (location_index % 2) != 0;

            for (Box<NDIM>::Iterator b(bc_coef_box); b; b++)
            {
                const Index<NDIM>& i_s_bdry = b();
                const double& a = acoef_data(i_s_bdry,0);
                const double& b = bcoef_data(i_s_bdry,0);
                const double& h = dx[bdry_normal_axis];

                // i_s_bdry: side index located on physical boundary
                //
                // i_c_intr: cell index located adjacent to physical boundary in
                // the patch interior
                Index<NDIM> i_c_intr = i_s_bdry;
                if (bdry_upper_side)
                {
                    i_c_intr(bdry_normal_axis) -= 1;
                }

                if (bdry_lower_side)
                {
                    const SideIndex<NDIM> ilower(i_c_intr, bdry_normal_axis, SideIndex<NDIM>::Lower);
                    diagonal(i_c_intr) += off_diagonal(ilower)*(-(a*h-2.0*b)/(a*h+2.0*b));
                    off_diagonal(ilower) = 0.0;
                }

                if (bdry_upper_side)
                {
                    const SideIndex<NDIM> iupper(i_c_intr, bdry_normal_axis, SideIndex<NDIM>::Upper);
                    diagonal(i_c_intr) += off_diagonal(iupper)*(-(a*h-2.0*b)/(a*h+2.0*b));
                    off_diagonal(iupper) = 0.0;
                }
            }
        }

        // Copy matrix entries to the hypre matrix structure.
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

            const SideIndex<NDIM> ixlower(i, SideIndex<NDIM>::X, SideIndex<NDIM>::Lower);
            mat_vals[0] = off_diagonal(ixlower);

            const SideIndex<NDIM> iylower(i, SideIndex<NDIM>::Y, SideIndex<NDIM>::Lower);
            mat_vals[1] = off_diagonal(iylower);
#if (NDIM == 3)
            SideIndex<NDIM> izlower(i, SideIndex<NDIM>::Z, SideIndex<NDIM>::Lower);
            mat_vals[2] = off_diagonal(izlower);
#endif
            const SideIndex<NDIM> ixupper(i, SideIndex<NDIM>::X, SideIndex<NDIM>::Upper);
            mat_vals[NDIM+0] = off_diagonal(ixupper);

            const SideIndex<NDIM> iyupper(i, SideIndex<NDIM>::Y, SideIndex<NDIM>::Upper);
            mat_vals[NDIM+1] = off_diagonal(iyupper);
#if (NDIM == 3)
            SideIndex<NDIM> izupper(i, SideIndex<NDIM>::Z, SideIndex<NDIM>::Upper);
            mat_vals[NDIM+2] = off_diagonal(izupper);
#endif
            mat_vals[2*NDIM] = diagonal(i);

            HYPRE_StructMatrixSetValues(d_matrix, i, stencil_sz, stencil_indices, &mat_vals[0]);
        }
    }

    // Assemble the hypre matrix.
    HYPRE_StructMatrixAssemble(d_matrix);
    return;
}// setMatrixCoefficients_aligned

void
CCPoissonHypreLevelSolver::setMatrixCoefficients_nonaligned()
{
    static const IntVector<NDIM> no_ghosts = 0;
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_num);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();

        // Compute all matrix coefficients.
        //
        // NOTE: Here we assume that no flux boundary conditions are imposed at
        // the physical domain.
        Pointer<CellData<NDIM,double> > C_data;
        if (!d_poisson_spec.cIsZero() && !d_poisson_spec.cIsConstant())
        {
            C_data = patch->getPatchData(d_poisson_spec.getCPatchDataId());
            if (C_data.isNull())
            {
                TBOX_ERROR(d_object_name << "::setMatrixCoefficients_nonaligned()\n"
                           << "  to solve (C u + div D grad u) = f with non-constant C,\n"
                           << "  C must be cell-centered double precision data" << std::endl);
            }
        }
        else
        {
            C_data = new CellData<NDIM,double>(patch_box, 1, no_ghosts);
            if (d_poisson_spec.cIsZero()) C_data->fill(0.0);
            else C_data->fill(d_poisson_spec.getCConstant());
        }

        Pointer<SideData<NDIM,double> > D_data;
        if (!d_poisson_spec.dIsConstant())
        {
            D_data = patch->getPatchData(d_poisson_spec.getDPatchDataId());
        }

        if (D_data.isNull())
        {
            TBOX_ERROR(d_object_name << "::setMatrixCoefficients_nonaligned()\n"
                       << "  to solve C u + div D grad u = f with non-constant D,\n"
                       << "  D must be side-centered double precision data" << std::endl);
        }

        // Setup the finite difference stencil.
        static const int stencil_sz = (NDIM == 2 ? 9 : 19);
        int stencil_indices[stencil_sz];
        for (int i = 0; i < stencil_sz; ++i)
        {
            stencil_indices[i] = i;
        }

        std::map<Index<NDIM>,int,IndexComp> stencil_index_map;
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
                        const Index<NDIM> i(x_offset,y_offset);
#endif
#if (NDIM == 3)
                        const Index<NDIM> i(x_offset,y_offset,z_offset);
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

            std::vector<double> mat_vals(stencil_sz,0.0);
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
                    const double& D_lower = (*D_data)(ilower,axis);
                    mat_vals[stencil_lower ] += D_lower/(h*h);
                    mat_vals[stencil_center] -= D_lower/(h*h);
                }
                {
                    Index<NDIM> i_stencil_upper(0);
                    ++i_stencil_upper[axis];
                    const int stencil_upper = stencil_index_map[i_stencil_upper];

                    const SideIndex<NDIM> iupper(i, axis, SideIndex<NDIM>::Upper);
                    const double& D_upper = (*D_data)(iupper,axis);
                    mat_vals[stencil_upper ] += D_upper/(h*h);
                    mat_vals[stencil_center] -= D_upper/(h*h);
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
                        patch->getPatchGeometry()->getTouchesRegularBoundary(norm_axis,lower) &&
                        i(norm_axis) == patch_box.lower()(norm_axis);

                    static const int upper = 1;
                    const bool skip_upper_side_fluxes =
                        patch->getPatchGeometry()->getTouchesRegularBoundary(norm_axis,upper) &&
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
                                i_stencil[ norm_axis] +=  norm_shift;
                                i_stencil[trans_axis] += trans_shift;
                                const int stencil_index = stencil_index_map[i_stencil];
                                if (trans_shift == 1)
                                {
                                    mat_vals[stencil_index] -= 0.25*(*D_data)(ilower,norm_axis)/(norm_h*trans_h);
                                }
                                else
                                {
                                    mat_vals[stencil_index] += 0.25*(*D_data)(ilower,norm_axis)/(norm_h*trans_h);
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
                                i_stencil[ norm_axis] +=  norm_shift;
                                i_stencil[trans_axis] += trans_shift;
                                const int stencil_index = stencil_index_map[i_stencil];
                                if (trans_shift == 1)
                                {
                                    mat_vals[stencil_index] += 0.25*(*D_data)(iupper,norm_axis)/(norm_h*trans_h);
                                }
                                else
                                {
                                    mat_vals[stencil_index] -= 0.25*(*D_data)(iupper,norm_axis)/(norm_h*trans_h);
                                }
                            }
                        }
                    }
                }
            }

            HYPRE_StructMatrixSetValues(d_matrix, i, stencil_sz, stencil_indices, &mat_vals[0]);
        }
    }

    // Assemble the hypre matrix.
    HYPRE_StructMatrixAssemble(d_matrix);
    return;
}// setMatrixCoefficients_nonaligned

void
CCPoissonHypreLevelSolver::setupHypreSolver()
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
CCPoissonHypreLevelSolver::solveSystem(
    const int x_idx,
    const int b_idx)
{
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_num);

    // Modify right-hand-side data to account for boundary conditions and copy
    // solution and right-hand-side data to hypre structures.
    const IntVector<NDIM> ghosts = 1;
    const IntVector<NDIM> no_ghosts = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();

        // Copy the solution data into the hypre vector, including ghost cell
        // values
        const Box<NDIM> x_ghost_box = Box<NDIM>::grow(patch_box, 1);
        Pointer<CellData<NDIM,double> > x_data = patch->getPatchData(x_idx);
        copyToHypre(d_sol_vec, x_data, x_ghost_box);

        // Modify the right-hand-side data to account for any inhomogeneous
        // boundary conditions and copy the right-hand-side into the hypre
        // vector.
        Pointer<CellData<NDIM,double> > b_data = patch->getPatchData(b_idx);
        if (!d_homogeneous_bc && pgeom->intersectsPhysicalBoundary())
        {
            b_data = new CellData<NDIM,double>(
                b_data->getBox(), b_data->getDepth(),
                b_data->getGhostCellWidth());
            b_data->copy(*(patch->getPatchData(b_idx)));

            OutersideData<NDIM,double> D_os_data(patch_box, 1);
            Pointer<OutersideData<NDIM,double> > D_os_data_ptr(&D_os_data,false);
            if (!d_poisson_spec.dIsConstant())
            {
                D_os_data.copy(*patch->getPatchData(d_poisson_spec.getDPatchDataId()));
            }
            else
            {
                D_os_data.fillAll(d_poisson_spec.getDConstant());
            }

            const Array<BoundaryBox<NDIM> > physical_codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);

            if (d_grid_aligned_anisotropy)
            {
                adjustBoundaryRhsEntries_aligned(
                    b_data, D_os_data_ptr, patch, physical_codim1_boxes, dx);
            }
            else
            {
                adjustBoundaryRhsEntries_nonaligned(
                    b_data, D_os_data_ptr, patch, physical_codim1_boxes, dx);
            }
        }
        copyToHypre(d_rhs_vec, b_data, patch_box);
    }

    // Assemble the hypre vectors.
    HYPRE_StructVectorAssemble(d_sol_vec);
    HYPRE_StructVectorAssemble(d_rhs_vec);

    // Solve the system.
    IBTK_TIMER_START(t_solve_system_hypre);

    if (d_solver_type == "PFMG")
    {
        HYPRE_StructPFMGSetMaxIter(d_solver, d_max_iterations);
        HYPRE_StructPFMGSetTol(d_solver, d_rel_residual_tol);
        if (d_initial_guess_nonzero)
        {
            HYPRE_StructPFMGSetNonZeroGuess(d_solver);
        }
        else
        {
            HYPRE_StructPFMGSetZeroGuess(d_solver);
        }
        HYPRE_StructPFMGSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_StructPFMGGetNumIterations(d_solver, &d_current_its);
        HYPRE_StructPFMGGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "SMG")
    {
        HYPRE_StructSMGSetMaxIter(d_solver, d_max_iterations);
        HYPRE_StructSMGSetTol(d_solver, d_rel_residual_tol);
        if (d_initial_guess_nonzero)
        {
            HYPRE_StructSMGSetNonZeroGuess(d_solver);
        }
        else
        {
            HYPRE_StructSMGSetZeroGuess(d_solver);
        }
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

    IBTK_TIMER_STOP(t_solve_system_hypre);

    // Pull the solution vector out of the hypre structures.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM,double> > x_data = patch->getPatchData(x_idx);
        copyFromHypre(x_data, d_sol_vec, patch_box);
    }
    return (d_current_residual_norm <= d_rel_residual_tol || d_current_residual_norm <= d_abs_residual_tol);
}// solveSystem

void
CCPoissonHypreLevelSolver::copyToHypre(
    HYPRE_StructVector vector,
    const Pointer<CellData<NDIM,double> > src_data,
    const Box<NDIM>& box)
{
    Index<NDIM> lower = box.lower();
    Index<NDIM> upper = box.upper();
    if (src_data->getGhostBox() == box)
    {
        HYPRE_StructVectorSetBoxValues(vector,lower,upper,src_data->getPointer(d_depth));
    }
    else
    {
        CellData<NDIM,double> hypre_data(box,1,0);
        hypre_data.copyDepth(0,*src_data,d_depth);
        HYPRE_StructVectorSetBoxValues(vector,lower,upper,hypre_data.getPointer());
    }
    return;
}// copyToHypre

void
CCPoissonHypreLevelSolver::copyFromHypre(
    Pointer<CellData<NDIM,double> > dst_data,
    HYPRE_StructVector vector,
    const Box<NDIM>& box)
{
    Index<NDIM> lower = box.lower();
    Index<NDIM> upper = box.upper();
    if (dst_data->getGhostBox() == box)
    {
        HYPRE_StructVectorGetBoxValues(vector,lower,upper,dst_data->getPointer(d_depth));
    }
    else
    {
        CellData<NDIM,double> hypre_data(box,1,0);
        HYPRE_StructVectorGetBoxValues(vector,lower,upper,hypre_data.getPointer());
        dst_data->copyDepth(d_depth,hypre_data,0);
    }
    return;
}// copyFromHypre

void
CCPoissonHypreLevelSolver::destroyHypreSolver()
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
CCPoissonHypreLevelSolver::deallocateHypreData()
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

void
CCPoissonHypreLevelSolver::adjustBoundaryRhsEntries_aligned(
    Pointer<CellData<NDIM,double> > rhs_data,
    const Pointer<OutersideData<NDIM,double> > D_data,
    const Pointer<Patch<NDIM> > patch,
    const Array<BoundaryBox<NDIM> >& codim1_boxes,
    const double* const dx)
{
    // Modify the rhs entries to account for inhomogeneous boundary conditions.
    //
    // Here, we follow the same linear extrapolation approach implemented in
    // class CartesianRobinBcHelper.  Namely, with u_i denoting
    // the interior cell, u_o denoting the ghost cell, and u_b and u_n denoting
    // the value and normal derivative of u at the boundary,
    //
    //     u_b = (u_i + u_o)/2   and   u_n = (u_o - u_i)/h
    //
    // Now, if
    //
    //     a*u_b + b*u_n = g
    //
    // then, with u_i = 0,
    //
    //     u_o = 2*h*g/(2*b + a*h)
    //
    // so that the boundary flux is
    //
    //     (u_i - u_o)/h = -2*g/(2*b + h*a)
    //
    // In this loop, we modify the rhs entries appropriately.
    const int n_bdry_boxes = codim1_boxes.size();
    for (int n = 0; n < n_bdry_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = codim1_boxes[n];
        const BoundaryBox<NDIM> trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
        const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

        ArrayData<NDIM,double> acoef_data(bc_coef_box, 1);
        ArrayData<NDIM,double> bcoef_data(bc_coef_box, 1);
        ArrayData<NDIM,double> gcoef_data(bc_coef_box, 1);

        Pointer<ArrayData<NDIM,double> > acoef_data_ptr(&acoef_data, false);
        Pointer<ArrayData<NDIM,double> > bcoef_data_ptr(&bcoef_data, false);
        Pointer<ArrayData<NDIM,double> > gcoef_data_ptr(&gcoef_data, false);

        ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coef);
        if (extended_bc_coef != NULL) extended_bc_coef->setHomogeneousBc(d_homogeneous_bc);
        d_bc_coef->setBcCoefs(
            acoef_data_ptr, bcoef_data_ptr, gcoef_data_ptr, NULL,
            *patch, trimmed_bdry_box, d_apply_time);

        const unsigned int location_index = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis =  location_index / 2;
        const bool bdry_upper_side = (location_index % 2) != 0;
        const int bdry_side = (bdry_upper_side ? 1 : 0);

        // i_s_bdry: side index located on physical boundary
        //
        // i_c_intr: cell index located adjacent to physical boundary in the
        // patch interior
        for (Box<NDIM>::Iterator b(bc_coef_box); b; b++)
        {
            const Index<NDIM>& i_s_bdry = b();
            const double& a = acoef_data(i_s_bdry,0);
            const double& b = bcoef_data(i_s_bdry,0);
            const double& g = gcoef_data(i_s_bdry,0);
            const double& h = dx[bdry_normal_axis];

            Index<NDIM> i_c_intr = i_s_bdry;
            if (bdry_upper_side)
            {
                i_c_intr(bdry_normal_axis) -= 1;
            }
            (*rhs_data)(i_c_intr) += ((D_data->getArrayData(bdry_normal_axis,bdry_side))(i_s_bdry,0)/h)*(-2.0*g)/(2.0*b+h*a);
        }
    }
    return;
}// adjustBoundaryRhsEntries_aligned

void
CCPoissonHypreLevelSolver::adjustBoundaryRhsEntries_nonaligned(
    Pointer<CellData<NDIM,double> > rhs_data,
    const Pointer<OutersideData<NDIM,double> > D_data,
    const Pointer<Patch<NDIM> > patch,
    const Array<BoundaryBox<NDIM> >& codim1_boxes,
    const double* const dx)
{
    // Modify the rhs entries to account for inhomogeneous boundary conditions.
    //
    // Here, we follow the same linear extrapolation approach implemented in
    // class CartesianRobinBcHelper.  Namely, with u_i denoting
    // the interior cell, u_o denoting the ghost cell, and u_b and u_n denoting
    // the value and normal derivative of u at the boundary,
    //
    //     u_b = (u_i + u_o)/2   and   u_n = (u_o - u_i)/h
    //
    // Now, if
    //
    //     a*u_b + b*u_n = g
    //
    // then, with u_i = 0,
    //
    //     u_o = 2*h*g/(2*b + a*h)
    //
    // so that the boundary flux is
    //
    //     (u_i - u_o)/h = -2*g/(2*b + h*a)
    //
    // In this loop, we modify the rhs entries appropriately.
    const int n_bdry_boxes = codim1_boxes.size();
    for (int n = 0; n < n_bdry_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = codim1_boxes[n];
        const BoundaryBox<NDIM> trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
        const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

        ArrayData<NDIM,double> acoef_data(bc_coef_box, 1);
        ArrayData<NDIM,double> bcoef_data(bc_coef_box, 1);
        ArrayData<NDIM,double> gcoef_data(bc_coef_box, 1);

        Pointer<ArrayData<NDIM,double> > acoef_data_ptr(&acoef_data, false);
        Pointer<ArrayData<NDIM,double> > bcoef_data_ptr(&bcoef_data, false);
        Pointer<ArrayData<NDIM,double> > gcoef_data_ptr(&gcoef_data, false);

        ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coef);
        if (extended_bc_coef != NULL) extended_bc_coef->setHomogeneousBc(d_homogeneous_bc);
        d_bc_coef->setBcCoefs(
            acoef_data_ptr, bcoef_data_ptr, gcoef_data_ptr, NULL,
            *patch, trimmed_bdry_box, d_apply_time);

        const unsigned int location_index = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis =  location_index / 2;
        const bool bdry_upper_side = (location_index % 2) != 0;
        const int bdry_side = (bdry_upper_side ? 1 : 0);

        // i_s_bdry: side index located on physical boundary
        //
        // i_c_intr: cell index located adjacent to physical boundary in the
        // patch interior
        for (Box<NDIM>::Iterator b(bc_coef_box); b; b++)
        {
            const Index<NDIM>& i_s_bdry = b();
            const double& a = acoef_data(i_s_bdry,0);
            const double& b = bcoef_data(i_s_bdry,0);
            const double& g = gcoef_data(i_s_bdry,0);
            const double& h = dx[bdry_normal_axis];

            Index<NDIM> i_c_intr = i_s_bdry;
            if (bdry_upper_side)
            {
                i_c_intr(bdry_normal_axis) -= 1;
            }
            (*rhs_data)(i_c_intr) += ((D_data->getArrayData(bdry_normal_axis,bdry_side))(i_s_bdry,bdry_normal_axis)/h)*(-2.0*g)/(2.0*b+h*a);
        }
    }
    return;
}// adjustBoundaryRhsEntries_nonaligned

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
