// Filename: SCPoissonHypreLevelSolver.C
// Created on 17 Sep 2008 by Boyce Griffith
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

#include "SCPoissonHypreLevelSolver.h"

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
#include <iterator>

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

inline Box<NDIM>
compute_tangential_extension(
    const Box<NDIM>& box,
    const int data_axis)
{
    Box<NDIM> extended_box = box;
    extended_box.upper()(data_axis) += 1;
    return extended_box;
}// compute_tangential_extension
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

SCPoissonHypreLevelSolver::SCPoissonHypreLevelSolver(
    const std::string& object_name,
    Pointer<Database> input_db)
    : d_object_name(object_name),
      d_is_initialized(false),
      d_hierarchy(),
      d_level_num(-1),
      d_poisson_spec(d_object_name+"::Poisson specs"),
      d_constant_coefficients(true),
      d_default_bc_coef(new LocationIndexRobinBcCoefs<NDIM>(
                            d_object_name+"::default_bc_coef", Pointer<Database>(NULL))),
      d_bc_coefs(),
      d_homogeneous_bc(true),
      d_apply_time(0.0),
      d_grid(NULL),
      d_stencil(),
      d_graph(NULL),
      d_matrix(NULL),
      d_rhs_vec(NULL),
      d_sol_vec(NULL),
      d_solver(NULL),
      d_precond(NULL),
      d_solver_type("Split"),
      d_precond_type("none"),
      d_split_solver_type("PFMG"),
      d_max_iterations(10),
      d_abs_residual_tol(0.0),
      d_rel_residual_tol(1.0e-6),
      d_initial_guess_nonzero(false),
      d_rel_change(0),
      d_num_pre_relax_steps(1),
      d_num_post_relax_steps(1),
      d_relax_type(RELAX_TYPE_WEIGHTED_JACOBI),
      d_skip_relax(1),
      d_two_norm(1),
      d_current_its(-1),
      d_current_residual_norm(-1.0),
      d_enable_logging(false)
{
    if (NDIM == 1 || NDIM > 3)
    {
        TBOX_ERROR(d_object_name << "::SCPoissonHypreLevelSolver()"
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

        if (d_solver_type == "SysPFMG" || d_precond_type == "SysPFMG")
        {
            d_num_pre_relax_steps = input_db->getIntegerWithDefault("num_pre_relax_steps", d_num_pre_relax_steps);
            d_num_post_relax_steps = input_db->getIntegerWithDefault("num_post_relax_steps", d_num_post_relax_steps);
            d_relax_type = input_db->getIntegerWithDefault("relax_type", d_relax_type);
            d_skip_relax = input_db->getIntegerWithDefault("skip_relax", d_skip_relax);
        }

        if (d_solver_type == "Split" || d_precond_type == "Split")
        {
            d_split_solver_type = input_db->getStringWithDefault("split_solver_type", d_split_solver_type);
        }

        if (d_solver_type == "PCG")
        {
            d_two_norm = input_db->getIntegerWithDefault("two_norm", d_two_norm);
        }
    }

    // Initialize the Poisson specifications.
    d_poisson_spec.setCZero();
    d_poisson_spec.setDConstant(-1.0);
    d_constant_coefficients = true;

    // Setup a default boundary condition object that specifies homogeneous
    // Dirichlet boundary conditions.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_default_bc_coef->setBoundaryValue(2*d  ,0.0);
        d_default_bc_coef->setBoundaryValue(2*d+1,0.0);
    }

    // Initialize the boundary conditions objects.
    setHomogeneousBc(d_homogeneous_bc);
    setPhysicalBcCoefs(blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>(static_cast<RobinBcCoefStrategy<NDIM>*>(NULL)));

    // Setup Timers.
    IBTK_DO_ONCE(
        t_solve_system            = TimerManager::getManager()->getTimer("IBTK::SCPoissonHypreLevelSolver::solveSystem()");
        t_solve_system_hypre      = TimerManager::getManager()->getTimer("IBTK::SCPoissonHypreLevelSolver::solveSystem()[hypre]");
        t_initialize_solver_state = TimerManager::getManager()->getTimer("IBTK::SCPoissonHypreLevelSolver::initializeSolverState()");
        t_deallocate_solver_state = TimerManager::getManager()->getTimer("IBTK::SCPoissonHypreLevelSolver::deallocateSolverState()");
                 );
    return;
}// SCPoissonHypreLevelSolver

SCPoissonHypreLevelSolver::~SCPoissonHypreLevelSolver()
{
    if (d_is_initialized) deallocateSolverState();
    delete d_default_bc_coef;
    return;
}// ~SCPoissonHypreLevelSolver

void
SCPoissonHypreLevelSolver::setPoissonSpecifications(
    const PoissonSpecifications& poisson_spec)
{
    d_poisson_spec = poisson_spec;
    if (d_poisson_spec.cIsVariable() || d_poisson_spec.dIsVariable())
    {
        d_constant_coefficients = false;
    }
    else
    {
        d_constant_coefficients = true;
    }
    return;
}// setPoissonSpecifications

void
SCPoissonHypreLevelSolver::setPhysicalBcCoefs(
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& bc_coefs)
{
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (bc_coefs[d] != NULL)
        {
            d_bc_coefs[d] = bc_coefs[d];
        }
        else
        {
            d_bc_coefs[d] = d_default_bc_coef;
        }
    }
    return;
}// setPhysicalBcCoefs

void
SCPoissonHypreLevelSolver::setHomogeneousBc(
    bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

void
SCPoissonHypreLevelSolver::setTime(
    const double time)
{
    d_apply_time = time;
    return;
}// setTime

bool
SCPoissonHypreLevelSolver::solveSystem(
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
SCPoissonHypreLevelSolver::initializeSolverState(
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
                   << "  coarsest_ln != finest_ln in SCPoissonHypreLevelSolver" << std::endl);
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
    if (d_constant_coefficients)
    {
        setMatrixCoefficients_constant_coefficients();
    }
    else
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                   << "  only constant-coefficient problems are currently supported" << std::endl);
    }
    setupHypreSolver();

    // Indicate that the solver is initialized.
    d_is_initialized = true;

    IBTK_TIMER_STOP(t_initialize_solver_state);
    return;
}// initializeSolverState

void
SCPoissonHypreLevelSolver::deallocateSolverState()
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
SCPoissonHypreLevelSolver::enableLogging(
    bool enabled)
{
    d_enable_logging = enabled;
    return;
}// enableLogging

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
SCPoissonHypreLevelSolver::allocateHypreData()
{
    // Get the MPI communicator.
#ifdef HAVE_MPI
    MPI_Comm communicator = SAMRAI_MPI::getCommunicator();
#else
    MPI_Comm communicator;
#endif

    // Setup the hypre grid and variables and assemble the grid.
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_num);
    Pointer<CartesianGridGeometry<NDIM> > grid_geometry = d_hierarchy->getGridGeometry();
    const IntVector<NDIM>& ratio = level->getRatio();
    const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift(ratio);

    HYPRE_SStructGridCreate(communicator, NDIM, NPARTS, &d_grid);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        const Box<NDIM>& patch_box = level->getPatch(p())->getBox();
        Index<NDIM> lower = patch_box.lower();
        Index<NDIM> upper = patch_box.upper();
        HYPRE_SStructGridSetExtents(d_grid, PART, lower, upper);
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
    HYPRE_SStructGridSetPeriodic(d_grid, PART, hypre_periodic_shift);

#if (NDIM == 2)
    HYPRE_SStructVariable vartypes[NVARS] = {HYPRE_SSTRUCT_VARIABLE_XFACE , HYPRE_SSTRUCT_VARIABLE_YFACE};
#endif
#if (NDIM == 3)
    HYPRE_SStructVariable vartypes[NVARS] = {HYPRE_SSTRUCT_VARIABLE_XFACE , HYPRE_SSTRUCT_VARIABLE_YFACE , HYPRE_SSTRUCT_VARIABLE_ZFACE};
#endif
    HYPRE_SStructGridSetVariables(d_grid, PART, NVARS, vartypes);

    HYPRE_SStructGridAssemble(d_grid);

    // Allocate stencil data and set stencil offsets.
    if (d_constant_coefficients)
    {
        static const int stencil_sz = 2*NDIM+1;
#if (NDIM == 2)
        int stencil_offsets[stencil_sz][2] = { { -1, 0 }, { 0, -1}, { +1, 0 }, { 0, +1}, { 0, 0 } };
#endif
#if (NDIM == 3)
        int stencil_offsets[stencil_sz][3] = { { -1, 0, 0 }, { 0, -1, 0 }, { 0, 0, -1 }, { +1, 0, 0 }, { 0, +1, 0 }, { 0, 0, +1 }, { 0, 0, 0 } };
#endif
        for (int var = 0; var < NVARS; ++var)
        {
            HYPRE_SStructStencilCreate(NDIM, stencil_sz, &d_stencil[var]);
            for (int s = 0; s < stencil_sz; ++s)
            {
                HYPRE_SStructStencilSetEntry(d_stencil[var], s, stencil_offsets[s], var);
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name << "::allocateHypreData()\n"
                   << "  only constant-coefficient problems are currently supported" << std::endl);
    }

    // Allocate the hypre graph.
    HYPRE_SStructGraphCreate(communicator, d_grid, &d_graph);
    for (int var = 0; var < NVARS; ++var)
    {
        HYPRE_SStructGraphSetStencil(d_graph, PART, var, d_stencil[var]);
    }
    HYPRE_SStructGraphAssemble(d_graph);

    // Allocate the hypre matrix.
    HYPRE_SStructMatrixCreate(communicator, d_graph, &d_matrix);
    HYPRE_SStructMatrixInitialize(d_matrix);

    // Allocate the hypre vectors.
    HYPRE_SStructVectorCreate(communicator, d_grid, &d_sol_vec);
    HYPRE_SStructVectorInitialize(d_sol_vec);

    HYPRE_SStructVectorCreate(communicator, d_grid, &d_rhs_vec);
    HYPRE_SStructVectorInitialize(d_rhs_vec);
    return;
}// allocateHypreData

void
SCPoissonHypreLevelSolver::setMatrixCoefficients_constant_coefficients()
{
    ArrayDataBasicOps<NDIM,double> array_ops;
    static const IntVector<NDIM> no_ghosts = 0;
    static const IntVector<NDIM> gcw_to_fill = 1;

    // The constant problem coefficients.
    const double C = (d_poisson_spec.cIsZero() ? 0.0 : d_poisson_spec.getCConstant());
    const double D = d_poisson_spec.getDConstant();

    // Setup the stencil information.
    static const int stencil_sz = 2*NDIM+1;
    std::vector<int> stencil_indices(stencil_sz);
    for (int k = 0; k < stencil_sz; ++k)
    {
        stencil_indices[k] = k;
    }

    // Set the matrix coefficients.
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_num);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Array<BoundaryBox<NDIM> > physical_codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
        const int n_physical_codim1_boxes = physical_codim1_boxes.size();
        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();
        const double* const patch_x_lower = pgeom->getXLower();
        const double* const patch_x_upper = pgeom->getXUpper();
        const IntVector<NDIM>& ratio_to_level_zero = pgeom->getRatio();
        Array<Array<bool> > touches_regular_bdry(NDIM), touches_periodic_bdry(NDIM);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            touches_regular_bdry [axis].resizeArray(2);
            touches_periodic_bdry[axis].resizeArray(2);
            for (int upperlower = 0; upperlower < 2; ++upperlower)
            {
                touches_regular_bdry [axis][upperlower] = pgeom->getTouchesRegularBoundary( axis,upperlower);
                touches_periodic_bdry[axis][upperlower] = pgeom->getTouchesPeriodicBoundary(axis,upperlower);
            }
        }
        for (int var = 0; var < NVARS; ++var)
        {
            const unsigned int axis = var;
            Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);

            // Compute all matrix coefficients, including those that touch the
            // physical boundary, but temporarily ignoring physical boundary
            // conditions.
            IntVector<NDIM> directions = 0;
            directions(axis) = 1;
            SideData<NDIM,double> mat_vals_data(patch_box, stencil_sz, no_ghosts, directions);

            std::vector<double> mat_vals(stencil_sz,0.0);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const double dx_sq = dx[d]*dx[d];
                mat_vals[     d] += D/dx_sq;          // lower off-diagonal
                mat_vals[NDIM+d] += D/dx_sq;          // upper off-diagonal
                mat_vals[2*NDIM] -= mat_vals[     d]; // diagonal
                mat_vals[2*NDIM] -= mat_vals[NDIM+d]; // diagonal
            }
            mat_vals[2*NDIM] += C; // diagonal

            for (int k = 0; k < stencil_sz; ++k)
            {
                mat_vals_data.fill(mat_vals[k],k);
            }

            // Set physical boundary conditions along boundaries which are NOT
            // aligned with the data axis.
            //
            // NOTE: It important to set these values first to avoid
            // problems at corners in the physical domain.  In particular,
            // since Dirichlet boundary conditions for values located on the
            // physical boundary override all other boundary conditions, we
            // set those values last.
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                const unsigned int location_index   = bdry_box.getLocationIndex();
                const unsigned int bdry_normal_axis = location_index/2;
                const bool is_lower        = location_index%2 == 0;
                if (bdry_normal_axis != axis)
                {
                    const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
                    const BoundaryBox<NDIM> trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
                    const Box<NDIM> bc_coef_box = compute_tangential_extension(PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box), axis);

                    ArrayData<NDIM,double> acoef_data(bc_coef_box, 1);
                    ArrayData<NDIM,double> bcoef_data(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM,double> > acoef_data_ptr(&acoef_data,false);
                    Pointer<ArrayData<NDIM,double> > bcoef_data_ptr(&bcoef_data,false);
                    Pointer<ArrayData<NDIM,double> > gcoef_data_ptr(NULL       ,false);

                    // Temporarily reset the patch geometry object associated
                    // with the patch so that boundary conditions are set at the
                    // correct spatial locations.
                    blitz::TinyVector<double,NDIM> shifted_patch_x_lower, shifted_patch_x_upper;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        shifted_patch_x_lower[d] = patch_x_lower[d];
                        shifted_patch_x_upper[d] = patch_x_upper[d];
                    }
                    shifted_patch_x_lower[axis] -= 0.5*dx[axis];
                    shifted_patch_x_upper[axis] -= 0.5*dx[axis];
                    patch->setPatchGeometry(
                        new CartesianPatchGeometry<NDIM>(
                            ratio_to_level_zero, touches_regular_bdry, touches_periodic_bdry,
                            dx, shifted_patch_x_lower.data(), shifted_patch_x_upper.data()));

                    // Set the boundary condition coefficients.
                    ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[axis]);
                    if (extended_bc_coef != NULL) extended_bc_coef->setHomogeneousBc(true);
                    d_bc_coefs[axis]->setBcCoefs(
                        acoef_data_ptr, bcoef_data_ptr, gcoef_data_ptr, NULL,
                        *patch, trimmed_bdry_box, d_apply_time);
                    if (extended_bc_coef != NULL) extended_bc_coef->setHomogeneousBc(d_homogeneous_bc);

                    // Restore the original patch geometry object.
                    patch->setPatchGeometry(pgeom);

                    for (Box<NDIM>::Iterator b(bc_coef_box); b; b++)
                    {
                        // Modify the diagonal and off-diagonal entries to
                        // account for homogeneous boundary conditions.
                        //
                        // Here, we follow the same linear extrapolation
                        // approach implemented in class
                        // CartesianRobinBcHelper.  Namely, with
                        // u_i denoting the interior cell, u_o denoting the
                        // ghost cell, and u_b and u_n denoting the value and
                        // normal derivative of u at the boundary,
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
                        const Index<NDIM>& i = b();
                        const double& a = acoef_data(i,0);
                        const double& b = bcoef_data(i,0);
                        const double& h = dx[bdry_normal_axis];

                        Index<NDIM> i_intr = i;
                        if (is_lower)
                        {
                            i_intr(bdry_normal_axis) += 0;
                        }
                        else
                        {
                            i_intr(bdry_normal_axis) -= 1;
                        }
                        const SideIndex<NDIM> i_s(i_intr, axis, SideIndex<NDIM>::Lower);

                        if (is_lower)
                        {
                            mat_vals_data(i_s,2*NDIM) += mat_vals_data(i_s,     bdry_normal_axis)*(-(a*h-2.0*b)/(a*h+2.0*b));
                            mat_vals_data(i_s,     bdry_normal_axis) = 0.0;
                        }
                        else
                        {
                            mat_vals_data(i_s,2*NDIM) += mat_vals_data(i_s,NDIM+bdry_normal_axis)*(-(a*h-2.0*b)/(a*h+2.0*b));
                            mat_vals_data(i_s,NDIM+bdry_normal_axis) = 0.0;
                        }
                    }
                }
            }

            // Set physical boundary conditions along boundaries which ARE
            // aligned with the data axis.
            //
            // NOTE: It important to set these values last to avoid problems at
            // corners in the physical domain.  In particular, since Dirichlet
            // boundary conditions for values located on the physical boundary
            // override all other boundary conditions, we set those values last.
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                const unsigned int location_index   = bdry_box.getLocationIndex();
                const unsigned int bdry_normal_axis = location_index/2;
                const bool is_lower        = location_index%2 == 0;
                if (bdry_normal_axis == axis)
                {
                    const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
                    const BoundaryBox<NDIM> trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
                    const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

                    ArrayData<NDIM,double> acoef_data(bc_coef_box, 1);
                    ArrayData<NDIM,double> bcoef_data(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM,double> > acoef_data_ptr(&acoef_data,false);
                    Pointer<ArrayData<NDIM,double> > bcoef_data_ptr(&bcoef_data,false);
                    Pointer<ArrayData<NDIM,double> > gcoef_data_ptr(NULL       ,false);

                    // Set the boundary condition coefficients.
                    ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[axis]);
                    if (extended_bc_coef != NULL) extended_bc_coef->setHomogeneousBc(true);
                    d_bc_coefs[axis]->setBcCoefs(
                        acoef_data_ptr, bcoef_data_ptr, gcoef_data_ptr, NULL,
                        *patch, trimmed_bdry_box, d_apply_time);
                    if (extended_bc_coef != NULL) extended_bc_coef->setHomogeneousBc(d_homogeneous_bc);

                    for (Box<NDIM>::Iterator b(bc_coef_box); b; b++)
                    {
                        const Index<NDIM>& i = b();
                        const double& a = acoef_data(i,0);
                        const double& b = bcoef_data(i,0);
                        const bool dirichlet_bc = MathUtilities<double>::equalEps(a,1.0);
                        const bool neumann_bc = MathUtilities<double>::equalEps(b,1.0);
#ifdef DEBUG_CHECK_ASSERTIONS
                        TBOX_ASSERT((dirichlet_bc || neumann_bc) && !(dirichlet_bc && neumann_bc));
#endif
                        const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                        if (dirichlet_bc)
                        {
                            // Along the physical boundary, Dirichlet boundary
                            // conditions are imposed exactly and set in the
                            // right-hand-side.
                            for (int k = 0; k < stencil_sz; ++k)
                            {
                                mat_vals_data(i_s,k) = 0.0;
                            }
                            mat_vals_data(i_s,2*NDIM) = 1.0;
                        }
                        else if (neumann_bc)
                        {
                            // Along Neumann boundaries, we can either use a
                            // symmetric (linear) or non-symmetric (quadratic)
                            // extrapolation.
                            //
                            // With u_i denoting the interior value, u_b
                            // denoting the boundary value, and u_o denoting the
                            // ghost value, the symmetric (linear) extrapolation
                            // is:
                            //
                            //     u_n = (u_o - u_b)/h
                            //
                            // so that:
                            //
                            //     b*u_n = 0 ===> u_o = u_b.
                            //
                            // The non-symmetric (quadratic) extrapolation is:
                            //
                            //     u_n = (u_o - u_i)/(2*h)
                            //
                            // so that:
                            //
                            //     b*u_n = 0 ===> u_o = u_i.
                            //
                            // \todo Add options to choose the boundary
                            // treatment.
                            if (is_lower)
                            {
                                mat_vals_data(i_s,NDIM+bdry_normal_axis) += mat_vals_data(i_s,     bdry_normal_axis);
                                mat_vals_data(i_s,     bdry_normal_axis)  = 0.0;
                            }
                            else
                            {
                                mat_vals_data(i_s,     bdry_normal_axis) += mat_vals_data(i_s,NDIM+bdry_normal_axis);
                                mat_vals_data(i_s,NDIM+bdry_normal_axis)  = 0.0;
                            }
                        }
                    }
                }
            }

            // Set the matrix coefficients in the hypre matrix.
            for (Box<NDIM>::Iterator b(side_box); b; b++)
            {
                // NOTE: In hypre, face-centered values are associated with the
                // cell index located on the "lower" side of the face.
                Index<NDIM> i = b();
                static const int lower = 0;
                const SideIndex<NDIM> i_s(i,axis,lower);
                i(axis) -= 1;
                for (int k = 0; k < stencil_sz; ++k)
                {
                    HYPRE_SStructMatrixSetValues(d_matrix, PART, i, var, 1, &stencil_indices[k], &mat_vals_data(i_s,k));
                }
            }
        }
    }

    // Assemble the hypre matrix.
    HYPRE_SStructMatrixAssemble(d_matrix);
    return;
}// setMatrixCoefficients_constant_coefficients

void
SCPoissonHypreLevelSolver::setupHypreSolver()
{
    // Get the MPI communicator.
#ifdef HAVE_MPI
    MPI_Comm communicator = SAMRAI_MPI::getCommunicator();
#else
    MPI_Comm communicator;
#endif

    // Determine the split solver type.
    int split_solver_type_id = -1;
    if (d_solver_type == "Split" || d_precond_type == "Split")
    {
        if (d_split_solver_type == "PFMG")
        {
            split_solver_type_id = HYPRE_PFMG;
        }
        else if (d_split_solver_type == "SMG")
        {
            split_solver_type_id = HYPRE_SMG;
        }
        else if (d_split_solver_type == "Jacobi")
        {
            split_solver_type_id = HYPRE_Jacobi;
        }
        else
        {
            TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                       << "  unknown split solver type: " << d_split_solver_type << std::endl);
        }
    }

    // When using a Krylov method, setup the preconditioner.
    if (d_solver_type == "PCG" || d_solver_type == "GMRES" || d_solver_type == "FlexGMRES" || d_solver_type == "LGMRES" || d_solver_type == "BiCGSTAB")
    {
        if (d_precond_type == "SysPFMG")
        {
            HYPRE_SStructSysPFMGCreate(communicator, &d_precond);
            HYPRE_SStructSysPFMGSetMaxIter(d_precond, 1);
            HYPRE_SStructSysPFMGSetTol(d_precond, 0.0);
            HYPRE_SStructSysPFMGSetZeroGuess(d_precond);
            HYPRE_SStructSysPFMGSetRelaxType(d_precond, d_relax_type);
            HYPRE_SStructSysPFMGSetNumPreRelax(d_precond, d_num_pre_relax_steps);
            HYPRE_SStructSysPFMGSetNumPostRelax(d_precond, d_num_post_relax_steps);
            HYPRE_SStructSysPFMGSetSkipRelax(d_precond, d_skip_relax);
        }
        else if (d_precond_type == "Split")
        {

            HYPRE_SStructSplitCreate(communicator, &d_precond);
            HYPRE_SStructSplitSetMaxIter(d_precond, 1);
            HYPRE_SStructSplitSetTol(d_precond, 0.0);
            HYPRE_SStructSplitSetZeroGuess(d_precond);
            HYPRE_SStructSplitSetStructSolver(d_precond, split_solver_type_id);
        }
    }

    // Setup the solver.
    if (d_solver_type == "SysPFMG")
    {
        HYPRE_SStructSysPFMGCreate(communicator, &d_solver);
        HYPRE_SStructSysPFMGSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructSysPFMGSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructSysPFMGSetRelChange(d_solver, d_rel_change);
        HYPRE_SStructSysPFMGSetRelaxType(d_solver, d_relax_type);
        HYPRE_SStructSysPFMGSetNumPreRelax(d_solver, d_num_pre_relax_steps);
        HYPRE_SStructSysPFMGSetNumPostRelax(d_solver, d_num_post_relax_steps);
        HYPRE_SStructSysPFMGSetSkipRelax(d_solver, d_skip_relax);
        if (d_initial_guess_nonzero)
        {
            HYPRE_SStructSysPFMGSetNonZeroGuess(d_solver);
        }
        else
        {
            HYPRE_SStructSysPFMGSetZeroGuess(d_solver);
        }
        HYPRE_SStructSysPFMGSetup(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
    }
    else if (d_solver_type == "Split")
    {
        HYPRE_SStructSplitCreate(communicator, &d_solver);
        HYPRE_SStructSplitSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructSplitSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructSplitSetStructSolver(d_solver, split_solver_type_id);
        if (d_initial_guess_nonzero)
        {
            HYPRE_SStructSplitSetNonZeroGuess(d_solver);
        }
        else
        {
            HYPRE_SStructSplitSetZeroGuess(d_solver);
        }
        HYPRE_SStructSplitSetup(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
    }
    else if (d_solver_type == "PCG")
    {
        HYPRE_SStructPCGCreate(communicator, &d_solver);
        HYPRE_SStructPCGSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructPCGSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructPCGSetAbsoluteTol(d_solver, d_abs_residual_tol);
        HYPRE_SStructPCGSetTwoNorm(d_solver, d_two_norm);
        HYPRE_SStructPCGSetRelChange(d_solver, d_rel_change);
        if (d_precond_type == "SysPFMG")
        {
            HYPRE_SStructPCGSetPrecond(d_solver,
                                       HYPRE_SStructSysPFMGSolve,
                                       HYPRE_SStructSysPFMGSetup,
                                       d_precond);
        }
        else if (d_precond_type == "Split")
        {
            HYPRE_SStructPCGSetPrecond(d_solver,
                                       HYPRE_SStructSplitSolve,
                                       HYPRE_SStructSplitSetup,
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
        HYPRE_SStructPCGSetup(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
    }
    else if (d_solver_type == "GMRES")
    {
        HYPRE_SStructGMRESCreate(communicator, &d_solver);
        HYPRE_SStructGMRESSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructGMRESSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructGMRESSetAbsoluteTol(d_solver, d_abs_residual_tol);
        if (d_precond_type == "SysPFMG")
        {
            HYPRE_SStructGMRESSetPrecond(d_solver,
                                         HYPRE_SStructSysPFMGSolve,
                                         HYPRE_SStructSysPFMGSetup,
                                         d_precond);
        }
        else if (d_precond_type == "Split")
        {
            HYPRE_SStructGMRESSetPrecond(d_solver,
                                         HYPRE_SStructSplitSolve,
                                         HYPRE_SStructSplitSetup,
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
        HYPRE_SStructGMRESSetup(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
    }
    else if (d_solver_type == "FlexGMRES")
    {
        HYPRE_SStructFlexGMRESCreate(communicator, &d_solver);
        HYPRE_SStructFlexGMRESSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructFlexGMRESSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructFlexGMRESSetAbsoluteTol(d_solver, d_abs_residual_tol);
        if (d_precond_type == "SysPFMG")
        {
            HYPRE_SStructFlexGMRESSetPrecond(d_solver,
                                             HYPRE_SStructSysPFMGSolve,
                                             HYPRE_SStructSysPFMGSetup,
                                             d_precond);
        }
        else if (d_precond_type == "Split")
        {
            HYPRE_SStructFlexGMRESSetPrecond(d_solver,
                                             HYPRE_SStructSplitSolve,
                                             HYPRE_SStructSplitSetup,
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
        HYPRE_SStructFlexGMRESSetup(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
    }
    else if (d_solver_type == "LGMRES")
    {
        HYPRE_SStructLGMRESCreate(communicator, &d_solver);
        HYPRE_SStructLGMRESSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructLGMRESSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructLGMRESSetAbsoluteTol(d_solver, d_abs_residual_tol);
        if (d_precond_type == "SysPFMG")
        {
            HYPRE_SStructLGMRESSetPrecond(d_solver,
                                          HYPRE_SStructSysPFMGSolve,
                                          HYPRE_SStructSysPFMGSetup,
                                          d_precond);
        }
        else if (d_precond_type == "Split")
        {
            HYPRE_SStructLGMRESSetPrecond(d_solver,
                                          HYPRE_SStructSplitSolve,
                                          HYPRE_SStructSplitSetup,
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
        HYPRE_SStructLGMRESSetup(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
    }
    else if (d_solver_type == "BiCGSTAB")
    {
        HYPRE_SStructBiCGSTABCreate(communicator, &d_solver);
        HYPRE_SStructBiCGSTABSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructBiCGSTABSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructBiCGSTABSetAbsoluteTol(d_solver, d_abs_residual_tol);
        if (d_precond_type == "SysPFMG")
        {
            HYPRE_SStructBiCGSTABSetPrecond(d_solver,
                                            HYPRE_SStructSysPFMGSolve,
                                            HYPRE_SStructSysPFMGSetup,
                                            d_precond);
        }
        else if (d_precond_type == "Split")
        {
            HYPRE_SStructBiCGSTABSetPrecond(d_solver,
                                            HYPRE_SStructSplitSolve,
                                            HYPRE_SStructSplitSetup,
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
        HYPRE_SStructBiCGSTABSetup(d_solver,d_matrix, d_rhs_vec, d_sol_vec);
    }
    else
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                   << "  unknown solver type: " << d_solver_type << std::endl);
    }
    return;
}// setupHypreSolver

bool
SCPoissonHypreLevelSolver::solveSystem(
    const int x_idx,
    const int b_idx)
{
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_level_num);

    // Modify right-hand-side data to account for boundary conditions and copy
    // solution and right-hand-side data to hypre structures.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

        // Copy the solution data into the hypre vector, including ghost cell
        // values
        const Box<NDIM> x_ghost_box = Box<NDIM>::grow(patch_box, 1);
        Pointer<SideData<NDIM,double> > x_data = patch->getPatchData(x_idx);
        copyToHypre(d_sol_vec, x_data, x_ghost_box);

        // Modify the right-hand-side data to account for any boundary
        // conditions and copy the right-hand-side into the hypre vector.
        Pointer<SideData<NDIM,double> > b_data = patch->getPatchData(b_idx);
        if (pgeom->intersectsPhysicalBoundary())
        {
            b_data = new SideData<NDIM,double>(
                b_data->getBox(), b_data->getDepth(),
                b_data->getGhostCellWidth());
            b_data->copy(*(patch->getPatchData(b_idx)));

            const Array<BoundaryBox<NDIM> > physical_codim1_boxes =
                PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);

            if (d_constant_coefficients)
            {
                adjustBoundaryRhsEntries_constant_coefficients(b_data, patch, physical_codim1_boxes);
            }
            else
            {
                TBOX_ERROR(d_object_name << "::solveSystem()\n"
                           << "  only constant-coefficient problems are currently supported" << std::endl);
            }
        }
        copyToHypre(d_rhs_vec, b_data, patch_box);
    }

    // Assemble the hypre vectors.
    HYPRE_SStructVectorAssemble(d_sol_vec);
    HYPRE_SStructVectorAssemble(d_rhs_vec);

    // Solve the system.
    IBTK_TIMER_START(t_solve_system_hypre);

    if (d_solver_type == "SysPFMG")
    {
        HYPRE_SStructSysPFMGSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructSysPFMGSetTol(d_solver, d_rel_residual_tol);
        if (d_initial_guess_nonzero)
        {
            HYPRE_SStructSysPFMGSetNonZeroGuess(d_solver);
        }
        else
        {
            HYPRE_SStructSysPFMGSetZeroGuess(d_solver);
        }
        HYPRE_SStructSysPFMGSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_SStructSysPFMGGetNumIterations(d_solver, &d_current_its);
        HYPRE_SStructSysPFMGGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "Split")
    {
        HYPRE_SStructSplitSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructSplitSetTol(d_solver, d_rel_residual_tol);
        if (d_initial_guess_nonzero)
        {
            HYPRE_SStructSplitSetNonZeroGuess(d_solver);
        }
        else
        {
            HYPRE_SStructSplitSetZeroGuess(d_solver);
        }
        HYPRE_SStructSplitSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_SStructSplitGetNumIterations(d_solver, &d_current_its);
        HYPRE_SStructSplitGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "PCG")
    {
        HYPRE_SStructPCGSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructPCGSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructPCGSetAbsoluteTol(d_solver, d_abs_residual_tol);
        HYPRE_SStructPCGSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_SStructPCGGetNumIterations(d_solver, &d_current_its);
        HYPRE_SStructPCGGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "GMRES")
    {
        HYPRE_SStructGMRESSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructGMRESSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructGMRESSetAbsoluteTol(d_solver, d_abs_residual_tol);
        HYPRE_SStructGMRESSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_SStructGMRESGetNumIterations(d_solver, &d_current_its);
        HYPRE_SStructGMRESGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "FlexGMRES")
    {
        HYPRE_SStructFlexGMRESSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructFlexGMRESSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructFlexGMRESSetAbsoluteTol(d_solver, d_abs_residual_tol);
        HYPRE_SStructFlexGMRESSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_SStructFlexGMRESGetNumIterations(d_solver, &d_current_its);
        HYPRE_SStructFlexGMRESGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "LGMRES")
    {
        HYPRE_SStructLGMRESSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructLGMRESSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructLGMRESSetAbsoluteTol(d_solver, d_abs_residual_tol);
        HYPRE_SStructLGMRESSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_SStructLGMRESGetNumIterations(d_solver, &d_current_its);
        HYPRE_SStructLGMRESGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "BiCGSTAB")
    {
        HYPRE_SStructBiCGSTABSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructBiCGSTABSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructBiCGSTABSetAbsoluteTol(d_solver, d_abs_residual_tol);
        HYPRE_SStructBiCGSTABSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_SStructBiCGSTABGetNumIterations(d_solver, &d_current_its);
        HYPRE_SStructBiCGSTABGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }

    IBTK_TIMER_STOP(t_solve_system_hypre);

    // Pull the solution vector out of the hypre structures.
    HYPRE_SStructVectorGather(d_sol_vec);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM,double> > x_data = patch->getPatchData(x_idx);
        copyFromHypre(x_data, d_sol_vec, patch_box);
    }
    return (d_current_residual_norm <= d_rel_residual_tol || d_current_residual_norm <= d_abs_residual_tol);
}// solveSystem

void
SCPoissonHypreLevelSolver::copyToHypre(
    HYPRE_SStructVector vector,
    const Pointer<SideData<NDIM,double> > src_data,
    const Box<NDIM>& box)
{
    const bool copy_data = src_data->getGhostBox() != box;
    Pointer<SideData<NDIM,double> > hypre_data =
        (copy_data ? Pointer<SideData<NDIM,double> >(new SideData<NDIM,double>(box,1,0)) : src_data);

    if (copy_data) hypre_data->copyOnBox(*src_data,box);

    for (int var = 0; var < NVARS; ++var)
    {
        const unsigned int axis = var;
        Index<NDIM> lower = box.lower();
        lower(axis) -= 1;
        Index<NDIM> upper = box.upper();
        HYPRE_SStructVectorSetBoxValues(vector,PART,lower,upper,var,hypre_data->getPointer(axis));
    }
    return;
}// copyToHypre

void
SCPoissonHypreLevelSolver::copyFromHypre(
    Pointer<SideData<NDIM,double> > dst_data,
    HYPRE_SStructVector vector,
    const Box<NDIM>& box)
{
    const bool copy_data = dst_data->getGhostBox() != box;
    Pointer<SideData<NDIM,double> > hypre_data =
        (copy_data ? Pointer<SideData<NDIM,double> >(new SideData<NDIM,double>(box,1,0)) : dst_data);

    for (int var = 0; var < NVARS; ++var)
    {
        const unsigned int axis = var;
        Index<NDIM> lower = box.lower();
        lower(axis) -= 1;
        Index<NDIM> upper = box.upper();
        HYPRE_SStructVectorGetBoxValues(vector,PART,lower,upper,var,hypre_data->getPointer(axis));
    }
    if (copy_data) dst_data->copyOnBox(*hypre_data,box);
    return;
}// copyFromHypre

void
SCPoissonHypreLevelSolver::destroyHypreSolver()
{
    // Destroy the solver.
    if (d_solver_type == "SysPFMG")
    {
        HYPRE_SStructSysPFMGDestroy(d_solver);
    }
    else if (d_solver_type == "Split")
    {
        HYPRE_SStructSplitDestroy(d_solver);
    }
    else if (d_solver_type == "PCG")
    {
        HYPRE_SStructPCGDestroy(d_solver);
    }
    else if (d_solver_type == "GMRES")
    {
        HYPRE_SStructGMRESDestroy(d_solver);
    }
    else if (d_solver_type == "FlexGMRES")
    {
        HYPRE_SStructFlexGMRESDestroy(d_solver);
    }
    else if (d_solver_type == "LGMRES")
    {
        HYPRE_SStructLGMRESDestroy(d_solver);
    }
    else if (d_solver_type == "BiCGSTAB")
    {
        HYPRE_SStructBiCGSTABDestroy(d_solver);
    }

    // When using a Krylov method, destroy the preconditioner.
    if (d_solver_type == "PCG" || d_solver_type == "GMRES" || d_solver_type == "FlexGMRES" || d_solver_type == "LGMRES" || d_solver_type == "BiCGSTAB")
    {
        if (d_precond_type == "SysPFMG")
        {
            HYPRE_SStructSysPFMGDestroy(d_precond);
        }
        else if (d_precond_type == "Split")
        {
            HYPRE_SStructSplitDestroy(d_precond);
        }
    }

    // Set the solver and preconditioner pointers to NULL.
    d_solver  = NULL;
    d_precond = NULL;
    return;
}// destroyHypreSolver

void
SCPoissonHypreLevelSolver::deallocateHypreData()
{
    if (d_graph  ) HYPRE_SStructGraphDestroy(d_graph);
    for (int var = 0; var < NVARS; ++var)
    {
        if (d_stencil[var]) HYPRE_SStructStencilDestroy(d_stencil[var]);
    }
    if (d_grid   ) HYPRE_SStructGridDestroy(d_grid);
    if (d_matrix ) HYPRE_SStructMatrixDestroy(d_matrix);
    if (d_sol_vec) HYPRE_SStructVectorDestroy(d_sol_vec);
    if (d_rhs_vec) HYPRE_SStructVectorDestroy(d_rhs_vec);
    d_grid    = NULL;
    for (int var = 0; var < NVARS; ++var)
    {
        d_stencil[var] = NULL;
    }
    d_matrix  = NULL;
    d_sol_vec = NULL;
    d_rhs_vec = NULL;
    return;
}// deallocateHypreData

void
SCPoissonHypreLevelSolver::adjustBoundaryRhsEntries_constant_coefficients(
    Pointer<SideData<NDIM,double> > rhs_data,
    const Pointer<Patch<NDIM> > patch,
    const Array<BoundaryBox<NDIM> >& physical_codim1_boxes)
{
    static const IntVector<NDIM> gcw_to_fill = 1;
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const Box<NDIM>& patch_box = patch->getBox();

    // Modify the rhs entries to account for homogeneous boundary conditions.
    //
    // NOTE: For values located on Dirichlet boundaries, we set the boundary
    // value explicitly in the right-hand-side.  Consequently, we must modify
    // the boundary values even for homogeneous boundary conditions.
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
        const unsigned int location_index   = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index/2;
//      const bool is_lower        = location_index%2 == 0;
        for (int var = 0; var < NVARS; ++var)
        {
            const unsigned int axis = var;
            if (bdry_normal_axis == axis)
            {
                const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
                const BoundaryBox<NDIM> trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
                const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

                ArrayData<NDIM,double> acoef_data(bc_coef_box, 1);
                ArrayData<NDIM,double> bcoef_data(bc_coef_box, 1);
                Pointer<ArrayData<NDIM,double> > acoef_data_ptr(&acoef_data,false);
                Pointer<ArrayData<NDIM,double> > bcoef_data_ptr(&bcoef_data,false);
                Pointer<ArrayData<NDIM,double> > gcoef_data_ptr(NULL       ,false);

                // Set the boundary condition coefficients.
                ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[axis]);
                if (extended_bc_coef != NULL) extended_bc_coef->setHomogeneousBc(d_homogeneous_bc);
                d_bc_coefs[axis]->setBcCoefs(
                    acoef_data_ptr, bcoef_data_ptr, gcoef_data_ptr, NULL,
                    *patch, trimmed_bdry_box, d_apply_time);

                for (Box<NDIM>::Iterator b(bc_coef_box); b; b++)
                {
                    const Index<NDIM>& i = b();
                    const SideIndex<NDIM> i_s_bdry(i, bdry_normal_axis, SideIndex<NDIM>::Lower);
                    const double& a = acoef_data(i,0);
                    const double& b = bcoef_data(i,0);
                    const bool dirichlet_bc = MathUtilities<double>::equalEps(a,1.0);
                    const bool neumann_bc = MathUtilities<double>::equalEps(b,1.0);
#ifdef DEBUG_CHECK_ASSERTIONS
                    TBOX_ASSERT((dirichlet_bc || neumann_bc) && !(dirichlet_bc && neumann_bc));
#endif
                    if (dirichlet_bc)
                    {
                        (*rhs_data)(i_s_bdry) = 0.0;
                    }
                    else if (neumann_bc)
                    {
                        // intentionally blank.
                    }
                }
            }
        }
    }

    // The remaining modifications are made only in the case of inhomogeneous
    // boundary conditions.
    if (d_homogeneous_bc) return;

    // Modify the rhs entries to account for inhomogeneous boundary conditions.
    TBOX_ERROR(d_object_name << "::adjustBoundaryRhsEntries_constant_coefficients()\n"
               << "  only homogeneous boundary conditions are currently supported" << std::endl);
    return;
}// adjustBoundaryRhsEntries_constant_coefficients

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
