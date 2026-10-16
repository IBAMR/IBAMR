// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/GeneralSolver.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PoissonUtilities.h>
#include <ibtk/SCPoissonHypreLevelSolver.h>
#include <ibtk/solver_utilities.h>

#include <tbox/Array.h>
#include <tbox/Database.h>
#include <tbox/PIO.h>
#include <tbox/Pointer.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

#include <ArrayDataBasicOps.h>
#include <ArrayDataNormOpsReal.h>
#include <BoundaryBox.h>
#include <Box.h>
#include <BoxList.h>
#include <CartesianPatchGeometry.h>
#include <CoarseFineBoundary.h>
#include <GridGeometry.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>
#include <SideVariable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <ibtk/namespaces.h> // IWYU pragma: keep

IBTK_DISABLE_EXTRA_WARNINGS
#include <HYPRE_sstruct_ls.h>
#include <HYPRE_sstruct_mv.h>
IBTK_ENABLE_EXTRA_WARNINGS

#include <mpi.h>

#include <algorithm>
#include <array>
#include <memory>
#include <numeric>
#include <ostream>
#include <string>
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

// Return the period to give hypre in the periodic direction d of a level: the period of the physical domain if the
// level extends across the periodic boundary, else zero. hypre treats a periodic grid as covering exactly one period
// and does not store the faces at the upper end of the index extent of the grid, which is wrong for a level that does
// not extend across the periodic boundary.
int
get_hypre_period(const PatchLevel<NDIM>& level, const std::string& object_name, const unsigned int d)
{
    const int period = level.getGridGeometry()->getPeriodicShift(level.getRatio())(d);
    if (period == 0)
    {
        return 0;
    }
    const Box<NDIM> domain = level.getPhysicalDomain().getBoundingBox();
    Box<NDIM> lower_layer = domain;
    lower_layer.upper(d) = domain.lower(d);
    Box<NDIM> upper_layer = domain;
    upper_layer.lower(d) = domain.upper(d);
    IntVector<NDIM> to_lower_layer(0);
    to_lower_layer(d) = domain.lower(d) - domain.upper(d);

    // The cells of the level in the first and last layer of cells of the domain in direction d, with the last layer
    // shifted by one period, so that the cells on the two sides of a periodic boundary coincide.
    BoxList<NDIM> lower_cells, upper_cells;
    for (int n = 0; n < level.getBoxes().getNumberOfBoxes(); ++n)
    {
        const Box<NDIM> lower_box = level.getBoxes()[n] * lower_layer;
        if (!lower_box.empty())
        {
            lower_cells.appendItem(lower_box);
        }
        const Box<NDIM> upper_box = level.getBoxes()[n] * upper_layer;
        if (!upper_box.empty())
        {
            upper_cells.appendItem(Box<NDIM>::shift(upper_box, to_lower_layer));
        }
    }
    if (lower_cells.isEmpty() || upper_cells.isEmpty())
    {
        return 0;
    }
    BoxList<NDIM> upper_only(upper_cells);
    upper_only.removeIntersections(lower_cells);
    BoxList<NDIM> lower_only(lower_cells);
    lower_only.removeIntersections(upper_cells);
    if (!upper_only.isEmpty() || !lower_only.isEmpty())
    {
        TBOX_ERROR(object_name << "::initializeSolverState():\n"
                               << "  patches of level " << level.getLevelNumber()
                               << " adjacent to the two periodic boundaries in direction " << d
                               << " are not periodic images of each other\n"
                               << "  the patches adjacent to the lower and upper boundary must be periodic images, or "
                                  "one of the two boundaries must have no adjacent patch\n");
    }
    return period;
}
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

SCPoissonHypreLevelSolver::SCPoissonHypreLevelSolver(const std::string& object_name,
                                                     Pointer<Database> input_db,
                                                     const std::string& /*default_options_prefix*/)
{
    if (NDIM == 1 || NDIM > 3)
    {
        TBOX_ERROR(d_object_name << "::SCPoissonHypreLevelSolver()"
                                 << "  hypre solvers are only provided for 2D and 3D problems" << std::endl);
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

        if (d_solver_type == "Split" || d_precond_type == "Split")
        {
            if (input_db->keyExists("split_solver_type"))
                d_split_solver_type = input_db->getString("split_solver_type");
        }

        if (d_solver_type == "PCG")
        {
            if (input_db->keyExists("two_norm")) d_two_norm = input_db->getInteger("two_norm");
        }
    }

    // Setup Timers.
    IBTK_DO_ONCE(t_solve_system =
                     TimerManager::getManager()->getTimer("IBTK::SCPoissonHypreLevelSolver::solveSystem()");
                 t_solve_system_hypre =
                     TimerManager::getManager()->getTimer("IBTK::SCPoissonHypreLevelSolver::solveSystem()[hypre]");
                 t_initialize_solver_state =
                     TimerManager::getManager()->getTimer("IBTK::SCPoissonHypreLevelSolver::initializeSolverState()");
                 t_deallocate_solver_state =
                     TimerManager::getManager()->getTimer("IBTK::SCPoissonHypreLevelSolver::deallocateSolverState()"););
    return;
} // SCPoissonHypreLevelSolver

SCPoissonHypreLevelSolver::~SCPoissonHypreLevelSolver()
{
    if (d_is_initialized) deallocateSolverState();
    return;
} // ~SCPoissonHypreLevelSolver

bool
SCPoissonHypreLevelSolver::solveSystem(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& b)
{
    IBTK_TIMER_START(t_solve_system);

    if (d_enable_logging) plog << d_object_name << "::solveSystem():" << std::endl;

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x, b);

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
SCPoissonHypreLevelSolver::initializeSolverState(const SAMRAIVectorReal<NDIM, double>& x,
                                                 const SAMRAIVectorReal<NDIM, double>& b)
{
    IBTK_TIMER_START(t_initialize_solver_state);

// Rudimentary error checking.
#if !defined(NDEBUG)
    if (x.getNumberOfComponents() != b.getNumberOfComponents())
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  vectors must have the same number of components" << std::endl);
    }

    const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy = x.getPatchHierarchy();
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
        if (!patch_hierarchy->getPatchLevel(ln))
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
    TBOX_ASSERT(d_level_num == x.getFinestLevelNumber());
    d_level = d_hierarchy->getPatchLevel(d_level_num);
    if (d_level_num > 0)
    {
        d_cf_boundary = new CoarseFineBoundary<NDIM>(*d_hierarchy, d_level_num, IntVector<NDIM>(1));
    }

    // Allocate and initialize the hypre data structures.
    allocateHypreData();
    setMatrixCoefficients();
    setupHypreSolver();

    // Indicate that the solver is initialized.
    d_is_initialized = true;

    IBTK_TIMER_STOP(t_initialize_solver_state);
    return;
} // initializeSolverState

void
SCPoissonHypreLevelSolver::deallocateSolverState()
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_deallocate_solver_state);

    // Deallocate the hypre data structures.
    destroyHypreSolver();
    deallocateHypreData();
    deallocateNullSpaceData();

    // Indicate that the solver is NOT initialized.
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_solver_state);
    return;
} // deallocateSolverState

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
SCPoissonHypreLevelSolver::allocateNullSpaceData()
{
    // Compute the control-volume weights, so that they are available when the weighted mean is removed.
    getNullSpaceMathOps()->getSideWeightPatchDescriptorIndex();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    if (d_rhs_idx == IBTK::invalid_index)
    {
        Pointer<VariableContext> context = var_db->getContext(d_object_name + "::CONTEXT");
        const std::string name = d_object_name + "::side_rhs";
        if (var_db->checkVariableExists(name))
        {
            Pointer<SideVariable<NDIM, double>> var = var_db->getVariable(name);
            d_rhs_idx = var_db->mapVariableAndContextToIndex(var, context);
        }
        else
        {
            Pointer<SideVariable<NDIM, double>> var = new SideVariable<NDIM, double>(name, 1);
            d_rhs_idx = var_db->registerVariableAndContext(var, context, IntVector<NDIM>(0));
        }
    }
    if (!d_level->checkAllocated(d_rhs_idx))
    {
        d_level->allocatePatchData(d_rhs_idx);
    }
    return;
} // allocateNullSpaceData

void
SCPoissonHypreLevelSolver::deallocateNullSpaceData()
{
    if (d_rhs_idx != IBTK::invalid_index && d_level->checkAllocated(d_rhs_idx))
    {
        d_level->deallocatePatchData(d_rhs_idx);
    }
    deallocateOwnMathOps();
    return;
} // deallocateNullSpaceData

Pointer<HierarchyMathOps>
SCPoissonHypreLevelSolver::getNullSpaceMathOps()
{
    // A given object provides the weights of this level only if its level range is this level.
    if (d_hier_math_ops && d_hier_math_ops->getCoarsestLevelNumber() == d_level_num &&
        d_hier_math_ops->getFinestLevelNumber() == d_level_num)
    {
        deallocateOwnMathOps();
        return d_hier_math_ops;
    }
    if (!d_level_math_ops)
    {
        d_level_math_ops =
            new HierarchyMathOps(d_object_name + "::hier_math_ops", d_hierarchy, d_level_num, d_level_num);
    }
    return d_level_math_ops;
} // getNullSpaceMathOps

void
SCPoissonHypreLevelSolver::deallocateOwnMathOps()
{
    if (d_level_math_ops)
    {
        const int wgt_idx = d_level_math_ops->getSideWeightPatchDescriptorIndex();
        if (d_level->checkAllocated(wgt_idx))
        {
            d_level->deallocatePatchData(wgt_idx);
        }
        d_level_math_ops.setNull();
    }
    return;
} // deallocateOwnMathOps

void
SCPoissonHypreLevelSolver::removeMeanFromRhs()
{
    // The data operations of the hierarchy act on all the axes together, so use those of the arrays of each axis.
    ArrayDataNormOpsReal<NDIM, double> norm_ops;
    ArrayDataBasicOps<NDIM, double> basic_ops;
    const int wgt_idx = getNullSpaceMathOps()->getSideWeightPatchDescriptorIndex();

    // Sum, over the sides normal to each axis, the weighted values and the weights.
    std::array<double, 2 * NDIM> sums;
    sums.fill(0.0);
    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = d_level->getPatch(p());
        Pointer<SideData<NDIM, double>> rhs_data = patch->getPatchData(d_rhs_idx);
        Pointer<SideData<NDIM, double>> wgt_data = patch->getPatchData(wgt_idx);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch->getBox(), axis);
            sums[axis] += norm_ops.integral(rhs_data->getArrayData(axis), wgt_data->getArrayData(axis), side_box);
            sums[NDIM + axis] +=
                norm_ops.sumControlVolumes(wgt_data->getArrayData(axis), wgt_data->getArrayData(axis), side_box);
        }
    }
    IBTK_MPI::sumReduction(sums.data(), static_cast<int>(sums.size()));

    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = d_level->getPatch(p());
        Pointer<SideData<NDIM, double>> rhs_data = patch->getPatchData(d_rhs_idx);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch->getBox(), axis);
            basic_ops.addScalar(
                rhs_data->getArrayData(axis), rhs_data->getArrayData(axis), -sums[axis] / sums[NDIM + axis], side_box);
        }
    }
    return;
} // removeMeanFromRhs

void
SCPoissonHypreLevelSolver::allocateHypreData()
{
    initialize_hypre();

    // Get the MPI communicator.
    MPI_Comm communicator = IBTK_MPI::getCommunicator();

    // Setup the hypre grid and variables and assemble the grid.
    HYPRE_SStructGridCreate(communicator, NDIM, NPARTS, &d_grid);
    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        const Box<NDIM>& patch_box = d_level->getPatch(p())->getBox();
        auto lower = hypre_array(patch_box.lower());
        auto upper = hypre_array(patch_box.upper());
        HYPRE_SStructGridSetExtents(d_grid, PART, lower.data(), upper.data());
    }

    std::array<HYPRE_Int, 3> hypre_periodic_shift;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        hypre_periodic_shift[d] = get_hypre_period(*d_level, d_object_name, d);
    }
    for (int d = NDIM; d < 3; ++d)
    {
        hypre_periodic_shift[d] = 0;
    }
    HYPRE_SStructGridSetPeriodic(d_grid, PART, hypre_periodic_shift.data());

#if (NDIM == 2)
    HYPRE_SStructVariable vartypes[NVARS] = { HYPRE_SSTRUCT_VARIABLE_XFACE, HYPRE_SSTRUCT_VARIABLE_YFACE };
#endif
#if (NDIM == 3)
    HYPRE_SStructVariable vartypes[NVARS] = { HYPRE_SSTRUCT_VARIABLE_XFACE,
                                              HYPRE_SSTRUCT_VARIABLE_YFACE,
                                              HYPRE_SSTRUCT_VARIABLE_ZFACE };
#endif
    HYPRE_SStructGridSetVariables(d_grid, PART, NVARS, vartypes);

    HYPRE_SStructGridAssemble(d_grid);

    // Allocate stencil data and set stencil offsets.
    static const int stencil_sz = 2 * NDIM + 1;
    d_stencil_offsets.resize(stencil_sz);
    std::fill(d_stencil_offsets.begin(), d_stencil_offsets.end(), hier::Index<NDIM>(0));
    for (unsigned int axis = 0, stencil_index = 1; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side, ++stencil_index)
        {
            d_stencil_offsets[stencil_index](axis) = (side == 0 ? -1 : +1);
        }
    }
    for (int var = 0; var < NVARS; ++var)
    {
        HYPRE_SStructStencilCreate(NDIM, stencil_sz, &d_stencil[var]);
        for (int s = 0; s < stencil_sz; ++s)
        {
            auto stencil_offset = hypre_array(d_stencil_offsets[s]);
            HYPRE_SStructStencilSetEntry(d_stencil[var], s, stencil_offset.data(), var);
            std::copy(stencil_offset.begin(), stencil_offset.end(), static_cast<int*>(d_stencil_offsets[s]));
        }
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
} // allocateHypreData

void
SCPoissonHypreLevelSolver::setMatrixCoefficients()
{
    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = d_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const auto stencil_size = d_stencil_offsets.size();
        SideData<NDIM, double> matrix_coefs(patch_box, stencil_size, IntVector<NDIM>(0));
        PoissonUtilities::computeMatrixCoefficients(
            matrix_coefs, patch, d_stencil_offsets, d_poisson_spec, d_bc_coefs, d_solution_time);
        // Drop the couplings to sides that are not sides of the level.
        clearOffLevelMatrixEntries(matrix_coefs, *d_level, d_stencil_offsets);

        // Copy matrix entries to the hypre matrix structure.
        std::vector<HYPRE_Int> stencil_indices(stencil_size);
        std::iota(stencil_indices.begin(), stencil_indices.end(), HYPRE_Int(0));
        std::vector<double> mat_vals(stencil_size, 0.0);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
            for (Box<NDIM>::Iterator b(side_box); b; b++)
            {
                SideIndex<NDIM> i(b(), axis, SideIndex<NDIM>::Lower);
                for (unsigned int k = 0; k < stencil_size; ++k)
                {
                    mat_vals[k] = matrix_coefs(i, k);
                }
                // NOTE: In SAMRAI, face-centered values are associated with the
                // cell index located on the "upper" side of the face, but in
                // hypre, face-centered values are associated with the cell
                // index located on the "lower" side of the face. Similarly,
                // in SAMRAI the index stores its axis, but here hypre expects
                // that to be a second argument (so slicing i is okay).
                i(axis) -= 1;
                auto hypre_i = hypre_array(i);
                HYPRE_SStructMatrixSetValues(d_matrix,
                                             PART,
                                             hypre_i.data(),
                                             axis,
                                             stencil_indices.size(),
                                             stencil_indices.data(),
                                             mat_vals.data());
            }
        }
    }

    // Assemble the hypre matrix.
    HYPRE_SStructMatrixAssemble(d_matrix);
    return;
} // setMatrixCoefficients

void
SCPoissonHypreLevelSolver::setupHypreSolver()
{
    // Get the MPI communicator.
    MPI_Comm communicator = IBTK_MPI::getCommunicator();

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
    if (d_solver_type == "PCG" || d_solver_type == "GMRES" || d_solver_type == "FlexGMRES" ||
        d_solver_type == "LGMRES" || d_solver_type == "BiCGSTAB")
    {
        if (d_precond_type == "Split")
        {
            HYPRE_SStructSplitCreate(communicator, &d_precond);
            HYPRE_SStructSplitSetMaxIter(d_precond, 1);
            HYPRE_SStructSplitSetTol(d_precond, 0.0);
            HYPRE_SStructSplitSetZeroGuess(d_precond);
            HYPRE_SStructSplitSetStructSolver(d_precond, split_solver_type_id);
        }
    }

    // Setup the solver.
    if (d_solver_type == "Split")
    {
        HYPRE_SStructSplitCreate(communicator, &d_solver);
        HYPRE_SStructSplitSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructSplitSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructSplitSetStructSolver(d_solver, split_solver_type_id);
        // hypre's zero-guess setting would make the solver treat every iterate
        // after the first as zero.  The initial guess is the solution vector.
        HYPRE_SStructSplitSetNonZeroGuess(d_solver);
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
        if (d_precond_type == "Split")
        {
            HYPRE_SStructPCGSetPrecond(d_solver, HYPRE_SStructSplitSolve, HYPRE_SStructSplitSetup, d_precond);
        }
        else if (d_precond_type == "none")
        {
            d_precond = nullptr;
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
        if (d_precond_type == "Split")
        {
            HYPRE_SStructGMRESSetPrecond(d_solver, HYPRE_SStructSplitSolve, HYPRE_SStructSplitSetup, d_precond);
        }
        else if (d_precond_type == "none")
        {
            d_precond = nullptr;
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
        if (d_precond_type == "Split")
        {
            HYPRE_SStructFlexGMRESSetPrecond(d_solver, HYPRE_SStructSplitSolve, HYPRE_SStructSplitSetup, d_precond);
        }
        else if (d_precond_type == "none")
        {
            d_precond = nullptr;
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
        if (d_precond_type == "Split")
        {
            HYPRE_SStructLGMRESSetPrecond(d_solver, HYPRE_SStructSplitSolve, HYPRE_SStructSplitSetup, d_precond);
        }
        else if (d_precond_type == "none")
        {
            d_precond = nullptr;
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
        if (d_precond_type == "Split")
        {
            HYPRE_SStructBiCGSTABSetPrecond(d_solver, HYPRE_SStructSplitSolve, HYPRE_SStructSplitSetup, d_precond);
        }
        else if (d_precond_type == "none")
        {
            d_precond = nullptr;
        }
        else
        {
            TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                     << "  unknown preconditioner type: " << d_precond_type << std::endl);
        }
        HYPRE_SStructBiCGSTABSetup(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
    }
    else
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  unknown solver type: " << d_solver_type << std::endl);
    }
    return;
} // setupHypreSolver

bool
SCPoissonHypreLevelSolver::solveSystem(const int x_idx, const int b_idx)
{
    const bool level_zero = (d_level_num == 0);

    // Return whether the right-hand-side data of the patch is modified to account for the boundary conditions.
    auto adjusts_rhs = [&](const Pointer<Patch<NDIM>>& patch) -> bool
    {
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const Array<BoundaryBox<NDIM>>& type_1_cf_bdry =
            level_zero ? Array<BoundaryBox<NDIM>>() :
                         d_cf_boundary->getBoundaries(patch->getPatchNumber(), /* boundary type */ 1);
        return pgeom->intersectsPhysicalBoundary() || type_1_cf_bdry.size() > 0;
    };

    // Modify rhs_data, a copy of the right-hand-side data of the patch, to account for any boundary conditions.
    auto adjust_rhs = [&](const Pointer<Patch<NDIM>>& patch, SideData<NDIM, double>& rhs_data) -> void
    {
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        Pointer<SideData<NDIM, double>> x_data = patch->getPatchData(x_idx);
        const Array<BoundaryBox<NDIM>>& type_1_cf_bdry =
            level_zero ? Array<BoundaryBox<NDIM>>() :
                         d_cf_boundary->getBoundaries(patch->getPatchNumber(), /* boundary type */ 1);
        const Array<BoundaryBox<NDIM>>& type_2_cf_bdry =
            level_zero ? Array<BoundaryBox<NDIM>>() :
                         d_cf_boundary->getBoundaries(patch->getPatchNumber(), /* boundary type */ 2);
        if (pgeom->intersectsPhysicalBoundary())
        {
            PoissonUtilities::adjustRHSAtPhysicalBoundary(
                rhs_data, patch, d_poisson_spec, d_bc_coefs, d_solution_time, d_homogeneous_bc);
        }
        if (type_1_cf_bdry.size() > 0)
        {
            PoissonUtilities::adjustRHSAtCoarseFineBoundary(rhs_data,
                                                            *x_data,
                                                            patch,
                                                            d_poisson_spec,
                                                            type_1_cf_bdry,
                                                            type_2_cf_bdry,
                                                            d_bc_coefs,
                                                            d_solution_time,
                                                            d_homogeneous_bc);
        }
    };

    // If the null space contains the constant vector, keep the right-hand side that hypre receives in scratch data and
    // remove its mean in each axis.
    if (d_nullspace_contains_constant_vec)
    {
        allocateNullSpaceData();
        for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = d_level->getPatch(p());
            Pointer<SideData<NDIM, double>> b_data = patch->getPatchData(b_idx);
            Pointer<SideData<NDIM, double>> rhs_data = patch->getPatchData(d_rhs_idx);
            rhs_data->copy(*b_data);
            if (adjusts_rhs(patch))
            {
                adjust_rhs(patch, *rhs_data);
            }
        }
        removeMeanFromRhs();
    }

    // Modify right-hand-side data to account for boundary conditions and copy
    // solution and right-hand-side data to hypre structures.
    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = d_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();

        // The initial guess is zero if initial_guess_nonzero is FALSE, without
        // reading the solution data.  Ghost values are never part of it.
        Pointer<SideData<NDIM, double>> x_data = patch->getPatchData(x_idx);
        if (d_initial_guess_nonzero)
        {
            copyToHypre(d_sol_vec, *x_data, patch_box);
        }
        else
        {
            SideData<NDIM, double> zero_data(patch_box, 1, IntVector<NDIM>(0));
            zero_data.fillAll(0.0);
            copyToHypre(d_sol_vec, zero_data, patch_box);
        }

        // Copy the right-hand-side into the hypre vector.
        Pointer<SideData<NDIM, double>> b_data = patch->getPatchData(b_idx);
        if (d_nullspace_contains_constant_vec)
        {
            Pointer<SideData<NDIM, double>> rhs_data = patch->getPatchData(d_rhs_idx);
            copyToHypre(d_rhs_vec, *rhs_data, patch_box);
        }
        else if (adjusts_rhs(patch))
        {
            SideData<NDIM, double> b_adj_data(patch_box, b_data->getDepth(), IntVector<NDIM>(0));
            b_adj_data.copyOnBox(*b_data, patch_box);
            adjust_rhs(patch, b_adj_data);
            copyToHypre(d_rhs_vec, b_adj_data, patch_box);
        }
        else
        {
            copyToHypre(d_rhs_vec, *b_data, patch_box);
        }
    }

    // Assemble the hypre vectors.
    HYPRE_SStructVectorAssemble(d_sol_vec);
    HYPRE_SStructVectorAssemble(d_rhs_vec);

    // Solve the system.
    IBTK_TIMER_START(t_solve_system_hypre);

    d_current_iterations = 0;
    // HYPRE_INT may be either long or int
    HYPRE_Int current_iterations = d_current_iterations;
    d_current_residual_norm = 0.0;

    if (d_solver_type == "Split")
    {
        HYPRE_SStructSplitSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructSplitSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructSplitSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_SStructSplitGetNumIterations(d_solver, &current_iterations);
        HYPRE_SStructSplitGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "PCG")
    {
        HYPRE_SStructPCGSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructPCGSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructPCGSetAbsoluteTol(d_solver, d_abs_residual_tol);
        HYPRE_SStructPCGSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_SStructPCGGetNumIterations(d_solver, &current_iterations);
        HYPRE_SStructPCGGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "GMRES")
    {
        HYPRE_SStructGMRESSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructGMRESSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructGMRESSetAbsoluteTol(d_solver, d_abs_residual_tol);
        HYPRE_SStructGMRESSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_SStructGMRESGetNumIterations(d_solver, &current_iterations);
        HYPRE_SStructGMRESGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "FlexGMRES")
    {
        HYPRE_SStructFlexGMRESSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructFlexGMRESSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructFlexGMRESSetAbsoluteTol(d_solver, d_abs_residual_tol);
        HYPRE_SStructFlexGMRESSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_SStructFlexGMRESGetNumIterations(d_solver, &current_iterations);
        HYPRE_SStructFlexGMRESGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "LGMRES")
    {
        HYPRE_SStructLGMRESSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructLGMRESSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructLGMRESSetAbsoluteTol(d_solver, d_abs_residual_tol);
        HYPRE_SStructLGMRESSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_SStructLGMRESGetNumIterations(d_solver, &current_iterations);
        HYPRE_SStructLGMRESGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }
    else if (d_solver_type == "BiCGSTAB")
    {
        HYPRE_SStructBiCGSTABSetMaxIter(d_solver, d_max_iterations);
        HYPRE_SStructBiCGSTABSetTol(d_solver, d_rel_residual_tol);
        HYPRE_SStructBiCGSTABSetAbsoluteTol(d_solver, d_abs_residual_tol);
        HYPRE_SStructBiCGSTABSolve(d_solver, d_matrix, d_rhs_vec, d_sol_vec);
        HYPRE_SStructBiCGSTABGetNumIterations(d_solver, &current_iterations);
        HYPRE_SStructBiCGSTABGetFinalRelativeResidualNorm(d_solver, &d_current_residual_norm);
    }

    d_current_iterations = current_iterations;
    IBTK_TIMER_STOP(t_solve_system_hypre);

    // This function reports convergence itself, so do not leave hypre's
    // nonconvergence error set.
    HYPRE_ClearError(HYPRE_ERROR_CONV);

    // Pull the solution vector out of the hypre structures.
    HYPRE_SStructVectorGather(d_sol_vec);
    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = d_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM, double>> x_data = patch->getPatchData(x_idx);
        copyFromHypre(*x_data, d_sol_vec, patch_box);
    }

    // During initialization we may call this function with zero vectors for
    // the RHS and solution - in that case we converge with zero iterations
    // and the relative error is NaN. If this is the case then return true.
    if (std::isnan(d_current_residual_norm) && d_current_iterations == 0) return true;

    return (d_current_residual_norm <= d_rel_residual_tol || d_current_residual_norm <= d_abs_residual_tol);
} // solveSystem

void
SCPoissonHypreLevelSolver::destroyHypreSolver()
{
    // Destroy the solver.
    if (d_solver_type == "Split")
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
    if (d_solver_type == "PCG" || d_solver_type == "GMRES" || d_solver_type == "FlexGMRES" ||
        d_solver_type == "LGMRES" || d_solver_type == "BiCGSTAB")
    {
        if (d_precond_type == "Split")
        {
            HYPRE_SStructSplitDestroy(d_precond);
        }
    }

    // Set the solver and preconditioner pointers to nullptr.
    d_solver = nullptr;
    d_precond = nullptr;
    return;
} // destroyHypreSolver

void
SCPoissonHypreLevelSolver::deallocateHypreData()
{
    if (d_graph) HYPRE_SStructGraphDestroy(d_graph);
    for (const auto& var : d_stencil)
    {
        if (var) HYPRE_SStructStencilDestroy(var);
    }
    if (d_grid) HYPRE_SStructGridDestroy(d_grid);
    if (d_matrix) HYPRE_SStructMatrixDestroy(d_matrix);
    if (d_sol_vec) HYPRE_SStructVectorDestroy(d_sol_vec);
    if (d_rhs_vec) HYPRE_SStructVectorDestroy(d_rhs_vec);
    d_grid = nullptr;
    for (auto& var : d_stencil)
    {
        var = nullptr;
    }
    d_matrix = nullptr;
    d_sol_vec = nullptr;
    d_rhs_vec = nullptr;
    return;
} // deallocateHypreData

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
