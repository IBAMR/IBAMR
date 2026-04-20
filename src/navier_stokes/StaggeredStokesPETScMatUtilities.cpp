// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2022 by the IBAMR developers
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

#include <ibamr/StaggeredStokesPETScMatUtilities.h>

#include <ibtk/ExtendedRobinBcCoefStrategy.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/SideSynchCopyFillPattern.h>
#include <ibtk/compiler_hints.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Array.h>
#include <tbox/MathUtilities.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <petsclog.h>
#include <petscmat.h>

#include <ArrayData.h>
#include <BoundaryBox.h>
#include <Box.h>
#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellGeometry.h>
#include <CellIndex.h>
#include <CoarseFineBoundary.h>
#include <Index.h>
#include <IntVector.h>
#include <MultiblockDataTranslator.h>
#include <Patch.h>
#include <PatchGeometry.h>
#include <PatchLevel.h>
#include <PoissonSpecifications.h>
#include <ProcessorMapping.h>
#include <RefineAlgorithm.h>
#include <RefineOperator.h>
#include <RefineSchedule.h>
#include <RobinBcCoefStrategy.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>
#include <SideIterator.h>
#include <SideVariable.h>
#include <Variable.h>
#include <VariableDatabase.h>
#include <VariableFillPattern.h>

#include <algorithm>
#include <array>
#include <memory>
#include <numeric>
#include <ostream>
#include <set>
#include <sstream>
#include <vector>

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
/*!
 * Return a box extended by one index in the tangential direction associated
 * with side-centered data on `data_axis`.
 */
inline Box<NDIM>
compute_tangential_extension(const Box<NDIM>& box, const int data_axis)
{
    Box<NDIM> extended_box = box;
    extended_box.upper()(data_axis) += 1;
    return extended_box;
} // compute_tangential_extension

std::array<int, NDIM>
get_seed_traversal_axis_order(const CouplingAwareASMSeedTraversalOrder seed_traversal_order)
{
    std::array<int, NDIM> axis_order(array_constant<int, NDIM>(0));
    if (NDIM == 2)
    {
        switch (seed_traversal_order)
        {
        case CouplingAwareASMSeedTraversalOrder::I_J:
            axis_order[0] = 0;
            axis_order[1] = 1;
            return axis_order;
        case CouplingAwareASMSeedTraversalOrder::J_I:
            axis_order[0] = 1;
            axis_order[1] = 0;
            return axis_order;
        case CouplingAwareASMSeedTraversalOrder::I_J_K:
        case CouplingAwareASMSeedTraversalOrder::J_K_I:
        case CouplingAwareASMSeedTraversalOrder::K_I_J:
            TBOX_ERROR("get_seed_traversal_axis_order():\n"
                       << "  3D-only seed traversal order " << enum_to_string(seed_traversal_order)
                       << " is invalid in 2D.\n");
        case CouplingAwareASMSeedTraversalOrder::UNKNOWN:
        default:
            TBOX_ERROR("get_seed_traversal_axis_order():\n"
                       << "  unsupported seed traversal order " << enum_to_string(seed_traversal_order) << ".\n");
        }
    }
    else if (NDIM == 3)
    {
        switch (seed_traversal_order)
        {
        case CouplingAwareASMSeedTraversalOrder::I_J_K:
            axis_order[0] = 0;
            axis_order[1] = 1;
            axis_order[2] = 2;
            return axis_order;
        case CouplingAwareASMSeedTraversalOrder::K_I_J:
            axis_order[0] = 2;
            axis_order[1] = 0;
            axis_order[2] = 1;
            return axis_order;
        case CouplingAwareASMSeedTraversalOrder::J_K_I:
            axis_order[0] = 1;
            axis_order[1] = 2;
            axis_order[2] = 0;
            return axis_order;
        case CouplingAwareASMSeedTraversalOrder::I_J:
        case CouplingAwareASMSeedTraversalOrder::J_I:
            TBOX_ERROR("get_seed_traversal_axis_order():\n"
                       << "  2D-only seed traversal order " << enum_to_string(seed_traversal_order)
                       << " is invalid in 3D.\n");
        case CouplingAwareASMSeedTraversalOrder::UNKNOWN:
        default:
            TBOX_ERROR("get_seed_traversal_axis_order():\n"
                       << "  unsupported seed traversal order " << enum_to_string(seed_traversal_order) << ".\n");
        }
    }
    else
    {
        TBOX_ERROR("get_seed_traversal_axis_order():\n"
                   << "  unsupported NDIM = " << NDIM << ".\n");
    }
    return axis_order;
}

/*!
 * Build STRICT-policy seed-pairing groups for velocity DOFs.
 *
 * Two velocity DOFs are paired when they belong to the same local per-cell
 * grouping used by STRICT closure construction.
 *
 * Output map setup:
 * - keys are velocity DOFs.
 * - each value is the set of paired velocity DOFs on other components.
 * - the map is cleared before insertion.
 */
void
finalize_sorted_unique_map_values(std::unordered_map<int, std::vector<int>>& map_data)
{
    for (auto& pair : map_data)
    {
        std::vector<int>& values = pair.second;
        std::sort(values.begin(), values.end());
        values.erase(std::unique(values.begin(), values.end()), values.end());
    }
}

void
build_coupling_aware_velocity_seed_pair_map(
    std::unordered_map<int, std::vector<int>>& velocity_dof_to_paired_seed_velocity_dofs,
    const int u_dof_index_idx,
    Pointer<PatchLevel<NDIM>> patch_level)
{
    velocity_dof_to_paired_seed_velocity_dofs.clear();
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM, int>> u_dof_data = patch->getPatchData(u_dof_index_idx);

        for (Box<NDIM>::Iterator b(patch_box); b; b++)
        {
            const CellIndex<NDIM>& ic = b();
            std::array<int, NDIM> axis_velocity_dofs(array_constant<int, NDIM>(IBTK::invalid_index));
            for (int axis = 0; axis < NDIM; ++axis)
            {
                axis_velocity_dofs[axis] = (*u_dof_data)(SideIndex<NDIM>(ic, axis, SideIndex<NDIM>::Lower));
            }

            for (int axis = 0; axis < NDIM; ++axis)
            {
                const int dof = axis_velocity_dofs[axis];
                if (dof < 0) continue;
                std::vector<int>& paired_set = velocity_dof_to_paired_seed_velocity_dofs[dof];
                for (int paired_axis = 0; paired_axis < NDIM; ++paired_axis)
                {
                    if (paired_axis == axis) continue;
                    const int paired_dof = axis_velocity_dofs[paired_axis];
                    if (paired_dof >= 0) paired_set.push_back(paired_dof);
                }
            }
        }
    }
    finalize_sorted_unique_map_values(velocity_dof_to_paired_seed_velocity_dofs);
}

/*!
 * Validate that locally owned DOFs are partitioned by `is_nonoverlap`.
 *
 * Input assumptions:
 * - `num_dofs_per_proc` is the global ownership partition for the current run.
 * - `is_nonoverlap` entries use global DOF indices.
 *
 * Throws with detailed diagnostics if any local DOF is missing or duplicated
 * across nonoverlap sets.
 */
void
validate_nonoverlap_subdomains(const std::vector<std::set<int>>& is_nonoverlap,
                               const std::vector<int>& num_dofs_per_proc,
                               const std::string& where)
{
    const int mpi_rank = IBTK_MPI::getRank();
    const int first_local_dof = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int one_past_local_dof = first_local_dof + num_dofs_per_proc[mpi_rank];

    // Each locally owned DOF must appear in exactly one nonoverlap IS.
    std::vector<int> local_nonoverlap_count(static_cast<std::size_t>(one_past_local_dof - first_local_dof), 0);
    for (const auto& subdomain_is : is_nonoverlap)
    {
        for (const int dof : subdomain_is)
        {
            if (dof < first_local_dof || dof >= one_past_local_dof) continue;
            ++local_nonoverlap_count[static_cast<std::size_t>(dof - first_local_dof)];
        }
    }

    int n_missing = 0, n_duplicated = 0;
    std::ostringstream missing_preview, duplicate_preview;
    int missing_preview_count = 0, duplicate_preview_count = 0;
    for (int local_offset = 0; local_offset < static_cast<int>(local_nonoverlap_count.size()); ++local_offset)
    {
        const int count = local_nonoverlap_count[static_cast<std::size_t>(local_offset)];
        const int dof = first_local_dof + local_offset;
        if (count == 0)
        {
            ++n_missing;
            if (missing_preview_count < 10)
            {
                if (missing_preview_count > 0) missing_preview << ", ";
                missing_preview << dof;
                ++missing_preview_count;
            }
        }
        else if (count > 1)
        {
            ++n_duplicated;
            if (duplicate_preview_count < 10)
            {
                if (duplicate_preview_count > 0) duplicate_preview << ", ";
                duplicate_preview << dof << "(x" << count << ")";
                ++duplicate_preview_count;
            }
        }
    }
    if (n_missing > 0 || n_duplicated > 0)
    {
        TBOX_ERROR(where << ":\n"
                         << "  invalid nonoverlap partition on rank " << mpi_rank << ".\n"
                         << "  locally owned DOFs = " << num_dofs_per_proc[mpi_rank] << ", missing = " << n_missing
                         << ", duplicated = " << n_duplicated << ".\n"
                         << "  first missing DOFs: [" << missing_preview.str() << "]\n"
                         << "  first duplicated DOFs: [" << duplicate_preview.str() << "]\n");
    }
}

/*!
 * Construct `is_nonoverlap` by assigning each locally owned DOF to the first
 * overlap subdomain that contains it.
 *
 * Input assumptions:
 * - `is_overlap` entries use global DOF indices.
 * - `num_dofs_per_proc` is the global ownership partition for the current run.
 *
 * Output set setup:
 * - `is_nonoverlap` is resized to `is_overlap.size()`.
 * - each locally owned DOF is inserted into exactly one set (first owner).
 * - resulting partition is validated.
 */
void
construct_nonoverlap_subdomains_from_overlap(const std::vector<std::set<int>>& is_overlap,
                                             std::vector<std::set<int>>& is_nonoverlap,
                                             const std::vector<int>& num_dofs_per_proc)
{
    is_nonoverlap.clear();
    is_nonoverlap.resize(is_overlap.size());

    const int mpi_rank = IBTK_MPI::getRank();
    const int first_local_dof = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int one_past_local_dof = first_local_dof + num_dofs_per_proc[mpi_rank];
    const std::size_t n_local_dofs = static_cast<std::size_t>(one_past_local_dof - first_local_dof);

    const std::size_t invalid_owner = is_overlap.size();
    std::vector<std::size_t> dof_owner(n_local_dofs, invalid_owner);
    for (std::size_t k = 0; k < is_overlap.size(); ++k)
    {
        for (const int dof : is_overlap[k])
        {
            if (dof < first_local_dof || dof >= one_past_local_dof) continue;
            const std::size_t local_offset = static_cast<std::size_t>(dof - first_local_dof);
            if (dof_owner[local_offset] == invalid_owner) dof_owner[local_offset] = k;
        }
    }
    for (std::size_t local_offset = 0; local_offset < n_local_dofs; ++local_offset)
    {
        const std::size_t owner = dof_owner[local_offset];
        if (owner == invalid_owner) continue;
        const int dof = first_local_dof + static_cast<int>(local_offset);
        is_nonoverlap[owner].insert(dof);
    }

    validate_nonoverlap_subdomains(is_nonoverlap, num_dofs_per_proc, "construct_nonoverlap_subdomains_from_overlap()");
}

/*!
 * Compute closure DOFs as the union of precomputed per-cell closures for
 * `involved_cell_dofs`.
 *
 * Input assumptions:
 * - `involved_cell_dofs` are pressure(cell) DOF indices.
 * - `cell_dof_to_closure_dofs` maps pressure(cell) DOFs to full closures.
 *
 * Output set setup:
 * - `closure_dofs` is cleared first.
 * - entries are the union of `cell_dof_to_closure_dofs[cell_dof]` over all
 *   `cell_dof` in `involved_cell_dofs`.
 */
void
find_cell_closure_dofs_from_map(std::set<int>& closure_dofs,
                                const std::set<int>& involved_cell_dofs,
                                const std::unordered_map<int, std::vector<int>>& cell_dof_to_closure_dofs)
{
    closure_dofs.clear();
    for (const int cell_dof : involved_cell_dofs)
    {
        const auto map_it = cell_dof_to_closure_dofs.find(cell_dof);
        if (map_it != cell_dof_to_closure_dofs.end())
        {
            closure_dofs.insert(map_it->second.begin(), map_it->second.end());
        }
    }
}

/*!
 * Validate that `A00_mat` has dimensions consistent with the coupling-aware
 * velocity/cell closure maps.
 */
void
validate_coupling_aware_a00_matrix(Mat A00_mat,
                                   const std::unordered_map<int, std::vector<int>>& velocity_dof_to_adjacent_cell_dofs,
                                   const std::unordered_map<int, std::vector<int>>& cell_dof_to_closure_dofs,
                                   const std::string& where)
{
    PetscInt nrows = -1, ncols = -1;
    int ierr = MatGetSize(A00_mat, &nrows, &ncols);
    IBTK_CHKERRQ(ierr);
    if (nrows != ncols)
    {
        TBOX_ERROR(where << ":\n"
                         << "  A00 block must be square when used for coupling-aware ASM.\n"
                         << "  got size (" << nrows << ", " << ncols << ").\n");
    }

    int local_max_velocity_dof = -1;
    int local_max_dof = -1;
    for (const auto& pair : velocity_dof_to_adjacent_cell_dofs)
    {
        local_max_velocity_dof = std::max(local_max_velocity_dof, pair.first);
        local_max_dof = std::max(local_max_dof, pair.first);
        for (const int cell_dof : pair.second) local_max_dof = std::max(local_max_dof, cell_dof);
    }
    for (const auto& pair : cell_dof_to_closure_dofs)
    {
        local_max_dof = std::max(local_max_dof, pair.first);
        for (const int closure_dof : pair.second) local_max_dof = std::max(local_max_dof, closure_dof);
    }
    const int global_max_velocity_dof = IBTK_MPI::maxReduction(local_max_velocity_dof);
    const int global_max_dof = IBTK_MPI::maxReduction(local_max_dof);
    const PetscInt expected_velocity_size = static_cast<PetscInt>(global_max_velocity_dof + 1);
    const PetscInt expected_a00_size = static_cast<PetscInt>(global_max_dof + 1);
    const bool matches_velocity_size = (expected_velocity_size > 0 && nrows == expected_velocity_size);
    const bool matches_full_size = (expected_a00_size > 0 && nrows == expected_a00_size);
    if (nrows > 0 && !matches_velocity_size && !matches_full_size)
    {
        TBOX_ERROR(where << ":\n"
                         << "  A00 size is inconsistent with the expected matrix spaces.\n"
                         << "  expected either velocity size " << expected_velocity_size << " x "
                         << expected_velocity_size << " (global max velocity DOF index " << global_max_velocity_dof
                         << ") or full size " << expected_a00_size << " x " << expected_a00_size
                         << " (global max DOF index " << global_max_dof << "), but got " << nrows << " x " << ncols
                         << ".\n");
    }
}

/*!
 * Build `initial_velocity_dofs` from seed components and one A00 row-neighbor
 * expansion.
 *
 * Output set setup:
 * - `initial_velocity_dofs` is reset to `seed_velocity_components`.
 * - for each locally owned seed row, row-neighbor velocity DOFs are inserted.
 */
void
build_initial_velocity_set_from_seed_components(std::set<int>& initial_velocity_dofs,
                                                const std::set<int>& seed_velocity_components,
                                                Mat A00_mat,
                                                const PetscInt first_local_row,
                                                const PetscInt one_past_local_row)
{
    initial_velocity_dofs = seed_velocity_components;
    for (const int velocity_dof : seed_velocity_components)
    {
        const PetscInt row = static_cast<PetscInt>(velocity_dof);
        if (row < first_local_row || row >= one_past_local_row) continue;
        PetscInt ncols = 0;
        const PetscInt* cols = nullptr;
        int ierr = MatGetRow(A00_mat, row, &ncols, &cols, nullptr);
        IBTK_CHKERRQ(ierr);
        for (PetscInt col_idx = 0; col_idx < ncols; ++col_idx)
        {
            initial_velocity_dofs.insert(static_cast<int>(cols[col_idx]));
        }
        ierr = MatRestoreRow(A00_mat, row, &ncols, &cols, nullptr);
        IBTK_CHKERRQ(ierr);
    }
}

/*!
 * Map a velocity set to involved pressure(cell) DOFs through the precomputed
 * velocity-to-adjacent-cell map.
 *
 * Output set setup:
 * - `involved_cell_dofs` is cleared first.
 * - entries are the union of adjacent cells for each velocity DOF in
 *   `velocity_dofs`.
 */
void
find_involved_cells_from_velocity_set(
    std::set<int>& involved_cell_dofs,
    const std::set<int>& velocity_dofs,
    const std::unordered_map<int, std::vector<int>>& velocity_dof_to_adjacent_cell_dofs)
{
    involved_cell_dofs.clear();
    for (const int velocity_dof : velocity_dofs)
    {
        const auto map_it = velocity_dof_to_adjacent_cell_dofs.find(velocity_dof);
        if (map_it != velocity_dof_to_adjacent_cell_dofs.end())
        {
            involved_cell_dofs.insert(map_it->second.begin(), map_it->second.end());
        }
    }
}

void
build_coupling_aware_velocity_dof_maps(std::unordered_map<int, std::vector<int>>& velocity_dof_to_adjacent_cell_dofs,
                                       std::unordered_map<int, int>& velocity_dof_to_component_axis,
                                       const int u_dof_index_idx,
                                       const int p_dof_index_idx,
                                       Pointer<PatchLevel<NDIM>> patch_level)
{
    velocity_dof_to_adjacent_cell_dofs.clear();
    velocity_dof_to_component_axis.clear();

    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM, int>> u_dof_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int>> p_dof_data = patch->getPatchData(p_dof_index_idx);

        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (SideIterator<NDIM> s(patch_box, axis); s; s++)
            {
                const SideIndex<NDIM>& i_s = s();
                const int velocity_dof = (*u_dof_data)(i_s);
                if (velocity_dof < 0) continue;

                velocity_dof_to_component_axis[velocity_dof] = axis;
                std::vector<int>& adjacent_cell_dofs = velocity_dof_to_adjacent_cell_dofs[velocity_dof];
                for (int side = 0; side <= 1; ++side)
                {
                    const int cell_dof = (*p_dof_data)(i_s.toCell(side));
                    if (cell_dof >= 0) adjacent_cell_dofs.push_back(cell_dof);
                }
            }
        }
    }
    finalize_sorted_unique_map_values(velocity_dof_to_adjacent_cell_dofs);
}

void
build_cell_dof_to_closure_map_from_velocity_maps(
    std::unordered_map<int, std::vector<int>>& cell_dof_to_closure_dofs,
    const std::unordered_map<int, std::vector<int>>& velocity_dof_to_adjacent_cell_dofs,
    const std::unordered_map<int, int>& velocity_dof_to_component_axis,
    const int u_dof_index_idx,
    const int p_dof_index_idx,
    Pointer<PatchLevel<NDIM>> patch_level)
{
    cell_dof_to_closure_dofs.clear();

    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM, int>> u_dof_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int>> p_dof_data = patch->getPatchData(p_dof_index_idx);

        for (Box<NDIM>::Iterator b(patch_box); b; b++)
        {
            const CellIndex<NDIM>& ic = b();
            const int cell_dof = (*p_dof_data)(ic);
            if (cell_dof < 0) continue;

            std::vector<int>& closure_dofs = cell_dof_to_closure_dofs[cell_dof];
            closure_dofs.push_back(cell_dof);
            auto add_velocity_to_closure = [&](const int velocity_dof)
            {
                if (velocity_dof < 0) return;
#if !defined(NDEBUG)
                const auto axis_it = velocity_dof_to_component_axis.find(velocity_dof);
                if (axis_it == velocity_dof_to_component_axis.end())
                {
                    TBOX_ERROR("build_cell_dof_to_closure_map_from_velocity_maps():\n"
                               << "  missing velocity axis entry for velocity DOF " << velocity_dof << ".\n");
                }
                const auto adjacent_it = velocity_dof_to_adjacent_cell_dofs.find(velocity_dof);
                if (adjacent_it == velocity_dof_to_adjacent_cell_dofs.end())
                {
                    TBOX_ERROR("build_cell_dof_to_closure_map_from_velocity_maps():\n"
                               << "  missing adjacent-cell entry for velocity DOF " << velocity_dof << ".\n");
                }
                if (!std::binary_search(adjacent_it->second.begin(), adjacent_it->second.end(), cell_dof))
                {
                    TBOX_ERROR("build_cell_dof_to_closure_map_from_velocity_maps():\n"
                               << "  inconsistent adjacency: velocity DOF " << velocity_dof
                               << " not adjacent to cell DOF " << cell_dof << ".\n");
                }
#endif
                closure_dofs.push_back(velocity_dof);
            };
            for (int axis = 0; axis < NDIM; ++axis)
            {
                const int lower_velocity_dof = (*u_dof_data)(SideIndex<NDIM>(ic, axis, SideIndex<NDIM>::Lower));
                const int upper_velocity_dof = (*u_dof_data)(SideIndex<NDIM>(ic, axis, SideIndex<NDIM>::Upper));
                add_velocity_to_closure(lower_velocity_dof);
                add_velocity_to_closure(upper_velocity_dof);
            }
        }
    }
    finalize_sorted_unique_map_values(cell_dof_to_closure_dofs);
}
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(
    Mat& mat,
    const PoissonSpecifications& u_problem_coefs,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    double data_time,
    const std::vector<int>& num_dofs_per_proc,
    int u_dof_index_idx,
    int p_dof_index_idx,
    Pointer<PatchLevel<NDIM>> patch_level)
{
    int ierr;
    if (mat)
    {
        ierr = MatDestroy(&mat);
        IBTK_CHKERRQ(ierr);
    }

    // Setup the finite difference stencils.
    static const int uu_stencil_sz = 2 * NDIM + 1;
    std::array<hier::Index<NDIM>, uu_stencil_sz> uu_stencil(
        array_constant<hier::Index<NDIM>, uu_stencil_sz>(hier::Index<NDIM>(0)));
    for (unsigned int axis = 0, uu_stencil_index = 1; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side, ++uu_stencil_index)
        {
            uu_stencil[uu_stencil_index](axis) = (side == 0 ? -1 : +1);
        }
    }
    static const int up_stencil_sz = 2;
    std::array<std::array<hier::Index<NDIM>, up_stencil_sz>, NDIM> up_stencil(
        array_constant<std::array<hier::Index<NDIM>, up_stencil_sz>, NDIM>(
            array_constant<hier::Index<NDIM>, up_stencil_sz>(hier::Index<NDIM>(0))));
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side)
        {
            up_stencil[axis][side](axis) = (side == 0 ? -1 : 0);
        }
    }
    static const int pu_stencil_sz = 2 * NDIM;
    std::array<hier::Index<NDIM>, pu_stencil_sz> pu_stencil(
        array_constant<hier::Index<NDIM>, pu_stencil_sz>(hier::Index<NDIM>(0)));
    for (unsigned int axis = 0, pu_stencil_index = 0; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side, ++pu_stencil_index)
        {
            pu_stencil[pu_stencil_index](axis) = (side == 0 ? 0 : +1);
        }
    }

    // Determine the index ranges.
    const int mpi_rank = IBTK_MPI::getRank();
    const int nlocal = num_dofs_per_proc[mpi_rank];
    const int ilower = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int iupper = ilower + nlocal;
    const int ntotal = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.end(), 0);

    // Determine the non-zero structure of the matrix.
    std::vector<int> d_nnz(nlocal, 0), o_nnz(nlocal, 0);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM, int>> u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int>> p_dof_index_data = patch->getPatchData(p_dof_index_idx);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
            {
                const CellIndex<NDIM>& ic = b();
                const SideIndex<NDIM> is(ic, axis, SideIndex<NDIM>::Lower);
                const int u_dof_index = (*u_dof_index_data)(is);
                if (UNLIKELY(ilower > u_dof_index || u_dof_index >= iupper)) continue;
                const int u_local_idx = u_dof_index - ilower;
                d_nnz[u_local_idx] += 1;
                for (unsigned int d = 0, uu_stencil_index = 1; d < NDIM; ++d)
                {
                    for (int side = 0; side <= 1; ++side, ++uu_stencil_index)
                    {
                        const int uu_dof_index = (*u_dof_index_data)(is + uu_stencil[uu_stencil_index]);
                        if (LIKELY(uu_dof_index >= ilower && uu_dof_index < iupper))
                        {
                            d_nnz[u_local_idx] += 1;
                        }
                        else
                        {
                            o_nnz[u_local_idx] += 1;
                        }
                    }
                }
                for (int side = 0, up_stencil_index = 0; side <= 1; ++side, ++up_stencil_index)
                {
                    const int up_dof_index = (*p_dof_index_data)(ic + up_stencil[axis][up_stencil_index]);
                    if (LIKELY(up_dof_index >= ilower && up_dof_index < iupper))
                    {
                        d_nnz[u_local_idx] += 1;
                    }
                    else
                    {
                        o_nnz[u_local_idx] += 1;
                    }
                }
                d_nnz[u_local_idx] = std::min(nlocal, d_nnz[u_local_idx]);
                o_nnz[u_local_idx] = std::min(ntotal - nlocal, o_nnz[u_local_idx]);
            }
        }
        for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
        {
            const CellIndex<NDIM>& ic = b();
            const int p_dof_index = (*p_dof_index_data)(ic);
            if (UNLIKELY(ilower > p_dof_index || p_dof_index >= iupper)) continue;
            const int p_local_idx = p_dof_index - ilower;
            d_nnz[p_local_idx] += 1;
            for (unsigned int axis = 0, pu_stencil_index = 0; axis < NDIM; ++axis)
            {
                for (int side = 0; side <= 1; ++side, ++pu_stencil_index)
                {
                    const int pu_dof_index = (*u_dof_index_data)(
                        SideIndex<NDIM>(ic + pu_stencil[pu_stencil_index], axis, SideIndex<NDIM>::Lower));
                    if (LIKELY(pu_dof_index >= ilower && pu_dof_index < iupper))
                    {
                        d_nnz[p_local_idx] += 1;
                    }
                    else
                    {
                        o_nnz[p_local_idx] += 1;
                    }
                }
            }
            d_nnz[p_local_idx] = std::min(nlocal, d_nnz[p_local_idx]);
            o_nnz[p_local_idx] = std::min(ntotal - nlocal, o_nnz[p_local_idx]);
        }
    }

    // Create an empty matrix.
    ierr = MatCreateAIJ(
        PETSC_COMM_WORLD, nlocal, nlocal, PETSC_DETERMINE, PETSC_DETERMINE, 0, d_nnz.data(), 0, o_nnz.data(), &mat);
    IBTK_CHKERRQ(ierr);

// Set some general matrix options.
#if !defined(NDEBUG)
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
#endif

    // Set the matrix coefficients.
    const double C = u_problem_coefs.getCConstant();
    const double D = u_problem_coefs.getDConstant();
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();

        const IntVector<NDIM> no_ghosts(0);
        SideData<NDIM, double> uu_matrix_coefs(patch_box, uu_stencil_sz, no_ghosts);
        SideData<NDIM, double> up_matrix_coefs(patch_box, up_stencil_sz, no_ghosts);
        CellData<NDIM, double> pu_matrix_coefs(patch_box, pu_stencil_sz, no_ghosts);

        // Compute all matrix coefficients, including those on the physical
        // boundary; however, do not yet take physical boundary conditions into
        // account.  Boundary conditions are handled subsequently.
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            std::vector<double> uu_mat_vals(uu_stencil_sz, 0.0);
            uu_mat_vals[0] = C; // diagonal
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const double dx_sq = dx[d] * dx[d];
                uu_mat_vals[0] -= 2 * D / dx_sq;    // diagonal
                uu_mat_vals[2 * d + 1] = D / dx_sq; // lower off-diagonal
                uu_mat_vals[2 * d + 2] = D / dx_sq; // upper off-diagonal
            }
            for (int uu_stencil_index = 0; uu_stencil_index < uu_stencil_sz; ++uu_stencil_index)
            {
                uu_matrix_coefs.fill(uu_mat_vals[uu_stencil_index], uu_stencil_index);
            }

            // grad p
            for (int d = 0; d < NDIM; ++d)
            {
                up_matrix_coefs.getArrayData(d).fill(-1.0 / dx[d], 0);
                up_matrix_coefs.getArrayData(d).fill(+1.0 / dx[d], 1);
            }

            // -div u
            std::vector<double> pu_mat_vals(pu_stencil_sz, 0.0);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                pu_matrix_coefs.fill(+1.0 / dx[d], 2 * d);
                pu_matrix_coefs.fill(-1.0 / dx[d], 2 * d + 1);
            }
        }

        // Data structures required to set physical boundary conditions.
        const Array<BoundaryBox<NDIM>> physical_codim1_boxes =
            PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
        const int n_physical_codim1_boxes = physical_codim1_boxes.size();
        const double* const patch_x_lower = pgeom->getXLower();
        const double* const patch_x_upper = pgeom->getXUpper();
        const IntVector<NDIM>& ratio_to_level_zero = pgeom->getRatio();
        Array<Array<bool>> touches_regular_bdry(NDIM), touches_periodic_bdry(NDIM);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            touches_regular_bdry[axis].resizeArray(2);
            touches_periodic_bdry[axis].resizeArray(2);
            for (int upperlower = 0; upperlower < 2; ++upperlower)
            {
                touches_regular_bdry[axis][upperlower] = pgeom->getTouchesRegularBoundary(axis, upperlower);
                touches_periodic_bdry[axis][upperlower] = pgeom->getTouchesPeriodicBoundary(axis, upperlower);
            }
        }

        // Modify matrix coefficients to account for physical boundary
        // conditions along boundaries which ARE NOT aligned with the data axis.
        //
        // NOTE: It important to set these values first to avoid problems at
        // corners in the physical domain.  In particular, since Dirichlet
        // boundary conditions for values located on the physical boundary
        // override all other boundary conditions, we set those values last.
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                const unsigned int location_index = bdry_box.getLocationIndex();
                const unsigned int bdry_normal_axis = location_index / 2;
                const bool is_lower = location_index % 2 == 0;

                if (bdry_normal_axis == axis) continue;

                const Box<NDIM> bc_fill_box =
                    pgeom->getBoundaryFillBox(bdry_box, patch_box, /* ghost_width_to_fill */ IntVector<NDIM>(1));
                const BoundaryBox<NDIM> trimmed_bdry_box =
                    PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
                const Box<NDIM> bc_coef_box = compute_tangential_extension(
                    PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box), axis);

                Pointer<ArrayData<NDIM, double>> acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double>> bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double>> gcoef_data;

                // Temporarily reset the patch geometry object associated with
                // the patch so that boundary conditions are set at the correct
                // spatial locations.
                std::array<double, NDIM> shifted_patch_x_lower, shifted_patch_x_upper;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    shifted_patch_x_lower[d] = patch_x_lower[d];
                    shifted_patch_x_upper[d] = patch_x_upper[d];
                }
                shifted_patch_x_lower[axis] -= 0.5 * dx[axis];
                shifted_patch_x_upper[axis] -= 0.5 * dx[axis];
                patch->setPatchGeometry(new CartesianPatchGeometry<NDIM>(ratio_to_level_zero,
                                                                         touches_regular_bdry,
                                                                         touches_periodic_bdry,
                                                                         dx,
                                                                         shifted_patch_x_lower.data(),
                                                                         shifted_patch_x_upper.data()));

                // Set the boundary condition coefficients.
                static const bool homogeneous_bc = true;
                auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(u_bc_coefs[axis]);
                if (extended_bc_coef)
                {
                    extended_bc_coef->clearTargetPatchDataIndex();
                    extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                }
                u_bc_coefs[axis]->setBcCoefs(
                    acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
                if (gcoef_data && homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

                // Restore the original patch geometry object.
                patch->setPatchGeometry(pgeom);

                // Modify the matrix coefficients to account for homogeneous
                // boundary conditions.
                for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
                {
                    const hier::Index<NDIM>& i = bc();
                    const double& a = (*acoef_data)(i, 0);
                    const double& b = (*bcoef_data)(i, 0);
                    const bool velocity_bc = (a == 1.0 || IBTK::rel_equal_eps(a, 1.0));
                    const bool traction_bc = (b == 1.0 || IBTK::rel_equal_eps(b, 1.0));
#if !defined(NDEBUG)
                    TBOX_ASSERT((velocity_bc || traction_bc) && !(velocity_bc && traction_bc));
#endif
                    hier::Index<NDIM> i_intr = i;
                    if (is_lower)
                    {
                        i_intr(bdry_normal_axis) += 0;
                    }
                    else
                    {
                        i_intr(bdry_normal_axis) -= 1;
                    }
                    const SideIndex<NDIM> i_s(i_intr, axis, SideIndex<NDIM>::Lower);

                    if (velocity_bc)
                    {
                        if (is_lower)
                        {
                            uu_matrix_coefs(i_s, 0) -= uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 1);
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 1) = 0.0;
                        }
                        else
                        {
                            uu_matrix_coefs(i_s, 0) -= uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 2);
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 2) = 0.0;
                        }
                    }
                    else if (traction_bc)
                    {
                        if (is_lower)
                        {
                            uu_matrix_coefs(i_s, 0) += uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 1);
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 1) = 0.0;
                        }
                        else
                        {
                            uu_matrix_coefs(i_s, 0) -= uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 2);
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 2) = 0.0;
                        }
                    }
                    else
                    {
                        TBOX_ERROR(
                            "StaggeredStokesPETScMatUtilities::"
                            "constructPatchLevelMACStokesOp(): Unknown BC type for "
                            "tangential velocity specified.");
                    }
                }
            }
        }

        // Modify matrix coefficients to account for physical boundary
        // conditions along boundaries which ARE aligned with the data axis.
        //
        // NOTE: It important to set these values last to avoid problems at corners
        // in the physical domain.  In particular, since Dirichlet boundary
        // conditions for values located on the physical boundary override all other
        // boundary conditions, we set those values last.
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                const unsigned int location_index = bdry_box.getLocationIndex();
                const unsigned int bdry_normal_axis = location_index / 2;
                const bool is_lower = location_index % 2 == 0;

                if (bdry_normal_axis != axis) continue;

                const Box<NDIM> bc_fill_box =
                    pgeom->getBoundaryFillBox(bdry_box, patch_box, /* ghost_width_to_fill */ IntVector<NDIM>(1));
                const BoundaryBox<NDIM> trimmed_bdry_box =
                    PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
                const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

                Pointer<ArrayData<NDIM, double>> acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double>> bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double>> gcoef_data;

                // Set the boundary condition coefficients.
                static const bool homogeneous_bc = true;
                auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(u_bc_coefs[axis]);
                if (extended_bc_coef)
                {
                    extended_bc_coef->clearTargetPatchDataIndex();
                    extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                }
                u_bc_coefs[axis]->setBcCoefs(
                    acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
                if (gcoef_data && homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

                // Modify the matrix coefficients to account for homogeneous
                // boundary conditions.
                for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
                {
                    const hier::Index<NDIM>& i = bc();
                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                    const double& a = (*acoef_data)(i, 0);
                    const double& b = (*bcoef_data)(i, 0);
                    const bool velocity_bc = (a == 1.0 || IBTK::rel_equal_eps(a, 1.0));
                    const bool traction_bc = (b == 1.0 || IBTK::rel_equal_eps(b, 1.0));
#if !defined(NDEBUG)
                    TBOX_ASSERT((velocity_bc || traction_bc) && !(velocity_bc && traction_bc));
#endif
                    if (velocity_bc)
                    {
                        uu_matrix_coefs(i_s, 0) = 1.0;
                        for (int k = 1; k < uu_stencil_sz; ++k)
                        {
                            uu_matrix_coefs(i_s, k) = 0.0;
                        }
                        for (int k = 0; k < up_stencil_sz; ++k)
                        {
                            up_matrix_coefs(i_s, k) = 0.0;
                        }
                    }
                    else if (traction_bc)
                    {
                        if (is_lower)
                        {
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 2) +=
                                uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 1);
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 1) = 0.0;
                        }
                        else
                        {
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 1) +=
                                uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 2);
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 2) = 0.0;
                        }
                    }
                    else
                    {
                        TBOX_ERROR(
                            "StaggeredStokesPETScMatUtilities::"
                            "constructPatchLevelMACStokesOp(): Unknown BC type for "
                            "normal velocity specified.");
                    }
                }
            }
        }

        // Set matrix coefficients.
        Pointer<SideData<NDIM, int>> u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int>> p_dof_index_data = patch->getPatchData(p_dof_index_idx);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
            {
                const CellIndex<NDIM>& ic = b();
                const SideIndex<NDIM> is(ic, axis, SideIndex<NDIM>::Lower);
                const int u_dof_index = (*u_dof_index_data)(is);
                if (UNLIKELY(ilower > u_dof_index || u_dof_index >= iupper)) continue;

                const int u_stencil_sz = uu_stencil_sz + up_stencil_sz;
                std::vector<double> u_mat_vals(u_stencil_sz);
                std::vector<int> u_mat_cols(u_stencil_sz);

                u_mat_vals[0] = uu_matrix_coefs(is, 0);
                u_mat_cols[0] = u_dof_index;
                for (unsigned int d = 0, uu_stencil_index = 1; d < NDIM; ++d)
                {
                    for (int side = 0; side <= 1; ++side, ++uu_stencil_index)
                    {
                        u_mat_vals[uu_stencil_index] = uu_matrix_coefs(is, uu_stencil_index);
                        u_mat_cols[uu_stencil_index] = (*u_dof_index_data)(is + uu_stencil[uu_stencil_index]);
                    }
                }
                for (int side = 0, up_stencil_index = 0; side <= 1; ++side, ++up_stencil_index)
                {
                    u_mat_vals[uu_stencil_sz + side] = up_matrix_coefs(is, up_stencil_index);
                    u_mat_cols[uu_stencil_sz + side] = (*p_dof_index_data)(ic + up_stencil[axis][up_stencil_index]);
                }

                ierr = MatSetValues(
                    mat, 1, &u_dof_index, u_stencil_sz, u_mat_cols.data(), u_mat_vals.data(), INSERT_VALUES);
                IBTK_CHKERRQ(ierr);
            }
        }

        for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
        {
            const CellIndex<NDIM>& ic = b();
            const int p_dof_index = (*p_dof_index_data)(ic);
            if (UNLIKELY(ilower > p_dof_index || p_dof_index >= iupper)) continue;

            const int p_stencil_sz = pu_stencil_sz + 1;
            std::vector<double> p_mat_vals(p_stencil_sz);
            std::vector<int> p_mat_cols(p_stencil_sz);

            for (unsigned int axis = 0, pu_stencil_index = 0; axis < NDIM; ++axis)
            {
                for (int side = 0; side <= 1; ++side, ++pu_stencil_index)
                {
                    p_mat_vals[pu_stencil_index] = pu_matrix_coefs(ic, pu_stencil_index);
                    p_mat_cols[pu_stencil_index] = (*u_dof_index_data)(
                        SideIndex<NDIM>(ic + pu_stencil[pu_stencil_index], axis, SideIndex<NDIM>::Lower));
                }
            }
            p_mat_vals[pu_stencil_sz] = 0.0;
            p_mat_cols[pu_stencil_sz] = p_dof_index;

            ierr =
                MatSetValues(mat, 1, &p_dof_index, p_stencil_sz, p_mat_cols.data(), p_mat_vals.data(), INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    return;
} // constructPatchLevelMACStokesOp

void
StaggeredStokesPETScMatUtilities::constructPatchLevelGeometricalASMSubdomains(
    std::vector<std::set<int>>& is_overlap,
    std::vector<std::set<int>>& is_nonoverlap,
    const std::vector<int>& num_dofs_per_proc,
    int u_dof_index_idx,
    Pointer<PatchLevel<NDIM>> patch_level,
    Pointer<CoarseFineBoundary<NDIM>> /*cf_boundary*/,
    int p_dof_index_idx,
    const IntVector<NDIM>& box_size,
    const IntVector<NDIM>& overlap_size)
{
    // Clear previously stored index sets.
    for (auto& k : is_overlap)
    {
        k.clear();
    }
    is_overlap.clear();
    for (auto& k : is_nonoverlap)
    {
        k.clear();
    }
    is_nonoverlap.clear();

    // Determine the subdomains associated with this processor.
    const int n_local_patches = patch_level->getProcessorMapping().getNumberOfLocalIndices();
    std::vector<std::vector<Box<NDIM>>> overlap_boxes(n_local_patches), nonoverlap_boxes(n_local_patches);
    int patch_counter = 0, subdomain_counter = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++, ++patch_counter)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        IndexUtilities::partitionPatchBox(
            overlap_boxes[patch_counter], nonoverlap_boxes[patch_counter], patch_box, box_size, overlap_size);
        subdomain_counter += overlap_boxes[patch_counter].size();
    }
    is_overlap.resize(subdomain_counter);
    is_nonoverlap.resize(subdomain_counter);

    // Fill in the IS'es.
    subdomain_counter = 0, patch_counter = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++, ++patch_counter)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        Pointer<SideData<NDIM, int>> u_dof_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int>> p_dof_data = patch->getPatchData(p_dof_index_idx);
#if !defined(NDEBUG)
        {
            const int u_data_depth = u_dof_data->getDepth();
            const int p_data_depth = p_dof_data->getDepth();
            TBOX_ASSERT(u_data_depth == 1);
            TBOX_ASSERT(p_data_depth == 1);
            TBOX_ASSERT(u_dof_data->getGhostCellWidth().min() >= overlap_size.max());
            TBOX_ASSERT(p_dof_data->getGhostCellWidth().min() >= overlap_size.max());
        }
#endif
        int n_patch_subdomains = static_cast<int>(overlap_boxes[patch_counter].size());
        for (int k = 0; k < n_patch_subdomains; ++k, ++subdomain_counter)
        {
            // The overlapping subdomains.
            const Box<NDIM>& overlap_sub_box = overlap_boxes[patch_counter][k];
            Box<NDIM> side_overlap_sub_box[NDIM];
            for (int axis = 0; axis < NDIM; ++axis)
            {
                side_overlap_sub_box[axis] = SideGeometry<NDIM>::toSideBox(overlap_sub_box, axis);
            }

            // Get the overlap DOFs.
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator b(side_overlap_sub_box[axis]); b; b++)
                {
                    const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
                    const int dof_idx = (*u_dof_data)(i_s);
                    if (dof_idx >= 0) is_overlap[subdomain_counter].insert(dof_idx);
                }
            }
            for (Box<NDIM>::Iterator b(overlap_sub_box); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                const int dof_idx = (*p_dof_data)(i);
                if (dof_idx >= 0) is_overlap[subdomain_counter].insert(dof_idx);
            }
        }
    }
    construct_nonoverlap_subdomains_from_overlap(is_overlap, is_nonoverlap, num_dofs_per_proc);
    return;
} // constructPatchLevelGeometricalASMSubdomains

void
StaggeredStokesPETScMatUtilities::ensurePatchLevelVelocityMapDataIsBuilt(PatchLevelCellClosureMapData& map_data,
                                                                         const int u_dof_index_idx,
                                                                         const int p_dof_index_idx,
                                                                         Pointer<PatchLevel<NDIM>> patch_level)
{
    if (map_data.velocity_maps_are_built) return;
    build_coupling_aware_velocity_dof_maps(map_data.velocity_dof_to_adjacent_cell_dofs,
                                           map_data.velocity_dof_to_component_axis,
                                           u_dof_index_idx,
                                           p_dof_index_idx,
                                           patch_level);
    map_data.velocity_maps_are_built = true;
}

void
StaggeredStokesPETScMatUtilities::ensurePatchLevelCellClosureMapIsBuilt(PatchLevelCellClosureMapData& map_data,
                                                                        const int u_dof_index_idx,
                                                                        const int p_dof_index_idx,
                                                                        Pointer<PatchLevel<NDIM>> patch_level)
{
    if (map_data.cell_closure_map_is_built) return;
    ensurePatchLevelVelocityMapDataIsBuilt(map_data, u_dof_index_idx, p_dof_index_idx, patch_level);
    build_cell_dof_to_closure_map_from_velocity_maps(map_data.cell_dof_to_closure_dofs,
                                                     map_data.velocity_dof_to_adjacent_cell_dofs,
                                                     map_data.velocity_dof_to_component_axis,
                                                     u_dof_index_idx,
                                                     p_dof_index_idx,
                                                     patch_level);
    map_data.cell_closure_map_is_built = true;
}

void
StaggeredStokesPETScMatUtilities::ensurePatchLevelVelocitySeedPairMapIsBuilt(PatchLevelCellClosureMapData& map_data,
                                                                             const int u_dof_index_idx,
                                                                             Pointer<PatchLevel<NDIM>> patch_level)
{
    if (map_data.velocity_seed_pair_map_is_built) return;
    build_coupling_aware_velocity_seed_pair_map(
        map_data.velocity_dof_to_paired_seed_velocity_dofs, u_dof_index_idx, patch_level);
    map_data.velocity_seed_pair_map_is_built = true;
}

void
StaggeredStokesPETScMatUtilities::computePatchLevelCouplingAwareASMSeedVelocityDofs(
    std::vector<int>& seed_velocity_dofs,
    const int u_dof_index_idx,
    Pointer<PatchLevel<NDIM>> patch_level,
    const PatchLevelCellClosureMapData& map_data,
    const int seed_velocity_axis,
    const int seed_velocity_stride,
    const CouplingAwareASMSeedTraversalOrder seed_traversal_order)
{
    if (seed_velocity_axis < 0 || seed_velocity_axis >= NDIM)
    {
        TBOX_ERROR("StaggeredStokesPETScMatUtilities::computePatchLevelCouplingAwareASMSeedVelocityDofs():\n"
                   << "  invalid seed_velocity_axis = " << seed_velocity_axis << ".\n");
    }
    if (seed_velocity_stride < 1)
    {
        TBOX_ERROR("StaggeredStokesPETScMatUtilities::computePatchLevelCouplingAwareASMSeedVelocityDofs():\n"
                   << "  invalid seed_velocity_stride = " << seed_velocity_stride << ".\n");
    }
    if (!map_data.velocity_maps_are_built)
    {
        TBOX_ERROR("StaggeredStokesPETScMatUtilities::computePatchLevelCouplingAwareASMSeedVelocityDofs():\n"
                   << "  velocity map data must be built before computing seed velocity DOFs.\n");
    }

    struct SeedVelocityRecord
    {
        int dof = -1;
        std::array<int, NDIM> logical_index{};
    };
    std::vector<SeedVelocityRecord> axis_velocity_records;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Box<NDIM> side_patch_box = SideGeometry<NDIM>::toSideBox(patch_box, seed_velocity_axis);
        Pointer<SideData<NDIM, int>> u_dof_data = patch->getPatchData(u_dof_index_idx);
        for (Box<NDIM>::Iterator b(side_patch_box); b; b++)
        {
            const SideIndex<NDIM> i_s(b(), seed_velocity_axis, SideIndex<NDIM>::Lower);
            const int dof = (*u_dof_data)(i_s);
            if (dof < 0) continue;
            const auto axis_it = map_data.velocity_dof_to_component_axis.find(dof);
            if (axis_it == map_data.velocity_dof_to_component_axis.end()) continue;
            if (axis_it->second != seed_velocity_axis) continue;
            SeedVelocityRecord rec;
            rec.dof = dof;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                rec.logical_index[d] = i_s(static_cast<int>(d));
            }
            axis_velocity_records.push_back(rec);
        }
    }

    const std::array<int, NDIM> axis_order = get_seed_traversal_axis_order(seed_traversal_order);
    std::sort(axis_velocity_records.begin(),
              axis_velocity_records.end(),
              [&axis_order](const SeedVelocityRecord& lhs, const SeedVelocityRecord& rhs)
              {
                  for (unsigned int q = 0; q < NDIM; ++q)
                  {
                      const int d = axis_order[q];
                      if (lhs.logical_index[d] < rhs.logical_index[d]) return true;
                      if (lhs.logical_index[d] > rhs.logical_index[d]) return false;
                  }
                  return lhs.dof < rhs.dof;
              });

    std::vector<int> axis_velocity_dofs;
    axis_velocity_dofs.reserve(axis_velocity_records.size());
    std::set<int> seen_dofs;
    for (const auto& rec : axis_velocity_records)
    {
        if (seen_dofs.insert(rec.dof).second) axis_velocity_dofs.push_back(rec.dof);
    }

    seed_velocity_dofs.clear();
    int axis_velocity_counter = 0;
    for (const int dof : axis_velocity_dofs)
    {
        if ((axis_velocity_counter % seed_velocity_stride) == 0) seed_velocity_dofs.push_back(dof);
        ++axis_velocity_counter;
    }
    return;
}

void
StaggeredStokesPETScMatUtilities::buildPatchLevelCellClosureMaps(PatchLevelCellClosureMapData& map_data,
                                                                 const int u_dof_index_idx,
                                                                 const int p_dof_index_idx,
                                                                 Pointer<PatchLevel<NDIM>> patch_level)
{
    ensurePatchLevelCellClosureMapIsBuilt(map_data, u_dof_index_idx, p_dof_index_idx, patch_level);
    ensurePatchLevelVelocitySeedPairMapIsBuilt(map_data, u_dof_index_idx, patch_level);
    return;
} // buildPatchLevelCellClosureMaps

void
StaggeredStokesPETScMatUtilities::findCoupledCellDofsFromA00(
    std::set<int>& involved_cell_dofs,
    Mat A00_mat,
    const std::set<int>& seed_velocity_dofs,
    const std::unordered_map<int, std::vector<int>>& velocity_dof_to_adjacent_cell_dofs,
    const std::unordered_map<int, std::vector<int>>& cell_dof_to_closure_dofs)
{
    involved_cell_dofs.clear();
    if (!A00_mat || seed_velocity_dofs.empty()) return;

    validate_coupling_aware_a00_matrix(A00_mat,
                                       velocity_dof_to_adjacent_cell_dofs,
                                       cell_dof_to_closure_dofs,
                                       "StaggeredStokesPETScMatUtilities::findCoupledCellDofsFromA00()");

    PetscInt first_local_row = -1, row_end = -1;
    int ierr = MatGetOwnershipRange(A00_mat, &first_local_row, &row_end);
    IBTK_CHKERRQ(ierr);
    std::set<int> initial_velocity_dofs;
    build_initial_velocity_set_from_seed_components(
        initial_velocity_dofs, seed_velocity_dofs, A00_mat, first_local_row, row_end);
    find_involved_cells_from_velocity_set(
        involved_cell_dofs, initial_velocity_dofs, velocity_dof_to_adjacent_cell_dofs);

    return;
} // findCoupledCellDofsFromA00

static void
construct_coupling_aware_asm_overlap_subdomains_with_cell_closure(
    std::vector<std::set<int>>& overlap_is,
    const std::vector<std::set<int>>& nonoverlap_is,
    Mat A00_mat,
    const std::unordered_map<int, std::vector<int>>& velocity_dof_to_adjacent_cell_dofs,
    const std::unordered_map<int, std::vector<int>>& cell_dof_to_closure_dofs,
    const std::unordered_map<int, int>& velocity_dof_to_component_axis,
    const int seed_velocity_axis,
    const CouplingAwareASMClosurePolicy closure_policy,
    const std::unordered_map<int, std::vector<int>>& velocity_dof_to_paired_seed_velocity_dofs)
{
    TBOX_ASSERT(overlap_is.size() == nonoverlap_is.size());
    if (!A00_mat || overlap_is.empty()) return;

    validate_coupling_aware_a00_matrix(A00_mat,
                                       velocity_dof_to_adjacent_cell_dofs,
                                       cell_dof_to_closure_dofs,
                                       "construct_coupling_aware_asm_overlap_subdomains_with_cell_closure()");

    PetscInt first_local_row = -1, row_end = -1;
    int ierr = MatGetOwnershipRange(A00_mat, &first_local_row, &row_end);
    IBTK_CHKERRQ(ierr);

    std::set<int> seed_velocity_dofs, involved_cell_dofs, closure_dofs, initial_velocity_dofs, initial_seed_components;
    for (std::size_t k = 0; k < overlap_is.size(); ++k)
    {
        seed_velocity_dofs.clear();
        involved_cell_dofs.clear();
        closure_dofs.clear();
        initial_velocity_dofs.clear();
        initial_seed_components.clear();

        for (const int dof : nonoverlap_is[k])
        {
            const auto axis_it = velocity_dof_to_component_axis.find(dof);
            if (axis_it != velocity_dof_to_component_axis.end() && axis_it->second == seed_velocity_axis)
            {
                seed_velocity_dofs.insert(dof);
            }
        }

        if (seed_velocity_dofs.empty())
        {
            overlap_is[k] = nonoverlap_is[k];
            continue;
        }

        initial_seed_components = seed_velocity_dofs;
        if (closure_policy == CouplingAwareASMClosurePolicy::STRICT)
        {
            for (const int seed_velocity_dof : seed_velocity_dofs)
            {
                const auto pair_it = velocity_dof_to_paired_seed_velocity_dofs.find(seed_velocity_dof);
                if (pair_it != velocity_dof_to_paired_seed_velocity_dofs.end())
                {
                    initial_seed_components.insert(pair_it->second.begin(), pair_it->second.end());
                }
            }
        }

        build_initial_velocity_set_from_seed_components(
            initial_velocity_dofs, initial_seed_components, A00_mat, first_local_row, row_end);
        find_involved_cells_from_velocity_set(
            involved_cell_dofs, initial_velocity_dofs, velocity_dof_to_adjacent_cell_dofs);
        if (closure_policy == CouplingAwareASMClosurePolicy::STRICT)
        {
            std::set<int> strict_involved_cells;
            for (const int cell_dof : involved_cell_dofs)
            {
                const auto closure_it = cell_dof_to_closure_dofs.find(cell_dof);
                if (closure_it == cell_dof_to_closure_dofs.end()) continue;

                bool valid_cell = true;
                for (const int dof : closure_it->second)
                {
                    if (velocity_dof_to_component_axis.find(dof) == velocity_dof_to_component_axis.end()) continue;
                    if (initial_velocity_dofs.find(dof) == initial_velocity_dofs.end())
                    {
                        valid_cell = false;
                        break;
                    }
                }
                if (valid_cell) strict_involved_cells.insert(cell_dof);
            }
            involved_cell_dofs = strict_involved_cells;
        }
        find_cell_closure_dofs_from_map(closure_dofs, involved_cell_dofs, cell_dof_to_closure_dofs);

        // Coupling-aware overlap mode: start from the nonoverlap subdomain and
        // then add the closure implied by A00 sparsity.
        overlap_is[k] = closure_dofs;
        if (closure_policy == CouplingAwareASMClosurePolicy::RELAXED)
        {
            overlap_is[k].insert(initial_velocity_dofs.begin(), initial_velocity_dofs.end());
        }
    }

    return;
} // construct_coupling_aware_asm_overlap_subdomains_with_cell_closure

void
StaggeredStokesPETScMatUtilities::constructPatchLevelCouplingAwareASMSubdomains(
    std::vector<std::set<int>>& is_overlap,
    std::vector<std::set<int>>& is_nonoverlap,
    const std::vector<int>& num_dofs_per_proc,
    const int u_dof_index_idx,
    Pointer<PatchLevel<NDIM>> patch_level,
    Pointer<CoarseFineBoundary<NDIM>> /*cf_boundary*/,
    Mat A00_mat,
    const PatchLevelCellClosureMapData& map_data,
    const int seed_velocity_axis,
    const int seed_velocity_stride,
    const CouplingAwareASMSeedTraversalOrder seed_traversal_order,
    const CouplingAwareASMClosurePolicy closure_policy)
{
    std::vector<int> seed_velocity_dofs;
    computePatchLevelCouplingAwareASMSeedVelocityDofs(seed_velocity_dofs,
                                                      u_dof_index_idx,
                                                      patch_level,
                                                      map_data,
                                                      seed_velocity_axis,
                                                      seed_velocity_stride,
                                                      seed_traversal_order);

    is_overlap.clear();
    is_nonoverlap.clear();
    is_overlap.resize(seed_velocity_dofs.size());
    is_nonoverlap.resize(seed_velocity_dofs.size());
    std::size_t subdomain_idx = 0;
    for (const int seed_velocity_dof : seed_velocity_dofs)
    {
        is_nonoverlap[subdomain_idx].insert(seed_velocity_dof);
        is_overlap[subdomain_idx] = is_nonoverlap[subdomain_idx];
        ++subdomain_idx;
    }

    construct_coupling_aware_asm_overlap_subdomains_with_cell_closure(
        is_overlap,
        is_nonoverlap,
        A00_mat,
        map_data.velocity_dof_to_adjacent_cell_dofs,
        map_data.cell_dof_to_closure_dofs,
        map_data.velocity_dof_to_component_axis,
        seed_velocity_axis,
        closure_policy,
        map_data.velocity_dof_to_paired_seed_velocity_dofs);

    // Paranoid coverage check: overlap subdomains must cover all locally owned
    // DOFs (equivalently all DOFs in serial).
    const int mpi_rank = IBTK_MPI::getRank();
    const int first_local_dof = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int one_past_local_dof = first_local_dof + num_dofs_per_proc[mpi_rank];
    std::vector<bool> local_dof_covered(static_cast<std::size_t>(one_past_local_dof - first_local_dof), false);
    for (const auto& subdomain_is : is_overlap)
    {
        for (const int dof : subdomain_is)
        {
            if (dof < first_local_dof || dof >= one_past_local_dof) continue;
            local_dof_covered[static_cast<std::size_t>(dof - first_local_dof)] = true;
        }
    }
    int n_missing_local_dofs = 0;
    for (const bool is_covered : local_dof_covered)
    {
        if (!is_covered) ++n_missing_local_dofs;
    }
    if (n_missing_local_dofs > 0)
    {
        std::ostringstream missing_preview;
        int preview_count = 0;
        for (int local_offset = 0; local_offset < static_cast<int>(local_dof_covered.size()) && preview_count < 10;
             ++local_offset)
        {
            if (!local_dof_covered[static_cast<std::size_t>(local_offset)])
            {
                if (preview_count > 0) missing_preview << ", ";
                missing_preview << (first_local_dof + local_offset);
                ++preview_count;
            }
        }
        TBOX_ERROR("StaggeredStokesPETScMatUtilities::constructPatchLevelCouplingAwareASMSubdomains():\n"
                   << "  overlap IS coverage is incomplete on rank " << mpi_rank << ".\n"
                   << "  missing " << n_missing_local_dofs << " of " << num_dofs_per_proc[mpi_rank]
                   << " locally owned DOFs.\n"
                   << "  first missing DOFs: [" << missing_preview.str() << "]\n");
    }

    // In RELAXED mode each closure should include one cell pressure DOF and its
    // incident velocity DOFs: 2*NDIM + 1 unknowns. STRICT mode may produce
    // smaller patches due to filtering and is allowed to do so.
    if (closure_policy == CouplingAwareASMClosurePolicy::RELAXED)
    {
        const std::size_t min_expected_overlap_size = static_cast<std::size_t>(2 * NDIM + 1);
        for (std::size_t k = 0; k < is_overlap.size(); ++k)
        {
            const std::size_t overlap_size = is_overlap[k].size();
            if (overlap_size < min_expected_overlap_size)
            {
                TBOX_ERROR("StaggeredStokesPETScMatUtilities::constructPatchLevelCouplingAwareASMSubdomains():\n"
                           << "  invalid coupling-aware overlap subdomain size for subdomain " << k << ".\n"
                           << "  got overlap size = " << overlap_size << ", expected at least "
                           << min_expected_overlap_size << " (= 2*NDIM + 1).\n");
            }
        }
    }

    construct_nonoverlap_subdomains_from_overlap(is_overlap, is_nonoverlap, num_dofs_per_proc);
    return;
} // constructPatchLevelCouplingAwareASMSubdomains

void
StaggeredStokesPETScMatUtilities::constructA00VelocitySubmatrix(
    Mat& A00_velocity_mat,
    Mat A00_mat,
    const std::vector<int>& num_dofs_per_proc,
    int u_dof_index_idx,
    int p_dof_index_idx,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level)
{
    if (A00_velocity_mat)
    {
        int ierr = MatDestroy(&A00_velocity_mat);
        IBTK_CHKERRQ(ierr);
    }
    if (!A00_mat) return;

    std::vector<std::set<int>> field_is;
    std::vector<std::string> field_names;
    constructPatchLevelFields(field_is, field_names, num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, patch_level);

    auto vel_it = std::find(field_names.begin(), field_names.end(), "velocity");
    if (vel_it == field_names.end())
    {
        TBOX_ERROR("StaggeredStokesPETScMatUtilities::constructA00VelocitySubmatrix():\n"
                   << "  unable to locate velocity field indices.\n");
    }
    const std::size_t vel_idx = static_cast<std::size_t>(std::distance(field_names.begin(), vel_it));
    std::vector<PetscInt> velocity_dofs(field_is[vel_idx].begin(), field_is[vel_idx].end());
    IS velocity_is = nullptr;
    int ierr = ISCreateGeneral(PETSC_COMM_WORLD,
                               static_cast<PetscInt>(velocity_dofs.size()),
                               velocity_dofs.empty() ? nullptr : velocity_dofs.data(),
                               PETSC_COPY_VALUES,
                               &velocity_is);
    IBTK_CHKERRQ(ierr);
    ierr = MatCreateSubMatrix(A00_mat, velocity_is, velocity_is, MAT_INITIAL_MATRIX, &A00_velocity_mat);
    IBTK_CHKERRQ(ierr);
    ierr = ISDestroy(&velocity_is);
    IBTK_CHKERRQ(ierr);
    return;
} // constructA00VelocitySubmatrix

void
StaggeredStokesPETScMatUtilities::constructPatchLevelFields(
    std::vector<std::set<int>>& is_field,
    std::vector<std::string>& is_field_name,
    const std::vector<int>& num_dofs_per_proc,
    int u_dof_index_idx,
    int p_dof_index_idx,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level)
{
    // Destroy the previously stored IS'es
    for (auto& k : is_field)
    {
        k.clear();
    }
    is_field.clear();
    is_field_name.clear();

    // Resize vectors
    is_field.resize(2);
    is_field_name.resize(2);

    // Name of the fields.
    static const int U_FIELD_IDX = 0;
    static const int P_FIELD_IDX = 1;
    is_field_name[U_FIELD_IDX] = "velocity";
    is_field_name[P_FIELD_IDX] = "pressure";

    // DOFs on this processor.
    const int mpi_rank = IBTK_MPI::getRank();
    const int n_local_dofs = num_dofs_per_proc[mpi_rank];

    const int first_local_dof = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int last_local_dof = first_local_dof + n_local_dofs;

    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Box<NDIM> side_patch_box[NDIM];
        for (int axis = 0; axis < NDIM; ++axis)
        {
            side_patch_box[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
        }

        Pointer<SideData<NDIM, int>> u_dof_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int>> p_dof_data = patch->getPatchData(p_dof_index_idx);
#if !defined(NDEBUG)
        const int u_data_depth = u_dof_data->getDepth();
        const int p_data_depth = p_dof_data->getDepth();
        TBOX_ASSERT(u_data_depth == 1);
        TBOX_ASSERT(p_data_depth == 1);
#endif

        // Get the local velocity DOFs.
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(side_patch_box[axis]); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                const int dof_idx = (*u_dof_data)(i_s);
                if (dof_idx >= first_local_dof && dof_idx < last_local_dof)
                {
                    is_field[0].insert(dof_idx);
                }
            }
        }

        // Get the local pressure DOFs.
        for (Box<NDIM>::Iterator b(patch_box); b; b++)
        {
            const CellIndex<NDIM>& i = b();
            const int dof_idx = (*p_dof_data)(i);
            if (dof_idx >= first_local_dof && dof_idx < last_local_dof)
            {
                is_field[1].insert(dof_idx);
            }
        }
    }

    return;
} // constructPatchLevelFields

void
StaggeredStokesPETScMatUtilities::constructProlongationOp(Mat& mat,
                                                          const std::string& u_op_type,
                                                          const std::string& p_op_type,
                                                          int u_dof_index_idx,
                                                          int p_dof_index_idx,
                                                          const std::vector<int>& num_fine_dofs_per_proc,
                                                          const std::vector<int>& num_coarse_dofs_per_proc,
                                                          Pointer<PatchLevel<NDIM>> fine_patch_level,
                                                          Pointer<PatchLevel<NDIM>> coarse_patch_level,
                                                          const AO& coarse_level_ao,
                                                          const int u_coarse_ao_offset,
                                                          const int p_coarse_ao_offset)
{
    int ierr;
    Mat p_prolong_mat = nullptr;
    PETScMatUtilities::constructProlongationOp(mat,
                                               u_op_type,
                                               u_dof_index_idx,
                                               num_fine_dofs_per_proc,
                                               num_coarse_dofs_per_proc,
                                               fine_patch_level,
                                               coarse_patch_level,
                                               coarse_level_ao,
                                               u_coarse_ao_offset);

    PETScMatUtilities::constructProlongationOp(p_prolong_mat,
                                               p_op_type,
                                               p_dof_index_idx,
                                               num_fine_dofs_per_proc,
                                               num_coarse_dofs_per_proc,
                                               fine_patch_level,
                                               coarse_patch_level,
                                               coarse_level_ao,
                                               p_coarse_ao_offset);

    // P{u,p} = (P_u + P_p){u,p}
    ierr = MatAXPY(mat, 1.0, p_prolong_mat, DIFFERENT_NONZERO_PATTERN);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&p_prolong_mat);
    IBTK_CHKERRQ(ierr);

} // constructPatchLevelProlongationOp

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
