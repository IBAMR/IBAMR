// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibamr/StaggeredStokesPETScMatUtilities.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>
#include <ibamr/private/StaggeredStokesPETScMatUtilities-inl.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/ibtk_utilities.h>

#include <petscmat.h>

#include <CellData.h>
#include <CellIndex.h>
#include <CellVariable.h>
#include <CoarseFineBoundary.h>
#include <Patch.h>
#include <PoissonSpecifications.h>
#include <SideVariable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <array>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
std::string
set_to_string(const std::set<int>& values)
{
    std::ostringstream stream;
    stream << "{";
    bool first = true;
    for (const int value : values)
    {
        if (!first) stream << ", ";
        stream << value;
        first = false;
    }
    stream << "}";
    return stream.str();
}

bool
check_set(const std::string& label, const std::set<int>& expected, const std::set<int>& actual)
{
    if (expected == actual) return true;
    std::cerr << "FAILED: " << label << "\n"
              << "  expected: " << set_to_string(expected) << "\n"
              << "  actual:   " << set_to_string(actual) << std::endl;
    return false;
}

std::set<int>
set_from_vector(const std::vector<int>& values)
{
    return std::set<int>(values.begin(), values.end());
}

bool
check_set(const std::string& label, const std::set<int>& expected, const std::vector<int>& actual)
{
    return check_set(label, expected, set_from_vector(actual));
}

std::set<int>
build_initial_velocity_dofs(const std::set<int>& seed_velocity_dofs,
                            Mat A00_mat,
                            const double relative_numerical_zero_tol)
{
    std::set<int> initial_velocity_dofs = seed_velocity_dofs;
    PetscInt first_local_row = -1;
    PetscInt row_end = -1;
    int ierr = MatGetOwnershipRange(A00_mat, &first_local_row, &row_end);
    IBTK_CHKERRQ(ierr);
    for (const int velocity_dof : seed_velocity_dofs)
    {
        const PetscInt row = static_cast<PetscInt>(velocity_dof);
        if (row < first_local_row || row >= row_end) continue;
        PetscInt ncols = 0;
        const PetscInt* cols = nullptr;
        const PetscScalar* vals = nullptr;
        ierr = MatGetRow(A00_mat, row, &ncols, &cols, &vals);
        IBTK_CHKERRQ(ierr);
        double row_max_abs = 0.0;
        for (PetscInt k = 0; k < ncols; ++k)
        {
            row_max_abs = std::max(row_max_abs, static_cast<double>(PetscAbsScalar(vals[k])));
        }
        const double numerical_zero_tol =
            std::max(static_cast<double>(ncols) * std::numeric_limits<double>::epsilon() * row_max_abs,
                     relative_numerical_zero_tol * row_max_abs);
        for (PetscInt k = 0; k < ncols; ++k)
        {
            const double value_abs = static_cast<double>(PetscAbsScalar(vals[k]));
            if (value_abs > numerical_zero_tol) initial_velocity_dofs.insert(static_cast<int>(cols[k]));
        }
        ierr = MatRestoreRow(A00_mat, row, &ncols, &cols, &vals);
        IBTK_CHKERRQ(ierr);
    }
    return initial_velocity_dofs;
}

std::set<int>
reference_extract_coupled_dofs_relaxed(
    const int seed_velocity_dof,
    Mat A00_mat,
    const std::unordered_map<int, std::vector<int>>& velocity_dof_to_adjacent_cell_dofs,
    const std::unordered_map<int, std::vector<int>>& cell_dof_to_closure_dofs,
    const double relative_numerical_zero_tol)
{
    const std::set<int> initial_velocity_dofs =
        build_initial_velocity_dofs(std::set<int>{ seed_velocity_dof }, A00_mat, relative_numerical_zero_tol);

    std::set<int> pressure_dofs;
    for (const int velocity_dof : initial_velocity_dofs)
    {
        const auto velocity_it = velocity_dof_to_adjacent_cell_dofs.find(velocity_dof);
        if (velocity_it != velocity_dof_to_adjacent_cell_dofs.end())
        {
            pressure_dofs.insert(velocity_it->second.begin(), velocity_it->second.end());
        }
    }

    std::set<int> coupled_dofs = initial_velocity_dofs;
    for (const int pressure_dof : pressure_dofs)
    {
        const auto closure_it = cell_dof_to_closure_dofs.find(pressure_dof);
        if (closure_it == cell_dof_to_closure_dofs.end()) continue;
        coupled_dofs.insert(closure_it->second.begin(), closure_it->second.end());
    }
    return coupled_dofs;
}

std::set<int>
reference_extract_coupled_dofs_strict(
    const int seed_velocity_dof,
    Mat A00_mat,
    const std::unordered_map<int, std::vector<int>>& velocity_dof_to_adjacent_cell_dofs,
    const std::unordered_map<int, std::vector<int>>& cell_dof_to_closure_dofs,
    const std::unordered_map<int, int>& velocity_dof_to_component_axis,
    const std::unordered_map<int, std::vector<int>>& velocity_dof_to_paired_seed_velocity_dofs,
    const double relative_numerical_zero_tol)
{
    std::set<int> seed_velocity_dofs = { seed_velocity_dof };
    const auto pair_it = velocity_dof_to_paired_seed_velocity_dofs.find(seed_velocity_dof);
    if (pair_it != velocity_dof_to_paired_seed_velocity_dofs.end())
    {
        seed_velocity_dofs.insert(pair_it->second.begin(), pair_it->second.end());
    }

    const std::set<int> initial_velocity_dofs =
        build_initial_velocity_dofs(seed_velocity_dofs, A00_mat, relative_numerical_zero_tol);

    std::set<int> candidate_pressure_dofs;
    for (const int velocity_dof : initial_velocity_dofs)
    {
        const auto velocity_it = velocity_dof_to_adjacent_cell_dofs.find(velocity_dof);
        if (velocity_it != velocity_dof_to_adjacent_cell_dofs.end())
        {
            candidate_pressure_dofs.insert(velocity_it->second.begin(), velocity_it->second.end());
        }
    }

    std::set<int> coupled_dofs;
    for (const int pressure_dof : candidate_pressure_dofs)
    {
        const auto closure_it = cell_dof_to_closure_dofs.find(pressure_dof);
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
        if (valid_cell) coupled_dofs.insert(closure_it->second.begin(), closure_it->second.end());
    }

    return coupled_dofs;
}

std::array<int, NDIM>
get_axis_order(const IBAMR::CouplingAwareASMSeedTraversalOrder traversal_order)
{
#if (NDIM == 2)
    if (traversal_order == IBAMR::CouplingAwareASMSeedTraversalOrder::I_J) return { 0, 1 };
    if (traversal_order == IBAMR::CouplingAwareASMSeedTraversalOrder::J_I) return { 1, 0 };
#else
    if (traversal_order == IBAMR::CouplingAwareASMSeedTraversalOrder::I_J_K) return { 0, 1, 2 };
    if (traversal_order == IBAMR::CouplingAwareASMSeedTraversalOrder::J_K_I) return { 1, 2, 0 };
    if (traversal_order == IBAMR::CouplingAwareASMSeedTraversalOrder::K_I_J) return { 2, 0, 1 };
#endif
    TBOX_ERROR("stokes_ib_cav_reference_parity: unsupported traversal order.\n");
    return {};
}

std::vector<int>
get_logically_ordered_pressure_dofs(const int p_dof_index_idx,
                                    Pointer<PatchLevel<NDIM>> level,
                                    const IBAMR::CouplingAwareASMSeedTraversalOrder traversal_order)
{
    struct PressureRecord
    {
        int dof = -1;
        std::array<int, NDIM> logical_index{};
    };
    std::vector<PressureRecord> records;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<CellData<NDIM, int>> p_dof_data = patch->getPatchData(p_dof_index_idx);
        for (Box<NDIM>::Iterator b(patch->getBox()); b; b++)
        {
            const CellIndex<NDIM>& cell = b();
            const int dof = (*p_dof_data)(cell);
            if (dof < 0) continue;
            PressureRecord record;
            record.dof = dof;
            for (unsigned int d = 0; d < NDIM; ++d) record.logical_index[d] = cell(static_cast<int>(d));
            records.push_back(record);
        }
    }
    const std::array<int, NDIM> axis_order = get_axis_order(traversal_order);
    std::sort(records.begin(),
              records.end(),
              [&axis_order](const PressureRecord& lhs, const PressureRecord& rhs)
              {
                  for (unsigned int q = 0; q < NDIM; ++q)
                  {
                      const int d = axis_order[q];
                      if (lhs.logical_index[d] < rhs.logical_index[d]) return true;
                      if (lhs.logical_index[d] > rhs.logical_index[d]) return false;
                  }
                  return lhs.dof < rhs.dof;
              });
    std::vector<int> pressure_dofs;
    std::set<int> seen;
    for (const PressureRecord& record : records)
    {
        if (seen.insert(record.dof).second) pressure_dofs.push_back(record.dof);
    }
    return pressure_dofs;
}

Mat
make_eulerian_elasticity_pattern(const int n_global_dofs, const std::vector<std::tuple<int, int, double>>& entries = {})
{
    Mat mat = nullptr;
    int ierr = MatCreateAIJ(PETSC_COMM_WORLD,
                            PETSC_DECIDE,
                            PETSC_DECIDE,
                            n_global_dofs,
                            n_global_dofs,
                            static_cast<PetscInt>(std::max<std::size_t>(1, entries.size())),
                            nullptr,
                            static_cast<PetscInt>(std::max<std::size_t>(1, entries.size())),
                            nullptr,
                            &mat);
    IBTK_CHKERRQ(ierr);
    for (const auto& entry : entries)
    {
        ierr = MatSetValue(mat, std::get<0>(entry), std::get<1>(entry), std::get<2>(entry), INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    return mat;
}

std::set<int>
reference_pressure_seeded_patch(const int pressure_seed_dof,
                                const std::set<int>& expanded_velocity_dofs,
                                const std::unordered_map<int, std::vector<int>>& velocity_dof_to_adjacent_cell_dofs,
                                const std::unordered_map<int, std::vector<int>>& cell_dof_to_closure_dofs,
                                const std::unordered_map<int, int>& velocity_dof_to_component_axis,
                                const IBAMR::CouplingAwareASMClosurePolicy closure_policy)
{
    std::set<int> candidate_pressure_dofs = { pressure_seed_dof };
    for (const int velocity_dof : expanded_velocity_dofs)
    {
        const auto velocity_it = velocity_dof_to_adjacent_cell_dofs.find(velocity_dof);
        if (velocity_it != velocity_dof_to_adjacent_cell_dofs.end())
        {
            candidate_pressure_dofs.insert(velocity_it->second.begin(), velocity_it->second.end());
        }
    }

    std::set<int> patch_dofs;
    for (const int pressure_dof : candidate_pressure_dofs)
    {
        const auto closure_it = cell_dof_to_closure_dofs.find(pressure_dof);
        if (closure_it == cell_dof_to_closure_dofs.end()) continue;
        if (closure_policy == IBAMR::CouplingAwareASMClosurePolicy::STRICT)
        {
            bool complete_velocity_stencil = true;
            for (const int dof : closure_it->second)
            {
                if (velocity_dof_to_component_axis.find(dof) == velocity_dof_to_component_axis.end()) continue;
                if (expanded_velocity_dofs.find(dof) == expanded_velocity_dofs.end())
                {
                    complete_velocity_stencil = false;
                    break;
                }
            }
            if (!complete_velocity_stencil) continue;
        }
        patch_dofs.insert(closure_it->second.begin(), closure_it->second.end());
    }
    if (closure_policy == IBAMR::CouplingAwareASMClosurePolicy::RELAXED)
    {
        patch_dofs.insert(expanded_velocity_dofs.begin(), expanded_velocity_dofs.end());
    }
    return patch_dofs;
}

bool
check_ordered_set(const std::string& label, const std::set<int>& expected, const std::vector<int>& actual)
{
    if (!std::is_sorted(actual.begin(), actual.end()))
    {
        std::cerr << "FAILED: " << label << " is not in deterministic increasing global-DOF order." << std::endl;
        return false;
    }
    return check_set(label, expected, actual);
}

} // namespace

int
main(int argc, char* argv[])
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app_initializer->getInputDatabase();

    int test_failures = 0;
    const auto hierarchy_tuple = setup_hierarchy<NDIM>(app_initializer);
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = std::get<0>(hierarchy_tuple);
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(0);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("stokes_ib_cav_reference_parity_ctx");
    Pointer<SideVariable<NDIM, int>> u_dof_index_var = new SideVariable<NDIM, int>("reference_parity_u_dof");
    Pointer<CellVariable<NDIM, int>> p_dof_index_var = new CellVariable<NDIM, int>("reference_parity_p_dof");
    const int u_dof_index_idx = var_db->registerVariableAndContext(u_dof_index_var, ctx, IntVector<NDIM>(1));
    const int p_dof_index_idx = var_db->registerVariableAndContext(p_dof_index_var, ctx, IntVector<NDIM>(1));

    level->allocatePatchData(u_dof_index_idx);
    level->allocatePatchData(p_dof_index_idx);

    std::vector<int> num_dofs_per_proc;
    IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
        num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, level);

    IBAMR::StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData map_data;
    IBAMR::StaggeredStokesPETScMatUtilities::buildPatchLevelCellClosureMaps(
        map_data, u_dof_index_idx, p_dof_index_idx, level);
    IBAMR::StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData pressure_map_data;
    IBAMR::StaggeredStokesPETScMatUtilities::buildPatchLevelCellClosureMaps(
        pressure_map_data, u_dof_index_idx, p_dof_index_idx, level);
    using MapAccess = IBAMR::StaggeredStokesPETScMatUtilitiesPrivateAccess;
    const auto& velocity_dof_to_adjacent_cell_dofs = MapAccess::getVelocityDofToAdjacentCellDofs(map_data);
    const auto& cell_dof_to_closure_dofs = MapAccess::getCellDofToClosureDofs(map_data);
    const auto& velocity_dof_to_component_axis = MapAccess::getVelocityDofToComponentAxis(map_data);
    bool velocity_relaxed_seed_pair_map_absent = !MapAccess::velocitySeedPairMapIsBuilt(map_data) &&
                                                 MapAccess::getVelocityDofToPairedSeedVelocityDofs(map_data).empty();

    std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);
    PoissonSpecifications u_problem_coefs("reference_parity_poisson");
    u_problem_coefs.setCConstant(1.0);
    u_problem_coefs.setDConstant(-1.0);
    const double relative_numerical_zero_tol = IBTK_RELATIVE_NUMERICAL_ZERO_TOL;
    Mat level_mat = nullptr;
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(
        level_mat, u_problem_coefs, u_bc_coefs, 0.0, num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, level);
    int ierr = 0;

#if (NDIM == 2)
    const auto primary_traversal_order = IBAMR::CouplingAwareASMSeedTraversalOrder::I_J;
    const auto alternate_traversal_order = IBAMR::CouplingAwareASMSeedTraversalOrder::J_I;
#else
    const auto primary_traversal_order = IBAMR::CouplingAwareASMSeedTraversalOrder::I_J_K;
    const auto alternate_traversal_order = IBAMR::CouplingAwareASMSeedTraversalOrder::K_I_J;
#endif

    std::vector<int> seed_velocity_dofs;
    IBAMR::StaggeredStokesPETScMatUtilitiesPrivateAccess::computePatchLevelCouplingAwareASMSeedVelocityDofs(
        seed_velocity_dofs, u_dof_index_idx, level, map_data, 0, 1, primary_traversal_order);

    std::vector<std::vector<int>> overlap_relaxed, nonoverlap_relaxed;
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelCouplingAwareASMSubdomains(
        overlap_relaxed,
        nonoverlap_relaxed,
        num_dofs_per_proc,
        u_dof_index_idx,
        level,
        Pointer<CoarseFineBoundary<NDIM>>(nullptr),
        level_mat,
        map_data,
        0,
        1,
        primary_traversal_order,
        IBAMR::CouplingAwareASMClosurePolicy::RELAXED,
        relative_numerical_zero_tol);
    velocity_relaxed_seed_pair_map_absent = velocity_relaxed_seed_pair_map_absent &&
                                            !MapAccess::velocitySeedPairMapIsBuilt(map_data) &&
                                            MapAccess::getVelocityDofToPairedSeedVelocityDofs(map_data).empty();
    if (!velocity_relaxed_seed_pair_map_absent) ++test_failures;

    std::vector<std::vector<int>> overlap_strict, nonoverlap_strict;
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelCouplingAwareASMSubdomains(
        overlap_strict,
        nonoverlap_strict,
        num_dofs_per_proc,
        u_dof_index_idx,
        level,
        Pointer<CoarseFineBoundary<NDIM>>(nullptr),
        level_mat,
        map_data,
        0,
        1,
        primary_traversal_order,
        IBAMR::CouplingAwareASMClosurePolicy::STRICT,
        relative_numerical_zero_tol);

    const auto& velocity_dof_to_paired_seed_velocity_dofs = MapAccess::getVelocityDofToPairedSeedVelocityDofs(map_data);
    std::size_t velocity_seed_pair_value_count = 0;
    bool velocity_strict_seed_pair_map_count_valid =
        MapAccess::velocitySeedPairMapIsBuilt(map_data) &&
        velocity_dof_to_paired_seed_velocity_dofs.size() == velocity_dof_to_component_axis.size();
    for (const auto& pair : velocity_dof_to_paired_seed_velocity_dofs)
    {
        velocity_seed_pair_value_count += pair.second.size();
        velocity_strict_seed_pair_map_count_valid =
            velocity_strict_seed_pair_map_count_valid && pair.second.size() == NDIM - 1;
    }
    if (!velocity_strict_seed_pair_map_count_valid) ++test_failures;

    if (overlap_relaxed.empty() || overlap_strict.empty()) ++test_failures;
    if (seed_velocity_dofs.empty()) ++test_failures;
    if (overlap_relaxed.size() != nonoverlap_relaxed.size()) ++test_failures;
    if (overlap_strict.size() != nonoverlap_strict.size()) ++test_failures;
    if (overlap_relaxed.size() != seed_velocity_dofs.size()) ++test_failures;
    if (overlap_strict.size() != seed_velocity_dofs.size()) ++test_failures;
    if (nonoverlap_relaxed.size() != seed_velocity_dofs.size()) ++test_failures;
    if (nonoverlap_strict.size() != seed_velocity_dofs.size()) ++test_failures;

    for (std::size_t k = 0; k < seed_velocity_dofs.size(); ++k)
    {
        const int seed_velocity_dof = seed_velocity_dofs[k];

        const std::set<int> expected_relaxed =
            reference_extract_coupled_dofs_relaxed(seed_velocity_dof,
                                                   level_mat,
                                                   velocity_dof_to_adjacent_cell_dofs,
                                                   cell_dof_to_closure_dofs,
                                                   relative_numerical_zero_tol);
        const std::set<int> expected_strict =
            reference_extract_coupled_dofs_strict(seed_velocity_dof,
                                                  level_mat,
                                                  velocity_dof_to_adjacent_cell_dofs,
                                                  cell_dof_to_closure_dofs,
                                                  velocity_dof_to_component_axis,
                                                  velocity_dof_to_paired_seed_velocity_dofs,
                                                  relative_numerical_zero_tol);

        if (!check_set("live relaxed CAV overlap for seed " + std::to_string(seed_velocity_dof),
                       expected_relaxed,
                       overlap_relaxed[k]))
        {
            ++test_failures;
        }
        if (!check_set("live strict CAV overlap for seed " + std::to_string(seed_velocity_dof),
                       expected_strict,
                       overlap_strict[k]))
        {
            ++test_failures;
        }
        if (!std::includes(
                expected_relaxed.begin(), expected_relaxed.end(), expected_strict.begin(), expected_strict.end()))
        {
            std::cerr
                << "FAILED: expected strict reference overlap to be contained in relaxed reference overlap for seed "
                << seed_velocity_dof << "." << std::endl;
            ++test_failures;
        }
    }

    const std::vector<int> primary_pressure_seed_dofs =
        get_logically_ordered_pressure_dofs(p_dof_index_idx, level, primary_traversal_order);
    const std::vector<int> alternate_pressure_seed_dofs =
        get_logically_ordered_pressure_dofs(p_dof_index_idx, level, alternate_traversal_order);
    const int n_global_dofs = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.end(), 0);

    Mat zero_eulerian_elasticity_mat = make_eulerian_elasticity_pattern(n_global_dofs);
    std::vector<std::vector<int>> pressure_patches_zero_relaxed;
    std::vector<int> constructed_pressure_seed_dofs;
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelPressureCellSeededCAVPatches(
        pressure_patches_zero_relaxed,
        constructed_pressure_seed_dofs,
        num_dofs_per_proc,
        p_dof_index_idx,
        level,
        zero_eulerian_elasticity_mat,
        pressure_map_data,
        1,
        primary_traversal_order,
        IBAMR::CouplingAwareASMClosurePolicy::RELAXED,
        relative_numerical_zero_tol);
    if (constructed_pressure_seed_dofs != primary_pressure_seed_dofs)
    {
        std::cerr << "FAILED: pressure-seed order does not equal the independent primary logical order." << std::endl;
        ++test_failures;
    }

    std::vector<std::vector<int>> pressure_patches_zero_strict;
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelPressureCellSeededCAVPatches(
        pressure_patches_zero_strict,
        constructed_pressure_seed_dofs,
        num_dofs_per_proc,
        p_dof_index_idx,
        level,
        zero_eulerian_elasticity_mat,
        pressure_map_data,
        1,
        primary_traversal_order,
        IBAMR::CouplingAwareASMClosurePolicy::STRICT,
        relative_numerical_zero_tol);
    if (constructed_pressure_seed_dofs != primary_pressure_seed_dofs)
    {
        std::cerr << "FAILED: STRICT pressure-seed order differs from the primary logical order." << std::endl;
        ++test_failures;
    }
    if (pressure_patches_zero_relaxed.size() != primary_pressure_seed_dofs.size() ||
        pressure_patches_zero_strict.size() != primary_pressure_seed_dofs.size())
    {
        std::cerr << "FAILED: zero-E_h pressure patch count differs from the pressure-cell count." << std::endl;
        ++test_failures;
    }
    for (std::size_t k = 0; k < std::min({ pressure_patches_zero_relaxed.size(),
                                           pressure_patches_zero_strict.size(),
                                           primary_pressure_seed_dofs.size() });
         ++k)
    {
        const int pressure_seed_dof = primary_pressure_seed_dofs[k];
        const std::set<int> expected_standard_patch = set_from_vector(cell_dof_to_closure_dofs.at(pressure_seed_dof));
        if (!check_ordered_set("zero-E_h RELAXED patch for pressure seed " + std::to_string(pressure_seed_dof),
                               expected_standard_patch,
                               pressure_patches_zero_relaxed[k]))
        {
            ++test_failures;
        }
        if (!check_ordered_set("zero-E_h STRICT patch for pressure seed " + std::to_string(pressure_seed_dof),
                               expected_standard_patch,
                               pressure_patches_zero_strict[k]))
        {
            ++test_failures;
        }
    }

    std::vector<std::vector<int>> pressure_patches_alternate;
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelPressureCellSeededCAVPatches(
        pressure_patches_alternate,
        constructed_pressure_seed_dofs,
        num_dofs_per_proc,
        p_dof_index_idx,
        level,
        zero_eulerian_elasticity_mat,
        pressure_map_data,
        1,
        alternate_traversal_order,
        IBAMR::CouplingAwareASMClosurePolicy::RELAXED,
        relative_numerical_zero_tol);
    if (constructed_pressure_seed_dofs != alternate_pressure_seed_dofs ||
        constructed_pressure_seed_dofs == primary_pressure_seed_dofs)
    {
        std::cerr << "FAILED: alternate pressure traversal order was not preserved." << std::endl;
        ++test_failures;
    }

    std::vector<int> expected_stride_two_seed_dofs;
    for (std::size_t k = 0; k < primary_pressure_seed_dofs.size(); k += 2)
    {
        expected_stride_two_seed_dofs.push_back(primary_pressure_seed_dofs[k]);
    }
    std::vector<std::vector<int>> pressure_patches_stride_two;
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelPressureCellSeededCAVPatches(
        pressure_patches_stride_two,
        constructed_pressure_seed_dofs,
        num_dofs_per_proc,
        p_dof_index_idx,
        level,
        zero_eulerian_elasticity_mat,
        pressure_map_data,
        2,
        primary_traversal_order,
        IBAMR::CouplingAwareASMClosurePolicy::RELAXED,
        relative_numerical_zero_tol);
    if (constructed_pressure_seed_dofs != expected_stride_two_seed_dofs ||
        pressure_patches_stride_two.size() != expected_stride_two_seed_dofs.size())
    {
        std::cerr << "FAILED: pressure seed stride was not applied after logical ordering and de-duplication."
                  << std::endl;
        ++test_failures;
    }
    for (std::size_t k = 0; k < std::min(pressure_patches_stride_two.size(), expected_stride_two_seed_dofs.size()); ++k)
    {
        const int pressure_seed_dof = expected_stride_two_seed_dofs[k];
        if (!check_ordered_set("stride-two standard patch for pressure seed " + std::to_string(pressure_seed_dof),
                               set_from_vector(cell_dof_to_closure_dofs.at(pressure_seed_dof)),
                               pressure_patches_stride_two[k]))
        {
            ++test_failures;
        }
    }

    if (!primary_pressure_seed_dofs.empty())
    {
        const int pressure_seed_dof = primary_pressure_seed_dofs.front();
        const std::set<int> standard_patch = set_from_vector(cell_dof_to_closure_dofs.at(pressure_seed_dof));
        std::set<int> standard_velocity_dofs;
        for (const int dof : standard_patch)
        {
            if (velocity_dof_to_component_axis.find(dof) != velocity_dof_to_component_axis.end())
            {
                standard_velocity_dofs.insert(dof);
            }
        }
        std::vector<int> external_velocity_dofs;
        for (const auto& entry : velocity_dof_to_component_axis)
        {
            if (standard_velocity_dofs.find(entry.first) == standard_velocity_dofs.end())
            {
                external_velocity_dofs.push_back(entry.first);
            }
        }
        std::sort(external_velocity_dofs.begin(), external_velocity_dofs.end());
        if (standard_velocity_dofs.empty() || external_velocity_dofs.size() < 2)
        {
            std::cerr << "FAILED: insufficient velocity DOFs for pressure-cell graph tests." << std::endl;
            ++test_failures;
        }
        else
        {
            const int standard_velocity_dof = *standard_velocity_dofs.begin();
            const int external_velocity_dof = external_velocity_dofs[0];
            std::set<int> expanded_velocity_dofs = standard_velocity_dofs;
            expanded_velocity_dofs.insert(external_velocity_dof);
            const std::set<int> expected_relaxed =
                reference_pressure_seeded_patch(pressure_seed_dof,
                                                expanded_velocity_dofs,
                                                velocity_dof_to_adjacent_cell_dofs,
                                                cell_dof_to_closure_dofs,
                                                velocity_dof_to_component_axis,
                                                IBAMR::CouplingAwareASMClosurePolicy::RELAXED);
            const std::set<int> expected_strict =
                reference_pressure_seeded_patch(pressure_seed_dof,
                                                expanded_velocity_dofs,
                                                velocity_dof_to_adjacent_cell_dofs,
                                                cell_dof_to_closure_dofs,
                                                velocity_dof_to_component_axis,
                                                IBAMR::CouplingAwareASMClosurePolicy::STRICT);

            std::vector<int> row_relaxed_patch;
            for (const bool use_column_edge : { false, true })
            {
                const int row = use_column_edge ? external_velocity_dof : standard_velocity_dof;
                const int column = use_column_edge ? standard_velocity_dof : external_velocity_dof;
                Mat one_edge_eulerian_elasticity_mat =
                    make_eulerian_elasticity_pattern(n_global_dofs, { { row, column, 1.0 } });
                std::vector<std::vector<int>> pressure_patches_relaxed;
                IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelPressureCellSeededCAVPatches(
                    pressure_patches_relaxed,
                    constructed_pressure_seed_dofs,
                    num_dofs_per_proc,
                    p_dof_index_idx,
                    level,
                    one_edge_eulerian_elasticity_mat,
                    pressure_map_data,
                    1,
                    primary_traversal_order,
                    IBAMR::CouplingAwareASMClosurePolicy::RELAXED,
                    relative_numerical_zero_tol);
                std::vector<std::vector<int>> pressure_patches_strict;
                IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelPressureCellSeededCAVPatches(
                    pressure_patches_strict,
                    constructed_pressure_seed_dofs,
                    num_dofs_per_proc,
                    p_dof_index_idx,
                    level,
                    one_edge_eulerian_elasticity_mat,
                    pressure_map_data,
                    1,
                    primary_traversal_order,
                    IBAMR::CouplingAwareASMClosurePolicy::STRICT,
                    relative_numerical_zero_tol);
                if (pressure_patches_relaxed.empty() ||
                    !check_ordered_set(std::string(use_column_edge ? "column-E_h" : "row-E_h") + " RELAXED patch",
                                       expected_relaxed,
                                       pressure_patches_relaxed.front()))
                {
                    ++test_failures;
                }
                if (pressure_patches_strict.empty() ||
                    !check_ordered_set(std::string(use_column_edge ? "column-E_h" : "row-E_h") + " STRICT patch",
                                       expected_strict,
                                       pressure_patches_strict.front()))
                {
                    ++test_failures;
                }
                if (!use_column_edge && !pressure_patches_relaxed.empty())
                {
                    row_relaxed_patch = pressure_patches_relaxed.front();
                }
                else if (use_column_edge && !pressure_patches_relaxed.empty() &&
                         pressure_patches_relaxed.front() != row_relaxed_patch)
                {
                    std::cerr << "FAILED: row and column E_h edges produce different undirected expansions."
                              << std::endl;
                    ++test_failures;
                }
                ierr = MatDestroy(&one_edge_eulerian_elasticity_mat);
                IBTK_CHKERRQ(ierr);
            }
            if (expected_relaxed == expected_strict ||
                expected_relaxed.find(external_velocity_dof) == expected_relaxed.end() ||
                expected_strict.find(external_velocity_dof) != expected_strict.end())
            {
                std::cerr << "FAILED: RELAXED and STRICT discriminator did not activate." << std::endl;
                ++test_failures;
            }

            const int ignored_velocity_dof = external_velocity_dofs[1];
            Mat threshold_eulerian_elasticity_mat = make_eulerian_elasticity_pattern(
                n_global_dofs,
                { { standard_velocity_dof, external_velocity_dof, 1.0 },
                  { standard_velocity_dof, ignored_velocity_dof, 0.1 * relative_numerical_zero_tol } });
            std::vector<std::vector<int>> pressure_patches_threshold;
            IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelPressureCellSeededCAVPatches(
                pressure_patches_threshold,
                constructed_pressure_seed_dofs,
                num_dofs_per_proc,
                p_dof_index_idx,
                level,
                threshold_eulerian_elasticity_mat,
                pressure_map_data,
                1,
                primary_traversal_order,
                IBAMR::CouplingAwareASMClosurePolicy::RELAXED,
                relative_numerical_zero_tol);
            if (pressure_patches_threshold.empty() || !check_ordered_set("relative-zero-filtered RELAXED patch",
                                                                         expected_relaxed,
                                                                         pressure_patches_threshold.front()))
            {
                ++test_failures;
            }
            ierr = MatDestroy(&threshold_eulerian_elasticity_mat);
            IBTK_CHKERRQ(ierr);

            int supported_pressure_dof = -1;
            std::set<int> supported_cell_velocity_dofs;
            for (const int candidate_pressure_dof : primary_pressure_seed_dofs)
            {
                if (candidate_pressure_dof == pressure_seed_dof) continue;
                std::set<int> candidate_velocity_dofs;
                for (const int dof : cell_dof_to_closure_dofs.at(candidate_pressure_dof))
                {
                    if (velocity_dof_to_component_axis.find(dof) != velocity_dof_to_component_axis.end())
                    {
                        candidate_velocity_dofs.insert(dof);
                    }
                }
                if (!std::includes(standard_velocity_dofs.begin(),
                                   standard_velocity_dofs.end(),
                                   candidate_velocity_dofs.begin(),
                                   candidate_velocity_dofs.end()))
                {
                    supported_pressure_dof = candidate_pressure_dof;
                    supported_cell_velocity_dofs = candidate_velocity_dofs;
                    break;
                }
            }
            std::vector<std::tuple<int, int, double>> full_support_entries;
            std::set<int> fully_expanded_velocity_dofs = standard_velocity_dofs;
            for (const int dof : supported_cell_velocity_dofs)
            {
                if (standard_velocity_dofs.find(dof) != standard_velocity_dofs.end()) continue;
                full_support_entries.emplace_back(standard_velocity_dof, dof, 1.0);
                fully_expanded_velocity_dofs.insert(dof);
            }
            if (supported_pressure_dof < 0 || full_support_entries.empty())
            {
                std::cerr << "FAILED: no candidate cell is available for the complete STRICT support test."
                          << std::endl;
                ++test_failures;
            }
            else
            {
                Mat full_support_eulerian_elasticity_mat =
                    make_eulerian_elasticity_pattern(n_global_dofs, full_support_entries);
                std::vector<std::vector<int>> pressure_patches_full_strict;
                IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelPressureCellSeededCAVPatches(
                    pressure_patches_full_strict,
                    constructed_pressure_seed_dofs,
                    num_dofs_per_proc,
                    p_dof_index_idx,
                    level,
                    full_support_eulerian_elasticity_mat,
                    pressure_map_data,
                    1,
                    primary_traversal_order,
                    IBAMR::CouplingAwareASMClosurePolicy::STRICT,
                    relative_numerical_zero_tol);
                const std::set<int> expected_full_strict =
                    reference_pressure_seeded_patch(pressure_seed_dof,
                                                    fully_expanded_velocity_dofs,
                                                    velocity_dof_to_adjacent_cell_dofs,
                                                    cell_dof_to_closure_dofs,
                                                    velocity_dof_to_component_axis,
                                                    IBAMR::CouplingAwareASMClosurePolicy::STRICT);
                if (pressure_patches_full_strict.empty() ||
                    !check_ordered_set(
                        "complete-support STRICT patch", expected_full_strict, pressure_patches_full_strict.front()) ||
                    expected_full_strict.find(supported_pressure_dof) == expected_full_strict.end())
                {
                    ++test_failures;
                }
                ierr = MatDestroy(&full_support_eulerian_elasticity_mat);
                IBTK_CHKERRQ(ierr);
            }
        }
    }

    ierr = MatDestroy(&zero_eulerian_elasticity_mat);
    IBTK_CHKERRQ(ierr);

    ierr = MatDestroy(&level_mat);
    IBTK_CHKERRQ(ierr);
    level->deallocatePatchData(u_dof_index_idx);
    level->deallocatePatchData(p_dof_index_idx);

    const bool pressure_seed_pair_map_absent =
        !MapAccess::velocitySeedPairMapIsBuilt(pressure_map_data) &&
        MapAccess::getVelocityDofToPairedSeedVelocityDofs(pressure_map_data).empty();
    if (!pressure_seed_pair_map_absent) ++test_failures;
    pout << "pressure_seed_pair_map_absent = " << (pressure_seed_pair_map_absent ? "true" : "false") << std::endl;
    pout << "velocity_relaxed_seed_pair_map_absent = " << (velocity_relaxed_seed_pair_map_absent ? "true" : "false")
         << std::endl;
    pout << "velocity_strict_seed_pair_map_count_valid = "
         << (velocity_strict_seed_pair_map_count_valid ? "true" : "false") << std::endl;
    pout << "velocity_strict_seed_pair_map_keys = " << velocity_dof_to_paired_seed_velocity_dofs.size() << std::endl;
    pout << "velocity_strict_seed_pair_map_values = " << velocity_seed_pair_value_count << std::endl;
    plog << "Input database:\n";
    input_db->printClassData(plog);
    pout << "test_failures = " << test_failures << std::endl;
    return test_failures > 0 ? 1 : 0;
}
