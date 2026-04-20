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

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/ibtk_utilities.h>

#include <petscmat.h>

#include <CellVariable.h>
#include <CoarseFineBoundary.h>
#include <PoissonSpecifications.h>
#include <SideVariable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <iostream>
#include <set>
#include <sstream>
#include <string>
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
build_initial_velocity_dofs(const std::set<int>& seed_velocity_dofs, Mat A00_mat)
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
        ierr = MatGetRow(A00_mat, row, &ncols, &cols, nullptr);
        IBTK_CHKERRQ(ierr);
        for (PetscInt k = 0; k < ncols; ++k) initial_velocity_dofs.insert(static_cast<int>(cols[k]));
        ierr = MatRestoreRow(A00_mat, row, &ncols, &cols, nullptr);
        IBTK_CHKERRQ(ierr);
    }
    return initial_velocity_dofs;
}

std::set<int>
reference_extract_coupled_dofs_relaxed(
    const int seed_velocity_dof,
    Mat A00_mat,
    const std::unordered_map<int, std::vector<int>>& velocity_dof_to_adjacent_cell_dofs,
    const std::unordered_map<int, std::vector<int>>& cell_dof_to_closure_dofs)
{
    const std::set<int> initial_velocity_dofs =
        build_initial_velocity_dofs(std::set<int>{ seed_velocity_dof }, A00_mat);

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
    const std::unordered_map<int, std::vector<int>>& velocity_dof_to_paired_seed_velocity_dofs)
{
    std::set<int> seed_velocity_dofs = { seed_velocity_dof };
    const auto pair_it = velocity_dof_to_paired_seed_velocity_dofs.find(seed_velocity_dof);
    if (pair_it != velocity_dof_to_paired_seed_velocity_dofs.end())
    {
        seed_velocity_dofs.insert(pair_it->second.begin(), pair_it->second.end());
    }

    const std::set<int> initial_velocity_dofs = build_initial_velocity_dofs(seed_velocity_dofs, A00_mat);

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

    std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);
    PoissonSpecifications u_problem_coefs("reference_parity_poisson");
    u_problem_coefs.setCConstant(1.0);
    u_problem_coefs.setDConstant(-1.0);
    Mat level_mat = nullptr;
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(
        level_mat, u_problem_coefs, u_bc_coefs, 0.0, num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, level);

    std::vector<int> seed_velocity_dofs;
    IBAMR::StaggeredStokesPETScMatUtilities::computePatchLevelCouplingAwareASMSeedVelocityDofs(
        seed_velocity_dofs, u_dof_index_idx, level, map_data, 0, 1, IBAMR::CouplingAwareASMSeedTraversalOrder::I_J);

    std::vector<std::set<int>> overlap_relaxed, nonoverlap_relaxed;
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
        IBAMR::CouplingAwareASMSeedTraversalOrder::I_J,
        IBAMR::CouplingAwareASMClosurePolicy::RELAXED);

    std::vector<std::set<int>> overlap_strict, nonoverlap_strict;
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
        IBAMR::CouplingAwareASMSeedTraversalOrder::I_J,
        IBAMR::CouplingAwareASMClosurePolicy::STRICT);

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
                                                   map_data.velocity_dof_to_adjacent_cell_dofs,
                                                   map_data.cell_dof_to_closure_dofs);
        const std::set<int> expected_strict =
            reference_extract_coupled_dofs_strict(seed_velocity_dof,
                                                  level_mat,
                                                  map_data.velocity_dof_to_adjacent_cell_dofs,
                                                  map_data.cell_dof_to_closure_dofs,
                                                  map_data.velocity_dof_to_component_axis,
                                                  map_data.velocity_dof_to_paired_seed_velocity_dofs);

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

    int ierr = MatDestroy(&level_mat);
    IBTK_CHKERRQ(ierr);
    level->deallocatePatchData(u_dof_index_idx);
    level->deallocatePatchData(p_dof_index_idx);

    plog << "Input database:\n";
    input_db->printClassData(plog);
    pout << "test_failures = " << test_failures << std::endl;
    return test_failures > 0 ? 1 : 0;
}
