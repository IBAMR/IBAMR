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

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>

#include <petscmat.h>

#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

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
matlab_extract_coupled_dofs_relaxed(const int seed_velocity_dof,
                                    Mat A00_mat,
                                    const std::unordered_map<int, std::vector<int>>& velocity_dof_to_adjacent_cell_dofs,
                                    const std::unordered_map<int, std::vector<int>>& cell_dof_to_closure_dofs,
                                    const std::map<int, int>& velocity_dof_to_component_axis)
{
    std::set<int> velocity_dofs = { seed_velocity_dof };

    PetscInt ncols = 0;
    const PetscInt* cols = nullptr;
    const PetscInt row = static_cast<PetscInt>(seed_velocity_dof);
    int ierr = MatGetRow(A00_mat, row, &ncols, &cols, nullptr);
    IBTK_CHKERRQ(ierr);
    for (PetscInt k = 0; k < ncols; ++k) velocity_dofs.insert(static_cast<int>(cols[k]));
    ierr = MatRestoreRow(A00_mat, row, &ncols, &cols, nullptr);
    IBTK_CHKERRQ(ierr);

    std::set<int> pressure_dofs;
    for (const int velocity_dof : velocity_dofs)
    {
        const auto vel_it = velocity_dof_to_adjacent_cell_dofs.find(velocity_dof);
        if (vel_it != velocity_dof_to_adjacent_cell_dofs.end())
        {
            pressure_dofs.insert(vel_it->second.begin(), vel_it->second.end());
        }
    }

    std::set<int> additional_velocity_dofs;
    for (const int pressure_dof : pressure_dofs)
    {
        const auto closure_it = cell_dof_to_closure_dofs.find(pressure_dof);
        if (closure_it == cell_dof_to_closure_dofs.end()) continue;
        for (const int dof : closure_it->second)
        {
            if (velocity_dof_to_component_axis.find(dof) != velocity_dof_to_component_axis.end())
            {
                additional_velocity_dofs.insert(dof);
            }
        }
    }

    std::set<int> coupled_dofs = velocity_dofs;
    coupled_dofs.insert(additional_velocity_dofs.begin(), additional_velocity_dofs.end());
    coupled_dofs.insert(pressure_dofs.begin(), pressure_dofs.end());
    return coupled_dofs;
}

std::set<int>
matlab_extract_coupled_dofs_strict(const int seed_x_velocity_dof,
                                   const int paired_y_velocity_dof,
                                   Mat A00_mat,
                                   const std::unordered_map<int, std::vector<int>>& velocity_dof_to_adjacent_cell_dofs,
                                   const std::unordered_map<int, std::vector<int>>& cell_dof_to_closure_dofs,
                                   const std::map<int, int>& velocity_dof_to_component_axis)
{
    std::set<int> initial_velocity_dofs = { seed_x_velocity_dof, paired_y_velocity_dof };

    for (const int seed_velocity_dof : std::set<int>{ seed_x_velocity_dof, paired_y_velocity_dof })
    {
        PetscInt ncols = 0;
        const PetscInt* cols = nullptr;
        const PetscInt row = static_cast<PetscInt>(seed_velocity_dof);
        int ierr = MatGetRow(A00_mat, row, &ncols, &cols, nullptr);
        IBTK_CHKERRQ(ierr);
        for (PetscInt k = 0; k < ncols; ++k) initial_velocity_dofs.insert(static_cast<int>(cols[k]));
        ierr = MatRestoreRow(A00_mat, row, &ncols, &cols, nullptr);
        IBTK_CHKERRQ(ierr);
    }

    std::set<int> candidate_pressure_dofs;
    for (const int velocity_dof : initial_velocity_dofs)
    {
        const auto vel_it = velocity_dof_to_adjacent_cell_dofs.find(velocity_dof);
        if (vel_it != velocity_dof_to_adjacent_cell_dofs.end())
        {
            candidate_pressure_dofs.insert(vel_it->second.begin(), vel_it->second.end());
        }
    }

    std::set<int> valid_pressure_dofs;
    for (const int pressure_dof : candidate_pressure_dofs)
    {
        const auto closure_it = cell_dof_to_closure_dofs.find(pressure_dof);
        if (closure_it == cell_dof_to_closure_dofs.end()) continue;

        std::set<int> pressure_coupled_velocities;
        for (const int dof : closure_it->second)
        {
            if (velocity_dof_to_component_axis.find(dof) != velocity_dof_to_component_axis.end())
            {
                pressure_coupled_velocities.insert(dof);
            }
        }

        bool is_valid = true;
        for (const int dof : pressure_coupled_velocities)
        {
            if (initial_velocity_dofs.find(dof) == initial_velocity_dofs.end())
            {
                is_valid = false;
                break;
            }
        }
        if (is_valid) valid_pressure_dofs.insert(pressure_dof);
    }

    std::set<int> final_selected_velocities;
    for (const int pressure_dof : valid_pressure_dofs)
    {
        const auto closure_it = cell_dof_to_closure_dofs.find(pressure_dof);
        if (closure_it == cell_dof_to_closure_dofs.end()) continue;
        for (const int dof : closure_it->second)
        {
            if (velocity_dof_to_component_axis.find(dof) != velocity_dof_to_component_axis.end())
            {
                final_selected_velocities.insert(dof);
            }
        }
    }

    std::set<int> coupled_dofs = final_selected_velocities;
    coupled_dofs.insert(valid_pressure_dofs.begin(), valid_pressure_dofs.end());
    return coupled_dofs;
}

void
construct_coupling_aware_overlap_subdomains_with_cell_closure(
    std::vector<std::set<int>>& overlap_is,
    const std::vector<std::set<int>>& nonoverlap_is,
    Mat A00_mat,
    const std::unordered_map<int, std::vector<int>>& velocity_dof_to_adjacent_cell_dofs,
    const std::unordered_map<int, std::vector<int>>& cell_dof_to_closure_dofs,
    const std::map<int, int>& velocity_dof_to_component_axis,
    const int seed_velocity_axis,
    const IBAMR::CouplingAwareASMClosurePolicy closure_policy,
    const std::unordered_map<int, std::vector<int>>& velocity_dof_to_paired_seed_velocity_dofs = {})
{
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
        if (closure_policy == IBAMR::CouplingAwareASMClosurePolicy::STRICT)
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
        initial_velocity_dofs = initial_seed_components;
        PetscInt first_local_row = -1, row_end = -1;
        int ierr = MatGetOwnershipRange(A00_mat, &first_local_row, &row_end);
        IBTK_CHKERRQ(ierr);
        for (const int velocity_dof : initial_seed_components)
        {
            const PetscInt row = static_cast<PetscInt>(velocity_dof);
            if (row < first_local_row || row >= row_end) continue;
            PetscInt ncols = 0;
            const PetscInt* cols = nullptr;
            ierr = MatGetRow(A00_mat, row, &ncols, &cols, nullptr);
            IBTK_CHKERRQ(ierr);
            for (PetscInt col_idx = 0; col_idx < ncols; ++col_idx)
            {
                initial_velocity_dofs.insert(static_cast<int>(cols[col_idx]));
            }
            ierr = MatRestoreRow(A00_mat, row, &ncols, &cols, nullptr);
            IBTK_CHKERRQ(ierr);
        }

        IBAMR::StaggeredStokesPETScMatUtilities::findCoupledCellDofsFromA00(involved_cell_dofs,
                                                                            A00_mat,
                                                                            initial_seed_components,
                                                                            velocity_dof_to_adjacent_cell_dofs,
                                                                            cell_dof_to_closure_dofs);
        if (closure_policy == IBAMR::CouplingAwareASMClosurePolicy::STRICT)
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
        for (const int cell_dof : involved_cell_dofs)
        {
            const auto closure_it = cell_dof_to_closure_dofs.find(cell_dof);
            if (closure_it == cell_dof_to_closure_dofs.end()) continue;
            closure_dofs.insert(closure_it->second.begin(), closure_it->second.end());
        }
        overlap_is[k] = closure_dofs;
        if (closure_policy == IBAMR::CouplingAwareASMClosurePolicy::RELAXED)
        {
            overlap_is[k].insert(initial_velocity_dofs.begin(), initial_velocity_dofs.end());
        }
    }
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app_initializer->getInputDatabase();
    int test_failures = 0;

    // Synthetic case intentionally chosen so strict and relaxed MATLAB
    // extractors differ.
    const std::unordered_map<int, std::vector<int>> velocity_dof_to_adjacent_cell_dofs = {
        { 0, { 10 } }, { 1, { 10, 11 } }, { 2, { 11 } }, { 3, { 10 } }, { 4, { 11 } }, { 5, { 11 } }
    };
    const std::unordered_map<int, std::vector<int>> cell_dof_to_closure_dofs = {
        { 10, { 0, 1, 3, 10 } }, { 11, { 1, 2, 4, 11 } }
    };
    const std::map<int, int> velocity_dof_to_component_axis = { { 0, 0 }, { 1, 0 }, { 2, 0 },
                                                                { 3, 1 }, { 4, 1 }, { 5, 1 } };

    Mat A00_mat = nullptr;
    int ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 6, 6, 4, nullptr, 4, nullptr, &A00_mat);
    IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(A00_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    IBTK_CHKERRQ(ierr);

    auto add_row = [&A00_mat](const int row_idx, const std::vector<int>& col_idxs)
    {
        std::vector<PetscInt> cols(col_idxs.begin(), col_idxs.end());
        std::vector<PetscScalar> vals(cols.size(), 1.0);
        const PetscInt row = static_cast<PetscInt>(row_idx);
        const int ierr_set_values =
            MatSetValues(A00_mat, 1, &row, static_cast<PetscInt>(cols.size()), cols.data(), vals.data(), INSERT_VALUES);
        IBTK_CHKERRQ(ierr_set_values);
    };
    add_row(0, { 0, 1, 3 });
    add_row(1, { 0, 1, 2 });
    add_row(2, { 1, 2, 4 });
    add_row(3, { 0, 3 });
    add_row(4, { 2, 4, 5 });
    add_row(5, { 4, 5 });

    ierr = MatAssemblyBegin(A00_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A00_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);

    const int seed_x_velocity_dof = 0;
    const int paired_y_velocity_dof = 3;

    const std::set<int> matlab_relaxed = matlab_extract_coupled_dofs_relaxed(seed_x_velocity_dof,
                                                                             A00_mat,
                                                                             velocity_dof_to_adjacent_cell_dofs,
                                                                             cell_dof_to_closure_dofs,
                                                                             velocity_dof_to_component_axis);
    const std::set<int> matlab_strict = matlab_extract_coupled_dofs_strict(seed_x_velocity_dof,
                                                                           paired_y_velocity_dof,
                                                                           A00_mat,
                                                                           velocity_dof_to_adjacent_cell_dofs,
                                                                           cell_dof_to_closure_dofs,
                                                                           velocity_dof_to_component_axis);

    if (!check_set("MATLAB relaxed expected", std::set<int>{ 0, 1, 2, 3, 4, 10, 11 }, matlab_relaxed)) ++test_failures;
    if (!check_set("MATLAB strict expected", std::set<int>{ 0, 1, 3, 10 }, matlab_strict)) ++test_failures;
    if (!(matlab_strict.size() < matlab_relaxed.size()))
    {
        std::cerr << "FAILED: expected strict MATLAB extractor set to be smaller than relaxed set." << std::endl;
        ++test_failures;
    }

    std::vector<std::set<int>> overlap_relaxed(1), overlap_strict(1), nonoverlap_is(1);
    nonoverlap_is[0] = { seed_x_velocity_dof };
    overlap_relaxed[0] = nonoverlap_is[0];
    overlap_strict[0] = nonoverlap_is[0];

    construct_coupling_aware_overlap_subdomains_with_cell_closure(overlap_relaxed,
                                                                  nonoverlap_is,
                                                                  A00_mat,
                                                                  velocity_dof_to_adjacent_cell_dofs,
                                                                  cell_dof_to_closure_dofs,
                                                                  velocity_dof_to_component_axis,
                                                                  /*seed_velocity_axis*/ 0,
                                                                  IBAMR::CouplingAwareASMClosurePolicy::RELAXED);

    construct_coupling_aware_overlap_subdomains_with_cell_closure(overlap_strict,
                                                                  nonoverlap_is,
                                                                  A00_mat,
                                                                  velocity_dof_to_adjacent_cell_dofs,
                                                                  cell_dof_to_closure_dofs,
                                                                  velocity_dof_to_component_axis,
                                                                  /*seed_velocity_axis*/ 0,
                                                                  IBAMR::CouplingAwareASMClosurePolicy::STRICT);

    if (!check_set("IBAMR relaxed overlap matches MATLAB relaxed extraction", matlab_relaxed, overlap_relaxed[0]))
        ++test_failures;
    if (!check_set("IBAMR strict overlap matches MATLAB strict extraction", matlab_strict, overlap_strict[0]))
        ++test_failures;

    ierr = MatDestroy(&A00_mat);
    IBTK_CHKERRQ(ierr);
    plog << "Input database:\n";
    input_db->printClassData(plog);
    pout << "test_failures = " << test_failures << std::endl;
    return test_failures > 0 ? 1 : 0;
}
