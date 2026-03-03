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

void
construct_coupling_aware_overlap_subdomains_with_cell_closure(
    std::vector<std::set<int>>& overlap_is,
    const std::vector<std::set<int>>& nonoverlap_is,
    Mat A00_mat,
    const std::map<int, std::set<int>>& velocity_dof_to_adjacent_cell_dofs,
    const std::map<int, std::set<int>>& cell_dof_to_closure_dofs,
    const std::map<int, int>& velocity_dof_to_component_axis,
    const int seed_velocity_axis,
    const IBAMR::CouplingAwareASMClosurePolicy closure_policy = IBAMR::CouplingAwareASMClosurePolicy::RELAXED,
    const std::map<int, std::set<int>>& velocity_dof_to_paired_seed_velocity_dofs = {})
{
    if (overlap_is.size() != nonoverlap_is.size())
    {
        std::cerr << "FAILED: overlap/nonoverlap size mismatch." << std::endl;
        return;
    }

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

    const std::map<int, std::set<int>> velocity_dof_to_adjacent_cell_dofs = {
        { 0, { 10 } }, { 1, { 10, 11 } }, { 2, { 11 } }, { 3, { 11, 12 } }, { 4, { 12 } }
    };
    const std::map<int, std::set<int>> cell_dof_to_closure_dofs = { { 10, { 0, 1, 10 } },
                                                                    { 11, { 1, 2, 3, 11 } },
                                                                    { 12, { 3, 4, 12 } } };
    const std::map<int, int> velocity_dof_to_component_axis = { { 0, 0 }, { 1, 1 }, { 2, 0 }, { 3, 1 }, { 4, 0 } };

    Mat SAJ_mat = nullptr;
    int ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 13, 13, 3, nullptr, 3, nullptr, &SAJ_mat);
    IBTK_CHKERRQ(ierr);

    auto add_row = [&SAJ_mat](const int row, const std::vector<int>& col_idxs)
    {
        std::vector<PetscInt> cols(col_idxs.begin(), col_idxs.end());
        std::vector<PetscScalar> vals(cols.size(), 1.0);
        const PetscInt row_i = static_cast<PetscInt>(row);
        const int ierr_set_values = MatSetValues(
            SAJ_mat, 1, &row_i, static_cast<PetscInt>(cols.size()), cols.data(), vals.data(), INSERT_VALUES);
        IBTK_CHKERRQ(ierr_set_values);
    };
    add_row(0, { 0, 1, 2 });
    add_row(1, { 0, 1 });
    add_row(2, { 0, 2, 3 });
    add_row(3, { 2, 3, 4 });
    add_row(4, { 3, 4 });

    ierr = MatAssemblyBegin(SAJ_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(SAJ_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);

    std::set<int> involved_cell_dofs;
    IBAMR::StaggeredStokesPETScMatUtilities::findCoupledCellDofsFromA00(
        involved_cell_dofs, SAJ_mat, std::set<int>{ 0 }, velocity_dof_to_adjacent_cell_dofs, cell_dof_to_closure_dofs);
    if (!check_set("coupled cell DOFs", std::set<int>{ 10, 11 }, involved_cell_dofs)) ++test_failures;

    std::set<int> closure_dofs;
    for (const int cell_dof : involved_cell_dofs)
    {
        const auto closure_it = cell_dof_to_closure_dofs.find(cell_dof);
        if (closure_it == cell_dof_to_closure_dofs.end()) continue;
        closure_dofs.insert(closure_it->second.begin(), closure_it->second.end());
    }
    if (!check_set("cell closure DOFs", std::set<int>{ 0, 1, 2, 3, 10, 11 }, closure_dofs)) ++test_failures;

    std::vector<std::set<int>> overlap_is(2), nonoverlap_is(2);
    nonoverlap_is[0] = { 0, 1, 10 };
    overlap_is[0] = { 0, 1, 10 };
    nonoverlap_is[1] = { 1, 12 };
    overlap_is[1] = { 1, 12 };
    const std::vector<std::set<int>> nonoverlap_is_expected = nonoverlap_is;
    construct_coupling_aware_overlap_subdomains_with_cell_closure(overlap_is,
                                                                  nonoverlap_is,
                                                                  SAJ_mat,
                                                                  velocity_dof_to_adjacent_cell_dofs,
                                                                  cell_dof_to_closure_dofs,
                                                                  velocity_dof_to_component_axis,
                                                                  /*seed_velocity_axis*/ 0);

    if (!check_set("constructed overlap IS[0]", std::set<int>{ 0, 1, 2, 3, 10, 11 }, overlap_is[0])) ++test_failures;
    if (!check_set("constructed overlap IS[1]", std::set<int>{ 1, 12 }, overlap_is[1])) ++test_failures;
    if (!check_set("nonoverlap IS[0] unchanged (axis 0)", nonoverlap_is_expected[0], nonoverlap_is[0])) ++test_failures;
    if (!check_set("nonoverlap IS[1] unchanged (axis 0)", nonoverlap_is_expected[1], nonoverlap_is[1])) ++test_failures;

    std::vector<std::set<int>> overlap_is_axis1 = nonoverlap_is;
    construct_coupling_aware_overlap_subdomains_with_cell_closure(overlap_is_axis1,
                                                                  nonoverlap_is,
                                                                  SAJ_mat,
                                                                  velocity_dof_to_adjacent_cell_dofs,
                                                                  cell_dof_to_closure_dofs,
                                                                  velocity_dof_to_component_axis,
                                                                  /*seed_velocity_axis*/ 1);

    if (!check_set("constructed overlap IS[0] axis 1", std::set<int>{ 0, 1, 2, 3, 10, 11 }, overlap_is_axis1[0]))
        ++test_failures;
    if (!check_set("constructed overlap IS[1] axis 1", std::set<int>{ 0, 1, 2, 3, 10, 11 }, overlap_is_axis1[1]))
        ++test_failures;
    if (!check_set("nonoverlap IS[0] unchanged (axis 1)", nonoverlap_is_expected[0], nonoverlap_is[0])) ++test_failures;
    if (!check_set("nonoverlap IS[1] unchanged (axis 1)", nonoverlap_is_expected[1], nonoverlap_is[1])) ++test_failures;

    ierr = MatDestroy(&SAJ_mat);
    IBTK_CHKERRQ(ierr);
    plog << "Input database:\n";
    input_db->printClassData(plog);
    pout << "test_failures = " << test_failures << std::endl;
    return test_failures > 0 ? 1 : 0;
}
