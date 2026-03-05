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

#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>

#include <petscmat.h>

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

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
matlab_extract_coupled_dofs(const int seed_velocity_dof,
                            Mat A00_mat,
                            const std::map<int, std::set<int>>& velocity_dof_to_adjacent_cell_dofs,
                            const std::map<int, std::set<int>>& cell_dof_to_closure_dofs,
                            const std::map<int, int>& velocity_dof_to_component_axis)
{
    std::set<int> initial_velocity_dofs = { seed_velocity_dof };

    PetscInt ncols = 0;
    const PetscInt* cols = nullptr;
    const PetscInt row = static_cast<PetscInt>(seed_velocity_dof);
    int ierr = MatGetRow(A00_mat, row, &ncols, &cols, nullptr);
    IBTK_CHKERRQ(ierr);
    for (PetscInt k = 0; k < ncols; ++k) initial_velocity_dofs.insert(static_cast<int>(cols[k]));
    ierr = MatRestoreRow(A00_mat, row, &ncols, &cols, nullptr);
    IBTK_CHKERRQ(ierr);

    std::set<int> pressure_dofs;
    for (const int velocity_dof : initial_velocity_dofs)
    {
        const auto vel_it = velocity_dof_to_adjacent_cell_dofs.find(velocity_dof);
        if (vel_it != velocity_dof_to_adjacent_cell_dofs.end())
        {
            pressure_dofs.insert(vel_it->second.begin(), vel_it->second.end());
        }
    }

    std::set<int> final_velocity_dofs;
    for (const int pressure_dof : pressure_dofs)
    {
        const auto closure_it = cell_dof_to_closure_dofs.find(pressure_dof);
        if (closure_it == cell_dof_to_closure_dofs.end()) continue;
        for (const int dof : closure_it->second)
        {
            if (velocity_dof_to_component_axis.find(dof) != velocity_dof_to_component_axis.end())
            {
                final_velocity_dofs.insert(dof);
            }
        }
    }

    std::set<int> coupled_dofs = final_velocity_dofs;
    coupled_dofs.insert(pressure_dofs.begin(), pressure_dofs.end());
    return coupled_dofs;
}

void
construct_coupling_aware_overlap_subdomains_with_cell_closure(
    std::vector<std::set<int>>& overlap_is,
    const std::vector<std::set<int>>& nonoverlap_is,
    Mat A00_mat,
    const std::map<int, std::set<int>>& velocity_dof_to_adjacent_cell_dofs,
    const std::map<int, std::set<int>>& cell_dof_to_closure_dofs,
    const std::map<int, int>& velocity_dof_to_component_axis,
    const int seed_velocity_axis)
{
    std::set<int> seed_velocity_dofs, involved_cell_dofs, closure_dofs, initial_velocity_dofs;
    for (std::size_t k = 0; k < overlap_is.size(); ++k)
    {
        seed_velocity_dofs.clear();
        involved_cell_dofs.clear();
        closure_dofs.clear();
        initial_velocity_dofs.clear();
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
        initial_velocity_dofs = seed_velocity_dofs;
        PetscInt first_local_row = -1, row_end = -1;
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
            for (PetscInt col_idx = 0; col_idx < ncols; ++col_idx)
            {
                initial_velocity_dofs.insert(static_cast<int>(cols[col_idx]));
            }
            ierr = MatRestoreRow(A00_mat, row, &ncols, &cols, nullptr);
            IBTK_CHKERRQ(ierr);
        }
        IBAMR::StaggeredStokesPETScMatUtilities::findCoupledCellDofsFromA00(involved_cell_dofs,
                                                                            A00_mat,
                                                                            seed_velocity_dofs,
                                                                            velocity_dof_to_adjacent_cell_dofs,
                                                                            cell_dof_to_closure_dofs);
        for (const int cell_dof : involved_cell_dofs)
        {
            const auto closure_it = cell_dof_to_closure_dofs.find(cell_dof);
            if (closure_it == cell_dof_to_closure_dofs.end()) continue;
            closure_dofs.insert(closure_it->second.begin(), closure_it->second.end());
        }
        overlap_is[k] = closure_dofs;
        overlap_is[k].insert(initial_velocity_dofs.begin(), initial_velocity_dofs.end());
    }
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    int test_failures = 0;

    // Synthetic velocity/cell connectivity for a tiny staggered layout.
    const std::map<int, std::set<int>> velocity_dof_to_adjacent_cell_dofs = {
        { 0, { 10 } }, { 1, { 10, 11 } }, { 2, { 11, 12 } }, { 3, { 12 } }, { 4, { 10 } }, { 5, { 11 } }, { 6, { 12 } }
    };
    const std::map<int, std::set<int>> cell_dof_to_closure_dofs = { { 10, { 0, 1, 4, 10 } },
                                                                    { 11, { 1, 2, 5, 11 } },
                                                                    { 12, { 2, 3, 6, 12 } } };
    const std::map<int, int> velocity_dof_to_component_axis = { { 0, 0 }, { 1, 0 }, { 2, 0 }, { 3, 0 },
                                                                { 4, 1 }, { 5, 1 }, { 6, 1 } };

    // Build a tiny A00 sparsity pattern (velocity-velocity block).
    Mat A00_mat = nullptr;
    int ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 7, 7, 4, nullptr, 4, nullptr, &A00_mat);
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
    add_row(0, { 0, 1, 4 });
    add_row(1, { 0, 1, 2, 5 });
    add_row(2, { 1, 2, 3, 6 });
    add_row(3, { 2, 3 });
    add_row(4, { 0, 4, 5 });
    add_row(5, { 1, 4, 5, 6 });
    add_row(6, { 2, 5, 6 });

    ierr = MatAssemblyBegin(A00_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A00_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);

    // Compare IBAMR overlap construction with MATLAB extract_coupled_dofs-style
    // construction for singleton x-velocity seeds (axis 0).
    for (const int seed_velocity_dof : std::set<int>{ 0, 1, 2, 3 })
    {
        std::vector<std::set<int>> overlap_is(1), nonoverlap_is(1);
        nonoverlap_is[0] = { seed_velocity_dof };
        overlap_is[0] = nonoverlap_is[0];

        construct_coupling_aware_overlap_subdomains_with_cell_closure(overlap_is,
                                                                      nonoverlap_is,
                                                                      A00_mat,
                                                                      velocity_dof_to_adjacent_cell_dofs,
                                                                      cell_dof_to_closure_dofs,
                                                                      velocity_dof_to_component_axis,
                                                                      /*seed_velocity_axis*/ 0);

        const std::set<int> matlab_ground_truth = matlab_extract_coupled_dofs(seed_velocity_dof,
                                                                              A00_mat,
                                                                              velocity_dof_to_adjacent_cell_dofs,
                                                                              cell_dof_to_closure_dofs,
                                                                              velocity_dof_to_component_axis);
        if (!check_set("IBAMR overlap matches MATLAB ground truth for seed " + std::to_string(seed_velocity_dof),
                       matlab_ground_truth,
                       overlap_is[0]))
        {
            ++test_failures;
        }
    }

    ierr = MatDestroy(&A00_mat);
    IBTK_CHKERRQ(ierr);
    std::ofstream out("output");
    out << "test_failures = " << test_failures << "\n";
    return test_failures > 0 ? 1 : 0;
}
