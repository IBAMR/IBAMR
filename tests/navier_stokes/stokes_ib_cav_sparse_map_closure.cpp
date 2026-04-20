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

} // namespace

int
main(int argc, char* argv[])
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app_initializer->getInputDatabase();
    int test_failures = 0;

    const std::unordered_map<int, std::vector<int>> velocity_dof_to_adjacent_cell_dofs = {
        { 0, { 10 } }, { 1, { 10, 11 } }, { 2, { 11 } }, { 3, { 11, 12 } }, { 4, { 12 } }
    };
    const std::unordered_map<int, std::vector<int>> cell_dof_to_closure_dofs = { { 10, { 0, 1, 10 } },
                                                                                 { 11, { 1, 2, 3, 11 } },
                                                                                 { 12, { 3, 4, 12 } } };
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

    ierr = MatDestroy(&SAJ_mat);
    IBTK_CHKERRQ(ierr);
    plog << "Input database:\n";
    input_db->printClassData(plog);
    pout << "test_failures = " << test_failures << std::endl;
    return test_failures > 0 ? 1 : 0;
}
