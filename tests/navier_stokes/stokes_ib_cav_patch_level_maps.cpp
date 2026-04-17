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

#include <CellData.h>
#include <CellGeometry.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>
#include <SideVariable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <iostream>
#include <map>
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
set_from_vector(const std::vector<int>& values)
{
    return std::set<int>(values.begin(), values.end());
}

bool
check_set(const std::string& label, const std::set<int>& expected, const std::vector<int>& actual)
{
    return check_set(label, expected, set_from_vector(actual));
}

void
construct_coupling_aware_overlap_subdomains_with_cell_closure(
    std::vector<std::set<int>>& overlap_is,
    const std::vector<std::set<int>>& nonoverlap_is,
    Mat A00_mat,
    const std::unordered_map<int, std::vector<int>>& velocity_dof_to_adjacent_cell_dofs,
    const std::unordered_map<int, std::vector<int>>& cell_dof_to_closure_dofs,
    const std::unordered_map<int, int>& velocity_dof_to_component_axis,
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
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app_initializer->getInputDatabase();

    int test_failures = 0;
    const auto hierarchy_tuple = setup_hierarchy<NDIM>(app_initializer);
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = std::get<0>(hierarchy_tuple);
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(0);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("stokes_ib_cav_closure_02_context");
    Pointer<SideVariable<NDIM, int>> u_dof_index_var = new SideVariable<NDIM, int>("stokes_ib_cav_closure_02_u_dof");
    Pointer<CellVariable<NDIM, int>> p_dof_index_var = new CellVariable<NDIM, int>("stokes_ib_cav_closure_02_p_dof");
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
    const auto& velocity_dof_to_adjacent_cell_dofs = map_data.velocity_dof_to_adjacent_cell_dofs;
    const auto& cell_dof_to_closure_dofs = map_data.cell_dof_to_closure_dofs;
    const auto& velocity_dof_to_component_axis = map_data.velocity_dof_to_component_axis;

    int seed_velocity_dof = IBTK::invalid_index;
    int lower_cell_dof = IBTK::invalid_index;
    int upper_cell_dof = IBTK::invalid_index;
    Pointer<Patch<NDIM>> seed_patch = nullptr;
    SideIndex<NDIM> seed_side_index(hier::Index<NDIM>(0), 0, SideIndex<NDIM>::Lower);
    for (PatchLevel<NDIM>::Iterator p(level); p && seed_velocity_dof == IBTK::invalid_index; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM, int>> u_dof_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int>> p_dof_data = patch->getPatchData(p_dof_index_idx);
        const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, /*axis*/ 0);
        for (Box<NDIM>::Iterator b(side_box); b && seed_velocity_dof == IBTK::invalid_index; b++)
        {
            const SideIndex<NDIM> i_s(b(), 0, SideIndex<NDIM>::Lower);
            const CellIndex<NDIM> lower_cell = i_s.toCell(0);
            const CellIndex<NDIM> upper_cell = i_s.toCell(1);
            if (!patch_box.contains(lower_cell) || !patch_box.contains(upper_cell)) continue;

            const int u_dof = (*u_dof_data)(i_s);
            const int p_lower = (*p_dof_data)(lower_cell);
            const int p_upper = (*p_dof_data)(upper_cell);
            if (u_dof < 0 || p_lower < 0 || p_upper < 0 || p_lower == p_upper) continue;

            seed_velocity_dof = u_dof;
            lower_cell_dof = p_lower;
            upper_cell_dof = p_upper;
            seed_patch = patch;
            seed_side_index = i_s;
        }
    }

    if (seed_velocity_dof == IBTK::invalid_index)
    {
        std::cerr << "FAILED: could not find an interior axis-0 velocity DOF with two adjacent cells." << std::endl;
        return 1;
    }

    {
        const auto velocity_it = velocity_dof_to_adjacent_cell_dofs.find(seed_velocity_dof);
        if (velocity_it == velocity_dof_to_adjacent_cell_dofs.end())
        {
            std::cerr << "FAILED: missing velocity->adjacent cell entry for seed velocity DOF." << std::endl;
            ++test_failures;
        }
        else if (!check_set("velocity->adjacent cells map",
                            std::set<int>{ lower_cell_dof, upper_cell_dof },
                            velocity_it->second))
        {
            ++test_failures;
        }
    }

    {
        Pointer<SideData<NDIM, int>> u_dof_data = seed_patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int>> p_dof_data = seed_patch->getPatchData(p_dof_index_idx);
        const CellIndex<NDIM> seed_cell = seed_side_index.toCell(0);
        const int seed_cell_dof = (*p_dof_data)(seed_cell);
        std::set<int> expected_cell_closure = { seed_cell_dof };
        for (int axis = 0; axis < NDIM; ++axis)
        {
            const int lower_u = (*u_dof_data)(SideIndex<NDIM>(seed_cell, axis, SideIndex<NDIM>::Lower));
            const int upper_u = (*u_dof_data)(SideIndex<NDIM>(seed_cell, axis, SideIndex<NDIM>::Upper));
            if (lower_u >= 0) expected_cell_closure.insert(lower_u);
            if (upper_u >= 0) expected_cell_closure.insert(upper_u);
        }
        const auto cell_it = cell_dof_to_closure_dofs.find(seed_cell_dof);
        if (cell_it == cell_dof_to_closure_dofs.end())
        {
            std::cerr << "FAILED: missing cell->closure entry for seed cell DOF." << std::endl;
            ++test_failures;
        }
        else if (!check_set("cell->closure map for representative cell", expected_cell_closure, cell_it->second))
        {
            ++test_failures;
        }
    }

    std::map<int, std::pair<Pointer<Patch<NDIM>>, CellIndex<NDIM>>> cell_dof_to_cell_index;
    int max_velocity_dof = IBTK::invalid_index;
    for (const auto& pair : velocity_dof_to_component_axis)
    {
        max_velocity_dof = std::max(max_velocity_dof, pair.first);
    }
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, int>> p_dof_data = patch->getPatchData(p_dof_index_idx);
        for (Box<NDIM>::Iterator b(patch_box); b; b++)
        {
            const CellIndex<NDIM> ic = b();
            const int cell_dof = (*p_dof_data)(ic);
            if (cell_dof >= 0) cell_dof_to_cell_index[cell_dof] = std::make_pair(patch, ic);
        }
    }

    auto make_cell_closure_from_patch_data =
        [u_dof_index_idx, p_dof_index_idx](const std::pair<Pointer<Patch<NDIM>>, CellIndex<NDIM>>& cell_info)
    {
        std::set<int> closure;
        Pointer<Patch<NDIM>> patch = cell_info.first;
        const CellIndex<NDIM>& cell = cell_info.second;
        Pointer<SideData<NDIM, int>> u_dof_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int>> p_dof_data = patch->getPatchData(p_dof_index_idx);
        const int cell_dof = (*p_dof_data)(cell);
        closure.insert(cell_dof);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            const int lower_u = (*u_dof_data)(SideIndex<NDIM>(cell, axis, SideIndex<NDIM>::Lower));
            const int upper_u = (*u_dof_data)(SideIndex<NDIM>(cell, axis, SideIndex<NDIM>::Upper));
            if (lower_u >= 0) closure.insert(lower_u);
            if (upper_u >= 0) closure.insert(upper_u);
        }
        return closure;
    };

    int max_dof = max_velocity_dof;
    for (const auto& pair : cell_dof_to_closure_dofs)
    {
        max_dof = std::max(max_dof, pair.first);
        for (const int dof : pair.second) max_dof = std::max(max_dof, dof);
    }

    Mat SAJ_mat = nullptr;
    const PetscInt n = static_cast<PetscInt>(max_dof + 1);
    int ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, n, n, 1, nullptr, 1, nullptr, &SAJ_mat);
    IBTK_CHKERRQ(ierr);
    const PetscInt row = static_cast<PetscInt>(seed_velocity_dof);
    const PetscInt col = static_cast<PetscInt>(seed_velocity_dof);
    const PetscScalar value = 1.0;
    ierr = MatSetValues(SAJ_mat, 1, &row, 1, &col, &value, INSERT_VALUES);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyBegin(SAJ_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(SAJ_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);

    std::vector<std::set<int>> overlap_is(1), nonoverlap_is(1);
    nonoverlap_is[0] = { seed_velocity_dof };
    overlap_is[0] = nonoverlap_is[0];
    const std::vector<std::set<int>> nonoverlap_expected = nonoverlap_is;

    construct_coupling_aware_overlap_subdomains_with_cell_closure(overlap_is,
                                                                  nonoverlap_is,
                                                                  SAJ_mat,
                                                                  velocity_dof_to_adjacent_cell_dofs,
                                                                  cell_dof_to_closure_dofs,
                                                                  velocity_dof_to_component_axis,
                                                                  /*seed_velocity_axis*/ 0);

    const std::set<int> lower_cell_closure = make_cell_closure_from_patch_data(cell_dof_to_cell_index[lower_cell_dof]);
    const std::set<int> upper_cell_closure = make_cell_closure_from_patch_data(cell_dof_to_cell_index[upper_cell_dof]);
    std::set<int> expected_overlap = nonoverlap_expected[0];
    expected_overlap.insert(lower_cell_closure.begin(), lower_cell_closure.end());
    expected_overlap.insert(upper_cell_closure.begin(), upper_cell_closure.end());

    if (!check_set("patch-level constructed overlap", expected_overlap, overlap_is[0])) ++test_failures;
    if (!check_set("patch-level nonoverlap unchanged", nonoverlap_expected[0], nonoverlap_is[0])) ++test_failures;

    ierr = MatDestroy(&SAJ_mat);
    IBTK_CHKERRQ(ierr);
    level->deallocatePatchData(u_dof_index_idx);
    level->deallocatePatchData(p_dof_index_idx);
    plog << "Input database:\n";
    input_db->printClassData(plog);
    pout << "test_failures = " << test_failures << std::endl;
    return test_failures > 0 ? 1 : 0;
}
