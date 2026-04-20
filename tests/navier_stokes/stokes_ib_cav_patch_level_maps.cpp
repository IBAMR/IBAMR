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

    level->deallocatePatchData(u_dof_index_idx);
    level->deallocatePatchData(p_dof_index_idx);
    plog << "Input database:\n";
    input_db->printClassData(plog);
    pout << "test_failures = " << test_failures << std::endl;
    return test_failures > 0 ? 1 : 0;
}
