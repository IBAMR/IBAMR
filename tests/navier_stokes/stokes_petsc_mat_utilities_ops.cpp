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

#include <CoarseFineBoundary.h>
#include <PoissonSpecifications.h>
#include <SideVariable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <fstream>
#include <set>
#include <string>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

int
main(int argc, char* argv[])
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv);

    int test_failures = 0;
    const auto hierarchy_tuple = setup_hierarchy<NDIM>(app_initializer);
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = std::get<0>(hierarchy_tuple);
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(0);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("stokes_petsc_mat_utilities_ops_ctx");
    Pointer<SideVariable<NDIM, int>> u_dof_index_var = new SideVariable<NDIM, int>("ops_u_dof");
    Pointer<CellVariable<NDIM, int>> p_dof_index_var = new CellVariable<NDIM, int>("ops_p_dof");
    const int u_dof_index_idx = var_db->registerVariableAndContext(u_dof_index_var, ctx, IntVector<NDIM>(1));
    const int p_dof_index_idx = var_db->registerVariableAndContext(p_dof_index_var, ctx, IntVector<NDIM>(1));

    level->allocatePatchData(u_dof_index_idx);
    level->allocatePatchData(p_dof_index_idx);

    std::vector<int> num_dofs_per_proc;
    IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
        num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, level);

    std::vector<std::set<int>> fields;
    std::vector<std::string> field_names;
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelFields(
        fields, field_names, num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, level);
    if (field_names.size() != 2 || fields.size() != 2)
    {
        ++test_failures;
    }

    IBAMR::StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData map_data;
    IBAMR::StaggeredStokesPETScMatUtilities::buildPatchLevelCellClosureMaps(
        map_data, u_dof_index_idx, p_dof_index_idx, level);

    if (map_data.velocity_dof_to_adjacent_cell_dofs.empty() || map_data.cell_dof_to_closure_dofs.empty() ||
        map_data.velocity_dof_to_component_axis.empty())
    {
        ++test_failures;
    }
    int axis0_velocity_count = 0;
    for (const auto& axis_pair : map_data.velocity_dof_to_component_axis)
    {
        if (axis_pair.second == 0) ++axis0_velocity_count;
    }
    if (axis0_velocity_count <= 0) ++test_failures;

    std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);
    PoissonSpecifications u_problem_coefs("stokes_ops_poisson");
    u_problem_coefs.setCConstant(1.0);
    u_problem_coefs.setDConstant(-1.0);
    Mat level_mat = nullptr;
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(
        level_mat, u_problem_coefs, u_bc_coefs, 0.0, num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, level);

    PetscInt nrows = 0, ncols = 0;
    int ierr = MatGetSize(level_mat, &nrows, &ncols);
    IBTK_CHKERRQ(ierr);
    if (nrows <= 0 || ncols <= 0 || nrows != ncols) ++test_failures;

    Mat velocity_submat = nullptr;
    IBAMR::StaggeredStokesPETScMatUtilities::constructA00VelocitySubmatrix(
        velocity_submat, level_mat, num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, level);
    PetscInt vel_rows = 0, vel_cols = 0;
    ierr = MatGetSize(velocity_submat, &vel_rows, &vel_cols);
    IBTK_CHKERRQ(ierr);
    if (vel_rows <= 0 || vel_rows != vel_cols) ++test_failures;

    std::vector<std::set<int>> overlap_geom, nonoverlap_geom;
    const IntVector<NDIM> box_size(2);
    const IntVector<NDIM> overlap_size(1);
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelGeometricalASMSubdomains(
        overlap_geom,
        nonoverlap_geom,
        box_size,
        overlap_size,
        num_dofs_per_proc,
        u_dof_index_idx,
        p_dof_index_idx,
        level,
        Pointer<CoarseFineBoundary<NDIM>>(nullptr));
    if (overlap_geom.empty() || nonoverlap_geom.empty()) ++test_failures;

    std::vector<std::set<int>> overlap_cav, nonoverlap_cav;
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelCouplingAwareASMSubdomains(
        overlap_cav,
        nonoverlap_cav,
        box_size,
        overlap_size,
        num_dofs_per_proc,
        u_dof_index_idx,
        level,
        Pointer<CoarseFineBoundary<NDIM>>(nullptr),
        level_mat,
        map_data,
        0,
        1,
        IBAMR::CouplingAwareASMClosurePolicy::RELAXED);
    if (overlap_cav.empty() || nonoverlap_cav.empty()) ++test_failures;
    if (static_cast<int>(overlap_cav.size()) != axis0_velocity_count) ++test_failures;

    std::vector<std::set<int>> overlap_cav_stride2, nonoverlap_cav_stride2;
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelCouplingAwareASMSubdomains(
        overlap_cav_stride2,
        nonoverlap_cav_stride2,
        box_size,
        overlap_size,
        num_dofs_per_proc,
        u_dof_index_idx,
        level,
        Pointer<CoarseFineBoundary<NDIM>>(nullptr),
        level_mat,
        map_data,
        0,
        2,
        IBAMR::CouplingAwareASMClosurePolicy::RELAXED);
    if (overlap_cav_stride2.empty() || nonoverlap_cav_stride2.empty()) ++test_failures;
    const int expected_stride2_count = (axis0_velocity_count + 1) / 2;
    if (static_cast<int>(overlap_cav_stride2.size()) != expected_stride2_count) ++test_failures;

    ierr = MatDestroy(&velocity_submat);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&level_mat);
    IBTK_CHKERRQ(ierr);

    level->deallocatePatchData(u_dof_index_idx);
    level->deallocatePatchData(p_dof_index_idx);

    std::ofstream out("output");
    out << "test_failures = " << test_failures << "\n";
    return test_failures > 0 ? 1 : 0;
}
