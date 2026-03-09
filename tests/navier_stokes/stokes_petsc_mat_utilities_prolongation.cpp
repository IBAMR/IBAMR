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

#include <petscmat.h>

#include <SideVariable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <numeric>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

int
main(int argc, char* argv[])
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app_initializer->getInputDatabase();

    int test_failures = 0;
    const auto hierarchy_tuple = setup_hierarchy<NDIM>(app_initializer);
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = std::get<0>(hierarchy_tuple);
    if (patch_hierarchy->getFinestLevelNumber() < 1)
    {
        plog << "Input database:\n";
        input_db->printClassData(plog);
        pout << "test_failures = 1\n";
        return 1;
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("stokes_petsc_mat_utilities_prolongation_ctx");
    Pointer<SideVariable<NDIM, int>> u_dof_index_var = new SideVariable<NDIM, int>("prolong_u_dof");
    Pointer<CellVariable<NDIM, int>> p_dof_index_var = new CellVariable<NDIM, int>("prolong_p_dof");
    const int u_dof_index_idx = var_db->registerVariableAndContext(u_dof_index_var, ctx, IntVector<NDIM>(1));
    const int p_dof_index_idx = var_db->registerVariableAndContext(p_dof_index_var, ctx, IntVector<NDIM>(1));

    const int coarse_ln = 0;
    const int fine_ln = 1;
    Pointer<PatchLevel<NDIM>> coarse_level = patch_hierarchy->getPatchLevel(coarse_ln);
    Pointer<PatchLevel<NDIM>> fine_level = patch_hierarchy->getPatchLevel(fine_ln);
    coarse_level->allocatePatchData(u_dof_index_idx);
    coarse_level->allocatePatchData(p_dof_index_idx);
    fine_level->allocatePatchData(u_dof_index_idx);
    fine_level->allocatePatchData(p_dof_index_idx);

    std::vector<int> num_coarse_dofs_per_proc, num_fine_dofs_per_proc;
    IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
        num_coarse_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, coarse_level);
    IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
        num_fine_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, fine_level);

    AO coarse_level_ao = nullptr;
    int u_coarse_ao_offset = 0, p_coarse_ao_offset = 0;
    IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelAO(coarse_level_ao,
                                                                   num_coarse_dofs_per_proc,
                                                                   u_dof_index_idx,
                                                                   p_dof_index_idx,
                                                                   coarse_level,
                                                                   u_coarse_ao_offset,
                                                                   p_coarse_ao_offset);

    Mat prolong_mat = nullptr;
    IBAMR::StaggeredStokesPETScMatUtilities::constructProlongationOp(prolong_mat,
                                                                     "RT0",
                                                                     "CONSERVATIVE",
                                                                     u_dof_index_idx,
                                                                     p_dof_index_idx,
                                                                     num_fine_dofs_per_proc,
                                                                     num_coarse_dofs_per_proc,
                                                                     fine_level,
                                                                     coarse_level,
                                                                     coarse_level_ao,
                                                                     u_coarse_ao_offset,
                                                                     p_coarse_ao_offset);

    PetscInt nrows = 0, ncols = 0;
    int ierr = MatGetSize(prolong_mat, &nrows, &ncols);
    IBTK_CHKERRQ(ierr);
    const int n_fine = std::accumulate(num_fine_dofs_per_proc.begin(), num_fine_dofs_per_proc.end(), 0);
    const int n_coarse = std::accumulate(num_coarse_dofs_per_proc.begin(), num_coarse_dofs_per_proc.end(), 0);
    if (nrows != static_cast<PetscInt>(n_fine) || ncols != static_cast<PetscInt>(n_coarse)) ++test_failures;

    ierr = MatDestroy(&prolong_mat);
    IBTK_CHKERRQ(ierr);
    ierr = AODestroy(&coarse_level_ao);
    IBTK_CHKERRQ(ierr);

    coarse_level->deallocatePatchData(u_dof_index_idx);
    coarse_level->deallocatePatchData(p_dof_index_idx);
    fine_level->deallocatePatchData(u_dof_index_idx);
    fine_level->deallocatePatchData(p_dof_index_idx);

    plog << "Input database:\n";
    input_db->printClassData(plog);
    pout << "test_failures = " << test_failures << "\n";
    return test_failures > 0 ? 1 : 0;
}
