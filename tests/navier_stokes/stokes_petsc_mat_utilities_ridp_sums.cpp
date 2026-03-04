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
#include <ibtk/PETScMatUtilities.h>

#include <petscmat.h>

#include <CellData.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>
#include <SideVariable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <cmath>
#include <fstream>
#include <numeric>
#include <set>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
void
collect_patch_level_dofs(std::set<int>& velocity_dofs,
                         std::set<int>& pressure_dofs,
                         const int u_dof_index_idx,
                         const int p_dof_index_idx,
                         Pointer<PatchLevel<NDIM>> level)
{
    velocity_dofs.clear();
    pressure_dofs.clear();
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM, int>> u_dof_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int>> p_dof_data = patch->getPatchData(p_dof_index_idx);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
            for (Box<NDIM>::Iterator b(side_box); b; b++)
            {
                const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
                const int dof = (*u_dof_data)(i_s);
                if (dof >= 0) velocity_dofs.insert(dof);
            }
        }
        for (Box<NDIM>::Iterator b(patch_box); b; b++)
        {
            const int dof = (*p_dof_data)(b());
            if (dof >= 0) pressure_dofs.insert(dof);
        }
    }
}

void
compute_sum_errors(const Vec sums,
                   const std::set<int>& velocity_dofs,
                   const std::set<int>& pressure_dofs,
                   double& max_vel_error,
                   double& max_pres_error,
                   int& bad_vel_count,
                   int& bad_pres_count,
                   const double tol)
{
    max_vel_error = 0.0;
    max_pres_error = 0.0;
    bad_vel_count = 0;
    bad_pres_count = 0;

    PetscInt ilower = 0, iupper = 0;
    int ierr = VecGetOwnershipRange(sums, &ilower, &iupper);
    IBTK_CHKERRQ(ierr);

    for (PetscInt i = ilower; i < iupper; ++i)
    {
        PetscScalar value = 0.0;
        ierr = VecGetValues(sums, 1, &i, &value);
        IBTK_CHKERRQ(ierr);
        const double v = PetscRealPart(value);
        const auto vel_it = velocity_dofs.find(static_cast<int>(i));
        const auto pres_it = pressure_dofs.find(static_cast<int>(i));
        if (vel_it != velocity_dofs.end())
        {
            const double err = std::abs(v - 1.0);
            max_vel_error = std::max(max_vel_error, err);
            if (err > tol) ++bad_vel_count;
        }
        else if (pres_it != pressure_dofs.end())
        {
            const double err = std::abs(v);
            max_pres_error = std::max(max_pres_error, err);
            if (err > tol) ++bad_pres_count;
        }
    }
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv);

    int test_failures = 0;
    const auto hierarchy_tuple = setup_hierarchy<NDIM>(app_initializer);
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = std::get<0>(hierarchy_tuple);
    if (patch_hierarchy->getFinestLevelNumber() < 1)
    {
        std::ofstream out("output");
        out << "test_failures = 1\n";
        return 1;
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("stokes_petsc_mat_utilities_ridp_sums_ctx");
    Pointer<SideVariable<NDIM, int>> u_dof_index_var = new SideVariable<NDIM, int>("ridp_u_dof");
    Pointer<CellVariable<NDIM, int>> p_dof_index_var = new CellVariable<NDIM, int>("ridp_p_dof");
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

    std::set<int> fine_velocity_dofs, fine_pressure_dofs;
    std::set<int> coarse_velocity_dofs, coarse_pressure_dofs;
    collect_patch_level_dofs(fine_velocity_dofs, fine_pressure_dofs, u_dof_index_idx, p_dof_index_idx, fine_level);
    collect_patch_level_dofs(
        coarse_velocity_dofs, coarse_pressure_dofs, u_dof_index_idx, p_dof_index_idx, coarse_level);

    Mat id_vel_mat = nullptr;
    ierr = MatCreateAIJ(PETSC_COMM_WORLD,
                        PETSC_DECIDE,
                        PETSC_DECIDE,
                        static_cast<PetscInt>(n_fine),
                        static_cast<PetscInt>(n_fine),
                        1,
                        nullptr,
                        0,
                        nullptr,
                        &id_vel_mat);
    IBTK_CHKERRQ(ierr);
    PetscInt fine_ilower = 0, fine_iupper = 0;
    ierr = MatGetOwnershipRange(id_vel_mat, &fine_ilower, &fine_iupper);
    IBTK_CHKERRQ(ierr);
    for (PetscInt i = fine_ilower; i < fine_iupper; ++i)
    {
        const double diag = fine_velocity_dofs.count(static_cast<int>(i)) ? 1.0 : 0.0;
        ierr = MatSetValue(id_vel_mat, i, i, diag, INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(id_vel_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(id_vel_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);

    Vec restriction_scale = nullptr;
    PETScMatUtilities::constructRestrictionScalingOp(prolong_mat, restriction_scale);

    Mat ridp_mat = nullptr;
    ierr = MatPtAP(id_vel_mat, prolong_mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &ridp_mat);
    IBTK_CHKERRQ(ierr);
    ierr = MatDiagonalScale(ridp_mat, restriction_scale, nullptr);
    IBTK_CHKERRQ(ierr);

    Vec ones = nullptr, row_sums = nullptr, col_sums = nullptr;
    ierr = MatCreateVecs(ridp_mat, &ones, &row_sums);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(ones, &col_sums);
    IBTK_CHKERRQ(ierr);
    ierr = VecSet(ones, 1.0);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(ridp_mat, ones, row_sums);
    IBTK_CHKERRQ(ierr);
    ierr = MatMultTranspose(ridp_mat, ones, col_sums);
    IBTK_CHKERRQ(ierr);

    const double tol = 1.0e-12;
    double row_vel_max_err = 0.0, row_pres_max_err = 0.0;
    double col_vel_max_err = 0.0, col_pres_max_err = 0.0;
    int row_bad_vel_count = 0, row_bad_pres_count = 0;
    int col_bad_vel_count = 0, col_bad_pres_count = 0;
    compute_sum_errors(row_sums,
                       coarse_velocity_dofs,
                       coarse_pressure_dofs,
                       row_vel_max_err,
                       row_pres_max_err,
                       row_bad_vel_count,
                       row_bad_pres_count,
                       tol);
    compute_sum_errors(col_sums,
                       coarse_velocity_dofs,
                       coarse_pressure_dofs,
                       col_vel_max_err,
                       col_pres_max_err,
                       col_bad_vel_count,
                       col_bad_pres_count,
                       tol);

    if (row_bad_vel_count > 0 || row_bad_pres_count > 0 || col_bad_vel_count > 0 || col_bad_pres_count > 0)
    {
        ++test_failures;
    }

    ierr = MatDestroy(&prolong_mat);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&id_vel_mat);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&ridp_mat);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&restriction_scale);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&ones);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&row_sums);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&col_sums);
    IBTK_CHKERRQ(ierr);
    ierr = AODestroy(&coarse_level_ao);
    IBTK_CHKERRQ(ierr);

    coarse_level->deallocatePatchData(u_dof_index_idx);
    coarse_level->deallocatePatchData(p_dof_index_idx);
    fine_level->deallocatePatchData(u_dof_index_idx);
    fine_level->deallocatePatchData(p_dof_index_idx);

    std::ofstream out("output");
    out << "n_fine_dofs = " << n_fine << "\n";
    out << "n_coarse_dofs = " << n_coarse << "\n";
    out << "n_fine_velocity_dofs = " << fine_velocity_dofs.size() << "\n";
    out << "n_fine_pressure_dofs = " << fine_pressure_dofs.size() << "\n";
    out << "n_coarse_velocity_dofs = " << coarse_velocity_dofs.size() << "\n";
    out << "n_coarse_pressure_dofs = " << coarse_pressure_dofs.size() << "\n";
    out << "row_bad_velocity_count = " << row_bad_vel_count << "\n";
    out << "row_bad_pressure_count = " << row_bad_pres_count << "\n";
    out << "col_bad_velocity_count = " << col_bad_vel_count << "\n";
    out << "col_bad_pressure_count = " << col_bad_pres_count << "\n";
    out << "row_velocity_max_error = " << row_vel_max_err << "\n";
    out << "row_pressure_max_error = " << row_pres_max_err << "\n";
    out << "col_velocity_max_error = " << col_vel_max_err << "\n";
    out << "col_pressure_max_error = " << col_pres_max_err << "\n";
    out << "test_failures = " << test_failures << "\n";
    return test_failures > 0 ? 1 : 0;
}
