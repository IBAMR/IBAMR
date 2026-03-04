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
#include <ibtk/CartSideDoubleRT0Refine.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>

#include <petscmat.h>

#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <HierarchySideDataOpsReal.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <RefineSchedule.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>
#include <SideVariable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
enum class ProfileType
{
    Affine,
    PiecewiseRT0,
    Nonlinear
};

ProfileType
string_to_profile_type(const std::string& profile)
{
    if (profile == "affine") return ProfileType::Affine;
    if (profile == "piecewise_rt0") return ProfileType::PiecewiseRT0;
    if (profile == "nonlinear") return ProfileType::Nonlinear;
    TBOX_ERROR("Unknown profile_type = " << profile << "\n");
    return ProfileType::Affine;
}

double
compute_side_coordinate(const SideIndex<NDIM>& i_s,
                        const int side_axis,
                        const int coord_axis,
                        const Box<NDIM>& patch_box,
                        const double* x_lower,
                        const double* dx)
{
    const int i = i_s(coord_axis) - patch_box.lower(coord_axis);
    const double offset = (coord_axis == side_axis) ? 0.0 : 0.5;
    return x_lower[coord_axis] + dx[coord_axis] * (static_cast<double>(i) + offset);
}

double
compute_piecewise_rt0_value(const SideIndex<NDIM>& i_s,
                            const int axis,
                            const Box<NDIM>& patch_box,
                            const double* x_lower,
                            const double* dx)
{
    int transverse_sum = 0;
    for (int d = 0; d < NDIM; ++d)
    {
        if (d == axis) continue;
        transverse_sum += (i_s(d) - patch_box.lower(d));
    }
    const int mod3 = ((transverse_sum % 3) + 3) % 3;
    const int mod5 = ((transverse_sum % 5) + 5) % 5;
    const double c0 = 0.2 * static_cast<double>(axis + 1) + 0.04 * static_cast<double>(mod3);
    const double c1 = 0.15 + 0.02 * static_cast<double>(mod5);
    const double x_axis = x_lower[axis] + dx[axis] * static_cast<double>(i_s(axis) - patch_box.lower(axis));
    return c0 + c1 * x_axis;
}

double
compute_nonlinear_value(const SideIndex<NDIM>& i_s,
                        const int axis,
                        const Box<NDIM>& patch_box,
                        const double* x_lower,
                        const double* dx)
{
    const double pi = 3.14159265358979323846;
    const double x_axis = compute_side_coordinate(i_s, axis, axis, patch_box, x_lower, dx);
    double transverse_sum = 0.0;
    double transverse_prod = 1.0;
    for (int d = 0; d < NDIM; ++d)
    {
        if (d == axis) continue;
        const double xd = compute_side_coordinate(i_s, axis, d, patch_box, x_lower, dx);
        transverse_sum += xd;
        transverse_prod *= (xd + 0.15 * static_cast<double>(d + 1));
    }
    return std::sin(2.0 * pi * x_axis) + 0.25 * std::cos(pi * transverse_sum) + 0.12 * x_axis * x_axis +
           0.04 * transverse_prod;
}

void
set_affine_side_field(Pointer<SideData<NDIM, double>> u_data,
                      Pointer<Patch<NDIM>> patch,
                      const std::array<double, NDIM>& coeffs)
{
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();
    for (int axis = 0; axis < NDIM; ++axis)
    {
        const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
        for (Box<NDIM>::Iterator b(side_box); b; b++)
        {
            const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
            const double x_axis = x_lower[axis] + dx[axis] * static_cast<double>(i_s(axis) - patch_box.lower(axis));
            (*u_data)(i_s) = coeffs[axis] * x_axis;
        }
    }
}

void
set_piecewise_rt0_side_field(Pointer<SideData<NDIM, double>> u_data, Pointer<Patch<NDIM>> patch)
{
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();
    for (int axis = 0; axis < NDIM; ++axis)
    {
        const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
        for (Box<NDIM>::Iterator b(side_box); b; b++)
        {
            const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
            (*u_data)(i_s) = compute_piecewise_rt0_value(i_s, axis, patch_box, x_lower, dx);
        }
    }
}

void
set_nonlinear_side_field(Pointer<SideData<NDIM, double>> u_data, Pointer<Patch<NDIM>> patch)
{
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();
    for (int axis = 0; axis < NDIM; ++axis)
    {
        const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
        for (Box<NDIM>::Iterator b(side_box); b; b++)
        {
            const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
            (*u_data)(i_s) = compute_nonlinear_value(i_s, axis, patch_box, x_lower, dx);
        }
    }
}

void
set_test_profile_side_field(Pointer<SideData<NDIM, double>> u_data,
                            Pointer<Patch<NDIM>> patch,
                            const ProfileType profile_type,
                            const std::array<double, NDIM>& coeffs)
{
    if (profile_type == ProfileType::Affine)
    {
        set_affine_side_field(u_data, patch, coeffs);
    }
    else
    {
        if (profile_type == ProfileType::PiecewiseRT0)
        {
            set_piecewise_rt0_side_field(u_data, patch);
        }
        else
        {
            set_nonlinear_side_field(u_data, patch);
        }
    }
}

void
check_nontrivial(const std::string& label,
                 const double samrai_max_norm,
                 const double petsc_max_norm,
                 const double tol,
                 int& test_failures)
{
    if (samrai_max_norm <= tol || petsc_max_norm <= tol)
    {
        ++test_failures;
        pout << "nontriviality check failed for " << label << ": SAMRAI max norm = " << samrai_max_norm
             << ", PETSc max norm = " << petsc_max_norm << ", tolerance = " << tol << "\n";
    }
}

void
check_max_norm_consistency(const std::string& label,
                           const double samrai_max_norm,
                           const double petsc_max_norm,
                           const double tol,
                           int& test_failures)
{
    if (std::abs(samrai_max_norm - petsc_max_norm) > tol)
    {
        ++test_failures;
        pout << "max-norm consistency check failed for " << label << ": SAMRAI max norm = " << samrai_max_norm
             << ", PETSc max norm = " << petsc_max_norm << ", tolerance = " << tol << "\n";
    }
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app_initializer->getInputDatabase();
    Pointer<Database> test_db =
        input_db->keyExists("test") ? input_db->getDatabase("test") : Pointer<Database>(input_db, false);
    const std::string test_mode = test_db->getStringWithDefault("test_mode", "size_check");
    const std::string profile_string = test_db->getStringWithDefault("profile_type", "affine");
    const ProfileType profile_type = string_to_profile_type(profile_string);
    const double equivalence_tol = test_db->getDoubleWithDefault("equivalence_tol", 1.0e-12);
    std::array<double, NDIM> affine_coeffs;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        affine_coeffs[axis] = test_db->getDoubleWithDefault("affine_coefficient_" + std::to_string(axis), 1.0);
    }

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
    HierarchySideDataOpsReal<NDIM, double> coarse_side_ops(patch_hierarchy, coarse_ln, coarse_ln);
    HierarchySideDataOpsReal<NDIM, double> fine_side_ops(patch_hierarchy, fine_ln, fine_ln);
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

    if (test_mode == "rt0_matrix_free_equivalence")
    {
        Pointer<SideVariable<NDIM, double>> u_coarse_var = new SideVariable<NDIM, double>("u_coarse");
        Pointer<SideVariable<NDIM, double>> u_fine_mat_var = new SideVariable<NDIM, double>("u_fine_mat");
        Pointer<SideVariable<NDIM, double>> u_fine_rt0_var = new SideVariable<NDIM, double>("u_fine_rt0");
        Pointer<CellVariable<NDIM, double>> p_coarse_var = new CellVariable<NDIM, double>("p_coarse");
        Pointer<CellVariable<NDIM, double>> p_fine_mat_var = new CellVariable<NDIM, double>("p_fine_mat");
        Pointer<CellVariable<NDIM, double>> p_fine_rt0_var = new CellVariable<NDIM, double>("p_fine_rt0");
        const int u_coarse_idx = var_db->registerVariableAndContext(u_coarse_var, ctx, IntVector<NDIM>(1));
        const int u_fine_mat_idx = var_db->registerVariableAndContext(u_fine_mat_var, ctx, IntVector<NDIM>(1));
        const int u_fine_rt0_idx = var_db->registerVariableAndContext(u_fine_rt0_var, ctx, IntVector<NDIM>(1));
        const int p_coarse_idx = var_db->registerVariableAndContext(p_coarse_var, ctx, IntVector<NDIM>(1));
        const int p_fine_mat_idx = var_db->registerVariableAndContext(p_fine_mat_var, ctx, IntVector<NDIM>(1));
        const int p_fine_rt0_idx = var_db->registerVariableAndContext(p_fine_rt0_var, ctx, IntVector<NDIM>(1));

        coarse_level->allocatePatchData(u_coarse_idx);
        coarse_level->allocatePatchData(p_coarse_idx);
        fine_level->allocatePatchData(u_fine_mat_idx);
        fine_level->allocatePatchData(p_fine_mat_idx);
        fine_level->allocatePatchData(u_fine_rt0_idx);
        fine_level->allocatePatchData(p_fine_rt0_idx);

        TBOX_ASSERT(coarse_level->getNumberOfPatches() == 1);
        TBOX_ASSERT(fine_level->getNumberOfPatches() == 1);
        Pointer<Patch<NDIM>> coarse_patch = coarse_level->getPatch(0);
        Pointer<Patch<NDIM>> fine_patch = fine_level->getPatch(0);
        Pointer<SideData<NDIM, double>> u_coarse_data = coarse_patch->getPatchData(u_coarse_idx);
        Pointer<SideData<NDIM, double>> u_fine_mat_data = fine_patch->getPatchData(u_fine_mat_idx);
        Pointer<SideData<NDIM, double>> u_fine_rt0_data = fine_patch->getPatchData(u_fine_rt0_idx);
        Pointer<CellData<NDIM, double>> p_coarse_data = coarse_patch->getPatchData(p_coarse_idx);
        Pointer<CellData<NDIM, double>> p_fine_mat_data = fine_patch->getPatchData(p_fine_mat_idx);
        Pointer<CellData<NDIM, double>> p_fine_rt0_data = fine_patch->getPatchData(p_fine_rt0_idx);
        u_coarse_data->fillAll(0.0);
        u_fine_mat_data->fillAll(0.0);
        u_fine_rt0_data->fillAll(0.0);
        p_coarse_data->fillAll(0.0);
        p_fine_mat_data->fillAll(0.0);
        p_fine_rt0_data->fillAll(0.0);
        set_test_profile_side_field(u_coarse_data, coarse_patch, profile_type, affine_coeffs);
        const double nontrivial_tol = 1.0e-14;
        const double norm_match_tol = 1.0e-12;

        Vec coarse_vec = nullptr, fine_vec = nullptr, fine_rt0_vec = nullptr, diff_vec = nullptr;
        ierr = MatCreateVecs(prolong_mat, &coarse_vec, &fine_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecDuplicate(fine_vec, &fine_rt0_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecDuplicate(fine_vec, &diff_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecSet(coarse_vec, 0.0);
        IBTK_CHKERRQ(ierr);
        ierr = VecSet(fine_vec, 0.0);
        IBTK_CHKERRQ(ierr);
        ierr = VecSet(fine_rt0_vec, 0.0);
        IBTK_CHKERRQ(ierr);
        ierr = VecSet(diff_vec, 0.0);
        IBTK_CHKERRQ(ierr);

        IBAMR::StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            coarse_vec, u_coarse_idx, u_dof_index_idx, p_coarse_idx, p_dof_index_idx, coarse_level);
        const double coarse_samrai_norm = coarse_side_ops.maxNorm(u_coarse_idx, IBTK::invalid_index);
        double coarse_vec_norm = 0.0;
        ierr = VecNorm(coarse_vec, NORM_INFINITY, &coarse_vec_norm);
        IBTK_CHKERRQ(ierr);
        check_nontrivial("input coarse field", coarse_samrai_norm, coarse_vec_norm, nontrivial_tol, test_failures);
        check_max_norm_consistency(
            "input coarse field", coarse_samrai_norm, coarse_vec_norm, norm_match_tol, test_failures);
        ierr = MatMult(prolong_mat, coarse_vec, fine_vec);
        IBTK_CHKERRQ(ierr);

        Pointer<RefineSchedule<NDIM>> fine_data_synch_sched =
            IBAMR::StaggeredStokesPETScVecUtilities::constructDataSynchSchedule(
                u_fine_mat_idx, p_fine_mat_idx, fine_level);
        Pointer<RefineSchedule<NDIM>> fine_ghost_fill_sched =
            IBAMR::StaggeredStokesPETScVecUtilities::constructGhostFillSchedule(
                u_fine_mat_idx, p_fine_mat_idx, fine_level);
        IBAMR::StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(fine_vec,
                                                                       u_fine_mat_idx,
                                                                       u_dof_index_idx,
                                                                       p_fine_mat_idx,
                                                                       p_dof_index_idx,
                                                                       fine_level,
                                                                       fine_data_synch_sched,
                                                                       fine_ghost_fill_sched);

        const IntVector<NDIM> ratio = fine_level->getRatioToCoarserLevel();
        IBTK::CartSideDoubleRT0Refine refine_op;
        refine_op.refine(*fine_patch, *coarse_patch, u_fine_rt0_idx, u_coarse_idx, fine_patch->getBox(), ratio);
        IBAMR::StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            fine_rt0_vec, u_fine_rt0_idx, u_dof_index_idx, p_fine_rt0_idx, p_dof_index_idx, fine_level);
        const double fine_mat_samrai_norm = fine_side_ops.maxNorm(u_fine_mat_idx, IBTK::invalid_index);
        const double fine_rt0_samrai_norm = fine_side_ops.maxNorm(u_fine_rt0_idx, IBTK::invalid_index);
        double fine_vec_norm = 0.0, fine_rt0_vec_norm = 0.0;
        ierr = VecNorm(fine_vec, NORM_INFINITY, &fine_vec_norm);
        IBTK_CHKERRQ(ierr);
        ierr = VecNorm(fine_rt0_vec, NORM_INFINITY, &fine_rt0_vec_norm);
        IBTK_CHKERRQ(ierr);
        check_nontrivial(
            "output fine matrix field", fine_mat_samrai_norm, fine_vec_norm, nontrivial_tol, test_failures);
        check_nontrivial(
            "output fine matrix-free field", fine_rt0_samrai_norm, fine_rt0_vec_norm, nontrivial_tol, test_failures);
        check_max_norm_consistency(
            "output fine matrix field", fine_mat_samrai_norm, fine_vec_norm, norm_match_tol, test_failures);
        check_max_norm_consistency(
            "output fine matrix-free field", fine_rt0_samrai_norm, fine_rt0_vec_norm, norm_match_tol, test_failures);

        ierr = VecWAXPY(diff_vec, -1.0, fine_rt0_vec, fine_vec);
        IBTK_CHKERRQ(ierr);
        double max_err = 0.0;
        ierr = VecNorm(diff_vec, NORM_INFINITY, &max_err);
        IBTK_CHKERRQ(ierr);
        if (max_err > equivalence_tol) ++test_failures;

        ierr = VecDestroy(&coarse_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&fine_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&fine_rt0_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&diff_vec);
        IBTK_CHKERRQ(ierr);

        fine_level->deallocatePatchData(u_fine_mat_idx);
        fine_level->deallocatePatchData(p_fine_mat_idx);
        fine_level->deallocatePatchData(u_fine_rt0_idx);
        fine_level->deallocatePatchData(p_fine_rt0_idx);
        coarse_level->deallocatePatchData(u_coarse_idx);
        coarse_level->deallocatePatchData(p_coarse_idx);
    }
    else if (test_mode != "size_check")
    {
        TBOX_ERROR("Unknown test_mode = " << test_mode << "\n");
    }

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
