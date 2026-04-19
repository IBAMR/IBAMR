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

/*
 * This test checks linear cell-centered pressure prolongation and restriction
 * parity between the PETSc transfer operators and the corresponding SAMRAI
 * transfer schedules on periodic analytic profiles.
 */

#include <ibamr/StaggeredStokesPETScMatUtilities.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartCellDoubleLinearCoarsen.h>
#include <ibtk/CartCellDoubleLinearRefine.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/ibtk_utilities.h>

#include <petscmat.h>

#include <boost/math/constants/constants.hpp>

#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <CellGeometry.h>
#include <CoarsenAlgorithm.h>
#include <CoarsenSchedule.h>
#include <HierarchyCellDataOpsReal.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RefineAlgorithm.h>
#include <RefineOperator.h>
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

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
enum class ProfileType
{
    TRIGONOMETRIC,
    NONLINEAR,
    UNKNOWN
};

ProfileType
string_to_profile_type(const std::string& profile)
{
    if (profile == "trigonometric") return ProfileType::TRIGONOMETRIC;
    if (profile == "nonlinear") return ProfileType::NONLINEAR;
    if (profile == "unknown") return ProfileType::UNKNOWN;
    TBOX_ERROR("Unknown profile_type = " << profile << "\n");
    return ProfileType::UNKNOWN;
}

double
compute_trigonometric_cell_value(const VectorNd& X, const std::array<double, NDIM>& coeffs, const double constant)
{
    const double pi = boost::math::constants::pi<double>();
    double val = constant;
    for (int d = 0; d < NDIM; ++d)
    {
        val += coeffs[d] * std::sin(2.0 * pi * (d + 1.0) * X[d]);
    }
    return val;
}

double
compute_nonlinear_cell_value(const VectorNd& X)
{
    const double pi = boost::math::constants::pi<double>();
    double phase_1 = 0.0;
    double phase_2 = 0.0;
    double prod = 1.0;
    for (int d = 0; d < NDIM; ++d)
    {
        phase_1 += (d + 1.0) * X[d];
        phase_2 += (d + 2.0) * X[d];
        prod *= std::sin(2.0 * pi * (d + 1.0) * X[d] + 0.15 * pi);
    }
    return 0.35 * std::sin(2.0 * pi * phase_1) + 0.2 * std::cos(2.0 * pi * phase_2) + 0.1 * prod;
}

void
set_pressure_profile(Pointer<CellData<NDIM, double>> p_data,
                     Pointer<Patch<NDIM>> patch,
                     const ProfileType profile_type,
                     const std::array<double, NDIM>& profile_coeffs,
                     const double profile_constant)
{
    const Box<NDIM>& fill_box = p_data->getGhostBox();

    for (Box<NDIM>::Iterator b(fill_box); b; b++)
    {
        const CellIndex<NDIM>& i = b();
        const VectorNd X = IndexUtilities::getCellCenter(*patch, i);
        double val = 0.0;
        switch (profile_type)
        {
        case ProfileType::TRIGONOMETRIC:
            val = compute_trigonometric_cell_value(X, profile_coeffs, profile_constant);
            break;
        case ProfileType::NONLINEAR:
            val = compute_nonlinear_cell_value(X);
            break;
        case ProfileType::UNKNOWN:
        default:
            TBOX_ERROR("Unknown profile_type.\n");
        }
        (*p_data)(i) = val;
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
    if (!IBTK::rel_equal_eps(samrai_max_norm, petsc_max_norm, tol))
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

    const std::string test_mode = test_db->getStringWithDefault("test_mode", "prolongation");
    const std::string profile_string = test_db->getStringWithDefault("profile_type", "unknown");
    const ProfileType profile_type = string_to_profile_type(profile_string);
    const double equivalence_tol = test_db->getDoubleWithDefault("equivalence_tol", 1.0e-12);
    double profile_constant = 0.0;
    std::array<double, NDIM> profile_coeffs = {};
    if (profile_type == ProfileType::UNKNOWN)
    {
        TBOX_ERROR("Unknown profile_type = " << profile_string << "\n");
    }
    else if (profile_type == ProfileType::TRIGONOMETRIC)
    {
        profile_constant = test_db->keyExists("profile_constant") ? test_db->getDouble("profile_constant") :
                           test_db->keyExists("affine_constant")  ? test_db->getDouble("affine_constant") :
                                                                    0.7;
        for (int axis = 0; axis < NDIM; ++axis)
        {
            const std::string coeff_key = "profile_coefficient_" + std::to_string(axis);
            const std::string affine_coeff_key = "affine_coefficient_" + std::to_string(axis);
            profile_coeffs[axis] = test_db->keyExists(coeff_key)        ? test_db->getDouble(coeff_key) :
                                   test_db->keyExists(affine_coeff_key) ? test_db->getDouble(affine_coeff_key) :
                                                                          1.0;
        }
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

    const int coarse_ln = 0;
    const int fine_ln = 1;
    Pointer<PatchLevel<NDIM>> coarse_level = patch_hierarchy->getPatchLevel(coarse_ln);
    Pointer<PatchLevel<NDIM>> fine_level = patch_hierarchy->getPatchLevel(fine_ln);
    HierarchyCellDataOpsReal<NDIM, double> coarse_cell_ops(patch_hierarchy, coarse_ln, coarse_ln);
    HierarchyCellDataOpsReal<NDIM, double> fine_cell_ops(patch_hierarchy, fine_ln, fine_ln);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("stokes_petsc_mat_utilities_pressure_transfer_ctx");

    Pointer<SideVariable<NDIM, int>> u_dof_index_var = new SideVariable<NDIM, int>("pressure_transfer_u_dof");
    Pointer<CellVariable<NDIM, int>> p_dof_index_var = new CellVariable<NDIM, int>("pressure_transfer_p_dof");
    const int u_dof_index_idx = var_db->registerVariableAndContext(u_dof_index_var, ctx, IntVector<NDIM>(1));
    const int p_dof_index_idx = var_db->registerVariableAndContext(p_dof_index_var, ctx, IntVector<NDIM>(1));

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
                                                                     "LINEAR",
                                                                     "LINEAR",
                                                                     u_dof_index_idx,
                                                                     p_dof_index_idx,
                                                                     num_fine_dofs_per_proc,
                                                                     num_coarse_dofs_per_proc,
                                                                     fine_level,
                                                                     coarse_level,
                                                                     coarse_level_ao,
                                                                     u_coarse_ao_offset,
                                                                     p_coarse_ao_offset);

    if (test_mode == "prolongation")
    {
        Pointer<SideVariable<NDIM, double>> u_coarse_var = new SideVariable<NDIM, double>("u_coarse");
        Pointer<SideVariable<NDIM, double>> u_fine_mat_var = new SideVariable<NDIM, double>("u_fine_mat");
        Pointer<SideVariable<NDIM, double>> u_fine_mf_var = new SideVariable<NDIM, double>("u_fine_mf");
        Pointer<CellVariable<NDIM, double>> p_coarse_var = new CellVariable<NDIM, double>("p_coarse");
        Pointer<CellVariable<NDIM, double>> p_fine_mat_var = new CellVariable<NDIM, double>("p_fine_mat");
        Pointer<CellVariable<NDIM, double>> p_fine_mf_var = new CellVariable<NDIM, double>("p_fine_mf");

        const int u_coarse_idx = var_db->registerVariableAndContext(u_coarse_var, ctx, IntVector<NDIM>(1));
        const int u_fine_mat_idx = var_db->registerVariableAndContext(u_fine_mat_var, ctx, IntVector<NDIM>(1));
        const int u_fine_mf_idx = var_db->registerVariableAndContext(u_fine_mf_var, ctx, IntVector<NDIM>(1));
        const int p_coarse_idx = var_db->registerVariableAndContext(p_coarse_var, ctx, IntVector<NDIM>(1));
        const int p_fine_mat_idx = var_db->registerVariableAndContext(p_fine_mat_var, ctx, IntVector<NDIM>(1));
        const int p_fine_mf_idx = var_db->registerVariableAndContext(p_fine_mf_var, ctx, IntVector<NDIM>(1));

        coarse_level->allocatePatchData(u_coarse_idx);
        coarse_level->allocatePatchData(p_coarse_idx);
        fine_level->allocatePatchData(u_fine_mat_idx);
        fine_level->allocatePatchData(p_fine_mat_idx);
        fine_level->allocatePatchData(u_fine_mf_idx);
        fine_level->allocatePatchData(p_fine_mf_idx);

        TBOX_ASSERT(coarse_level->getNumberOfPatches() == 1);
        TBOX_ASSERT(fine_level->getNumberOfPatches() == 1);
        Pointer<Patch<NDIM>> coarse_patch = coarse_level->getPatch(0);
        Pointer<Patch<NDIM>> fine_patch = fine_level->getPatch(0);

        Pointer<SideData<NDIM, double>> u_coarse_data = coarse_patch->getPatchData(u_coarse_idx);
        Pointer<SideData<NDIM, double>> u_fine_mat_data = fine_patch->getPatchData(u_fine_mat_idx);
        Pointer<SideData<NDIM, double>> u_fine_mf_data = fine_patch->getPatchData(u_fine_mf_idx);
        Pointer<CellData<NDIM, double>> p_coarse_data = coarse_patch->getPatchData(p_coarse_idx);
        Pointer<CellData<NDIM, double>> p_fine_mat_data = fine_patch->getPatchData(p_fine_mat_idx);
        Pointer<CellData<NDIM, double>> p_fine_mf_data = fine_patch->getPatchData(p_fine_mf_idx);

        u_coarse_data->fillAll(0.0);
        u_fine_mat_data->fillAll(0.0);
        u_fine_mf_data->fillAll(0.0);
        p_fine_mat_data->fillAll(0.0);
        p_fine_mf_data->fillAll(0.0);
        set_pressure_profile(p_coarse_data, coarse_patch, profile_type, profile_coeffs, profile_constant);
        const double nontrivial_tol = 1.0e-14;
        const double norm_match_tol = 1.0e-12;

        Vec coarse_vec = nullptr, fine_vec_mat = nullptr, fine_vec_mf = nullptr, diff_vec = nullptr;
        int ierr = MatCreateVecs(prolong_mat, &coarse_vec, &fine_vec_mat);
        IBTK_CHKERRQ(ierr);
        ierr = VecDuplicate(fine_vec_mat, &fine_vec_mf);
        IBTK_CHKERRQ(ierr);
        ierr = VecDuplicate(fine_vec_mat, &diff_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecSet(fine_vec_mf, 0.0);
        IBTK_CHKERRQ(ierr);

        IBAMR::StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            coarse_vec, u_coarse_idx, u_dof_index_idx, p_coarse_idx, p_dof_index_idx, coarse_level);
        const double coarse_samrai_norm = coarse_cell_ops.maxNorm(p_coarse_idx, IBTK::invalid_index);
        double coarse_vec_norm = 0.0;
        ierr = VecNorm(coarse_vec, NORM_INFINITY, &coarse_vec_norm);
        IBTK_CHKERRQ(ierr);
        check_nontrivial(
            "input coarse pressure field", coarse_samrai_norm, coarse_vec_norm, nontrivial_tol, test_failures);
        check_max_norm_consistency(
            "input coarse pressure field", coarse_samrai_norm, coarse_vec_norm, norm_match_tol, test_failures);
        ierr = MatMult(prolong_mat, coarse_vec, fine_vec_mat);
        IBTK_CHKERRQ(ierr);

        Pointer<RefineSchedule<NDIM>> fine_data_synch_sched =
            IBAMR::StaggeredStokesPETScVecUtilities::constructDataSynchSchedule(
                u_fine_mat_idx, p_fine_mat_idx, fine_level);
        Pointer<RefineSchedule<NDIM>> fine_ghost_fill_sched =
            IBAMR::StaggeredStokesPETScVecUtilities::constructGhostFillSchedule(
                u_fine_mat_idx, p_fine_mat_idx, fine_level);
        IBAMR::StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(fine_vec_mat,
                                                                       u_fine_mat_idx,
                                                                       u_dof_index_idx,
                                                                       p_fine_mat_idx,
                                                                       p_dof_index_idx,
                                                                       fine_level,
                                                                       fine_data_synch_sched,
                                                                       fine_ghost_fill_sched);

        Pointer<RefineOperator<NDIM>> p_refine_op = new IBTK::CartCellDoubleLinearRefine();
        Pointer<RefineAlgorithm<NDIM>> p_refine_alg = new RefineAlgorithm<NDIM>();
        p_refine_alg->registerRefine(p_fine_mf_idx, p_coarse_idx, p_fine_mf_idx, p_refine_op, nullptr);
        Pointer<RefineSchedule<NDIM>> p_refine_sched =
            p_refine_alg->createSchedule(fine_level, Pointer<PatchLevel<NDIM>>(), coarse_ln, patch_hierarchy, nullptr);
        p_refine_sched->fillData(0.0);

        IBAMR::StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            fine_vec_mf, u_fine_mf_idx, u_dof_index_idx, p_fine_mf_idx, p_dof_index_idx, fine_level);
        const double fine_mat_samrai_norm = fine_cell_ops.maxNorm(p_fine_mat_idx, IBTK::invalid_index);
        const double fine_mf_samrai_norm = fine_cell_ops.maxNorm(p_fine_mf_idx, IBTK::invalid_index);
        double fine_vec_mat_norm = 0.0, fine_vec_mf_norm = 0.0;
        ierr = VecNorm(fine_vec_mat, NORM_INFINITY, &fine_vec_mat_norm);
        IBTK_CHKERRQ(ierr);
        ierr = VecNorm(fine_vec_mf, NORM_INFINITY, &fine_vec_mf_norm);
        IBTK_CHKERRQ(ierr);
        check_nontrivial("output fine matrix pressure field",
                         fine_mat_samrai_norm,
                         fine_vec_mat_norm,
                         nontrivial_tol,
                         test_failures);
        check_nontrivial("output fine matrix-free pressure field",
                         fine_mf_samrai_norm,
                         fine_vec_mf_norm,
                         nontrivial_tol,
                         test_failures);
        check_max_norm_consistency("output fine matrix pressure field",
                                   fine_mat_samrai_norm,
                                   fine_vec_mat_norm,
                                   norm_match_tol,
                                   test_failures);
        check_max_norm_consistency("output fine matrix-free pressure field",
                                   fine_mf_samrai_norm,
                                   fine_vec_mf_norm,
                                   norm_match_tol,
                                   test_failures);

        ierr = VecWAXPY(diff_vec, -1.0, fine_vec_mf, fine_vec_mat);
        IBTK_CHKERRQ(ierr);
        double max_err = 0.0;
        ierr = VecNorm(diff_vec, NORM_INFINITY, &max_err);
        IBTK_CHKERRQ(ierr);
        if (max_err > equivalence_tol) ++test_failures;

        ierr = VecDestroy(&coarse_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&fine_vec_mat);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&fine_vec_mf);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&diff_vec);
        IBTK_CHKERRQ(ierr);

        coarse_level->deallocatePatchData(u_coarse_idx);
        coarse_level->deallocatePatchData(p_coarse_idx);
        fine_level->deallocatePatchData(u_fine_mat_idx);
        fine_level->deallocatePatchData(p_fine_mat_idx);
        fine_level->deallocatePatchData(u_fine_mf_idx);
        fine_level->deallocatePatchData(p_fine_mf_idx);
    }
    else if (test_mode == "restriction")
    {
        Pointer<SideVariable<NDIM, double>> u_fine_var = new SideVariable<NDIM, double>("u_fine");
        Pointer<SideVariable<NDIM, double>> u_coarse_mat_var = new SideVariable<NDIM, double>("u_coarse_mat");
        Pointer<SideVariable<NDIM, double>> u_coarse_mf_var = new SideVariable<NDIM, double>("u_coarse_mf");
        Pointer<CellVariable<NDIM, double>> p_fine_var = new CellVariable<NDIM, double>("p_fine");
        Pointer<CellVariable<NDIM, double>> p_coarse_mat_var = new CellVariable<NDIM, double>("p_coarse_mat");
        Pointer<CellVariable<NDIM, double>> p_coarse_mf_var = new CellVariable<NDIM, double>("p_coarse_mf");

        const int u_fine_idx = var_db->registerVariableAndContext(u_fine_var, ctx, IntVector<NDIM>(1));
        const int u_coarse_mat_idx = var_db->registerVariableAndContext(u_coarse_mat_var, ctx, IntVector<NDIM>(1));
        const int u_coarse_mf_idx = var_db->registerVariableAndContext(u_coarse_mf_var, ctx, IntVector<NDIM>(1));
        const int p_fine_idx = var_db->registerVariableAndContext(p_fine_var, ctx, IntVector<NDIM>(1));
        const int p_coarse_mat_idx = var_db->registerVariableAndContext(p_coarse_mat_var, ctx, IntVector<NDIM>(1));
        const int p_coarse_mf_idx = var_db->registerVariableAndContext(p_coarse_mf_var, ctx, IntVector<NDIM>(1));

        fine_level->allocatePatchData(u_fine_idx);
        fine_level->allocatePatchData(p_fine_idx);
        coarse_level->allocatePatchData(u_coarse_mat_idx);
        coarse_level->allocatePatchData(p_coarse_mat_idx);
        coarse_level->allocatePatchData(u_coarse_mf_idx);
        coarse_level->allocatePatchData(p_coarse_mf_idx);

        TBOX_ASSERT(coarse_level->getNumberOfPatches() == 1);
        TBOX_ASSERT(fine_level->getNumberOfPatches() == 1);
        Pointer<Patch<NDIM>> coarse_patch = coarse_level->getPatch(0);
        Pointer<Patch<NDIM>> fine_patch = fine_level->getPatch(0);

        Pointer<SideData<NDIM, double>> u_fine_data = fine_patch->getPatchData(u_fine_idx);
        Pointer<SideData<NDIM, double>> u_coarse_mat_data = coarse_patch->getPatchData(u_coarse_mat_idx);
        Pointer<SideData<NDIM, double>> u_coarse_mf_data = coarse_patch->getPatchData(u_coarse_mf_idx);
        Pointer<CellData<NDIM, double>> p_fine_data = fine_patch->getPatchData(p_fine_idx);
        Pointer<CellData<NDIM, double>> p_coarse_mat_data = coarse_patch->getPatchData(p_coarse_mat_idx);
        Pointer<CellData<NDIM, double>> p_coarse_mf_data = coarse_patch->getPatchData(p_coarse_mf_idx);

        u_fine_data->fillAll(0.0);
        u_coarse_mat_data->fillAll(0.0);
        u_coarse_mf_data->fillAll(0.0);
        p_coarse_mat_data->fillAll(0.0);
        p_coarse_mf_data->fillAll(0.0);
        set_pressure_profile(p_fine_data, fine_patch, profile_type, profile_coeffs, profile_constant);
        const double nontrivial_tol = 1.0e-14;
        const double norm_match_tol = 1.0e-12;

        Vec fine_vec = nullptr, coarse_vec_mat = nullptr, coarse_vec_mf = nullptr, diff_vec = nullptr;
        int ierr = MatCreateVecs(prolong_mat, &coarse_vec_mat, &fine_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecDuplicate(coarse_vec_mat, &coarse_vec_mf);
        IBTK_CHKERRQ(ierr);
        ierr = VecDuplicate(coarse_vec_mat, &diff_vec);
        IBTK_CHKERRQ(ierr);

        IBAMR::StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            fine_vec, u_fine_idx, u_dof_index_idx, p_fine_idx, p_dof_index_idx, fine_level);
        const double fine_samrai_norm = fine_cell_ops.maxNorm(p_fine_idx, IBTK::invalid_index);
        double fine_vec_norm = 0.0;
        ierr = VecNorm(fine_vec, NORM_INFINITY, &fine_vec_norm);
        IBTK_CHKERRQ(ierr);
        check_nontrivial("input fine pressure field", fine_samrai_norm, fine_vec_norm, nontrivial_tol, test_failures);
        check_max_norm_consistency(
            "input fine pressure field", fine_samrai_norm, fine_vec_norm, norm_match_tol, test_failures);

        Vec restriction_scale = nullptr;
        PETScMatUtilities::constructRestrictionScalingOp(prolong_mat, restriction_scale);
        Mat restrict_mat = nullptr;
        ierr = MatTranspose(prolong_mat, MAT_INITIAL_MATRIX, &restrict_mat);
        IBTK_CHKERRQ(ierr);
        ierr = MatDiagonalScale(restrict_mat, restriction_scale, nullptr);
        IBTK_CHKERRQ(ierr);
        ierr = MatMult(restrict_mat, fine_vec, coarse_vec_mat);
        IBTK_CHKERRQ(ierr);

        Pointer<RefineSchedule<NDIM>> coarse_data_synch_sched =
            IBAMR::StaggeredStokesPETScVecUtilities::constructDataSynchSchedule(
                u_coarse_mat_idx, p_coarse_mat_idx, coarse_level);
        Pointer<RefineSchedule<NDIM>> coarse_ghost_fill_sched =
            IBAMR::StaggeredStokesPETScVecUtilities::constructGhostFillSchedule(
                u_coarse_mat_idx, p_coarse_mat_idx, coarse_level);
        IBAMR::StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(coarse_vec_mat,
                                                                       u_coarse_mat_idx,
                                                                       u_dof_index_idx,
                                                                       p_coarse_mat_idx,
                                                                       p_dof_index_idx,
                                                                       coarse_level,
                                                                       coarse_data_synch_sched,
                                                                       coarse_ghost_fill_sched);

        Pointer<CoarsenOperator<NDIM>> p_coarsen_op = new IBTK::CartCellDoubleLinearCoarsen();
        Pointer<CoarsenAlgorithm<NDIM>> p_coarsen_alg = new CoarsenAlgorithm<NDIM>();
        p_coarsen_alg->registerCoarsen(p_coarse_mf_idx, p_fine_idx, p_coarsen_op, IntVector<NDIM>(0));
        Pointer<CoarsenSchedule<NDIM>> p_coarsen_sched = p_coarsen_alg->createSchedule(coarse_level, fine_level);
        p_coarsen_sched->coarsenData();

        IBAMR::StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            coarse_vec_mf, u_coarse_mf_idx, u_dof_index_idx, p_coarse_mf_idx, p_dof_index_idx, coarse_level);
        const double coarse_mat_samrai_norm = coarse_cell_ops.maxNorm(p_coarse_mat_idx, IBTK::invalid_index);
        const double coarse_mf_samrai_norm = coarse_cell_ops.maxNorm(p_coarse_mf_idx, IBTK::invalid_index);
        double coarse_vec_mat_norm = 0.0, coarse_vec_mf_norm = 0.0;
        ierr = VecNorm(coarse_vec_mat, NORM_INFINITY, &coarse_vec_mat_norm);
        IBTK_CHKERRQ(ierr);
        ierr = VecNorm(coarse_vec_mf, NORM_INFINITY, &coarse_vec_mf_norm);
        IBTK_CHKERRQ(ierr);
        check_nontrivial("output coarse matrix pressure field",
                         coarse_mat_samrai_norm,
                         coarse_vec_mat_norm,
                         nontrivial_tol,
                         test_failures);
        check_nontrivial("output coarse matrix-free pressure field",
                         coarse_mf_samrai_norm,
                         coarse_vec_mf_norm,
                         nontrivial_tol,
                         test_failures);
        check_max_norm_consistency("output coarse matrix pressure field",
                                   coarse_mat_samrai_norm,
                                   coarse_vec_mat_norm,
                                   norm_match_tol,
                                   test_failures);
        check_max_norm_consistency("output coarse matrix-free pressure field",
                                   coarse_mf_samrai_norm,
                                   coarse_vec_mf_norm,
                                   norm_match_tol,
                                   test_failures);

        ierr = VecWAXPY(diff_vec, -1.0, coarse_vec_mf, coarse_vec_mat);
        IBTK_CHKERRQ(ierr);
        double max_err = 0.0;
        ierr = VecNorm(diff_vec, NORM_INFINITY, &max_err);
        IBTK_CHKERRQ(ierr);
        if (max_err > equivalence_tol) ++test_failures;

        ierr = MatDestroy(&restrict_mat);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&restriction_scale);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&fine_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&coarse_vec_mat);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&coarse_vec_mf);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&diff_vec);
        IBTK_CHKERRQ(ierr);

        fine_level->deallocatePatchData(u_fine_idx);
        fine_level->deallocatePatchData(p_fine_idx);
        coarse_level->deallocatePatchData(u_coarse_mat_idx);
        coarse_level->deallocatePatchData(p_coarse_mat_idx);
        coarse_level->deallocatePatchData(u_coarse_mf_idx);
        coarse_level->deallocatePatchData(p_coarse_mf_idx);
    }
    else
    {
        TBOX_ERROR("Unknown test_mode = " << test_mode << "\n");
    }

    int ierr = MatDestroy(&prolong_mat);
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
