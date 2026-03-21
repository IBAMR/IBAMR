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

#include <ibamr/StaggeredStokesPETScLevelSolver.h>
#include <ibamr/StaggeredStokesPETScMatUtilities.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>

#include <tbox/MemoryDatabase.h>

#include <petscksp.h>
#include <petscsys.h>

#include <CellData.h>
#include <CellVariable.h>
#include <PoissonSpecifications.h>
#include <SAMRAIVectorReal.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
void
fill_nontrivial_rhs(const Pointer<PatchLevel<NDIM>>& level, const int f_u_idx, const int f_p_idx)
{
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> f_u_data = patch->getPatchData(f_u_idx);
        Pointer<CellData<NDIM, double>> f_p_data = patch->getPatchData(f_p_idx);
        const Box<NDIM> patch_box = patch->getBox();
        for (int axis = 0; axis < NDIM; ++axis)
        {
            const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
            for (Box<NDIM>::Iterator b(side_box); b; b++)
            {
                const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
                const double val =
                    0.05 * (axis + 1) + 0.01 * static_cast<double>(i_s(0)) - 0.02 * static_cast<double>(i_s(1));
                (*f_u_data)(i_s) = val;
            }
        }
        for (Box<NDIM>::Iterator b(patch_box); b; b++)
        {
            const CellIndex<NDIM>& i_c = b();
            (*f_p_data)(i_c) = 0.03 * static_cast<double>(i_c(0)) + 0.04 * static_cast<double>(i_c(1)) - 0.1;
        }
    }
}

void
zero_solution_fields(const Pointer<PatchLevel<NDIM>>& level, const int u_idx, const int p_idx)
{
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
        u_data->fillAll(0.0);
        p_data->fillAll(0.0);
    }
}

void
apply_matlab_cav_sweep(Vec y,
                       Vec b,
                       Mat A,
                       const std::vector<std::vector<int>>& subdomain_dofs,
                       const double alpha,
                       IS pressure_is)
{
    int ierr = VecZeroEntries(y);
    IBTK_CHKERRQ(ierr);
    Vec r = nullptr;
    ierr = VecDuplicate(y, &r);
    IBTK_CHKERRQ(ierr);

    for (const auto& subdomain_dof_list : subdomain_dofs)
    {
        IS overlap_subdomain = nullptr;
        ierr = ISCreateGeneral(PETSC_COMM_SELF,
                               static_cast<PetscInt>(subdomain_dof_list.size()),
                               subdomain_dof_list.data(),
                               PETSC_COPY_VALUES,
                               &overlap_subdomain);
        IBTK_CHKERRQ(ierr);
        ierr = MatMult(A, y, r);
        IBTK_CHKERRQ(ierr);
        ierr = VecAYPX(r, -1.0, b);
        IBTK_CHKERRQ(ierr);

        Vec r_sub_view = nullptr;
        ierr = VecGetSubVector(r, overlap_subdomain, &r_sub_view);
        IBTK_CHKERRQ(ierr);
        Vec r_sub = nullptr;
        ierr = VecDuplicate(r_sub_view, &r_sub);
        IBTK_CHKERRQ(ierr);
        ierr = VecCopy(r_sub_view, r_sub);
        IBTK_CHKERRQ(ierr);
        ierr = VecRestoreSubVector(r, overlap_subdomain, &r_sub_view);
        IBTK_CHKERRQ(ierr);

        Mat A_sub = nullptr;
        ierr = MatCreateSubMatrix(A, overlap_subdomain, overlap_subdomain, MAT_INITIAL_MATRIX, &A_sub);
        IBTK_CHKERRQ(ierr);

        Vec delta_sub = nullptr;
        ierr = VecDuplicate(r_sub, &delta_sub);
        IBTK_CHKERRQ(ierr);
        KSP sub_ksp = nullptr;
        ierr = KSPCreate(PETSC_COMM_SELF, &sub_ksp);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSetReusePreconditioner(sub_ksp, PETSC_TRUE);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSetType(sub_ksp, KSPPREONLY);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSetOperators(sub_ksp, A_sub, A_sub);
        IBTK_CHKERRQ(ierr);
        PC sub_pc = nullptr;
        ierr = KSPGetPC(sub_ksp, &sub_pc);
        IBTK_CHKERRQ(ierr);
        ierr = PCSetType(sub_pc, PCLU);
        IBTK_CHKERRQ(ierr);
        ierr = PCFactorReorderForNonzeroDiagonal(sub_pc, std::numeric_limits<double>::epsilon());
        IBTK_CHKERRQ(ierr);
        ierr = KSPSolve(sub_ksp, r_sub, delta_sub);
        IBTK_CHKERRQ(ierr);
        if (alpha != 1.0)
        {
            ierr = VecScale(delta_sub, alpha);
            IBTK_CHKERRQ(ierr);
        }

        Vec y_sub = nullptr;
        ierr = VecGetSubVector(y, overlap_subdomain, &y_sub);
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPY(y_sub, 1.0, delta_sub);
        IBTK_CHKERRQ(ierr);
        ierr = VecRestoreSubVector(y, overlap_subdomain, &y_sub);
        IBTK_CHKERRQ(ierr);

        ierr = KSPDestroy(&sub_ksp);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&delta_sub);
        IBTK_CHKERRQ(ierr);
        ierr = MatDestroy(&A_sub);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&r_sub);
        IBTK_CHKERRQ(ierr);
        ierr = ISDestroy(&overlap_subdomain);
        IBTK_CHKERRQ(ierr);
    }

    Vec p_sub = nullptr;
    ierr = VecGetSubVector(y, pressure_is, &p_sub);
    IBTK_CHKERRQ(ierr);
    PetscScalar p_sum = 0.0;
    ierr = VecSum(p_sub, &p_sum);
    IBTK_CHKERRQ(ierr);
    PetscInt n_p = 0;
    ierr = VecGetSize(p_sub, &n_p);
    IBTK_CHKERRQ(ierr);
    if (n_p > 0)
    {
        const PetscScalar p_mean = p_sum / static_cast<PetscScalar>(n_p);
        ierr = VecShift(p_sub, -p_mean);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecRestoreSubVector(y, pressure_is, &p_sub);
    IBTK_CHKERRQ(ierr);

    ierr = VecDestroy(&r);
    IBTK_CHKERRQ(ierr);
}

double
vec_norm_inf(Vec x)
{
    PetscReal norm = 0.0;
    int ierr = VecNorm(x, NORM_INFINITY, &norm);
    IBTK_CHKERRQ(ierr);
    return static_cast<double>(norm);
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app_initializer->getInputDatabase();
    Pointer<Database> test_db = input_db->keyExists("test") ? input_db->getDatabase("test") : input_db;

    const std::string closure_policy = test_db->getStringWithDefault("coupling_aware_asm_closure_policy", "RELAXED");
    const int seed_axis = test_db->getIntegerWithDefault("coupling_aware_asm_seed_axis", 0);
    const int seed_stride = test_db->getIntegerWithDefault("coupling_aware_asm_seed_stride", 1);
    const double tol = test_db->getDoubleWithDefault("parity_tol", 1.0e-11);
    const std::string shell_pc_type = test_db->getStringWithDefault("shell_pc_type", "multiplicative");
    const bool test_all_eigen_reference_solver_types =
        test_db->getBoolWithDefault("test_all_eigen_reference_solver_types", false);
    const bool verify_reference_parity = test_db->getBoolWithDefault("verify_reference_parity", true);

    std::vector<std::string> solver_types;
    const bool is_eigen_reference_case = shell_pc_type.find("eigen-reference") != std::string::npos;
    if (is_eigen_reference_case && test_all_eigen_reference_solver_types)
    {
        solver_types = { "llt",
                         "ldlt",
                         "partial-piv-lu",
                         "full-piv-lu",
                         "householder-qr",
                         "col-piv-householder-qr",
                         "complete-orthogonal-decomposition",
                         "full-piv-householder-qr",
                         "jacobi-svd",
                         "bdc-svd" };
    }
    else if (test_db->keyExists("eigen_subdomain_solver_type"))
    {
        solver_types = { test_db->getString("eigen_subdomain_solver_type") };
    }
    else
    {
        solver_types = { "" };
    }

    int test_failures = 0;

    const auto hierarchy_tuple = setup_hierarchy<NDIM>(app_initializer);
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = std::get<0>(hierarchy_tuple);
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(0);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("stokes_petsc_level_solver_shell_multiplicative_semantics_ctx");
    Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("semantics_u");
    Pointer<CellVariable<NDIM, double>> p_var = new CellVariable<NDIM, double>("semantics_p");
    Pointer<SideVariable<NDIM, double>> f_u_var = new SideVariable<NDIM, double>("semantics_f_u");
    Pointer<CellVariable<NDIM, double>> f_p_var = new CellVariable<NDIM, double>("semantics_f_p");
    Pointer<SideVariable<NDIM, int>> u_dof_index_var = new SideVariable<NDIM, int>("semantics_u_dof");
    Pointer<CellVariable<NDIM, int>> p_dof_index_var = new CellVariable<NDIM, int>("semantics_p_dof");

    const IntVector<NDIM> one_ghost(1);
    const int u_idx = var_db->registerVariableAndContext(u_var, ctx, one_ghost);
    const int p_idx = var_db->registerVariableAndContext(p_var, ctx, one_ghost);
    const int f_u_idx = var_db->registerVariableAndContext(f_u_var, ctx, one_ghost);
    const int f_p_idx = var_db->registerVariableAndContext(f_p_var, ctx, one_ghost);
    const int u_dof_index_idx = var_db->registerVariableAndContext(u_dof_index_var, ctx, one_ghost);
    const int p_dof_index_idx = var_db->registerVariableAndContext(p_dof_index_var, ctx, one_ghost);

    for (const int data_idx : { u_idx, p_idx, f_u_idx, f_p_idx, u_dof_index_idx, p_dof_index_idx })
    {
        level->allocatePatchData(data_idx);
    }

    SAMRAIVectorReal<NDIM, double> x_vec("x", patch_hierarchy, 0, 0);
    SAMRAIVectorReal<NDIM, double> b_vec("b", patch_hierarchy, 0, 0);
    x_vec.addComponent(u_var, u_idx);
    x_vec.addComponent(p_var, p_idx);
    b_vec.addComponent(f_u_var, f_u_idx);
    b_vec.addComponent(f_p_var, f_p_idx);

    zero_solution_fields(level, u_idx, p_idx);
    fill_nontrivial_rhs(level, f_u_idx, f_p_idx);

    std::vector<int> num_dofs_per_proc;
    IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
        num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, level);

    std::vector<std::set<int>> field_is;
    std::vector<std::string> field_names;
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelFields(
        field_is, field_names, num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, level);
    const auto pressure_name_it = std::find(field_names.begin(), field_names.end(), "pressure");
    if (pressure_name_it == field_names.end())
    {
        TBOX_ERROR("stokes_petsc_level_solver_shell_multiplicative_semantics: pressure field not found.\n");
    }
    const std::size_t pressure_idx = static_cast<std::size_t>(std::distance(field_names.begin(), pressure_name_it));
    std::vector<PetscInt> pressure_dofs(field_is[pressure_idx].begin(), field_is[pressure_idx].end());
    IS pressure_is = nullptr;
    int ierr = ISCreateGeneral(PETSC_COMM_WORLD,
                               static_cast<PetscInt>(pressure_dofs.size()),
                               pressure_dofs.empty() ? nullptr : pressure_dofs.data(),
                               PETSC_COPY_VALUES,
                               &pressure_is);
    IBTK_CHKERRQ(ierr);

    const int rank = IBTK_MPI::getRank();
    Vec b_petsc = nullptr;
    Vec x_petsc = nullptr;
    Vec x_expected = nullptr;
    Vec x_diff = nullptr;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, num_dofs_per_proc[rank], PETSC_DETERMINE, &b_petsc);
    IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, num_dofs_per_proc[rank], PETSC_DETERMINE, &x_petsc);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(b_petsc, &x_expected);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(b_petsc, &x_diff);
    IBTK_CHKERRQ(ierr);
    IBAMR::StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        b_petsc, f_u_idx, u_dof_index_idx, f_p_idx, p_dof_index_idx, level);

    for (const std::string& solver_type : solver_types)
    {
        zero_solution_fields(level, u_idx, p_idx);

        Pointer<MemoryDatabase> solver_db = new MemoryDatabase("solver_db");
        solver_db->putString("ksp_type", "preonly");
        solver_db->putString("pc_type", "shell");
        solver_db->putString("shell_pc_type", shell_pc_type);
        solver_db->putInteger("max_iterations", 1);
        solver_db->putBool("initial_guess_nonzero", false);
        solver_db->putString("asm_subdomain_construction_mode", "COUPLING_AWARE");
        solver_db->putString("coupling_aware_asm_closure_policy", closure_policy);
        solver_db->putInteger("coupling_aware_asm_seed_axis", seed_axis);
        solver_db->putInteger("coupling_aware_asm_seed_stride", seed_stride);
        if (!solver_type.empty()) solver_db->putString("eigen_subdomain_solver_type", solver_type);
        if (test_db->keyExists("eigen_subdomain_solver_threshold"))
        {
            solver_db->putDouble("eigen_subdomain_solver_threshold", test_db->getDouble("eigen_subdomain_solver_threshold"));
        }

        Pointer<IBAMR::StaggeredStokesPETScLevelSolver> solver = new IBAMR::StaggeredStokesPETScLevelSolver(
            "solver_shell_multiplicative_semantics", solver_db, "stokes_shell_sem_");
        PoissonSpecifications problem_coefs("stokes_shell_sem_poisson");
        problem_coefs.setCConstant(1.0);
        problem_coefs.setDConstant(-1.0);
        solver->setVelocityPoissonSpecifications(problem_coefs);
        solver->initializeSolverState(x_vec, b_vec);

        std::vector<std::vector<int>>* overlap_is = nullptr;
        std::vector<std::vector<int>>* nonoverlap_is = nullptr;
        solver->getASMSubdomains(&nonoverlap_is, &overlap_is);

        const KSP& petsc_ksp = solver->getPETScKSP();
        Mat A_mat = nullptr;
        Mat pc_mat = nullptr;
        ierr = KSPGetOperators(petsc_ksp, &A_mat, &pc_mat);
        IBTK_CHKERRQ(ierr);

        apply_matlab_cav_sweep(x_expected, b_petsc, A_mat, *overlap_is, 1.0, pressure_is);
        const double expected_inf_norm = vec_norm_inf(x_expected);
        if (!(expected_inf_norm > 0.0)) ++test_failures;

        const bool converged = solver->solveSystem(x_vec, b_vec);
        if (!converged) ++test_failures;
        IBAMR::StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            x_petsc, u_idx, u_dof_index_idx, p_idx, p_dof_index_idx, level);
        const double actual_inf_norm = vec_norm_inf(x_petsc);
        if (!(actual_inf_norm > 0.0)) ++test_failures;

        ierr = VecCopy(x_petsc, x_diff);
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPY(x_diff, -1.0, x_expected);
        IBTK_CHKERRQ(ierr);
        const double error_inf_norm = vec_norm_inf(x_diff);

        if (verify_reference_parity)
        {
            if (!(error_inf_norm <= tol)) ++test_failures;
        }
        else
        {
            if (!std::isfinite(error_inf_norm)) ++test_failures;
        }

        solver->deallocateSolverState();

        const std::string solver_label = solver_type.empty() ? std::string("default") : solver_type;
        pout << "solver_type = " << solver_label << "\n";
        if (verify_reference_parity)
        {
            pout << "expected_inf_norm = " << expected_inf_norm << "\n";
            pout << "actual_inf_norm = " << actual_inf_norm << "\n";
            pout << "error_inf_norm = " << error_inf_norm << "\n";
        }
    }

    ierr = VecDestroy(&x_diff);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&x_expected);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&x_petsc);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&b_petsc);
    IBTK_CHKERRQ(ierr);
    ierr = ISDestroy(&pressure_is);
    IBTK_CHKERRQ(ierr);

    for (const int data_idx : { u_idx, p_idx, f_u_idx, f_p_idx, u_dof_index_idx, p_dof_index_idx })
    {
        level->deallocatePatchData(data_idx);
    }

    pout << "shell_pc_type = " << shell_pc_type << "\n";
    pout << "coupling_aware_asm_closure_policy = " << closure_policy << "\n";
    pout << "parity_tol = " << tol << "\n";
    pout << "test_failures = " << test_failures << "\n";
    return test_failures > 0 ? 1 : 0;
}
