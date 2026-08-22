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
#include <array>
#include <limits>
#include <numeric>
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
apply_reference_cav_sweep(Vec y,
                          Vec b,
                          Mat A,
                          const std::vector<std::vector<int>>& overlap_subdomain_dofs,
                          const double alpha,
                          IS pressure_is)
{
    int ierr = VecZeroEntries(y);
    IBTK_CHKERRQ(ierr);
    Vec r = nullptr;
    ierr = VecDuplicate(y, &r);
    IBTK_CHKERRQ(ierr);

    for (const auto& overlap_dof_list : overlap_subdomain_dofs)
    {
        IS overlap_subdomain = nullptr;
        ierr = ISCreateGeneral(PETSC_COMM_SELF,
                               static_cast<PetscInt>(overlap_dof_list.size()),
                               overlap_dof_list.data(),
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
    const std::string patch_seed_type = test_db->keyExists("coupling_aware_asm_patch_seed_type") ?
                                            test_db->getString("coupling_aware_asm_patch_seed_type") :
                                            "VELOCITY_COMPONENT";
    const std::string shell_pc_type = test_db->getStringWithDefault("shell_pc_type", "multiplicative");
    const int seed_axis = test_db->getIntegerWithDefault("coupling_aware_asm_seed_axis", 0);
    const int seed_stride = test_db->getIntegerWithDefault("coupling_aware_asm_seed_stride", 1);
    const std::string default_seed_traversal_order = NDIM == 2 ? "I_J" : "I_J_K";
    const std::string seed_traversal_order = test_db->keyExists("coupling_aware_asm_seed_traversal_order") ?
                                                 test_db->getString("coupling_aware_asm_seed_traversal_order") :
                                                 default_seed_traversal_order;
    const double alpha = test_db->getDoubleWithDefault("shell_pc_relaxation_factor", 1.0);
    const double tol = test_db->getDoubleWithDefault("parity_tol", 1.0e-11);
    const bool require_parity = test_db->getBoolWithDefault("require_parity", true);
    const bool verify_reinitialize =
        test_db->keyExists("verify_reinitialize") && test_db->getBool("verify_reinitialize");

    int test_failures = 0;

    const auto hierarchy_tuple = setup_hierarchy<NDIM>(app_initializer);
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = std::get<0>(hierarchy_tuple);
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(0);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("stokes_petsc_level_solver_cav_shell_semantics_ctx");
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

    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
        u_data->fillAll(0.0);
        p_data->fillAll(0.0);
    }
    fill_nontrivial_rhs(level, f_u_idx, f_p_idx);

    std::vector<int> num_dofs_per_proc;
    IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
        num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, level);

    std::vector<std::set<int>> field_is;
    std::vector<std::string> field_names;
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelFields(
        field_is, field_names, num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, level);
    const auto pressure_name_it = std::find(field_names.begin(), field_names.end(), "pressure");
    const auto velocity_name_it = std::find(field_names.begin(), field_names.end(), "velocity");
    if (pressure_name_it == field_names.end() || velocity_name_it == field_names.end())
    {
        TBOX_ERROR("stokes_petsc_level_solver_cav_shell_semantics: velocity or pressure field not found.\n");
    }
    const std::size_t pressure_idx = static_cast<std::size_t>(std::distance(field_names.begin(), pressure_name_it));
    const std::size_t velocity_idx = static_cast<std::size_t>(std::distance(field_names.begin(), velocity_name_it));
    std::vector<PetscInt> pressure_dofs(field_is[pressure_idx].begin(), field_is[pressure_idx].end());
    std::vector<PetscInt> velocity_dofs(field_is[velocity_idx].begin(), field_is[velocity_idx].end());
    struct PressureSeedRecord
    {
        std::array<int, NDIM> logical_index{};
        int dof = -1;
    };
    std::vector<PressureSeedRecord> pressure_seed_records;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<CellData<NDIM, int>> p_dof_data = patch->getPatchData(p_dof_index_idx);
        for (Box<NDIM>::Iterator b(patch->getBox()); b; b++)
        {
            PressureSeedRecord record;
            for (unsigned int d = 0; d < NDIM; ++d) record.logical_index[d] = b()(static_cast<int>(d));
            record.dof = (*p_dof_data)(b());
            pressure_seed_records.push_back(record);
        }
    }
    std::sort(pressure_seed_records.begin(),
              pressure_seed_records.end(),
              [](const PressureSeedRecord& lhs, const PressureSeedRecord& rhs)
              {
                  for (unsigned int d = 0; d < NDIM; ++d)
                  {
                      if (lhs.logical_index[d] < rhs.logical_index[d]) return true;
                      if (lhs.logical_index[d] > rhs.logical_index[d]) return false;
                  }
                  return lhs.dof < rhs.dof;
              });
    std::vector<int> expected_pressure_seed_dofs;
    std::set<int> seen_pressure_dofs;
    for (const PressureSeedRecord& record : pressure_seed_records)
    {
        if (seen_pressure_dofs.insert(record.dof).second) expected_pressure_seed_dofs.push_back(record.dof);
    }
    IS pressure_is = nullptr;
    int ierr = ISCreateGeneral(PETSC_COMM_WORLD,
                               static_cast<PetscInt>(pressure_dofs.size()),
                               pressure_dofs.empty() ? nullptr : pressure_dofs.data(),
                               PETSC_COPY_VALUES,
                               &pressure_is);
    IBTK_CHKERRQ(ierr);

    Pointer<MemoryDatabase> solver_db = new MemoryDatabase("solver_db");
    solver_db->putString("ksp_type", "preonly");
    solver_db->putString("pc_type", "shell");
    solver_db->putString("shell_pc_type", shell_pc_type);
    solver_db->putDouble("shell_pc_relaxation_factor", alpha);
    solver_db->putInteger("max_iterations", 1);
    solver_db->putBool("initial_guess_nonzero", false);
    solver_db->putString("asm_subdomain_construction_mode", "COUPLING_AWARE");
    solver_db->putString("coupling_aware_asm_closure_policy", closure_policy);
    solver_db->putString("coupling_aware_asm_patch_seed_type", patch_seed_type);
    solver_db->putInteger("coupling_aware_asm_seed_axis", seed_axis);
    solver_db->putInteger("coupling_aware_asm_seed_stride", seed_stride);
    solver_db->putString("coupling_aware_asm_seed_traversal_order", seed_traversal_order);
    if (test_db->keyExists("blas_lapack_subdomain_solver_type"))
    {
        solver_db->putString("blas_lapack_subdomain_solver_type",
                             test_db->getString("blas_lapack_subdomain_solver_type"));
    }

    Pointer<IBAMR::StaggeredStokesPETScLevelSolver> solver =
        new IBAMR::StaggeredStokesPETScLevelSolver("solver_cav_shell_semantics", solver_db, "stokes_shell_sem_");
    PoissonSpecifications problem_coefs("stokes_shell_sem_poisson");
    problem_coefs.setCConstant(1.0);
    problem_coefs.setDConstant(-1.0);
    solver->setVelocityPoissonSpecifications(problem_coefs);
    Mat construction_mat = nullptr;
    if (patch_seed_type == "PRESSURE_CELL")
    {
        const int n_global_dofs = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.end(), 0);
        ierr = MatCreateAIJ(PETSC_COMM_WORLD,
                            PETSC_DECIDE,
                            PETSC_DECIDE,
                            n_global_dofs,
                            n_global_dofs,
                            2,
                            nullptr,
                            1,
                            nullptr,
                            &construction_mat);
        IBTK_CHKERRQ(ierr);
        if (velocity_dofs.size() >= 2)
        {
            const PetscInt coupled_velocity_dofs[2] = { velocity_dofs.front(), velocity_dofs.back() };
            const PetscScalar coupling_values[4] = { 0.0, 0.25, 0.25, 0.0 };
            ierr = MatSetValues(
                construction_mat, 2, coupled_velocity_dofs, 2, coupled_velocity_dofs, coupling_values, INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
        ierr = MatAssemblyBegin(construction_mat, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyEnd(construction_mat, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);
        solver->setCouplingAwareASMConstructionMat(construction_mat);
    }
    solver->initializeSolverState(x_vec, b_vec);

    const auto& overlap_is = solver->getASMSubdomains();
    const auto& nonoverlap_is = solver->getASMNonoverlapSubdomains();
    const std::size_t pressure_seed_patch_count = overlap_is.size();
    bool pressure_seed_order_valid = patch_seed_type != "PRESSURE_CELL";
    bool pressure_seed_partition_absent = patch_seed_type != "PRESSURE_CELL";
    if (patch_seed_type == "PRESSURE_CELL")
    {
        const auto& pressure_seeds = solver->getCouplingAwareASMPressureSeedDOFs();
        pressure_seed_order_valid = pressure_seeds.size() == expected_pressure_seed_dofs.size() &&
                                    pressure_seeds.size() == overlap_is.size() &&
                                    pressure_seeds == expected_pressure_seed_dofs;
        for (std::size_t k = 0; k < std::min(pressure_seeds.size(), overlap_is.size()); ++k)
        {
            pressure_seed_order_valid =
                pressure_seed_order_valid &&
                std::binary_search(overlap_is[k].begin(), overlap_is[k].end(), pressure_seeds[k]);
        }
        pressure_seed_partition_absent =
            nonoverlap_is.size() == overlap_is.size() &&
            std::all_of(nonoverlap_is.begin(), nonoverlap_is.end(), [](const auto& dofs) { return dofs.empty(); });
        if (!pressure_seed_order_valid || !pressure_seed_partition_absent) ++test_failures;
    }

    const KSP& petsc_ksp = solver->getPETScKSP();
    Mat A_mat = nullptr;
    Mat pc_mat = nullptr;
    ierr = KSPGetOperators(petsc_ksp, &A_mat, &pc_mat);
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

    apply_reference_cav_sweep(x_expected, b_petsc, A_mat, overlap_is, alpha, pressure_is);
    const double expected_inf_norm = vec_norm_inf(x_expected);
    if (require_parity && !(expected_inf_norm > 0.0)) ++test_failures;

    const bool converged = solver->solveSystem(x_vec, b_vec);
    if (require_parity && !converged) ++test_failures;
    IBAMR::StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        x_petsc, u_idx, u_dof_index_idx, p_idx, p_dof_index_idx, level);
    const double actual_inf_norm = vec_norm_inf(x_petsc);
    if (require_parity && !(actual_inf_norm > 0.0)) ++test_failures;

    ierr = VecCopy(x_petsc, x_diff);
    IBTK_CHKERRQ(ierr);
    ierr = VecAXPY(x_diff, -1.0, x_expected);
    IBTK_CHKERRQ(ierr);
    const double error_inf_norm = vec_norm_inf(x_diff);
    if (require_parity && !(error_inf_norm <= tol)) ++test_failures;

    solver->deallocateSolverState();
    if (patch_seed_type == "PRESSURE_CELL" && !solver->getCouplingAwareASMPressureSeedDOFs().empty()) ++test_failures;

    bool reinitialize_pressure_seed_order_valid = true;
    double reinitialize_error_inf_norm = std::numeric_limits<double>::quiet_NaN();
    if (verify_reinitialize)
    {
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
            u_data->fillAll(0.0);
            p_data->fillAll(0.0);
        }
        // The construction matrix is borrowed only for one initialized state,
        // so its live owner must resupply it before reinitialization.
        if (construction_mat) solver->setCouplingAwareASMConstructionMat(construction_mat);
        solver->initializeSolverState(x_vec, b_vec);
        if (patch_seed_type == "PRESSURE_CELL")
        {
            reinitialize_pressure_seed_order_valid =
                solver->getCouplingAwareASMPressureSeedDOFs() == expected_pressure_seed_dofs;
            if (!reinitialize_pressure_seed_order_valid) ++test_failures;
        }
        const bool reinitialize_converged = solver->solveSystem(x_vec, b_vec);
        if (!reinitialize_converged) ++test_failures;
        IBAMR::StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            x_petsc, u_idx, u_dof_index_idx, p_idx, p_dof_index_idx, level);
        ierr = VecCopy(x_petsc, x_diff);
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPY(x_diff, -1.0, x_expected);
        IBTK_CHKERRQ(ierr);
        reinitialize_error_inf_norm = vec_norm_inf(x_diff);
        if (!(reinitialize_error_inf_norm <= tol)) ++test_failures;
        solver->deallocateSolverState();
        if (patch_seed_type == "PRESSURE_CELL" && !solver->getCouplingAwareASMPressureSeedDOFs().empty())
            ++test_failures;
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

    ierr = MatDestroy(&construction_mat);
    IBTK_CHKERRQ(ierr);

    for (const int data_idx : { u_idx, p_idx, f_u_idx, f_p_idx, u_dof_index_idx, p_dof_index_idx })
    {
        level->deallocatePatchData(data_idx);
    }

    plog << "Input database:\n";
    input_db->printClassData(plog);
    pout << "coupling_aware_asm_closure_policy = " << closure_policy << "\n";
    if (patch_seed_type == "PRESSURE_CELL")
    {
        pout << "coupling_aware_asm_patch_seed_type = " << patch_seed_type << "\n";
        pout << "pressure_seed_patch_count = " << pressure_seed_patch_count << "\n";
        pout << "pressure_seed_order_valid = " << (pressure_seed_order_valid ? "true" : "false") << "\n";
        pout << "pressure_seed_partition_absent = " << (pressure_seed_partition_absent ? "true" : "false") << "\n";
    }
    pout << "shell_pc_type = " << shell_pc_type << "\n";
    pout << "shell_pc_relaxation_factor = " << alpha << "\n";
    pout << "require_parity = " << require_parity << "\n";
    pout << "expected_inf_norm = " << expected_inf_norm << "\n";
    pout << "actual_inf_norm = " << actual_inf_norm << "\n";
    pout << "error_inf_norm = " << error_inf_norm << "\n";
    if (verify_reinitialize)
    {
        pout << "reinitialize_pressure_seed_order_valid = "
             << (reinitialize_pressure_seed_order_valid ? "true" : "false") << "\n";
        pout << "reinitialize_error_inf_norm = "
             << (reinitialize_error_inf_norm <= tol ? 0.0 : reinitialize_error_inf_norm) << "\n";
    }
    pout << "parity_tol = " << tol << "\n";
    pout << "test_failures = " << test_failures << "\n";
    return test_failures > 0 ? 1 : 0;
}
