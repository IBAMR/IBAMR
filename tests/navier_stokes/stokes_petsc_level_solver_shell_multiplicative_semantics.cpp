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
#include <ibtk/private/PETScLevelSolverShellBackend.h>

#include <tbox/MemoryDatabase.h>

#include <petscksp.h>
#include <petscsys.h>

#include <CellData.h>
#include <CellVariable.h>
#include <LocationIndexRobinBcCoefs.h>
#include <PoissonSpecifications.h>
#include <RobinBcCoefStrategy.h>
#include <SAMRAIVectorReal.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "../tests.h"
#include "cav_raw_operator_comparator.h"

#include <ibtk/app_namespaces.h>

namespace
{
struct RowAccessMatrixContext
{
    std::vector<std::vector<PetscInt>> columns;
    std::vector<std::vector<PetscScalar>> values;
    std::size_t get_row_calls = 0;
};

PetscErrorCode
get_row_from_test_matrix(Mat mat,
                         const PetscInt row,
                         PetscInt* n_columns,
                         const PetscInt** columns,
                         const PetscScalar** values)
{
    RowAccessMatrixContext* context = nullptr;
    PetscErrorCode ierr = MatShellGetContext(mat, &context);
    CHKERRQ(ierr);
    ++context->get_row_calls;
    *n_columns = static_cast<PetscInt>(context->columns[static_cast<std::size_t>(row)].size());
    *columns = context->columns[static_cast<std::size_t>(row)].data();
    if (values) *values = context->values[static_cast<std::size_t>(row)].data();
    return PETSC_SUCCESS;
}

PetscErrorCode
restore_row_from_test_matrix(Mat /*mat*/,
                             PetscInt /*row*/,
                             PetscInt* /*n_columns*/,
                             const PetscInt** /*columns*/,
                             const PetscScalar** /*values*/)
{
    return PETSC_SUCCESS;
}

class RecordingShellBackend : public IBTK::PETScLevelSolverShellBackend
{
public:
    const std::string& getName() const override
    {
        static const std::string name = "Recording";
        return name;
    }

    void initializeSolverState(const IBTK::PETScLevelSolverShellBackendState& solver_state) override
    {
        initializeCorrectionCompositionState(solver_state);
        finalizeCorrectionCompositionState();
    }

    void deallocateSolverState() override
    {
        deallocateCorrectionCompositionState();
    }

    const std::vector<std::string>& getEvents() const
    {
        return d_events;
    }

    const std::vector<PetscScalar>& getRhsSamples() const
    {
        return d_rhs_samples;
    }

    void recordEvent(std::string event)
    {
        d_events.push_back(std::move(event));
    }

protected:
    std::size_t getNumberOfSubdomains() const override
    {
        return 3;
    }

    void initializeSubdomainSweep(Vec /*x*/, Vec /*y*/) override
    {
        d_events.push_back("initialize");
    }

    void beginSubdomainRhs(const std::size_t subdomain_num, Vec x, Vec /*y*/) override
    {
        d_events.push_back("begin" + std::to_string(subdomain_num));
        const PetscInt dof = (*d_solver_state.subdomain_dofs)[subdomain_num].front();
        PetscScalar value = 0.0;
        const int ierr = VecGetValues(getSubdomainResidualSource(x), 1, &dof, &value);
        IBTK_CHKERRQ(ierr);
        d_rhs_samples.push_back(value);
    }

    void endSubdomainRhs(const std::size_t subdomain_num, Vec /*x*/, Vec /*y*/) override
    {
        d_events.push_back("end" + std::to_string(subdomain_num));
    }

    void solveSubdomain(const std::size_t subdomain_num) override
    {
        d_events.push_back("solve" + std::to_string(subdomain_num));
    }

    void observeSubdomainSolve(const std::size_t subdomain_num, Vec /*x*/, Vec /*y*/) override
    {
        d_events.push_back("observe" + std::to_string(subdomain_num));
    }

    void accumulateSubdomainCorrection(const std::size_t subdomain_num, Vec y) override
    {
        d_events.push_back("accumulate" + std::to_string(subdomain_num));
        const auto& dofs = (*d_solver_state.subdomain_dofs)[subdomain_num];
        std::vector<PetscScalar> values(dofs.size(), 1.0);
        int ierr = VecSetValues(y, static_cast<PetscInt>(dofs.size()), dofs.data(), values.data(), ADD_VALUES);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyBegin(y);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyEnd(y);
        IBTK_CHKERRQ(ierr);
    }

    const std::vector<int>& getSubdomainCorrectionDofs(const std::size_t subdomain_num) const override
    {
        return (*d_solver_state.subdomain_dofs)[subdomain_num];
    }

    void copySubdomainCorrection(const std::size_t subdomain_num, PetscScalar* correction_values) override
    {
        d_events.push_back("update" + std::to_string(subdomain_num));
        std::fill_n(correction_values, (*d_solver_state.subdomain_dofs)[subdomain_num].size(), 1.0);
    }

    void finalizeSubdomainSweep(Vec /*x*/, Vec /*y*/) override
    {
        d_events.push_back("finalize");
    }

private:
    std::vector<std::string> d_events;
    std::vector<PetscScalar> d_rhs_samples;
};

bool
verify_shared_composition_driver(const bool use_multiplicative, std::size_t& event_count)
{
    RecordingShellBackend backend;
    Vec x = nullptr, y = nullptr;
    int ierr = VecCreateSeq(PETSC_COMM_SELF, 1, &x);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(x, &y);
    IBTK_CHKERRQ(ierr);
    ierr = VecSet(x, 1.0);
    IBTK_CHKERRQ(ierr);
    Mat mat = nullptr;
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, 1, 1, 1, nullptr, &mat);
    IBTK_CHKERRQ(ierr);
    const PetscInt row = 0;
    const PetscScalar value = 1.0;
    ierr = MatSetValue(mat, row, row, value, INSERT_VALUES);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);

    const std::vector<std::vector<int>> subdomains = { { 0 }, { 0 }, { 0 } };
    IBTK::PETScLevelSolverShellBackendState solver_state;
    solver_state.use_multiplicative = use_multiplicative;
    solver_state.petsc_mat = mat;
    solver_state.petsc_x = x;
    solver_state.petsc_b = x;
    solver_state.subdomain_dofs = &subdomains;
    solver_state.nonoverlap_subdomain_dofs = &subdomains;
    std::function<void(int, Mat, Vec, Vec, Vec)> observer = [](int, Mat, Vec, Vec, Vec) {};
    std::function<bool(int)> predicate = [](const int subdomain_num) { return subdomain_num != 1; };
    solver_state.subdomain_solve_observer = &observer;
    solver_state.subdomain_solve_observer_predicate = &predicate;
    solver_state.postprocess_result = [&backend](Vec) { backend.recordEvent("postprocess"); };
    backend.initializeSolverState(solver_state);

    backend.apply(x, y);
    backend.deallocateSolverState();
    ierr = MatDestroy(&mat);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&y);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    IBTK_CHKERRQ(ierr);

    const std::vector<std::string> additive_events = {
        "initialize", "begin0",      "begin1", "begin2", "end0",     "solve0",      "observe0", "accumulate0", "end1",
        "solve1",     "accumulate1", "end2",   "solve2", "observe2", "accumulate2", "finalize", "postprocess"
    };
    const std::vector<std::string> multiplicative_events = { "initialize",  "begin0",      "end0",        "solve0",
                                                             "observe0",    "accumulate0", "update0",     "begin1",
                                                             "end1",        "solve1",      "accumulate1", "update1",
                                                             "begin2",      "end2",        "solve2",      "observe2",
                                                             "accumulate2", "finalize",    "postprocess" };
    event_count = backend.getEvents().size();
    return backend.getEvents() == (use_multiplicative ? multiplicative_events : additive_events);
}

bool
verify_affected_row_residual_update(std::size_t& row_visits,
                                    bool& stagewise_residual_valid,
                                    const bool invalidate_pattern = false)
{
    RowAccessMatrixContext context;
    context.columns = { { 0, 1 }, { 0, 1, 2 }, { 1, 2, 3 }, { 2, 3 } };
    context.values = { { 9.0, 9.0 }, { 9.0, 9.0, 9.0 }, { 9.0, 9.0, 9.0 }, { 9.0, 9.0 } };
    const std::vector<std::vector<PetscScalar>> live_values = {
        { 2.0, -1.0 }, { -1.0, 2.0, -1.0 }, { -1.0, 2.0, -1.0 }, { -1.0, 2.0 }
    };
    Mat mat = nullptr;
    int ierr = MatCreateShell(PETSC_COMM_SELF, 4, 4, 4, 4, &context, &mat);
    IBTK_CHKERRQ(ierr);
    ierr = MatShellSetOperation(mat, MATOP_GET_ROW, reinterpret_cast<void (*)(void)>(get_row_from_test_matrix));
    IBTK_CHKERRQ(ierr);
    ierr = MatShellSetOperation(mat, MATOP_RESTORE_ROW, reinterpret_cast<void (*)(void)>(restore_row_from_test_matrix));
    IBTK_CHKERRQ(ierr);

    Vec x = nullptr, y = nullptr;
    ierr = VecCreateSeq(PETSC_COMM_SELF, 4, &x);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(x, &y);
    IBTK_CHKERRQ(ierr);
    const PetscInt indices[] = { 0, 1, 2, 3 };
    const PetscScalar x_values[] = { 10.0, 20.0, 30.0, 40.0 };
    ierr = VecSetValues(x, 4, indices, x_values, INSERT_VALUES);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyBegin(x);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x);
    IBTK_CHKERRQ(ierr);

    const std::vector<std::vector<int>> subdomains = { { 0 }, { 1 }, { 2 } };
    IBTK::PETScLevelSolverShellBackendState solver_state;
    solver_state.object_name = "affected_row_test";
    solver_state.options_prefix = "affected_row_";
    solver_state.use_multiplicative = true;
    solver_state.petsc_mat = mat;
    solver_state.petsc_x = x;
    solver_state.petsc_b = x;
    solver_state.subdomain_dofs = &subdomains;
    solver_state.nonoverlap_subdomain_dofs = &subdomains;
    solver_state.postprocess_result = [](Vec) {};

    RecordingShellBackend backend;
    backend.initializeSolverState(solver_state);
    // The setup values are deliberately wrong: the isolated composition
    // check succeeds only when application reads values from the borrowed
    // matrix instead of a copied update matrix. Production operators still
    // require complete solver reinitialization after reassembly.
    context.values = live_values;
    if (invalidate_pattern) std::swap(context.columns[0][0], context.columns[0][1]);
    context.get_row_calls = 0;
    backend.apply(x, y);
    row_visits = context.get_row_calls;

    const std::vector<PetscScalar> expected_rhs = { 10.0, 21.0, 31.0 };
    stagewise_residual_valid = backend.getRhsSamples() == expected_rhs;
    const PetscScalar* y_values = nullptr;
    ierr = VecGetArrayRead(y, &y_values);
    IBTK_CHKERRQ(ierr);
    const bool correction_valid = y_values[0] == 1.0 && y_values[1] == 1.0 && y_values[2] == 1.0 && y_values[3] == 0.0;
    ierr = VecRestoreArrayRead(y, &y_values);
    IBTK_CHKERRQ(ierr);

    backend.deallocateSolverState();
    ierr = VecDestroy(&y);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&mat);
    IBTK_CHKERRQ(ierr);
    return row_visits == 5 && stagewise_residual_valid && correction_valid;
}

class TestableStaggeredStokesPETScLevelSolver : public IBAMR::StaggeredStokesPETScLevelSolver
{
public:
    using IBAMR::StaggeredStokesPETScLevelSolver::StaggeredStokesPETScLevelSolver;

    void copyAdjustedRHS(Vec petsc_b, SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& b)
    {
        setupKSPVecs(d_petsc_x, petsc_b, x, b);
    }

    bool usesMultiplicativeShellSmootherForTest() const
    {
        return usesMultiplicativeShellSmoother();
    }

    bool usesRestrictShellSmootherPartitionForTest() const
    {
        return usesRestrictShellSmootherPartition();
    }
};

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
normalize_pressure_gauge(Vec y, IS pressure_is)
{
    Vec p_sub = nullptr;
    int ierr = VecGetSubVector(y, pressure_is, &p_sub);
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
}

void
apply_reference_shell_action(Vec y,
                             Vec b,
                             Mat A,
                             const std::vector<std::vector<int>>& overlap_subdomain_dofs,
                             const std::vector<std::vector<int>>& nonoverlap_subdomain_dofs,
                             const bool use_multiplicative,
                             const bool use_restrict_partition,
                             const double alpha,
                             IS pressure_is)
{
    int ierr = VecZeroEntries(y);
    IBTK_CHKERRQ(ierr);
    Vec r = nullptr;
    ierr = VecDuplicate(y, &r);
    IBTK_CHKERRQ(ierr);

    for (std::size_t subdomain_num = 0; subdomain_num < overlap_subdomain_dofs.size(); ++subdomain_num)
    {
        const auto& overlap_dof_list = overlap_subdomain_dofs[subdomain_num];
        const auto& update_dof_list =
            use_restrict_partition ? nonoverlap_subdomain_dofs[subdomain_num] : overlap_subdomain_dofs[subdomain_num];
        IS overlap_subdomain = nullptr;
        ierr = ISCreateGeneral(PETSC_COMM_SELF,
                               static_cast<PetscInt>(overlap_dof_list.size()),
                               overlap_dof_list.data(),
                               PETSC_COPY_VALUES,
                               &overlap_subdomain);
        IBTK_CHKERRQ(ierr);
        if (use_multiplicative)
        {
            ierr = MatMult(A, y, r);
            IBTK_CHKERRQ(ierr);
            ierr = VecAYPX(r, -1.0, b);
            IBTK_CHKERRQ(ierr);
        }
        else
        {
            ierr = VecCopy(b, r);
            IBTK_CHKERRQ(ierr);
        }

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

        const PetscScalar* delta_sub_arr = nullptr;
        ierr = VecGetArrayRead(delta_sub, &delta_sub_arr);
        IBTK_CHKERRQ(ierr);

        PetscScalar* y_arr = nullptr;
        ierr = VecGetArray(y, &y_arr);
        IBTK_CHKERRQ(ierr);
        for (const int dof : update_dof_list)
        {
            const auto pos_it = std::find(overlap_dof_list.begin(), overlap_dof_list.end(), dof);
            TBOX_ASSERT(pos_it != overlap_dof_list.end());
            const std::size_t local_pos = static_cast<std::size_t>(std::distance(overlap_dof_list.begin(), pos_it));
            y_arr[dof] += delta_sub_arr[local_pos];
        }
        ierr = VecRestoreArray(y, &y_arr);
        IBTK_CHKERRQ(ierr);
        ierr = VecRestoreArrayRead(delta_sub, &delta_sub_arr);
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

    // Match StaggeredStokesPETScLevelSolver: only multiplicative shell
    // corrections receive an explicit zero-mean pressure gauge.
    if (use_multiplicative) normalize_pressure_gauge(y, pressure_is);

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
    Pointer<Logger::Appender> abort_append(new TestAppender());
    Logger::getInstance()->setAbortAppender(abort_append);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app_initializer->getInputDatabase();
    Pointer<Database> test_db = input_db->keyExists("test") ? input_db->getDatabase("test") : input_db;

    const std::string closure_policy = test_db->getStringWithDefault("coupling_aware_asm_closure_policy", "RELAXED");
    const int seed_axis = test_db->getIntegerWithDefault("coupling_aware_asm_seed_axis", 0);
    const int seed_stride = test_db->getIntegerWithDefault("coupling_aware_asm_seed_stride", 1);
    const double tol = test_db->getDoubleWithDefault("parity_tol", 1.0e-11);
    const bool has_shell_pc_type = test_db->keyExists("shell_pc_type");
    const std::string shell_pc_type = has_shell_pc_type ? test_db->getString("shell_pc_type") : "";
    const bool test_all_eigen_reference_solver_types =
        test_db->getBoolWithDefault("test_all_eigen_reference_solver_types", false);
    const bool test_all_blas_lapack_solver_types =
        test_db->getBoolWithDefault("test_all_blas_lapack_solver_types", false);
    const bool verify_reference_parity = test_db->getBoolWithDefault("verify_reference_parity", true);
    const bool report_actual_inf_norm = test_db->getBoolWithDefault("report_actual_inf_norm", false);
    const bool verify_rhs_boundary_adjustment = test_db->getBoolWithDefault("verify_rhs_boundary_adjustment", false);
    const bool verify_reinitialize = test_db->getBoolWithDefault("verify_reinitialize", false);
    const bool verify_local_solve_observer = test_db->getBoolWithDefault("verify_local_solve_observer", false);
    const bool verify_stagewise_original_residual =
        test_db->getBoolWithDefault("verify_stagewise_original_residual", false);
    const bool verify_shared_composition = test_db->getBoolWithDefault("verify_shared_composition_driver", false);
    const bool verify_affected_row_invalidation_failure =
        test_db->getBoolWithDefault("verify_affected_row_invalidation_failure", false);
    const bool verify_live_operator_state_view = test_db->getBoolWithDefault("verify_live_operator_state_view", false);
    const bool verify_supplied_operator_lifetime =
        test_db->getBoolWithDefault("verify_supplied_operator_lifetime", false);
    const bool verify_raw_export_comparator_contract =
        test_db->getBoolWithDefault("verify_raw_export_comparator_contract", false);
    const std::string configured_pc_type = test_db->getStringWithDefault("configured_pc_type", "shell");
    const std::string expected_shell_pc_name = test_db->getStringWithDefault("expected_shell_pc_name", "");
    const double velocity_poisson_c = test_db->getDoubleWithDefault("velocity_poisson_c", 1.0);
    const double velocity_poisson_d = test_db->getDoubleWithDefault("velocity_poisson_d", -1.0);

    std::vector<std::string> solver_types;
    const bool is_eigen_reference_case = shell_pc_type.find("eigen-reference") != std::string::npos;
    const bool is_blas_lapack_case = shell_pc_type.find("blas-lapack") != std::string::npos;
    if (is_eigen_reference_case && test_all_eigen_reference_solver_types)
    {
        solver_types = { "partial-piv-lu",
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
    else if (is_blas_lapack_case && test_all_blas_lapack_solver_types)
    {
        solver_types = { "", "svd", "lu", "symmetric-indefinite", "qr" };
    }
    else if (test_db->keyExists("blas_lapack_subdomain_solver_type"))
    {
        solver_types = { test_db->getString("blas_lapack_subdomain_solver_type") };
    }
    else
    {
        solver_types = { "" };
    }

    int test_failures = 0;
    if (verify_affected_row_invalidation_failure)
    {
        std::size_t ignored_row_visits = 0;
        bool ignored_stagewise_residual_valid = false;
        verify_affected_row_residual_update(ignored_row_visits, ignored_stagewise_residual_valid, true);
        ++test_failures;
    }
    bool shared_additive_order_valid = true, shared_multiplicative_order_valid = true;
    bool affected_row_update_valid = true, affected_row_stagewise_residual_valid = true;
    std::size_t shared_additive_event_count = 0, shared_multiplicative_event_count = 0;
    std::size_t affected_row_update_row_visits = 0;
    if (verify_shared_composition)
    {
        shared_additive_order_valid = verify_shared_composition_driver(false, shared_additive_event_count);
        shared_multiplicative_order_valid = verify_shared_composition_driver(true, shared_multiplicative_event_count);
        affected_row_update_valid =
            verify_affected_row_residual_update(affected_row_update_row_visits, affected_row_stagewise_residual_valid);
        if (!shared_additive_order_valid || !shared_multiplicative_order_valid || !affected_row_update_valid)
            ++test_failures;
    }

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
        solver_db->putString("pc_type", configured_pc_type);
        if (has_shell_pc_type) solver_db->putString("shell_pc_type", shell_pc_type);
        solver_db->putInteger("max_iterations", 1);
        solver_db->putBool("initial_guess_nonzero", false);
        solver_db->putString("asm_subdomain_construction_mode", "COUPLING_AWARE");
        solver_db->putString("coupling_aware_asm_closure_policy", closure_policy);
        solver_db->putInteger("coupling_aware_asm_seed_axis", seed_axis);
        solver_db->putInteger("coupling_aware_asm_seed_stride", seed_stride);
        if (!solver_type.empty())
        {
            if (is_eigen_reference_case)
            {
                solver_db->putString("eigen_subdomain_solver_type", solver_type);
            }
            else if (is_blas_lapack_case)
            {
                solver_db->putString("blas_lapack_subdomain_solver_type", solver_type);
            }
        }
        if (test_db->keyExists("eigen_subdomain_solver_threshold"))
        {
            solver_db->putDouble("eigen_subdomain_solver_threshold",
                                 test_db->getDouble("eigen_subdomain_solver_threshold"));
        }
        if (test_db->keyExists("blas_lapack_subdomain_solver_rcond"))
        {
            solver_db->putDouble("blas_lapack_subdomain_solver_rcond",
                                 test_db->getDouble("blas_lapack_subdomain_solver_rcond"));
        }

        Pointer<TestableStaggeredStokesPETScLevelSolver> solver = new TestableStaggeredStokesPETScLevelSolver(
            "solver_shell_multiplicative_semantics", solver_db, "stokes_shell_sem_");
        bool live_view_preinit_empty = true;
        bool live_view_operator_identity = true;
        bool live_view_ownership_valid = true;
        bool live_view_block_partition_valid = true;
        bool live_view_provenance_valid = true;
        bool live_view_nullspace_gauge_valid = true;
        bool live_view_post_deallocate_empty = true;
        bool supplied_operator_identity = !verify_supplied_operator_lifetime;
        bool supplied_operator_valid_after_deallocate = !verify_supplied_operator_lifetime;
        bool raw_export_round_trip = true;
        bool raw_export_global_dof_mapping = true;
        bool raw_export_pressure_row_minus_div = true;
        bool raw_export_borrowed_matrix_round_trip = true;
        bool raw_export_borrowed_vector_round_trip = true;
        bool raw_export_manifest_round_trip = true;
        bool raw_export_manifest_mutation_detected = true;
        bool raw_export_invalid_provenance_rejected = true;
        double raw_export_matrix_max_abs_error = 0.0;
        if (verify_live_operator_state_view || verify_raw_export_comparator_contract ||
            verify_supplied_operator_lifetime)
        {
            const auto view = solver->getLiveOperatorStateView();
            live_view_preinit_empty = !view.initialized && !view.operator_mat && !view.num_dofs_per_proc &&
                                      !view.locally_owned_velocity_dofs && !view.locally_owned_pressure_dofs;
            if (!live_view_preinit_empty) ++test_failures;
            solver->setComponentsHaveNullSpace(false, true);
        }
        const bool configured_use_multiplicative = solver->usesMultiplicativeShellSmootherForTest();
        if (is_eigen_reference_case && !configured_use_multiplicative)
        {
            pout << "unsupported shell_pc_type = " << shell_pc_type << "\n";
            pout << "reason = eigen-reference backend only supports multiplicative mode\n";
            pout << "suggested shell_pc_type = multiplicative-eigen-reference\n";
            return 1;
        }
        PoissonSpecifications problem_coefs("stokes_shell_sem_poisson");
        problem_coefs.setCConstant(velocity_poisson_c);
        problem_coefs.setDConstant(velocity_poisson_d);
        solver->setVelocityPoissonSpecifications(problem_coefs);
        std::vector<std::unique_ptr<LocationIndexRobinBcCoefs<NDIM>>> u_bc_coef_storage;
        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);
        if (verify_rhs_boundary_adjustment)
        {
            u_bc_coef_storage.reserve(NDIM);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                u_bc_coef_storage.push_back(std::make_unique<LocationIndexRobinBcCoefs<NDIM>>(
                    "rhs_adjust_u_bc_" + std::to_string(axis), Pointer<Database>(nullptr)));
                for (int location = 0; location < 2 * NDIM; ++location)
                {
                    u_bc_coef_storage.back()->setBoundaryValue(
                        location, 0.125 * static_cast<double>((axis + 1) * (location + 1)));
                }
                u_bc_coefs[axis] = u_bc_coef_storage.back().get();
            }
            Pointer<IBAMR::StaggeredStokesPhysicalBoundaryHelper> bc_helper =
                new IBAMR::StaggeredStokesPhysicalBoundaryHelper();
            bc_helper->cacheBcCoefData(u_bc_coefs, 0.0, patch_hierarchy);
            solver->setPhysicalBcCoefs(u_bc_coefs, nullptr);
            solver->setPhysicalBoundaryHelper(bc_helper);
        }
        Mat supplied_operator = nullptr;
        if (verify_supplied_operator_lifetime)
        {
            IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(supplied_operator,
                                                                                    problem_coefs,
                                                                                    u_bc_coefs,
                                                                                    0.0,
                                                                                    num_dofs_per_proc,
                                                                                    u_dof_index_idx,
                                                                                    p_dof_index_idx,
                                                                                    level);
            solver->setOperatorMat(supplied_operator);
        }
        solver->initializeSolverState(x_vec, b_vec);
        const bool use_multiplicative = solver->usesMultiplicativeShellSmootherForTest();
        const bool use_restrict_partition = solver->usesRestrictShellSmootherPartitionForTest();
        const KSP& petsc_ksp = solver->getPETScKSP();
        if (verify_live_operator_state_view || verify_raw_export_comparator_contract ||
            verify_supplied_operator_lifetime)
        {
            const auto view = solver->getLiveOperatorStateView();
            Mat ksp_operator = nullptr;
            Mat ksp_preconditioner = nullptr;
            ierr = KSPGetOperators(petsc_ksp, &ksp_operator, &ksp_preconditioner);
            IBTK_CHKERRQ(ierr);
            live_view_operator_identity = view.initialized && view.operator_mat == ksp_operator;

            PetscInt global_rows = 0, global_cols = 0;
            ierr = MatGetSize(view.operator_mat, &global_rows, &global_cols);
            IBTK_CHKERRQ(ierr);
            const int owned_dofs = view.num_dofs_per_proc ? (*view.num_dofs_per_proc)[rank] : -1;
            const int partition_dofs = view.locally_owned_velocity_dofs && view.locally_owned_pressure_dofs ?
                                           static_cast<int>(view.locally_owned_velocity_dofs->size() +
                                                            view.locally_owned_pressure_dofs->size()) :
                                           -1;
            live_view_ownership_valid =
                view.num_dofs_per_proc &&
                std::accumulate(view.num_dofs_per_proc->begin(), view.num_dofs_per_proc->end(), 0) == global_rows &&
                global_rows == global_cols;
            live_view_block_partition_valid =
                view.locally_owned_velocity_dofs && view.locally_owned_pressure_dofs &&
                std::is_sorted(view.locally_owned_velocity_dofs->begin(), view.locally_owned_velocity_dofs->end()) &&
                std::is_sorted(view.locally_owned_pressure_dofs->begin(), view.locally_owned_pressure_dofs->end()) &&
                owned_dofs == partition_dofs &&
                std::all_of(view.locally_owned_velocity_dofs->begin(),
                            view.locally_owned_velocity_dofs->end(),
                            [&](const PetscInt dof) { return solver->isVelocityDOF(static_cast<int>(dof)); }) &&
                std::all_of(view.locally_owned_pressure_dofs->begin(),
                            view.locally_owned_pressure_dofs->end(),
                            [&](const PetscInt dof) { return solver->isPressureDOF(static_cast<int>(dof)); });
            live_view_provenance_valid = view.level_number == 0 &&
                                         view.operator_was_provided == verify_supplied_operator_lifetime &&
                                         !view.includes_augmented_operator;
            supplied_operator_identity = !verify_supplied_operator_lifetime || view.operator_mat == supplied_operator;
            live_view_nullspace_gauge_valid = !view.velocity_nullspace_declared && view.pressure_nullspace_declared &&
                                              view.operator_nullspace_attached && view.zero_mean_pressure_correction;
            if (!live_view_operator_identity || !live_view_ownership_valid || !live_view_block_partition_valid ||
                !live_view_provenance_valid || !live_view_nullspace_gauge_valid || !supplied_operator_identity)
            {
                ++test_failures;
            }

            if (verify_raw_export_comparator_contract)
            {
                using namespace IBAMR::TestSupport;
                const CAVRawOperatorBundle raw_bundle =
                    captureCAVRawOperatorBundle(view, level, u_dof_index_idx, p_dof_index_idx, b_petsc);
                const std::string raw_prefix = "cav_raw_live_operator";
                writeCAVRawOperatorBundle(raw_bundle, raw_prefix);
                const CAVRawOperatorBundle reread_bundle = readCAVRawOperatorBundle(raw_prefix);
                CAVRawMappingSpec identity_mapping;
                identity_mapping.coordinate_tolerance = 0.0;
                identity_mapping.comparison_tolerance = 0.0;
                const CAVRawComparison raw_comparison =
                    compareCAVRawOperatorBundles(reread_bundle, raw_bundle, identity_mapping);
                raw_export_round_trip = raw_comparison.matched;
                raw_export_matrix_max_abs_error = raw_comparison.matrix_max_abs_error;

                const std::string borrowed_matrix_filename = "cav_raw_live_borrowed_matrix.mtx";
                writeCAVRawMatrixMarket(view.operator_mat, borrowed_matrix_filename);
                const CAVRawMatrixMarketData borrowed_matrix = readCAVRawMatrixMarket(borrowed_matrix_filename);
                raw_export_borrowed_matrix_round_trip =
                    borrowed_matrix.nrows == raw_bundle.nrows && borrowed_matrix.ncols == raw_bundle.ncols &&
                    borrowed_matrix.entries.size() == raw_bundle.matrix_entries.size() &&
                    std::equal(borrowed_matrix.entries.begin(),
                               borrowed_matrix.entries.end(),
                               raw_bundle.matrix_entries.begin(),
                               [](const CAVRawMatrixEntry& lhs, const CAVRawMatrixEntry& rhs)
                               { return lhs.row == rhs.row && lhs.column == rhs.column && lhs.value == rhs.value; });

                const std::string borrowed_vector_filename = "cav_raw_live_borrowed_vector.mtx";
                writeCAVRawVectorMarket(b_petsc, borrowed_vector_filename);
                raw_export_borrowed_vector_round_trip =
                    readCAVRawVectorMarket(borrowed_vector_filename) == raw_bundle.equation_values;

                CAVLiveExportManifest manifest;
                manifest.candidate_sha = "0123456789abcdef0123456789abcdef01234567";
                manifest.candidate_dirty = false;
                manifest.oracle_sha = "5b77344db6746269f8c77695c99e9043907ba74b";
                manifest.case_id = "native-raw-export-contract";
                manifest.dimension = NDIM;
                manifest.mpi_ranks = IBTK_MPI::getNodes();
                manifest.pressure_equation = "minus-div";
                manifest.pressure_equation_row_multiplier_to_oracle = -1.0;
                manifest.pressure_gauge = "zero-mean-correction";
                manifest.patch_seed_type = "VELOCITY_COMPONENT";
                manifest.closure_policy = closure_policy;
                manifest.seed_stride = seed_stride;
                manifest.traversal_order = "I_J";
                manifest.composition = use_multiplicative ? "multiplicative" : "additive";
                manifest.local_solver_backend = shell_pc_type;
                const std::string manifest_filename = "cav_raw_live_manifest.txt";
                writeCAVLiveExportManifest(manifest, manifest_filename);
                const CAVLiveExportManifest reread_manifest = readCAVLiveExportManifest(manifest_filename);
                raw_export_manifest_round_trip = sameCAVLiveExportManifest(manifest, reread_manifest);
                CAVLiveExportManifest mutated_manifest = reread_manifest;
                mutated_manifest.pressure_equation_row_multiplier_to_oracle = 1.0;
                raw_export_manifest_mutation_detected = !sameCAVLiveExportManifest(manifest, mutated_manifest);
                CAVLiveExportManifest invalid_manifest = reread_manifest;
                invalid_manifest.candidate_sha = "not-a-sha";
                try
                {
                    writeCAVLiveExportManifest(invalid_manifest, "cav_raw_invalid_manifest.txt");
                    raw_export_invalid_provenance_rejected = false;
                }
                catch (const std::runtime_error&)
                {
                }

                std::vector<PetscInt> exported_velocity_dofs;
                std::vector<PetscInt> exported_pressure_dofs;
                for (const CAVRawDofRecord& record : raw_bundle.dofs)
                {
                    if (record.kind == CAVRawDofKind::VELOCITY)
                        exported_velocity_dofs.push_back(record.dof);
                    else
                        exported_pressure_dofs.push_back(record.dof);
                }
                raw_export_global_dof_mapping = raw_bundle.nrows == 48 && raw_bundle.ncols == 48 &&
                                                raw_bundle.dofs.size() == 48 && view.locally_owned_velocity_dofs &&
                                                view.locally_owned_pressure_dofs &&
                                                exported_velocity_dofs == *view.locally_owned_velocity_dofs &&
                                                exported_pressure_dofs == *view.locally_owned_pressure_dofs;

                std::array<int, NDIM> pressure_index = {};
                pressure_index[0] = 1;
                pressure_index[1] = 1;
                std::array<int, NDIM> lower_x_index = pressure_index;
                std::array<int, NDIM> upper_x_index = pressure_index;
                std::array<int, NDIM> lower_y_index = pressure_index;
                std::array<int, NDIM> upper_y_index = pressure_index;
                ++upper_x_index[0];
                ++upper_y_index[1];
                const PetscInt pressure_dof = findCAVRawDof(raw_bundle, CAVRawDofKind::PRESSURE, -1, pressure_index);
                const PetscInt lower_x_dof = findCAVRawDof(raw_bundle, CAVRawDofKind::VELOCITY, 0, lower_x_index);
                const PetscInt upper_x_dof = findCAVRawDof(raw_bundle, CAVRawDofKind::VELOCITY, 0, upper_x_index);
                const PetscInt lower_y_dof = findCAVRawDof(raw_bundle, CAVRawDofKind::VELOCITY, 1, lower_y_index);
                const PetscInt upper_y_dof = findCAVRawDof(raw_bundle, CAVRawDofKind::VELOCITY, 1, upper_y_index);
                // The Stokes saddle-point operator uses -Div in its pressure equation. For h=1/4, the lower
                // x/y faces therefore have coefficient +4 and the upper faces coefficient -4.
                raw_export_pressure_row_minus_div =
                    pressure_dof >= 0 && lower_x_dof >= 0 && upper_x_dof >= 0 && lower_y_dof >= 0 && upper_y_dof >= 0 &&
                    countCAVRawMatrixRowNonzeros(raw_bundle, pressure_dof) == 4 &&
                    getCAVRawMatrixValue(raw_bundle, pressure_dof, lower_x_dof) == 4.0 &&
                    getCAVRawMatrixValue(raw_bundle, pressure_dof, upper_x_dof) == -4.0 &&
                    getCAVRawMatrixValue(raw_bundle, pressure_dof, lower_y_dof) == 4.0 &&
                    getCAVRawMatrixValue(raw_bundle, pressure_dof, upper_y_dof) == -4.0;

                if (!raw_export_round_trip || !raw_export_global_dof_mapping || !raw_export_pressure_row_minus_div ||
                    !raw_export_borrowed_matrix_round_trip || !raw_export_borrowed_vector_round_trip ||
                    !raw_export_manifest_round_trip || !raw_export_manifest_mutation_detected ||
                    !raw_export_invalid_provenance_rejected)
                {
                    ++test_failures;
                }
            }
        }
        if (!expected_shell_pc_name.empty())
        {
            PC shell_pc = nullptr;
            ierr = KSPGetPC(petsc_ksp, &shell_pc);
            IBTK_CHKERRQ(ierr);
            const char* actual_shell_pc_name = nullptr;
            ierr = PCShellGetName(shell_pc, &actual_shell_pc_name);
            IBTK_CHKERRQ(ierr);
            const std::string actual_shell_pc_name_str = actual_shell_pc_name ? actual_shell_pc_name : "";
            if (actual_shell_pc_name_str != expected_shell_pc_name) ++test_failures;
            pout << "expected_shell_pc_name = " << expected_shell_pc_name << "\n";
            pout << "actual_shell_pc_name = " << actual_shell_pc_name_str << "\n";
        }

        if (verify_rhs_boundary_adjustment)
        {
            Vec adjusted_b = nullptr;
            ierr = VecDuplicate(b_petsc, &adjusted_b);
            IBTK_CHKERRQ(ierr);
            solver->copyAdjustedRHS(adjusted_b, x_vec, b_vec);
            ierr = VecAXPY(adjusted_b, -1.0, b_petsc);
            IBTK_CHKERRQ(ierr);
            const double rhs_adjustment_inf_norm = vec_norm_inf(adjusted_b);
            const bool rhs_adjustment_detected = rhs_adjustment_inf_norm > 0.0;
            if (!rhs_adjustment_detected) ++test_failures;
            ierr = VecDestroy(&adjusted_b);
            IBTK_CHKERRQ(ierr);
            solver->deallocateSolverState();
            if (supplied_operator)
            {
                ierr = MatDestroy(&supplied_operator);
                IBTK_CHKERRQ(ierr);
            }

            const std::string solver_label = solver_type.empty() ? std::string("default") : solver_type;
            pout << "solver_type = " << solver_label << "\n";
            pout << "rhs_boundary_adjustment_detected = " << (rhs_adjustment_detected ? "true" : "false") << "\n";
            continue;
        }

        const auto& overlap_is = solver->getASMSubdomains();
        const auto& nonoverlap_is = solver->getASMNonoverlapSubdomains();

        Mat A_mat = nullptr;
        Mat pc_mat = nullptr;
        ierr = KSPGetOperators(petsc_ksp, &A_mat, &pc_mat);
        IBTK_CHKERRQ(ierr);

        int observer_call_count = 0;
        bool observer_data_valid = !verify_local_solve_observer;
        bool observer_disabled_silent = !verify_local_solve_observer;
        bool local_trace_round_trip = true;
        bool local_trace_selection_valid = !verify_local_solve_observer;
        bool local_trace_mutation_detected = !verify_local_solve_observer;
        bool local_trace_all_round_trip = true;
        bool local_trace_all_selection_valid = !verify_local_solve_observer;
        bool local_trace_all_omission_detected = !verify_local_solve_observer;
        int observer_call_count_after_clear = 0;
        std::vector<IBAMR::TestSupport::CAVLocalSolveTraceRecord> local_trace_records;
        std::vector<IBAMR::TestSupport::CAVLocalSolveTraceRecord> local_trace_all_records;
        const auto write_local_trace = [](std::vector<IBAMR::TestSupport::CAVLocalSolveTraceRecord>& records,
                                          const int sweep,
                                          const int ordinal,
                                          const std::string& prefix,
                                          Mat local_matrix,
                                          Vec local_rhs,
                                          Vec local_solution,
                                          Vec current_global_source)
        {
            const std::string trace_stem =
                prefix + "_sweep" + std::to_string(sweep) + "_patch" + std::to_string(ordinal);
            const IBAMR::TestSupport::CAVLocalSolveTraceRecord record = { sweep, ordinal, trace_stem };
            const bool round_trip = IBAMR::TestSupport::writeCAVLocalSolveTraceArtifacts(
                record, local_matrix, local_rhs, local_solution, current_global_source);
            records.push_back(record);
            return round_trip;
        };
        int stagewise_observer_call_count = 0;
        bool stagewise_patch_order_valid = !verify_stagewise_original_residual;
        double stagewise_local_rhs_max_error = 0.0;
        double stagewise_original_residual_max_error = 0.0;
        double stagewise_final_state_error_inf_norm = 0.0;
        Vec stagewise_iterate = nullptr;
        Vec stagewise_fresh_residual = nullptr;
        Vec stagewise_difference = nullptr;
        if (verify_stagewise_original_residual)
        {
            if (!use_multiplicative || !is_blas_lapack_case)
            {
                TBOX_ERROR(
                    "Stagewise original-residual verification requires a multiplicative BLAS/LAPACK "
                    "shell backend.\n");
            }
            ierr = VecDuplicate(b_petsc, &stagewise_iterate);
            IBTK_CHKERRQ(ierr);
            ierr = VecDuplicate(b_petsc, &stagewise_fresh_residual);
            IBTK_CHKERRQ(ierr);
            ierr = VecDuplicate(b_petsc, &stagewise_difference);
            IBTK_CHKERRQ(ierr);
            ierr = VecZeroEntries(stagewise_iterate);
            IBTK_CHKERRQ(ierr);
            stagewise_patch_order_valid = true;
            // Maintain an independent full-space iterate so every observed
            // pre-update residual can be checked against a fresh b - A x.
            solver->setShellSubdomainSolveObserver(
                [&](const int ordinal,
                    Mat /*local_matrix*/,
                    Vec local_rhs,
                    Vec local_solution,
                    Vec current_global_source)
                {
                    const std::size_t expected_ordinal = static_cast<std::size_t>(stagewise_observer_call_count++);
                    if (ordinal < 0 || static_cast<std::size_t>(ordinal) != expected_ordinal ||
                        expected_ordinal >= overlap_is.size())
                    {
                        stagewise_patch_order_valid = false;
                        return;
                    }

                    int observer_ierr = MatMult(A_mat, stagewise_iterate, stagewise_fresh_residual);
                    IBTK_CHKERRQ(observer_ierr);
                    observer_ierr = VecAYPX(stagewise_fresh_residual, -1.0, b_petsc);
                    IBTK_CHKERRQ(observer_ierr);
                    observer_ierr = VecCopy(current_global_source, stagewise_difference);
                    IBTK_CHKERRQ(observer_ierr);
                    observer_ierr = VecAXPY(stagewise_difference, -1.0, stagewise_fresh_residual);
                    IBTK_CHKERRQ(observer_ierr);
                    stagewise_original_residual_max_error =
                        std::max(stagewise_original_residual_max_error, vec_norm_inf(stagewise_difference));

                    const PetscScalar* source_values = nullptr;
                    const PetscScalar* rhs_values = nullptr;
                    observer_ierr = VecGetArrayRead(current_global_source, &source_values);
                    IBTK_CHKERRQ(observer_ierr);
                    observer_ierr = VecGetArrayRead(local_rhs, &rhs_values);
                    IBTK_CHKERRQ(observer_ierr);
                    const auto& overlap_dofs = overlap_is[expected_ordinal];
                    for (std::size_t local_dof = 0; local_dof < overlap_dofs.size(); ++local_dof)
                    {
                        stagewise_local_rhs_max_error =
                            std::max(stagewise_local_rhs_max_error,
                                     static_cast<double>(PetscAbsScalar(rhs_values[local_dof] -
                                                                        source_values[overlap_dofs[local_dof]])));
                    }
                    observer_ierr = VecRestoreArrayRead(local_rhs, &rhs_values);
                    IBTK_CHKERRQ(observer_ierr);
                    observer_ierr = VecRestoreArrayRead(current_global_source, &source_values);
                    IBTK_CHKERRQ(observer_ierr);

                    const auto& update_dofs = use_restrict_partition ? nonoverlap_is[expected_ordinal] : overlap_dofs;
                    const PetscScalar* solution_values = nullptr;
                    PetscScalar* iterate_values = nullptr;
                    observer_ierr = VecGetArrayRead(local_solution, &solution_values);
                    IBTK_CHKERRQ(observer_ierr);
                    observer_ierr = VecGetArray(stagewise_iterate, &iterate_values);
                    IBTK_CHKERRQ(observer_ierr);
                    for (const int update_dof : update_dofs)
                    {
                        const auto position = std::find(overlap_dofs.begin(), overlap_dofs.end(), update_dof);
                        if (position == overlap_dofs.end())
                        {
                            stagewise_patch_order_valid = false;
                            continue;
                        }
                        const std::size_t local_position =
                            static_cast<std::size_t>(std::distance(overlap_dofs.begin(), position));
                        iterate_values[update_dof] += solution_values[local_position];
                    }
                    observer_ierr = VecRestoreArray(stagewise_iterate, &iterate_values);
                    IBTK_CHKERRQ(observer_ierr);
                    observer_ierr = VecRestoreArrayRead(local_solution, &solution_values);
                    IBTK_CHKERRQ(observer_ierr);
                });
        }
        else if (verify_local_solve_observer)
        {
            solver->setShellSubdomainSolveObserver(
                [&](const int ordinal, Mat local_matrix, Vec local_rhs, Vec local_solution, Vec current_global_source)
                {
                    ++observer_call_count;
                    local_trace_round_trip = local_trace_round_trip && write_local_trace(local_trace_records,
                                                                                         0,
                                                                                         ordinal,
                                                                                         "cav_local_solve",
                                                                                         local_matrix,
                                                                                         local_rhs,
                                                                                         local_solution,
                                                                                         current_global_source);
                    PetscInt rows = 0, cols = 0, rhs_size = 0, solution_size = 0, global_source_size = 0;
                    int observer_ierr = MatGetSize(local_matrix, &rows, &cols);
                    IBTK_CHKERRQ(observer_ierr);
                    observer_ierr = VecGetSize(local_rhs, &rhs_size);
                    IBTK_CHKERRQ(observer_ierr);
                    observer_ierr = VecGetSize(local_solution, &solution_size);
                    IBTK_CHKERRQ(observer_ierr);
                    observer_ierr = VecGetSize(current_global_source, &global_source_size);
                    IBTK_CHKERRQ(observer_ierr);
                    Vec defect = nullptr;
                    observer_ierr = VecDuplicate(local_rhs, &defect);
                    IBTK_CHKERRQ(observer_ierr);
                    observer_ierr = MatMult(local_matrix, local_solution, defect);
                    IBTK_CHKERRQ(observer_ierr);
                    observer_ierr = VecAYPX(defect, -1.0, local_rhs);
                    IBTK_CHKERRQ(observer_ierr);
                    PetscReal defect_inf = 0.0, matrix_inf = 0.0, solution_inf = 0.0, rhs_inf = 0.0;
                    observer_ierr = VecNorm(defect, NORM_INFINITY, &defect_inf);
                    IBTK_CHKERRQ(observer_ierr);
                    observer_ierr = MatNorm(local_matrix, NORM_INFINITY, &matrix_inf);
                    IBTK_CHKERRQ(observer_ierr);
                    observer_ierr = VecNorm(local_solution, NORM_INFINITY, &solution_inf);
                    IBTK_CHKERRQ(observer_ierr);
                    observer_ierr = VecNorm(local_rhs, NORM_INFINITY, &rhs_inf);
                    IBTK_CHKERRQ(observer_ierr);
                    const double backward_error = defect_inf / (matrix_inf * solution_inf + rhs_inf + 1.0e-30);
                    const PetscInt expected_size = static_cast<PetscInt>(overlap_is.front().size());
                    PetscInt expected_global_size = 0;
                    observer_ierr = VecGetSize(b_petsc, &expected_global_size);
                    IBTK_CHKERRQ(observer_ierr);
                    observer_data_valid = ordinal == 0 && rows == expected_size && cols == expected_size &&
                                          rhs_size == expected_size && solution_size == expected_size &&
                                          global_source_size == expected_global_size && std::isfinite(backward_error) &&
                                          backward_error <= tol;
                    observer_ierr = VecDestroy(&defect);
                    IBTK_CHKERRQ(observer_ierr);
                },
                [](const int ordinal) { return ordinal == 0; });
        }

        apply_reference_shell_action(x_expected,
                                     b_petsc,
                                     A_mat,
                                     overlap_is,
                                     nonoverlap_is,
                                     use_multiplicative,
                                     use_restrict_partition,
                                     1.0,
                                     pressure_is);
        const double expected_inf_norm = vec_norm_inf(x_expected);
        if (!(expected_inf_norm > 0.0)) ++test_failures;

        const bool converged = solver->solveSystem(x_vec, b_vec);
        if (!converged) ++test_failures;
        if (verify_local_solve_observer)
        {
            IBAMR::TestSupport::writeCAVLocalSolveTraceIndex(local_trace_records, "cav_local_solve_trace.txt");
            const auto reread_trace_records =
                IBAMR::TestSupport::readCAVLocalSolveTraceIndex("cav_local_solve_trace.txt");
            local_trace_round_trip = local_trace_round_trip && IBAMR::TestSupport::sameCAVLocalSolveTraceIndex(
                                                                   local_trace_records, reread_trace_records);
            local_trace_selection_valid = local_trace_records.size() == 1 && local_trace_records.front().sweep == 0 &&
                                          local_trace_records.front().patch_ordinal == 0;
            auto mutated_trace_records = reread_trace_records;
            if (!mutated_trace_records.empty()) ++mutated_trace_records.front().patch_ordinal;
            local_trace_mutation_detected =
                !IBAMR::TestSupport::sameCAVLocalSolveTraceIndex(local_trace_records, mutated_trace_records);
            if (observer_call_count != 1 || !observer_data_valid || !local_trace_round_trip ||
                !local_trace_selection_valid || !local_trace_mutation_detected)
                ++test_failures;
            solver->setShellSubdomainSolveObserver({});
            observer_call_count_after_clear = observer_call_count;
        }
        if (verify_stagewise_original_residual)
        {
            if (stagewise_observer_call_count != static_cast<int>(overlap_is.size()) || !stagewise_patch_order_valid ||
                !(stagewise_local_rhs_max_error <= tol) || !(stagewise_original_residual_max_error <= tol))
            {
                ++test_failures;
            }
            solver->setShellSubdomainSolveObserver({});
        }
        IBAMR::StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            x_petsc, u_idx, u_dof_index_idx, p_idx, p_dof_index_idx, level);
        const double actual_inf_norm = vec_norm_inf(x_petsc);
        if (!(actual_inf_norm > 0.0)) ++test_failures;
        if (verify_stagewise_original_residual)
        {
            // The shell normalizes only the pressure state after completing
            // the sweep; it does not project intermediate residuals.
            normalize_pressure_gauge(stagewise_iterate, pressure_is);
            ierr = VecCopy(x_petsc, stagewise_difference);
            IBTK_CHKERRQ(ierr);
            ierr = VecAXPY(stagewise_difference, -1.0, stagewise_iterate);
            IBTK_CHKERRQ(ierr);
            stagewise_final_state_error_inf_norm = vec_norm_inf(stagewise_difference);
            if (!(stagewise_final_state_error_inf_norm <= tol)) ++test_failures;
            ierr = VecDestroy(&stagewise_difference);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&stagewise_fresh_residual);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&stagewise_iterate);
            IBTK_CHKERRQ(ierr);
        }

        ierr = VecCopy(x_petsc, x_diff);
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPY(x_diff, -1.0, x_expected);
        IBTK_CHKERRQ(ierr);
        const double error_inf_norm = vec_norm_inf(x_diff);
        const double reported_error_inf_norm = verify_reference_parity && error_inf_norm <= tol ? 0.0 : error_inf_norm;

        if (verify_reference_parity)
        {
            if (!(error_inf_norm <= tol)) ++test_failures;
        }
        else
        {
            if (!std::isfinite(error_inf_norm)) ++test_failures;
        }

        solver->deallocateSolverState();

        if (verify_supplied_operator_lifetime)
        {
            double supplied_operator_norm = 0.0;
            ierr = MatNorm(supplied_operator, NORM_INFINITY, &supplied_operator_norm);
            IBTK_CHKERRQ(ierr);
            supplied_operator_valid_after_deallocate =
                std::isfinite(supplied_operator_norm) && supplied_operator_norm > 0.0;
            if (!supplied_operator_valid_after_deallocate) ++test_failures;
        }

        double reinitialize_error_inf_norm = std::numeric_limits<double>::quiet_NaN();
        if (verify_reinitialize)
        {
            zero_solution_fields(level, u_idx, p_idx);
            solver->initializeSolverState(x_vec, b_vec);
            const bool reinitialize_converged = solver->solveSystem(x_vec, b_vec);
            if (!reinitialize_converged) ++test_failures;
            if (verify_local_solve_observer)
            {
                observer_disabled_silent = observer_call_count == observer_call_count_after_clear;
                if (!observer_disabled_silent) ++test_failures;
            }
            IBAMR::StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
                x_petsc, u_idx, u_dof_index_idx, p_idx, p_dof_index_idx, level);
            ierr = VecCopy(x_petsc, x_diff);
            IBTK_CHKERRQ(ierr);
            ierr = VecAXPY(x_diff, -1.0, x_expected);
            IBTK_CHKERRQ(ierr);
            reinitialize_error_inf_norm = vec_norm_inf(x_diff);
            if (!(reinitialize_error_inf_norm <= tol)) ++test_failures;
            if (verify_local_solve_observer)
            {
                zero_solution_fields(level, u_idx, p_idx);
                solver->setShellSubdomainSolveObserver(
                    [&](const int ordinal,
                        Mat local_matrix,
                        Vec local_rhs,
                        Vec local_solution,
                        Vec current_global_source)
                    {
                        local_trace_all_round_trip =
                            local_trace_all_round_trip && write_local_trace(local_trace_all_records,
                                                                            0,
                                                                            ordinal,
                                                                            "cav_local_solve_all",
                                                                            local_matrix,
                                                                            local_rhs,
                                                                            local_solution,
                                                                            current_global_source);
                    });
                const bool all_trace_converged = solver->solveSystem(x_vec, b_vec);
                solver->setShellSubdomainSolveObserver({});
                IBAMR::TestSupport::writeCAVLocalSolveTraceIndex(local_trace_all_records,
                                                                 "cav_local_solve_all_trace.txt");
                const auto reread_all_records =
                    IBAMR::TestSupport::readCAVLocalSolveTraceIndex("cav_local_solve_all_trace.txt");
                local_trace_all_round_trip =
                    local_trace_all_round_trip &&
                    IBAMR::TestSupport::sameCAVLocalSolveTraceIndex(local_trace_all_records, reread_all_records);
                local_trace_all_selection_valid = local_trace_all_records.size() == overlap_is.size();
                for (std::size_t k = 0; local_trace_all_selection_valid && k < local_trace_all_records.size(); ++k)
                {
                    local_trace_all_selection_valid = local_trace_all_records[k].sweep == 0 &&
                                                      local_trace_all_records[k].patch_ordinal == static_cast<int>(k);
                }
                auto omitted_all_records = reread_all_records;
                if (!omitted_all_records.empty()) omitted_all_records.pop_back();
                local_trace_all_omission_detected =
                    !IBAMR::TestSupport::sameCAVLocalSolveTraceIndex(local_trace_all_records, omitted_all_records);
                IBAMR::StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
                    x_petsc, u_idx, u_dof_index_idx, p_idx, p_dof_index_idx, level);
                ierr = VecCopy(x_petsc, x_diff);
                IBTK_CHKERRQ(ierr);
                ierr = VecAXPY(x_diff, -1.0, x_expected);
                IBTK_CHKERRQ(ierr);
                const bool all_trace_solution_valid = vec_norm_inf(x_diff) <= tol;
                if (!all_trace_converged || !local_trace_all_round_trip || !local_trace_all_selection_valid ||
                    !local_trace_all_omission_detected || !all_trace_solution_valid)
                    ++test_failures;
            }
            solver->deallocateSolverState();
        }

        if (supplied_operator)
        {
            ierr = MatDestroy(&supplied_operator);
            IBTK_CHKERRQ(ierr);
        }

        if (verify_live_operator_state_view)
        {
            const auto view = solver->getLiveOperatorStateView();
            live_view_post_deallocate_empty = !view.initialized && !view.operator_mat && !view.num_dofs_per_proc &&
                                              !view.locally_owned_velocity_dofs && !view.locally_owned_pressure_dofs;
            if (!live_view_post_deallocate_empty) ++test_failures;
        }
        if (verify_supplied_operator_lifetime)
        {
            pout << "supplied_operator_identity = " << (supplied_operator_identity ? "true" : "false") << "\n";
            pout << "supplied_operator_valid_after_deallocate = "
                 << (supplied_operator_valid_after_deallocate ? "true" : "false") << "\n";
        }

        const std::string solver_label = solver_type.empty() ? std::string("default") : solver_type;
        pout << "solver_type = " << solver_label << "\n";
        if (verify_reference_parity)
        {
            pout << "expected_inf_norm = " << expected_inf_norm << "\n";
            pout << "actual_inf_norm = " << actual_inf_norm << "\n";
            pout << "error_inf_norm = " << reported_error_inf_norm << "\n";
        }
        else if (report_actual_inf_norm)
        {
            pout << "actual_inf_norm = " << actual_inf_norm << "\n";
        }
        if (verify_reinitialize)
        {
            pout << "reinitialize_error_inf_norm = "
                 << (reinitialize_error_inf_norm <= tol ? 0.0 : reinitialize_error_inf_norm) << "\n";
        }
        if (verify_local_solve_observer)
        {
            pout << "observer_call_count = " << observer_call_count << "\n";
            pout << "observer_data_valid = " << (observer_data_valid ? "true" : "false") << "\n";
            pout << "observer_disabled_silent = " << (observer_disabled_silent ? "true" : "false") << "\n";
            pout << "local_trace_round_trip = " << (local_trace_round_trip ? "true" : "false") << "\n";
            pout << "local_trace_selection_valid = " << (local_trace_selection_valid ? "true" : "false") << "\n";
            pout << "local_trace_mutation_detected = " << (local_trace_mutation_detected ? "true" : "false") << "\n";
            pout << "local_trace_all_round_trip = " << (local_trace_all_round_trip ? "true" : "false") << "\n";
            pout << "local_trace_all_selection_valid = " << (local_trace_all_selection_valid ? "true" : "false")
                 << "\n";
            pout << "local_trace_all_omission_detected = " << (local_trace_all_omission_detected ? "true" : "false")
                 << "\n";
        }
        if (verify_stagewise_original_residual)
        {
            pout << "stagewise_observer_call_count = " << stagewise_observer_call_count << "\n";
            pout << "stagewise_patch_order_valid = " << (stagewise_patch_order_valid ? "true" : "false") << "\n";
            pout << "stagewise_local_rhs_max_error = "
                 << (stagewise_local_rhs_max_error <= tol ? 0.0 : stagewise_local_rhs_max_error) << "\n";
            pout << "stagewise_original_residual_max_error = "
                 << (stagewise_original_residual_max_error <= tol ? 0.0 : stagewise_original_residual_max_error)
                 << "\n";
            pout << "stagewise_final_state_error_inf_norm = "
                 << (stagewise_final_state_error_inf_norm <= tol ? 0.0 : stagewise_final_state_error_inf_norm) << "\n";
        }
        if (verify_live_operator_state_view)
        {
            pout << "live_view_preinit_empty = " << (live_view_preinit_empty ? "true" : "false") << "\n";
            pout << "live_view_operator_identity = " << (live_view_operator_identity ? "true" : "false") << "\n";
            pout << "live_view_ownership_valid = " << (live_view_ownership_valid ? "true" : "false") << "\n";
            pout << "live_view_block_partition_valid = " << (live_view_block_partition_valid ? "true" : "false")
                 << "\n";
            pout << "live_view_provenance_valid = " << (live_view_provenance_valid ? "true" : "false") << "\n";
            pout << "live_view_nullspace_gauge_valid = " << (live_view_nullspace_gauge_valid ? "true" : "false")
                 << "\n";
            pout << "live_view_post_deallocate_empty = " << (live_view_post_deallocate_empty ? "true" : "false")
                 << "\n";
        }
        if (verify_raw_export_comparator_contract)
        {
            pout << "raw_export_round_trip = " << (raw_export_round_trip ? "true" : "false") << "\n";
            pout << "raw_export_global_dof_mapping = " << (raw_export_global_dof_mapping ? "true" : "false") << "\n";
            pout << "raw_export_pressure_row_minus_div = " << (raw_export_pressure_row_minus_div ? "true" : "false")
                 << "\n";
            pout << "raw_export_borrowed_matrix_round_trip = "
                 << (raw_export_borrowed_matrix_round_trip ? "true" : "false") << "\n";
            pout << "raw_export_borrowed_vector_round_trip = "
                 << (raw_export_borrowed_vector_round_trip ? "true" : "false") << "\n";
            pout << "raw_export_manifest_round_trip = " << (raw_export_manifest_round_trip ? "true" : "false") << "\n";
            pout << "raw_export_manifest_mutation_detected = "
                 << (raw_export_manifest_mutation_detected ? "true" : "false") << "\n";
            pout << "raw_export_invalid_provenance_rejected = "
                 << (raw_export_invalid_provenance_rejected ? "true" : "false") << "\n";
            pout << "raw_export_matrix_max_abs_error = " << raw_export_matrix_max_abs_error << "\n";
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

    if (verify_shared_composition)
    {
        pout << "shared_additive_event_count = " << shared_additive_event_count << "\n";
        pout << "shared_additive_order_valid = " << (shared_additive_order_valid ? "true" : "false") << "\n";
        pout << "shared_multiplicative_event_count = " << shared_multiplicative_event_count << "\n";
        pout << "shared_multiplicative_order_valid = " << (shared_multiplicative_order_valid ? "true" : "false")
             << "\n";
        pout << "affected_row_update_row_visits = " << affected_row_update_row_visits << "\n";
        pout << "affected_row_stagewise_residual_valid = " << (affected_row_stagewise_residual_valid ? "true" : "false")
             << "\n";
        pout << "affected_row_update_valid = " << (affected_row_update_valid ? "true" : "false") << "\n";
    }
    pout << "shell_pc_type = " << (has_shell_pc_type ? shell_pc_type : std::string("<default>")) << "\n";
    pout << "coupling_aware_asm_closure_policy = " << closure_policy << "\n";
    pout << "parity_tol = " << tol << "\n";
    pout << "test_failures = " << test_failures << "\n";
    return test_failures > 0 ? 1 : 0;
}
