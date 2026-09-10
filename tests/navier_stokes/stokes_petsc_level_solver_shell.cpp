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
#include <ibamr/StaggeredStokesPETScVecUtilities.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/private/PETScLevelSolverShellBackend.h>

#include <tbox/MemoryDatabase.h>

#include <CellData.h>
#include <CellVariable.h>
#include <LocationIndexRobinBcCoefs.h>
#include <PoissonSpecifications.h>
#include <SAMRAIVectorReal.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <memory>
#include <numeric>

#include "../tests.h"

#include <ibamr/app_namespaces.h>

namespace
{
double
norm_inf(Vec x)
{
    PetscReal value = 0.0;
    const int ierr = VecNorm(x, NORM_INFINITY, &value);
    IBTK_CHKERRQ(ierr);
    return value;
}

struct RowMatrix
{
    std::vector<std::vector<PetscInt>> columns{ { 0, 1 }, { 0, 1, 2 }, { 1, 2, 3 }, { 2, 3 } };
    std::vector<std::vector<PetscScalar>> values{ { 9, 9 }, { 9, 9, 9 }, { 9, 9, 9 }, { 9, 9 } };
    int row_reads = 0, multiplies = 0;
};

PetscErrorCode
get_test_row(Mat mat, PetscInt row, PetscInt* n, const PetscInt** columns, const PetscScalar** values)
{
    RowMatrix* context = nullptr;
    int ierr = MatShellGetContext(mat, &context);
    IBTK_CHKERRQ(ierr);
    ++context->row_reads;
    *n = static_cast<PetscInt>(context->columns[row].size());
    *columns = context->columns[row].data();
    if (values)
    {
        *values = context->values[row].data();
    }
    return 0;
}

PetscErrorCode
restore_test_row(Mat, PetscInt, PetscInt*, const PetscInt**, const PetscScalar**)
{
    return 0;
}

PetscErrorCode
multiply_test_matrix(Mat mat, Vec x, Vec y)
{
    RowMatrix* context = nullptr;
    int ierr = MatShellGetContext(mat, &context);
    IBTK_CHKERRQ(ierr);
    ++context->multiplies;
    const PetscScalar* input = nullptr;
    PetscScalar* output = nullptr;
    ierr = VecGetArrayRead(x, &input);
    IBTK_CHKERRQ(ierr);
    ierr = VecGetArray(y, &output);
    IBTK_CHKERRQ(ierr);
    for (std::size_t row = 0; row < context->columns.size(); ++row)
    {
        output[row] = 0.0;
        for (std::size_t j = 0; j < context->columns[row].size(); ++j)
        {
            output[row] += context->values[row][j] * input[context->columns[row][j]];
        }
    }
    ierr = VecRestoreArray(y, &output);
    IBTK_CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(x, &input);
    IBTK_CHKERRQ(ierr);
    return 0;
}

// Prescribed local corrections isolate the real shared composer's residual action.
class StageBackend : public PETScLevelSolverShellBackend
{
public:
    /*! \brief Record residual samples in test-owned storage. */
    StageBackend(std::vector<PetscScalar>& samples, std::vector<std::size_t>& visits, std::size_t stages)
        : d_samples(samples), d_visits(visits), d_stages(stages)
    {
    }
    ~StageBackend() override
    {
        deallocateSolverState();
    }
    void
    initializeSolverState(Mat mat,
                          Vec x,
                          Vec b,
                          const std::vector<IS>&,
                          const std::vector<IS>&,
                          const std::string&,
                          bool multiplicative = false,
                          PETScLevelSolverShellTraversal traversal = PETScLevelSolverShellTraversal::FORWARD) override
    {
        deallocateSolverState();
        initializeComposition(mat, x, b, multiplicative, traversal);
        finalizeComposition();
    }
    void deallocateSolverState() override
    {
        deallocateComposition();
    }

protected:
    std::size_t getNumberOfSubdomains() const override
    {
        return d_stages;
    }
    void beginSubdomainRhs(std::size_t i, Vec source) override
    {
        PetscScalar value = 0.0;
        int ierr = VecGetValues(source, 1, d_dofs[i].data(), &value);
        IBTK_CHKERRQ(ierr);
        d_samples.push_back(value);
        d_visits.push_back(i);
    }
    void endSubdomainRhs(std::size_t, Vec) override
    {
    }
    void solveSubdomain(std::size_t) override
    {
    }
    void accumulateSubdomainCorrection(std::size_t i, Vec y) override
    {
        int ierr = VecSetValue(y, d_dofs[i][0], 1.0, ADD_VALUES);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyBegin(y);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyEnd(y);
        IBTK_CHKERRQ(ierr);
    }
    const std::vector<PetscInt>& getSubdomainCorrectionDofs(std::size_t i) const override
    {
        return d_dofs[i];
    }
    void copySubdomainCorrection(std::size_t, PetscScalar* values) override
    {
        values[0] = 1.0;
    }

private:
    std::vector<PetscScalar>& d_samples;
    std::vector<std::size_t>& d_visits;
    const std::size_t d_stages;
    const std::vector<std::vector<PetscInt>> d_dofs{ { 0 }, { 1 }, { 2 } };
};

int
check_stages(const bool fallback,
             const bool invalidate,
             const PETScLevelSolverShellTraversal traversal = PETScLevelSolverShellTraversal::FORWARD,
             const std::size_t stages = 3)
{
    RowMatrix context;
    Mat mat = nullptr;
    Vec x = nullptr, y = nullptr;
    int ierr = MatCreateShell(PETSC_COMM_SELF, 4, 4, 4, 4, &context, &mat);
    IBTK_CHKERRQ(ierr);
    ierr = MatShellSetOperation(mat, MATOP_MULT, reinterpret_cast<void (*)(void)>(multiply_test_matrix));
    IBTK_CHKERRQ(ierr);
    if (!fallback)
    {
        ierr = MatShellSetOperation(mat, MATOP_GET_ROW, reinterpret_cast<void (*)(void)>(get_test_row));
        IBTK_CHKERRQ(ierr);
        ierr = MatShellSetOperation(mat, MATOP_RESTORE_ROW, reinterpret_cast<void (*)(void)>(restore_test_row));
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecCreateSeq(PETSC_COMM_SELF, 4, &x);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(x, &y);
    IBTK_CHKERRQ(ierr);
    const PetscInt indices[] = { 0, 1, 2, 3 };
    const PetscScalar rhs[] = { 10, 20, 30, 40 };
    ierr = VecSetValues(x, 4, indices, rhs, INSERT_VALUES);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyBegin(x);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x);
    IBTK_CHKERRQ(ierr);
    std::vector<PetscScalar> samples;
    std::vector<std::size_t> visits;
    StageBackend backend(samples, visits, stages);
    std::vector<std::size_t> expected_visits;
    std::vector<PetscScalar> expected_samples;
    std::vector<PetscScalar> expected_correction(4, 0.0);
    int expected_reads = 0;
    if (stages == 1)
    {
        expected_visits = { 0 };
        expected_samples = { 10 };
        expected_correction[0] = 1;
    }
    else if (stages == 3)
    {
        if (traversal == PETScLevelSolverShellTraversal::FORWARD)
        {
            expected_visits = { 0, 1, 2 };
            expected_samples = { 10, 21, 31 };
            expected_correction = { 1, 1, 1, 0 };
            expected_reads = 5;
        }
        else if (traversal == PETScLevelSolverShellTraversal::REVERSE)
        {
            expected_visits = { 2, 1, 0 };
            expected_samples = { 30, 21, 11 };
            expected_correction = { 1, 1, 1, 0 };
            expected_reads = 6;
        }
        else
        {
            expected_visits = { 0, 1, 2, 1, 0 };
            expected_samples = { 10, 21, 31, 20, 10 };
            expected_correction = { 2, 2, 1, 0 };
            expected_reads = 11;
        }
    }
    const int expected_multiplies = expected_visits.empty() ? 0 : static_cast<int>(expected_visits.size()) - 1;
    int failures = 0;
    for (int cycle = 0; cycle < 2; ++cycle)
    {
        backend.initializeSolverState(mat, x, x, {}, {}, "", true, traversal);
        // Application must read values at cached offsets, not an update-matrix copy.
        context.values = { { 2, -1 }, { -1, 2, -1 }, { -1, 2, -1 }, { -1, 2 } };
        if (invalidate)
        {
            std::swap(context.columns[0][0], context.columns[0][1]);
        }
        for (int application = 0; application < 2; ++application)
        {
            samples.clear();
            visits.clear();
            context.row_reads = context.multiplies = 0;
            backend.apply(x, y);
            if (invalidate)
            {
                backend.deallocateSolverState();
                ierr = MatDestroy(&mat);
                IBTK_CHKERRQ(ierr);
                ierr = VecDestroy(&x);
                IBTK_CHKERRQ(ierr);
                ierr = VecDestroy(&y);
                IBTK_CHKERRQ(ierr);
                return 0;
            }
            if (visits != expected_visits || samples != expected_samples ||
                context.row_reads != (fallback ? 0 : expected_reads) ||
                context.multiplies != (fallback ? expected_multiplies : 0))
            {
                ++failures;
            }
            const PetscScalar* values = nullptr;
            ierr = VecGetArrayRead(y, &values);
            IBTK_CHKERRQ(ierr);
            if (!std::equal(expected_correction.begin(), expected_correction.end(), values))
            {
                ++failures;
            }
            if (stages == 3 && cycle == 0 && application == 0)
            {
                plog << "stage_rhs =";
                for (PetscScalar sample : samples)
                {
                    plog << ' ' << sample;
                }
                plog << "\ncorrection = " << values[0] << ' ' << values[1] << ' ' << values[2] << ' ' << values[3]
                     << "\nrow_reads = " << context.row_reads << "\nmatmult_calls = " << context.multiplies << '\n';
            }
            ierr = VecRestoreArrayRead(y, &values);
            IBTK_CHKERRQ(ierr);
        }
        backend.deallocateSolverState();
        backend.deallocateSolverState();
    }
    ierr = MatDestroy(&mat);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&y);
    IBTK_CHKERRQ(ierr);
    return failures;
}

int
check_hand_solve(const bool blas,
                 const PETScLevelSolverShellTraversal traversal = PETScLevelSolverShellTraversal::FORWARD,
                 const bool traversal_case = false)
{
    const int rank = IBTK_MPI::getRank();
    Mat mat = nullptr;
    Vec x = nullptr, y = nullptr, expected = nullptr;
    int ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 3, 3, 3, nullptr, 3, nullptr, &mat);
    IBTK_CHKERRQ(ierr);
    PetscInt lo = 0, hi = 0;
    ierr = MatGetOwnershipRange(mat, &lo, &hi);
    IBTK_CHKERRQ(ierr);
    for (PetscInt row = lo; row < hi; ++row)
    {
        for (PetscInt col = std::max<PetscInt>(0, row - 1); col <= std::min<PetscInt>(2, row + 1); ++col)
        {
            ierr = MatSetValue(mat, row, col, row == col ? 2.0 : -1.0, INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
    }
    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatCreateVecs(mat, &x, &y);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(y, &expected);
    IBTK_CHKERRQ(ierr);
    for (PetscInt row = lo; row < hi; ++row)
    {
        ierr = VecSetValue(x, row, row + 1.0, INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
        const bool parallel = traversal_case && IBTK_MPI::getNodes() == 2;
        const std::vector<PetscScalar> target =
            traversal == PETScLevelSolverShellTraversal::REVERSE ?
                (parallel ? std::vector<PetscScalar>{ 40.0 / 9.0, 41.0 / 9.0, 8.0 / 3.0 } :
                            std::vector<PetscScalar>{ 20.0 / 9.0, 31.0 / 9.0, 8.0 / 3.0 }) :
            traversal == PETScLevelSolverShellTraversal::SYMMETRIC ?
                (parallel ? std::vector<PetscScalar>{ 64.0 / 27.0, 107.0 / 27.0, 32.0 / 9.0 } :
                            std::vector<PetscScalar>{ 64.0 / 27.0, 101.0 / 27.0, 28.0 / 9.0 }) :
                (parallel ? std::vector<PetscScalar>{ 8.0 / 3.0, 37.0 / 9.0, 32.0 / 9.0 } :
                            std::vector<PetscScalar>{ 4.0 / 3.0, 29.0 / 9.0, 28.0 / 9.0 });
        ierr = VecSetValue(expected, row, target[row], INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(x);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyBegin(expected);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(expected);
    IBTK_CHKERRQ(ierr);
    // Rank 1 contributes only to stage zero; traversal cases duplicate rank zero
    // there to distinguish the summed parallel correction from the serial result.
    const std::vector<std::vector<PetscInt>> dofs =
        rank == 0 ? std::vector<std::vector<PetscInt>>{ { 0, 1 }, { 1, 2 } } :
                    std::vector<std::vector<PetscInt>>{ traversal_case ? std::vector<PetscInt>{ 0, 1 } :
                                                                         std::vector<PetscInt>{ 1, 2 } };
    std::vector<IS> overlap(dofs.size()), partition(dofs.size());
    for (std::size_t i = 0; i < dofs.size(); ++i)
    {
        ierr = ISCreateGeneral(PETSC_COMM_SELF, 2, dofs[i].data(), PETSC_COPY_VALUES, &overlap[i]);
        IBTK_CHKERRQ(ierr);
        ierr = ISCreateGeneral(PETSC_COMM_SELF, 0, nullptr, PETSC_COPY_VALUES, &partition[i]);
        IBTK_CHKERRQ(ierr);
    }
    Pointer<MemoryDatabase> db = new MemoryDatabase("hand");
    db->putString("blas_lapack_subdomain_solver_type", "lu");
    std::unique_ptr<PETScLevelSolverShellBackend> backend =
        PETScLevelSolverShellBackendManager::get_manager().allocateBackend(blas ? "blas-lapack" : "petsc", db);
    int failures = 0;
    for (int cycle = 0; cycle < 2; ++cycle)
    {
        backend->initializeSolverState(mat, x, y, overlap, partition, "r09_hand", true, traversal);
        for (int application = 0; application < 2; ++application)
        {
            backend->apply(x, y);
            ierr = VecAXPY(y, -1.0, expected);
            IBTK_CHKERRQ(ierr);
            const double error = norm_inf(y);
            if (!std::isfinite(error) || error > 1.0e-12)
            {
                ++failures;
            }
            if (cycle == 0 && application == 0)
            {
                plog << "hand_error = " << error << '\n';
            }
        }
        backend->deallocateSolverState();
        backend->deallocateSolverState();
    }
    for (IS& is : overlap)
    {
        ierr = ISDestroy(&is);
        IBTK_CHKERRQ(ierr);
    }
    for (IS& is : partition)
    {
        ierr = ISDestroy(&is);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatDestroy(&mat);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&y);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&expected);
    IBTK_CHKERRQ(ierr);
    return failures;
}

// Gather the RHS independently of the backend's restriction/prolongation
// scatters. Additive writes use the partition; multiplicative visits add full corrections.
void
reference_action(Mat mat,
                 Vec rhs,
                 Vec result,
                 const std::vector<IS>& overlap,
                 const std::vector<IS>& partition,
                 const bool multiplicative,
                 const PETScLevelSolverShellTraversal traversal)
{
    Vec residual = nullptr, gathered = nullptr;
    VecScatter gather = nullptr;
    int ierr = VecDuplicate(rhs, &residual);
    IBTK_CHKERRQ(ierr);
    ierr = VecScatterCreateToAll(rhs, &gather, &gathered);
    IBTK_CHKERRQ(ierr);
    ierr = VecSet(result, 0.0);
    IBTK_CHKERRQ(ierr);
    Mat* submat = nullptr;
    ierr = MatCreateSubMatrices(
        mat, static_cast<PetscInt>(overlap.size()), overlap.data(), overlap.data(), MAT_INITIAL_MATRIX, &submat);
    IBTK_CHKERRQ(ierr);
    if (!multiplicative)
    {
        ierr = VecScatterBegin(gather, rhs, gathered, INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
        ierr = VecScatterEnd(gather, rhs, gathered, INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
    }
    const int n_stages = IBTK_MPI::maxReduction(static_cast<int>(overlap.size()));
    std::vector<int> visits(n_stages);
    std::iota(visits.begin(), visits.end(), 0);
    if (multiplicative && traversal == PETScLevelSolverShellTraversal::REVERSE)
    {
        std::reverse(visits.begin(), visits.end());
    }
    else if (multiplicative && traversal == PETScLevelSolverShellTraversal::SYMMETRIC)
    {
        for (int i = n_stages - 2; i >= 0; --i)
        {
            visits.push_back(i);
        }
    }
    for (const int i : visits)
    {
        if (multiplicative)
        {
            // Recompute the global original residual independently before each stage.
            ierr = MatMult(mat, result, residual);
            IBTK_CHKERRQ(ierr);
            ierr = VecAYPX(residual, -1.0, rhs);
            IBTK_CHKERRQ(ierr);
            ierr = VecScatterBegin(gather, residual, gathered, INSERT_VALUES, SCATTER_FORWARD);
            IBTK_CHKERRQ(ierr);
            ierr = VecScatterEnd(gather, residual, gathered, INSERT_VALUES, SCATTER_FORWARD);
            IBTK_CHKERRQ(ierr);
        }
        if (i < static_cast<int>(overlap.size()))
        {
            PetscInt n = 0, m = 0;
            const PetscInt* indices = nullptr;
            const PetscInt* owned = nullptr;
            ierr = ISGetLocalSize(overlap[i], &n);
            IBTK_CHKERRQ(ierr);
            ierr = ISGetLocalSize(partition[i], &m);
            IBTK_CHKERRQ(ierr);
            ierr = ISGetIndices(overlap[i], &indices);
            IBTK_CHKERRQ(ierr);
            ierr = ISGetIndices(partition[i], &owned);
            IBTK_CHKERRQ(ierr);
            Vec local_rhs = nullptr, local_solution = nullptr;
            ierr = MatCreateVecs(submat[i], &local_solution, &local_rhs);
            IBTK_CHKERRQ(ierr);
            const PetscScalar* global_values = nullptr;
            PetscScalar* local_values = nullptr;
            ierr = VecGetArrayRead(gathered, &global_values);
            IBTK_CHKERRQ(ierr);
            ierr = VecGetArray(local_rhs, &local_values);
            IBTK_CHKERRQ(ierr);
            for (PetscInt j = 0; j < n; ++j)
            {
                local_values[j] = global_values[indices[j]];
            }
            ierr = VecRestoreArray(local_rhs, &local_values);
            IBTK_CHKERRQ(ierr);
            ierr = VecRestoreArrayRead(gathered, &global_values);
            IBTK_CHKERRQ(ierr);
            KSP ksp = nullptr;
            PC pc = nullptr;
            ierr = KSPCreate(PETSC_COMM_SELF, &ksp);
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetOperators(ksp, submat[i], submat[i]);
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetType(ksp, KSPPREONLY);
            IBTK_CHKERRQ(ierr);
            ierr = KSPGetPC(ksp, &pc);
            IBTK_CHKERRQ(ierr);
            ierr = PCSetType(pc, PCSVD);
            IBTK_CHKERRQ(ierr);
            ierr = KSPSolve(ksp, local_rhs, local_solution);
            IBTK_CHKERRQ(ierr);
            const PetscScalar* solution_values = nullptr;
            ierr = VecGetArrayRead(local_solution, &solution_values);
            IBTK_CHKERRQ(ierr);
            for (PetscInt j = 0; j < (multiplicative ? n : m); ++j)
            {
                const PetscInt position =
                    multiplicative ? j :
                                     static_cast<PetscInt>(std::lower_bound(indices, indices + n, owned[j]) - indices);
                TBOX_ASSERT(position < n && (multiplicative || indices[position] == owned[j]));
                ierr = VecSetValue(result,
                                   multiplicative ? indices[j] : owned[j],
                                   solution_values[position],
                                   multiplicative ? ADD_VALUES : INSERT_VALUES);
                IBTK_CHKERRQ(ierr);
            }
            ierr = VecRestoreArrayRead(local_solution, &solution_values);
            IBTK_CHKERRQ(ierr);
            ierr = ISRestoreIndices(partition[i], &owned);
            IBTK_CHKERRQ(ierr);
            ierr = ISRestoreIndices(overlap[i], &indices);
            IBTK_CHKERRQ(ierr);
            ierr = KSPDestroy(&ksp);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&local_rhs);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&local_solution);
            IBTK_CHKERRQ(ierr);
        }
        if (multiplicative)
        {
            ierr = VecAssemblyBegin(result);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyEnd(result);
            IBTK_CHKERRQ(ierr);
        }
    }
    ierr = VecAssemblyBegin(result);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(result);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroyMatrices(static_cast<PetscInt>(overlap.size()), &submat);
    IBTK_CHKERRQ(ierr);
    ierr = VecScatterDestroy(&gather);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&gathered);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&residual);
    IBTK_CHKERRQ(ierr);
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    Pointer<Logger::Appender> appender = new TestAppender();
    Logger::getInstance()->setAbortAppender(appender);
    Pointer<AppInitializer> app = new AppInitializer(argc, argv, "output");
    Pointer<Database> input = app->getInputDatabase();
    Pointer<Database> test = input->getDatabase("test");
    const bool boundary = test->getBoolWithDefault("boundary", false);
    const bool lifetime = test->getBoolWithDefault("lifetime", false);
    const bool invalid = test->getBoolWithDefault("invalid", false);
    const std::string shell_type = test->getStringWithDefault("shell_pc_type", "multiplicative");
    const bool multiplicative = shell_type.rfind("multiplicative", 0) == 0;
    std::string traversal_name = test->getStringWithDefault("shell_pc_subdomain_traversal", "FORWARD");
    std::transform(traversal_name.begin(),
                   traversal_name.end(),
                   traversal_name.begin(),
                   [](const unsigned char c) { return static_cast<char>(std::toupper(c)); });
    const PETScLevelSolverShellTraversal traversal =
        traversal_name == "REVERSE"   ? PETScLevelSolverShellTraversal::REVERSE :
        traversal_name == "SYMMETRIC" ? PETScLevelSolverShellTraversal::SYMMETRIC :
                                        PETScLevelSolverShellTraversal::FORWARD;
    if (test->getBoolWithDefault("traversal_stages", false))
    {
        int failures = 0;
        for (PETScLevelSolverShellTraversal order : { PETScLevelSolverShellTraversal::FORWARD,
                                                      PETScLevelSolverShellTraversal::REVERSE,
                                                      PETScLevelSolverShellTraversal::SYMMETRIC })
        {
            for (const bool fallback : { false, true })
            {
                failures += check_stages(fallback, false, order);
                failures += check_stages(fallback, false, order, 0);
                failures += check_stages(fallback, false, order, 1);
            }
        }
        plog << "failures = " << failures << '\n';
        return failures ? 1 : 0;
    }
    if (test->getBoolWithDefault("stages", false))
    {
        return check_stages(test->getBoolWithDefault("fallback", false), test->getBoolWithDefault("invalidate", false));
    }
    const int hand_failures = test->getBoolWithDefault("hand_solve", false) ?
                                  check_hand_solve(shell_type == "multiplicative-blas-lapack",
                                                   traversal,
                                                   test->keyExists("shell_pc_subdomain_traversal")) :
                                  0;
    const bool diagonal_operator = test->getBoolWithDefault("diagonal_operator", false);
    const double small_diagonal = test->getDoubleWithDefault("small_diagonal", 1.0e-12);
    const bool all_blas_modes = test->getBoolWithDefault("all_blas_modes", false);
    const std::vector<std::string> solver_types =
        all_blas_modes ?
            std::vector<std::string>{ "", "svd", "lu", "symmetric-indefinite", "qr" } :
            std::vector<std::string>{ test->getStringWithDefault("blas_lapack_subdomain_solver_type", "") };
    const auto hierarchy_data = setup_hierarchy<NDIM>(app);
    Pointer<PatchHierarchy<NDIM>> hierarchy = std::get<0>(hierarchy_data);
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
    VariableDatabase<NDIM>* variables = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = variables->getContext("shell_test");
    Pointer<SideVariable<NDIM, double>> u = new SideVariable<NDIM, double>("u");
    Pointer<CellVariable<NDIM, double>> p = new CellVariable<NDIM, double>("p");
    Pointer<SideVariable<NDIM, double>> f = new SideVariable<NDIM, double>("f");
    Pointer<CellVariable<NDIM, double>> h = new CellVariable<NDIM, double>("h");
    Pointer<SideVariable<NDIM, int>> u_dof = new SideVariable<NDIM, int>("u_dof");
    Pointer<CellVariable<NDIM, int>> p_dof = new CellVariable<NDIM, int>("p_dof");
    const int ui = variables->registerVariableAndContext(u, context, IntVector<NDIM>(1));
    const int pi = variables->registerVariableAndContext(p, context, IntVector<NDIM>(1));
    const int fi = variables->registerVariableAndContext(f, context, IntVector<NDIM>(1));
    const int hi = variables->registerVariableAndContext(h, context, IntVector<NDIM>(1));
    const int udi = variables->registerVariableAndContext(u_dof, context, IntVector<NDIM>(1));
    const int pdi = variables->registerVariableAndContext(p_dof, context, IntVector<NDIM>(1));
    for (int index : { ui, pi, fi, hi, udi, pdi })
    {
        level->allocatePatchData(index);
    }
    SAMRAIVectorReal<NDIM, double> x("x", hierarchy, 0, 0), b("b", hierarchy, 0, 0);
    x.addComponent(u, ui);
    x.addComponent(p, pi);
    b.addComponent(f, fi);
    b.addComponent(h, hi);
    x.setToScalar(0.0);
    b.setToScalar(0.0);
    const double wavenumber = 2.0 * std::acos(-1.0) / input->getInteger("N");
    for (PatchLevel<NDIM>::Iterator patch_number(level); patch_number; patch_number++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(patch_number());
        Pointer<SideData<NDIM, double>> force = patch->getPatchData(fi);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            const Box<NDIM> box = SideGeometry<NDIM>::toSideBox(patch->getBox(), axis);
            for (Box<NDIM>::Iterator it(box); it; it++)
            {
                (*force)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower)) =
                    boundary ? axis + 1.0 : std::sin(wavenumber * it()(axis)) + 0.25 * (axis + 1.0);
            }
        }
    }
    std::vector<int> dofs;
    StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(dofs, udi, pdi, level);
    Vec rhs = nullptr, expected = nullptr, actual = nullptr;
    int ierr = VecCreateMPI(PETSC_COMM_WORLD, dofs[IBTK_MPI::getRank()], PETSC_DETERMINE, &rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(rhs, &expected);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(rhs, &actual);
    IBTK_CHKERRQ(ierr);
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(rhs, fi, udi, hi, pdi, level);
    Pointer<MemoryDatabase> db = new MemoryDatabase("solver");
    db->putString("ksp_type", "preonly");
    db->putString("pc_type", test->getStringWithDefault("pc_type", "shell"));
    if (test->keyExists("shell_pc_type"))
    {
        db->putString("shell_pc_type", shell_type);
    }
    db->putBool("initial_guess_nonzero", false);
    db->putInteger("max_iterations", 1);
    int box_size[NDIM];
    std::fill_n(box_size, NDIM, boundary ? input->getInteger("N") : 4);
    db->putIntegerArray("subdomain_box_size", box_size, NDIM);
    int failures = hand_failures;
    plog << std::setprecision(12);
    for (const std::string& solver_type : solver_types)
    {
        if (!solver_type.empty())
        {
            db->putString("blas_lapack_subdomain_solver_type", solver_type);
        }
        if (test->keyExists("blas_lapack_subdomain_solver_rcond"))
        {
            db->putDouble("blas_lapack_subdomain_solver_rcond", test->getDouble("blas_lapack_subdomain_solver_rcond"));
        }
        StaggeredStokesPETScLevelSolver solver("shell_solver", db, "shell_");
        PoissonSpecifications coefficients("coefficients");
        coefficients.setCConstant(1.0);
        coefficients.setDConstant(-1.0);
        solver.setVelocityPoissonSpecifications(coefficients);
        solver.setComponentsHaveNullSpace(false, !diagonal_operator);
        std::vector<std::unique_ptr<LocationIndexRobinBcCoefs<NDIM>>> bc_storage;
        std::vector<RobinBcCoefStrategy<NDIM>*> bcs(NDIM, nullptr);
        if (boundary)
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                bc_storage.push_back(std::make_unique<LocationIndexRobinBcCoefs<NDIM>>("bc", nullptr));
                for (int face = 0; face < 2 * NDIM; ++face)
                {
                    bc_storage.back()->setBoundaryValue(face, axis + 1.0);
                }
                bcs[axis] = bc_storage.back().get();
            }
            Pointer<StaggeredStokesPhysicalBoundaryHelper> helper = new StaggeredStokesPhysicalBoundaryHelper();
            helper->cacheBcCoefData(bcs, 0.0, hierarchy);
            solver.setPhysicalBcCoefs(bcs, nullptr);
            solver.setPhysicalBoundaryHelper(helper);
            solver.setHomogeneousBc(false);
        }
        Mat supplied = nullptr;
        if (diagonal_operator)
        {
            // A diagonal operator has an analytic truncated pseudoinverse: entries
            // below the requested cutoff contribute zero, and the others divide by two.
            PetscInt n = 0;
            ierr = VecGetSize(rhs, &n);
            IBTK_CHKERRQ(ierr);
            ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, 1, nullptr, &supplied);
            IBTK_CHKERRQ(ierr);
            for (PetscInt j = 0; j < n; ++j)
            {
                ierr = MatSetValue(supplied, j, j, j % 2 == 0 ? 2.0 : small_diagonal, INSERT_VALUES);
                IBTK_CHKERRQ(ierr);
            }
            ierr = MatAssemblyBegin(supplied, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);
            ierr = MatAssemblyEnd(supplied, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);
            solver.setOperatorMat(supplied);
        }
        Vec default_forward = nullptr;
        if (test->getBoolWithDefault("compare_forward", false))
        {
            ierr = VecDuplicate(rhs, &default_forward);
            IBTK_CHKERRQ(ierr);
            solver.initializeSolverState(x, b);
            if (!solver.solveSystem(x, b))
            {
                ++failures;
            }
            StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(default_forward, ui, udi, pi, pdi, level);
            solver.deallocateSolverState();
        }
        if (test->keyExists("shell_pc_subdomain_traversal"))
        {
            db->putString("shell_pc_subdomain_traversal", test->getString("shell_pc_subdomain_traversal"));
        }
        solver.initializeSolverState(x, b);
        if (invalid)
        {
            solver.deallocateSolverState();
            return 0;
        }
        if (lifetime && !diagonal_operator)
        {
            Mat assembled = nullptr;
            ierr = KSPGetOperators(solver.getPETScKSP(), &assembled, nullptr);
            IBTK_CHKERRQ(ierr);
            ierr = MatDuplicate(assembled, MAT_COPY_VALUES, &supplied);
            IBTK_CHKERRQ(ierr);
            solver.deallocateSolverState();
            solver.setOperatorMat(supplied);
            solver.setOperatorMat(supplied);
            solver.initializeSolverState(x, b);
        }
        if (all_blas_modes)
        {
            plog << "solver_type = " << (solver_type.empty() ? "default" : solver_type) << '\n';
        }
        for (int cycle = 0; cycle < (lifetime ? 2 : 1); ++cycle)
        {
            Mat mat = nullptr;
            PC pc = nullptr;
            ierr = KSPGetOperators(solver.getPETScKSP(), &mat, nullptr);
            IBTK_CHKERRQ(ierr);
            ierr = KSPGetPC(solver.getPETScKSP(), &pc);
            IBTK_CHKERRQ(ierr);
            PCType pc_type = nullptr;
            ierr = PCGetType(pc, &pc_type);
            IBTK_CHKERRQ(ierr);
            if (std::string(pc_type) != "shell" || (lifetime && mat != supplied))
            {
                ++failures;
            }
            if (diagonal_operator)
            {
                const PetscScalar* rhs_values = nullptr;
                PetscScalar* expected_values = nullptr;
                PetscInt n = 0;
                ierr = VecGetSize(rhs, &n);
                IBTK_CHKERRQ(ierr);
                ierr = VecGetArrayRead(rhs, &rhs_values);
                IBTK_CHKERRQ(ierr);
                ierr = VecGetArray(expected, &expected_values);
                IBTK_CHKERRQ(ierr);
                for (PetscInt j = 0; j < n; ++j)
                {
                    expected_values[j] = j % 2 == 0 ? rhs_values[j] / (cycle == 0 ? 2.0 : 4.0) : 0.0;
                }
                ierr = VecRestoreArray(expected, &expected_values);
                IBTK_CHKERRQ(ierr);
                ierr = VecRestoreArrayRead(rhs, &rhs_values);
                IBTK_CHKERRQ(ierr);
            }
            else if (!boundary)
            {
                std::vector<IS>* overlap = nullptr;
                std::vector<IS>* partition = nullptr;
                solver.getASMSubdomains(&partition, &overlap);
                reference_action(mat, rhs, expected, *overlap, *partition, multiplicative, traversal);
                // Left-preconditioned PETSc KSP removes the operator nullspace after PCApply.
                MatNullSpace nullspace = nullptr;
                ierr = MatGetNullSpace(mat, &nullspace);
                IBTK_CHKERRQ(ierr);
                if (nullspace)
                {
                    ierr = MatNullSpaceRemove(nullspace, expected);
                    IBTK_CHKERRQ(ierr);
                }
            }
            x.setToScalar(0.0);
            if (!solver.solveSystem(x, b))
            {
                ++failures;
            }
            StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(actual, ui, udi, pi, pdi, level);
            const double action_norm = norm_inf(actual);
            if (default_forward)
            {
                ierr = VecAXPY(default_forward, -1.0, actual);
                IBTK_CHKERRQ(ierr);
                const double default_error = norm_inf(default_forward);
                if (!std::isfinite(default_error) || default_error > 1.0e-12)
                {
                    ++failures;
                }
                plog << "default_forward_error = " << default_error << '\n';
                ierr = VecCopy(actual, default_forward);
                IBTK_CHKERRQ(ierr);
            }
            if (boundary)
            {
                // Constant velocity (1,2), pressure zero solves -Laplace(u)+u+grad(p)=u.
                // Its Dirichlet values are prescribed on every face, independently of RHS packing.
                StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(expected, fi, udi, hi, pdi, level);
            }
            ierr = VecAXPY(actual, -1.0, expected);
            IBTK_CHKERRQ(ierr);
            const double error = norm_inf(actual);
            if (!std::isfinite(error) || error > 1.0e-9 || action_norm <= 0.0)
            {
                ++failures;
            }
            plog << "action_norm = " << action_norm << "\nerror = " << error << '\n';
            if (lifetime)
            {
                if (!solver.solveSystem(x, b))
                {
                    ++failures;
                }
                StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(actual, ui, udi, pi, pdi, level);
                ierr = VecAXPY(actual, -1.0, expected);
                IBTK_CHKERRQ(ierr);
                const double repeated_error = norm_inf(actual);
                if (!std::isfinite(repeated_error) || repeated_error > 1.0e-9)
                {
                    ++failures;
                }
                plog << "repeated_error = " << repeated_error << '\n';
            }
            solver.deallocateSolverState();
            if (lifetime)
            {
                PetscReal matrix_norm = 0.0;
                ierr = MatNorm(supplied, NORM_INFINITY, &matrix_norm);
                IBTK_CHKERRQ(ierr);
                if (!(matrix_norm > 0.0))
                {
                    ++failures;
                }
                if (cycle == 0)
                {
                    if (shell_type == "additive-blas-lapack")
                    {
                        ierr = MatScale(supplied, 2.0);
                        IBTK_CHKERRQ(ierr);
                    }
                    solver.initializeSolverState(x, b);
                }
            }
        }
        ierr = VecDestroy(&default_forward);
        IBTK_CHKERRQ(ierr);
        ierr = MatDestroy(&supplied);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecDestroy(&rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&expected);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&actual);
    IBTK_CHKERRQ(ierr);
    plog << "failures = " << failures << '\n';
    return failures ? 1 : 0;
}
