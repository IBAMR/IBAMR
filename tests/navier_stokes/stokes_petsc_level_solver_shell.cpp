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

#include <ibamr/StaggeredStokesEigenSchurComplementSubdomainSolver.h>
#include <ibamr/StaggeredStokesPETScLevelSolver.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>

#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PETScLevelSolverSubdomainSolver.h>

#include <tbox/MemoryDatabase.h>

#include <CellData.h>
#include <CellVariable.h>
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
#include <optional>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

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

// The action of the shell preconditioner, computed independently of the solver. It gathers the
// right-hand side on every rank and solves each local matrix. The additive shell writes the
// nonoverlapping part of each scaled solution. The multiplicative shell visits the subdomains stage by
// stage, with the subdomain i of each rank at stage i, computes the residual of the original system
// from the current result, and adds the whole scaled solution of each subdomain.
void
reference_action(Mat mat,
                 Vec rhs,
                 Vec result,
                 const std::vector<IS>& overlap,
                 const std::vector<IS>& partition,
                 const bool multiplicative,
                 const double scale,
                 const std::string& traversal)
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
    // The subdomains that this rank visits at each stage of the multiplicative shell.
    std::vector<std::size_t> visits;
    for (std::size_t i = 0; i < overlap.size(); ++i)
    {
        visits.push_back(traversal == "REVERSE" ? overlap.size() - 1 - i : i);
    }
    if (traversal == "SYMMETRIC")
    {
        for (std::size_t i = 1; i < overlap.size(); ++i)
        {
            visits.push_back(overlap.size() - 1 - i);
        }
    }
    const std::size_t n_stages = multiplicative ?
                                     static_cast<std::size_t>(IBTK_MPI::maxReduction(static_cast<int>(visits.size()))) :
                                     overlap.size();
    for (std::size_t stage = 0; stage < n_stages; ++stage)
    {
        const std::size_t i = multiplicative && stage < visits.size() ? visits[stage] : stage;
        if (multiplicative)
        {
            ierr = MatMult(mat, result, residual);
            IBTK_CHKERRQ(ierr);
            ierr = VecAYPX(residual, -1.0, rhs);
            IBTK_CHKERRQ(ierr);
            ierr = VecScatterBegin(gather, residual, gathered, INSERT_VALUES, SCATTER_FORWARD);
            IBTK_CHKERRQ(ierr);
            ierr = VecScatterEnd(gather, residual, gathered, INSERT_VALUES, SCATTER_FORWARD);
            IBTK_CHKERRQ(ierr);
        }
        if (multiplicative ? stage < visits.size() : i < overlap.size())
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
            ierr = VecScale(local_solution, scale);
            IBTK_CHKERRQ(ierr);
            const PetscScalar* solution_values = nullptr;
            ierr = VecGetArrayRead(local_solution, &solution_values);
            IBTK_CHKERRQ(ierr);
            if (multiplicative)
            {
                ierr = VecSetValues(result, n, indices, solution_values, ADD_VALUES);
                IBTK_CHKERRQ(ierr);
            }
            else
            {
                for (PetscInt j = 0; j < m; ++j)
                {
                    const PetscInt position =
                        static_cast<PetscInt>(std::lower_bound(indices, indices + n, owned[j]) - indices);
                    TBOX_ASSERT(position < n && indices[position] == owned[j]);
                    ierr = VecSetValue(result, owned[j], solution_values[position], INSERT_VALUES);
                    IBTK_CHKERRQ(ierr);
                }
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

// The factor by which the application-supplied subdomain solver scales the solution of the built-in one.
constexpr double SUBDOMAIN_SOLVER_SCALE = 0.5;

// A value no correction of the test problem approaches, used to detect entries a subdomain solver fails to overwrite.
constexpr double OUTPUT_SENTINEL = 1.0e30;

// Counts the use of an application-supplied subdomain solver. The test owns the counters, since the level solver
// owns the solver.
struct SubdomainSolverCounters
{
    int constructions = 0, initializations = 0, deallocations = 0, calls = 0, subdomain_solves = 0,
        contract_violations = 0;
};

// An application-supplied subdomain solver that scales the solution of the built-in one, counts its use and checks
// that the built-in one respects the input and output contract.
class ScaledCountingSubdomainSolver
{
public:
    explicit ScaledCountingSubdomainSolver(SubdomainSolverCounters& counters) : d_counters(counters)
    {
        ++d_counters.constructions;
    }
    ScaledCountingSubdomainSolver(const ScaledCountingSubdomainSolver&) = delete;
    ScaledCountingSubdomainSolver(ScaledCountingSubdomainSolver&&) = delete;
    ScaledCountingSubdomainSolver& operator=(const ScaledCountingSubdomainSolver&) = delete;
    ScaledCountingSubdomainSolver& operator=(ScaledCountingSubdomainSolver&&) = delete;
    void initializeSolverState(const std::vector<Mat>& matrices,
                               const std::vector<IS>& subdomains,
                               const std::string& options_prefix)
    {
        ++d_counters.initializations;
        d_counters.contract_violations += matrices.size() != subdomains.size();
        d_offsets.assign(matrices.size() + 1, 0);
        for (std::size_t i = 0; i < matrices.size(); ++i)
        {
            PetscInt order = 0;
            int ierr = MatGetSize(matrices[i], &order, nullptr);
            IBTK_CHKERRQ(ierr);
            PetscInt subdomain_size = 0;
            if (i < subdomains.size())
            {
                ierr = ISGetLocalSize(subdomains[i], &subdomain_size);
                IBTK_CHKERRQ(ierr);
            }
            d_counters.contract_violations += subdomain_size != order;
            d_offsets[i + 1] = d_offsets[i] + order;
        }
        d_subdomain_solver.initializeSolverState(matrices, subdomains, options_prefix);
    }
    void deallocateSolverState()
    {
        ++d_counters.deallocations;
        d_subdomain_solver.deallocateSolverState();
    }
    void solve(const std::size_t first, const std::size_t last, Vec b, Vec x)
    {
        ++d_counters.calls;
        d_counters.subdomain_solves += static_cast<int>(last - first);
        // The subdomain solver must leave its input unchanged, leave the solutions of other subdomains
        // unchanged, and overwrite every entry of the solutions that it computes.
        Vec b_before = nullptr, x_before = nullptr;
        int ierr = VecDuplicate(b, &b_before);
        IBTK_CHKERRQ(ierr);
        ierr = VecCopy(b, b_before);
        IBTK_CHKERRQ(ierr);
        ierr = VecDuplicate(x, &x_before);
        IBTK_CHKERRQ(ierr);
        ierr = VecCopy(x, x_before);
        IBTK_CHKERRQ(ierr);
        const PetscInt begin = d_offsets[first], end = d_offsets[last];
        PetscScalar* x_values = nullptr;
        ierr = VecGetArray(x, &x_values);
        IBTK_CHKERRQ(ierr);
        std::fill(x_values + begin, x_values + end, OUTPUT_SENTINEL);
        ierr = VecRestoreArray(x, &x_values);
        IBTK_CHKERRQ(ierr);
        d_subdomain_solver.solve(first, last, b, x);
        ierr = VecAXPY(b_before, -1.0, b);
        IBTK_CHKERRQ(ierr);
        PetscReal b_change = 0.0;
        ierr = VecNorm(b_before, NORM_INFINITY, &b_change);
        IBTK_CHKERRQ(ierr);
        const PetscScalar *x_now = nullptr, *x_old = nullptr;
        ierr = VecGetArrayRead(x, &x_now);
        IBTK_CHKERRQ(ierr);
        ierr = VecGetArrayRead(x_before, &x_old);
        IBTK_CHKERRQ(ierr);
        PetscInt n = 0;
        ierr = VecGetSize(x, &n);
        IBTK_CHKERRQ(ierr);
        bool violated = b_change != 0.0;
        for (PetscInt k = 0; k < n; ++k)
        {
            violated = violated || (begin <= k && k < end ? !(PetscAbsScalar(x_now[k]) < 0.5 * OUTPUT_SENTINEL) :
                                                            x_now[k] != x_old[k]);
        }
        d_counters.contract_violations += violated;
        ierr = VecRestoreArrayRead(x_before, &x_old);
        IBTK_CHKERRQ(ierr);
        ierr = VecRestoreArrayRead(x, &x_now);
        IBTK_CHKERRQ(ierr);
        ierr = VecGetArray(x, &x_values);
        IBTK_CHKERRQ(ierr);
        for (PetscInt k = begin; k < end; ++k)
        {
            x_values[k] *= SUBDOMAIN_SOLVER_SCALE;
        }
        ierr = VecRestoreArray(x, &x_values);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&b_before);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&x_before);
        IBTK_CHKERRQ(ierr);
    }

private:
    SubdomainSolverCounters& d_counters;
    PETScLevelSolverSubdomainSolver d_subdomain_solver = make_petsc_subdomain_solver();
    std::vector<PetscInt> d_offsets;
};

struct MissingSolve
{
    void initializeSolverState(const std::vector<Mat>&, const std::vector<IS>&, const std::string&);
    void deallocateSolverState();
};

// An implementation needs neither a base class nor copy, move, or default construction.
static_assert(PETScLevelSolverSubdomainSolverImplementation<ScaledCountingSubdomainSolver>);
static_assert(!std::is_copy_constructible_v<ScaledCountingSubdomainSolver>);
static_assert(!std::is_move_constructible_v<ScaledCountingSubdomainSolver>);
static_assert(!PETScLevelSolverSubdomainSolverImplementation<MissingSolve>);
static_assert(!std::is_constructible_v<PETScLevelSolverSubdomainSolver, std::in_place_type_t<MissingSolve>>);
static_assert(
    !std::is_constructible_v<PETScLevelSolverSubdomainSolver, std::in_place_type_t<ScaledCountingSubdomainSolver>>);
static_assert(!std::is_copy_constructible_v<PETScLevelSolverSubdomainSolver>);

// Exposes whether the solver communicates.
class CommunicationProbe : public StaggeredStokesPETScLevelSolver
{
public:
    using StaggeredStokesPETScLevelSolver::StaggeredStokesPETScLevelSolver;
    bool restrictionCommunicates() const
    {
        return d_restriction_communicates;
    }
    int numberOfStages() const
    {
        return d_n_stages;
    }
    int haloCommunicatingStages() const
    {
        return static_cast<int>(std::count(d_halo_communicates.begin(), d_halo_communicates.end(), true));
    }
    int correctionCommunicatingStages() const
    {
        return static_cast<int>(std::count(d_correction_communicates.begin(), d_correction_communicates.end(), true));
    }
    // Remove a DOF from the nonoverlapping set of the first subdomain, so that the sets no longer partition the DOFs.
    bool break_partition = false;

protected:
    void generateASMSubdomains(std::vector<std::set<int>>& overlap_is,
                               std::vector<std::set<int>>& nonoverlap_is) override
    {
        StaggeredStokesPETScLevelSolver::generateASMSubdomains(overlap_is, nonoverlap_is);
        if (break_partition && !nonoverlap_is.empty() && !nonoverlap_is.front().empty())
        {
            nonoverlap_is.front().erase(nonoverlap_is.front().begin());
        }
    }
};

// Solve small dense systems directly with the Eigen subdomain solvers, whose solutions are known.
int
check_eigen_local(Pointer<Database> test)
{
    const std::string mode = test->getString("eigen_local");
    const bool spd = mode == "all";
    const bool rank_test = mode == "rank";
    // A coupled (non-diagonal) rank-deficient matrix: unlike a diagonal matrix, truncating its
    // rank actually eliminates a nonzero off-diagonal coupling term, so this case distinguishes a
    // solver that applies the configured threshold before factoring from one that applies it after.
    const bool rank_coupled_test = mode == "rank_coupled";
    const std::vector<std::string> types =
        mode == "all"     ? std::vector<std::string>{ "LLT",
                                                      "LDLT",
                                                      "PARTIAL_PIV_LU",
                                                      "FULL_PIV_LU",
                                                      "HOUSEHOLDER_QR",
                                                      "COL_PIV_HOUSEHOLDER_QR",
                                                      "COMPLETE_ORTHOGONAL_DECOMPOSITION",
                                                      "FULL_PIV_HOUSEHOLDER_QR",
                                                      "JACOBI_SVD",
                                                      "BDC_SVD" } :
        rank_coupled_test ? std::vector<std::string>{ "COMPLETE_ORTHOGONAL_DECOMPOSITION" } :
        rank_test         ? std::vector<std::string>{ "COMPLETE_ORTHOGONAL_DECOMPOSITION", "JACOBI_SVD", "BDC_SVD" } :
                    std::vector<std::string>{ "FULL_PIV_HOUSEHOLDER_QR", "COL_PIV_HOUSEHOLDER_QR", "PARTIAL_PIV_LU" };
    const double spd_entries[3][3] = { { 4.0, 1.0, 1.0 }, { 1.0, 3.0, -1.0 }, { 1.0, -1.0, 5.0 } };
    const double nonsymmetric_entries[3][3] = { { 4.0, 2.0, 1.0 }, { 1.0, 3.0, -1.0 }, { 0.0, -1.0, 5.0 } };
    Mat mat = nullptr;
    int ierr = MatCreateSeqDense(PETSC_COMM_SELF, 3, 3, nullptr, &mat);
    IBTK_CHKERRQ(ierr);
    const double target[3] = { 1.0, -2.0, 3.0 };
    double rhs_values[3] = { 0.0, 0.0, 0.0 }, expected[3] = { 0.0, 0.0, 0.0 };
    for (PetscInt i = 0; i < 3; ++i)
    {
        for (PetscInt j = 0; j < 3; ++j)
        {
            // With a rank threshold of 0.1, the second and third pivots of the diagonal matrix are dropped.
            // rank_coupled_test additionally couples the first two rows, so truncating to rank one also
            // eliminates that coupling term rather than leaving it untouched.
            const double value = rank_coupled_test ? (i == 0 && j == 1 ? 1.0 :
                                                      i == j           ? (i == 0 ? 2.0 :
                                                                          i == 1 ? 0.01 :
                                                                                   0.0) :
                                                                         0.0) :
                                 rank_test         ? (i == j ? (i == 0 ? 2.0 :
                                                                i == 1 ? 0.01 :
                                                                         0.0) :
                                                               0.0) :
                                                     (spd ? spd_entries[i][j] : nonsymmetric_entries[i][j]);
            ierr = MatSetValue(mat, i, j, value, INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
            rhs_values[i] += (rank_test || rank_coupled_test) ? 0.0 : value * target[j];
        }
        expected[i] = (rank_test || rank_coupled_test) ? 0.0 : target[i];
    }
    if (rank_test)
    {
        std::fill_n(rhs_values, 3, 1.0);
        rhs_values[0] = 2.0;
        expected[0] = 1.0;
    }
    if (rank_coupled_test)
    {
        // The unique rank-one pseudoinverse action: truncating the trailing 0.01 pivot also removes
        // its contribution to the coupled first row, which a solver reusing rank-two factors would not do.
        rhs_values[0] = 1.0;
        expected[0] = 0.4;
        expected[1] = 0.2;
    }
    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    // The subdomain has global DOFs 10, 11, 12. The Schur cases use the velocity and pressure DOFs
    // {10, 12} and {11}, {} and {10, 11, 12}, or {10, 11, 12} and {}, and an invalid case that omits DOF 11.
    const std::set<int> velocity = mode == "schur_pressure_only" ? std::set<int>{} :
                                   mode == "schur_velocity_only" ? std::set<int>{ 10, 11, 12 } :
                                                                   std::set<int>{ 10, 12 };
    const std::set<int> pressure = mode == "schur_pressure_only"  ? std::set<int>{ 10, 11, 12 } :
                                   mode == "schur_velocity_only"  ? std::set<int>{} :
                                   mode == "schur_invalid_fields" ? std::set<int>{} :
                                                                    std::set<int>{ 11 };
    IS subdomain = nullptr;
    const PetscInt dofs[3] = { 10, 11, 12 };
    ierr = ISCreateGeneral(PETSC_COMM_SELF, 3, dofs, PETSC_COPY_VALUES, &subdomain);
    IBTK_CHKERRQ(ierr);
    Vec b = nullptr, x = nullptr;
    ierr = VecCreateSeq(PETSC_COMM_SELF, 3, &b);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(b, &x);
    IBTK_CHKERRQ(ierr);
    for (PetscInt i = 0; i < 3; ++i)
    {
        ierr = VecSetValue(b, i, rhs_values[i], INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
    }
    double max_error = 0.0;
    const bool schur = mode.rfind("schur", 0) == 0;
    for (const std::string& type : schur ? std::vector<std::string>{ "" } : types)
    {
        for (const std::string& name :
             schur ? std::vector<std::string>{ "schur" } : std::vector<std::string>{ "eigen", "eigen-pseudoinverse" })
        {
            Pointer<MemoryDatabase> db = new MemoryDatabase("eigen");
            db->putString("eigen_subdomain_solver_type", type);
            db->putString("eigen_subdomain_pseudoinverse_type", type);
            db->putDouble("eigen_subdomain_solver_threshold", (rank_test || rank_coupled_test) ? 0.1 : -1.0);
            db->putDouble("eigen_subdomain_pseudoinverse_threshold", (rank_test || rank_coupled_test) ? 0.1 : -1.0);
            std::optional<PETScLevelSolverSubdomainSolver> subdomain_solver;
            if (schur)
            {
                subdomain_solver.emplace(std::in_place_type<StaggeredStokesEigenSchurComplementSubdomainSolver>,
                                         db,
                                         [&]()
                                         {
                                             Vec indicators = nullptr;
                                             int provider_ierr = VecCreateSeq(PETSC_COMM_SELF, 13, &indicators);
                                             IBTK_CHKERRQ(provider_ierr);
                                             provider_ierr = VecSet(indicators, -1.0);
                                             IBTK_CHKERRQ(provider_ierr);
                                             for (const int dof : velocity)
                                             {
                                                 provider_ierr = VecSetValue(indicators, dof, 0.0, INSERT_VALUES);
                                                 IBTK_CHKERRQ(provider_ierr);
                                             }
                                             for (const int dof : pressure)
                                             {
                                                 provider_ierr = VecSetValue(indicators, dof, 1.0, INSERT_VALUES);
                                                 IBTK_CHKERRQ(provider_ierr);
                                             }
                                             return indicators;
                                         });
            }
            else
            {
                subdomain_solver.emplace(name == "eigen" ? make_eigen_subdomain_solver(db) :
                                                           make_eigen_pseudoinverse_subdomain_solver(db));
            }
            // Reinitialization reuses the subdomain solver.
            for (int cycle = 0; cycle < 2; ++cycle)
            {
                subdomain_solver->initializeSolverState({ mat }, { subdomain }, "eigen_");
                ierr = VecSet(x, OUTPUT_SENTINEL);
                IBTK_CHKERRQ(ierr);
                subdomain_solver->solve(0, 1, b, x);
                subdomain_solver->deallocateSolverState();
                const PetscScalar* values = nullptr;
                ierr = VecGetArrayRead(x, &values);
                IBTK_CHKERRQ(ierr);
                for (PetscInt i = 0; i < 3; ++i)
                {
                    // std::max would discard a NaN.
                    const double entry_error = std::abs(values[i] - expected[i]);
                    max_error = std::isfinite(entry_error) ? std::max(max_error, entry_error) :
                                                             std::numeric_limits<double>::infinity();
                }
                ierr = VecRestoreArrayRead(x, &values);
                IBTK_CHKERRQ(ierr);
            }
        }
    }
    if (!(max_error < 1.0e-10))
    {
        TBOX_ERROR("Failed check: max_error < 1.0e-10 (max_error = " << max_error << ").\n");
    }
    plog << "max_error = " << max_error << '\n';
    ierr = MatDestroy(&mat);
    IBTK_CHKERRQ(ierr);
    ierr = ISDestroy(&subdomain);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&b);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    IBTK_CHKERRQ(ierr);
    return 0;
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
    if (test->keyExists("eigen_local"))
    {
        return check_eigen_local(test);
    }
    const bool lifetime = test->getBoolWithDefault("lifetime", false);
    const bool application_subdomain_solver = test->getBoolWithDefault("application_subdomain_solver", false);
    const bool null_subdomain_solver = test->getBoolWithDefault("null_subdomain_solver", false);
    const bool subdomain_solver_unused = test->getBoolWithDefault("subdomain_solver_unused", false);
    const bool replace_initialized = test->getBoolWithDefault("replace_initialized", false);
    const bool uneven_subdomains = test->getBoolWithDefault("uneven_subdomains", false);
    const std::string shell_type = test->getString("shell_pc_type");
    std::string traversal = test->getStringWithDefault("shell_pc_subdomain_traversal", "FORWARD");
    const std::string subdomain_solver_type = test->getStringWithDefault("subdomain_solver", "petsc");
    const bool multiplicative = shell_type == "multiplicative";
    const bool diagonal_operator = test->getBoolWithDefault("diagonal_operator", false);
    const double small_diagonal = test->getDoubleWithDefault("small_diagonal", 1.0e-12);
    // "rank_one", "upper" or "symmetric" supplies a dense operator times operator_scale.
    const std::string dense_operator = test->getStringWithDefault("dense_operator", "");
    const double operator_scale = test->getDoubleWithDefault("operator_scale", 1.0);
    const bool supplied_operator = diagonal_operator || !dense_operator.empty();
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
                    std::sin(wavenumber * it()(axis)) + 0.25 * (axis + 1.0);
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
    db->putString("shell_pc_type", shell_type);
    db->putString("shell_pc_subdomain_traversal", traversal);
    std::transform(
        traversal.begin(), traversal.end(), traversal.begin(), [](const unsigned char c) { return std::toupper(c); });
    db->putString("subdomain_solver", subdomain_solver_type);
    db->putBool("initial_guess_nonzero", false);
    db->putBool("check_subdomain_coverage", true);
    db->putInteger("max_iterations", 1);
    int box_size[NDIM];
    std::fill_n(box_size, NDIM, 4);
    db->putIntegerArray("subdomain_box_size", box_size, NDIM);
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
        for (const std::string key : { "eigen_subdomain_solver_type", "eigen_subdomain_pseudoinverse_type" })
        {
            if (test->keyExists(key))
            {
                db->putString(key, test->getString(key));
            }
        }
        CommunicationProbe solver("shell_solver", db, "shell_");
        PoissonSpecifications coefficients("coefficients");
        coefficients.setCConstant(1.0);
        coefficients.setDConstant(-1.0);
        solver.break_partition = test->getBoolWithDefault("break_partition", false);
        solver.setVelocityPoissonSpecifications(coefficients);
        solver.setComponentsHaveNullSpace(false, !supplied_operator);
        SubdomainSolverCounters application_counters;
        SubdomainSolverCounters* counts = nullptr;
        if (application_subdomain_solver)
        {
            PETScLevelSolverSubdomainSolver subdomain_solver(std::in_place_type<ScaledCountingSubdomainSolver>,
                                                             application_counters);
            counts = &application_counters;
            // Moving a handle transfers the implementation, which is constructed once and never moved.
            PETScLevelSolverSubdomainSolver installed(std::move(subdomain_solver));
            if (subdomain_solver || !installed || application_counters.constructions != 1)
            {
                TBOX_ERROR("Failed check: moving a subdomain solver did not transfer ownership.\n");
            }
            solver.setSubdomainSolver(std::move(installed));
        }
        if (null_subdomain_solver)
        {
            // A handle that has been moved from is empty.
            PETScLevelSolverSubdomainSolver spent = make_petsc_subdomain_solver();
            PETScLevelSolverSubdomainSolver taken(std::move(spent));
            solver.setSubdomainSolver(std::move(spent));
            return 0;
        }
        Mat supplied = nullptr;
        if (diagonal_operator)
        {
            // The operator is diagonal, with 2 on even rows and small_diagonal on odd rows. For the SVD
            // solver, the expected action is that of its truncated pseudoinverse: entries below the
            // requested cutoff contribute zero, and the others divide by two.
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
        else if (!dense_operator.empty())
        {
            // Every principal submatrix of "rank_one" has rank one, of "upper" is nonsymmetric with
            // relative asymmetry of order one, and of "symmetric" is symmetric positive definite.
            PetscInt n = 0;
            ierr = VecGetSize(rhs, &n);
            IBTK_CHKERRQ(ierr);
            ierr = MatCreateSeqDense(PETSC_COMM_SELF, n, n, nullptr, &supplied);
            IBTK_CHKERRQ(ierr);
            for (PetscInt i = 0; i < n; ++i)
            {
                for (PetscInt j = 0; j < n; ++j)
                {
                    double value = 0.0;
                    if (dense_operator == "rank_one")
                    {
                        value = (1.0 + i % 3) * (1.0 + j % 3);
                    }
                    else if (dense_operator == "upper")
                    {
                        value = i <= j ? 1.0 : 0.0;
                    }
                    else
                    {
                        value = i == j ? 2.0 : 1.0;
                    }
                    ierr = MatSetValue(supplied, i, j, operator_scale * value, INSERT_VALUES);
                    IBTK_CHKERRQ(ierr);
                }
            }
            ierr = MatAssemblyBegin(supplied, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);
            ierr = MatAssemblyEnd(supplied, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);
            solver.setOperatorMat(supplied);
        }
        solver.initializeSolverState(x, b);
        if (test->getBoolWithDefault("report_communication", false))
        {
            plog << "restriction communicates = " << solver.restrictionCommunicates() << '\n';
            if (multiplicative)
            {
                plog << "stages = " << solver.numberOfStages()
                     << "\nhalo communicating stages = " << solver.haloCommunicatingStages()
                     << "\ncorrection communicating stages = " << solver.correctionCommunicatingStages() << '\n';
            }
        }
        if (replace_initialized)
        {
            SubdomainSolverCounters replacement_counters;
            solver.setSubdomainSolver(PETScLevelSolverSubdomainSolver(std::in_place_type<ScaledCountingSubdomainSolver>,
                                                                      replacement_counters));
            return 0;
        }
        if (subdomain_solver_unused)
        {
            // The subdomain solver is initialized only for a shell preconditioner.
            x.setToScalar(0.0);
            if (!solver.solveSystem(x, b))
            {
                TBOX_ERROR("Failed check: solver.solveSystem(x, b).\n");
            }
            solver.deallocateSolverState();
            plog << "subdomain solver initializations = " << counts->initializations
                 << "\nsubdomain solver deallocations = " << counts->deallocations
                 << "\nsubdomain solver calls = " << counts->calls << '\n';
            if (counts->initializations != 0 || counts->deallocations != 0 || counts->calls != 0)
            {
                TBOX_ERROR("The subdomain solver was used although the preconditioner is not a shell.\n");
            }
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
                TBOX_ERROR("Failed check: std::string(pc_type) != 'shell' || (lifetime && mat != supplied).\n");
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
            else
            {
                std::vector<IS>* overlap = nullptr;
                std::vector<IS>* partition = nullptr;
                solver.getASMSubdomains(&partition, &overlap);
                if (IBTK_MPI::getNodes() > 1)
                {
                    // Ranks with different numbers of subdomains must still all take part in the collective scatters.
                    const int local_count = static_cast<int>(overlap->size());
                    const int min_count = IBTK_MPI::minReduction(local_count);
                    const int max_count = IBTK_MPI::maxReduction(local_count);
                    plog << "local subdomains: min = " << min_count << ", max = " << max_count << '\n';
                    if (uneven_subdomains && min_count == max_count)
                    {
                        TBOX_ERROR("Failed check: uneven_subdomains && min_count == max_count.\n");
                    }
                }
                // The supplied subdomain solver, not the built-in one, determines the action.
                reference_action(mat,
                                 rhs,
                                 expected,
                                 *overlap,
                                 *partition,
                                 multiplicative,
                                 counts ? SUBDOMAIN_SOLVER_SCALE : 1.0,
                                 traversal);
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
                TBOX_ERROR("Failed check: !solver.solveSystem(x, b).\n");
            }
            StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(actual, ui, udi, pi, pdi, level);
            const double action_norm = norm_inf(actual);
            ierr = VecAXPY(actual, -1.0, expected);
            IBTK_CHKERRQ(ierr);
            const double error = norm_inf(actual);
            if (!std::isfinite(error) || error > 1.0e-9 * action_norm || action_norm <= 0.0)
            {
                TBOX_ERROR(
                    "Failed check: !std::isfinite(error) || error > 1.0e-9 * action_norm || action_norm <= 0.0.\n");
            }
            plog << "action_norm = " << action_norm << "\nerror = " << error / action_norm << '\n';
            if (lifetime)
            {
                if (!solver.solveSystem(x, b))
                {
                    TBOX_ERROR("Failed check: !solver.solveSystem(x, b).\n");
                }
                StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(actual, ui, udi, pi, pdi, level);
                ierr = VecAXPY(actual, -1.0, expected);
                IBTK_CHKERRQ(ierr);
                const double repeated_error = norm_inf(actual);
                if (!std::isfinite(repeated_error) || repeated_error > 1.0e-9)
                {
                    TBOX_ERROR("Failed check: !std::isfinite(repeated_error) || repeated_error > 1.0e-9.\n");
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
                    TBOX_ERROR("Failed check: !(matrix_norm > 0.0).\n");
                }
                if (cycle == 0)
                {
                    if (subdomain_solver_type == "blas-lapack")
                    {
                        ierr = MatScale(supplied, 2.0);
                        IBTK_CHKERRQ(ierr);
                    }
                    solver.initializeSolverState(x, b);
                }
            }
        }
        if (counts)
        {
            // The subdomain solver is initialized and released once per solver initialization.
            plog << "subdomain solver initializations = " << counts->initializations
                 << "\nsubdomain solver deallocations = " << counts->deallocations
                 << "\nsubdomain solver calls = " << counts->calls
                 << "\nsubdomain solver subdomain solves = " << counts->subdomain_solves
                 << "\nsubdomain solver contract violations = " << counts->contract_violations << '\n';
            if (counts->contract_violations != 0)
            {
                TBOX_ERROR("The subdomain solver contract was violated " << counts->contract_violations << " times.\n");
            }
            if (counts->constructions != 1 || counts->calls < 1 || counts->initializations != counts->deallocations)
            {
                TBOX_ERROR(
                    "Failed check: the subdomain solver was not constructed once, called, and released as often as "
                    "initialized.\n");
            }
        }
        ierr = MatDestroy(&supplied);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecDestroy(&rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&expected);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&actual);
    IBTK_CHKERRQ(ierr);
    return 0;
}
