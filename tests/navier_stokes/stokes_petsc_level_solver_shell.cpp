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
#include <ibamr/StaggeredStokesPETScMatUtilities.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>

#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScLevelSolverSubdomainSolver.h>
#include <ibtk/ibtk_utilities.h>

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
#include <array>
#include <cmath>
#include <deque>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

#include "../tests.h"
#include "coupling_aware_asm_test_utilities.h"

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

// The string setting key of db, or value if db is null or lacks the key. Unlike getStringWithDefault(), this does not
// add the key to db, which the solver reads.
std::string
read_setting(Pointer<Database> db, const std::string& key, const std::string& value)
{
    return db && db->keyExists(key) ? db->getString(key) : value;
}

// The action of subdomain relaxation, computed independently of the solver from the level operator and direct
// solves of the subdomain matrices. Every rank gathers the whole input r. ADDITIVE composition solves each
// subdomain of this rank with r(P_i). MULTIPLICATIVE composition, with one group for each rank, keeps the
// correction z of this rank's group as a vector of all of the DOFs and visits the subdomains in the order of the
// traversal, solving each with r(P_i) - A(P_i, :) z before adding its scaled solution to z. FULL output adds every
// entry of the scaled solutions (ADDITIVE) or of z on the subdomains of this rank (MULTIPLICATIVE) to the result, and
// OWNED output only the entries of the nonoverlapping sets.
void
reference_action(Mat mat,
                 Vec rhs,
                 Vec result,
                 const std::vector<IS>& overlap,
                 const std::vector<IS>& partition,
                 const bool multiplicative,
                 const bool owned_output,
                 const double scale,
                 const std::string& traversal)
{
    Vec gathered = nullptr, correction = nullptr;
    VecScatter gather = nullptr;
    int ierr = VecScatterCreateToAll(rhs, &gather, &gathered);
    IBTK_CHKERRQ(ierr);
    ierr = VecScatterBegin(gather, rhs, gathered, INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);
    ierr = VecScatterEnd(gather, rhs, gathered, INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(gathered, &correction);
    IBTK_CHKERRQ(ierr);
    ierr = VecSet(correction, 0.0);
    IBTK_CHKERRQ(ierr);
    ierr = VecSet(result, 0.0);
    IBTK_CHKERRQ(ierr);
    const PetscInt n_subdomains = static_cast<PetscInt>(overlap.size());
    PetscInt n_dofs = 0;
    ierr = VecGetSize(rhs, &n_dofs);
    IBTK_CHKERRQ(ierr);
    IS all_columns = nullptr;
    ierr = ISCreateStride(PETSC_COMM_SELF, n_dofs, 0, 1, &all_columns);
    IBTK_CHKERRQ(ierr);
    const std::vector<IS> all_columns_of_subdomains(overlap.size(), all_columns);
    Mat *submat = nullptr, *rows = nullptr;
    ierr = MatCreateSubMatrices(mat, n_subdomains, overlap.data(), overlap.data(), MAT_INITIAL_MATRIX, &submat);
    IBTK_CHKERRQ(ierr);
    ierr = MatCreateSubMatrices(
        mat, n_subdomains, overlap.data(), all_columns_of_subdomains.data(), MAT_INITIAL_MATRIX, &rows);
    IBTK_CHKERRQ(ierr);
    std::vector<PetscInt> visits;
    for (PetscInt i = 0; i < n_subdomains; ++i)
    {
        visits.push_back(traversal == "REVERSE" && multiplicative ? n_subdomains - 1 - i : i);
    }
    if (traversal == "SYMMETRIC" && multiplicative)
    {
        for (PetscInt i = n_subdomains - 2; i >= 0; --i)
        {
            visits.push_back(i);
        }
    }
    std::set<PetscInt> solved_dofs, owned_dofs;
    for (const PetscInt i : visits)
    {
        PetscInt n = 0;
        const PetscInt* indices = nullptr;
        ierr = ISGetLocalSize(overlap[i], &n);
        IBTK_CHKERRQ(ierr);
        ierr = ISGetIndices(overlap[i], &indices);
        IBTK_CHKERRQ(ierr);
        Vec local_rhs = nullptr, local_solution = nullptr;
        ierr = MatCreateVecs(submat[i], &local_solution, &local_rhs);
        IBTK_CHKERRQ(ierr);
        if (multiplicative)
        {
            // The residual that the earlier solves of the group leave.
            ierr = MatMult(rows[i], correction, local_rhs);
            IBTK_CHKERRQ(ierr);
            ierr = VecScale(local_rhs, -1.0);
            IBTK_CHKERRQ(ierr);
        }
        else
        {
            ierr = VecSet(local_rhs, 0.0);
            IBTK_CHKERRQ(ierr);
        }
        const PetscScalar* global_values = nullptr;
        PetscScalar* local_values = nullptr;
        ierr = VecGetArrayRead(gathered, &global_values);
        IBTK_CHKERRQ(ierr);
        ierr = VecGetArray(local_rhs, &local_values);
        IBTK_CHKERRQ(ierr);
        for (PetscInt j = 0; j < n; ++j)
        {
            local_values[j] += global_values[indices[j]];
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
            ierr = VecSetValues(correction, n, indices, solution_values, ADD_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyBegin(correction);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyEnd(correction);
            IBTK_CHKERRQ(ierr);
            solved_dofs.insert(indices, indices + n);
        }
        else if (!owned_output)
        {
            ierr = VecSetValues(result, n, indices, solution_values, ADD_VALUES);
            IBTK_CHKERRQ(ierr);
        }
        if (owned_output)
        {
            PetscInt m = 0;
            const PetscInt* owned = nullptr;
            ierr = ISGetLocalSize(partition[i], &m);
            IBTK_CHKERRQ(ierr);
            ierr = ISGetIndices(partition[i], &owned);
            IBTK_CHKERRQ(ierr);
            owned_dofs.insert(owned, owned + m);
            if (!multiplicative)
            {
                for (PetscInt j = 0; j < m; ++j)
                {
                    const PetscInt position =
                        static_cast<PetscInt>(std::lower_bound(indices, indices + n, owned[j]) - indices);
                    TBOX_ASSERT(position < n && indices[position] == owned[j]);
                    ierr = VecSetValue(result, owned[j], solution_values[position], ADD_VALUES);
                    IBTK_CHKERRQ(ierr);
                }
            }
            ierr = ISRestoreIndices(partition[i], &owned);
            IBTK_CHKERRQ(ierr);
        }
        ierr = VecRestoreArrayRead(local_solution, &solution_values);
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
        // Only the completed correction of the group is added to the result.
        const std::set<PetscInt>& output_dofs = owned_output ? owned_dofs : solved_dofs;
        const PetscScalar* correction_values = nullptr;
        ierr = VecGetArrayRead(correction, &correction_values);
        IBTK_CHKERRQ(ierr);
        for (const PetscInt dof : output_dofs)
        {
            ierr = VecSetValue(result, dof, correction_values[dof], ADD_VALUES);
            IBTK_CHKERRQ(ierr);
        }
        ierr = VecRestoreArrayRead(correction, &correction_values);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(result);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(result);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroyMatrices(n_subdomains, &submat);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroyMatrices(n_subdomains, &rows);
    IBTK_CHKERRQ(ierr);
    ierr = ISDestroy(&all_columns);
    IBTK_CHKERRQ(ierr);
    ierr = VecScatterDestroy(&gather);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&gathered);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&correction);
    IBTK_CHKERRQ(ierr);
}

// The factor by which the application-supplied subdomain solver scales the solution of the built-in one.
constexpr double SUBDOMAIN_SOLVER_SCALE = 0.5;

// A value no correction of the test problem approaches, used to detect entries that a subdomain solver or the
// preconditioner fails to overwrite.
constexpr double OUTPUT_SENTINEL = 1.0e30;

// Counts the use of an application-supplied subdomain solver. The test owns the counters, since the level solver
// owns the solver.
struct SubdomainSolverCounters
{
    int constructions = 0, initializations = 0, deallocations = 0, calls = 0, subdomain_solves = 0,
        contract_violations = 0;
    bool operator==(const SubdomainSolverCounters&) const = default;
};

// An application-supplied subdomain solver that scales the solution of the built-in one, counts its use and checks
// that the built-in one respects the input and output contract.
class ScaledCountingSubdomainSolver
{
public:
    explicit ScaledCountingSubdomainSolver(SubdomainSolverCounters& counters,
                                           const double scale = SUBDOMAIN_SOLVER_SCALE)
        : d_counters(counters), d_scale(scale)
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
            x_values[k] *= d_scale;
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
    const double d_scale;
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
            // The subdomain_solver database of each factory holds only its own settings.
            Pointer<MemoryDatabase> db = new MemoryDatabase("subdomain_solver");
            if (!schur)
            {
                const std::string key = name == "eigen" ? "eigen_subdomain_solver" : "eigen_subdomain_pseudoinverse";
                db->putString(key + "_type", type);
                db->putDouble(key + "_threshold", (rank_test || rank_coupled_test) ? 0.1 : -1.0);
            }
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
// The application fixture has one distant face edge initially and a complete
// remote velocity stencil after reinitialization. Enumerate logical cell supports
// independently of the production geometry owner and its global-index closure.
std::vector<std::set<int>>
cav_application_patches(const CAFields& fields, const int n, const bool strict, const int cycle, Mat elasticity)
{
    const CACell origin{};
    CACell remote{};
    remote.fill(2);
    const int source = fields.at(origin)[0];
    std::set<int> targets{ fields.at(remote)[0] };
    if (cycle == 1)
    {
        targets = ca_cell_stencil(fields, remote, n);
        targets.erase(fields.at(remote)[NDIM]);
    }
    int ierr = MatZeroEntries(elasticity);
    IBTK_CHKERRQ(ierr);
    for (const int target : targets)
    {
        // A symmetric positive semidefinite spring couples the two velocities.
        ierr = MatSetValue(elasticity, source, source, 0.25, ADD_VALUES);
        IBTK_CHKERRQ(ierr);
        ierr = MatSetValue(elasticity, target, target, 0.25, ADD_VALUES);
        IBTK_CHKERRQ(ierr);
        ierr = MatSetValue(elasticity, source, target, -0.25, ADD_VALUES);
        IBTK_CHKERRQ(ierr);
        ierr = MatSetValue(elasticity, target, source, -0.25, ADD_VALUES);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(elasticity, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(elasticity, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    std::vector<std::set<int>> patches;
    for (const auto& seed : fields)
    {
        std::set<int> patch = ca_cell_stencil(fields, seed.first, n);
        std::set<int> velocities = patch;
        velocities.erase(seed.second[NDIM]);
        const std::set<int> original = velocities;
        if (original.count(source))
        {
            velocities.insert(targets.begin(), targets.end());
        }
        if (std::any_of(targets.begin(), targets.end(), [&](int dof) { return original.count(dof); }))
        {
            velocities.insert(source);
        }
        if (velocities != original)
        {
            for (const auto& candidate : fields)
            {
                const std::set<int> stencil = ca_cell_stencil(fields, candidate.first, n);
                int incident = 0;
                for (const int dof : stencil)
                {
                    incident += velocities.count(dof);
                }
                if (strict ? incident == 2 * NDIM : incident > 0)
                {
                    patch.insert(stencil.begin(), stencil.end());
                }
            }
        }
        patches.push_back(std::move(patch));
    }
    return patches;
}

// Index sets of the given sets of DOFs, which the caller destroys.
std::vector<IS>
make_index_sets(const std::vector<std::set<int>>& sets)
{
    std::vector<IS> result;
    for (const std::set<int>& dofs : sets)
    {
        const std::vector<PetscInt> indices(dofs.begin(), dofs.end());
        IS is = nullptr;
        const int ierr = ISCreateGeneral(
            PETSC_COMM_SELF, static_cast<PetscInt>(indices.size()), indices.data(), PETSC_COPY_VALUES, &is);
        IBTK_CHKERRQ(ierr);
        result.push_back(is);
    }
    return result;
}

void
destroy_index_sets(std::vector<IS>& sets)
{
    for (IS& is : sets)
    {
        const int ierr = ISDestroy(&is);
        IBTK_CHKERRQ(ierr);
    }
    sets.clear();
}

// Solve the CAV patch with the given DOFs of mat for the entries of residual at those DOFs with an SVD, as the
// subdomain solvers of the CAV cases do, and add the solution to correction at the same DOFs.
void
add_cav_patch_solution(Mat mat, const std::set<int>& dofs, Vec residual, Vec correction)
{
    const std::vector<PetscInt> indices(dofs.begin(), dofs.end());
    const PetscInt n = static_cast<PetscInt>(indices.size());
    IS is = nullptr;
    int ierr = ISCreateGeneral(PETSC_COMM_SELF, n, indices.data(), PETSC_COPY_VALUES, &is);
    IBTK_CHKERRQ(ierr);
    Mat* submat = nullptr;
    ierr = MatCreateSubMatrices(mat, 1, &is, &is, MAT_INITIAL_MATRIX, &submat);
    IBTK_CHKERRQ(ierr);
    Vec local_rhs = nullptr, local_solution = nullptr;
    ierr = MatCreateVecs(submat[0], &local_solution, &local_rhs);
    IBTK_CHKERRQ(ierr);
    PetscScalar* local_values = nullptr;
    ierr = VecGetArray(local_rhs, &local_values);
    IBTK_CHKERRQ(ierr);
    ierr = VecGetValues(residual, n, indices.data(), local_values);
    IBTK_CHKERRQ(ierr);
    ierr = VecRestoreArray(local_rhs, &local_values);
    IBTK_CHKERRQ(ierr);
    KSP ksp = nullptr;
    PC pc = nullptr;
    ierr = KSPCreate(PETSC_COMM_SELF, &ksp);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, submat[0], submat[0]);
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
    ierr = VecSetValues(correction, n, indices.data(), solution_values, ADD_VALUES);
    IBTK_CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(local_solution, &solution_values);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyBegin(correction);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(correction);
    IBTK_CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&local_rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&local_solution);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroySubMatrices(1, &submat);
    IBTK_CHKERRQ(ierr);
    ierr = ISDestroy(&is);
    IBTK_CHKERRQ(ierr);
}

// The SAMRAI_PATCH groups, built from logical cells, the boxes of the SAMRAI patches, and the independently
// enumerated CAV patches, whose seeds are the cells in order. Each box makes one group, which owns the DOFs of its
// cells and solves the CAV patches whose seeds are within standard_width cells of a periodic image of the box, and
// the CAV patches that are larger than the standard Vanka patch of their seed whose seeds are within ib_width cells.
// Each group visits its CAV patches in the order of the traversal.
void
reference_patch_groups(const std::string& traversal,
                       const int standard_width,
                       const int ib_width,
                       const std::vector<Box<NDIM>>& boxes,
                       const CAFields& fields,
                       const std::vector<std::set<int>>& patches,
                       const int n,
                       std::vector<std::vector<int>>& visits,
                       std::vector<std::set<int>>& owned)
{
    std::vector<CACell> seeds;
    for (const auto& field : fields)
    {
        seeds.push_back(field.first);
    }
    const auto near = [&](const CACell& cell, const Box<NDIM>& box, const int width)
    {
        for (int d = 0; d < NDIM; ++d)
        {
            bool image_near = false;
            for (const int shift : { -n, 0, n })
            {
                image_near =
                    image_near || (box.lower(d) - width <= cell[d] + shift && cell[d] + shift <= box.upper(d) + width);
            }
            if (!image_near)
            {
                return false;
            }
        }
        return true;
    };
    visits.clear();
    owned.clear();
    for (const Box<NDIM>& box : boxes)
    {
        std::vector<int> members;
        for (std::size_t k = 0; k < seeds.size(); ++k)
        {
            const bool expanded = patches[k] != ca_cell_stencil(fields, seeds[k], n);
            if (near(seeds[k], box, standard_width) || (expanded && near(seeds[k], box, ib_width)))
            {
                members.push_back(static_cast<int>(k));
            }
        }
        owned.emplace_back();
        for (Box<NDIM>::Iterator b(box); b; b++)
        {
            CACell cell{};
            for (int d = 0; d < NDIM; ++d)
            {
                cell[d] = b()(d);
            }
            owned.back().insert(fields.at(cell).begin(), fields.at(cell).end());
        }
        visits.push_back(members);
        if (traversal == "REVERSE")
        {
            std::reverse(visits.back().begin(), visits.back().end());
        }
        else if (traversal == "SYMMETRIC")
        {
            visits.back().insert(visits.back().end(), members.rbegin() + 1, members.rend());
        }
    }
}

// The action of groups computed independently of the solver: each group starts from a zero correction, visits its
// CAV patches in order, solves each for the residual of the original system with the correction, and adds the whole
// solution to the correction. With FULL output the completed correction is added to the result on the CAV patches of
// the group, and with OWNED output it is written at the DOFs that the group owns, which must partition the DOFs. The
// groups are applied in the given order.
void
reference_group_action(Mat mat,
                       Vec rhs,
                       Vec result,
                       const std::vector<std::set<int>>& patches,
                       const std::vector<std::vector<int>>& visits,
                       const std::vector<std::set<int>>& owned,
                       const bool owned_output,
                       const std::vector<std::size_t>& group_order)
{
    Vec correction = nullptr, residual = nullptr;
    int ierr = VecDuplicate(rhs, &correction);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(rhs, &residual);
    IBTK_CHKERRQ(ierr);
    ierr = VecSet(result, 0.0);
    IBTK_CHKERRQ(ierr);
    PetscInt n_dofs = 0;
    ierr = VecGetSize(rhs, &n_dofs);
    IBTK_CHKERRQ(ierr);
    std::vector<int> owners(n_dofs, 0);
    for (const std::size_t g : group_order)
    {
        ierr = VecSet(correction, 0.0);
        IBTK_CHKERRQ(ierr);
        for (const int k : visits[g])
        {
            ierr = MatMult(mat, correction, residual);
            IBTK_CHKERRQ(ierr);
            ierr = VecAYPX(residual, -1.0, rhs);
            IBTK_CHKERRQ(ierr);
            add_cav_patch_solution(mat, patches[k], residual, correction);
        }
        std::set<int> output_dofs;
        if (owned_output)
        {
            output_dofs = owned[g];
            for (const int dof : output_dofs)
            {
                ++owners[dof];
            }
        }
        else
        {
            for (const int k : visits[g])
            {
                output_dofs.insert(patches[k].begin(), patches[k].end());
            }
        }
        for (const int dof : output_dofs)
        {
            PetscScalar value = 0.0;
            ierr = VecGetValues(correction, 1, &dof, &value);
            IBTK_CHKERRQ(ierr);
            ierr = VecSetValue(result, dof, value, ADD_VALUES);
            IBTK_CHKERRQ(ierr);
        }
        ierr = VecAssemblyBegin(result);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyEnd(result);
        IBTK_CHKERRQ(ierr);
    }
    if (owned_output && std::any_of(owners.begin(), owners.end(), [](const int count) { return count != 1; }))
    {
        TBOX_ERROR("Failed check: the reference groups do not partition the DOFs.\n");
    }
    ierr = VecDestroy(&correction);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&residual);
    IBTK_CHKERRQ(ierr);
}

// Apply subdomain relaxation on the CAV patches of the application problem with each composition and output, and
// compare its action with the reference action on the independently enumerated patches, each of which owns the DOFs
// of its seed cell. The solver is also applied again, and rebuilt with a changed construction matrix, which changes
// the patches and the operator.
int
check_cav_composition(Pointer<PatchLevel<NDIM>> level,
                      SAMRAIVectorReal<NDIM, double>& x,
                      SAMRAIVectorReal<NDIM, double>& b,
                      Vec rhs,
                      const CAFields& fields,
                      const int n)
{
    PetscInt n_dofs = 0;
    int ierr = VecGetSize(rhs, &n_dofs);
    IBTK_CHKERRQ(ierr);
    std::vector<std::set<int>> owned;
    for (const auto& field : fields)
    {
        owned.emplace_back(field.second.begin(), field.second.end());
    }
    Vec action = nullptr, expected = nullptr, repeated = nullptr;
    for (Vec* v : { &action, &expected, &repeated })
    {
        ierr = VecDuplicate(rhs, v);
        IBTK_CHKERRQ(ierr);
    }
    std::map<std::string, Vec> actions;
    for (const bool multiplicative : { false, true })
    {
        for (const bool owned_output : { false, true })
        {
            const std::string name =
                std::string(multiplicative ? "MULTIPLICATIVE" : "ADDITIVE") + (owned_output ? " OWNED" : " FULL");
            Mat elasticity = nullptr;
            ierr = MatCreateSeqDense(PETSC_COMM_WORLD, n_dofs, n_dofs, nullptr, &elasticity);
            IBTK_CHKERRQ(ierr);
            ierr = MatAssemblyBegin(elasticity, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);
            ierr = MatAssemblyEnd(elasticity, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);
            Pointer<MemoryDatabase> db = new MemoryDatabase("composition_solver");
            db->putString("ksp_type", "preonly");
            db->putString("pc_type", "shell");
            db->putBool("initial_guess_nonzero", false);
            db->putBool("check_subdomain_coverage", true);
            db->putString("subdomain_construction", "COUPLING_AWARE");
            Pointer<Database> relaxation_db = db->putDatabase("subdomain_relaxation");
            relaxation_db->putString("composition", multiplicative ? "MULTIPLICATIVE" : "ADDITIVE");
            if (multiplicative)
            {
                relaxation_db->putString("grouping", "RANK");
            }
            relaxation_db->putString("output", owned_output ? "OWNED" : "FULL");
            StaggeredStokesPETScLevelSolver solver("composition_solver", db, "shell_");
            PoissonSpecifications coefficients("coefficients");
            coefficients.setCConstant(1.0);
            coefficients.setDConstant(-1.0);
            solver.setVelocityPoissonSpecifications(coefficients);
            solver.setComponentsHaveNullSpace(false, true);
            for (int cycle = 0; cycle < 2; ++cycle)
            {
                const std::vector<std::set<int>> patches = cav_application_patches(fields, n, false, cycle, elasticity);
                solver.setCouplingAwareASMConstructionMat(elasticity);
                solver.setAugmentedOperatorMat(elasticity);
                solver.initializeSolverState(x, b);
                std::vector<IS>* overlap = nullptr;
                std::vector<IS>* partition = nullptr;
                solver.getASMSubdomains(&partition, &overlap);
                if (ca_read_sets(*overlap) != patches || ca_read_sets(*partition) != owned)
                {
                    TBOX_ERROR(
                        "Failed check: the solver's CAV patches or their owned DOFs differ from the expected "
                        "ones.\n");
                }
                Mat level_mat = nullptr;
                PC pc = nullptr;
                ierr = KSPGetOperators(solver.getPETScKSP(), &level_mat, nullptr);
                IBTK_CHKERRQ(ierr);
                ierr = KSPGetPC(solver.getPETScKSP(), &pc);
                IBTK_CHKERRQ(ierr);
                std::vector<IS> patch_sets = make_index_sets(patches), owned_sets = make_index_sets(owned);
                reference_action(
                    level_mat, rhs, expected, patch_sets, owned_sets, multiplicative, owned_output, 1.0, "FORWARD");
                destroy_index_sets(patch_sets);
                destroy_index_sets(owned_sets);
                // The output holds other values before each application.
                for (Vec v : { action, repeated })
                {
                    ierr = VecSet(v, OUTPUT_SENTINEL);
                    IBTK_CHKERRQ(ierr);
                    ierr = PCApply(pc, rhs, v);
                    IBTK_CHKERRQ(ierr);
                }
                ierr = VecAXPY(repeated, -1.0, action);
                IBTK_CHKERRQ(ierr);
                const double repeated_difference = norm_inf(repeated);
                ierr = VecAXPY(action, -1.0, expected);
                IBTK_CHKERRQ(ierr);
                const double error = norm_inf(action) / norm_inf(expected);
                if (!(error <= 1.0e-9) || repeated_difference != 0.0)
                {
                    TBOX_ERROR("Failed check: " << name << " error = " << error
                                                << ", repeated difference = " << repeated_difference << ".\n");
                }
                plog << name << (cycle == 0 ? "" : " rebuilt") << " error = " << error << '\n';
                if (cycle == 0)
                {
                    ierr = VecDuplicate(expected, &actions[name]);
                    IBTK_CHKERRQ(ierr);
                    ierr = VecCopy(expected, actions[name]);
                    IBTK_CHKERRQ(ierr);
                }
                solver.deallocateSolverState();
                solver.setAugmentedOperatorMat(nullptr);
            }
            ierr = MatDestroy(&elasticity);
            IBTK_CHKERRQ(ierr);
        }
    }
    // SAMRAI_PATCH grouping makes one group for each of the four SAMRAI patches of the level, and the spring of the CAV
    // patches couples DOFs of two of them. Each case is compared with the independent reference of its groups.
    std::vector<Box<NDIM>> boxes;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        boxes.push_back(level->getPatch(p())->getBox());
    }
    CACell remote{};
    remote.fill(2);
    const auto box_of = [&](const CACell& cell)
    {
        hier::Index<NDIM> index;
        for (int d = 0; d < NDIM; ++d)
        {
            index(d) = cell[d];
        }
        return std::find_if(boxes.begin(), boxes.end(), [&](const Box<NDIM>& box) { return box.contains(index); }) -
               boxes.begin();
    };
    if (boxes.size() != 4 || box_of(CACell{}) == box_of(remote))
    {
        TBOX_ERROR("Failed check: the level needs four SAMRAI patches, and the spring must couple two of them.\n");
    }
    const auto group_order = [](const std::size_t n_groups, const bool reverse)
    {
        std::vector<std::size_t> order(n_groups);
        for (std::size_t g = 0; g < n_groups; ++g)
        {
            order[g] = reverse ? n_groups - 1 - g : g;
        }
        return order;
    };
    struct PatchCase
    {
        std::string traversal;
        int ib_width;
        bool owned_output;
    };
    for (const PatchCase& patch_case : std::vector<PatchCase>{ { "FORWARD", 0, true },
                                                               { "REVERSE", 0, true },
                                                               { "SYMMETRIC", 0, true },
                                                               { "FORWARD", 0, false },
                                                               { "SYMMETRIC", 1, true },
                                                               { "SYMMETRIC", 1, false } })
    {
        const std::string name = std::string("SAMRAI_PATCH ") + (patch_case.owned_output ? "OWNED " : "FULL ") +
                                 patch_case.traversal + " width " + std::to_string(patch_case.ib_width);
        Mat elasticity = nullptr;
        ierr = MatCreateSeqDense(PETSC_COMM_WORLD, n_dofs, n_dofs, nullptr, &elasticity);
        IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyBegin(elasticity, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyEnd(elasticity, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);
        Pointer<MemoryDatabase> db = new MemoryDatabase("patch_solver");
        db->putString("ksp_type", "preonly");
        db->putString("pc_type", "shell");
        db->putBool("initial_guess_nonzero", false);
        db->putBool("check_subdomain_coverage", true);
        db->putString("subdomain_construction", "COUPLING_AWARE");
        Pointer<Database> ca_db = db->putDatabase("coupling_aware_subdomains");
        ca_db->putInteger("group_standard_seed_ghost_width", 0);
        ca_db->putInteger("group_ib_seed_ghost_width", patch_case.ib_width);
        Pointer<Database> relaxation_db = db->putDatabase("subdomain_relaxation");
        relaxation_db->putString("composition", "MULTIPLICATIVE");
        relaxation_db->putString("grouping", "SAMRAI_PATCH");
        relaxation_db->putString("output", patch_case.owned_output ? "OWNED" : "FULL");
        relaxation_db->putString("traversal", patch_case.traversal);
        StaggeredStokesPETScLevelSolver solver("patch_solver", db, "shell_");
        PoissonSpecifications coefficients("coefficients");
        coefficients.setCConstant(1.0);
        coefficients.setDConstant(-1.0);
        solver.setVelocityPoissonSpecifications(coefficients);
        solver.setComponentsHaveNullSpace(false, true);
        for (int cycle = 0; cycle < 2; ++cycle)
        {
            const std::vector<std::set<int>> patches = cav_application_patches(fields, n, false, cycle, elasticity);
            solver.setCouplingAwareASMConstructionMat(elasticity);
            solver.setAugmentedOperatorMat(elasticity);
            solver.initializeSolverState(x, b);
            std::vector<IS>* overlap = nullptr;
            std::vector<IS>* partition = nullptr;
            solver.getASMSubdomains(&partition, &overlap);
            if (ca_read_sets(*overlap) != patches)
            {
                TBOX_ERROR("Failed check: the solver's CAV patches differ from the expected ones.\n");
            }
            Mat level_mat = nullptr;
            PC pc = nullptr;
            ierr = KSPGetOperators(solver.getPETScKSP(), &level_mat, nullptr);
            IBTK_CHKERRQ(ierr);
            ierr = KSPGetPC(solver.getPETScKSP(), &pc);
            IBTK_CHKERRQ(ierr);
            std::vector<std::vector<int>> visits;
            std::vector<std::set<int>> group_owned;
            reference_patch_groups(
                patch_case.traversal, 0, patch_case.ib_width, boxes, fields, patches, n, visits, group_owned);
            reference_group_action(level_mat,
                                   rhs,
                                   expected,
                                   patches,
                                   visits,
                                   group_owned,
                                   patch_case.owned_output,
                                   group_order(visits.size(), false));
            for (Vec v : { action, repeated })
            {
                ierr = VecSet(v, OUTPUT_SENTINEL);
                IBTK_CHKERRQ(ierr);
                ierr = PCApply(pc, rhs, v);
                IBTK_CHKERRQ(ierr);
            }
            ierr = VecAXPY(repeated, -1.0, action);
            IBTK_CHKERRQ(ierr);
            const double repeated_difference = norm_inf(repeated);
            ierr = VecAXPY(action, -1.0, expected);
            IBTK_CHKERRQ(ierr);
            const double error = norm_inf(action) / norm_inf(expected);
            if (!(error <= 1.0e-9) || repeated_difference != 0.0)
            {
                TBOX_ERROR("Failed check: " << name << " error = " << error
                                            << ", repeated difference = " << repeated_difference << ".\n");
            }
            plog << name << (cycle == 0 ? "" : " rebuilt") << " error = " << error << '\n';
            if (cycle == 0)
            {
                ierr = VecDuplicate(expected, &actions[name]);
                IBTK_CHKERRQ(ierr);
                ierr = VecCopy(expected, actions[name]);
                IBTK_CHKERRQ(ierr);
                if (patch_case.traversal == "FORWARD" && patch_case.ib_width == 0 && patch_case.owned_output)
                {
                    // The groups are independent, so a reference that applies them in the reverse order agrees; the
                    // solver itself is not reordered. Some CAV patch of a group also contains the pressure DOF of a
                    // cell that the group does not own, so the group's correction has entries outside its owned DOFs.
                    reference_group_action(level_mat,
                                           rhs,
                                           action,
                                           patches,
                                           visits,
                                           group_owned,
                                           patch_case.owned_output,
                                           group_order(visits.size(), true));
                    ierr = VecAXPY(action, -1.0, expected);
                    IBTK_CHKERRQ(ierr);
                    plog << "reordered reference difference = " << norm_inf(action) / norm_inf(expected) << '\n';
                    bool reaches_beyond = false;
                    for (std::size_t g = 0; g < visits.size(); ++g)
                    {
                        for (const int k : visits[g])
                        {
                            for (const auto& field : fields)
                            {
                                reaches_beyond = reaches_beyond || (patches[k].count(field.second[NDIM]) &&
                                                                    !group_owned[g].count(field.second[NDIM]));
                            }
                        }
                    }
                    if (!reaches_beyond)
                    {
                        TBOX_ERROR("Failed check: no CAV patch of a group reaches beyond the cells of the group.\n");
                    }
                }
                if (patch_case.ib_width > 0 && patch_case.owned_output)
                {
                    // With a positive IB width, IB-expanded CAV patches are solved by several groups.
                    std::size_t occurrences = 0;
                    for (const std::vector<int>& group : visits)
                    {
                        occurrences += (group.size() + 1) / 2;
                    }
                    if (!(occurrences > patches.size()))
                    {
                        TBOX_ERROR("Failed check: no CAV patch is solved by several groups.\n");
                    }
                    plog << "group memberships = " << occurrences << " of " << patches.size() << " CAV patches\n";
                }
            }
            solver.deallocateSolverState();
            solver.setAugmentedOperatorMat(nullptr);
        }
        ierr = MatDestroy(&elasticity);
        IBTK_CHKERRQ(ierr);
    }
    // With one rank the two outputs of the multiplicative composition are the correction of its one group, and the
    // other actions differ.
    for (const auto& pair : { std::make_pair("ADDITIVE FULL", "ADDITIVE OWNED"),
                              std::make_pair("MULTIPLICATIVE FULL", "MULTIPLICATIVE OWNED"),
                              std::make_pair("ADDITIVE OWNED", "MULTIPLICATIVE OWNED"),
                              std::make_pair("SAMRAI_PATCH FULL FORWARD width 0", "SAMRAI_PATCH OWNED FORWARD width 0"),
                              std::make_pair("SAMRAI_PATCH OWNED FORWARD width 0", "MULTIPLICATIVE OWNED"),
                              std::make_pair("SAMRAI_PATCH OWNED FORWARD width 0", "ADDITIVE OWNED") })
    {
        ierr = VecWAXPY(action, -1.0, actions[pair.second], actions[pair.first]);
        IBTK_CHKERRQ(ierr);
        plog << pair.first << " - " << pair.second << " = " << norm_inf(action) << '\n';
    }
    for (auto& entry : actions)
    {
        ierr = VecDestroy(&entry.second);
        IBTK_CHKERRQ(ierr);
    }
    for (Vec* v : { &action, &expected, &repeated })
    {
        ierr = VecDestroy(v);
        IBTK_CHKERRQ(ierr);
    }
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
    const bool cav = test->keyExists("cav_application");
    const std::string cav_scenario = cav ? test->getString("cav_application") : "";
    const bool lifetime = cav || test->getBoolWithDefault("lifetime", false);
    const bool application_subdomain_solver = test->getBoolWithDefault("application_subdomain_solver", false);
    const bool null_subdomain_solver = test->getBoolWithDefault("null_subdomain_solver", false);
    const bool subdomain_solver_unused = test->getBoolWithDefault("subdomain_solver_unused", false);
    const bool replace_initialized = test->getBoolWithDefault("replace_initialized", false);
    const bool named_subdomain_solver = test->getBoolWithDefault("named_subdomain_solver", false);
    const bool uneven_subdomains = test->getBoolWithDefault("uneven_subdomains", false);
    const bool diagonal_operator = test->getBoolWithDefault("diagonal_operator", false);
    const double small_diagonal = test->getDoubleWithDefault("small_diagonal", 1.0e-12);
    // "rank_one", "upper" or "symmetric" supplies a dense operator times operator_scale.
    const std::string dense_operator = test->getStringWithDefault("dense_operator", "");
    const double operator_scale = test->getDoubleWithDefault("operator_scale", 1.0);
    const bool supplied_operator = diagonal_operator || !dense_operator.empty();
    const bool all_blas_modes = test->getBoolWithDefault("all_blas_modes", false);
    const std::vector<std::string> solver_types =
        all_blas_modes ? std::vector<std::string>{ "", "svd", "lu", "symmetric-indefinite", "qr" } :
                         std::vector<std::string>{ "" };
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
    CAFields cav_fields;
    std::vector<int> dofs;
    StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(dofs, udi, pdi, level);
    if (cav)
    {
        for (PatchLevel<NDIM>::Iterator patch_number(level); patch_number; patch_number++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(patch_number());
            Pointer<SideData<NDIM, int>> velocity = patch->getPatchData(udi);
            Pointer<CellData<NDIM, int>> pressure = patch->getPatchData(pdi);
            Pointer<CellData<NDIM, double>> pressure_rhs = patch->getPatchData(hi);
            for (Box<NDIM>::Iterator it(patch->getBox()); it; it++)
            {
                CACell cell{};
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    cell[axis] = it()(axis);
                }
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    cav_fields[cell][axis] = (*velocity)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower));
                }
                cav_fields[cell][NDIM] = (*pressure)(it());
                (*pressure_rhs)(it()) = 0.125 * std::sin(wavenumber * it()(0));
            }
        }
    }
    Vec rhs = nullptr, expected = nullptr, actual = nullptr;
    int ierr = VecCreateMPI(PETSC_COMM_WORLD, dofs[IBTK_MPI::getRank()], PETSC_DETERMINE, &rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(rhs, &expected);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(rhs, &actual);
    IBTK_CHKERRQ(ierr);
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(rhs, fi, udi, hi, pdi, level);
    // The input supplies the settings of each case, and the settings that every case shares are added here.
    Pointer<Database> db = input->getDatabase("level_solver");
    db->putString("ksp_type", "preonly");
    if (!db->keyExists("pc_type"))
    {
        db->putString("pc_type", "shell");
    }
    db->putBool("initial_guess_nonzero", false);
    if (!db->keyExists("check_subdomain_coverage"))
    {
        db->putBool("check_subdomain_coverage", true);
    }
    db->putInteger("max_iterations", 1);
    // The CAV cases select pressure-cell patches through the default seed type, and the patches that the test
    // expects follow the closure policy of the input.
    const bool cav_strict =
        read_setting(db->keyExists("coupling_aware_subdomains") ? db->getDatabase("coupling_aware_subdomains") :
                                                                  Pointer<Database>(),
                     "closure_policy",
                     "RELAXED") == "STRICT";
    int box_size[NDIM];
    std::fill_n(box_size, NDIM, 4);
    db->putIntegerArray("subdomain_box_size", box_size, NDIM);
    if (cav_scenario == "composition")
    {
        plog << std::setprecision(12);
        const int status = check_cav_composition(level, x, b, rhs, cav_fields, input->getInteger("N"));
        for (Vec* vector : { &rhs, &expected, &actual })
        {
            ierr = VecDestroy(vector);
            IBTK_CHKERRQ(ierr);
        }
        return status;
    }
    // The settings of the reference action. The cases that stop with an error before the reference is used need not
    // supply them.
    Pointer<Database> relaxation_db =
        db->keyExists("subdomain_relaxation") ? db->getDatabase("subdomain_relaxation") : Pointer<Database>();
    Pointer<Database> solver_db = relaxation_db && relaxation_db->keyExists("subdomain_solver") ?
                                      relaxation_db->getDatabase("subdomain_solver") :
                                      Pointer<Database>();
    const bool multiplicative = read_setting(relaxation_db, "composition", "") == "MULTIPLICATIVE";
    const bool owned_output = read_setting(relaxation_db, "output", "") == "OWNED";
    const std::string traversal = read_setting(relaxation_db, "traversal", "FORWARD");
    const std::string subdomain_solver_type = read_setting(solver_db, "type", "petsc");
    // A named application subdomain solver, which the type of the subdomain_solver database selects. The factory
    // reads the scale of its solutions from that database, and each call creates an independent solver with
    // counters of its own.
    std::deque<SubdomainSolverCounters> factory_counters;
    PETScLevelSolver::SubdomainSolverFactories factories;
    if (named_subdomain_solver)
    {
        factories.emplace_back(test->getStringWithDefault("factory_name", "application"),
                               [&factory_counters](Pointer<Database> factory_db)
                               {
                                   check_database_keys("application factory", factory_db, { "type", "solution_scale" });
                                   factory_counters.emplace_back();
                                   return PETScLevelSolverSubdomainSolver(
                                       std::in_place_type<ScaledCountingSubdomainSolver>,
                                       factory_counters.back(),
                                       factory_db->getDouble("solution_scale"));
                               });
    }
    const double application_scale = named_subdomain_solver && solver_db->keyExists("solution_scale") ?
                                         solver_db->getDouble("solution_scale") :
                                         SUBDOMAIN_SOLVER_SCALE;
    plog << std::setprecision(12);
    for (const std::string& solver_type : solver_types)
    {
        if (!solver_type.empty())
        {
            solver_db->putString("blas_lapack_subdomain_solver_type", solver_type);
        }
        if (test->keyExists("pc_type_option"))
        {
            // Override the preconditioner type through the PETSc options database, which initialization applies.
            ierr = PetscOptionsSetValue(nullptr, "-shell_pc_type", test->getString("pc_type_option").c_str());
            IBTK_CHKERRQ(ierr);
        }
        // The level solver owns its subdomain solver, so the counters must outlive it.
        SubdomainSolverCounters application_counters, replacement_counters;
        CommunicationProbe solver("shell_solver", db, "shell_", factories);
        PoissonSpecifications coefficients("coefficients");
        coefficients.setCConstant(1.0);
        coefficients.setDConstant(-1.0);
        solver.break_partition = test->getBoolWithDefault("break_partition", false);
        solver.setVelocityPoissonSpecifications(coefficients);
        solver.setComponentsHaveNullSpace(false, !supplied_operator);
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
        Mat elasticity = nullptr, fixed_augmentation = nullptr;
        std::vector<std::set<int>> cav_expected, previous_patches;
        if (cav)
        {
            ierr = MatCreateSeqDense(PETSC_COMM_WORLD, dofs.front(), dofs.front(), nullptr, &elasticity);
            IBTK_CHKERRQ(ierr);
            ierr = MatAssemblyBegin(elasticity, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);
            ierr = MatAssemblyEnd(elasticity, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);
            cav_expected = cav_application_patches(cav_fields, input->getInteger("N"), cav_strict, 0, elasticity);
            if (cav_scenario != "missing_matrix")
            {
                PetscInt before = 0, after = 0;
                ierr = PetscObjectGetReference(reinterpret_cast<PetscObject>(elasticity), &before);
                IBTK_CHKERRQ(ierr);
                solver.setCouplingAwareASMConstructionMat(elasticity);
                ierr = PetscObjectGetReference(reinterpret_cast<PetscObject>(elasticity), &after);
                IBTK_CHKERRQ(ierr);
                if (before != after)
                {
                    TBOX_ERROR("Supplying the construction matrix changed its reference count.\n");
                }
            }
            // The solve matrix is augmented with the elasticity matrix, or, to show that the construction matrix
            // is a separate input, with a copy that stays fixed while the construction matrix changes.
            if (test->getBoolWithDefault("fixed_augmentation", false))
            {
                ierr = MatDuplicate(elasticity, MAT_COPY_VALUES, &fixed_augmentation);
                IBTK_CHKERRQ(ierr);
            }
            solver.setAugmentedOperatorMat(fixed_augmentation ? fixed_augmentation : elasticity);
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
        if (named_subdomain_solver)
        {
            // The factory runs when a level solver first sets up its shell preconditioner, not at each
            // initialization, and another level solver has an independent subdomain solver.
            if (factory_counters.size() != 1)
            {
                TBOX_ERROR("Failed check: the factory did not run when the shell preconditioner was set up.\n");
            }
            counts = &factory_counters.front();
            CommunicationProbe other_solver("other_shell_solver", db, "shell_", factories);
            other_solver.setVelocityPoissonSpecifications(coefficients);
            other_solver.setComponentsHaveNullSpace(false, !supplied_operator);
            other_solver.initializeSolverState(x, b);
            other_solver.deallocateSolverState();
            if (factory_counters.size() != 2)
            {
                TBOX_ERROR("Failed check: each level solver needs its own subdomain solver from the factory.\n");
            }
        }
        if (test->getBoolWithDefault("report_communication", false))
        {
            plog << "restriction communicates = " << solver.restrictionCommunicates() << '\n';
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
        if (cav_scenario == "initialized_setter")
        {
            solver.setCouplingAwareASMConstructionMat(elasticity);
        }
        if (cav && cav_scenario != "apply")
        {
            solver.deallocateSolverState();
            ierr = MatDestroy(&elasticity);
            IBTK_CHKERRQ(ierr);
            ierr = MatDestroy(&fixed_augmentation);
            IBTK_CHKERRQ(ierr);
            return 0;
        }
        if (lifetime && !diagonal_operator && !cav)
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
        // An application subdomain solver is also replaced by another one after the cycles.
        const int reinitialization_cycles = lifetime ? 2 : 1;
        SubdomainSolverCounters original_counters;
        for (int cycle = 0; cycle < reinitialization_cycles + (application_subdomain_solver ? 1 : 0); ++cycle)
        {
            const bool replacement = cycle == reinitialization_cycles;
            const char* const label = replacement ? "replacement " : "";
            if (replacement)
            {
                original_counters = *counts;
                PETScLevelSolverSubdomainSolver replaced(std::in_place_type<ScaledCountingSubdomainSolver>,
                                                         replacement_counters);
                solver.setSubdomainSolver(std::move(replaced));
                if (replaced)
                {
                    TBOX_ERROR("Failed check: replacing a subdomain solver did not transfer ownership.\n");
                }
                solver.initializeSolverState(x, b);
            }
            Mat mat = nullptr;
            PC pc = nullptr;
            ierr = KSPGetOperators(solver.getPETScKSP(), &mat, nullptr);
            IBTK_CHKERRQ(ierr);
            ierr = KSPGetPC(solver.getPETScKSP(), &pc);
            IBTK_CHKERRQ(ierr);
            PCType pc_type = nullptr;
            ierr = PCGetType(pc, &pc_type);
            IBTK_CHKERRQ(ierr);
            const char* shell_name = nullptr;
            ierr = PCShellGetName(pc, &shell_name);
            IBTK_CHKERRQ(ierr);
            if (std::string(pc_type) != "shell" || std::string(shell_name) != "subdomain_relaxation" ||
                (lifetime && !cav && mat != supplied))
            {
                TBOX_ERROR(
                    "Failed check: the preconditioner is not the subdomain_relaxation shell, or the operator "
                    "is not the supplied one.\n");
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
                if (cav)
                {
                    // Each patch owns the DOFs of its seed cell.
                    std::vector<std::set<int>> seed_cell_dofs;
                    for (const auto& field : cav_fields)
                    {
                        seed_cell_dofs.emplace_back(field.second.begin(), field.second.end());
                    }
                    const std::vector<std::set<int>> patches = ca_read_sets(*overlap);
                    if (ca_read_sets(*partition) != seed_cell_dofs || patches != cav_expected ||
                        (cycle == 1 && patches == previous_patches))
                    {
                        TBOX_ERROR(
                            "Failed check: the owned DOFs are not those of the seed cells, patches != cav_expected, "
                            "or (cycle == 1 && patches == previous_patches).\n");
                    }
                    previous_patches = patches;
                    plog << "pressure_patches = " << patches.size() << "\npartition_size = " << partition->size()
                         << '\n';
                }
                // The supplied subdomain solver, not the built-in one, determines the action.
                const double scale = counts ? application_scale : 1.0;
                reference_action(
                    mat, rhs, expected, *overlap, *partition, multiplicative, owned_output, scale, traversal);
                if (cycle == 0 && solver_type == solver_types.front())
                {
                    // The other compositions and outputs give different actions, except that with one rank both
                    // outputs of the multiplicative composition are the correction of the one group.
                    for (const bool other_multiplicative : { false, true })
                    {
                        for (const bool other_owned_output : { false, true })
                        {
                            // OWNED output needs the nonoverlapping sets, which CAV patches do not have.
                            if ((other_multiplicative == multiplicative && other_owned_output == owned_output) ||
                                (other_owned_output && partition->size() != overlap->size()))
                            {
                                continue;
                            }
                            reference_action(mat,
                                             rhs,
                                             actual,
                                             *overlap,
                                             *partition,
                                             other_multiplicative,
                                             other_owned_output,
                                             scale,
                                             traversal);
                            ierr = VecAXPY(actual, -1.0, expected);
                            IBTK_CHKERRQ(ierr);
                            plog << "distance from the " << (other_multiplicative ? "MULTIPLICATIVE " : "ADDITIVE ")
                                 << (other_owned_output ? "OWNED" : "FULL") << " action = " << norm_inf(actual) << '\n';
                        }
                    }
                }
            }
            // Applying the preconditioner twice to an output that holds other values gives the expected action
            // both times.
            double apply_error = 0.0;
            for (int application = 0; application < 2; ++application)
            {
                ierr = VecSet(actual, OUTPUT_SENTINEL);
                IBTK_CHKERRQ(ierr);
                ierr = PCApply(pc, rhs, actual);
                IBTK_CHKERRQ(ierr);
                ierr = VecAXPY(actual, -1.0, expected);
                IBTK_CHKERRQ(ierr);
                apply_error = std::max(apply_error, norm_inf(actual));
            }
            if (!std::isfinite(apply_error) || apply_error > 1.0e-9 * norm_inf(expected))
            {
                TBOX_ERROR("Failed check: applying the preconditioner does not give the expected action.\n");
            }
            plog << label << "apply_error = " << apply_error / norm_inf(expected) << '\n';
            // Left-preconditioned PETSc KSP removes the operator nullspace after PCApply.
            MatNullSpace nullspace = nullptr;
            ierr = MatGetNullSpace(mat, &nullspace);
            IBTK_CHKERRQ(ierr);
            if (nullspace)
            {
                ierr = MatNullSpaceRemove(nullspace, expected);
                IBTK_CHKERRQ(ierr);
            }
            x.setToScalar(0.0);
            if (!solver.solveSystem(x, b))
            {
                TBOX_ERROR("Failed check: !solver.solveSystem(x, b).\n");
            }
            StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(actual, ui, udi, pi, pdi, level);
            const double action_norm = norm_inf(actual);
            if (cav)
            {
                PetscScalar pressure_sum = 0.0;
                const PetscScalar* values = nullptr;
                ierr = VecGetArrayRead(actual, &values);
                IBTK_CHKERRQ(ierr);
                for (const auto& cell : cav_fields)
                {
                    pressure_sum += values[cell.second[NDIM]];
                }
                ierr = VecRestoreArrayRead(actual, &values);
                IBTK_CHKERRQ(ierr);
                const double pressure_mean = static_cast<double>(PetscRealPart(pressure_sum)) / cav_fields.size();
                if (!std::isfinite(pressure_mean) || std::abs(pressure_mean) > 1.0e-9)
                {
                    TBOX_ERROR("Failed check: !std::isfinite(pressure_mean) || std::abs(pressure_mean) > 1.0e-9.\n");
                }
                plog << "pressure_mean = " << pressure_mean << '\n';
            }
            ierr = VecAXPY(actual, -1.0, expected);
            IBTK_CHKERRQ(ierr);
            const double error = norm_inf(actual);
            if (!std::isfinite(error) || error > 1.0e-9 * action_norm || action_norm <= 0.0)
            {
                TBOX_ERROR(
                    "Failed check: !std::isfinite(error) || error > 1.0e-9 * action_norm || action_norm <= 0.0.\n");
            }
            plog << label << "action_norm = " << action_norm << '\n'
                 << label << "error = " << error / action_norm << '\n';
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
                plog << label << "repeated_error = " << repeated_error << '\n';
            }
            std::vector<IS>* cav_overlap = nullptr;
            std::vector<IS>* cav_partition = nullptr;
            if (cav)
            {
                // Query while initialized; the returned containers outlive solver state.
                solver.getASMSubdomains(&cav_partition, &cav_overlap);
            }
            solver.deallocateSolverState();
            if (replacement)
            {
                // The replaced solver is not used again, and the new one is initialized, used, and released once.
                if (*counts != original_counters || replacement_counters.constructions != 1 ||
                    replacement_counters.initializations != 1 || replacement_counters.deallocations != 1 ||
                    replacement_counters.calls < 1 || replacement_counters.subdomain_solves < 1 ||
                    replacement_counters.contract_violations != 0)
                {
                    TBOX_ERROR(
                        "Failed check: the replacement subdomain solver was not the only one used after "
                        "replacement.\n");
                }
            }
            if (cav)
            {
                if (!cav_overlap->empty() || !cav_partition->empty())
                {
                    TBOX_ERROR("Failed check: !cav_overlap->empty() || !cav_partition->empty().\n");
                }
                {
                    PetscReal matrix_norm = 0.0;
                    ierr = MatNorm(elasticity, NORM_INFINITY, &matrix_norm);
                    IBTK_CHKERRQ(ierr);
                    if (matrix_norm != 0.5 * (cycle == 0 ? 1 : 2 * NDIM))
                    {
                        TBOX_ERROR("Failed check: matrix_norm != 0.5 * (cycle == 0 ? 1 : 2 * NDIM).\n");
                    }
                }
                if (cycle == 0)
                {
                    solver.setAugmentedOperatorMat(nullptr);
                    cav_expected =
                        cav_application_patches(cav_fields, input->getInteger("N"), cav_strict, 1, elasticity);
                    solver.setCouplingAwareASMConstructionMat(elasticity);
                    solver.setAugmentedOperatorMat(fixed_augmentation ? fixed_augmentation : elasticity);
                    solver.initializeSolverState(x, b);
                }
            }
            else if (lifetime)
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
        if (named_subdomain_solver && (factory_counters[1].constructions != 1 ||
                                       factory_counters[1].initializations != 1 || factory_counters[1].calls != 0))
        {
            TBOX_ERROR("Failed check: the other level solver's subdomain solver was used.\n");
        }
        ierr = MatDestroy(&elasticity);
        IBTK_CHKERRQ(ierr);
        ierr = MatDestroy(&fixed_augmentation);
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
    return 0;
}
