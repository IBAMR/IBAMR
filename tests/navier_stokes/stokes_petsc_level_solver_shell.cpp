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

#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScLevelSolverSubdomainSolver.h>

#include <CellData.h>
#include <CellVariable.h>
#include <PoissonSpecifications.h>
#include <SAMRAIVectorReal.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <memory>
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

// The action of subdomain relaxation, computed independently of the solver from the level operator and direct
// solves of the subdomain matrices. Every rank gathers the whole input r. ADDITIVE composition solves each
// subdomain of this rank with r(P_i). MULTIPLICATIVE composition, with one group for each rank, keeps the
// correction z of this rank's group as a vector of all of the DOFs and solves each subdomain in turn with
// r(P_i) - A(P_i, :) z before adding its scaled solution to z. FULL output adds every entry of the scaled
// solutions (ADDITIVE) or of z on the subdomains of this rank (MULTIPLICATIVE) to the result, and OWNED output
// only the entries of the nonoverlapping sets.
void
reference_action(Mat mat,
                 Vec rhs,
                 Vec result,
                 const std::vector<IS>& overlap,
                 const std::vector<IS>& partition,
                 const bool multiplicative,
                 const bool owned_output,
                 const double scale)
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
    std::set<PetscInt> solved_dofs, owned_dofs;
    for (PetscInt i = 0; i < n_subdomains; ++i)
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
    const bool lifetime = test->getBoolWithDefault("lifetime", false);
    const bool application_subdomain_solver = test->getBoolWithDefault("application_subdomain_solver", false);
    const bool null_subdomain_solver = test->getBoolWithDefault("null_subdomain_solver", false);
    const bool subdomain_solver_unused = test->getBoolWithDefault("subdomain_solver_unused", false);
    const bool replace_initialized = test->getBoolWithDefault("replace_initialized", false);
    const bool uneven_subdomains = test->getBoolWithDefault("uneven_subdomains", false);
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
    // The input supplies the settings of each case, and the settings that every case shares are added here.
    Pointer<Database> db = input->getDatabase("level_solver");
    db->putString("ksp_type", "preonly");
    if (!db->keyExists("pc_type"))
    {
        db->putString("pc_type", "shell");
    }
    db->putBool("initial_guess_nonzero", false);
    db->putBool("check_subdomain_coverage", true);
    db->putInteger("max_iterations", 1);
    int box_size[NDIM];
    std::fill_n(box_size, NDIM, 4);
    db->putIntegerArray("subdomain_box_size", box_size, NDIM);
    // The level solver owns its subdomain solver, so the counters must outlive it.
    SubdomainSolverCounters application_counters, replacement_counters;
    SubdomainSolverCounters* counts = nullptr;
    CommunicationProbe solver("shell_solver", db, "shell_");
    PoissonSpecifications coefficients("coefficients");
    coefficients.setCConstant(1.0);
    coefficients.setDConstant(-1.0);
    solver.break_partition = test->getBoolWithDefault("break_partition", false);
    solver.setVelocityPoissonSpecifications(coefficients);
    solver.setComponentsHaveNullSpace(false, true);
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
    solver.initializeSolverState(x, b);
    if (test->getBoolWithDefault("report_communication", false))
    {
        plog << "restriction communicates = " << solver.restrictionCommunicates() << '\n';
    }
    if (replace_initialized)
    {
        solver.setSubdomainSolver(
            PETScLevelSolverSubdomainSolver(std::in_place_type<ScaledCountingSubdomainSolver>, replacement_counters));
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
    Mat supplied = nullptr;
    if (lifetime)
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
    Pointer<Database> relaxation_db = db->getDatabase("subdomain_relaxation");
    const bool multiplicative = relaxation_db->getString("composition") == "MULTIPLICATIVE";
    const bool owned_output = relaxation_db->getString("output") == "OWNED";
    plog << std::setprecision(12);
    // An application subdomain solver is also replaced by another one after the cycles.
    const int reinitialization_cycles = lifetime ? 2 : 1;
    SubdomainSolverCounters original_counters;
    for (int cycle = 0; cycle < reinitialization_cycles + (counts ? 1 : 0); ++cycle)
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
            (lifetime && mat != supplied))
        {
            TBOX_ERROR(
                "Failed check: the preconditioner is not the subdomain_relaxation shell, or the operator is "
                "not the supplied one.\n");
        }
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
        const double scale = counts ? SUBDOMAIN_SOLVER_SCALE : 1.0;
        reference_action(mat, rhs, expected, *overlap, *partition, multiplicative, owned_output, scale);
        if (cycle == 0)
        {
            // The other compositions and outputs give different actions, except that with one rank both outputs of
            // the multiplicative composition are the correction of the one group.
            for (const bool other_multiplicative : { false, true })
            {
                for (const bool other_owned_output : { false, true })
                {
                    if (other_multiplicative == multiplicative && other_owned_output == owned_output)
                    {
                        continue;
                    }
                    reference_action(
                        mat, rhs, actual, *overlap, *partition, other_multiplicative, other_owned_output, scale);
                    ierr = VecAXPY(actual, -1.0, expected);
                    IBTK_CHKERRQ(ierr);
                    plog << "distance from the " << (other_multiplicative ? "MULTIPLICATIVE " : "ADDITIVE ")
                         << (other_owned_output ? "OWNED" : "FULL") << " action = " << norm_inf(actual) << '\n';
                }
            }
        }
        // Applying the preconditioner twice to an output that holds other values gives the reference action both
        // times.
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
        if (!std::isfinite(apply_error) || apply_error > 1.0e-9)
        {
            TBOX_ERROR("Failed check: applying the preconditioner does not give the reference action.\n");
        }
        plog << label << "apply_error = " << apply_error << '\n';
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
        ierr = VecAXPY(actual, -1.0, expected);
        IBTK_CHKERRQ(ierr);
        const double error = norm_inf(actual);
        if (!std::isfinite(error) || error > 1.0e-9 || action_norm <= 0.0)
        {
            TBOX_ERROR("Failed check: !std::isfinite(error) || error > 1.0e-9 || action_norm <= 0.0.\n");
        }
        plog << label << "action_norm = " << action_norm << '\n' << label << "error = " << error << '\n';
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
                    "Failed check: the replacement subdomain solver was not the only one used after replacement.\n");
            }
        }
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
                solver.initializeSolverState(x, b);
            }
        }
    }
    ierr = MatDestroy(&supplied);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&expected);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&actual);
    IBTK_CHKERRQ(ierr);
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
    return 0;
}
