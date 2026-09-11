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
#include <cmath>
#include <iomanip>
#include <memory>

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

// Gather the RHS independently of the backend's restriction/prolongation
// scatters, solve each local matrix, and write only its partition subset.
void
reference_action(Mat mat,
                 Vec rhs,
                 Vec result,
                 const std::vector<IS>& overlap,
                 const std::vector<IS>& partition,
                 const bool legacy)
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
    if (!legacy)
    {
        ierr = VecScatterBegin(gather, rhs, gathered, INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
        ierr = VecScatterEnd(gather, rhs, gathered, INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
    }
    for (std::size_t i = 0; i < overlap.size(); ++i)
    {
        if (legacy)
        {
            // The preservation case is serial: reproduce its ordered restricted writes.
            ierr = MatMult(mat, result, residual);
            IBTK_CHKERRQ(ierr);
            ierr = VecAYPX(residual, -1.0, rhs);
            IBTK_CHKERRQ(ierr);
            ierr = VecScatterBegin(gather, residual, gathered, INSERT_VALUES, SCATTER_FORWARD);
            IBTK_CHKERRQ(ierr);
            ierr = VecScatterEnd(gather, residual, gathered, INSERT_VALUES, SCATTER_FORWARD);
            IBTK_CHKERRQ(ierr);
        }
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
        for (PetscInt j = 0; j < m; ++j)
        {
            const PetscInt position = static_cast<PetscInt>(std::lower_bound(indices, indices + n, owned[j]) - indices);
            TBOX_ASSERT(position < n && indices[position] == owned[j]);
            ierr = VecSetValue(result, owned[j], solution_values[position], INSERT_VALUES);
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
        if (legacy)
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
    const std::string shell_type = test->getString("shell_pc_type");
    const bool legacy = shell_type == "multiplicative";
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
    db->putString("shell_pc_type", shell_type);
    db->putBool("initial_guess_nonzero", false);
    db->putInteger("max_iterations", 1);
    int box_size[NDIM];
    std::fill_n(box_size, NDIM, boundary ? input->getInteger("N") : 4);
    db->putIntegerArray("subdomain_box_size", box_size, NDIM);
    int failures = 0;
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
                reference_action(mat, rhs, expected, *overlap, *partition, legacy);
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
