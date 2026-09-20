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

// Lifetime of the shell preconditioner state of PETScLevelSolver and its subclasses on one periodic level: the
// number of local subdomains that a solver generates from std::set data, the release of the per-subdomain work
// vectors on deallocation, and repeated initialization. A shell preconditioner without subdomains would leave the
// unpreconditioned Krylov solve unable to converge within the iteration limit, so the solve checks would fail.

#include <ibamr/StaggeredStokesPETScLevelSolver.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CCPoissonPETScLevelSolver.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/SCPoissonPETScLevelSolver.h>

#include <PoissonSpecifications.h>
#include <RobinBcCoefStrategy.h>

#include "level_solver_test_utilities.h"

#include <ibamr/app_namespaces.h>

using namespace level_solver_test;

namespace
{
// Solve with the installed KSP, for the residual of a consistent right-hand side. The Stokes operator of a periodic
// level is singular, so the solution is compared only through its residual.
void
require_level_solve(PETScLevelSolver& solver, const std::string& name)
{
    Mat matrix;
    PetscErrorCode ierr = KSPGetOperators(solver.getPETScKSP(), &matrix, nullptr);
    IBTK_CHKERRQ(ierr);
    Vec exact, rhs, solution, residual;
    ierr = MatCreateVecs(matrix, &exact, &rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(exact, &solution);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(exact, &residual);
    IBTK_CHKERRQ(ierr);
    PetscInt first, last;
    ierr = VecGetOwnershipRange(exact, &first, &last);
    IBTK_CHKERRQ(ierr);
    PetscScalar* values;
    ierr = VecGetArray(exact, &values);
    IBTK_CHKERRQ(ierr);
    for (PetscInt i = first; i < last; ++i) values[i - first] = std::sin(0.13 * i) + 0.5;
    ierr = VecRestoreArray(exact, &values);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(matrix, exact, rhs);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSolve(solver.getPETScKSP(), rhs, solution);
    IBTK_CHKERRQ(ierr);
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(solver.getPETScKSP(), &reason);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(matrix, solution, residual);
    IBTK_CHKERRQ(ierr);
    ierr = VecAXPY(residual, -1.0, rhs);
    IBTK_CHKERRQ(ierr);
    PetscReal residual_norm, rhs_norm;
    ierr = VecNorm(residual, NORM_2, &residual_norm);
    IBTK_CHKERRQ(ierr);
    ierr = VecNorm(rhs, NORM_2, &rhs_norm);
    IBTK_CHKERRQ(ierr);
    for (Vec* v : { &exact, &rhs, &solution, &residual })
    {
        ierr = VecDestroy(v);
        IBTK_CHKERRQ(ierr);
    }
    if (reason <= 0 || !std::isfinite(residual_norm) || residual_norm > 1.0e-9 * rhs_norm)
    {
        TBOX_ERROR(name << ": the level solve did not converge; reason = " << reason
                        << ", relative residual = " << residual_norm / rhs_norm << ".\n");
    }
}

// Initialize and deallocate twice. Each time the solver must have 8 subdomains, solve, and hold no shell vector
// after deallocation apart from the references retained here. Return the total size of its overlapping subdomains.
template <class Solver>
PetscInt
check_shell_state(LevelSolverProbe<Solver>& solver,
                  Pointer<HierarchyVector> x,
                  Pointer<HierarchyVector> b,
                  const std::string& name,
                  const bool singular,
                  const bool multiplicative)
{
    PetscInt overlap_total = 0;
    for (int cycle = 0; cycle < 2; ++cycle)
    {
        solver.initializeSolverState(*x, *b);
        std::vector<IS>*nonoverlap, *overlap;
        solver.getASMSubdomains(&nonoverlap, &overlap);
        if (overlap->size() != 8 || nonoverlap->size() != 8)
        {
            TBOX_ERROR(name << ": expected 8 subdomains, found " << overlap->size() << " overlapping and "
                            << nonoverlap->size() << " nonoverlapping.\n");
        }
        overlap_total = 0;
        for (IS is : *overlap)
        {
            PetscInt n;
            PetscErrorCode ierr = ISGetSize(is, &n);
            IBTK_CHKERRQ(ierr);
            overlap_total += n;
        }
        if (singular)
        {
            require_level_solve(solver, name);
        }
        else if (!check_level_solve(solver))
        {
            TBOX_ERROR(name << ": the level solve did not converge to the exact solution.\n");
        }
        std::vector<Vec> retained = solver.retainShellVectors();
        // The packed right-hand sides and solutions, and for the multiplicative shell also the packed residuals,
        // an empty vector, a view of the right-hand side, residual and solution of each subdomain, and a vector
        // of the halo of each stage, of which there is one for each subdomain.
        const size_t expected_vectors = 2 + (multiplicative ? 2 + 4 * overlap->size() : 0);
        if (retained.size() != expected_vectors)
        {
            TBOX_ERROR(name << ": expected " << expected_vectors << " shell work vectors, found " << retained.size()
                            << ".\n");
        }
        solver.deallocateSolverState();
        if (!solver.shellStorageEmpty())
        {
            TBOX_ERROR(name << ": shell storage remains after deallocation.\n");
        }
        for (Vec& v : retained)
        {
            PetscInt references;
            PetscErrorCode ierr = PetscObjectGetReference(reinterpret_cast<PetscObject>(v), &references);
            IBTK_CHKERRQ(ierr);
            if (references != 1)
            {
                TBOX_ERROR(name << ": a shell work vector has " << references << " references after deallocation.\n");
            }
            ierr = VecDestroy(&v);
            IBTK_CHKERRQ(ierr);
        }
    }
    return overlap_total;
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app = new AppInitializer(argc, argv, "output");
    LevelFixture fixture(app->getComponentDatabase("CartesianGeometry"));
    Pointer<HierarchyVector> cc_x = new HierarchyVector("cc_x", fixture.hierarchy, 0, 0);
    Pointer<HierarchyVector> sc_x = new HierarchyVector("sc_x", fixture.hierarchy, 0, 0);
    cc_x->addComponent(fixture.x->getComponentVariable(1), fixture.x->getComponentDescriptorIndex(1));
    sc_x->addComponent(fixture.x->getComponentVariable(0), fixture.x->getComponentDescriptorIndex(0));
    Pointer<HierarchyVector> cc_b = cc_x->cloneVector("cc_b"), sc_b = sc_x->cloneVector("sc_b");
    cc_b->allocateVectorData();
    sc_b->allocateVectorData();
    PoissonSpecifications coefs("state_coefs");
    coefs.setCConstant(2.0);
    coefs.setDConstant(0.0);
    PoissonSpecifications stokes_coefs("state_stokes_coefs");
    stokes_coefs.setCConstant(2.0);
    stokes_coefs.setDConstant(-0.01);
    PetscInt previous[3] = { 0, 0, 0 };
    for (const int width : { 0, 2 })
    {
        Pointer<Database> db = level_solver_database("shell", width);
        // Exercise both shell compositions across these lifetimes.
        db->putString("shell_pc_type", width == 0 ? "multiplicative" : "additive");
        LevelSolverProbe<CCPoissonPETScLevelSolver> cc("state_cc", db);
        LevelSolverProbe<SCPoissonPETScLevelSolver> sc("state_sc", db);
        LevelSolverProbe<StaggeredStokesPETScLevelSolver> stokes("state_stokes", db);
        cc.setPoissonSpecifications(coefs);
        sc.setPoissonSpecifications(coefs);
        cc.setPhysicalBcCoef(nullptr);
        sc.setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, nullptr));
        stokes.setVelocityPoissonSpecifications(stokes_coefs);
        stokes.setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, nullptr), nullptr);
        cc.setTimeInterval(0.0, 1.0);
        sc.setTimeInterval(0.0, 1.0);
        stokes.setTimeInterval(0.0, 1.0);
        stokes.setSolutionTime(1.0);
        const PetscInt totals[3] = { check_shell_state(cc, cc_x, cc_b, "Poisson (cell)", false, width == 0),
                                     check_shell_state(sc, sc_x, sc_b, "Poisson (side)", false, width == 0),
                                     check_shell_state(stokes, fixture.x, fixture.b, "Stokes", true, width == 0) };
        const char* names[3] = { "cell", "side", "Stokes" };
        for (int k = 0; k < 3; ++k)
        {
            if (width == 2 && totals[k] <= previous[k])
            {
                TBOX_ERROR(names[k] << ": overlap did not enlarge the subdomains.\n");
            }
            pout << names[k] << " subdomain size total, overlap " << width << " = " << totals[k] << '\n';
            previous[k] = totals[k];
        }
    }
    free_vector_components(*cc_b);
    free_vector_components(*sc_b);
    return 0;
}
