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

#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScLevelSolverSubdomainSolver.h>

#include <limits>

namespace IBTK
{
namespace
{
class PETScLevelSolverPetscSubdomainSolver
{
public:
    PETScLevelSolverPetscSubdomainSolver() = default;
    PETScLevelSolverPetscSubdomainSolver(const PETScLevelSolverPetscSubdomainSolver&) = delete;
    PETScLevelSolverPetscSubdomainSolver& operator=(const PETScLevelSolverPetscSubdomainSolver&) = delete;
    ~PETScLevelSolverPetscSubdomainSolver();
    void initializeSolverState(const std::vector<Mat>& matrices,
                               const std::vector<IS>& subdomains,
                               const std::string& options_prefix);
    void deallocateSolverState();
    void solve(std::size_t first, std::size_t last, Vec b, Vec x);

private:
    std::vector<KSP> d_sub_ksp;
    // Vectors that are placed on the packed arrays of one subdomain at a time, so that solving does not copy.
    std::vector<Vec> d_rhs_views, d_solution_views;
    std::vector<PetscInt> d_offsets;
};

PETScLevelSolverPetscSubdomainSolver::~PETScLevelSolverPetscSubdomainSolver()
{
    deallocateSolverState();
}

void
PETScLevelSolverPetscSubdomainSolver::initializeSolverState(const std::vector<Mat>& matrices,
                                                            const std::vector<IS>& /*subdomains*/,
                                                            const std::string& options_prefix)
{
    deallocateSolverState();
    d_sub_ksp.resize(matrices.size());
    d_rhs_views.resize(matrices.size());
    d_solution_views.resize(matrices.size());
    d_offsets.assign(matrices.size() + 1, 0);
    for (std::size_t i = 0; i < matrices.size(); ++i)
    {
        PetscInt order = 0;
        int ierr = MatGetSize(matrices[i], &order, nullptr);
        IBTK_CHKERRQ(ierr);
        d_offsets[i + 1] = d_offsets[i] + order;
        for (Vec* view : { &d_rhs_views[i], &d_solution_views[i] })
        {
            ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, order, nullptr, view);
            IBTK_CHKERRQ(ierr);
        }
        KSP& sub_ksp = d_sub_ksp[i];
        ierr = KSPCreate(PETSC_COMM_SELF, &sub_ksp);
        IBTK_CHKERRQ(ierr);
        const std::string sub_prefix = options_prefix + "_sub";
        ierr = KSPSetOptionsPrefix(sub_ksp, sub_prefix.c_str());
        IBTK_CHKERRQ(ierr);
        ierr = KSPSetOperators(sub_ksp, matrices[i], matrices[i]);
        IBTK_CHKERRQ(ierr);

        // Set default configuration.
        ierr = KSPSetReusePreconditioner(sub_ksp, PETSC_TRUE);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSetType(sub_ksp, KSPPREONLY);
        IBTK_CHKERRQ(ierr);
        PC sub_pc;
        ierr = KSPGetPC(sub_ksp, &sub_pc);
        IBTK_CHKERRQ(ierr);
        ierr = PCSetType(sub_pc, PCLU);
        IBTK_CHKERRQ(ierr);
        ierr = PCFactorReorderForNonzeroDiagonal(sub_pc, std::numeric_limits<double>::epsilon());
        IBTK_CHKERRQ(ierr);

        // Set from options.
        ierr = KSPSetFromOptions(sub_ksp);
        IBTK_CHKERRQ(ierr);

        // Always use a zero initial guess.
        ierr = KSPSetInitialGuessNonzero(sub_ksp, PETSC_FALSE);
        IBTK_CHKERRQ(ierr);

        // Factor now. PETSc records a failed factorization in the KSP instead of returning an error, so check for it
        // here, as the other subdomain solvers report a singular block at initialization.
        ierr = KSPSetUp(sub_ksp);
        IBTK_CHKERRQ(ierr);
        KSPConvergedReason reason;
        ierr = KSPGetConvergedReason(sub_ksp, &reason);
        IBTK_CHKERRQ(ierr);
        if (reason == KSP_DIVERGED_PC_FAILED)
        {
            TBOX_ERROR("PETScLevelSolverPetscSubdomainSolver::initializeSolverState():\n"
                       << "  options prefix \"" << options_prefix
                       << "\": setting up the preconditioner failed for subdomain " << i
                       << ", for example because the subdomain matrix is singular.\n");
        }
    }
}

void
PETScLevelSolverPetscSubdomainSolver::deallocateSolverState()
{
    for (KSP& sub_ksp : d_sub_ksp)
    {
        const int ierr = KSPDestroy(&sub_ksp);
        IBTK_CHKERRQ(ierr);
    }
    d_sub_ksp.clear();
    for (std::vector<Vec>* views : { &d_rhs_views, &d_solution_views })
    {
        for (Vec& view : *views)
        {
            const int ierr = VecDestroy(&view);
            IBTK_CHKERRQ(ierr);
        }
        views->clear();
    }
    d_offsets.clear();
}

void
PETScLevelSolverPetscSubdomainSolver::solve(const std::size_t first, const std::size_t last, Vec b, Vec x)
{
    const PetscScalar* rhs = nullptr;
    PetscScalar* solution = nullptr;
    int ierr = VecGetArrayRead(b, &rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecGetArray(x, &solution);
    IBTK_CHKERRQ(ierr);
    for (std::size_t i = first; i < last; ++i)
    {
        ierr = VecPlaceArray(d_rhs_views[i], rhs + d_offsets[i]);
        IBTK_CHKERRQ(ierr);
        ierr = VecPlaceArray(d_solution_views[i], solution + d_offsets[i]);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSolve(d_sub_ksp[i], d_rhs_views[i], d_solution_views[i]);
        IBTK_CHKERRQ(ierr);
        ierr = VecResetArray(d_rhs_views[i]);
        IBTK_CHKERRQ(ierr);
        ierr = VecResetArray(d_solution_views[i]);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecRestoreArray(x, &solution);
    IBTK_CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(b, &rhs);
    IBTK_CHKERRQ(ierr);
}
} // namespace

PETScLevelSolverSubdomainSolver
make_petsc_subdomain_solver()
{
    return PETScLevelSolverSubdomainSolver(std::in_place_type<PETScLevelSolverPetscSubdomainSolver>);
}
} // namespace IBTK
