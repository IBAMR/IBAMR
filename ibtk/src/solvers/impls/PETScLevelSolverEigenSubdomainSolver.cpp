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

#include <ibtk/EigenDenseSolver.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScLevelSolverSubdomainSolver.h>

#include <tbox/Database.h>
#include <tbox/Utilities.h>

#include <Eigen/Core>

#include <cmath>

namespace IBTK
{
namespace
{
/*! \brief Eigen factorizations, or precomputed solve matrices, for PETScLevelSolverSubdomainSolver.
 *
 * Solves read and write the packed vectors of the level solver through Eigen maps, without copies. Types
 * whose Eigen solve allocates a temporary vector form their solve matrix at setup, so that applying the
 * solver never allocates.
 */
class PETScLevelSolverEigenSubdomainSolver
{
public:
    /*! \brief Read the settings, using defaults for a null database. */
    PETScLevelSolverEigenSubdomainSolver(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db, bool precompute);
    void initializeSolverState(const std::vector<Mat>& matrices,
                               const std::vector<IS>& subdomains,
                               const std::string& options_prefix);
    void deallocateSolverState();
    void solve(std::size_t first, std::size_t last, Vec b, Vec x);

private:
    const bool d_precompute;
    std::string d_solver_type;
    double d_threshold = -1.0;
    std::vector<PetscInt> d_offsets;
    std::vector<std::unique_ptr<EigenDenseSolver>> d_solvers;
    std::vector<Eigen::MatrixXd> d_solve_matrices;
};

PETScLevelSolverEigenSubdomainSolver::PETScLevelSolverEigenSubdomainSolver(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
    const bool precompute)
    : d_precompute(precompute), d_solver_type("COL_PIV_HOUSEHOLDER_QR")
{
    const std::string key = precompute ? "eigen_subdomain_pseudoinverse" : "eigen_subdomain_solver";
    if (input_db)
    {
        d_solver_type = input_db->getStringWithDefault(key + "_type", d_solver_type);
        d_threshold = input_db->getDoubleWithDefault(key + "_threshold", -1.0);
    }
    validate_eigen_solver_type(d_solver_type);
    if (!std::isfinite(d_threshold))
    {
        TBOX_ERROR("PETScLevelSolverEigenSubdomainSolver::PETScLevelSolverEigenSubdomainSolver():\n"
                   << "  " << key << "_threshold must be finite.\n");
    }
}

void
PETScLevelSolverEigenSubdomainSolver::initializeSolverState(const std::vector<Mat>& matrices,
                                                            const std::vector<IS>& /*subdomains*/,
                                                            const std::string& /*options_prefix*/)
{
    deallocateSolverState();
    d_offsets.assign(matrices.size() + 1, 0);
    d_solve_matrices.resize(matrices.size());
    d_solvers.resize(matrices.size());
    for (std::size_t i = 0; i < matrices.size(); ++i)
    {
        PetscInt order = 0;
        const int ierr = MatGetSize(matrices[i], &order, nullptr);
        IBTK_CHKERRQ(ierr);
        d_offsets[i + 1] = d_offsets[i] + order;
        if (order == 0)
        {
            continue;
        }
        std::unique_ptr<EigenDenseSolver> solver =
            make_eigen_dense_solver(d_solver_type, make_eigen_matrix(matrices[i]), d_threshold);
        if (d_precompute || solver->solveAllocates())
        {
            d_solve_matrices[i] = solver->getSolveMatrix();
        }
        else
        {
            d_solvers[i] = std::move(solver);
        }
    }
}

void
PETScLevelSolverEigenSubdomainSolver::deallocateSolverState()
{
    d_offsets.clear();
    d_solvers.clear();
    d_solve_matrices.clear();
}

void
PETScLevelSolverEigenSubdomainSolver::solve(const std::size_t first, const std::size_t last, Vec b, Vec x)
{
    const PetscScalar* rhs = nullptr;
    PetscScalar* solution = nullptr;
    int ierr = VecGetArrayRead(b, &rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecGetArray(x, &solution);
    IBTK_CHKERRQ(ierr);
    for (std::size_t i = first; i < last; ++i)
    {
        const Eigen::Index order = d_offsets[i + 1] - d_offsets[i];
        if (order == 0)
        {
            continue;
        }
        if (d_solvers[i])
        {
            d_solvers[i]->solve(rhs + d_offsets[i], solution + d_offsets[i]);
        }
        else
        {
            Eigen::Map<Eigen::VectorXd>(solution + d_offsets[i], order).noalias() =
                d_solve_matrices[i] * Eigen::Map<const Eigen::VectorXd>(rhs + d_offsets[i], order);
        }
    }
    ierr = VecRestoreArray(x, &solution);
    IBTK_CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(b, &rhs);
    IBTK_CHKERRQ(ierr);
}
} // namespace

PETScLevelSolverSubdomainSolver
make_eigen_subdomain_solver(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
{
    return PETScLevelSolverSubdomainSolver(std::in_place_type<PETScLevelSolverEigenSubdomainSolver>, input_db, false);
}

PETScLevelSolverSubdomainSolver
make_eigen_pseudoinverse_subdomain_solver(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
{
    return PETScLevelSolverSubdomainSolver(std::in_place_type<PETScLevelSolverEigenSubdomainSolver>, input_db, true);
}
} // namespace IBTK
