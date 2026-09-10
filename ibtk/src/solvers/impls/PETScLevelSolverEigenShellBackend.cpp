// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/private/PETScLevelSolverEigenShellBackend.h>

#include <tbox/Utilities.h>

#include <cmath>

namespace IBTK
{
PETScLevelSolverEigenShellBackend::PETScLevelSolverEigenShellBackend(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
    const bool precompute)
    : d_precompute(precompute), d_solver_type(precompute ? "COL_PIV_HOUSEHOLDER_QR" : "PARTIAL_PIV_LU")
{
    const std::string key = precompute ? "eigen_subdomain_pseudoinverse" : "eigen_subdomain_solver";
    if (input_db)
    {
        d_solver_type = input_db->getStringWithDefault(key + "_type", d_solver_type);
        d_threshold = input_db->getDoubleWithDefault(key + "_threshold", -1.0);
    }
    validate_eigen_local_solver_type(d_solver_type);
    if (!std::isfinite(d_threshold))
    {
        TBOX_ERROR("Eigen local solver threshold must be finite.\n");
    }
}
PETScLevelSolverEigenShellBackend::~PETScLevelSolverEigenShellBackend()
{
    deallocateSolverState();
}
void
PETScLevelSolverEigenShellBackend::initializeSolverState(Mat mat,
                                                         Vec x,
                                                         Vec b,
                                                         const std::vector<IS>& overlap,
                                                         const std::vector<IS>& nonoverlap,
                                                         const std::string& /*options_prefix*/,
                                                         const bool use_multiplicative,
                                                         const PETScLevelSolverShellTraversal traversal)
{
    deallocateSolverState();
    initializeSubdomains(mat, x, b, overlap, nonoverlap, use_multiplicative, traversal);
    if (d_precompute)
    {
        d_solve_matrices.resize(overlap.size());
    }
    else
    {
        d_solvers.resize(overlap.size());
    }
    for (std::size_t i = 0; i < overlap.size(); ++i)
    {
        if (d_subdomains[i].dofs.empty())
        {
            continue;
        }
        std::unique_ptr<EigenLocalSolver> solver =
            make_eigen_local_solver(d_solver_type, extractLocalOperator(mat, i), d_threshold);
        if (d_precompute)
        {
            d_solve_matrices[i] = solver->getSolveMatrix();
        }
        else
        {
            d_solvers[i] = std::move(solver);
        }
    }
    finalizeComposition();
}
void
PETScLevelSolverEigenShellBackend::deallocateSolverState()
{
    clearSubdomains();
    d_solvers.clear();
    d_solve_matrices.clear();
}
void
PETScLevelSolverEigenShellBackend::solveSubdomain(const std::size_t i)
{
    if (d_subdomains[i].dofs.empty())
    {
        return;
    }
    if (d_precompute)
    {
        d_subdomains[i].solution.noalias() = d_solve_matrices[i] * d_subdomains[i].rhs;
    }
    else
    {
        d_solvers[i]->solve(d_subdomains[i].solution, d_subdomains[i].rhs);
    }
}
} // namespace IBTK
