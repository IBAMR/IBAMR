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

#include <ibamr/private/StaggeredStokesEigenSchurComplementShellBackend.h>

#include <tbox/Utilities.h>

#include <cmath>

namespace IBAMR
{
namespace
{
Eigen::MatrixXd
extract_block(const Eigen::MatrixXd& matrix,
              const std::vector<Eigen::Index>& rows,
              const std::vector<Eigen::Index>& columns)
{
    Eigen::MatrixXd block(rows.size(), columns.size());
    for (std::size_t j = 0; j < columns.size(); ++j)
    {
        for (std::size_t i = 0; i < rows.size(); ++i)
        {
            block(i, j) = matrix(rows[i], columns[j]);
        }
    }
    return block;
}
} // namespace
StaggeredStokesEigenSchurComplementShellBackend::StaggeredStokesEigenSchurComplementShellBackend(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
{
    if (input_db)
    {
        d_a00_type = input_db->getStringWithDefault("a00_solver_type", d_a00_type);
        d_schur_type = input_db->getStringWithDefault("schur_solver_type", d_schur_type);
        d_a00_threshold = input_db->getDoubleWithDefault("a00_solver_threshold", -1.0);
        d_schur_threshold = input_db->getDoubleWithDefault("schur_solver_threshold", -1.0);
    }
    IBTK::validate_eigen_local_solver_type(d_a00_type);
    IBTK::validate_eigen_local_solver_type(d_schur_type);
    if (!std::isfinite(d_a00_threshold) || !std::isfinite(d_schur_threshold))
    {
        TBOX_ERROR("Eigen Schur solver thresholds must be finite.\n");
    }
}
StaggeredStokesEigenSchurComplementShellBackend::~StaggeredStokesEigenSchurComplementShellBackend()
{
    deallocateSolverState();
}
void
StaggeredStokesEigenSchurComplementShellBackend::initializeSolverState(
    Mat /*mat*/,
    Vec /*x*/,
    Vec /*b*/,
    const std::vector<IS>& /*overlap*/,
    const std::vector<IS>& /*nonoverlap*/,
    const std::string& /*options_prefix*/,
    bool /*use_multiplicative*/,
    IBTK::PETScLevelSolverShellTraversal /*traversal*/)
{
    TBOX_ERROR("Eigen Schur shell backend requires Stokes velocity and pressure field IDs.\n");
}
void
StaggeredStokesEigenSchurComplementShellBackend::initializeSolverState(
    Mat mat,
    Vec x,
    Vec b,
    const std::vector<IS>& overlap,
    const std::vector<IS>& nonoverlap,
    const std::set<int>& velocity_dofs,
    const std::set<int>& pressure_dofs,
    const std::string& /*options_prefix*/,
    const bool use_multiplicative,
    const IBTK::PETScLevelSolverShellTraversal traversal)
{
    deallocateSolverState();
    initializeSubdomains(mat, x, b, overlap, nonoverlap, use_multiplicative, traversal);
    d_blocks.resize(overlap.size());
    for (std::size_t i = 0; i < overlap.size(); ++i)
    {
        for (std::size_t j = 0; j < d_subdomains[i].dofs.size(); ++j)
        {
            const bool velocity = velocity_dofs.count(d_subdomains[i].dofs[j]) != 0;
            const bool pressure = pressure_dofs.count(d_subdomains[i].dofs[j]) != 0;
            if (velocity == pressure)
            {
                TBOX_ERROR("Eigen Schur overlap DOF must belong to exactly one Stokes field.\n");
            }
            if (velocity)
            {
                d_blocks[i].velocity.push_back(j);
            }
            else
            {
                d_blocks[i].pressure.push_back(j);
            }
        }
        const Eigen::MatrixXd matrix = extractLocalOperator(mat, i);
        const Eigen::Index nv = d_blocks[i].velocity.size(), np = d_blocks[i].pressure.size();
        d_blocks[i].velocity_rhs.resize(nv);
        d_blocks[i].velocity_solution.resize(nv);
        d_blocks[i].pressure_rhs.resize(np);
        d_blocks[i].pressure_solution.resize(np);
        if (nv > 0)
        {
            d_blocks[i].a00_solver = IBTK::make_eigen_local_solver(
                d_a00_type, extract_block(matrix, d_blocks[i].velocity, d_blocks[i].velocity), d_a00_threshold);
        }
        if (np > 0)
        {
            Eigen::MatrixXd schur = extract_block(matrix, d_blocks[i].pressure, d_blocks[i].pressure);
            if (nv > 0)
            {
                d_blocks[i].a10 = extract_block(matrix, d_blocks[i].pressure, d_blocks[i].velocity);
                d_blocks[i].a00_inv_a01 = d_blocks[i].a00_solver->solveMatrix(
                    extract_block(matrix, d_blocks[i].velocity, d_blocks[i].pressure));
                schur.noalias() -= d_blocks[i].a10 * d_blocks[i].a00_inv_a01;
            }
            d_blocks[i].schur_solve =
                IBTK::make_eigen_local_solver(d_schur_type, schur, d_schur_threshold)->getSolveMatrix();
        }
    }
    finalizeComposition();
}
void
StaggeredStokesEigenSchurComplementShellBackend::deallocateSolverState()
{
    clearSubdomains();
    d_blocks.clear();
}
void
StaggeredStokesEigenSchurComplementShellBackend::solveSubdomain(const std::size_t i)
{
    for (std::size_t j = 0; j < d_blocks[i].velocity.size(); ++j)
    {
        d_blocks[i].velocity_rhs[j] = d_subdomains[i].rhs[d_blocks[i].velocity[j]];
    }
    for (std::size_t j = 0; j < d_blocks[i].pressure.size(); ++j)
    {
        d_blocks[i].pressure_rhs[j] = d_subdomains[i].rhs[d_blocks[i].pressure[j]];
    }
    if (!d_blocks[i].velocity.empty())
    {
        d_blocks[i].a00_solver->solve(d_blocks[i].velocity_solution, d_blocks[i].velocity_rhs);
    }
    if (!d_blocks[i].pressure.empty())
    {
        if (!d_blocks[i].velocity.empty())
        {
            d_blocks[i].pressure_rhs.noalias() -= d_blocks[i].a10 * d_blocks[i].velocity_solution;
        }
        d_blocks[i].pressure_solution.noalias() = d_blocks[i].schur_solve * d_blocks[i].pressure_rhs;
        if (!d_blocks[i].velocity.empty())
        {
            d_blocks[i].velocity_solution.noalias() -= d_blocks[i].a00_inv_a01 * d_blocks[i].pressure_solution;
        }
    }
    for (std::size_t j = 0; j < d_blocks[i].velocity.size(); ++j)
    {
        d_subdomains[i].solution[d_blocks[i].velocity[j]] = d_blocks[i].velocity_solution[j];
    }
    for (std::size_t j = 0; j < d_blocks[i].pressure.size(); ++j)
    {
        d_subdomains[i].solution[d_blocks[i].pressure[j]] = d_blocks[i].pressure_solution[j];
    }
}
} // namespace IBAMR
