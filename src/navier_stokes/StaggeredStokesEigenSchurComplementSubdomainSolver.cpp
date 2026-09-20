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

#include <ibtk/IBTK_CHKERRQ.h>

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

StaggeredStokesEigenSchurComplementSubdomainSolver::StaggeredStokesEigenSchurComplementSubdomainSolver(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
    FieldProvider fields)
    : d_fields(std::move(fields))
{
    if (input_db)
    {
        d_a00_type = input_db->getStringWithDefault("a00_solver_type", d_a00_type);
        d_schur_type = input_db->getStringWithDefault("schur_solver_type", d_schur_type);
        d_a00_threshold = input_db->getDoubleWithDefault("a00_solver_threshold", -1.0);
        d_schur_threshold = input_db->getDoubleWithDefault("schur_solver_threshold", -1.0);
    }
    IBTK::validate_eigen_solver_type(d_a00_type);
    IBTK::validate_eigen_solver_type(d_schur_type);
    if (!std::isfinite(d_a00_threshold))
    {
        TBOX_ERROR(
            "StaggeredStokesEigenSchurComplementSubdomainSolver::StaggeredStokesEigenSchurComplementSubdomainSolver():"
            "\n"
            << "  a00_solver_threshold must be finite.\n");
    }
    if (!std::isfinite(d_schur_threshold))
    {
        TBOX_ERROR(
            "StaggeredStokesEigenSchurComplementSubdomainSolver::StaggeredStokesEigenSchurComplementSubdomainSolver():"
            "\n"
            << "  schur_solver_threshold must be finite.\n");
    }
}

void
StaggeredStokesEigenSchurComplementSubdomainSolver::initializeSolverState(const std::vector<Mat>& matrices,
                                                                          const std::vector<IS>& subdomains,
                                                                          const std::string& /*options_prefix*/)
{
    deallocateSolverState();
    d_offsets.assign(matrices.size() + 1, 0);
    d_blocks.resize(matrices.size());
    std::vector<PetscInt> gathered_dofs;
    for (std::size_t i = 0; i < matrices.size(); ++i)
    {
        PetscInt order = 0;
        int ierr = ISGetLocalSize(subdomains[i], &order);
        IBTK_CHKERRQ(ierr);
        d_offsets[i + 1] = d_offsets[i] + order;
        const PetscInt* dofs = nullptr;
        ierr = ISGetIndices(subdomains[i], &dofs);
        IBTK_CHKERRQ(ierr);
        gathered_dofs.insert(gathered_dofs.end(), dofs, dofs + order);
        ierr = ISRestoreIndices(subdomains[i], &dofs);
        IBTK_CHKERRQ(ierr);
    }

    // Gather the field of every DOF of every subdomain, which may be owned by another rank.
    Vec fields = d_fields();
    const PetscInt n_gathered = d_offsets.back();
    Vec gathered = nullptr;
    IS gathered_is = nullptr, positions_is = nullptr;
    VecScatter scatter = nullptr;
    int ierr = VecCreateSeq(PETSC_COMM_SELF, n_gathered, &gathered);
    IBTK_CHKERRQ(ierr);
    ierr = ISCreateGeneral(PETSC_COMM_SELF, n_gathered, gathered_dofs.data(), PETSC_COPY_VALUES, &gathered_is);
    IBTK_CHKERRQ(ierr);
    ierr = ISCreateStride(PETSC_COMM_SELF, n_gathered, 0, 1, &positions_is);
    IBTK_CHKERRQ(ierr);
    ierr = VecScatterCreate(fields, gathered_is, gathered, positions_is, &scatter);
    IBTK_CHKERRQ(ierr);
    ierr = VecScatterBegin(scatter, fields, gathered, INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);
    ierr = VecScatterEnd(scatter, fields, gathered, INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);
    const PetscScalar* field_values = nullptr;
    ierr = VecGetArrayRead(gathered, &field_values);
    IBTK_CHKERRQ(ierr);
    for (std::size_t i = 0; i < matrices.size(); ++i)
    {
        Block& block = d_blocks[i];
        for (PetscInt j = 0; j < d_offsets[i + 1] - d_offsets[i]; ++j)
        {
            const PetscReal field = PetscRealPart(field_values[d_offsets[i] + j]);
            if (field != 0.0 && field != 1.0)
            {
                TBOX_ERROR("StaggeredStokesEigenSchurComplementSubdomainSolver::initializeSolverState():\n"
                           << "  every subdomain DOF must belong to exactly one Stokes field.\n");
            }
            (field == 0.0 ? block.velocity : block.pressure).push_back(j);
        }
    }
    ierr = VecRestoreArrayRead(gathered, &field_values);
    IBTK_CHKERRQ(ierr);
    ierr = VecScatterDestroy(&scatter);
    IBTK_CHKERRQ(ierr);
    ierr = ISDestroy(&gathered_is);
    IBTK_CHKERRQ(ierr);
    ierr = ISDestroy(&positions_is);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&gathered);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&fields);
    IBTK_CHKERRQ(ierr);
    for (std::size_t i = 0; i < matrices.size(); ++i)
    {
        Block& block = d_blocks[i];
        const Eigen::MatrixXd matrix = IBTK::make_eigen_matrix(matrices[i]);
        const Eigen::Index nv = block.velocity.size(), np = block.pressure.size();
        block.velocity_rhs.resize(nv);
        block.velocity_solution.resize(nv);
        block.pressure_rhs.resize(np);
        block.pressure_solution.resize(np);
        std::unique_ptr<IBTK::EigenDenseSolver> a00_solver;
        if (nv > 0)
        {
            a00_solver = IBTK::make_eigen_dense_solver(
                d_a00_type, extract_block(matrix, block.velocity, block.velocity), d_a00_threshold);
            block.a00_solve = a00_solver->getSolveMatrix();
        }
        if (np > 0)
        {
            Eigen::MatrixXd schur = extract_block(matrix, block.pressure, block.pressure);
            if (nv > 0)
            {
                block.a10 = extract_block(matrix, block.pressure, block.velocity);
                block.a00_inv_a01 = a00_solver->solveMatrix(extract_block(matrix, block.velocity, block.pressure));
                schur.noalias() -= block.a10 * block.a00_inv_a01;
            }
            block.schur_solve = IBTK::make_eigen_dense_solver(d_schur_type, schur, d_schur_threshold)->getSolveMatrix();
        }
    }
}

void
StaggeredStokesEigenSchurComplementSubdomainSolver::deallocateSolverState()
{
    d_offsets.clear();
    d_blocks.clear();
}

void
StaggeredStokesEigenSchurComplementSubdomainSolver::solve(const std::size_t first, const std::size_t last, Vec b, Vec x)
{
    const PetscScalar* rhs = nullptr;
    PetscScalar* solution = nullptr;
    int ierr = VecGetArrayRead(b, &rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecGetArray(x, &solution);
    IBTK_CHKERRQ(ierr);
    for (std::size_t i = first; i < last; ++i)
    {
        Block& block = d_blocks[i];
        const PetscScalar* subdomain_rhs = rhs + d_offsets[i];
        PetscScalar* subdomain_solution = solution + d_offsets[i];
        for (std::size_t j = 0; j < block.velocity.size(); ++j)
        {
            block.velocity_rhs[j] = subdomain_rhs[block.velocity[j]];
        }
        for (std::size_t j = 0; j < block.pressure.size(); ++j)
        {
            block.pressure_rhs[j] = subdomain_rhs[block.pressure[j]];
        }
        if (!block.velocity.empty())
        {
            block.velocity_solution.noalias() = block.a00_solve * block.velocity_rhs;
        }
        if (!block.pressure.empty())
        {
            if (!block.velocity.empty())
            {
                block.pressure_rhs.noalias() -= block.a10 * block.velocity_solution;
            }
            block.pressure_solution.noalias() = block.schur_solve * block.pressure_rhs;
            if (!block.velocity.empty())
            {
                block.velocity_solution.noalias() -= block.a00_inv_a01 * block.pressure_solution;
            }
        }
        for (std::size_t j = 0; j < block.velocity.size(); ++j)
        {
            subdomain_solution[block.velocity[j]] = block.velocity_solution[j];
        }
        for (std::size_t j = 0; j < block.pressure.size(); ++j)
        {
            subdomain_solution[block.pressure[j]] = block.pressure_solution[j];
        }
    }
    ierr = VecRestoreArray(x, &solution);
    IBTK_CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(b, &rhs);
    IBTK_CHKERRQ(ierr);
}
} // namespace IBAMR
