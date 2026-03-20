// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2025 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/IBTK_MPI.h>
#include <ibtk/PETScLevelSolver.h>
#include <ibtk/private/PETScLevelSolverEigenPseudoinverseShellBackend.h>

#include <tbox/Database.h>

namespace IBTK
{
PETScLevelSolverEigenPseudoinverseShellBackend::PETScLevelSolverEigenPseudoinverseShellBackend(PETScLevelSolver& solver)
    : PETScLevelSolverEigenShellBackendBase(solver)
{
}

void
PETScLevelSolverEigenPseudoinverseShellBackend::configure(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
{
    d_solver_type = EigenSubdomainSolverType::COL_PIV_HOUSEHOLDER_QR;
    d_solver_threshold = -1.0;
    if (!input_db) return;

    if (input_db->keyExists("eigen_subdomain_pseudoinverse_type"))
    {
        d_solver_type = parseEigenSubdomainSolverType(input_db->getString("eigen_subdomain_pseudoinverse_type"),
                                                      "PETScLevelSolverEigenPseudoinverseShellBackend::configure()");
    }
    if (input_db->keyExists("eigen_subdomain_pseudoinverse_threshold"))
    {
        d_solver_threshold = input_db->getDouble("eigen_subdomain_pseudoinverse_threshold");
    }
}

void
PETScLevelSolverEigenPseudoinverseShellBackend::initialize()
{
    initializeCommonDataWithLocalOperatorHook(
        [this](const Eigen::MatrixXd& local_operator, const std::size_t subdomain_num)
        {
            auto& cache = getCommonSubdomains()[subdomain_num];
            cache.local_pseudoinverse = buildSubdomainPseudoinverse(local_operator);
        });
}

void
PETScLevelSolverEigenPseudoinverseShellBackend::deallocate()
{
    clearCommonData();
}

void
PETScLevelSolverEigenPseudoinverseShellBackend::initializeSubdomainPseudoinverse(const Eigen::MatrixXd& local_operator,
                                                                                 const std::size_t subdomain_num)
{
    getCommonSubdomains()[subdomain_num].local_pseudoinverse = buildSubdomainPseudoinverse(local_operator);
}

Eigen::VectorXd
PETScLevelSolverEigenPseudoinverseShellBackend::solveSubdomainSystem(const Eigen::VectorXd& rhs,
                                                                     const std::size_t subdomain_num) const
{
    return getCommonSubdomains()[subdomain_num].local_pseudoinverse * rhs;
}

Eigen::MatrixXd
PETScLevelSolverEigenPseudoinverseShellBackend::buildSubdomainPseudoinverse(const Eigen::MatrixXd& local_operator) const
{
    Eigen::MatrixXd pseudoinverse;
    dispatchEigenSolverType(
        getSolverType(),
        [this, &local_operator, &pseudoinverse](auto solver_tag)
        {
            using SolverType = typename decltype(solver_tag)::type;
            if constexpr (std::is_same_v<SolverType, Eigen::LLT<Eigen::MatrixXd>>)
            {
                pseudoinverse = buildLLTSolveMatrix(local_operator);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::LDLT<Eigen::MatrixXd>>)
            {
                pseudoinverse = buildLDLTSolveMatrix(local_operator);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::PartialPivLU<Eigen::MatrixXd>>)
            {
                pseudoinverse = buildPartialPivLUSolveMatrix(local_operator);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::FullPivLU<Eigen::MatrixXd>>)
            {
                pseudoinverse = buildFullPivLUSolveMatrix(local_operator, getSolverThreshold());
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::HouseholderQR<Eigen::MatrixXd>>)
            {
                pseudoinverse = buildHouseholderQRSolveMatrix(local_operator);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::ColPivHouseholderQR<Eigen::MatrixXd>>)
            {
                pseudoinverse = buildQRSolveMatrix<SolverType>(local_operator, getSolverThreshold());
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd>>)
            {
                pseudoinverse = buildCompleteOrthogonalDecompositionPseudoinverse(local_operator, getSolverThreshold());
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::FullPivHouseholderQR<Eigen::MatrixXd>>)
            {
                pseudoinverse = buildQRSolveMatrix<SolverType>(local_operator, getSolverThreshold());
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::JacobiSVD<Eigen::MatrixXd>> ||
                               std::is_same_v<SolverType, Eigen::BDCSVD<Eigen::MatrixXd>>)
            {
                pseudoinverse = buildSVDPseudoinverse<SolverType>(local_operator, getSolverThreshold());
            }
        });
    return pseudoinverse;
}

PETScLevelSolverEigenPseudoinverseShellBackend::EigenSubdomainSolverType
PETScLevelSolverEigenPseudoinverseShellBackend::getSolverType() const
{
    TBOX_ASSERT(d_solver.d_shell_smoother_backend == PETScLevelSolver::ShellSmootherBackend::EIGEN_PSEUDOINVERSE);
    return d_solver_type;
}

double
PETScLevelSolverEigenPseudoinverseShellBackend::getSolverThreshold() const
{
    return d_solver_threshold;
}

void
PETScLevelSolverEigenPseudoinverseShellBackend::applyAdditive(Vec x, Vec y)
{
    TBOX_ASSERT(IBTK_MPI::getNodes() == 1);
    auto& common_subdomains = getCommonSubdomains();
    const Eigen::Index n = getNumDofs();
    TBOX_ASSERT(n > 0);

    {
        ConstPetscVecArrayMap x_array(x, n);
        PetscVecArrayMap y_array(y, n);
        const auto x_map = x_array.getMap();
        auto y_map = y_array.getMap();
        y_map.setZero();
        for (std::size_t subdomain_num = 0; subdomain_num < common_subdomains.size(); ++subdomain_num)
        {
            auto& cache = common_subdomains[subdomain_num];
            std::size_t rhs_idx = 0;
            for (const int dof : cache.overlap_dofs)
            {
                cache.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = x_map[dof];
            }
            cache.delta_workspace.noalias() = cache.local_pseudoinverse * cache.rhs_workspace;
            std::size_t update_idx = 0;
            for (const int dof : cache.update_dofs)
            {
                y_map[dof] +=
                    cache.delta_workspace[static_cast<Eigen::Index>(cache.update_local_positions[update_idx++])];
            }
        }
    }
    d_solver.postprocessShellResult(y);
}

void
PETScLevelSolverEigenPseudoinverseShellBackend::applyMultiplicative(Vec x, Vec y)
{
    TBOX_ASSERT(IBTK_MPI::getNodes() == 1);
    auto& common_subdomains = getCommonSubdomains();
    const Eigen::Index n = getNumDofs();
    TBOX_ASSERT(n > 0);

    {
        ConstPetscVecArrayMap x_array(x, n);
        PetscVecArrayMap y_array(y, n);
        const auto x_map = x_array.getMap();
        auto y_map = y_array.getMap();
        Eigen::VectorXd residual(n);
        y_map.setZero();
        residual = x_map;
        const std::size_t n_subdomains = common_subdomains.size();
        for (std::size_t subdomain_num = 0; subdomain_num + 1 < n_subdomains; ++subdomain_num)
        {
            auto& cache = common_subdomains[subdomain_num];
            std::size_t rhs_idx = 0;
            for (const int dof : cache.overlap_dofs)
            {
                cache.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = residual[dof];
            }
            cache.delta_workspace.noalias() = cache.local_pseudoinverse * cache.rhs_workspace;

            std::size_t update_idx = 0;
            for (const int dof : cache.update_dofs)
            {
                y_map[dof] +=
                    cache.delta_workspace[static_cast<Eigen::Index>(cache.update_local_positions[update_idx++])];
            }
            if (cache.active_residual_update_mat.rows() > 0)
            {
                std::size_t residual_input_idx = 0;
                for (const int local_pos : cache.update_local_positions)
                {
                    cache.residual_input_workspace[static_cast<Eigen::Index>(residual_input_idx++)] =
                        cache.delta_workspace[static_cast<Eigen::Index>(local_pos)];
                }
                cache.residual_delta_workspace.noalias() =
                    cache.active_residual_update_mat * cache.residual_input_workspace;
                std::size_t row_idx = 0;
                for (const int row : cache.active_residual_update_rows)
                {
                    residual[row] -= cache.residual_delta_workspace[static_cast<Eigen::Index>(row_idx++)];
                }
            }
        }
        if (n_subdomains > 0)
        {
            auto& cache = common_subdomains.back();
            std::size_t rhs_idx = 0;
            for (const int dof : cache.overlap_dofs)
            {
                cache.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = residual[dof];
            }
            cache.delta_workspace.noalias() = cache.local_pseudoinverse * cache.rhs_workspace;

            std::size_t update_idx = 0;
            for (const int dof : cache.update_dofs)
            {
                y_map[dof] +=
                    cache.delta_workspace[static_cast<Eigen::Index>(cache.update_local_positions[update_idx++])];
            }
        }
    }
    d_solver.postprocessShellResult(y);
}
} // namespace IBTK
