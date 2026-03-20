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

#ifndef included_IBTK_private_PETScLevelSolverEigenShellBackendBase_inl
#define included_IBTK_private_PETScLevelSolverEigenShellBackendBase_inl

namespace IBTK
{
template <class SolverType>
inline void
PETScLevelSolverEigenShellBackendBase::initializeEigenSolver(SolverType& solver,
                                                             const Eigen::MatrixXd& local_operator,
                                                             const double threshold)
{
    if constexpr (eigen_solver_uses_svd_compute<SolverType>::value)
    {
        solver.compute(local_operator, Eigen::ComputeThinU | Eigen::ComputeThinV);
    }
    else
    {
        solver.compute(local_operator);
    }

    if constexpr (eigen_solver_supports_threshold<SolverType>::value)
    {
        if (threshold >= 0.0) solver.setThreshold(threshold);
    }
}

template <class SolverType, class RhsType>
inline auto
PETScLevelSolverEigenShellBackendBase::solveEigenSystem(const SolverType& solver, const RhsType& rhs)
    -> decltype(solver.solve(rhs))
{
    return solver.solve(rhs);
}

template <class SolverType>
inline Eigen::MatrixXd
PETScLevelSolverEigenShellBackendBase::buildQRSolveMatrix(const Eigen::MatrixXd& matrix, const double threshold)
{
    SolverType solver(matrix);
    if (threshold >= 0.0) solver.setThreshold(threshold);
    return solver.solve(Eigen::MatrixXd::Identity(matrix.rows(), matrix.cols()));
}

inline Eigen::SparseMatrix<double, Eigen::RowMajor>
PETScLevelSolverEigenShellBackendBase::copyPETScMatToEigenSparse(Mat mat)
{
    PetscInt m = 0, n = 0;
    int ierr = MatGetSize(mat, &m, &n);
    IBTK_CHKERRQ(ierr);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<std::size_t>(m) * 16);
    for (PetscInt row = 0; row < m; ++row)
    {
        PetscInt row_nnz = 0;
        const PetscInt* cols = nullptr;
        const PetscScalar* vals = nullptr;
        ierr = MatGetRow(mat, row, &row_nnz, &cols, &vals);
        IBTK_CHKERRQ(ierr);
        for (PetscInt k = 0; k < row_nnz; ++k)
        {
            triplets.emplace_back(static_cast<int>(row), static_cast<int>(cols[k]), PetscRealPart(vals[k]));
        }
        ierr = MatRestoreRow(mat, row, &row_nnz, &cols, &vals);
        IBTK_CHKERRQ(ierr);
    }
    Eigen::SparseMatrix<double, Eigen::RowMajor> A(static_cast<int>(m), static_cast<int>(n));
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();
    return A;
}

template <class Handler>
inline void
PETScLevelSolverEigenShellBackendBase::dispatchEigenSolverType(const EigenSubdomainSolverType solver_type,
                                                               Handler&& handler) const
{
    PETScLevelSolverEigenShell::dispatchSubdomainSolverType(
        d_solver.d_object_name, d_solver.d_options_prefix, solver_type, std::forward<Handler>(handler));
}

inline PETScLevelSolverEigenShellBackendBase::EigenSubdomainSolverType
PETScLevelSolverEigenShellBackendBase::parseEigenSubdomainSolverType(const std::string& type, const char* caller) const
{
    return PETScLevelSolverEigenShell::parseSubdomainSolverType(
        type, d_solver.d_object_name, d_solver.d_options_prefix, caller);
}

template <class SVDType>
inline Eigen::MatrixXd
PETScLevelSolverEigenShellBackendBase::buildSVDPseudoinverse(const Eigen::MatrixXd& matrix, const double threshold)
{
    SVDType svd(matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
    if (threshold >= 0.0) svd.setThreshold(threshold);

    const auto& singular_values = svd.singularValues();
    Eigen::VectorXd inv_singular_values = Eigen::VectorXd::Zero(singular_values.size());
    if (singular_values.size() > 0)
    {
        const double cutoff = singular_values[0] * svd.threshold();
        for (Eigen::Index k = 0; k < singular_values.size(); ++k)
        {
            if (singular_values[k] > cutoff) inv_singular_values[k] = 1.0 / singular_values[k];
        }
    }

    return svd.matrixV() * inv_singular_values.asDiagonal() * svd.matrixU().adjoint();
}

template <class InitializeSubdomainSolver>
inline void
PETScLevelSolverEigenShellBackendBase::initializeCommonDataWithLocalOperatorHook(
    InitializeSubdomainSolver initialize_subdomain_solver)
{
    const bool use_multiplicative =
        d_solver.d_shell_smoother_composition == PETScLevelSolver::ShellSmootherComposition::MULTIPLICATIVE;
    const Eigen::SparseMatrix<double, Eigen::RowMajor> eigen_level_mat =
        copyPETScMatToEigenSparse(d_solver.d_petsc_mat);
    Eigen::SparseMatrix<double> eigen_level_mat_transpose;
    if (use_multiplicative) eigen_level_mat_transpose = eigen_level_mat.transpose();
    d_n_dofs = eigen_level_mat.rows();
    d_common_subdomains.clear();
    d_common_subdomains.resize(d_solver.d_subdomain_dofs.size());
    for (std::size_t subdomain_num = 0; subdomain_num < d_common_subdomains.size(); ++subdomain_num)
    {
        auto& cache = d_common_subdomains[subdomain_num];
        cache = CommonSubdomainCache();
        cache.overlap_dofs = d_solver.d_subdomain_dofs[subdomain_num];
        const Eigen::Index overlap_size = static_cast<Eigen::Index>(cache.overlap_dofs.size());
        std::unordered_map<int, int> overlap_col_map;
        overlap_col_map.reserve(static_cast<std::size_t>(overlap_size));
        for (Eigen::Index local_col = 0; local_col < overlap_size; ++local_col)
        {
            overlap_col_map[cache.overlap_dofs[static_cast<std::size_t>(local_col)]] = static_cast<int>(local_col);
        }

        cache.nonoverlap_dofs = d_solver.d_nonoverlap_subdomain_dofs[subdomain_num];
        cache.nonoverlap_local_positions.resize(cache.nonoverlap_dofs.size());
        for (std::size_t local_col = 0; local_col < cache.nonoverlap_dofs.size(); ++local_col)
        {
            const int dof = cache.nonoverlap_dofs[local_col];
            const auto overlap_pos_it = overlap_col_map.find(dof);
#if !defined(NDEBUG)
            TBOX_ASSERT(overlap_pos_it != overlap_col_map.end());
#endif
            cache.nonoverlap_local_positions[local_col] = overlap_pos_it->second;
        }

        Eigen::MatrixXd local_operator = Eigen::MatrixXd::Zero(overlap_size, overlap_size);
        for (Eigen::Index local_row = 0; local_row < overlap_size; ++local_row)
        {
            const int global_row = cache.overlap_dofs[static_cast<std::size_t>(local_row)];
            for (auto it = Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator(eigen_level_mat, global_row); it;
                 ++it)
            {
                const auto col_it = overlap_col_map.find(static_cast<int>(it.col()));
                if (col_it != overlap_col_map.end()) local_operator(local_row, col_it->second) = it.value();
            }
        }

        initialize_subdomain_solver(local_operator, subdomain_num);

        switch (d_solver.d_shell_smoother_partition)
        {
        case PETScLevelSolver::ShellSmootherPartition::BASIC:
            cache.update_dofs = cache.overlap_dofs;
            cache.update_local_positions.resize(cache.overlap_dofs.size());
            for (std::size_t local_pos = 0; local_pos < cache.update_local_positions.size(); ++local_pos)
            {
                cache.update_local_positions[local_pos] = static_cast<int>(local_pos);
            }
            break;
        case PETScLevelSolver::ShellSmootherPartition::RESTRICT:
            cache.update_dofs = cache.nonoverlap_dofs;
            cache.update_local_positions = cache.nonoverlap_local_positions;
            break;
        }
        cache.rhs_workspace.resize(overlap_size);
        cache.delta_workspace.resize(overlap_size);

        if (use_multiplicative)
        {
            const std::vector<int>& active_update_dofs = cache.update_dofs;
            const int active_num_cols = static_cast<int>(active_update_dofs.size());
            std::vector<Eigen::Triplet<double>> triplets;
            std::unordered_map<int, int> row_map;
            row_map.reserve(static_cast<std::size_t>(active_num_cols));
            for (int local_col = 0; local_col < active_num_cols; ++local_col)
            {
                const int global_col = active_update_dofs[static_cast<std::size_t>(local_col)];
                for (auto it = Eigen::SparseMatrix<double>::InnerIterator(eigen_level_mat_transpose, global_col); it;
                     ++it)
                {
                    const int row = static_cast<int>(it.row());
                    const auto row_it = row_map.find(row);
                    int local_row = -1;
                    if (row_it == row_map.end())
                    {
                        local_row = static_cast<int>(cache.active_residual_update_rows.size());
                        cache.active_residual_update_rows.push_back(row);
                        row_map.emplace(row, local_row);
                    }
                    else
                    {
                        local_row = row_it->second;
                    }
                    triplets.emplace_back(local_row, local_col, it.value());
                }
            }

            cache.active_residual_update_mat.resize(static_cast<int>(cache.active_residual_update_rows.size()),
                                                    active_num_cols);
            cache.active_residual_update_mat.setFromTriplets(triplets.begin(), triplets.end());
            cache.residual_input_workspace.resize(active_num_cols);
            cache.residual_delta_workspace.resize(static_cast<Eigen::Index>(cache.active_residual_update_rows.size()));
        }
    }
}
} // namespace IBTK

#endif
