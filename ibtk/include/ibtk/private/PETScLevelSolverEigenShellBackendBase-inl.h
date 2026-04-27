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

#include <ibtk/IBTK_MPI.h>

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
        d_solver_state.object_name, d_solver_state.options_prefix, solver_type, std::forward<Handler>(handler));
}

inline PETScLevelSolverEigenShellBackendBase::EigenSubdomainSolverType
PETScLevelSolverEigenShellBackendBase::parseEigenSubdomainSolverType(const std::string& type, const char* caller) const
{
    return PETScLevelSolverEigenShell::parseSubdomainSolverType(
        type, d_solver_state.object_name, d_solver_state.options_prefix, caller);
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

inline void
PETScLevelSolverEigenShellBackendBase::checkSerialEigenShellBackend(
    const PETScLevelSolverShellBackendState& solver_state,
    const char* caller)
{
    if (IBTK_MPI::getNodes() != 1)
    {
        TBOX_ERROR(solver_state.object_name << " " << solver_state.options_prefix << " " << caller << ":\n"
                                            << "  Eigen shell smoother backends are currently serial-only.\n");
    }
}

inline Eigen::MatrixXd
PETScLevelSolverEigenShellBackendBase::buildDenseEigenSolveMatrix(const EigenSubdomainSolverType solver_type,
                                                                  const Eigen::MatrixXd& matrix,
                                                                  const double threshold,
                                                                  const char* caller) const
{
    Eigen::MatrixXd solve_matrix;
    dispatchEigenSolverType(
        solver_type,
        [this, &matrix, threshold, caller, &solve_matrix](auto solver_tag)
        {
            using SolverType = typename decltype(solver_tag)::type;
            if constexpr (std::is_same_v<SolverType, Eigen::LLT<Eigen::MatrixXd>>)
            {
                solve_matrix = buildLLTSolveMatrix(matrix);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::LDLT<Eigen::MatrixXd>>)
            {
                solve_matrix = buildLDLTSolveMatrix(matrix);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::PartialPivLU<Eigen::MatrixXd>>)
            {
                solve_matrix = buildPartialPivLUSolveMatrix(matrix);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::FullPivLU<Eigen::MatrixXd>>)
            {
                solve_matrix = buildFullPivLUSolveMatrix(matrix, threshold);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::HouseholderQR<Eigen::MatrixXd>>)
            {
                solve_matrix = buildHouseholderQRSolveMatrix(matrix);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::ColPivHouseholderQR<Eigen::MatrixXd>>)
            {
                solve_matrix = buildQRSolveMatrix<SolverType>(matrix, threshold);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd>>)
            {
                solve_matrix = buildCompleteOrthogonalDecompositionPseudoinverse(matrix, threshold);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::FullPivHouseholderQR<Eigen::MatrixXd>>)
            {
                solve_matrix = buildQRSolveMatrix<SolverType>(matrix, threshold);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::JacobiSVD<Eigen::MatrixXd>> ||
                               std::is_same_v<SolverType, Eigen::BDCSVD<Eigen::MatrixXd>>)
            {
                solve_matrix = buildSVDPseudoinverse<SolverType>(matrix, threshold);
            }
            else
            {
                TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix << " " << caller << ":\n"
                                                      << "Unsupported Eigen subdomain solver type.\n");
            }
        });
    return solve_matrix;
}

template <class InitializeSubdomainSolver>
inline void
PETScLevelSolverEigenShellBackendBase::initializeCommonDataWithLocalOperatorHook(
    InitializeSubdomainSolver initialize_subdomain_solver)
{
    const bool use_multiplicative = d_solver_state.use_multiplicative;
    const bool use_restrict_partition = d_solver_state.use_restrict_partition;
    const Eigen::SparseMatrix<double, Eigen::RowMajor> eigen_level_mat =
        copyPETScMatToEigenSparse(d_solver_state.petsc_mat);
    Eigen::SparseMatrix<double> eigen_level_mat_transpose;
    if (use_multiplicative) eigen_level_mat_transpose = eigen_level_mat.transpose();
    d_n_dofs = eigen_level_mat.rows();
    d_common_subdomains.clear();
    d_common_subdomains.resize(d_solver_state.subdomain_dofs->size());
    for (std::size_t subdomain_num = 0; subdomain_num < d_common_subdomains.size(); ++subdomain_num)
    {
        auto& cache = d_common_subdomains[subdomain_num];
        cache = CommonSubdomainCache();
        cache.overlap_dofs = (*d_solver_state.subdomain_dofs)[subdomain_num];
        const Eigen::Index overlap_size = static_cast<Eigen::Index>(cache.overlap_dofs.size());
        std::unordered_map<int, int> overlap_col_map;
        overlap_col_map.reserve(static_cast<std::size_t>(overlap_size));
        for (Eigen::Index local_col = 0; local_col < overlap_size; ++local_col)
        {
            overlap_col_map[cache.overlap_dofs[static_cast<std::size_t>(local_col)]] = static_cast<int>(local_col);
        }

        cache.nonoverlap_dofs = (*d_solver_state.nonoverlap_subdomain_dofs)[subdomain_num];
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

        if (use_restrict_partition)
        {
            cache.update_dofs = cache.nonoverlap_dofs;
            cache.update_local_positions = cache.nonoverlap_local_positions;
        }
        else
        {
            cache.update_dofs = cache.overlap_dofs;
            cache.update_local_positions.resize(cache.overlap_dofs.size());
            for (std::size_t local_pos = 0; local_pos < cache.update_local_positions.size(); ++local_pos)
            {
                cache.update_local_positions[local_pos] = static_cast<int>(local_pos);
            }
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

inline PETScLevelSolverEigenShellBackendBase::SubdomainSweepView
PETScLevelSolverEigenShellBackendBase::getCommonSubdomainSweepView(CommonSubdomainCache& cache)
{
    return { cache.overlap_dofs,
             cache.update_dofs,
             cache.update_local_positions,
             cache.active_residual_update_rows,
             cache.active_residual_update_mat,
             cache.rhs_workspace,
             cache.delta_workspace,
             cache.residual_input_workspace,
             cache.residual_delta_workspace };
}

template <class GetSubdomainSweepView, class SolveSubdomain>
inline void
PETScLevelSolverEigenShellBackendBase::applyAdditiveSubdomainSweep(Vec x,
                                                                   Vec y,
                                                                   const std::size_t n_subdomains,
                                                                   GetSubdomainSweepView get_subdomain_sweep_view,
                                                                   SolveSubdomain solve_subdomain)
{
    TBOX_ASSERT(IBTK_MPI::getNodes() == 1);
    const Eigen::Index n = getNumDofs();
    TBOX_ASSERT(n > 0);

    {
        ConstPetscVecArrayMap x_array(x, n);
        PetscVecArrayMap y_array(y, n);
        const auto x_map = x_array.getMap();
        auto y_map = y_array.getMap();
        y_map.setZero();
        for (std::size_t subdomain_num = 0; subdomain_num < n_subdomains; ++subdomain_num)
        {
            SubdomainSweepView view = get_subdomain_sweep_view(subdomain_num);
            std::size_t rhs_idx = 0;
            for (const int dof : view.overlap_dofs)
            {
                view.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = x_map[dof];
            }
            solve_subdomain(view, subdomain_num);
            std::size_t update_idx = 0;
            for (const int dof : view.update_dofs)
            {
                y_map[dof] +=
                    view.delta_workspace[static_cast<Eigen::Index>(view.update_local_positions[update_idx++])];
            }
        }
    }
    getSolverState().postprocess_result(y);
}

template <class GetSubdomainSweepView, class SolveSubdomain>
inline void
PETScLevelSolverEigenShellBackendBase::applyMultiplicativeSubdomainSweep(Vec x,
                                                                         Vec y,
                                                                         const std::size_t n_subdomains,
                                                                         GetSubdomainSweepView get_subdomain_sweep_view,
                                                                         SolveSubdomain solve_subdomain)
{
    TBOX_ASSERT(IBTK_MPI::getNodes() == 1);
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
        for (std::size_t subdomain_num = 0; subdomain_num < n_subdomains; ++subdomain_num)
        {
            SubdomainSweepView view = get_subdomain_sweep_view(subdomain_num);
            std::size_t rhs_idx = 0;
            for (const int dof : view.overlap_dofs)
            {
                view.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = residual[dof];
            }
            solve_subdomain(view, subdomain_num);

            std::size_t update_idx = 0;
            for (const int dof : view.update_dofs)
            {
                y_map[dof] +=
                    view.delta_workspace[static_cast<Eigen::Index>(view.update_local_positions[update_idx++])];
            }
            if (subdomain_num + 1 < n_subdomains && view.active_residual_update_mat.rows() > 0)
            {
                std::size_t residual_input_idx = 0;
                for (const int local_pos : view.update_local_positions)
                {
                    view.residual_input_workspace[static_cast<Eigen::Index>(residual_input_idx++)] =
                        view.delta_workspace[static_cast<Eigen::Index>(local_pos)];
                }
                view.residual_delta_workspace.noalias() =
                    view.active_residual_update_mat * view.residual_input_workspace;
                std::size_t row_idx = 0;
                for (const int row : view.active_residual_update_rows)
                {
                    residual[row] -= view.residual_delta_workspace[static_cast<Eigen::Index>(row_idx++)];
                }
            }
        }
    }
    getSolverState().postprocess_result(y);
}
} // namespace IBTK

#endif
