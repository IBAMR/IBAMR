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
    const bool use_restrict_partition = d_solver_state.use_restrict_partition;
    const Eigen::SparseMatrix<double, Eigen::RowMajor> eigen_level_mat =
        copyPETScMatToEigenSparse(d_solver_state.petsc_mat);
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
    }
    finalizeCorrectionCompositionState();
}

inline std::size_t
PETScLevelSolverEigenShellBackendBase::getNumberOfSubdomains() const
{
    return d_common_subdomains.size();
}

inline void
PETScLevelSolverEigenShellBackendBase::initializeSubdomainSweep(Vec x, Vec y)
{
    TBOX_ASSERT(IBTK_MPI::getNodes() == 1);
    TBOX_ASSERT(getNumDofs() > 0);
    const int ierr = VecGetArray(y, &d_sweep_y_values);
    IBTK_CHKERRQ(ierr);
}

inline void
PETScLevelSolverEigenShellBackendBase::beginSubdomainRhs(const std::size_t subdomain_num, Vec x, Vec /*y*/)
{
    auto& cache = d_common_subdomains[subdomain_num];
    Vec source = getSubdomainResidualSource(x);
    const PetscScalar* source_values = nullptr;
    int ierr = VecGetArrayRead(source, &source_values);
    IBTK_CHKERRQ(ierr);
    const auto source_map =
        Eigen::Map<const Eigen::VectorXd>(reinterpret_cast<const double*>(source_values), getNumDofs());
    std::size_t rhs_idx = 0;
    for (const int dof : cache.overlap_dofs)
    {
        cache.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = source_map[dof];
    }
    ierr = VecRestoreArrayRead(source, &source_values);
    IBTK_CHKERRQ(ierr);
}

inline void
PETScLevelSolverEigenShellBackendBase::endSubdomainRhs(std::size_t /*subdomain_num*/, Vec /*x*/, Vec /*y*/)
{
}

inline void
PETScLevelSolverEigenShellBackendBase::accumulateSubdomainCorrection(const std::size_t subdomain_num, Vec /*y*/)
{
    auto& cache = d_common_subdomains[subdomain_num];
    auto y_map = Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double*>(d_sweep_y_values), getNumDofs());
    std::size_t update_idx = 0;
    for (const int dof : cache.update_dofs)
    {
        y_map[dof] += cache.delta_workspace[static_cast<Eigen::Index>(cache.update_local_positions[update_idx++])];
    }
}

inline const std::vector<int>&
PETScLevelSolverEigenShellBackendBase::getSubdomainCorrectionDofs(const std::size_t subdomain_num) const
{
    return d_common_subdomains[subdomain_num].update_dofs;
}

inline void
PETScLevelSolverEigenShellBackendBase::copySubdomainCorrection(const std::size_t subdomain_num,
                                                               PetscScalar* correction_values)
{
    const auto& cache = d_common_subdomains[subdomain_num];
    for (std::size_t update_num = 0; update_num < cache.update_local_positions.size(); ++update_num)
        correction_values[update_num] = cache.delta_workspace[cache.update_local_positions[update_num]];
}

inline void
PETScLevelSolverEigenShellBackendBase::finalizeSubdomainSweep(Vec x, Vec y)
{
    int ierr = VecRestoreArray(y, &d_sweep_y_values);
    IBTK_CHKERRQ(ierr);
    d_sweep_y_values = nullptr;
}
} // namespace IBTK

#endif
