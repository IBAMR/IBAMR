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

#ifndef included_IBTK_private_PETScLevelSolverEigenShellBackend_inl
#define included_IBTK_private_PETScLevelSolverEigenShellBackend_inl

#include <ibtk/IBTK_MPI.h>

namespace IBTK
{
template <class SolverType>
inline PETScLevelSolverEigenShellBackend::TypedSolveStorage<SolverType>&
PETScLevelSolverEigenShellBackend::getSolveStorage()
{
    auto* storage = static_cast<TypedSolveStorage<SolverType>*>(d_solve_storage.get());
    TBOX_ASSERT(storage);
    return *storage;
}

template <class SolverType>
inline const PETScLevelSolverEigenShellBackend::TypedSolveStorage<SolverType>&
PETScLevelSolverEigenShellBackend::getSolveStorage() const
{
    const auto* storage = static_cast<const TypedSolveStorage<SolverType>*>(d_solve_storage.get());
    TBOX_ASSERT(storage);
    return *storage;
}

template <class SolverType>
inline void
PETScLevelSolverEigenShellBackend::initializeSolveStorage()
{
    auto storage = std::make_unique<TypedSolveStorage<SolverType>>();
    storage->solvers.resize(d_solver.d_subdomain_dofs.size());
    d_solve_storage = std::move(storage);

    initializeCommonDataWithLocalOperatorHook(
        [this](const Eigen::MatrixXd& local_operator, const std::size_t subdomain_num)
        {
            initializeEigenSolver(
                getSolveStorage<SolverType>().solvers[subdomain_num], local_operator, getSolverThreshold());
        });
}

template <class SolverType>
inline void
PETScLevelSolverEigenShellBackend::initializeSubdomainSolver(const Eigen::MatrixXd& local_operator,
                                                             const std::size_t subdomain_num)
{
    initializeEigenSolver(getSolveStorage<SolverType>().solvers[subdomain_num], local_operator, getSolverThreshold());
}

template <class SolverType>
inline Eigen::VectorXd
PETScLevelSolverEigenShellBackend::solveSubdomainSystem(const Eigen::VectorXd& rhs,
                                                        const std::size_t subdomain_num) const
{
    return solveEigenSystem(getSolveStorage<SolverType>().solvers[subdomain_num], rhs);
}

template <class SolverType>
inline void
PETScLevelSolverEigenShellBackend::applyAdditiveImpl(Vec x, Vec y)
{
    TBOX_ASSERT(IBTK_MPI::getNodes() == 1);
    auto& common_subdomains = getCommonSubdomains();
    const auto& typed_storage = getSolveStorage<SolverType>();
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
            auto& solver = typed_storage.solvers[subdomain_num];
            std::size_t rhs_idx = 0;
            for (const int dof : cache.overlap_dofs)
            {
                cache.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = x_map[dof];
            }
            cache.delta_workspace = solveEigenSystem(solver, cache.rhs_workspace);
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

template <class SolverType>
inline void
PETScLevelSolverEigenShellBackend::applyMultiplicativeImpl(Vec x, Vec y)
{
    TBOX_ASSERT(IBTK_MPI::getNodes() == 1);
    auto& common_subdomains = getCommonSubdomains();
    const auto& typed_storage = getSolveStorage<SolverType>();
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
            auto& solver = typed_storage.solvers[subdomain_num];
            std::size_t rhs_idx = 0;
            for (const int dof : cache.overlap_dofs)
            {
                cache.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = residual[dof];
            }
            cache.delta_workspace = solveEigenSystem(solver, cache.rhs_workspace);

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
            auto& solver = typed_storage.solvers.back();
            std::size_t rhs_idx = 0;
            for (const int dof : cache.overlap_dofs)
            {
                cache.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = residual[dof];
            }
            cache.delta_workspace = solveEigenSystem(solver, cache.rhs_workspace);

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

#endif
