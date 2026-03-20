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

#ifndef included_IBAMR_private_StaggeredStokesEigenSchurComplementShellBackend_inl
#define included_IBAMR_private_StaggeredStokesEigenSchurComplementShellBackend_inl

#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>

#include <type_traits>

namespace IBAMR
{
namespace StaggeredStokesEigenSchurComplementShellBackendDetail
{
inline void
extract_subvector_into(Eigen::VectorXd& subvector, const Eigen::VectorXd& vector, const std::vector<int>& positions)
{
    TBOX_ASSERT(subvector.size() == static_cast<Eigen::Index>(positions.size()));
    for (std::size_t k = 0; k < positions.size(); ++k)
    {
        subvector(static_cast<Eigen::Index>(k)) = vector(positions[k]);
    }
}

inline void
insert_subvector(Eigen::VectorXd& vector, const std::vector<int>& positions, const Eigen::VectorXd& subvector)
{
    for (std::size_t k = 0; k < positions.size(); ++k)
    {
        vector(positions[k]) = subvector(static_cast<Eigen::Index>(k));
    }
}

template <class SolverType>
inline void
initialize_custom_a00_solver(SolverType& solver,
                             const Eigen::MatrixXd& matrix,
                             const double threshold,
                             const std::size_t subdomain_num)
{
    if constexpr (std::is_same_v<SolverType, Eigen::LLT<Eigen::MatrixXd>>)
    {
        solver.compute(matrix);
        if (solver.info() != Eigen::Success)
        {
            TBOX_ERROR("initialize_custom_a00_solver():\n"
                       << "  LLT factorization failed for the local A00 block on subdomain " << subdomain_num << ".\n");
        }
    }
    else if constexpr (std::is_same_v<SolverType, Eigen::LDLT<Eigen::MatrixXd>>)
    {
        solver.compute(matrix);
        if (solver.info() != Eigen::Success)
        {
            TBOX_ERROR("initialize_custom_a00_solver():\n"
                       << "  LDLT factorization failed for the local A00 block on subdomain " << subdomain_num
                       << ".\n");
        }
    }
    else
    {
        if constexpr (std::is_same_v<SolverType, Eigen::JacobiSVD<Eigen::MatrixXd>> ||
                      std::is_same_v<SolverType, Eigen::BDCSVD<Eigen::MatrixXd>>)
        {
            solver.compute(matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
            if (threshold >= 0.0) solver.setThreshold(threshold);
        }
        else
        {
            solver.compute(matrix);
            if constexpr (std::is_same_v<SolverType, Eigen::FullPivLU<Eigen::MatrixXd>> ||
                          std::is_same_v<SolverType, Eigen::ColPivHouseholderQR<Eigen::MatrixXd>> ||
                          std::is_same_v<SolverType, Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd>> ||
                          std::is_same_v<SolverType, Eigen::FullPivHouseholderQR<Eigen::MatrixXd>>)
            {
                if (threshold >= 0.0) solver.setThreshold(threshold);
            }
        }
    }
}

template <class SolverType, class RhsType>
inline auto
solve_custom_a00(const SolverType& solver, const RhsType& rhs) -> decltype(solver.solve(rhs))
{
    return solver.solve(rhs);
}
} // namespace StaggeredStokesEigenSchurComplementShellBackendDetail

template <class SolverType>
inline StaggeredStokesEigenSchurComplementShellBackend::CustomEigenA00TypedSolveStorage<SolverType>&
StaggeredStokesEigenSchurComplementShellBackend::getCustomEigenA00SolveStorage()
{
    auto* storage = static_cast<CustomEigenA00TypedSolveStorage<SolverType>*>(d_a00_solver_storage.get());
    TBOX_ASSERT(storage);
    return *storage;
}

template <class SolverType>
inline const StaggeredStokesEigenSchurComplementShellBackend::CustomEigenA00TypedSolveStorage<SolverType>&
StaggeredStokesEigenSchurComplementShellBackend::getCustomEigenA00SolveStorage() const
{
    const auto* storage = static_cast<const CustomEigenA00TypedSolveStorage<SolverType>*>(d_a00_solver_storage.get());
    TBOX_ASSERT(storage);
    return *storage;
}

template <class SolverType>
inline void
StaggeredStokesEigenSchurComplementShellBackend::initializeCustomEigenA00SolveStorage(const std::size_t n_subdomains)
{
    auto storage = std::make_unique<CustomEigenA00TypedSolveStorage<SolverType>>();
    storage->solvers.resize(n_subdomains);
    d_a00_solver_storage = std::move(storage);
}

template <class SolverType>
inline void
StaggeredStokesEigenSchurComplementShellBackend::solveCustomEigenSubdomain(CustomEigenSchurSubdomainCache& custom_cache,
                                                                           const SolverType& a00_solver) const
{
    using namespace StaggeredStokesEigenSchurComplementShellBackendDetail;
    extract_subvector_into(
        custom_cache.velocity_rhs_workspace, custom_cache.rhs_workspace, custom_cache.velocity_positions);
    extract_subvector_into(
        custom_cache.pressure_rhs_workspace, custom_cache.rhs_workspace, custom_cache.pressure_positions);
    custom_cache.delta_workspace.setZero();

    if (custom_cache.pressure_positions.empty())
    {
        custom_cache.velocity_solution_workspace =
            solve_custom_a00(a00_solver, custom_cache.velocity_rhs_workspace).eval();
        insert_subvector(
            custom_cache.delta_workspace, custom_cache.velocity_positions, custom_cache.velocity_solution_workspace);
    }
    else if (custom_cache.velocity_positions.empty())
    {
        custom_cache.pressure_solution_workspace.noalias() =
            custom_cache.schur_solve_matrix * custom_cache.pressure_rhs_workspace;
        insert_subvector(
            custom_cache.delta_workspace, custom_cache.pressure_positions, custom_cache.pressure_solution_workspace);
    }
    else
    {
        custom_cache.velocity_solution_workspace =
            solve_custom_a00(a00_solver, custom_cache.velocity_rhs_workspace).eval();
        custom_cache.pressure_rhs_workspace -= custom_cache.A10 * custom_cache.velocity_solution_workspace;
        custom_cache.pressure_solution_workspace.noalias() =
            custom_cache.schur_solve_matrix * custom_cache.pressure_rhs_workspace;
        custom_cache.velocity_solution_workspace -= custom_cache.A00_inv_A01 * custom_cache.pressure_solution_workspace;
        insert_subvector(
            custom_cache.delta_workspace, custom_cache.velocity_positions, custom_cache.velocity_solution_workspace);
        insert_subvector(
            custom_cache.delta_workspace, custom_cache.pressure_positions, custom_cache.pressure_solution_workspace);
    }
}

template <class SolverType>
inline void
StaggeredStokesEigenSchurComplementShellBackend::applyAdditiveImpl(Vec x, Vec y)
{
    TBOX_ASSERT(IBTK::IBTK_MPI::getNodes() == 1);
    auto& typed_storage = getCustomEigenA00SolveStorage<SolverType>();
    PetscInt n_local = 0;
    int ierr = VecGetLocalSize(x, &n_local);
    IBTK_CHKERRQ(ierr);
    const Eigen::Index n = static_cast<Eigen::Index>(n_local);
    TBOX_ASSERT(n > 0);

    {
        ConstPetscVecArrayMap x_array(x, n);
        PetscVecArrayMap y_array(y, n);
        const auto x_map = x_array.getMap();
        auto y_map = y_array.getMap();
        y_map.setZero();
        for (std::size_t subdomain_num = 0; subdomain_num < d_subdomain_caches.size(); ++subdomain_num)
        {
            auto& custom_cache = d_subdomain_caches[subdomain_num];
            const auto& a00_solver = typed_storage.solvers[subdomain_num];
            const auto& overlap_dofs = *custom_cache.overlap_dofs;
            const auto& update_dofs = *custom_cache.update_dofs;
            const auto& update_local_positions = *custom_cache.update_local_positions;
            std::size_t rhs_idx = 0;
            for (const int dof : overlap_dofs)
            {
                custom_cache.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = x_map[dof];
            }

            solveCustomEigenSubdomain(custom_cache, a00_solver);

            std::size_t update_idx = 0;
            for (const int dof : update_dofs)
            {
                y_map[dof] +=
                    custom_cache.delta_workspace[static_cast<Eigen::Index>(update_local_positions[update_idx++])];
            }
        }
    }
    d_stokes_solver.postprocessShellResult(y);
}

template <class SolverType>
inline void
StaggeredStokesEigenSchurComplementShellBackend::applyMultiplicativeImpl(Vec x, Vec y)
{
    TBOX_ASSERT(IBTK::IBTK_MPI::getNodes() == 1);
    auto& typed_storage = getCustomEigenA00SolveStorage<SolverType>();
    PetscInt n_local = 0;
    int ierr = VecGetLocalSize(x, &n_local);
    IBTK_CHKERRQ(ierr);
    const Eigen::Index n = static_cast<Eigen::Index>(n_local);
    TBOX_ASSERT(n > 0);

    {
        ConstPetscVecArrayMap x_array(x, n);
        PetscVecArrayMap y_array(y, n);
        const auto x_map = x_array.getMap();
        auto y_map = y_array.getMap();
        Eigen::VectorXd residual(n);
        y_map.setZero();
        residual = x_map;
        const std::size_t n_subdomains = d_subdomain_caches.size();
        for (std::size_t subdomain_num = 0; subdomain_num + 1 < n_subdomains; ++subdomain_num)
        {
            auto& custom_cache = d_subdomain_caches[subdomain_num];
            const auto& a00_solver = typed_storage.solvers[subdomain_num];
            const auto& overlap_dofs = *custom_cache.overlap_dofs;
            const auto& update_dofs = *custom_cache.update_dofs;
            const auto& update_local_positions = *custom_cache.update_local_positions;
            const auto& active_residual_update_rows = *custom_cache.active_residual_update_rows;
            const auto& active_residual_update_mat = *custom_cache.active_residual_update_mat;
            std::size_t rhs_idx = 0;
            for (const int dof : overlap_dofs)
            {
                custom_cache.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = residual[dof];
            }

            solveCustomEigenSubdomain(custom_cache, a00_solver);

            std::size_t update_idx = 0;
            for (const int dof : update_dofs)
            {
                y_map[dof] +=
                    custom_cache.delta_workspace[static_cast<Eigen::Index>(update_local_positions[update_idx++])];
            }
            if (active_residual_update_mat.rows() > 0)
            {
                std::size_t residual_input_idx = 0;
                for (const int local_pos : update_local_positions)
                {
                    custom_cache.residual_input_workspace[static_cast<Eigen::Index>(residual_input_idx++)] =
                        custom_cache.delta_workspace[static_cast<Eigen::Index>(local_pos)];
                }
                custom_cache.residual_delta_workspace.noalias() =
                    active_residual_update_mat * custom_cache.residual_input_workspace;
                std::size_t row_idx = 0;
                for (const int row : active_residual_update_rows)
                {
                    residual[row] -= custom_cache.residual_delta_workspace[static_cast<Eigen::Index>(row_idx++)];
                }
            }
        }
        if (n_subdomains > 0)
        {
            auto& custom_cache = d_subdomain_caches.back();
            const auto& a00_solver = typed_storage.solvers.back();
            const auto& overlap_dofs = *custom_cache.overlap_dofs;
            const auto& update_dofs = *custom_cache.update_dofs;
            const auto& update_local_positions = *custom_cache.update_local_positions;
            std::size_t rhs_idx = 0;
            for (const int dof : overlap_dofs)
            {
                custom_cache.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = residual[dof];
            }

            solveCustomEigenSubdomain(custom_cache, a00_solver);

            std::size_t update_idx = 0;
            for (const int dof : update_dofs)
            {
                y_map[dof] +=
                    custom_cache.delta_workspace[static_cast<Eigen::Index>(update_local_positions[update_idx++])];
            }
        }
    }
    d_stokes_solver.postprocessShellResult(y);
}
} // namespace IBAMR

#endif
