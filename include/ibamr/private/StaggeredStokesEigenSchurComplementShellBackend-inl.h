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

#ifndef included_IBAMR_private_StaggeredStokesEigenSchurComplementShellBackend_inl
#define included_IBAMR_private_StaggeredStokesEigenSchurComplementShellBackend_inl

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

inline IBTK::PETScLevelSolverEigenShellBackendBase::SubdomainSweepView
StaggeredStokesEigenSchurComplementShellBackend::getCustomSubdomainSweepView(
    CustomEigenSchurSubdomainCache& custom_cache) const
{
    TBOX_ASSERT(custom_cache.overlap_dofs);
    TBOX_ASSERT(custom_cache.update_dofs);
    TBOX_ASSERT(custom_cache.update_local_positions);
    TBOX_ASSERT(custom_cache.active_residual_update_rows);
    TBOX_ASSERT(custom_cache.active_residual_update_mat);
    return { *custom_cache.overlap_dofs,
             *custom_cache.update_dofs,
             *custom_cache.update_local_positions,
             *custom_cache.active_residual_update_rows,
             *custom_cache.active_residual_update_mat,
             custom_cache.rhs_workspace,
             custom_cache.delta_workspace,
             custom_cache.residual_input_workspace,
             custom_cache.residual_delta_workspace };
}

template <class SolverType>
inline void
StaggeredStokesEigenSchurComplementShellBackend::initializeCustomEigenA00Solver(SolverType& solver,
                                                                                const Eigen::MatrixXd& matrix,
                                                                                const std::size_t subdomain_num) const
{
    initializeEigenSolver(solver, matrix, d_a00_solver_threshold);
    if constexpr (std::is_same_v<SolverType, Eigen::LLT<Eigen::MatrixXd>> ||
                  std::is_same_v<SolverType, Eigen::LDLT<Eigen::MatrixXd>>)
    {
        const char* const factorization = std::is_same_v<SolverType, Eigen::LLT<Eigen::MatrixXd>> ? "LLT" : "LDLT";
        if (solver.info() != Eigen::Success)
        {
            TBOX_ERROR("StaggeredStokesEigenSchurComplementShellBackend::initializeCustomEigenA00Solver():\n"
                       << "  " << factorization << " factorization failed for the local A00 block on subdomain "
                       << subdomain_num << ".\n");
        }
    }
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
StaggeredStokesEigenSchurComplementShellBackend::applyAdditive(Vec x, Vec y)
{
    auto& typed_storage = getCustomEigenA00SolveStorage<SolverType>();
    applyAdditiveSubdomainSweep(
        x,
        y,
        d_subdomain_caches.size(),
        [this](const std::size_t subdomain_num)
        { return getCustomSubdomainSweepView(d_subdomain_caches[subdomain_num]); },
        [this, &typed_storage](SubdomainSweepView& /*view*/, const std::size_t subdomain_num)
        { solveCustomEigenSubdomain(d_subdomain_caches[subdomain_num], typed_storage.solvers[subdomain_num]); });
}

template <class SolverType>
inline void
StaggeredStokesEigenSchurComplementShellBackend::applyMultiplicative(Vec x, Vec y)
{
    auto& typed_storage = getCustomEigenA00SolveStorage<SolverType>();
    applyMultiplicativeSubdomainSweep(
        x,
        y,
        d_subdomain_caches.size(),
        [this](const std::size_t subdomain_num)
        { return getCustomSubdomainSweepView(d_subdomain_caches[subdomain_num]); },
        [this, &typed_storage](SubdomainSweepView& /*view*/, const std::size_t subdomain_num)
        { solveCustomEigenSubdomain(d_subdomain_caches[subdomain_num], typed_storage.solvers[subdomain_num]); });
}
} // namespace IBAMR

#endif
