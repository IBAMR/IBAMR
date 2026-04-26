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

#ifndef included_IBTK_private_PETScLevelSolverEigenFactorizedLocalSolveShellBackend_inl
#define included_IBTK_private_PETScLevelSolverEigenFactorizedLocalSolveShellBackend_inl

namespace IBTK
{
template <class SolverType>
inline PETScLevelSolverEigenFactorizedLocalSolveShellBackend::TypedSolveStorage<SolverType>&
PETScLevelSolverEigenFactorizedLocalSolveShellBackend::getSolveStorage()
{
    auto* storage = static_cast<TypedSolveStorage<SolverType>*>(d_solve_storage.get());
    TBOX_ASSERT(storage);
    return *storage;
}

template <class SolverType>
inline const PETScLevelSolverEigenFactorizedLocalSolveShellBackend::TypedSolveStorage<SolverType>&
PETScLevelSolverEigenFactorizedLocalSolveShellBackend::getSolveStorage() const
{
    const auto* storage = static_cast<const TypedSolveStorage<SolverType>*>(d_solve_storage.get());
    TBOX_ASSERT(storage);
    return *storage;
}

template <class SolverType>
inline void
PETScLevelSolverEigenFactorizedLocalSolveShellBackend::initializeSolveStorage()
{
    auto storage = std::make_unique<TypedSolveStorage<SolverType>>();
    storage->solvers.resize(d_solver_state.subdomain_dofs->size());
    d_solve_storage = std::move(storage);
}

template <class SolverType>
inline void
PETScLevelSolverEigenFactorizedLocalSolveShellBackend::initializeTypedSubdomainSolver(
    const Eigen::MatrixXd& local_operator,
    const std::size_t subdomain_num)
{
    initializeEigenSolver(getSolveStorage<SolverType>().solvers[subdomain_num], local_operator, getSolverThreshold());
}

template <class SolverType>
inline Eigen::VectorXd
PETScLevelSolverEigenFactorizedLocalSolveShellBackend::solveTypedSubdomainSystem(const Eigen::VectorXd& rhs,
                                                                                 const std::size_t subdomain_num) const
{
    return solveEigenSystem(getSolveStorage<SolverType>().solvers[subdomain_num], rhs);
}
} // namespace IBTK

#endif
