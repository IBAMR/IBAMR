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

#ifndef included_IBTK_private_PETScLevelSolverEigenLocalSolveShellBackend_inl
#define included_IBTK_private_PETScLevelSolverEigenLocalSolveShellBackend_inl

namespace IBTK
{
inline void
PETScLevelSolverEigenLocalSolveShellBackend::initializeSolverState(
    const PETScLevelSolverShellBackendState& solver_state)
{
    checkSerialEigenShellBackend(solver_state, "PETScLevelSolverEigenLocalSolveShellBackend::initializeSolverState()");
    initializeCorrectionCompositionState(solver_state);
    configureFromInputDatabase(getSolverState().input_db);
    initializeAdditionalSolverState();
    initializeCommonDataWithLocalOperatorHook(
        [this](const Eigen::MatrixXd& local_operator, const std::size_t subdomain_num)
        { initializeLocalSubdomainSolver(local_operator, subdomain_num); });
}

inline void
PETScLevelSolverEigenLocalSolveShellBackend::deallocateSolverState()
{
    deallocateCorrectionCompositionState();
    deallocateAdditionalSolverState();
    clearCommonData();
}

inline void
PETScLevelSolverEigenLocalSolveShellBackend::initializeAdditionalSolverState()
{
    return;
}

inline void
PETScLevelSolverEigenLocalSolveShellBackend::deallocateAdditionalSolverState()
{
    return;
}

inline void
PETScLevelSolverEigenLocalSolveShellBackend::solveSubdomain(const std::size_t subdomain_num)
{
    auto& common_subdomains = getCommonSubdomains();
    auto& cache = common_subdomains[subdomain_num];
    cache.delta_workspace = solveLocalSubdomainSystem(cache.rhs_workspace, subdomain_num);
}
} // namespace IBTK

#endif
