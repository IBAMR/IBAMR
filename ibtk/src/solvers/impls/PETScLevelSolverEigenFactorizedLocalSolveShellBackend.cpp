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

#include <ibtk/private/PETScLevelSolverEigenFactorizedLocalSolveShellBackend.h>

#include <tbox/Database.h>

namespace IBTK
{
void
PETScLevelSolverEigenFactorizedLocalSolveShellBackend::configureFromInputDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
{
    d_solver_type = EigenSubdomainSolverType::PARTIAL_PIV_LU;
    d_solver_threshold = -1.0;
    if (!input_db) return;

    if (input_db->keyExists("eigen_subdomain_solver_type"))
    {
        d_solver_type = parseEigenSubdomainSolverType(
            input_db->getString("eigen_subdomain_solver_type"),
            "PETScLevelSolverEigenFactorizedLocalSolveShellBackend::configureFromInputDatabase()");
    }
    if (input_db->keyExists("eigen_subdomain_solver_threshold"))
    {
        d_solver_threshold = input_db->getDouble("eigen_subdomain_solver_threshold");
    }
    return;
}

void
PETScLevelSolverEigenFactorizedLocalSolveShellBackend::initializeLocalSubdomainSolver(
    const Eigen::MatrixXd& local_operator,
    const std::size_t subdomain_num)
{
    dispatchEigenSolverType(getSolverType(),
                            [this, &local_operator, subdomain_num](auto solver_tag)
                            {
                                using SolverType = typename decltype(solver_tag)::type;
                                initializeTypedSubdomainSolver<SolverType>(local_operator, subdomain_num);
                            });
    return;
}

Eigen::VectorXd
PETScLevelSolverEigenFactorizedLocalSolveShellBackend::solveLocalSubdomainSystem(const Eigen::VectorXd& rhs,
                                                                                 const std::size_t subdomain_num) const
{
    Eigen::VectorXd solution;
    dispatchEigenSolverType(getSolverType(),
                            [this, &rhs, subdomain_num, &solution](auto solver_tag)
                            {
                                using SolverType = typename decltype(solver_tag)::type;
                                solution = solveTypedSubdomainSystem<SolverType>(rhs, subdomain_num);
                            });
    return solution;
}

void
PETScLevelSolverEigenFactorizedLocalSolveShellBackend::initializeAdditionalSolverState()
{
    dispatchEigenSolverType(getSolverType(),
                            [this](auto solver_tag)
                            {
                                using SolverType = typename decltype(solver_tag)::type;
                                initializeSolveStorage<SolverType>();
                            });
    return;
}

void
PETScLevelSolverEigenFactorizedLocalSolveShellBackend::deallocateAdditionalSolverState()
{
    d_solve_storage.reset();
    return;
}

PETScLevelSolverEigenFactorizedLocalSolveShellBackend::EigenSubdomainSolverType
PETScLevelSolverEigenFactorizedLocalSolveShellBackend::getSolverType() const
{
    return d_solver_type;
}

double
PETScLevelSolverEigenFactorizedLocalSolveShellBackend::getSolverThreshold() const
{
    return d_solver_threshold;
}
} // namespace IBTK
