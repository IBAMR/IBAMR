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

#include <ibtk/PETScLevelSolver.h>
#include <ibtk/private/PETScLevelSolverEigenShellBackend.h>

#include <tbox/Database.h>

namespace IBTK
{
PETScLevelSolverEigenShellBackend::PETScLevelSolverEigenShellBackend(PETScLevelSolver& solver)
    : PETScLevelSolverEigenShellBackendBase(solver)
{
}

const std::string&
PETScLevelSolverEigenShellBackend::getTypeKey() const
{
    return d_type_key;
}

const char*
PETScLevelSolverEigenShellBackend::getPCNameSuffixAdditive() const
{
    return "PC_AdditiveEigen";
}

const char*
PETScLevelSolverEigenShellBackend::getPCNameSuffixMultiplicative() const
{
    return "PC_MultiplicativeEigen";
}

void
PETScLevelSolverEigenShellBackend::configure(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
{
    d_solver_type = EigenSubdomainSolverType::PARTIAL_PIV_LU;
    d_solver_threshold = -1.0;
    if (!input_db) return;

    if (input_db->keyExists("eigen_subdomain_solver_type"))
    {
        d_solver_type = parseEigenSubdomainSolverType(input_db->getString("eigen_subdomain_solver_type"),
                                                      "PETScLevelSolverEigenShellBackend::configure()");
    }
    if (input_db->keyExists("eigen_subdomain_solver_threshold"))
    {
        d_solver_threshold = input_db->getDouble("eigen_subdomain_solver_threshold");
    }
}

void
PETScLevelSolverEigenShellBackend::initialize()
{
    dispatchEigenSolverType(getSolverType(),
                            [this](auto solver_tag)
                            {
                                using SolverType = typename decltype(solver_tag)::type;
                                initializeSolveStorage<SolverType>();
                            });
}

void
PETScLevelSolverEigenShellBackend::deallocate()
{
    clearCommonData();
    d_solve_storage.reset();
}

void
PETScLevelSolverEigenShellBackend::initializeSubdomainSolver(const Eigen::MatrixXd& local_operator,
                                                             const std::size_t subdomain_num)
{
    dispatchEigenSolverType(getSolverType(),
                            [this, &local_operator, subdomain_num](auto solver_tag)
                            {
                                using SolverType = typename decltype(solver_tag)::type;
                                initializeSubdomainSolver<SolverType>(local_operator, subdomain_num);
                            });
}

Eigen::VectorXd
PETScLevelSolverEigenShellBackend::solveSubdomainSystem(const Eigen::VectorXd& rhs,
                                                        const std::size_t subdomain_num) const
{
    Eigen::VectorXd solution;
    dispatchEigenSolverType(getSolverType(),
                            [this, &rhs, subdomain_num, &solution](auto solver_tag)
                            {
                                using SolverType = typename decltype(solver_tag)::type;
                                solution = solveSubdomainSystem<SolverType>(rhs, subdomain_num);
                            });
    return solution;
}

void
PETScLevelSolverEigenShellBackend::applyAdditive(Vec x, Vec y)
{
    dispatchEigenSolverType(getSolverType(),
                            [this, x, y](auto solver_tag)
                            {
                                using SolverType = typename decltype(solver_tag)::type;
                                applyAdditiveImpl<SolverType>(x, y);
                            });
}

void
PETScLevelSolverEigenShellBackend::applyMultiplicative(Vec x, Vec y)
{
    dispatchEigenSolverType(getSolverType(),
                            [this, x, y](auto solver_tag)
                            {
                                using SolverType = typename decltype(solver_tag)::type;
                                applyMultiplicativeImpl<SolverType>(x, y);
                            });
}

PETScLevelSolverEigenShellBackend::EigenSubdomainSolverType
PETScLevelSolverEigenShellBackend::getSolverType() const
{
    return d_solver_type;
}

double
PETScLevelSolverEigenShellBackend::getSolverThreshold() const
{
    return d_solver_threshold;
}
} // namespace IBTK
