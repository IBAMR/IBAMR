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

#include <ibtk/PETScLevelSolver.h>
#include <ibtk/private/PETScLevelSolverEigenPseudoinverseShellBackend.h>

#include <tbox/Database.h>

namespace IBTK
{
namespace
{
std::unique_ptr<PETScLevelSolverShellBackend>
allocate_eigen_pseudoinverse_shell_backend(PETScLevelSolver& /*solver*/)
{
    return std::make_unique<PETScLevelSolverEigenPseudoinverseShellBackend>();
}

class PETScLevelSolverEigenPseudoinverseShellBackendRegistrar
{
public:
    PETScLevelSolverEigenPseudoinverseShellBackendRegistrar()
    {
        PETScLevelSolverShellBackendManager::getManager()->registerShellBackendFactoryFunction(
            "eigen-pseudoinverse", allocate_eigen_pseudoinverse_shell_backend);
    }
};

static PETScLevelSolverEigenPseudoinverseShellBackendRegistrar eigen_pseudoinverse_shell_backend_registrar;
} // namespace

const std::string PETScLevelSolverEigenPseudoinverseShellBackend::s_backend_name = "EigenPseudoinverse";

void
PETScLevelSolverEigenPseudoinverseShellBackend::configureFromInputDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
{
    d_solver_type = EigenSubdomainSolverType::COL_PIV_HOUSEHOLDER_QR;
    d_solver_threshold = -1.0;
    if (!input_db) return;

    if (input_db->keyExists("eigen_subdomain_pseudoinverse_type"))
    {
        d_solver_type = parseEigenSubdomainSolverType(
            input_db->getString("eigen_subdomain_pseudoinverse_type"),
            "PETScLevelSolverEigenPseudoinverseShellBackend::PETScLevelSolverEigenPseudoinverseShellBackend()");
    }
    if (input_db->keyExists("eigen_subdomain_pseudoinverse_threshold"))
    {
        d_solver_threshold = input_db->getDouble("eigen_subdomain_pseudoinverse_threshold");
    }
}

const std::string&
PETScLevelSolverEigenPseudoinverseShellBackend::getName() const
{
    return s_backend_name;
}

void
PETScLevelSolverEigenPseudoinverseShellBackend::initializeLocalSubdomainSolver(const Eigen::MatrixXd& local_operator,
                                                                               const std::size_t subdomain_num)
{
    getCommonSubdomains()[subdomain_num].local_pseudoinverse = buildSubdomainPseudoinverse(local_operator);
}

Eigen::VectorXd
PETScLevelSolverEigenPseudoinverseShellBackend::solveLocalSubdomainSystem(const Eigen::VectorXd& rhs,
                                                                          const std::size_t subdomain_num) const
{
    return getCommonSubdomains()[subdomain_num].local_pseudoinverse * rhs;
}

Eigen::MatrixXd
PETScLevelSolverEigenPseudoinverseShellBackend::buildSubdomainPseudoinverse(const Eigen::MatrixXd& local_operator) const
{
    return buildDenseEigenSolveMatrix(getSolverType(),
                                      local_operator,
                                      getSolverThreshold(),
                                      "PETScLevelSolverEigenPseudoinverseShellBackend::buildSubdomainPseudoinverse()");
}

PETScLevelSolverEigenPseudoinverseShellBackend::EigenSubdomainSolverType
PETScLevelSolverEigenPseudoinverseShellBackend::getSolverType() const
{
    return d_solver_type;
}

double
PETScLevelSolverEigenPseudoinverseShellBackend::getSolverThreshold() const
{
    return d_solver_threshold;
}
} // namespace IBTK
