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

#include <ibtk/IBTK_MPI.h>
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
    Eigen::MatrixXd pseudoinverse;
    dispatchEigenSolverType(
        getSolverType(),
        [this, &local_operator, &pseudoinverse](auto solver_tag)
        {
            using SolverType = typename decltype(solver_tag)::type;
            if constexpr (std::is_same_v<SolverType, Eigen::LLT<Eigen::MatrixXd>>)
            {
                pseudoinverse = buildLLTSolveMatrix(local_operator);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::LDLT<Eigen::MatrixXd>>)
            {
                pseudoinverse = buildLDLTSolveMatrix(local_operator);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::PartialPivLU<Eigen::MatrixXd>>)
            {
                pseudoinverse = buildPartialPivLUSolveMatrix(local_operator);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::FullPivLU<Eigen::MatrixXd>>)
            {
                pseudoinverse = buildFullPivLUSolveMatrix(local_operator, getSolverThreshold());
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::HouseholderQR<Eigen::MatrixXd>>)
            {
                pseudoinverse = buildHouseholderQRSolveMatrix(local_operator);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::ColPivHouseholderQR<Eigen::MatrixXd>>)
            {
                pseudoinverse = buildQRSolveMatrix<SolverType>(local_operator, getSolverThreshold());
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd>>)
            {
                pseudoinverse = buildCompleteOrthogonalDecompositionPseudoinverse(local_operator, getSolverThreshold());
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::FullPivHouseholderQR<Eigen::MatrixXd>>)
            {
                pseudoinverse = buildQRSolveMatrix<SolverType>(local_operator, getSolverThreshold());
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::JacobiSVD<Eigen::MatrixXd>> ||
                               std::is_same_v<SolverType, Eigen::BDCSVD<Eigen::MatrixXd>>)
            {
                pseudoinverse = buildSVDPseudoinverse<SolverType>(local_operator, getSolverThreshold());
            }
            else
            {
                TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix
                                                      << " PETScLevelSolverEigenPseudoinverseShellBackend::"
                                                      << "buildSubdomainPseudoinverse()\n"
                                                      << "Unsupported Eigen subdomain solver type.\n");
            }
        });
    return pseudoinverse;
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
