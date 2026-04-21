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

#ifndef included_IBTK_private_PETScLevelSolverEigenPseudoinverseShellBackend
#define included_IBTK_private_PETScLevelSolverEigenPseudoinverseShellBackend

#include <ibtk/private/PETScLevelSolverEigenLocalSolveShellBackend.h>

namespace IBTK
{
/*!
 * \brief Eigen pseudoinverse shell smoother backend for PETScLevelSolver.
 *
 * This backend explicitly forms a dense pseudoinverse for each local
 * subdomain operator so that shell applications can reuse that matrix during
 * additive or multiplicative sweeps.
 */
class PETScLevelSolverEigenPseudoinverseShellBackend : public PETScLevelSolverEigenLocalSolveShellBackend
{
public:
    static const std::string s_backend_name;

    PETScLevelSolverEigenPseudoinverseShellBackend() = default;

    const std::string& getName() const override;

private:
    using EigenSubdomainSolverType = PETScLevelSolverEigenLocalSolveShellBackend::EigenSubdomainSolverType;

    EigenSubdomainSolverType getSolverType() const;
    double getSolverThreshold() const;
    void configureFromInputDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db) override;
    void initializeLocalSubdomainSolver(const Eigen::MatrixXd& local_operator, std::size_t subdomain_num) override;
    Eigen::VectorXd solveLocalSubdomainSystem(const Eigen::VectorXd& rhs, std::size_t subdomain_num) const override;

    Eigen::MatrixXd buildSubdomainPseudoinverse(const Eigen::MatrixXd& local_operator) const;
    EigenSubdomainSolverType d_solver_type = EigenSubdomainSolverType::COL_PIV_HOUSEHOLDER_QR;
    double d_solver_threshold = -1.0;
};
} // namespace IBTK

#endif
