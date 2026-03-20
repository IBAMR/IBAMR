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

#ifndef included_IBTK_private_PETScLevelSolverEigenPseudoinverseShellBackend
#define included_IBTK_private_PETScLevelSolverEigenPseudoinverseShellBackend

#include <ibtk/private/PETScLevelSolverEigenShellBackendBase.h>

namespace IBTK
{
class PETScLevelSolverEigenPseudoinverseShellBackend : public PETScLevelSolverEigenShellBackendBase
{
public:
    explicit PETScLevelSolverEigenPseudoinverseShellBackend(PETScLevelSolver& solver);

    void configure(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);
    void initialize();
    void deallocate();
    void initializeSubdomainPseudoinverse(const Eigen::MatrixXd& local_operator, std::size_t subdomain_num);
    Eigen::VectorXd solveSubdomainSystem(const Eigen::VectorXd& rhs, std::size_t subdomain_num) const;
    void applyAdditive(Vec x, Vec y);
    void applyMultiplicative(Vec x, Vec y);

private:
    using EigenSubdomainSolverType = PETScLevelSolverEigenShellBackendBase::EigenSubdomainSolverType;

    EigenSubdomainSolverType getSolverType() const;
    double getSolverThreshold() const;

    Eigen::MatrixXd buildSubdomainPseudoinverse(const Eigen::MatrixXd& local_operator) const;
    EigenSubdomainSolverType d_solver_type = EigenSubdomainSolverType::COL_PIV_HOUSEHOLDER_QR;
    double d_solver_threshold = -1.0;
};
} // namespace IBTK

#endif
