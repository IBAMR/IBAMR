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

    const std::string& getTypeKey() const override;
    void configure(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db) override;
    const char* getPCNameSuffixAdditive() const override;
    const char* getPCNameSuffixMultiplicative() const override;
    void initialize() override;
    void deallocate() override;
    void initializeSubdomainPseudoinverse(const Eigen::MatrixXd& local_operator, std::size_t subdomain_num);
    Eigen::VectorXd solveSubdomainSystem(const Eigen::VectorXd& rhs, std::size_t subdomain_num) const;
    void applyAdditive(Vec x, Vec y) override;
    void applyMultiplicative(Vec x, Vec y) override;

private:
    using EigenSubdomainSolverType = PETScLevelSolverEigenShellBackendBase::EigenSubdomainSolverType;

    EigenSubdomainSolverType getSolverType() const;
    double getSolverThreshold() const;

    Eigen::MatrixXd buildSubdomainPseudoinverse(const Eigen::MatrixXd& local_operator) const;
    std::string d_type_key = "eigen-pseudoinverse";
    EigenSubdomainSolverType d_solver_type = EigenSubdomainSolverType::COL_PIV_HOUSEHOLDER_QR;
    double d_solver_threshold = -1.0;
};
} // namespace IBTK

#endif
