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

#ifndef included_IBTK_private_PETScLevelSolverEigenReferenceShellBackend
#define included_IBTK_private_PETScLevelSolverEigenReferenceShellBackend

#include <ibtk/private/PETScLevelSolverEigenFactorizedLocalSolveShellBackend.h>

namespace IBTK
{
/*!
 * \brief Reference multiplicative Eigen shell smoother backend for
 * PETScLevelSolver.
 *
 * This backend recomputes the global residual on each multiplicative step to
 * provide a simple reference implementation for validating the optimized
 * shell smoother backends.
 */
class PETScLevelSolverEigenReferenceShellBackend : public PETScLevelSolverEigenFactorizedLocalSolveShellBackend
{
public:
    static const std::string s_backend_name;

    PETScLevelSolverEigenReferenceShellBackend() = default;

    const std::string& getName() const override;

private:
    void initializeAdditionalSolverState() override;
    void deallocateAdditionalSolverState() override;
    void initializeSubdomainSweep(Vec x, Vec y) override;
    void beginSubdomainRhs(std::size_t subdomain_num, Vec x, Vec y) override;
    void updateSubdomainResidual(std::size_t subdomain_num, Vec x, Vec y) override;

    Eigen::SparseMatrix<double, Eigen::RowMajor> d_level_mat;
};
} // namespace IBTK

#endif
