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
 * This backend retains the reference Eigen local factorization choices while
 * using the same live-operator correction-composition path as every other
 * backend.
 */
class PETScLevelSolverEigenReferenceShellBackend : public PETScLevelSolverEigenFactorizedLocalSolveShellBackend
{
public:
    static const std::string s_backend_name;

    PETScLevelSolverEigenReferenceShellBackend() = default;

    const std::string& getName() const override;

private:
    void initializeSubdomainSweep(Vec x, Vec y) override;
};
} // namespace IBTK

#endif
