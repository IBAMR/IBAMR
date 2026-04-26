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

#ifndef included_IBTK_private_PETScLevelSolverEigenShellBackend
#define included_IBTK_private_PETScLevelSolverEigenShellBackend

#include <ibtk/private/PETScLevelSolverEigenFactorizedLocalSolveShellBackend.h>

namespace IBTK
{
/*!
 * \brief Eigen direct-solve shell smoother backend for PETScLevelSolver.
 *
 * This backend factors each local subdomain operator with a configurable
 * Eigen direct solver and applies the resulting local solves in either
 * additive or multiplicative shell mode.
 */
class PETScLevelSolverEigenShellBackend : public PETScLevelSolverEigenFactorizedLocalSolveShellBackend
{
public:
    static const std::string s_backend_name;

    PETScLevelSolverEigenShellBackend() = default;

    const std::string& getName() const override;
};
} // namespace IBTK

#endif
