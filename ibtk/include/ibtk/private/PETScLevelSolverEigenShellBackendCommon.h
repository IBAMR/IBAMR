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

#ifndef included_IBTK_private_PETScLevelSolverEigenShellBackendCommon
#define included_IBTK_private_PETScLevelSolverEigenShellBackendCommon

#include <tbox/Utilities.h>

#include <petscvec.h>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/SVD>

#include <algorithm>
#include <cctype>
#include <string>
#include <utility>

namespace IBTK
{
namespace PETScLevelSolverEigenShell
{
enum class SubdomainSolverType
{
    LLT,
    LDLT,
    PARTIAL_PIV_LU,
    FULL_PIV_LU,
    HOUSEHOLDER_QR,
    COL_PIV_HOUSEHOLDER_QR,
    COMPLETE_ORTHOGONAL_DECOMPOSITION,
    FULL_PIV_HOUSEHOLDER_QR,
    JACOBI_SVD,
    BDC_SVD
};

template <class SolverType>
struct SolverTag
{
    using type = SolverType;
};

SubdomainSolverType parseSubdomainSolverType(const std::string& type,
                                             const std::string& object_name,
                                             const std::string& options_prefix,
                                             const char* caller);

template <class Handler>
void dispatchSubdomainSolverType(const std::string& object_name,
                                 const std::string& options_prefix,
                                 SubdomainSolverType solver_type,
                                 Handler&& handler);
} // namespace PETScLevelSolverEigenShell
} // namespace IBTK

#include <ibtk/private/PETScLevelSolverEigenShellBackendCommon-inl.h>

#endif
