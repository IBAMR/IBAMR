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

#ifndef included_IBTK_private_PETScLevelSolverEigenShellBackendCommon_inl
#define included_IBTK_private_PETScLevelSolverEigenShellBackendCommon_inl

namespace IBTK
{
namespace PETScLevelSolverEigenShell
{
inline SubdomainSolverType
parseSubdomainSolverType(const std::string& type,
                         const std::string& object_name,
                         const std::string& options_prefix,
                         const char* caller)
{
    std::string normalized = type;
    std::transform(normalized.begin(),
                   normalized.end(),
                   normalized.begin(),
                   [](const unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (normalized == "llt")
    {
        return SubdomainSolverType::LLT;
    }
    if (normalized == "ldlt")
    {
        return SubdomainSolverType::LDLT;
    }
    if (normalized == "partialpivlu" || normalized == "partial-piv-lu" || normalized == "partial_piv_lu")
    {
        return SubdomainSolverType::PARTIAL_PIV_LU;
    }
    if (normalized == "fullpivlu" || normalized == "full-piv-lu" || normalized == "full_piv_lu")
    {
        return SubdomainSolverType::FULL_PIV_LU;
    }
    if (normalized == "householderqr" || normalized == "householder-qr" || normalized == "householder_qr")
    {
        return SubdomainSolverType::HOUSEHOLDER_QR;
    }
    if (normalized == "colpivhouseholderqr" || normalized == "col-piv-householder-qr" ||
        normalized == "col_piv_householder_qr")
    {
        return SubdomainSolverType::COL_PIV_HOUSEHOLDER_QR;
    }
    if (normalized == "completeorthogonaldecomposition" || normalized == "complete-orthogonal-decomposition" ||
        normalized == "complete_orthogonal_decomposition" || normalized == "cod")
    {
        return SubdomainSolverType::COMPLETE_ORTHOGONAL_DECOMPOSITION;
    }
    if (normalized == "fullpivhouseholderqr" || normalized == "full-piv-householder-qr" ||
        normalized == "full_piv_householder_qr")
    {
        return SubdomainSolverType::FULL_PIV_HOUSEHOLDER_QR;
    }
    if (normalized == "jacobisvd" || normalized == "jacobi-svd" || normalized == "jacobi_svd")
    {
        return SubdomainSolverType::JACOBI_SVD;
    }
    if (normalized == "bdcsvd" || normalized == "bdc-svd" || normalized == "bdc_svd")
    {
        return SubdomainSolverType::BDC_SVD;
    }

    TBOX_ERROR(object_name << " " << options_prefix << " " << caller << "\n"
                           << "Unknown Eigen subdomain solver type: " << type << std::endl);
    return SubdomainSolverType::JACOBI_SVD;
}

template <class Handler>
inline void
dispatchSubdomainSolverType(const std::string& object_name,
                            const std::string& options_prefix,
                            const SubdomainSolverType solver_type,
                            Handler&& handler)
{
    switch (solver_type)
    {
    case SubdomainSolverType::LLT:
        std::forward<Handler>(handler)(SolverTag<Eigen::LLT<Eigen::MatrixXd>>{});
        return;
    case SubdomainSolverType::LDLT:
        std::forward<Handler>(handler)(SolverTag<Eigen::LDLT<Eigen::MatrixXd>>{});
        return;
    case SubdomainSolverType::PARTIAL_PIV_LU:
        std::forward<Handler>(handler)(SolverTag<Eigen::PartialPivLU<Eigen::MatrixXd>>{});
        return;
    case SubdomainSolverType::FULL_PIV_LU:
        std::forward<Handler>(handler)(SolverTag<Eigen::FullPivLU<Eigen::MatrixXd>>{});
        return;
    case SubdomainSolverType::HOUSEHOLDER_QR:
        std::forward<Handler>(handler)(SolverTag<Eigen::HouseholderQR<Eigen::MatrixXd>>{});
        return;
    case SubdomainSolverType::COL_PIV_HOUSEHOLDER_QR:
        std::forward<Handler>(handler)(SolverTag<Eigen::ColPivHouseholderQR<Eigen::MatrixXd>>{});
        return;
    case SubdomainSolverType::COMPLETE_ORTHOGONAL_DECOMPOSITION:
        std::forward<Handler>(handler)(SolverTag<Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd>>{});
        return;
    case SubdomainSolverType::FULL_PIV_HOUSEHOLDER_QR:
        std::forward<Handler>(handler)(SolverTag<Eigen::FullPivHouseholderQR<Eigen::MatrixXd>>{});
        return;
    case SubdomainSolverType::JACOBI_SVD:
        std::forward<Handler>(handler)(SolverTag<Eigen::JacobiSVD<Eigen::MatrixXd>>{});
        return;
    case SubdomainSolverType::BDC_SVD:
        std::forward<Handler>(handler)(SolverTag<Eigen::BDCSVD<Eigen::MatrixXd>>{});
        return;
    }
}
} // namespace PETScLevelSolverEigenShell
} // namespace IBTK

#endif
