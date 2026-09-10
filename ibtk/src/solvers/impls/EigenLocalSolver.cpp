// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <EigenLocalSolverFactorization.h>

#include <algorithm>
#include <cctype>
#include <vector>

namespace IBTK
{
namespace
{
std::string
normalize_solver_type(const std::string& type)
{
    std::string normalized;
    for (const unsigned char c : type)
    {
        if (c != '-' && c != '_')
        {
            normalized.push_back(static_cast<char>(std::tolower(c)));
        }
    }
    if (normalized == "cod")
    {
        normalized = "completeorthogonaldecomposition";
    }
    return normalized;
}
} // namespace
void
validate_eigen_local_solver_type(const std::string& type)
{
    const std::string normalized = normalize_solver_type(type);
    const std::vector<std::string> names = { "llt",
                                             "ldlt",
                                             "partialpivlu",
                                             "fullpivlu",
                                             "householderqr",
                                             "colpivhouseholderqr",
                                             "completeorthogonaldecomposition",
                                             "fullpivhouseholderqr",
                                             "jacobisvd",
                                             "bdcsvd" };
    if (std::find(names.begin(), names.end(), normalized) == names.end())
    {
        TBOX_ERROR("Unknown Eigen subdomain solver type: " << type << "\n");
    }
}
std::unique_ptr<EigenLocalSolver>
make_eigen_local_solver(const std::string& type, const Eigen::MatrixXd& matrix, const double threshold)
{
    const std::string normalized = normalize_solver_type(type);
    if (normalized == "llt")
    {
        return std::make_unique<EigenLocalSolverDetail::Factorization<Eigen::LLT<Eigen::MatrixXd>>>(matrix, threshold);
    }
    if (normalized == "ldlt")
    {
        return std::make_unique<EigenLocalSolverDetail::Factorization<Eigen::LDLT<Eigen::MatrixXd>>>(matrix, threshold);
    }
    if (normalized == "partialpivlu")
    {
        return std::make_unique<EigenLocalSolverDetail::Factorization<Eigen::PartialPivLU<Eigen::MatrixXd>>>(matrix,
                                                                                                             threshold);
    }
    if (normalized == "fullpivlu")
    {
        return std::make_unique<EigenLocalSolverDetail::Factorization<Eigen::FullPivLU<Eigen::MatrixXd>>>(matrix,
                                                                                                          threshold);
    }
    if (normalized == "householderqr")
    {
        return std::make_unique<EigenLocalSolverDetail::Factorization<Eigen::HouseholderQR<Eigen::MatrixXd>>>(
            matrix, threshold);
    }
    if (normalized == "colpivhouseholderqr")
    {
        return std::make_unique<EigenLocalSolverDetail::Factorization<Eigen::ColPivHouseholderQR<Eigen::MatrixXd>>>(
            matrix, threshold);
    }
    if (normalized == "completeorthogonaldecomposition")
    {
        return std::make_unique<
            EigenLocalSolverDetail::Factorization<Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd>>>(matrix,
                                                                                                            threshold);
    }
    if (normalized == "fullpivhouseholderqr")
    {
        return std::make_unique<EigenLocalSolverDetail::Factorization<Eigen::FullPivHouseholderQR<Eigen::MatrixXd>>>(
            matrix, threshold);
    }
    if (normalized == "jacobisvd")
    {
        return std::make_unique<EigenLocalSolverDetail::Factorization<EigenLocalSolverDetail::JacobiSolver>>(matrix,
                                                                                                             threshold);
    }
    if (normalized == "bdcsvd")
    {
        return std::make_unique<EigenLocalSolverDetail::Factorization<EigenLocalSolverDetail::BDCSolver>>(matrix,
                                                                                                          threshold);
    }
    TBOX_ERROR("Unknown Eigen subdomain solver type: " << type << "\n");
    return nullptr;
}
} // namespace IBTK
