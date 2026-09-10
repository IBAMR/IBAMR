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

#include <ibtk/private/EigenLocalSolver.h>

#include <tbox/Utilities.h>

#include <Eigen/Dense>

#include <algorithm>
#include <cctype>
#include <type_traits>
#include <vector>

namespace IBTK
{
namespace
{
#if EIGEN_VERSION_AT_LEAST(5, 0, 0)
using JacobiSolver = Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV>;
using BDCSolver = Eigen::BDCSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV>;
#else
using JacobiSolver = Eigen::JacobiSVD<Eigen::MatrixXd>;
using BDCSolver = Eigen::BDCSVD<Eigen::MatrixXd>;
#endif

template <class Solver>
class Factorization : public EigenLocalSolver
{
public:
    /*! \brief Compute the factorization and configure its rank threshold. */
    Factorization(const Eigen::MatrixXd& matrix, double threshold);
    /*! \copydoc EigenLocalSolver::solve */
    void solve(Eigen::VectorXd& solution, const Eigen::VectorXd& rhs) const override;
    /*! \copydoc EigenLocalSolver::solveMatrix */
    Eigen::MatrixXd solveMatrix(const Eigen::MatrixXd& rhs) const override;
    /*! \copydoc EigenLocalSolver::getSolveMatrix */
    Eigen::MatrixXd getSolveMatrix() const override;

private:
    Solver d_solver;
};

template <class Solver>
Factorization<Solver>::Factorization(const Eigen::MatrixXd& matrix, const double threshold)
{
    if constexpr (std::is_same_v<Solver, JacobiSolver> || std::is_same_v<Solver, BDCSolver>)
    {
#if EIGEN_VERSION_AT_LEAST(5, 0, 0)
        d_solver.compute(matrix);
#else
        d_solver.compute(matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
#endif
    }
    else
    {
        d_solver.compute(matrix);
    }
    if constexpr (!std::is_same_v<Solver, Eigen::LLT<Eigen::MatrixXd>> &&
                  !std::is_same_v<Solver, Eigen::LDLT<Eigen::MatrixXd>> &&
                  !std::is_same_v<Solver, Eigen::PartialPivLU<Eigen::MatrixXd>> &&
                  !std::is_same_v<Solver, Eigen::HouseholderQR<Eigen::MatrixXd>>)
    {
        if (threshold >= 0.0)
        {
            d_solver.setThreshold(threshold);
        }
    }
    if constexpr (std::is_same_v<Solver, Eigen::LLT<Eigen::MatrixXd>> ||
                  std::is_same_v<Solver, Eigen::LDLT<Eigen::MatrixXd>>)
    {
        if (d_solver.info() != Eigen::Success)
        {
            TBOX_ERROR("Eigen local LLT/LDLT factorization failed.\n");
        }
    }
}
template <class Solver>
void
Factorization<Solver>::solve(Eigen::VectorXd& solution, const Eigen::VectorXd& rhs) const
{
    solution = d_solver.solve(rhs);
}
template <class Solver>
Eigen::MatrixXd
Factorization<Solver>::solveMatrix(const Eigen::MatrixXd& rhs) const
{
    return d_solver.solve(rhs);
}
template <class Solver>
Eigen::MatrixXd
Factorization<Solver>::getSolveMatrix() const
{
    if constexpr (std::is_same_v<Solver, Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd>>)
    {
        return d_solver.pseudoInverse();
    }
    else if constexpr (std::is_same_v<Solver, JacobiSolver> || std::is_same_v<Solver, BDCSolver>)
    {
        Eigen::VectorXd inverse = Eigen::VectorXd::Zero(d_solver.singularValues().size());
        const double cutoff = d_solver.singularValues()[0] * d_solver.threshold();
        for (Eigen::Index k = 0; k < inverse.size(); ++k)
        {
            if (d_solver.singularValues()[k] > cutoff)
            {
                inverse[k] = 1.0 / d_solver.singularValues()[k];
            }
        }
        return d_solver.matrixV() * inverse.asDiagonal() * d_solver.matrixU().adjoint();
    }
    else
    {
        return d_solver.solve(Eigen::MatrixXd::Identity(d_solver.rows(), d_solver.cols()));
    }
}

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
        return std::make_unique<Factorization<Eigen::LLT<Eigen::MatrixXd>>>(matrix, threshold);
    }
    if (normalized == "ldlt")
    {
        return std::make_unique<Factorization<Eigen::LDLT<Eigen::MatrixXd>>>(matrix, threshold);
    }
    if (normalized == "partialpivlu")
    {
        return std::make_unique<Factorization<Eigen::PartialPivLU<Eigen::MatrixXd>>>(matrix, threshold);
    }
    if (normalized == "fullpivlu")
    {
        return std::make_unique<Factorization<Eigen::FullPivLU<Eigen::MatrixXd>>>(matrix, threshold);
    }
    if (normalized == "householderqr")
    {
        return std::make_unique<Factorization<Eigen::HouseholderQR<Eigen::MatrixXd>>>(matrix, threshold);
    }
    if (normalized == "colpivhouseholderqr")
    {
        return std::make_unique<Factorization<Eigen::ColPivHouseholderQR<Eigen::MatrixXd>>>(matrix, threshold);
    }
    if (normalized == "completeorthogonaldecomposition")
    {
        return std::make_unique<Factorization<Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd>>>(matrix,
                                                                                                        threshold);
    }
    if (normalized == "fullpivhouseholderqr")
    {
        return std::make_unique<Factorization<Eigen::FullPivHouseholderQR<Eigen::MatrixXd>>>(matrix, threshold);
    }
    if (normalized == "jacobisvd")
    {
        return std::make_unique<Factorization<JacobiSolver>>(matrix, threshold);
    }
    if (normalized == "bdcsvd")
    {
        return std::make_unique<Factorization<BDCSolver>>(matrix, threshold);
    }
    TBOX_ERROR("Unknown Eigen subdomain solver type: " << type << "\n");
    return nullptr;
}
} // namespace IBTK
