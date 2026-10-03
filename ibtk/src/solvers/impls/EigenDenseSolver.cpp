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

#include <ibtk/EigenDenseSolver.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/string_utilities.h>

#include <tbox/Utilities.h>

#include <Eigen/Dense>

#include <algorithm>
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
class Factorization : public EigenDenseSolver
{
public:
    /*! \brief Compute the factorization and configure its rank threshold. */
    Factorization(const Eigen::MatrixXd& matrix, double threshold);
    /*! \copydoc EigenDenseSolver::solve */
    void solve(const double* rhs, double* solution) const override;
    /*! \copydoc EigenDenseSolver::solveAllocates */
    bool solveAllocates() const override;
    /*! \copydoc EigenDenseSolver::solveMatrix */
    Eigen::MatrixXd solveMatrix(const Eigen::MatrixXd& rhs) const override;
    /*! \copydoc EigenDenseSolver::getSolveMatrix */
    Eigen::MatrixXd getSolveMatrix() const override;

private:
    Solver d_solver;
};

template <class Solver>
Factorization<Solver>::Factorization(const Eigen::MatrixXd& matrix, const double threshold)
{
    // The threshold must be set before compute(): CompleteOrthogonalDecomposition builds its
    // rank-dependent Z transformation during compute() itself, using whatever threshold is
    // configured at that time. Reusing those factors later with a different notion of rank
    // (e.g. from a threshold set afterward) does not recompute the eliminated coupling, so the
    // action of solve() silently differs from the true rank-truncated pseudoinverse action.
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
    if constexpr (std::is_same_v<Solver, Eigen::LLT<Eigen::MatrixXd>> ||
                  std::is_same_v<Solver, Eigen::LDLT<Eigen::MatrixXd>>)
    {
        if (d_solver.info() != Eigen::Success)
        {
            TBOX_ERROR("make_eigen_dense_solver():\n"
                       << "  the LLT/LDLT factorization failed.\n");
        }
    }
}
template <class Solver>
void
Factorization<Solver>::solve(const double* rhs, double* solution) const
{
    const Eigen::Index n = d_solver.rows();
    Eigen::Map<Eigen::VectorXd>(solution, n) = d_solver.solve(Eigen::Map<const Eigen::VectorXd>(rhs, n));
}
template <class Solver>
bool
Factorization<Solver>::solveAllocates() const
{
    return !(std::is_same_v<Solver, Eigen::LLT<Eigen::MatrixXd>> ||
             std::is_same_v<Solver, Eigen::LDLT<Eigen::MatrixXd>> ||
             std::is_same_v<Solver, Eigen::PartialPivLU<Eigen::MatrixXd>>);
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

} // namespace
void
validate_eigen_solver_type(const std::string& type)
{
    // Construct a solver of the smallest matrix, so that the names that make_eigen_dense_solver() accepts are
    // defined in one place.
    make_eigen_dense_solver(type, Eigen::MatrixXd::Identity(1, 1), -1.0);
}

std::unique_ptr<EigenDenseSolver>
make_eigen_dense_solver(const std::string& type, const Eigen::MatrixXd& matrix, const double threshold)
{
    if (equals_ignore_case(type, "LLT"))
    {
        return std::make_unique<Factorization<Eigen::LLT<Eigen::MatrixXd>>>(matrix, threshold);
    }
    if (equals_ignore_case(type, "LDLT"))
    {
        return std::make_unique<Factorization<Eigen::LDLT<Eigen::MatrixXd>>>(matrix, threshold);
    }
    if (equals_ignore_case(type, "PARTIAL_PIV_LU"))
    {
        return std::make_unique<Factorization<Eigen::PartialPivLU<Eigen::MatrixXd>>>(matrix, threshold);
    }
    if (equals_ignore_case(type, "FULL_PIV_LU"))
    {
        return std::make_unique<Factorization<Eigen::FullPivLU<Eigen::MatrixXd>>>(matrix, threshold);
    }
    if (equals_ignore_case(type, "HOUSEHOLDER_QR"))
    {
        return std::make_unique<Factorization<Eigen::HouseholderQR<Eigen::MatrixXd>>>(matrix, threshold);
    }
    if (equals_ignore_case(type, "COL_PIV_HOUSEHOLDER_QR"))
    {
        return std::make_unique<Factorization<Eigen::ColPivHouseholderQR<Eigen::MatrixXd>>>(matrix, threshold);
    }
    if (equals_ignore_case(type, "COMPLETE_ORTHOGONAL_DECOMPOSITION"))
    {
        return std::make_unique<Factorization<Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd>>>(matrix,
                                                                                                        threshold);
    }
    if (equals_ignore_case(type, "FULL_PIV_HOUSEHOLDER_QR"))
    {
        return std::make_unique<Factorization<Eigen::FullPivHouseholderQR<Eigen::MatrixXd>>>(matrix, threshold);
    }
    if (equals_ignore_case(type, "JACOBI_SVD"))
    {
        return std::make_unique<Factorization<JacobiSolver>>(matrix, threshold);
    }
    if (equals_ignore_case(type, "BDC_SVD"))
    {
        return std::make_unique<Factorization<BDCSolver>>(matrix, threshold);
    }
    TBOX_ERROR("make_eigen_dense_solver():\n"
               << "  unknown Eigen subdomain solver type = " << type
               << "; valid values are LLT, LDLT, PARTIAL_PIV_LU, FULL_PIV_LU, HOUSEHOLDER_QR, COL_PIV_HOUSEHOLDER_QR, "
                  "COMPLETE_ORTHOGONAL_DECOMPOSITION, FULL_PIV_HOUSEHOLDER_QR, JACOBI_SVD, and BDC_SVD.\n");
    return nullptr;
}

Eigen::MatrixXd
make_eigen_matrix(Mat matrix)
{
    PetscInt rows = 0, columns = 0;
    int ierr = MatGetSize(matrix, &rows, &columns);
    IBTK_CHKERRQ(ierr);
    // Only the stored entries of a sparse matrix are visited.
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(rows, columns);
    for (PetscInt row = 0; row < rows; ++row)
    {
        PetscInt count = 0;
        const PetscInt* column_indices = nullptr;
        const PetscScalar* values = nullptr;
        ierr = MatGetRow(matrix, row, &count, &column_indices, &values);
        IBTK_CHKERRQ(ierr);
        for (PetscInt k = 0; k < count; ++k)
        {
            result(row, column_indices[k]) = values[k];
        }
        ierr = MatRestoreRow(matrix, row, &count, &column_indices, &values);
        IBTK_CHKERRQ(ierr);
    }
    return result;
}
} // namespace IBTK
