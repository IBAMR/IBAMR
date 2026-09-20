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

#ifndef included_IBTK_EigenDenseSolver
#define included_IBTK_EigenDenseSolver

#include <ibtk/config.h>

#include <petscmat.h>

#include <Eigen/Core>

#include <memory>
#include <string>

namespace IBTK
{
/*! \brief Dense factorization of a square matrix, shared by the Eigen subdomain solvers. */
class EigenDenseSolver
{
public:
    /*! \brief Release the factorization. */
    virtual ~EigenDenseSolver() = default;
    /*! \brief Solve for one right-hand side, reading and writing arrays of the order of the matrix. */
    virtual void solve(const double* rhs, double* solution) const = 0;
    /*!
     * \brief Return whether solve() allocates a temporary vector, as Eigen's QR, SVD and full-pivoting LU
     * solves do. LLT, LDLT and partial-pivoting LU solve in the output.
     */
    virtual bool solveAllocates() const = 0;
    /*! \brief Solve multiple right-hand sides during setup. */
    virtual Eigen::MatrixXd solveMatrix(const Eigen::MatrixXd& rhs) const = 0;
    /*! \brief Form the configured solve matrix, with mode-specific rank handling. */
    virtual Eigen::MatrixXd getSolveMatrix() const = 0;
};

/*! \brief Require one of the canonical solver names, in any case. */
void validate_eigen_solver_type(const std::string& type);

/*!
 * \brief Construct a named Eigen factorization; negative thresholds use Eigen defaults.
 *
 * The matrix must be nonempty and square. LLT and LDLT require their usual
 * definiteness or symmetry assumptions. COMPLETE_ORTHOGONAL_DECOMPOSITION and SVD
 * solve matrices are Moore-Penrose pseudoinverses; other modes solve against the
 * identity and retain their own pivot and rank policy.
 */
std::unique_ptr<EigenDenseSolver>
make_eigen_dense_solver(const std::string& type, const Eigen::MatrixXd& matrix, double threshold);

/*! \brief Copy a sequential PETSc matrix into a dense Eigen matrix. */
Eigen::MatrixXd make_eigen_matrix(Mat matrix);
} // namespace IBTK
#endif
