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

#ifndef included_IBTK_private_EigenLocalSolver
#define included_IBTK_private_EigenLocalSolver

#include <ibtk/config.h>

#include <Eigen/Core>

#include <memory>
#include <string>

namespace IBTK
{
/*! \brief Dense local factorization shared by Eigen shell backends. */
class EigenLocalSolver
{
public:
    /*! \brief Release the factorization. */
    virtual ~EigenLocalSolver() = default;
    /*! \brief Solve into an already sized vector. */
    virtual void solve(Eigen::VectorXd& solution, const Eigen::VectorXd& rhs) const = 0;
    /*! \brief Solve multiple right-hand sides during setup. */
    virtual Eigen::MatrixXd solveMatrix(const Eigen::MatrixXd& rhs) const = 0;
    /*! \brief Form the configured solve matrix, with mode-specific rank handling. */
    virtual Eigen::MatrixXd getSolveMatrix() const = 0;
};
/*! \brief Validate a solver name, preserving Eigen shell spelling aliases. */
void validate_eigen_local_solver_type(const std::string& type);
/*! \brief Construct a named Eigen factorization; negative thresholds use Eigen defaults.
 *
 * The matrix must be nonempty and square.
 * LLT and LDLT require their usual definiteness/symmetry assumptions. COD and
 * SVD solve matrices are Moore-Penrose pseudoinverses; other modes solve against
 * the identity and retain their own pivot/rank policy.
 */
std::unique_ptr<EigenLocalSolver>
make_eigen_local_solver(const std::string& type, const Eigen::MatrixXd& matrix, double threshold);
} // namespace IBTK
#endif
