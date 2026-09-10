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

#ifndef included_IBTK_EigenLocalSolverFactorization
#define included_IBTK_EigenLocalSolverFactorization

#include <ibtk/config.h>

#include <ibtk/private/EigenLocalSolver.h>

#include <tbox/Utilities.h>

#include <Eigen/Dense>

#include <type_traits>

namespace IBTK
{
namespace EigenLocalSolverDetail
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
} // namespace EigenLocalSolverDetail
} // namespace IBTK
#include <EigenLocalSolverFactorization-inl.h>
#endif
