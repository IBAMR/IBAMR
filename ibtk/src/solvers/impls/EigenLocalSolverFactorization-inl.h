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

#ifndef included_IBTK_EigenLocalSolverFactorization_inl
#define included_IBTK_EigenLocalSolverFactorization_inl
namespace IBTK
{
namespace EigenLocalSolverDetail
{
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
} // namespace EigenLocalSolverDetail
} // namespace IBTK
#endif
