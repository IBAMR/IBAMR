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

#ifndef included_IBTK_private_PETScLevelSolverEigenShellBackendBase
#define included_IBTK_private_PETScLevelSolverEigenShellBackendBase

#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScLevelSolver.h>
#include <ibtk/private/PETScLevelSolverEigenShellBackendCommon.h>
#include <ibtk/private/PETScLevelSolverShellBackend.h>

#include <petscvec.h>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <Eigen/Sparse>

#include <type_traits>
#include <unordered_map>

namespace IBTK
{
/*!
 * \brief Common Eigen-based shell smoother support for PETScLevelSolver.
 *
 * This base class caches subdomain metadata shared by the concrete Eigen
 * backends and provides helper utilities for translating PETSc matrices and
 * vectors into the Eigen data structures used by those backends.
 */
class PETScLevelSolverEigenShellBackendBase : public PETScLevelSolverShellBackend
{
public:
    PETScLevelSolverEigenShellBackendBase() = default;

protected:
    using EigenSubdomainSolverType = PETScLevelSolverEigenShell::SubdomainSolverType;

    struct CommonSubdomainCache
    {
        std::vector<int> overlap_dofs;
        std::vector<int> nonoverlap_dofs;
        std::vector<int> nonoverlap_local_positions;
        std::vector<int> update_dofs;
        std::vector<int> update_local_positions;
        Eigen::MatrixXd local_pseudoinverse;
        Eigen::VectorXd rhs_workspace;
        Eigen::VectorXd delta_workspace;
    };

    template <class Handler>
    void dispatchEigenSolverType(EigenSubdomainSolverType solver_type, Handler&& handler) const;

    EigenSubdomainSolverType parseEigenSubdomainSolverType(const std::string& type, const char* caller) const;

    template <class SolverType>
    static void initializeEigenSolver(SolverType& solver, const Eigen::MatrixXd& local_operator, double threshold);

    template <class SolverType, class RhsType>
    static auto solveEigenSystem(const SolverType& solver, const RhsType& rhs) -> decltype(solver.solve(rhs));

    template <class SolverType>
    static Eigen::MatrixXd buildQRSolveMatrix(const Eigen::MatrixXd& matrix, double threshold);

    static Eigen::SparseMatrix<double, Eigen::RowMajor> copyPETScMatToEigenSparse(Mat mat);

    static Eigen::MatrixXd buildLLTSolveMatrix(const Eigen::MatrixXd& matrix)
    {
        Eigen::LLT<Eigen::MatrixXd> solver(matrix);
        return solver.solve(Eigen::MatrixXd::Identity(matrix.rows(), matrix.cols()));
    }

    static Eigen::MatrixXd buildLDLTSolveMatrix(const Eigen::MatrixXd& matrix)
    {
        Eigen::LDLT<Eigen::MatrixXd> solver(matrix);
        return solver.solve(Eigen::MatrixXd::Identity(matrix.rows(), matrix.cols()));
    }

    static Eigen::MatrixXd buildHouseholderQRSolveMatrix(const Eigen::MatrixXd& matrix)
    {
        Eigen::HouseholderQR<Eigen::MatrixXd> solver(matrix);
        return solver.solve(Eigen::MatrixXd::Identity(matrix.rows(), matrix.cols()));
    }

    static Eigen::MatrixXd buildPartialPivLUSolveMatrix(const Eigen::MatrixXd& matrix)
    {
        Eigen::PartialPivLU<Eigen::MatrixXd> solver(matrix);
        return solver.solve(Eigen::MatrixXd::Identity(matrix.rows(), matrix.cols()));
    }

    static Eigen::MatrixXd buildFullPivLUSolveMatrix(const Eigen::MatrixXd& matrix, double threshold)
    {
        Eigen::FullPivLU<Eigen::MatrixXd> solver(matrix);
        if (threshold >= 0.0) solver.setThreshold(threshold);
        return solver.solve(Eigen::MatrixXd::Identity(matrix.rows(), matrix.cols()));
    }

    template <class SVDType>
    static Eigen::MatrixXd buildSVDPseudoinverse(const Eigen::MatrixXd& matrix, double threshold);

    static void checkSerialEigenShellBackend(const PETScLevelSolverShellBackendState& solver_state, const char* caller);

    Eigen::MatrixXd buildDenseEigenSolveMatrix(EigenSubdomainSolverType solver_type,
                                               const Eigen::MatrixXd& matrix,
                                               double threshold,
                                               const char* caller) const;

    static Eigen::MatrixXd buildCompleteOrthogonalDecompositionPseudoinverse(const Eigen::MatrixXd& matrix,
                                                                             double threshold)
    {
        Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> solver(matrix);
        if (threshold >= 0.0) solver.setThreshold(threshold);
        return solver.pseudoInverse();
    }

    void clearCommonData()
    {
        d_n_dofs = 0;
        d_common_subdomains.clear();
        d_sweep_y_values = nullptr;
    }

    Eigen::Index getNumDofs() const
    {
        return d_n_dofs;
    }

    std::vector<CommonSubdomainCache>& getCommonSubdomains()
    {
        return d_common_subdomains;
    }

    const std::vector<CommonSubdomainCache>& getCommonSubdomains() const
    {
        return d_common_subdomains;
    }

    template <class InitializeSubdomainSolver>
    void initializeCommonDataWithLocalOperatorHook(InitializeSubdomainSolver initialize_subdomain_solver);

    std::size_t getNumberOfSubdomains() const override;
    void initializeSubdomainSweep(Vec x, Vec y) override;
    void beginSubdomainRhs(std::size_t subdomain_num, Vec x, Vec y) override;
    void endSubdomainRhs(std::size_t subdomain_num, Vec x, Vec y) override;
    void accumulateSubdomainCorrection(std::size_t subdomain_num, Vec y) override;
    const std::vector<int>& getSubdomainCorrectionDofs(std::size_t subdomain_num) const override;
    void copySubdomainCorrection(std::size_t subdomain_num, PetscScalar* correction_values) override;
    void finalizeSubdomainSweep(Vec x, Vec y) override;

    const PETScLevelSolverShellBackendState& getSolverState() const
    {
        return d_solver_state;
    }

    Eigen::Index d_n_dofs = 0;
    std::vector<CommonSubdomainCache> d_common_subdomains;
    PetscScalar* d_sweep_y_values = nullptr;

private:
    template <class SolverType>
    struct eigen_solver_supports_threshold : std::false_type
    {
    };

    template <class SolverType>
    struct eigen_solver_uses_svd_compute : std::false_type
    {
    };
};

template <>
struct PETScLevelSolverEigenShellBackendBase::eigen_solver_supports_threshold<Eigen::FullPivLU<Eigen::MatrixXd>>
    : std::true_type
{
};

template <>
struct PETScLevelSolverEigenShellBackendBase::eigen_solver_supports_threshold<
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd>> : std::true_type
{
};

template <>
struct PETScLevelSolverEigenShellBackendBase::eigen_solver_supports_threshold<
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd>> : std::true_type
{
};

template <>
struct PETScLevelSolverEigenShellBackendBase::eigen_solver_supports_threshold<
    Eigen::FullPivHouseholderQR<Eigen::MatrixXd>> : std::true_type
{
};

template <>
struct PETScLevelSolverEigenShellBackendBase::eigen_solver_supports_threshold<Eigen::JacobiSVD<Eigen::MatrixXd>>
    : std::true_type
{
};

template <>
struct PETScLevelSolverEigenShellBackendBase::eigen_solver_supports_threshold<Eigen::BDCSVD<Eigen::MatrixXd>>
    : std::true_type
{
};

template <>
struct PETScLevelSolverEigenShellBackendBase::eigen_solver_uses_svd_compute<Eigen::JacobiSVD<Eigen::MatrixXd>>
    : std::true_type
{
};

template <>
struct PETScLevelSolverEigenShellBackendBase::eigen_solver_uses_svd_compute<Eigen::BDCSVD<Eigen::MatrixXd>>
    : std::true_type
{
};

} // namespace IBTK

#include <ibtk/private/PETScLevelSolverEigenShellBackendBase-inl.h>

#endif
