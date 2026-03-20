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

#include <petscvec.h>

#include <type_traits>
#include <unordered_map>

namespace IBTK
{
class PETScLevelSolverEigenShellBackendBase
{
public:
    explicit PETScLevelSolverEigenShellBackendBase(PETScLevelSolver& solver) : d_solver(solver)
    {
    }

protected:
    using EigenSubdomainSolverType = PETScLevelSolverEigenShell::SubdomainSolverType;

    struct CommonSubdomainCache
    {
        std::vector<int> overlap_dofs;
        std::vector<int> nonoverlap_dofs;
        std::vector<int> nonoverlap_local_positions;
        std::vector<int> update_dofs;
        std::vector<int> update_local_positions;
        std::vector<int> active_residual_update_rows;
        Eigen::SparseMatrix<double, Eigen::RowMajor> active_residual_update_mat;
        Eigen::MatrixXd local_pseudoinverse;
        Eigen::VectorXd rhs_workspace;
        Eigen::VectorXd delta_workspace;
        Eigen::VectorXd residual_input_workspace;
        Eigen::VectorXd residual_delta_workspace;
    };

    class ConstPetscVecArrayMap
    {
    public:
        ConstPetscVecArrayMap(Vec vec, Eigen::Index size) : d_vec(vec), d_size(size)
        {
            const int ierr = VecGetArrayRead(d_vec, &d_array);
            IBTK_CHKERRQ(ierr);
        }

        ~ConstPetscVecArrayMap()
        {
            const int ierr = VecRestoreArrayRead(d_vec, &d_array);
            IBTK_CHKERRQ(ierr);
        }

        Eigen::Map<const Eigen::VectorXd> getMap() const
        {
            return Eigen::Map<const Eigen::VectorXd>(reinterpret_cast<const double*>(d_array), d_size);
        }

    private:
        Vec d_vec = nullptr;
        Eigen::Index d_size = 0;
        const PetscScalar* d_array = nullptr;
    };

    class PetscVecArrayMap
    {
    public:
        PetscVecArrayMap(Vec vec, Eigen::Index size) : d_vec(vec), d_size(size)
        {
            const int ierr = VecGetArray(d_vec, &d_array);
            IBTK_CHKERRQ(ierr);
        }

        ~PetscVecArrayMap()
        {
            const int ierr = VecRestoreArray(d_vec, &d_array);
            IBTK_CHKERRQ(ierr);
        }

        Eigen::Map<Eigen::VectorXd> getMap() const
        {
            return Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double*>(d_array), d_size);
        }

    private:
        Vec d_vec = nullptr;
        Eigen::Index d_size = 0;
        PetscScalar* d_array = nullptr;
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

    PETScLevelSolver& d_solver;
    Eigen::Index d_n_dofs = 0;
    std::vector<CommonSubdomainCache> d_common_subdomains;

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
