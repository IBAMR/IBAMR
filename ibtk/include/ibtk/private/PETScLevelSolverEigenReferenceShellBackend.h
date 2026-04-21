// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_private_PETScLevelSolverEigenReferenceShellBackend
#define included_IBTK_private_PETScLevelSolverEigenReferenceShellBackend

#include <ibtk/private/PETScLevelSolverEigenLocalSolveShellBackend.h>

namespace IBTK
{
/*!
 * \brief Reference multiplicative Eigen shell smoother backend for
 * PETScLevelSolver.
 *
 * This backend recomputes the global residual on each multiplicative step to
 * provide a simple reference implementation for validating the optimized
 * shell smoother backends.
 */
class PETScLevelSolverEigenReferenceShellBackend : public PETScLevelSolverEigenLocalSolveShellBackend
{
public:
    static const std::string s_backend_name;

    PETScLevelSolverEigenReferenceShellBackend() = default;

    using EigenSubdomainSolverType = PETScLevelSolverEigenLocalSolveShellBackend::EigenSubdomainSolverType;

    const std::string& getName() const override;
    void apply(Vec x, Vec y) override;

private:
    struct SolveStorageBase
    {
        virtual ~SolveStorageBase() = default;
    };

    template <class SolverType>
    struct TypedSolveStorage : public SolveStorageBase
    {
        std::vector<SolverType> solvers;
    };

    template <class SolverType>
    TypedSolveStorage<SolverType>& getSolveStorage();

    template <class SolverType>
    const TypedSolveStorage<SolverType>& getSolveStorage() const;

    template <class SolverType>
    void initializeSolveStorage();

    template <class SolverType>
    void initializeTypedSubdomainSolver(const Eigen::MatrixXd& local_operator, std::size_t subdomain_num);

    template <class SolverType>
    Eigen::VectorXd solveTypedSubdomainSystem(const Eigen::VectorXd& rhs, std::size_t subdomain_num) const;

    EigenSubdomainSolverType getSolverType() const;
    double getSolverThreshold() const;
    void configureFromInputDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db) override;
    void initializeLocalSubdomainSolver(const Eigen::MatrixXd& local_operator, std::size_t subdomain_num) override;
    Eigen::VectorXd solveLocalSubdomainSystem(const Eigen::VectorXd& rhs, std::size_t subdomain_num) const override;
    void initializeAdditionalSolverState() override;
    void deallocateAdditionalSolverState() override;

    std::unique_ptr<SolveStorageBase> d_solve_storage;
    Eigen::SparseMatrix<double, Eigen::RowMajor> d_level_mat;
    EigenSubdomainSolverType d_solver_type = EigenSubdomainSolverType::PARTIAL_PIV_LU;
    double d_solver_threshold = -1.0;
};
} // namespace IBTK

#include <ibtk/private/PETScLevelSolverEigenReferenceShellBackend-inl.h>

#endif
