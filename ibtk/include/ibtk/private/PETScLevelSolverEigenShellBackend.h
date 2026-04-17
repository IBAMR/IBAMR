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

#ifndef included_IBTK_private_PETScLevelSolverEigenShellBackend
#define included_IBTK_private_PETScLevelSolverEigenShellBackend

#include <ibtk/private/PETScLevelSolverEigenShellBackendBase.h>

namespace IBTK
{
class PETScLevelSolverEigenShellBackend : public PETScLevelSolverEigenShellBackendBase
{
public:
    explicit PETScLevelSolverEigenShellBackend(PETScLevelSolver& solver);

    using EigenSubdomainSolverType = PETScLevelSolverEigenShellBackendBase::EigenSubdomainSolverType;

    const std::string& getTypeKey() const override;
    void configure(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db) override;
    const char* getPCNameSuffixAdditive() const override;
    const char* getPCNameSuffixMultiplicative() const override;
    void initialize() override;
    void deallocate() override;
    void initializeSubdomainSolver(const Eigen::MatrixXd& local_operator, std::size_t subdomain_num);
    Eigen::VectorXd solveSubdomainSystem(const Eigen::VectorXd& rhs, std::size_t subdomain_num) const;
    void applyAdditive(Vec x, Vec y) override;
    void applyMultiplicative(Vec x, Vec y) override;

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
    void initializeSubdomainSolver(const Eigen::MatrixXd& local_operator, std::size_t subdomain_num);

    template <class SolverType>
    Eigen::VectorXd solveSubdomainSystem(const Eigen::VectorXd& rhs, std::size_t subdomain_num) const;

    template <class SolverType>
    void applyAdditiveImpl(Vec x, Vec y);

    template <class SolverType>
    void applyMultiplicativeImpl(Vec x, Vec y);

    EigenSubdomainSolverType getSolverType() const;
    double getSolverThreshold() const;

    std::unique_ptr<SolveStorageBase> d_solve_storage;
    std::string d_type_key = "eigen";
    EigenSubdomainSolverType d_solver_type = EigenSubdomainSolverType::PARTIAL_PIV_LU;
    double d_solver_threshold = -1.0;
};
} // namespace IBTK

#include <ibtk/private/PETScLevelSolverEigenShellBackend-inl.h>

#endif
