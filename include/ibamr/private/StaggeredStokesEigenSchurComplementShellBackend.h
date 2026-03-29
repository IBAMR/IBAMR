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

#ifndef included_IBAMR_private_StaggeredStokesEigenSchurComplementShellBackend
#define included_IBAMR_private_StaggeredStokesEigenSchurComplementShellBackend

#include <ibamr/StaggeredStokesPETScLevelSolver.h>

#include <ibtk/private/PETScLevelSolverEigenShellBackendBase.h>

#include <tbox/Pointer.h>

#include <petscvec.h>

#include <memory>
#include <string>
#include <vector>

namespace IBAMR
{
class StaggeredStokesEigenSchurComplementShellBackend : public IBTK::PETScLevelSolverEigenShellBackendBase
{
public:
    explicit StaggeredStokesEigenSchurComplementShellBackend(StaggeredStokesPETScLevelSolver& solver,
                                                             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    const std::string& getTypeKey() const override;

    void configure(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db) override;

    const char* getPCNameSuffixAdditive() const override;

    const char* getPCNameSuffixMultiplicative() const override;

    void reset();

    void initialize() override;

    void deallocate() override;

    void applyAdditive(Vec x, Vec y) override;

    void applyMultiplicative(Vec x, Vec y) override;

private:
    using EigenSubdomainSolverType = PETScLevelSolverEigenShellBackendBase::EigenSubdomainSolverType;

    struct CustomEigenSchurSubdomainCache
    {
        int overlap_size = 0;
        const std::vector<int>* overlap_dofs = nullptr;
        const std::vector<int>* update_dofs = nullptr;
        const std::vector<int>* update_local_positions = nullptr;
        const std::vector<int>* active_residual_update_rows = nullptr;
        std::vector<int> velocity_positions;
        std::vector<int> pressure_positions;
        const Eigen::SparseMatrix<double, Eigen::RowMajor>* active_residual_update_mat = nullptr;
        Eigen::MatrixXd A00;
        Eigen::MatrixXd A01;
        Eigen::MatrixXd A10;
        Eigen::MatrixXd A11;
        Eigen::MatrixXd A00_inv_A01;
        Eigen::MatrixXd schur;
        Eigen::MatrixXd schur_solve_matrix;
        Eigen::VectorXd rhs_workspace;
        Eigen::VectorXd velocity_rhs_workspace;
        Eigen::VectorXd pressure_rhs_workspace;
        Eigen::VectorXd delta_workspace;
        Eigen::VectorXd residual_input_workspace;
        Eigen::VectorXd residual_delta_workspace;
        Eigen::VectorXd velocity_solution_workspace;
        Eigen::VectorXd pressure_solution_workspace;
    };

    struct CustomEigenA00SolverStorageBase
    {
        virtual ~CustomEigenA00SolverStorageBase() = default;
    };

    template <class SolverType>
    struct CustomEigenA00TypedSolveStorage : public CustomEigenA00SolverStorageBase
    {
        std::vector<SolverType> solvers;
    };

    template <class SolverType>
    CustomEigenA00TypedSolveStorage<SolverType>& getCustomEigenA00SolveStorage();

    template <class SolverType>
    const CustomEigenA00TypedSolveStorage<SolverType>& getCustomEigenA00SolveStorage() const;

    template <class SolverType>
    void initializeCustomEigenA00SolveStorage(std::size_t n_subdomains);

    template <class SolverType>
    void solveCustomEigenSubdomain(CustomEigenSchurSubdomainCache& custom_cache, const SolverType& a00_solver) const;

    template <class SolverType>
    void applyAdditiveImpl(Vec x, Vec y);

    template <class SolverType>
    void applyMultiplicativeImpl(Vec x, Vec y);

    EigenSubdomainSolverType parseBuiltinSolverType(const std::string& type) const;

    Eigen::MatrixXd buildSchurSolveMatrix(const Eigen::MatrixXd& schur) const;

    StaggeredStokesPETScLevelSolver& d_stokes_solver;
    std::vector<CustomEigenSchurSubdomainCache> d_subdomain_caches;
    std::unique_ptr<CustomEigenA00SolverStorageBase> d_a00_solver_storage;
    std::string d_type_key = "eigen-schur-complement";
    EigenSubdomainSolverType d_a00_solver_type = EigenSubdomainSolverType::FULL_PIV_HOUSEHOLDER_QR;
    double d_a00_solver_threshold = -1.0;
    EigenSubdomainSolverType d_schur_solver_type = EigenSubdomainSolverType::FULL_PIV_HOUSEHOLDER_QR;
    double d_schur_solver_threshold = -1.0;
};
} // namespace IBAMR

#include <ibamr/private/StaggeredStokesEigenSchurComplementShellBackend-inl.h>

#endif
