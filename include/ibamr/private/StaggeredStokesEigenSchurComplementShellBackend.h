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
/*!
 * \brief Eigen shell smoother backend that applies subdomain Schur-complement
 * solves for staggered Stokes operators.
 *
 * This backend splits each cached subdomain into velocity and pressure blocks,
 * factors the local A00 operator, and supplies the corresponding one-patch
 * Schur solve to the backend-independent correction composer.
 */
class StaggeredStokesEigenSchurComplementShellBackend : public IBTK::PETScLevelSolverEigenShellBackendBase
{
public:
    static const std::string s_backend_name;

    explicit StaggeredStokesEigenSchurComplementShellBackend(StaggeredStokesPETScLevelSolver& solver);

    const std::string& getName() const override;

    void initializeSolverState(const IBTK::PETScLevelSolverShellBackendState& solver_state) override;

    void deallocateSolverState() override;

private:
    using EigenSubdomainSolverType = PETScLevelSolverEigenShellBackendBase::EigenSubdomainSolverType;

    struct CustomEigenSchurSubdomainCache
    {
        std::vector<int> velocity_positions;
        std::vector<int> pressure_positions;
        Eigen::MatrixXd A00;
        Eigen::MatrixXd A01;
        Eigen::MatrixXd A10;
        Eigen::MatrixXd A11;
        Eigen::MatrixXd A00_inv_A01;
        Eigen::MatrixXd schur;
        Eigen::MatrixXd schur_solve_matrix;
        Eigen::VectorXd velocity_rhs_workspace;
        Eigen::VectorXd pressure_rhs_workspace;
        Eigen::VectorXd velocity_solution_workspace;
        Eigen::VectorXd pressure_solution_workspace;
    };

    struct CustomEigenA00SolverStorageBase
    {
        virtual ~CustomEigenA00SolverStorageBase() = default;

        virtual void solveSubdomain(const StaggeredStokesEigenSchurComplementShellBackend& backend,
                                    CustomEigenSchurSubdomainCache& custom_cache,
                                    const Eigen::VectorXd& rhs,
                                    Eigen::VectorXd& delta,
                                    std::size_t subdomain_num) const = 0;
    };

    template <class SolverType>
    struct CustomEigenA00TypedSolveStorage : public CustomEigenA00SolverStorageBase
    {
        void solveSubdomain(const StaggeredStokesEigenSchurComplementShellBackend& backend,
                            CustomEigenSchurSubdomainCache& custom_cache,
                            const Eigen::VectorXd& rhs,
                            Eigen::VectorXd& delta,
                            std::size_t subdomain_num) const override
        {
            backend.solveCustomEigenSubdomain(custom_cache, rhs, delta, solvers[subdomain_num]);
        }

        std::vector<SolverType> solvers;
    };

    template <class SolverType>
    CustomEigenA00TypedSolveStorage<SolverType>& getCustomEigenA00SolveStorage();

    template <class SolverType>
    const CustomEigenA00TypedSolveStorage<SolverType>& getCustomEigenA00SolveStorage() const;

    template <class SolverType>
    void initializeCustomEigenA00SolveStorage(std::size_t n_subdomains);

    template <class SolverType>
    void
    initializeCustomEigenA00Solver(SolverType& solver, const Eigen::MatrixXd& matrix, std::size_t subdomain_num) const;

    template <class SolverType>
    void solveCustomEigenSubdomain(CustomEigenSchurSubdomainCache& custom_cache,
                                   const Eigen::VectorXd& rhs,
                                   Eigen::VectorXd& delta,
                                   const SolverType& a00_solver) const;

    void solveSubdomain(std::size_t subdomain_num) override;

    EigenSubdomainSolverType parseBuiltinSolverType(const std::string& type) const;
    void configureFromInputDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    Eigen::MatrixXd buildSchurSolveMatrix(const Eigen::MatrixXd& schur) const;

    StaggeredStokesPETScLevelSolver& d_stokes_solver;
    std::vector<CustomEigenSchurSubdomainCache> d_subdomain_caches;
    std::unique_ptr<CustomEigenA00SolverStorageBase> d_a00_solver_storage;
    EigenSubdomainSolverType d_a00_solver_type = EigenSubdomainSolverType::FULL_PIV_HOUSEHOLDER_QR;
    double d_a00_solver_threshold = -1.0;
    EigenSubdomainSolverType d_schur_solver_type = EigenSubdomainSolverType::FULL_PIV_HOUSEHOLDER_QR;
    double d_schur_solver_threshold = -1.0;
};

} // namespace IBAMR

#include <ibamr/private/StaggeredStokesEigenSchurComplementShellBackend-inl.h>

#endif
