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

#ifndef included_IBTK_private_PETScLevelSolverBlasLapackShellBackend
#define included_IBTK_private_PETScLevelSolverBlasLapackShellBackend

#include <ibtk/private/PETScLevelSolverShellBackend.h>

#include <petscblaslapack.h>

#include <memory>
#include <string>
#include <vector>

namespace IBTK
{
/*!
 * \brief Serial real-scalar BLAS/LAPACK shell smoother backend.
 *
 * The backend borrows the full PETSc operator in global numbering and stores
 * one dense, column-major factor or solve representation per subdomain. The
 * unfactored local matrices are not retained. Multiplicative residual updates
 * use the borrowed operator and reusable global vectors, so the application
 * path does not copy the full operator or build update matrices.
 *
 * The configurable backend name uses SVD by default and also supports LU,
 * symmetric-indefinite, and QR local solvers. The legacy fixed-LU name shares
 * this implementation. Cholesky is rejected because these Stokes subdomain
 * matrices are indefinite.
 */
class PETScLevelSolverBlasLapackShellBackend : public PETScLevelSolverShellBackend
{
public:
    PETScLevelSolverBlasLapackShellBackend(std::string backend_name, bool configurable_solver_type);

    const std::string& getName() const override;
    void initializeSolverState(const PETScLevelSolverShellBackendState& solver_state) override;
    void deallocateSolverState() override;

private:
    enum class SubdomainSolverType
    {
        SVD,
        LU,
        SYMMETRIC_INDEFINITE,
        QR
    };

    struct SubdomainData
    {
        const std::vector<int>* overlap_dofs = nullptr;
        const std::vector<int>* update_dofs = nullptr;
        std::vector<PetscBLASInt> update_local_positions;
        PetscBLASInt local_size = 0;
        std::vector<PetscScalar> solve_data;
        std::vector<PetscBLASInt> pivots;
        std::vector<PetscScalar> rhs_workspace;
        std::vector<PetscScalar> solution_workspace;
    };

    struct Data
    {
        Vec residual = nullptr;
        Vec patch_correction = nullptr;
        Vec residual_update = nullptr;
        PetscInt n_dofs = 0;
        std::vector<SubdomainData> subdomains;
    };

    void configureFromInputDatabase();
    void initializeSubdomainSolver(SubdomainData& subdomain_data, std::size_t subdomain_num);
    void initializeQRSolver(SubdomainData& subdomain_data, std::size_t subdomain_num);
    void initializeSVDSolver(SubdomainData& subdomain_data, std::size_t subdomain_num);
    void verifySymmetricSubdomainMatrix(const SubdomainData& subdomain_data, std::size_t subdomain_num) const;
    void solveSubdomainSystem(SubdomainData& subdomain_data, std::size_t subdomain_num) const;
    const char* getSolverTypeName() const;
    void observeSubdomain(std::size_t subdomain_num, SubdomainData& subdomain_data, Vec current_global_source) const;

    std::size_t getNumberOfSubdomains() const override;
    void initializeSubdomainSweep(Vec x, Vec y) override;
    void beginSubdomainRhs(std::size_t subdomain_num, Vec x, Vec y) override;
    void endSubdomainRhs(std::size_t subdomain_num, Vec x, Vec y) override;
    void solveSubdomain(std::size_t subdomain_num) override;
    void observeSubdomainSolve(std::size_t subdomain_num, Vec x, Vec y) override;
    void accumulateSubdomainCorrection(std::size_t subdomain_num, Vec y) override;
    void updateSubdomainResidual(std::size_t subdomain_num, Vec x, Vec y) override;

    std::unique_ptr<Data> d_data;
    std::string d_backend_name;
    bool d_configurable_solver_type = false;
    SubdomainSolverType d_subdomain_solver_type = SubdomainSolverType::LU;
    PetscReal d_subdomain_solver_rcond = -1.0;
};
} // namespace IBTK

#endif
