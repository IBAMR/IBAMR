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

#ifndef included_IBTK_private_PETScLevelSolverBlasLapackLUShellBackend
#define included_IBTK_private_PETScLevelSolverBlasLapackLUShellBackend

#include <ibtk/private/PETScLevelSolverShellBackend.h>

#include <petscblaslapack.h>

#include <memory>
#include <string>
#include <vector>

namespace IBTK
{
/*!
 * \brief Serial real-scalar BLAS/LAPACK LU shell smoother backend.
 *
 * The backend borrows the full PETSc operator in global numbering and stores
 * one dense, column-major LU factorization per subdomain. The unfactored local
 * matrices are overwritten during setup and are not retained. Multiplicative
 * residual updates use the borrowed operator and reusable global vectors, so
 * the application path does not copy the full operator or build update
 * matrices.
 */
class PETScLevelSolverBlasLapackLUShellBackend : public PETScLevelSolverShellBackend
{
public:
    static const std::string s_backend_name;

    PETScLevelSolverBlasLapackLUShellBackend() = default;

    const std::string& getName() const override;
    void initializeSolverState(const PETScLevelSolverShellBackendState& solver_state) override;
    void deallocateSolverState() override;
    void apply(Vec x, Vec y) override;

private:
    struct SubdomainData
    {
        const std::vector<int>* overlap_dofs = nullptr;
        const std::vector<int>* update_dofs = nullptr;
        std::vector<PetscBLASInt> update_local_positions;
        PetscBLASInt local_size = 0;
        std::vector<PetscScalar> lu_factor;
        std::vector<PetscBLASInt> pivots;
        std::vector<PetscScalar> rhs_workspace;
    };

    struct Data
    {
        Vec residual = nullptr;
        Vec patch_correction = nullptr;
        Vec residual_update = nullptr;
        PetscInt n_dofs = 0;
        std::vector<SubdomainData> subdomains;
    };

    void solveSubdomain(SubdomainData& subdomain_data, std::size_t subdomain_num) const;
    bool shouldObserveSubdomain(std::size_t subdomain_num) const;
    void observeSubdomain(std::size_t subdomain_num,
                          SubdomainData& subdomain_data,
                          std::vector<PetscScalar>& local_rhs,
                          Vec current_global_source) const;
    void applyAdditive(Vec x, Vec y);
    void applyMultiplicative(Vec x, Vec y);

    PETScLevelSolverShellBackendState d_solver_state;
    std::unique_ptr<Data> d_data;
};
} // namespace IBTK

#endif
