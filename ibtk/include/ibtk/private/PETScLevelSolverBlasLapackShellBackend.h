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

#include <ibtk/config.h>

#include <ibtk/private/PETScLevelSolverShellBackend.h>

#include <petscblaslapack.h>

namespace IBTK
{
/*! \brief Serial real-scalar BLAS/LAPACK solves for PETScLevelSolverShellBackend.
 *
 * Local solvers use the settings documented in PETScLevelSolver. Each subdomain
 * retains a factorization or inverse/pseudoinverse and reusable solve vectors.
 */
class PETScLevelSolverBlasLapackShellBackend : public PETScLevelSolverShellBackend
{
public:
    /*! \brief Read local solver settings, using defaults for a null database. */
    explicit PETScLevelSolverBlasLapackShellBackend(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);
    /*! \copydoc PETScLevelSolverShellBackend::initializeSolverState */
    void initializeSolverState(Mat mat,
                               Vec x,
                               Vec b,
                               const std::vector<IS>& overlap,
                               const std::vector<IS>& nonoverlap,
                               const std::string& options_prefix) override;
    /*! \copydoc PETScLevelSolverShellBackend::deallocateSolverState */
    void deallocateSolverState() override;
    /*! \copydoc PETScLevelSolverShellBackend::apply */
    void apply(Vec x, Vec y) override;

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
        std::vector<PetscInt> overlap_dofs;
        std::vector<PetscBLASInt> update_local_positions;
        PetscBLASInt local_size = 0;
        std::vector<PetscScalar> solve_data;
        std::vector<PetscBLASInt> pivots;
        std::vector<PetscScalar> rhs_workspace;
        std::vector<PetscScalar> solution_workspace;
    };

    /*! \brief Factor the subdomain matrix or construct its solve matrix. */
    void initializeSubdomainSolver(SubdomainData& subdomain_data, std::size_t subdomain_num);
    /*! \brief Construct an inverse from the full-rank QR factorization. */
    void initializeQRSolver(SubdomainData& subdomain_data, std::size_t subdomain_num);
    /*! \brief Construct the SVD pseudoinverse using the configured rank tolerance. */
    void initializeSVDSolver(SubdomainData& subdomain_data, std::size_t subdomain_num);
    /*! \brief Require symmetry for the symmetric-indefinite factorization. */
    void verifySymmetricSubdomainMatrix(const SubdomainData& subdomain_data, std::size_t subdomain_num) const;
    /*! \brief Overwrite the local RHS with its correction. */
    void solveSubdomainSystem(SubdomainData& subdomain_data, std::size_t subdomain_num) const;

    std::vector<SubdomainData> d_subdomains;
    std::string d_options_prefix;
    SubdomainSolverType d_subdomain_solver_type = SubdomainSolverType::SVD;
    PetscReal d_subdomain_solver_rcond = -1.0;
};
} // namespace IBTK
#endif
