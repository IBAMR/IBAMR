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

#include <ibtk/PETScLevelSolver.h>

#include <petscblaslapack.h>

#include <memory>
#include <string>
#include <vector>

namespace IBTK
{
class PETScLevelSolverBlasLapackShellBackend : public PETScLevelSolverShellBackend
{
public:
    explicit PETScLevelSolverBlasLapackShellBackend(PETScLevelSolver& solver);

    const std::string& getTypeKey() const override;
    void configure(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db) override;
    const char* getPCNameSuffixAdditive() const override;
    const char* getPCNameSuffixMultiplicative() const override;
    void initialize() override;
    void deallocate() override;
    void applyAdditive(Vec x, Vec y) override;
    void applyMultiplicative(Vec x, Vec y) override;

private:
    struct SubdomainData
    {
        std::vector<int> overlap_dofs;
        std::vector<int> update_dofs;
        std::vector<int> update_local_positions;
        std::vector<PetscScalar> local_operator_lu;
        std::vector<PetscBLASInt> pivot;
        PetscBLASInt local_size = 0;
        PetscBLASInt local_lda = 1;

        std::vector<int> active_residual_update_rows;
        std::vector<PetscScalar> active_residual_update_mat;
        PetscBLASInt active_residual_update_num_rows = 0;
        PetscBLASInt active_residual_update_num_cols = 0;
        PetscBLASInt active_residual_update_lda = 1;

        std::vector<PetscScalar> rhs_workspace;
        std::vector<PetscScalar> delta_workspace;
        std::vector<PetscScalar> residual_input_workspace;
        std::vector<PetscScalar> residual_delta_workspace;
    };

    void parseSolverType(const std::string& type_name);
    void initializeSubdomainData(SubdomainData& subdomain_data,
                                 const std::vector<int>& overlap_dofs,
                                 const std::vector<int>& nonoverlap_dofs,
                                 bool use_restrict_partition,
                                 bool use_multiplicative,
                                 Mat level_mat) const;
    void factorizeSubdomainMatrix(SubdomainData& subdomain_data) const;
    void solveSubdomainSystem(SubdomainData& subdomain_data) const;
    void updateResidual(SubdomainData& subdomain_data, std::vector<PetscScalar>& residual) const;

    PETScLevelSolverBackendContext& d_context;
    std::vector<SubdomainData> d_subdomains;
    PetscBLASInt d_n_dofs = 0;
    std::string d_type_key = "blas-lapack";
    std::string d_subdomain_solver_type = "lu";
};
} // namespace IBTK

#endif
