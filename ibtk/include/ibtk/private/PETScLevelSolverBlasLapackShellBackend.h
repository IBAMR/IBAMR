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
    enum class SubdomainSolverType
    {
        SVD,
        LU,
        CHOLESKY,
        SYMMETRIC_INDEFINITE,
        QR
    };

    struct SolverDataBase
    {
        virtual ~SolverDataBase() = default;
    };

    struct LUSolverData : public SolverDataBase
    {
        std::vector<PetscScalar> factor;
        std::vector<PetscBLASInt> pivot;
    };

    struct CholeskySolverData : public SolverDataBase
    {
        std::vector<PetscScalar> factor;
    };

    struct SymmetricIndefiniteSolverData : public SolverDataBase
    {
        std::vector<PetscScalar> factor;
        std::vector<PetscBLASInt> pivot;
        std::vector<PetscScalar> work;
    };

    struct DenseSolveMatrixSolverData : public SolverDataBase
    {
        std::vector<PetscScalar> solve_matrix;
    };

    struct SVDSolverData : public DenseSolveMatrixSolverData
    {
        std::vector<PetscReal> singular_values;
        PetscBLASInt effective_rank = 0;
    };

    struct SubdomainData
    {
        std::size_t subdomain_num = 0;
        std::vector<int> overlap_dofs;
        std::vector<int> update_dofs;
        std::vector<int> update_local_positions;
        PetscBLASInt local_size = 0;
        PetscBLASInt local_lda = 1;
        std::unique_ptr<SolverDataBase> solver_data;

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
    template <class SolverDataType>
    SolverDataType& getSolverData(SubdomainData& subdomain_data) const;

    template <class SolverDataType>
    const SolverDataType& getSolverData(const SubdomainData& subdomain_data) const;

    void initializeSubdomainSolveData(SubdomainData& subdomain_data,
                                      const std::vector<PetscScalar>& local_operator) const;
    void initializeSubdomainData(SubdomainData& subdomain_data,
                                 std::vector<PetscScalar>& local_operator,
                                 std::size_t subdomain_num,
                                 const std::vector<int>& overlap_dofs,
                                 const std::vector<int>& nonoverlap_dofs,
                                 bool use_restrict_partition,
                                 bool use_multiplicative,
                                 Mat level_mat) const;
    void buildQRSolveMatrix(SubdomainData& subdomain_data, const std::vector<PetscScalar>& local_operator) const;
    void buildSVDSolveMatrix(SubdomainData& subdomain_data, const std::vector<PetscScalar>& local_operator) const;
    void verifySymmetricSubdomainMatrix(const SubdomainData& subdomain_data,
                                        const std::vector<PetscScalar>& local_operator) const;
    void solveSubdomainSystem(SubdomainData& subdomain_data) const;
    void updateResidual(SubdomainData& subdomain_data, std::vector<PetscScalar>& residual) const;
    const char* getSolverTypeName() const;

    PETScLevelSolverBackendContext& d_context;
    std::vector<SubdomainData> d_subdomains;
    PetscBLASInt d_n_dofs = 0;
    std::string d_type_key = "blas-lapack";
    SubdomainSolverType d_subdomain_solver_type = SubdomainSolverType::SVD;
    PetscReal d_subdomain_solver_rcond = -1.0;
};
} // namespace IBTK

#endif
