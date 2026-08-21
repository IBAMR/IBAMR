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

#ifndef included_IBTK_private_PETScLevelSolverPetscShellBackend
#define included_IBTK_private_PETScLevelSolverPetscShellBackend

#include <ibtk/private/PETScLevelSolverShellBackend.h>

#include <petscksp.h>

#include <memory>
#include <vector>

namespace IBTK
{
class PETScLevelSolver;

/*!
 * \brief PETSc-native shell smoother backend for PETScLevelSolver.
 *
 * This backend mirrors the traditional PETSc ASM-style implementation by
 * constructing PETSc submatrices, scatters, and sub-KSP objects for the
 * cached subdomain description supplied by PETScLevelSolver.
 */
class PETScLevelSolverPetscShellBackend : public PETScLevelSolverShellBackend
{
public:
    static const std::string s_backend_name;

    PETScLevelSolverPetscShellBackend() = default;

    const std::string& getName() const override;
    void initializeSolverState(const PETScLevelSolverShellBackendState& solver_state) override;
    void deallocateSolverState() override;

private:
    struct Data
    {
        InsertMode prolongation_insert_mode = INSERT_VALUES;
        std::vector<IS> global_overlap_is, global_nonoverlap_is;
        std::vector<IS> local_overlap_is, local_nonoverlap_is;
        std::vector<VecScatter> restriction, prolongation;
        std::vector<KSP> sub_ksp;
        Mat* sub_mat = nullptr;
        std::vector<Vec> sub_x, sub_y;
    };

    std::size_t getNumberOfSubdomains() const override;
    void initializeSubdomainSweep(Vec x, Vec y) override;
    void beginSubdomainRhs(std::size_t subdomain_num, Vec x, Vec y) override;
    void endSubdomainRhs(std::size_t subdomain_num, Vec x, Vec y) override;
    void solveSubdomain(std::size_t subdomain_num) override;
    void accumulateSubdomainCorrection(std::size_t subdomain_num, Vec y) override;
    const std::vector<int>& getSubdomainCorrectionDofs(std::size_t subdomain_num) const override;
    void copySubdomainCorrection(std::size_t subdomain_num, PetscScalar* correction_values) override;

    void beginAccumulateCorrection(int subdomain_num, Vec sub_y, Vec y);
    void endAccumulateCorrection(int subdomain_num, Vec sub_y, Vec y);

    std::unique_ptr<Data> d_data;
};
} // namespace IBTK

#endif
