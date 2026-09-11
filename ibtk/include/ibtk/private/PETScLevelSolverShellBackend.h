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

#ifndef included_IBTK_private_PETScLevelSolverShellBackend
#define included_IBTK_private_PETScLevelSolverShellBackend

#include <ibtk/config.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <petscksp.h>

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace IBTK
{
/*! \brief Visit order for multiplicative shell composition.
 *
 * FORWARD visits 0,...,N-1; REVERSE visits N-1,...,0; SYMMETRIC visits
 * 0,...,N-1,N-2,...,0, with the turnaround visited once. Symmetric order
 * alone does not guarantee a symmetric positive definite preconditioner.
 */
enum class PETScLevelSolverShellTraversal
{
    FORWARD,
    REVERSE,
    SYMMETRIC
};

/*! \brief Shared additive and multiplicative shell composition.
 *
 * initializeSolverState() receives sorted overlapping index sets and matching
 * nonoverlapping subsets. Additive application solves independent overlapping
 * systems and writes the subsets, which must partition locally owned DOFs.
 * Multiplicative application ignores the subsets, adds each full overlapping
 * correction, and updates
 * the original residual before the next subdomain. Parallel ranks advance through
 * the same number of stages, summing all corrections in each stage. Traversal
 * defaults to FORWARD; additive composition requires FORWARD. See
 * PETScLevelSolverShellTraversal for the multiplicative visit orders.
 *
 * Inputs are borrowed during initialization; implementations retain objects they
 * need afterwards. Reinitialize after changing the operator or subdomains.
 * apply() overwrites its distinct output without modifying its input and requires
 * initialized state. Operations are collective over the level communicator.
 * deallocateSolverState() may be called repeatedly, including during destruction.
 */
class PETScLevelSolverShellBackend
{
public:
    /*! \brief Construct uninitialized composition state. */
    PETScLevelSolverShellBackend() = default;
    /*! \brief Disable copying of owned PETSc objects. */
    PETScLevelSolverShellBackend(const PETScLevelSolverShellBackend&) = delete;
    /*! \brief Disable assignment of owned PETSc objects. */
    PETScLevelSolverShellBackend& operator=(const PETScLevelSolverShellBackend&) = delete;
    /*! \brief Release composition state. */
    virtual ~PETScLevelSolverShellBackend();
    /*! \brief Initialize local solves and communication for the supplied layout. */
    virtual void
    initializeSolverState(Mat mat,
                          Vec x,
                          Vec b,
                          const std::vector<IS>& overlap,
                          const std::vector<IS>& nonoverlap,
                          const std::string& options_prefix,
                          bool use_multiplicative = false,
                          PETScLevelSolverShellTraversal traversal = PETScLevelSolverShellTraversal::FORWARD) = 0;
    /*! \brief Release the initialized state. */
    virtual void deallocateSolverState() = 0;
    /*! \brief Apply the shared correction sequence. */
    void apply(Vec x, Vec y);

protected:
    /*! \brief Retain the operator and allocate residual scratch before local setup.
     *
     * Finalize composition after local correction views exist. Release composition
     * before destroying those views. Row-access application reads current values
     * at cached offsets and rejects changed columns at those offsets; this is not
     * a detector for arbitrary sparsity changes. Reassembly requires reinitialization.
     */
    void initializeComposition(Mat mat,
                               Vec x,
                               Vec b,
                               bool use_multiplicative,
                               PETScLevelSolverShellTraversal traversal = PETScLevelSolverShellTraversal::FORWARD);
    /*! \brief Construct residual-update metadata after local setup. */
    void finalizeComposition();
    /*! \brief Release composition scratch and references. */
    void deallocateComposition();
    /*! \brief Return the collective number of stages, including empty local stages. */
    virtual std::size_t getNumberOfSubdomains() const = 0;
    /*! \brief Begin gathering this stage's RHS from the supplied residual or RHS. */
    virtual void beginSubdomainRhs(std::size_t i, Vec source) = 0;
    /*! \brief Complete the matching RHS gather. */
    virtual void endSubdomainRhs(std::size_t i, Vec source) = 0;
    /*! \brief Solve the gathered local system. */
    virtual void solveSubdomain(std::size_t i) = 0;
    /*! \brief Accumulate the stage correction using the initialized composition. */
    virtual void accumulateSubdomainCorrection(std::size_t i, Vec y) = 0;
    /*! \brief Return sorted, unique global correction indices, empty for an empty stage.
     *
     * Keep indices unchanged until teardown. They describe every value accumulated
     * by a multiplicative stage, including nonlocal overlapping contributions.
     */
    virtual const std::vector<PetscInt>& getSubdomainCorrectionDofs(std::size_t i) const = 0;
    /*! \brief Copy correction values in getSubdomainCorrectionDofs() order. */
    virtual void copySubdomainCorrection(std::size_t i, PetscScalar* values) = 0;

private:
    struct UpdateEntry
    {
        PetscInt matrix_entry;
        PetscInt correction_entry;
    };
    struct UpdateRow
    {
        PetscInt global_row;
        std::size_t entry_begin;
        std::size_t entry_end;
    };
    struct SubdomainUpdate
    {
        std::size_t row_begin;
        std::size_t row_end;
    };
    /*! \brief Cache flat affected-row offsets for serial row-access matrices. */
    void initializeAffectedRows();
    /*! \brief Subtract the actual accumulated stage action from the residual. */
    void updateResidual(std::size_t i);

    PETScLevelSolverShellTraversal d_traversal = PETScLevelSolverShellTraversal::FORWARD;
    Mat d_mat = nullptr;
    Vec d_residual = nullptr, d_correction = nullptr, d_action = nullptr;
    bool d_multiplicative = false, d_initialized = false, d_use_rows = false;
    std::vector<SubdomainUpdate> d_updates;
    std::vector<UpdateRow> d_rows;
    std::vector<UpdateEntry> d_entries;
    std::vector<PetscScalar> d_values;
};

/*! \brief Factories selected by the backend suffix of shell_pc_type.
 *
 * The built-in "petsc" factory is also selected by bare composition names. Register other
 * factories before initializing a level solver. Keys are case-sensitive and
 * must be nonempty; registering an existing key replaces its factory.
 */
class PETScLevelSolverShellBackendManager
{
public:
    using Factory =
        std::unique_ptr<PETScLevelSolverShellBackend> (*)(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);
    /*! \brief Return the process-local registry. */
    static PETScLevelSolverShellBackendManager& get_manager();
    /*! \brief Register or replace a factory with a nonnull function. */
    void registerFactory(const std::string& key, Factory factory);
    /*! \brief Construct the named backend; unknown keys are fatal errors. */
    std::unique_ptr<PETScLevelSolverShellBackend>
    allocateBackend(const std::string& key, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db) const;

private:
    /*! \brief Register the built-in backends. */
    PETScLevelSolverShellBackendManager();
    std::map<std::string, Factory> d_factories;
};
} // namespace IBTK
#endif
