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

#include <tbox/Pointer.h>

#include <petscis.h>
#include <petscmat.h>
#include <petscvec.h>

#include <cstddef>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace SAMRAI
{
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

namespace IBTK
{
class PETScLevelSolver;

/*!
 * \brief Hierarchy-dependent PETSc and subdomain data supplied to shell
 * smoother backends.
 *
 * PETScLevelSolver populates this structure after assembling the level
 * operator and cached subdomain description. Shell backends then use this
 * state to initialize any backend-specific data structures and to apply
 * additive or multiplicative sweeps.
 */
struct PETScLevelSolverShellBackendState
{
    std::string object_name;
    std::string options_prefix;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db;
    bool use_multiplicative = false;
    bool use_restrict_partition = false;
    Mat petsc_mat = nullptr;
    Vec petsc_x = nullptr;
    Vec petsc_b = nullptr;
    const std::vector<std::vector<int>>* subdomain_dofs = nullptr;
    const std::vector<std::vector<int>>* nonoverlap_subdomain_dofs = nullptr;
    const std::function<void(int, Mat, Vec, Vec, Vec)>* subdomain_solve_observer = nullptr;
    const std::function<bool(int)>* subdomain_solve_observer_predicate = nullptr;
    std::function<void(Vec)> postprocess_result;
};

/*!
 * \brief Abstract interface implemented by shell preconditioner backends used
 * by PETScLevelSolver.
 *
 * A backend owns the algebra and representation required to solve one local
 * subdomain problem. This base class owns correction-composition ordering so
 * that concrete backends do not duplicate additive or multiplicative sweep
 * semantics.
 */
class PETScLevelSolverShellBackend
{
public:
    virtual ~PETScLevelSolverShellBackend() = default;

    /*!
     * \brief Return the backend name appended after the standardized PETSc
     * shell composition prefix.
     */
    virtual const std::string& getName() const = 0;

    /*!
     * \brief Initialize backend-specific data for one solver state.
     */
    virtual void initializeSolverState(const PETScLevelSolverShellBackendState& solver_state) = 0;

    /*!
     * \brief Release hierarchy-dependent backend data.
     */
    virtual void deallocateSolverState() = 0;

    /*!
     * \brief Apply the configured shell smoother action.
     *
     * The common composer controls subdomain order, observer placement,
     * correction accumulation, residual-update placement, and final
     * postprocessing. Concrete backends implement only the representation-
     * specific one-subdomain operations below.
     */
    void apply(Vec x, Vec y);

protected:
    /*!
     * \brief Initialize and finalize the representation shared by all
     * correction-composition backends.
     *
     * Concrete backends call initializeCorrectionCompositionState() before
     * constructing their local solver data and
     * finalizeCorrectionCompositionState() after the correction DOF views
     * returned below are valid. The views and the live level matrix must
     * remain valid until deallocateCorrectionCompositionState() is called.
     * As with the backend's local solve data, operator reassembly requires
     * solver-state reinitialization; borrowing live values avoids an update-
     * matrix copy but does not relax that lifecycle invariant.
     */
    void initializeCorrectionCompositionState(const PETScLevelSolverShellBackendState& solver_state);
    void finalizeCorrectionCompositionState();
    void deallocateCorrectionCompositionState();

    /*!
     * \brief Return the original right-hand side for additive composition or
     * the current original residual for multiplicative composition.
     */
    Vec getSubdomainResidualSource(Vec original_rhs) const;

    virtual std::size_t getNumberOfSubdomains() const = 0;
    virtual void initializeSubdomainSweep(Vec x, Vec y) = 0;
    virtual void beginSubdomainRhs(std::size_t subdomain_num, Vec x, Vec y) = 0;
    virtual void endSubdomainRhs(std::size_t subdomain_num, Vec x, Vec y) = 0;
    virtual void solveSubdomain(std::size_t subdomain_num) = 0;
    virtual void observeSubdomainSolve(std::size_t subdomain_num, Vec x, Vec y);
    virtual void accumulateSubdomainCorrection(std::size_t subdomain_num, Vec y) = 0;
    virtual const std::vector<int>& getSubdomainCorrectionDofs(std::size_t subdomain_num) const = 0;
    virtual void copySubdomainCorrection(std::size_t subdomain_num, PetscScalar* correction_values) = 0;
    virtual void finalizeSubdomainSweep(Vec x, Vec y);

    bool shouldObserveSubdomainSolve(std::size_t subdomain_num) const;

    PETScLevelSolverShellBackendState d_solver_state;

private:
    struct ResidualUpdateEntry
    {
        PetscInt matrix_entry = -1;
        PetscInt correction_entry = -1;
    };

    struct ResidualUpdateRow
    {
        PetscInt global_row = -1;
        std::vector<ResidualUpdateEntry> entries;
    };

    struct SubdomainResidualUpdate
    {
        const std::vector<int>* correction_dofs = nullptr;
        std::vector<ResidualUpdateRow> rows;
    };

    void initializeAffectedRowResidualUpdates();
    void initializeFullResidualUpdateFallback();
    void updateSubdomainResidual(std::size_t subdomain_num);
    void updateAffectedResidualRows(std::size_t subdomain_num);
    void updateFullResidualFallback(std::size_t subdomain_num);

    Vec d_sweep_residual = nullptr;
    Vec d_full_correction = nullptr;
    Vec d_full_residual_update = nullptr;
    PetscInt d_local_row_begin = 0;
    PetscInt d_local_row_end = 0;
    std::vector<SubdomainResidualUpdate> d_subdomain_residual_updates;
    std::vector<PetscScalar> d_correction_values;
    bool d_use_affected_row_updates = false;
    bool d_correction_composition_initialized = false;
};

/*!
 * \brief Registry of PETScLevelSolver shell backend factory functions.
 *
 * Concrete backend implementations register factory functions with this
 * manager so PETScLevelSolver can construct the configured backend without
 * depending on concrete backend class definitions in its public header.
 */
class PETScLevelSolverShellBackendManager
{
public:
    using ShellBackendMaker = std::unique_ptr<PETScLevelSolverShellBackend> (*)(PETScLevelSolver& solver);

    static PETScLevelSolverShellBackendManager* getManager();

    static void freeManager();

    std::unique_ptr<PETScLevelSolverShellBackend> allocateShellBackend(const std::string& type_key,
                                                                       PETScLevelSolver& solver) const;

    bool hasShellBackendFactoryFunction(const std::string& type_key) const;

    void registerShellBackendFactoryFunction(const std::string& type_key, ShellBackendMaker backend_maker);

    std::vector<std::string> getRegisteredShellBackendTypes() const;

protected:
    PETScLevelSolverShellBackendManager();
    ~PETScLevelSolverShellBackendManager() = default;

private:
    PETScLevelSolverShellBackendManager(const PETScLevelSolverShellBackendManager& from) = delete;
    PETScLevelSolverShellBackendManager& operator=(const PETScLevelSolverShellBackendManager& that) = delete;

    static PETScLevelSolverShellBackendManager* s_shell_backend_manager_instance;
    static bool s_registered_callback;
    static unsigned char s_shutdown_priority;

    std::map<std::string, ShellBackendMaker> d_shell_backend_maker_map;
};

void destroy_petsc_index_sets(std::vector<IS>& index_sets);

void build_petsc_subdomain_index_sets(std::vector<IS>& subdomain_is,
                                      std::vector<IS>& nonoverlap_subdomain_is,
                                      const std::vector<std::vector<int>>& subdomain_dofs,
                                      const std::vector<std::vector<int>>& nonoverlap_subdomain_dofs);
} // namespace IBTK

#endif
