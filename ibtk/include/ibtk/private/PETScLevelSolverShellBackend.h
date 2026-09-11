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
/*! \brief Additive shell preconditioner implemented on geometric ASM subdomains.
 *
 * initializeSolverState() receives sorted overlapping index sets and matching
 * nonoverlapping subsets that partition the locally owned DOFs. The matrix and
 * vectors share the level solver's parallel layout. Inputs are borrowed during
 * initialization; implementations retain any objects needed by apply(). Rebuild
 * state after changing the matrix or subdomains. deallocateSolverState() releases
 * owned state and may be called repeatedly. apply() overwrites its distinct output
 * vector with the restricted additive correction, without modifying its input.
 * Call apply() only while solver state is initialized. Solver-state operations
 * and application are collective over the level communicator.
 */
class PETScLevelSolverShellBackend
{
public:
    /*! \brief Release the backend. */
    virtual ~PETScLevelSolverShellBackend() = default;
    /*! \brief Initialize local solves and communication for the supplied layout. */
    virtual void initializeSolverState(Mat mat,
                                       Vec x,
                                       Vec b,
                                       const std::vector<IS>& overlap,
                                       const std::vector<IS>& nonoverlap,
                                       const std::string& options_prefix) = 0;
    /*! \brief Release the initialized state. */
    virtual void deallocateSolverState() = 0;
    /*! \brief Apply the restricted additive correction. */
    virtual void apply(Vec x, Vec y) = 0;
};

/*! \brief Factories selected by the suffix of shell_pc_type = "additive-KEY".
 *
 * The built-in "petsc" factory is also selected by "additive". Register other
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
