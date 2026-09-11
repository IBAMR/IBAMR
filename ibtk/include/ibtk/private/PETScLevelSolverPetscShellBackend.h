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

#ifndef included_IBTK_private_PETScLevelSolverPetscShellBackend
#define included_IBTK_private_PETScLevelSolverPetscShellBackend

#include <ibtk/config.h>

#include <ibtk/private/PETScLevelSolverShellBackend.h>

namespace IBTK
{
/*! \brief PETSc submatrix solves for PETScLevelSolverShellBackend.
 *
 * Local solvers default to preonly/LU, with the level options prefix followed
 * by "_sub". PETSc options may override the local solver configuration.
 */
class PETScLevelSolverPetscShellBackend : public PETScLevelSolverShellBackend
{
public:
    /*! \brief Construct an uninitialized backend. */
    PETScLevelSolverPetscShellBackend() = default;
    /*! \brief Disable copying of owned PETSc objects. */
    PETScLevelSolverPetscShellBackend(const PETScLevelSolverPetscShellBackend&) = delete;
    /*! \brief Disable assignment of owned PETSc objects. */
    PETScLevelSolverPetscShellBackend& operator=(const PETScLevelSolverPetscShellBackend&) = delete;

    /*! \brief Release local PETSc objects. */
    ~PETScLevelSolverPetscShellBackend() override;
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
    std::vector<KSP> d_sub_ksp;
    std::vector<Vec> d_sub_x, d_sub_y;
    std::vector<VecScatter> d_restriction, d_prolongation;
};
} // namespace IBTK
#endif
