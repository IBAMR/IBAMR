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

#ifndef included_IBTK_private_PETScLevelSolverEigenShellBackend
#define included_IBTK_private_PETScLevelSolverEigenShellBackend

#include <ibtk/config.h>

#include <ibtk/private/EigenLocalSolver.h>
#include <ibtk/private/PETScLevelSolverEigenShellBackendBase.h>

namespace IBTK
{
/*! \brief Serial factorized or precomputed Eigen local solves.
 *
 * Settings and rank policies are documented in PETScLevelSolver. Retains only
 * local factorizations or solve matrices and reusable subdomain vectors.
 */
class PETScLevelSolverEigenShellBackend : public PETScLevelSolverEigenShellBackendBase
{
public:
    /*! \brief Read the factorized or precomputed solve settings. */
    explicit PETScLevelSolverEigenShellBackend(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                               bool precompute = false);
    /*! \brief Release local and composition state. */
    ~PETScLevelSolverEigenShellBackend() override;
    /*! \copydoc PETScLevelSolverShellBackend::initializeSolverState */
    void
    initializeSolverState(Mat mat,
                          Vec x,
                          Vec b,
                          const std::vector<IS>& overlap,
                          const std::vector<IS>& nonoverlap,
                          const std::string& options_prefix,
                          bool use_multiplicative = false,
                          PETScLevelSolverShellTraversal traversal = PETScLevelSolverShellTraversal::FORWARD) override;
    /*! \copydoc PETScLevelSolverShellBackend::deallocateSolverState */
    void deallocateSolverState() override;

protected:
    /*! \copydoc PETScLevelSolverShellBackend::solveSubdomain */
    void solveSubdomain(std::size_t i) override;

private:
    const bool d_precompute;
    std::string d_solver_type;
    double d_threshold = -1.0;
    std::vector<std::unique_ptr<EigenLocalSolver>> d_solvers;
    std::vector<Eigen::MatrixXd> d_solve_matrices;
};
} // namespace IBTK
#endif
