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

#ifndef included_IBAMR_private_StaggeredStokesEigenSchurComplementShellBackend
#define included_IBAMR_private_StaggeredStokesEigenSchurComplementShellBackend

#include <ibamr/config.h>

#include <ibtk/private/EigenLocalSolver.h>
#include <ibtk/private/PETScLevelSolverEigenShellBackendBase.h>

#include <set>

namespace IBAMR
{
/*! \brief Serial local Stokes block solves using S = A11 - A10 A00^{-1} A01.
 *
 * The field-aware initializer requires disjoint velocity and pressure sets
 * covering every overlap DOF. It consumes these sets synchronously and keeps
 * only local positions, factorizations and solve workspaces. Field numbering
 * need not be contiguous. Settings are documented in StaggeredStokesPETScLevelSolver.
 */
class StaggeredStokesEigenSchurComplementShellBackend : public IBTK::PETScLevelSolverEigenShellBackendBase
{
public:
    /*! \brief Read the A00 and Schur solver settings. */
    explicit StaggeredStokesEigenSchurComplementShellBackend(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);
    /*! \brief Release local and composition state. */
    ~StaggeredStokesEigenSchurComplementShellBackend() override;
    /*! \brief Reject initialization without Stokes field IDs. */
    void initializeSolverState(
        Mat mat,
        Vec x,
        Vec b,
        const std::vector<IS>& overlap,
        const std::vector<IS>& nonoverlap,
        const std::string& options_prefix,
        bool use_multiplicative = false,
        IBTK::PETScLevelSolverShellTraversal traversal = IBTK::PETScLevelSolverShellTraversal::FORWARD) override;
    /*! \brief Initialize local blocks from the supplied velocity and pressure IDs. */
    void initializeSolverState(Mat mat,
                               Vec x,
                               Vec b,
                               const std::vector<IS>& overlap,
                               const std::vector<IS>& nonoverlap,
                               const std::set<int>& velocity_dofs,
                               const std::set<int>& pressure_dofs,
                               const std::string& options_prefix,
                               bool use_multiplicative,
                               IBTK::PETScLevelSolverShellTraversal traversal);
    /*! \copydoc IBTK::PETScLevelSolverShellBackend::deallocateSolverState */
    void deallocateSolverState() override;

protected:
    /*! \copydoc IBTK::PETScLevelSolverShellBackend::solveSubdomain */
    void solveSubdomain(std::size_t i) override;

private:
    struct SchurSubdomain
    {
        std::vector<Eigen::Index> velocity, pressure;
        std::unique_ptr<IBTK::EigenLocalSolver> a00_solver;
        Eigen::MatrixXd a10, a00_inv_a01, schur_solve;
        Eigen::VectorXd velocity_rhs, pressure_rhs, velocity_solution, pressure_solution;
    };
    std::string d_a00_type = "FULL_PIV_HOUSEHOLDER_QR", d_schur_type = "FULL_PIV_HOUSEHOLDER_QR";
    double d_a00_threshold = -1.0, d_schur_threshold = -1.0;
    std::vector<SchurSubdomain> d_blocks;
};
} // namespace IBAMR
#endif
