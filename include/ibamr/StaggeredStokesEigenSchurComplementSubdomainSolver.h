// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBAMR_StaggeredStokesEigenSchurComplementSubdomainSolver
#define included_IBAMR_StaggeredStokesEigenSchurComplementSubdomainSolver

#include <ibamr/config.h>

#include <ibtk/EigenDenseSolver.h>
#include <ibtk/PETScLevelSolverSubdomainSolver.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <Eigen/Core>

#include <functional>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace IBAMR
{
/*! \brief Local Stokes block solves using the Schur complement S = A11 - A10 A00^{-1} A01.
 *
 * The velocity and pressure DOFs of the level are obtained from a callback each time the solver state is
 * initialized, as a global vector that is 0 at velocity DOFs and 1 at pressure DOFs, and gathered on every
 * DOF of every subdomain, including those owned by other ranks. Every such DOF must have one of these
 * values. Field numbering need not be contiguous. The solver keeps only the positions of the fields in
 * each subdomain and solve matrices, so that applying the solver does not allocate. Gathering the
 * indicators communicates, so initializeSolverState() must be called on every rank, including ranks
 * without subdomains, as PETScLevelSolver does; solve() does not communicate. The settings are
 * documented in StaggeredStokesPETScLevelSolver.
 */
class StaggeredStokesEigenSchurComplementSubdomainSolver
{
public:
    /*! \brief Create the global vector of field indicators, which the solver destroys. */
    using FieldProvider = std::function<Vec()>;

    /*! \brief Read the A00 and Schur solver settings. */
    StaggeredStokesEigenSchurComplementSubdomainSolver(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                                       FieldProvider fields);
    /*! \copydoc IBTK::PETScLevelSolverSubdomainSolver::initializeSolverState */
    void initializeSolverState(const std::vector<Mat>& matrices,
                               const std::vector<IS>& subdomains,
                               const std::string& options_prefix);
    /*! \copydoc IBTK::PETScLevelSolverSubdomainSolver::deallocateSolverState */
    void deallocateSolverState();
    /*! \copydoc IBTK::PETScLevelSolverSubdomainSolver::solve */
    void solve(std::size_t first, std::size_t last, Vec b, Vec x);

private:
    struct Block
    {
        // Positions of the velocity and pressure DOFs in the subdomain.
        std::vector<Eigen::Index> velocity, pressure;
        Eigen::MatrixXd a10, a00_solve, a00_inv_a01, schur_solve;
        Eigen::VectorXd velocity_rhs, pressure_rhs, velocity_solution, pressure_solution;
    };

    FieldProvider d_fields;
    std::string d_a00_type = "FULL_PIV_HOUSEHOLDER_QR", d_schur_type = "FULL_PIV_HOUSEHOLDER_QR";
    double d_a00_threshold = -1.0, d_schur_threshold = -1.0;
    std::vector<PetscInt> d_offsets;
    std::vector<Block> d_blocks;
};
} // namespace IBAMR
#endif
