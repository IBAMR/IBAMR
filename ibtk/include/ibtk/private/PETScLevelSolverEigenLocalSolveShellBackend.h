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

#ifndef included_IBTK_private_PETScLevelSolverEigenLocalSolveShellBackend
#define included_IBTK_private_PETScLevelSolverEigenLocalSolveShellBackend

#include <ibtk/private/PETScLevelSolverEigenShellBackendBase.h>

namespace IBTK
{
/*!
 * \brief Common Eigen shell smoother support for overlap-based local solves.
 *
 * This class owns the standard solver-state lifecycle and adapts an Eigen
 * local solve to the backend-independent correction composer. Concrete
 * backends only provide local subdomain setup and solve operations.
 */
class PETScLevelSolverEigenLocalSolveShellBackend : public PETScLevelSolverEigenShellBackendBase
{
public:
    void initializeSolverState(const PETScLevelSolverShellBackendState& solver_state) override final;

    void deallocateSolverState() override final;

protected:
    virtual void configureFromInputDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db) = 0;

    virtual void initializeLocalSubdomainSolver(const Eigen::MatrixXd& local_operator, std::size_t subdomain_num) = 0;

    virtual Eigen::VectorXd solveLocalSubdomainSystem(const Eigen::VectorXd& rhs, std::size_t subdomain_num) const = 0;

    virtual void initializeAdditionalSolverState();

    virtual void deallocateAdditionalSolverState();

private:
    void solveSubdomain(std::size_t subdomain_num) override;
};
} // namespace IBTK

#include <ibtk/private/PETScLevelSolverEigenLocalSolveShellBackend-inl.h>

#endif
