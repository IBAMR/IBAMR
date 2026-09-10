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

#ifndef included_IBTK_private_PETScLevelSolverEigenShellBackendBase
#define included_IBTK_private_PETScLevelSolverEigenShellBackendBase

#include <ibtk/config.h>

#include <ibtk/private/PETScLevelSolverShellBackend.h>

#include <Eigen/Core>

namespace IBTK
{
/*! \brief Serial Eigen subdomain gathers and corrections for the shared composer. */
class PETScLevelSolverEigenShellBackendBase : public PETScLevelSolverShellBackend
{
protected:
    struct Subdomain
    {
        std::vector<PetscInt> dofs;
        std::vector<Eigen::Index> restricted_positions;
        Eigen::VectorXd rhs, solution;
    };
    /*! \brief Copy subdomain indices and allocate local workspaces. */
    void initializeSubdomains(Mat mat,
                              Vec x,
                              Vec b,
                              const std::vector<IS>& overlap,
                              const std::vector<IS>& nonoverlap,
                              bool multiplicative,
                              PETScLevelSolverShellTraversal traversal);
    /*! \brief Extract one temporary dense local operator during setup. */
    Eigen::MatrixXd extractLocalOperator(Mat mat, std::size_t i) const;
    /*! \brief Release common local and composition state. */
    void clearSubdomains();
    /*! \copydoc PETScLevelSolverShellBackend::getNumberOfSubdomains */
    std::size_t getNumberOfSubdomains() const override;
    /*! \copydoc PETScLevelSolverShellBackend::beginSubdomainRhs */
    void beginSubdomainRhs(std::size_t i, Vec source) override;
    /*! \copydoc PETScLevelSolverShellBackend::endSubdomainRhs */
    void endSubdomainRhs(std::size_t i, Vec source) override;
    /*! \copydoc PETScLevelSolverShellBackend::accumulateSubdomainCorrection */
    void accumulateSubdomainCorrection(std::size_t i, Vec y) override;
    /*! \copydoc PETScLevelSolverShellBackend::getSubdomainCorrectionDofs */
    const std::vector<PetscInt>& getSubdomainCorrectionDofs(std::size_t i) const override;
    /*! \copydoc PETScLevelSolverShellBackend::copySubdomainCorrection */
    void copySubdomainCorrection(std::size_t i, PetscScalar* values) override;
    std::vector<Subdomain> d_subdomains;

private:
    bool d_multiplicative = false;
};
} // namespace IBTK
#endif
