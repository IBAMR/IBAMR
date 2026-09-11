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

#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/private/PETScLevelSolverEigenShellBackendBase.h>

#include <tbox/Utilities.h>

#include <algorithm>

namespace IBTK
{
void
PETScLevelSolverEigenShellBackendBase::initializeSubdomains(Mat mat,
                                                            Vec x,
                                                            Vec b,
                                                            const std::vector<IS>& overlap,
                                                            const std::vector<IS>& nonoverlap,
                                                            const bool multiplicative,
                                                            const PETScLevelSolverShellTraversal traversal)
{
    if (IBTK_MPI::getNodes() != 1)
    {
        TBOX_ERROR("Eigen shell backends require one MPI rank.\n");
    }
    PetscInt rows = 0, cols = 0, local_rows = 0;
    int ierr = MatGetSize(mat, &rows, &cols);
    IBTK_CHKERRQ(ierr);
    ierr = MatGetLocalSize(mat, &local_rows, nullptr);
    IBTK_CHKERRQ(ierr);
    if (rows != cols || rows != local_rows)
    {
        TBOX_ERROR("Eigen shell backends require a serial square operator.\n");
    }
    initializeComposition(mat, x, b, multiplicative, traversal);
    d_multiplicative = multiplicative;
    TBOX_ASSERT(multiplicative || overlap.size() == nonoverlap.size());
    d_subdomains.resize(overlap.size());
    for (std::size_t i = 0; i < overlap.size(); ++i)
    {
        PetscInt n = 0;
        const PetscInt* ids = nullptr;
        ierr = ISGetLocalSize(overlap[i], &n);
        IBTK_CHKERRQ(ierr);
        ierr = ISGetIndices(overlap[i], &ids);
        IBTK_CHKERRQ(ierr);
        if (n > 0)
        {
            d_subdomains[i].dofs.assign(ids, ids + n);
        }
        ierr = ISRestoreIndices(overlap[i], &ids);
        IBTK_CHKERRQ(ierr);
        d_subdomains[i].rhs.resize(n);
        d_subdomains[i].solution.resize(n);
        if (!multiplicative)
        {
            PetscInt m = 0;
            ierr = ISGetLocalSize(nonoverlap[i], &m);
            IBTK_CHKERRQ(ierr);
            ierr = ISGetIndices(nonoverlap[i], &ids);
            IBTK_CHKERRQ(ierr);
            for (PetscInt j = 0; j < m; ++j)
            {
                const std::vector<PetscInt>::const_iterator position =
                    std::lower_bound(d_subdomains[i].dofs.cbegin(), d_subdomains[i].dofs.cend(), ids[j]);
                TBOX_ASSERT(position != d_subdomains[i].dofs.cend() && *position == ids[j]);
                d_subdomains[i].restricted_positions.push_back(position - d_subdomains[i].dofs.cbegin());
            }
            ierr = ISRestoreIndices(nonoverlap[i], &ids);
            IBTK_CHKERRQ(ierr);
        }
    }
}
Eigen::MatrixXd
PETScLevelSolverEigenShellBackendBase::extractLocalOperator(Mat mat, const std::size_t i) const
{
    const PetscInt n = static_cast<PetscInt>(d_subdomains[i].dofs.size());
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(n, n);
    for (PetscInt j = 0; j < n; ++j)
    {
        PetscInt count = 0;
        const PetscInt* columns = nullptr;
        const PetscScalar* values = nullptr;
        int ierr = MatGetRow(mat, d_subdomains[i].dofs[j], &count, &columns, &values);
        IBTK_CHKERRQ(ierr);
        for (PetscInt k = 0; k < count; ++k)
        {
            const std::vector<PetscInt>::const_iterator position =
                std::lower_bound(d_subdomains[i].dofs.cbegin(), d_subdomains[i].dofs.cend(), columns[k]);
            if (position != d_subdomains[i].dofs.cend() && *position == columns[k])
            {
                matrix(j, position - d_subdomains[i].dofs.cbegin()) = PetscRealPart(values[k]);
            }
        }
        ierr = MatRestoreRow(mat, d_subdomains[i].dofs[j], &count, &columns, &values);
        IBTK_CHKERRQ(ierr);
    }
    return matrix;
}
void
PETScLevelSolverEigenShellBackendBase::clearSubdomains()
{
    deallocateComposition();
    d_subdomains.clear();
}
std::size_t
PETScLevelSolverEigenShellBackendBase::getNumberOfSubdomains() const
{
    return d_subdomains.size();
}
void
PETScLevelSolverEigenShellBackendBase::beginSubdomainRhs(const std::size_t i, Vec source)
{
    const PetscScalar* values = nullptr;
    int ierr = VecGetArrayRead(source, &values);
    IBTK_CHKERRQ(ierr);
    for (std::size_t j = 0; j < d_subdomains[i].dofs.size(); ++j)
    {
        d_subdomains[i].rhs[j] = PetscRealPart(values[d_subdomains[i].dofs[j]]);
    }
    ierr = VecRestoreArrayRead(source, &values);
    IBTK_CHKERRQ(ierr);
}
void
PETScLevelSolverEigenShellBackendBase::endSubdomainRhs(std::size_t /*i*/, Vec /*source*/)
{
}
void
PETScLevelSolverEigenShellBackendBase::accumulateSubdomainCorrection(const std::size_t i, Vec y)
{
    PetscScalar* values = nullptr;
    int ierr = VecGetArray(y, &values);
    IBTK_CHKERRQ(ierr);
    if (d_multiplicative)
    {
        for (std::size_t j = 0; j < d_subdomains[i].dofs.size(); ++j)
        {
            values[d_subdomains[i].dofs[j]] += d_subdomains[i].solution[j];
        }
    }
    else
    {
        for (const Eigen::Index j : d_subdomains[i].restricted_positions)
        {
            values[d_subdomains[i].dofs[j]] += d_subdomains[i].solution[j];
        }
    }
    ierr = VecRestoreArray(y, &values);
    IBTK_CHKERRQ(ierr);
}
const std::vector<PetscInt>&
PETScLevelSolverEigenShellBackendBase::getSubdomainCorrectionDofs(const std::size_t i) const
{
    return d_subdomains[i].dofs;
}
void
PETScLevelSolverEigenShellBackendBase::copySubdomainCorrection(const std::size_t i, PetscScalar* values)
{
    for (Eigen::Index j = 0; j < d_subdomains[i].solution.size(); ++j)
    {
        values[j] = d_subdomains[i].solution[j];
    }
}
} // namespace IBTK
