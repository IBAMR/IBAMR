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
#include <ibtk/private/PETScLevelSolverShellBackend.h>

#include <tbox/Utilities.h>

#include <algorithm>
#include <tuple>

namespace IBTK
{
PETScLevelSolverShellBackend::~PETScLevelSolverShellBackend()
{
    deallocateComposition();
}

void
PETScLevelSolverShellBackend::apply(Vec x, Vec y)
{
    TBOX_ASSERT(d_initialized);
    int ierr = VecZeroEntries(y);
    IBTK_CHKERRQ(ierr);
    const std::size_t n_subdomains = getNumberOfSubdomains();
    if (d_multiplicative)
    {
        ierr = VecCopy(x, d_residual);
        IBTK_CHKERRQ(ierr);
        for (std::size_t i = 0; i < n_subdomains; ++i)
        {
            beginSubdomainRhs(i, d_residual);
            endSubdomainRhs(i, d_residual);
            solveSubdomain(i);
            accumulateSubdomainCorrection(i, y);
            if (i + 1 < n_subdomains)
            {
                updateResidual(i);
            }
        }
    }
    else
    {
        for (std::size_t i = 0; i < n_subdomains; ++i)
        {
            beginSubdomainRhs(i, x);
        }
        for (std::size_t i = 0; i < n_subdomains; ++i)
        {
            endSubdomainRhs(i, x);
            solveSubdomain(i);
            accumulateSubdomainCorrection(i, y);
        }
    }
}

void
PETScLevelSolverShellBackend::initializeComposition(Mat mat, Vec x, Vec b, const bool use_multiplicative)
{
    TBOX_ASSERT(!d_initialized && !d_mat);
    d_multiplicative = use_multiplicative;
    if (!d_multiplicative)
    {
        return;
    }
    d_mat = mat;
    int ierr = PetscObjectReference(reinterpret_cast<PetscObject>(d_mat));
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(b, &d_residual);
    IBTK_CHKERRQ(ierr);
    PetscBool has_rows = PETSC_FALSE;
    ierr = MatHasOperation(mat, MATOP_GET_ROW, &has_rows);
    IBTK_CHKERRQ(ierr);
    int ranks = 0;
    ierr = MPI_Comm_size(PetscObjectComm(reinterpret_cast<PetscObject>(mat)), &ranks);
    IBTK_CHKERRQ(ierr);
    d_use_rows = ranks == 1 && has_rows;
    if (!d_use_rows)
    {
        ierr = VecDuplicate(x, &d_correction);
        IBTK_CHKERRQ(ierr);
        ierr = VecDuplicate(b, &d_action);
        IBTK_CHKERRQ(ierr);
    }
}

void
PETScLevelSolverShellBackend::finalizeComposition()
{
    if (d_multiplicative)
    {
        std::size_t max_size = 0;
        for (std::size_t i = 0; i < getNumberOfSubdomains(); ++i)
        {
            max_size = std::max(max_size, getSubdomainCorrectionDofs(i).size());
        }
        d_values.resize(max_size);
        if (d_use_rows)
        {
            initializeAffectedRows();
        }
    }
    d_initialized = true;
}

void
PETScLevelSolverShellBackend::deallocateComposition()
{
    int ierr = VecDestroy(&d_residual);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&d_correction);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&d_action);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&d_mat);
    IBTK_CHKERRQ(ierr);
    d_updates = std::vector<SubdomainUpdate>();
    d_rows = std::vector<UpdateRow>();
    d_entries = std::vector<UpdateEntry>();
    d_values = std::vector<PetscScalar>();
    d_initialized = d_multiplicative = d_use_rows = false;
}

void
PETScLevelSolverShellBackend::initializeAffectedRows()
{
    struct Location
    {
        PetscInt column, row, entry;
    };
    struct CorrectionLocation
    {
        PetscInt row, entry, correction;
    };
    PetscInt n = 0, m = 0;
    int ierr = MatGetSize(d_mat, &n, &m);
    IBTK_CHKERRQ(ierr);
    if (n != m)
    {
        TBOX_ERROR("Multiplicative shell composition requires a square operator.\n");
    }
    std::vector<Location> locations;
    for (PetscInt row = 0; row < n; ++row)
    {
        PetscInt count = 0;
        const PetscInt* columns = nullptr;
        ierr = MatGetRow(d_mat, row, &count, &columns, nullptr);
        IBTK_CHKERRQ(ierr);
        for (PetscInt j = 0; j < count; ++j)
        {
            locations.push_back({ columns[j], row, j });
        }
        ierr = MatRestoreRow(d_mat, row, &count, &columns, nullptr);
        IBTK_CHKERRQ(ierr);
    }
    std::sort(locations.begin(),
              locations.end(),
              [](const Location& a, const Location& b)
              { return std::tie(a.column, a.row, a.entry) < std::tie(b.column, b.row, b.entry); });
    std::vector<CorrectionLocation> selected;
    for (std::size_t i = 0; i < getNumberOfSubdomains(); ++i)
    {
        selected.clear();
        const std::vector<PetscInt>& dofs = getSubdomainCorrectionDofs(i);
        for (std::size_t j = 0; j < dofs.size(); ++j)
        {
            if (dofs[j] < 0 || dofs[j] >= n || (j > 0 && dofs[j] <= dofs[j - 1]))
            {
                TBOX_ERROR("Multiplicative correction DOFs must be sorted, unique and in range.\n");
            }
            std::vector<Location>::const_iterator it =
                std::lower_bound(locations.cbegin(),
                                 locations.cend(),
                                 dofs[j],
                                 [](const Location& location, PetscInt column) { return location.column < column; });
            for (; it != locations.cend() && it->column == dofs[j]; ++it)
            {
                selected.push_back({ it->row, it->entry, static_cast<PetscInt>(j) });
            }
        }
        std::sort(selected.begin(),
                  selected.end(),
                  [](const CorrectionLocation& a, const CorrectionLocation& b)
                  { return std::tie(a.row, a.entry) < std::tie(b.row, b.entry); });
        const std::size_t first = d_rows.size();
        for (std::size_t j = 0; j < selected.size();)
        {
            const PetscInt row = selected[j].row;
            const std::size_t begin = d_entries.size();
            do
            {
                d_entries.push_back({ selected[j].entry, selected[j].correction });
                ++j;
            } while (j < selected.size() && selected[j].row == row);
            d_rows.push_back({ row, begin, d_entries.size() });
        }
        d_updates.push_back({ first, d_rows.size() });
    }
}

void
PETScLevelSolverShellBackend::updateResidual(const std::size_t i)
{
    copySubdomainCorrection(i, d_values.data());
    const std::vector<PetscInt>& dofs = getSubdomainCorrectionDofs(i);
    if (!d_use_rows)
    {
        int ierr = VecZeroEntries(d_correction);
        IBTK_CHKERRQ(ierr);
        // Match full overlap accumulation, including contributions from other ranks.
        ierr = VecSetValues(d_correction, static_cast<PetscInt>(dofs.size()), dofs.data(), d_values.data(), ADD_VALUES);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyBegin(d_correction);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyEnd(d_correction);
        IBTK_CHKERRQ(ierr);
        ierr = MatMult(d_mat, d_correction, d_action);
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPY(d_residual, -1.0, d_action);
        IBTK_CHKERRQ(ierr);
        return;
    }
    PetscScalar* residual = nullptr;
    int ierr = VecGetArray(d_residual, &residual);
    IBTK_CHKERRQ(ierr);
    for (std::size_t j = d_updates[i].row_begin; j < d_updates[i].row_end; ++j)
    {
        PetscInt count = 0;
        const PetscInt* columns = nullptr;
        const PetscScalar* values = nullptr;
        ierr = MatGetRow(d_mat, d_rows[j].global_row, &count, &columns, &values);
        IBTK_CHKERRQ(ierr);
        PetscScalar action = 0.0;
        for (std::size_t k = d_rows[j].entry_begin; k < d_rows[j].entry_end; ++k)
        {
            if (d_entries[k].matrix_entry >= count ||
                columns[d_entries[k].matrix_entry] != dofs[d_entries[k].correction_entry])
            {
                TBOX_ERROR("level-operator sparsity changed without shell-backend reinitialization.\n");
            }
            action += values[d_entries[k].matrix_entry] * d_values[d_entries[k].correction_entry];
        }
        residual[d_rows[j].global_row] -= action;
        ierr = MatRestoreRow(d_mat, d_rows[j].global_row, &count, &columns, &values);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecRestoreArray(d_residual, &residual);
    IBTK_CHKERRQ(ierr);
}
} // namespace IBTK
