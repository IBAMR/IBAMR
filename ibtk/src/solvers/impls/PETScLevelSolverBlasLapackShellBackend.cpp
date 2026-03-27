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

#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PETScLevelSolver.h>
#include <ibtk/private/PETScLevelSolverBlasLapackShellBackend.h>

#include <tbox/Database.h>

#include <petscmat.h>
#include <petscsys.h>
#include <petscvec.h>

#include <algorithm>
#include <cctype>
#include <string>
#include <unordered_map>
#include <vector>

namespace IBTK
{
namespace
{
std::string
to_lower(std::string value)
{
    std::transform(value.begin(),
                   value.end(),
                   value.begin(),
                   [](unsigned char c) -> unsigned char { return static_cast<unsigned char>(std::tolower(c)); });
    return value;
}

std::vector<PetscScalar>
copy_vec_to_std_vector(Vec vec, PetscInt n_entries)
{
    const PetscScalar* vec_arr = nullptr;
    int ierr = VecGetArrayRead(vec, &vec_arr);
    IBTK_CHKERRQ(ierr);
    std::vector<PetscScalar> values(static_cast<std::size_t>(n_entries), 0.0);
    std::copy_n(vec_arr, n_entries, values.begin());
    ierr = VecRestoreArrayRead(vec, &vec_arr);
    IBTK_CHKERRQ(ierr);
    return values;
}
} // namespace

PETScLevelSolverBlasLapackShellBackend::PETScLevelSolverBlasLapackShellBackend(PETScLevelSolver& solver)
    : d_context(solver)
{
}

const std::string&
PETScLevelSolverBlasLapackShellBackend::getTypeKey() const
{
    return d_type_key;
}

void
PETScLevelSolverBlasLapackShellBackend::configure(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
{
    d_subdomain_solver_type = "lu";
    if (!input_db) return;

    if (input_db->keyExists("blas_lapack_subdomain_solver_type"))
    {
        parseSolverType(input_db->getString("blas_lapack_subdomain_solver_type"));
    }
}

const char*
PETScLevelSolverBlasLapackShellBackend::getPCNameSuffixAdditive() const
{
    return "PC_AdditiveBlasLapack";
}

const char*
PETScLevelSolverBlasLapackShellBackend::getPCNameSuffixMultiplicative() const
{
    return "PC_MultiplicativeBlasLapack";
}

void
PETScLevelSolverBlasLapackShellBackend::initialize()
{
    TBOX_ASSERT(IBTK_MPI::getNodes() == 1);
#if defined(PETSC_USE_COMPLEX)
    TBOX_ERROR("PETScLevelSolverBlasLapackShellBackend::initialize():\n"
               << "  complex PETSc scalar builds are not supported by the blas-lapack shell backend.\n");
#endif

    PetscInt n_dofs = 0;
    int ierr = VecGetSize(d_context.getPETScXForBackend(), &n_dofs);
    IBTK_CHKERRQ(ierr);
    ierr = PetscBLASIntCast(n_dofs, &d_n_dofs);
    IBTK_CHKERRQ(ierr);

    const bool use_restrict_partition = d_context.useRestrictPartitionForBackend();
    const bool use_multiplicative = d_context.isShellMultiplicativeForBackend();
    Mat level_mat = d_context.getPETScMatForBackend();

    const auto& overlap_subdomain_dofs = d_context.getSubdomainDOFsForBackend();
    const auto& nonoverlap_subdomain_dofs = d_context.getNonoverlapSubdomainDOFsForBackend();
    d_subdomains.clear();
    d_subdomains.resize(overlap_subdomain_dofs.size());
    for (std::size_t subdomain_num = 0; subdomain_num < d_subdomains.size(); ++subdomain_num)
    {
        initializeSubdomainData(d_subdomains[subdomain_num],
                                overlap_subdomain_dofs[subdomain_num],
                                nonoverlap_subdomain_dofs[subdomain_num],
                                use_restrict_partition,
                                use_multiplicative,
                                level_mat);
        factorizeSubdomainMatrix(d_subdomains[subdomain_num]);
    }
}

void
PETScLevelSolverBlasLapackShellBackend::deallocate()
{
    d_subdomains.clear();
    d_n_dofs = 0;
}

void
PETScLevelSolverBlasLapackShellBackend::applyAdditive(Vec x, Vec y)
{
    TBOX_ASSERT(IBTK_MPI::getNodes() == 1);
    PetscScalar* y_arr = nullptr;
    const PetscScalar* x_arr = nullptr;
    int ierr = VecGetArray(y, &y_arr);
    IBTK_CHKERRQ(ierr);
    ierr = VecGetArrayRead(x, &x_arr);
    IBTK_CHKERRQ(ierr);

    std::fill_n(y_arr, d_n_dofs, 0.0);
    for (auto& subdomain_data : d_subdomains)
    {
        std::size_t rhs_idx = 0;
        for (const int dof : subdomain_data.overlap_dofs)
        {
            subdomain_data.rhs_workspace[rhs_idx++] = x_arr[dof];
        }
        subdomain_data.delta_workspace = subdomain_data.rhs_workspace;
        solveSubdomainSystem(subdomain_data);
        std::size_t update_idx = 0;
        for (const int dof : subdomain_data.update_dofs)
        {
            y_arr[dof] +=
                subdomain_data
                    .delta_workspace[static_cast<std::size_t>(subdomain_data.update_local_positions[update_idx++])];
        }
    }

    ierr = VecRestoreArrayRead(x, &x_arr);
    IBTK_CHKERRQ(ierr);
    ierr = VecRestoreArray(y, &y_arr);
    IBTK_CHKERRQ(ierr);
    d_context.postprocessShellResultForBackend(y);
}

void
PETScLevelSolverBlasLapackShellBackend::applyMultiplicative(Vec x, Vec y)
{
    TBOX_ASSERT(IBTK_MPI::getNodes() == 1);
    PetscScalar* y_arr = nullptr;
    int ierr = VecGetArray(y, &y_arr);
    IBTK_CHKERRQ(ierr);
    std::fill_n(y_arr, d_n_dofs, 0.0);
    ierr = VecRestoreArray(y, &y_arr);
    IBTK_CHKERRQ(ierr);

    std::vector<PetscScalar> residual = copy_vec_to_std_vector(x, d_n_dofs);
    ierr = VecGetArray(y, &y_arr);
    IBTK_CHKERRQ(ierr);
    for (std::size_t subdomain_num = 0; subdomain_num < d_subdomains.size(); ++subdomain_num)
    {
        auto& subdomain_data = d_subdomains[subdomain_num];
        std::size_t rhs_idx = 0;
        for (const int dof : subdomain_data.overlap_dofs)
        {
            subdomain_data.rhs_workspace[rhs_idx++] = residual[static_cast<std::size_t>(dof)];
        }
        subdomain_data.delta_workspace = subdomain_data.rhs_workspace;
        solveSubdomainSystem(subdomain_data);

        std::size_t update_idx = 0;
        for (const int dof : subdomain_data.update_dofs)
        {
            y_arr[dof] +=
                subdomain_data
                    .delta_workspace[static_cast<std::size_t>(subdomain_data.update_local_positions[update_idx++])];
        }

        if (subdomain_num + 1 < d_subdomains.size())
        {
            updateResidual(subdomain_data, residual);
        }
    }
    ierr = VecRestoreArray(y, &y_arr);
    IBTK_CHKERRQ(ierr);
    d_context.postprocessShellResultForBackend(y);
}

void
PETScLevelSolverBlasLapackShellBackend::parseSolverType(const std::string& type_name)
{
    d_subdomain_solver_type = to_lower(type_name);
    if (d_subdomain_solver_type != "lu")
    {
        TBOX_ERROR(d_context.getObjectNameForBackend()
                   << " " << d_context.getOptionsPrefixForBackend()
                   << " PETScLevelSolverBlasLapackShellBackend::configure()\n"
                   << "Unsupported blas_lapack_subdomain_solver_type = " << type_name << "\n"
                   << "Supported values: lu" << std::endl);
    }
}

void
PETScLevelSolverBlasLapackShellBackend::initializeSubdomainData(SubdomainData& subdomain_data,
                                                                const std::vector<int>& overlap_dofs,
                                                                const std::vector<int>& nonoverlap_dofs,
                                                                const bool use_restrict_partition,
                                                                const bool use_multiplicative,
                                                                Mat level_mat) const
{
    subdomain_data = SubdomainData();
    subdomain_data.overlap_dofs = overlap_dofs;
    const PetscBLASInt overlap_size = static_cast<PetscBLASInt>(overlap_dofs.size());
    subdomain_data.local_size = overlap_size;
    subdomain_data.local_lda = std::max<PetscBLASInt>(1, overlap_size);
    subdomain_data.local_operator_lu.resize(
        static_cast<std::size_t>(subdomain_data.local_lda) * static_cast<std::size_t>(overlap_size), 0.0);
    subdomain_data.pivot.resize(static_cast<std::size_t>(std::max<PetscBLASInt>(1, overlap_size)), 0);
    subdomain_data.rhs_workspace.resize(static_cast<std::size_t>(overlap_size), 0.0);
    subdomain_data.delta_workspace.resize(static_cast<std::size_t>(overlap_size), 0.0);

    std::unordered_map<int, PetscBLASInt> overlap_col_map;
    overlap_col_map.reserve(overlap_dofs.size());
    for (PetscBLASInt local_col = 0; local_col < overlap_size; ++local_col)
    {
        overlap_col_map[overlap_dofs[static_cast<std::size_t>(local_col)]] = local_col;
    }

    for (PetscBLASInt local_row = 0; local_row < overlap_size; ++local_row)
    {
        const int global_row = overlap_dofs[static_cast<std::size_t>(local_row)];
        PetscInt row_nnz = 0;
        const PetscInt* cols = nullptr;
        const PetscScalar* vals = nullptr;
        int ierr = MatGetRow(level_mat, global_row, &row_nnz, &cols, &vals);
        IBTK_CHKERRQ(ierr);
        for (PetscInt entry = 0; entry < row_nnz; ++entry)
        {
            const auto col_it = overlap_col_map.find(static_cast<int>(cols[entry]));
            if (col_it != overlap_col_map.end())
            {
                subdomain_data.local_operator_lu[static_cast<std::size_t>(local_row) +
                                                 static_cast<std::size_t>(col_it->second) *
                                                     static_cast<std::size_t>(subdomain_data.local_lda)] = vals[entry];
            }
        }
        ierr = MatRestoreRow(level_mat, global_row, &row_nnz, &cols, &vals);
        IBTK_CHKERRQ(ierr);
    }

    if (use_restrict_partition)
    {
        subdomain_data.update_dofs = nonoverlap_dofs;
        subdomain_data.update_local_positions.resize(nonoverlap_dofs.size(), 0);
        for (std::size_t local_idx = 0; local_idx < nonoverlap_dofs.size(); ++local_idx)
        {
            const auto pos_it = overlap_col_map.find(nonoverlap_dofs[local_idx]);
            TBOX_ASSERT(pos_it != overlap_col_map.end());
            subdomain_data.update_local_positions[local_idx] = static_cast<int>(pos_it->second);
        }
    }
    else
    {
        subdomain_data.update_dofs = overlap_dofs;
        subdomain_data.update_local_positions.resize(overlap_dofs.size(), 0);
        for (std::size_t local_idx = 0; local_idx < overlap_dofs.size(); ++local_idx)
        {
            subdomain_data.update_local_positions[local_idx] = static_cast<int>(local_idx);
        }
    }

    if (!use_multiplicative) return;

    std::vector<int> active_rows;
    std::unordered_map<int, bool> update_dof_map;
    update_dof_map.reserve(subdomain_data.update_dofs.size());
    for (const int dof : subdomain_data.update_dofs)
    {
        update_dof_map[dof] = true;
    }

    PetscInt n_rows = 0, n_cols = 0;
    int ierr = MatGetSize(level_mat, &n_rows, &n_cols);
    IBTK_CHKERRQ(ierr);
    for (PetscInt row = 0; row < n_rows; ++row)
    {
        PetscInt row_nnz = 0;
        const PetscInt* cols = nullptr;
        const PetscScalar* vals = nullptr;
        ierr = MatGetRow(level_mat, row, &row_nnz, &cols, &vals);
        IBTK_CHKERRQ(ierr);
        bool row_is_active = false;
        for (PetscInt entry = 0; entry < row_nnz && !row_is_active; ++entry)
        {
            row_is_active = update_dof_map.find(static_cast<int>(cols[entry])) != update_dof_map.end();
        }
        ierr = MatRestoreRow(level_mat, row, &row_nnz, &cols, &vals);
        IBTK_CHKERRQ(ierr);
        if (row_is_active) active_rows.push_back(static_cast<int>(row));
    }

    subdomain_data.active_residual_update_rows = std::move(active_rows);
    subdomain_data.active_residual_update_num_rows =
        static_cast<PetscBLASInt>(subdomain_data.active_residual_update_rows.size());
    subdomain_data.active_residual_update_num_cols = static_cast<PetscBLASInt>(subdomain_data.update_dofs.size());
    subdomain_data.active_residual_update_lda =
        std::max<PetscBLASInt>(1, subdomain_data.active_residual_update_num_rows);
    subdomain_data.active_residual_update_mat.resize(
        static_cast<std::size_t>(subdomain_data.active_residual_update_lda) *
            static_cast<std::size_t>(subdomain_data.active_residual_update_num_cols),
        0.0);
    subdomain_data.residual_input_workspace.resize(static_cast<std::size_t>(subdomain_data.update_dofs.size()), 0.0);
    subdomain_data.residual_delta_workspace.resize(
        static_cast<std::size_t>(subdomain_data.active_residual_update_num_rows), 0.0);

    for (PetscBLASInt local_row = 0; local_row < subdomain_data.active_residual_update_num_rows; ++local_row)
    {
        const int global_row = subdomain_data.active_residual_update_rows[static_cast<std::size_t>(local_row)];
        PetscInt row_nnz = 0;
        const PetscInt* cols = nullptr;
        const PetscScalar* vals = nullptr;
        ierr = MatGetRow(level_mat, global_row, &row_nnz, &cols, &vals);
        IBTK_CHKERRQ(ierr);
        for (PetscInt entry = 0; entry < row_nnz; ++entry)
        {
            const auto col_it = std::find(
                subdomain_data.update_dofs.begin(), subdomain_data.update_dofs.end(), static_cast<int>(cols[entry]));
            if (col_it != subdomain_data.update_dofs.end())
            {
                const PetscBLASInt local_col =
                    static_cast<PetscBLASInt>(std::distance(subdomain_data.update_dofs.begin(), col_it));
                subdomain_data.active_residual_update_mat[static_cast<std::size_t>(local_row) +
                                                          static_cast<std::size_t>(local_col) *
                                                              static_cast<std::size_t>(
                                                                  subdomain_data.active_residual_update_lda)] =
                    vals[entry];
            }
        }
        ierr = MatRestoreRow(level_mat, global_row, &row_nnz, &cols, &vals);
        IBTK_CHKERRQ(ierr);
    }
}

void
PETScLevelSolverBlasLapackShellBackend::factorizeSubdomainMatrix(SubdomainData& subdomain_data) const
{
    if (subdomain_data.local_size == 0) return;

    PetscBLASInt info = 0;
    LAPACKgetrf_(&subdomain_data.local_size,
                 &subdomain_data.local_size,
                 subdomain_data.local_operator_lu.data(),
                 &subdomain_data.local_lda,
                 subdomain_data.pivot.data(),
                 &info);
    if (info != 0)
    {
        TBOX_ERROR(d_context.getObjectNameForBackend()
                   << " " << d_context.getOptionsPrefixForBackend()
                   << " PETScLevelSolverBlasLapackShellBackend::factorizeSubdomainMatrix()\n"
                   << "LAPACKgetrf failed with info = " << info << ".\n");
    }
}

void
PETScLevelSolverBlasLapackShellBackend::solveSubdomainSystem(SubdomainData& subdomain_data) const
{
    if (subdomain_data.local_size == 0) return;

    const PetscBLASInt nrhs = 1;
    PetscBLASInt info = 0;
    const char trans = 'N';
    LAPACKgetrs_(&trans,
                 &subdomain_data.local_size,
                 &nrhs,
                 subdomain_data.local_operator_lu.data(),
                 &subdomain_data.local_lda,
                 subdomain_data.pivot.data(),
                 subdomain_data.delta_workspace.data(),
                 &subdomain_data.local_size,
                 &info);
    if (info != 0)
    {
        TBOX_ERROR(d_context.getObjectNameForBackend()
                   << " " << d_context.getOptionsPrefixForBackend()
                   << " PETScLevelSolverBlasLapackShellBackend::solveSubdomainSystem()\n"
                   << "LAPACKgetrs failed with info = " << info << ".\n");
    }
}

void
PETScLevelSolverBlasLapackShellBackend::updateResidual(SubdomainData& subdomain_data,
                                                       std::vector<PetscScalar>& residual) const
{
    if (subdomain_data.active_residual_update_num_rows == 0 || subdomain_data.active_residual_update_num_cols == 0)
    {
        return;
    }

    std::size_t input_idx = 0;
    for (const int local_pos : subdomain_data.update_local_positions)
    {
        subdomain_data.residual_input_workspace[input_idx++] =
            subdomain_data.delta_workspace[static_cast<std::size_t>(local_pos)];
    }

    const char trans = 'N';
    const PetscBLASInt inc = 1;
    const PetscScalar alpha = 1.0;
    const PetscScalar beta = 0.0;
    BLASgemv_(&trans,
              &subdomain_data.active_residual_update_num_rows,
              &subdomain_data.active_residual_update_num_cols,
              &alpha,
              subdomain_data.active_residual_update_mat.data(),
              &subdomain_data.active_residual_update_lda,
              subdomain_data.residual_input_workspace.data(),
              &inc,
              &beta,
              subdomain_data.residual_delta_workspace.data(),
              &inc);

    for (std::size_t row_idx = 0; row_idx < subdomain_data.active_residual_update_rows.size(); ++row_idx)
    {
        const int row = subdomain_data.active_residual_update_rows[row_idx];
        residual[static_cast<std::size_t>(row)] -= subdomain_data.residual_delta_workspace[row_idx];
    }
}
} // namespace IBTK
