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
#include <ibtk/private/PETScLevelSolverBlasLapackLUShellBackend.h>

#include <tbox/Utilities.h>

#include <cstddef>
#include <memory>
#include <vector>

namespace IBTK
{
namespace
{
std::unique_ptr<PETScLevelSolverShellBackend>
allocate_blas_lapack_lu_shell_backend(PETScLevelSolver& /*solver*/)
{
    return std::make_unique<PETScLevelSolverBlasLapackLUShellBackend>();
}

class PETScLevelSolverBlasLapackLUShellBackendRegistrar
{
public:
    PETScLevelSolverBlasLapackLUShellBackendRegistrar()
    {
        PETScLevelSolverShellBackendManager::getManager()->registerShellBackendFactoryFunction(
            "blas-lapack-lu", allocate_blas_lapack_lu_shell_backend);
    }
};

static PETScLevelSolverBlasLapackLUShellBackendRegistrar blas_lapack_lu_shell_backend_registrar;

std::size_t
column_major_index(const PetscBLASInt row, const PetscBLASInt col, const PetscBLASInt leading_dimension)
{
    return static_cast<std::size_t>(row) +
           static_cast<std::size_t>(col) * static_cast<std::size_t>(leading_dimension);
}
} // namespace

const std::string PETScLevelSolverBlasLapackLUShellBackend::s_backend_name = "BlasLapackLU";

const std::string&
PETScLevelSolverBlasLapackLUShellBackend::getName() const
{
    return s_backend_name;
}

void
PETScLevelSolverBlasLapackLUShellBackend::initializeSolverState(
    const PETScLevelSolverShellBackendState& solver_state)
{
    if (IBTK_MPI::getNodes() != 1)
    {
        TBOX_ERROR(solver_state.object_name << " " << solver_state.options_prefix
                                            << " BLAS/LAPACK LU shell backend currently requires one MPI rank.\n");
    }
#if defined(PETSC_USE_COMPLEX)
    TBOX_ERROR("PETScLevelSolverBlasLapackLUShellBackend requires real PETSc scalars.\n");
#endif

    TBOX_ASSERT(!d_data);
    d_solver_state = solver_state;
    d_data = std::make_unique<Data>();
    auto& data = *d_data;
    PetscInt n_columns = 0;
    int ierr = MatGetSize(d_solver_state.petsc_mat, &data.n_dofs, &n_columns);
    IBTK_CHKERRQ(ierr);
    if (data.n_dofs != n_columns)
    {
        TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix
                                              << " BLAS/LAPACK LU shell backend requires a square operator.\n");
    }
    PetscInt local_rows = 0;
    ierr = MatGetLocalSize(d_solver_state.petsc_mat, &local_rows, nullptr);
    IBTK_CHKERRQ(ierr);
    if (local_rows != data.n_dofs)
    {
        TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix
                                              << " BLAS/LAPACK LU shell backend currently requires one MPI rank.\n");
    }
    if (d_solver_state.use_multiplicative)
    {
        ierr = VecDuplicate(d_solver_state.petsc_x, &data.residual);
        IBTK_CHKERRQ(ierr);
        ierr = VecDuplicate(d_solver_state.petsc_x, &data.patch_correction);
        IBTK_CHKERRQ(ierr);
        ierr = VecDuplicate(d_solver_state.petsc_b, &data.residual_update);
        IBTK_CHKERRQ(ierr);
    }

    const auto& overlap_subdomains = *d_solver_state.subdomain_dofs;
    const auto& nonoverlap_subdomains = *d_solver_state.nonoverlap_subdomain_dofs;
    TBOX_ASSERT(overlap_subdomains.size() == nonoverlap_subdomains.size());
    data.subdomains.resize(overlap_subdomains.size());
    std::vector<PetscBLASInt> local_position(static_cast<std::size_t>(data.n_dofs), -1);

    for (std::size_t subdomain_num = 0; subdomain_num < data.subdomains.size(); ++subdomain_num)
    {
        auto& subdomain_data = data.subdomains[subdomain_num];
        subdomain_data.overlap_dofs = &overlap_subdomains[subdomain_num];
        subdomain_data.update_dofs = d_solver_state.use_restrict_partition ? &nonoverlap_subdomains[subdomain_num] :
                                                                          &overlap_subdomains[subdomain_num];
        const auto& overlap_dofs = *subdomain_data.overlap_dofs;
        ierr = PetscBLASIntCast(static_cast<PetscInt>(overlap_dofs.size()), &subdomain_data.local_size);
        IBTK_CHKERRQ(ierr);
        const PetscBLASInt local_size = subdomain_data.local_size;
        subdomain_data.lu_factor.assign(static_cast<std::size_t>(local_size) * static_cast<std::size_t>(local_size),
                                        0.0);
        subdomain_data.pivots.resize(static_cast<std::size_t>(local_size));
        subdomain_data.rhs_workspace.resize(static_cast<std::size_t>(local_size));

        for (PetscBLASInt local_dof = 0; local_dof < local_size; ++local_dof)
        {
            const int global_dof = overlap_dofs[static_cast<std::size_t>(local_dof)];
            TBOX_ASSERT(global_dof >= 0 && global_dof < data.n_dofs);
            local_position[static_cast<std::size_t>(global_dof)] = local_dof;
        }
        for (PetscBLASInt local_row = 0; local_row < local_size; ++local_row)
        {
            const PetscInt global_row = overlap_dofs[static_cast<std::size_t>(local_row)];
            PetscInt row_nnz = 0;
            const PetscInt* columns = nullptr;
            const PetscScalar* values = nullptr;
            ierr = MatGetRow(d_solver_state.petsc_mat, global_row, &row_nnz, &columns, &values);
            IBTK_CHKERRQ(ierr);
            for (PetscInt entry = 0; entry < row_nnz; ++entry)
            {
                const PetscInt global_col = columns[entry];
                if (global_col < 0 || global_col >= data.n_dofs) continue;
                const PetscBLASInt local_col = local_position[static_cast<std::size_t>(global_col)];
                if (local_col >= 0)
                {
                    subdomain_data.lu_factor[column_major_index(local_row, local_col, local_size)] = values[entry];
                }
            }
            ierr = MatRestoreRow(d_solver_state.petsc_mat, global_row, &row_nnz, &columns, &values);
            IBTK_CHKERRQ(ierr);
        }

        subdomain_data.update_local_positions.reserve(subdomain_data.update_dofs->size());
        for (const int update_dof : *subdomain_data.update_dofs)
        {
            const PetscBLASInt position = local_position[static_cast<std::size_t>(update_dof)];
            TBOX_ASSERT(position >= 0);
            subdomain_data.update_local_positions.push_back(position);
        }
        for (const int global_dof : overlap_dofs)
        {
            local_position[static_cast<std::size_t>(global_dof)] = -1;
        }

        if (local_size == 0) continue;
        PetscBLASInt info = 0;
        LAPACKgetrf_(&subdomain_data.local_size,
                     &subdomain_data.local_size,
                     subdomain_data.lu_factor.data(),
                     &subdomain_data.local_size,
                     subdomain_data.pivots.data(),
                     &info);
        if (info != 0)
        {
            TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix
                                                  << " BLAS/LAPACK LU factorization failed for subdomain "
                                                  << subdomain_num << " with info = " << info << ".");
        }
    }
}

void
PETScLevelSolverBlasLapackLUShellBackend::deallocateSolverState()
{
    if (d_data)
    {
        int ierr;
        if (d_data->residual)
        {
            ierr = VecDestroy(&d_data->residual);
            IBTK_CHKERRQ(ierr);
        }
        if (d_data->patch_correction)
        {
            ierr = VecDestroy(&d_data->patch_correction);
            IBTK_CHKERRQ(ierr);
        }
        if (d_data->residual_update)
        {
            ierr = VecDestroy(&d_data->residual_update);
            IBTK_CHKERRQ(ierr);
        }
    }
    d_data.reset();
    d_solver_state = PETScLevelSolverShellBackendState();
}

void
PETScLevelSolverBlasLapackLUShellBackend::solveSubdomain(SubdomainData& subdomain_data,
                                                         const std::size_t subdomain_num) const
{
    if (subdomain_data.local_size == 0) return;
    const char trans = 'N';
    const PetscBLASInt nrhs = 1;
    PetscBLASInt info = 0;
    LAPACKgetrs_(&trans,
                 &subdomain_data.local_size,
                 &nrhs,
                 subdomain_data.lu_factor.data(),
                 &subdomain_data.local_size,
                 subdomain_data.pivots.data(),
                 subdomain_data.rhs_workspace.data(),
                 &subdomain_data.local_size,
                 &info);
    if (info != 0)
    {
        TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix
                                              << " BLAS/LAPACK LU solve failed for subdomain " << subdomain_num
                                              << " with info = " << info << ".\n");
    }
}

bool
PETScLevelSolverBlasLapackLUShellBackend::shouldObserveSubdomain(const std::size_t subdomain_num) const
{
    if (!d_solver_state.subdomain_solve_observer || !*d_solver_state.subdomain_solve_observer) return false;
    return !d_solver_state.subdomain_solve_observer_predicate ||
           !*d_solver_state.subdomain_solve_observer_predicate ||
           (*d_solver_state.subdomain_solve_observer_predicate)(static_cast<int>(subdomain_num));
}

void
PETScLevelSolverBlasLapackLUShellBackend::observeSubdomain(const std::size_t subdomain_num,
                                                          SubdomainData& subdomain_data,
                                                          std::vector<PetscScalar>& local_rhs,
                                                          Vec current_global_source) const
{
    const auto& overlap_dofs = *subdomain_data.overlap_dofs;
    const PetscInt local_size = static_cast<PetscInt>(overlap_dofs.size());
    std::vector<PetscInt> petsc_dofs(overlap_dofs.begin(), overlap_dofs.end());
    IS subdomain_is = nullptr;
    int ierr = ISCreateGeneral(PetscObjectComm(reinterpret_cast<PetscObject>(d_solver_state.petsc_mat)),
                               local_size,
                               petsc_dofs.data(),
                               PETSC_COPY_VALUES,
                               &subdomain_is);
    IBTK_CHKERRQ(ierr);
    Mat local_matrix = nullptr;
    ierr = MatCreateSubMatrix(
        d_solver_state.petsc_mat, subdomain_is, subdomain_is, MAT_INITIAL_MATRIX, &local_matrix);
    IBTK_CHKERRQ(ierr);
    Vec rhs = nullptr, solution = nullptr;
    ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, local_size, local_rhs.data(), &rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecCreateSeqWithArray(
        PETSC_COMM_SELF, 1, local_size, subdomain_data.rhs_workspace.data(), &solution);
    IBTK_CHKERRQ(ierr);
    (*d_solver_state.subdomain_solve_observer)(
        static_cast<int>(subdomain_num), local_matrix, rhs, solution, current_global_source);
    ierr = VecDestroy(&solution);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&rhs);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&local_matrix);
    IBTK_CHKERRQ(ierr);
    ierr = ISDestroy(&subdomain_is);
    IBTK_CHKERRQ(ierr);
}

void
PETScLevelSolverBlasLapackLUShellBackend::applyAdditive(Vec x, Vec y)
{
    auto& data = *d_data;
    int ierr = VecZeroEntries(y);
    IBTK_CHKERRQ(ierr);
    for (std::size_t subdomain_num = 0; subdomain_num < data.subdomains.size(); ++subdomain_num)
    {
        auto& subdomain_data = data.subdomains[subdomain_num];
        const PetscScalar* x_values = nullptr;
        ierr = VecGetArrayRead(x, &x_values);
        IBTK_CHKERRQ(ierr);
        for (std::size_t local_dof = 0; local_dof < subdomain_data.overlap_dofs->size(); ++local_dof)
        {
            subdomain_data.rhs_workspace[local_dof] = x_values[(*subdomain_data.overlap_dofs)[local_dof]];
        }
        ierr = VecRestoreArrayRead(x, &x_values);
        IBTK_CHKERRQ(ierr);
        std::vector<PetscScalar> observed_rhs;
        if (shouldObserveSubdomain(subdomain_num)) observed_rhs = subdomain_data.rhs_workspace;
        solveSubdomain(subdomain_data, subdomain_num);
        if (!observed_rhs.empty()) observeSubdomain(subdomain_num, subdomain_data, observed_rhs, x);
        PetscScalar* y_values = nullptr;
        ierr = VecGetArray(y, &y_values);
        IBTK_CHKERRQ(ierr);
        for (std::size_t update_num = 0; update_num < subdomain_data.update_dofs->size(); ++update_num)
        {
            y_values[(*subdomain_data.update_dofs)[update_num]] +=
                subdomain_data.rhs_workspace[static_cast<std::size_t>(
                    subdomain_data.update_local_positions[update_num])];
        }
        ierr = VecRestoreArray(y, &y_values);
        IBTK_CHKERRQ(ierr);
    }
}

void
PETScLevelSolverBlasLapackLUShellBackend::applyMultiplicative(Vec x, Vec y)
{
    auto& data = *d_data;
    int ierr = VecZeroEntries(y);
    IBTK_CHKERRQ(ierr);
    ierr = VecCopy(x, data.residual);
    IBTK_CHKERRQ(ierr);
    for (std::size_t subdomain_num = 0; subdomain_num < data.subdomains.size(); ++subdomain_num)
    {
        auto& subdomain_data = data.subdomains[subdomain_num];
        const PetscScalar* residual_values = nullptr;
        ierr = VecGetArrayRead(data.residual, &residual_values);
        IBTK_CHKERRQ(ierr);
        for (std::size_t local_dof = 0; local_dof < subdomain_data.overlap_dofs->size(); ++local_dof)
        {
            subdomain_data.rhs_workspace[local_dof] =
                residual_values[(*subdomain_data.overlap_dofs)[local_dof]];
        }
        ierr = VecRestoreArrayRead(data.residual, &residual_values);
        IBTK_CHKERRQ(ierr);

        std::vector<PetscScalar> observed_rhs;
        if (shouldObserveSubdomain(subdomain_num)) observed_rhs = subdomain_data.rhs_workspace;
        solveSubdomain(subdomain_data, subdomain_num);
        if (!observed_rhs.empty()) observeSubdomain(subdomain_num, subdomain_data, observed_rhs, data.residual);
        PetscScalar* y_values = nullptr;
        ierr = VecGetArray(y, &y_values);
        IBTK_CHKERRQ(ierr);
        for (std::size_t update_num = 0; update_num < subdomain_data.update_dofs->size(); ++update_num)
        {
            y_values[(*subdomain_data.update_dofs)[update_num]] +=
                subdomain_data.rhs_workspace[static_cast<std::size_t>(
                    subdomain_data.update_local_positions[update_num])];
        }
        ierr = VecRestoreArray(y, &y_values);
        IBTK_CHKERRQ(ierr);

        if (subdomain_num + 1 < data.subdomains.size())
        {
            ierr = VecZeroEntries(data.patch_correction);
            IBTK_CHKERRQ(ierr);
            PetscScalar* correction_values = nullptr;
            ierr = VecGetArray(data.patch_correction, &correction_values);
            IBTK_CHKERRQ(ierr);
            for (std::size_t update_num = 0; update_num < subdomain_data.update_dofs->size(); ++update_num)
            {
                correction_values[(*subdomain_data.update_dofs)[update_num]] =
                    subdomain_data.rhs_workspace[static_cast<std::size_t>(
                        subdomain_data.update_local_positions[update_num])];
            }
            ierr = VecRestoreArray(data.patch_correction, &correction_values);
            IBTK_CHKERRQ(ierr);
            ierr = MatMult(d_solver_state.petsc_mat, data.patch_correction, data.residual_update);
            IBTK_CHKERRQ(ierr);
            ierr = VecAXPY(data.residual, -1.0, data.residual_update);
            IBTK_CHKERRQ(ierr);
        }
    }
}

void
PETScLevelSolverBlasLapackLUShellBackend::apply(Vec x, Vec y)
{
    TBOX_ASSERT(d_data);
    if (d_solver_state.use_multiplicative)
        applyMultiplicative(x, y);
    else
        applyAdditive(x, y);
    d_solver_state.postprocess_result(y);
}
} // namespace IBTK
