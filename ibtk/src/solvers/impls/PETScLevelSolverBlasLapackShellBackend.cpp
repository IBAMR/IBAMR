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
#include <tbox/Utilities.h>

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace IBTK
{
namespace
{
std::unique_ptr<PETScLevelSolverShellBackend>
allocate_blas_lapack_shell_backend(PETScLevelSolver& /*solver*/)
{
    return std::make_unique<PETScLevelSolverBlasLapackShellBackend>("BlasLapack", true);
}

std::unique_ptr<PETScLevelSolverShellBackend>
allocate_blas_lapack_lu_shell_backend(PETScLevelSolver& /*solver*/)
{
    return std::make_unique<PETScLevelSolverBlasLapackShellBackend>("BlasLapackLU", false);
}

class PETScLevelSolverBlasLapackShellBackendRegistrar
{
public:
    PETScLevelSolverBlasLapackShellBackendRegistrar()
    {
        PETScLevelSolverShellBackendManager::getManager()->registerShellBackendFactoryFunction(
            "blas-lapack", allocate_blas_lapack_shell_backend);
        PETScLevelSolverShellBackendManager::getManager()->registerShellBackendFactoryFunction(
            "blas-lapack-lu", allocate_blas_lapack_lu_shell_backend);
    }
};

static PETScLevelSolverBlasLapackShellBackendRegistrar blas_lapack_shell_backend_registrar;

std::string
to_lower(std::string value)
{
    std::transform(value.begin(),
                   value.end(),
                   value.begin(),
                   [](const unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

std::size_t
column_major_index(const PetscBLASInt row, const PetscBLASInt col, const PetscBLASInt leading_dimension)
{
    return static_cast<std::size_t>(row) + static_cast<std::size_t>(col) * static_cast<std::size_t>(leading_dimension);
}

PetscBLASInt
query_lapack_workspace_size(const PetscScalar query_value)
{
    return std::max<PetscBLASInt>(1, static_cast<PetscBLASInt>(PetscRealPart(query_value)));
}
} // namespace

PETScLevelSolverBlasLapackShellBackend::PETScLevelSolverBlasLapackShellBackend(std::string backend_name,
                                                                               const bool configurable_solver_type)
    : d_backend_name(std::move(backend_name)), d_configurable_solver_type(configurable_solver_type)
{
}

const std::string&
PETScLevelSolverBlasLapackShellBackend::getName() const
{
    return d_backend_name;
}

void
PETScLevelSolverBlasLapackShellBackend::configureFromInputDatabase()
{
    d_subdomain_solver_rcond = -1.0;
    if (!d_configurable_solver_type)
    {
        d_subdomain_solver_type = SubdomainSolverType::LU;
        if (d_solver_state.input_db && d_solver_state.input_db->keyExists("blas_lapack_subdomain_solver_type") &&
            to_lower(d_solver_state.input_db->getString("blas_lapack_subdomain_solver_type")) != "lu")
        {
            TBOX_ERROR(d_solver_state.object_name
                       << " " << d_solver_state.options_prefix
                       << " shell backend 'blas-lapack-lu' has a fixed LU local solver; "
                       << "select shell backend 'blas-lapack' to configure another BLAS/LAPACK mode.");
        }
        return;
    }

    d_subdomain_solver_type = SubdomainSolverType::SVD;
    if (!d_solver_state.input_db) return;
    if (d_solver_state.input_db->keyExists("blas_lapack_subdomain_solver_type"))
    {
        const std::string solver_type =
            to_lower(d_solver_state.input_db->getString("blas_lapack_subdomain_solver_type"));
        if (solver_type == "svd")
            d_subdomain_solver_type = SubdomainSolverType::SVD;
        else if (solver_type == "lu")
            d_subdomain_solver_type = SubdomainSolverType::LU;
        else if (solver_type == "symmetric-indefinite")
            d_subdomain_solver_type = SubdomainSolverType::SYMMETRIC_INDEFINITE;
        else if (solver_type == "qr")
            d_subdomain_solver_type = SubdomainSolverType::QR;
        else if (solver_type == "cholesky")
            TBOX_ERROR(d_solver_state.object_name
                       << " " << d_solver_state.options_prefix
                       << " Cholesky is not supported for indefinite Stokes subdomain matrices; "
                       << "use 'symmetric-indefinite' when a symmetry-specific factorization is required.");
        else
            TBOX_ERROR(d_solver_state.object_name
                       << " " << d_solver_state.options_prefix << " unsupported blas_lapack_subdomain_solver_type = "
                       << solver_type << "; supported values are 'svd', 'lu', 'symmetric-indefinite', and 'qr'.");
    }
    if (d_solver_state.input_db->keyExists("blas_lapack_subdomain_solver_rcond"))
    {
        d_subdomain_solver_rcond = d_solver_state.input_db->getDouble("blas_lapack_subdomain_solver_rcond");
    }
}

void
PETScLevelSolverBlasLapackShellBackend::initializeSolverState(const PETScLevelSolverShellBackendState& solver_state)
{
    if (IBTK_MPI::getNodes() != 1)
    {
        TBOX_ERROR(solver_state.object_name << " " << solver_state.options_prefix
                                            << " BLAS/LAPACK shell backend currently requires one MPI rank.\n");
    }
#if defined(PETSC_USE_COMPLEX)
    TBOX_ERROR("PETScLevelSolverBlasLapackShellBackend requires real PETSc scalars.\n");
#endif

    TBOX_ASSERT(!d_data);
    d_solver_state = solver_state;
    configureFromInputDatabase();
    d_data = std::make_unique<Data>();
    auto& data = *d_data;
    PetscInt n_columns = 0;
    int ierr = MatGetSize(d_solver_state.petsc_mat, &data.n_dofs, &n_columns);
    IBTK_CHKERRQ(ierr);
    if (data.n_dofs != n_columns)
    {
        TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix
                                              << " BLAS/LAPACK shell backend requires a square operator.\n");
    }
    PetscInt local_rows = 0;
    ierr = MatGetLocalSize(d_solver_state.petsc_mat, &local_rows, nullptr);
    IBTK_CHKERRQ(ierr);
    if (local_rows != data.n_dofs)
    {
        TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix
                                              << " BLAS/LAPACK shell backend currently requires one MPI rank.\n");
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
        subdomain_data.solve_data.assign(static_cast<std::size_t>(local_size) * static_cast<std::size_t>(local_size),
                                         0.0);
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
                    subdomain_data.solve_data[column_major_index(local_row, local_col, local_size)] = values[entry];
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

        initializeSubdomainSolver(subdomain_data, subdomain_num);
    }
}

void
PETScLevelSolverBlasLapackShellBackend::initializeSubdomainSolver(SubdomainData& subdomain_data,
                                                                  const std::size_t subdomain_num)
{
    if (subdomain_data.local_size == 0) return;
    PetscBLASInt info = 0;
    switch (d_subdomain_solver_type)
    {
    case SubdomainSolverType::LU:
        subdomain_data.pivots.resize(static_cast<std::size_t>(subdomain_data.local_size));
        LAPACKgetrf_(&subdomain_data.local_size,
                     &subdomain_data.local_size,
                     subdomain_data.solve_data.data(),
                     &subdomain_data.local_size,
                     subdomain_data.pivots.data(),
                     &info);
        if (info != 0)
        {
            TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix
                                                  << " BLAS/LAPACK LU factorization failed for subdomain "
                                                  << subdomain_num << " with info = " << info << ".");
        }
        return;
    case SubdomainSolverType::SYMMETRIC_INDEFINITE:
    {
        verifySymmetricSubdomainMatrix(subdomain_data, subdomain_num);
        // SYTRF overwrites the only persistent matrix-sized object. Its
        // workspace is setup-only and is released after factorization.
        subdomain_data.pivots.resize(static_cast<std::size_t>(subdomain_data.local_size));
        const char uplo = 'L';
        PetscBLASInt lwork = -1;
        PetscScalar work_query = 0.0;
        LAPACKsytrf_(&uplo,
                     &subdomain_data.local_size,
                     subdomain_data.solve_data.data(),
                     &subdomain_data.local_size,
                     subdomain_data.pivots.data(),
                     &work_query,
                     &lwork,
                     &info);
        if (info != 0)
        {
            TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix
                                                  << " BLAS/LAPACK symmetric-indefinite workspace query failed for "
                                                  << "subdomain " << subdomain_num << " with info = " << info << ".");
        }
        lwork = query_lapack_workspace_size(work_query);
        std::vector<PetscScalar> work(static_cast<std::size_t>(lwork));
        LAPACKsytrf_(&uplo,
                     &subdomain_data.local_size,
                     subdomain_data.solve_data.data(),
                     &subdomain_data.local_size,
                     subdomain_data.pivots.data(),
                     work.data(),
                     &lwork,
                     &info);
        if (info != 0)
        {
            TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix
                                                  << " BLAS/LAPACK symmetric-indefinite factorization failed for "
                                                  << "subdomain " << subdomain_num << " with info = " << info << ".");
        }
        return;
    }
    case SubdomainSolverType::QR:
        initializeQRSolver(subdomain_data, subdomain_num);
        return;
    case SubdomainSolverType::SVD:
        initializeSVDSolver(subdomain_data, subdomain_num);
        return;
    }
}

void
PETScLevelSolverBlasLapackShellBackend::initializeQRSolver(SubdomainData& subdomain_data,
                                                           const std::size_t subdomain_num)
{
    const PetscBLASInt n = subdomain_data.local_size;
    // For A = Q R, form A^{-1} = R^{-1} Q^T directly in the persistent
    // solve-data buffer. The QR factor and LAPACK work arrays are setup-only.
    std::vector<PetscScalar> qr_factor = std::move(subdomain_data.solve_data);
    std::vector<PetscScalar> tau(static_cast<std::size_t>(n));
    PetscBLASInt info = 0;
    PetscBLASInt lwork = -1;
    PetscScalar work_query = 0.0;
    LAPACKgeqrf_(&n, &n, qr_factor.data(), &n, tau.data(), &work_query, &lwork, &info);
    if (info != 0)
    {
        TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix
                                              << " BLAS/LAPACK QR workspace query failed for subdomain "
                                              << subdomain_num << " with info = " << info << ".");
    }
    lwork = query_lapack_workspace_size(work_query);
    std::vector<PetscScalar> work(static_cast<std::size_t>(lwork));
    LAPACKgeqrf_(&n, &n, qr_factor.data(), &n, tau.data(), work.data(), &lwork, &info);
    if (info != 0)
    {
        TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix
                                              << " BLAS/LAPACK QR factorization failed for subdomain " << subdomain_num
                                              << " with info = " << info << ".");
    }

    PetscReal max_abs_diagonal = 0.0;
    for (PetscBLASInt k = 0; k < n; ++k)
    {
        max_abs_diagonal = std::max(max_abs_diagonal, PetscAbsScalar(qr_factor[column_major_index(k, k, n)]));
    }
    if (max_abs_diagonal == 0.0)
    {
        TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix
                                              << " BLAS/LAPACK QR factorization detected a zero-rank subdomain "
                                              << subdomain_num << ".");
    }
    if (d_subdomain_solver_rcond >= 0.0)
    {
        // Unpivoted QR is accepted only as a full-rank solver. This threshold
        // checks the computed R diagonal; SVD is the rank-deficient mode.
        for (PetscBLASInt k = 0; k < n; ++k)
        {
            if (PetscAbsScalar(qr_factor[column_major_index(k, k, n)]) <= d_subdomain_solver_rcond * max_abs_diagonal)
            {
                TBOX_ERROR(d_solver_state.object_name
                           << " " << d_solver_state.options_prefix
                           << " BLAS/LAPACK QR factorization detected a rank-deficient subdomain " << subdomain_num
                           << " with blas_lapack_subdomain_solver_rcond = " << d_subdomain_solver_rcond << ".");
            }
        }
    }

    subdomain_data.solve_data.assign(static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 0.0);
    for (PetscBLASInt k = 0; k < n; ++k)
    {
        subdomain_data.solve_data[column_major_index(k, k, n)] = 1.0;
    }
    const char side = 'L';
    const char transpose = 'T';
    lwork = -1;
    work_query = 0.0;
    LAPACKormqr_(&side,
                 &transpose,
                 &n,
                 &n,
                 &n,
                 qr_factor.data(),
                 &n,
                 tau.data(),
                 subdomain_data.solve_data.data(),
                 &n,
                 &work_query,
                 &lwork,
                 &info);
    if (info != 0)
    {
        TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix
                                              << " BLAS/LAPACK QR-application workspace query failed for subdomain "
                                              << subdomain_num << " with info = " << info << ".");
    }
    lwork = query_lapack_workspace_size(work_query);
    work.assign(static_cast<std::size_t>(lwork), 0.0);
    LAPACKormqr_(&side,
                 &transpose,
                 &n,
                 &n,
                 &n,
                 qr_factor.data(),
                 &n,
                 tau.data(),
                 subdomain_data.solve_data.data(),
                 &n,
                 work.data(),
                 &lwork,
                 &info);
    if (info != 0)
    {
        TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix
                                              << " BLAS/LAPACK QR application failed for subdomain " << subdomain_num
                                              << " with info = " << info << ".");
    }
    const char uplo = 'U';
    const char no_transpose = 'N';
    const char nonunit = 'N';
    LAPACKtrtrs_(
        &uplo, &no_transpose, &nonunit, &n, &n, qr_factor.data(), &n, subdomain_data.solve_data.data(), &n, &info);
    if (info != 0)
    {
        TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix
                                              << " BLAS/LAPACK triangular QR solve failed for subdomain "
                                              << subdomain_num << " with info = " << info << ".");
    }
    subdomain_data.solution_workspace.resize(static_cast<std::size_t>(n));
}

void
PETScLevelSolverBlasLapackShellBackend::initializeSVDSolver(SubdomainData& subdomain_data,
                                                            const std::size_t subdomain_num)
{
#if defined(PETSC_USE_COMPLEX)
    TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix
                                          << " BLAS/LAPACK SVD setup requires real PETSc scalars.");
#else
    const PetscBLASInt n = subdomain_data.local_size;
    // Solve A X = I with GELSS once during setup. The overwritten right-hand
    // side is A^+ and becomes the only persistent matrix-sized solve object.
    std::vector<PetscScalar> svd_factor = std::move(subdomain_data.solve_data);
    subdomain_data.solve_data.assign(static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 0.0);
    for (PetscBLASInt k = 0; k < n; ++k)
    {
        subdomain_data.solve_data[column_major_index(k, k, n)] = 1.0;
    }
    std::vector<PetscReal> singular_values(static_cast<std::size_t>(n));
    const PetscReal rcond = d_subdomain_solver_rcond >= 0.0 ? d_subdomain_solver_rcond : -1.0;
    PetscBLASInt effective_rank = 0;
    PetscBLASInt info = 0;
    PetscBLASInt lwork = -1;
    PetscScalar work_query = 0.0;
    LAPACKgelss_(&n,
                 &n,
                 &n,
                 svd_factor.data(),
                 &n,
                 subdomain_data.solve_data.data(),
                 &n,
                 singular_values.data(),
                 &rcond,
                 &effective_rank,
                 &work_query,
                 &lwork,
                 &info);
    if (info != 0)
    {
        TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix
                                              << " BLAS/LAPACK SVD workspace query failed for subdomain "
                                              << subdomain_num << " with info = " << info << ".");
    }
    lwork = query_lapack_workspace_size(work_query);
    std::vector<PetscScalar> work(static_cast<std::size_t>(lwork));
    LAPACKgelss_(&n,
                 &n,
                 &n,
                 svd_factor.data(),
                 &n,
                 subdomain_data.solve_data.data(),
                 &n,
                 singular_values.data(),
                 &rcond,
                 &effective_rank,
                 work.data(),
                 &lwork,
                 &info);
    if (info != 0)
    {
        TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix
                                              << " BLAS/LAPACK SVD solve-matrix construction failed for subdomain "
                                              << subdomain_num << " with info = " << info
                                              << " and effective rank = " << effective_rank << ".");
    }
    subdomain_data.solution_workspace.resize(static_cast<std::size_t>(n));
#endif
}

void
PETScLevelSolverBlasLapackShellBackend::verifySymmetricSubdomainMatrix(const SubdomainData& subdomain_data,
                                                                       const std::size_t subdomain_num) const
{
    PetscReal max_entry = 0.0;
    PetscReal max_asymmetry = 0.0;
    for (PetscBLASInt row = 0; row < subdomain_data.local_size; ++row)
    {
        for (PetscBLASInt col = row; col < subdomain_data.local_size; ++col)
        {
            const PetscScalar upper =
                subdomain_data.solve_data[column_major_index(row, col, subdomain_data.local_size)];
            const PetscScalar lower =
                subdomain_data.solve_data[column_major_index(col, row, subdomain_data.local_size)];
            max_entry = std::max(max_entry, PetscAbsScalar(upper));
            max_entry = std::max(max_entry, PetscAbsScalar(lower));
            max_asymmetry = std::max(max_asymmetry, PetscAbsScalar(upper - lower));
        }
    }
    const PetscReal tolerance = 100.0 * std::numeric_limits<PetscReal>::epsilon() * std::max<PetscReal>(1.0, max_entry);
    if (max_asymmetry > tolerance)
    {
        TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix << " subdomain " << subdomain_num
                                              << " is not symmetric enough for the " << getSolverTypeName()
                                              << " solver; maximum asymmetry = " << max_asymmetry << ".");
    }
}

const char*
PETScLevelSolverBlasLapackShellBackend::getSolverTypeName() const
{
    switch (d_subdomain_solver_type)
    {
    case SubdomainSolverType::SVD:
        return "svd";
    case SubdomainSolverType::LU:
        return "lu";
    case SubdomainSolverType::SYMMETRIC_INDEFINITE:
        return "symmetric-indefinite";
    case SubdomainSolverType::QR:
        return "qr";
    }
    return "unknown";
}

void
PETScLevelSolverBlasLapackShellBackend::deallocateSolverState()
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
PETScLevelSolverBlasLapackShellBackend::solveSubdomainSystem(SubdomainData& subdomain_data,
                                                             const std::size_t subdomain_num) const
{
    if (subdomain_data.local_size == 0) return;
    const PetscBLASInt nrhs = 1;
    PetscBLASInt info = 0;
    switch (d_subdomain_solver_type)
    {
    case SubdomainSolverType::LU:
    {
        const char transpose = 'N';
        LAPACKgetrs_(&transpose,
                     &subdomain_data.local_size,
                     &nrhs,
                     subdomain_data.solve_data.data(),
                     &subdomain_data.local_size,
                     subdomain_data.pivots.data(),
                     subdomain_data.rhs_workspace.data(),
                     &subdomain_data.local_size,
                     &info);
        if (info != 0)
        {
            TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix
                                                  << " BLAS/LAPACK LU solve failed for subdomain " << subdomain_num
                                                  << " with info = " << info << ".");
        }
        return;
    }
    case SubdomainSolverType::SYMMETRIC_INDEFINITE:
    {
        const char uplo = 'L';
        LAPACKsytrs_(&uplo,
                     &subdomain_data.local_size,
                     &nrhs,
                     subdomain_data.solve_data.data(),
                     &subdomain_data.local_size,
                     subdomain_data.pivots.data(),
                     subdomain_data.rhs_workspace.data(),
                     &subdomain_data.local_size,
                     &info);
        if (info != 0)
        {
            TBOX_ERROR(d_solver_state.object_name << " " << d_solver_state.options_prefix
                                                  << " BLAS/LAPACK symmetric-indefinite solve failed for subdomain "
                                                  << subdomain_num << " with info = " << info << ".");
        }
        return;
    }
    case SubdomainSolverType::QR:
    case SubdomainSolverType::SVD:
    {
        const char transpose = 'N';
        const PetscBLASInt increment = 1;
        const PetscScalar alpha = 1.0;
        const PetscScalar beta = 0.0;
        BLASgemv_(&transpose,
                  &subdomain_data.local_size,
                  &subdomain_data.local_size,
                  &alpha,
                  subdomain_data.solve_data.data(),
                  &subdomain_data.local_size,
                  subdomain_data.rhs_workspace.data(),
                  &increment,
                  &beta,
                  subdomain_data.solution_workspace.data(),
                  &increment);
        // Swapping the equally sized workspaces avoids a per-patch solution
        // copy while keeping the next matrix-vector product nonaliasing.
        subdomain_data.rhs_workspace.swap(subdomain_data.solution_workspace);
        return;
    }
    }
}

void
PETScLevelSolverBlasLapackShellBackend::observeSubdomain(const std::size_t subdomain_num,
                                                         SubdomainData& subdomain_data,
                                                         Vec current_global_source) const
{
    const auto& overlap_dofs = *subdomain_data.overlap_dofs;
    const PetscInt local_size = static_cast<PetscInt>(overlap_dofs.size());
    std::vector<PetscScalar> local_rhs(static_cast<std::size_t>(local_size));
    const PetscScalar* source_values = nullptr;
    int ierr = VecGetArrayRead(current_global_source, &source_values);
    IBTK_CHKERRQ(ierr);
    for (std::size_t local_dof = 0; local_dof < overlap_dofs.size(); ++local_dof)
    {
        local_rhs[local_dof] = source_values[overlap_dofs[local_dof]];
    }
    ierr = VecRestoreArrayRead(current_global_source, &source_values);
    IBTK_CHKERRQ(ierr);

    std::vector<PetscInt> petsc_dofs(overlap_dofs.begin(), overlap_dofs.end());
    IS subdomain_is = nullptr;
    ierr = ISCreateGeneral(PetscObjectComm(reinterpret_cast<PetscObject>(d_solver_state.petsc_mat)),
                           local_size,
                           petsc_dofs.data(),
                           PETSC_COPY_VALUES,
                           &subdomain_is);
    IBTK_CHKERRQ(ierr);
    Mat local_matrix = nullptr;
    ierr = MatCreateSubMatrix(d_solver_state.petsc_mat, subdomain_is, subdomain_is, MAT_INITIAL_MATRIX, &local_matrix);
    IBTK_CHKERRQ(ierr);
    Vec rhs = nullptr, solution = nullptr;
    ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, local_size, local_rhs.data(), &rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, local_size, subdomain_data.rhs_workspace.data(), &solution);
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

std::size_t
PETScLevelSolverBlasLapackShellBackend::getNumberOfSubdomains() const
{
    return d_data->subdomains.size();
}

void
PETScLevelSolverBlasLapackShellBackend::initializeSubdomainSweep(Vec x, Vec /*y*/)
{
    TBOX_ASSERT(d_data);
    if (d_solver_state.use_multiplicative)
    {
        const int ierr = VecCopy(x, d_data->residual);
        IBTK_CHKERRQ(ierr);
    }
}

void
PETScLevelSolverBlasLapackShellBackend::beginSubdomainRhs(const std::size_t subdomain_num, Vec x, Vec /*y*/)
{
    auto& subdomain_data = d_data->subdomains[subdomain_num];
    Vec source = d_solver_state.use_multiplicative ? d_data->residual : x;
    const PetscScalar* source_values = nullptr;
    int ierr = VecGetArrayRead(source, &source_values);
    IBTK_CHKERRQ(ierr);
    for (std::size_t local_dof = 0; local_dof < subdomain_data.overlap_dofs->size(); ++local_dof)
    {
        subdomain_data.rhs_workspace[local_dof] = source_values[(*subdomain_data.overlap_dofs)[local_dof]];
    }
    ierr = VecRestoreArrayRead(source, &source_values);
    IBTK_CHKERRQ(ierr);
}

void
PETScLevelSolverBlasLapackShellBackend::endSubdomainRhs(std::size_t /*subdomain_num*/, Vec /*x*/, Vec /*y*/)
{
}

void
PETScLevelSolverBlasLapackShellBackend::solveSubdomain(const std::size_t subdomain_num)
{
    solveSubdomainSystem(d_data->subdomains[subdomain_num], subdomain_num);
}

void
PETScLevelSolverBlasLapackShellBackend::observeSubdomainSolve(const std::size_t subdomain_num, Vec x, Vec /*y*/)
{
    Vec source = d_solver_state.use_multiplicative ? d_data->residual : x;
    observeSubdomain(subdomain_num, d_data->subdomains[subdomain_num], source);
}

void
PETScLevelSolverBlasLapackShellBackend::accumulateSubdomainCorrection(const std::size_t subdomain_num, Vec y)
{
    auto& subdomain_data = d_data->subdomains[subdomain_num];
    PetscScalar* y_values = nullptr;
    int ierr = VecGetArray(y, &y_values);
    IBTK_CHKERRQ(ierr);
    for (std::size_t update_num = 0; update_num < subdomain_data.update_dofs->size(); ++update_num)
    {
        y_values[(*subdomain_data.update_dofs)[update_num]] +=
            subdomain_data.rhs_workspace[static_cast<std::size_t>(subdomain_data.update_local_positions[update_num])];
    }
    ierr = VecRestoreArray(y, &y_values);
    IBTK_CHKERRQ(ierr);
}

void
PETScLevelSolverBlasLapackShellBackend::updateSubdomainResidual(const std::size_t subdomain_num, Vec /*x*/, Vec /*y*/)
{
    auto& data = *d_data;
    auto& subdomain_data = data.subdomains[subdomain_num];
    int ierr = VecZeroEntries(data.patch_correction);
    IBTK_CHKERRQ(ierr);
    PetscScalar* correction_values = nullptr;
    ierr = VecGetArray(data.patch_correction, &correction_values);
    IBTK_CHKERRQ(ierr);
    for (std::size_t update_num = 0; update_num < subdomain_data.update_dofs->size(); ++update_num)
    {
        correction_values[(*subdomain_data.update_dofs)[update_num]] =
            subdomain_data.rhs_workspace[static_cast<std::size_t>(subdomain_data.update_local_positions[update_num])];
    }
    ierr = VecRestoreArray(data.patch_correction, &correction_values);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(d_solver_state.petsc_mat, data.patch_correction, data.residual_update);
    IBTK_CHKERRQ(ierr);
    ierr = VecAXPY(data.residual, -1.0, data.residual_update);
    IBTK_CHKERRQ(ierr);
}
} // namespace IBTK
