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
#include <ibtk/private/PETScLevelSolverBlasLapackShellBackend.h>

#include <tbox/Database.h>
#include <tbox/Utilities.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>

namespace IBTK
{
namespace
{
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

PETScLevelSolverBlasLapackShellBackend::PETScLevelSolverBlasLapackShellBackend(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
{
    if (!input_db)
    {
        return;
    }
    std::string solver_type = input_db->getStringWithDefault("blas_lapack_subdomain_solver_type", "svd");
    std::transform(solver_type.begin(),
                   solver_type.end(),
                   solver_type.begin(),
                   [](const unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (solver_type == "svd")
    {
        d_subdomain_solver_type = SubdomainSolverType::SVD;
    }
    else if (solver_type == "lu")
    {
        d_subdomain_solver_type = SubdomainSolverType::LU;
    }
    else if (solver_type == "symmetric-indefinite")
    {
        d_subdomain_solver_type = SubdomainSolverType::SYMMETRIC_INDEFINITE;
    }
    else if (solver_type == "qr")
    {
        d_subdomain_solver_type = SubdomainSolverType::QR;
    }
    else if (solver_type == "cholesky")
    {
        TBOX_ERROR("BLAS/LAPACK Cholesky is not supported for indefinite Stokes subdomain matrices; "
                   << "use 'symmetric-indefinite' when a symmetry-specific factorization is required.\n");
    }
    else
    {
        TBOX_ERROR("Unsupported blas_lapack_subdomain_solver_type = "
                   << solver_type << "; supported values are 'svd', 'lu', 'symmetric-indefinite', and 'qr'.\n");
    }
    d_subdomain_solver_rcond = input_db->getDoubleWithDefault("blas_lapack_subdomain_solver_rcond", -1.0);
    if (!std::isfinite(d_subdomain_solver_rcond))
    {
        TBOX_ERROR("blas_lapack_subdomain_solver_rcond must be finite.\n");
    }
}

PETScLevelSolverBlasLapackShellBackend::~PETScLevelSolverBlasLapackShellBackend()
{
    deallocateSolverState();
}

void
PETScLevelSolverBlasLapackShellBackend::initializeSolverState(Mat mat,
                                                              Vec x,
                                                              Vec b,
                                                              const std::vector<IS>& overlap,
                                                              const std::vector<IS>& nonoverlap,
                                                              const std::string& options_prefix,
                                                              const bool use_multiplicative,
                                                              const PETScLevelSolverShellTraversal traversal)
{
    deallocateSolverState();
    initializeComposition(mat, x, b, use_multiplicative, traversal);
    d_options_prefix = options_prefix;
    if (IBTK_MPI::getNodes() != 1)
    {
        TBOX_ERROR(d_options_prefix << " BLAS/LAPACK shell backend requires one MPI rank.\n");
    }
#if defined(PETSC_USE_COMPLEX)
    TBOX_ERROR("PETScLevelSolverBlasLapackShellBackend requires real PETSc scalars.\n");
#endif
    PetscInt rows = 0, columns = 0, local_rows = 0;
    int ierr = MatGetSize(mat, &rows, &columns);
    IBTK_CHKERRQ(ierr);
    ierr = MatGetLocalSize(mat, &local_rows, nullptr);
    IBTK_CHKERRQ(ierr);
    if (rows != columns || local_rows != rows)
    {
        TBOX_ERROR(d_options_prefix << " BLAS/LAPACK shell backend requires a serial square operator.\n");
    }
    TBOX_ASSERT(overlap.size() == nonoverlap.size());
    d_subdomains.resize(overlap.size());
    for (std::size_t i = 0; i < overlap.size(); ++i)
    {
        PetscInt n = 0, m = 0;
        ierr = ISGetLocalSize(overlap[i], &n);
        IBTK_CHKERRQ(ierr);
        ierr = ISGetLocalSize(nonoverlap[i], &m);
        IBTK_CHKERRQ(ierr);
        ierr = PetscBLASIntCast(n, &d_subdomains[i].local_size);
        IBTK_CHKERRQ(ierr);
        const PetscInt* indices = nullptr;
        ierr = ISGetIndices(overlap[i], &indices);
        IBTK_CHKERRQ(ierr);
        if (n > 0)
        {
            d_subdomains[i].overlap_dofs.assign(indices, indices + n);
        }
        ierr = ISRestoreIndices(overlap[i], &indices);
        IBTK_CHKERRQ(ierr);
        ierr = ISGetIndices(nonoverlap[i], &indices);
        IBTK_CHKERRQ(ierr);
        for (PetscInt j = 0; j < m; ++j)
        {
            const std::vector<PetscInt>::const_iterator position = std::lower_bound(
                d_subdomains[i].overlap_dofs.cbegin(), d_subdomains[i].overlap_dofs.cend(), indices[j]);
            TBOX_ASSERT(position != d_subdomains[i].overlap_dofs.cend() && *position == indices[j]);
            d_subdomains[i].update_local_positions.push_back(
                static_cast<PetscBLASInt>(position - d_subdomains[i].overlap_dofs.cbegin()));
        }
        ierr = ISRestoreIndices(nonoverlap[i], &indices);
        IBTK_CHKERRQ(ierr);
        if (use_multiplicative)
        {
            d_subdomains[i].update_local_positions.resize(n);
            std::iota(d_subdomains[i].update_local_positions.begin(), d_subdomains[i].update_local_positions.end(), 0);
        }
        d_subdomains[i].solve_data.resize(static_cast<std::size_t>(n) * static_cast<std::size_t>(n));
        d_subdomains[i].rhs_workspace.resize(static_cast<std::size_t>(n));
        // Reuse the RHS buffer for one extracted row during setup.
        for (PetscBLASInt row = 0; row < n; ++row)
        {
            ierr = MatGetValues(mat,
                                1,
                                &d_subdomains[i].overlap_dofs[row],
                                n,
                                d_subdomains[i].overlap_dofs.data(),
                                d_subdomains[i].rhs_workspace.data());
            IBTK_CHKERRQ(ierr);
            for (PetscBLASInt col = 0; col < n; ++col)
            {
                d_subdomains[i].solve_data[column_major_index(row, col, d_subdomains[i].local_size)] =
                    d_subdomains[i].rhs_workspace[col];
            }
        }
        initializeSubdomainSolver(d_subdomains[i], i);
    }
    finalizeComposition();
}

void
PETScLevelSolverBlasLapackShellBackend::deallocateSolverState()
{
    deallocateComposition();
    d_subdomains.clear();
}

std::size_t
PETScLevelSolverBlasLapackShellBackend::getNumberOfSubdomains() const
{
    return d_subdomains.size();
}

void
PETScLevelSolverBlasLapackShellBackend::beginSubdomainRhs(const std::size_t i, Vec source)
{
    const PetscScalar* rhs = nullptr;
    int ierr = VecGetArrayRead(source, &rhs);
    IBTK_CHKERRQ(ierr);
    for (PetscBLASInt j = 0; j < d_subdomains[i].local_size; ++j)
    {
        d_subdomains[i].rhs_workspace[j] = rhs[d_subdomains[i].overlap_dofs[j]];
    }
    ierr = VecRestoreArrayRead(source, &rhs);
    IBTK_CHKERRQ(ierr);
}

void
PETScLevelSolverBlasLapackShellBackend::endSubdomainRhs(std::size_t /*i*/, Vec /*source*/)
{
}

void
PETScLevelSolverBlasLapackShellBackend::solveSubdomain(const std::size_t i)
{
    solveSubdomainSystem(d_subdomains[i], i);
}

void
PETScLevelSolverBlasLapackShellBackend::accumulateSubdomainCorrection(const std::size_t i, Vec y)
{
    PetscScalar* correction = nullptr;
    int ierr = VecGetArray(y, &correction);
    IBTK_CHKERRQ(ierr);
    for (const PetscBLASInt position : d_subdomains[i].update_local_positions)
    {
        correction[d_subdomains[i].overlap_dofs[position]] += d_subdomains[i].rhs_workspace[position];
    }
    ierr = VecRestoreArray(y, &correction);
    IBTK_CHKERRQ(ierr);
}

const std::vector<PetscInt>&
PETScLevelSolverBlasLapackShellBackend::getSubdomainCorrectionDofs(const std::size_t i) const
{
    return d_subdomains[i].overlap_dofs;
}

void
PETScLevelSolverBlasLapackShellBackend::copySubdomainCorrection(const std::size_t i, PetscScalar* values)
{
    std::copy(d_subdomains[i].rhs_workspace.begin(), d_subdomains[i].rhs_workspace.end(), values);
}

void
PETScLevelSolverBlasLapackShellBackend::initializeSubdomainSolver(SubdomainData& subdomain_data,
                                                                  const std::size_t subdomain_num)
{
    if (subdomain_data.local_size == 0)
    {
        return;
    }
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
            TBOX_ERROR(d_options_prefix << " BLAS/LAPACK LU factorization failed for subdomain " << subdomain_num
                                        << " with info = " << info << ".");
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
            TBOX_ERROR(d_options_prefix << " BLAS/LAPACK symmetric-indefinite workspace query failed for "
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
            TBOX_ERROR(d_options_prefix << " BLAS/LAPACK symmetric-indefinite factorization failed for "
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
    // Older PETSc LAPACK declarations require mutable dimension pointers.
    PetscBLASInt n = subdomain_data.local_size;
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
        TBOX_ERROR(d_options_prefix << " BLAS/LAPACK QR workspace query failed for subdomain " << subdomain_num
                                    << " with info = " << info << ".");
    }
    lwork = query_lapack_workspace_size(work_query);
    std::vector<PetscScalar> work(static_cast<std::size_t>(lwork));
    LAPACKgeqrf_(&n, &n, qr_factor.data(), &n, tau.data(), work.data(), &lwork, &info);
    if (info != 0)
    {
        TBOX_ERROR(d_options_prefix << " BLAS/LAPACK QR factorization failed for subdomain " << subdomain_num
                                    << " with info = " << info << ".");
    }

    PetscReal max_abs_diagonal = 0.0;
    for (PetscBLASInt k = 0; k < n; ++k)
    {
        max_abs_diagonal = std::max(max_abs_diagonal, PetscAbsScalar(qr_factor[column_major_index(k, k, n)]));
    }
    if (max_abs_diagonal == 0.0)
    {
        TBOX_ERROR(d_options_prefix << " BLAS/LAPACK QR factorization detected a zero-rank subdomain " << subdomain_num
                                    << ".");
    }
    if (d_subdomain_solver_rcond >= 0.0)
    {
        // Unpivoted QR is accepted only as a full-rank solver. This threshold
        // checks the computed R diagonal; SVD is the rank-deficient mode.
        for (PetscBLASInt k = 0; k < n; ++k)
        {
            if (PetscAbsScalar(qr_factor[column_major_index(k, k, n)]) <= d_subdomain_solver_rcond * max_abs_diagonal)
            {
                TBOX_ERROR(d_options_prefix
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
        TBOX_ERROR(d_options_prefix << " BLAS/LAPACK QR-application workspace query failed for subdomain "
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
        TBOX_ERROR(d_options_prefix << " BLAS/LAPACK QR application failed for subdomain " << subdomain_num
                                    << " with info = " << info << ".");
    }
    const char uplo = 'U';
    const char no_transpose = 'N';
    const char nonunit = 'N';
    LAPACKtrtrs_(
        &uplo, &no_transpose, &nonunit, &n, &n, qr_factor.data(), &n, subdomain_data.solve_data.data(), &n, &info);
    if (info != 0)
    {
        TBOX_ERROR(d_options_prefix << " BLAS/LAPACK triangular QR solve failed for subdomain " << subdomain_num
                                    << " with info = " << info << ".");
    }
    subdomain_data.solution_workspace.resize(static_cast<std::size_t>(n));
}

void
PETScLevelSolverBlasLapackShellBackend::initializeSVDSolver(SubdomainData& subdomain_data,
                                                            const std::size_t subdomain_num)
{
#if defined(PETSC_USE_COMPLEX)
    TBOX_ERROR(d_options_prefix << " BLAS/LAPACK SVD setup requires real PETSc scalars.");
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
        TBOX_ERROR(d_options_prefix << " BLAS/LAPACK SVD workspace query failed for subdomain " << subdomain_num
                                    << " with info = " << info << ".");
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
        TBOX_ERROR(d_options_prefix << " BLAS/LAPACK SVD solve-matrix construction failed for subdomain "
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
        TBOX_ERROR(d_options_prefix << " subdomain " << subdomain_num << " is not symmetric enough for the "
                                    << "symmetric-indefinite solver; maximum asymmetry = " << max_asymmetry << ".");
    }
}

void
PETScLevelSolverBlasLapackShellBackend::solveSubdomainSystem(SubdomainData& subdomain_data,
                                                             const std::size_t subdomain_num) const
{
    if (subdomain_data.local_size == 0)
    {
        return;
    }
    PetscBLASInt nrhs = 1;
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
            TBOX_ERROR(d_options_prefix << " BLAS/LAPACK LU solve failed for subdomain " << subdomain_num
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
            TBOX_ERROR(d_options_prefix << " BLAS/LAPACK symmetric-indefinite solve failed for subdomain "
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

} // namespace IBTK
