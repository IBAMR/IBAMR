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
#include <ibtk/PETScLevelSolverSubdomainSolver.h>
#include <ibtk/string_utilities.h>

#include <tbox/Database.h>
#include <tbox/Utilities.h>

#include <petscblaslapack.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

#if defined(PETSC_USE_COMPLEX)
#error "The BLAS/LAPACK subdomain solver requires real PETSc scalars."
#endif

namespace IBTK
{
namespace
{
/*! \brief Real-scalar BLAS/LAPACK solves for PETScLevelSolverSubdomainSolver.
 *
 * The settings are those of make_blas_lapack_subdomain_solver(). Each subdomain retains a
 * factorization or inverse/pseudoinverse.
 */
class PETScLevelSolverBlasLapackSubdomainSolver
{
public:
    /*! \brief Read subdomain solver settings, using defaults for a null database. */
    explicit PETScLevelSolverBlasLapackSubdomainSolver(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);
    void initializeSolverState(const std::vector<Mat>& matrices,
                               const std::vector<IS>& subdomains,
                               const std::string& options_prefix);
    void deallocateSolverState();
    void solve(std::size_t first, std::size_t last, Vec b, Vec x);

private:
    enum class SubdomainSolverType
    {
        SVD,
        LU,
        SYMMETRIC_INDEFINITE,
        QR
    };
    struct SubdomainData
    {
        PetscBLASInt local_size = 0;
        std::vector<PetscScalar> solve_data;
        std::vector<PetscBLASInt> pivots;
    };

    /*! \brief Factor the subdomain matrix or construct its solve matrix. */
    void initializeSubdomainSolver(SubdomainData& subdomain_data, std::size_t subdomain_num);
    /*! \brief Construct an inverse from the full-rank QR factorization. */
    void initializeQRSolver(SubdomainData& subdomain_data, std::size_t subdomain_num);
    /*! \brief Construct the SVD pseudoinverse using the configured rank tolerance. */
    void initializeSVDSolver(SubdomainData& subdomain_data, std::size_t subdomain_num);
    /*! \brief Require symmetry for the symmetric-indefinite factorization. */
    void verifySymmetricSubdomainMatrix(const SubdomainData& subdomain_data, std::size_t subdomain_num) const;
    /*! \brief Solve for one subdomain, reading the right-hand side and writing the solution in place. */
    void solveSubdomainSystem(const SubdomainData& subdomain_data,
                              std::size_t subdomain_num,
                              const PetscScalar* rhs,
                              PetscScalar* solution) const;

    std::vector<SubdomainData> d_subdomains;
    std::vector<PetscInt> d_offsets;
    std::string d_options_prefix;
    SubdomainSolverType d_subdomain_solver_type = SubdomainSolverType::SVD;
    PetscReal d_subdomain_solver_rcond = -1.0;
};

std::string
error_header(const char* const method, const std::string& options_prefix)
{
    return std::string("PETScLevelSolverBlasLapackSubdomainSolver::") + method + "():\n  options prefix \"" +
           options_prefix + "\": ";
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

PETScLevelSolverBlasLapackSubdomainSolver::PETScLevelSolverBlasLapackSubdomainSolver(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
{
    if (!input_db)
    {
        return;
    }
    const std::string solver_type = input_db->getStringWithDefault("blas_lapack_subdomain_solver_type", "svd");
    if (equals_ignore_case(solver_type, "svd"))
    {
        d_subdomain_solver_type = SubdomainSolverType::SVD;
    }
    else if (equals_ignore_case(solver_type, "lu"))
    {
        d_subdomain_solver_type = SubdomainSolverType::LU;
    }
    else if (equals_ignore_case(solver_type, "symmetric-indefinite"))
    {
        d_subdomain_solver_type = SubdomainSolverType::SYMMETRIC_INDEFINITE;
    }
    else if (equals_ignore_case(solver_type, "qr"))
    {
        d_subdomain_solver_type = SubdomainSolverType::QR;
    }
    else
    {
        TBOX_ERROR("PETScLevelSolverBlasLapackSubdomainSolver::PETScLevelSolverBlasLapackSubdomainSolver():\n"
                   << "  unsupported blas_lapack_subdomain_solver_type = " << solver_type
                   << "; supported values are \"svd\", \"lu\", \"symmetric-indefinite\", and \"qr\".\n");
    }
    d_subdomain_solver_rcond = input_db->getDoubleWithDefault("blas_lapack_subdomain_solver_rcond", -1.0);
    if (!std::isfinite(d_subdomain_solver_rcond))
    {
        TBOX_ERROR("PETScLevelSolverBlasLapackSubdomainSolver::PETScLevelSolverBlasLapackSubdomainSolver():\n"
                   << "  blas_lapack_subdomain_solver_rcond must be finite.\n");
    }
}

void
PETScLevelSolverBlasLapackSubdomainSolver::initializeSolverState(const std::vector<Mat>& matrices,
                                                                 const std::vector<IS>& /*subdomains*/,
                                                                 const std::string& options_prefix)
{
    deallocateSolverState();
    d_options_prefix = options_prefix;
    d_subdomains.resize(matrices.size());
    d_offsets.assign(matrices.size() + 1, 0);
    for (std::size_t i = 0; i < matrices.size(); ++i)
    {
        PetscInt rows = 0, columns = 0;
        int ierr = MatGetSize(matrices[i], &rows, &columns);
        IBTK_CHKERRQ(ierr);
        if (rows != columns)
        {
            TBOX_ERROR(error_header("initializeSolverState", d_options_prefix)
                       << "BLAS/LAPACK subdomain solver requires square subdomain matrices.\n");
        }
        SubdomainData& subdomain_data = d_subdomains[i];
        ierr = PetscBLASIntCast(rows, &subdomain_data.local_size);
        IBTK_CHKERRQ(ierr);
        subdomain_data.solve_data.resize(static_cast<std::size_t>(rows) * static_cast<std::size_t>(rows));
        d_offsets[i + 1] = d_offsets[i] + rows;
        // Copy the stored entries of the sparse matrix into the zeroed dense matrix.
        std::fill(subdomain_data.solve_data.begin(), subdomain_data.solve_data.end(), PetscScalar(0.0));
        for (PetscInt row = 0; row < rows; ++row)
        {
            PetscInt count = 0;
            const PetscInt* column_indices = nullptr;
            const PetscScalar* values = nullptr;
            ierr = MatGetRow(matrices[i], row, &count, &column_indices, &values);
            IBTK_CHKERRQ(ierr);
            for (PetscInt k = 0; k < count; ++k)
            {
                subdomain_data.solve_data[column_major_index(static_cast<PetscBLASInt>(row),
                                                             static_cast<PetscBLASInt>(column_indices[k]),
                                                             subdomain_data.local_size)] = values[k];
            }
            ierr = MatRestoreRow(matrices[i], row, &count, &column_indices, &values);
            IBTK_CHKERRQ(ierr);
        }
        initializeSubdomainSolver(subdomain_data, i);
    }
}

void
PETScLevelSolverBlasLapackSubdomainSolver::deallocateSolverState()
{
    d_subdomains.clear();
    d_offsets.clear();
}

void
PETScLevelSolverBlasLapackSubdomainSolver::solve(const std::size_t first, const std::size_t last, Vec b, Vec x)
{
    const PetscScalar* rhs = nullptr;
    PetscScalar* solution = nullptr;
    int ierr = VecGetArrayRead(b, &rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecGetArray(x, &solution);
    IBTK_CHKERRQ(ierr);
    for (std::size_t i = first; i < last; ++i)
    {
        solveSubdomainSystem(d_subdomains[i], i, rhs + d_offsets[i], solution + d_offsets[i]);
    }
    ierr = VecRestoreArray(x, &solution);
    IBTK_CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(b, &rhs);
    IBTK_CHKERRQ(ierr);
}

void
PETScLevelSolverBlasLapackSubdomainSolver::initializeSubdomainSolver(SubdomainData& subdomain_data,
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
            TBOX_ERROR(error_header("initializeSubdomainSolver", d_options_prefix)
                       << "BLAS/LAPACK LU factorization failed for subdomain " << subdomain_num
                       << " with info = " << info << ".\n");
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
            TBOX_ERROR(error_header("initializeSubdomainSolver", d_options_prefix)
                       << "BLAS/LAPACK symmetric-indefinite workspace query failed for "
                       << "subdomain " << subdomain_num << " with info = " << info << ".\n");
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
            TBOX_ERROR(error_header("initializeSubdomainSolver", d_options_prefix)
                       << "BLAS/LAPACK symmetric-indefinite factorization failed for "
                       << "subdomain " << subdomain_num << " with info = " << info << ".\n");
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
PETScLevelSolverBlasLapackSubdomainSolver::initializeQRSolver(SubdomainData& subdomain_data,
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
        TBOX_ERROR(error_header("initializeQRSolver", d_options_prefix)
                   << "BLAS/LAPACK QR workspace query failed for subdomain " << subdomain_num << " with info = " << info
                   << ".\n");
    }
    lwork = query_lapack_workspace_size(work_query);
    std::vector<PetscScalar> work(static_cast<std::size_t>(lwork));
    LAPACKgeqrf_(&n, &n, qr_factor.data(), &n, tau.data(), work.data(), &lwork, &info);
    if (info != 0)
    {
        TBOX_ERROR(error_header("initializeQRSolver", d_options_prefix)
                   << "BLAS/LAPACK QR factorization failed for subdomain " << subdomain_num << " with info = " << info
                   << ".\n");
    }

    PetscReal max_abs_diagonal = 0.0;
    for (PetscBLASInt k = 0; k < n; ++k)
    {
        max_abs_diagonal = std::max(max_abs_diagonal, PetscAbsScalar(qr_factor[column_major_index(k, k, n)]));
    }
    if (max_abs_diagonal == 0.0)
    {
        TBOX_ERROR(error_header("initializeQRSolver", d_options_prefix)
                   << "BLAS/LAPACK QR factorization detected a zero-rank subdomain " << subdomain_num << ".\n");
    }
    // Unpivoted QR is accepted only as a full-rank solver. Rounding can turn an exactly singular matrix
    // into an R with tiny but nonzero diagonal entries, so the computed diagonal is always compared
    // with a threshold: the requested rcond, or n times machine epsilon when rcond is negative.
    // SVD is the rank-deficient mode.
    const PetscReal rcond = d_subdomain_solver_rcond >= 0.0 ?
                                d_subdomain_solver_rcond :
                                static_cast<PetscReal>(n) * std::numeric_limits<PetscReal>::epsilon();
    for (PetscBLASInt k = 0; k < n; ++k)
    {
        if (PetscAbsScalar(qr_factor[column_major_index(k, k, n)]) <= rcond * max_abs_diagonal)
        {
            TBOX_ERROR(error_header("initializeQRSolver", d_options_prefix)
                       << "BLAS/LAPACK QR factorization detected a rank-deficient subdomain " << subdomain_num
                       << " with a relative threshold of " << rcond << ".\n");
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
        TBOX_ERROR(error_header("initializeQRSolver", d_options_prefix)
                   << "BLAS/LAPACK QR-application workspace query failed for subdomain " << subdomain_num
                   << " with info = " << info << ".\n");
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
        TBOX_ERROR(error_header("initializeQRSolver", d_options_prefix)
                   << "BLAS/LAPACK QR application failed for subdomain " << subdomain_num << " with info = " << info
                   << ".\n");
    }
    const char uplo = 'U';
    const char no_transpose = 'N';
    const char nonunit = 'N';
    LAPACKtrtrs_(
        &uplo, &no_transpose, &nonunit, &n, &n, qr_factor.data(), &n, subdomain_data.solve_data.data(), &n, &info);
    if (info != 0)
    {
        TBOX_ERROR(error_header("initializeQRSolver", d_options_prefix)
                   << "BLAS/LAPACK triangular QR solve failed for subdomain " << subdomain_num
                   << " with info = " << info << ".\n");
    }
}

void
PETScLevelSolverBlasLapackSubdomainSolver::initializeSVDSolver(SubdomainData& subdomain_data,
                                                               const std::size_t subdomain_num)
{
    const PetscBLASInt n = subdomain_data.local_size;
    // Use GELSS once during setup to compute the minimum-norm least-squares solutions with
    // right-hand side I. The overwritten right-hand side stores the pseudoinverse with the
    // configured singular-value cutoff, and becomes the only persistent matrix-sized solve object.
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
        TBOX_ERROR(error_header("initializeSVDSolver", d_options_prefix)
                   << "BLAS/LAPACK SVD workspace query failed for subdomain " << subdomain_num
                   << " with info = " << info << ".\n");
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
        TBOX_ERROR(error_header("initializeSVDSolver", d_options_prefix)
                   << "BLAS/LAPACK SVD solve-matrix construction failed for subdomain " << subdomain_num
                   << " with info = " << info << " and effective rank = " << effective_rank << ".\n");
    }
}

void
PETScLevelSolverBlasLapackSubdomainSolver::verifySymmetricSubdomainMatrix(const SubdomainData& subdomain_data,
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
    // The tolerance scales with the matrix, so that multiplying it by a constant does not change the
    // outcome. A zero matrix is symmetric; the factorization reports that it is singular.
    if (max_entry == 0.0)
    {
        return;
    }
    const PetscReal tolerance = 100.0 * std::numeric_limits<PetscReal>::epsilon() * max_entry;
    if (max_asymmetry > tolerance)
    {
        TBOX_ERROR(error_header("verifySymmetricSubdomainMatrix", d_options_prefix)
                   << "subdomain " << subdomain_num << " is not symmetric enough for the "
                   << "symmetric-indefinite solver; maximum asymmetry is " << max_asymmetry / max_entry
                   << " times the largest entry, above the limit of " << tolerance / max_entry << ".\n");
    }
}

void
PETScLevelSolverBlasLapackSubdomainSolver::solveSubdomainSystem(const SubdomainData& subdomain_data,
                                                                const std::size_t subdomain_num,
                                                                const PetscScalar* rhs,
                                                                PetscScalar* solution) const
{
    if (subdomain_data.local_size == 0)
    {
        return;
    }
    PetscBLASInt n = subdomain_data.local_size;
    PetscBLASInt nrhs = 1;
    PetscBLASInt info = 0;
    // Some PETSc releases declare the LAPACK solve arguments non-const although they are only read.
    auto* factors = const_cast<PetscScalar*>(subdomain_data.solve_data.data());
    auto* pivots = const_cast<PetscBLASInt*>(subdomain_data.pivots.data());
    // The factorization solvers overwrite their right-hand side, so they solve in the output.
    if (d_subdomain_solver_type == SubdomainSolverType::LU ||
        d_subdomain_solver_type == SubdomainSolverType::SYMMETRIC_INDEFINITE)
    {
        std::copy(rhs, rhs + subdomain_data.local_size, solution);
    }
    switch (d_subdomain_solver_type)
    {
    case SubdomainSolverType::LU:
    {
        const char transpose = 'N';
        LAPACKgetrs_(&transpose, &n, &nrhs, factors, &n, pivots, solution, &n, &info);
        if (info != 0)
        {
            TBOX_ERROR(error_header("solveSubdomainSystem", d_options_prefix)
                       << "BLAS/LAPACK LU solve failed for subdomain " << subdomain_num << " with info = " << info
                       << ".\n");
        }
        return;
    }
    case SubdomainSolverType::SYMMETRIC_INDEFINITE:
    {
        const char uplo = 'L';
        LAPACKsytrs_(&uplo, &n, &nrhs, factors, &n, pivots, solution, &n, &info);
        if (info != 0)
        {
            TBOX_ERROR(error_header("solveSubdomainSystem", d_options_prefix)
                       << "BLAS/LAPACK symmetric-indefinite solve failed for subdomain " << subdomain_num
                       << " with info = " << info << ".\n");
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
                  rhs,
                  &increment,
                  &beta,
                  solution,
                  &increment);
        return;
    }
    }
}

} // namespace

PETScLevelSolverSubdomainSolver
make_blas_lapack_subdomain_solver(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
{
    return PETScLevelSolverSubdomainSolver(std::in_place_type<PETScLevelSolverBlasLapackSubdomainSolver>, input_db);
}
} // namespace IBTK
