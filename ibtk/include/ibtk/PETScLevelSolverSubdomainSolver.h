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

#ifndef included_IBTK_PETScLevelSolverSubdomainSolver
#define included_IBTK_PETScLevelSolverSubdomainSolver

#include <ibtk/config.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <petscksp.h>

#include <concepts>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace IBTK
{
/*!
 * \brief Requirements of an implementation of the subdomain solver of a PETScLevelSolver
 * shell preconditioner.
 *
 * A type T models this concept when a mutable T provides the operations of
 * PETScLevelSolverSubdomainSolver, with no result. T need not derive from an IBAMR class
 * and need not be copyable, movable, or default constructible. The concept checks only that
 * the operations can be called; the requirements on what they do are those documented for
 * PETScLevelSolverSubdomainSolver.
 */
template <class T>
concept PETScLevelSolverSubdomainSolverImplementation = requires(T & solver,
                                                                 const std::vector<Mat>& matrices,
                                                                 const std::vector<IS>& subdomains,
                                                                 const std::string& options_prefix,
                                                                 std::size_t first,
                                                                 std::size_t last,
                                                                 Vec b,
                                                                 Vec x)
{
    {
        solver.initializeSolverState(matrices, subdomains, options_prefix)
    } -> std::same_as<void>;
    {
        solver.deallocateSolverState()
    } -> std::same_as<void>;
    {
        solver.solve(first, last, b, x)
    } -> std::same_as<void>;
};

/*!
 * \brief Solver for the subdomain problems of a PETScLevelSolver shell preconditioner.
 *
 * With A the level operator, PETScLevelSolver divides the level into overlapping
 * subdomains and forms the sequential matrix A_i = A(O_i, O_i) for each
 * subdomain i of this rank, in which O_i is the sorted index set of the
 * subdomain. PETScLevelSolver composes the subdomain corrections, additively or
 * multiplicatively, and performs all communication. A subdomain solver only applies
 * (an approximation to) the inverse of A_i to a sequential vector.
 *
 * This class is a move-only handle that owns an implementation, which is any type that
 * satisfies PETScLevelSolverSubdomainSolverImplementation. It is constructed in place, so
 * an implementation need not be copyable or movable. For example,
 * \code
 * class MySolver
 * {
 * public:
 *     explicit MySolver(double tolerance);
 *     void initializeSolverState(const std::vector<Mat>& matrices,
 *                                const std::vector<IS>& subdomains,
 *                                const std::string& options_prefix);
 *     void deallocateSolverState();
 *     void solve(std::size_t first, std::size_t last, Vec b, Vec x);
 * };
 *
 * level_solver.setSubdomainSolver(
 *     IBTK::PETScLevelSolverSubdomainSolver(std::in_place_type<MySolver>, 1.0e-10));
 * \endcode
 *
 * initializeSolverState() receives the matrices of this rank's subdomains and their
 * index sets, which are borrowed during initialization; implementations retain any
 * objects they need afterwards. deallocateSolverState() releases owned state, is called
 * only after initializeSolverState(), and may be called repeatedly. solve() may be
 * called only while solver state is initialized. A rank may have no subdomains.
 * PETScLevelSolver calls initializeSolverState() and deallocateSolverState() on every
 * rank, so they may communicate, and every rank must take part in that communication.
 * solve() is not collective and must not communicate.
 *
 * The right-hand sides and solutions of the subdomains are packed, in order, into
 * two sequential vectors, so that solve() needs no copies to gather or scatter
 * them and can solve many subdomains in one call: the entries of subdomain i
 * follow those of subdomains 0, ..., i - 1, and their number is the order of A_i.
 *
 * Destroying a handle destroys its implementation without calling
 * deallocateSolverState(), which can communicate and can depend on the state of the
 * owning PETScLevelSolver. The implementation must release its own resources when it is
 * destroyed.
 *
 * A handle that has been moved from is empty: it converts to false and must not be used.
 */
class PETScLevelSolverSubdomainSolver
{
public:
    /*!
     * \brief Construct an implementation of type Implementation in place from args.
     */
    template <class Implementation, class... Args>
    requires PETScLevelSolverSubdomainSolverImplementation<Implementation>&&
        std::constructible_from<Implementation, Args...> explicit PETScLevelSolverSubdomainSolver(
            std::in_place_type_t<Implementation>,
            Args&&... args);

    /*!
     * \brief Copying is not supported.
     */
    PETScLevelSolverSubdomainSolver(const PETScLevelSolverSubdomainSolver&) = delete;

    /*!
     * \brief Copying is not supported.
     */
    PETScLevelSolverSubdomainSolver& operator=(const PETScLevelSolverSubdomainSolver&) = delete;

    /*!
     * \brief Transfer ownership of the implementation of other, which becomes empty.
     */
    PETScLevelSolverSubdomainSolver(PETScLevelSolverSubdomainSolver&& other) noexcept;

    /*!
     * \brief Replace the implementation with that of other, which becomes empty.
     */
    PETScLevelSolverSubdomainSolver& operator=(PETScLevelSolverSubdomainSolver&& other) noexcept;

    /*!
     * \brief Destroy the implementation.
     */
    ~PETScLevelSolverSubdomainSolver();

    /*!
     * \brief Return whether this handle owns an implementation.
     */
    explicit operator bool() const;

    /*!
     * \brief Initialize the solvers of the subdomain matrices.
     *
     * The rows and columns of matrices[i] are those of the sorted index set subdomains[i], which
     * contains global DOF indices; the two vectors have one entry for each subdomain of this rank.
     * The options prefix is that of the level solver.
     */
    void initializeSolverState(const std::vector<Mat>& matrices,
                               const std::vector<IS>& subdomains,
                               const std::string& options_prefix);

    /*!
     * \brief Release the initialized state.
     */
    void deallocateSolverState();

    /*!
     * \brief Solve A_i x_i = b_i for the subdomains i = first, ..., last - 1.
     *
     * The vectors b and x are distinct sequential vectors in host memory that hold the
     * packed right-hand sides and solutions of all subdomains, so 0 <= first <= last <= the
     * number of subdomains, and first == last is a valid call that does nothing. This
     * method overwrites the entries of x of the subdomains that it solves, whatever their
     * initial values, and does not modify b or the other entries of x. The entries of b of
     * subdomains outside first, ..., last - 1 are unspecified and must not be used. The
     * subdomains of one call are independent of each other, so an implementation may solve
     * them in any order or concurrently.
     */
    void solve(std::size_t first, std::size_t last, Vec b, Vec x);

private:
    //! Type-erased operations of an implementation.
    class Adapter
    {
    public:
        virtual ~Adapter() = default;
        virtual void initializeSolverState(const std::vector<Mat>& matrices,
                                           const std::vector<IS>& subdomains,
                                           const std::string& options_prefix) = 0;
        virtual void deallocateSolverState() = 0;
        virtual void solve(std::size_t first, std::size_t last, Vec b, Vec x) = 0;
    };

    //! Adapter that holds an implementation.
    template <class Implementation>
    class ImplementationAdapter;

    std::unique_ptr<Adapter> d_adapter;
};

/*!
 * \brief Return the built-in subdomain solver, which uses PETSc.
 *
 * There is one KSP for each subdomain, with the level options prefix followed by
 * "_sub". It defaults to preonly with LU, and PETSc options may override this
 * configuration.
 */
PETScLevelSolverSubdomainSolver make_petsc_subdomain_solver();

/*!
 * \brief Return the subdomain solver that uses BLAS/LAPACK routines. Each subdomain
 * retains a factorization or inverse/pseudoinverse.
 *
 * This subdomain solver requires real PETSc scalars. The settings are read from
 * input_db, and defaults are used for a null database.
 * blas_lapack_subdomain_solver_type is "svd" (default), "lu",
 * "symmetric-indefinite", or "qr". LU, symmetric-indefinite and QR fail on a
 * singular subdomain matrix (QR: on one that is rank-deficient by its threshold),
 * whereas SVD never fails and forms a pseudoinverse.
 * blas_lapack_subdomain_solver_rcond must be finite, defaults to -1.0, and
 * affects only SVD and QR. Let epsilon be the machine epsilon of PetscReal and
 * n the subdomain size. A negative value selects LAPACK's SVD default cutoff
 * and, for QR, the threshold n * epsilon; a nonnegative value is the relative
 * SVD cutoff or the QR threshold. QR is a full-rank solver: it fails if any
 * |R_ii| is at most the threshold times max_j |R_jj|, with R the triangular
 * factor. The symmetric-indefinite solver requires the subdomain matrix A to
 * satisfy max_ij |A_ij - A_ji| <= 100 * epsilon * max_ij |A_ij|.
 */
PETScLevelSolverSubdomainSolver
make_blas_lapack_subdomain_solver(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

/*!
 * \brief Return the subdomain solver that factors each subdomain matrix with Eigen.
 *
 * The factorization is selected by eigen_subdomain_solver_type (default COL_PIV_HOUSEHOLDER_QR,
 * for the same reason make_blas_lapack_subdomain_solver()'s default is SVD: a singular subdomain gets a
 * bounded pseudo-solution from the truncated pivots instead of the Inf or NaN that PARTIAL_PIV_LU or
 * HOUSEHOLDER_QR would return, at the cost of allocating a solve matrix at setup instead of solving
 * in place; PARTIAL_PIV_LU is available for callers who already know their subdomains are nonsingular
 * and want to skip that cost) and eigen_subdomain_solver_threshold
 * (default -1). Available types are LLT, LDLT,
 * PARTIAL_PIV_LU, FULL_PIV_LU, HOUSEHOLDER_QR, COL_PIV_HOUSEHOLDER_QR,
 * COMPLETE_ORTHOGONAL_DECOMPOSITION, FULL_PIV_HOUSEHOLDER_QR, JACOBI_SVD and
 * BDC_SVD. Names are case-insensitive. LLT requires positive definiteness and LDLT requires symmetry.
 * Thresholds must be finite. A nonnegative value sets Eigen's relative rank threshold
 * where the type supports one, and a negative value retains its default. Eigen's LLT,
 * LDLT and PARTIAL_PIV_LU solve in place. The other types allocate a temporary vector
 * in each solve, so they form the solve matrix at setup instead, as
 * make_eigen_pseudoinverse_subdomain_solver() does, and an application is a matrix-vector
 * product that does not allocate.
 *
 * Only LLT and LDLT report a failed factorization. The other types do not detect a
 * singular subdomain matrix: PARTIAL_PIV_LU and HOUSEHOLDER_QR then return infinite or
 * NaN values. COL_PIV_HOUSEHOLDER_QR uses Eigen's own internal pivot cutoff in solve():
 * eigen_subdomain_solver_threshold does not control truncation for this type.
 * FULL_PIV_HOUSEHOLDER_QR, COMPLETE_ORTHOGONAL_DECOMPOSITION, JACOBI_SVD and BDC_SVD do
 * honor the configured threshold in solve(), treating pivots or singular values below it
 * as zero.
 */
PETScLevelSolverSubdomainSolver make_eigen_subdomain_solver(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

/*!
 * \brief Return the subdomain solver that forms the solve matrix of each subdomain with Eigen.
 *
 * The type is eigen_subdomain_pseudoinverse_type (default COL_PIV_HOUSEHOLDER_QR) and the
 * threshold is eigen_subdomain_pseudoinverse_threshold (default -1), with the values that
 * make_eigen_subdomain_solver() accepts. COMPLETE_ORTHOGONAL_DECOMPOSITION and the SVD types form Moore-Penrose
 * pseudoinverses. The other types solve against the identity and retain their pivot and rank policy.
 */
PETScLevelSolverSubdomainSolver
make_eigen_pseudoinverse_subdomain_solver(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);
} // namespace IBTK

#include <ibtk/private/PETScLevelSolverSubdomainSolver-inl.h>

#endif
