// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_FischerGuess
#define included_IBTK_FischerGuess

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#ifdef IBTK_HAVE_LIBMESH

#include <libmesh/numeric_vector.h>

IBTK_DISABLE_EXTRA_WARNINGS
#include <Eigen/Core>
IBTK_ENABLE_EXTRA_WARNINGS

namespace IBTK
{
/**
 * Class implementing a modified version of Fischer's first algorithm from the
 * 1998 manuscript "Projection techniques for iterative solution of A x = b with
 * successive right-hand sides".
 *
 * In essence, the caller submits pairs of solution and corresponding RHS
 * vectors from solving the same matrix equation several times (i.e., at
 * different time steps). This class uses those pairs to compute an estimated
 * solution vector for a new RHS. In practice, this estimate is a good enough
 * that it cuts the number of required solver iterations by about 50%, while
 * computing the initial guess is relatively inexpensive.
 *
 * The algorithm used below varies from the one Fischer proposed in a few major
 * ways:
 * <ol>
 *   <li>The RHS vectors are not orthogonalized. Instead, the least-squares
 *     equation is solved with the SVD, which enables the projection matrix to be
 *     singular (or nearly so).</li>
 *   <li>Since this class does not orthogonalize, it can compute the projection
 *     matrix incrementally. This means (unlike the PETSc implementation) we can
 *     always use the full set of available vectors to compute a guess.</li>
 *   <li>The projection is done with the standard $R^n$ inner product, which may
 *     not be optimal, but is much faster than doing sparse matrix-vector
 *     multiplications.</li>
 * </ol>
 *
 * Since the systems we solve typically require low iteration counts, these
 * tricks to get a less optimal guess at a lower cost are beneficial.
 */
class FischerGuess
{
public:
    /**
     * Constructor. Sets the number of stored vectors to five.
     */
    FischerGuess();

    /**
     * Constructor.
     *
     * @param n_vectors The number of stored vectors.
     */
    FischerGuess(int n_vectors);

    /**
     * Add a new solution and RHS pair to the stored collection. If the
     * collection is full then the oldest vectors are removed.
     */
    void submit(const libMesh::NumericVector<double>& solution, const libMesh::NumericVector<double>& rhs);

    /**
     * Given a RHS vector, use the stored collection of vectors to
     * compute an estimate of its corresponding solution vector.
     */
    void guess(libMesh::NumericVector<double>& solution, const libMesh::NumericVector<double>& rhs) const;

protected:
    int d_n_max_vectors = 5;

    int d_n_stored_vectors = 0;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> d_correlation_matrix;

    std::vector<std::unique_ptr<libMesh::NumericVector<double> > > d_solutions;
    std::vector<std::unique_ptr<libMesh::NumericVector<double> > > d_rhs;
};
} // namespace IBTK

#endif //#ifdef IBTK_HAVE_LIBMESH
#endif
