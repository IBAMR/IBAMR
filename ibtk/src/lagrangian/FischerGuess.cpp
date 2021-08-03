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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/FischerGuess.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/TimerManager.h>

#include <Eigen/Dense>

#include <ibtk/app_namespaces.h>

namespace IBTK
{
namespace
{
static Timer* t_submit;
static Timer* t_guess;
} // namespace

FischerGuess::FischerGuess()
{
    using namespace SAMRAI::tbox;
    IBTK_DO_ONCE(t_submit = TimerManager::getManager()->getTimer("IBTK::FischerGuess::submit()");
                 t_guess = TimerManager::getManager()->getTimer("IBTK::FischerGuess::guess()"););
}

FischerGuess::FischerGuess(const int n_vectors) : FischerGuess()
{
    d_n_max_vectors = n_vectors;
}

void
FischerGuess::submit(const libMesh::NumericVector<double>& solution, const libMesh::NumericVector<double>& rhs)
{
    if (d_n_max_vectors == 0) return;
    IBTK_TIMER_START(t_submit);
    // update our list of vectors:
    if (d_n_stored_vectors == d_n_max_vectors)
    {
        // oldest entry is first
        d_solutions.erase(d_solutions.begin());
        d_solutions.push_back(solution.clone());
        d_rhs.erase(d_rhs.begin());
        d_rhs.push_back(rhs.clone());

        // shift the computed dot products up and to the left:
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat_copy(d_correlation_matrix);
        for (int i = 1; i < d_n_max_vectors; ++i)
        {
            for (int j = 1; j < d_n_max_vectors; ++j)
            {
                d_correlation_matrix(i - 1, j - 1) = mat_copy(i, j);
            }
        }
    }
    else
    {
        ++d_n_stored_vectors;
        d_solutions.push_back(solution.clone());
        d_rhs.push_back(rhs.clone());

        // Save the prior dot products in a different way:
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat_copy;

        mat_copy = d_correlation_matrix;
        d_correlation_matrix.resize(d_n_stored_vectors, d_n_stored_vectors);
        for (int i = 0; i < d_n_stored_vectors - 1; ++i)
        {
            for (int j = 0; j < d_n_stored_vectors - 1; ++j)
            {
                d_correlation_matrix(i, j) = mat_copy(i, j);
            }
        }
    }

    // Compute the last row and then copy it into the last column.
    for (int j = 0; j < d_n_stored_vectors; ++j)
    {
        const double inner = d_rhs.back()->dot(*d_rhs[j]);
        d_correlation_matrix(d_n_stored_vectors - 1, j) = inner;
        d_correlation_matrix(j, d_n_stored_vectors - 1) = inner;
    }
    IBTK_TIMER_STOP(t_submit);
}

void
FischerGuess::guess(libMesh::NumericVector<double>& solution, const libMesh::NumericVector<double>& rhs) const
{
    if (d_n_stored_vectors == 0)
    {
        return;
    }

    IBTK_TIMER_START(t_guess);
    Eigen::VectorXd coef_rhs(d_n_stored_vectors);
    for (int i = 0; i < d_n_stored_vectors; ++i)
    {
        coef_rhs(i) = d_rhs[i]->dot(rhs);
    }
    // TODO - if we were especially clever we could use these dot
    // products to update d_correlation_matrix for next time.

    // TODO - at some point in the future (if and when we require Eigen 3.3 or
    // newer) we can use the slightly better bcdSvd algorithm here. This isn't
    // worth worrying about yet since the SVD is O(n^3) and typically n = 5.
    Eigen::VectorXd coefs(d_n_stored_vectors);
    // Should the SVD fail for any reason just use the last solution as a guess.
    try
    {
        // SVD doesn't make sense with invalid input - throw something and
        // immediately catch it
        if (!d_correlation_matrix.allFinite()) throw int();
        if (!coef_rhs.allFinite()) throw int();
        Eigen::JacobiSVD<decltype(d_correlation_matrix)> svd(d_correlation_matrix, ComputeThinU | ComputeThinV);
        coefs = svd.solve(coef_rhs);
    }
    catch (...)
    {
        TBOX_WARNING(
            "FischerGuess::guess()\n"
            "  Unable to compute the SVD of the correlation matrix.\n"
            "  This is not a fatal error, but usually indicates that\n"
            "  the finite element solution vectors are not valid,\n"
            "  e.g., they may contain infinities or NaNs.");
        coefs.fill(0.0);
        coefs(d_n_stored_vectors - 1) = 1.0;
    }

    solution = 0.0;
    for (int i = 0; i < d_n_stored_vectors; ++i)
    {
        solution.add(coefs(i), *d_solutions[i]);
    }
    IBTK_TIMER_STOP(t_guess);
}
} // namespace IBTK
