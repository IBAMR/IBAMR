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

#ifndef included_IBTK_ib_kernels
#define included_IBTK_ib_kernels

#include <ibtk/config.h>

#include <ibtk/IBKernelConcepts.h>

#include <cstddef>

namespace IBTK
{
/*!
 * \brief Supplied one-dimensional IB kernel evaluators.
 *
 * \see IBKernelEvaluatorScalar for the coordinate and coefficient conventions.
 */
namespace IBKernels
{
/*!
 * \brief Evaluate the N-point cardinal B-spline kernel of degree N - 1.
 *
 * N must be positive.
 */
template <std::size_t N>
requires(N > 0) class BSpline
{
public:
    /*! \brief Return the N stencil weights. */
    IBKernels::Weights<double, N> operator()(double r) const;
};

/*!
 * \brief Evaluate the 3-point IB kernel.
 */
class IB3
{
public:
    /*! \brief Return the three stencil weights. */
    IBKernels::Weights<double, 3> operator()(double r) const;
};

/*!
 * \brief Evaluate the 4-point IB kernel.
 */
class IB4
{
public:
    /*! \brief Return the four stencil weights. */
    IBKernels::Weights<double, 4> operator()(double r) const;
};

/*!
 * \brief Evaluate the 5-point IB kernel.
 */
class IB5
{
public:
    /*! \brief Return the five stencil weights. */
    IBKernels::Weights<double, 5> operator()(double r) const;
};

/*!
 * \brief Evaluate the 6-point IB kernel.
 */
class IB6
{
public:
    /*! \brief Return the six stencil weights. */
    IBKernels::Weights<double, 6> operator()(double r) const;
};

} // namespace IBKernels
} // namespace IBTK

#include <ibtk/private/ib_kernels-inl.h>

#endif
