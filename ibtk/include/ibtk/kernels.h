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

#ifndef included_IBTK_kernels
#define included_IBTK_kernels

#include <ibtk/config.h>

#include <array>
#include <cstddef>

namespace IBTK
{
/*!
 * \brief One-dimensional IB kernel functions returning consecutive weights.
 *
 * For a normalized grid coordinate q and width N, matrix construction selects
 * first index floor(q + 1/2) - (N-1)/2 for odd N. For even N it selects
 * floor(q) - N/2 + 1 in the side-normal direction and ceil(q) - N/2 in
 * tangential directions. Here q includes the grid-data centering. The argument r is q
 * minus this first index; returned entry i is phi(r-i), for 0 <= i < N.
 *
 * For an odd width N, (N - 2)/2 <= r < N/2; for an even width N,
 * N/2 - 1 <= r <= N/2. No grid-spacing factors are included in the weights.
 */
namespace Kernels
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
    //! Consecutive weights in the kernel stencil.
    using Weights = std::array<double, N>;

    /*! \brief Return the stencil weights described in \ref Kernels. */
    Weights operator()(double r) const;
};

/*!
 * \brief Evaluate the 3-point IB kernel.
 */
class IB3
{
public:
    //! Consecutive weights in the kernel stencil.
    using Weights = std::array<double, 3>;

    /*! \brief Return the stencil weights described in \ref Kernels. */
    Weights operator()(double r) const;
};

/*!
 * \brief Evaluate the 4-point IB kernel.
 */
class IB4
{
public:
    //! Consecutive weights in the kernel stencil.
    using Weights = std::array<double, 4>;

    /*! \brief Return the stencil weights described in \ref Kernels. */
    Weights operator()(double r) const;
};

/*!
 * \brief Evaluate the 5-point IB kernel.
 */
class IB5
{
public:
    //! Consecutive weights in the kernel stencil.
    using Weights = std::array<double, 5>;

    /*! \brief Return the stencil weights described in \ref Kernels. */
    Weights operator()(double r) const;
};

/*!
 * \brief Evaluate the 6-point IB kernel.
 */
class IB6
{
public:
    //! Consecutive weights in the kernel stencil.
    using Weights = std::array<double, 6>;

    /*! \brief Return the stencil weights described in \ref Kernels. */
    Weights operator()(double r) const;
};

} // namespace Kernels
} // namespace IBTK

#include <ibtk/private/kernels-inl.h>

#endif
