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

#ifndef included_IBTK_ib_kernel_evaluators
#define included_IBTK_ib_kernel_evaluators

#include <ibtk/config.h>

#include <ibtk/IBKernelConcepts.h>

#include <cstddef>

namespace IBTK
{
/*!
 * \brief Supplied one-dimensional IB kernel evaluators.
 *
 * Arithmetic uses Output's real coefficient type, including conversion of r.
 * Supplied evaluators require default-initializable, movable Output storage
 * with writable indexed coefficients and the exact stencil extent.
 *
 * \see IBKernelEvaluatorScalar for the coordinate and coefficient conventions.
 */
namespace IBKernelEvaluators
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
    template <detail::IBKernelWritableWeights<N> Output, std::floating_point Input>
    Output evaluate(Input r) const;

    /*! \brief Return the number of stencil points. */
    static constexpr std::size_t get_stencil_width();
};

/*!
 * \brief Evaluate the 3-point IB kernel.
 */
class IB3
{
public:
    /*! \brief Return the three stencil weights. */
    template <detail::IBKernelWritableWeights<3> Output, std::floating_point Input>
    Output evaluate(Input r) const;

    /*! \brief Return the number of stencil points. */
    static constexpr std::size_t get_stencil_width();
};

/*!
 * \brief Evaluate the 4-point IB kernel.
 */
class IB4
{
public:
    /*! \brief Return the four stencil weights. */
    template <detail::IBKernelWritableWeights<4> Output, std::floating_point Input>
    Output evaluate(Input r) const;

    /*! \brief Return the number of stencil points. */
    static constexpr std::size_t get_stencil_width();
};

/*!
 * \brief Evaluate the 5-point IB kernel.
 */
class IB5
{
public:
    /*! \brief Return the five stencil weights. */
    template <detail::IBKernelWritableWeights<5> Output, std::floating_point Input>
    Output evaluate(Input r) const;

    /*! \brief Return the number of stencil points. */
    static constexpr std::size_t get_stencil_width();
};

/*!
 * \brief Evaluate the 6-point IB kernel.
 */
class IB6
{
public:
    /*! \brief Return the six stencil weights. */
    template <detail::IBKernelWritableWeights<6> Output, std::floating_point Input>
    Output evaluate(Input r) const;

    /*! \brief Return the number of stencil points. */
    static constexpr std::size_t get_stencil_width();
};

} // namespace IBKernelEvaluators
} // namespace IBTK

#include <ibtk/private/ib_kernel_evaluators-inl.h>

#endif
