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

#ifndef included_IBTK_IBKernelEvaluators
#define included_IBTK_IBKernelEvaluators

#include <ibtk/config.h>

#include <array>

namespace IBTK
{
/*!
 * \defgroup IBKernelEvaluators IB kernel evaluators
 * \brief Evaluate consecutive weights of a one-dimensional IB kernel.
 *
 * The argument is the displacement from the lower stencil point to the IB
 * point, divided by the grid spacing. Entry i is the weight at displacement
 * r - i. The result contains the natural number of stencil points, without
 * padding. Geometry, grid spacing factors, and tensor products are not part
 * of these evaluators.
 *
 * For an odd width N, (N - 2)/2 <= r < N/2; for an even width N,
 * N/2 - 1 <= r <= N/2. Stencil placement is the caller's responsibility.
 */
//\{
/*!
 * \brief Evaluate the 1-point B-spline kernel.
 */
class IBKernelEvaluatorBSpline1
{
public:
    //! Consecutive weights in the kernel stencil.
    using Weights = std::array<double, 1>;

    /*! \brief Return the stencil weights described in \ref IBKernelEvaluators. */
    Weights operator()(double r) const;
};

/*!
 * \brief Evaluate the 2-point B-spline kernel.
 */
class IBKernelEvaluatorBSpline2
{
public:
    //! Consecutive weights in the kernel stencil.
    using Weights = std::array<double, 2>;

    /*! \brief Return the stencil weights described in \ref IBKernelEvaluators. */
    Weights operator()(double r) const;
};

/*!
 * \brief Evaluate the 3-point B-spline kernel.
 */
class IBKernelEvaluatorBSpline3
{
public:
    //! Consecutive weights in the kernel stencil.
    using Weights = std::array<double, 3>;

    /*! \brief Return the stencil weights described in \ref IBKernelEvaluators. */
    Weights operator()(double r) const;
};

/*!
 * \brief Evaluate the 4-point B-spline kernel.
 */
class IBKernelEvaluatorBSpline4
{
public:
    //! Consecutive weights in the kernel stencil.
    using Weights = std::array<double, 4>;

    /*! \brief Return the stencil weights described in \ref IBKernelEvaluators. */
    Weights operator()(double r) const;
};

/*!
 * \brief Evaluate the 5-point B-spline kernel.
 */
class IBKernelEvaluatorBSpline5
{
public:
    //! Consecutive weights in the kernel stencil.
    using Weights = std::array<double, 5>;

    /*! \brief Return the stencil weights described in \ref IBKernelEvaluators. */
    Weights operator()(double r) const;
};

/*!
 * \brief Evaluate the 6-point B-spline kernel.
 */
class IBKernelEvaluatorBSpline6
{
public:
    //! Consecutive weights in the kernel stencil.
    using Weights = std::array<double, 6>;

    /*! \brief Return the stencil weights described in \ref IBKernelEvaluators. */
    Weights operator()(double r) const;
};

/*!
 * \brief Evaluate the 3-point IB kernel.
 */
class IBKernelEvaluatorIB3
{
public:
    //! Consecutive weights in the kernel stencil.
    using Weights = std::array<double, 3>;

    /*! \brief Return the stencil weights described in \ref IBKernelEvaluators. */
    Weights operator()(double r) const;
};

/*!
 * \brief Evaluate the 4-point IB kernel.
 */
class IBKernelEvaluatorIB4
{
public:
    //! Consecutive weights in the kernel stencil.
    using Weights = std::array<double, 4>;

    /*! \brief Return the stencil weights described in \ref IBKernelEvaluators. */
    Weights operator()(double r) const;
};

/*!
 * \brief Evaluate the 5-point IB kernel.
 */
class IBKernelEvaluatorIB5
{
public:
    //! Consecutive weights in the kernel stencil.
    using Weights = std::array<double, 5>;

    /*! \brief Return the stencil weights described in \ref IBKernelEvaluators. */
    Weights operator()(double r) const;
};

/*!
 * \brief Evaluate the 6-point IB kernel.
 */
class IBKernelEvaluatorIB6
{
public:
    //! Consecutive weights in the kernel stencil.
    using Weights = std::array<double, 6>;

    /*! \brief Return the stencil weights described in \ref IBKernelEvaluators. */
    Weights operator()(double r) const;
};

//\}
} // namespace IBTK

#include <ibtk/private/IBKernelEvaluators-inl.h>

#endif
