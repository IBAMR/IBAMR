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

#ifndef included_IBTK_IBKernelTensorProductEvaluator
#define included_IBTK_IBKernelTensorProductEvaluator

#include <ibtk/config.h>

#include <array>
#include <tuple>
#include <type_traits>

namespace IBTK
{
/*!
 * \brief Evaluate a tensor product of one-dimensional IB kernels.
 *
 * Each scalar evaluator must accept a double through a const call operator
 * and return std::array<double, N> with its natural stencil width. A single evaluator
 * applies in every direction; with two evaluators, the first applies along
 * the selected axis and the second in the remaining directions.
 * Evaluators are stored by value. Dim must be two or three, and Axis identifies
 * x (0), y (1), or z (2), with Axis < Dim. For side-centered data, Axis is the
 * face-normal direction.
 *
 * The input contains displacements from the lower stencil point in grid
 * units, following \ref IBKernelEvaluators. The caller must account for the
 * grid-data centering when computing these displacements. Output coefficients
 * are ordered with coordinate zero varying fastest. No grid-spacing factors
 * are applied.
 *
 * For example:
 * \code
 * IBKernelTensorProductEvaluator evaluator{IBKernelEvaluatorIB4{}, IBKernelEvaluatorIB3{}};
 * const auto weights = evaluator.evaluate<1>(std::array<double, 3>{1.0, 1.5, 1.0});
 * // weights contains 3 * 4 * 3 coefficients.
 * \endcode
 */
template <class NormalEvaluator, class TangentialEvaluator = NormalEvaluator>
class IBKernelTensorProductEvaluator
{
public:
    /*! \brief Construct an isotropic product of a copyable scalar evaluator. */
    explicit IBKernelTensorProductEvaluator(NormalEvaluator evaluator);

    /*! \brief Construct a product of normal and tangential kernel evaluators. */
    IBKernelTensorProductEvaluator(NormalEvaluator normal_evaluator, TangentialEvaluator tangential_evaluator);

    /*! \brief Return the stencil width in each coordinate direction. */
    template <int Axis, std::size_t Dim>
    static constexpr std::array<int, Dim> get_stencil_widths();

    /*!
     * \brief Return a tuple of one-dimensional weight arrays in coordinate order.
     *
     * Each scalar evaluator is called once per corresponding direction.
     * Consumers may form products while accumulating values without storing
     * the full tensor product.
     */
    template <int Axis, std::size_t Dim>
    auto evaluateFactors(const std::array<double, Dim>& r) const;

    /*! \brief Return the tensor-product coefficients with coordinate zero varying fastest. */
    template <int Axis, std::size_t Dim>
    auto evaluate(const std::array<double, Dim>& r) const;

private:
    //! Normal kernel weights.
    using NormalWeights = std::invoke_result_t<const NormalEvaluator&, double>;
    //! Tangential kernel weights.
    using TangentialWeights = std::invoke_result_t<const TangentialEvaluator&, double>;

    //! Owned normal kernel evaluator.
    NormalEvaluator d_normal;
    //! Owned tangential kernel evaluator.
    TangentialEvaluator d_tangential;
};
} // namespace IBTK

#include <ibtk/private/IBKernelTensorProductEvaluator-inl.h>

#endif
