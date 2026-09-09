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
 * and return std::array<double, N>, where N > 0 is its natural stencil width.
 * A single evaluator applies in every direction; with two evaluators, the
 * first applies along the selected axis and the second in the remaining
 * directions. Evaluators are stored by value. Axis identifies x (0), y (1),
 * or z (2), with Axis < NDIM. For side-centered data, Axis is the face-normal
 * direction.
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
 * std::array<double, NDIM> r;
 * r.fill(1.0);
 * r[1] = 1.5;
 * const std::array<double, NDIM == 2 ? 12 : 36> weights = evaluator.evaluate<1>(r);
 * // weights contains 3 * 4 coefficients in 2D, or 3 * 4 * 3 in 3D.
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
    template <int Axis>
    static constexpr std::array<int, NDIM> get_stencil_widths();

    /*!
     * \brief Return the tensor-product coefficients with coordinate zero varying fastest.
     *
     * Each scalar evaluator is called once per corresponding direction.
     */
    template <int Axis>
    auto evaluate(const std::array<double, NDIM>& r) const;

private:
    /*! \brief Return a tuple of one-dimensional weight arrays in coordinate order. */
    template <int Axis>
    auto evaluateFactors(const std::array<double, NDIM>& r) const;

    //! Normal kernel weights.
    using NormalWeights = std::invoke_result_t<const NormalEvaluator&, double>;
    //! Tangential kernel weights.
    using TangentialWeights = std::invoke_result_t<const TangentialEvaluator&, double>;

    static_assert(std::tuple_size<NormalWeights>::value > 0 && std::tuple_size<TangentialWeights>::value > 0,
                  "Kernel stencils must be nonempty");
    static_assert(std::is_same_v<NormalWeights, std::array<double, std::tuple_size<NormalWeights>::value>> &&
                      std::is_same_v<TangentialWeights, std::array<double, std::tuple_size<TangentialWeights>::value>>,
                  "Kernel evaluators must return std::array<double, N>");

    //! Owned normal kernel evaluator.
    NormalEvaluator d_normal;
    //! Owned tangential kernel evaluator.
    TangentialEvaluator d_tangential;
};
} // namespace IBTK

#include <ibtk/private/IBKernelTensorProductEvaluator-inl.h>

#endif
