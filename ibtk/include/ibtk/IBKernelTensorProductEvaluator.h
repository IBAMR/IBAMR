// ---------------------------------------------------------------------
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
// This file is part of IBAMR and is distributed under the 3-clause BSD
// license. The full text of the license can be found in COPYRIGHT.
// ---------------------------------------------------------------------

#ifndef included_IBTK_IBKernelTensorProductEvaluator
#define included_IBTK_IBKernelTensorProductEvaluator
#include <ibtk/config.h>

#include <ibtk/KernelConcepts.h>

#include <array>
#include <type_traits>

namespace IBTK
{
/*!
 * \brief Evaluate normal/tangential products of scalar kernels.
 *
 * Axis selects the normal coordinate: x (0), y (1), or z (2).
 * Other coordinates use the tangential kernel. Each r[d] follows the
 * ScalarKernel coordinate convention. Returned coefficients have coordinate
 * zero varying fastest. Evaluators are owned by value.
 */
template <ScalarKernel Normal, ScalarKernel Tangential = Normal>
class IBKernelTensorProductEvaluator
{
public:
    //! Tensor-product coefficients.
    using Weights =
        std::array<std::common_type_t<typename KernelWeightTraits<typename Normal::Weights>::value_type,
                                      typename KernelWeightTraits<typename Tangential::Weights>::value_type>,
                   KernelWeightTraits<typename Normal::Weights>::extent *
                       KernelWeightTraits<typename Tangential::Weights>::extent
#if (NDIM == 3)
                       * KernelWeightTraits<typename Tangential::Weights>::extent
#endif
                   >;

    /*! \brief Copy a scalar evaluator into both factor roles. */
    explicit IBKernelTensorProductEvaluator(
        Normal evaluator) requires std::same_as<Normal, Tangential>&& std::copy_constructible<Normal>;

    /*! \brief Own the supplied normal and tangential evaluators. */
    IBKernelTensorProductEvaluator(Normal normal, Tangential tangential);

    /*! \brief Return coordinate widths with Axis selecting the normal factor. */
    template <int Axis>
    requires(Axis >= 0 && Axis < NDIM) static constexpr std::array<std::size_t, NDIM> get_stencil_widths();

    /*! \brief Return the product of the scalar weights at r. */
    template <int Axis>
    requires(Axis >= 0 && Axis < NDIM) Weights evaluate(const std::array<double, NDIM>& r) const;

private:
    /*! \brief Evaluate the scalar factor for coordinate Direction. */
    template <int Axis, int Direction>
    auto evaluateDirection(double r) const;

    //! Owned scalar factors.
    [[no_unique_address]] Normal d_normal;
    [[no_unique_address]] Tangential d_tangential;
};
} // namespace IBTK
#include <ibtk/private/IBKernelTensorProductEvaluator-inl.h>
#endif
