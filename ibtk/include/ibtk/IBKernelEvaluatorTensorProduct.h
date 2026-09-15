// ---------------------------------------------------------------------
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
// This file is part of IBAMR and is distributed under the 3-clause BSD
// license. The full text of the license can be found in COPYRIGHT.
// ---------------------------------------------------------------------

#ifndef included_IBTK_IBKernelEvaluatorTensorProduct
#define included_IBTK_IBKernelEvaluatorTensorProduct
#include <ibtk/config.h>

#include <ibtk/IBKernelConcepts.h>

#include <array>
#include <type_traits>

namespace IBTK
{
/*!
 * \brief Evaluate a tensor product of one-dimensional IB kernels.
 *
 * Axis selects the coordinate using Normal; other coordinates use Tangential.
 * The scalar evaluators are stored by value. Both factors and their product
 * use the coefficient type selected by Output.
 *
 * \see IBKernelEvaluatorCartesian
 */
template <detail::IBKernelScalarShape Normal, detail::IBKernelScalarShape Tangential = Normal>
class IBKernelEvaluatorTensorProduct
{
public:
    /*! \brief Use copies of the same scalar evaluator in every coordinate. */
    explicit IBKernelEvaluatorTensorProduct(Normal evaluator)
        requires(std::same_as<Normal, Tangential>&& std::constructible_from<Normal, Normal&>&&
                     std::constructible_from<Normal, Normal&&>);

    /*! \brief Use separate normal and tangential evaluators. */
    IBKernelEvaluatorTensorProduct(Normal normal, Tangential tangential)
        requires(std::constructible_from<Normal, Normal&&>&& std::constructible_from<Tangential, Tangential&&>);

    /*! \brief Return coordinate widths with Axis selecting the normal factor. */
    template <int Axis>
    static constexpr std::array<std::size_t, NDIM> get_stencil_widths() requires(Axis >= 0 && Axis < NDIM);

    /*! \brief Return the product of the scalar weights at r. */
    template <int Axis, IBKernelWeights Output, std::floating_point Input>
    requires(
        detail::IBKernelCartesianShape<IBKernelEvaluatorTensorProduct, Axis>&&
            detail::IBKernelWritableWeights<Output,
                                            detail::ib_kernel_stencil_size<IBKernelEvaluatorTensorProduct, Axis>()>&&
                IBKernelEvaluatorScalar<Normal, Input, typename IBKernelWeightsTraits<Output>::value_type>&&
                    IBKernelEvaluatorScalar<Tangential, Input, typename IBKernelWeightsTraits<Output>::value_type>)
        Output evaluate(const std::array<Input, NDIM>& r) const;

private:
    /*! \brief Evaluate the scalar factor for coordinate Direction. */
    template <int Axis, int Direction, class Coefficient, std::floating_point Input>
    auto evaluateDirection(const Input& r) const;

    //! Owned scalar factors.
    [[no_unique_address]] Normal d_normal;
    [[no_unique_address]] Tangential d_tangential;
};
} // namespace IBTK
#include <ibtk/private/IBKernelEvaluatorTensorProduct-inl.h>
#endif
