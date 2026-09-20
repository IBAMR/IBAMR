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

#ifndef included_IBTK_IBKernelEvaluatorTensorProduct
#define included_IBTK_IBKernelEvaluatorTensorProduct
#include <ibtk/config.h>

#include <ibtk/IBKernelConcepts.h>

#include <array>
#include <type_traits>

namespace IBTK
{
namespace detail
{
/*! \brief The requirements of IBKernelEvaluatorTensorProduct::evaluate(). */
template <class Product, class NormalEvaluator, class TransverseEvaluator, int Axis, class Output, class Input>
concept IBKernelTensorProductEvaluable =
    Axis >= 0 && Axis < NDIM && IBKernelWritableWeights<Output, ib_kernel_stencil_size<Product, Axis>()> &&
    IBKernelEvaluatorScalar<NormalEvaluator, Input, ib_kernel_weights_value_t<Output>> &&
    IBKernelEvaluatorScalar<TransverseEvaluator, Input, ib_kernel_weights_value_t<Output>>;
} // namespace detail

/*!
 * \brief Tensor product of one-dimensional IB kernel evaluators.
 *
 * With one evaluator, that evaluator is used in every coordinate. With two,
 * the normal evaluator is used along Axis and the transverse evaluator in the
 * other NDIM - 1 coordinates. Evaluators are stored by value. The weights of
 * each factor and of the product have the coefficient type of Output.
 *
 * \see IBKernelEvaluatorCartesian
 */
template <IBKernelScalarStencil NormalEvaluator, IBKernelScalarStencil TransverseEvaluator = NormalEvaluator>
class IBKernelEvaluatorTensorProduct
{
public:
    /*! \brief Use evaluator in every coordinate. */
    explicit IBKernelEvaluatorTensorProduct(NormalEvaluator evaluator)
        requires(std::same_as<NormalEvaluator, TransverseEvaluator>&&
                     std::constructible_from<NormalEvaluator, NormalEvaluator&>&&
                         std::constructible_from<NormalEvaluator, NormalEvaluator&&>);

    /*! \brief Use normal_evaluator along Axis and transverse_evaluator in the other coordinates. */
    IBKernelEvaluatorTensorProduct(NormalEvaluator normal_evaluator, TransverseEvaluator transverse_evaluator)
        requires(std::constructible_from<NormalEvaluator, NormalEvaluator&&>&&
                     std::constructible_from<TransverseEvaluator, TransverseEvaluator&&>);

    /*!
     * \brief Return the stencil width in each coordinate: the normal
     * evaluator's at index Axis, the transverse evaluator's elsewhere.
     */
    template <int Axis>
    static constexpr std::array<std::size_t, NDIM> get_stencil_widths() requires(Axis >= 0 && Axis < NDIM);

    /*!
     * \brief Return the product of the one-dimensional weights at r: the
     * normal evaluator's at r[Axis], the transverse evaluator's at every other
     * r[d].
     */
    template <int Axis, IBKernelWeights Output, std::floating_point Input>
    requires detail::IBKernelTensorProductEvaluable<IBKernelEvaluatorTensorProduct,
                                                    NormalEvaluator,
                                                    TransverseEvaluator,
                                                    Axis,
                                                    Output,
                                                    Input>
        Output evaluate(const std::array<Input, NDIM>& r) const;

private:
    /*! \brief Evaluate the one-dimensional factor for coordinate Direction. */
    template <int Axis, int Direction, class Coefficient, std::floating_point Input>
    std::conditional_t<Axis == Direction,
                       IBKernelEvaluators::Weights<Coefficient, NormalEvaluator::get_stencil_width()>,
                       IBKernelEvaluators::Weights<Coefficient, TransverseEvaluator::get_stencil_width()>>
    evaluateDirection(const Input& r) const;

    [[no_unique_address]] NormalEvaluator d_normal_evaluator;
    [[no_unique_address]] TransverseEvaluator d_transverse_evaluator;
};
} // namespace IBTK
#include <ibtk/private/IBKernelEvaluatorTensorProduct-inl.h>
#endif
