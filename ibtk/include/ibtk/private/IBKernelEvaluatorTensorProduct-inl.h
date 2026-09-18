// ---------------------------------------------------------------------
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
// This file is part of IBAMR and is distributed under the 3-clause BSD
// license. The full text of the license can be found in COPYRIGHT.
// ---------------------------------------------------------------------

#ifndef included_IBTK_IBKernelEvaluatorTensorProduct_inl
#define included_IBTK_IBKernelEvaluatorTensorProduct_inl
#include <ibtk/config.h>

#include <ibtk/IBKernelEvaluatorTensorProduct.h>

#include <utility>

namespace IBTK
{
template <IBKernelEvaluatorScalar Normal, IBKernelEvaluatorScalar Tangential>
IBKernelEvaluatorTensorProduct<Normal, Tangential>::IBKernelEvaluatorTensorProduct(
    Normal evaluator) requires std::same_as<Normal, Tangential> : d_normal(evaluator),
                                                                  d_tangential(std::move(evaluator))
{
}

template <IBKernelEvaluatorScalar Normal, IBKernelEvaluatorScalar Tangential>
IBKernelEvaluatorTensorProduct<Normal, Tangential>::IBKernelEvaluatorTensorProduct(Normal normal, Tangential tangential)
    requires(std::constructible_from<Normal, Normal&&>&& std::constructible_from<Tangential, Tangential&&>)
    : d_normal(std::move(normal)), d_tangential(std::move(tangential))
{
}

template <IBKernelEvaluatorScalar Normal, IBKernelEvaluatorScalar Tangential>
template <int Axis>
constexpr std::array<std::size_t, NDIM>
IBKernelEvaluatorTensorProduct<Normal, Tangential>::get_stencil_widths() requires(Axis >= 0 && Axis < NDIM)
{
    std::array<std::size_t, NDIM> widths;
    widths.fill(Tangential::get_stencil_width());
    widths[Axis] = Normal::get_stencil_width();
    return widths;
}

template <IBKernelEvaluatorScalar Normal, IBKernelEvaluatorScalar Tangential>
template <int Axis, int Direction, class Coefficient, std::floating_point Input>
std::conditional_t<Axis == Direction,
                   IBKernelEvaluators::Weights<Coefficient, Normal::get_stencil_width()>,
                   IBKernelEvaluators::Weights<Coefficient, Tangential::get_stencil_width()>>
IBKernelEvaluatorTensorProduct<Normal, Tangential>::evaluateDirection(const Input& r) const
{
    if constexpr (Axis == Direction)
    {
        return d_normal.template evaluate<IBKernelEvaluators::Weights<Coefficient, Normal::get_stencil_width()>>(r);
    }
    else
    {
        return d_tangential
            .template evaluate<IBKernelEvaluators::Weights<Coefficient, Tangential::get_stencil_width()>>(r);
    }
}

template <IBKernelEvaluatorScalar Normal, IBKernelEvaluatorScalar Tangential>
template <int Axis, IBKernelWeights Output, std::floating_point Input>
requires(detail::IBKernelWritableWeights<
         Output,
         detail::ib_kernel_stencil_size<IBKernelEvaluatorTensorProduct<Normal, Tangential>, Axis>()>&&
             IBKernelEvaluatorScalar<Normal, Input, ib_kernel_weights_value_t<Output>>&&
                 IBKernelEvaluatorScalar<Tangential, Input, ib_kernel_weights_value_t<Output>>) Output
    IBKernelEvaluatorTensorProduct<Normal, Tangential>::evaluate(const std::array<Input, NDIM>& r) const
{
    constexpr std::array<std::size_t, NDIM> widths = get_stencil_widths<Axis>();
    using Coefficient = ib_kernel_weights_value_t<Output>;
    using NormalWeights = IBKernelEvaluators::Weights<Coefficient, Normal::get_stencil_width()>;
    using TangentialWeights = IBKernelEvaluators::Weights<Coefficient, Tangential::get_stencil_width()>;
    const std::conditional_t<Axis == 0, NormalWeights, TangentialWeights> wx =
        evaluateDirection<Axis, 0, Coefficient>(r[0]);
    const std::conditional_t<Axis == 1, NormalWeights, TangentialWeights> wy =
        evaluateDirection<Axis, 1, Coefficient>(r[1]);
    Output weights{};
#if (NDIM == 3)
    const std::conditional_t<Axis == 2, NormalWeights, TangentialWeights> wz =
        evaluateDirection<Axis, 2, Coefficient>(r[2]);
    for (std::size_t k = 0; k < widths[2]; ++k)
    {
#endif
        for (std::size_t j = 0; j < widths[1]; ++j)
        {
            for (std::size_t i = 0; i < widths[0]; ++i)
            {
#if (NDIM == 3)
                weights[i + widths[0] * (j + widths[1] * k)] = wx[i] * wy[j] * wz[k];
#else
            weights[i + widths[0] * j] = wx[i] * wy[j];
#endif
            }
        }
#if (NDIM == 3)
    }
#endif
    return weights;
}
} // namespace IBTK
#endif
