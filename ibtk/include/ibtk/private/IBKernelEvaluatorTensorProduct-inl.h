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
template <detail::IBKernelScalarShape Normal, detail::IBKernelScalarShape Tangential>
IBKernelEvaluatorTensorProduct<Normal, Tangential>::IBKernelEvaluatorTensorProduct(Normal evaluator)
    requires(std::same_as<Normal, Tangential>&& std::constructible_from<Normal, Normal&>&&
                 std::constructible_from<Normal, Normal&&>)
    : d_normal(evaluator), d_tangential(std::move(evaluator))
{
}

template <detail::IBKernelScalarShape Normal, detail::IBKernelScalarShape Tangential>
IBKernelEvaluatorTensorProduct<Normal, Tangential>::IBKernelEvaluatorTensorProduct(Normal normal, Tangential tangential)
    requires(std::constructible_from<Normal, Normal&&>&& std::constructible_from<Tangential, Tangential&&>)
    : d_normal(std::move(normal)), d_tangential(std::move(tangential))
{
}

template <detail::IBKernelScalarShape Normal, detail::IBKernelScalarShape Tangential>
template <int Axis>
constexpr std::array<std::size_t, NDIM>
IBKernelEvaluatorTensorProduct<Normal, Tangential>::get_stencil_widths() requires(Axis >= 0 && Axis < NDIM)
{
    std::array<std::size_t, NDIM> widths;
    widths.fill(Tangential::get_stencil_width());
    widths[Axis] = Normal::get_stencil_width();
    return widths;
}

template <detail::IBKernelScalarShape Normal, detail::IBKernelScalarShape Tangential>
template <int Axis, int Direction, class Coefficient, std::floating_point Input>
auto
IBKernelEvaluatorTensorProduct<Normal, Tangential>::evaluateDirection(const Input& r) const
{
    if constexpr (Axis == Direction)
    {
        return d_normal.template evaluate<IBKernels::Weights<Coefficient, Normal::get_stencil_width()>>(r);
    }
    else
    {
        return d_tangential.template evaluate<IBKernels::Weights<Coefficient, Tangential::get_stencil_width()>>(r);
    }
}

template <detail::IBKernelScalarShape Normal, detail::IBKernelScalarShape Tangential>
template <int Axis, IBKernelWeights Output, std::floating_point Input>
requires(detail::IBKernelCartesianShape<IBKernelEvaluatorTensorProduct<Normal, Tangential>, Axis>&&
             detail::IBKernelWritableWeights<
                 Output,
                 detail::ib_kernel_stencil_size<IBKernelEvaluatorTensorProduct<Normal, Tangential>, Axis>()>&&
                 IBKernelEvaluatorScalar<Normal, Input, typename IBKernelWeightsTraits<Output>::value_type>&&
                     IBKernelEvaluatorScalar<Tangential, Input, typename IBKernelWeightsTraits<Output>::value_type>)
    Output IBKernelEvaluatorTensorProduct<Normal, Tangential>::evaluate(const std::array<Input, NDIM>& r) const
{
    constexpr std::array<std::size_t, NDIM> widths = get_stencil_widths<Axis>();
    using Coefficient = typename IBKernelWeightsTraits<Output>::value_type;
    const auto wx = evaluateDirection<Axis, 0, Coefficient>(r[0]);
    const auto wy = evaluateDirection<Axis, 1, Coefficient>(r[1]);
    Output weights;
#if (NDIM == 3)
    const auto wz = evaluateDirection<Axis, 2, Coefficient>(r[2]);
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
