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

#ifndef included_IBTK_IBKernelEvaluatorTensorProduct_inl
#define included_IBTK_IBKernelEvaluatorTensorProduct_inl
#include <ibtk/config.h>

#include <ibtk/IBKernelEvaluatorTensorProduct.h>

#include <utility>

namespace IBTK
{
template <IBKernelScalarStencil NormalEvaluator, IBKernelScalarStencil TransverseEvaluator>
IBKernelEvaluatorTensorProduct<NormalEvaluator, TransverseEvaluator>::IBKernelEvaluatorTensorProduct(
    NormalEvaluator evaluator) requires(std::same_as<NormalEvaluator, TransverseEvaluator>&&
                                            std::constructible_from<NormalEvaluator, NormalEvaluator&>&&
                                                std::constructible_from<NormalEvaluator, NormalEvaluator&&>)
    : d_normal_evaluator(evaluator), d_transverse_evaluator(std::move(evaluator))
{
}

template <IBKernelScalarStencil NormalEvaluator, IBKernelScalarStencil TransverseEvaluator>
IBKernelEvaluatorTensorProduct<NormalEvaluator, TransverseEvaluator>::IBKernelEvaluatorTensorProduct(
    NormalEvaluator normal_evaluator,
    TransverseEvaluator transverse_evaluator)
    requires(std::constructible_from<NormalEvaluator, NormalEvaluator&&>&&
                 std::constructible_from<TransverseEvaluator, TransverseEvaluator&&>)
    : d_normal_evaluator(std::move(normal_evaluator)), d_transverse_evaluator(std::move(transverse_evaluator))
{
}

template <IBKernelScalarStencil NormalEvaluator, IBKernelScalarStencil TransverseEvaluator>
template <int Axis>
constexpr std::array<std::size_t, NDIM>
IBKernelEvaluatorTensorProduct<NormalEvaluator, TransverseEvaluator>::get_stencil_widths()
    requires(Axis >= 0 && Axis < NDIM)
{
    std::array<std::size_t, NDIM> widths;
    widths.fill(TransverseEvaluator::get_stencil_width());
    widths[Axis] = NormalEvaluator::get_stencil_width();
    return widths;
}

template <IBKernelScalarStencil NormalEvaluator, IBKernelScalarStencil TransverseEvaluator>
template <int Axis, int Direction, class Coefficient, std::floating_point Input>
std::conditional_t<Axis == Direction,
                   IBKernelEvaluators::Weights<Coefficient, NormalEvaluator::get_stencil_width()>,
                   IBKernelEvaluators::Weights<Coefficient, TransverseEvaluator::get_stencil_width()>>
IBKernelEvaluatorTensorProduct<NormalEvaluator, TransverseEvaluator>::evaluateDirection(const Input& r) const
{
    if constexpr (Axis == Direction)
    {
        return d_normal_evaluator
            .template evaluate<IBKernelEvaluators::Weights<Coefficient, NormalEvaluator::get_stencil_width()>>(r);
    }
    else
    {
        return d_transverse_evaluator
            .template evaluate<IBKernelEvaluators::Weights<Coefficient, TransverseEvaluator::get_stencil_width()>>(r);
    }
}

template <IBKernelScalarStencil NormalEvaluator, IBKernelScalarStencil TransverseEvaluator>
template <int Axis, IBKernelWeights Output, std::floating_point Input>
requires detail::IBKernelTensorProductEvaluable<IBKernelEvaluatorTensorProduct<NormalEvaluator, TransverseEvaluator>,
                                                NormalEvaluator,
                                                TransverseEvaluator,
                                                Axis,
                                                Output,
                                                Input>
    Output
    IBKernelEvaluatorTensorProduct<NormalEvaluator, TransverseEvaluator>::evaluate(
        const std::array<Input, NDIM>& r) const
{
    constexpr std::array<std::size_t, NDIM> widths = get_stencil_widths<Axis>();
    using Coefficient = ib_kernel_weights_value_t<Output>;
    using NormalWeights = IBKernelEvaluators::Weights<Coefficient, NormalEvaluator::get_stencil_width()>;
    using TransverseWeights = IBKernelEvaluators::Weights<Coefficient, TransverseEvaluator::get_stencil_width()>;
    const std::conditional_t<Axis == 0, NormalWeights, TransverseWeights> wx =
        evaluateDirection<Axis, 0, Coefficient>(r[0]);
    const std::conditional_t<Axis == 1, NormalWeights, TransverseWeights> wy =
        evaluateDirection<Axis, 1, Coefficient>(r[1]);
    Output weights{};
#if (NDIM == 3)
    const std::conditional_t<Axis == 2, NormalWeights, TransverseWeights> wz =
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
