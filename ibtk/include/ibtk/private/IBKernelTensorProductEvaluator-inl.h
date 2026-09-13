// ---------------------------------------------------------------------
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
// This file is part of IBAMR and is distributed under the 3-clause BSD
// license. The full text of the license can be found in COPYRIGHT.
// ---------------------------------------------------------------------

#ifndef included_IBTK_IBKernelTensorProductEvaluator_inl
#define included_IBTK_IBKernelTensorProductEvaluator_inl
#include <ibtk/config.h>

#include <ibtk/IBKernelTensorProductEvaluator.h>

#include <utility>

namespace IBTK
{
template <ScalarKernel Normal, ScalarKernel Tangential>
IBKernelTensorProductEvaluator<Normal, Tangential>::IBKernelTensorProductEvaluator(
    Normal evaluator) requires std::same_as<Normal, Tangential>&& std::copy_constructible<Normal>
    : d_normal(evaluator), d_tangential(std::move(evaluator))
{
}

template <ScalarKernel Normal, ScalarKernel Tangential>
IBKernelTensorProductEvaluator<Normal, Tangential>::IBKernelTensorProductEvaluator(Normal normal, Tangential tangential)
    : d_normal(std::move(normal)), d_tangential(std::move(tangential))
{
}

template <ScalarKernel Normal, ScalarKernel Tangential>
template <int Axis>
requires(Axis >= 0 &&
         Axis <
             NDIM) constexpr std::array<std::size_t,
                                        NDIM> IBKernelTensorProductEvaluator<Normal, Tangential>::get_stencil_widths()
{
    std::array<std::size_t, NDIM> widths;
    widths.fill(KernelWeightTraits<typename Tangential::Weights>::extent);
    widths[Axis] = KernelWeightTraits<typename Normal::Weights>::extent;
    return widths;
}

template <ScalarKernel Normal, ScalarKernel Tangential>
template <int Axis, int Direction>
auto
IBKernelTensorProductEvaluator<Normal, Tangential>::evaluateDirection(double r) const
{
    if constexpr (Axis == Direction)
    {
        return d_normal(r);
    }
    else
    {
        return d_tangential(r);
    }
}

template <ScalarKernel Normal, ScalarKernel Tangential>
template <int Axis>
requires(Axis >= 0 && Axis < NDIM) typename IBKernelTensorProductEvaluator<Normal, Tangential>::Weights
    IBKernelTensorProductEvaluator<Normal, Tangential>::evaluate(const std::array<double, NDIM>& r) const
{
    constexpr std::array<std::size_t, NDIM> widths = get_stencil_widths<Axis>();
    const auto wx = evaluateDirection<Axis, 0>(r[0]);
    const auto wy = evaluateDirection<Axis, 1>(r[1]);
    Weights weights;
#if (NDIM == 3)
    const auto wz = evaluateDirection<Axis, 2>(r[2]);
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
