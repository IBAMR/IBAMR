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

#ifndef included_IBTK_IBKernelTensorProductEvaluator_inl
#define included_IBTK_IBKernelTensorProductEvaluator_inl

#include <ibtk/config.h>

#include <ibtk/IBKernelTensorProductEvaluator.h>

#include <utility>

namespace IBTK
{
template <class NormalEvaluator, class TangentialEvaluator>
inline IBKernelTensorProductEvaluator<NormalEvaluator, TangentialEvaluator>::IBKernelTensorProductEvaluator(
    NormalEvaluator evaluator)
    : d_normal(evaluator), d_tangential(std::move(evaluator))
{
    static_assert(std::is_same_v<NormalEvaluator, TangentialEvaluator>, "An isotropic product uses one evaluator type");
}

template <class NormalEvaluator, class TangentialEvaluator>
inline IBKernelTensorProductEvaluator<NormalEvaluator, TangentialEvaluator>::IBKernelTensorProductEvaluator(
    NormalEvaluator normal_evaluator,
    TangentialEvaluator tangential_evaluator)
    : d_normal(std::move(normal_evaluator)), d_tangential(std::move(tangential_evaluator))
{
}

template <class NormalEvaluator, class TangentialEvaluator>
template <int Axis>
inline constexpr std::array<int, NDIM>
IBKernelTensorProductEvaluator<NormalEvaluator, TangentialEvaluator>::get_stencil_widths()
{
    static_assert(Axis >= 0 && Axis < NDIM, "Invalid tensor-product axis");
    constexpr int normal_width = std::tuple_size<NormalWeights>::value;
    constexpr int tangential_width = std::tuple_size<TangentialWeights>::value;
    std::array<int, NDIM> widths = {};
    for (std::size_t d = 0; d < NDIM; ++d) widths[d] = d == Axis ? normal_width : tangential_width;
    return widths;
}

template <class NormalEvaluator, class TangentialEvaluator>
template <int Axis>
inline auto
IBKernelTensorProductEvaluator<NormalEvaluator, TangentialEvaluator>::evaluateFactors(
    const std::array<double, NDIM>& r) const
{
    static_assert(Axis >= 0 && Axis < NDIM, "Invalid tensor-product axis");
    const auto evaluate_direction = [&](auto direction)
    {
        constexpr int d = decltype(direction)::value;
        if constexpr (d == Axis)
            return d_normal(r[d]);
        else
            return d_tangential(r[d]);
    };
    if constexpr (NDIM == 2)
        return std::make_tuple(evaluate_direction(std::integral_constant<int, 0>{}),
                               evaluate_direction(std::integral_constant<int, 1>{}));
    else
        return std::make_tuple(evaluate_direction(std::integral_constant<int, 0>{}),
                               evaluate_direction(std::integral_constant<int, 1>{}),
                               evaluate_direction(std::integral_constant<int, 2>{}));
}

template <class NormalEvaluator, class TangentialEvaluator>
template <int Axis>
inline auto
IBKernelTensorProductEvaluator<NormalEvaluator, TangentialEvaluator>::evaluate(const std::array<double, NDIM>& r) const
{
    constexpr auto widths = get_stencil_widths<Axis>();
    const auto factors = evaluateFactors<Axis>(r);
    const auto& wx = std::get<0>(factors);
    const auto& wy = std::get<1>(factors);
    constexpr int nx = widths[0], ny = widths[1];
    if constexpr (NDIM == 2)
    {
        std::array<double, nx * ny> weights;
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) weights[i + nx * j] = wx[i] * wy[j];
        return weights;
    }
    else
    {
        const auto& wz = std::get<2>(factors);
        constexpr int nz = widths[2];
        std::array<double, nx * ny * nz> weights;
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) weights[i + nx * (j + ny * k)] = wx[i] * wy[j] * wz[k];
        return weights;
    }
}
} // namespace IBTK

#endif
