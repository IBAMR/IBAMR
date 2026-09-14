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
template <IBKernelEvaluatorScalar Normal, detail::IBKernelEvaluatorScalarPair<Normal> Tangential>
IBKernelEvaluatorTensorProduct<Normal, Tangential>::IBKernelEvaluatorTensorProduct(Normal evaluator)
    requires(std::same_as<Normal, Tangential>&& std::constructible_from<Normal, Normal&>&&
                 std::constructible_from<Normal, Normal&&>)
    : d_normal(evaluator), d_tangential(std::move(evaluator))
{
}

template <IBKernelEvaluatorScalar Normal, detail::IBKernelEvaluatorScalarPair<Normal> Tangential>
IBKernelEvaluatorTensorProduct<Normal, Tangential>::IBKernelEvaluatorTensorProduct(Normal normal, Tangential tangential)
    requires(std::constructible_from<Normal, Normal&&>&& std::constructible_from<Tangential, Tangential&&>)
    : d_normal(std::move(normal)), d_tangential(std::move(tangential))
{
}

template <IBKernelEvaluatorScalar Normal, detail::IBKernelEvaluatorScalarPair<Normal> Tangential>
template <int Axis>
constexpr std::array<std::size_t, NDIM>
IBKernelEvaluatorTensorProduct<Normal, Tangential>::get_stencil_widths() requires(Axis >= 0 && Axis < NDIM)
{
    std::array<std::size_t, NDIM> widths;
    widths.fill(IBKernelWeightsTraits<std::invoke_result_t<const Tangential&, const double&>>::extent);
    widths[Axis] = IBKernelWeightsTraits<std::invoke_result_t<const Normal&, const double&>>::extent;
    return widths;
}

template <IBKernelEvaluatorScalar Normal, detail::IBKernelEvaluatorScalarPair<Normal> Tangential>
template <int Axis, int Direction>
auto
IBKernelEvaluatorTensorProduct<Normal, Tangential>::evaluateDirection(const double& r) const
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

template <IBKernelEvaluatorScalar Normal, detail::IBKernelEvaluatorScalarPair<Normal> Tangential>
template <int Axis>
auto
IBKernelEvaluatorTensorProduct<Normal, Tangential>::evaluate(const std::array<double, NDIM>& r) const
    requires(Axis >= 0 && Axis < NDIM)
{
    constexpr std::array<std::size_t, NDIM> widths = get_stencil_widths<Axis>();
    const auto wx = evaluateDirection<Axis, 0>(r[0]);
    const auto wy = evaluateDirection<Axis, 1>(r[1]);
    using Value = std::remove_cv_t<
        typename IBKernelWeightsTraits<std::invoke_result_t<const Normal&, const double&>>::value_type>;
    constexpr std::size_t count = widths[0] * widths[1]
#if (NDIM == 3)
                                  * widths[2]
#endif
        ;
    IBKernels::Weights<Value, count> weights;
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
                weights[i + widths[0] * (j + widths[1] * k)] =
                    static_cast<Value>(wx[i]) * static_cast<Value>(wy[j]) * static_cast<Value>(wz[k]);
#else
            weights[i + widths[0] * j] = static_cast<Value>(wx[i]) * static_cast<Value>(wy[j]);
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
