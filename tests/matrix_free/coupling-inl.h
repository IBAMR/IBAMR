// Copyright (c) 2026 by the IBAMR developers
// This file is part of IBAMR and is distributed under the 3-clause BSD license.
#ifndef included_matrix_free_coupling_inl
#define included_matrix_free_coupling_inl

#include <ibtk/config.h>

#include <cmath>
#include <utility>

#include "coupling.h"

namespace MatrixFreeTest
{
template <class Visitor>
void
for_each_bspline(const Visitor& visit)
{
    [&]<std::size_t... K>(std::index_sequence<K...>)
    {
        (visit("BSPLINE_" + std::to_string(K + 2),
               IBTK::IBKernelEvaluatorTensorProduct{ IBTK::IBKernels::BSpline<K + 2>{} }),
         ...);
        (visit("COMPOSITE_BSPLINE_" + std::to_string(K + 2) + std::to_string(K + 1),
               IBTK::IBKernelEvaluatorTensorProduct{ IBTK::IBKernels::BSpline<K + 2>{},
                                                     IBTK::IBKernels::BSpline<K + 1>{} }),
         ...);
        (visit("COMPOSITE_BSPLINE_" + std::to_string(K + 1) + std::to_string(K + 2),
               IBTK::IBKernelEvaluatorTensorProduct{ IBTK::IBKernels::BSpline<K + 1>{},
                                                     IBTK::IBKernels::BSpline<K + 2>{} }),
         ...);
    }(std::make_index_sequence<5>{});
}

constexpr std::size_t
CosineKernel::get_stencil_width()
{
    return 4;
}

template <class Output, class Input>
Output
CosineKernel::evaluate(const Input r) const
{
    static_assert(IBTK::IBKernelWeightsTraits<Output>::extent == get_stencil_width());
    using Coefficient = typename IBTK::IBKernelWeightsTraits<Output>::value_type;
    const Coefficient x = r;
    const Coefficient half_pi = std::acos(Coefficient{ -1 }) / 2;
    Output weights{};
    for (std::size_t i = 0; i < get_stencil_width(); ++i)
    {
        weights[i] = Coefficient{ 0.25 } * (1 + std::cos(half_pi * (x - i)));
    }
    return weights;
}

template <int Axis>
constexpr std::array<std::size_t, NDIM>
CartesianCosineKernel::get_stencil_widths()
{
    std::array<std::size_t, NDIM> widths;
    widths.fill(4);
    return widths;
}

template <int Axis, class Output, class Input>
Output
CartesianCosineKernel::evaluate(const std::array<Input, NDIM>& r) const
{
    return IBTK::IBKernelEvaluatorTensorProduct{ CosineKernel{} }.template evaluate<Axis, Output>(r);
}

template <bool Spread, class Coefficient, IBTK::Experimental::TensorProductMode Mode, class Evaluator>
void
couple(const Evaluator& evaluator,
       const SAMRAI::hier::Patch<NDIM>& patch,
       SAMRAI::pdat::SideData<NDIM, double>& field,
       const std::span<const double> positions,
       const std::span<const int> indices,
       const std::span<const double> shifts,
       const std::conditional_t<Spread, const double*, double*> values)
{
    const IBTK::Experimental::SideCoupling coupling(patch, field);
    [&]<std::size_t... Axis>(std::index_sequence<Axis...>)
    {
        if constexpr (Spread)
        {
            (coupling.template spreadAxis<Axis, Coefficient, Mode>(
                 evaluator, field.getPointer(Axis), positions, indices, shifts, values + Axis, NDIM),
             ...);
        }
        else
        {
            (coupling.template interpolateAxis<Axis, Coefficient, Mode>(
                 evaluator, field.getPointer(Axis), positions, indices, shifts, values + Axis, NDIM),
             ...);
        }
    }(std::make_index_sequence<NDIM>{});
}
} // namespace MatrixFreeTest
#endif
