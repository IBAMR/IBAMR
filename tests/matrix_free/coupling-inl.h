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

template <bool Spread, class Coefficient, class Evaluator>
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
            (coupling.template spreadAxis<Axis, Coefficient>(
                 evaluator, field.getPointer(Axis), positions, indices, shifts, values + Axis, NDIM),
             ...);
        }
        else
        {
            (coupling.template interpolateAxis<Axis, Coefficient>(
                 evaluator, field.getPointer(Axis), positions, indices, shifts, values + Axis, NDIM),
             ...);
        }
    }(std::make_index_sequence<NDIM>{});
}
} // namespace MatrixFreeTest
#endif
