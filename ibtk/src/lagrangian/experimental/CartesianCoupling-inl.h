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

#ifndef included_IBTK_Experimental_CartesianCoupling_inl
#define included_IBTK_Experimental_CartesianCoupling_inl

#include <ibtk/config.h>

#include <tbox/Utilities.h>

#include <CartesianCoupling.h>
#include <CartesianPatchGeometry.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

namespace IBTK::Experimental
{
template <DataCentering C>
inline CartesianCoupling<C>::CartesianCoupling(const SAMRAI::hier::Patch<NDIM>& patch,
                                               const typename CartesianCentering<C>::template Data<double>& field)
{
    if (field.getDepth() < 1 || field.getBox() != patch.getBox())
    {
        TBOX_ERROR("CartesianCoupling requires a positive-depth field on the patch box.\n");
    }
    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM>> geometry = patch.getPatchGeometry();
    for (int d = 0; d < NDIM; ++d)
    {
        d_dx[d] = geometry->getDx()[d];
        d_x_lower[d] = geometry->getXLower()[d];
        d_patch_lower[d] = patch.getBox().lower()(d);
        d_inverse_volume /= d_dx[d];
    }
    for (int axis = 0; axis < s_n_arrays; ++axis)
    {
        d_has_axis[axis] = CartesianCentering<C>::template has_axis<double>(field, axis);
        if (!d_has_axis[axis])
        {
            continue;
        }
        const SAMRAI::pdat::ArrayData<NDIM, double>& array = [&]() -> const SAMRAI::pdat::ArrayData<NDIM, double>&
        {
            if constexpr (CartesianCentering<C>::is_staggered())
            {
                return field.getArrayData(axis);
            }
            else
            {
                return field.getArrayData();
            }
        }();
        const SAMRAI::hier::Box<NDIM>& box = array.getBox();
        const VectorNd offset = CartesianCentering<C>::offset(axis);
        std::ptrdiff_t stride = 1;
        for (int s = 0; s < NDIM; ++s)
        {
            const int d = C == DataCentering::FACE ? (axis + s) % NDIM : s;
            d_offset[axis][d] = offset[d];
            d_lower[axis][d] = box.lower()(s);
            d_upper[axis][d] = box.upper()(s);
            d_stride[axis][s] = stride;
            stride *= box.numberCells(s);
        }
    }
}

template <DataCentering C>
template <int Axis, class Coefficient, TensorProductMode Mode, int KernelAxis, class Evaluator>
requires IBKernelEvaluatorCartesian<Evaluator, double, Coefficient> inline void
CartesianCoupling<C>::interpolateAxis(const Evaluator& evaluator,
                                      const double* const field,
                                      const std::span<const double> positions,
                                      const std::span<const int> indices,
                                      const std::span<const double> shifts,
                                      double* const values,
                                      const std::ptrdiff_t marker_stride) const
{
    applyAxis<Axis, false, Coefficient, Mode, KernelAxis>(
        evaluator, field, positions, indices, shifts, values, marker_stride);
}

template <DataCentering C>
template <int Axis, class Coefficient, TensorProductMode Mode, int KernelAxis, class Evaluator>
requires IBKernelEvaluatorCartesian<Evaluator, double, Coefficient> inline void
CartesianCoupling<C>::spreadAxis(const Evaluator& evaluator,
                                 double* const field,
                                 const std::span<const double> positions,
                                 const std::span<const int> indices,
                                 const std::span<const double> shifts,
                                 const double* const values,
                                 const std::ptrdiff_t marker_stride) const
{
    applyAxis<Axis, true, Coefficient, Mode, KernelAxis>(
        evaluator, field, positions, indices, shifts, values, marker_stride);
}

template <DataCentering C>
template <int Axis, bool Spread, class Coefficient, TensorProductMode Mode, int KernelAxis, class Evaluator>
inline void
CartesianCoupling<C>::applyAxis(const Evaluator& evaluator,
                                const std::conditional_t<Spread, double*, const double*> field,
                                const std::span<const double> positions,
                                const std::span<const int> indices,
                                const std::span<const double> shifts,
                                const std::conditional_t<Spread, const double*, double*> values,
                                const std::ptrdiff_t marker_stride) const
{
    static_assert(Axis >= 0 && Axis < s_n_arrays);
    static_assert(KernelAxis >= 0 && KernelAxis < NDIM);
    if constexpr (C == DataCentering::SIDE)
    {
        if (!d_has_axis[Axis])
        {
            TBOX_ERROR("CartesianCoupling: requested direction is not allocated.\n");
        }
    }
    // Face arrays store their normal direction first; evaluators retain Cartesian order.
    constexpr int p0 = C == DataCentering::FACE ? Axis : 0;
    constexpr int p1 = C == DataCentering::FACE ? (Axis + 1) % NDIM : 1;
#if (NDIM == 3)
    constexpr int p2 = C == DataCentering::FACE ? (Axis + 2) % NDIM : 2;
#endif
    constexpr std::array<std::size_t, NDIM> widths = Evaluator::template get_stencil_widths<KernelAxis>();
    static_assert(std::all_of(widths.begin(),
                              widths.end(),
                              [](std::size_t n)
                              { return n <= static_cast<std::size_t>(std::numeric_limits<int>::max()); }));
    using Weights = IBKernelEvaluators::Weights<Coefficient, detail::ib_kernel_stencil_size<Evaluator, KernelAxis>()>;
    using Factors = std::tuple<IBKernelEvaluators::Weights<Coefficient, widths[0]>,
                               IBKernelEvaluators::Weights<Coefficient, widths[1]>
#if (NDIM == 3)
                               ,
                               IBKernelEvaluators::Weights<Coefficient, widths[2]>
#endif
                               >;
    constexpr bool factorized =
        Mode != TensorProductMode::EXPANDED && requires(const Evaluator& kernel, const std::array<double, NDIM>& r)
    {
        {
            kernel.template evaluateFactors<KernelAxis, Coefficient>(r)
        } -> std::same_as<Factors>;
    };
    constexpr bool contracted = factorized && Mode == TensorProductMode::CONTRACTED;
    constexpr std::array<std::size_t, NDIM> weight_stride
    {
        1, widths[0]
#if (NDIM == 3)
            ,
            widths[0] * widths[1]
#endif
    };
#if !defined(NDEBUG)
    TBOX_ASSERT(marker_stride > 0);
    TBOX_ASSERT(shifts.empty() || shifts.size() == NDIM * indices.size());
#endif
    for (std::size_t point = 0; point < indices.size(); ++point)
    {
        const int marker = indices[point];
#if !defined(NDEBUG)
        TBOX_ASSERT(marker >= 0 && NDIM * static_cast<std::size_t>(marker + 1) <= positions.size());
#endif
        std::array<double, NDIM> r;
        std::array<int, NDIM> lower, first, last;
        bool full_stencil = true;
        for (int d = 0; d < NDIM; ++d)
        {
            const double x = positions[NDIM * marker + d] + (shifts.empty() ? 0.0 : shifts[NDIM * point + d]);
            const double grid_position = (x - d_x_lower[d]) / d_dx[d] - d_offset[Axis][d];
            const int width = static_cast<int>(widths[d]);
            // Odd stencils center on the nearest grid point; even stencils straddle it.
            const int local_lower =
                static_cast<int>(std::floor(grid_position + (width % 2 ? 0.5 : 0.0))) - width / 2 + (width % 2 ? 0 : 1);
            lower[d] = d_patch_lower[d] + local_lower;
            r[d] = grid_position - local_lower;
            first[d] = std::max(0, d_lower[Axis][d] - lower[d]);
            last[d] = std::min(width, d_upper[Axis][d] - lower[d] + 1);
            full_stencil = full_stencil && first[d] == 0 && last[d] == width;
        }
        const auto weights = [&]
        {
            if constexpr (factorized)
            {
                return evaluator.template evaluateFactors<KernelAxis, Coefficient>(r);
            }
            else
            {
                return evaluator.template evaluate<KernelAxis, Weights>(r);
            }
        }();
        double value = 0.0;
        if constexpr (Spread)
        {
            value = values[marker_stride * marker] * d_inverse_volume;
        }
        // Fixed bounds let the compiler specialize complete stencils; ghosts use the same path.
        const auto apply_stencil = [&]<bool Clipped>()
        {
#if (NDIM == 3)
            for (int k = Clipped ? first[p2] : 0; k < (Clipped ? last[p2] : static_cast<int>(widths[p2])); ++k)
            {
                double plane_value = 0.0;
                if constexpr (contracted && Spread)
                {
                    plane_value = value * std::get<p2>(weights)[k];
                }
#endif
                for (int j = Clipped ? first[p1] : 0; j < (Clipped ? last[p1] : static_cast<int>(widths[p1])); ++j)
                {
                    std::ptrdiff_t offset =
                        lower[p0] - d_lower[Axis][p0] + (lower[p1] + j - d_lower[Axis][p1]) * d_stride[Axis][1];
                    std::size_t weight_offset = weight_stride[p1] * j;
#if (NDIM == 3)
                    offset += (lower[p2] + k - d_lower[Axis][p2]) * d_stride[Axis][2];
                    weight_offset += weight_stride[p2] * k;
#endif
                    double row_value = 0.0;
                    if constexpr (contracted && Spread)
                    {
#if (NDIM == 3)
                        row_value = plane_value * std::get<p1>(weights)[j];
#else
                    row_value = value * std::get<p1>(weights)[j];
#endif
                    }
                    for (int i = Clipped ? first[p0] : 0; i < (Clipped ? last[p0] : static_cast<int>(widths[p0])); ++i)
                    {
                        if constexpr (contracted)
                        {
                            if constexpr (Spread)
                            {
                                field[offset + i] += std::get<p0>(weights)[i] * row_value;
                            }
                            else
                            {
                                row_value += std::get<p0>(weights)[i] * field[offset + i];
                            }
                        }
                        else
                        {
                            const Coefficient weight = [&]() -> Coefficient
                            {
                                if constexpr (factorized)
                                {
#if (NDIM == 3)
                                    return std::get<p0>(weights)[i] * std::get<p1>(weights)[j] *
                                           std::get<p2>(weights)[k];
#else
                                return std::get<p0>(weights)[i] * std::get<p1>(weights)[j];
#endif
                                }
                                else
                                {
                                    return weights[weight_offset + weight_stride[p0] * i];
                                }
                            }();
                            if constexpr (Spread)
                            {
                                field[offset + i] += weight * value;
                            }
                            else
                            {
                                value += weight * field[offset + i];
                            }
                        }
                    }
                    if constexpr (contracted && !Spread)
                    {
#if (NDIM == 3)
                        plane_value += std::get<p1>(weights)[j] * row_value;
#else
                    value += std::get<p1>(weights)[j] * row_value;
#endif
                    }
                }
#if (NDIM == 3)
                if constexpr (contracted && !Spread)
                {
                    value += std::get<p2>(weights)[k] * plane_value;
                }
            }
#endif
        };
        if (full_stencil)
        {
            apply_stencil.template operator()<false>();
        }
        else
        {
            apply_stencil.template operator()<true>();
        }
        if constexpr (!Spread)
        {
            values[marker_stride * marker] = value;
        }
    }
}
} // namespace IBTK::Experimental
#endif
