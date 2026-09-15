// Copyright (c) 2026 by the IBAMR developers
// This file is part of IBAMR and is distributed under the 3-clause BSD license.
#ifndef included_IBTK_Experimental_SideCoupling_inl
#define included_IBTK_Experimental_SideCoupling_inl

#include <ibtk/config.h>

#include <tbox/Utilities.h>

#include <CartesianPatchGeometry.h>
#include <SideCoupling.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace IBTK::Experimental
{
inline SideCoupling::SideCoupling(const SAMRAI::hier::Patch<NDIM>& patch,
                                  const SAMRAI::pdat::SideData<NDIM, double>& field)
{
    if (field.getDepth() != 1 || field.getDirectionVector() != SAMRAI::hier::IntVector<NDIM>(1) ||
        field.getBox() != patch.getBox())
    {
        TBOX_ERROR("SideCoupling requires an all-direction, depth-one field on the patch box.\n");
    }
    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM>> geometry = patch.getPatchGeometry();
    for (int d = 0; d < NDIM; ++d)
    {
        d_dx[d] = geometry->getDx()[d];
        d_x_lower[d] = geometry->getXLower()[d];
        d_patch_lower[d] = patch.getBox().lower()(d);
        d_inverse_volume /= d_dx[d];
    }
    for (int axis = 0; axis < NDIM; ++axis)
    {
        const SAMRAI::hier::Box<NDIM>& box = field.getArrayData(axis).getBox();
        std::ptrdiff_t stride = 1;
        for (int d = 0; d < NDIM; ++d)
        {
            d_lower[axis][d] = box.lower()(d);
            d_upper[axis][d] = box.upper()(d);
            d_stride[axis][d] = stride;
            stride *= box.numberCells(d);
        }
    }
}

template <int Axis, class Coefficient, class Evaluator>
requires IBKernelEvaluatorCartesian<Evaluator, double, Coefficient> inline void
SideCoupling::interpolateAxis(const Evaluator& evaluator,
                              const double* const field,
                              const std::span<const double> positions,
                              const std::span<const int> indices,
                              const std::span<const double> shifts,
                              double* const values,
                              const std::ptrdiff_t marker_stride) const
{
    applyAxis<Axis, false, Coefficient>(evaluator, field, positions, indices, shifts, values, marker_stride);
}

template <int Axis, class Coefficient, class Evaluator>
requires IBKernelEvaluatorCartesian<Evaluator, double, Coefficient> inline void
SideCoupling::spreadAxis(const Evaluator& evaluator,
                         double* const field,
                         const std::span<const double> positions,
                         const std::span<const int> indices,
                         const std::span<const double> shifts,
                         const double* const values,
                         const std::ptrdiff_t marker_stride) const
{
    applyAxis<Axis, true, Coefficient>(evaluator, field, positions, indices, shifts, values, marker_stride);
}

template <int Axis, bool Spread, class Coefficient, class Evaluator>
inline void
SideCoupling::applyAxis(const Evaluator& evaluator,
                        const std::conditional_t<Spread, double*, const double*> field,
                        const std::span<const double> positions,
                        const std::span<const int> indices,
                        const std::span<const double> shifts,
                        const std::conditional_t<Spread, const double*, double*> values,
                        const std::ptrdiff_t marker_stride) const
{
    static_assert(Axis >= 0 && Axis < NDIM);
    constexpr std::array<std::size_t, NDIM> widths = Evaluator::template get_stencil_widths<Axis>();
    static_assert(std::all_of(widths.begin(),
                              widths.end(),
                              [](std::size_t n)
                              { return n <= static_cast<std::size_t>(std::numeric_limits<int>::max()); }));
    using Weights = IBKernels::Weights<Coefficient, detail::ib_kernel_stencil_size<Evaluator, Axis>()>;
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
            const double grid_position = (x - d_x_lower[d]) / d_dx[d] - (d == Axis ? 0.0 : 0.5);
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
        const Weights weights = evaluator.template evaluate<Axis, Weights>(r);
        double value = 0.0;
        if constexpr (Spread)
        {
            value = values[marker_stride * marker] * d_inverse_volume;
        }
        // Fixed bounds let the compiler specialize complete stencils; ghosts use the same path.
        const auto apply_stencil = [&]<bool Clipped>()
        {
#if (NDIM == 3)
            for (int k = Clipped ? first[2] : 0; k < (Clipped ? last[2] : static_cast<int>(widths[2])); ++k)
            {
#endif
                for (int j = Clipped ? first[1] : 0; j < (Clipped ? last[1] : static_cast<int>(widths[1])); ++j)
                {
                    std::ptrdiff_t offset =
                        lower[0] - d_lower[Axis][0] + (lower[1] + j - d_lower[Axis][1]) * d_stride[Axis][1];
                    std::size_t weight_offset = widths[0] * j;
#if (NDIM == 3)
                    offset += (lower[2] + k - d_lower[Axis][2]) * d_stride[Axis][2];
                    weight_offset += widths[0] * widths[1] * k;
#endif
                    for (int i = Clipped ? first[0] : 0; i < (Clipped ? last[0] : static_cast<int>(widths[0])); ++i)
                    {
                        if constexpr (Spread)
                        {
                            field[offset + i] += weights[weight_offset + i] * value;
                        }
                        else
                        {
                            value += weights[weight_offset + i] * field[offset + i];
                        }
                    }
                }
#if (NDIM == 3)
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
