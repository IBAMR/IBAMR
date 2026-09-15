// Copyright (c) 2026 by the IBAMR developers
// This file is part of IBAMR and is distributed under the 3-clause BSD license.
#ifndef included_IBTK_Experimental_SideCoupling
#define included_IBTK_Experimental_SideCoupling

#include <ibtk/config.h>

#include <ibtk/IBKernelConcepts.h>

#include <Patch.h>
#include <SideData.h>

#include <array>
#include <cstddef>
#include <span>
#include <type_traits>

namespace IBTK::Experimental
{
/*!
 * \brief Patch-local, matrix-free coupling of depth-one side vectors.
 *
 * This experimental, uninstalled entry point copies patch geometry once and
 * retains no field or marker storage. Each axis uses the corresponding SAMRAI
 * SideData pointer with the construction-time box and ghost widths.
 * Positions have NDIM consecutive doubles per marker. indices selects marker
 * records; optional shifts has NDIM entries per selected index, in list order.
 * Indices and shifted coordinates must be valid and representable as grid indices.
 * Stencils are clipped to the allocated ghost box without renormalization.
 *
 * Interpolation overwrites selected values; spreading adds marker values divided
 * by cell volume. marker_stride permits either scalar component buffers or
 * interleaved vectors (pointer + axis, stride NDIM). Unselected values are untouched.
 * Field, marker values, positions, shifts, and indices must not overlap whenever
 * either range is written. Concurrent spreads require disjoint destinations or
 * caller-provided synchronization. Concepts do not establish these conditions.
 */
class SideCoupling
{
public:
    /*! \brief Copy geometry for an all-direction, depth-one side field. */
    SideCoupling(const SAMRAI::hier::Patch<NDIM>& patch, const SAMRAI::pdat::SideData<NDIM, double>& field);

    /*! \brief Interpolate one component with the selected coefficient arithmetic. */
    template <int Axis, class Coefficient = double, class Evaluator>
    requires IBKernelEvaluatorCartesian<Evaluator, double, Coefficient> void
    interpolateAxis(const Evaluator& evaluator,
                    const double* field,
                    std::span<const double> positions,
                    std::span<const int> indices,
                    std::span<const double> shifts,
                    double* values,
                    std::ptrdiff_t marker_stride = 1) const;

    /*! \brief Add one component to the field, spreading values rather than densities. */
    template <int Axis, class Coefficient = double, class Evaluator>
    requires IBKernelEvaluatorCartesian<Evaluator, double, Coefficient> void
    spreadAxis(const Evaluator& evaluator,
               double* field,
               std::span<const double> positions,
               std::span<const int> indices,
               std::span<const double> shifts,
               const double* values,
               std::ptrdiff_t marker_stride = 1) const;

private:
    /*! \brief Apply a component stencil with coordinate zero varying fastest. */
    template <int Axis, bool Spread, class Coefficient, class Evaluator>
    void applyAxis(const Evaluator& evaluator,
                   std::conditional_t<Spread, double*, const double*> field,
                   std::span<const double> positions,
                   std::span<const int> indices,
                   std::span<const double> shifts,
                   std::conditional_t<Spread, const double*, double*> values,
                   std::ptrdiff_t marker_stride) const;

    std::array<double, NDIM> d_dx, d_x_lower;
    std::array<int, NDIM> d_patch_lower;
    std::array<std::array<int, NDIM>, NDIM> d_lower, d_upper;
    std::array<std::array<std::ptrdiff_t, NDIM>, NDIM> d_stride;
    double d_inverse_volume = 1.0;
};
} // namespace IBTK::Experimental

#include <SideCoupling-inl.h>
#endif
