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

#ifndef included_IBTK_Experimental_CartesianCoupling
#define included_IBTK_Experimental_CartesianCoupling

#include <ibtk/config.h>

#include <ibtk/CartesianCentering.h>
#include <ibtk/IBKernelConcepts.h>

#include <Patch.h>

#include <array>
#include <cstddef>
#include <span>
#include <type_traits>

namespace IBTK::Experimental
{
/*! \brief Compile-time choices for comparing tensor-product applications. */
enum class TensorProductMode
{
    EXPANDED,
    FACTORIZED,
    CONTRACTED
};

/*!
 * \brief Patch-local, matrix-free coupling of Cartesian field components.
 *
 * This experimental, uninstalled entry point copies patch geometry once and
 * retains no field or marker storage. C selects cell, node, side, face, or edge
 * centering. Field pointers must address the start of one depth plane with the
 * construction-time allocation bounds. For staggered data, Axis selects an
 * allocated normal (side/face) or tangent (edge) direction; cell/node use Axis=0.
 * Select depth with getPointer(depth) or getPointer(Axis, depth) on the data.
 * KernelAxis independently selects the evaluator's distinguished Cartesian
 * direction, defaulting to Axis; it does not select a depth plane.
 * Positions have NDIM consecutive doubles per marker. indices selects marker
 * records; optional shifts has NDIM entries per selected index, in list order.
 * Indices and shifted coordinates must be valid and representable as grid indices.
 * Stencils are clipped to the allocated ghost box without renormalization.
 *
 * Interpolation overwrites selected values; spreading adds marker values divided
 * by cell volume. marker_stride permits either scalar component buffers or
 * interleaved values (pointer + component, stride number_of_components).
 * Unselected values are untouched.
 * Field, marker values, positions, shifts, and indices must not overlap whenever
 * either range is written. Concurrent spreads require disjoint destinations or
 * caller-provided synchronization. Concepts do not establish these conditions.
 *
 * FACTORIZED forms coefficient products on demand in storage traversal order.
 * CONTRACTED reuses plane/row products and sums with double field arithmetic,
 * using the selected precision for the one-dimensional factors. Evaluators without an
 * owning evaluateFactors() result use EXPANDED in every mode.
 */
template <DataCentering C>
class CartesianCoupling
{
public:
    /*! \brief Copy geometry and allocation bounds for the field's available directions. */
    CartesianCoupling(const SAMRAI::hier::Patch<NDIM>& patch,
                      const typename CartesianCentering<C>::template Data<double>& field);

    /*! \brief Interpolate one component with the selected coefficient arithmetic. */
    template <int Axis,
              class Coefficient = double,
              TensorProductMode Mode = TensorProductMode::CONTRACTED,
              int KernelAxis = Axis,
              class Evaluator>
    requires IBKernelEvaluatorCartesian<Evaluator, double, Coefficient> void
    interpolateAxis(const Evaluator& evaluator,
                    const double* field,
                    std::span<const double> positions,
                    std::span<const int> indices,
                    std::span<const double> shifts,
                    double* values,
                    std::ptrdiff_t marker_stride = 1) const;

    /*! \brief Add one component to the field, spreading values rather than densities. */
    template <int Axis,
              class Coefficient = double,
              TensorProductMode Mode = TensorProductMode::CONTRACTED,
              int KernelAxis = Axis,
              class Evaluator>
    requires IBKernelEvaluatorCartesian<Evaluator, double, Coefficient> void
    spreadAxis(const Evaluator& evaluator,
               double* field,
               std::span<const double> positions,
               std::span<const int> indices,
               std::span<const double> shifts,
               const double* values,
               std::ptrdiff_t marker_stride = 1) const;

private:
    /*! \brief Apply a component stencil with the contiguous storage direction varying fastest. */
    template <int Axis, bool Spread, class Coefficient, TensorProductMode Mode, int KernelAxis, class Evaluator>
    void applyAxis(const Evaluator& evaluator,
                   std::conditional_t<Spread, double*, const double*> field,
                   std::span<const double> positions,
                   std::span<const int> indices,
                   std::span<const double> shifts,
                   std::conditional_t<Spread, const double*, double*> values,
                   std::ptrdiff_t marker_stride) const;

    static constexpr int s_n_arrays = CartesianCentering<C>::is_staggered() ? NDIM : 1;
    std::array<bool, s_n_arrays> d_has_axis{};
    std::array<std::array<double, NDIM>, s_n_arrays> d_offset{};
    std::array<double, NDIM> d_dx, d_x_lower;
    std::array<int, NDIM> d_patch_lower;
    std::array<std::array<int, NDIM>, s_n_arrays> d_lower{}, d_upper{};
    std::array<std::array<std::ptrdiff_t, NDIM>, s_n_arrays> d_stride{};
    double d_inverse_volume = 1.0;
};
} // namespace IBTK::Experimental

#include <CartesianCoupling-inl.h>
#endif
