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

#ifndef included_IBTK_IBOperator
#define included_IBTK_IBOperator
#include <ibtk/config.h>

#include <ibtk/CartesianCentering.h>
#include <ibtk/IBKernelConcepts.h>
#include <ibtk/IBKernelTensorProduct.h>

#include <CartesianGridGeometry.h>
#include <Patch.h>
#include <PatchData.h>

#include <concepts>
#include <cstddef>
#include <memory>
#include <span>
#include <type_traits>

namespace IBTK
{
/*!
 * \brief Immutable, evaluator-owning matrix-free IB coupling on Cartesian patches.
 *
 * Copies share a const evaluator; the handle retains no patch or marker storage.
 * Field axis, depth, kernel orientation, and marker stride are independent.
 * Positions contain NDIM coordinates per record; indices select those records.
 * Optional shifts contain NDIM physical displacements per selection-list entry.
 * Values select a component by values[marker_stride * indices[l]]. The caller
 * supplies sufficient storage and finite, representable coordinates. Writable
 * outputs must not overlap inputs. Concurrent spreads need disjoint destinations
 * or caller synchronization.
 *
 * Gather overwrites selected values (the last repeated index wins). Spread adds
 * every selected value divided by cell volume. Quadrature weights, ghost fills,
 * boundary accumulation, and hierarchy transfers belong to the caller. Stencils
 * clip to allocated storage without renormalization. An empty selection accesses
 * no marker values; an empty stencil gathers zero and spreads nothing.
 *
 * Supply grid_geometry for hierarchy calls and comparison with assembled matrices.
 * Its single-box domain and the patch refinement ratio define the common lattice
 * origin and nearest-endpoint cell mapping. Without it, the patch bounds define
 * that lattice; floating-point placement can differ between these conventions.
 * Negative SAMRAI ratios denote coarsening and multiply the reference spacing
 * by their magnitude. Assembled-matrix parity is specified for refined levels;
 * assembled operators may impose additional ratio restrictions.
 * The geometry is borrowed only during the call.
 *
 * Cell/node fields require field_axis zero. Side/face axes are normal directions;
 * edge axes are tangent directions. The selected side direction must be allocated.
 * Odd-width stencils center on the nearest grid point. At even-width exact ties,
 * zero geometric offsets select the right window and half-cell offsets the left.
 */
class IBOperator
{
public:
    /*! \brief Own a custom evaluator by copying an lvalue or moving an rvalue. */
    template <class Evaluator>
    requires(IBKernelEvaluatorCartesian<std::remove_cvref_t<Evaluator>, double, double>&&
                 std::constructible_from<std::remove_cvref_t<Evaluator>,
                                         Evaluator&&>) explicit IBOperator(Evaluator&& evaluator);

    /*! \brief Own the built-in evaluator for kernel; reject unsupported kernels. */
    explicit IBOperator(const IBKernelTensorProduct& kernel);

    /*! \brief Return whether named construction supports this kernel. */
    static bool is_built_in(const IBKernelTensorProduct& kernel);

    /*! \brief Return max_stencil_width / 2 + 1, including one-cell marker movement. */
    int getMinimumGhostWidth() const;

    /*! \brief Overwrite one marker component from a field depth plane. */
    template <DataCentering C>
    void interpolate(const SAMRAI::hier::Patch<NDIM>& patch,
                     const typename CartesianCentering<C>::template Data<double>& field,
                     int field_axis,
                     int depth,
                     int kernel_axis,
                     std::span<const double> positions,
                     std::span<const int> indices,
                     std::span<const double> shifts,
                     double* values,
                     std::ptrdiff_t marker_stride = 1,
                     const SAMRAI::geom::CartesianGridGeometry<NDIM>* grid_geometry = nullptr) const;

    /*! \brief Add one marker component to a field depth plane. */
    template <DataCentering C>
    void spread(const SAMRAI::hier::Patch<NDIM>& patch,
                typename CartesianCentering<C>::template Data<double>& field,
                int field_axis,
                int depth,
                int kernel_axis,
                std::span<const double> positions,
                std::span<const int> indices,
                std::span<const double> shifts,
                const double* values,
                std::ptrdiff_t marker_stride = 1,
                const SAMRAI::geom::CartesianGridGeometry<NDIM>* grid_geometry = nullptr) const;

private:
    /*! \brief Dispatch component batches to immutable evaluator storage. */
    class Concept
    {
    public:
        virtual ~Concept() = default;
        /*! \copydoc IBOperator::getMinimumGhostWidth */
        virtual int getMinimumGhostWidth() const = 0;
        /*! \brief Interpolate a batch whose concrete centering has been checked by the typed entry point. */
        virtual void interpolate(DataCentering centering,
                                 const SAMRAI::hier::Patch<NDIM>& patch,
                                 const SAMRAI::hier::PatchData<NDIM>& field,
                                 int field_axis,
                                 int depth,
                                 int kernel_axis,
                                 std::span<const double> positions,
                                 std::span<const int> indices,
                                 std::span<const double> shifts,
                                 double* values,
                                 std::ptrdiff_t marker_stride,
                                 const SAMRAI::geom::CartesianGridGeometry<NDIM>* grid_geometry) const = 0;
        /*! \brief Spread a batch whose concrete centering has been checked by the typed entry point. */
        virtual void spread(DataCentering centering,
                            const SAMRAI::hier::Patch<NDIM>& patch,
                            SAMRAI::hier::PatchData<NDIM>& field,
                            int field_axis,
                            int depth,
                            int kernel_axis,
                            std::span<const double> positions,
                            std::span<const int> indices,
                            std::span<const double> shifts,
                            const double* values,
                            std::ptrdiff_t marker_stride,
                            const SAMRAI::geom::CartesianGridGeometry<NDIM>* grid_geometry) const = 0;
    };
    template <class Evaluator>
    class Model;
    std::shared_ptr<const Concept> d_operations;
};
} // namespace IBTK
#include <ibtk/private/IBOperator-inl.h>
#endif
