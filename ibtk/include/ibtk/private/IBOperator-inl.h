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

#ifndef included_IBTK_IBOperator_inl
#define included_IBTK_IBOperator_inl
#include <ibtk/config.h>

#include <ibtk/IBOperator.h>
#include <ibtk/private/CartesianCoupling.h>

#include <algorithm>
#include <utility>

namespace IBTK
{
namespace detail
{
/*! \brief Select a spatial axis once per component batch. */
template <int Count, class Function>
void
dispatch_coupling_axis(const int axis, const Function& function)
{
    switch (axis)
    {
    case 0:
        function.template operator()<0>();
        break;
    case 1:
        if constexpr (Count > 1)
        {
            function.template operator()<1>();
            break;
        }
        else
        {
            TBOX_ERROR("IBOperator: invalid field axis.\n");
        }
#if NDIM == 3
    case 2:
        if constexpr (Count > 2)
        {
            function.template operator()<2>();
            break;
        }
        else
        {
            TBOX_ERROR("IBOperator: invalid field axis.\n");
        }
#endif
    default:
        TBOX_ERROR("IBOperator: invalid axis.\n");
    }
}

/*! \brief Apply the shared numerical core with a borrowed concrete evaluator. */
template <DataCentering C, bool Spread, class Evaluator>
void
apply_cartesian_coupling(const Evaluator& evaluator,
                         const SAMRAI::hier::Patch<NDIM>& patch,
                         std::conditional_t<Spread,
                                            typename CartesianCentering<C>::template Data<double>&,
                                            const typename CartesianCentering<C>::template Data<double>&> field,
                         const int field_axis,
                         const int depth,
                         const int kernel_axis,
                         const std::span<const double> positions,
                         const std::span<const int> indices,
                         const std::span<const double> shifts,
                         const std::conditional_t<Spread, const double*, double*> values,
                         const std::ptrdiff_t marker_stride,
                         const SAMRAI::geom::CartesianGridGeometry<NDIM>* const grid_geometry)
{
    if (depth < 0 || depth >= field.getDepth() || marker_stride <= 0 || positions.size() % NDIM != 0 ||
        (!shifts.empty() && shifts.size() != NDIM * indices.size()))
    {
        TBOX_ERROR("IBOperator: invalid depth, stride, or marker array shape.\n");
    }
    for (const int index : indices)
    {
        if (index < 0 || static_cast<std::size_t>(index) >= positions.size() / NDIM)
        {
            TBOX_ERROR("IBOperator: marker index outside coordinate storage.\n");
        }
    }
    const CartesianCoupling<C> coupling(patch, field, grid_geometry);
    dispatch_coupling_axis<CartesianCentering<C>::is_staggered() ? NDIM : 1>(
        field_axis,
        [&]<int Axis>()
        {
            if (!CartesianCentering<C>::template has_axis<double>(field, Axis))
            {
                TBOX_ERROR("IBOperator: requested direction is not allocated.\n");
            }
            const auto pointer = [&]()
            {
                if constexpr (CartesianCentering<C>::is_staggered())
                {
                    return field.getPointer(Axis, depth);
                }
                else
                {
                    return field.getPointer(depth);
                }
            }();
            dispatch_coupling_axis<NDIM>(
                kernel_axis,
                [&]<int KernelAxis>()
                {
                    if constexpr (Spread)
                    {
                        coupling.template spreadAxis<Axis, double, TensorProductMode::CONTRACTED, KernelAxis>(
                            evaluator, pointer, positions, indices, shifts, values, marker_stride);
                    }
                    else
                    {
                        coupling.template interpolateAxis<Axis, double, TensorProductMode::CONTRACTED, KernelAxis>(
                            evaluator, pointer, positions, indices, shifts, values, marker_stride);
                    }
                });
        });
}
} // namespace detail

template <class Evaluator>
class IBOperator::Model final : public IBOperator::Concept
{
public:
    template <class Argument>
    requires std::constructible_from<Evaluator, Argument&&> explicit Model(Argument&& evaluator)
        : d_evaluator(std::forward<Argument>(evaluator))
    {
    }
    int getMinimumGhostWidth() const override
    {
        std::size_t width = 0;
        for (const std::array<std::size_t, NDIM>& widths : {
                 Evaluator::template get_stencil_widths<0>(), Evaluator::template get_stencil_widths<1>()
#if (NDIM == 3)
                                                                  ,
                     Evaluator::template get_stencil_widths<2>()
#endif
             })
        {
            width = std::max(width, *std::max_element(widths.begin(), widths.end()));
        }
        return static_cast<int>(width / 2 + 1);
    }

    void interpolate(DataCentering centering,
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
                     const SAMRAI::geom::CartesianGridGeometry<NDIM>* grid_geometry) const override
    {
        dispatch_data_centering(
            centering,
            [&]<DataCentering C>()
            {
                detail::apply_cartesian_coupling<C, false>(
                    d_evaluator,
                    patch,
                    static_cast<const typename CartesianCentering<C>::template Data<double>&>(field),
                    field_axis,
                    depth,
                    kernel_axis,
                    positions,
                    indices,
                    shifts,
                    values,
                    marker_stride,
                    grid_geometry);
            });
    }
    void spread(DataCentering centering,
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
                const SAMRAI::geom::CartesianGridGeometry<NDIM>* grid_geometry) const override
    {
        dispatch_data_centering(centering,
                                [&]<DataCentering C>()
                                {
                                    detail::apply_cartesian_coupling<C, true>(
                                        d_evaluator,
                                        patch,
                                        static_cast<typename CartesianCentering<C>::template Data<double>&>(field),
                                        field_axis,
                                        depth,
                                        kernel_axis,
                                        positions,
                                        indices,
                                        shifts,
                                        values,
                                        marker_stride,
                                        grid_geometry);
                                });
    }

private:
    const Evaluator d_evaluator;
};

template <class Evaluator>
requires(IBKernelEvaluatorCartesian<std::remove_cvref_t<Evaluator>, double, double>&&
             std::constructible_from<std::remove_cvref_t<Evaluator>, Evaluator&&>)
    IBOperator::IBOperator(Evaluator&& evaluator)
    : d_operations(std::make_shared<Model<std::remove_cvref_t<Evaluator>>>(std::forward<Evaluator>(evaluator)))
{
}

inline int
IBOperator::getMinimumGhostWidth() const
{
    return d_operations->getMinimumGhostWidth();
}
template <DataCentering C>
inline void
IBOperator::interpolate(const SAMRAI::hier::Patch<NDIM>& patch,
                        const typename CartesianCentering<C>::template Data<double>& field,
                        int field_axis,
                        int depth,
                        int kernel_axis,
                        std::span<const double> positions,
                        std::span<const int> indices,
                        std::span<const double> shifts,
                        double* values,
                        std::ptrdiff_t marker_stride,
                        const SAMRAI::geom::CartesianGridGeometry<NDIM>* grid_geometry) const
{
    d_operations->interpolate(C,
                              patch,
                              field,
                              field_axis,
                              depth,
                              kernel_axis,
                              positions,
                              indices,
                              shifts,
                              values,
                              marker_stride,
                              grid_geometry);
}
template <DataCentering C>
inline void
IBOperator::spread(const SAMRAI::hier::Patch<NDIM>& patch,
                   typename CartesianCentering<C>::template Data<double>& field,
                   int field_axis,
                   int depth,
                   int kernel_axis,
                   std::span<const double> positions,
                   std::span<const int> indices,
                   std::span<const double> shifts,
                   const double* values,
                   std::ptrdiff_t marker_stride,
                   const SAMRAI::geom::CartesianGridGeometry<NDIM>* grid_geometry) const
{
    d_operations->spread(C,
                         patch,
                         field,
                         field_axis,
                         depth,
                         kernel_axis,
                         positions,
                         indices,
                         shifts,
                         values,
                         marker_stride,
                         grid_geometry);
}
} // namespace IBTK
#endif
