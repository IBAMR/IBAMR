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

#include <ibtk/CartesianCentering.h>
#include <ibtk/IBOperator.h>

#include <tbox/Array.h>

#include <CartesianCoupling.h>
#include <CartesianPatchGeometry.h>
#include <CellIterator.h>
#include <EdgeIterator.h>
#include <FaceIterator.h>
#include <NodeIterator.h>
#include <PatchDescriptor.h>
#include <SideIterator.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <span>
#include <string>
#include <utility>

#include "cartesian_coupling.h"
#include "coupling.h"
#include "fixture.h"

#include <ibtk/app_namespaces.h>

namespace
{
using IBTK::CartesianCentering;
using IBTK::DataCentering;
using Mode = IBTK::Experimental::TensorProductMode;

/*! \brief Construct a noncubic patch with shifted indices and unequal grid spacings. */
Pointer<Patch<NDIM>>
make_cartesian_patch()
{
    hier::Index<NDIM> lower, upper;
    std::array<double, NDIM> dx, x_lower, x_upper;
    tbox::Array<tbox::Array<bool>> boundaries(NDIM);
    for (int d = 0; d < NDIM; ++d)
    {
        const int cells = 7 + 2 * d;
        lower(d) = -5 + 7 * d;
        upper(d) = lower(d) + cells - 1;
        dx[d] = 0.125 + 0.03125 * d;
        x_lower[d] = -0.75 + 0.125 * d;
        x_upper[d] = x_lower[d] + cells * dx[d];
        boundaries[d].resizeArray(2);
        boundaries[d][0] = boundaries[d][1] = false;
    }
    Pointer<Patch<NDIM>> patch = new Patch<NDIM>(Box<NDIM>(lower, upper), new PatchDescriptor<NDIM>());
    patch->setPatchGeometry(new CartesianPatchGeometry<NDIM>(
        IntVector<NDIM>(1), boundaries, boundaries, dx.data(), x_lower.data(), x_upper.data()));
    return patch;
}

/*! \brief Visit SAMRAI indices and independently reconstruct their Cartesian coordinates. */
template <DataCentering C, class Function>
void
visit_grid(const Box<NDIM>& ghost_box, const int axis, const Function& function)
{
    if constexpr (C == DataCentering::CELL)
    {
        for (CellIterator<NDIM> i(ghost_box); i; i++)
        {
            function(i(), static_cast<const hier::Index<NDIM>&>(i()));
        }
    }
    else if constexpr (C == DataCentering::NODE)
    {
        for (NodeIterator<NDIM> i(ghost_box); i; i++)
        {
            function(i(), static_cast<const hier::Index<NDIM>&>(i()));
        }
    }
    else if constexpr (C == DataCentering::SIDE)
    {
        for (SideIterator<NDIM> i(ghost_box, axis); i; i++)
        {
            function(i(), static_cast<const hier::Index<NDIM>&>(i()));
        }
    }
    else if constexpr (C == DataCentering::FACE)
    {
        for (FaceIterator<NDIM> i(ghost_box, axis); i; i++)
        {
            hier::Index<NDIM> cartesian;
            for (int d = 0; d < NDIM; ++d)
            {
                cartesian((axis + d) % NDIM) = i()(d);
            }
            function(i(), cartesian);
        }
    }
    else
    {
        for (EdgeIterator<NDIM> i(ghost_box, axis); i; i++)
        {
            function(i(), static_cast<const hier::Index<NDIM>&>(i()));
        }
    }
}

/*! \brief Return the physical grid offset without using the coupling's centering helper. */
template <DataCentering C>
double
grid_offset(const int axis, const int direction)
{
    if constexpr (C == DataCentering::CELL)
    {
        return 0.5;
    }
    else if constexpr (C == DataCentering::NODE)
    {
        return 0.0;
    }
    else if constexpr (C == DataCentering::EDGE)
    {
        return axis == direction ? 0.5 : 0.0;
    }
    else
    {
        return axis == direction ? 0.0 : 0.5;
    }
}

/*! \brief Use a nonseparable field to expose index and component permutations. */
double
field_value(const hier::Index<NDIM>& index, const int axis)
{
    double value = 0.5 + 0.3 * axis;
    for (int d = 0; d < NDIM; ++d)
    {
        value += std::sin(0.11 * (d + 1) * index(d));
    }
    return value + 0.007 * index(0) * index(NDIM - 1);
}

/*! \brief Compare one chosen array and depth against dense scalar gather and spread. */
template <DataCentering C, int Axis, int KernelAxis, Mode Application = Mode::CONTRACTED, class Evaluator>
std::array<double, 5>
check_axis(const std::string& name, const Evaluator& evaluator, const bool partial_side = false)
{
    using Data = typename CartesianCentering<C>::template Data<double>;
    constexpr bool staggered = C == DataCentering::SIDE || C == DataCentering::FACE || C == DataCentering::EDGE;
    const Pointer<Patch<NDIM>> patch = make_cartesian_patch();
    const Pointer<CartesianPatchGeometry<NDIM>> geometry = patch->getPatchGeometry();
    IntVector<NDIM> ghosts, directions(1);
    for (int d = 0; d < NDIM; ++d)
    {
        ghosts(d) = d;
    }
    if (partial_side)
    {
        directions = IntVector<NDIM>(0);
        directions(Axis) = 1;
    }
    const auto make_data = [&]() -> std::unique_ptr<Data>
    {
        if constexpr (C == DataCentering::SIDE)
        {
            return std::make_unique<Data>(patch->getBox(), 3, ghosts, directions);
        }
        else
        {
            return std::make_unique<Data>(patch->getBox(), 3, ghosts);
        }
    };
    const std::unique_ptr<Data> field = make_data(), spread = make_data();
    constexpr double initial_spread = 0.125;
    spread->fillAll(initial_spread);
    for (int axis = 0; axis < (staggered ? NDIM : 1); ++axis)
    {
        if constexpr (C == DataCentering::SIDE)
        {
            if (!directions(axis))
            {
                continue;
            }
        }
        visit_grid<C>(field->getGhostBox(),
                      axis,
                      [&](const auto& index, const hier::Index<NDIM>& cartesian)
                      {
                          (*field)(index, 0) = -4.25;
                          (*field)(index, 1) = field_value(cartesian, axis);
                          (*field)(index, 2) = -9.5;
                      });
    }

    // Markers zero and two are unselected. Marker three has an empty clipped stencil.
    const std::array<int, 4> indices{ 4, 1, 5, 3 };
    std::array<double, 6 * NDIM> positions;
    std::array<double, 4 * NDIM> shifts;
    for (int d = 0; d < NDIM; ++d)
    {
        for (int marker = 0; marker < 6; ++marker)
        {
            double coordinate = 3.125 + 0.1875 * d;
            if (marker == 1)
            {
                coordinate = -0.125;
            }
            else if (marker == 3)
            {
                coordinate = -12.0;
            }
            else if (marker == 5)
            {
                coordinate = patch->getBox().numberCells(d) - 0.125;
            }
            positions[NDIM * marker + d] = geometry->getXLower()[d] + coordinate * geometry->getDx()[d];
        }
        for (std::size_t point = 0; point < indices.size(); ++point)
        {
            shifts[NDIM * point + d] = (0.0625 * static_cast<double>(point + 1) - 0.03125 * d) * geometry->getDx()[d];
        }
    }
    const std::array<double, 6> force{ 0.3, -0.5, 0.7, 1.1, 0.9, -0.4 };
    std::array<double, 6> values, reference;
    values.fill(-17.0);
    reference.fill(-17.0);
    for (const int marker : indices)
    {
        reference[marker] = 0.0;
    }
    const auto selected_pointer = [](Data& data) -> double*
    {
        if constexpr (staggered)
        {
            return data.getPointer(Axis, 1);
        }
        else
        {
            return data.getPointer(1);
        }
    };
    const IBTK::Experimental::CartesianCoupling<C> coupling(*patch, *field);
    coupling.template interpolateAxis<Axis, double, Application, KernelAxis>(
        evaluator, selected_pointer(*field), positions, indices, shifts, values.data(), 1);
    coupling.template spreadAxis<Axis, double, Application, KernelAxis>(
        evaluator, selected_pointer(*spread), positions, indices, shifts, force.data(), 1);
    const IBTK::IBOperator op = name == "COSINE_4" ? IBTK::IBOperator(MatrixFreeTest::CartesianCosineKernel{}) :
                                                     IBTK::IBOperator(IBTK::IBKernelTensorProduct(name));
    const std::unique_ptr<Data> handle_spread = make_data();
    handle_spread->fillAll(initial_spread);
    std::array<double, 18> handle_values, handle_force;
    handle_values.fill(-17.0);
    handle_force.fill(-99.0);
    for (int k = 0; k < 6; ++k)
    {
        handle_force[3 * k + 1] = force[k];
    }
    op.template interpolate<C>(
        *patch, *field, Axis, 1, KernelAxis, positions, indices, shifts, handle_values.data() + 1, 3);
    op.template spread<C>(
        *patch, *handle_spread, Axis, 1, KernelAxis, positions, indices, shifts, handle_force.data() + 1, 3);
    op.template interpolate<C>(*patch, *field, Axis, 1, KernelAxis, {}, {}, {}, nullptr);
    op.template spread<C>(*patch, *handle_spread, Axis, 1, KernelAxis, {}, {}, {}, nullptr);
    for (int k = 0; k < 6; ++k)
    {
        TBOX_ASSERT(std::isfinite(handle_values[3 * k + 1]));
        TBOX_ASSERT(std::abs(handle_values[3 * k + 1] - values[k]) < 1.0e-12);
        TBOX_ASSERT(handle_values[3 * k] == -17.0 && handle_values[3 * k + 2] == -17.0);
        values[k] = handle_values[3 * k + 1];
    }
    visit_grid<C>(field->getGhostBox(),
                  Axis,
                  [&](const auto& index, const hier::Index<NDIM>&)
                  {
                      TBOX_ASSERT(std::isfinite((*handle_spread)(index, 1)));
                      TBOX_ASSERT(std::abs((*handle_spread)(index, 1) - (*spread)(index, 1)) < 1.0e-10);
                      (*spread)(index, 1) = (*handle_spread)(index, 1);
                  });
    double volume = 1.0;
    for (int d = 0; d < NDIM; ++d)
    {
        volume *= geometry->getDx()[d];
    }
    std::array<double, 5> result{};
    double grid_pairing = 0.0, marker_pairing = 0.0;
    visit_grid<C>(
        field->getGhostBox(),
        Axis,
        [&](const auto& index, const hier::Index<NDIM>& cartesian)
        {
            double expected_spread = 0.0;
            for (std::size_t point = 0; point < indices.size(); ++point)
            {
                const int marker = indices[point];
                double weight = 1.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    const double grid_x =
                        geometry->getXLower()[d] +
                        (cartesian(d) - patch->getBox().lower()(d) + grid_offset<C>(Axis, d)) * geometry->getDx()[d];
                    weight *= MatrixFreeTest::reference_weight(
                        name,
                        KernelAxis,
                        d,
                        (positions[NDIM * marker + d] + shifts[NDIM * point + d] - grid_x) / geometry->getDx()[d]);
                }
                reference[marker] += weight * (*field)(index, 1);
                expected_spread += weight * force[marker];
            }
            const double actual_spread = volume * ((*spread)(index, 1) - initial_spread);
            if (!std::isfinite(actual_spread) || !std::isfinite(expected_spread))
            {
                TBOX_ERROR("Nonfinite Cartesian spread for " << name << '\n');
            }
            result[1] = std::max(result[1], std::abs(actual_spread - expected_spread));
            grid_pairing += actual_spread * (*field)(index, 1);
            result[4] += actual_spread * actual_spread;
        });
    for (int marker = 0; marker < 6; ++marker)
    {
        if (!std::isfinite(values[marker]) || !std::isfinite(reference[marker]))
        {
            TBOX_ERROR("Nonfinite Cartesian gather for " << name << '\n');
        }
        result[0] = std::max(result[0], std::abs(values[marker] - reference[marker]));
    }
    if (values[0] != -17.0 || values[2] != -17.0 || values[3] != 0.0)
    {
        TBOX_ERROR("Cartesian gather changed an unselected marker or did not overwrite an empty stencil.\n");
    }
    for (const int marker : indices)
    {
        marker_pairing += values[marker] * force[marker];
    }
    result[2] = std::abs(grid_pairing - marker_pairing);
    result[3] = values[4];
    if (!std::isfinite(result[2]) || !std::isfinite(result[4]))
    {
        TBOX_ERROR("Nonfinite Cartesian coupling norm for " << name << '\n');
    }
    for (int axis = 0; axis < (staggered ? NDIM : 1); ++axis)
    {
        if constexpr (C == DataCentering::SIDE)
        {
            if (!directions(axis))
            {
                continue;
            }
        }
        visit_grid<C>(field->getGhostBox(),
                      axis,
                      [&](const auto& index, const hier::Index<NDIM>& cartesian)
                      {
                          for (int depth = 0; depth < 3; ++depth)
                          {
                              TBOX_ASSERT(std::isfinite((*handle_spread)(index, depth)));
                              TBOX_ASSERT(std::abs((*handle_spread)(index, depth) - (*spread)(index, depth)) < 1.0e-10);
                          }
                          if ((*field)(index, 0) != -4.25 || (*field)(index, 1) != field_value(cartesian, axis) ||
                              (*field)(index, 2) != -9.5 || (*spread)(index, 0) != initial_spread ||
                              (*spread)(index, 2) != initial_spread ||
                              (axis != Axis && (*spread)(index, 1) != initial_spread))
                          {
                              TBOX_ERROR("Cartesian coupling modified an input field or unselected depth/axis.\n");
                          }
                      });
    }
    return result;
}

/*! \brief Combine maximum errors and representative numerical values across cases. */
void
accumulate_result(std::array<double, 5>& result, const std::array<double, 5>& next)
{
    for (int i = 0; i < 3; ++i)
    {
        result[i] = std::max(result[i], next[i]);
    }
    result[3] += next[3];
    result[4] += next[4];
}

/*! \brief Print compact numerical evidence from the dense reference comparisons. */
void
print_result(const std::string& label, const std::array<double, 5>& result)
{
    plog << label << " reference_errors " << result[0] << ' ' << result[1] << " adjoint_error " << result[2]
         << " selected_sum " << result[3] << " spread_norm " << std::sqrt(result[4]) << '\n';
}

/*! \brief Exercise every spline and axis, then compare modes with independent kernel orientation. */
template <DataCentering C>
void
check_centering()
{
    constexpr bool staggered = C == DataCentering::SIDE || C == DataCentering::FACE || C == DataCentering::EDGE;
    std::array<double, 5> spline_result{};
    MatrixFreeTest::for_each_bspline(
        [&](const std::string& name, const auto& evaluator)
        {
            [&]<std::size_t... Axis>(std::index_sequence<Axis...>) {
                (accumulate_result(spline_result, check_axis < C, staggered ? Axis : 0, Axis > (name, evaluator)), ...);
            }(std::make_index_sequence<NDIM>{});
        });
    print_result(IBTK::enum_to_string(C) + " BS2-6_CBS_adjacent", spline_result);
    const IBTK::IBKernelEvaluatorTensorProduct evaluator{ IBTK::IBKernelEvaluators::BSpline<3>{},
                                                          IBTK::IBKernelEvaluators::BSpline<2>{} };
    const auto check_mode = [&]<Mode Application>(const std::string& label)
    {
        std::array<double, 5> result{};
        [&]<std::size_t... Axis>(std::index_sequence<Axis...>)
        {
            (accumulate_result(result,
                               check_axis < C,
                               staggered ? Axis : 0,
                               (Axis + 1) % NDIM,
                               Application > ("COMPOSITE_BSPLINE_32", evaluator)),
             ...);
        }(std::make_index_sequence<NDIM>{});
        print_result(IBTK::enum_to_string(C) + " CBS32_" + label, result);
    };
    check_mode.template operator()<Mode::EXPANDED>("expanded");
    check_mode.template operator()<Mode::FACTORIZED>("factorized");
    check_mode.template operator()<Mode::CONTRACTED>("contracted");
}
} // namespace

namespace MatrixFreeTest
{
void
check_cartesian_coupling()
{
    check_centering<DataCentering::CELL>();
    check_centering<DataCentering::NODE>();
    check_centering<DataCentering::SIDE>();
    check_centering<DataCentering::FACE>();
    check_centering<DataCentering::EDGE>();
    std::array<double, 5> fallback{};
    [&]<std::size_t... Axis>(std::index_sequence<Axis...>)
    {
        (accumulate_result(fallback, check_axis<DataCentering::FACE, Axis, Axis>("COSINE_4", CartesianCosineKernel{})),
         ...);
    }(std::make_index_sequence<NDIM>{});
    print_result("FACE custom_fallback", fallback);
    const IBTK::IBKernelEvaluatorTensorProduct evaluator{ IBTK::IBKernelEvaluators::BSpline<3>{},
                                                          IBTK::IBKernelEvaluators::BSpline<2>{} };
    print_result("SIDE partial_directions",
                 check_axis<DataCentering::SIDE, NDIM - 1, 0>("COMPOSITE_BSPLINE_32", evaluator, true));
}
} // namespace MatrixFreeTest
