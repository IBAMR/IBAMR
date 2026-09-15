// Copyright (c) 2026 by the IBAMR developers
// This file is part of IBAMR and is distributed under the 3-clause BSD license.
#include <ibtk/ibtk_utilities.h>

#include <tbox/Array.h>

#include <CartesianPatchGeometry.h>
#include <PatchDescriptor.h>
#include <SideGeometry.h>
#include <SideIterator.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <random>

#include "fixture.h"

#include <ibtk/app_namespaces.h>

namespace
{
using MatrixFreeTest::FortranCompositeInterpolate;
using MatrixFreeTest::FortranCompositeSpread;
using MatrixFreeTest::FortranInterpolate;
using MatrixFreeTest::FortranSpread;
extern "C"
{
    double IBTK_FC_FUNC_(lagrangian_ib_5_delta, LAGRANGIAN_IB_5_DELTA)(const double&);
#if (NDIM == 2)
    FortranInterpolate IBTK_FC_FUNC_(lagrangian_ib_4_interp2d, LAGRANGIAN_IB_4_INTERP2D);
    FortranSpread IBTK_FC_FUNC_(lagrangian_ib_4_spread2d, LAGRANGIAN_IB_4_SPREAD2D);
    FortranInterpolate IBTK_FC_FUNC_(lagrangian_ib_5_interp2d, LAGRANGIAN_IB_5_INTERP2D);
    FortranSpread IBTK_FC_FUNC_(lagrangian_ib_5_spread2d, LAGRANGIAN_IB_5_SPREAD2D);
    FortranInterpolate IBTK_FC_FUNC_(lagrangian_piecewise_linear_interp2d, LAGRANGIAN_PIECEWISE_LINEAR_INTERP2D);
    FortranSpread IBTK_FC_FUNC_(lagrangian_piecewise_linear_spread2d, LAGRANGIAN_PIECEWISE_LINEAR_SPREAD2D);
    FortranInterpolate IBTK_FC_FUNC_(lagrangian_bspline_3_interp2d, LAGRANGIAN_BSPLINE_3_INTERP2D);
    FortranSpread IBTK_FC_FUNC_(lagrangian_bspline_3_spread2d, LAGRANGIAN_BSPLINE_3_SPREAD2D);
    FortranInterpolate IBTK_FC_FUNC_(lagrangian_bspline_4_interp2d, LAGRANGIAN_BSPLINE_4_INTERP2D);
    FortranSpread IBTK_FC_FUNC_(lagrangian_bspline_4_spread2d, LAGRANGIAN_BSPLINE_4_SPREAD2D);
    FortranInterpolate IBTK_FC_FUNC_(lagrangian_bspline_5_interp2d, LAGRANGIAN_BSPLINE_5_INTERP2D);
    FortranSpread IBTK_FC_FUNC_(lagrangian_bspline_5_spread2d, LAGRANGIAN_BSPLINE_5_SPREAD2D);
    FortranInterpolate IBTK_FC_FUNC_(lagrangian_bspline_6_interp2d, LAGRANGIAN_BSPLINE_6_INTERP2D);
    FortranSpread IBTK_FC_FUNC_(lagrangian_bspline_6_spread2d, LAGRANGIAN_BSPLINE_6_SPREAD2D);
    FortranCompositeInterpolate IBTK_FC_FUNC_(lagrangian_discontinuous_linear_interp2d,
                                              LAGRANGIAN_DISCONTINUOUS_LINEAR_INTERP2D);
    FortranCompositeSpread IBTK_FC_FUNC_(lagrangian_discontinuous_linear_spread2d,
                                         LAGRANGIAN_DISCONTINUOUS_LINEAR_SPREAD2D);
    FortranCompositeInterpolate IBTK_FC_FUNC_(lagrangian_composite_bspline_23_interp2d,
                                              LAGRANGIAN_COMPOSITE_BSPLINE_23_INTERP2D);
    FortranCompositeSpread IBTK_FC_FUNC_(lagrangian_composite_bspline_23_spread2d,
                                         LAGRANGIAN_COMPOSITE_BSPLINE_23_SPREAD2D);
    FortranCompositeInterpolate IBTK_FC_FUNC_(lagrangian_composite_bspline_32_interp2d,
                                              LAGRANGIAN_COMPOSITE_BSPLINE_32_INTERP2D);
    FortranCompositeSpread IBTK_FC_FUNC_(lagrangian_composite_bspline_32_spread2d,
                                         LAGRANGIAN_COMPOSITE_BSPLINE_32_SPREAD2D);
    FortranCompositeInterpolate IBTK_FC_FUNC_(lagrangian_composite_bspline_34_interp2d,
                                              LAGRANGIAN_COMPOSITE_BSPLINE_34_INTERP2D);
    FortranCompositeSpread IBTK_FC_FUNC_(lagrangian_composite_bspline_34_spread2d,
                                         LAGRANGIAN_COMPOSITE_BSPLINE_34_SPREAD2D);
    FortranCompositeInterpolate IBTK_FC_FUNC_(lagrangian_composite_bspline_43_interp2d,
                                              LAGRANGIAN_COMPOSITE_BSPLINE_43_INTERP2D);
    FortranCompositeSpread IBTK_FC_FUNC_(lagrangian_composite_bspline_43_spread2d,
                                         LAGRANGIAN_COMPOSITE_BSPLINE_43_SPREAD2D);
    FortranCompositeInterpolate IBTK_FC_FUNC_(lagrangian_composite_bspline_45_interp2d,
                                              LAGRANGIAN_COMPOSITE_BSPLINE_45_INTERP2D);
    FortranCompositeSpread IBTK_FC_FUNC_(lagrangian_composite_bspline_45_spread2d,
                                         LAGRANGIAN_COMPOSITE_BSPLINE_45_SPREAD2D);
    FortranCompositeInterpolate IBTK_FC_FUNC_(lagrangian_composite_bspline_54_interp2d,
                                              LAGRANGIAN_COMPOSITE_BSPLINE_54_INTERP2D);
    FortranCompositeSpread IBTK_FC_FUNC_(lagrangian_composite_bspline_54_spread2d,
                                         LAGRANGIAN_COMPOSITE_BSPLINE_54_SPREAD2D);
    FortranCompositeInterpolate IBTK_FC_FUNC_(lagrangian_composite_bspline_56_interp2d,
                                              LAGRANGIAN_COMPOSITE_BSPLINE_56_INTERP2D);
    FortranCompositeSpread IBTK_FC_FUNC_(lagrangian_composite_bspline_56_spread2d,
                                         LAGRANGIAN_COMPOSITE_BSPLINE_56_SPREAD2D);
    FortranCompositeInterpolate IBTK_FC_FUNC_(lagrangian_composite_bspline_65_interp2d,
                                              LAGRANGIAN_COMPOSITE_BSPLINE_65_INTERP2D);
    FortranCompositeSpread IBTK_FC_FUNC_(lagrangian_composite_bspline_65_spread2d,
                                         LAGRANGIAN_COMPOSITE_BSPLINE_65_SPREAD2D);
#endif
#if (NDIM == 3)
    FortranInterpolate IBTK_FC_FUNC_(lagrangian_ib_4_interp3d, LAGRANGIAN_IB_4_INTERP3D);
    FortranSpread IBTK_FC_FUNC_(lagrangian_ib_4_spread3d, LAGRANGIAN_IB_4_SPREAD3D);
    FortranInterpolate IBTK_FC_FUNC_(lagrangian_ib_5_interp3d, LAGRANGIAN_IB_5_INTERP3D);
    FortranSpread IBTK_FC_FUNC_(lagrangian_ib_5_spread3d, LAGRANGIAN_IB_5_SPREAD3D);
    FortranInterpolate IBTK_FC_FUNC_(lagrangian_piecewise_linear_interp3d, LAGRANGIAN_PIECEWISE_LINEAR_INTERP3D);
    FortranSpread IBTK_FC_FUNC_(lagrangian_piecewise_linear_spread3d, LAGRANGIAN_PIECEWISE_LINEAR_SPREAD3D);
    FortranInterpolate IBTK_FC_FUNC_(lagrangian_bspline_3_interp3d, LAGRANGIAN_BSPLINE_3_INTERP3D);
    FortranSpread IBTK_FC_FUNC_(lagrangian_bspline_3_spread3d, LAGRANGIAN_BSPLINE_3_SPREAD3D);
    FortranInterpolate IBTK_FC_FUNC_(lagrangian_bspline_4_interp3d, LAGRANGIAN_BSPLINE_4_INTERP3D);
    FortranSpread IBTK_FC_FUNC_(lagrangian_bspline_4_spread3d, LAGRANGIAN_BSPLINE_4_SPREAD3D);
    FortranInterpolate IBTK_FC_FUNC_(lagrangian_bspline_5_interp3d, LAGRANGIAN_BSPLINE_5_INTERP3D);
    FortranSpread IBTK_FC_FUNC_(lagrangian_bspline_5_spread3d, LAGRANGIAN_BSPLINE_5_SPREAD3D);
    FortranInterpolate IBTK_FC_FUNC_(lagrangian_bspline_6_interp3d, LAGRANGIAN_BSPLINE_6_INTERP3D);
    FortranSpread IBTK_FC_FUNC_(lagrangian_bspline_6_spread3d, LAGRANGIAN_BSPLINE_6_SPREAD3D);
    FortranCompositeInterpolate IBTK_FC_FUNC_(lagrangian_discontinuous_linear_interp3d,
                                              LAGRANGIAN_DISCONTINUOUS_LINEAR_INTERP3D);
    FortranCompositeSpread IBTK_FC_FUNC_(lagrangian_discontinuous_linear_spread3d,
                                         LAGRANGIAN_DISCONTINUOUS_LINEAR_SPREAD3D);
    FortranCompositeInterpolate IBTK_FC_FUNC_(lagrangian_composite_bspline_23_interp3d,
                                              LAGRANGIAN_COMPOSITE_BSPLINE_23_INTERP3D);
    FortranCompositeSpread IBTK_FC_FUNC_(lagrangian_composite_bspline_23_spread3d,
                                         LAGRANGIAN_COMPOSITE_BSPLINE_23_SPREAD3D);
    FortranCompositeInterpolate IBTK_FC_FUNC_(lagrangian_composite_bspline_32_interp3d,
                                              LAGRANGIAN_COMPOSITE_BSPLINE_32_INTERP3D);
    FortranCompositeSpread IBTK_FC_FUNC_(lagrangian_composite_bspline_32_spread3d,
                                         LAGRANGIAN_COMPOSITE_BSPLINE_32_SPREAD3D);
    FortranCompositeInterpolate IBTK_FC_FUNC_(lagrangian_composite_bspline_34_interp3d,
                                              LAGRANGIAN_COMPOSITE_BSPLINE_34_INTERP3D);
    FortranCompositeSpread IBTK_FC_FUNC_(lagrangian_composite_bspline_34_spread3d,
                                         LAGRANGIAN_COMPOSITE_BSPLINE_34_SPREAD3D);
    FortranCompositeInterpolate IBTK_FC_FUNC_(lagrangian_composite_bspline_43_interp3d,
                                              LAGRANGIAN_COMPOSITE_BSPLINE_43_INTERP3D);
    FortranCompositeSpread IBTK_FC_FUNC_(lagrangian_composite_bspline_43_spread3d,
                                         LAGRANGIAN_COMPOSITE_BSPLINE_43_SPREAD3D);
    FortranCompositeInterpolate IBTK_FC_FUNC_(lagrangian_composite_bspline_45_interp3d,
                                              LAGRANGIAN_COMPOSITE_BSPLINE_45_INTERP3D);
    FortranCompositeSpread IBTK_FC_FUNC_(lagrangian_composite_bspline_45_spread3d,
                                         LAGRANGIAN_COMPOSITE_BSPLINE_45_SPREAD3D);
    FortranCompositeInterpolate IBTK_FC_FUNC_(lagrangian_composite_bspline_54_interp3d,
                                              LAGRANGIAN_COMPOSITE_BSPLINE_54_INTERP3D);
    FortranCompositeSpread IBTK_FC_FUNC_(lagrangian_composite_bspline_54_spread3d,
                                         LAGRANGIAN_COMPOSITE_BSPLINE_54_SPREAD3D);
    FortranCompositeInterpolate IBTK_FC_FUNC_(lagrangian_composite_bspline_56_interp3d,
                                              LAGRANGIAN_COMPOSITE_BSPLINE_56_INTERP3D);
    FortranCompositeSpread IBTK_FC_FUNC_(lagrangian_composite_bspline_56_spread3d,
                                         LAGRANGIAN_COMPOSITE_BSPLINE_56_SPREAD3D);
    FortranCompositeInterpolate IBTK_FC_FUNC_(lagrangian_composite_bspline_65_interp3d,
                                              LAGRANGIAN_COMPOSITE_BSPLINE_65_INTERP3D);
    FortranCompositeSpread IBTK_FC_FUNC_(lagrangian_composite_bspline_65_spread3d,
                                         LAGRANGIAN_COMPOSITE_BSPLINE_65_SPREAD3D);
#endif
}
} // namespace

namespace MatrixFreeTest
{
double
reference_weight(const std::string& kernel, const int axis, const int direction, const double distance)
{
    const double x = std::abs(distance);
    if (kernel == "IB_5")
    {
        // Geometry is evaluated independently of the scalar Fortran delta.
        return IBTK_FC_FUNC_(lagrangian_ib_5_delta, LAGRANGIAN_IB_5_DELTA)(x);
    }
    if (kernel == "IB_4")
    {
        if (x < 1.0)
        {
            return (3.0 - 2.0 * x + std::sqrt(1.0 + 4.0 * x - 4.0 * x * x)) / 8.0;
        }
        if (x < 2.0)
        {
            return (5.0 - 2.0 * x - std::sqrt(-7.0 + 12.0 * x - 4.0 * x * x)) / 8.0;
        }
        return 0.0;
    }
    if (kernel == "COSINE_4")
    {
        return x < 2.0 ? 0.25 * (1.0 + std::cos(0.5 * std::acos(-1.0) * x)) : 0.0;
    }
    const int width = kernel.starts_with("COMPOSITE_BSPLINE_") ?
                          kernel[kernel.size() - (axis == direction ? 2 : 1)] - '0' :
                          kernel.back() - '0';
    if (width == 1)
    {
        // Nearest-grid ties select the upper index: distance is in [-1/2, 1/2).
        return distance >= -0.5 && distance < 0.5 ? 1.0 : 0.0;
    }
    if (x >= 0.5 * width)
    {
        return 0.0;
    }
    // Truncated-power cardinal spline formula, independent of the evaluator recurrence.
    const long double t = static_cast<long double>(x) + 0.5L * width;
    long double sum = 0.0L, binomial = 1.0L, factorial = 1.0L;
    for (int k = 0; k <= width; ++k)
    {
        if (t > k)
        {
            sum += (k % 2 == 0 ? binomial : -binomial) * std::pow(t - k, width - 1);
        }
        binomial *= static_cast<long double>(width - k) / (k + 1);
        if (k > 0 && k < width)
        {
            factorial *= k;
        }
    }
    return static_cast<double>(sum / factorial);
}

Pointer<Patch<NDIM>>
make_patch(const int cells)
{
    hier::Index<NDIM> lower, upper;
    std::array<double, NDIM> dx, x_lower, x_upper;
    tbox::Array<tbox::Array<bool>> boundaries(NDIM);
    for (int d = 0; d < NDIM; ++d)
    {
        lower(d) = -3 + 7 * d;
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

void
fill_field(SideData<NDIM, double>& field)
{
    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (SideIterator<NDIM> i(field.getGhostBox(), axis); i; i++)
        {
            double value = 0.5 + 0.3 * axis;
            for (int d = 0; d < NDIM; ++d)
            {
                value += std::sin(0.11 * (d + 1) * i()(d));
            }
            field(i()) = value;
        }
    }
}

std::vector<double>
make_positions(const Patch<NDIM>& patch, const int count, const bool shuffle)
{
    const Pointer<CartesianPatchGeometry<NDIM>> geometry = patch.getPatchGeometry();
    const Box<NDIM>& box = patch.getBox();
    std::vector<int> order(count);
    std::iota(order.begin(), order.end(), 0);
    if (shuffle)
    {
        std::mt19937 generator(1729);
        std::shuffle(order.begin(), order.end(), generator);
    }
    std::vector<double> positions(NDIM * count);
    std::size_t total_cells = 1;
    for (int d = 0; d < NDIM; ++d)
    {
        total_cells *= box.numberCells(d);
    }
    for (int point = 0; point < count; ++point)
    {
        std::size_t cell = static_cast<std::size_t>(order[point]) * total_cells / count;
        for (int d = 0; d < NDIM; ++d)
        {
            const int width = box.numberCells(d);
            // Binary fractions give exact center/face ties in the regression.
            const double fraction = 0.125 + 0.125 * ((order[point] + 2 * d) % 7);
            positions[NDIM * point + d] = geometry->getXLower()[d] + (cell % width + fraction) * geometry->getDx()[d];
            cell /= width;
        }
    }
    return positions;
}

FortranInterpolate*
get_fortran_interpolate(const std::string& kernel)
{
#if (NDIM == 2)
    if (kernel == "IB_4")
    {
        return &IBTK_FC_FUNC_(lagrangian_ib_4_interp2d, LAGRANGIAN_IB_4_INTERP2D);
    }
    if (kernel == "IB_5")
    {
        return &IBTK_FC_FUNC_(lagrangian_ib_5_interp2d, LAGRANGIAN_IB_5_INTERP2D);
    }
    if (kernel == "BSPLINE_2")
    {
        return &IBTK_FC_FUNC_(lagrangian_piecewise_linear_interp2d, LAGRANGIAN_PIECEWISE_LINEAR_INTERP2D);
    }
    if (kernel == "BSPLINE_3")
    {
        return &IBTK_FC_FUNC_(lagrangian_bspline_3_interp2d, LAGRANGIAN_BSPLINE_3_INTERP2D);
    }
    if (kernel == "BSPLINE_4")
    {
        return &IBTK_FC_FUNC_(lagrangian_bspline_4_interp2d, LAGRANGIAN_BSPLINE_4_INTERP2D);
    }
    if (kernel == "BSPLINE_5")
    {
        return &IBTK_FC_FUNC_(lagrangian_bspline_5_interp2d, LAGRANGIAN_BSPLINE_5_INTERP2D);
    }
    if (kernel == "BSPLINE_6")
    {
        return &IBTK_FC_FUNC_(lagrangian_bspline_6_interp2d, LAGRANGIAN_BSPLINE_6_INTERP2D);
    }
#endif
#if (NDIM == 3)
    if (kernel == "IB_4")
    {
        return &IBTK_FC_FUNC_(lagrangian_ib_4_interp3d, LAGRANGIAN_IB_4_INTERP3D);
    }
    if (kernel == "IB_5")
    {
        return &IBTK_FC_FUNC_(lagrangian_ib_5_interp3d, LAGRANGIAN_IB_5_INTERP3D);
    }
    if (kernel == "BSPLINE_2")
    {
        return &IBTK_FC_FUNC_(lagrangian_piecewise_linear_interp3d, LAGRANGIAN_PIECEWISE_LINEAR_INTERP3D);
    }
    if (kernel == "BSPLINE_3")
    {
        return &IBTK_FC_FUNC_(lagrangian_bspline_3_interp3d, LAGRANGIAN_BSPLINE_3_INTERP3D);
    }
    if (kernel == "BSPLINE_4")
    {
        return &IBTK_FC_FUNC_(lagrangian_bspline_4_interp3d, LAGRANGIAN_BSPLINE_4_INTERP3D);
    }
    if (kernel == "BSPLINE_5")
    {
        return &IBTK_FC_FUNC_(lagrangian_bspline_5_interp3d, LAGRANGIAN_BSPLINE_5_INTERP3D);
    }
    if (kernel == "BSPLINE_6")
    {
        return &IBTK_FC_FUNC_(lagrangian_bspline_6_interp3d, LAGRANGIAN_BSPLINE_6_INTERP3D);
    }
#endif
    return nullptr;
}

FortranSpread*
get_fortran_spread(const std::string& kernel)
{
#if (NDIM == 2)
    if (kernel == "IB_4")
    {
        return &IBTK_FC_FUNC_(lagrangian_ib_4_spread2d, LAGRANGIAN_IB_4_SPREAD2D);
    }
    if (kernel == "IB_5")
    {
        return &IBTK_FC_FUNC_(lagrangian_ib_5_spread2d, LAGRANGIAN_IB_5_SPREAD2D);
    }
    if (kernel == "BSPLINE_2")
    {
        return &IBTK_FC_FUNC_(lagrangian_piecewise_linear_spread2d, LAGRANGIAN_PIECEWISE_LINEAR_SPREAD2D);
    }
    if (kernel == "BSPLINE_3")
    {
        return &IBTK_FC_FUNC_(lagrangian_bspline_3_spread2d, LAGRANGIAN_BSPLINE_3_SPREAD2D);
    }
    if (kernel == "BSPLINE_4")
    {
        return &IBTK_FC_FUNC_(lagrangian_bspline_4_spread2d, LAGRANGIAN_BSPLINE_4_SPREAD2D);
    }
    if (kernel == "BSPLINE_5")
    {
        return &IBTK_FC_FUNC_(lagrangian_bspline_5_spread2d, LAGRANGIAN_BSPLINE_5_SPREAD2D);
    }
    if (kernel == "BSPLINE_6")
    {
        return &IBTK_FC_FUNC_(lagrangian_bspline_6_spread2d, LAGRANGIAN_BSPLINE_6_SPREAD2D);
    }
#endif
#if (NDIM == 3)
    if (kernel == "IB_4")
    {
        return &IBTK_FC_FUNC_(lagrangian_ib_4_spread3d, LAGRANGIAN_IB_4_SPREAD3D);
    }
    if (kernel == "IB_5")
    {
        return &IBTK_FC_FUNC_(lagrangian_ib_5_spread3d, LAGRANGIAN_IB_5_SPREAD3D);
    }
    if (kernel == "BSPLINE_2")
    {
        return &IBTK_FC_FUNC_(lagrangian_piecewise_linear_spread3d, LAGRANGIAN_PIECEWISE_LINEAR_SPREAD3D);
    }
    if (kernel == "BSPLINE_3")
    {
        return &IBTK_FC_FUNC_(lagrangian_bspline_3_spread3d, LAGRANGIAN_BSPLINE_3_SPREAD3D);
    }
    if (kernel == "BSPLINE_4")
    {
        return &IBTK_FC_FUNC_(lagrangian_bspline_4_spread3d, LAGRANGIAN_BSPLINE_4_SPREAD3D);
    }
    if (kernel == "BSPLINE_5")
    {
        return &IBTK_FC_FUNC_(lagrangian_bspline_5_spread3d, LAGRANGIAN_BSPLINE_5_SPREAD3D);
    }
    if (kernel == "BSPLINE_6")
    {
        return &IBTK_FC_FUNC_(lagrangian_bspline_6_spread3d, LAGRANGIAN_BSPLINE_6_SPREAD3D);
    }
#endif
    return nullptr;
}

FortranCompositeInterpolate*
get_fortran_composite_interpolate(const std::string& kernel)
{
#if (NDIM == 2)
    if (kernel == "COMPOSITE_BSPLINE_21")
    {
        return &IBTK_FC_FUNC_(lagrangian_discontinuous_linear_interp2d, LAGRANGIAN_DISCONTINUOUS_LINEAR_INTERP2D);
    }
    if (kernel == "COMPOSITE_BSPLINE_23")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_23_interp2d, LAGRANGIAN_COMPOSITE_BSPLINE_23_INTERP2D);
    }
    if (kernel == "COMPOSITE_BSPLINE_32")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_32_interp2d, LAGRANGIAN_COMPOSITE_BSPLINE_32_INTERP2D);
    }
    if (kernel == "COMPOSITE_BSPLINE_34")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_34_interp2d, LAGRANGIAN_COMPOSITE_BSPLINE_34_INTERP2D);
    }
    if (kernel == "COMPOSITE_BSPLINE_43")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_43_interp2d, LAGRANGIAN_COMPOSITE_BSPLINE_43_INTERP2D);
    }
    if (kernel == "COMPOSITE_BSPLINE_45")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_45_interp2d, LAGRANGIAN_COMPOSITE_BSPLINE_45_INTERP2D);
    }
    if (kernel == "COMPOSITE_BSPLINE_54")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_54_interp2d, LAGRANGIAN_COMPOSITE_BSPLINE_54_INTERP2D);
    }
    if (kernel == "COMPOSITE_BSPLINE_56")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_56_interp2d, LAGRANGIAN_COMPOSITE_BSPLINE_56_INTERP2D);
    }
    if (kernel == "COMPOSITE_BSPLINE_65")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_65_interp2d, LAGRANGIAN_COMPOSITE_BSPLINE_65_INTERP2D);
    }
#endif
#if (NDIM == 3)
    if (kernel == "COMPOSITE_BSPLINE_21")
    {
        return &IBTK_FC_FUNC_(lagrangian_discontinuous_linear_interp3d, LAGRANGIAN_DISCONTINUOUS_LINEAR_INTERP3D);
    }
    if (kernel == "COMPOSITE_BSPLINE_23")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_23_interp3d, LAGRANGIAN_COMPOSITE_BSPLINE_23_INTERP3D);
    }
    if (kernel == "COMPOSITE_BSPLINE_32")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_32_interp3d, LAGRANGIAN_COMPOSITE_BSPLINE_32_INTERP3D);
    }
    if (kernel == "COMPOSITE_BSPLINE_34")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_34_interp3d, LAGRANGIAN_COMPOSITE_BSPLINE_34_INTERP3D);
    }
    if (kernel == "COMPOSITE_BSPLINE_43")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_43_interp3d, LAGRANGIAN_COMPOSITE_BSPLINE_43_INTERP3D);
    }
    if (kernel == "COMPOSITE_BSPLINE_45")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_45_interp3d, LAGRANGIAN_COMPOSITE_BSPLINE_45_INTERP3D);
    }
    if (kernel == "COMPOSITE_BSPLINE_54")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_54_interp3d, LAGRANGIAN_COMPOSITE_BSPLINE_54_INTERP3D);
    }
    if (kernel == "COMPOSITE_BSPLINE_56")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_56_interp3d, LAGRANGIAN_COMPOSITE_BSPLINE_56_INTERP3D);
    }
    if (kernel == "COMPOSITE_BSPLINE_65")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_65_interp3d, LAGRANGIAN_COMPOSITE_BSPLINE_65_INTERP3D);
    }
#endif
    return nullptr;
}

FortranCompositeSpread*
get_fortran_composite_spread(const std::string& kernel)
{
#if (NDIM == 2)
    if (kernel == "COMPOSITE_BSPLINE_21")
    {
        return &IBTK_FC_FUNC_(lagrangian_discontinuous_linear_spread2d, LAGRANGIAN_DISCONTINUOUS_LINEAR_SPREAD2D);
    }
    if (kernel == "COMPOSITE_BSPLINE_23")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_23_spread2d, LAGRANGIAN_COMPOSITE_BSPLINE_23_SPREAD2D);
    }
    if (kernel == "COMPOSITE_BSPLINE_32")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_32_spread2d, LAGRANGIAN_COMPOSITE_BSPLINE_32_SPREAD2D);
    }
    if (kernel == "COMPOSITE_BSPLINE_34")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_34_spread2d, LAGRANGIAN_COMPOSITE_BSPLINE_34_SPREAD2D);
    }
    if (kernel == "COMPOSITE_BSPLINE_43")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_43_spread2d, LAGRANGIAN_COMPOSITE_BSPLINE_43_SPREAD2D);
    }
    if (kernel == "COMPOSITE_BSPLINE_45")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_45_spread2d, LAGRANGIAN_COMPOSITE_BSPLINE_45_SPREAD2D);
    }
    if (kernel == "COMPOSITE_BSPLINE_54")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_54_spread2d, LAGRANGIAN_COMPOSITE_BSPLINE_54_SPREAD2D);
    }
    if (kernel == "COMPOSITE_BSPLINE_56")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_56_spread2d, LAGRANGIAN_COMPOSITE_BSPLINE_56_SPREAD2D);
    }
    if (kernel == "COMPOSITE_BSPLINE_65")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_65_spread2d, LAGRANGIAN_COMPOSITE_BSPLINE_65_SPREAD2D);
    }
#endif
#if (NDIM == 3)
    if (kernel == "COMPOSITE_BSPLINE_21")
    {
        return &IBTK_FC_FUNC_(lagrangian_discontinuous_linear_spread3d, LAGRANGIAN_DISCONTINUOUS_LINEAR_SPREAD3D);
    }
    if (kernel == "COMPOSITE_BSPLINE_23")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_23_spread3d, LAGRANGIAN_COMPOSITE_BSPLINE_23_SPREAD3D);
    }
    if (kernel == "COMPOSITE_BSPLINE_32")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_32_spread3d, LAGRANGIAN_COMPOSITE_BSPLINE_32_SPREAD3D);
    }
    if (kernel == "COMPOSITE_BSPLINE_34")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_34_spread3d, LAGRANGIAN_COMPOSITE_BSPLINE_34_SPREAD3D);
    }
    if (kernel == "COMPOSITE_BSPLINE_43")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_43_spread3d, LAGRANGIAN_COMPOSITE_BSPLINE_43_SPREAD3D);
    }
    if (kernel == "COMPOSITE_BSPLINE_45")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_45_spread3d, LAGRANGIAN_COMPOSITE_BSPLINE_45_SPREAD3D);
    }
    if (kernel == "COMPOSITE_BSPLINE_54")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_54_spread3d, LAGRANGIAN_COMPOSITE_BSPLINE_54_SPREAD3D);
    }
    if (kernel == "COMPOSITE_BSPLINE_56")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_56_spread3d, LAGRANGIAN_COMPOSITE_BSPLINE_56_SPREAD3D);
    }
    if (kernel == "COMPOSITE_BSPLINE_65")
    {
        return &IBTK_FC_FUNC_(lagrangian_composite_bspline_65_spread3d, LAGRANGIAN_COMPOSITE_BSPLINE_65_SPREAD3D);
    }
#endif
    return nullptr;
}

} // namespace MatrixFreeTest
