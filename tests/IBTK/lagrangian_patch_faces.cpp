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

// Test that interpolation and spreading are independent of the division of the domain into patches, for points
// exactly on cell faces.

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/LEInteractor.h>

#include <tbox/Pointer.h>

#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchDescriptor.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
constexpr int N = 8; // cells along each axis of the unit domain
constexpr int GHOSTS = 6;

// The field depends only on the cell index. Every layout therefore has the same value in each cell, ghost cells
// included.
double
field(const hier::Index<NDIM>& i)
{
    double value = 1.0;
    for (int d = 0; d < NDIM; ++d) value += std::sin(1.3 * i(d) + 0.7 * d) + 0.1 * i(d) * i(d) * (d + 1);
    return value;
}

// A patch over the box, in the unit domain.
Pointer<Patch<NDIM>>
make_patch(const Box<NDIM>& box)
{
    Pointer<PatchDescriptor<NDIM>> descriptor = new PatchDescriptor<NDIM>();
    Pointer<Patch<NDIM>> patch = new Patch<NDIM>(box, descriptor);
    tbox::Array<tbox::Array<bool>> regular(NDIM), periodic(NDIM);
    double dx[NDIM], x_lower[NDIM], x_upper[NDIM];
    for (int d = 0; d < NDIM; ++d)
    {
        regular[d].resizeArray(2);
        periodic[d].resizeArray(2);
        for (int side = 0; side < 2; ++side) regular[d][side] = periodic[d][side] = false;
        dx[d] = 1.0 / N;
        x_lower[d] = dx[d] * box.lower(d);
        x_upper[d] = dx[d] * (box.upper(d) + 1);
    }
    patch->setPatchGeometry(
        new CartesianPatchGeometry<NDIM>(IntVector<NDIM>(1), regular, periodic, dx, x_lower, x_upper));
    return patch;
}

// The patches of a layout: one patch for the domain, or two along each axis.
std::vector<Box<NDIM>>
layout(const bool divided)
{
    std::vector<Box<NDIM>> boxes;
    const int parts = divided ? (1 << NDIM) : 1;
    for (int p = 0; p < parts; ++p)
    {
        hier::Index<NDIM> lower, upper;
        for (int d = 0; d < NDIM; ++d)
        {
            const int half = divided ? ((p >> d) & 1) : 0;
            lower(d) = half * (N / 2);
            upper(d) = divided ? lower(d) + N / 2 - 1 : N - 1;
        }
        boxes.push_back(Box<NDIM>(lower, upper));
    }
    return boxes;
}

// Points on cell faces, including the faces shared by patches of the divided layout (at 4 cells from the lower end), at
// cell centers, and between.
std::vector<double>
points()
{
    std::vector<double> X;
    static const double along[] = { 1.0, 2.0, 3.0, 4.0, 4.5, 5.0, 5.25, 6.0 };
    const int num_along = static_cast<int>(sizeof(along) / sizeof(along[0]));
    int num_points = 1;
    for (int d = 0; d < NDIM; ++d) num_points *= num_along;
    for (int k = 0; k < num_points; ++k)
    {
        int remaining = k;
        for (int d = 0; d < NDIM; ++d)
        {
            X.push_back(along[remaining % num_along] / N);
            remaining /= num_along;
        }
    }
    return X;
}

// The values at the points of a layout, with each point interpolated on the patch whose box contains its cell.
std::vector<double>
interpolate(const std::vector<Box<NDIM>>& boxes, const std::vector<double>& X, const std::string& kernel)
{
    const int num_points = static_cast<int>(X.size() / NDIM);
    std::vector<double> result(num_points, std::numeric_limits<double>::quiet_NaN());
    std::vector<int> count(num_points, 0);
    for (const Box<NDIM>& box : boxes)
    {
        Pointer<Patch<NDIM>> patch = make_patch(box);
        Pointer<CellData<NDIM, double>> q = new CellData<NDIM, double>(box, 1, IntVector<NDIM>(GHOSTS));
        for (Box<NDIM>::Iterator b(q->getGhostBox()); b; b++) (*q)(CellIndex<NDIM>(b()), 0) = field(b());
        std::vector<double> Q(num_points, std::numeric_limits<double>::quiet_NaN());
        LEInteractor::interpolate(
            Q.data(), num_points, 1, X.data(), static_cast<int>(X.size()), NDIM, q, patch, box, kernel);
        for (int m = 0; m < num_points; ++m)
        {
            hier::Index<NDIM> i;
            for (int d = 0; d < NDIM; ++d) i(d) = static_cast<int>(std::floor(X[NDIM * m + d] * N));
            if (box.contains(i))
            {
                result[m] = Q[m];
                ++count[m];
            }
        }
    }
    for (int m = 0; m < num_points; ++m)
    {
        if (count[m] != 1)
            TBOX_ERROR(kernel << ": point " << m << " lies in " << count[m] << " patches of the layout, not one\n");
    }
    return result;
}

// The values spread into each cell, summed over the patches of a layout.
std::map<std::array<int, NDIM>, double>
spread(const std::vector<Box<NDIM>>& boxes, const std::vector<double>& X, const std::string& kernel)
{
    const int num_points = static_cast<int>(X.size() / NDIM);
    std::vector<double> F(num_points);
    for (int m = 0; m < num_points; ++m) F[m] = 1.0 + 0.01 * m;
    std::map<std::array<int, NDIM>, double> total;
    for (const Box<NDIM>& box : boxes)
    {
        Pointer<Patch<NDIM>> patch = make_patch(box);
        Pointer<CellData<NDIM, double>> q = new CellData<NDIM, double>(box, 1, IntVector<NDIM>(GHOSTS));
        q->fillAll(0.0);
        LEInteractor::spread(
            q, F.data(), num_points, 1, X.data(), static_cast<int>(X.size()), NDIM, patch, box, kernel);
        for (Box<NDIM>::Iterator b(q->getGhostBox()); b; b++)
        {
            std::array<int, NDIM> key;
            for (int d = 0; d < NDIM; ++d) key[d] = b()(d);
            total[key] += (*q)(CellIndex<NDIM>(b()), 0);
        }
    }
    return total;
}
} // namespace

int
main(int argc, char** argv)
{
    IBTK::IBTKInit ibtk_init(argc, argv);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv);
    {
        const std::vector<double> X = points();
        const std::vector<Box<NDIM>> one = layout(false), divided = layout(true);
        const std::vector<std::string> kernels = { "PIECEWISE_CONSTANT",
                                                   "DISCONTINUOUS_LINEAR",
                                                   "PIECEWISE_LINEAR",
                                                   "BSPLINE_3",
                                                   "BSPLINE_4",
                                                   "IB_3",
                                                   "IB_4",
                                                   "COMPOSITE_BSPLINE_45",
                                                   "COMPOSITE_BSPLINE_54" };
        for (const std::string& kernel : kernels)
        {
            const std::vector<double> Q1 = interpolate(one, X, kernel), Q2 = interpolate(divided, X, kernel);
            for (std::size_t m = 0; m < Q1.size(); ++m)
            {
                if (!(std::abs(Q1[m] - Q2[m]) <= 1.0e-12 * std::max(1.0, std::abs(Q1[m]))))
                {
                    std::ostringstream point;
                    for (int d = 0; d < NDIM; ++d) point << (d ? ", " : "") << X[NDIM * m + d] * N;
                    TBOX_ERROR(kernel << ": interpolation at the point (" << point.str()
                                      << ") in cells from the lower end gives " << Q1[m] << " with one patch and "
                                      << Q2[m] << " with the domain divided\n");
                }
            }
            const auto S1 = spread(one, X, kernel), S2 = spread(divided, X, kernel);
            double worst = 0.0;
            std::array<int, NDIM> worst_key{};
            for (const auto& entry : S1)
            {
                const auto it = S2.find(entry.first);
                const double difference = std::abs(entry.second - (it == S2.end() ? 0.0 : it->second));
                if (difference > worst)
                {
                    worst = difference;
                    worst_key = entry.first;
                }
            }
            for (const auto& entry : S2)
            {
                if (S1.find(entry.first) == S1.end() && std::abs(entry.second) > worst)
                {
                    worst = std::abs(entry.second);
                    worst_key = entry.first;
                }
            }
            if (!(worst <= 1.0e-10))
            {
                std::ostringstream cell;
                for (int d = 0; d < NDIM; ++d) cell << (d ? ", " : "") << worst_key[d];
                TBOX_ERROR(kernel << ": the spread values of the layouts differ by " << worst << " in the cell ("
                                  << cell.str() << ")\n");
            }
            pout << kernel << ": passed\n";
        }
    }
    return 0;
}
