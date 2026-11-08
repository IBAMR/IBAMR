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

// Test that the LEInteractor overloads that take index data select each point once and give the same results for
// every division of the domain into patches, for every kernel and for cell-, node-, side-, and (in three dimensions)
// edge-centered data. The points lie on cell faces and centers, within six ulps of them, at random positions, and up to
// 0.9 cells beyond the domain, as between redistributions: beyond a periodic boundary they are found through their
// periodic images. On three layouts of patches, the program checks that each point is selected on exactly one patch,
// that the weights of each point sum to 1, that interpolation gives the same values on every layout, and that
// spreading into the ghost regions gives the same grid values on every layout. It also checks that interpolation
// with the continuous kernels agrees with the overloads that take positions only, at the points inside the domain.

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/LNode.h>
#include <ibtk/LNodeSet.h>
#include <ibtk/LNodeSetData.h>

#include <tbox/Pointer.h>

#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIndex.h>
#include <EdgeData.h>
#include <EdgeGeometry.h>
#include <EdgeIndex.h>
#include <NodeData.h>
#include <NodeGeometry.h>
#include <NodeIndex.h>
#include <PatchLevel.h>
#include <ProcessorMapping.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
// The four-point IB kernel as a user-defined kernel with an even stencil size.
double
four_point_kernel(const double r)
{
    const double abs_r = std::abs(r);
    if (abs_r < 1.0)
    {
        return 0.125 * (3.0 - 2.0 * abs_r + std::sqrt(1.0 + 4.0 * abs_r - 4.0 * abs_r * abs_r));
    }
    if (abs_r < 2.0)
    {
        return 0.125 * (5.0 - 2.0 * abs_r - std::sqrt(-7.0 + 12.0 * abs_r - 4.0 * abs_r * abs_r));
    }
    return 0.0;
}

struct Domain
{
    std::array<double, NDIM> x_lo{}, x_up{};
    std::array<bool, NDIM> periodic{};
    int N = 0;
};

// The positions along one axis: each face and center and the six doubles on either side of it, eight random
// positions, and positions beyond each end of the domain.
std::vector<double>
axis_positions(const Domain& domain, const int d)
{
    const double x_lo = domain.x_lo[d], x_up = domain.x_up[d];
    const double h = (x_up - x_lo) / domain.N;
    std::vector<double> positions;
    const auto add_with_neighbors = [&positions](const double x)
    {
        positions.push_back(x);
        double below = x, above = x;
        for (int k = 0; k < 6; ++k)
        {
            below = std::nextafter(below, -std::numeric_limits<double>::infinity());
            above = std::nextafter(above, std::numeric_limits<double>::infinity());
            positions.push_back(below);
            positions.push_back(above);
        }
    };
    // std::fma rounds once, so the positions do not depend on whether the compiler contracts multiply-adds.
    for (int j = 0; j <= domain.N; ++j)
    {
        add_with_neighbors(std::fma(j, h, x_lo));
        if (j < domain.N)
        {
            add_with_neighbors(std::fma(j + 0.5, h, x_lo));
        }
    }
    for (int k = 1; k <= 8; ++k)
    {
        const double fraction = std::fma(k, 0.6180339887498949, -std::floor(k * 0.6180339887498949));
        positions.push_back(std::fma(fraction, x_up - x_lo, x_lo));
    }
    for (const double cells : { 0.3, 0.9 })
    {
        positions.push_back(std::fma(-cells, h, x_lo));
        positions.push_back(std::fma(cells, h, x_up));
    }
    positions.push_back(x_lo - 1.0e-16);
    positions.push_back(std::nextafter(x_up, std::numeric_limits<double>::infinity()));
    std::vector<double> kept;
    for (const double x : positions)
    {
        if (x > x_lo - h && x < x_up + h)
        {
            kept.push_back(x);
        }
    }
    std::sort(kept.begin(), kept.end());
    kept.erase(std::unique(kept.begin(), kept.end()), kept.end());
    return kept;
}

// The points: for each axis, every position along it, combined with positions along the other axes in a fixed order.
std::vector<double>
build_points(const Domain& domain)
{
    std::array<std::vector<double>, NDIM> along;
    for (int d = 0; d < NDIM; ++d)
    {
        along[d] = axis_positions(domain, d);
    }
    std::vector<double> X;
    for (int lead = 0; lead < NDIM; ++lead)
    {
        for (std::size_t k = 0; k < along[lead].size(); ++k)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                const std::size_t n = along[d].size();
                X.push_back(d == lead ? along[d][k] : along[d][(7 * k + 13 * d + 5) % n]);
            }
        }
    }
    return X;
}

bool
inside_domain(const Domain& domain, const double* X)
{
    for (int d = 0; d < NDIM; ++d)
    {
        if (X[d] < domain.x_lo[d] || X[d] >= domain.x_up[d])
        {
            return false;
        }
    }
    return true;
}

// A level of patches whose boxes divide each axis at the given cells.
Pointer<PatchLevel<NDIM>>
make_level(Pointer<CartesianGridGeometry<NDIM>> grid_geom, const Domain& domain, const std::vector<int>& splits)
{
    std::vector<std::pair<int, int>> intervals;
    int lower = 0;
    for (const int split : splits)
    {
        intervals.emplace_back(lower, split - 1);
        lower = split;
    }
    intervals.emplace_back(lower, domain.N - 1);
    int num_patches = 1;
    for (int d = 0; d < NDIM; ++d)
    {
        num_patches *= static_cast<int>(intervals.size());
    }
    BoxArray<NDIM> boxes(num_patches);
    ProcessorMapping mapping(num_patches);
    for (int p = 0; p < num_patches; ++p)
    {
        hier::Index<NDIM> box_lower(0), box_upper(0);
        int remaining = p;
        for (int d = 0; d < NDIM; ++d)
        {
            const std::pair<int, int>& interval = intervals[remaining % intervals.size()];
            remaining /= static_cast<int>(intervals.size());
            box_lower(d) = interval.first;
            box_upper(d) = interval.second;
        }
        boxes[p] = Box<NDIM>(box_lower, box_upper);
        mapping.setProcessorAssignment(p, 0);
    }
    return new PatchLevel<NDIM>(
        boxes, mapping, IntVector<NDIM>(1), grid_geom, VariableDatabase<NDIM>::getDatabase()->getPatchDescriptor());
}

// Index data of the patch with each point in the cell of its position, clamped to the domain, and in each periodic
// image of that cell that lies in the ghost box, as after a redistribution. The local PETSc index of a point is its
// number.
Pointer<LNodeSetData>
build_index_data(Pointer<Patch<NDIM>> patch,
                 Pointer<CartesianGridGeometry<NDIM>> grid_geom,
                 const Domain& domain,
                 const std::vector<double>& X,
                 const int ghost_width)
{
    Pointer<LNodeSetData> idx_data = new LNodeSetData(patch->getBox(), IntVector<NDIM>(ghost_width));
    const Box<NDIM>& ghost_box = idx_data->getGhostBox();
    const int num_points = static_cast<int>(X.size() / NDIM);
    int num_images = 1;
    for (int d = 0; d < NDIM; ++d)
    {
        num_images *= 3;
    }
    for (int k = 0; k < num_points; ++k)
    {
        hier::Index<NDIM> cell(0);
        for (int d = 0; d < NDIM; ++d)
        {
            const double h = (domain.x_up[d] - domain.x_lo[d]) / domain.N;
            const int c = static_cast<int>(std::floor((X[NDIM * k + d] - domain.x_lo[d]) / h));
            cell(d) = std::min(std::max(c, 0), domain.N - 1);
        }
        for (int image = 0; image < num_images; ++image)
        {
            CellIndex<NDIM> i(cell);
            bool valid = true;
            int remaining = image;
            for (int d = 0; d < NDIM; ++d)
            {
                const int shift = remaining % 3 - 1;
                remaining /= 3;
                if (shift != 0 && !domain.periodic[d])
                {
                    valid = false;
                }
                i(d) += shift * domain.N;
            }
            if (!valid || !ghost_box.contains(i))
            {
                continue;
            }
            if (!idx_data->isElement(i))
            {
                idx_data->appendItemPointer(i, new LNodeSet());
            }
            idx_data->getItem(i)->push_back(new LNode(k, k, k));
        }
    }
    idx_data->cacheLocalIndices(patch, grid_geom, IntVector<NDIM>(1));
    return idx_data;
}

struct Layout
{
    std::string name;
    Pointer<PatchLevel<NDIM>> level;
    std::vector<Pointer<LNodeSetData>> idx_data;
};

// The data index wrapped into the domain along each periodic axis.
hier::Index<NDIM>
wrap_index(const hier::Index<NDIM>& i, const Domain& domain)
{
    hier::Index<NDIM> wrapped = i;
    for (int d = 0; d < NDIM; ++d)
    {
        if (domain.periodic[d])
        {
            wrapped(d) = ((i(d) % domain.N) + domain.N) % domain.N;
        }
    }
    return wrapped;
}

// A field of the data index that varies from cell to cell and is periodic along the periodic axes, or 1.
double
field(const hier::Index<NDIM>& i, const Domain& domain, const int component, const bool one)
{
    if (one)
    {
        return 1.0;
    }
    const hier::Index<NDIM> w = wrap_index(i, domain);
    double value = 2.0 + 0.1 * component;
    for (int d = 0; d < NDIM; ++d)
    {
        value += std::sin(1.7 * w(d) + 0.9 * d + 0.3 * component) + 0.05 * w(d);
    }
    return value;
}

enum class Centering
{
    CELL,
    NODE,
    SIDE,
    EDGE
};

std::string
centering_name(const Centering centering)
{
    switch (centering)
    {
    case Centering::CELL:
        return "cell";
    case Centering::NODE:
        return "node";
    case Centering::SIDE:
        return "side";
    case Centering::EDGE:
        return "edge";
    }
    return "";
}

int
centering_depth(const Centering centering)
{
    return centering == Centering::SIDE || centering == Centering::EDGE ? NDIM : 1;
}

// Patch data of one centering, with one array per component for side and edge data.
struct PatchData
{
    Pointer<CellData<NDIM, double>> cell;
    Pointer<NodeData<NDIM, double>> node;
    Pointer<SideData<NDIM, double>> side;
#if (NDIM == 3)
    Pointer<EdgeData<NDIM, double>> edge;
#endif
};

PatchData
make_data(const Box<NDIM>& box, const Centering centering, const int ghost_width)
{
    PatchData data;
    const IntVector<NDIM> ghosts(ghost_width);
    switch (centering)
    {
    case Centering::CELL:
        data.cell = new CellData<NDIM, double>(box, 1, ghosts);
        data.cell->fillAll(0.0);
        break;
    case Centering::NODE:
        data.node = new NodeData<NDIM, double>(box, 1, ghosts);
        data.node->fillAll(0.0);
        break;
    case Centering::SIDE:
        data.side = new SideData<NDIM, double>(box, 1, ghosts);
        data.side->fillAll(0.0);
        break;
    case Centering::EDGE:
#if (NDIM == 3)
        data.edge = new EdgeData<NDIM, double>(box, 1, ghosts);
        data.edge->fillAll(0.0);
#endif
        break;
    }
    return data;
}

// Call fcn(index, component, value) for each data point of the patch data, in the data box if interior and in the
// ghost box otherwise.
template <class Function>
void
for_each_data_point(PatchData& data, const Centering centering, const bool interior, Function fcn)
{
    switch (centering)
    {
    case Centering::CELL:
    {
        const Box<NDIM> box = interior ? data.cell->getBox() : data.cell->getGhostBox();
        for (Box<NDIM>::Iterator b(box); b; b++)
        {
            fcn(b(), 0, (*data.cell)(CellIndex<NDIM>(b()), 0));
        }
        break;
    }
    case Centering::NODE:
    {
        const Box<NDIM> box = NodeGeometry<NDIM>::toNodeBox(interior ? data.node->getBox() : data.node->getGhostBox());
        for (Box<NDIM>::Iterator b(box); b; b++)
        {
            fcn(b(), 0, (*data.node)(NodeIndex<NDIM>(b(), IntVector<NDIM>(0)), 0));
        }
        break;
    }
    case Centering::SIDE:
    {
        for (int axis = 0; axis < NDIM; ++axis)
        {
            const Box<NDIM> box =
                SideGeometry<NDIM>::toSideBox(interior ? data.side->getBox() : data.side->getGhostBox(), axis);
            for (Box<NDIM>::Iterator b(box); b; b++)
            {
                fcn(b(), axis, (*data.side)(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Lower)));
            }
        }
        break;
    }
    case Centering::EDGE:
    {
#if (NDIM == 3)
        for (int axis = 0; axis < NDIM; ++axis)
        {
            const Box<NDIM> box =
                EdgeGeometry<NDIM>::toEdgeBox(interior ? data.edge->getBox() : data.edge->getGhostBox(), axis);
            for (Box<NDIM>::Iterator b(box); b; b++)
            {
                fcn(b(), axis, (*data.edge)(EdgeIndex<NDIM>(b(), axis, EdgeIndex<NDIM>::LowerLeft)));
            }
        }
#endif
        break;
    }
    }
}

// Interpolate at every point on every patch of the layout, through the index data or, if positions, through the
// positions only. Each point must be selected on exactly one patch, except that the positions select no point outside
// the domain.
std::vector<double>
interpolate(const Layout& layout,
            const Domain& domain,
            const std::vector<double>& X,
            const std::string& kernel,
            const Centering centering,
            const int ghost_width,
            const bool positions,
            const bool one)
{
    const int num_points = static_cast<int>(X.size() / NDIM);
    const int depth = centering_depth(centering);
    std::vector<double> result(depth * num_points, std::numeric_limits<double>::quiet_NaN());
    std::vector<int> count(num_points, 0);
    int patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(layout.level); p; p++, ++patch_num)
    {
        Pointer<Patch<NDIM>> patch = layout.level->getPatch(p());
        const Box<NDIM>& box = patch->getBox();
        PatchData data = make_data(box, centering, ghost_width);
        for_each_data_point(data,
                            centering,
                            /*interior*/ false,
                            [&](const hier::Index<NDIM>& i, const int component, double& value)
                            { value = field(i, domain, component, one); });
        std::vector<double> Q(depth * num_points, std::numeric_limits<double>::quiet_NaN());
        const Pointer<LNodeSetData>& idx_data = layout.idx_data[patch_num];
        const int Q_size = static_cast<int>(Q.size()), X_size = static_cast<int>(X.size());
        switch (centering)
        {
        case Centering::CELL:
            if (positions)
            {
                LEInteractor::interpolate(
                    Q.data(), Q_size, depth, X.data(), X_size, NDIM, data.cell, patch, box, kernel);
            }
            else
            {
                LEInteractor::interpolate(Q.data(), depth, X.data(), NDIM, idx_data, data.cell, patch, box, kernel);
            }
            break;
        case Centering::NODE:
            if (positions)
            {
                LEInteractor::interpolate(
                    Q.data(), Q_size, depth, X.data(), X_size, NDIM, data.node, patch, box, kernel);
            }
            else
            {
                LEInteractor::interpolate(Q.data(), depth, X.data(), NDIM, idx_data, data.node, patch, box, kernel);
            }
            break;
        case Centering::SIDE:
            if (positions)
            {
                LEInteractor::interpolate(
                    Q.data(), Q_size, depth, X.data(), X_size, NDIM, data.side, patch, box, kernel);
            }
            else
            {
                LEInteractor::interpolate(Q.data(), depth, X.data(), NDIM, idx_data, data.side, patch, box, kernel);
            }
            break;
        case Centering::EDGE:
#if (NDIM == 3)
            if (positions)
            {
                LEInteractor::interpolate(
                    Q.data(), Q_size, depth, X.data(), X_size, NDIM, data.edge, patch, box, kernel);
            }
            else
            {
                LEInteractor::interpolate(Q.data(), depth, X.data(), NDIM, idx_data, data.edge, patch, box, kernel);
            }
#endif
            break;
        }

        // A point is selected on this patch if the patch wrote its value.
        for (int k = 0; k < num_points; ++k)
        {
            if (std::isnan(Q[depth * k]))
            {
                continue;
            }
            ++count[k];
            std::copy(&Q[depth * k], &Q[depth * k] + depth, &result[depth * k]);
        }
    }
    for (int k = 0; k < num_points; ++k)
    {
        const int expected = positions && !inside_domain(domain, &X[NDIM * k]) ? 0 : 1;
        if (count[k] != expected)
        {
            TBOX_ERROR(kernel << ", " << centering_name(centering) << " data, layout " << layout.name << ": the "
                              << (positions ? "positions select" : "index data select") << " point " << k << " on "
                              << count[k] << " patches, not " << expected << "\n");
        }
    }
    return result;
}

// Spread into the ghost region of every patch of the layout through the index data, and return the values in the
// data boxes of the patches, by data index wrapped into the domain and by component. A data point that two patches
// hold, or that a periodic boundary repeats, must have the same value on each.
std::map<std::array<int, NDIM + 1>, double>
spread(const Layout& layout,
       const Domain& domain,
       const std::vector<double>& X,
       const std::string& kernel,
       const Centering centering,
       const int ghost_width)
{
    const int num_points = static_cast<int>(X.size() / NDIM);
    const int depth = centering_depth(centering);
    std::vector<double> F(depth * num_points);
    for (int k = 0; k < depth * num_points; ++k)
    {
        F[k] = 1.0 + 0.01 * (k % 97);
    }
    std::map<std::array<int, NDIM + 1>, double> values;
    double largest = 0.0, mismatch = 0.0;
    int patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(layout.level); p; p++, ++patch_num)
    {
        Pointer<Patch<NDIM>> patch = layout.level->getPatch(p());
        PatchData data = make_data(patch->getBox(), centering, ghost_width);
        const Pointer<LNodeSetData>& idx_data = layout.idx_data[patch_num];
        const Box<NDIM> ghost_box = Box<NDIM>::grow(patch->getBox(), IntVector<NDIM>(ghost_width));
        switch (centering)
        {
        case Centering::CELL:
            LEInteractor::spread(data.cell, F.data(), depth, X.data(), NDIM, idx_data, patch, ghost_box, kernel);
            break;
        case Centering::NODE:
            LEInteractor::spread(data.node, F.data(), depth, X.data(), NDIM, idx_data, patch, ghost_box, kernel);
            break;
        case Centering::SIDE:
            LEInteractor::spread(data.side, F.data(), depth, X.data(), NDIM, idx_data, patch, ghost_box, kernel);
            break;
        case Centering::EDGE:
#if (NDIM == 3)
            LEInteractor::spread(data.edge, F.data(), depth, X.data(), NDIM, idx_data, patch, ghost_box, kernel);
#endif
            break;
        }
        for_each_data_point(data,
                            centering,
                            /*interior*/ true,
                            [&](const hier::Index<NDIM>& i, const int component, const double& value)
                            {
                                const hier::Index<NDIM> w = wrap_index(i, domain);
                                std::array<int, NDIM + 1> key{};
                                for (int d = 0; d < NDIM; ++d)
                                {
                                    key[d] = w(d);
                                }
                                key[NDIM] = component;
                                largest = std::max(largest, std::abs(value));
                                const std::map<std::array<int, NDIM + 1>, double>::iterator it = values.find(key);
                                if (it == values.end())
                                {
                                    values[key] = value;
                                }
                                else
                                {
                                    mismatch = std::max(mismatch, std::abs(it->second - value));
                                }
                            });
    }
    if (!(mismatch <= 1.0e-12 * largest))
    {
        TBOX_ERROR(kernel << ", " << centering_name(centering) << " data, layout " << layout.name
                          << ": spreading gives values that differ by " << mismatch
                          << " at the same data point, relative to the largest value " << largest << "\n");
    }
    return values;
}

// The largest difference between the values of two interpolations, relative to the values and at least 1.
double
max_difference(const std::vector<double>& a, const std::vector<double>& b)
{
    double worst = 0.0;
    for (std::size_t k = 0; k < a.size(); ++k)
    {
        if (std::isnan(a[k]) || std::isnan(b[k]))
        {
            continue;
        }
        worst = std::max(worst, std::abs(a[k] - b[k]) / std::max(1.0, std::abs(a[k])));
    }
    return worst;
}

void
check(const std::vector<Layout>& layouts,
      const Domain& domain,
      const std::vector<double>& X,
      const std::string& kernel,
      const Centering centering,
      const int ghost_width)
{
    const std::string name = kernel + ", " + centering_name(centering) + " data";
    const int num_points = static_cast<int>(X.size() / NDIM);
    const int depth = centering_depth(centering);

    // The weights of every point sum to 1.
    const std::vector<double> ones =
        interpolate(layouts[0], domain, X, kernel, centering, ghost_width, /*positions*/ false, /*one*/ true);
    for (std::size_t k = 0; k < ones.size(); ++k)
    {
        if (!(std::abs(ones[k] - 1.0) <= 1.0e-12))
        {
            TBOX_ERROR(name << ": interpolation of the field 1 gives " << ones[k] << " at point " << k / depth << "\n");
        }
    }

    // Interpolation gives the same values on every layout and, for the continuous kernels, the same values as the
    // positions at the points inside the domain. At a point midway between two data points, the discontinuous kernels
    // may take either one.
    const std::vector<double> Q =
        interpolate(layouts[0], domain, X, kernel, centering, ghost_width, /*positions*/ false, /*one*/ false);
    for (std::size_t l = 1; l < layouts.size(); ++l)
    {
        const std::vector<double> Q_l =
            interpolate(layouts[l], domain, X, kernel, centering, ghost_width, /*positions*/ false, /*one*/ false);
        const double difference = max_difference(Q, Q_l);
        if (!(difference <= 1.0e-12))
        {
            TBOX_ERROR(name << ": interpolation on the layouts " << layouts[0].name << " and " << layouts[l].name
                            << " differs by " << difference << "\n");
        }
    }
    if (kernel != "PIECEWISE_CONSTANT" && kernel != "DISCONTINUOUS_LINEAR")
    {
        const std::vector<double> Q_positions =
            interpolate(layouts[0], domain, X, kernel, centering, ghost_width, /*positions*/ true, /*one*/ false);
        const double positions_difference = max_difference(Q, Q_positions);
        if (!(positions_difference <= 1.0e-12))
        {
            TBOX_ERROR(name << ": interpolation through the index data and the positions differs by "
                            << positions_difference << "\n");
        }
    }

    // Spreading gives the same grid values on every layout.
    const std::map<std::array<int, NDIM + 1>, double> S = spread(layouts[0], domain, X, kernel, centering, ghost_width);
    double largest = 0.0;
    for (const std::pair<const std::array<int, NDIM + 1>, double>& entry : S)
    {
        largest = std::max(largest, std::abs(entry.second));
    }
    for (std::size_t l = 1; l < layouts.size(); ++l)
    {
        const std::map<std::array<int, NDIM + 1>, double> S_l =
            spread(layouts[l], domain, X, kernel, centering, ghost_width);
        double difference = S.size() == S_l.size() ? 0.0 : largest;
        for (const std::pair<const std::array<int, NDIM + 1>, double>& entry : S)
        {
            const std::map<std::array<int, NDIM + 1>, double>::const_iterator it = S_l.find(entry.first);
            difference = std::max(difference, std::abs(entry.second - (it == S_l.end() ? 0.0 : it->second)));
        }
        if (!(difference <= 1.0e-12 * largest))
        {
            TBOX_ERROR(name << ": spreading on the layouts " << layouts[0].name << " and " << layouts[l].name
                            << " differs by " << difference << ", relative to the largest value " << largest << "\n");
        }
    }

    // Print the sum of the interpolated values, which does not depend on the layout.
    double sum = 0.0;
    for (int k = 0; k < depth * num_points; ++k)
    {
        sum += Q[k];
    }
    pout << name << ": passed, sum of the interpolated values = " << std::setprecision(10) << sum << '\n';
}
} // namespace

int
main(int argc, char** argv)
{
    IBTK::IBTKInit ibtk_init(argc, argv);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv);
    {
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        LEInteractor::s_kernel_fcn = &four_point_kernel;
        LEInteractor::s_kernel_fcn_stencil_size = 4;
        const std::vector<std::string> kernels = { "PIECEWISE_CONSTANT",
                                                   "DISCONTINUOUS_LINEAR",
                                                   "PIECEWISE_LINEAR",
                                                   "PIECEWISE_CUBIC",
                                                   "IB_3",
                                                   "IB_4",
                                                   "IB_4_W8",
                                                   "IB_5",
                                                   "IB_6",
                                                   "BSPLINE_3",
                                                   "BSPLINE_4",
                                                   "BSPLINE_5",
                                                   "BSPLINE_6",
                                                   "COMPOSITE_BSPLINE_23",
                                                   "COMPOSITE_BSPLINE_32",
                                                   "COMPOSITE_BSPLINE_34",
                                                   "COMPOSITE_BSPLINE_43",
                                                   "COMPOSITE_BSPLINE_45",
                                                   "COMPOSITE_BSPLINE_54",
                                                   "COMPOSITE_BSPLINE_56",
                                                   "COMPOSITE_BSPLINE_65",
                                                   "USER_DEFINED" };
        int ghost_width = 0;
        for (const std::string& kernel : kernels)
        {
            ghost_width = std::max(ghost_width, LEInteractor::getMinimumGhostWidth(kernel));
        }

        Pointer<CartesianGridGeometry<NDIM>> grid_geom = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Domain domain;
        domain.N = input_db->getInteger("N");
        for (int d = 0; d < NDIM; ++d)
        {
            domain.x_lo[d] = grid_geom->getXLower()[d];
            domain.x_up[d] = grid_geom->getXUpper()[d];
            domain.periodic[d] = grid_geom->getPeriodicShift()(d) != 0;
        }
        const std::vector<double> X = build_points(domain);
        int num_outside = 0;
        for (std::size_t k = 0; k < X.size() / NDIM; ++k)
        {
            num_outside += inside_domain(domain, &X[NDIM * k]) ? 0 : 1;
        }
        pout << "points = " << X.size() / NDIM << ", of which outside the domain = " << num_outside << '\n';

        // One patch; patches that meet at the middle face; and patches that meet one ghost width from the lower end.
        // SAMRAI requires each patch to be at least the ghost width wide.
        if (domain.N / 2 < ghost_width || domain.N - ghost_width < ghost_width)
        {
            TBOX_ERROR("N = " << domain.N << " is too small for patches of at least " << ghost_width << " cells\n");
        }
        std::vector<Layout> layouts = { { "one patch", nullptr, {} },
                                        { "halves", nullptr, {} },
                                        { "one ghost width", nullptr, {} } };
        const std::vector<std::vector<int>> splits = { {}, { domain.N / 2 }, { ghost_width } };
        for (std::size_t l = 0; l < layouts.size(); ++l)
        {
            layouts[l].level = make_level(grid_geom, domain, splits[l]);
            for (PatchLevel<NDIM>::Iterator p(layouts[l].level); p; p++)
            {
                layouts[l].idx_data.push_back(
                    build_index_data(layouts[l].level->getPatch(p()), grid_geom, domain, X, ghost_width));
            }
        }

        std::vector<Centering> centerings = { Centering::CELL, Centering::NODE, Centering::SIDE };
        if (NDIM == 3)
        {
            centerings.push_back(Centering::EDGE);
        }
        for (const std::string& kernel : kernels)
        {
            for (const Centering centering : centerings)
            {
                check(layouts, domain, X, kernel, centering, ghost_width);
            }
        }
    }
    return 0;
}
