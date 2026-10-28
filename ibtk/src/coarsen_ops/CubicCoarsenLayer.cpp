// ---------------------------------------------------------------------
//
// Copyright (c) 2026 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <tbox/Utilities.h>

#include <Index.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

#include "./CubicCoarsenLayer.h"

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Extrapolation weights to the position adjacent to a row of n values (nearest first).
const double EXTRAPOLATION_WEIGHTS[4][4] = { { 1.0, 0.0, 0.0, 0.0 },
                                             { 2.0, -1.0, 0.0, 0.0 },
                                             { 3.0, -3.0, 1.0, 0.0 },
                                             { 4.0, -6.0, 4.0, -1.0 } };

constexpr int MAX_BLOCK_SIZE = 4;
// Number of positions of the largest block: MAX_BLOCK_SIZE ^ NDIM.
constexpr int MAX_BLOCK_POSITIONS = (NDIM == 2) ? 16 : 64;
// Number of nonempty subsets of the axes: 2 ^ NDIM - 1.
constexpr int MAX_SUBSETS = (NDIM == 2) ? 3 : 7;

// Element strides of an array along each axis (zero beyond NDIM) and between depths.
struct Layout
{
    std::array<std::ptrdiff_t, 3> stride = { { 0, 0, 0 } };
    std::ptrdiff_t depth_stride = 0;
    int depth = 0;
};

Layout
get_layout(const ArrayData<NDIM, double>& data)
{
    Layout layout;
    std::ptrdiff_t s = 1;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        layout.stride[d] = s;
        s *= data.getBox().numberCells(d);
    }
    layout.depth_stride = s;
    layout.depth = data.getDepth();
    return layout;
}

// Return the number of axes in the subset of n_axes axes that the number set encodes: axis k belongs to the subset if
// digit k of set in base 2 is 1.
int
set_size(const int set, const int n_axes)
{
    int size = 0;
    int rest = set;
    for (int k = 0; k < n_axes; ++k)
    {
        size += rest % 2;
        rest /= 2;
    }
    return size;
}

// Return true if the box contains a NaN.
[[maybe_unused]] bool
contains_nan(const ArrayData<NDIM, double>& data, const Box<NDIM>& box)
{
    for (Box<NDIM>::Iterator b(box); b; b++)
    {
        for (int depth = 0; depth < data.getDepth(); ++depth)
        {
            if (std::isnan(data(*b, depth)))
            {
                return true;
            }
        }
    }
    return false;
}

// The extrapolation block of the positions of a region along a subset of the axes beyond the box: the offsets from a
// position to the values of the block and their weights, in the order in which they are summed.
struct Block
{
    bool all_inside = false; // whether every value of the block lies in the box
    bool built = false;      // whether the offsets and weights have been computed
    int n_axes = 0;
    std::array<int, NDIM> n_values;
    std::array<std::ptrdiff_t, NDIM> steps;
    std::ptrdiff_t first = 0; // offset to the value nearest to the position
    int n_positions = 0;
    std::array<std::ptrdiff_t, MAX_BLOCK_POSITIONS> offsets;
    std::array<double, MAX_BLOCK_POSITIONS> weights;
};

// Compute the offsets and weights of a block by walking through it with an odometer over the offsets
// (1, ..., n_values) along its axes.
void
build_block(Block& block)
{
    std::array<int, NDIM> offset;
    offset.fill(1);
    block.n_positions = 0;
    while (true)
    {
        std::ptrdiff_t o = 0;
        double w = 1.0;
        for (int k = 0; k < block.n_axes; ++k)
        {
            o += offset[k] * block.steps[k];
            w *= EXTRAPOLATION_WEIGHTS[block.n_values[k] - 1][offset[k] - 1];
        }
        block.offsets[block.n_positions] = o;
        block.weights[block.n_positions] = w;
        ++block.n_positions;

        int k = 0;
        while (k < block.n_axes && offset[k] == block.n_values[k])
        {
            offset[k] = 1;
            ++k;
        }
        if (k == block.n_axes)
        {
            break;
        }
        ++offset[k];
    }
    block.built = true;
    return;
}

// Extrapolate along one axis, from n values that start one step from p.
template <int N>
inline double
extrapolate(const double* const p, const std::ptrdiff_t step)
{
    double value = 0.0;
    for (int k = 0; k < N; ++k)
    {
        value += EXTRAPOLATION_WEIGHTS[N - 1][k] * p[(k + 1) * step];
    }
    return value;
}

// Replace the NaN values of a region that lies beyond the box along one axis only.
template <int N>
void
extrapolate_face(double* const base,
                 const Layout& layout,
                 const std::array<int, 3>& lo,
                 const std::array<int, 3>& hi,
                 const std::ptrdiff_t step,
                 [[maybe_unused]] const ArrayData<NDIM, double>& data,
                 [[maybe_unused]] const Box<NDIM>& box)
{
    for (int k = lo[2]; k <= hi[2]; ++k)
    {
        for (int j = lo[1]; j <= hi[1]; ++j)
        {
            double* row = base + k * layout.stride[2] + j * layout.stride[1];
            for (int i = lo[0]; i <= hi[0]; ++i)
            {
                double* const p = row + i;
                if (std::isnan(*p))
                {
                    for (int depth = 0; depth < layout.depth; ++depth)
                    {
                        double* const q = p + depth * layout.depth_stride;
                        *q = extrapolate<N>(q, step);
                        TBOX_ASSERT(!std::isnan(*q) || contains_nan(data, box));
                    }
                }
            }
        }
    }
    return;
}

// Replace the NaN values of a region that lies beyond the box along the axes s_axes (ascending). A position takes the
// extrapolation along the smallest set of these axes (the sets of one size in the order of their numbers) for which
// every value of the block is a valid one, and in particular every value of a block along all of them is.
void
extrapolate_region(double* const base,
                   const Layout& layout,
                   const std::array<int, 3>& lo,
                   const std::array<int, 3>& hi,
                   const std::array<int, NDIM>& s_axes,
                   const std::array<int, NDIM>& s_dirs,
                   const int n_s,
                   const Box<NDIM>& box,
                   [[maybe_unused]] const ArrayData<NDIM, double>& data)
{
    int n_subsets = 1;
    for (int k = 0; k < n_s; ++k)
    {
        n_subsets *= 2;
    }

    // Describe the blocks of the nonempty subsets, by size and then by number. Build a block when a position needs it.
    std::array<Block, MAX_SUBSETS> blocks;
    int n_blocks = 0;
    for (int size = 1; size <= n_s; ++size)
    {
        for (int set = 1; set < n_subsets; ++set)
        {
            if (set_size(set, n_s) != size)
            {
                continue;
            }
            Block& block = blocks[n_blocks++];
            block.all_inside = (size == n_s);
            int rest = set;
            for (int k = 0; k < n_s; ++k)
            {
                if (rest % 2 == 1)
                {
                    const int d = s_axes[k];
                    block.n_values[block.n_axes] = std::min(MAX_BLOCK_SIZE, box.numberCells(d));
                    block.steps[block.n_axes] = s_dirs[k] * layout.stride[d];
                    block.first += block.steps[block.n_axes];
                    ++block.n_axes;
                }
                rest /= 2;
            }
        }
    }

    for (int k = lo[2]; k <= hi[2]; ++k)
    {
        for (int j = lo[1]; j <= hi[1]; ++j)
        {
            double* row = base + k * layout.stride[2] + j * layout.stride[1];
            for (int i = lo[0]; i <= hi[0]; ++i)
            {
                double* const p = row + i;
                if (!std::isnan(*p))
                {
                    continue;
                }
                bool found = false;
                for (int b = 0; b < n_blocks && !found; ++b)
                {
                    Block& block = blocks[b];
                    if (!block.all_inside && std::isnan(p[block.first]))
                    {
                        continue;
                    }
                    if (!block.built)
                    {
                        build_block(block);
                    }
                    bool valid = true;
                    if (!block.all_inside)
                    {
                        for (int m = 0; m < block.n_positions && valid; ++m)
                        {
                            valid = !std::isnan(p[block.offsets[m]]);
                        }
                    }
                    if (!valid)
                    {
                        continue;
                    }
                    for (int depth = 0; depth < layout.depth; ++depth)
                    {
                        double* const q = p + depth * layout.depth_stride;
                        double value = 0.0;
                        for (int m = 0; m < block.n_positions; ++m)
                        {
                            value += block.weights[m] * q[block.offsets[m]];
                        }
                        *q = value;
                        TBOX_ASSERT(!std::isnan(value) || contains_nan(data, box));
                    }
                    found = true;
                }
                TBOX_ASSERT(found);
                (void)found;
            }
        }
    }
    return;
}
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
mark_cubic_coarsen_layer(ArrayData<NDIM, double>& data, const Box<NDIM>& box)
{
    TBOX_ASSERT(data.getBox().contains(box));
    const double marker = std::numeric_limits<double>::quiet_NaN();
    Box<NDIM> rest = data.getBox();
    for (int d = NDIM - 1; d >= 0; --d)
    {
        if (rest.lower(d) < box.lower(d))
        {
            Box<NDIM> slab = rest;
            slab.upper(d) = box.lower(d) - 1;
            data.fillAll(marker, slab);
        }
        if (rest.upper(d) > box.upper(d))
        {
            Box<NDIM> slab = rest;
            slab.lower(d) = box.upper(d) + 1;
            data.fillAll(marker, slab);
        }
        rest.lower(d) = box.lower(d);
        rest.upper(d) = box.upper(d);
    }
    return;
}

void
extrapolate_cubic_coarsen_layer(ArrayData<NDIM, double>& data, const Box<NDIM>& box, const IntVector<NDIM>& width)
{
    std::array<int, NDIM> axes;
    int n_axes = 0;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(width(d) == 0 || width(d) == 1);
        if (width(d) == 1)
        {
            axes[n_axes++] = d;
        }
    }
    if (n_axes == 0)
    {
        return;
    }
    TBOX_ASSERT(data.getBox().contains(Box<NDIM>::grow(box, width)));

    const Layout layout = get_layout(data);
    // base_zero points to the element with index zero along every axis; the region loops use indices directly.
    double* base_zero = data.getPointer(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        base_zero -= data.getBox().lower(d) * layout.stride[d];
    }

    int n_sets = 1;
    for (int k = 0; k < n_axes; ++k)
    {
        n_sets *= 2;
    }

    // Visit the regions from the highest number of axes beyond the box to the lowest, so that an extrapolation never
    // reads a value that an extrapolation of this call has written.
    for (int size = n_axes; size >= 1; --size)
    {
        for (int set = 1; set < n_sets; ++set)
        {
            if (set_size(set, n_axes) != size)
            {
                continue;
            }
            std::array<int, NDIM> s_axes;
            int rest = set;
            int n_s = 0;
            for (int k = 0; k < n_axes; ++k)
            {
                if (rest % 2 == 1)
                {
                    s_axes[n_s++] = axes[k];
                }
                rest /= 2;
            }
            int n_signs = 1;
            for (int k = 0; k < size; ++k)
            {
                n_signs *= 2;
            }
            for (int signs = 0; signs < n_signs; ++signs)
            {
                std::array<int, 3> lo = { { 0, 0, 0 } };
                std::array<int, 3> hi = { { 0, 0, 0 } };
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    lo[d] = box.lower(d);
                    hi[d] = box.upper(d);
                }
                std::array<int, NDIM> s_dirs;
                int sign_rest = signs;
                for (int k = 0; k < size; ++k)
                {
                    const int d = s_axes[k];
                    if (sign_rest % 2 == 0)
                    {
                        lo[d] = hi[d] = box.lower(d) - 1;
                        s_dirs[k] = 1;
                    }
                    else
                    {
                        lo[d] = hi[d] = box.upper(d) + 1;
                        s_dirs[k] = -1;
                    }
                    sign_rest /= 2;
                }
                if (size == 1)
                {
                    const int d = s_axes[0];
                    const std::ptrdiff_t step = s_dirs[0] * layout.stride[d];
                    switch (std::min(MAX_BLOCK_SIZE, box.numberCells(d)))
                    {
                    case 1:
                        extrapolate_face<1>(base_zero, layout, lo, hi, step, data, box);
                        break;
                    case 2:
                        extrapolate_face<2>(base_zero, layout, lo, hi, step, data, box);
                        break;
                    case 3:
                        extrapolate_face<3>(base_zero, layout, lo, hi, step, data, box);
                        break;
                    default:
                        extrapolate_face<4>(base_zero, layout, lo, hi, step, data, box);
                        break;
                    }
                }
                else
                {
                    extrapolate_region(base_zero, layout, lo, hi, s_axes, s_dirs, size, box, data);
                }
            }
        }
    }
    return;
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
