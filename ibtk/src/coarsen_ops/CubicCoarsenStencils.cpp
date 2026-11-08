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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <tbox/Array.h>

#include <Box.h>
#include <BoxArray.h>
#include <GridGeometry.h>

#include <limits>

#include "./CubicCoarsenStencils.h"

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// A one-dimensional stencil, in fine offsets from the lower of the two fine cells nearest to the coarse center, and
// the range of coarse neighbor offsets that it involves along a RATIO_2 axis.
struct Stencil1D
{
    int size;
    int offsets[4];
    double weights[4];
    int lowest_neighbor;
    int highest_neighbor;
    int num_averaged;
    int num_one_sided;
};

// The centered cubic, the one-sided quadratics extending upward and downward, the two-point average, and the fine
// side of a NORMAL axis.
constexpr int CENTERED = 0;
constexpr int AVERAGE = 3;
constexpr int IDENTITY = 4;
constexpr Stencil1D STENCILS[] = { { 4, { -1, 0, 1, 2 }, { -0.0625, 0.5625, 0.5625, -0.0625 }, -1, 1, 0, 0 },
                                   { 3, { 0, 1, 2, 0 }, { 0.375, 0.75, -0.125, 0.0 }, 0, 1, 0, 1 },
                                   { 3, { -1, 0, 1, 0 }, { -0.125, 0.75, 0.375, 0.0 }, -1, 0, 0, 1 },
                                   { 2, { 0, 1, 0, 0 }, { 0.5, 0.5, 0.0, 0.0 }, 0, 0, 1, 0 },
                                   { 1, { 0, 0, 0, 0 }, { 1.0, 0.0, 0.0, 0.0 }, 0, 0, 0, 0 } };

// The number of coarse neighbors of a position, itself included.
constexpr int NUM_NEIGHBORS = (NDIM == 2) ? 9 : 27;

// Return the index into the neighbor coverage of the coarse offset with components -1, 0, or 1.
int
neighbor_index(const std::array<int, NDIM>& offset)
{
    int index = 0;
    int stride = 1;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        index += (offset[d] + 1) * stride;
        stride *= 3;
    }
    return index;
}
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

LevelCoverage::LevelCoverage(Pointer<PatchLevel<NDIM>> level)
    : d_box_tree(level->getBoxTree()),
      d_periodic_shift(level->getGridGeometry()->getPeriodicShift(level->getRatio())),
      d_domain_lower(level->getPhysicalDomain().getBoundingBox().lower())
{
    // intentionally blank
    return;
} // LevelCoverage

bool
LevelCoverage::covers(Index<NDIM> cell) const
{
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (d_periodic_shift(d) > 0)
        {
            const int offset = cell(d) - d_domain_lower(d);
            cell(d) = d_domain_lower(d) + ((offset % d_periodic_shift(d)) + d_periodic_shift(d)) % d_periodic_shift(d);
        }
    }
    tbox::Array<int> indices;
    d_box_tree->findOverlapIndices(indices, Box<NDIM>(cell, cell));
    return indices.size() > 0;
} // covers

CubicCoarsenStencil::CubicCoarsenStencil(const std::array<CubicCoarsenAxis, NDIM>& axes,
                                         const std::function<bool(const IntVector<NDIM>&)>& covered)
{
    // Find which coarse neighbors the fine level covers, along the RATIO_2 axes.
    std::array<bool, NUM_NEIGHBORS> neighbor_covered;
    neighbor_covered.fill(true);
    for (int k = 0; k < NUM_NEIGHBORS; ++k)
    {
        IntVector<NDIM> offset(0);
        bool is_neighbor = false;
        bool along_ratio_2 = true;
        int rest = k;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            offset(d) = rest % 3 - 1;
            rest /= 3;
            if (offset(d) != 0)
            {
                is_neighbor = true;
                along_ratio_2 = along_ratio_2 && axes[d] == CubicCoarsenAxis::RATIO_2;
            }
        }
        if (is_neighbor && along_ratio_2)
        {
            neighbor_covered[k] = covered(offset);
        }
    }

    // Return whether a product uses only covered coarse neighbors.
    const auto uses_only_covered_neighbors = [&](const std::array<int, NDIM>& choice)
    {
        std::array<int, NDIM> lower, upper, offset;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            const bool ratio_2 = axes[d] == CubicCoarsenAxis::RATIO_2;
            lower[d] = ratio_2 ? STENCILS[choice[d]].lowest_neighbor : 0;
            upper[d] = ratio_2 ? STENCILS[choice[d]].highest_neighbor : 0;
            offset[d] = lower[d];
        }
        while (true)
        {
            if (!neighbor_covered[neighbor_index(offset)])
            {
                return false;
            }
            unsigned int d = 0;
            while (d < NDIM && offset[d] == upper[d])
            {
                offset[d] = lower[d];
                ++d;
            }
            if (d == NDIM)
            {
                return true;
            }
            ++offset[d];
        }
    };

    // Visit every product of the stencils of the RATIO_2 axes and select the valid ones with the lowest score: the
    // number of averaged axes first, then the number of one-sided axes.
    std::array<int, NDIM> choice;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        choice[d] = axes[d] == CubicCoarsenAxis::NORMAL ? IDENTITY : CENTERED;
    }
    int best_score = std::numeric_limits<int>::max();
    while (true)
    {
        int score = 0;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (axes[d] == CubicCoarsenAxis::RATIO_2)
            {
                score += (NDIM + 1) * STENCILS[choice[d]].num_averaged + STENCILS[choice[d]].num_one_sided;
            }
        }
        if (score <= best_score && uses_only_covered_neighbors(choice))
        {
            if (score < best_score)
            {
                best_score = score;
                d_num_products = 0;
            }
            d_products[d_num_products++] = choice;
        }
        unsigned int d = 0;
        while (d < NDIM && (axes[d] != CubicCoarsenAxis::RATIO_2 || choice[d] == AVERAGE))
        {
            if (axes[d] == CubicCoarsenAxis::RATIO_2)
            {
                choice[d] = CENTERED;
            }
            ++d;
        }
        if (d == NDIM)
        {
            break;
        }
        ++choice[d];
    }
    return;
} // CubicCoarsenStencil

double
CubicCoarsenStencil::apply(const ArrayData<NDIM, double>& fine_data,
                           const Index<NDIM>& fine_base_index,
                           const int depth) const
{
    double sum = 0.0;
    for (int p = 0; p < d_num_products; ++p)
    {
        const std::array<int, NDIM>& choice = d_products[p];
        std::array<int, NDIM> k;
        k.fill(0);
        double value = 0.0;
        while (true)
        {
            Index<NDIM> i = fine_base_index;
            double weight = 1.0;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                i(d) += STENCILS[choice[d]].offsets[k[d]];
                weight *= STENCILS[choice[d]].weights[k[d]];
            }
            value += weight * fine_data(i, depth);
            unsigned int d = 0;
            while (d < NDIM && k[d] + 1 == STENCILS[choice[d]].size)
            {
                k[d] = 0;
                ++d;
            }
            if (d == NDIM)
            {
                break;
            }
            ++k[d];
        }
        sum += value;
    }
    return sum / static_cast<double>(d_num_products);
} // apply

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
