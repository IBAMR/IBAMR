// ---------------------------------------------------------------------
//
// Copyright (c) 2023 - 2023 by the IBAMR developers
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

#include <ibtk/IBTK_MPI.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/MarkerPatchHierarchy.h>

#include <CartesianGridGeometry.h>

#include <ibtk/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
MarkerPatch::MarkerPatch(const Box<NDIM>& patch_box,
                         const std::vector<Box<NDIM> >& nonoverlapping_patch_boxes,
                         const Pointer<CartesianGridGeometry<NDIM> >& grid_geom,
                         const IntVector<NDIM>& ratio)
    : d_patch_box(patch_box), d_nonoverlapping_patch_boxes(nonoverlapping_patch_boxes)
{
    for (const Box<NDIM>& box : nonoverlapping_patch_boxes)
    {
        TBOX_ASSERT(patch_box.contains(box));
    }
    d_domain_box = Box<NDIM>::refine(grid_geom->getPhysicalDomain()[0], ratio);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_x_lo[d] = grid_geom->getXLower()[d];
        d_x_up[d] = grid_geom->getXUpper()[d];
        d_dx[d] = grid_geom->getDx()[d] / ratio(d);
    }
}

void
MarkerPatch::insert(const int& index, const IBTK::Point& position, const IBTK::Vector& velocity)
{
    TBOX_ASSERT(index >= 0);
    const auto p = std::lower_bound(d_indices.begin(), d_indices.end(), index) - d_indices.begin();
    d_indices.insert(d_indices.begin() + p, index);
    d_positions.insert(d_positions.begin() + NDIM * p, position.data(), position.data() + position.size());
    d_velocities.insert(d_velocities.begin() + NDIM * p, velocity.data(), velocity.data() + velocity.size());
}

bool
MarkerPatch::contains(const IBTK::Point& position) const
{
    const auto index = IndexUtilities::getCellIndex(
        position.data(), d_x_lo.data(), d_x_up.data(), d_dx.data(), d_domain_box.lower(), d_domain_box.upper());
    if (!d_patch_box.contains(index)) return false;

    for (const auto& box : d_nonoverlapping_patch_boxes)
        if (box.contains(index)) return true;
    return false;
}

std::tuple<std::vector<int>, EigenAlignedVector<IBTK::Point>, EigenAlignedVector<IBTK::Vector> >
MarkerPatch::prune()
{
    std::vector<int> indices;
    EigenAlignedVector<IBTK::Point> positions;
    EigenAlignedVector<IBTK::Vector> velocities;
    // Prune from the back to the front to avoid moving stored markers
    for (int index = size() - 1; index >= 0; index--)
    {
        IBTK::Point point;
        IBTK::Point velocity;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            point[d] = d_positions[index * NDIM + d];
            velocity[d] = d_velocities[index * NDIM + d];
        }
        if (!contains(point))
        {
            // emplace so that things remain sorted
            indices.insert(indices.begin(), d_indices[index]);
            d_indices.erase(d_indices.begin() + index);
            auto p0 = d_positions.begin() + NDIM * index;
            auto p1 = d_positions.begin() + NDIM * (index + 1);
            positions.emplace(positions.begin(), &*p0);
            d_positions.erase(p0, p1);
            auto v0 = d_velocities.begin() + NDIM * index;
            auto v1 = d_velocities.begin() + NDIM * (index + 1);
            velocities.emplace(velocities.begin(), &*v0);
            d_velocities.erase(v0, v1);
        }
    }

    return std::make_tuple(std::move(indices), std::move(positions), std::move(velocities));
}

std::tuple<int, IBTK::Point, IBTK::Vector>
MarkerPatch::operator[](const unsigned int local_index) const
{
#ifndef NDEBUG
    TBOX_ASSERT(local_index < size());
#endif
    IBTK::Point position(&d_positions[local_index * NDIM]);
    IBTK::Point velocity(&d_velocities[local_index * NDIM]);
    return std::make_tuple(d_indices[local_index], position, velocity);
}

std::size_t
MarkerPatch::size() const
{
#ifndef NDEBUG
    TBOX_ASSERT(d_positions.size() == d_indices.size() * NDIM);
    TBOX_ASSERT(d_velocities.size() == d_indices.size() * NDIM);
#endif
    return d_indices.size();
}
} // namespace IBTK
