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
namespace
{
std::vector<std::vector<hier::Box<NDIM> > >
compute_nonoverlapping_patch_boxes(const Pointer<BasePatchLevel<NDIM> >& c_level,
                                   const Pointer<BasePatchLevel<NDIM> >& f_level)
{
    const Pointer<PatchLevel<NDIM> > coarse_level = c_level;
    const Pointer<PatchLevel<NDIM> > fine_level = f_level;
    TBOX_ASSERT(coarse_level);
    TBOX_ASSERT(fine_level);
    TBOX_ASSERT(coarse_level->getLevelNumber() + 1 == fine_level->getLevelNumber());

    const IntVector<NDIM> ratio = fine_level->getRatioToCoarserLevel();

    // Get all (including those not on this processor) fine-level boxes:
    BoxList<NDIM> finer_box_list;
    long combined_size = 0;
    for (int i = 0; i < fine_level->getNumberOfPatches(); ++i)
    {
        Box<NDIM> patch_box = fine_level->getBoxForPatch(i);
        patch_box.coarsen(ratio);
        combined_size += patch_box.size();
        finer_box_list.addItem(patch_box);
    }
    finer_box_list.simplifyBoxes();

    // Remove said boxes from each coarse-level patch:
    const auto rank = IBTK_MPI::getRank();
    std::vector<std::vector<Box<NDIM> > > result;
    long coarse_size = 0;
    for (int i = 0; i < coarse_level->getNumberOfPatches(); ++i)
    {
        BoxList<NDIM> coarse_box_list;
        coarse_box_list.addItem(coarse_level->getBoxForPatch(i));
        coarse_size += coarse_box_list.getFirstItem().size();
        coarse_box_list.removeIntersections(finer_box_list);

        const bool patch_is_local = rank == coarse_level->getMappingForPatch(i);
        if (patch_is_local) result.emplace_back();
        typename tbox::List<Box<NDIM> >::Iterator it(coarse_box_list);
        while (it)
        {
            if (patch_is_local) result.back().push_back(*it);
            combined_size += (*it).size();
            it++;
        }
    }

    TBOX_ASSERT(coarse_size == combined_size);

    return result;
}
} // namespace
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

MarkerPatchHierarchy::MarkerPatchHierarchy(const std::string& name,
                                           Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                           const EigenAlignedVector<IBTK::Point>& positions,
                                           const EigenAlignedVector<IBTK::Point>& velocities)
    : d_object_name(name), d_hierarchy(patch_hierarchy)
{
    reinit(positions, velocities);
}

void
MarkerPatchHierarchy::reinit(const EigenAlignedVector<IBTK::Point>& positions,
                             const EigenAlignedVector<IBTK::Point>& velocities)
{
    TBOX_ASSERT(positions.size() == velocities.size());
    const auto rank = IBTK_MPI::getRank();
    const Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    unsigned int num_emplaced_markers = 0;
    std::vector<bool> marker_emplaced(positions.size());

    auto insert_markers = [&](MarkerPatch& marker_patch) {
        for (unsigned int k = 0; k < positions.size(); ++k)
        {
            if (!marker_emplaced[k] && marker_patch.contains(positions[k]))
            {
                marker_emplaced[k] = true;
                marker_patch.insert(k, positions[k], velocities[k]);
                ++num_emplaced_markers;
            }
        }
    };

    d_marker_patches.clear();
    d_marker_patches.resize(d_hierarchy->getFinestLevelNumber() + 1);

    // Assign particles to levels in a top-down way:
    //
    // 1. Compute a set of boxes for the patch which do not intersect any
    //    finer patch level.
    // 2. Emplace particles in the vector corresponding to the local index of a
    //    patch in the present level.
    for (int ln = d_hierarchy->getFinestLevelNumber(); ln >= 0; --ln)
    {
        Pointer<PatchLevel<NDIM> > current_level = d_hierarchy->getPatchLevel(ln);
        Pointer<PatchLevel<NDIM> > finer_level =
            ln == d_hierarchy->getFinestLevelNumber() ? nullptr : d_hierarchy->getPatchLevel(ln + 1);
        const IntVector<NDIM>& ratio = current_level->getRatio();

        // If there is no finer level then each Patch has exactly one
        // nonoverlapping box
        if (!finer_level)
        {
            for (int i = 0; i < current_level->getNumberOfPatches(); ++i)
            {
                if (rank == current_level->getMappingForPatch(i))
                {
                    const Box<NDIM> box = current_level->getPatch(i)->getBox();
                    d_marker_patches[ln].emplace_back(box, std::vector<Box<NDIM> >{ box }, grid_geom, ratio);
                    insert_markers(d_marker_patches[ln].back());
                }
            }
        }
        // otherwise we need to subtract off the boxes on the finer level first.
        else
        {
            const std::vector<std::vector<hier::Box<NDIM> > > nonoverlapping_patch_boxes =
                compute_nonoverlapping_patch_boxes(current_level, finer_level);
            unsigned int local_num = 0;
            for (int i = 0; i < current_level->getNumberOfPatches(); ++i)
            {
                if (rank == current_level->getMappingForPatch(i))
                {
                    d_marker_patches[ln].emplace_back(
                        current_level->getPatch(i)->getBox(), nonoverlapping_patch_boxes[local_num], grid_geom, ratio);
                    insert_markers(d_marker_patches[ln].back());
                    ++local_num;
                }
            }
        }
    }
    num_emplaced_markers = IBTK_MPI::sumReduction(num_emplaced_markers);
    TBOX_ASSERT(num_emplaced_markers == positions.size());
}

const MarkerPatch&
MarkerPatchHierarchy::getMarkerPatch(const int ln, const int local_patch_num) const
{
    TBOX_ASSERT(ln < static_cast<int>(d_marker_patches.size()));
    TBOX_ASSERT(local_patch_num < static_cast<int>(d_marker_patches[ln].size()));
    return d_marker_patches[ln][local_patch_num];
}

} // namespace IBTK
