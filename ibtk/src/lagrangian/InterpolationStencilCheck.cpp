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

#include <ibtk/IndexUtilities.h>
#include <ibtk/InterpolationStencilCheck.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/LIndexSetData.h>
#include <ibtk/LNode.h>
#include <ibtk/PhysicalBoundaryUtilities.h>

#include <tbox/Array.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <Box.h>
#include <BoxArray.h>
#include <BoxList.h>
#include <BoxTree.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <EdgeData.h>
#include <GridGeometry.h>
#include <IntVector.h>
#include <NodeData.h>
#include <Patch.h>
#include <PatchData.h>
#include <PatchDescriptor.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SideData.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <sstream>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// The possible staggerings of data, as the cells of a unit box: the entry along an axis is 1 for data on the lower
// faces of the cells along that axis, and 0 for data at the cell centers. The position in the box indexes the lists.
const Box<NDIM> STAGGERINGS(hier::Index<NDIM>(0), hier::Index<NDIM>(1));

// The staggering of each component of q_data.
std::vector<IntVector<NDIM>>
staggerings(const PatchData<NDIM>& q_data)
{
    std::vector<IntVector<NDIM>> result;
    if (dynamic_cast<const CellData<NDIM, double>*>(&q_data))
    {
        result.emplace_back(0);
    }
    else if (dynamic_cast<const NodeData<NDIM, double>*>(&q_data))
    {
        result.emplace_back(1);
    }
    else if (dynamic_cast<const SideData<NDIM, double>*>(&q_data))
    {
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            result.emplace_back(0);
            result.back()(axis) = 1;
        }
    }
    else if (dynamic_cast<const EdgeData<NDIM, double>*>(&q_data))
    {
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            result.emplace_back(1);
            result.back()(axis) = 0;
        }
    }
    else
    {
        TBOX_ERROR("InterpolationStencilCheck::check():\n"
                   << "  unsupported patch data type\n"
                   << "  use cell-, node-, side-, or edge-centered data with double precision values\n");
    }
    return result;
}

// The nonzero offsets from zero to the staggering. Data of the staggering at index j touch the cell j and the cells
// j - offset. Because no caller passes node- or edge-centered data, the lists for those data are unused.
std::vector<IntVector<NDIM>>
touching_offsets(const IntVector<NDIM>& staggering)
{
    std::vector<IntVector<NDIM>> offsets;
    for (Box<NDIM>::Iterator b(Box<NDIM>(hier::Index<NDIM>(0), hier::Index<NDIM>(staggering))); b; b++)
    {
        const IntVector<NDIM> offset(b());
        if (offset != IntVector<NDIM>(0)) offsets.push_back(offset);
    }
    return offsets;
}

// The number of nodes along each axis that the kernel function uses. This number is less than the stencil size for the
// kernels of odd order.
int
num_nodes_used(const std::string& kernel_fcn)
{
    if (kernel_fcn == "IB_3" || kernel_fcn == "BSPLINE_3" || kernel_fcn == "COMPOSITE_BSPLINE_23" ||
        kernel_fcn == "COMPOSITE_BSPLINE_32")
    {
        return 3;
    }
    if (kernel_fcn == "IB_5" || kernel_fcn == "BSPLINE_5")
    {
        return 5;
    }
    return LEInteractor::getStencilSize(kernel_fcn);
}

// The nodes of the stencil of the point X, which are the n nearest to X along each axis. A node lies at the center of
// the cell with the same index along an axis, or on its lower face along an axis with staggering 1.
Box<NDIM>
stencil_box(const std::array<double, NDIM>& X, const Patch<NDIM>& patch, const IntVector<NDIM>& staggering, const int n)
{
    const Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    Box<NDIM> box;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        const double node_offset = staggering(d) ? 0.0 : 0.5;
        const double position =
            (X[d] - pgeom->getXLower()[d]) / pgeom->getDx()[d] + patch.getBox().lower(d) - node_offset;
        box.lower(d) = static_cast<int>(std::floor(position - 0.5 * n)) + 1;
        box.upper(d) = box.lower(d) + n - 1;
    }
    return box;
}

// Emit an error if the stencil includes a flagged cell. The list outside contains the flagged cells outside the
// physical domain beside a boundary that the patch does not touch.
void
check_stencil(const std::string& caller,
              const int ln,
              const Box<NDIM>& patch_box,
              const std::array<double, NDIM>& X,
              const Box<NDIM>& stencil,
              const BoxList<NDIM>& flagged,
              const BoxList<NDIM>& outside)
{
    for (const BoxList<NDIM>* list : { &outside, &flagged })
    {
        for (BoxList<NDIM>::Iterator b(*list); b; b++)
        {
            const Box<NDIM> overlap = b() * stencil;
            if (overlap.empty()) continue;
            std::ostringstream position;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                position << (d == 0 ? "(" : ", ") << X[d];
            }
            if (list == &outside)
            {
                TBOX_ERROR(caller << "():\n"
                                  << "  the stencil of the point at " << position.str() << ") on level " << ln
                                  << " includes the cell " << overlap.lower()
                                  << " outside the physical domain, beside a boundary that the patch with box "
                                  << patch_box << " does not touch\n"
                                  << "  patches must touch a physical boundary or be at least as far from it as the "
                                     "ghost cell width\n");
            }
            TBOX_ERROR(caller << "():\n"
                              << "  the stencil of the point at " << position.str() << ") on level " << ln
                              << " includes the cell " << overlap.lower() << ", where the data are unset\n"
                              << "  enlarge the tag buffer or the refined region of level " << ln << "\n");
        }
    }
}
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
InterpolationStencilCheck::reset(const Pointer<PatchHierarchy<NDIM>> hierarchy, const bool fill_sets_uncovered)
{
    d_max_ghost_width = VariableDatabase<NDIM>::getDatabase()->getPatchDescriptor()->getMaxGhostWidth();
    const Pointer<GridGeometry<NDIM>> grid_geom = hierarchy->getGridGeometry();
    d_flagged_boxes.clear();
    d_flagged_boxes.resize(hierarchy->getNumberOfLevels());
    d_levels.clear();
    for (int ln = 0; ln < hierarchy->getNumberOfLevels(); ++ln)
    {
        d_levels.push_back(hierarchy->getPatchLevel(ln));
    }
    for (int ln = 0; ln < hierarchy->getNumberOfLevels(); ++ln)
    {
        // The cells of the level and its periodic images, and those of the domain.
        const Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        const IntVector<NDIM>& ratio = level->getRatio();
        const BoxList<NDIM> domain(PhysicalBoundaryUtilities::getPeriodicPhysicalDomain(grid_geom, ratio));
        BoxList<NDIM> covered(level->getBoxes());
        const std::vector<IntVector<NDIM>> image_shifts =
            PhysicalBoundaryUtilities::getPeriodicImageShifts(grid_geom->getPeriodicShift(ratio));
        for (const IntVector<NDIM>& shift : image_shifts)
        {
            BoxList<NDIM> image(level->getBoxes());
            image.shift(shift);
            covered.catenateItems(image);
        }
        const BoxArray<NDIM> covered_boxes(covered);
        const BoxTree<NDIM> covered_tree(covered_boxes);

        d_flagged_boxes[ln].resize(level->getBoxes().size());
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            // The cells in the ghost region of the patch that are in the domain and not on the level.
            BoxList<NDIM> uncovered;
            if (!fill_sets_uncovered)
            {
                uncovered.addItem(Box<NDIM>::grow(level->getBoxes()[p()], d_max_ghost_width));
                uncovered.intersectBoxes(domain);
                covered_tree.removeIntersections(uncovered);
            }
            // Ghost cells outside the domain beside a boundary that the patch does not touch. Nothing fills these
            // cells, because SAMRAI gives boundary boxes only to patches that touch a boundary.
            BoxList<NDIM> unfilled;
            {
                const Box<NDIM>& patch_box = level->getBoxes()[p()];
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    if (grid_geom->getPeriodicShift(ratio)(d) != 0) continue;
                    for (const int side : { -1, +1 })
                    {
                        Box<NDIM> slab = patch_box;
                        slab.lower(d) = slab.upper(d) = (side < 0 ? patch_box.lower(d) : patch_box.upper(d)) + side;
                        BoxList<NDIM> slab_outside_domain(slab);
                        slab_outside_domain.removeIntersections(domain);
                        if (!slab_outside_domain.isEmpty()) continue;
                        // The ghost cells beyond the domain along axis d on this side, extended along the other axes
                        // where they reach the end of the patch.
                        Box<NDIM> exterior_ghost_box = patch_box;
                        if (side < 0)
                        {
                            exterior_ghost_box.upper(d) = patch_box.lower(d) - 1;
                            exterior_ghost_box.lower(d) = patch_box.lower(d) - d_max_ghost_width(d);
                        }
                        else
                        {
                            exterior_ghost_box.lower(d) = patch_box.upper(d) + 1;
                            exterior_ghost_box.upper(d) = patch_box.upper(d) + d_max_ghost_width(d);
                        }
                        if (exterior_ghost_box.empty()) continue;
                        BoxList<NDIM> outside(exterior_ghost_box);
                        outside.removeIntersections(domain);
                        for (BoxList<NDIM>::Iterator o(outside); o; o++)
                        {
                            Box<NDIM> ghost = o();
                            for (unsigned int e = 0; e < NDIM; ++e)
                            {
                                if (e == d) continue;
                                if (ghost.lower(e) == patch_box.lower(e)) ghost.lower(e) -= d_max_ghost_width(e);
                                if (ghost.upper(e) == patch_box.upper(e)) ghost.upper(e) += d_max_ghost_width(e);
                            }
                            unfilled.addItem(ghost);
                        }
                    }
                }
            }
            BoxList<NDIM> unfilled_copy(unfilled);
            uncovered.catenateItems(unfilled_copy);
            if (uncovered.isEmpty()) continue;

            // Data at an index are covered if a cell that touches them is covered. Along the axes of staggering, those
            // cells have that index or the next lower index.
            std::vector<BoxList<NDIM>>& flagged_boxes = d_flagged_boxes[ln][p()];
            for (Box<NDIM>::Iterator s(STAGGERINGS); s; s++)
            {
                flagged_boxes.push_back(uncovered);
                for (const IntVector<NDIM>& offset : touching_offsets(IntVector<NDIM>(s())))
                {
                    BoxList<NDIM> shifted_uncovered(uncovered);
                    shifted_uncovered.shift(offset);
                    flagged_boxes.back().intersectBoxes(shifted_uncovered);
                }
            }
            // The flagged data that touch an unfilled cell, for the message.
            for (Box<NDIM>::Iterator s(STAGGERINGS); s; s++)
            {
                BoxList<NDIM> touching_unfilled(unfilled);
                for (const IntVector<NDIM>& offset : touching_offsets(IntVector<NDIM>(s())))
                {
                    BoxList<NDIM> shifted_unfilled(unfilled);
                    shifted_unfilled.shift(offset);
                    touching_unfilled.catenateItems(shifted_unfilled);
                }
                BoxList<NDIM> flagged_unfilled(flagged_boxes[STAGGERINGS.offset(s())]);
                flagged_unfilled.intersectBoxes(touching_unfilled);
                flagged_boxes.push_back(flagged_unfilled);
            }
        }
    }
    return;
} // reset

template <class T>
void
InterpolationStencilCheck::check(const std::string& caller,
                                 const Pointer<Patch<NDIM>> patch,
                                 const PatchData<NDIM>& q_data,
                                 const Box<NDIM>& interp_box,
                                 const double* const X_data,
                                 const Pointer<LIndexSetData<T>> idx_data,
                                 const std::string& kernel_fcn) const
{
    if (!hasFlaggedCells(*patch)) return;
    std::vector<int> local_indices;
    std::vector<double> periodic_shifts;
    LEInteractor::buildLocalIndices(local_indices, periodic_shifts, interp_box, patch, X_data, NDIM, idx_data);
    checkPoints(caller, *patch, q_data, kernel_fcn, X_data, local_indices, periodic_shifts);
    return;
} // check

void
InterpolationStencilCheck::check(const std::string& caller,
                                 const Pointer<Patch<NDIM>> patch,
                                 const PatchData<NDIM>& q_data,
                                 const Box<NDIM>& interp_box,
                                 const double* const X_data,
                                 const int X_size,
                                 const std::string& kernel_fcn) const
{
    if (!hasFlaggedCells(*patch)) return;
    const Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
    std::vector<int> local_indices;
    for (int k = 0; k < X_size / NDIM; ++k)
    {
        const hier::Index<NDIM> i = IndexUtilities::getCellIndex(&X_data[NDIM * k],
                                                                 pgeom->getXLower(),
                                                                 pgeom->getXUpper(),
                                                                 pgeom->getDx(),
                                                                 patch->getBox().lower(),
                                                                 patch->getBox().upper());
        if (interp_box.contains(i)) local_indices.push_back(k);
    }
    checkPoints(caller,
                *patch,
                q_data,
                kernel_fcn,
                X_data,
                local_indices,
                std::vector<double>(NDIM * local_indices.size(), 0.0));
    return;
} // check

/////////////////////////////// PRIVATE //////////////////////////////////////

bool
InterpolationStencilCheck::hasFlaggedCells(const Patch<NDIM>& patch) const
{
    const int ln = patch.getPatchLevelNumber();
    const int pn = patch.getPatchNumber();
    if (ln >= static_cast<int>(d_levels.size()) || pn >= d_levels[ln]->getNumberOfPatches() ||
        d_levels[ln]->getPatch(pn).getPointer() != &patch)
    {
        TBOX_ERROR("InterpolationStencilCheck::check():\n"
                   << "  patch " << pn << " of level " << ln << " is not a patch that the last call to reset() found\n"
                   << "  reset the check whenever the patch hierarchy changes\n");
    }
    return !d_flagged_boxes[ln][pn].empty();
}

void
InterpolationStencilCheck::checkPoints(const std::string& caller,
                                       const Patch<NDIM>& patch,
                                       const PatchData<NDIM>& q_data,
                                       const std::string& kernel_fcn,
                                       const double* const X_data,
                                       const std::vector<int>& local_indices,
                                       const std::vector<double>& periodic_shifts) const
{
    const int ln = patch.getPatchLevelNumber();
    const std::vector<BoxList<NDIM>>& flagged_boxes = d_flagged_boxes[ln][patch.getPatchNumber()];
    const int n = num_nodes_used(kernel_fcn);
    for (const IntVector<NDIM>& staggering : staggerings(q_data))
    {
        // The stencil is trimmed to the ghost region of the data.
        Box<NDIM> data_box = q_data.getGhostBox();
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            data_box.upper(d) += staggering(d);
        }
        for (std::size_t k = 0; k < local_indices.size(); ++k)
        {
            std::array<double, NDIM> X;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X[d] = X_data[NDIM * local_indices[k] + d] + periodic_shifts[NDIM * k + d];
            }
            check_stencil(caller,
                          ln,
                          patch.getBox(),
                          X,
                          stencil_box(X, patch, staggering, n) * data_box,
                          flagged_boxes[STAGGERINGS.offset(hier::Index<NDIM>(staggering))],
                          flagged_boxes[STAGGERINGS.size() + STAGGERINGS.offset(hier::Index<NDIM>(staggering))]);
        }
    }
    return;
} // checkPoints

/////////////////////////////// TEMPLATES ////////////////////////////////////

template void InterpolationStencilCheck::check(const std::string& caller,
                                               Pointer<Patch<NDIM>> patch,
                                               const PatchData<NDIM>& q_data,
                                               const Box<NDIM>& interp_box,
                                               const double* X_data,
                                               Pointer<LIndexSetData<LNode>> idx_data,
                                               const std::string& kernel_fcn) const;

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK
