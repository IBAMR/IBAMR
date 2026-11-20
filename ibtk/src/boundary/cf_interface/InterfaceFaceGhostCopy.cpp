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

#include <ibtk/InterfaceFaceGhostCopy.h>
#include <ibtk/RegionFillPattern.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <BoxArray.h>
#include <BoxGeometryFillPattern.h>
#include <BoxList.h>
#include <BoxTree.h>
#include <CartesianGridGeometry.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <ProcessorMapping.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <SideData.h>
#include <VariableDatabase.h>
#include <VariableFillPattern.h>

#include <algorithm>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

InterfaceFaceGhostCopy::InterfaceFaceGhostCopy(Pointer<PatchHierarchy<NDIM>> hierarchy,
                                               const int coarsest_ln,
                                               const int finest_ln,
                                               const int dst_data_idx,
                                               const int src_data_idx,
                                               Pointer<VariableFillPattern<NDIM>> fill_pattern,
                                               const bool protect_faces)
    : d_hierarchy(hierarchy),
      d_data_idx(dst_data_idx),
      d_src_idx(src_data_idx),
      d_fill_pattern(fill_pattern),
      d_protect(protect_faces),
      d_levels(finest_ln + 1)
{
    // Ghost faces exist only if the data have ghost cells.
    const bool has_ghosts = VariableDatabase<NDIM>::getDatabase()
                                ->getPatchDescriptor()
                                ->getPatchDataFactory(d_data_idx)
                                ->getGhostCellWidth()
                                .max() > 0;
    for (int ln = std::max(1, coarsest_ln); ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        Pointer<BoxTree<NDIM>> tree = level->getBoxTree();
        const BoxList<NDIM> domain_cells(level->getPhysicalDomain());

        // Offsets of the periodic images.
        const IntVector<NDIM> periodic_shift = d_hierarchy->getGridGeometry()->getPeriodicShift(level->getRatio());
        std::vector<IntVector<NDIM>> offsets;
        {
            int num = 1;
            for (int d = 0; d < NDIM; ++d) num *= 3;
            for (int k = 0; k < num; ++k)
            {
                IntVector<NDIM> offset(0);
                int rem = k;
                bool allowed = true;
                for (int d = 0; d < NDIM; ++d)
                {
                    const int m = rem % 3 - 1;
                    rem /= 3;
                    if (m != 0 && periodic_shift(d) == 0) allowed = false;
                    offset(d) = m * periodic_shift(d);
                }
                if (allowed) offsets.push_back(offset);
            }
        }

        // The cells that are in the domain or have a periodic image in it: a refine operator fills the ghost cells
        // across a periodic boundary from the coarse level, and that overwrites the faces next to them.
        BoxList<NDIM> image_domain;
        for (const auto& offset : offsets)
        {
            BoxList<NDIM> image = domain_cells;
            image.shift(-offset);
            image_domain.unionBoxes(image);
        }

        const BoxArray<NDIM>& boxes = level->getBoxes();
        const ProcessorMapping& mapping = level->getProcessorMapping();
        Level& L = d_levels[ln];
        L.regions.assign(NDIM, BoxList<NDIM>());
        BoxList<NDIM> carrier_cells;
        for (int k = 0; k < boxes.getNumberOfBoxes(); ++k)
        {
            const Box<NDIM>& b = boxes[k];
            std::vector<BoxList<NDIM>> faces(NDIM);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (const int side : { -1, +1 })
                {
                    // The cells next to the patch that are not in the fine level and are in the domain.
                    Box<NDIM> adjacent = b;
                    adjacent.lower(axis) = adjacent.upper(axis) = (side < 0 ? b.lower(axis) - 1 : b.upper(axis) + 1);
                    BoxList<NDIM> list(adjacent);
                    for (const auto& offset : offsets)
                    {
                        list.shift(-offset);
                        tree->removeIntersections(list);
                        list.shift(offset);
                    }
                    list.intersectBoxes(image_domain);
                    // The cells on both sides of the faces, inside the domain.
                    for (const auto& offset : offsets)
                    {
                        BoxList<NDIM> outside = list, inside = list;
                        IntVector<NDIM> toward(0);
                        toward(axis) = (side < 0 ? 1 : -1);
                        inside.shift(toward);
                        outside.shift(offset);
                        outside.intersectBoxes(domain_cells);
                        carrier_cells.unionBoxes(outside);
                        inside.shift(offset);
                        inside.intersectBoxes(domain_cells);
                        carrier_cells.unionBoxes(inside);
                    }
                    if (side < 0)
                    {
                        IntVector<NDIM> up(0);
                        up(axis) = 1;
                        list.shift(up);
                    }
                    faces[axis].unionBoxes(list);
                }
                // The faces of all patches, with periodic images, which other patches contain in their ghost cells.
                for (const auto& offset : offsets)
                {
                    BoxList<NDIM> image = faces[axis];
                    image.shift(offset);
                    L.regions[axis].unionBoxes(image);
                }
            }
            if (d_protect && mapping.isMappingLocal(k)) L.faces[k] = faces;
        }
        // The level of the cells. A cell belongs to the process of a patch next to it.
        const int num_cells = carrier_cells.getNumberOfItems();
        if (num_cells == 0 || !has_ghosts) continue;
        BoxArray<NDIM> carrier_boxes(num_cells);
        ProcessorMapping carrier_mapping(num_cells);
        int n = 0;
        for (typename BoxList<NDIM>::Iterator it(carrier_cells); it; it++, ++n)
        {
            carrier_boxes[n] = *it;
            tbox::Array<int> indices;
            tree->findOverlapIndices(indices, Box<NDIM>::grow(*it, IntVector<NDIM>(1)));
            int patch = indices.size() > 0 ? indices[0] : 0;
            for (int j = 1; j < indices.size(); ++j) patch = std::min(patch, indices[j]);
            carrier_mapping.setProcessorAssignment(n, mapping.getProcessorAssignment(patch));
        }
        L.carrier = new PatchLevel<NDIM>(carrier_boxes,
                                         carrier_mapping,
                                         level->getRatio(),
                                         d_hierarchy->getGridGeometry(),
                                         VariableDatabase<NDIM>::getDatabase()->getPatchDescriptor(),
                                         nullptr,
                                         true);
        L.carrier->allocatePatchData(d_data_idx);
        setUpSchedule(ln);
    }
}

void
InterfaceFaceGhostCopy::setUpSchedule(const int ln)
{
    Level& L = d_levels[ln];
    if (L.carrier.isNull()) return;
    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
    // The faces go from the patches to the carrier without a pattern, and from the carrier to the ghost cells within
    // the pattern of the component.
    Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> pattern_in =
        new RegionFillPattern(new SAMRAI::xfer::BoxGeometryFillPattern<NDIM>(), L.regions);
    Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> inner = d_fill_pattern;
    if (inner.isNull()) inner = new SAMRAI::xfer::BoxGeometryFillPattern<NDIM>();
    Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> pattern_out = new RegionFillPattern(inner, L.regions);
    L.alg_in = new RefineAlgorithm<NDIM>();
    L.alg_in->registerRefine(d_data_idx, d_src_idx, d_data_idx, Pointer<RefineOperator<NDIM>>(nullptr), pattern_in);
    L.sched_in = L.alg_in->createSchedule(L.carrier, level);
    L.alg_out = new RefineAlgorithm<NDIM>();
    L.alg_out->registerRefine(d_data_idx, d_data_idx, d_data_idx, Pointer<RefineOperator<NDIM>>(nullptr), pattern_out);
    L.sched_out = L.alg_out->createSchedule(level, L.carrier);
}

void
InterfaceFaceGhostCopy::resetPatchDataIndices(const int dst_data_idx,
                                              const int src_data_idx,
                                              Pointer<VariableFillPattern<NDIM>> fill_pattern)
{
    for (int ln = 1; ln < static_cast<int>(d_levels.size()); ++ln)
    {
        Level& L = d_levels[ln];
        if (L.carrier.isNull()) continue;
        L.carrier->deallocatePatchData(d_data_idx);
        L.carrier->allocatePatchData(dst_data_idx);
    }
    d_data_idx = dst_data_idx;
    d_src_idx = src_data_idx;
    d_fill_pattern = fill_pattern;
    for (int ln = 1; ln < static_cast<int>(d_levels.size()); ++ln) setUpSchedule(ln);
}

void
InterfaceFaceGhostCopy::save(const int ln)
{
    Level& L = d_levels[ln];
    L.saved_values.clear();
    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        auto it = L.faces.find(p());
        if (it == L.faces.end()) continue;
        Pointer<SideData<NDIM, double>> data = level->getPatch(p())->getPatchData(d_data_idx);
        std::vector<std::vector<double>>& saved = L.saved_values[p()];
        saved.assign(NDIM, std::vector<double>());
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (typename BoxList<NDIM>::Iterator bl(it->second[axis]); bl; bl++)
            {
                for (BoxIterator<NDIM> bi(*bl); bi; bi++)
                {
                    for (int depth = 0; depth < data->getDepth(); ++depth)
                    {
                        saved[axis].push_back(data->getArrayData(axis)(bi(), depth));
                    }
                }
            }
        }
    }
}

void
InterfaceFaceGhostCopy::restore(const Patch<NDIM>& patch)
{
    const int ln = patch.getPatchLevelNumber();
    if (ln < 1 || ln >= static_cast<int>(d_levels.size())) return;
    Level& L = d_levels[ln];
    auto it = L.faces.find(patch.getPatchNumber());
    auto saved_it = L.saved_values.find(patch.getPatchNumber());
    if (it == L.faces.end() || saved_it == L.saved_values.end()) return;
    Pointer<SideData<NDIM, double>> data = patch.getPatchData(d_data_idx);
    const std::vector<std::vector<double>>& saved = saved_it->second;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        std::size_t k = 0;
        for (typename BoxList<NDIM>::Iterator bl(it->second[axis]); bl; bl++)
        {
            for (BoxIterator<NDIM> bi(*bl); bi; bi++)
            {
                for (int depth = 0; depth < data->getDepth(); ++depth)
                {
                    data->getArrayData(axis)(bi(), depth) = saved[axis][k++];
                }
            }
        }
    }
}

void
InterfaceFaceGhostCopy::copyToGhostCells(const int ln, const double fill_time)
{
    Level& L = d_levels[ln];
    L.saved_values.clear();
    if (L.sched_in)
    {
        L.sched_in->fillData(fill_time);
        L.sched_out->fillData(fill_time);
    }
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
