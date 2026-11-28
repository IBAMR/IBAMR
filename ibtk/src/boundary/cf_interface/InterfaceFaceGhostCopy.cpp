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

#include <ibtk/ibtk_utilities.h>
#include <ibtk/private/ArrayDataTransaction.h>
#include <ibtk/private/InterfaceFaceGhostCopy.h>
#include <ibtk/private/LevelBoxIndex.h>
#include <ibtk/private/periodic_image_offsets.h>

#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <BoxArray.h>
#include <BoxList.h>
#include <CartesianGridGeometry.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <ProcessorMapping.h>
#include <SideData.h>
#include <SideDataFactory.h>
#include <Variable.h>
#include <VariableDatabase.h>
#include <VariableFillPattern.h>

#include <set>
#include <utility>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
namespace detail
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// The faces of each axis on the boundary of the patch with box b along the coarse-fine interface: the faces between
// the patch and the cells next to it that are not in the level and are in the domain or have a periodic image in it.
std::vector<BoxList<NDIM>>
interface_faces(const Box<NDIM>& b,
                LevelBoxIndex& index,
                const std::vector<IntVector<NDIM>>& offsets,
                const BoxList<NDIM>& image_domain)
{
    std::vector<BoxList<NDIM>> faces(NDIM);
    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (const int side : { -1, +1 })
        {
            Box<NDIM> adjacent = b;
            adjacent.lower(axis) = adjacent.upper(axis) = (side < 0 ? b.lower(axis) - 1 : b.upper(axis) + 1);
            BoxList<NDIM> list(adjacent);
            for (const auto& offset : offsets)
            {
                list.shift(-offset);
                index.removeIntersections(list);
                list.shift(offset);
            }
            list.intersectBoxes(image_domain);
            if (side < 0)
            {
                IntVector<NDIM> up(0);
                up(axis) = 1;
                list.shift(up);
            }
            faces[axis].unionBoxes(list);
        }
    }
    return faces;
} // interface_faces
} // namespace

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
      d_protect(protect_faces),
      d_levels(finest_ln + 1)
{
    // Ghost faces exist only if the data have ghost cells.
    const IntVector<NDIM> ghost_width = VariableDatabase<NDIM>::getDatabase()
                                            ->getPatchDescriptor()
                                            ->getPatchDataFactory(d_data_idx)
                                            ->getGhostCellWidth();
    const bool has_ghosts = ghost_width.max() > 0;
    d_shape = fill_pattern_shape(fill_pattern, ghost_width);
    for (int ln = std::max(1, coarsest_ln); ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        const BoxList<NDIM> domain_cells(level->getPhysicalDomain());
        const std::vector<IntVector<NDIM>> offsets =
            periodic_image_offsets(d_hierarchy->getGridGeometry()->getPeriodicShift(level->getRatio()));

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
        LevelBoxIndex index(*level);
        Level& L = d_levels[ln];
        std::map<int, std::vector<BoxList<NDIM>>> faces;
        auto faces_of = [&](const int k) -> const std::vector<BoxList<NDIM>>&
        {
            auto it = faces.find(k);
            if (it == faces.end())
            {
                it = faces.emplace(k, interface_faces(boxes[k], index, offsets, image_domain)).first;
            }
            return it->second;
        };
        if (d_protect)
        {
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                L.faces[p()] = faces_of(p());
            }
        }
        if (!has_ghosts)
        {
            continue;
        }

        // The pairs of a destination and a source patch, one of them local, whose boxes or periodic images are closer
        // than the ghost width plus one: a face is on the boundary of its patch, and the ghost cells of the other
        // patch can contain it.
        std::set<std::pair<int, int>> pairs;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            for (const auto& offset : offsets)
            {
                Box<NDIM> query = Box<NDIM>::grow(boxes[p()], ghost_width + 1);
                query.shift(offset);
                tbox::Array<int> indices;
                index.findOverlapIndices(indices, query);
                for (int j = 0; j < indices.size(); ++j)
                {
                    pairs.emplace(p(), indices[j]);
                    pairs.emplace(indices[j], p());
                }
            }
        }

        // The faces of the source within the fill pattern of the destination, shifted to the ghost cells of the
        // destination. The faces of a patch are in its interior, which the copy leaves alone.
        std::map<int, std::vector<BoxList<NDIM>>> fill_indices;
        for (const auto& pair : pairs)
        {
            const int dst = pair.first, src = pair.second;
            auto fill_it = fill_indices.find(dst);
            if (fill_it == fill_indices.end())
            {
                fill_it =
                    fill_indices.emplace(dst, fill_pattern_indices(d_data_idx, fill_pattern, boxes[dst], ln)).first;
            }
            const std::vector<BoxList<NDIM>>& src_faces = faces_of(src);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (const auto& offset : offsets)
                {
                    if (dst == src && offset == IntVector<NDIM>(0))
                    {
                        continue;
                    }
                    BoxList<NDIM> list = src_faces[axis];
                    list.shift(offset);
                    list.intersectBoxes(fill_it->second[axis]);
                    for (typename BoxList<NDIM>::Iterator it(list); it; it++)
                    {
                        L.copies.push_back({ dst, src, axis, offset, *it });
                    }
                }
            }
        }
    }
} // InterfaceFaceGhostCopy

void
InterfaceFaceGhostCopy::resetPatchDataIndices(const int dst_data_idx,
                                              const int src_data_idx,
                                              Pointer<VariableFillPattern<NDIM>> fill_pattern)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const IntVector<NDIM> ghost_width =
        var_db->getPatchDescriptor()->getPatchDataFactory(dst_data_idx)->getGhostCellWidth();
    if (fill_pattern_shape(fill_pattern, ghost_width) != d_shape)
    {
        Pointer<Variable<NDIM>> var;
        var_db->mapIndexToVariable(dst_data_idx, var);
        TBOX_ERROR("InterfaceFaceGhostCopy::resetPatchDataIndices():\n"
                   << "  the fill width or fill pattern of variable " << var->getName() << " changes\n"
                   << "  call initializeOperatorState() again to change the fill width or fill pattern\n");
    }
    d_data_idx = dst_data_idx;
    d_src_idx = src_data_idx;
} // resetPatchDataIndices

bool
InterfaceFaceGhostCopy::protectsFaces() const
{
    return d_protect;
} // protectsFaces

void
InterfaceFaceGhostCopy::save(const int ln)
{
    Level& L = d_levels[ln];
    L.saved_values.clear();
    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        auto it = L.faces.find(p());
        if (it == L.faces.end())
        {
            continue;
        }
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
} // save

void
InterfaceFaceGhostCopy::restore(const Patch<NDIM>& patch)
{
    const int ln = patch.getPatchLevelNumber();
    if (ln < 1 || ln >= static_cast<int>(d_levels.size()))
    {
        return;
    }
    Level& L = d_levels[ln];
    auto it = L.faces.find(patch.getPatchNumber());
    auto saved_it = L.saved_values.find(patch.getPatchNumber());
    if (it == L.faces.end() || saved_it == L.saved_values.end())
    {
        return;
    }
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
} // restore

void
InterfaceFaceGhostCopy::appendTransactions(const int ln, Schedule& schedule)
{
    if (ln < 1 || ln >= static_cast<int>(d_levels.size()))
    {
        return;
    }
    const Level& L = d_levels[ln];
    if (L.copies.empty())
    {
        return;
    }
    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
    const ProcessorMapping& mapping = level->getProcessorMapping();
    Pointer<SideDataFactory<NDIM, double>> factory =
        VariableDatabase<NDIM>::getDatabase()->getPatchDescriptor()->getPatchDataFactory(d_data_idx);
    const int depth = factory->getDefaultDepth();
    for (const Copy& c : L.copies)
    {
        // The arrays are looked up when the transaction runs, since the patch data indices can change.
        auto src = [this, level, c]() -> ArrayData<NDIM, double>&
        {
            Pointer<SideData<NDIM, double>> data = level->getPatch(c.src_patch)->getPatchData(d_src_idx);
            return data->getArrayData(c.axis);
        };
        auto dst = [this, level, c]() -> ArrayData<NDIM, double>&
        {
            Pointer<SideData<NDIM, double>> data = level->getPatch(c.dst_patch)->getPatchData(d_data_idx);
            return data->getArrayData(c.axis);
        };
        schedule.appendTransaction(new ArrayDataTransaction(mapping.getProcessorAssignment(c.src_patch),
                                                            mapping.getProcessorAssignment(c.dst_patch),
                                                            c.box,
                                                            c.offset,
                                                            depth,
                                                            src,
                                                            dst));
    }
} // appendTransactions

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
