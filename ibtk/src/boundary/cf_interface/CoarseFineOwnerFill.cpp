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
#include <ibtk/private/CoarseFineOwnerFill.h>
#include <ibtk/private/LevelBoxIndex.h>
#include <ibtk/private/data_kind.h>
#include <ibtk/private/periodic_image_offsets.h>

#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <BoxArray.h>
#include <BoxList.h>
#include <BoxTree.h>
#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <CellDataFactory.h>
#include <CellVariable.h>
#include <Patch.h>
#include <PatchDataFactory.h>
#include <PatchDescriptor.h>
#include <PatchLevel.h>
#include <ProcessorMapping.h>
#include <SideData.h>
#include <SideDataFactory.h>
#include <SideVariable.h>
#include <Variable.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <limits>
#include <set>
#include <string>
#include <utility>

#include <ibtk/namespaces.h> // IWYU pragma: keep

// FORTRAN ROUTINES
#if (NDIM == 2)
#define CC_QUAD_TANGENTIAL_INTERPOLATION_FC                                                                            \
    IBTK_FC_FUNC(ccquadtangentialinterpolation2d, CCQUADTANGENTIALINTERPOLATION2D)
#define CC_QUAD_NORMAL_INTERPOLATION_FC IBTK_FC_FUNC(ccquadnormalinterpolation2d, CCQUADNORMALINTERPOLATION2D)
#define SC_QUAD_TANGENTIAL_INTERPOLATION_FC                                                                            \
    IBTK_FC_FUNC(scquadtangentialinterpolation2d, SCQUADTANGENTIALINTERPOLATION2D)
#define SC_QUAD_NORMAL_INTERPOLATION_FC IBTK_FC_FUNC(scquadnormalinterpolation2d, SCQUADNORMALINTERPOLATION2D)
#endif
#if (NDIM == 3)
#define CC_QUAD_TANGENTIAL_INTERPOLATION_FC                                                                            \
    IBTK_FC_FUNC(ccquadtangentialinterpolation3d, CCQUADTANGENTIALINTERPOLATION3D)
#define CC_QUAD_NORMAL_INTERPOLATION_FC IBTK_FC_FUNC(ccquadnormalinterpolation3d, CCQUADNORMALINTERPOLATION3D)
#define SC_QUAD_TANGENTIAL_INTERPOLATION_FC                                                                            \
    IBTK_FC_FUNC(scquadtangentialinterpolation3d, SCQUADTANGENTIALINTERPOLATION3D)
#define SC_QUAD_NORMAL_INTERPOLATION_FC IBTK_FC_FUNC(scquadnormalinterpolation3d, SCQUADNORMALINTERPOLATION3D)
#endif

// Function interfaces
extern "C"
{
    void CC_QUAD_TANGENTIAL_INTERPOLATION_FC(double* U_fine,
                                             const int& U_fine_gcw,
                                             const double* U_coarse,
                                             const int& U_crse_gcw,
                                             const int& ilowerf0,
                                             const int& iupperf0,
                                             const int& ilowerf1,
                                             const int& iupperf1,
#if (NDIM == 3)
                                             const int& ilowerf2,
                                             const int& iupperf2,
#endif
                                             const int& ilowerc0,
                                             const int& iupperc0,
                                             const int& ilowerc1,
                                             const int& iupperc1,
#if (NDIM == 3)
                                             const int& ilowerc2,
                                             const int& iupperc2,
#endif
                                             const int& loc_index,
                                             const int* ratio_to_coarser,
                                             const int* blower,
                                             const int* bupper);

    void CC_QUAD_NORMAL_INTERPOLATION_FC(double* U,
                                         const int& U_gcw,
                                         const int& ilower0,
                                         const int& iupper0,
                                         const int& ilower1,
                                         const int& iupper1,
#if (NDIM == 3)
                                         const int& ilower2,
                                         const int& iupper2,
#endif
                                         const int& loc_index,
                                         const int* ratio_to_coarser,
                                         const int* blower,
                                         const int* bupper);

    void SC_QUAD_TANGENTIAL_INTERPOLATION_FC(double* U_fine0,
                                             double* U_fine1,
#if (NDIM == 3)
                                             double* U_fine2,
#endif
                                             const int& U_fine_gcw,
                                             const double* U_crse0,
                                             const double* U_crse1,
#if (NDIM == 3)
                                             const double* U_crse2,
#endif
                                             const int& U_crse_gcw,
                                             const int* sc_indicator0,
                                             const int* sc_indicator1,
#if (NDIM == 3)
                                             const int* sc_indicator2,
#endif
                                             const int& sc_indicator_gcw,
                                             const int& ilowerf0,
                                             const int& iupperf0,
                                             const int& ilowerf1,
                                             const int& iupperf1,
#if (NDIM == 3)
                                             const int& ilowerf2,
                                             const int& iupperf2,
#endif
                                             const int& ilowerc0,
                                             const int& iupperc0,
                                             const int& ilowerc1,
                                             const int& iupperc1,
#if (NDIM == 3)
                                             const int& ilowerc2,
                                             const int& iupperc2,
#endif
                                             const int& loc_index,
                                             const int* ratio_to_coarser,
                                             const int* blower,
                                             const int* bupper);

    void SC_QUAD_NORMAL_INTERPOLATION_FC(double* U0,
                                         double* U1,
#if (NDIM == 3)
                                         double* U2,
#endif
                                         const int& U_gcw,
                                         const double* W0,
                                         const double* W1,
#if (NDIM == 3)
                                         const double* W2,
#endif
                                         const int& W_gcw,
                                         const int* sc_indicator0,
                                         const int* sc_indicator1,
#if (NDIM == 3)
                                         const int* sc_indicator2,
#endif
                                         const int& sc_indicator_gcw,
                                         const int& ilower0,
                                         const int& iupper0,
                                         const int& ilower1,
                                         const int& iupper1,
#if (NDIM == 3)
                                         const int& ilower2,
                                         const int& iupper2,
#endif
                                         const int& loc_index,
                                         const int* ratio_to_coarser,
                                         const int* blower,
                                         const int* bupper);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
namespace detail
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// The kernels support at most this many layers.
constexpr int MAX_DEPTH = 8;

// The cells (comp = -1) or the faces normal to axis comp that belong to a box of cells.
Box<NDIM>
index_box(const Box<NDIM>& cell_box, const int comp)
{
    Box<NDIM> box = cell_box;
    if (comp >= 0)
    {
        box.upper(comp) += 1;
    }
    return box;
} // index_box

int
location_index(const int d, const int s)
{
    return 2 * d + (s > 0 ? 0 : 1);
} // location_index

// Whether the first nonzero component of offset is positive.
bool
lexicographically_positive(const IntVector<NDIM>& offset)
{
    for (int d = 0; d < NDIM; ++d)
    {
        if (offset(d) != 0)
        {
            return offset(d) > 0;
        }
    }
    return false;
} // lexicographically_positive

// The fill box of the kernels for a coarse-fine value along (d, s) on the cells (faces) phi of component comp, which a
// patch with the box patch_box computes. The faces that the kernels fill are all faces of the patch.
Box<NDIM>
kernel_box(const Box<NDIM>& phi, const int comp, const int d, const int s, const Box<NDIM>& patch_box)
{
    Box<NDIM> box = phi;
    if (comp < 0)
    {
        return box;
    }
    if (comp == d)
    {
        // The normal component sits on the face of the ghost cell that is away from the fine level.
        if (s < 0)
        {
            box.shift(d, -1);
        }
    }
    else if (phi.upper(comp) > phi.lower(comp))
    {
        // The tangential components fill one more face than the fill box has cells.
        box.upper(comp) = phi.upper(comp) - 1;
    }
    else
    {
        box.lower(comp) = box.upper(comp) =
            (phi.lower(comp) <= patch_box.upper(comp) ? phi.lower(comp) : phi.lower(comp) - 1);
    }
    return box;
} // kernel_box

// The cells or faces of component comp that the kernels fill on the fill box kernel.
Box<NDIM>
filled_index_box(const Box<NDIM>& kernel, const int comp, const int d, const int s)
{
    Box<NDIM> box = kernel;
    if (comp < 0)
    {
        return box;
    }
    if (comp == d)
    {
        if (s < 0)
        {
            box.shift(d, +1);
        }
    }
    else
    {
        box.upper(comp) += 1;
    }
    return box;
} // filled_index_box

// Copy on the part of box that lies in both arrays.
void
copy_on(ArrayData<NDIM, double>& dst, const ArrayData<NDIM, double>& src, const Box<NDIM>& box)
{
    const Box<NDIM> overlap = box * dst.getBox() * src.getBox();
    if (!overlap.empty())
    {
        dst.copy(src, overlap);
    }
} // copy_on

double*
data_pointer(Pointer<PatchData<NDIM>> data, const int comp, const int depth)
{
    if (comp < 0)
    {
        Pointer<CellData<NDIM, double>> cc = data;
        return cc->getPointer(depth);
    }
    Pointer<SideData<NDIM, double>> sc = data;
    return sc->getPointer(comp, depth);
} // data_pointer

ArrayData<NDIM, double>&
array_data(Pointer<PatchData<NDIM>> data, const int comp)
{
    if (comp < 0)
    {
        Pointer<CellData<NDIM, double>> cc = data;
        return cc->getArrayData();
    }
    Pointer<SideData<NDIM, double>> sc = data;
    return sc->getArrayData(comp);
} // array_data

int
ghost_width(Pointer<PatchData<NDIM>> data)
{
    return data->getGhostCellWidth().max();
} // ghost_width

// A slab of cells (faces) of component ci along the direction with location index code, which a patch computes: the
// cells (faces) of the slab and those of them that get values, in the index space of the patch.
struct TaskGeometry
{
    std::size_t ci = 0;
    int code = 0;
    Box<NDIM> index_box;
    BoxList<NDIM> value_indices;
};

// The level data that compute_tasks() uses.
struct LevelGeometry
{
    int ln = 0;
    const BoxArray<NDIM>* boxes = nullptr;
    detail::LevelBoxIndex* index = nullptr;
    IntVector<NDIM> ratio = IntVector<NDIM>(0);
    std::vector<IntVector<NDIM>> offsets;
    // For each component, the cells (faces) that are in the domain or have a periodic image in it.
    std::vector<BoxList<NDIM>> image_domain;
};

// The tasks of the patch P of the level, in an order that is the same on every process.
std::vector<TaskGeometry>
compute_tasks(const LevelGeometry& g, const std::vector<int>& comps, const IntVector<NDIM>& fill_width, const int P)
{
    const BoxArray<NDIM>& boxes = *g.boxes;

    // Remove the cells (faces) of the fine level, with periodic images, from a list.
    auto remove_fine = [&](BoxList<NDIM>& list, const int comp)
    {
        for (const auto& offset : g.offsets)
        {
            for (int pass = 0; pass < (comp < 0 ? 1 : 2); ++pass)
            {
                // A face is in F if the cell above it or the cell below it is.
                IntVector<NDIM> shift = offset;
                if (pass == 1)
                {
                    shift(comp) += 1;
                }
                list.shift(-shift);
                g.index->removeIntersections(list);
                list.shift(shift);
            }
        }
    };

    std::vector<TaskGeometry> tasks;
    for (std::size_t ci = 0; ci < comps.size(); ++ci)
    {
        const int comp = comps[ci];
        const Box<NDIM> patch_box_index = index_box(boxes[P], comp);
        for (int d = 0; d < NDIM; ++d)
        {
            for (const int s : { +1, -1 })
            {
                // The cells (faces) of the first layer of the patch on the side of c + s e_d, and the cells (faces) c.
                Box<NDIM> layer = patch_box_index;
                if (comp == d)
                {
                    layer.lower(d) = layer.upper(d) = (s > 0 ? patch_box_index.lower(d) : patch_box_index.upper(d));
                }
                else
                {
                    layer.lower(d) = layer.upper(d) = (s > 0 ? boxes[P].lower(d) : boxes[P].upper(d));
                }
                if (comp != d && boxes[P].upper(d) - boxes[P].lower(d) < 1)
                {
                    TBOX_ERROR("CoarseFineOwnerFill::computeGeometry():\n"
                               << "  patch " << boxes[P] << " of level " << g.ln << " is one cell thick along axis "
                               << d << "\n  the coarse-fine interpolation requires patches of at least two cells\n");
                }
                // The cells (faces) c at depths 1 to the maximum depth in front of the patch, each with no cell (face)
                // of F between it and the patch.
                const int code = location_index(d, s);
                const int max_depth = std::min({ fill_width(d), g.ratio(d), MAX_DEPTH });
                IntVector<NDIM> toward(0);
                toward(d) = -s;
                std::vector<BoxList<NDIM>> depth_boxes(1);
                {
                    BoxList<NDIM> chain(layer);
                    for (int k = 1; k <= max_depth; ++k)
                    {
                        chain.shift(toward);
                        remove_fine(chain, comp);
                        chain.intersectBoxes(g.image_domain[ci]);
                        if (chain.isEmpty())
                        {
                            break;
                        }
                        depth_boxes.push_back(chain);
                    }
                }
                const int num_depths = static_cast<int>(depth_boxes.size()) - 1;
                if (num_depths == 0)
                {
                    continue;
                }
                for (int k = 1; k <= num_depths; ++k)
                {
                    BoxList<NDIM>& phi = depth_boxes[k];
                    if (comp < 0)
                    {
                        continue;
                    }
                    // A face on the boundary of two patches, or of a patch and a periodic image of a patch, has one
                    // owner: the patch of smaller index.
                    std::set<int> smaller;
                    for (typename BoxList<NDIM>::Iterator it(phi); it; it++)
                    {
                        Box<NDIM> b = *it;
                        b.shift(d, k * s);
                        tbox::Array<int> indices;
                        g.index->findOverlapIndices(indices, b);
                        b.shift(comp, -1);
                        tbox::Array<int> indices2;
                        g.index->findOverlapIndices(indices2, b);
                        for (int j = 0; j < indices.size(); ++j)
                        {
                            if (indices[j] < P)
                            {
                                smaller.insert(indices[j]);
                            }
                        }
                        for (int j = 0; j < indices2.size(); ++j)
                        {
                            if (indices2[j] < P)
                            {
                                smaller.insert(indices2[j]);
                            }
                        }
                    }
                    // The index finds a patch through any of its periodic images, so every image of its faces is
                    // removed.
                    for (const int other : smaller)
                    {
                        for (const auto& offset : g.offsets)
                        {
                            Box<NDIM> taken = index_box(boxes[other], comp);
                            taken.shift(offset);
                            taken.shift(d, -k * s);
                            phi.removeIntersections(taken);
                        }
                    }
                    // A patch that spans a period meets its own image. Of two faces that are periodic images of each
                    // other, the patch computes the one at the lower end.
                    for (const auto& offset : g.offsets)
                    {
                        if (lexicographically_positive(offset))
                        {
                            Box<NDIM> taken = patch_box_index;
                            taken.shift(offset);
                            taken.shift(d, -k * s);
                            phi.removeIntersections(taken);
                        }
                    }
                }
                // The slabs of cells (faces) that the kernels fill: the columns of exactly m depths, for each m.
                for (int m = 1; m <= num_depths; ++m)
                {
                    IntVector<NDIM> to_first(0), to_next(0);
                    to_first(d) = (m - 1) * s;
                    to_next(d) = m * s;
                    BoxList<NDIM> columns = depth_boxes[m];
                    columns.shift(to_first);
                    if (m < num_depths)
                    {
                        BoxList<NDIM> deeper = depth_boxes[m + 1];
                        deeper.shift(to_next);
                        columns.removeIntersections(deeper);
                    }
                    for (typename BoxList<NDIM>::Iterator it(columns); it; it++)
                    {
                        Box<NDIM> slab = *it;
                        if (s > 0)
                        {
                            slab.lower(d) -= (m - 1);
                        }
                        else
                        {
                            slab.upper(d) += (m - 1);
                        }
                        TaskGeometry task{ ci, code, slab, BoxList<NDIM>() };
                        for (int k = 1; k <= num_depths; ++k)
                        {
                            BoxList<NDIM> at_depth = depth_boxes[k];
                            at_depth.intersectBoxes(slab);
                            for (typename BoxList<NDIM>::Iterator pit(at_depth); pit; pit++)
                            {
                                task.value_indices.addItem(*pit);
                            }
                        }
                        tasks.push_back(task);
                    }
                }
            }
        }
    }
    return tasks;
} // compute_tasks
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

CoarseFineOwnerFill::CoarseFineOwnerFill(Pointer<PatchHierarchy<NDIM>> hierarchy,
                                         const int coarsest_ln,
                                         const int finest_ln,
                                         const int dst_data_idx,
                                         const int src_data_idx,
                                         Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern)
    : d_hierarchy(hierarchy),
      d_coarsest_ln(coarsest_ln),
      d_finest_ln(finest_ln),
      d_dst_idx(dst_data_idx),
      d_src_idx(src_data_idx),
      d_fill_pattern(fill_pattern)
{
    d_is_cell = get_data_kind(d_src_idx).centering == DataCentering::CELL;
    d_depth = get_data_kind(d_dst_idx).depth;
    if (d_is_cell)
    {
        d_comps.push_back(-1);
    }
    else
    {
        for (int axis = 0; axis < NDIM; ++axis)
        {
            d_comps.push_back(axis);
        }
    }

    setGhostWidth();

    d_levels = std::vector<Level>(d_finest_ln + 1);
    for (int ln = std::max(1, d_coarsest_ln); ln <= d_finest_ln; ++ln)
    {
        setUpLevel(ln);
    }
} // CoarseFineOwnerFill

void
CoarseFineOwnerFill::resetPatchDataIndices(const int dst_data_idx,
                                           const int src_data_idx,
                                           Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern)
{
    const DataKind kind = { get_data_kind(src_data_idx).centering, get_data_kind(dst_data_idx).depth };
    const DataKind expected = { d_is_cell ? DataCentering::CELL : DataCentering::SIDE, d_depth };
    if (kind != expected)
    {
        Pointer<Variable<NDIM>> var;
        VariableDatabase<NDIM>::getDatabase()->mapIndexToVariable(dst_data_idx, var);
        TBOX_ERROR("CoarseFineOwnerFill::resetPatchDataIndices():\n"
                   << "  variable " << var->getName() << " has " << data_kind_to_string(kind)
                   << ", but the fill was set up for " << data_kind_to_string(expected) << "\n"
                   << "  call initializeOperatorState() again to change the data\n");
    }
    d_dst_idx = dst_data_idx;
    d_src_idx = src_data_idx;
    const detail::FillPatternShape old_shape = detail::fill_pattern_shape(d_fill_pattern, d_ghost_width);
    d_fill_pattern = fill_pattern;
    const IntVector<NDIM> old_fill_width = d_fill_width;
    setGhostWidth();
    if (d_fill_width != old_fill_width)
    {
        Pointer<Variable<NDIM>> var;
        VariableDatabase<NDIM>::getDatabase()->mapIndexToVariable(d_dst_idx, var);
        TBOX_ERROR("CoarseFineOwnerFill::resetPatchDataIndices():\n"
                   << "  the fill width of variable " << var->getName() << " changes from " << old_fill_width << " to "
                   << d_fill_width << "\n"
                   << "  call initializeOperatorState() again to change the fill width\n");
    }
    if (detail::fill_pattern_shape(d_fill_pattern, d_ghost_width) != old_shape)
    {
        Pointer<Variable<NDIM>> var;
        VariableDatabase<NDIM>::getDatabase()->mapIndexToVariable(d_dst_idx, var);
        TBOX_ERROR("CoarseFineOwnerFill::resetPatchDataIndices():\n"
                   << "  the fill pattern of variable " << var->getName() << " changes shape\n"
                   << "  call initializeOperatorState() again to change the fill pattern\n");
    }
} // resetPatchDataIndices

void
CoarseFineOwnerFill::resetPatchDataIndex(const int data_idx)
{
    if (d_dst_idx == data_idx && d_src_idx == data_idx)
    {
        return;
    }
    resetPatchDataIndices(data_idx, data_idx, d_fill_pattern);
} // resetPatchDataIndex

void
CoarseFineOwnerFill::zeroTangentialValues(const int ln)
{
    for (auto& entry : d_levels[ln].work)
    {
        for (const auto& data : entry.second.tangential)
        {
            if (!data)
            {
                continue;
            }
            for (const int comp : d_comps)
            {
                array_data(data, comp).fillAll(0.0);
            }
        }
    }
} // zeroTangentialValues

void
CoarseFineOwnerFill::copyTangentialValues(const CoarseFineOwnerFill& from, const int ln)
{
    for (auto& entry : d_levels[ln].work)
    {
        const PatchWork& source = from.d_levels[ln].work.at(entry.first);
        for (std::size_t code = 0; code < entry.second.tangential.size(); ++code)
        {
            if (!entry.second.tangential[code])
            {
                continue;
            }
            for (const int comp : d_comps)
            {
                ArrayData<NDIM, double>& dst = array_data(entry.second.tangential[code], comp);
                copy_on(dst, array_data(source.tangential[code], comp), dst.getBox());
            }
        }
    }
} // copyTangentialValues

bool
CoarseFineOwnerFill::hasValues(const int ln) const
{
    return ln >= 0 && ln < static_cast<int>(d_levels.size()) && d_levels[ln].active;
} // hasValues

void
CoarseFineOwnerFill::computeValues(const int ln)
{
    if (!hasValues(ln))
    {
        return;
    }
    Level& L = d_levels[ln];
    for (auto& entry : L.work)
    {
        Pointer<Patch<NDIM>> patch = L.level->getPatch(entry.first);
        computeNormal(entry.second, *patch, L.level->getRatioToCoarserLevel());
    }
} // computeValues

void
CoarseFineOwnerFill::appendTransactions(const int ln, Schedule& schedule)
{
    if (!hasValues(ln))
    {
        return;
    }
    Level& L = d_levels[ln];
    const ProcessorMapping& mapping = L.level->getProcessorMapping();
    for (const Transfer& t : L.transfers)
    {
        // The arrays do not depend on the patch data indices.
        auto work_it = L.work.find(t.src_patch);
        auto receiver_it = L.receivers.find(t.dst_patch);
        ArrayData<NDIM, double>* src = work_it == L.work.end() ? nullptr : work_it->second.tasks[t.task].owned.get();
        ArrayData<NDIM, double>* dst =
            receiver_it == L.receivers.end() ? nullptr : receiver_it->second.slots[t.ci][t.code].get();
        schedule.appendTransaction(new ArrayDataTransaction(
            mapping.getProcessorAssignment(t.src_patch),
            mapping.getProcessorAssignment(t.dst_patch),
            t.box,
            t.offset,
            d_depth,
            [src]() -> ArrayData<NDIM, double>& { return *src; },
            [dst]() -> ArrayData<NDIM, double>& { return *dst; }));
    }
} // appendTransactions

void
CoarseFineOwnerFill::copyValues(const int ln)
{
    if (!hasValues(ln))
    {
        return;
    }
    Level& L = d_levels[ln];
    // The value is the average of the values of the directions, summed in the order of the directions.
    for (auto& entry : L.receivers)
    {
        const Receiver& r = entry.second;
        Pointer<PatchData<NDIM>> dst_data = L.level->getPatch(entry.first)->getPatchData(d_dst_idx);
        for (std::size_t ci = 0; ci < d_comps.size(); ++ci)
        {
            ArrayData<NDIM, double>& result = array_data(dst_data, d_comps[ci]);
            for (const auto& piece : r.pieces[ci])
            {
                int count = 0;
                std::array<const ArrayData<NDIM, double>*, 2 * NDIM> slots{};
                for (int code = 0; code < 2 * NDIM; ++code)
                {
                    if (piece.directions[code])
                    {
                        slots[count++] = r.slots[ci][code].get();
                    }
                }
                const double n = static_cast<double>(count);
                for (BoxIterator<NDIM> bi(piece.box); bi; bi++)
                {
                    const hier::Index<NDIM>& i = bi();
                    for (int depth = 0; depth < d_depth; ++depth)
                    {
                        double sum = 0.0;
                        for (int k = 0; k < count; ++k)
                        {
                            sum += (*slots[k])(i, depth);
                        }
                        result(i, depth) = sum / n;
                    }
                }
            }
        }
    }
} // copyValues

void
CoarseFineOwnerFill::setPhysicalBoundaryConditions(Patch<NDIM>& /*patch*/,
                                                   const double /*fill_time*/,
                                                   const IntVector<NDIM>& /*ghost_width_to_fill*/)
{
    // intentionally blank
} // setPhysicalBoundaryConditions

IntVector<NDIM>
CoarseFineOwnerFill::getRefineOpStencilWidth() const
{
    return IntVector<NDIM>(1);
} // getRefineOpStencilWidth

void
CoarseFineOwnerFill::preprocessRefine(Patch<NDIM>& /*fine*/,
                                      const Patch<NDIM>& /*coarse*/,
                                      const Box<NDIM>& /*fine_box*/,
                                      const IntVector<NDIM>& /*ratio*/)
{
    // intentionally blank
} // preprocessRefine

void
CoarseFineOwnerFill::postprocessRefine(Patch<NDIM>& fine,
                                       const Patch<NDIM>& coarse,
                                       const Box<NDIM>& fine_box,
                                       const IntVector<NDIM>& ratio)
{
    // The schedule also uses this object for temporary coarser levels, whose patches are not in the hierarchy.
    if (!fine.inHierarchy())
    {
        return;
    }
    const int ln = fine.getPatchLevelNumber();
    if (ln < 1 || ln >= static_cast<int>(d_levels.size()) || !d_levels[ln].active)
    {
        return;
    }
    Level& L = d_levels[ln];
    auto it = L.work.find(fine.getPatchNumber());
    if (it == L.work.end() || &fine != L.level->getPatch(fine.getPatchNumber()).getPointer())
    {
        return;
    }
    computeTangential(ln, it->second, coarse, fine_box, ratio);
} // postprocessRefine

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CoarseFineOwnerFill::setGhostWidth()
{
    Pointer<PatchDataFactory<NDIM>> dst_factory =
        VariableDatabase<NDIM>::getDatabase()->getPatchDescriptor()->getPatchDataFactory(d_dst_idx);
    d_ghost_width = dst_factory->getGhostCellWidth();
    // The fill pattern can restrict the width of the fill to less than the ghost width of the data.
    d_fill_width = detail::fill_pattern_shape(d_fill_pattern, d_ghost_width).width;
    d_work_gcw = std::max(2, std::min(d_fill_width.max(), MAX_DEPTH) + 1);
} // setGhostWidth

void
CoarseFineOwnerFill::checkFillWidth(const int ln) const
{
    const IntVector<NDIM> ratio = d_levels[ln].level->getRatioToCoarserLevel();
    int axis = -1;
    for (int d = 0; d < NDIM && axis < 0; ++d)
    {
        if (d_fill_width(d) > ratio(d))
        {
            axis = d;
        }
    }
    if (axis < 0)
    {
        return;
    }
    Pointer<Variable<NDIM>> var;
    VariableDatabase<NDIM>::getDatabase()->mapIndexToVariable(d_dst_idx, var);

    TBOX_ERROR("CoarseFineOwnerFill::checkFillWidth():\n"
               << "  the quadratic coarse-fine interpolation of variable " << var->getName() << " fills "
               << d_fill_width(axis) << " ghost cell layers along axis " << axis
               << ", but the refinement ratio of level " << ln << " is " << ratio(axis) << "\n"
               << "  the interpolation reaches at most one coarse cell from the coarse-fine interface\n"
               << "  use a ghost width or a fill pattern width of at most the refinement ratio\n");
} // checkFillWidth

void
CoarseFineOwnerFill::setUpLevel(const int ln)
{
    Level& L = d_levels[ln];
    L.level = d_hierarchy->getPatchLevel(ln);
    checkFillWidth(ln);
    const IntVector<NDIM> ratio = L.level->getRatioToCoarserLevel();
    for (int d = 0; d < NDIM; ++d)
    {
        if (ratio(d) < 2)
        {
            TBOX_ERROR(
                "CoarseFineOwnerFill::computeGeometry():\n"
                << "  the coarse-fine interpolation requires a refinement ratio of at least 2 but the ratio of level "
                << ln << " is " << ratio(d) << " along axis " << d << "\n");
        }
    }

    const BoxArray<NDIM>& boxes = L.level->getBoxes();
    detail::LevelBoxIndex index(*L.level);
    LevelGeometry g;
    g.ln = ln;
    g.boxes = &boxes;
    g.index = &index;
    g.ratio = ratio;
    g.offsets = periodic_image_offsets(d_hierarchy->getGridGeometry()->getPeriodicShift(L.level->getRatio()));
    const BoxList<NDIM> domain_cells(L.level->getPhysicalDomain());
    for (const int comp : d_comps)
    {
        // Each cell (face) is listed once, including the faces that two domain boxes share and those on a periodic
        // boundary, which are in the domain and in an image of it.
        BoxList<NDIM> domain_list;
        for (typename BoxList<NDIM>::Iterator it(domain_cells); it; it++)
        {
            BoxList<NDIM> faces(index_box(*it, comp));
            faces.removeIntersections(domain_list);
            domain_list.unionBoxes(faces);
        }
        BoxList<NDIM> image_domain = domain_list;
        for (const auto& offset : g.offsets)
        {
            BoxList<NDIM> image = domain_list;
            image.shift(-offset);
            image.removeIntersections(image_domain);
            image_domain.unionBoxes(image);
        }
        g.image_domain.push_back(image_domain);
    }

    // The pairs of a destination and a source patch, one of them local, whose boxes or periodic images are close
    // enough for the source to own values in the ghost cells of the destination: the values are within the fill width
    // of the source, and the fill width is at most the ghost width.
    IntVector<NDIM> reach = d_ghost_width;
    reach *= 2;
    reach += IntVector<NDIM>(1);
    std::set<std::pair<int, int>> pairs;
    for (PatchLevel<NDIM>::Iterator p(L.level); p; p++)
    {
        for (const auto& offset : g.offsets)
        {
            Box<NDIM> query = Box<NDIM>::grow(boxes[p()], reach);
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

    // The tasks of the local patches and of the patches that send values to them or get values from them.
    std::map<int, std::vector<TaskGeometry>> tasks;
    auto tasks_of = [&](const int P) -> const std::vector<TaskGeometry>&
    {
        auto it = tasks.find(P);
        if (it == tasks.end())
        {
            it = tasks.emplace(P, compute_tasks(g, d_comps, d_fill_width, P)).first;
        }
        return it->second;
    };

    // The values of each task within the fill pattern of the destination, shifted to its ghost cells.
    L.transfers.clear();
    std::map<int, std::vector<BoxList<NDIM>>> fill_indices;
    for (const auto& pair : pairs)
    {
        const int dst = pair.first, src = pair.second;
        const std::vector<TaskGeometry>& src_tasks = tasks_of(src);
        if (src_tasks.empty())
        {
            continue;
        }
        auto fill_it = fill_indices.find(dst);
        if (fill_it == fill_indices.end())
        {
            fill_it = fill_indices.emplace(dst, detail::fill_pattern_indices(d_dst_idx, d_fill_pattern, boxes[dst], ln))
                          .first;
        }
        for (std::size_t ti = 0; ti < src_tasks.size(); ++ti)
        {
            const TaskGeometry& t = src_tasks[ti];
            for (const auto& offset : g.offsets)
            {
                BoxList<NDIM> list = t.value_indices;
                list.shift(offset);
                list.intersectBoxes(fill_it->second[t.ci]);
                for (typename BoxList<NDIM>::Iterator it(list); it; it++)
                {
                    L.transfers.push_back({ dst, src, static_cast<int>(ti), t.ci, t.code, offset, *it });
                }
            }
        }
    }

    // The tasks that the local patches compute.
    L.work.clear();
    for (PatchLevel<NDIM>::Iterator p(L.level); p; p++)
    {
        const std::vector<TaskGeometry>& patch_tasks = tasks_of(p());
        if (patch_tasks.empty())
        {
            continue;
        }
        PatchWork& pw = L.work[p()];
        pw.box = boxes[p()];
        for (const TaskGeometry& t : patch_tasks)
        {
            const int comp = d_comps[t.ci];
            const int d = t.code / 2, s = (t.code % 2 == 0 ? +1 : -1);
            const Box<NDIM> kernel = kernel_box(t.index_box, comp, d, s, boxes[p()]);
            auto owned = std::make_unique<ArrayData<NDIM, double>>(index_box(kernel, comp), d_depth);
            owned->fillAll(0.0);
            pw.tasks.push_back({ comp, d, s, t.index_box, kernel, std::move(owned) });
        }
    }

    setUpReceivers(ln);
    L.active = !L.work.empty() || !L.receivers.empty() || !L.transfers.empty();
    setUpWork(ln);
} // setUpLevel

void
CoarseFineOwnerFill::setUpReceivers(const int ln)
{
    Level& L = d_levels[ln];
    L.receivers.clear();
    const ProcessorMapping& mapping = L.level->getProcessorMapping();

    // The disjoint cells (faces) of each component and direction that each local patch gets. Each value has one source.
    std::map<int, std::vector<std::array<BoxList<NDIM>, 2 * NDIM>>> regions;
    for (const Transfer& t : L.transfers)
    {
        if (!mapping.isMappingLocal(t.dst_patch))
        {
            continue;
        }
        std::vector<std::array<BoxList<NDIM>, 2 * NDIM>>& r = regions[t.dst_patch];
        if (r.empty())
        {
            r.resize(d_comps.size());
        }
        BoxList<NDIM> overlap(t.box);
        overlap.intersectBoxes(r[t.ci][t.code]);
        if (!overlap.isEmpty())
        {
            TBOX_ERROR("CoarseFineOwnerFill::setUpReceivers():\n"
                       << "  patch " << t.dst_patch << " of level " << ln << " gets the coarse-fine values in "
                       << overlap.getFirstItem() << " (component " << d_comps[t.ci] << ", direction " << t.code
                       << ") twice, once from patch " << t.src_patch << "\n"
                       << "  each coarse-fine value must have one owner\n");
        }
        r[t.ci][t.code].unionBoxes(BoxList<NDIM>(t.box));
    }

    for (const auto& entry : regions)
    {
        Receiver& receiver = L.receivers[entry.first];
        receiver.slots.resize(d_comps.size());
        receiver.pieces.resize(d_comps.size());
        for (std::size_t ci = 0; ci < d_comps.size(); ++ci)
        {
            std::vector<Piece> pieces;
            for (int code = 0; code < 2 * NDIM; ++code)
            {
                const BoxList<NDIM>& region = entry.second[ci][code];
                if (region.isEmpty())
                {
                    continue;
                }
                receiver.slots[ci][code] = std::make_unique<ArrayData<NDIM, double>>(region.getBoundingBox(), d_depth);
                receiver.slots[ci][code]->fillAll(std::numeric_limits<double>::quiet_NaN());
                std::array<bool, 2 * NDIM> only_code = {};
                only_code[code] = true;
                std::vector<Piece> next;
                BoxList<NDIM> fresh = region;
                for (const auto& piece : pieces)
                {
                    BoxList<NDIM> inside(piece.box);
                    inside.intersectBoxes(region);
                    BoxList<NDIM> outside(piece.box);
                    outside.removeIntersections(region);
                    std::array<bool, 2 * NDIM> with_code = piece.directions;
                    with_code[code] = true;
                    for (typename BoxList<NDIM>::Iterator it(inside); it; it++)
                    {
                        next.push_back({ *it, with_code });
                    }
                    for (typename BoxList<NDIM>::Iterator it(outside); it; it++)
                    {
                        next.push_back({ *it, piece.directions });
                    }
                    fresh.removeIntersections(piece.box);
                }
                for (typename BoxList<NDIM>::Iterator it(fresh); it; it++)
                {
                    next.push_back({ *it, only_code });
                }
                pieces = next;
            }
            receiver.pieces[ci] = pieces;
        }
    }
} // setUpReceivers

void
CoarseFineOwnerFill::setUpWork(const int ln)
{
    Level& L = d_levels[ln];
    for (auto& entry : L.work)
    {
        PatchWork& pw = entry.second;
        for (const auto& task : pw.tasks)
        {
            const int code = location_index(task.d, task.s);
            if (pw.tangential[code])
            {
                continue;
            }
            if (d_is_cell)
            {
                pw.tangential[code] = new CellData<NDIM, double>(pw.box, d_depth, IntVector<NDIM>(d_work_gcw));
            }
            else
            {
                pw.tangential[code] = new SideData<NDIM, double>(pw.box, d_depth, IntVector<NDIM>(d_work_gcw));
            }
        }
        if (d_is_cell)
        {
            Pointer<CellData<NDIM, double>> work =
                new CellData<NDIM, double>(pw.box, d_depth, IntVector<NDIM>(d_work_gcw));
            work->fillAll(0.0);
            pw.work = work;
        }
        else
        {
            Pointer<SideData<NDIM, double>> work =
                new SideData<NDIM, double>(pw.box, d_depth, IntVector<NDIM>(d_work_gcw));
            Pointer<SideData<NDIM, double>> result =
                new SideData<NDIM, double>(pw.box, d_depth, IntVector<NDIM>(d_work_gcw));
            work->fillAll(0.0);
            result->fillAll(0.0);
            pw.work = work;
            pw.result = result;
            pw.indicator = new SideData<NDIM, int>(pw.box, 1, IntVector<NDIM>(d_work_gcw));
            pw.indicator->fillAll(0);
        }
    }
    zeroTangentialValues(ln);
} // setUpWork

void
CoarseFineOwnerFill::computeTangential(const int ln,
                                       PatchWork& pw,
                                       const Patch<NDIM>& coarse,
                                       const Box<NDIM>& fine_box,
                                       const IntVector<NDIM>& ratio)
{
#if defined(NDEBUG)
    NULL_USE(ln);
#endif
    Pointer<PatchData<NDIM>> cdata = coarse.getPatchData(d_dst_idx);
    const Box<NDIM>& cbox = coarse.getBox();
    const Box<NDIM>& fbox = pw.box;
    const int crse_gcw = ghost_width(cdata);
    for (const auto& task : pw.tasks)
    {
        // The coarse patch covers the cells of the fine box only.
        const Box<NDIM> kernel = task.kernel_box * fine_box;
        if (kernel.empty())
        {
            continue;
        }
#if !defined(NDEBUG)
        {
            // The kernel uses the coarse cells of the kernel box and their neighbors along the tangential directions.
            // Every one of them is in the coarse level or outside the physical domain, where the physical boundary
            // fill sets it. The coarse patch of the schedule sets no other cell.
            Pointer<PatchLevel<NDIM>> coarse_level = d_hierarchy->getPatchLevel(ln - 1);
            Box<NDIM> used = Box<NDIM>::coarsen(kernel, ratio);
            for (int d = 0; d < NDIM; ++d)
            {
                if (d != task.d)
                {
                    used.grow(d, 1);
                }
            }
            const IntVector<NDIM> shift = d_hierarchy->getGridGeometry()->getPeriodicShift(coarse_level->getRatio());
            const BoxList<NDIM> coarse_domain(coarse_level->getPhysicalDomain());
            BoxList<NDIM> unset;
            for (const auto& offset : periodic_image_offsets(shift))
            {
                Box<NDIM> image = used;
                image.shift(offset);
                BoxList<NDIM> in_domain(image);
                in_domain.intersectBoxes(coarse_domain);
                unset.unionBoxes(in_domain);
            }
            coarse_level->getBoxTree()->removeIntersections(unset);
            if (!unset.isEmpty())
            {
                TBOX_ERROR("CoarseFineOwnerFill::computeTangential():\n"
                           << "  the tangential kernel of level " << ln << " uses coarse cells "
                           << *typename BoxList<NDIM>::Iterator(unset)
                           << " that are in the physical domain but not in level " << ln - 1 << "\n");
            }
        }
#endif
        const int loc = location_index(task.d, task.s);
        Pointer<PatchData<NDIM>> T = pw.tangential[loc];
        for (int depth = 0; depth < d_depth; ++depth)
        {
            if (d_is_cell)
            {
                CC_QUAD_TANGENTIAL_INTERPOLATION_FC(data_pointer(T, -1, depth),
                                                    d_work_gcw,
                                                    data_pointer(cdata, -1, depth),
                                                    crse_gcw,
                                                    fbox.lower(0),
                                                    fbox.upper(0),
                                                    fbox.lower(1),
                                                    fbox.upper(1),
#if (NDIM == 3)
                                                    fbox.lower(2),
                                                    fbox.upper(2),
#endif
                                                    cbox.lower(0),
                                                    cbox.upper(0),
                                                    cbox.lower(1),
                                                    cbox.upper(1),
#if (NDIM == 3)
                                                    cbox.lower(2),
                                                    cbox.upper(2),
#endif
                                                    loc,
                                                    ratio,
                                                    kernel.lower(),
                                                    kernel.upper());
            }
            else
            {
                Pointer<SideData<NDIM, int>> indicator = pw.indicator;
                SC_QUAD_TANGENTIAL_INTERPOLATION_FC(data_pointer(T, 0, depth),
                                                    data_pointer(T, 1, depth),
#if (NDIM == 3)
                                                    data_pointer(T, 2, depth),
#endif
                                                    d_work_gcw,
                                                    data_pointer(cdata, 0, depth),
                                                    data_pointer(cdata, 1, depth),
#if (NDIM == 3)
                                                    data_pointer(cdata, 2, depth),
#endif
                                                    crse_gcw,
                                                    indicator->getPointer(0),
                                                    indicator->getPointer(1),
#if (NDIM == 3)
                                                    indicator->getPointer(2),
#endif
                                                    d_work_gcw,
                                                    fbox.lower(0),
                                                    fbox.upper(0),
                                                    fbox.lower(1),
                                                    fbox.upper(1),
#if (NDIM == 3)
                                                    fbox.lower(2),
                                                    fbox.upper(2),
#endif
                                                    cbox.lower(0),
                                                    cbox.upper(0),
                                                    cbox.lower(1),
                                                    cbox.upper(1),
#if (NDIM == 3)
                                                    cbox.lower(2),
                                                    cbox.upper(2),
#endif
                                                    loc,
                                                    ratio,
                                                    kernel.lower(),
                                                    kernel.upper());
            }
        }
    }
} // computeTangential

void
CoarseFineOwnerFill::computeNormal(PatchWork& pw, Patch<NDIM>& patch, const IntVector<NDIM>& ratio)
{
    Pointer<PatchData<NDIM>> source = patch.getPatchData(d_src_idx);
    const Box<NDIM>& fbox = pw.box;
    for (const auto& task : pw.tasks)
    {
        const int loc = location_index(task.d, task.s);
        const Box<NDIM> filled = filled_index_box(task.kernel_box, task.comp, task.d, task.s);
        Pointer<PatchData<NDIM>> T = pw.tangential[loc];
        // The tangential values and the fine values that the formula uses, where the kernel uses them.
        ArrayData<NDIM, double>& W = array_data(pw.work, task.comp);
        copy_on(W, array_data(T, task.comp), filled);
        // The kernels use the first two layers of the patch next to the slab.
        Box<NDIM> first = task.kernel_box;
        if (task.s > 0)
        {
            first.lower(task.d) = first.upper(task.d);
        }
        else
        {
            first.upper(task.d) = first.lower(task.d);
        }
        for (const int n : { 1, 2 })
        {
            Box<NDIM> stencil = filled_index_box(first, task.comp, task.d, task.s);
            stencil.shift(task.d, n * task.s);
            copy_on(W, array_data(source, task.comp), stencil);
        }
        for (int depth = 0; depth < d_depth; ++depth)
        {
            if (d_is_cell)
            {
                CC_QUAD_NORMAL_INTERPOLATION_FC(data_pointer(pw.work, -1, depth),
                                                d_work_gcw,
                                                fbox.lower(0),
                                                fbox.upper(0),
                                                fbox.lower(1),
                                                fbox.upper(1),
#if (NDIM == 3)
                                                fbox.lower(2),
                                                fbox.upper(2),
#endif
                                                loc,
                                                ratio,
                                                task.kernel_box.lower(),
                                                task.kernel_box.upper());
            }
            else
            {
                Pointer<SideData<NDIM, int>> indicator = pw.indicator;
                SC_QUAD_NORMAL_INTERPOLATION_FC(data_pointer(pw.result, 0, depth),
                                                data_pointer(pw.result, 1, depth),
#if (NDIM == 3)
                                                data_pointer(pw.result, 2, depth),
#endif
                                                d_work_gcw,
                                                data_pointer(pw.work, 0, depth),
                                                data_pointer(pw.work, 1, depth),
#if (NDIM == 3)
                                                data_pointer(pw.work, 2, depth),
#endif
                                                d_work_gcw,
                                                indicator->getPointer(0),
                                                indicator->getPointer(1),
#if (NDIM == 3)
                                                indicator->getPointer(2),
#endif
                                                d_work_gcw,
                                                fbox.lower(0),
                                                fbox.upper(0),
                                                fbox.lower(1),
                                                fbox.upper(1),
#if (NDIM == 3)
                                                fbox.lower(2),
                                                fbox.upper(2),
#endif
                                                loc,
                                                ratio,
                                                task.kernel_box.lower(),
                                                task.kernel_box.upper());
            }
        }
        // The values go to the owned values of the task.
        copy_on(*task.owned, d_is_cell ? array_data(pw.work, -1) : array_data(pw.result, task.comp), filled);
    }
} // computeNormal

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
