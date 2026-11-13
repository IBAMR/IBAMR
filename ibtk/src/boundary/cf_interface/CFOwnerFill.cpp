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

#include <ibtk/CFOwnerFill.h>
#include <ibtk/CellNoCornersFillPattern.h>
#include <ibtk/RegionFillPattern.h>
#include <ibtk/SideNoCornersFillPattern.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <BoxArray.h>
#include <BoxGeometryFillPattern.h>
#include <BoxList.h>
#include <BoxTree.h>
#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <CellDataFactory.h>
#include <CellOverlap.h>
#include <CellVariable.h>
#include <Patch.h>
#include <PatchDataFactory.h>
#include <PatchDescriptor.h>
#include <PatchLevel.h>
#include <ProcessorMapping.h>
#include <SideData.h>
#include <SideDataFactory.h>
#include <SideOverlap.h>
#include <SideVariable.h>
#include <Variable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <limits>
#include <set>
#include <string>

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
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// The kernels may reach one cell beyond a patch.
// The kernels support at most this many layers.
constexpr int MAX_DEPTH = 8;

// The cells (comp = -1) or the faces normal to axis comp that belong to a box of cells.
Box<NDIM>
index_box(const Box<NDIM>& cell_box, const int comp)
{
    Box<NDIM> box = cell_box;
    if (comp >= 0) box.upper(comp) += 1;
    return box;
}

int
location_index(const int d, const int s)
{
    return 2 * d + (s > 0 ? 0 : 1);
}

// The fill box of the kernels for a coarse-fine value along (d, s) on the cells (faces) phi of component comp, which a
// patch with the box patch_box computes. The faces that the kernels fill are all faces of the patch.
Box<NDIM>
kernel_box(const Box<NDIM>& phi, const int comp, const int d, const int s, const Box<NDIM>& patch_box)
{
    Box<NDIM> box = phi;
    if (comp < 0) return box;
    if (comp == d)
    {
        // The normal component sits on the face of the ghost cell that is away from the fine level.
        if (s < 0) box.shift(d, -1);
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
}

// The cells or faces of component comp that the kernels fill on the fill box kernel.
Box<NDIM>
filled_index_box(const Box<NDIM>& kernel, const int comp, const int d, const int s)
{
    Box<NDIM> box = kernel;
    if (comp < 0) return box;
    if (comp == d)
    {
        if (s < 0) box.shift(d, +1);
    }
    else
    {
        box.upper(comp) += 1;
    }
    return box;
}

// Copy on the part of box that lies in both arrays.
void
copy_on(ArrayData<NDIM, double>& dst, const ArrayData<NDIM, double>& src, const Box<NDIM>& box)
{
    const Box<NDIM> overlap = box * dst.getBox() * src.getBox();
    if (!overlap.empty()) dst.copy(src, overlap);
}

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
}

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
}

int
ghost_width(Pointer<PatchData<NDIM>> data)
{
    return data->getGhostCellWidth().max();
}

// The centering of the data with index idx: 0 for cell data, 1 for side data, and -1 otherwise.
int
data_centering(const int idx)
{
    Pointer<Variable<NDIM>> var;
    VariableDatabase<NDIM>::getDatabase()->mapIndexToVariable(idx, var);
    Pointer<CellVariable<NDIM, double>> cc_var = var;
    Pointer<SideVariable<NDIM, double>> sc_var = var;
    return cc_var ? 0 : (sc_var ? 1 : -1);
}

// The depth of the cell or side data with index idx, as the patch descriptor has it.
int
data_depth(const int idx, const bool is_cell)
{
    Pointer<PatchDataFactory<NDIM>> factory =
        VariableDatabase<NDIM>::getDatabase()->getPatchDescriptor()->getPatchDataFactory(idx);
    Pointer<CellDataFactory<NDIM, double>> cc_factory = factory;
    Pointer<SideDataFactory<NDIM, double>> sc_factory = factory;
    return is_cell ? cc_factory->getDefaultDepth() : sc_factory->getDefaultDepth();
}

std::string
data_to_string(const int centering, const int depth)
{
    return std::string(centering == 0 ? "cell" : (centering == 1 ? "side" : "other")) + "-centered data of depth " +
           std::to_string(depth);
}
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

CFOwnerFill::CFOwnerFill(Pointer<PatchHierarchy<NDIM>> hierarchy,
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
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_is_cell = data_centering(d_src_idx) == 0;
    d_depth = data_depth(d_dst_idx, d_is_cell);
    if (d_is_cell)
    {
        d_comps.push_back(-1);
    }
    else
    {
        for (int axis = 0; axis < NDIM; ++axis) d_comps.push_back(axis);
    }

    setGhostWidth();

    // The temporary levels have ghost width zero, so they do not change the maximum ghost width. The patch data
    // indices are shared by all objects of the same centering and depth.
    static std::map<std::pair<bool, int>, std::array<int, 2 * NDIM + 2>> shared_indices;
    const std::string kind = std::string(d_is_cell ? "cc" : "sc") + std::to_string(d_depth);
    const std::array<std::string, 2 * NDIM + 2> names = [&]
    {
        std::array<std::string, 2 * NDIM + 2> n;
        n[0] = "IBTK::CFOwnerFill::owned_" + kind;
        for (int code = 0; code < 2 * NDIM; ++code)
            n[1 + code] = "IBTK::CFOwnerFill::slot" + std::to_string(code) + "_" + kind;
        n[2 * NDIM + 1] = "IBTK::CFOwnerFill::result_" + kind;
        return n;
    }();
    auto& indices = shared_indices[std::make_pair(d_is_cell, d_depth)];
    bool valid = indices[0] != 0 || indices[1] != 0;
    for (std::size_t k = 0; k < names.size() && valid; ++k)
    {
        Pointer<Variable<NDIM>> v;
        var_db->mapIndexToVariable(indices[k], v);
        valid = v && v->getName() == names[k];
    }
    if (!valid)
    {
        Pointer<VariableContext> ctx = var_db->getContext("IBTK::CFOwnerFill::CONTEXT");
        for (std::size_t k = 0; k < names.size(); ++k)
        {
            Pointer<Variable<NDIM>> v = var_db->checkVariableExists(names[k]) ? var_db->getVariable(names[k]) : nullptr;
            if (!v && d_is_cell)
                v = new CellVariable<NDIM, double>(names[k], d_depth);
            else if (!v)
                v = new SideVariable<NDIM, double>(names[k], d_depth);
            indices[k] = var_db->registerVariableAndContext(v, ctx, IntVector<NDIM>(0));
        }
    }
    d_owned_idx = indices[0];
    for (int code = 0; code < 2 * NDIM; ++code) d_slot_idx[code] = indices[1 + code];
    d_result_idx = indices[2 * NDIM + 1];

    d_levels.resize(d_finest_ln + 1);
    for (int ln = std::max(1, d_coarsest_ln); ln <= d_finest_ln; ++ln)
    {
        setUpLevel(ln);
    }
}

CFOwnerFill::~CFOwnerFill()
{
    clearLevels();
}

void
CFOwnerFill::setGhostWidth()
{
    Pointer<PatchDataFactory<NDIM>> dst_factory =
        VariableDatabase<NDIM>::getDatabase()->getPatchDescriptor()->getPatchDataFactory(d_dst_idx);
    d_ghost_width = dst_factory->getGhostCellWidth();
    // The fill pattern can restrict the width of the fill to less than the ghost width of the data.
    d_fill_width = d_ghost_width;
    IntVector<NDIM> pattern_width(std::numeric_limits<int>::max());
    if (Pointer<CellNoCornersFillPattern> cc_pattern = d_fill_pattern)
        pattern_width = cc_pattern->getStencilWidth();
    else if (Pointer<SideNoCornersFillPattern> sc_pattern = d_fill_pattern)
        pattern_width = sc_pattern->getStencilWidth();
    for (int d = 0; d < NDIM; ++d) d_fill_width(d) = std::min(d_ghost_width(d), pattern_width(d));
    d_work_gcw = std::max(2, std::min(d_fill_width.max(), MAX_DEPTH) + 1);
}

void
CFOwnerFill::clearLevels()
{
    for (auto& L : d_levels)
    {
        L.copy_sched.setNull();
        L.copy_alg.setNull();
        for (auto& o : L.owned)
        {
            o.sched.setNull();
            o.alg.setNull();
            o.level->deallocatePatchData(d_owned_idx);
            o.level.setNull();
        }
        if (L.value_level)
        {
            for (int code = 0; code < 2 * NDIM; ++code) L.value_level->deallocatePatchData(d_slot_idx[code]);
            L.value_level->deallocatePatchData(d_result_idx);
        }
        L.work.clear();
        L.value_level.setNull();
    }
    d_levels.clear();
}

void
CFOwnerFill::resetPatchDataIndices(const int dst_data_idx,
                                   const int src_data_idx,
                                   Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern)
{
    const int centering = data_centering(src_data_idx);
    const int depth = centering < 0 ? 0 : data_depth(dst_data_idx, centering == 0);
    if (centering != (d_is_cell ? 0 : 1) || depth != d_depth)
    {
        Pointer<Variable<NDIM>> var;
        VariableDatabase<NDIM>::getDatabase()->mapIndexToVariable(dst_data_idx, var);
        TBOX_ERROR("CFOwnerFill::resetPatchDataIndices():\n"
                   << "  variable " << var->getName() << " has " << data_to_string(centering, depth)
                   << ", but the fill was set up for " << data_to_string(d_is_cell ? 0 : 1, d_depth) << "\n"
                   << "  call initializeOperatorState() again to change the data\n");
    }
    d_dst_idx = dst_data_idx;
    d_src_idx = src_data_idx;
    d_fill_pattern = fill_pattern;
    const IntVector<NDIM> old_fill_width = d_fill_width;
    setGhostWidth();
    if (d_fill_width != old_fill_width)
    {
        Pointer<Variable<NDIM>> var;
        VariableDatabase<NDIM>::getDatabase()->mapIndexToVariable(d_dst_idx, var);
        TBOX_ERROR("CFOwnerFill::resetPatchDataIndices():\n"
                   << "  the fill width of variable " << var->getName() << " changes from " << old_fill_width << " to "
                   << d_fill_width << "\n"
                   << "  call initializeOperatorState() again to change the fill width\n");
    }
    for (int ln = std::max(1, d_coarsest_ln); ln <= d_finest_ln; ++ln)
    {
        if (d_levels[ln].active) setUpSchedules(ln);
    }
}

bool
CFOwnerFill::hasValues(const int ln) const
{
    return ln >= 0 && ln < static_cast<int>(d_levels.size()) && d_levels[ln].active;
}

void
CFOwnerFill::copyValues(const int ln, const double fill_time)
{
    d_levels[ln].copy_sched->fillData(fill_time);
}

void
CFOwnerFill::setPhysicalBoundaryConditions(Patch<NDIM>& /*patch*/,
                                           const double /*fill_time*/,
                                           const IntVector<NDIM>& /*ghost_width_to_fill*/)
{
    // intentionally blank
}

IntVector<NDIM>
CFOwnerFill::getRefineOpStencilWidth() const
{
    return IntVector<NDIM>(1);
}

void
CFOwnerFill::preprocessRefine(Patch<NDIM>& /*fine*/,
                              const Patch<NDIM>& /*coarse*/,
                              const Box<NDIM>& /*fine_box*/,
                              const IntVector<NDIM>& /*ratio*/)
{
    // intentionally blank
}

void
CFOwnerFill::postprocessRefine(Patch<NDIM>& fine,
                               const Patch<NDIM>& coarse,
                               const Box<NDIM>& fine_box,
                               const IntVector<NDIM>& ratio)
{
    // The schedule also uses this object for temporary coarser levels, whose patches are not in the hierarchy.
    if (!fine.inHierarchy()) return;
    const int ln = fine.getPatchLevelNumber();
    if (ln < 1 || ln >= static_cast<int>(d_levels.size()) || !d_levels[ln].active) return;
    Level& L = d_levels[ln];
    auto it = L.work.find(fine.getPatchNumber());
    if (it == L.work.end() || &fine != L.level->getPatch(fine.getPatchNumber()).getPointer()) return;
    computeTangential(ln, it->second, coarse, fine_box, ratio);
}

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CFOwnerFill::setUpLevel(const int ln)
{
    Level& L = d_levels[ln];
    L.level = d_hierarchy->getPatchLevel(ln);
    checkFillWidth(ln);
    BoxList<NDIM> value_boxes;
    computeGeometry(ln, value_boxes);
    if (!L.active) return;
    setUpValueLevel(ln, value_boxes);
    setUpDirect(ln);
    setUpPieces(ln);
    setUpWork(ln);
    setUpSchedules(ln);
}

void
CFOwnerFill::checkFillWidth(const int ln) const
{
    const IntVector<NDIM> ratio = d_levels[ln].level->getRatioToCoarserLevel();
    int axis = -1;
    for (int d = 0; d < NDIM && axis < 0; ++d)
    {
        if (d_fill_width(d) > ratio(d)) axis = d;
    }
    if (axis < 0) return;
    Pointer<Variable<NDIM>> var;
    VariableDatabase<NDIM>::getDatabase()->mapIndexToVariable(d_dst_idx, var);

    TBOX_ERROR("CFOwnerFill::checkFillWidth():\n"
               << "  the quadratic coarse-fine interpolation of variable " << var->getName() << " fills "
               << d_fill_width(axis) << " ghost cell layers along axis " << axis
               << ", but the refinement ratio of level " << ln << " is " << ratio(axis) << "\n"
               << "  the interpolation reaches at most one coarse cell from the coarse-fine interface\n"
               << "  use a ghost width or a fill pattern width of at most the refinement ratio\n");
}

void
CFOwnerFill::setUpSchedules(const int ln)
{
    Level& L = d_levels[ln];
    const bool first = !L.copy_sched;
    if (first)
    {
        // Send the values of each direction and component to the level of values.
        for (auto& o : L.owned)
        {
            // A schedule sets up only the values that cross processes or periodic boundaries.
            std::vector<BoxList<NDIM>> regions(d_is_cell ? 1 : NDIM);
            const std::size_t ci = d_is_cell ? 0 : static_cast<std::size_t>(o.comp);
            regions[ci] = L.cross_regions[ci][o.code];
            if (regions[ci].isEmpty()) continue;
            Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> pattern =
                new RegionFillPattern(new SAMRAI::xfer::BoxGeometryFillPattern<NDIM>(), regions);
            o.alg = new RefineAlgorithm<NDIM>();
            o.alg->registerRefine(
                d_slot_idx[o.code], d_owned_idx, d_slot_idx[o.code], Pointer<RefineOperator<NDIM>>(nullptr), pattern);
            o.sched = o.alg->createSchedule(L.value_level, o.level);
        }
    }

    // Copy the averages into the ghost cells of the level.
    Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> inner = d_fill_pattern;
    if (!inner) inner = new SAMRAI::xfer::BoxGeometryFillPattern<NDIM>();
    Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> pattern = new RegionFillPattern(inner, L.all_regions);
    L.copy_alg = new RefineAlgorithm<NDIM>();
    L.copy_alg->registerRefine(d_dst_idx, d_result_idx, d_dst_idx, Pointer<RefineOperator<NDIM>>(nullptr), pattern);
    if (first)
    {
        L.copy_sched = L.copy_alg->createSchedule(L.level, L.value_level);
    }
    else
    {
        L.copy_alg->resetSchedule(L.copy_sched);
    }
}

void
CFOwnerFill::computeGeometry(const int ln, BoxList<NDIM>& value_boxes)
{
    Level& L = d_levels[ln];
    const BoxArray<NDIM>& boxes = L.level->getBoxes();
    const IntVector<NDIM> ratio = L.level->getRatioToCoarserLevel();
    for (int d = 0; d < NDIM; ++d)
    {
        if (ratio(d) < 2)
        {
            TBOX_ERROR(
                "CFOwnerFill::computeGeometry():\n"
                << "  the coarse-fine interpolation requires a refinement ratio of at least 2 but the ratio of level "
                << ln << " is " << ratio(d) << " along axis " << d << "\n");
        }
    }
    Pointer<BoxTree<NDIM>> tree = L.level->getBoxTree();
    const BoxList<NDIM> domain_cells(L.level->getPhysicalDomain());

    // Offsets of the periodic images.
    const IntVector<NDIM> periodic_shift = d_hierarchy->getGridGeometry()->getPeriodicShift(L.level->getRatio());
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

    // Remove the cells (faces) of the fine level, with periodic images, from a list.
    auto remove_fine = [&](BoxList<NDIM>& list, const int comp)
    {
        for (const auto& offset : offsets)
        {
            for (int pass = 0; pass < (comp < 0 ? 1 : 2); ++pass)
            {
                // A face is in F if the cell above it or the cell below it is.
                IntVector<NDIM> shift = offset;
                if (pass == 1) shift(comp) += 1;
                list.shift(-shift);
                tree->removeIntersections(list);
                list.shift(shift);
            }
        }
    };

    L.regions.assign(d_comps.size(), std::array<BoxList<NDIM>, 2 * NDIM>());
    L.all_regions.assign(d_comps.size(), BoxList<NDIM>());
    // The tasks of every patch, for the levels of owned values.
    struct GlobalTask
    {
        int patch;
        std::size_t ci;
        int code;
        Box<NDIM> index_box;
        BoxList<NDIM> local_parts;   // The cells (faces) that the slab fills, in the index space of the domain.
        BoxList<NDIM> wrapped_parts; // The same cells (faces) beyond a periodic boundary, as cells (faces) inside it.
    };
    std::vector<GlobalTask> tasks;
    for (std::size_t ci = 0; ci < d_comps.size(); ++ci)
    {
        const int comp = d_comps[ci];
        BoxList<NDIM> domain_list;
        for (typename BoxList<NDIM>::Iterator it(domain_cells); it; it++) domain_list.addItem(index_box(*it, comp));
        // The cells (faces) that are in the domain or have a periodic image in it.
        BoxList<NDIM> image_domain;
        for (const auto& offset : offsets)
        {
            BoxList<NDIM> image = domain_list;
            image.shift(-offset);
            image_domain.unionBoxes(image);
        }
        for (int P = 0; P < boxes.getNumberOfBoxes(); ++P)
        {
            const Box<NDIM> patch_box_index = index_box(boxes[P], comp);
            for (int d = 0; d < NDIM; ++d)
            {
                for (const int s : { +1, -1 })
                {
                    // The cells (faces) of the first layer of the patch on the side of c + s e_d, and the cells (faces)
                    // c.
                    Box<NDIM> layer = patch_box_index;
                    if (comp == d)
                        layer.lower(d) = layer.upper(d) = (s > 0 ? patch_box_index.lower(d) : patch_box_index.upper(d));
                    else
                        layer.lower(d) = layer.upper(d) = (s > 0 ? boxes[P].lower(d) : boxes[P].upper(d));
                    if (comp != d && boxes[P].upper(d) - boxes[P].lower(d) < 1)
                    {
                        TBOX_ERROR("CFOwnerFill::computeGeometry():\n"
                                   << "  patch " << boxes[P] << " of level " << ln << " is one cell thick along axis "
                                   << d
                                   << "\n  the coarse-fine interpolation requires patches of at least two cells\n");
                    }
                    // The cells (faces) c at depths 1 to the maximum depth in front of the patch, each with no cell
                    // (face) of F between it and the patch.
                    const int code = location_index(d, s);
                    const int max_depth = std::min({ d_fill_width(d), ratio(d), MAX_DEPTH });
                    IntVector<NDIM> toward(0);
                    toward(d) = -s;
                    std::vector<BoxList<NDIM>> depth_boxes(1);
                    {
                        BoxList<NDIM> chain(layer);
                        for (int k = 1; k <= max_depth; ++k)
                        {
                            chain.shift(toward);
                            remove_fine(chain, comp);
                            chain.intersectBoxes(image_domain);
                            if (chain.isEmpty()) break;
                            depth_boxes.push_back(chain);
                        }
                    }
                    const int num_depths = static_cast<int>(depth_boxes.size()) - 1;
                    if (num_depths == 0) continue;
                    for (int k = 1; k <= num_depths; ++k)
                    {
                        BoxList<NDIM>& phi = depth_boxes[k];
                        if (comp < 0) continue;
                        // A face on the boundary of two patches has two owners: the patch of smaller index.
                        std::set<int> smaller;
                        for (typename BoxList<NDIM>::Iterator it(phi); it; it++)
                        {
                            Box<NDIM> g = *it;
                            g.shift(d, k * s);
                            tbox::Array<int> indices;
                            tree->findOverlapIndices(indices, g);
                            g.shift(comp, -1);
                            tbox::Array<int> indices2;
                            tree->findOverlapIndices(indices2, g);
                            for (int j = 0; j < indices.size(); ++j)
                                if (indices[j] < P) smaller.insert(indices[j]);
                            for (int j = 0; j < indices2.size(); ++j)
                                if (indices2[j] < P) smaller.insert(indices2[j]);
                        }
                        for (const int other : smaller)
                        {
                            Box<NDIM> taken = index_box(boxes[other], comp);
                            taken.shift(d, -k * s);
                            phi.removeIntersections(taken);
                        }
                    }
                    // A ghost cell (face) beyond a periodic boundary gets the value of its periodic image inside the
                    // domain.
                    for (int k = 1; k <= num_depths; ++k)
                    {
                        for (const auto& offset : offsets)
                        {
                            BoxList<NDIM> part = depth_boxes[k];
                            part.shift(offset);
                            part.intersectBoxes(domain_list);
                            for (typename BoxList<NDIM>::Iterator it(part); it; it++) L.regions[ci][code].addItem(*it);
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
                                slab.lower(d) -= (m - 1);
                            else
                                slab.upper(d) += (m - 1);
                            GlobalTask task{ P, ci, code, slab, BoxList<NDIM>(), BoxList<NDIM>() };
                            for (int k = 1; k <= num_depths; ++k)
                            {
                                BoxList<NDIM> at_depth = depth_boxes[k];
                                at_depth.intersectBoxes(slab);
                                if (at_depth.isEmpty()) continue;
                                for (const auto& offset : offsets)
                                {
                                    BoxList<NDIM> part = at_depth;
                                    part.shift(offset);
                                    part.intersectBoxes(domain_list);
                                    BoxList<NDIM>& parts =
                                        (offset == IntVector<NDIM>(0)) ? task.local_parts : task.wrapped_parts;
                                    for (typename BoxList<NDIM>::Iterator pit(part); pit; pit++) parts.addItem(*pit);
                                }
                            }
                            tasks.push_back(task);
                        }
                    }
                }
            }
        }
    }

    // Levels of owned values, one for each direction and component; the patches that a rank computes.
    const ProcessorMapping& fine_mapping = L.level->getProcessorMapping();
    L.owned.clear();
    L.work.clear();
    L.task_info.clear();
    for (std::size_t ci = 0; ci < d_comps.size(); ++ci)
    {
        for (int code = 0; code < 2 * NDIM; ++code)
        {
            int count = 0;
            for (const auto& t : tasks) count += (t.ci == ci && t.code == code);
            if (count == 0) continue;
            OwnedLevel o;
            o.code = code;
            o.comp = d_comps[ci];
            BoxArray<NDIM> owned_boxes(count);
            ProcessorMapping mapping(count);
            int k = 0;
            const int d = code / 2, s = (code % 2 == 0 ? +1 : -1);
            const int owned_level = static_cast<int>(L.owned.size());
            for (const auto& t : tasks)
            {
                if (t.ci != ci || t.code != code) continue;
                const Box<NDIM> kernel = kernel_box(t.index_box, o.comp, d, s, boxes[t.patch]);
                owned_boxes[k] = kernel;
                mapping.setProcessorAssignment(k, fine_mapping.getProcessorAssignment(t.patch));
                L.task_info.push_back({ ci,
                                        code,
                                        owned_level,
                                        k,
                                        fine_mapping.getProcessorAssignment(t.patch),
                                        t.local_parts,
                                        t.wrapped_parts });
                if (fine_mapping.isMappingLocal(t.patch))
                {
                    PatchWork& pw = L.work[t.patch];
                    pw.box = boxes[t.patch];
                    pw.tasks.push_back({ o.comp, d, s, t.index_box, kernel, owned_level, k, nullptr });
                }
                ++k;
            }
            o.level = new PatchLevel<NDIM>(owned_boxes,
                                           mapping,
                                           L.level->getRatio(),
                                           d_hierarchy->getGridGeometry(),
                                           VariableDatabase<NDIM>::getDatabase()->getPatchDescriptor(),
                                           nullptr,
                                           true);
            o.level->allocatePatchData(d_owned_idx);
            for (auto& entry : L.work)
            {
                for (auto& task : entry.second.tasks)
                {
                    if (task.owned_level != owned_level) continue;
                    task.owned_array =
                        &array_data(o.level->getPatch(task.owned_patch)->getPatchData(d_owned_idx), task.comp);
                }
            }
            L.owned.push_back(o);
        }
    }

    // The cells (faces) that get values, and the strips of cells that contain them. A face belongs to the cell above
    // it, or below it if that cell is outside the domain.
    BoxList<NDIM> union_list;
    for (std::size_t ci = 0; ci < d_comps.size(); ++ci)
    {
        const int comp = d_comps[ci];
        for (int code = 0; code < 2 * NDIM; ++code)
        {
            BoxList<NDIM> cells;
            for (typename BoxList<NDIM>::Iterator it(L.regions[ci][code]); it; it++)
            {
                cells.addItem(*it);
                if (comp >= 0)
                {
                    Box<NDIM> below = *it;
                    below.shift(comp, -1);
                    cells.addItem(below);
                }
                for (const auto& offset : offsets)
                {
                    Box<NDIM> image = *it;
                    image.shift(offset);
                    L.all_regions[ci].addItem(image);
                }
            }
            cells.intersectBoxes(domain_cells);
            cells.simplifyBoxes();
            if (!union_list.isEmpty())
            {
                const BoxArray<NDIM> union_array(union_list);
                BoxTree<NDIM> union_tree(union_array);
                union_tree.removeIntersections(cells);
            }
            union_list.unionBoxes(cells);
        }
    }
    union_list.coalesceBoxes();
    value_boxes = union_list;
    L.active = !union_list.isEmpty();
}

void
CFOwnerFill::setUpValueLevel(const int ln, const BoxList<NDIM>& value_boxes)
{
    Level& L = d_levels[ln];
    const ProcessorMapping& fine_mapping = L.level->getProcessorMapping();
    Pointer<BoxTree<NDIM>> tree = L.level->getBoxTree();
    const int num_boxes = value_boxes.getNumberOfItems();
    BoxArray<NDIM> boxes(num_boxes);
    ProcessorMapping mapping(num_boxes);
    int k = 0;
    for (typename BoxList<NDIM>::Iterator it(value_boxes); it; it++, ++k)
    {
        boxes[k] = *it;
        // A patch of values belongs to the process of the fine patch of smallest index next to it.
        tbox::Array<int> indices;
        tree->findOverlapIndices(indices, Box<NDIM>::grow(*it, IntVector<NDIM>(1)));
        int owner_patch = indices.size() > 0 ? indices[0] : 0;
        for (int j = 1; j < indices.size(); ++j) owner_patch = std::min(owner_patch, indices[j]);
        mapping.setProcessorAssignment(k, fine_mapping.getProcessorAssignment(owner_patch));
    }
    L.value_level = new PatchLevel<NDIM>(boxes,
                                         mapping,
                                         L.level->getRatio(),
                                         d_hierarchy->getGridGeometry(),
                                         VariableDatabase<NDIM>::getDatabase()->getPatchDescriptor(),
                                         nullptr,
                                         true);
    for (int code = 0; code < 2 * NDIM; ++code) L.value_level->allocatePatchData(d_slot_idx[code]);
    L.value_level->allocatePatchData(d_result_idx);
}

void
CFOwnerFill::setUpDirect(const int ln)
{
    Level& L = d_levels[ln];
    const BoxArray<NDIM>& value_boxes = L.value_level->getBoxes();
    const ProcessorMapping& value_mapping = L.value_level->getProcessorMapping();
    BoxTree<NDIM> value_tree(value_boxes);
    L.cross_regions.assign(d_comps.size(), std::array<BoxList<NDIM>, 2 * NDIM>());
    L.direct.clear();
    for (const TaskInfo& t : L.task_info)
    {
        const int comp = d_comps[t.ci];
        // The values that cross a periodic boundary go through a schedule.
        for (typename BoxList<NDIM>::Iterator it(t.wrapped_parts); it; it++) L.cross_regions[t.ci][t.code].addItem(*it);
        for (typename BoxList<NDIM>::Iterator it(t.local_parts); it; it++)
        {
            Box<NDIM> query = *it;
            if (comp >= 0) query.lower(comp) -= 1;
            tbox::Array<int> indices;
            value_tree.findOverlapIndices(indices, query);
            for (int k = 0; k < indices.size(); ++k)
            {
                const int j = indices[k];
                const Box<NDIM> piece = *it * index_box(value_boxes[j], comp);
                if (piece.empty()) continue;
                if (value_mapping.getProcessorAssignment(j) != t.rank)
                {
                    L.cross_regions[t.ci][t.code].addItem(piece);
                }
                else if (value_mapping.isMappingLocal(j))
                {
                    // Both patches are on this process: copy directly.
                    const OwnedLevel& o = L.owned[t.owned_level];
                    L.direct.push_back(
                        { &array_data(L.value_level->getPatch(j)->getPatchData(d_slot_idx[t.code]), comp),
                          &array_data(o.level->getPatch(t.owned_patch)->getPatchData(d_owned_idx), comp),
                          piece });
                }
            }
        }
    }
    L.task_info.clear();
}

void
CFOwnerFill::setUpPieces(const int ln)
{
    Level& L = d_levels[ln];
    L.value_patches.clear();
    const BoxArray<NDIM>& value_boxes = L.value_level->getBoxes();
    const int num_value_boxes = value_boxes.getNumberOfBoxes();
    BoxTree<NDIM> value_tree(value_boxes);

    // The part of each region in each value patch, found with the tree.
    std::vector<std::vector<std::array<BoxList<NDIM>, 2 * NDIM>>> local(num_value_boxes);
    for (std::size_t ci = 0; ci < d_comps.size(); ++ci)
    {
        for (int code = 0; code < 2 * NDIM; ++code)
        {
            for (typename BoxList<NDIM>::Iterator it(L.regions[ci][code]); it; it++)
            {
                // A face belongs to the patch whose box contains the cell above it or the cell below it.
                Box<NDIM> query = *it;
                if (d_comps[ci] >= 0) query.lower(d_comps[ci]) -= 1;
                tbox::Array<int> indices;
                value_tree.findOverlapIndices(indices, query);
                for (int k = 0; k < indices.size(); ++k)
                {
                    const int j = indices[k];
                    if (!L.value_level->getProcessorMapping().isMappingLocal(j)) continue;
                    const Box<NDIM> part = *it * index_box(value_boxes[j], d_comps[ci]);
                    if (part.empty()) continue;
                    if (local[j].empty()) local[j].resize(d_comps.size());
                    local[j][ci][code].addItem(part);
                }
            }
        }
    }

    for (PatchLevel<NDIM>::Iterator p(L.value_level); p; p++)
    {
        ValuePatch vp;
        vp.pieces.resize(d_comps.size());
        vp.result.resize(d_comps.size());
        vp.slots.resize(d_comps.size());
        Pointer<Patch<NDIM>> patch = L.value_level->getPatch(p());
        for (std::size_t ci = 0; ci < d_comps.size(); ++ci)
        {
            vp.result[ci] = &array_data(patch->getPatchData(d_result_idx), d_comps[ci]);
            for (int code = 0; code < 2 * NDIM; ++code)
                vp.slots[ci][code] = &array_data(patch->getPatchData(d_slot_idx[code]), d_comps[ci]);
            if (local[p()].empty()) continue;
            std::vector<Piece> pieces;
            for (int code = 0; code < 2 * NDIM; ++code)
            {
                const BoxList<NDIM>& region = local[p()][ci][code];
                if (region.isEmpty()) continue;
                const int bit = 1 << code;
                std::vector<Piece> next;
                BoxList<NDIM> fresh = region;
                for (const auto& piece : pieces)
                {
                    BoxList<NDIM> inside(piece.box);
                    inside.intersectBoxes(region);
                    BoxList<NDIM> outside(piece.box);
                    outside.removeIntersections(region);
                    for (typename BoxList<NDIM>::Iterator it(inside); it; it++)
                        next.push_back({ *it, piece.mask | bit });
                    for (typename BoxList<NDIM>::Iterator it(outside); it; it++) next.push_back({ *it, piece.mask });
                    fresh.removeIntersections(piece.box);
                }
                for (typename BoxList<NDIM>::Iterator it(fresh); it; it++) next.push_back({ *it, bit });
                pieces = next;
            }
            vp.pieces[ci] = pieces;
        }
        L.value_patches.push_back(vp);
    }
}

void
CFOwnerFill::setUpWork(const int ln)
{
    Level& L = d_levels[ln];
    for (auto& entry : L.work)
    {
        PatchWork& pw = entry.second;
        for (const auto& task : pw.tasks)
        {
            const int code = location_index(task.d, task.s);
            if (pw.tangential[code]) continue;
            if (d_is_cell)
                pw.tangential[code] = new CellData<NDIM, double>(pw.box, d_depth, IntVector<NDIM>(d_work_gcw));
            else
                pw.tangential[code] = new SideData<NDIM, double>(pw.box, d_depth, IntVector<NDIM>(d_work_gcw));
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
}

void
CFOwnerFill::computeTangential(const int ln,
                               PatchWork& pw,
                               const Patch<NDIM>& coarse,
                               const Box<NDIM>& fine_box,
                               const IntVector<NDIM>& ratio)
{
    Pointer<PatchData<NDIM>> cdata = coarse.getPatchData(d_dst_idx);
    const Box<NDIM>& cbox = coarse.getBox();
    const Box<NDIM>& fbox = pw.box;
    const int crse_gcw = ghost_width(cdata);
    for (const auto& task : pw.tasks)
    {
        // The coarse patch covers the cells of the fine box only.
        const Box<NDIM> kernel = task.kernel_box * fine_box;
        if (kernel.empty()) continue;
#if !defined(NDEBUG)
        {
            // The kernel uses the coarse cells of the kernel box and their neighbors along the tangential directions.
            // Every one of them is in the coarse level or outside the physical domain, where the physical boundary
            // fill sets it. The coarse patch of the schedule sets no other cell.
            Pointer<PatchLevel<NDIM>> coarse_level = d_hierarchy->getPatchLevel(ln - 1);
            Box<NDIM> used = Box<NDIM>::coarsen(kernel, ratio);
            for (int d = 0; d < NDIM; ++d)
            {
                if (d != task.d) used.grow(d, 1);
            }
            const IntVector<NDIM> shift = d_hierarchy->getGridGeometry()->getPeriodicShift(coarse_level->getRatio());
            const BoxList<NDIM> coarse_domain(coarse_level->getPhysicalDomain());
            BoxList<NDIM> unset;
            for (int k = 0; k < (NDIM == 2 ? 9 : 27); ++k)
            {
                IntVector<NDIM> offset(0);
                int rem = k;
                bool allowed = true;
                for (int d = 0; d < NDIM; ++d)
                {
                    const int m = rem % 3 - 1;
                    rem /= 3;
                    if (m != 0 && shift(d) == 0) allowed = false;
                    offset(d) = m * shift(d);
                }
                if (!allowed) continue;
                Box<NDIM> image = used;
                image.shift(offset);
                BoxList<NDIM> in_domain(image);
                in_domain.intersectBoxes(coarse_domain);
                unset.unionBoxes(in_domain);
            }
            coarse_level->getBoxTree()->removeIntersections(unset);
            if (!unset.isEmpty())
            {
                TBOX_ERROR("CFOwnerFill::computeTangential():\n"
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
}

void
CFOwnerFill::computeNormal(PatchWork& pw, Patch<NDIM>& patch, const IntVector<NDIM>& ratio)
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
            first.lower(task.d) = first.upper(task.d);
        else
            first.upper(task.d) = first.lower(task.d);
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
        // The values go to the patch of owned values.
        copy_on(*task.owned_array, d_is_cell ? array_data(pw.work, -1) : array_data(pw.result, task.comp), filled);
    }
}

void
CFOwnerFill::computeValues(const int ln)
{
    Level& L = d_levels[ln];
    for (PatchLevel<NDIM>::Iterator p(L.level); p; p++)
    {
        auto it = L.work.find(p());
        if (it == L.work.end()) continue;
        Pointer<Patch<NDIM>> patch = L.level->getPatch(p());
        computeNormal(it->second, *patch, L.level->getRatioToCoarserLevel());
    }
    for (const DirectCopy& c : L.direct) c.dst->copy(*c.src, c.box);
    for (auto& o : L.owned)
    {
        if (o.sched) o.sched->fillData(0.0);
    }

    // The value is the average of the values of the directions, summed in the order of the directions.
    for (const ValuePatch& vp : L.value_patches)
    {
        for (std::size_t ci = 0; ci < d_comps.size(); ++ci)
        {
            ArrayData<NDIM, double>& result = *vp.result[ci];
            for (const auto& piece : vp.pieces[ci])
            {
                int count = 0;
                std::array<const ArrayData<NDIM, double>*, 2 * NDIM> slots;
                for (int code = 0; code < 2 * NDIM; ++code)
                {
                    if (piece.mask & (1 << code)) slots[count++] = vp.slots[ci][code];
                }
                const double n = static_cast<double>(count);
                for (BoxIterator<NDIM> bi(piece.box); bi; bi++)
                {
                    const hier::Index<NDIM>& i = bi();
                    for (int depth = 0; depth < d_depth; ++depth)
                    {
                        double sum = 0.0;
                        for (int k = 0; k < count; ++k) sum += (*slots[k])(i, depth);
                        result(i, depth) = sum / n;
                    }
                }
            }
        }
    }
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
