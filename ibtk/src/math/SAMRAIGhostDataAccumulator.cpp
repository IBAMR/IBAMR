// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2026 by the IBAMR developers
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
#include <ibtk/RobinPhysBdryPatchStrategy.h>
#include <ibtk/SAMRAIGhostDataAccumulator.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/AbstractStream.h>
#include <tbox/Array.h>
#include <tbox/Pointer.h>
#include <tbox/SAMRAI_MPI.h>
#include <tbox/Schedule.h>
#include <tbox/Transaction.h>
#include <tbox/Utilities.h>

#include <ArrayData.h>
#include <BasePatchHierarchy.h>
#include <Box.h>
#include <BoxArray.h>
#include <BoxList.h>
#include <BoxOverlap.h>
#include <BoxTree.h>
#include <CellData.h>
#include <CellDataFactory.h>
#include <CellOverlap.h>
#include <CellVariable.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchData.h>
#include <PatchLevel.h>
#include <RefineAlgorithm.h>
#include <RefineOperator.h>
#include <RefineSchedule.h>
#include <RefineTransactionFactory.h>
#include <SideData.h>
#include <SideDataFactory.h>
#include <SideGeometry.h>
#include <SideOverlap.h>
#include <SideVariable.h>
#include <Variable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <array>
#include <memory>
#include <vector>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////
namespace
{
static Timer* t_constructor;
static Timer* t_accumulate_ghost_data;

// Centerings of the data that can be accumulated. A transaction treats the
// values of the patch data as a list of arrays, one for each component of the
// centering.
enum class Centering
{
    CELL,
    SIDE
};

int
num_components(const Centering centering)
{
    return centering == Centering::CELL ? 1 : NDIM;
}

ArrayData<NDIM, double>&
get_array_data(const Centering centering, PatchData<NDIM>& data, const int component)
{
    if (centering == Centering::CELL)
    {
        return static_cast<CellData<NDIM, double>&>(data).getArrayData();
    }
    return static_cast<SideData<NDIM, double>&>(data).getArrayData(component);
}

const BoxList<NDIM>&
get_overlap_boxes(const Centering centering, const BoxOverlap<NDIM>& overlap, const int component)
{
    if (centering == Centering::CELL)
    {
        return static_cast<const CellOverlap<NDIM>&>(overlap).getDestinationBoxList();
    }
    return static_cast<const SideOverlap<NDIM>&>(overlap).getDestinationBoxList(component);
}

const IntVector<NDIM>&
get_overlap_offset(const Centering centering, const BoxOverlap<NDIM>& overlap)
{
    if (centering == Centering::CELL)
    {
        return static_cast<const CellOverlap<NDIM>&>(overlap).getSourceOffset();
    }
    return static_cast<const SideOverlap<NDIM>&>(overlap).getSourceOffset();
}

// State that a schedule and its transactions share. It is set before each use of
// the schedule.
struct ScheduleState
{
    // Patch data index of the data that receive the sums.
    int target_idx = invalid_index;

    // Flags that tell which patches of the level may hold nonzero values.
    bool enabled = false;
    std::vector<char> active; // indexed by patch number
};

/*!
 * Transaction of a schedule that sums values into the interiors of patches.
 *
 * The schedule is created so that its destination patches are the patches that
 * hold values in their ghost regions, in the patch data with index d_values_idx,
 * and its source patches are the patches that own the values. The refine
 * schedule then finds, for each destination patch p, the source patches q whose
 * data boxes intersect the ghost data box of p, but a refine transaction moves
 * data from q to p. This transaction moves the values in the opposite direction,
 * from p to q, and adds them to the patch data of q with the index
 * ScheduleState::target_idx. The class tbox::Schedule decides from the source
 * and destination processors of a transaction which processor packs and which
 * unpacks the data, so those are reported in the opposite order.
 *
 * The schedule is registered with the same patch data index for its source and
 * destination, so the refine schedule does not create a transaction that moves
 * data of a patch to itself without a shift; the periodic images of a patch are
 * different locations of its ghost region and do get a transaction.
 *
 * The overlap boxes are in the index space of p. The boxes in the index space of
 * q are obtained with the negative of the source offset of the overlap. The same
 * transaction moves the values that the refine schedule does not move (see
 * compute_gap_transfers()), with overlaps that are built from the patch boxes.
 */
class GhostSumTransaction : public tbox::Transaction
{
public:
    GhostSumTransaction(const Centering centering,
                        const int depth,
                        const int values_idx,
                        Pointer<PatchLevel<NDIM>> level,
                        const BoxOverlap<NDIM>& overlap,
                        const int p,
                        const int q,
                        std::shared_ptr<const ScheduleState> state)
        : d_centering(centering),
          d_values_idx(values_idx),
          d_level(level),
          d_p(p),
          d_q(q),
          d_offset(get_overlap_offset(centering, overlap)),
          d_state(state)
    {
        d_num_values = 0;
        for (int component = 0; component < num_components(centering); ++component)
        {
            d_boxes_p[component] = get_overlap_boxes(centering, overlap, component);
            d_boxes_q[component] = d_boxes_p[component];
            d_boxes_q[component].shift(-d_offset);
            d_num_values += d_boxes_p[component].getTotalSizeOfBoxes();
        }
        d_num_values *= depth;
    }

    ~GhostSumTransaction() override = default;

    bool canEstimateIncomingMessageSize() override
    {
        return true;
    }

    int computeIncomingMessageSize() override
    {
        return skip() ? 0 : static_cast<int>(d_num_values * sizeof(double));
    }

    int computeOutgoingMessageSize() override
    {
        return computeIncomingMessageSize();
    }

    int getSourceProcessor() override
    {
        return d_level->getMappingForPatch(d_p);
    }

    int getDestinationProcessor() override
    {
        return d_level->getMappingForPatch(d_q);
    }

    void packStream(tbox::AbstractStream& stream) override
    {
        if (skip())
        {
            return;
        }
        Pointer<PatchData<NDIM>> data = d_level->getPatch(d_p)->getPatchData(d_values_idx);
        for (int component = 0; component < num_components(d_centering); ++component)
        {
            if (d_boxes_p[component].getNumberOfItems() == 0)
            {
                continue;
            }
            get_array_data(d_centering, *data, component).packStream(stream, d_boxes_p[component], IntVector<NDIM>(0));
        }
    }

    void unpackStream(tbox::AbstractStream& stream) override
    {
        if (skip())
        {
            return;
        }
        Pointer<PatchData<NDIM>> data = d_level->getPatch(d_q)->getPatchData(d_state->target_idx);
        for (int component = 0; component < num_components(d_centering); ++component)
        {
            if (d_boxes_q[component].getNumberOfItems() == 0)
            {
                continue;
            }
            get_array_data(d_centering, *data, component)
                .unpackStreamAndSum(stream, d_boxes_q[component], IntVector<NDIM>(0));
        }
    }

    void copyLocalData() override
    {
        if (skip())
        {
            return;
        }
        Pointer<PatchData<NDIM>> values = d_level->getPatch(d_p)->getPatchData(d_values_idx);
        Pointer<PatchData<NDIM>> target = d_level->getPatch(d_q)->getPatchData(d_state->target_idx);
        for (int component = 0; component < num_components(d_centering); ++component)
        {
            if (d_boxes_q[component].getNumberOfItems() == 0)
            {
                continue;
            }
            get_array_data(d_centering, *target, component)
                .sum(get_array_data(d_centering, *values, component), d_boxes_q[component], -d_offset);
        }
    }

    void printClassData(std::ostream& stream) const override
    {
        stream << "IBTK ghost sum transaction: values on patch " << d_p << ", target patch " << d_q << ", offset "
               << d_offset << "\n";
    }

private:
    bool skip() const
    {
        return d_state->enabled && !d_state->active[d_p];
    }

    const Centering d_centering;
    const int d_values_idx;
    Pointer<PatchLevel<NDIM>> d_level;
    const int d_p, d_q;
    const IntVector<NDIM> d_offset;
    std::shared_ptr<const ScheduleState> d_state;
    long d_num_values = 0;
    BoxList<NDIM> d_boxes_p[NDIM], d_boxes_q[NDIM];
};

class GhostSumTransactionFactory : public RefineTransactionFactory<NDIM>
{
public:
    GhostSumTransactionFactory(const Centering centering,
                               const int depth,
                               const int values_idx,
                               std::shared_ptr<const ScheduleState> state)
        : d_centering(centering), d_depth(depth), d_values_idx(values_idx), d_state(state)
    {
    }

    // The patch data indices are fixed by the constructor and by the state, so
    // the items of the refine algorithm are not needed.
    void setRefineItems(const RefineClasses<NDIM>::Data** /*items*/, const int /*num_items*/) override
    {
    }

    void unsetRefineItems() override
    {
    }

    Pointer<tbox::Transaction> allocate(Pointer<PatchLevel<NDIM>> dst_level,
                                        Pointer<PatchLevel<NDIM>> src_level,
                                        Pointer<BoxOverlap<NDIM>> overlap,
                                        const int dst_patch_id,
                                        const int src_patch_id,
                                        const int /*ritem_id*/,
                                        const Box<NDIM>& /*box*/,
                                        const bool /*use_time_interpolation*/,
                                        Pointer<tbox::Arena> /*pool*/) const override
    {
        TBOX_ASSERT(dst_level == src_level);
        return Pointer<tbox::Transaction>(new GhostSumTransaction(
            d_centering, d_depth, d_values_idx, dst_level, *overlap, dst_patch_id, src_patch_id, d_state));
    }

private:
    const Centering d_centering;
    const int d_depth;
    const int d_values_idx;
    std::shared_ptr<const ScheduleState> d_state;
};

using SideBoxLists = std::array<BoxList<NDIM>, NDIM>;

// Values that patch p holds at sides of its ghost side box and that the refine schedule of the level does not move to
// the patch q that owns the sides.
struct GapTransfer
{
    int p = -1;
    int q = -1;
    IntVector<NDIM> shift = IntVector<NDIM>(0); // q is shifted by this to be near p
    SideBoxLists boxes;                         // in the index space of p
};

/*!
 * Find the sides of the ghost side box of a patch p that lie in the side box of another patch q (or of a periodic
 * image of q), but that the refine schedule of the level does not move to q. The result has an entry for each pair
 * of patches and periodic shift with such sides, for the pairs in which p or q is owned by this process. The entries
 * are sorted by patch numbers and shift, so that the processes that own p and q list their common entries in the same
 * order.
 *
 * The schedule moves values from a patch p to a patch q, or to a periodic image of q that SAMRAI gives to q (only
 * patches that touch a periodic boundary have images), if the cell box of q intersects the ghost cell box of p, and
 * then moves the values at the locations in the intersection of the ghost side box of p and the side box of the
 * intersection of the two cell boxes. A side in the ghost side box of p and the side box of q that is not in that set
 * is on the outermost layer of the ghost side box of p, and q is separated from the ghost cell box of p by a zero
 * width. If the ghost cell width is zero, the schedule makes the intersection with the cell box of q after growing the
 * ghost cell box of p by one. The sides of every periodic image of q, also an image that SAMRAI does not give to q,
 * are considered here; the schedule moves nothing to such an image.
 */
std::vector<GapTransfer>
compute_gap_transfers(Pointer<PatchLevel<NDIM>> level, const IntVector<NDIM>& gcw)
{
    const BoxArray<NDIM>& patch_boxes = level->getBoxes();
    Pointer<BoxTree<NDIM>> box_tree = level->getBoxTree();
    const IntVector<NDIM> no_shift(0);
    const IntVector<NDIM> one(1);
    const bool zero_gcw = (gcw == no_shift);
    const int rank = IBTK_MPI::getRank();
    std::vector<GapTransfer> transfers;

    // All shifts that map a patch to its periodic images, and the test whether the schedule uses a shift for a patch.
    const IntVector<NDIM> period = level->getGridGeometry()->getPeriodicShift(level->getRatio());
    std::vector<IntVector<NDIM>> all_shifts(1, no_shift);
    for (int d = 0; d < NDIM; ++d)
    {
        if (period(d) == 0)
        {
            continue;
        }
        const std::size_t n = all_shifts.size();
        for (std::size_t k = 0; k < n; ++k)
        {
            for (const int sign : { -1, 1 })
            {
                IntVector<NDIM> shift = all_shifts[k];
                shift(d) += sign * period(d);
                all_shifts.push_back(shift);
            }
        }
    }
    auto is_schedule_shift = [&](const int q, const IntVector<NDIM>& shift)
    {
        if (shift == no_shift)
        {
            return true;
        }
        for (tbox::List<IntVector<NDIM>>::Iterator s(level->getShiftsForPatch(q)); s; s++)
        {
            if (s() == shift)
            {
                return true;
            }
        }
        return false;
    };

    auto find_sides = [&](const int p, const int q, const IntVector<NDIM>& shift)
    {
        if (q == p && shift == no_shift)
        {
            return;
        }
        const Box<NDIM> ghost_box = Box<NDIM>::grow(patch_boxes[p], gcw);
        const Box<NDIM> wide_ghost_box = Box<NDIM>::grow(ghost_box, one);
        const Box<NDIM> shifted = Box<NDIM>::shift(patch_boxes[q], shift);
        if ((wide_ghost_box * shifted).empty())
        {
            return;
        }

        Box<NDIM> moved;
        if (is_schedule_shift(q, shift))
        {
            moved = ghost_box * shifted;
            if (moved.empty() && zero_gcw)
            {
                moved = wide_ghost_box * shifted;
            }
        }
        GapTransfer transfer;
        transfer.p = p;
        transfer.q = q;
        transfer.shift = shift;
        bool found = false;
        for (int axis = 0; axis < NDIM; ++axis)
        {
            const Box<NDIM> ghost_side_box = SideGeometry<NDIM>::toSideBox(ghost_box, axis);
            const Box<NDIM> shared = ghost_side_box * SideGeometry<NDIM>::toSideBox(shifted, axis);
            if (shared.empty())
            {
                continue;
            }
            BoxList<NDIM> missing(shared);
            if (!moved.empty())
            {
                missing.removeIntersections(ghost_side_box * SideGeometry<NDIM>::toSideBox(moved, axis));
            }
            for (BoxList<NDIM>::Iterator b(missing); b; b++)
            {
                transfer.boxes[axis].appendItem(b());
                found = true;
            }
        }
        if (found)
        {
            transfers.push_back(transfer);
        }
    };

    // The box tree holds the patches together with their periodic images, so it finds every patch with an image that
    // intersects the box. The indices are sorted and a patch may appear more than once.
    tbox::Array<int> candidates;
    auto for_each_candidate = [&](const Box<NDIM>& region, const auto& function)
    {
        box_tree->findOverlapIndices(candidates, region);
        int previous = -1;
        for (int k = 0; k < candidates.size(); ++k)
        {
            if (candidates[k] != previous)
            {
                previous = candidates[k];
                function(previous);
            }
        }
    };

    // Pairs in which this process owns p.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        const Box<NDIM> region = Box<NDIM>::grow(patch_boxes[p()], gcw + one);
        for (const IntVector<NDIM>& shift : all_shifts)
        {
            for_each_candidate(Box<NDIM>::shift(region, -shift), [&](const int q) { find_sides(p(), q, shift); });
        }
    }

    // Pairs in which this process owns q but not p.
    for (PatchLevel<NDIM>::Iterator q(level); q; q++)
    {
        for (const IntVector<NDIM>& shift : all_shifts)
        {
            for_each_candidate(Box<NDIM>::grow(Box<NDIM>::shift(patch_boxes[q()], shift), gcw + one),
                               [&](const int p)
                               {
                                   if (level->getMappingForPatch(p) != rank)
                                   {
                                       find_sides(p, q(), shift);
                                   }
                               });
        }
    }

    std::sort(transfers.begin(),
              transfers.end(),
              [](const GapTransfer& a, const GapTransfer& b)
              {
                  if (a.p != b.p)
                  {
                      return a.p < b.p;
                  }
                  if (a.q != b.q)
                  {
                      return a.q < b.q;
                  }
                  for (int d = 0; d < NDIM; ++d)
                  {
                      if (a.shift(d) != b.shift(d))
                      {
                          return a.shift(d) < b.shift(d);
                      }
                  }
                  return false;
              });
    return transfers;
}

/*!
 * Emit an error if the patches of the level are too narrow for the ghost cell
 * width, in a direction in which the level has patches next to each other or is
 * periodic.
 *
 * The periodic images of a patch are given only to patches that touch a periodic
 * boundary, and the transpose of the physical boundary fill acts only on patches
 * that touch the physical boundary. A patch that is narrower than the ghost
 * cell width has a ghost region that extends past its neighbor to patches that
 * do not touch the boundary, so values would be lost.
 */
void
check_patch_widths(Pointer<PatchLevel<NDIM>> level, const IntVector<NDIM>& gcw)
{
    const BoxArray<NDIM>& patch_boxes = level->getBoxes();
    const IntVector<NDIM> period = level->getGridGeometry()->getPeriodicShift(level->getRatio());
    for (int d = 0; d < NDIM; ++d)
    {
        // Patches next to each other in direction d have different extents in
        // that direction.
        bool check = (period(d) != 0);
        for (int p = 1; p < patch_boxes.getNumberOfBoxes() && !check; ++p)
        {
            check = patch_boxes[p].lower(d) != patch_boxes[0].lower(d) ||
                    patch_boxes[p].upper(d) != patch_boxes[0].upper(d);
        }
        if (!check)
        {
            continue;
        }
        for (int p = 0; p < patch_boxes.getNumberOfBoxes(); ++p)
        {
            if (patch_boxes[p].numberCells(d) < gcw(d))
            {
                TBOX_ERROR("SAMRAIGhostDataAccumulator:\n"
                           << "  patch " << p << " of level " << level->getLevelNumber() << " (box " << patch_boxes[p]
                           << ") has width " << patch_boxes[p].numberCells(d) << " in direction " << d
                           << ", which is less than the ghost cell width " << gcw(d) << ".\n"
                           << "  Patches must be at least as wide as the ghost cell width in every direction in "
                              "which the level has patches next to each other or is periodic.\n");
            }
        }
    }
}
} // namespace

/////////////////////////////// NESTED TYPES /////////////////////////////////

struct SAMRAIGhostDataAccumulator::LevelInfo
{
    // The level that holds the scratch patch data.
    Pointer<PatchLevel<NDIM>> level;

    // Set when the level is first used.
    bool initialized = false;

    // The schedule of the level, built when the level is first used.
    Pointer<RefineAlgorithm<NDIM>> algorithm;
    Pointer<RefineSchedule<NDIM>> schedule;
    std::shared_ptr<ScheduleState> state;

    // For side-centered data, the schedule that moves the values that the schedule above does not move. It is null if
    // this process sends and receives no such values.
    Pointer<tbox::Schedule> gap_schedule;
};

/////////////////////////////// PUBLIC ///////////////////////////////////////

SAMRAIGhostDataAccumulator::SAMRAIGhostDataAccumulator(Pointer<BasePatchHierarchy<NDIM>> patch_hierarchy,
                                                       Pointer<Variable<NDIM>> var,
                                                       const IntVector<NDIM> gcw,
                                                       const int coarsest_ln,
                                                       const int finest_ln)
    : d_hierarchy(patch_hierarchy), d_var(var), d_gcw(gcw), d_coarsest_ln(coarsest_ln), d_finest_ln(finest_ln)
{
    auto set_timer = [&](const char* name) { return TimerManager::getManager()->getTimer(name); };
    t_constructor = set_timer("IBTK::SAMRAIGhostDataAccumulator::SAMRAIGhostDataAccumulator()");
    t_accumulate_ghost_data = set_timer("IBTK::SAMRAIGhostDataAccumulator::accumulateGhostData()");

    IBTK_TIMER_START(t_constructor);
    Pointer<CellVariable<NDIM, double>> cc_var = var;
    Pointer<SideVariable<NDIM, double>> sc_var = var;
    if (!cc_var && !sc_var)
    {
        TBOX_ERROR("SAMRAIGhostDataAccumulator::SAMRAIGhostDataAccumulator():\n"
                   << "  only cell-centered and side-centered variables with double precision values are supported.\n");
    }
    d_cc_data = cc_var;

    // Register a patch data index for a copy of the values that are summed. The
    // index is private to this object.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("SAMRAIGhostDataAccumulator");
    const int shared_idx = var_db->registerVariableAndContext(var, context, d_gcw);
    d_scratch_idx = var_db->registerClonedPatchDataIndex(var, shared_idx);
    var_db->removePatchDataIndex(shared_idx);

    d_level_info.resize(d_finest_ln + 1);
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_scratch_idx);
        d_level_info[ln] = std::make_unique<LevelInfo>();
        d_level_info[ln]->level = level;
    }
    IBTK_TIMER_STOP(t_constructor);
}

SAMRAIGhostDataAccumulator::~SAMRAIGhostDataAccumulator()
{
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        d_level_info[ln]->level->deallocatePatchData(d_scratch_idx);
    }
    VariableDatabase<NDIM>::getDatabase()->removePatchDataIndex(d_scratch_idx);
}

void
SAMRAIGhostDataAccumulator::accumulateGhostData(const int idx)
{
    accumulateGhostData(idx, nullptr, 0.0);
}

void
SAMRAIGhostDataAccumulator::accumulateGhostData(const int idx,
                                                RobinPhysBdryPatchStrategy* const bdry_op,
                                                const double fill_time,
                                                const std::vector<std::vector<int>>* const active_patch_nums)
{
    IBTK_TIMER_START(t_accumulate_ghost_data);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM>> var;
    var_db->mapIndexToVariable(idx, var);
    TBOX_ASSERT(var == d_var);

    const int num_levels = d_finest_ln - d_coarsest_ln + 1;

    // Determine the patches that may hold nonzero values. Every process needs
    // to know these for all patches since transactions are skipped consistently
    // by the sending and the receiving process.
    std::vector<int> active_flags;
    std::vector<int> level_offsets(num_levels + 1, 0);
    if (active_patch_nums)
    {
        TBOX_ASSERT(static_cast<int>(active_patch_nums->size()) == num_levels);
        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            level_offsets[ln - d_coarsest_ln + 1] =
                level_offsets[ln - d_coarsest_ln] + d_level_info[ln]->level->getNumberOfPatches();
        }
        active_flags.assign(level_offsets[num_levels], 0);
        for (int l = 0; l < num_levels; ++l)
        {
            for (const int patch_num : (*active_patch_nums)[l])
            {
                active_flags[level_offsets[l] + patch_num] = 1;
            }
        }
        IBTK_MPI::sumReduction(active_flags.data(), static_cast<int>(active_flags.size()));
    }

    if (bdry_op)
    {
        bdry_op->setPatchDataIndex(idx);
    }
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        const int l = ln - d_coarsest_ln;
        LevelInfo& level_info = *d_level_info[ln];
        Pointer<PatchLevel<NDIM>> level = level_info.level;

        // A level without active patches has nothing to accumulate.
        if (active_patch_nums && std::none_of(active_flags.begin() + level_offsets[l],
                                              active_flags.begin() + level_offsets[l + 1],
                                              [](const int flag) { return flag > 0; }))
        {
            continue;
        }
        initializeLevel(ln);

        ScheduleState& state = *level_info.state;
        state.target_idx = idx;
        state.enabled = static_cast<bool>(active_patch_nums);
        if (active_patch_nums)
        {
            state.active.assign(level->getNumberOfPatches(), 0);
            for (int n = 0; n < level->getNumberOfPatches(); ++n)
            {
                state.active[n] = active_flags[level_offsets[l] + n] > 0;
            }
        }

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<PatchData<NDIM>> data = patch->getPatchData(idx);
            TBOX_ASSERT(d_gcw == data->getGhostCellWidth());
            if (state.enabled && !state.active[patch->getPatchNumber()])
            {
                continue;
            }

            // 1. Transpose of the physical boundary fill.
            if (bdry_op)
            {
                bdry_op->accumulateFromPhysicalBoundaryData(*patch, fill_time, d_gcw);
            }

            // 2. Copy the values that are summed into the interiors.
            patch->getPatchData(d_scratch_idx)->copy(*data);
        }

        // 3. Sum the values in ghost regions into the interiors of the patches that own them.
        level_info.schedule->fillData(fill_time, /*do_physical_boundary_fill*/ false);
        if (level_info.gap_schedule)
        {
            level_info.gap_schedule->communicate();
        }
    }
    IBTK_TIMER_STOP(t_accumulate_ghost_data);
}

/////////////////////////////// PRIVATE //////////////////////////////////////

void
SAMRAIGhostDataAccumulator::initializeLevel(const int ln)
{
    LevelInfo& level_info = *d_level_info[ln];
    if (level_info.initialized)
    {
        return;
    }
    Pointer<PatchLevel<NDIM>> level = level_info.level;

    check_patch_widths(level, d_gcw);

    int depth = 0;
    if (d_cc_data)
    {
        Pointer<CellDataFactory<NDIM, double>> factory =
            level->getPatchDescriptor()->getPatchDataFactory(d_scratch_idx);
        depth = factory->getDefaultDepth();
    }
    else
    {
        Pointer<SideDataFactory<NDIM, double>> factory =
            level->getPatchDescriptor()->getPatchDataFactory(d_scratch_idx);
        depth = factory->getDefaultDepth();
    }

    // The destination of the schedule is the copy of the data, which holds the
    // values in the ghost regions. The source is the copy as well, so that the
    // schedule depends only on the level; the data that receive the sums are
    // given by the state of the schedule.
    level_info.state = std::make_shared<ScheduleState>();
    level_info.algorithm = new RefineAlgorithm<NDIM>();
    level_info.algorithm->registerRefine(
        d_scratch_idx, d_scratch_idx, d_scratch_idx, Pointer<RefineOperator<NDIM>>(nullptr));
    Pointer<RefineTransactionFactory<NDIM>> factory = new GhostSumTransactionFactory(
        d_cc_data ? Centering::CELL : Centering::SIDE, depth, d_scratch_idx, level_info.state);
    level_info.schedule = level_info.algorithm->createSchedule("DEFAULT_FILL", level, nullptr, factory);

    // The schedule does not move the values at the sides on the outermost layer of the ghost side box of a patch p to a
    // patch q whose cells are all outside the ghost cell box of p. They are moved by a schedule of the same
    // transactions.
    if (!d_cc_data)
    {
        const std::vector<GapTransfer> transfers = compute_gap_transfers(level, d_gcw);
        if (!transfers.empty())
        {
            level_info.gap_schedule = new tbox::Schedule();
            for (const GapTransfer& transfer : transfers)
            {
                const SideOverlap<NDIM> overlap(transfer.boxes.data(), transfer.shift);
                level_info.gap_schedule->appendTransaction(Pointer<tbox::Transaction>(new GhostSumTransaction(
                    Centering::SIDE, depth, d_scratch_idx, level, overlap, transfer.p, transfer.q, level_info.state)));
            }
        }
    }
    level_info.initialized = true;
}

////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
