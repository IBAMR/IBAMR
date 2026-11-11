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

// Centerings of the data that can be accumulated. The values of a patch data are a list of arrays, one for each
// component.
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

// Values that a schedule and its transactions share, set before each use.
struct ScheduleState
{
    // Patch data index of the data to which the sums are added.
    int target_idx = invalid_index;

    // Flags indicating which patches of the level may contain nonzero values.
    bool restrict_to_active_patches = false;
    std::vector<char> is_active_patch; // indexed by patch number
};

/*!
 * Transaction that sums the values in the ghost region of patch p into the patch data with index
 * ScheduleState::target_idx of patch q. This transaction reverses the source and destination processors of a refine
 * transaction.
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
        return d_state->restrict_to_active_patches && !d_state->is_active_patch[d_p];
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

    // The patch data indices come from the constructor and the state, not from the refine items.
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

// Sides of the ghost side box of patch p that patch q owns. The refine schedule of the level does not transfer their
// values.
struct GapTransfer
{
    int p = -1;
    int q = -1;
    IntVector<NDIM> shift = IntVector<NDIM>(0); // q is shifted by this to be near p
    SideBoxLists boxes;                         // in the index space of p
};

/*!
 * Return the sides of the ghost side box of patch p that lie in the side box of patch q or of a periodic image of q,
 * for all pairs of patches. The refine schedule does not transfer these sides. The result is sorted by patch numbers
 * and shift.
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

    // The shifts that map a patch to its periodic images, and a test of whether the schedule uses a shift for a patch.
    const IntVector<NDIM> periodic_shift = level->getGridGeometry()->getPeriodicShift(level->getRatio());
    std::vector<IntVector<NDIM>> all_shifts(1, no_shift);
    for (int d = 0; d < NDIM; ++d)
    {
        if (periodic_shift(d) == 0)
        {
            continue;
        }
        const std::size_t n = all_shifts.size();
        for (std::size_t k = 0; k < n; ++k)
        {
            for (const int sign : { -1, 1 })
            {
                IntVector<NDIM> shift = all_shifts[k];
                shift(d) += sign * periodic_shift(d);
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
        const Box<NDIM> shifted_q_box = Box<NDIM>::shift(patch_boxes[q], shift);
        if ((wide_ghost_box * shifted_q_box).empty())
        {
            return;
        }

        Box<NDIM> scheduled_box;
        if (is_schedule_shift(q, shift))
        {
            scheduled_box = ghost_box * shifted_q_box;
            if (scheduled_box.empty() && zero_gcw)
            {
                scheduled_box = wide_ghost_box * shifted_q_box;
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
            const Box<NDIM> shared_side_box = ghost_side_box * SideGeometry<NDIM>::toSideBox(shifted_q_box, axis);
            if (shared_side_box.empty())
            {
                continue;
            }
            BoxList<NDIM> missing_sides(shared_side_box);
            if (!scheduled_box.empty())
            {
                missing_sides.removeIntersections(ghost_side_box * SideGeometry<NDIM>::toSideBox(scheduled_box, axis));
            }
            for (BoxList<NDIM>::Iterator b(missing_sides); b; b++)
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

    // The box tree contains the patches together with their periodic images. The overlap indices are sorted and unique.
    tbox::Array<int> candidates;
    auto for_each_candidate = [&](const Box<NDIM>& region, const auto& function)
    {
        box_tree->findOverlapIndices(candidates, region);
        for (int k = 0; k < candidates.size(); ++k)
        {
            function(candidates[k]);
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
 * Emit an error if patches are narrower than the ghost cell width in a direction in which the level has adjacent
 * patches or is periodic. Patches are adjacent if their boxes grown by the ghost cell width in that direction
 * intersect.
 */
void
check_patch_widths(Pointer<PatchLevel<NDIM>> level, const IntVector<NDIM>& gcw)
{
    const BoxArray<NDIM>& patch_boxes = level->getBoxes();
    const IntVector<NDIM> periodic_shift = level->getGridGeometry()->getPeriodicShift(level->getRatio());
    Pointer<BoxTree<NDIM>> box_tree = level->getBoxTree();
    tbox::Array<int> candidates;
    for (int d = 0; d < NDIM; ++d)
    {
        bool check = (periodic_shift(d) != 0);
        IntVector<NDIM> growth(0);
        growth(d) = gcw(d);
        for (int p = 0; p < patch_boxes.getNumberOfBoxes() && !check; ++p)
        {
            box_tree->findOverlapIndices(candidates, Box<NDIM>::grow(patch_boxes[p], growth));
            for (int k = 0; k < candidates.size() && !check; ++k)
            {
                check = (candidates[k] != p);
            }
        }
        if (!check)
        {
            continue;
        }
        for (int p = 0; p < patch_boxes.getNumberOfBoxes(); ++p)
        {
            if (patch_boxes[p].numberCells(d) < gcw(d))
            {
                TBOX_ERROR("SAMRAIGhostDataAccumulator::check_patch_widths():\n"
                           << "  patch " << p << " of level " << level->getLevelNumber() << " (box " << patch_boxes[p]
                           << ") is " << patch_boxes[p].numberCells(d) << " cells wide along axis " << d
                           << ", less than the ghost cell width " << gcw(d) << "\n"
                           << "  patches must be at least as wide as the ghost cell width along every axis in which "
                              "the level has adjacent patches or is periodic\n");
            }
        }
    }
}

/*!
 * Emit an error if a patch lies within the ghost cell width of a periodic boundary without touching it, because SAMRAI
 * lists periodic images only for patches that touch the boundary.
 */
void
check_periodic_images(Pointer<PatchLevel<NDIM>> level, const IntVector<NDIM>& gcw)
{
    const BoxArray<NDIM>& patch_boxes = level->getBoxes();
    const IntVector<NDIM> periodic_shift = level->getGridGeometry()->getPeriodicShift(level->getRatio());
    const BoxArray<NDIM> domain = level->getPhysicalDomain();

    // The periodic extent is that of the bounding box of the domain.
    Box<NDIM> domain_bounds = domain[0];
    for (int n = 1; n < domain.size(); ++n)
    {
        domain_bounds += domain[n];
    }
    for (int d = 0; d < NDIM; ++d)
    {
        if (periodic_shift(d) == 0)
        {
            continue;
        }
        const int lower = domain_bounds.lower(d);
        const int upper = domain_bounds.upper(d);
        for (int p = 0; p < patch_boxes.getNumberOfBoxes(); ++p)
        {
            const int distance = std::min(patch_boxes[p].lower(d) == lower ? gcw(d) : patch_boxes[p].lower(d) - lower,
                                          patch_boxes[p].upper(d) == upper ? gcw(d) : upper - patch_boxes[p].upper(d));
            if (distance < gcw(d))
            {
                TBOX_ERROR("SAMRAIGhostDataAccumulator::check_periodic_images():\n"
                           << "  patch " << p << " of level " << level->getLevelNumber() << " (box " << patch_boxes[p]
                           << ") is " << distance << " cells from the periodic boundary along axis " << d
                           << ", less than the ghost cell width " << gcw(d) << "\n"
                           << "  patches must touch a periodic boundary or be at least as far from it as the ghost "
                              "cell width\n");
            }
        }
    }
}

/*!
 * Emit an error if a patch lies within the ghost cell width of a physical boundary without touching it, because SAMRAI
 * gives boundary boxes only to patches that touch the boundary.
 */
void
check_physical_boundaries(Pointer<PatchLevel<NDIM>> level, const IntVector<NDIM>& gcw)
{
    const BoxArray<NDIM>& patch_boxes = level->getBoxes();
    const IntVector<NDIM> periodic_shift = level->getGridGeometry()->getPeriodicShift(level->getRatio());
    const BoxList<NDIM> domain(level->getPhysicalDomain());

    // Whether the cells of the box are in the physical domain, which may have several boxes.
    const auto in_domain = [&domain](const Box<NDIM>& box)
    {
        BoxList<NDIM> outside(box);
        outside.removeIntersections(domain);
        return outside.isEmpty();
    };

    for (int d = 0; d < NDIM; ++d)
    {
        if (periodic_shift(d) != 0 || gcw(d) < 2)
        {
            continue;
        }
        for (int p = 0; p < patch_boxes.getNumberOfBoxes(); ++p)
        {
            // The ghost region of a side is the slab of width gcw next to the patch. A patch touches the boundary
            // on a side if the slab of width one extends outside the domain.
            for (const int side : { 0, 1 })
            {
                const auto slab = [&](const int width)
                {
                    Box<NDIM> box = patch_boxes[p];
                    if (side == 0)
                    {
                        box.upper(d) = box.lower(d) - 1;
                        box.lower(d) -= width;
                    }
                    else
                    {
                        box.lower(d) = box.upper(d) + 1;
                        box.upper(d) += width;
                    }
                    return box;
                };
                if (!in_domain(slab(1)) || in_domain(slab(gcw(d))))
                {
                    continue;
                }
                int distance = 1;
                while (in_domain(slab(distance + 1)))
                {
                    ++distance;
                }
                TBOX_ERROR("SAMRAIGhostDataAccumulator::check_physical_boundaries():\n"
                           << "  patch " << p << " of level " << level->getLevelNumber() << " (box " << patch_boxes[p]
                           << ") is " << distance << " cells from the physical boundary along axis " << d
                           << ", less than the ghost cell width " << gcw(d) << "\n"
                           << "  patches must touch a physical boundary or be at least as far from it as the ghost "
                              "cell width\n");
            }
        }
    }
}
} // namespace

/////////////////////////////// NESTED TYPES /////////////////////////////////

struct SAMRAIGhostDataAccumulator::LevelSchedules
{
    // The level on which the scratch patch data are allocated.
    Pointer<PatchLevel<NDIM>> level;

    // Set when the level is first used.
    bool initialized = false;

    // The schedule of the level, built when the level is first used.
    Pointer<RefineSchedule<NDIM>> schedule;
    std::shared_ptr<ScheduleState> state;

    // For side-centered data, the schedule of the values that the schedule above omits. Null if this process neither
    // sends nor receives such values.
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
                   << "  unsupported variable type\n"
                   << "  use a cell-centered or side-centered variable with double precision values\n");
    }
    d_cc_data = cc_var;

    // Register a patch data index, private to this object, for a copy of the values that are summed.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("SAMRAIGhostDataAccumulator");
    const int shared_idx = var_db->registerVariableAndContext(var, context, d_gcw);
    d_scratch_idx = var_db->registerClonedPatchDataIndex(var, shared_idx);
    var_db->removePatchDataIndex(shared_idx);

    d_level_schedules.resize(d_finest_ln + 1);
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_scratch_idx);
        d_level_schedules[ln] = std::make_unique<LevelSchedules>();
        d_level_schedules[ln]->level = level;
    }
    IBTK_TIMER_STOP(t_constructor);
}

SAMRAIGhostDataAccumulator::~SAMRAIGhostDataAccumulator()
{
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        d_level_schedules[ln]->level->deallocatePatchData(d_scratch_idx);
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

    // The caller destroys the object when a level is replaced. Check the levels before the patch numbers of the caller
    // are used.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        if (Pointer<PatchLevel<NDIM>>(d_hierarchy->getPatchLevel(ln)) != d_level_schedules[ln]->level)
        {
            TBOX_ERROR("SAMRAIGhostDataAccumulator::accumulateGhostData():\n"
                       << "  level " << ln
                       << " of the patch hierarchy is not the level for which the accumulator was built\n"
                       << "  construct a new accumulator after the levels are regridded or redistributed\n");
        }
    }

    // Determine the patches that may contain nonzero values, on all processes, so that the sending and the receiving
    // process skip the same transactions.
    std::vector<int> active_flags;
    std::vector<int> level_offsets(num_levels + 1, 0);
    if (active_patch_nums)
    {
        TBOX_ASSERT(static_cast<int>(active_patch_nums->size()) == num_levels);
        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            level_offsets[ln - d_coarsest_ln + 1] =
                level_offsets[ln - d_coarsest_ln] + d_level_schedules[ln]->level->getNumberOfPatches();
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
        LevelSchedules& level_schedules = *d_level_schedules[ln];
        Pointer<PatchLevel<NDIM>> level = level_schedules.level;

        // A level without active patches has nothing to accumulate.
        if (active_patch_nums && std::none_of(active_flags.begin() + level_offsets[l],
                                              active_flags.begin() + level_offsets[l + 1],
                                              [](const int flag) { return flag > 0; }))
        {
            continue;
        }
        initializeLevel(ln);

        ScheduleState& state = *level_schedules.state;
        state.target_idx = idx;
        state.restrict_to_active_patches = static_cast<bool>(active_patch_nums);
        if (active_patch_nums)
        {
            state.is_active_patch.assign(level->getNumberOfPatches(), 0);
            for (int n = 0; n < level->getNumberOfPatches(); ++n)
            {
                state.is_active_patch[n] = active_flags[level_offsets[l] + n] > 0;
            }
        }

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<PatchData<NDIM>> data = patch->getPatchData(idx);
            TBOX_ASSERT(d_gcw == data->getGhostCellWidth());
            if (state.restrict_to_active_patches && !state.is_active_patch[patch->getPatchNumber()])
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
        level_schedules.schedule->fillData(fill_time, /*do_physical_boundary_fill*/ false);
        if (level_schedules.gap_schedule)
        {
            level_schedules.gap_schedule->communicate();
        }
    }
    IBTK_TIMER_STOP(t_accumulate_ghost_data);
}

/////////////////////////////// PRIVATE //////////////////////////////////////

void
SAMRAIGhostDataAccumulator::initializeLevel(const int ln)
{
    LevelSchedules& level_schedules = *d_level_schedules[ln];
    if (level_schedules.initialized)
    {
        return;
    }
    Pointer<PatchLevel<NDIM>> level = level_schedules.level;

    check_patch_widths(level, d_gcw);
    check_physical_boundaries(level, d_gcw);
    if (d_cc_data)
    {
        check_periodic_images(level, d_gcw);
    }

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

    // Because the schedule uses the scratch copy as source and destination, it depends only on the level.
    // The state contains the patch data index of the data to which the sums are added.
    level_schedules.state = std::make_shared<ScheduleState>();
    Pointer<RefineAlgorithm<NDIM>> algorithm = new RefineAlgorithm<NDIM>();
    algorithm->registerRefine(d_scratch_idx, d_scratch_idx, d_scratch_idx, Pointer<RefineOperator<NDIM>>(nullptr));
    Pointer<RefineTransactionFactory<NDIM>> factory = new GhostSumTransactionFactory(
        d_cc_data ? Centering::CELL : Centering::SIDE, depth, d_scratch_idx, level_schedules.state);
    level_schedules.schedule = algorithm->createSchedule("DEFAULT_FILL", level, nullptr, factory);

    // The refine schedule omits the values on the outer edge of the ghost side box of patch p if patch q has no cell in
    // the ghost cell box of p. A second schedule of the same transactions sends them.
    if (!d_cc_data)
    {
        const std::vector<GapTransfer> transfers = compute_gap_transfers(level, d_gcw);
        if (!transfers.empty())
        {
            level_schedules.gap_schedule = new tbox::Schedule();
            for (const GapTransfer& transfer : transfers)
            {
                const SideOverlap<NDIM> overlap(transfer.boxes.data(), transfer.shift);
                level_schedules.gap_schedule->appendTransaction(
                    Pointer<tbox::Transaction>(new GhostSumTransaction(Centering::SIDE,
                                                                       depth,
                                                                       d_scratch_idx,
                                                                       level,
                                                                       overlap,
                                                                       transfer.p,
                                                                       transfer.q,
                                                                       level_schedules.state)));
            }
        }
    }
    level_schedules.initialized = true;
}

////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
