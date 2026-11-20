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

#include <ibtk/AppInitializer.h>
#include <ibtk/EdgeDataSynchronization.h>
#include <ibtk/FaceDataSynchronization.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/NodeDataSynchronization.h>
#include <ibtk/SideDataSynchronization.h>

#include <tbox/Database.h>
#include <tbox/PIO.h>
#include <tbox/Pointer.h>

#include <ArrayData.h>
#include <Box.h>
#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <CoarsenAlgorithm.h>
#include <CoarsenSchedule.h>
#include <EdgeData.h>
#include <EdgeGeometry.h>
#include <EdgeVariable.h>
#include <FaceData.h>
#include <FaceGeometry.h>
#include <FaceVariable.h>
#include <Index.h>
#include <IntVector.h>
#include <NodeData.h>
#include <NodeGeometry.h>
#include <NodeVariable.h>
#include <OuteredgeVariable.h>
#include <OuterfaceVariable.h>
#include <OuternodeVariable.h>
#include <OutersideVariable.h>
#include <Patch.h>
#include <PatchData.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <ProcessorMapping.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <Variable.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

// Report which copy of each patch border value an operation uses. The test
// counts values by the source cells that touch them
//   level border: one patch and cells outside the source level
//   two patches:  two patches and no others
//   junctions:    three or more of these
//   elsewhere:    one patch only, or none
// and by origin; the owner is the patch that contains the highest touching
// cell, comparing the last coordinate first.

namespace
{
using namespace SAMRAI;

using ArrayVector = std::vector<pdat::ArrayData<NDIM, double>*>;
// (array, index) or (patch, periodic shift).
using Key = std::array<int, NDIM + 1>;

// Values that are not marks; marks are nonnegative.
const double COARSE = -1.0;
const double NOT_SET = -2.0;
// Larger than any patch number and than any index on level 0 plus 2.
const double MARK_BASE = 64.0;
const int RATIO = 2;
const std::array<std::string, 2> PRIORITY_NAMES = { "coarse-priority", "fine-priority" };

// Variable 0 gives priority to the coarse values on a coarse-fine interface and
// variable 1 to the fine values.
struct Centering
{
    std::string name, refine_op, coarsen_op;
    int n_arrays;
    bool has_outer_coarsen_op;
    std::array<tbox::Pointer<hier::Variable<NDIM>>, 2> vars;
    std::array<int, 2> idxs, ghost_idxs;
    tbox::Pointer<hier::Variable<NDIM>> outer_var;
    int outer_idx;
    std::function<ArrayVector(hier::PatchData<NDIM>&)> arrays;
    // The index box of array n for data on the given box of cells.
    std::function<hier::Box<NDIM>(const hier::Box<NDIM>&, int)> array_box;
    std::function<void(int, tbox::Pointer<hier::PatchHierarchy<NDIM>>)> synchronize;
};

template <class DataType>
pdat::ArrayData<NDIM, double>&
get_array(DataType& data, const int axis)
{
    return data.getArrayData(axis);
}

pdat::ArrayData<NDIM, double>&
get_array(pdat::NodeData<NDIM, double>& data, const int /*n*/)
{
    return data.getArrayData();
}

// Register the data of a centering; this must precede making the hierarchy.
template <template <int, class> class VariableType,
          template <int, class>
          class OuterVariableType,
          template <int, class>
          class DataType,
          class SynchronizationType>
Centering
make_centering(const std::string& name,
               const std::string& refine_op,
               const std::string& coarsen_op,
               const int n_arrays,
               const bool has_outer_coarsen_op,
               std::function<hier::Box<NDIM>(const hier::Box<NDIM>&, int)> array_box)
{
    Centering centering;
    centering.name = name;
    centering.refine_op = refine_op;
    centering.coarsen_op = coarsen_op;
    centering.n_arrays = n_arrays;
    centering.has_outer_coarsen_op = has_outer_coarsen_op;
    auto* var_db = hier::VariableDatabase<NDIM>::getDatabase();
    const hier::IntVector<NDIM> no_ghosts(0), ghosts(1);
    for (int priority = 0; priority <= 1; ++priority)
    {
        tbox::Pointer<hier::Variable<NDIM>> var =
            new VariableType<NDIM, double>(name + "::" + PRIORITY_NAMES[priority], 1, priority == 1);
        centering.vars[priority] = var;
        centering.idxs[priority] = var_db->registerVariableAndContext(var, var_db->getContext("no ghosts"), no_ghosts);
        centering.ghost_idxs[priority] = var_db->registerVariableAndContext(var, var_db->getContext("ghosts"), ghosts);
    }
    centering.outer_var = new OuterVariableType<NDIM, double>(name + "::outer", 1);
    centering.outer_idx =
        var_db->registerVariableAndContext(centering.outer_var, var_db->getContext("no ghosts"), no_ghosts);
    centering.arrays = [n_arrays](hier::PatchData<NDIM>& data)
    {
        ArrayVector arrays;
        for (int n = 0; n < n_arrays; ++n)
        {
            arrays.push_back(&get_array(dynamic_cast<DataType<NDIM, double>&>(data), n));
        }
        return arrays;
    };
    centering.array_box = array_box;
    centering.synchronize = [coarsen_op](const int idx, tbox::Pointer<hier::PatchHierarchy<NDIM>> hierarchy)
    {
        SynchronizationType synch_op;
        synch_op.initializeOperatorState(
            typename SynchronizationType::SynchronizationTransactionComponent(idx, coarsen_op), hierarchy);
        synch_op.synchronizeData(0.0);
    };
    return centering;
}

Key
make_key(const int n, const hier::IntVector<NDIM>& i)
{
    Key key;
    key[0] = n;
    for (int d = 0; d < NDIM; ++d)
    {
        key[d + 1] = i(d);
    }
    return key;
}

// Identifies a value's array, patch, and position, rounded down to level 0 so
// that the fine values that a coarsen operator combines have one mark.
double
mark(const int patch, const Key& key, const hier::IntVector<NDIM>& ratio)
{
    double position = key[0];
    for (int d = 0; d < NDIM; ++d)
    {
        position = MARK_BASE * position + std::floor(static_cast<double>(key[d + 1]) / ratio(d)) + 2.0;
    }
    return MARK_BASE * position + patch;
}

// Includes ghost values.
void
fill_data(const Centering& centering, tbox::Pointer<hier::PatchLevel<NDIM>> level, const int idx, const double value)
{
    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        for (auto* array : centering.arrays(*level->getPatch(p())->getPatchData(idx)))
        {
            array->fillAll(value);
        }
    }
}

void
mark_data(const Centering& centering, tbox::Pointer<hier::PatchLevel<NDIM>> level, const int idx)
{
    fill_data(centering, level, idx, NOT_SET);
    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        tbox::Pointer<hier::Patch<NDIM>> patch = level->getPatch(p());
        const ArrayVector arrays = centering.arrays(*patch->getPatchData(idx));
        for (int n = 0; n < centering.n_arrays; ++n)
        {
            for (hier::Box<NDIM>::Iterator i(centering.array_box(patch->getBox(), n)); i; i++)
            {
                (*arrays[n])(i(), 0) = mark(p(), make_key(n, i()), level->getRatio());
            }
        }
    }
}

// The cells that touch a value, each as (patch, periodic shift), or -1 and
// zeros for a cell that is in no patch.
struct Touching
{
    std::set<Key> sources;
    // The mark of the copy that belongs to the highest cell that is in a patch.
    double owner_mark = std::numeric_limits<double>::quiet_NaN();
};

// Find the source-level cells that touch each value on a level.
std::map<Key, Touching>
find_touching_cells(const Centering& centering,
                    tbox::Pointer<hier::PatchLevel<NDIM>> level,
                    tbox::Pointer<hier::PatchLevel<NDIM>> src_level)
{
    hier::BoxArray<NDIM> src_boxes = src_level->getBoxes();
    src_boxes.coarsen(src_level->getRatio() / level->getRatio());
    const hier::Box<NDIM> domain = level->getPhysicalDomain()[0];
    const hier::IntVector<NDIM> periodic_shift = level->getGridGeometry()->getPeriodicShift(level->getRatio());
    std::map<Key, Touching> touching;
    // The iterator visits the cells in increasing order, comparing the last
    // coordinate first, so the last cell that touches a value is the highest.
    // A cell outside a periodic domain is the image of the cell that the
    // periodic shift takes to it.
    for (hier::Box<NDIM>::Iterator c(hier::Box<NDIM>::grow(domain, hier::IntVector<NDIM>(1))); c; c++)
    {
        hier::IntVector<NDIM> shift(0);
        for (int d = 0; d < NDIM; ++d)
        {
            shift(d) = c()(d) < domain.lower(d) ? -periodic_shift(d) : c()(d) > domain.upper(d) ? periodic_shift(d) : 0;
        }
        int patch = -1;
        for (int k = 0; k < src_boxes.size(); ++k)
        {
            if (src_boxes[k].contains(c() - shift))
            {
                patch = k;
            }
        }
        // The values of an image cell are the images of the original cell's values.
        const hier::Box<NDIM> cell(c(), c()), src_cell(c() - shift, c() - shift);
        for (int n = 0; n < centering.n_arrays; ++n)
        {
            hier::Box<NDIM>::Iterator src_i(centering.array_box(src_cell, n));
            for (hier::Box<NDIM>::Iterator i(centering.array_box(cell, n)); i; i++, src_i++)
            {
                Touching& value_touching = touching[make_key(n, i())];
                value_touching.sources.insert(make_key(patch, patch < 0 ? hier::IntVector<NDIM>(0) : shift));
                if (patch >= 0)
                {
                    value_touching.owner_mark = mark(patch, make_key(n, src_i()), level->getRatio());
                }
            }
        }
    }
    return touching;
}

// Print the nonzero counts of a level's values by touching source cells and by
// origin. Count values off the source patch borders only on request.
void
report(const std::string& operation,
       const Centering& centering,
       tbox::Pointer<hier::PatchLevel<NDIM>> level,
       const int idx,
       tbox::Pointer<hier::PatchLevel<NDIM>> src_level,
       const bool examine_elsewhere = false)
{
    const std::array<std::string, 4> group_names = { "level border", "two patches", "junctions", "elsewhere" };
    const std::array<std::string, 5> origin_names = { "owner", "own copy", "coarse", "not set", "other" };
    const std::map<Key, Touching> touching = find_touching_cells(centering, level, src_level);
    std::array<std::array<int, 5>, 4> counts{};
    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        tbox::Pointer<hier::Patch<NDIM>> patch = level->getPatch(p());
        const ArrayVector arrays = centering.arrays(*patch->getPatchData(idx));
        for (int n = 0; n < centering.n_arrays; ++n)
        {
            for (hier::Box<NDIM>::Iterator i(centering.array_box(patch->getBox(), n)); i; i++)
            {
                // Number groups and origins in the order of their names.
                // The source of a cell outside the source level sorts first,
                // and a value that matches no candidate has the last origin.
                const Key key = make_key(n, i());
                const Touching& value_touching = touching.at(key);
                const std::size_t n_sources = value_touching.sources.size();
                const bool on_level_border = value_touching.sources.begin()->front() < 0;
                const int group = n_sources == 1 ? 3 : n_sources > 2 ? 2 : on_level_border ? 0 : 1;
                const std::array<double, 4> candidates = {
                    value_touching.owner_mark, mark(p(), key, level->getRatio()), COARSE, NOT_SET
                };
                const std::ptrdiff_t origin =
                    std::find(candidates.begin(), candidates.end(), (*arrays[n])(i(), 0)) - candidates.begin();
                if (group != 3 || examine_elsewhere)
                {
                    ++counts[group][origin];
                }
            }
        }
    }

    tbox::plog << operation << ':';
    std::string separator = " ";
    for (std::size_t group = 0; group < counts.size(); ++group)
    {
        IBTK::IBTK_MPI::sumReduction(counts[group].data(), static_cast<int>(counts[group].size()));
        std::string origins;
        for (std::size_t origin = 0; origin < counts[group].size(); ++origin)
        {
            if (counts[group][origin] > 0)
            {
                origins += (origins.empty() ? "" : " / ") + origin_names[origin] + " (" +
                           std::to_string(counts[group][origin]) + ")";
            }
        }
        if (!origins.empty())
        {
            tbox::plog << separator << group_names[group] << " = " << origins;
            separator = "; ";
        }
    }
    tbox::plog << '\n';
}

// Assign patches to the processors in turn.
tbox::Pointer<hier::PatchHierarchy<NDIM>>
make_hierarchy(const std::vector<Centering>& centerings,
               tbox::Pointer<geom::CartesianGridGeometry<NDIM>> grid_geometry,
               const std::vector<hier::BoxArray<NDIM>>& boxes)
{
    tbox::Pointer<hier::PatchHierarchy<NDIM>> hierarchy =
        new hier::PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    for (int ln = 0; ln < static_cast<int>(boxes.size()); ++ln)
    {
        hier::ProcessorMapping mapping(boxes[ln].size());
        for (int k = 0; k < boxes[ln].size(); ++k)
        {
            mapping.setProcessorAssignment(k, k % IBTK::IBTK_MPI::getNodes());
        }
        hierarchy->makeNewPatchLevel(ln, hier::IntVector<NDIM>(ln == 0 ? 1 : RATIO), boxes[ln], mapping);
        tbox::Pointer<hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (const Centering& centering : centerings)
        {
            level->allocatePatchData(centering.outer_idx, 0.0);
            for (int priority = 0; priority <= 1; ++priority)
            {
                level->allocatePatchData(centering.idxs[priority], 0.0);
                level->allocatePatchData(centering.ghost_idxs[priority], 0.0);
            }
        }
    }
    return hierarchy;
}

// Fill level 1 and coarsen it to level 0 with each kind of schedule and the
// classes built on them. Level 1 is an L of three patches, numbered so that a
// rule based on patch numbers, or one that compares the first coordinate
// first, picks a different owner than the rule tested.
void
test_two_levels(const std::vector<Centering>& centerings,
                tbox::Pointer<geom::CartesianGridGeometry<NDIM>> grid_geometry)
{
    auto fine_box = [](const int i, const int j)
    {
        hier::Index<NDIM> lower(8);
        lower(0) = i;
        lower(1) = j;
        return hier::Box<NDIM>(lower, lower + 7);
    };
    std::vector<hier::BoxArray<NDIM>> boxes = { grid_geometry->getPhysicalDomain(), hier::BoxArray<NDIM>(3) };
    boxes[1][0] = fine_box(8, 16);
    boxes[1][1] = fine_box(8, 8);
    boxes[1][2] = fine_box(16, 8);
    tbox::Pointer<hier::PatchHierarchy<NDIM>> hierarchy = make_hierarchy(centerings, grid_geometry, boxes);
    tbox::Pointer<hier::PatchLevel<NDIM>> coarse_level = hierarchy->getPatchLevel(0);
    tbox::Pointer<hier::PatchLevel<NDIM>> fine_level = hierarchy->getPatchLevel(1);

    for (const Centering& centering : centerings)
    {
        auto prepare = [&](const int idx)
        {
            fill_data(centering, coarse_level, idx, COARSE);
            mark_data(centering, fine_level, idx);
        };

        // A fill takes the priority from its destination and a coarsening from
        // its source, so use every pair of priorities.
        tbox::plog << centering.name << " data\n";
        for (int src_priority = 0; src_priority <= 1; ++src_priority)
        {
            for (int dst_priority = 0; dst_priority <= 1; ++dst_priority)
            {
                const bool in_place = src_priority == dst_priority;
                const int dst_idx = centering.ghost_idxs[dst_priority];
                const int src_idx = in_place ? dst_idx : centering.idxs[src_priority];
                const std::string data =
                    PRIORITY_NAMES[src_priority] +
                    (in_place ? " data in place" : " into " + PRIORITY_NAMES[dst_priority] + " data");
                if (src_priority == 0 || dst_priority == 0)
                {
                    prepare(src_idx);
                    prepare(dst_idx);
                    xfer::RefineAlgorithm<NDIM> refine_alg;
                    refine_alg.registerRefine(
                        dst_idx,
                        src_idx,
                        dst_idx,
                        grid_geometry->lookupRefineOperator(centering.vars[0], centering.refine_op));
                    refine_alg.createSchedule(fine_level, 0, hierarchy)->fillData(0.0);
                    report("  refine schedule, " + data, centering, fine_level, dst_idx, fine_level);
                }

                prepare(src_idx);
                prepare(dst_idx);
                xfer::CoarsenAlgorithm<NDIM> coarsen_alg;
                coarsen_alg.registerCoarsen(
                    dst_idx, src_idx, grid_geometry->lookupCoarsenOperator(centering.vars[0], centering.coarsen_op));
                coarsen_alg.createSchedule(coarse_level, fine_level)->coarsenData();
                report("  coarsen schedule, " + data, centering, coarse_level, dst_idx, fine_level);
            }
        }

        const int idx = centering.ghost_idxs[1];
        prepare(idx);
        centering.synchronize(idx, hierarchy);
        report("  " + centering.name + "DataSynchronization, level 1", centering, fine_level, idx, fine_level);
        report("  " + centering.name + "DataSynchronization, level 0", centering, coarse_level, idx, fine_level);

        prepare(idx);
        using ITC = IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        IBTK::HierarchyGhostCellInterpolation ghost_fill_op;
        ghost_fill_op.initializeOperatorState(ITC(idx, "NONE", false, "NONE", "NONE", false), hierarchy);
        ghost_fill_op.fillData(0.0);
        report("  HierarchyGhostCellInterpolation without a refine operator", centering, fine_level, idx, fine_level);
    }
}

// Copy and coarsen from outer data when the source patches are one cell thick,
// so that the outer data of a patch hold both of its faces in that direction.
// Border values must come from the owner, and the operation must leave all
// other values alone.
void
test_outer_data(const std::vector<Centering>& centerings,
                tbox::Pointer<geom::CartesianGridGeometry<NDIM>> grid_geometry)
{
    std::vector<hier::BoxArray<NDIM>> boxes = { hier::BoxArray<NDIM>(16), hier::BoxArray<NDIM>(8) };
    for (int k = 0; k < 16; ++k)
    {
        hier::Index<NDIM> lower(0), upper(15);
        lower(0) = upper(0) = k;
        boxes[0][(5 * k + 3) % 16] = hier::Box<NDIM>(lower, upper);
    }
    for (int k = 0; k < 8; ++k)
    {
        hier::Index<NDIM> lower(8), upper(23);
        lower(0) = 8 + 2 * k;
        upper(0) = lower(0) + 1;
        boxes[1][(5 * k + 3) % 8] = hier::Box<NDIM>(lower, upper);
    }
    tbox::Pointer<hier::PatchHierarchy<NDIM>> hierarchy = make_hierarchy(centerings, grid_geometry, boxes);
    tbox::Pointer<hier::PatchLevel<NDIM>> dst_level = hierarchy->getPatchLevel(0);

    for (const Centering& centering : centerings)
    {
        // Copy when the source is level 0 and coarsen when the source is level 1.
        // The outer data get their marks from a copy of ordinary data.
        const int dst_idx = centering.idxs[0], tmp_idx = centering.idxs[1], outer_idx = centering.outer_idx;
        for (int src_ln = 0; src_ln <= (centering.has_outer_coarsen_op ? 1 : 0); ++src_ln)
        {
            tbox::Pointer<hier::PatchLevel<NDIM>> src_level = hierarchy->getPatchLevel(src_ln);
            fill_data(centering, dst_level, dst_idx, NOT_SET);
            mark_data(centering, src_level, tmp_idx);
            for (hier::PatchLevel<NDIM>::Iterator p(src_level); p; p++)
            {
                tbox::Pointer<hier::Patch<NDIM>> patch = src_level->getPatch(p());
                patch->getPatchData(outer_idx)->copy(*patch->getPatchData(tmp_idx));
            }

            if (src_ln == 0)
            {
                xfer::RefineAlgorithm<NDIM> copy_alg;
                copy_alg.registerRefine(dst_idx, outer_idx, dst_idx, nullptr);
                copy_alg.createSchedule(dst_level)->fillData(0.0);
            }
            else
            {
                xfer::CoarsenAlgorithm<NDIM> coarsen_alg;
                coarsen_alg.registerCoarsen(
                    dst_idx,
                    outer_idx,
                    grid_geometry->lookupCoarsenOperator(centering.outer_var, centering.coarsen_op));
                coarsen_alg.createSchedule(dst_level, src_level)->coarsenData();
            }
            const std::string operation = src_ln == 0 ? " data, copy" : " data, coarsening";
            report(centering.name + operation + " from outer data", centering, dst_level, dst_idx, src_level, true);
        }
    }
}

// Copy and synchronize on a single patch covering a periodic domain, which is
// its own neighbor: upper-border values must come from the lowest-index copy.
void
test_periodic_domain(const std::vector<Centering>& centerings,
                     tbox::Pointer<geom::CartesianGridGeometry<NDIM>> grid_geometry)
{
    tbox::Pointer<hier::PatchHierarchy<NDIM>> hierarchy =
        make_hierarchy(centerings, grid_geometry, { grid_geometry->getPhysicalDomain() });
    tbox::Pointer<hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
    for (const Centering& centering : centerings)
    {
        const int idx = centering.ghost_idxs[1];
        mark_data(centering, level, idx);
        xfer::RefineAlgorithm<NDIM> copy_alg;
        copy_alg.registerRefine(idx, idx, idx, nullptr);
        copy_alg.createSchedule(level)->fillData(0.0);
        report(centering.name + " data, refine schedule without a refine operator", centering, level, idx, level);

        mark_data(centering, level, idx);
        centering.synchronize(idx, hierarchy);
        report(centering.name + "DataSynchronization", centering, level, idx, level);
    }
}
} // namespace

int
main(int argc, char** argv)
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    tbox::Pointer<IBTK::AppInitializer> app_initializer = new IBTK::AppInitializer(argc, argv);
    tbox::Pointer<geom::CartesianGridGeometry<NDIM>> grid_geometry = new geom::CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));

    using pdat::EdgeGeometry;
    using pdat::FaceGeometry;
    using pdat::NodeGeometry;
    using pdat::SideGeometry;
    const std::vector<Centering> centerings = {
        make_centering<pdat::SideVariable, pdat::OutersideVariable, pdat::SideData, IBTK::SideDataSynchronization>(
            "Side", "CONSTANT_REFINE", "CONSERVATIVE_COARSEN", NDIM, true, SideGeometry<NDIM>::toSideBox),
        make_centering<pdat::FaceVariable, pdat::OuterfaceVariable, pdat::FaceData, IBTK::FaceDataSynchronization>(
            "Face", "CONSTANT_REFINE", "CONSERVATIVE_COARSEN", NDIM, true, FaceGeometry<NDIM>::toFaceBox),
        make_centering<pdat::EdgeVariable, pdat::OuteredgeVariable, pdat::EdgeData, IBTK::EdgeDataSynchronization>(
            "Edge", "CONSTANT_REFINE", "CONSERVATIVE_COARSEN", NDIM, false, EdgeGeometry<NDIM>::toEdgeBox),
        make_centering<pdat::NodeVariable, pdat::OuternodeVariable, pdat::NodeData, IBTK::NodeDataSynchronization>(
            "Node",
            "LINEAR_REFINE",
            "CONSTANT_COARSEN",
            1,
            true,
            [](const hier::Box<NDIM>& box, const int /*n*/) { return NodeGeometry<NDIM>::toNodeBox(box); })
    };

    tbox::Pointer<tbox::Database> input_db = app_initializer->getInputDatabase();
    if (input_db->getBoolWithDefault("TEST_OUTER_DATA", false))
    {
        test_outer_data(centerings, grid_geometry);
    }
    else if (input_db->getBoolWithDefault("TEST_PERIODIC_DOMAIN", false))
    {
        test_periodic_domain(centerings, grid_geometry);
    }
    else
    {
        test_two_levels(centerings, grid_geometry);
    }
}
