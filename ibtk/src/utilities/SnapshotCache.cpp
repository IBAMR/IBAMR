// ---------------------------------------------------------------------
//
// Copyright (c) 2022 - 2026 by the IBAMR developers
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
#include <ibtk/SnapshotCache.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/RestartManager.h>

#include <CartesianGridGeometry.h>
#include <HierarchyDataOpsManager.h>
#include <RefineAlgorithm.h>
#include <RefineOperator.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <iterator>
#include <limits>

#include <ibtk/app_namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
SnapshotCache::SnapshotCache(std::string object_name,
                             Pointer<Variable<NDIM>> var,
                             Pointer<Database> input_db,
                             Pointer<GridGeometry<NDIM>> grid_geom,
                             bool register_for_restart)
    : d_object_name(std::move(object_name)), d_snapshot_var(var)
{
    TBOX_ASSERT(grid_geom);

    if (input_db)
        if (input_db->keyExists("gcw")) input_db->getIntegerArray("gcw", &d_gcw[0], NDIM);

    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_ctx = var_db->getContext(d_object_name + "::context");
    d_snapshot_idx = var_db->registerVariableAndContext(d_snapshot_var, d_ctx, d_gcw);

    auto restart_manager = RestartManager::getManager();
    if (register_for_restart)
    {
        restart_manager->registerRestartItem(d_object_name, this);
        var_db->registerPatchDataForRestart(d_snapshot_idx);
        d_registered_for_restart = true;
    }

    bool from_restart = restart_manager->isFromRestart();
    if (from_restart) getFromRestart(grid_geom);
}

SnapshotCache::~SnapshotCache()
{
    clearSnapshots();
}

void
SnapshotCache::clearSnapshots()
{
    // Deallocate all patch data.
    for (auto& snapshot : d_snapshots)
    {
        Pointer<PatchHierarchy<NDIM>> hierarchy = snapshot.second;
        for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            if (level->checkAllocated(d_snapshot_idx)) level->deallocatePatchData(d_snapshot_idx);
        }
    }

    // Clear snapshots
    d_schedules.clear();
    d_current_levels = detail::HierarchyLevelsIdentity();
    d_snapshots.clear();
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    var_db->removePatchDataIndex(d_snapshot_idx);
    d_snapshot_var = nullptr;
    d_snapshot_idx = IBTK::invalid_index;
}

SnapshotCache::Schedules&
SnapshotCache::resetSchedules(const PatchHierarchy<NDIM>& current_hierarchy,
                              const Pointer<PatchHierarchy<NDIM>>& snapshot_hierarchy)
{
    if (d_current_levels.changed(current_hierarchy, 0, current_hierarchy.getFinestLevelNumber()))
    {
        d_schedules.clear();
    }
    auto it = std::find_if(d_schedules.begin(),
                           d_schedules.end(),
                           [&snapshot_hierarchy](const Schedules& schedules) -> bool
                           { return schedules.snapshot_hierarchy.getPointer() == snapshot_hierarchy.getPointer(); });
    if (it == d_schedules.end())
    {
        d_schedules.emplace_back();
        it = std::prev(d_schedules.end());
        it->snapshot_hierarchy = snapshot_hierarchy;
    }
    if (it->snapshot_levels.changed(*snapshot_hierarchy, 0, snapshot_hierarchy->getFinestLevelNumber()))
    {
        it->update_idx = IBTK::invalid_index;
        it->update.clear();
        it->fill.clear();
    }
    return *it;
}

SnapshotCache::value_type
SnapshotCache::getSnapshot(double time, double tol)
{
    auto it =
        std::find_if(d_snapshots.begin(),
                     d_snapshots.end(),
                     [time, tol](const value_type& t) -> bool { return IBTK::abs_equal_eps(t.first, time, tol); });
    if (it == d_snapshots.end())
        return std::make_pair(std::numeric_limits<double>::quiet_NaN(), nullptr);
    else
        return *it;
}

void
SnapshotCache::storeSnapshot(const int u_idx, const double time, Pointer<PatchHierarchy<NDIM>> hierarchy)
{
    // Create a snapshot of the patch hierarchy.
    Pointer<PatchHierarchy<NDIM>> snapshot_hierarchy = hierarchy->makeRefinedPatchHierarchy(
        d_object_name + "::SnapshotHierarchy_" + std::to_string(d_snapshots.size()),
        1 /*ratio*/,
        false /*register_for_restart*/);
    // The schedules that copy the data also serve updateSnapshot().
    Schedules& schedules = resetSchedules(*hierarchy, snapshot_hierarchy);
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    schedules.update_idx = u_idx;
    schedules.update_factory = var_db->getPatchDescriptor()->getPatchDataFactory(u_idx);
    // Allocate patch data and copy data to new patch hierarchy.
    int coarsest_ln = 0;
    int finest_ln = snapshot_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> snapshot_level = snapshot_hierarchy->getPatchLevel(ln);
        if (!snapshot_level->checkAllocated(d_snapshot_idx)) snapshot_level->allocatePatchData(d_snapshot_idx, time);
        // We've allocated data, now copy it
        // Note these levels should cover the same space, so we shouldn't need to have a refine operator.
        Pointer<PatchLevel<NDIM>> old_level = hierarchy->getPatchLevel(ln);
        // First fill in from the snapshot hierarchy
        Pointer<RefineAlgorithm<NDIM>> refine_alg = new RefineAlgorithm<NDIM>();
        Pointer<RefineOperator<NDIM>> refine_op = nullptr;
        refine_alg->registerRefine(d_snapshot_idx, u_idx, u_idx, refine_op);
        schedules.update.push_back(refine_alg->createSchedule(snapshot_level, old_level));
        schedules.update.back()->fillData(time);
    }

    // Store the index, time, hierarchy, and variable.
    // First determine where to store this index
    auto it = std::find_if(d_snapshots.begin(),
                           d_snapshots.end(),
                           [time](const value_type& snapshot) -> bool { return time < snapshot.first; });
    d_snapshots.insert(it, std::make_pair(time, snapshot_hierarchy));
}

void
SnapshotCache::updateSnapshot(const int u_idx,
                              const double time,
                              Pointer<PatchHierarchy<NDIM>> current_hierarchy,
                              const double tol)
{
    // The time must match a stored snapshot.
    const std::pair<double, Pointer<PatchHierarchy<NDIM>>>& snapshot = getSnapshot(time, tol);
    if (!snapshot.second || !IBTK::abs_equal_eps(snapshot.first, time, tol))
        TBOX_ERROR("Snapshot at time: " << time << " with tolerance " << tol << " does not exist!\n");

    const Pointer<PatchHierarchy<NDIM>>& snapshot_hierarchy = snapshot.second;
    int coarsest_ln = 0;
    int finest_ln = current_hierarchy->getFinestLevelNumber();
#ifndef NDEBUG
    TBOX_ASSERT(snapshot_hierarchy->getFinestLevelNumber() == finest_ln);
#endif
    Schedules& schedules = resetSchedules(*current_hierarchy, snapshot_hierarchy);
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<PatchDataFactory<NDIM>> u_factory = var_db->getPatchDescriptor()->getPatchDataFactory(u_idx);
    if (schedules.update_idx != u_idx || schedules.update_factory != u_factory)
    {
        schedules.update.clear();
    }
    if (schedules.update.empty())
    {
        // A refined patch hierarchy with refinement ratio 1 cannot copy the snapshot_hierarchy's name, so the
        // snapshot_hierarchy's levels are removed and replaced with copies of the current_hierarchy's levels.
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> current_level = current_hierarchy->getPatchLevel(ln);
            snapshot_hierarchy->removePatchLevel(ln);
            snapshot_hierarchy->makeNewPatchLevel(
                ln, current_level->getRatio(), current_level->getBoxes(), current_level->getProcessorMapping());
            Pointer<PatchLevel<NDIM>> new_level = snapshot_hierarchy->getPatchLevel(ln);

            // Allocate data on the new level.
            new_level->allocatePatchData(d_snapshot_idx, time);
            Pointer<RefineAlgorithm<NDIM>> refine_alg = new RefineAlgorithm<NDIM>();
            Pointer<RefineOperator<NDIM>> refine_op = nullptr;
            refine_alg->registerRefine(d_snapshot_idx, u_idx, u_idx, refine_op);
            schedules.update.push_back(refine_alg->createSchedule(new_level, current_level));
        }
        schedules.update_idx = u_idx;
        schedules.update_factory = u_factory;
        schedules.snapshot_levels.changed(*snapshot_hierarchy, coarsest_ln, finest_ln);
        schedules.fill.clear();
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        snapshot_hierarchy->getPatchLevel(ln)->setTime(time, d_snapshot_idx);
        schedules.update[ln]->fillData(time);
    }
}

void
SnapshotCache::fillSnapshotOnHierarchy(const int u_idx,
                                       const double time,
                                       Pointer<PatchHierarchy<NDIM>> current_hierarchy,
                                       const std::string& snapshot_refine_type,
                                       const double tol)
{
    // The time must match a stored snapshot.
    const std::pair<double, Pointer<PatchHierarchy<NDIM>>>& snapshot = getSnapshot(time, tol);
    if (!snapshot.second || !IBTK::abs_equal_eps(snapshot.first, time, tol))
        TBOX_ERROR("Snapshot at time: " << time << " with tolerance " << tol << " does not exist!\n");

    // The fill needs the snapshot index on the snapshot hierarchy and a scratch index on both hierarchies; the snapshot
    // index serves as the scratch index.
    int snapshot_idx = d_snapshot_idx;
    int scr_idx = snapshot_idx;

    const Pointer<PatchHierarchy<NDIM>>& snapshot_hierarchy = snapshot.second;
    // The current hierarchy and the snapshot hierarchy must have the same number of levels.
    TBOX_ASSERT(current_hierarchy->getFinestLevelNumber() == snapshot_hierarchy->getFinestLevelNumber());
    // Transfer data from the snapshot to the current hierarchy.
    int coarsest_ln = 0;
    int finest_ln = snapshot_hierarchy->getFinestLevelNumber();
    Schedules& schedules = resetSchedules(*current_hierarchy, snapshot_hierarchy);
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<PatchDataFactory<NDIM>> u_factory = var_db->getPatchDescriptor()->getPatchDataFactory(u_idx);
    auto fill = std::find_if(schedules.fill.begin(),
                             schedules.fill.end(),
                             [&](const Schedules::Fill& f) -> bool
                             { return f.fill_idx == u_idx && f.refine_type == snapshot_refine_type; });
    if (fill == schedules.fill.end())
    {
        schedules.fill.emplace_back();
        fill = std::prev(schedules.fill.end());
        fill->fill_idx = u_idx;
        fill->refine_type = snapshot_refine_type;
    }
    if (fill->fill_factory != u_factory || fill->schedules.empty())
    {
        fill->fill_factory = u_factory;
        fill->schedules.clear();
        Pointer<CartesianGridGeometry<NDIM>> grid_geom = current_hierarchy->getGridGeometry();
        Pointer<RefineOperator<NDIM>> refine_op = grid_geom->lookupRefineOperator(d_snapshot_var, snapshot_refine_type);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<RefineAlgorithm<NDIM>> refine_alg = new RefineAlgorithm<NDIM>();
            refine_alg->registerRefine(u_idx, snapshot_idx, scr_idx, refine_op);
            fill->schedules.push_back(refine_alg->createSchedule(current_hierarchy->getPatchLevel(ln),
                                                                 snapshot_hierarchy->getPatchLevel(ln),
                                                                 ln - 1,
                                                                 current_hierarchy));
        }
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> cur_level = current_hierarchy->getPatchLevel(ln);
        Pointer<PatchLevel<NDIM>> snp_level = snapshot_hierarchy->getPatchLevel(ln);

        // RefineSchedule requires source and destination data at the same time, so set the snapshot data to the time of
        // the current data, from the ranks that own patches on this level.
        double local_time = std::numeric_limits<double>::lowest();
        const PatchLevel<NDIM>::Iterator first_patch(cur_level);
        if (first_patch)
        {
            local_time = cur_level->getPatch(*first_patch)->getPatchData(u_idx)->getTime();
        }
        const double allocated_time = IBTK_MPI::maxReduction(local_time);
        // The scratch data on the current hierarchy use the snapshot index.
        cur_level->allocatePatchData(scr_idx, allocated_time);
        snp_level->setTime(allocated_time, snapshot_idx);

        fill->schedules[ln]->fillData(allocated_time);
    }

    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = current_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(scr_idx);
    }
}

void
SnapshotCache::putToDatabase(Pointer<Database> db)
{
    // Write data to the database.
    db->putInteger("num_snapshots_stored", d_snapshots.size());
    // Loop through snapshots and store the hierarchies
    int i = 0;
    for (const auto& snapshot : d_snapshots)
    {
        const double time = snapshot.first;
        db->putDouble("time_" + std::to_string(i), time);
        Pointer<Database> snapshot_db = db->putDatabase("snapshot_hierarchy_" + std::to_string(i));
        snapshot.second->putToDatabase(snapshot_db);
        ++i;
    }
}

void
SnapshotCache::getFromRestart(Pointer<GridGeometry<NDIM>> grid_geom)
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
        db = restart_db->getDatabase(d_object_name);
    else
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to " << d_object_name
                                 << " not found in restart file.\n");
    // Get data from the database
    unsigned int num_snapshots_stored = db->getInteger("num_snapshots_stored");
    d_snapshots.reserve(num_snapshots_stored);
    ComponentSelector comp_selector;
    comp_selector.setFlag(d_snapshot_idx);
    // Because we wrote these when they were sorted, we can retrieve them in the same order.
    // NOTE: this will give the hierarchies different names upon restarts. As an alternative, instead of writing
    // databases in the stored order, we can write the number in the names of the hierarchy. This would require keeping
    // track of the order we store hierarchies and writing those labels to the restart database.
    for (unsigned int i = 0; i < num_snapshots_stored; ++i)
    {
        const double time = db->getDouble("time_" + std::to_string(i));
        Pointer<PatchHierarchy<NDIM>> hierarchy =
            new PatchHierarchy<NDIM>(d_object_name + "::SnapshotHierarchy_" + std::to_string(i), grid_geom, false);
        hierarchy->getFromDatabase(db->getDatabase("snapshot_hierarchy_" + std::to_string(i)), comp_selector);
        d_snapshots.push_back(std::make_pair(time, hierarchy));
    }
}
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
