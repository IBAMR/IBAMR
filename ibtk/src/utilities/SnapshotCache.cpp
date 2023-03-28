// ---------------------------------------------------------------------
//
// Copyright (c) 2022 - 2022 by the IBAMR developers
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

#include <ibtk/SnapshotCache.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/RestartManager.h>

#include <CartesianGridGeometry.h>
#include <HierarchyDataOpsManager.h>
#include <RefineAlgorithm.h>
#include <RefineOperator.h>
#include <VariableDatabase.h>

#include <algorithm>

#include <ibtk/app_namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
SnapshotCache::SnapshotCache(std::string object_name,
                             Pointer<Variable<NDIM> > var,
                             Pointer<Database> input_db,
                             Pointer<GridGeometry<NDIM> > grid_geom,
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
        Pointer<PatchHierarchy<NDIM> > hierarchy = snapshot.second;
        for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            if (level->checkAllocated(d_snapshot_idx)) level->deallocatePatchData(d_snapshot_idx);
        }
    }

    // Clear snapshots
    d_snapshots.clear();
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    var_db->removePatchDataIndex(d_snapshot_idx);
    d_snapshot_var = nullptr;
    d_snapshot_idx = IBTK::invalid_index;
}

SnapshotCache::value_type
SnapshotCache::getSnapshot(double time, double tol)
{
    auto it = std::find_if(d_snapshots.begin(), d_snapshots.end(), [time, tol](const value_type& t) -> bool {
        return IBTK::abs_equal_eps(t.first, time, tol);
    });
    if (it == d_snapshots.end())
        return std::make_pair(std::numeric_limits<double>::quiet_NaN(), nullptr);
    else
        return *it;
}

void
SnapshotCache::storeSnapshot(const int u_idx, const double time, Pointer<PatchHierarchy<NDIM> > hierarchy)
{
    // Create a snapshot of the patch hierarchy.
    Pointer<PatchHierarchy<NDIM> > snapshot_hierarchy = hierarchy->makeRefinedPatchHierarchy(
        d_object_name + "::SnapshotHierarchy_" + std::to_string(d_snapshots.size()),
        1 /*ratio*/,
        false /*register_for_restart*/);
    // Allocate patch data and copy data to new patch hierarchy.
    int coarsest_ln = 0;
    int finest_ln = snapshot_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > snapshot_level = snapshot_hierarchy->getPatchLevel(ln);
        if (!snapshot_level->checkAllocated(d_snapshot_idx)) snapshot_level->allocatePatchData(d_snapshot_idx, time);
        // We've allocated data, now copy it
        // Note these levels should cover the same space, so we shouldn't need to have a refine operator.
        Pointer<PatchLevel<NDIM> > old_level = hierarchy->getPatchLevel(ln);
        {
            // First fill in from the snapshot hierarchy
            Pointer<RefineAlgorithm<NDIM> > refine_alg = new RefineAlgorithm<NDIM>();
            Pointer<RefineOperator<NDIM> > refine_op = nullptr;
            refine_alg->registerRefine(d_snapshot_idx, u_idx, u_idx, refine_op);
            Pointer<RefineSchedule<NDIM> > schedule = refine_alg->createSchedule(snapshot_level, old_level);
            schedule->fillData(time);
        }
    }

    // Store the index, time, hierarchy, and variable.
    // First determine where to store this index
    auto it = std::find_if(d_snapshots.begin(), d_snapshots.end(), [time](const value_type& snapshot) -> bool {
        return time < snapshot.first;
    });
    d_snapshots.insert(it, std::make_pair(time, snapshot_hierarchy));
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
SnapshotCache::getFromRestart(Pointer<GridGeometry<NDIM> > grid_geom)
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
        Pointer<PatchHierarchy<NDIM> > hierarchy =
            new PatchHierarchy<NDIM>(d_object_name + "::SnapshotHierarchy_" + std::to_string(i), grid_geom, false);
        hierarchy->getFromDatabase(db->getDatabase("snapshot_hierarchy_" + std::to_string(i)), comp_selector);
        d_snapshots.push_back(std::make_pair(time, hierarchy));
    }
}
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
