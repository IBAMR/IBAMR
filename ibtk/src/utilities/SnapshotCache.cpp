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

#include "ibtk/SnapshotCache.h"

#include "CartesianGridGeometry.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "VariableDatabase.h"

#include <algorithm>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
template <class VariableType>
SnapshotCache<VariableType>::SnapshotCache(std::string object_name, Pointer<Database> input_db)
    : d_object_name(std::move(object_name))
{
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_ctx = var_db->getContext(d_object_name + "::context");

    if (input_db)
    {
        d_snapshot_refine_type = input_db->getStringWithDefault("refine_type", d_snapshot_refine_type);
        if (input_db->keyExists("gcw")) input_db->getIntegerArray("gcw", &d_gcw[0], NDIM);
        d_depth = input_db->getIntegerWithDefault("depth", d_depth);
    }
}

template <class VariableType>
SnapshotCache<VariableType>::~SnapshotCache()
{
    clearSnapshots();
}

template <class VariableType>
void
SnapshotCache<VariableType>::clearSnapshots()
{
    // Clear the map.
    d_idx_time_map.clear();
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    // Deallocate all patch data.
    for (const auto& idx : d_snapshot_idxs)
    {
        const Pointer<PatchHierarchy<NDIM> >& hierarchy = d_snapshot_hierarchies[idx];
        int coarsest_ln = 0;
        int finest_ln = hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            if (level->checkAllocated(idx)) level->deallocatePatchData(idx);
        }
        // Delete variable from variable database
        var_db->removePatchDataIndex(idx);
    }

    // Clear snapshots
    d_snapshot_hierarchies.clear();
    d_snapshot_idxs.clear();
    d_snapshot_vars.clear();
    d_num_snapshots_stored = 0;
}

template <class VariableType>
void
SnapshotCache<VariableType>::setSnapshot(const int u_idx,
                                         const double time,
                                         Pointer<PatchHierarchy<NDIM> > current_hierarchy)
{
    ++d_num_snapshots_stored;
    // Push back a new variable
    Pointer<VariableType> snapshot_var =
        new VariableType(d_object_name + "::SnapshotVar_" + std::to_string(d_num_snapshots_stored), d_depth);
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    int snapshot_idx = var_db->registerVariableAndContext(snapshot_var, d_ctx, d_gcw);
    d_snapshot_idxs.push_back(snapshot_idx);
    d_snapshot_vars[snapshot_idx] = snapshot_var;

    // Now fill in snapshot
    fillSnapshot(u_idx, current_hierarchy, snapshot_idx, time);
}

template <class VariableType>
void
SnapshotCache<VariableType>::updateSnapshot(const int u_idx,
                                            const double time,
                                            Pointer<PatchHierarchy<NDIM> > current_hierarchy,
                                            const double tol)
{
    // Make sure we are on a stored snapshot.
    auto it =
        std::find_if(d_idx_time_map.begin(),
                     d_idx_time_map.end(),
                     [time, tol](const std::pair<int, double>& t) -> bool { return std::abs(t.second - time) < tol; });
    if (it == d_idx_time_map.end())
        TBOX_ERROR("Snapshot at time: " << time << " with tolerance " << tol << " does not exist!\n");

    // We have the snapshot index. Clear the current snapshot and copy the new one.
    const std::pair<int, double>& snapshot_idx_t_pair = *it;
    const int snapshot_idx = snapshot_idx_t_pair.first;
    const Pointer<PatchHierarchy<NDIM> >& snapshot_hierarchy = d_snapshot_hierarchies.at(snapshot_idx);

    // First clear the current snapshot
    const int coarsest_ln = 0;
    const int finest_ln = snapshot_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = snapshot_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(snapshot_idx)) level->deallocatePatchData(snapshot_idx);
    }

    // Now fill in snapshot
    fillSnapshot(u_idx, current_hierarchy, snapshot_idx, time);
}

template <class VariableType>
void
SnapshotCache<VariableType>::getSnapshot(const int u_idx,
                                         const double time,
                                         Pointer<PatchHierarchy<NDIM> > current_hierarchy,
                                         const double tol)
{
    getSnapshot(u_idx, time, current_hierarchy, d_snapshot_refine_type, tol);
}

template <class VariableType>
void
SnapshotCache<VariableType>::getSnapshot(const int u_idx,
                                         const double time,
                                         Pointer<PatchHierarchy<NDIM> > current_hierarchy,
                                         const std::string& snapshot_refine_type,
                                         const double tol)
{
    // Make sure we are on a stored snapshot.
    auto it =
        std::find_if(d_idx_time_map.begin(),
                     d_idx_time_map.end(),
                     [time, tol](const std::pair<int, double>& t) -> bool { return std::abs(t.second - time) < tol; });
    if (it == d_idx_time_map.end())
        TBOX_ERROR("Snapshot at time: " << time << " with tolerance " << tol << " does not exist!\n");

    // We have our snapshot index. Now copy the data to u_idx.
    const std::pair<int, double>& snapshot_idx_t_pair = *it;
    const int snapshot_idx = snapshot_idx_t_pair.first;
    const Pointer<PatchHierarchy<NDIM> >& snapshot_hierarchy = d_snapshot_hierarchies.at(snapshot_idx);
    const Pointer<VariableType>& snapshot_var = d_snapshot_vars.at(snapshot_idx);
    // Make sure the current hierarchy and snapshot hierarchy have the same number of levels
    TBOX_ASSERT(current_hierarchy->getFinestLevelNumber() == snapshot_hierarchy->getFinestLevelNumber());
    // Now transfer data from the snapshot to the current hierarchy
    int coarsest_ln = 0;
    int finest_ln = snapshot_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > cur_level = current_hierarchy->getPatchLevel(ln);
        Pointer<PatchLevel<NDIM> > snp_level = snapshot_hierarchy->getPatchLevel(ln);

        // RefineSchedule requires that patch data be stored at the same points.
        // We reset the snapshot patch data time to that of the current time.
        double allocated_time =
            cur_level->getPatch(*PatchLevel<NDIM>::Iterator(cur_level))->getPatchData(u_idx)->getTime();
        // Note we need scratch allocated on the current hierarchy. Use the snapshot for that.
        cur_level->allocatePatchData(snapshot_idx, allocated_time);
        for (PatchLevel<NDIM>::Iterator p(snp_level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = snp_level->getPatch(p());
            Pointer<PatchData<NDIM> > data = patch->getPatchData(snapshot_idx);
            data->setTime(allocated_time);
        }

        // Now copy the data.
        Pointer<RefineAlgorithm<NDIM> > refine_alg = new RefineAlgorithm<NDIM>();
        Pointer<CartesianGridGeometry<NDIM> > grid_geom = current_hierarchy->getGridGeometry();
        Pointer<RefineOperator<NDIM> > refine_op = grid_geom->lookupRefineOperator(snapshot_var, snapshot_refine_type);
        refine_alg->registerRefine(u_idx, snapshot_idx, snapshot_idx, refine_op);
        Pointer<RefineSchedule<NDIM> > schedule =
            refine_alg->createSchedule(cur_level, snp_level, ln - 1, current_hierarchy);
        schedule->fillData(allocated_time);
    }

    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = current_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(snapshot_idx);
    }
}

template <class VariableType>
void
SnapshotCache<VariableType>::fillSnapshot(const int u_idx,
                                          Pointer<PatchHierarchy<NDIM> > hierarchy,
                                          const int snapshot_idx,
                                          const double time)
{
    // Create a snapshot of the patch hierarchy.
    Pointer<PatchHierarchy<NDIM> > snapshot_hierarchy = hierarchy->makeRefinedPatchHierarchy(
        d_object_name + "::SnapshotHierarchy_" + std::to_string(d_num_snapshots_stored),
        1 /*ratio*/,
        false /*register_for_restart*/);
    // Allocate patch data and copy data to new patch hierarchy.
    int coarsest_ln = 0;
    int finest_ln = snapshot_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > snapshot_level = snapshot_hierarchy->getPatchLevel(ln);
        if (!snapshot_level->checkAllocated(snapshot_idx)) snapshot_level->allocatePatchData(snapshot_idx, time);
        // We've allocated data, now copy it
        // Note these levels should cover the same space, so we shouldn't need to have a refine operator.
        Pointer<PatchLevel<NDIM> > old_level = hierarchy->getPatchLevel(ln);
        {
            // First fill in from the snapshot hierarchy
            Pointer<RefineAlgorithm<NDIM> > refine_alg = new RefineAlgorithm<NDIM>();
            Pointer<RefineOperator<NDIM> > refine_op = nullptr;
            refine_alg->registerRefine(snapshot_idx, u_idx, u_idx, refine_op);
            Pointer<RefineSchedule<NDIM> > schedule = refine_alg->createSchedule(snapshot_level, old_level);
            schedule->fillData(time);
        }
    }

    // Store the index, time, hierarchy, and variable.
    d_idx_time_map[snapshot_idx] = time;
    d_snapshot_hierarchies[snapshot_idx] = snapshot_hierarchy;
}

// Instantiate the viable templates
template class SnapshotCache<SAMRAI::pdat::CellVariable<NDIM, double> >;
template class SnapshotCache<SAMRAI::pdat::SideVariable<NDIM, double> >;
template class SnapshotCache<SAMRAI::pdat::NodeVariable<NDIM, double> >;
template class SnapshotCache<SAMRAI::pdat::EdgeVariable<NDIM, double> >;
template class SnapshotCache<SAMRAI::pdat::FaceVariable<NDIM, double> >;
//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
