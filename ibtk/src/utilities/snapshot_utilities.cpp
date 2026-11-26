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

#include <ibtk/ibtk_utilities.h>
#include <ibtk/snapshot_utilities.h>

#include <HierarchyDataOpsManager.h>

#include <algorithm>

#include <ibtk/app_namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
void
update_snapshot(SnapshotCache& cache,
                int u_idx,
                double time,
                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> current_hierarchy,
                double tol)
{
    cache.updateSnapshot(u_idx, time, current_hierarchy, tol);
}

void
fill_snapshot_on_hierarchy(SnapshotCache& cache,
                           const int u_idx,
                           const double time,
                           Pointer<PatchHierarchy<NDIM>> current_hierarchy,
                           const std::string& snapshot_refine_type,
                           const double tol)
{
    cache.fillSnapshotOnHierarchy(u_idx, time, current_hierarchy, snapshot_refine_type, tol);
}

void
fill_snapshot_at_time(SnapshotCache& cache,
                      const int u_idx,
                      double time,
                      const int scr_idx,
                      Pointer<PatchHierarchy<NDIM>> hierarchy,
                      const std::string& refine_type,
                      Pointer<HierarchyDataOpsReal<NDIM, double>> hier_data_ops,
                      const double period)
{
    // If there's only one snapshot, return it
    if (cache.getNumSnapshots() == 1)
    {
        // Get the actual time value from the snapshot
        time = (*cache.begin()).first;
        fill_snapshot_on_hierarchy(cache, u_idx, time, hierarchy, refine_type);
        return;
    }

    // Determine the correct snapshot index
    double snapshot_time_low = 0.0, snapshot_time_up = 0.0;
    double t_low = 0.0, t_up = 0.0;
    auto it_up =
        std::upper_bound(cache.begin(),
                         cache.end(),
                         time,
                         [](const double a, const SnapshotCache::value_type& b) -> bool { return a < b.first; });
    if (period == period && it_up == cache.end())
    {
        // Snapshot is storing periodic values, and the time value is between the last element and the first.
        snapshot_time_up = (*cache.begin()).first;
        t_up = snapshot_time_up + period;
        snapshot_time_low = (*std::next(cache.end(), -1)).first;
        t_low = snapshot_time_low;
    }
    else if (it_up == cache.end())
    {
        TBOX_ERROR("Could not find time points around time " << time << " in snapshot cache!\n");
    }
    else
    {
        // everything is normal
        snapshot_time_up = (*it_up).first;
        t_up = snapshot_time_up;
        snapshot_time_low = (*std::next(it_up, -1)).first;
        t_low = snapshot_time_low;
    }
    // Create a valid hier_data_ops object
    if (!hier_data_ops)
    {
        auto hier_math_ops = HierarchyDataOpsManager<NDIM>::getManager();
        hier_data_ops = hier_math_ops->getOperationsDouble(cache.getVariable(), hierarchy, true);
    }

    // Now grab the indices from the snapshot.
    // Allocate patch data. We need the time point at which u_idx is allocated
    fill_snapshot_on_hierarchy(cache, u_idx, snapshot_time_low, hierarchy, refine_type, 1.0e-8);
    fill_snapshot_on_hierarchy(cache, scr_idx, snapshot_time_up, hierarchy, refine_type, 1.0e-8);
    // We have the upper and lower indices, now we can interpolate in time.
    hier_data_ops->linearSum(u_idx, (t_up - time) / (t_up - t_low), u_idx, (time - t_low) / (t_up - t_low), scr_idx);
}
//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
