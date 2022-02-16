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

#include "ibamr/config.h"

#include "ibamr/HierarchyTimeInterpolator.h"

#include "HierarchyDataOpsManager.h"

#include "ibamr/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
double
map_to_period(const double t_start, const double t_end, double time)
{
    const double period = t_end - t_start;
#ifndef NDEBUG
    // We only deal with positive period lengths.
    TBOX_ASSERT(period > 0.0);
#endif

    while (time < t_end)
    {
        // Subtract period until we're less than t_end.
        time -= period;
    }
    while (time > t_start)
    {
        // Add period until we're greater than t_start
        time += period;
    }
#ifndef NDEBUG
    TBOX_ASSERT(time <= t_end && time >= t_start);
#endif
    return time;
}

void
allocate_patch_data(int idx, Pointer<PatchHierarchy<NDIM> > hierarchy, int coarsest_ln, int finest_ln)
{
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(idx)) level->allocatePatchData(idx);
    }
}

void
deallocate_patch_data(int idx, Pointer<PatchHierarchy<NDIM> > hierarchy, int coarsest_ln, int finest_ln)
{
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(idx)) level->deallocatePatchData(idx);
    }
}

template <class VariableType>
HierarchyTimeInterpolator<VariableType>::HierarchyTimeInterpolator(std::string object_name,
                                                                   Pointer<Database> input_db,
                                                                   Pointer<PatchHierarchy<NDIM> > hierarchy)
    : d_object_name(std::move(object_name))
{
    // Get depth
    d_depth = input_db->getIntegerWithDefault("depth", d_depth);
    // Fill in period information
    d_t_start = input_db->getDouble("t_start");
    d_t_end = input_db->getDouble("t_end");
    d_t_period = d_t_end - d_t_start;
    // Snapshot points are equally spaced between start and end
    int num_snapshots = input_db->getInteger("num_snapshots");
    double dt = d_t_period / num_snapshots;
    for (int i = 0; i < num_snapshots; ++i) d_snapshot_time_pts.insert(d_t_start + static_cast<double>(i) * dt);

    // Set refine type
    d_mean_refine_type = input_db->getStringWithDefault("refine_type", d_mean_refine_type);
    commonConstructor(input_db, hierarchy);
}

template <class VariableType>
HierarchyTimeInterpolator<VariableType>::HierarchyTimeInterpolator(std::string object_name,
                                                                   Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                                   std::set<double> snapshot_time_points,
                                                                   const double t_start,
                                                                   const double t_end,
                                                                   const int depth)
    : d_object_name(std::move(object_name)),
      d_depth(depth),
      d_t_start(t_start),
      d_t_end(t_end),
      d_t_period(t_end - t_start),
      d_snapshot_time_pts(std::move(snapshot_time_points))
{
    commonConstructor(nullptr, hierarchy);
}

template <class VariableType>
void
HierarchyTimeInterpolator<VariableType>::commonConstructor(Pointer<Database> input_db,
                                                           Pointer<PatchHierarchy<NDIM> > hierarchy)
{
    // Register the scratch variable
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_scratch_var = new VariableType(d_object_name + "::Variable", d_depth);
    d_scratch_idx =
        var_db->registerVariableAndContext(d_scratch_var, var_db->getContext(d_object_name + "::ctx"), 1 /*ghosts*/);

    // Create the cache
    d_snapshot_cache = std::unique_ptr<SnapshotCache<VariableType> >(
        new SnapshotCache<VariableType>(d_object_name + "::SnapshotCache", input_db));

    // Create the hierarchy data ops
    auto hier_math_ops = HierarchyDataOpsManager<NDIM>::getManager();
    d_hier_data_ops = hier_math_ops->getOperationsDouble(d_scratch_var, hierarchy);
}

template <class VariableType>
HierarchyTimeInterpolator<VariableType>::~HierarchyTimeInterpolator()
{
    clearSnapshots();
}

template <class VariableType>
void
HierarchyTimeInterpolator<VariableType>::clearSnapshots()
{
    // Reset the snapshots.
    d_snapshot_cache->clearSnapshots();

    // Reset the number of updates.
    d_idx_num_updates_map.clear();
}

template <class VariableType>
bool
HierarchyTimeInterpolator<VariableType>::updateTimeAveragedSnapshot(const int u_idx,
                                                                    double time,
                                                                    Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                                    const std::string& refine_type,
                                                                    const int wgt_idx,
                                                                    const double tol)
{
    if (d_t_period == 0.0)
    {
        // Ensure time is the only time point we are storing
        time = *(d_snapshot_time_pts.begin());
    }
    else
    {
        time = map_to_period(d_t_start, d_t_end, time);
    }
    time = getTimePt(time, tol);
    // If this is the first snapshot, record it.
    if (d_idx_num_updates_map[time] == 0)
    {
        d_snapshot_cache->setSnapshot(u_idx, time, hierarchy);
        d_idx_steady_state_map[time] = false;
        return false;
    }
    // Otherwise, we need to update the mean.
    // Fill the scratch index with the current mean
    allocate_patch_data(d_scratch_idx, hierarchy, 0, hierarchy->getFinestLevelNumber());
    d_snapshot_cache->getSnapshot(d_scratch_idx, time, hierarchy, refine_type, tol);
    // The mean is updated via
    // u_avg = u_avg + (1/N)*(u - u_avg)
    // Note first mean is already calculated, so we increment steady state idx.
    const int N = ++d_idx_num_updates_map[time];
    const double weight = 1.0 / static_cast<double>(N);
    d_hier_data_ops->linearSum(d_scratch_idx, weight, u_idx, 1.0 - weight, d_scratch_idx);
    // Update snapshot
    d_snapshot_cache->updateSnapshot(d_scratch_idx, time, hierarchy, tol);

    // Determine if we are at a steady state
    d_hier_data_ops->subtract(d_scratch_idx, d_scratch_idx, u_idx);
    const double L2_norm = d_hier_data_ops->L2Norm(d_scratch_idx, wgt_idx);
    plog << "At time " << time << ", the L^2 norm of the change in u_mean is " << L2_norm << "\n";
    deallocate_patch_data(d_scratch_idx, hierarchy, 0, hierarchy->getFinestLevelNumber());
    if (L2_norm <= d_periodic_thresh)
    {
        d_idx_steady_state_map[time] = true;
        return true;
    }
    return false;
}

template <class VariableType>
void
HierarchyTimeInterpolator<VariableType>::fillSnapshotAtTime(const int u_idx,
                                                            double time,
                                                            Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                            const std::string& refine_type)
{
    // If the period is 0.0, then just return the mean
    if (d_t_period == 0.0)
    {
        d_snapshot_cache->getSnapshot(u_idx, time, hierarchy, "LINEAR_REFINE", 1.0e-8);
        return;
    }
    // Ensure that time is between our periodic time values
    time = map_to_period(d_t_start, d_t_end, time);
    // Determine the correct snapshot index
    auto it_low = d_snapshot_time_pts.lower_bound(time);
    auto it_up = d_snapshot_time_pts.upper_bound(time);
    // Wrap around to beginning if we need to.
    // This can happend if time is between the last and first stored value (Assuming periodicity).
    double t_low = *it_low;
    double t_up = it_up == d_snapshot_time_pts.end() ? *(d_snapshot_time_pts.begin()) + d_t_period : *it_up;

    // Now grab the indices from the snapshot.
    allocate_patch_data(d_scratch_idx, hierarchy, 0, hierarchy->getFinestLevelNumber());
    d_snapshot_cache->getSnapshot(u_idx, t_low, hierarchy, refine_type, 1.0e-8);
    d_snapshot_cache->getSnapshot(d_scratch_idx, t_up, hierarchy, refine_type, 1.0e-8);
    // We have the upper and lower indices, now we can interpolate in time.
    d_hier_data_ops->linearSum(u_idx, (time - t_low) / d_t_period, u_idx, -(time - t_up) / d_t_period, d_scratch_idx);
    deallocate_patch_data(d_scratch_idx, hierarchy, 0, hierarchy->getFinestLevelNumber());
}

template <class VariableType>
double
HierarchyTimeInterpolator<VariableType>::getTimePt(const double time, const double tol)
{
    double t_low = *(d_snapshot_time_pts.lower_bound(time));
    double t_up = *(d_snapshot_time_pts.upper_bound(time));
    if (std::abs(t_low - time) < tol)
        return t_low;
    else if (std::abs(t_up - time) < tol)
        return t_up;
    else
        TBOX_ERROR("Time point is not within the given tolerance!");
    return 0.0;
}
// Instantiate the viable templates
template class HierarchyTimeInterpolator<SAMRAI::pdat::CellVariable<NDIM, double> >;
template class HierarchyTimeInterpolator<SAMRAI::pdat::SideVariable<NDIM, double> >;
template class HierarchyTimeInterpolator<SAMRAI::pdat::NodeVariable<NDIM, double> >;
template class HierarchyTimeInterpolator<SAMRAI::pdat::EdgeVariable<NDIM, double> >;
template class HierarchyTimeInterpolator<SAMRAI::pdat::FaceVariable<NDIM, double> >;
//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
