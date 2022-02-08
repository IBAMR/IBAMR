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

HierarchyTimeInterpolator::INSStaggeredMeanFlowCalculator(std::string object_name,
                                                          Pointer<Database> input_db,
                                                          Pointer<PatchHierarchy<NDIM> > patch_hierarchy)
    : d_object_name(std::move(object_name)), d_hierarchy(patch_hierarchy)
{
    // blank for now
}

HierarchyTimeInterpolator::~INSStaggeredMeanFlowCalculator()
{
    clearSnapshots();
}

void
HierarchyTimeInterpolator::clearSnapshots()
{
    // Clear the map.
    d_idx_time_map.clear();

    // Deallocate all patch data.
    int coarsest_ln = 0;
    int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (const auto& snapshot_idx : d_snapshot_idxs)
            if (level->checkAllocated(snapshot_idx)) level->deallocatePatchData(snapshot_idx);
    }
    d_num_snapshots_stored = 0;
}

void
HierarchyTimeInterpolator::setVelocitySnapshot(const int u_idx,
                                               const double time,
                                               Pointer<PatchHierarchy<NDIM> > hierarchy)
{
    // Ensure we aren't over-storing
    if (d_num_snapshots_stored >= d_max_snapshots)
        TBOX_ERROR(d_object_name + ": max number of snapshots already stored!\n");
    // This currently only works with a single patch hierarchy
    TBOX_ASSERT(d_hierarchy == hierarchy);

    // Allocate patch data and copy velocity data.
    int coarsest_ln = 0;
    int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const int snapshot_idx = d_snapshot_idxs[d_num_snapshots_stored];
        if (!level->checkAllocated(snapshot_idx)) level->allocatePatchData(snapshot_idx);
        // We've allocated data, now copy it
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM, double> > snapshot_data = patch->getPatchData(snapshot_idx);
            Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_idx);
            snapshot_data->copy(*u_data);
        }
    }

    // Store this index into our map
    d_idx_time_map[time] = d_num_snapshots_stored++;
}

void
HierarchyTimeInterpolator::fillVelocityAtTime(const int u_idx, double time, Pointer<PatchHierarchy<NDIM> > hierarchy)
{
    // This currently only works with a single patch hierarchy
    TBOX_ASSERT(d_hierarchy == hierarchy);
    // Ensure that time is between our periodic time values
    time = map_to_period(d_t_start, d_t_end, time);
    // Determine the correct snapshot index
    auto it_low = d_idx_time_map.lower_bound(time);
    auto it_up = d_idx_time_map.upper_bound(time);
    // Wrap around to beginning if we need to.
    // This can happend if time is between the last and first stored value (Assuming periodicity).
    it_up = it_up == d_idx_time_map.end() ? d_idx_time_map.begin() : it_up;

    // TODO We should check that the times are in the right order. Note that it_up can be the beginning of the map...
    double t_low = (*it_low).first;
    int idx_low = (*it_low).second;
    double t_up = (*it_up).first;
    int idx_up = (*it_up).second;
    // We have the upper and lower indices, now we can interpolate in time.
    d_hier_sc_data_ops->linearSum(u_idx, (time - t_low) / d_t_period, idx_low, -(time - t_up) / d_t_period, idx_up);
}

IBTK::VectorNd
HierarchyTimeInterpolator::computeTimeAvgVelocity()
{
    // blank for now
    return IBTK::VectorNd::Zero();
}

double
HierarchyTimeInterpolator::computeTKE()
{
    // blank for now
    return 0.0;
}
//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
