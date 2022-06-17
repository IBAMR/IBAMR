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

#include <ibtk/ibtk_utilities.h>
#include <ibtk/snapshot_utilities.h>

#include <CartesianGridGeometry.h>
#include <HierarchyDataOpsManager.h>
#include <RefineAlgorithm.h>
#include <RefineOperator.h>

#include <algorithm>

#include <ibtk/app_namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
void
update_snapshot(SnapshotCache& cache,
                int u_idx,
                double time,
                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > current_hierarchy,
                double tol)
{
    // Make sure we are on a stored snapshot.
    const std::pair<double, Pointer<PatchHierarchy<NDIM> > >& snapshot = cache.getSnapshot(time, tol);
    if (!snapshot.second || !IBTK::abs_equal_eps(snapshot.first, time, tol))
        TBOX_ERROR("Snapshot at time: " << time << " with tolerance " << tol << " does not exist!\n");

    // We need the snapshot index
    int snapshot_idx = cache.getPatchIndex();

    const Pointer<PatchHierarchy<NDIM> >& snapshot_hierarchy = snapshot.second;

    // Note we can not do the usual trick of making a refined patch hierarchy with refinement ratio 1 because we can't
    // copy the snapshot_hierarchy's name. Instead, we delete all the snapshot_hierarchy's levels and replace them with
    // copies of the current_hierarchy.
    int coarsest_ln = 0;
    int finest_ln = current_hierarchy->getFinestLevelNumber();
#ifndef NDEBUG
    TBOX_ASSERT(snapshot_hierarchy->getFinestLevelNumber() == finest_ln);
#endif
    Pointer<GridGeometry<NDIM> > new_geometry =
        current_hierarchy->getGridGeometry()->makeRefinedGridGeometry("SnapshotGeometry", 1, false);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        snapshot_hierarchy->removePatchLevel(ln);
        Pointer<PatchLevel<NDIM> > new_level = new PatchLevel<NDIM>();
        Pointer<PatchLevel<NDIM> > current_level = current_hierarchy->getPatchLevel(ln);
        new_level->setRefinedPatchLevel(current_level, 1, new_geometry);

        // Now allocate data on the new level
        new_level->allocatePatchData(snapshot_idx, time);
        Pointer<RefineAlgorithm<NDIM> > refine_alg = new RefineAlgorithm<NDIM>();
        Pointer<RefineOperator<NDIM> > refine_op = nullptr;
        refine_alg->registerRefine(snapshot_idx, u_idx, u_idx, refine_op);
        Pointer<RefineSchedule<NDIM> > schedule = refine_alg->createSchedule(new_level, current_level);
        schedule->fillData(time);
    }
}

void
fill_snapshot_on_hierarchy(SnapshotCache& cache,
                           const int u_idx,
                           const double time,
                           Pointer<PatchHierarchy<NDIM> > current_hierarchy,
                           const std::string& snapshot_refine_type,
                           const double tol)
{
    // Make sure we are on a stored snapshot.
    const std::pair<double, Pointer<PatchHierarchy<NDIM> > >& snapshot = cache.getSnapshot(time, tol);
    if (!snapshot.second || !IBTK::abs_equal_eps(snapshot.first, time, tol))
        TBOX_ERROR("Snapshot at time: " << time << " with tolerance " << tol << " does not exist!\n");

    // We need the snapshot index on the snapshot hierarchy and a scratch index on both hierarchies. We use the snapshot
    // index as the scratch index, but maybe we should clone the index?
    int snapshot_idx = cache.getPatchIndex();
    int scr_idx = snapshot_idx;

    // We have our snapshot index. Now copy the data to u_idx.
    const Pointer<PatchHierarchy<NDIM> >& snapshot_hierarchy = snapshot.second;
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
        cur_level->allocatePatchData(scr_idx, allocated_time);
        snp_level->setTime(allocated_time, snapshot_idx);

        // Now copy the data.
        Pointer<RefineAlgorithm<NDIM> > refine_alg = new RefineAlgorithm<NDIM>();
        Pointer<CartesianGridGeometry<NDIM> > grid_geom = current_hierarchy->getGridGeometry();
        Pointer<RefineOperator<NDIM> > refine_op =
            grid_geom->lookupRefineOperator(cache.getVariable(), snapshot_refine_type);
        refine_alg->registerRefine(u_idx, snapshot_idx, scr_idx, refine_op);
        Pointer<RefineSchedule<NDIM> > schedule =
            refine_alg->createSchedule(cur_level, snp_level, ln - 1, current_hierarchy);
        schedule->fillData(allocated_time);
    }

    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = current_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(scr_idx);
    }
}

void
fill_snapshot_at_time(SnapshotCache& cache,
                      const int u_idx,
                      double time,
                      const int scr_idx,
                      Pointer<PatchHierarchy<NDIM> > hierarchy,
                      const std::string& refine_type,
                      Pointer<HierarchyDataOpsReal<NDIM, double> > hier_data_ops,
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
    auto it_up = std::upper_bound(
        cache.begin(), cache.end(), time, [](const double a, const SnapshotCache::value_type& b) -> bool {
            return a < b.first;
        });
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
        auto var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<Variable<NDIM> > var;
        var_db->mapIndexToVariable(u_idx, var);
        Pointer<HierarchyDataOpsReal<NDIM, double> > hier_data_ops =
            hier_math_ops->getOperationsDouble(cache.getVariable(), hierarchy, true);
    }

    // Now grab the indices from the snapshot.
    // Allocate patch data. We need the time point at which u_idx is allocated
    fill_snapshot_on_hierarchy(cache, u_idx, snapshot_time_low, hierarchy, refine_type, 1.0e-8);
    fill_snapshot_on_hierarchy(cache, scr_idx, snapshot_time_up, hierarchy, refine_type, 1.0e-8);
    // We have the upper and lower indices, now we can interpolate in time.
    hier_data_ops->linearSum(u_idx, (time - t_low) / (t_up - t_low), u_idx, -(time - t_up) / (t_up - t_low), scr_idx);
}
//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
