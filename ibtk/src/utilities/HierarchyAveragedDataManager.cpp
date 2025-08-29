// ---------------------------------------------------------------------
//
// Copyright (c) 2022 - 2025 by the IBAMR developers
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
#include <ibtk/HierarchyAveragedDataManager.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/snapshot_utilities.h>

#include <HierarchyDataOpsManager.h>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
double
map_to_period(const double t_start, const double t_end, double time)
{
    const double period = t_end - t_start;
    if (IBTK::abs_equal_eps(period, 0.0)) return t_start;
#ifndef NDEBUG
    // We only deal with positive period lengths.
    TBOX_ASSERT(period > 0.0);
#endif
    // Shift so that t_start is 0.0
    time -= t_start;
    // If time is BEFORE t_start, we will get a negative value here
    time = std::abs(std::fmod(time, period));
    time += t_start;
    return time;
}

namespace
{
void
allocate_patch_data(const int idx,
                    const double time,
                    Pointer<PatchHierarchy<NDIM> > hierarchy,
                    const int coarsest_ln,
                    const int finest_ln)
{
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(idx)) level->allocatePatchData(idx, time);
    }
}

void
deallocate_patch_data(const int idx,
                      Pointer<PatchHierarchy<NDIM> > hierarchy,
                      const int coarsest_ln,
                      const int finest_ln)
{
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(idx)) level->deallocatePatchData(idx);
    }
}
} // namespace

HierarchyAveragedDataManager::HierarchyAveragedDataManager(std::string object_name,
                                                           Pointer<Variable<NDIM> > var,
                                                           Pointer<Database> input_db,
                                                           Pointer<GridGeometry<NDIM> > grid_geom,
                                                           bool register_for_restart)
    : d_object_name(std::move(object_name)),
      d_var(var),
      d_snapshot_cache(d_object_name + "::SnapshotCache", var, input_db, grid_geom, register_for_restart)
{
    // Fill in period information
    d_period_start = input_db->getDouble("period_start");
    d_period_end = input_db->getDouble("period_end");
    d_period_length = d_period_end - d_period_start;
    d_periodic_thresh = input_db->getDouble("threshold");
    // Snapshot points are equally spaced between start and end
    int num_snapshots = input_db->getInteger("num_snapshots");
    double dt = d_period_length / num_snapshots;
    for (int i = 0; i < num_snapshots; ++i) d_snapshot_time_pts.insert(d_period_start + static_cast<double>(i) * dt);

    commonConstructor(input_db);

    auto restart_manager = RestartManager::getManager();
    if (register_for_restart)
    {
        restart_manager->registerRestartItem(d_object_name, this);
        auto var_db = VariableDatabase<NDIM>::getDatabase();
        var_db->registerPatchDataForRestart(d_scratch_idx);
    }

    bool from_restart = restart_manager->isFromRestart();
    if (from_restart) getFromRestart();
}

HierarchyAveragedDataManager::HierarchyAveragedDataManager(std::string object_name,
                                                           Pointer<Variable<NDIM> > var,
                                                           Pointer<Database> input_db,
                                                           std::set<double> snapshot_time_points,
                                                           const double period_start,
                                                           const double period_end,
                                                           const double threshold,
                                                           Pointer<GridGeometry<NDIM> > grid_geom,
                                                           bool register_for_restart)
    : d_object_name(std::move(object_name)),
      d_var(var),
      d_period_start(period_start),
      d_period_end(period_end),
      d_period_length(period_end - period_start),
      d_periodic_thresh(threshold),
      d_snapshot_time_pts(std::move(snapshot_time_points)),
      d_snapshot_cache(d_object_name + "::SnapshotCache", var, nullptr, grid_geom, register_for_restart)
{
    // Bad things can happen if the snapshot time points contain NaN's.
    for (const auto& time_point : snapshot_time_points)
    {
        if (std::isnan(time_point))
            TBOX_ERROR(d_object_name + "::HierarchyAveragedDataManager(): NaN found in the time point!\n");
    }
    commonConstructor(input_db);
}

void
HierarchyAveragedDataManager::commonConstructor(Pointer<Database> input_db)
{
    // Get some information from the database
    d_enable_logging = input_db->getBool("enable_logging");
    d_output_data = input_db->getBool("output_data");
    std::string dir_dump_name = input_db->getStringWithDefault("dir_dump_name", "./");
    d_mean_refine_type = input_db->getStringWithDefault("refine_type", d_mean_refine_type);

    // Register the scratch variable. Note we need the scratch variable to have sufficient ghost cell width for whatever
    // refinement operation must be done. We should read this from the input database.
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    d_scratch_idx =
        var_db->registerVariableAndContext(d_var, var_db->getContext(d_object_name + "::ctx"), 1 /*ghosts*/);

    // Determine depth for drawing.
    int depth = 1;
    Pointer<CellVariable<NDIM, double> > c_var = d_var;
    Pointer<NodeVariable<NDIM, double> > n_var = d_var;
    Pointer<SideVariable<NDIM, double> > s_var = d_var;
    Pointer<FaceVariable<NDIM, double> > f_var = d_var;
    Pointer<EdgeVariable<NDIM, double> > e_var = d_var;
    if (c_var)
    {
        Pointer<CellDataFactory<NDIM, double> > data_factory = c_var->getPatchDataFactory();
        depth = data_factory->getDefaultDepth();
    }
    else if (n_var)
    {
        Pointer<NodeDataFactory<NDIM, double> > data_factory = n_var->getPatchDataFactory();
        depth = data_factory->getDefaultDepth();
    }
    else if (s_var || f_var)
    {
        depth = NDIM;
    }
    else if (e_var)
    {
        Pointer<EdgeDataFactory<NDIM, double> > data_factory = e_var->getPatchDataFactory();
        depth = data_factory->getDefaultDepth();
    }
    d_mean_var = new CellVariable<NDIM, double>(d_object_name + "::MeanVar", depth);
    d_dev_var = new CellVariable<NDIM, double>(d_object_name + "::DevVar", depth);
    if (d_output_data)
    {
        d_mean_idx = var_db->registerVariableAndContext(d_mean_var, var_db->getContext(d_object_name + "::ctx"), 1);
        d_dev_idx = var_db->registerVariableAndContext(d_dev_var, var_db->getContext(d_object_name + "::ctx"), 1);
        d_visit_data_writer = std::make_unique<VisItDataWriter<NDIM> >(d_object_name + "::VisitWriter", dir_dump_name);
        d_visit_data_writer->registerPlotQuantity("mean_flow_field", "VECTOR", d_mean_idx);
        d_visit_data_writer->registerPlotQuantity("deviation", "VECTOR", d_dev_idx);
    }
}

void
HierarchyAveragedDataManager::clearSnapshots()
{
    // Reset the snapshots.
    d_snapshot_cache.clearSnapshots();

    // Reset the number of updates.
    d_idx_num_updates_map.clear();
}

bool
HierarchyAveragedDataManager::updateTimeAveragedSnapshot(const int u_idx,
                                                         double time,
                                                         Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                         const std::string& refine_type,
                                                         const int wgt_idx,
                                                         const double tol)
{
    // Create the hierarchy data ops
    auto hier_math_ops = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<HierarchyDataOpsReal<NDIM, double> > hier_data_ops =
        hier_math_ops->getOperationsDouble(d_var, hierarchy, true /*get_unique*/);

    if (d_period_length == 0.0)
    {
        // Ensure time is the only time point we are storing
        time = *(d_snapshot_time_pts.begin());
    }
    else
    {
        time = map_to_period(d_period_start, d_period_end, time);
    }
    time = getTimePoint(time, tol);
    // If this is the first snapshot, record it.
    if (d_idx_num_updates_map[time] == 0)
    {
        // Need to copy the data to the scratch index.
        allocate_patch_data(d_scratch_idx, time, hierarchy, 0, hierarchy->getFinestLevelNumber());
        hier_data_ops->copyData(d_scratch_idx, u_idx);
        d_snapshot_cache.storeSnapshot(d_scratch_idx, time, hierarchy);
        d_idx_steady_state_map[time] = false;
        d_idx_num_updates_map[time] = 1;
        return false;
    }

    // Otherwise, we need to update the mean.
    // Fill the scratch index with the current snapshot
    allocate_patch_data(d_scratch_idx, time, hierarchy, 0, hierarchy->getFinestLevelNumber());
    fill_snapshot_on_hierarchy(d_snapshot_cache, d_scratch_idx, time, hierarchy, refine_type, tol);
    // The mean is updated via
    // u_avg^N = u_avg^(N-1) + (1/N)*(u - u_avg^(N-1))
    // Note first mean is already calculated, so we increment steady state idx.
    const double N = static_cast<double>(d_idx_num_updates_map[time]++);
    hier_data_ops->resetLevels(0, hierarchy->getFinestLevelNumber());
    hier_data_ops->linearSum(d_scratch_idx, 1.0 / (N + 1.0), u_idx, N / (N + 1.0), d_scratch_idx);
    // Update snapshot with new mean
    update_snapshot(d_snapshot_cache, d_scratch_idx, time, hierarchy, tol);

    // Output data if necessary.
    if (d_output_data)
    {
        pout << "Printing out averaged data\n";
        allocate_patch_data(d_mean_idx, time, hierarchy, 0, hierarchy->getFinestLevelNumber());
        allocate_patch_data(d_dev_idx, time, hierarchy, 0, hierarchy->getFinestLevelNumber());
        // Copy data to visit index
        Pointer<HierarchyGhostCellInterpolation> ghost_fill = new HierarchyGhostCellInterpolation();
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> itc = { ITC(
            d_scratch_idx, "CONSERVATIVE_LINEAR_REFINE", false, "CONSERVATIVE_COARSEN", "LINEAR") };
        ghost_fill->initializeOperatorState(itc, hierarchy);
        HierarchyMathOps hier_math_ops("HierarchyMathOps", hierarchy);
        hier_math_ops.resetLevels(0, hierarchy->getFinestLevelNumber());
        hier_math_ops.setPatchHierarchy(hierarchy);
        Pointer<SideVariable<NDIM, double> > sc_var = d_var;
        Pointer<NodeVariable<NDIM, double> > nc_var = d_var;
        Pointer<EdgeVariable<NDIM, double> > ec_var = d_var;
        Pointer<FaceVariable<NDIM, double> > fc_var = d_var;
        Pointer<CellVariable<NDIM, double> > cc_var = d_var;
        if (sc_var) hier_math_ops.interp(d_mean_idx, d_mean_var, d_scratch_idx, sc_var, ghost_fill, time, false);
        if (nc_var) hier_math_ops.interp(d_mean_idx, d_mean_var, d_scratch_idx, nc_var, ghost_fill, time, false);
        if (ec_var) hier_math_ops.interp(d_mean_idx, d_mean_var, d_scratch_idx, ec_var, ghost_fill, time, false);
        if (fc_var) hier_math_ops.interp(d_mean_idx, d_mean_var, d_scratch_idx, fc_var, ghost_fill, time, false);
        if (cc_var) hier_data_ops->copyData(d_mean_idx, d_scratch_idx);
        // Determine our convergence criteria: ||1/(N + 1) * (u_avg - u)||
        hier_data_ops->linearSum(d_scratch_idx, 1.0 / (N + 1.0), u_idx, -1.0 / (N + 1.0), d_scratch_idx);
        if (sc_var) hier_math_ops.interp(d_dev_idx, d_dev_var, d_scratch_idx, sc_var, ghost_fill, time, false);
        if (nc_var) hier_math_ops.interp(d_dev_idx, d_dev_var, d_scratch_idx, nc_var, ghost_fill, time, false);
        if (ec_var) hier_math_ops.interp(d_dev_idx, d_dev_var, d_scratch_idx, ec_var, ghost_fill, time, false);
        if (fc_var) hier_math_ops.interp(d_dev_idx, d_dev_var, d_scratch_idx, fc_var, ghost_fill, time, false);
        if (cc_var) hier_data_ops->copyData(d_dev_idx, d_scratch_idx);
        d_visit_data_writer->writePlotData(hierarchy, d_visit_ts++, time);
        deallocate_patch_data(d_mean_idx, hierarchy, 0, hierarchy->getFinestLevelNumber());
        deallocate_patch_data(d_dev_idx, hierarchy, 0, hierarchy->getFinestLevelNumber());
    }
    else
    {
        // Determine our convergence criteria: ||1/(N + 1) * (u_avg - u)||
        hier_data_ops->linearSum(d_scratch_idx, 1.0 / (N + 1.0), u_idx, -1.0 / (N + 1.0), d_scratch_idx);
    }

    const double L2_norm = hier_data_ops->L2Norm(d_scratch_idx, wgt_idx);
    if (d_enable_logging)
    {
        // Print statistics to the log file.
        const double L1_norm = hier_data_ops->L1Norm(d_scratch_idx, wgt_idx);
        const double max_norm = hier_data_ops->maxNorm(d_scratch_idx, wgt_idx);
        plog << "At time " << time << ", the L^1 norm of the change in u_mean is " << L1_norm << "\n";
        plog << "At time " << time << ", the L^2 norm of the change in u_mean is " << L2_norm << "\n";
        plog << "At time " << time << ", the max norm of the change in u_mean is " << max_norm << "\n";
    }

    deallocate_patch_data(d_scratch_idx, hierarchy, 0, hierarchy->getFinestLevelNumber());
    if (L2_norm <= d_periodic_thresh)
    {
        d_idx_steady_state_map[time] = true;
        return true;
    }
    return false;
}

double
HierarchyAveragedDataManager::getTimePoint(double time, const double tol)
{
    time = map_to_period(d_period_start, d_period_end, time);
    auto it_up = d_snapshot_time_pts.upper_bound(time);
    double t_low, t_up;
    if (it_up == d_snapshot_time_pts.begin())
    {
        t_up = *it_up;
        // Set the "lower" bound to be the other end of the period.
        t_low = *d_snapshot_time_pts.rbegin();
    }
    else
    {
        auto it_low = std::next(it_up, -1);
        t_low = *(it_low);
        t_up = *(it_up);
    }

    if (std::abs(t_low - time) <= tol)
        return t_low;
    else if (std::abs(t_up - time) <= tol)
        return t_up;
    else
        TBOX_ERROR("Time point: " << time << " is not within the given tolerance " << tol << "!\n");
    return 0.0;
}

void
HierarchyAveragedDataManager::putToDatabase(Pointer<Database> db)
{
    db->putInteger("num_snapshots", d_snapshot_time_pts.size());
    int i = 0;
    for (const auto& time : d_snapshot_time_pts)
    {
        db->putDouble("time_" + std::to_string(i), time);
        db->putBool("at_steady_state_" + std::to_string(i), d_idx_steady_state_map.at(time));
        db->putInteger("num_updates_" + std::to_string(i), d_idx_num_updates_map.at(time));
    }

    d_snapshot_cache.putToDatabase(db);
    if (d_output_data) db->putInteger("visit_ts", d_visit_ts);
    db->putDouble("period_start", d_period_start);
    db->putDouble("period_end", d_period_end);
}

void
HierarchyAveragedDataManager::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
        db = restart_db->getDatabase(d_object_name);
    else
        TBOX_ERROR(d_object_name << ": Restart database corresponding to " << d_object_name
                                 << " not found in restart file.\n");

    // Now grab the data from the database
    int num_snapshots = db->getInteger("num_snapshots");
    for (int i = 0; i < num_snapshots; ++i)
    {
        double time = db->getDouble("time_" + std::to_string(i));
        d_snapshot_time_pts.insert(time);
        d_idx_steady_state_map.insert(std::make_pair(time, db->getBool("at_steady_state_" + std::to_string(i))));
        d_idx_num_updates_map.insert(std::make_pair(time, db->getInteger("num_updates_" + std::to_string(i))));
    }

    if (d_output_data) d_visit_ts = db->getInteger("visit_ts");
    d_period_start = db->getDouble("period_start");
    d_period_end = db->getDouble("period_end");
    d_period_length = d_period_end - d_period_start;
}
//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
