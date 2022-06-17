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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_snapshot_utilities
#define included_IBTK_snapshot_utilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/SnapshotCache.h>

#include <tbox/Pointer.h>

#include <HierarchyDataOpsReal.h>
#include <PatchHierarchy.h>

#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * Update a snapshot at the given time. This function copies the data stored in u_idx into the snapshot index in the
 * cache object. The snapshotted hierarchy is also updated to the current hierarchy.
 *
 * This outputs an error if a snapshot with that time point within the specified tolerance can not be found.
 *
 * Calls to this function may use a different patch index than the one used in setSnapshot(), but the underlying
 * data layout must be consistent. This data layout is specified by the variable used to construct the cache.
 */
void update_snapshot(SnapshotCache& cache,
                     int u_idx,
                     double time,
                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > current_hierarchy,
                     double tol = 1.0e-8);

/*!
 * Transfer the snapshot at the specified time value within a given tolerance to the patch index u_idx on the
 * supplied patch hierarchy. The patch index u_idx must contain sufficient ghost cell width to perform the
 * operations used by the refinement operator.
 *
 * Note this can require a significant amount of communication if the supplied patch hierarchy has a different
 * configuration of patches than the snapshot patch hierarchy.
 *
 * This function does not synchronize the data on the current hierarchy (i.e. no coarsening is performed). It
 * transfers the snapshot using the refinement operator.
 */
void fill_snapshot_on_hierarchy(SnapshotCache& cache,
                                int u_idx,
                                double time,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > current_hierarchy,
                                const std::string& snapshot_refine_type,
                                double tol = 1.0e-8);

/*!
 * This function interpolates in time between two snapshots. We find the two closest time points stored in the snapshot
 * cache and linear interpolate in time between them. Each snapshot will be transferred onto the provided hierarchy. A
 * scratch patch index with the same patch layout as u_idx must be provided, and data for both must be already
 * allocated. The refine_type must be a valid refinement operation for the data layout.
 *
 * This assumes that the time is between two snapshots stored in the cache. An error will occur if two time points can
 * not be found. If there is a single snapshot, this function returns that snapshot.
 *
 * The hier_data_ops object, if provided, must match the same variable type used in for the snapshots.
 *
 * If the optional period parameter is provided, this function will treat the first snapshot time point t_1 as also the
 * last snapshot time point t_1 + period.
 */
void
fill_snapshot_at_time(SnapshotCache& snapshot_cache,
                      int u_idx,
                      double t,
                      int scr_idx,
                      SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                      const std::string& refine_type,
                      SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyDataOpsReal<NDIM, double> > hier_data_ops = nullptr,
                      double period = std::numeric_limits<double>::quiet_NaN());

} // namespace IBTK
#endif
