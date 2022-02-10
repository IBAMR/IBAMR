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

#ifndef included_IBAMR_INSStaggeredMeanFlowCalculator
#define included_IBAMR_INSStaggeredMeanFlowCalculator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/SnapshotCache.h>

#include "CellVariable.h"
#include "IntVector.h"
#include "SideVariable.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class HierarchyTimeInterpolator provides a method of determining and storing average fields over periodic
 * intervals. This class can also be used to determine intermediate velocity fields when at a periodic steady state.
 */
class HierarchyTimeInterpolator
{
public:
    /*!
     * The constructor for class HierarchyTimeInterpolator sets some default values and determines data centering. The
     * expected period and number of snapshots must be set in the input database. This class assumes the snapshots will
     * be taken at equidistant points along the periodic interval.
     */
    HierarchyTimeInterpolator(std::string object_name,
                              SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                              SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                              SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy);

    /*!
     * The constructor for class HierarchyTimeInterpolator sets some default values and determines data centering. In
     * this constructor, the times at which snapshots are taken are set, as well as the period length.
     */
    HierarchyTimeInterpolator(std::string object_name,
                              std::set<double> snapshot_time_points,
                              double period,
                              SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                              SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy);

    /*!
     * \brief The destructor for class HierarchyTimeInterpolator deallocates patch data as needed.
     */
    ~HierarchyTimeInterpolator();

    /*!
     * Clear the snapshots taken with this class.
     */
    void clearSnapshots();

    /*!
     * Fill a patch index at a particular time point. This will perform a time interpolation of the snapshots stored in
     * the SnapshotCache object. Note if there is only one snapshot (i.e. the period is 0 and we are assumed to be at a
     * steady state), a copy will be performed
     */
    void
    fillSnapshotAtTime(int u_idx, double time, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * Update the time averaged snapshot. The update is computed through the relation
     * ``` u_avg = u_avg + 1/N * (u - u_avg)```
     * in which N is the number of updates to the snapshot. If this the first update after clearMean() is called, we
     * simply copy data.
     *
     * This function returns true if the snapshot is at a steady state, or that ||u - u_avg|| < threshold.
     */
    bool updateTimeAveragedSnapshot(int u_idx,
                                    double time,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * Return true if all the tracked mean fields are at a steady state.
     */
    bool isAtPeriodicSteadyState();

    /*!
     * Return the time points this class is storing.
     */
    const std::set<double>& getSnapshotTimePts()
    {
        return d_snapshot_time_pts;
    }

    /*!
     * Returns true if a flow snapshot is being stored within the specified tolerance.
     */
    bool timePtIsStored(double time, double tol);

    /*!
     * Returns the index of the snapshot, if that time point is stored. Otherwise returns -1 (the largest unsigned
     * integer).
     */
    unsigned int getTimePtIndex(double time, double tol);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    HierarchyTimeInterpolator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    HierarchyTimeInterpolator(const HierarchyTimeInterpolator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    HierarchyTimeInterpolator& operator=(const HierarchyTimeInterpolator& that) = delete;

    std::string d_object_name;

    /*
     * Data for tracking mean flow quantities.
     */
    unsigned int d_num_averaging_cycles = 0;
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_scratch_var;
    int d_scratch_idx;
    std::string d_u_mean_coarsen_type = "CONSERVATIVE_COARSEN";
    std::string d_u_mean_refine_type = "CONSERVATIVE_LINEAR_REFINE";

    /*
     * Interval to keep flow quantities
     */
    double d_t_start = std::numeric_limits<double>::quiet_NaN(), d_t_end = std::numeric_limits<double>::quiet_NaN(),
           d_t_period = std::numeric_limits<double>::quiet_NaN(),
           d_periodic_thresh = std::numeric_limits<double>::quiet_NaN();

    /*
     * Map between time point index and the number of updates. We use this value to update average snapshots.
     */
    std::map<int, unsigned int> d_idx_num_updates_map;

    /*
     * Set of time points that snapshots are determined. These values are between d_t_start and d_t_end.
     */
    std::set<double> d_snapshot_time_pts;

    std::unique_ptr<SnapshotCache> d_hier_time_interpolator;
};

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_HierarchyTimeInterpolator
