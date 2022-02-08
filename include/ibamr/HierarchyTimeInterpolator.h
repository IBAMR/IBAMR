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

#ifndef included_IBAMR_HierarchyTimeInterpolator
#define included_IBAMR_HierarchyTimeInterpolator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/INSStaggeredHierarchyIntegrator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"

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
 * \brief Class HierarchyTimeInterpolator provides a method of storing flow snapshots to compute flow statistics
 * over periodic intervals. This class can also be used to determine intermediate velocity fields when at a periodic
 * steady state.
 */
class HierarchyTimeInterpolator
{
public:
    /*!
     * The constructor for class HierarchyTimeInterpolator sets some
     * default values, reads in configuration information from input databases, and registers some variable/context
     * pairs with the VariableDatabase.
     */
    HierarchyTimeInterpolator(std::string object_name,
                              SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                              SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy);

    /*!
     * \brief The destructor for class HierarchyTimeInterpolator deallocates patch data as needed.
     */
    ~HierarchyTimeInterpolator();

    /*!
     * Clear the snapshots stored in this object. This deallocates all patch indices stored by this object.
     */
    void clearSnapshots();

    /*!
     * Create a velocity snapshot. The data from this patch index is copied into an array of indices controlled by this
     * object. A mapping between the time point and the patch index is stored. This causes an error if more snapshots
     * are created than initially allocated. If you want to overwrite the snapshot, you should clear the data first.
     * Data is not allocated in this class until this function is called.
     */
    void
    setVelocitySnapshot(int u_idx, double time, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * Fill a velocity at a particular time point. We first find the two snapshots closest to the provided time values.
     * We then linearly interpolate in time the velocity and fill in the provided patch index. If no snapshots of the
     * velocity have occured, this emits an error. This assumes the flow is periodic, and that the period has been
     * correctly set up in the input file.
     */
    void
    fillVelocityAtTime(int u_idx, double time, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * Compute the time averaged velocity over the snapshots stored in this class. We use a trapezoidal rule to compute
     * the integral.
     */
    void computeTimeAvgVelocity(int avg_idx);

    /*!
     * Compute the turbulence kinetic energy over the snapshots stored in this class. We use a trapezoidal rule to
     * compute the integral.
     */
    void computeTKE(int tke_idx);

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

    /*!
     * Note we store one patch hierarchy, which assumes all velocity data is stored on the same hierarchy.
     * Realistically, each snapshot also comes with a snapshot of the hierarchy. One solution is to project each
     * velocity at each regrid step. This can be expensive and potentially do a lot of unecessary regrids. An
     * alternative is store snapshots of the patch hierarchy, and project the velocity as needed.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;

    /*
     * Data for tracking mean flow quantities.
     */
    unsigned int d_averaging_interval = 0, d_averaging_period = 0, d_num_averaging_cycles = 0, d_max_snapshots = 0,
                 d_num_snapshots_stored = 0;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > > d_snapshot_vars;
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_ctx;
    std::string d_snapshot_coarsen_type = "CONSERVATIVE_COARSEN";
    std::string d_snapshot_refine_type = "CONSERVATIVE_LINEAR_REFINE";
    std::vector<int> d_snapshot_idxs;

    /*
     * Interval to keep flow quantities
     */
    double d_t_start = std::numeric_limits<double>::quiet_NaN(), d_t_end = std::numeric_limits<double>::quiet_NaN(),
           d_t_period = std::numeric_limits<double>::quiet_NaN();

    /*
     * Map between time point and corresponding snapshot index.
     *
     * NOTE: this is potentially problematic, as due to finite precision we can have trouble finding keys.
     */
    std::map<double, int> d_idx_time_map;

    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM, double> > d_hier_sc_data_ops;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_HierarchyTimeInterpolator
