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

#ifndef included_IBAMR_SnapshotCache
#define included_IBAMR_SnapshotCache

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "HierarchyDataOpsReal.h"
#include "IntVector.h"
#include "Variable.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <map>
#include <set>
#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class SnapshotCache provides a method of storing snapshots. This class can also be used to determine
 * intermediate velocity fields when at a periodic steady state.
 */
class SnapshotCache
{
public:
    /*!
     * The constructor for class SnapshotCache sets some default values, reads in configuration information from input
     * databases, and registers some variable/context pairs with the VariableDatabase.
     */
    SnapshotCache(std::string object_name,
                  SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                  SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                  SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy);

    /*!
     * \brief The destructor for class SnapshotCache deallocates patch data as needed.
     */
    ~SnapshotCache();

    /*!
     * Clear the snapshots stored in this object. This deallocates all patch indices stored by this object.
     */
    void clearSnapshots();

    /*!
     * Create a snapshot. The data from this patch index is copied into an array of indices controlled by this
     * object. A mapping between the time point and the patch index is stored. This causes an error if more snapshots
     * are created than initially allocated. If you want to overwrite the snapshot, you should call updateSnapshot().
     * Data is not allocated in this class until this function is called.
     */
    void setSnapshot(int u_idx, double time, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * Update a snapshot. The data from this patch index is copied into an array of indices controlled by this
     * object. This outputs an error if a snapshot with that time point within the specified tolerance can not be found.
     */
    void updateSnapshot(int u_idx,
                        double time,
                        double tol,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * Transfer the snapshot at the time value within a given tolerance to the patch index u_idx on the supplied patch
     * hierarchy.
     *
     * Note this can cause a significant amount of communication if the supplied patch hierarchy has a different
     * configuration of patches than the snapshot patch hierarchy.
     */
    void getSnapshot(int u_idx,
                     double time,
                     double tol,
                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    SnapshotCache() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    SnapshotCache(const SnapshotCache& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    SnapshotCache& operator=(const SnapshotCache& that) = delete;

    std::string d_object_name;

    /*
     * Data for tracking snapshots.
     */
    unsigned int d_max_snapshots = 0, d_num_snapshots_stored = 0;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > > d_snapshot_vars;
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_ctx;
    std::string d_snapshot_coarsen_type = "CONSERVATIVE_COARSEN";
    std::string d_snapshot_refine_type = "CONSERVATIVE_LINEAR_REFINE";
    std::vector<int> d_snapshot_idxs;

    /*
     * Map between the snapshot patch index and the corresponding time point.
     */
    std::map<int, double> d_idx_time_map;

    /*!
     * Note we store a copy of the patch hierarchy at which the snapshot is stored. As needed, we will transfer data
     * from the snapshot hierarchy to the provided hierarchy.
     */
    std::map<int, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > > d_hierarchy_vec;

    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyDataOpsReal<NDIM, double> > d_hier_data_ops;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_SnapshotCache
