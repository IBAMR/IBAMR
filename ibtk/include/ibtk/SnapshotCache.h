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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_SnapshotCache
#define included_IBTK_SnapshotCache

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/ibtk_utilities.h>
#include <ibtk/private/HierarchyLevelsIdentity.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Serializable.h>

#include <CellVariable.h>
#include <EdgeVariable.h>
#include <FaceVariable.h>
#include <HierarchyDataOpsReal.h>
#include <IntVector.h>
#include <NodeVariable.h>
#include <PatchDataFactory.h>
#include <PatchHierarchy.h>
#include <RefineSchedule.h>
#include <SideVariable.h>

#include <map>
#include <set>
#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class SnapshotCache provides a method of storing snapshots of patch data on a snapshot of the patch hierarchy.
 *
 * This class can be used to store one type of variable, which is specified in the constructor. When retrieving or
 * storing data with this class, the same patch data layout must be used, e.g. you must consistently store/retrieve
 * CellData<NDIM, double> with depth of 2.
 *
 * This class stores snapshots in a list ordered by increasing time values. Iterators to the stored snapshots are
 * provided via the begin(), end(), and getSnapshot() methods.
 *
 * Communication schedules for updating and filling snapshots are reused until a level of the current hierarchy or of a
 * snapshot hierarchy is replaced. Replaced levels of the current hierarchy remain in memory until the next call to
 * updateSnapshot() or fillSnapshotOnHierarchy(), or to clearSnapshots().
 */
class SnapshotCache : public SAMRAI::tbox::Serializable
{
public:
    using value_type = std::pair<double, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>>;
    using iterator = std::vector<value_type>::iterator;
    using const_iterator = std::vector<value_type>::const_iterator;

    /*!
     * The constructor for class SnapshotCache sets some default values, reads in configuration information from input
     * databases, and registers a variable/context pair with the VariableDatabase. Can optionally set up this class for
     * restarts.
     *
     * The optional input database is searched for the following key:
     * - 'gcw' : Integer array that is the ghost cell width with which the internal variable is set up.
     *
     * If the input database is an invalid pointer, the ghost cell width is 1. Note that this value should be set as the
     * maximum ghost cell width needed to perform any needed refinement operation upon snapshot retrieval.
     */
    SnapshotCache(std::string object_name,
                  SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                  SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                  SAMRAI::tbox::Pointer<SAMRAI::hier::GridGeometry<NDIM>> grid_geom,
                  bool register_for_restart = true);

    /*!
     * \brief The destructor for class SnapshotCache deallocates patch data as needed.
     */
    ~SnapshotCache();

    /*!
     * Clear the snapshots stored in this object. This deallocates all patch indices stored by this object.
     */
    void clearSnapshots();

    /*!
     * Create a snapshot. This creates a snapshot of the current hierarchy, and copies that data from this patch index
     * into a patch index maintained by this object. If you want to overwrite a snapshot, you should call
     * getSnapshot() instead and copy the data manually. Data is not allocated in this class until this function is
     * called.
     */
    void storeSnapshot(int u_idx, double time, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy);

    /*!
     * Update the snapshot at the given time to the data in u_idx on the current hierarchy, and replace the snapshot
     * hierarchy by the current hierarchy. This outputs an error if no snapshot is stored within the tolerance. Index
     * u_idx may differ from the index used in storeSnapshot(), but its data layout must match that of the variable of
     * this object. Free function update_snapshot() calls this function.
     */
    void updateSnapshot(int u_idx,
                        double time,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> current_hierarchy,
                        double tol = 1.0e-8);

    /*!
     * Transfer the snapshot at the given time to u_idx on the current hierarchy, using the named refinement operator
     * and without synchronizing the levels. Index u_idx must have a ghost cell width sufficient for the operator. This
     * outputs an error if no snapshot is stored within the tolerance. Free function fill_snapshot_on_hierarchy() calls
     * this function.
     */
    void fillSnapshotOnHierarchy(int u_idx,
                                 double time,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> current_hierarchy,
                                 const std::string& snapshot_refine_type,
                                 double tol = 1.0e-8);

    /*!
     * Get a copy of the snapshot time and patch hierarchy pair. If a snapshot is not present within the provided
     * tolerance, a pair of NaN and nullptr is returned.
     */
    value_type getSnapshot(double time, double tol = 1.0e-8);

    /*!
     * Get a const iterator to the first element of the snapshots. Dereferencing the iterator returns a pair of the
     * snapshot time and patch hierarchy.
     *
     * \note This iterator is invalid if the number of snapshots changes.
     */
    //@{
    inline iterator begin()
    {
        return d_snapshots.begin();
    }
    //@}

    /*!
     * Get an iterator to element following the last element of the snapshots. Dereferencing this element results in
     * undefined behavior.
     *
     * \note This iterator is invalid if the number of snapshots changes.
     */
    //@{
    inline iterator end()
    {
        return d_snapshots.end();
    }
    //@}

    /*!
     * Get the number of stored snapshots.
     */
    inline size_t getNumSnapshots()
    {
        return d_snapshots.size();
    }

    /*!
     * Write the snapshot cache to a database. This writes both the time points and the snapshotted patch hierarchies.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    /*!
     * Get the patch index owned by this class.
     */
    inline int getPatchIndex()
    {
        return d_snapshot_idx;
    }

    /*!
     * Get the variable owned by this class.
     */
    inline SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> getVariable()
    {
        return d_snapshot_var;
    }

private:
    /*!
     * Read stored information from the restart files. Note the GridGeometry object is required in order to construct
     * the snapshotted patch hierarchies.
     */
    void getFromRestart(SAMRAI::tbox::Pointer<SAMRAI::hier::GridGeometry<NDIM>> grid_geom);

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

    /*!
     * Communication schedules of one snapshot hierarchy. The schedules contain the levels of the current hierarchy and
     * of the snapshot hierarchy, so they are dropped when either changes.
     */
    struct Schedules
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> snapshot_hierarchy;
        detail::HierarchyLevelsIdentity snapshot_levels;

        // One schedule per level that copies index update_idx on the current hierarchy to the snapshot levels.
        int update_idx = IBTK::invalid_index;
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchDataFactory<NDIM>> update_factory;
        std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM>>> update;

        // One set of schedules per level that transfer the snapshot to index fill_idx on the current hierarchy.
        struct Fill
        {
            int fill_idx = IBTK::invalid_index;
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchDataFactory<NDIM>> fill_factory;
            std::string refine_type;
            std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM>>> schedules;
        };
        std::vector<Fill> fill;
    };

    /*!
     * Drop all schedules if the current hierarchy differs from the one they were built on, and return the schedules of
     * the given snapshot hierarchy, which are dropped if the levels of that hierarchy changed.
     */
    Schedules& resetSchedules(const SAMRAI::hier::PatchHierarchy<NDIM>& current_hierarchy,
                              const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& snapshot_hierarchy);

    std::string d_object_name;

    /*
     * Data for tracking snapshots.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> d_snapshot_var;
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_ctx;
    int d_snapshot_idx = IBTK::invalid_index;
    SAMRAI::hier::IntVector<NDIM> d_gcw = 1;

    /*!
     * Vector of snapshots consisting of time points and hierarchies. Each hierarchy allocates data only with the patch
     * data index stored by this class.
     */
    std::vector<value_type> d_snapshots;

    /*!
     * Schedules for the snapshots, and the levels of the current hierarchy on which they are built.
     */
    //\{
    detail::HierarchyLevelsIdentity d_current_levels;
    std::vector<Schedules> d_schedules;
    //\}

    bool d_registered_for_restart = false;
};
} // namespace IBTK
//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_SnapshotCache
