// ---------------------------------------------------------------------
//
// Copyright (c) 2023 - 2026 by the IBAMR developers
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

#ifndef included_IBTK_MarkerPatchHierarchy
#define included_IBTK_MarkerPatchHierarchy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/IBTK_MPI.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/private/BoxArrayIndex.h>
#include <ibtk/private/HierarchyLayoutKey.h>

#include <tbox/Database.h>

#include <Box.h>
#include <CartesianGridGeometry.h>
#include <PatchHierarchy.h>

#include <array>
#include <deque>
#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
class MarkerPatchHierarchy;

/*!
 * Structure containing all relevant data to store markers on a Patch.
 */
class MarkerPatch
{
public:
    /*!
     * Constructor.
     *
     * @param[in] patch_box Box, defined over the patch's level's index space,
     *            for the present Patch.
     * @param[in] nonoverlapping_patch_boxes The subset of @p patch_box which is
     *            uniquely indexed by this Patch.
     * @param[in] grid_geom The grid geometry.
     * @param[in] ratio The ratio of the Path's level to the coarsest level.
     */
    MarkerPatch(const SAMRAI::hier::Box<NDIM>& patch_box,
                const std::vector<SAMRAI::hier::Box<NDIM>>& nonoverlapping_patch_boxes,
                const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM>>& grid_geom,
                const SAMRAI::hier::IntVector<NDIM>& ratio);

    /**
     * Add a marker point.
     */
    void insert(const int& index, const IBTK::Point& position, const IBTK::Vector& velocity);

    /*!
     * Return whether or not the point should be uniquely assigned to this Patch.
     */
    bool contains(const IBTK::Point& position) const;

    /**
     * Remove all marker points which do not lie in cells uniquely owned by
     * this patch.
     */
    std::tuple<std::vector<int>, std::vector<IBTK::Point>, std::vector<IBTK::Vector>> prune();

    /*!
     * Return the @p index-th marker point stored by the present Patch.
     */
    std::tuple<int, IBTK::Point, IBTK::Vector> operator[](const unsigned int index) const;

    /*!
     * Return the number of marker points presently stored by this object.
     *
     * @note This may contain points outside the cells uniquely owned by this
     * Patch if prune() is not called first.
     */
    std::size_t size() const;

private:
    // Marker data. Stored in a packed format so we can directly call
    // Fortran routines on it.
    std::vector<int> d_indices;
    std::vector<double> d_positions;
    std::vector<double> d_velocities;

    // Box for the patch itself.
    SAMRAI::hier::Box<NDIM> d_patch_box;

    // The unique subset of index space assigned to this Patch.
    std::vector<SAMRAI::hier::Box<NDIM>> d_nonoverlapping_patch_boxes;

    // Data for computing cell indices. Passed along to IndexUtilities::getCellIndex();
    SAMRAI::hier::Box<NDIM> d_domain_box;
    std::array<double, NDIM> d_x_lo;
    std::array<double, NDIM> d_x_up;
    std::array<double, NDIM> d_dx;

    // Allow MarkerPatchHierarchy to directly manipulate the arrays owned by this
    // class for timestepping.
    friend class MarkerPatchHierarchy;
};

/*!
 * Implementation of marker points, also sometimes called particles. These
 * points are akin to IB points except they are advected by the fluid
 * and do not apply any force upon it.
 *
 * Marker velocities are independent of the patch layout only if the caller
 * fills the ghost cells of the velocity with the refine operator alone,
 * without coarse-fine interpolation.
 *
 * Markers are indexed in the index space of the bounding box of the physical
 * domain. A marker that is outside the bounding box or in a notch of a
 * physical domain of several boxes (inside the bounding box, but in no box of
 * the domain) is outside the domain: it is stored with zero velocity and is
 * not advected.
 *
 * A marker stored on a coarse level that moves beyond a wall or into a notch
 * next to a finer patch can lie several cells of the finer level from the box
 * of the patch that supplies its velocity. The interpolation reaches one cell
 * of that patch beyond its box. For a marker farther away, the program stops
 * in a debug build, and the marker has zero velocity in an optimized build.
 *
 * Markers are redistributed by exchanging them only with the processes whose
 * patches are near the patches of a process. The exchange requires that, in
 * each time step and at each intermediate stage position of the step, a marker
 * moves at most one cell of the level that stores it. If a marker moves to a
 * patch of a process that is not a neighbor, the program stops in a debug
 * build. In an optimized build, the marker has zero velocity at that stage of
 * the step, or it is lost when markers are redistributed.
 *
 * Each process stores the markers that leave the domain from its own patches.
 * Call reinit() after every regrid of the patch hierarchy: the program stops
 * at a time step or a call to setVelocities() if the layout of the hierarchy
 * has changed since the last call to reinit().
 *
 * The constructor, the destructor, reinit(), getFromDatabase(), the
 * time-stepping functions, collectAllMarkers(), and writeH5Part() are
 * collective. The object duplicates the MPI communicator, so it cannot be
 * copied.
 */
class MarkerPatchHierarchy : public SAMRAI::tbox::Serializable
{
public:
    /**
     * Constructor.
     *
     * @p positions and @p velocities should be the positions and velocities
     * of the markers. These vectors should contain the complete (i.e., the
     * same on each processor) list of markers and are implicitly numbered by
     * their array index.
     */
    MarkerPatchHierarchy(const std::string& name,
                         SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> patch_hierarchy,
                         const std::vector<IBTK::Point>& positions,
                         const std::vector<IBTK::Point>& velocities,
                         const bool register_for_restart = true);

    virtual ~MarkerPatchHierarchy();

    MarkerPatchHierarchy(const MarkerPatchHierarchy&) = delete;

    MarkerPatchHierarchy& operator=(const MarkerPatchHierarchy&) = delete;

    /*!
     * Replace the markers by the given markers and rebuild the data that
     * depend on the layout of the patch hierarchy.
     */
    void reinit(const std::vector<IBTK::Point>& positions, const std::vector<IBTK::Point>& velocities);

    /**
     * Save the present state of the object (i.e., the marker point data) to a
     * SAMRAI database.
     *
     * @note Since the constructor requires a pointer to the PatchHierarchy this
     * function relies on the constructor to set up all geometric and index
     * space data and only saves the markers themselves.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    /*!
     * Write the particles to a single H5Part file.
     *
     * H5Part is a simple HDF5 schema for particles which VisIt understands. You
     * should either use the file suffix .h5part or use the VisIt H5Part reader
     * explicitly. To visualize, select 'Mesh -> particles'.
     */
    void writeH5Part(const std::string& filename,
                     const int time_step,
                     const double simulation_time,
                     const bool write_velocities = true) const;

    /**
     * Load the marker points from a database.
     */
    virtual void getFromDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /**
     * Get the MarkerPatch associated with the present level and local patch number.
     */
    const MarkerPatch& getMarkerPatch(const int ln, const int local_patch_num) const;

    /**
     * Get the total number of marker points stored across all processors.
     */
    std::size_t getNumberOfMarkers() const;

    /**
     * Collect all markers on all processors in a single array.
     */
    std::pair<std::vector<IBTK::Point>, std::vector<IBTK::Vector>> collectAllMarkers() const;

    /**
     * Set marker velocities with some given velocity field @p u_idx and
     * interaction kernel @p kernel.
     */
    void setVelocities(const int u_idx, const std::string& kernel);

    /**
     * Advect the markers with their present velocity values with the forward
     * Euler method. After the markers are moved their velocities are reset to
     * values interpolated from @p u_new_idx via the IB kernel @p kernel.
     *
     * @note This function assumes that @p u_new_idx has up-to-date ghost
     * values with sufficient width for @p kernel. See LEInteractor for kernel
     * names.
     */
    void forwardEulerStep(const double dt, const int u_new_idx, const std::string& kernel);

    /**
     * Advect the markers with their present velocity values with the explicit
     * backward Euler method (i.e., a forward Euler prediction step followed by
     * a backward Euler correction step). After the markers are moved their
     * velocities are reset to values interpolated from @p u_new_idx via the IB
     * kernel @p kernel.
     *
     * @note This function assumes that @p u_new_idx has up-to-date ghost
     * values with sufficient width for @p kernel. See LEInteractor for kernel
     * names.
     */
    void backwardEulerStep(const double dt, const int u_new_idx, const std::string& kernel);

    /**
     * Advect the markers by some given velocity field at half and new times
     * and the interaction kernel @p kernel.
     *
     * This time integrator is the 'explicit midpoint' method, which uses the
     * present velocity to predict a midpoint position and then uses the
     * midpoint velocity to compute the new position.
     *
     * @note This function assumes that @p u_half_idx and @p u_new_idx have
     * up-to-date ghost values with sufficient width for @p kernel. See
     * LEInteractor for kernel names.
     */
    void midpointStep(const double dt, const int u_half_idx, const int u_new_idx, const std::string& kernel);

    /**
     * Advect the markers by some given velocity field at the end of a timestep
     * and the interaction kernel @p kernel.
     *
     * This time integrator is the 'explicit trapezoid' method, which uses the
     * present velocity to predict a new position and then uses the
     * trapezoid rule to correct that position.
     *
     * @note This function assumes that @p u_new_idx has up-to-date ghost values
     * with sufficient width for @p kernel. See LEInteractor for kernel names.
     */
    void trapezoidalStep(const double dt, const int u_new_idx, const std::string& kernel);

protected:
    /**
     * Return the marker patches of the patches that this process owns, from the
     * finest level to the coarsest, with their patches.
     */
    std::vector<std::pair<MarkerPatch*, SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>>>> getLocalMarkerPatches();

    /**
     * Interpolate the velocity @p u_idx at @p positions, the positions of the markers of @p local_marker_patches in
     * order, into @p velocities, with the patch whose box contains each position (beyond a physical boundary, or in a
     * notch of the domain, the patch that contains the nearest point of the domain). A marker must move less than one
     * cell width from the patch that stores it.
     */
    void interpolateInCurrentPatches(
        int u_idx,
        const std::string& kernel,
        const std::vector<std::pair<MarkerPatch*, SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>>>>&
            local_marker_patches,
        const std::vector<std::vector<double>>& positions,
        std::vector<std::vector<double>>& velocities);

    /**
     * Redistribute any particles which may have moved outside of their
     * patches, taking into account physical boundaries (markers are not
     * allowed to escape) and periodicity.
     */
    void pruneAndRedistribute();

    std::string d_object_name;

    bool d_register_for_restart;

    std::size_t d_num_markers;

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;

    std::vector<std::deque<MarkerPatch>> d_marker_patches;

    MarkerPatch d_markers_outside_domain;

private:
    /*!
     * Owner of a cell: the level, the rank, and the patch number on that level.
     */
    struct Owner
    {
        int level = -1;
        int rank = -1;
        int patch = -1;
    };

    /*!
     * The cells that each patch of a level owns, for all patches of the level.
     */
    struct LevelOwnership
    {
        // The parts of the patch boxes that no finer level covers, and the rank and patch number of each box.
        detail::BoxArrayIndex boxes;
        std::vector<int> box_rank;
        std::vector<int> box_patch;

        // The boxes of each patch, by patch number.
        std::vector<std::vector<SAMRAI::hier::Box<NDIM>>> patch_boxes;

        // The number of the marker patch of each patch of this process, by patch number, or -1.
        std::vector<int> local_patch;

        // The index frame of the level, as in MarkerPatch.
        SAMRAI::hier::Box<NDIM> domain_box;
        std::array<double, NDIM> dx;
    };

    /*!
     * Build the ownership of cells and the neighbors of this process for the current layout.
     */
    void buildExchangeData();

    /*!
     * Return the owner of the cell that contains X, or an owner of level -1 if no patch owns it. X is used in the
     * index frame of each level. @p found is scratch space.
     */
    Owner findOwner(const double* X, std::vector<int>& found) const;

    /*!
     * Stop with an error if the layout of the patch hierarchy has changed since the last call to reinit().
     */
    void checkHierarchyIsCurrent(const char* function) const;

    // The layout of the levels that the data below were built for.
    detail::HierarchyLayoutKey d_layout_key;

    // Data built for the layout: the ownership of each level and the bounds of the domain.
    std::vector<LevelOwnership> d_ownership;
    std::array<double, NDIM> d_x_lo, d_x_up;

    // The ranks other than this one whose patches are within two cells of a patch of this process, in ascending order.
    std::vector<int> d_neighbor_ranks;

    // A communicator of this object, so that its messages cannot match other receives.
    MPI_Comm d_comm = MPI_COMM_NULL;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_MarkerPatchHierarchy
