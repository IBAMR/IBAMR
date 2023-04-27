// ---------------------------------------------------------------------
//
// Copyright (c) 2023 - 2023 by the IBAMR developers
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

#include <ibtk/ibtk_utilities.h>

#include <tbox/Database.h>

#include <Box.h>
#include <CartesianGridGeometry.h>
#include <PatchHierarchy.h>

#include <deque>
#include <iosfwd>
#include <string>
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
     *
     * @param[in] patch_box Box, defined over the patch's level's index space,
     *            for the present Patch.
     * @param[in] nonoverlapping_patch_boxes The subset of @p patch_box which is
     *            uniquely indexed by this Patch.
     * @param[in] grid_geom The grid geometry.
     * @param[in] ratio The ratio of the Path's level to the coarsest level.
     */
    MarkerPatch(const SAMRAI::hier::Box<NDIM>& patch_box,
                const std::vector<SAMRAI::hier::Box<NDIM> >& nonoverlapping_patch_boxes,
                const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> >& grid_geom,
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
    std::tuple<std::vector<int>, EigenAlignedVector<IBTK::Point>, EigenAlignedVector<IBTK::Vector> > prune();

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
    std::vector<SAMRAI::hier::Box<NDIM> > d_nonoverlapping_patch_boxes;

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
 * points are akin to IB points except they are simply advected by the fluid
 * and do not apply any force upon it.
 */
class MarkerPatchHierarchy : public SAMRAI::tbox::DescribedClass
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
                         SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                         const EigenAlignedVector<IBTK::Point>& positions,
                         const EigenAlignedVector<IBTK::Point>& velocities);

    void reinit(const EigenAlignedVector<IBTK::Point>& positions, const EigenAlignedVector<IBTK::Point>& velocities);

    /**
     * Get the MarkerPatch associated with the present level and local patch number.
     */
    const MarkerPatch& getMarkerPatch(const int ln, const int local_patch_num) const;

    /**
     * Get the total number of marker points stored across all processors.
     */
    std::size_t getNumberOfMarkers() const;

#if 0
    /**
     * Save the present state of the object to a SAMRAI database.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;
#endif

    /**
     * Set marker velocities with some given velocity field @p u_idx and
     * interaction kernel @p kernel.
     */
    void setVelocities(const int u_idx, const std::string& kernel);

    /**
     * Advect the markers by some given velocity field at half and new times
     * and the interaction kernel @p kernel.
     *
     * This time integrator is the 'explicit midpoint' method, which uses the
     * present velocity to predict a midpoint position and then uses the
     * midpoint velocity velocity to compute the new position.
     *
     * @note This function assumes that @p u_half_idx and @p u_new_idx have
     * up-to-date ghost values with sufficient width for @p kernel. See
     * LEInteractor for kernel names.
     */
    void midpointStep(const double dt, const int u_half_idx, const int u_new_idx, const std::string& kernel);

protected:
    /**
     * Redistribute any particles which may have moved outside of their
     * patches, taking into account physical boundaries (markers are not
     * allowed to escape) and periodicity.
     */
    void pruneAndRedistribute();

    std::string d_object_name;

    std::size_t d_num_markers;

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;

    std::vector<std::deque<MarkerPatch> > d_marker_patches;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_MarkerPatchHierarchy
