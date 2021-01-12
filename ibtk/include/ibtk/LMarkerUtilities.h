// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#ifndef included_IBTK_LMarkerUtilities
#define included_IBTK_LMarkerUtilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/LMarkerSetData.h"
#include "ibtk/ibtk_utilities.h"

#include "tbox/Pointer.h"

#include <string>
#include <vector>

namespace SAMRAI
{
namespace geom
{
template <int DIM>
class CartesianGridGeometry;
} // namespace geom
namespace hier
{
template <int DIM>
class BasePatchLevel;
template <int DIM>
class PatchHierarchy;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LMarkerUtilities is a utility class that defines useful
 * functions for dealing with Lagrangian marker particles.
 */
class LMarkerUtilities
{
public:
    /*!
     * Read the initial positions of a collection of Lagrangian markers from a
     * text file.  Returns the number of markers read.
     */
    static unsigned int
    readMarkerPositions(std::vector<Point>& mark_init_posns,
                        const std::string& mark_input_file_name,
                        SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom);

    /*!
     * Advect all markers by the specified advection velocity using forward
     * Euler.
     */
    static void eulerStep(int mark_current_idx,
                          int mark_new_idx,
                          int u_current_idx,
                          double dt,
                          const std::string& weighting_fcn,
                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                          int coarsest_ln = -1,
                          int finest_ln = -1);

    /*!
     * Advect all markers by the specified advection velocity using the explicit
     * midpoint rule.
     *
     * \note This function requires an initial call to eulerStep to compute the
     * predicted marker positions.
     */
    static void midpointStep(int mark_current_idx,
                             int mark_new_idx,
                             int u_half_idx,
                             double dt,
                             const std::string& weighting_fcn,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                             int coarsest_ln = -1,
                             int finest_ln = -1);

    /*!
     * Advect all markers by the specified advection velocity using the explicit
     * trapezoidal rule.
     *
     * \note This function requires an initial call to eulerStep to compute the
     * current marker velocities and the predicted marker positions.
     */
    static void trapezoidalStep(int mark_current_idx,
                                int mark_new_idx,
                                int u_new_idx,
                                double dt,
                                const std::string& weighting_fcn,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                int coarsest_ln = -1,
                                int finest_ln = -1);

    /*!
     * Collect all marker data onto the coarsest level of the patch hierarchy
     * (to prepare for regridding the patch hierarchy).
     */
    static void collectMarkersOnPatchHierarchy(int mark_idx,
                                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * Initialize marker data on the specified level of the patch hierarchy by
     * refining markers from coarser levels in the patch hierarchy or by copying
     * markers from old_level.
     */
    static void initializeMarkersOnLevel(int mark_idx,
                                         const std::vector<Point>& mark_init_posns,
                                         SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                         int level_number,
                                         bool initial_time,
                                         SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level);

    /*!
     * Prune marker data in refined regions of the specified levels of the patch
     * hierarchy.
     */
    static void pruneInvalidMarkers(int mark_idx,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                    int coarsest_ln = -1,
                                    int finest_ln = -1);

    /*!
     * Count the markers.
     */
    static unsigned int countMarkers(int mark_idx,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                     int coarsest_ln = -1,
                                     int finest_ln = -1);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     */
    LMarkerUtilities() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    LMarkerUtilities(const LMarkerUtilities& from) = delete;

    /*!
     * \brief Unimplemented destructor.
     */
    ~LMarkerUtilities() = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LMarkerUtilities& operator=(const LMarkerUtilities& that) = delete;

    /*!
     * Determine the number of markers in a patch.
     */
    static unsigned int countMarkersOnPatch(SAMRAI::tbox::Pointer<LMarkerSetData> mark_data);

    /*!
     * Collect marker positions into a single vector.
     */
    static void collectMarkerPositionsOnPatch(std::vector<double>& X_mark,
                                              SAMRAI::tbox::Pointer<LMarkerSetData> mark_data);

    /*!
     * Reset marker positions from a single vector.
     */
    static void resetMarkerPositionsOnPatch(const std::vector<double>& X_mark,
                                            SAMRAI::tbox::Pointer<LMarkerSetData> mark_data);

    /*!
     * Collect marker velocities into a single vector.
     */
    static void collectMarkerVelocitiesOnPatch(std::vector<double>& U_mark,
                                               SAMRAI::tbox::Pointer<LMarkerSetData> mark_data);

    /*!
     * Reset marker velocities from a single vector.
     */
    static void resetMarkerVelocitiesOnPatch(const std::vector<double>& U_mark,
                                             SAMRAI::tbox::Pointer<LMarkerSetData> mark_data);

    /*!
     * Prevent markers from leaving the computational domain through physical
     * boundaries.
     */
    static void preventMarkerEscape(std::vector<double>& X_mark,
                                    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LMarkerUtilities
