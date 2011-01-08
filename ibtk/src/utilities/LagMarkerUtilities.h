// Filename: LagMarkerUtilities.h
// Created on 28 Apr 2010 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_LagMarkerUtilities
#define included_LagMarkerUtilities

/////////////////////////////// INCLUDES /////////////////////////////////////

// C++ STDLIB INCLUDES
#include <vector>

// PETSC INCLUDES
#include <petsc.h>

// SAMRAI INCLUDES
#include <PatchHierarchy.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LagMarkerUtilities is a utility class that defines useful
 * functions for dealing with Lagrangian marker particles.
 */
class LagMarkerUtilities
{
public:

    /*!
     * Read the initial positions of the markers from a file.  Returns the
     * number of markers read.
     */
    static int
    readMarkerPositions(
        std::vector<double>& mark_init_posns,
        const std::string& mark_input_file_name,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * Advect all markers by the specified advection velocity.  Uses a
     * single-step explicit midpoint rule.
     */
    static void
    advectMarkers(
        const int mark_current_idx,
        const int mark_new_idx,
        const int u_idx,
        const double dt,
        const std::string& weighting_fcn,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int coarsest_ln=-1,
        const int finest_ln=-1);

    /*!
     * Collect all marker data onto the coarsest level of the patch hierarchy
     * (to prepare for regridding the patch hierarchy).
     */
    static void
    collectMarkersOnPatchHierarchy(
        const int mark_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * Initialize marker data on the specified level of the patch hierarchy by
     * refining markers from coarser levels in the patch hierarchy or by copying
     * markers from old_level.
     */
    static void
    initializeMarkersOnLevel(
        const int mark_idx,
        const std::vector<double>& mark_init_posns,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const bool initial_time,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level);

    /*!
     * Prune marker data in refined regions of the specified levels of the patch
     * hierarchy.
     */
    static void
    pruneDuplicateMarkers(
        const int mark_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int coarsest_ln=-1,
        const int finest_ln=-1);

    /*!
     * Count the markers.
     */
    static int
    countMarkers(
        const int mark_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int coarsest_ln=-1,
        const int finest_ln=-1);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     */
    LagMarkerUtilities();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    LagMarkerUtilities(
        const LagMarkerUtilities& from);

    /*!
     * \brief Unimplemented destructor.
     */
    ~LagMarkerUtilities();

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LagMarkerUtilities& operator=(
        const LagMarkerUtilities& that);
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/LagMarkerUtilities.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LagMarkerUtilities
