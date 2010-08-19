// Filename: AdvDiffHypPatchOps.h
// Created on 19 Mar 2004 by Boyce Griffith
//
// Copyright (c) 2002-2010 Boyce Griffith
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef included_AdvDiffHypPatchOps
#define included_AdvDiffHypPatchOps

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/AdvectHypPatchOps.h>
#include <ibamr/GodunovAdvector.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <PatchLevel.h>
#include <tbox/Database.h>
#include <tbox/Pointer.h>

// C++ STDLIB INCLUDES
#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class AdvDiffHypPatchOps is a version of class AdvectHypPatchOps that
 * has been specialized for use with a linearly-implicit time integrator for the
 * advection-diffusion equation.
 *
 * \see AdvDiffHierarchyIntegrator
 * \see AdvectHypPatchOps
 */
class AdvDiffHypPatchOps
    : public AdvectHypPatchOps
{
public:
    /*!
     * The constructor for AdvDiffHypPatchOps sets default parameters for the
     * patch strategy.  The constructor also registers this object for restart
     * with the restart manager using the object name when so requested.
     *
     * After default values are set, this routine calls getFromRestart() if
     * execution from a restart file is specified.  Finally, getFromInput() is
     * called to read values from the given input database (potentially
     * overriding those found in the restart file).
     */
    AdvDiffHypPatchOps(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        SAMRAI::tbox::Pointer<GodunovAdvector> godunov_advector,
        SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom,
        bool register_for_restart=true);

    /*!
     * The destructor for AdvDiffHypPatchOps unregisters the patch strategy
     * object with the restart manager when so registered.
     */
    virtual
    ~AdvDiffHypPatchOps();

    /*!
     * Return the name of the patch operations object.
     */
    const std::string&
    getName() const;

    ///
    ///  The following routines:
    ///
    ///      conservativeDifferenceOnPatch(),
    ///      preprocessAdvanceLevelState(),
    ///      postprocessAdvanceLevelState()
    ///
    ///  are redefined from the AdvectHypPatchOps base class.
    ///

    /*!
     * Update solution variables by performing a conservative difference using
     * the fluxes calculated in computeFluxesOnPatch().
     */
    virtual void
    conservativeDifferenceOnPatch(
        SAMRAI::hier::Patch<NDIM>& patch,
        const double time,
        const double dt,
        bool at_synchronization);

#if 0
    // NOTE: Should the following method be removed?

    /*!
     * Compute the values of any time-dependent source terms for use by the
     * explicit predictor.
     *
     * This routine is called after patch boundary data is filled (i.e., ghosts)
     * and before computeFluxesOnPatch().
     *
     * Note that when this routine is called, the scratch data is filled on all
     * patches (i.e., ghost cells) and that data is the same as the current
     * level data on all patch interiors.  That is, both scratch and current
     * data correspond to current_time.
     */
    virtual void
    preprocessAdvanceLevelState(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >& level,
        double current_time,
        double dt,
        bool first_step,
        bool last_step,
        bool regrid_advance);
#endif

    /*!
     * Add source terms to the updated solution.
     *
     * This routine is called after conservativeDifferenceOnPatch() is called
     * and before computeStableDtOnPatch().
     *
     * Note that when this routine is called, the scratch data is filled on all
     * patches (i.e., ghost cells) and that data is the same as the new level
     * data on all patch interiors.  That is, both scratch and new data
     * correspond to current_time + dt on patch interiors.  The current data and
     * ghost values correspond to the current_time.
     */
    virtual void
    postprocessAdvanceLevelState(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >& level,
        double current_time,
        double dt,
        bool first_step,
        bool last_step,
        bool regrid_advance);


private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    AdvDiffHypPatchOps();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AdvDiffHypPatchOps(
        const AdvDiffHypPatchOps& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AdvDiffHypPatchOps&
    operator=(
        const AdvDiffHypPatchOps& that);
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/AdvDiffHypPatchOps.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_AdvDiffHypPatchOps
