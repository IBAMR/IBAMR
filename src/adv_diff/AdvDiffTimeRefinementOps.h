#ifndef included_AdvDiffTimeRefinementOps
#define included_AdvDiffTimeRefinementOps

// Filename: AdvDiffTimeRefinementOps.h
// Last modified: <06.Dec.2006 17:55:49 griffith@box221.cims.nyu.edu>
// Created on 06 Dec 2006 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/GodunovAdvector.h>
#include <ibamr/GodunovHypPatchOps.h>

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
 * \brief A specialized version of the GodunovHypPatchOps class for
 * use with an implicit, time-refined time integrator for the
 * advection-diffusion equation.
 *
 * \see GodunovHypPatchOps
 */
class AdvDiffTimeRefinementOps
    : public GodunovHypPatchOps
{
public:
    /*!
     * The constructor for AdvDiffTimeRefinementOps sets default
     * parameters for the patch strategy.  The constructor also
     * registers this object for restart with the restart manager
     * using the object name when so requested.
     *
     * After default values are set, this routine calls
     * getFromRestart() if execution from a restart file is specified.
     * Finally, getFromInput() is called to read values from the given
     * input database (potentially overriding those found in the
     * restart file).
     */
    AdvDiffTimeRefinementOps(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        SAMRAI::tbox::Pointer<GodunovAdvector> godunov_advector,
        SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom,
        bool register_for_restart=true);

    /*!
     * The destructor for AdvDiffTimeRefinementOps unregisters the patch
     * strategy object with the restart manager when so registered.
     */
    virtual ~AdvDiffTimeRefinementOps();

    ///
    ///  The following routines:
    ///
    ///      conservativeDifferenceOnPatch(),
    ///      preprocessAdvanceLevelState(),
    ///      postprocessAdvanceLevelState()
    ///
    ///  are redefined from the GodunovHypPatchOps base class.
    ///

    /*!
     * Update solution variables by performing a conservative
     * difference using the fluxes calculated in
     * computeFluxesOnPatch().
     */
    virtual void conservativeDifferenceOnPatch(
        SAMRAI::hier::Patch<NDIM>& patch,
        const double time,
        const double dt,
        bool at_synchronization);

    /*!
     * Compute the values of any time-dependent source terms for use
     * by the explicit predictor.
     *
     * This routine is called after patch boundary data is filled
     * (i.e., ghosts) and before computeFluxesOnPatch().
     *
     * Note that when this routine is called, the scratch data is
     * filled on all patches (i.e., ghost cells) and that data is the
     * same as the current level data on all patch interiors.  That
     * is, both scratch and current data correspond to current_time.
     */
    virtual void preprocessAdvanceLevelState(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >& level,
        double current_time,
        double dt,
        bool first_step,
        bool last_step,
        bool regrid_advance);

    /*!
     * Add source terms to the updated solution.
     *
     * This routine is called after conservativeDifferenceOnPatch() is
     * called and before computeStableDtOnPatch().
     *
     * Note that when this routine is called, the scratch data is
     * filled on all patches (i.e., ghost cells) and that data is the
     * same as the new level data on all patch interiors.  That is,
     * both scratch and new data correspond to current_time + dt on
     * patch interiors.  The current data and ghost values correspond
     * to the current_time.
     */
    virtual void postprocessAdvanceLevelState(
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
     * NOTE: This constructor is not implemented and should not be
     * used.
     */
    AdvDiffTimeRefinementOps();

    /*!
     * \brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    AdvDiffTimeRefinementOps(
        const AdvDiffTimeRefinementOps& from);

    /*!
     * \brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AdvDiffTimeRefinementOps& operator=(
        const AdvDiffTimeRefinementOps& that);
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/AdvDiffTimeRefinementOps.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_AdvDiffTimeRefinementOps
