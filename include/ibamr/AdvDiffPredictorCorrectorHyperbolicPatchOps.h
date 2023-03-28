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

#ifndef included_IBAMR_AdvDiffPredictorCorrectorHyperbolicPatchOps
#define included_IBAMR_AdvDiffPredictorCorrectorHyperbolicPatchOps

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/AdvectorPredictorCorrectorHyperbolicPatchOps.h"

#include "tbox/Pointer.h"

#include <string>

namespace IBAMR
{
class AdvectorExplicitPredictorPatchOps;
} // namespace IBAMR
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
class Patch;
template <int DIM>
class PatchLevel;
} // namespace hier
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class AdvDiffPredictorCorrectorHyperbolicPatchOps is a specialization of class
 * AdvectorPredictorCorrectorHyperbolicPatchOps for use with a linearly-implicit time
 *integrator for
 * the advection-diffusion equation.
 *
 * \see AdvDiffGodunovHierarchyIntegrator
 * \see AdvectorPredictorCorrectorHyperbolicPatchOps
 */
class AdvDiffPredictorCorrectorHyperbolicPatchOps : public AdvectorPredictorCorrectorHyperbolicPatchOps
{
public:
    /*!
     * The constructor for AdvDiffPredictorCorrectorHyperbolicPatchOps sets default parameters
     *for
     * the patch strategy.  The constructor also registers this object for
     * restart with the restart manager using the object name when so requested.
     *
     * After default values are set, this routine calls getFromRestart() if
     * execution from a restart file is specified.  Finally, getFromInput() is
     * called to read values from the given input database (potentially
     * overriding those found in the restart file).
     */
    AdvDiffPredictorCorrectorHyperbolicPatchOps(
        std::string object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        SAMRAI::tbox::Pointer<AdvectorExplicitPredictorPatchOps> explicit_predictor,
        SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom,
        bool register_for_restart = true);

    /*!
     * The destructor for AdvDiffPredictorCorrectorHyperbolicPatchOps unregisters the patch
     * strategy object with the restart manager when so registered.
     */
    ~AdvDiffPredictorCorrectorHyperbolicPatchOps() = default;

    /*!
     * Update solution variables by performing a conservative difference using
     * the fluxes calculated in computeFluxesOnPatch().
     */
    void conservativeDifferenceOnPatch(SAMRAI::hier::Patch<NDIM>& patch,
                                       double time,
                                       double dt,
                                       bool at_synchronization) override;

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
    void preprocessAdvanceLevelState(const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >& level,
                                     double current_time,
                                     double dt,
                                     bool first_step,
                                     bool last_step,
                                     bool regrid_advance) override;

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
    void postprocessAdvanceLevelState(const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >& level,
                                      double current_time,
                                      double dt,
                                      bool first_step,
                                      bool last_step,
                                      bool regrid_advance) override;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    AdvDiffPredictorCorrectorHyperbolicPatchOps() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AdvDiffPredictorCorrectorHyperbolicPatchOps(const AdvDiffPredictorCorrectorHyperbolicPatchOps& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AdvDiffPredictorCorrectorHyperbolicPatchOps&
    operator=(const AdvDiffPredictorCorrectorHyperbolicPatchOps& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_AdvDiffPredictorCorrectorHyperbolicPatchOps
