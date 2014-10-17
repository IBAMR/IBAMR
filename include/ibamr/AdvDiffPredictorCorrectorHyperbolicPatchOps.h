// Filename: AdvDiffPredictorCorrectorHyperbolicPatchOps.h
// Created on 19 Mar 2004 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

#ifndef included_AdvDiffPredictorCorrectorHyperbolicPatchOps
#define included_AdvDiffPredictorCorrectorHyperbolicPatchOps

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>

#include "ibamr/AdvectorPredictorCorrectorHyperbolicPatchOps.h"
#include "tbox/Pointer.h"

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
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        SAMRAI::tbox::Pointer<AdvectorExplicitPredictorPatchOps> explicit_predictor,
        SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom,
        bool register_for_restart = true);

    /*!
     * The destructor for AdvDiffPredictorCorrectorHyperbolicPatchOps unregisters the patch
     * strategy object with the restart manager when so registered.
     */
    ~AdvDiffPredictorCorrectorHyperbolicPatchOps();

    /*!
     * Update solution variables by performing a conservative difference using
     * the fluxes calculated in computeFluxesOnPatch().
     */
    void
    conservativeDifferenceOnPatch(SAMRAI::hier::Patch<NDIM>& patch, double time, double dt, bool at_synchronization);

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
                                     bool regrid_advance);

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
                                      bool regrid_advance);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    AdvDiffPredictorCorrectorHyperbolicPatchOps();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AdvDiffPredictorCorrectorHyperbolicPatchOps(const AdvDiffPredictorCorrectorHyperbolicPatchOps& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AdvDiffPredictorCorrectorHyperbolicPatchOps& operator=(const AdvDiffPredictorCorrectorHyperbolicPatchOps& that);
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_AdvDiffPredictorCorrectorHyperbolicPatchOps
