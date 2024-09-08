// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2023 by the IBAMR developers
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

#ifndef included_IBAMR_SOAcousticStreamingBcCoefs
#define included_IBAMR_SOAcousticStreamingBcCoefs

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibtk/ibtk_utilities.h"

#include "CartesianGridGeometry.h"
#include "IntVector.h"
#include "RobinBcCoefStrategy.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BoundaryBox;
template <int DIM>
class Patch;
template <int DIM>
class Variable;
} // namespace hier
namespace pdat
{
template <int DIM, class TYPE>
class ArrayData;
} // namespace pdat
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!  \brief Class SOAcousticStreamingBcCoefs is an implementation of the strategy
 * class SAMRAI::solv::RobinBcCoefStrategy that allows for the run-time
 * specification of (possibly spatially- and temporally-varying) velocity boundary
 * conditions for the second-order acoustic streaming system.
 */
class SOAcousticStreamingBcCoefs : public SAMRAI::solv::RobinBcCoefStrategy<NDIM>
{
public:
    /*
     * \brief Set velocity component
     */
    void setSOVelocityComponent(int U2_comp);

    /*
     * \brief Set first-order solution patch data indices.
     */
    void setFOVelocityPressurePatchDataIndices(int U1_real_idx, int U1_imag_idx, int p1_real_idx, int p1_imag_idx);

    /*
     * \brief Set (zeroth-order) density patch data index.
     */
    void setDensityPatchDataIndex(int rho_idx);

    /*
     * \brief Set sound speed.
     */
    void setSoundSpeed(double sound_speed);

    /*
     * \brief Set acoustic angular frequency.
     */
    void setAcousticAngularFrequency(double omega);

    /*
     * \brief Set if we are using Stokes drift as the velocity boundary condition for the second order system.
     */
    void useStokesDriftVelocityForm(bool use_stokes_drift_vel);

    /*!
     * \name Implementation of SAMRAI::solv::RobinBcCoefStrategy interface.
     */
    //\{

    /*!
     * \brief Function to fill arrays of Robin boundary condition coefficients
     * at a patch boundary.
     *
     * \see SAMRAI::solv::RobinBcCoefStrategy::setBcCoefs()
     *
     * \param acoef_data  Boundary coefficient data.
     *        The array will have been defined to include index range
     *        for corresponding to the boundary box \a bdry_box and
     *        appropriate for the alignment of the given variable.  If
     *        this is a null pointer, then the calling function is not
     *        interested in a, and you can disregard it.
     * \param bcoef_data  Boundary coefficient data.
     *        This array is exactly like \a acoef_data, except that it
     *        is to be filled with the b coefficient.
     * \param gcoef_data  Boundary coefficient data.
     *        This array is exactly like \a acoef_data, except that it
     *        is to be filled with the g coefficient.
     * \param variable    Variable to set the coefficients for.
     *        If implemented for multiple variables, this parameter
     *        can be used to determine which variable's coefficients
     *        are being sought.
     * \param patch       Patch requiring bc coefficients.
     * \param bdry_box    Boundary box showing where on the boundary the coefficient data is
     *needed.
     * \param fill_time   Solution time corresponding to filling, for use when coefficients are
     *time-dependent.
     */
    void setBcCoefs(SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double> >& acoef_data,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double> >& bcoef_data,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double> >& gcoef_data,
                    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& variable,
                    const SAMRAI::hier::Patch<NDIM>& patch,
                    const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
                    double fill_time = 0.0) const override;

    /*
     * \brief Return how many cells past the edge or corner of the patch the
     * object can fill.
     *
     * The "extension" used here is the number of cells that a boundary box
     * extends past the patch in the direction parallel to the boundary.
     *
     * Note that the inability to fill the sufficient number of cells past the
     * edge or corner of the patch may preclude the child class from being used
     * in data refinement operations that require the extra data, such as linear
     * refinement.
     *
     * The boundary box that setBcCoefs() is required to fill should not extend
     * past the limits returned by this function.
     */
    SAMRAI::hier::IntVector<NDIM> numberOfExtensionsFillable() const override;

    //\}

private:
    /*!
     * \brief Second order velocity component
     */
    int d_U2_comp = -1;

    /*!
     * \brief Patch data indices of first-order velocity and pressure components.
     */
    int d_U1_real_idx = IBTK::invalid_index, d_U1_imag_idx = IBTK::invalid_index, d_p1_real_idx = IBTK::invalid_index,
        d_p1_imag_idx = IBTK::invalid_index;

    /*
     * \brief Patch data index of zeroth order density
     */
    int d_rho_idx = IBTK::invalid_index;

    /*
     * \brief Sound speed and acoustic angular frequency
     */
    double d_sound_speed = std::numeric_limits<double>::signaling_NaN(),
           d_acoustic_freq = std::numeric_limits<double>::signaling_NaN();

    /*
     * \brief If we are using Stokes drift as the velocity boundary condition.
     */
    bool d_use_stokes_drift_velocity = true;

    /*
     * \brief Set u2 from first-order mass flux.
     */
    void setBcCoefsFromFOMassFlux(SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double> >& acoef_data,
                                  SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double> >& bcoef_data,
                                  SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double> >& gcoef_data,
                                  const SAMRAI::hier::Patch<NDIM>& patch,
                                  const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
                                  double fill_time) const;

    /*
     * \brief Set u2 from Stokes drift.
     */
    void setBcCoefsFromStokesDrift(SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double> >& acoef_data,
                                   SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double> >& bcoef_data,
                                   SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double> >& gcoef_data,
                                   const SAMRAI::hier::Patch<NDIM>& patch,
                                   const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
                                   double fill_time) const;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_SOAcousticStreamingBcCoefs
