// Filename: INSProjectionBcCoef.h
// Created on 22 Feb 2007 by Boyce Griffith
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

#ifndef included_INSProjectionBcCoef
#define included_INSProjectionBcCoef

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/ibamr_enums.h>

// IBTK INCLUDES
#include <ibtk/ExtendedRobinBcCoefStrategy.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSProjectionBcCoef is a concrete
 * SAMRAI::solv::RobinBcCoefStrategy that is used to specify boundary conditions
 * for the solution of the discrete Poisson problem that must be solved in the
 * implementation of a projection method.
 */
class INSProjectionBcCoef
    : public virtual IBTK::ExtendedRobinBcCoefStrategy
{
public:
    /*!
     * \brief Constructor.
     *
     * \param P_idx            Patch data descriptor index for the cell-centered pressure field.
     * \param P_bc_coef        Boundary condition specification object corresponding to the pressure.
     * \param projection_type  The type of projection ("pressure_increment" or "pressure_update").
     * \param u_idx            Patch data descriptor index for the face-centered intermediate velocity field.
     * \param u_bc_coefs       Vector of boundary condition specification objects corresponding to the components of the velocity.
     * \param homogeneous_bc   Whether to employ homogeneous (as opposed to inhomogeneous) boundary conditions
     *
     * \note Precisely NDIM boundary condition objects must be provided to the
     * class constructor.
     */
    INSProjectionBcCoef(
        const int P_idx,
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* const P_bc_coef,
        const ProjectionMethodType& projection_type,
        const int u_idx,
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
        const bool homogeneous_bc=false);

    /*!
     * \brief Destructor.
     */
    virtual
    ~INSProjectionBcCoef();

    /*!
     * \brief Reset the problem coefficients required to specify the boundary
     * conditions for the scalar function phi.
     */
    void
    setProblemCoefs(
        const double rho,
        const double dt);

    /*!
     * \brief Reset the patch data descriptor index for the cell-centered
     * pressure.
     *
     * \param P_idx  Patch data descriptor index for the cell-centered pressure field.
     */
    void
    setCurrentPressurePatchDataIndex(
        const int P_idx);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions for the pressure.
     *
     * \param P_bc_coef  Boundary condition specification object corresponding to the pressure.
     */
    void
    setPressurePhysicalBcCoef(
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* const P_bc_coef);

    /*!
     * \brief Set the type of projection to perform.
     *
     * \param projection_type  The type of projection (PRESSURE_INCREMENT or PRESSURE_UPDATE).
     */
    void
    setProjectionType(
        const ProjectionMethodType& projection_type);

    /*!
     * \brief Reset the patch data descriptor index for the face-centered
     * intermediate velocity.
     *
     * \param u_idx  Patch data descriptor index for the face-centered intermediate velocity field.
     */
    void
    setIntermediateVelocityPatchDataIndex(
        const int u_idx);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions for the velocity.
     *
     * \param u_bc_coefs  Vector of boundary condition specification objects corresponding to the components of the velocity.
     */
    void
    setVelocityPhysicalBcCoefs(
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs);

    /*!
     * \name Implementation of IBTK::ExtendedRobinBcCoefStrategy interface.
     */
    //\{

    /*!
     * \brief Set the target data index.
     */
    virtual void
    setTargetPatchDataIndex(
        const int target_idx);

    /*!
     * \brief Set whether the class is filling homogeneous or inhomogeneous
     * boundary conditions.
     */
    virtual void
    setHomogeneousBc(
        const bool homogeneous_bc);

    //\}

    /*!
     * \name Implementation of SAMRAI::solv::RobinBcCoefStrategy interface.
     */
    //\{

    /*!
     * \brief Function to fill arrays of Robin boundary condition coefficients
     * at a patch boundary.
     *
     * \note In the original SAMRAI::solv::RobinBcCoefStrategy interface, it was
     * assumed that \f$ b = (1-a) \f$.  In the new interface, \f$a\f$ and
     * \f$b\f$ are independent.
     *
     * \see SAMRAI::solv::RobinBcCoefStrategy::setBcCoefs()
     *
     * \param acoef_data  Boundary coefficient data.
     *        The array will have been defined to include index range for
     *        corresponding to the boundary box \a bdry_box and appropriate for
     *        the alignment of the given variable.  If this is a null pointer,
     *        then the calling function is not interested in a, and you can
     *        disregard it.
     * \param bcoef_data  Boundary coefficient data.
     *        This array is exactly like \a acoef_data, except that it is to be
     *        filled with the b coefficient.
     * \param gcoef_data  Boundary coefficient data.
     *        This array is exactly like \a acoef_data, except that it is to be
     *        filled with the g coefficient.
     * \param variable    Variable to set the coefficients for.
     *        If implemented for multiple variables, this parameter can be used
     *        to determine which variable's coefficients are being sought.
     * \param patch       Patch requiring bc coefficients.
     * \param bdry_box    Boundary box showing where on the boundary the coefficient data is needed.
     * \param fill_time   Solution time corresponding to filling, for use when coefficients are time-dependent.
     */
    virtual void
    setBcCoefs(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& acoef_data,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& bcoef_data,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& gcoef_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& variable,
        const SAMRAI::hier::Patch<NDIM>& patch,
        const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
        double fill_time=0.0) const;

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
    virtual SAMRAI::hier::IntVector<NDIM>
    numberOfExtensionsFillable() const;

    //\}

protected:

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSProjectionBcCoef();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSProjectionBcCoef(
        const INSProjectionBcCoef& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSProjectionBcCoef&
    operator=(
        const INSProjectionBcCoef& that);

    /*!
     * \brief Implementation of boundary condition filling function.
     */
    void
    setBcCoefs_private(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& acoef_data,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& bcoef_data,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& gcoef_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& variable,
        const SAMRAI::hier::Patch<NDIM>& patch,
        const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
        double fill_time) const;

    /*
     * The patch data index corresponding to the current pressure.
     */
    int d_P_idx;

    /*
     * The physical boundary conditions for the pressure (only used at "open"
     * boundaries where the pressure is specified).
     */
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_P_bc_coef;

    /*
     * The type of projection.
     */
    ProjectionMethodType d_projection_type;

    /*
     * The patch data index corresponding to the intermediate velocity.
     */
    int d_u_idx;

    /*
     * The boundary condition specification objects for the updated velocity.
     */
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_u_bc_coefs;

    /*
     * Whether to use homogeneous boundary conditions.
     */
    bool d_homogeneous_bc;

    /*
     * Fluid density.
     */
    double d_rho;

    /*
     * Timestep size.
     */
    double d_dt;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/INSProjectionBcCoef.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSProjectionBcCoef
