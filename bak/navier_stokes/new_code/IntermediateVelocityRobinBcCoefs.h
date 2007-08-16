#ifndef included_IntermediateVelocityRobinBcCoefs
#define included_IntermediateVelocityRobinBcCoefs

// Filename: IntermediateVelocityRobinBcCoefs.h
// Last modified: <16.May.2007 21:13:09 griffith@box221.cims.nyu.edu>
// Created on 09 Mar 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// STOOLS INCLUDES
#include <stools/ExtendedRobinBcCoefStrategy.h>

// SAMRAI INCLUDES
#include <CellVariable.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IntermediateVelocityRobinBcCoefs is an implementation of the
 * strategy class STOOLS::ExtendedRobinBcCoefStrategy that is used to specify
 * physical boundary conditions for a single component of the intermediate
 * velocity in a projection method.
 */
class IntermediateVelocityRobinBcCoefs
    : public STOOLS::ExtendedRobinBcCoefStrategy
{
public:
    /*!
     * \brief Constructor
     */
    IntermediateVelocityRobinBcCoefs(
        const int velocity_depth,
        const bool using_pressure_increment_form,
        const SAMRAI::solv::RobinBcCoefStrategy<NDIM>* const bc_coef);

    /*!
     * \brief Destructor.
     */
    virtual
    ~IntermediateVelocityRobinBcCoefs();

    /*!
     * \brief Specify whether or not to modify the boundary conditions.
     */
    void
    correctBcCoefs(
        const bool correct_bc_coefs);

    /*!
     * \brief Specify the current time interval.
     */
    void
    setTimeInterval(
        const double current_time,
        const double new_time);

    /*!
     * \brief Specify the (uniform) mass density.
     */
    void
    setRho(
        const double rho);

    /*!
     * \brief Specify the pressure variable to be used to update the boundary
     * conditions.
     */
    void
    setPressureVariable(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_var,
        const int P_src_idx,
        const int P_dst_idx);

    /*!
     * \brief Specify the phi variable to be used to update the boundary
     * conditions.
     */
    void
    setPhiVariable(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Phi_var,
        const int Phi_src_idx,
        const int Phi_dst_idx);

    /*!
     * \name Extended SAMRAI::solv::RobinBcCoefStrategy interface.
     */
    //\{

    /*!
     * \brief Function to specify the variables and patch data indices whose
     * values must be filled in order to specify boundary condition data.
     */
    virtual void
    getBcFillVars(
        std::vector<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > >& bc_fill_vars,
        std::vector<int>& bc_fill_src_idxs,
        std::vector<int>& bc_fill_dst_idxs) const;

    //\}

    /*!
     * \name Implementation of SAMRAI::solv::RobinBcCoefStrategy interface.
     */
    //\{

    /*!
     * \brief Function to fill arrays of Robin boundary condition coefficients
     * at a patch boundary.  (New interface.)
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
     *
     * \note An unrecoverable exception will occur if this method is called when
     * IBAMR is compiled with SAMRAI version 2.1.
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

    /*!
     * \brief Function to fill arrays of Robin boundary condition coefficients
     * at a patch boundary.  (Old interface.)
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
     * \param gcoef_data  Boundary coefficient data.
     *        This array is exactly like \a acoef_data, except that it is to be
     *        filled with the g coefficient.
     * \param variable    Variable to set the coefficients for.
     *        If implemented for multiple variables, this parameter can be used
     *        to determine which variable's coefficients are being sought.
     * \param patch       Patch requiring bc coefficients.
     * \param bdry_box    Boundary box showing where on the boundary the coefficient data is needed.
     * \param fill_time   Solution time corresponding to filling, for use when coefficients are time-dependent.
     *
     * \note An unrecoverable exception will occur if this method is called when
     * IBAMR is compiled with SAMRAI versions after version 2.1.
     */
    virtual void
    setBcCoefs(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& acoef_data,
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

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IntermediateVelocityRobinBcCoefs();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IntermediateVelocityRobinBcCoefs(
        const IntermediateVelocityRobinBcCoefs& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IntermediateVelocityRobinBcCoefs&
    operator=(
        const IntermediateVelocityRobinBcCoefs& that);

    /*!
     * \brief Correct the boundary coefficients using p(n-1/2) and phi(n-1/2).
     */
    void
    computeCorrectedBcCoefs(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& acoef_data,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& bcoef_data,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& gcoef_data,
        const SAMRAI::hier::Patch<NDIM>& patch,
        const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box) const;

    /*
     * The component of the velocity to correct.
     */
    const int d_velocity_depth;

    /*
     * Whether the projection algorithm is using pressure increment form.
     */
    const bool d_using_pressure_increment_form;

    /*
     * We use this object to actually set the boundary values.
     */
    const SAMRAI::solv::RobinBcCoefStrategy<NDIM>* const d_bc_coef;

    /*
     * Whether the class should correct the boundary condition coefficients.
     */
    bool d_correct_bc_coefs;

    /*
     * The current and new simulation times.
     */
    double d_current_time, d_new_time;

    /*
     * The uniform mass density.
     */
    double d_rho;

    /*
     * Patch data descriptor indices indicating the data to be employed to set
     * the boundary conditions.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_P_var;
    int d_P_src_idx, d_P_dst_idx;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_Phi_var;
    int d_Phi_src_idx, d_Phi_dst_idx;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IntermediateVelocityRobinBcCoefs.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IntermediateVelocityRobinBcCoefs
