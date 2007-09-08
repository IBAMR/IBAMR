#ifndef included_INSIntermediateVelocityBcCoef
#define included_INSIntermediateVelocityBcCoef

// Filename: INSIntermediateVelocityBcCoef.h
// Last modified: <07.Sep.2007 22:18:37 griffith@box221.cims.nyu.edu>
// Created on 30 Aug 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// STOOLS INCLUDES
#include <stools/ExtendedRobinBcCoefStrategy.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSIntermediateVelocityBcCoef is a concrete
 * SAMRAI::solv::RobinBcCoefStrategy that is used to specify boundary conditions
 * for the intermediate velocity field that is computed in a projection method.
 */
class INSIntermediateVelocityBcCoef
    : public virtual STOOLS::ExtendedRobinBcCoefStrategy
{
public:
    /*!
     * \brief Constructor.
     *
     * \param comp_idx        Component of the velocity which this boundary condition specification is to operate on
     * \param Phi_idx         Patch data descriptor index for the cell-centered scalar function Phi
     * \param u_bc_coefs      Vector of boundary condition specification objects corresponding to the components of the velocity.
     * \param homogeneous_bc  Whether to employ homogeneous (as opposed to inhomogeneous) boundary conditions
     *
     * \note Precisely NDIM boundary condition objects must be provided to the
     * class constructor.
     */
    INSIntermediateVelocityBcCoef(
        const int comp_idx,
        const int Phi_idx,
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
        const bool homogneous_bc=false);

    /*!
     * \brief Destructor.
     */
    virtual
    ~INSIntermediateVelocityBcCoef();

    /*!
     * \brief Reset the boundary condition specification object so that it
     * specifies the "true" velocity boundary conditions.
     */
    void
    useTrueVelocityBcCoefs();

    /*!
     * \brief Reset the boundary condition specification object so that it
     * specifies the "intermediate" velocity boundary conditions.
     */
    void
    useIntermediateVelocityBcCoefs(
        const double current_time,
        const double new_time,
        const double rho,
        const bool velocity_correction);

    /*!
     * \brief Reset the patch data descriptor index for the cell-centered scalar
     * function Phi.
     *
     * \param Phi_idx  Patch data descriptor index for the cell-centered scalar function Phi
     */
    void
    setPhiPatchDataIndex(
        const int Phi_idx);

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
     * \name Implementation of STOOLS::ExtendedRobinBcCoefStrategy interface.
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

    /*!
     * \return Set of patch data descriptor indices whose values must be
     * provided in order to set the homogeneous Robin boundary conditions.
     */
    virtual const std::set<int>&
    getHomogeneousBcFillDataIndices() const;

    /*!
     * \return Set of patch data descriptor indices whose values must be
     * provided in order to set the inhomogeneous Robin boundary conditions.
     */
    virtual const std::set<int>&
    getInhomogeneousBcFillDataIndices() const;

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

protected:

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSIntermediateVelocityBcCoef();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSIntermediateVelocityBcCoef(
        const INSIntermediateVelocityBcCoef& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSIntermediateVelocityBcCoef&
    operator=(
        const INSIntermediateVelocityBcCoef& that);

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
     * Component of the velocity which this boundary condition specification is to operate on.
     */
    const int d_comp_idx;

    /*
     * The patch data index corresponding to the current value of U_star.
     */
    int d_target_idx;

    /*
     * The patch data index corresponding to the estimated value of Phi.
     */
    int d_Phi_idx;

    /*
     * Sets of patch data indices required to fill homogeneous and inhomogeneous
     * boundary conditions.
     */
    std::set<int> d_homo_patch_data_idxs, d_inhomo_patch_data_idxs;

    /*
     * The boundary condition specification objects for the updated velocity.
     */
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_u_bc_coefs;

    /*
     * Whether to use homogeneous boundary conditions.
     */
    bool d_homogeneous_bc;

    /*
     * Current time, new time, and timestep size.
     */
    double d_current_time, d_new_time;

    /*
     * Fluid density rho.
     */
    double d_rho;

    /*
     * Whether we are using the "true" or "intermediate" velocity bc coefs.
     */
    bool d_using_intermediate_velocity_bc_coefs;

    /*
     * Whether we are doing a velocity correction or not.
     */
    bool d_velocity_correction;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/INSIntermediateVelocityBcCoef.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSIntermediateVelocityBcCoef
