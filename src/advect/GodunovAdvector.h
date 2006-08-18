#ifndef included_GodunovAdvector
#define included_GodunovAdvector

// Filename: GodunovAdvector.h
// Last modified: <17.Aug.2006 20:11:32 boyce@bigboy.nyconnect.com>
// Created on 14 Feb 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <CellData.h>
#include <FaceData.h>
#include <Patch.h>
#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Serializable.h>

// C++ STDLIB INCLUDES
#include <ostream>
#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * The GodunovAdvector class implements the predictors required to use
 * an explicit predictor-corrector method to solve the
 * non-conservative advection equation
 *
 *      dQ/dt + u * grad(Q) = Psi
 * 
 * where Q is a cell centered quantity, u is a specified face centered
 * advection velocity, and Psi is an optional source term.  These
 * routines can also be used to solve the conservative advection
 * equation
 * 
 *      dQ/dt + div(u Q) = Psi
 *
 * The class implements an upwind (Godunov) explicit predictor which
 * can be used to generate time and face cenetered values or fluxes.
 * These predicted values can be used in a second order accurate
 * predictor-corrector method for solving the advection equation (and
 * related problems, including the advection-diffusion equation and
 * the equations of incompressible flow).
 * 
 * Note that the predicted fluxes are computed using the advective
 * (i.e. non-conservative) form of the advection equation.
 * Consequently, when the advection velocity u is not discretely
 * divergence free, the appropriate non-conservative form of the
 * source term must be supplied to the predictor in order to obtain a
 * consistent method.
 */
class GodunovAdvector
    : public SAMRAI::tbox::Serializable
{
public:
    /*!
     * Enumerated type for the different supported fluxes/limiters
     * available.
     */
    enum FLUX_TYPE
    {
        CTU_ONLY = -1,
        SECOND_ORDER = 0,
        FOURTH_ORDER = 1,
        MC_LIMITED = 2,
        MUSCL_LIMITED = 3
    };
    
    /*!
     * The constructor for GodunovAdvector sets default parameters for
     * the advection predictor.  The constructor also registers this
     * object for restart with the restart manager using the object
     * name.
     *
     * After default values are set, this routine calls
     * getFromRestart() if execution from a restart file is specified.
     * Finally, getFromInput() is called to read values from the given
     * input database (potentially overriding those found in the
     * restart file).
     */
    GodunovAdvector(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        bool register_for_restart=true);
    
    /*!
     * The destructor for GodunovAdvector unregisters the predictor
     * object with the restart manager when so registered.
     */
    ~GodunovAdvector();
    
    ///
    ///  The following routines:
    ///
    ///      computeStableDtOnPatch(),
    ///      computeAdvectiveDerivative()
    ///      computeFlux(),
    ///      predictValue(),
    ///      predictValueWithSourceTerm(),
    ///      predictNormalVelocity(),
    ///      predictNormalVelocityWithSourceTerm()
    ///      
    ///  provide the explicit predictors required to solve the
    ///  advection equation using a predictor-corrector methodology.
    ///
    
    /*!
     * Compute a stable time increment for the patch using an explicit
     * CFL condition and return the computed dt.
     */
    double computeStableDtOnPatch(
        const SAMRAI::pdat::FaceData<NDIM,double>& u_ADV,
        const SAMRAI::hier::Patch<NDIM>& patch) const;
 
    /*!
     * Compute the advective derivative N = [u_ADV * grad(q)], using
     * the specified advection velocity and predicted face centered
     * values.
     */
    void computeAdvectiveDerivative(
        SAMRAI::pdat::CellData<NDIM,double>& N,
        const SAMRAI::pdat::FaceData<NDIM,double>& u_ADV,
        const SAMRAI::pdat::FaceData<NDIM,double>& q,
        const SAMRAI::hier::Patch<NDIM>& patch) const;
    
    /*!
     * Computes the time integral of the advective fluxes
     * corresponding to a face centered value and a face centered
     * advective velocity.
     *
     * The advective fluxes are defined to be
     *
     *      f(i+1/2,j) = dt*u(i+1/2,j)*q(i+1/2,j)
     *      f(i,j+1/2) = dt*v(i,j+1/2)*q(i,j+1/2)
     *
     * in two spatial dimensions, and similarly for other dimensions.
     */
    void computeFlux(
        SAMRAI::pdat::FaceData<NDIM,double>& flux,
        const SAMRAI::pdat::FaceData<NDIM,double>& u_ADV,
        const SAMRAI::pdat::FaceData<NDIM,double>& q_half,
        const SAMRAI::hier::Patch<NDIM>& patch,
        const double dt) const;
    
    /*!
     * A Godunov predictor used to predict face and time centered
     * values from cell centered values using a Taylor expansion about
     * each cell center.
     * 
     * The predictor assumes that Q satisfies an equation of the form
     * 
     *      dQ/dt + u * grad(Q) = 0
     * 
     * i.e. Q satisfies a non-conservative advection equation.
     *
     * @see predictValueWithSourceTerm
     */
    void predictValue(
        SAMRAI::pdat::FaceData<NDIM,double>& q_half,
        const SAMRAI::pdat::FaceData<NDIM,double>& u_ADV,
        const SAMRAI::pdat::CellData<NDIM,double>& Q,
        const SAMRAI::hier::Patch<NDIM>& patch,
        const double dt) const;
    
    /*!
     * A Godunov predictor used to predict face and time centered
     * values from cell centered values using a Taylor expansion about
     * each cell center.
     * 
     * The predictor assumes that Q satisfies an equation of the form
     * 
     *      dQ/dt + u * grad(Q) = Psi
     * 
     * i.e. Q satisfies a non-conservative advection equation.
     */
    void predictValueWithSourceTerm(
        SAMRAI::pdat::FaceData<NDIM,double>& q_half,
        const SAMRAI::pdat::FaceData<NDIM,double>& u_ADV,
        const SAMRAI::pdat::CellData<NDIM,double>& Q,
        const SAMRAI::pdat::CellData<NDIM,double>& Psi,
        const SAMRAI::hier::Patch<NDIM>& patch,
        const double dt) const;
    
    /*!
     * A Godunov predictor used to predict face and time centered
     * normal velocities from cell centered values using a Taylor
     * expansion about each cell center.
     * 
     * The predictor assumes that V satisfies an equation of the form
     * 
     *      dV/dt + u * grad V = 0
     * 
     * i.e. V satisfies a non-conservative advection equation.
     *
     * Note that if the predicted velocities are to be conservatively
     * differenced (i.e. used in a discretization of the conservative
     * form of the equation), a consistent method is obtained only
     * when the proper source term is included should the advection
     * velocity not be discretely divergence free.
     *
     * @see predictNormalVelocityWithSourceTerm
     */
    void predictNormalVelocity(
        SAMRAI::pdat::FaceData<NDIM,double>& v_half,
        const SAMRAI::pdat::FaceData<NDIM,double>& u_ADV,
        const SAMRAI::pdat::CellData<NDIM,double>& V,
        const SAMRAI::hier::Patch<NDIM>& patch,
        const double dt) const;
    
    /*!
     * A Godunov predictor used to predict face and time centered
     * normal velocities from cell centered values using a Taylor
     * expansion about each cell center.
     * 
     * The predictor assumes that V satisfies an equation of the form
     * 
     *      dV/dt + u * grad V = Psi
     * 
     * i.e. V satisfies a non-conservative advection equation.
     *
     * Note that if the predicted velocities are to be conservatively
     * differenced (i.e. used in a discretization of the conservative
     * form of the equation), a consistent method is obtained only
     * when the proper source term is included should the advection
     * velocity not be discretely divergence free.
     */
    void predictNormalVelocityWithSourceTerm(
        SAMRAI::pdat::FaceData<NDIM,double>& V_half,
        const SAMRAI::pdat::FaceData<NDIM,double>& u_ADV,
        const SAMRAI::pdat::CellData<NDIM,double>& V,
        const SAMRAI::pdat::CellData<NDIM,double>& Psi,
        const SAMRAI::hier::Patch<NDIM>& patch,
        const double dt) const;
    
    /*!
     * Subtracts off the face centered gradient of a scalar from a
     * predicted velocity to approximately enforce incompressibility.
     */
    void enforceIncompressibility(
        SAMRAI::pdat::FaceData<NDIM,double>& q_half,
        const SAMRAI::pdat::FaceData<NDIM,double>& u_ADV,
        const SAMRAI::pdat::FaceData<NDIM,double>& grad,
        const SAMRAI::hier::Patch<NDIM>& patch,
        const double dt) const;

    ///
    ///  The following routines:
    ///
    ///      putToDatabase()
    ///
    ///  are concrete implementations of functions declared in the
    ///  SAMRAI::tbox::Serializable abstract base class.
    ///
    
    /*!
     * Write state of GodunovAdvector object to the given database for
     * restart.
     *
     * This routine is a concrete implementation of the function
     * declared in the SAMRAI::tbox::Serializable abstract base class.
     */
    virtual void putToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    ///
    ///  The following routines:
    ///
    ///      printClassData()
    ///
    ///  are provided for your viewing pleasure.
    ///
    
    /*!
     * Print all data members for GodunovAdvector class.
     */
    virtual void printClassData(
        std::ostream& os) const;
    
private:
    /*!
     * @brief Default constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     */
    GodunovAdvector();
    
    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     * 
     * @param from The value to copy to this object.
     */
    GodunovAdvector(
        const GodunovAdvector& from);
    
    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     * 
     * @param that The value to assign to this object.
     * 
     * @return A reference to this object.
     */
    GodunovAdvector& operator=(
        const GodunovAdvector& that);

    /*
     * Private functions used to compute the predicted values/fluxes.
     */
    void predict(
        SAMRAI::pdat::FaceData<NDIM,double>& q_half,
        const SAMRAI::pdat::FaceData<NDIM,double>& u_ADV,
        const SAMRAI::pdat::CellData<NDIM,double>& Q,
        const SAMRAI::hier::Patch<NDIM>& patch,
        const double dt) const;

    void predictWithSourceTerm(
        SAMRAI::pdat::FaceData<NDIM,double>& q_half,
        const SAMRAI::pdat::FaceData<NDIM,double>& u_ADV,
        const SAMRAI::pdat::CellData<NDIM,double>& Q,
        const SAMRAI::pdat::CellData<NDIM,double>& Psi,
        const SAMRAI::hier::Patch<NDIM>& patch,
        const double dt) const;
    
    /*
     * These private member functions read data from input and
     * restart.  When beginning a run from a restart file, all data
     * members are read from the restart file.  If the boolean flag is
     * true when reading from input, some restart values may be
     * overridden by those in the input file.
     *
     * An assertion results if the database pointer is null.
     */
    void getFromInput(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
        bool is_from_restart);
    
    void getFromRestart();
    
    /*
     * The object name is used as a handle to databases stored in
     * restart files and for error reporting purposes.  The boolean is
     * used to control restart file writing operations.
     */
    std::string d_object_name;
    bool d_registered_for_restart;

    /*
     *  Parameters for numerical method:
     *
     *    d_limiter_type ........ specifies the type of slope limiting
     *                            used in computing numerical fluxes
     *    d_use_full_ctu ........ specifies whether full corner transport
     *                            upwinding is used for 3D computations
     */
    int d_limiter_type;
#if (NDIM == 3)
    bool d_use_full_ctu;
#endif
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#ifndef DEBUG_NO_INLINE
//#include "GodunovAdvector.I"
//#endif

#endif //#ifndef included_GodunovAdvector
