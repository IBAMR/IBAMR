#ifndef included_PInit
#define included_PInit

// Filename: PInit.h
// Last modified: <16.Feb.2007 16:34:23 griffith@box221.cims.nyu.edu>
// Created on 19 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// STOOLS INCLUDES
#include <stools/SetDataStrategy.h>

// NAMESPACE
using namespace SAMRAI;
using namespace STOOLS;
using namespace std;

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Method to initialize the value of the hydrostatic pressure P
 * (used for convergence analysis only).
 */
class PInit
    : public SetDataStrategy
{
public:
    /*!
     * \brief Constructor.
     *
     * \param object_name the object name
     * \param nu the kinematic viscosity (\f$\frac{\mu}{\rho}\f$)
     */
    PInit(
        const string& object_name,
        const double nu);

    /*!
     * \brief Destructor.
     */
    virtual ~PInit();

    /*!
     * Indicates whether the concrete SetDataStrategy object is time
     * dependent.
     */
    virtual bool isTimeDependent() const { return true; }

    /*!
     * Set the data on the patch interior to the exact answer.
     */
    virtual void setDataOnPatch(
        const int data_idx,
        tbox::Pointer<hier::Variable<NDIM> > var,
        hier::Patch<NDIM>& patch,
        const double data_time,
        const bool initial_time=false);

protected:

private:
    /*!
     * \brief Default constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     */
    PInit();

    /*!
     * \brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    PInit(
        const PInit& from);

    /*!
     * \brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PInit& operator=(
        const PInit& that);

    /*
     * The kinematic viscosity.
     */
    double d_nu;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "PInit.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PInit
