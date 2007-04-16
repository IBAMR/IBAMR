#ifndef included_UInit
#define included_UInit

// Filename: UInit.h
// Last modified: <16.Apr.2007 01:34:19 boyce@trasnaform2.local>
// Created on 19 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// STOOLS INCLUDES
#include <stools/SetDataStrategy.h>

// SAMRAI INCLUDES
#include <tbox/Database.h>

// NAMESPACE
using namespace SAMRAI;
using namespace STOOLS;
using namespace std;

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Method to initialize the value of the fluid velocity U.
 */
class UInit
    : public SetDataStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    UInit(
        const string& object_name,
        tbox::Pointer<tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    virtual
    ~UInit();

    /*!
     * Indicates whether the concrete SetDataStrategy object is time
     * dependent.
     */
    virtual bool
    isTimeDependent() const
        { return true; }

    /*!
     * Set the data on the patch interior to the exact answer.
     */
    virtual void
    setDataOnPatch(
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
     * \note This constructor is not implemented and should not be used.
     */
    UInit();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    UInit(
        const UInit& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    UInit&
    operator=(
        const UInit& that);

    /*
     * The object name is used as a handle to databases stored in restart files
     * and for error reporting purposes.
     */
    string d_object_name;

    /*
     * Shear layer paramters.
     */
    double d_rho, d_delta;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "UInit.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_UInit
