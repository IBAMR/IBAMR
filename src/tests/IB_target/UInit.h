#ifndef included_UInit
#define included_UInit

// Filename: UInit.h
// Last modified: <16.Feb.2007 16:33:14 griffith@box221.cims.nyu.edu>
// Created on 24 Oct 2006 by Boyce Griffith (boyce@bigboy.nyconnect.com)

/////////////////////////////// INCLUDES /////////////////////////////////////

// STOOLS INCLUDES
#include <stools/SetDataStrategy.h>

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
        const string& object_name);

    /*!
     * \brief Destructor.
     */
    virtual ~UInit();

    /*!
     * Indicates whether the concrete SetDataStrategy object is time
     * dependent.
     */
    virtual bool isTimeDependent() const { return false; }

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
    UInit();

    /*!
     * \brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    UInit(
        const UInit& from);

    /*!
     * \brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    UInit& operator=(
        const UInit& that);

    /*
     * The object name is used as a handle to databases stored in
     * restart files and for error reporting purposes.
     */
    string d_object_name;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "UInit.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_UInit
