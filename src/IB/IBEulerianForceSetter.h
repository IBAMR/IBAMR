#ifndef included_IBEulerianForceSetter
#define included_IBEulerianForceSetter

// Filename: IBEulerianForceSetter.h
// Created on 28 Sep 2004 by Boyce Griffith (boyce@mstu1.cims.nyu.edu)
// Last modified: <11.Jan.2007 16:33:24 griffith@box221.cims.nyu.edu>

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/SetDataStrategy.h>

// SAMRAI INCLUDES
#include <Patch.h>
#include <Variable.h>
#include <tbox/Pointer.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * @brief Simple class to specify the body force on the Eulerian grid
 * during an IB computation.
 */
class IBEulerianForceSetter
    : public SetDataStrategy
{
public:
    /*!
     * @brief Constructor.
     */
    IBEulerianForceSetter(
        const std::string& object_name,
        const int F_idx);

    /*!
     * @brief Destructor.
     */
    ~IBEulerianForceSetter();

    //@{ @name Methods to set the data.

    /*!
     * This concrete SetDataStrategy class is time dependent.
     */
    bool isTimeDependent() const;

    /*!
     * Set the data on the patch interior.
     */
    void setDataOnPatch(
        const int data_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
        SAMRAI::hier::Patch<NDIM>& patch,
        const double data_time,
        const bool initial_time=false);

    //@}

private:
    /*!
     * @brief Default constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     */
    IBEulerianForceSetter();

    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * @param from The value to copy to this object.
     */
    IBEulerianForceSetter(
        const IBEulerianForceSetter& from);

    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * @param that The value to assign to this object.
     *
     * @return A reference to this object.
     */
    IBEulerianForceSetter& operator=(
        const IBEulerianForceSetter& that);

    /*
     * Patch data descriptor index for the force data.
     */
    const int d_F_idx;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBEulerianForceSetter.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBEulerianForceSetter
