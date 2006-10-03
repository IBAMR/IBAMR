#ifndef included_IBEulerianSourceSetter
#define included_IBEulerianSourceSetter

// Filename: IBEulerianForceSetter.h
// Created on 18 Jun 2005 by Boyce Griffith (boyce@bigboy.verizon.net)
// Last modified: <03.Oct.2006 10:24:32 boyce@boyce-griffiths-powerbook-g4-15.local>

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
 * @brief Simple class to specify sources and sinks on the Eulerian
 * grid during an IB computation.
 */
class IBEulerianSourceSetter
    : public SetDataStrategy
{
public:
    /*!
     * @brief Constructor.
     */
    IBEulerianSourceSetter(
        const std::string& object_name,
        const int Q_idx);

    /*!
     * @brief Destructor.
     */
    ~IBEulerianSourceSetter();

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
    IBEulerianSourceSetter();

    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * @param from The value to copy to this object.
     */
    IBEulerianSourceSetter(
        const IBEulerianSourceSetter& from);

    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * @param that The value to assign to this object.
     *
     * @return A reference to this object.
     */
    IBEulerianSourceSetter& operator=(
        const IBEulerianSourceSetter& that);

    /*
     * Patch data descriptor index for the source/sink data.
     */
    const int d_Q_idx;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "IBEulerianSourceSetter.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBEulerianSourceSetter
