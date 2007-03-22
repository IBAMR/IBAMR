#ifndef included_IBEulerianSourceSetter
#define included_IBEulerianSourceSetter

// Filename: IBEulerianForceSetter.h
// Created on 18 Jun 2005 by Boyce Griffith (boyce@bigboy.verizon.net)
// Last modified: <21.Mar.2007 20:19:52 griffith@box221.cims.nyu.edu>

/////////////////////////////// INCLUDES /////////////////////////////////////

// STOOLS INCLUDES
#include <stools/SetDataStrategy.h>

// SAMRAI INCLUDES
#include <Patch.h>
#include <Variable.h>
#include <tbox/Pointer.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Simple class to specify sources and sinks on the Eulerian
 * grid during an IB computation.
 */
class IBEulerianSourceSetter
    : public STOOLS::SetDataStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    IBEulerianSourceSetter(
        const std::string& object_name,
        const int Q_idx);

    /*!
     * \brief Virtual destructor.
     */
    virtual ~IBEulerianSourceSetter();

    /*!
     * \name Methods to set the data.
     */
    //\{

    /*!
     * \note This concrete STOOLS::SetDataStrategy is time-dependent.
     */
    virtual bool isTimeDependent() const;

    /*!
     * Set the data on the patch interior.
     */
    virtual void setDataOnPatch(
        const int data_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
        SAMRAI::hier::Patch<NDIM>& patch,
        const double data_time,
        const bool initial_time=false);

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     */
    IBEulerianSourceSetter();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    IBEulerianSourceSetter(
        const IBEulerianSourceSetter& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
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

//#include <ibamr/IBEulerianSourceSetter.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBEulerianSourceSetter
