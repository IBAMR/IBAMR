#ifndef included_SpringForceSpecFactory
#define included_SpringForceSpecFactory

// Filename: SpringForceSpecFactory.h
// Created on 14 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)
// Last modified: <03.Oct.2006 09:51:36 boyce@boyce-griffiths-powerbook-g4-15.local>

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/Stashable.h>
#include <ibamr/StashableFactory.h>

// SAMRAI INCLUDES
#include <IntVector.h>
#include <tbox/AbstractStream.h>
#include <tbox/Pointer.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * @brief Description of class.
 */
class SpringForceSpecFactory
    : public StashableFactory
{
public:
    /*!
     * @brief Default constructor.
     */
    SpringForceSpecFactory();

    /*!
     * @brief Destructor.
     */
    ~SpringForceSpecFactory();

    /*!
     * @brief Return the unique identifier used to specify the
     * StashableFactory object used by the StashableManager to extract
     * Stashable objects from data streams.
     */
    int getStashableID() const;

    /*!
     * @brief Set the unique identifier used to specify the
     * StashableFactory object used by the StashableManager to extract
     * Stashable objects from data streams.
     */
    void setStashableID(
        const int stashable_id);

    /*!
     * @brief Build a Stashable object by unpacking data from the
     * input stream.
     */
    SAMRAI::tbox::Pointer<Stashable> unpackStream(
        SAMRAI::tbox::AbstractStream& stream,
        const SAMRAI::hier::IntVector<NDIM>& offset);

private:
    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * @param from The value to copy to this object.
     */
    SpringForceSpecFactory(
        const SpringForceSpecFactory& from);

    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * @param that The value to assign to this object.
     *
     * @return A reference to this object.
     */
    SpringForceSpecFactory& operator=(
        const SpringForceSpecFactory& that);

    /*
     * The stashable ID for this object type.
     */
    static int s_stashable_id;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "SpringForceSpecFactory.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_SpringForceSpecFactory
