#ifndef included_Stashable
#define included_Stashable

// Filename: Stashable.h
// Created on 14 Jun 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)
// Last modified: <04.Oct.2006 19:49:45 boyce@boyce-griffiths-powerbook-g4-15.local>

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <tbox/AbstractStream.h>
#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * @brief Interface to facilitate object packing for
 * tbox::AbstractStream based communication.
 *
 * NOTE: The StashableManager should be used for all communications
 * and storage operations.  Do not directly use Stashable objects.
 *
 * @see StashableFactory
 * @see StashableManager
 */
class Stashable
    : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * @brief Default empty constructor.
     */
    Stashable();

    /*!
     * @brief Virtual destructor.
     */
    virtual ~Stashable();

    /*!
     * @brief Return the unique identifier used to specify the
     * StashableFactory object used by the StashableManager to extract
     * Stashable objects from data streams.
     */
    virtual int getStashableID() const = 0;

    /*!
     * @brief Return an upper bound on the amount of space required to
     * pack the object to a buffer.
     */
    virtual size_t getDataStreamSize() const = 0;

    /*!
     * @brief Pack data into the output stream.
     */
    virtual void packStream(
        SAMRAI::tbox::AbstractStream& stream) = 0;

private:
    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * @param from The value to copy to this object.
     */
    Stashable(
        const Stashable& from);

    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * @param that The value to assign to this object.
     *
     * @return A reference to this object.
     */
    Stashable& operator=(
        const Stashable& that);
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/Stashable.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_Stashable
