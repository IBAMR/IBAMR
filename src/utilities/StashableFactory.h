#ifndef included_StashableFactory
#define included_StashableFactory

// Filename: StashableFactory.h
// Created on 14 Jun 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)
// Last modified: <16.Nov.2006 00:16:44 boyce@bigboy.nyconnect.com>

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/Stashable.h>

// SAMRAI INCLUDES
#include <IntVector.h>
#include <tbox/AbstractStream.h>
#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Interface to facilitate object unpacking for
 * tbox::AbstractStream based communication.
 *
 * \note Class StashableManager should be used for all communications
 * and storage operations.  Do not directly use StashableFactory
 * objects.
 *
 * \see Stashable
 * \see StashableManager
 */
class StashableFactory
    : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default empty constructor.
     */
    StashableFactory();

    /*!
     * \brief Virtual destructor.
     */
    virtual ~StashableFactory();

    /*!
     * \brief Return the unique identifier used to specify the
     * StashableFactory object used by the StashableManager to extract
     * Stashable objects from data streams.
     */
    virtual int getStashableID() const = 0;

    /*!
     * \brief Set the unique identifier used to specify the
     * StashableFactory object used by the StashableManager to extract
     * Stashable objects from data streams.
     */
    virtual void setStashableID(
        const int stashable_id) = 0;

    /*!
     * \brief Build a Stashable object by unpacking data from the
     * input stream.
     */
    virtual SAMRAI::tbox::Pointer<Stashable> unpackStream(
        SAMRAI::tbox::AbstractStream& stream,
        const SAMRAI::hier::IntVector<NDIM>& offset) = 0;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    StashableFactory(
        const StashableFactory& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StashableFactory& operator=(
        const StashableFactory& that);
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/StashableFactory.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_StashableFactory
