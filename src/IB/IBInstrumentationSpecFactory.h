#ifndef included_IBInstrumentationSpecFactory
#define included_IBInstrumentationSpecFactory

// Filename: IBInstrumentationSpecFactory.h
// Last modified: <11.Jun.2007 17:41:49 griffith@box221.cims.nyu.edu>
// Created on 11 Jun 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

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
 * \brief Class IBInstrumentationSpecFactory is a factory class to rebuild
 * IBInstrumentationSpec objects from SAMRAI::tbox::AbstractStream data streams.
 */
class IBInstrumentationSpecFactory
    : public StashableFactory
{
public:
    /*!
     * \brief Default constructor.
     */
    IBInstrumentationSpecFactory();

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBInstrumentationSpecFactory();

    /*!
     * \brief Return the unique identifier used to specify the StashableFactory
     * object used by the StashableManager to extract Stashable objects from
     * data streams.
     */
    virtual int
    getStashableID() const;

    /*!
     * \brief Set the unique identifier used to specify the StashableFactory
     * object used by the StashableManager to extract Stashable objects from
     * data streams.
     */
    virtual void
    setStashableID(
        const int stashable_id);

    /*!
     * \brief Build a Stashable object by unpacking data from the input stream.
     */
    virtual SAMRAI::tbox::Pointer<Stashable>
    unpackStream(
        SAMRAI::tbox::AbstractStream& stream,
        const SAMRAI::hier::IntVector<NDIM>& offset);

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBInstrumentationSpecFactory(
        const IBInstrumentationSpecFactory& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBInstrumentationSpecFactory&
    operator=(
        const IBInstrumentationSpecFactory& that);

    /*
     * The stashable ID for this object type.
     */
    static int s_stashable_id;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBInstrumentationSpecFactory.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBInstrumentationSpecFactory
