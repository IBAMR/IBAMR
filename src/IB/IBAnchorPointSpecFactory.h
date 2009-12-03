#ifndef included_IBAnchorPointSpecFactory
#define included_IBAnchorPointSpecFactory

// Filename: IBAnchorPointSpecFactory.h
// Last modified: <18.Aug.2008 13:35:58 boyce@dm-linux.maths.gla.ac.uk>
// Created on 18 Aug 2008 by Boyce Griffith (boyce@dm-linux.maths.gla.ac.uk)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/Stashable.h>
#include <ibtk/StashableFactory.h>

// SAMRAI INCLUDES
#include <IntVector.h>
#include <tbox/AbstractStream.h>
#include <tbox/Pointer.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBAnchorPointSpecFactory is a factory class to rebuild
 * IBAnchorPointSpec objects from SAMRAI::tbox::AbstractStream data streams.
 */
class IBAnchorPointSpecFactory
    : public IBTK::StashableFactory
{
public:
    /*!
     * \brief Default constructor.
     */
    IBAnchorPointSpecFactory();

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBAnchorPointSpecFactory();

    /*!
     * \brief Return the unique identifier used to specify the
     * IBTK::StashableFactory object used by the IBTK::StashableManager to
     * extract Stashable objects from data streams.
     */
    virtual int
    getStashableID() const;

    /*!
     * \brief Set the unique identifier used to specify the
     * IBTK::StashableFactory object used by the IBTK::StashableManager to
     * extract Stashable objects from data streams.
     */
    virtual void
    setStashableID(
        const int stashable_id);

    /*!
     * \brief Build a IBTK::Stashable object by unpacking data from the input
     * stream.
     */
    virtual SAMRAI::tbox::Pointer<IBTK::Stashable>
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
    IBAnchorPointSpecFactory(
        const IBAnchorPointSpecFactory& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBAnchorPointSpecFactory&
    operator=(
        const IBAnchorPointSpecFactory& that);

    /*
     * The stashable ID for this object type.
     */
    static int s_stashable_id;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBAnchorPointSpecFactory.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAnchorPointSpecFactory
