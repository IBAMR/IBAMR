#ifndef included_IBMovingTargetPointForceSpecFactory
#define included_IBMovingTargetPointForceSpecFactory

// Filename: IBMovingTargetPointForceSpecFactory.h
// Last modified: <14.Aug.2008 13:02:36 boyce@dm-linux.maths.gla.ac.uk>
// Created on 14 Aug 2008 by Boyce Griffith (boyce@dm-linux.maths.gla.ac.uk)

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
 * \brief Class IBMovingTargetPointForceSpecFactory is a factory class to
 * rebuild IBMovingTargetPointForceSpec objects from
 * SAMRAI::tbox::AbstractStream data streams.
 */
class IBMovingTargetPointForceSpecFactory
    : public IBTK::StashableFactory
{
public:
    /*!
     * \brief Default constructor.
     */
    IBMovingTargetPointForceSpecFactory();

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBMovingTargetPointForceSpecFactory();

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
    IBMovingTargetPointForceSpecFactory(
        const IBMovingTargetPointForceSpecFactory& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBMovingTargetPointForceSpecFactory&
    operator=(
        const IBMovingTargetPointForceSpecFactory& that);

    /*
     * The stashable ID for this object type.
     */
    static int s_stashable_id;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBMovingTargetPointForceSpecFactory.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBMovingTargetPointForceSpecFactory
