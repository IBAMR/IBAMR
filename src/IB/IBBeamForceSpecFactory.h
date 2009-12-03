#ifndef included_IBBeamForceSpecFactory
#define included_IBBeamForceSpecFactory

// Filename: IBBeamForceSpecFactory.h
// Last modified: <12.Mar.2008 22:38:43 griffith@box221.cims.nyu.edu>
// Created on 22 Mar 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

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
 * \brief Class IBBeamForceSpecFactory is a factory class to rebuild
 * IBBeamForceSpec objects from SAMRAI::tbox::AbstractStream data streams.
 */
class IBBeamForceSpecFactory
    : public IBTK::StashableFactory
{
public:
    /*!
     * \brief Default constructor.
     */
    IBBeamForceSpecFactory();

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBBeamForceSpecFactory();

    /*!
     * \brief Return the unique identifier used to specify the IBTK::StashableFactory
     * object used by the IBTK::StashableManager to extract Stashable objects from
     * data streams.
     */
    virtual int
    getStashableID() const;

    /*!
     * \brief Set the unique identifier used to specify the IBTK::StashableFactory
     * object used by the IBTK::StashableManager to extract Stashable objects from
     * data streams.
     */
    virtual void
    setStashableID(
        const int stashable_id);

    /*!
     * \brief Build a IBTK::Stashable object by unpacking data from the input stream.
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
    IBBeamForceSpecFactory(
        const IBBeamForceSpecFactory& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBBeamForceSpecFactory&
    operator=(
        const IBBeamForceSpecFactory& that);

    /*
     * The stashable ID for this object type.
     */
    static int s_stashable_id;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBBeamForceSpecFactory.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBBeamForceSpecFactory
