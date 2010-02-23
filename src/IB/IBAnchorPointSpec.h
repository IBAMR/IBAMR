#ifndef included_IBAnchorPointSpec
#define included_IBAnchorPointSpec

// Filename: IBAnchorPointSpec.h
// Last modified: <15.Dec.2009 19:15:20 griffith@boyce-griffiths-mac-pro.local>
// Created on 18 Aug 2008 by Boyce Griffith (boyce@dm-linux.maths.gla.ac.uk)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/DynamicArena.h>
#include <ibtk/Stashable.h>

// SAMRAI INCLUDES
#include <tbox/AbstractStream.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBAnchorPointSpec is used to indicate that a particular node of
 * the curvilinear mesh is anchored in place.
 *
 * \note Anchored curvilinear mesh nodes are fixed in space and are not allowed
 * to spread force to the Cartesian grid.
 */
class IBAnchorPointSpec
    : public IBTK::Stashable
{
public:
    /*!
     * \brief Register this class and its factory class with the singleton
     * IBTK::StashableManager object.  This method must be called before any
     * IBAnchorPointSpec objects are created.
     *
     * \note This method is collective on all MPI processes.  This is done to
     * ensure that all processes employ the same stashable ID for the
     * IBAnchorPointSpec class.
     */
    static void
    registerWithStashableManager();

    /*!
     * \brief Returns a boolean indicating whether the class has been registered
     * with the singleton IBTK::StashableManager object.
     */
    static bool
    getIsRegisteredWithStashableManager();

    /*!
     * \brief Operator new.
     */
    static void*
    operator new(
        std::size_t size);

    /*!
     * \brief Operator delete.
     */
    static void
    operator delete(
        void* ptr,
        std::size_t size);

    /*!
     * \brief Default constructor.
     */
    IBAnchorPointSpec(
        const int node_idx=-1);

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBAnchorPointSpec();

    /*!
     * \return A const reference to the node index.
     */
    const int&
    getNodeIndex() const;

    /*!
     * \return A non-const reference to the node index.
     */
    int&
    getNodeIndex();

    /*!
     * \brief Return the unique identifier used to specify the IBTK::StashableFactory
     * object used by the IBTK::StashableManager to extract Stashable objects from
     * data streams.
     */
    virtual int
    getStashableID() const;

    /*!
     * \brief Return an upper bound on the amount of space required to pack the
     * object to a buffer.
     */
    virtual size_t
    getDataStreamSize() const;

    /*!
     * \brief Pack data into the output stream.
     */
    virtual void
    packStream(
        SAMRAI::tbox::AbstractStream& stream);

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBAnchorPointSpec(
        const IBAnchorPointSpec& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBAnchorPointSpec&
    operator=(
        const IBAnchorPointSpec& that);

    /*!
     * Indicates whether the factory has been registered with the
     * IBTK::StashableManager.
     */
    static bool s_registered_factory;

    /*!
     * The stashable ID for this object type.
     */
    static int s_stashable_id;

    /*!
     * Memory arena for allocating objects.
     */
    static IBTK::DynamicArena s_arena;

    /*!
     * The Lagrangian index of the anchored curvilinear mesh node.
     */
    int d_node_idx;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibamr/IBAnchorPointSpec.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAnchorPointSpec
