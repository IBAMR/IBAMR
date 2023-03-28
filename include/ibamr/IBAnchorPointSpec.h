// ---------------------------------------------------------------------
//
// Copyright (c) 2008 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_IBAnchorPointSpec
#define included_IBAMR_IBAnchorPointSpec

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibtk/Streamable.h"
#include "ibtk/StreamableFactory.h"

#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class IntVector;
} // namespace hier
namespace tbox
{
class AbstractStream;
} // namespace tbox
} // namespace SAMRAI

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
class IBAnchorPointSpec : public IBTK::Streamable
{
public:
    /*!
     * \brief Register this class and its factory class with the singleton
     * IBTK::StreamableManager object.  This method must be called before any
     * IBAnchorPointSpec objects are created.
     *
     * \note This method is collective on all MPI processes.  This is done to
     * ensure that all processes employ the same class ID for the
     * IBAnchorPointSpec class.
     */
    static void registerWithStreamableManager();

    /*!
     * \brief Returns a boolean indicating whether the class has been registered
     * with the singleton IBTK::StreamableManager object.
     */
    static bool getIsRegisteredWithStreamableManager();

    /*!
     * The unique class ID for this object type assigned by the
     * IBTK::StreamableManager.
     */
    static int STREAMABLE_CLASS_ID;

    /*!
     * \brief Default constructor.
     */
    IBAnchorPointSpec(int node_idx = IBTK::invalid_index);

    /*!
     * \brief Destructor.
     */
    ~IBAnchorPointSpec();

    /*!
     * \return A const reference to the node index.
     */
    const int& getNodeIndex() const;

    /*!
     * \return A non-const reference to the node index.
     */
    int& getNodeIndex();

    /*!
     * \brief Return the unique identifier used to specify the
     * IBTK::StreamableFactory object used by the IBTK::StreamableManager to
     * extract Streamable objects from data streams.
     */
    int getStreamableClassID() const override;

    /*!
     * \brief Return an upper bound on the amount of space required to pack the
     * object to a buffer.
     */
    size_t getDataStreamSize() const override;

    /*!
     * \brief Pack data into the output stream.
     */
    void packStream(SAMRAI::tbox::AbstractStream& stream) override;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBAnchorPointSpec(const IBAnchorPointSpec& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBAnchorPointSpec& operator=(const IBAnchorPointSpec& that) = delete;

    /*!
     * The Lagrangian index of the anchored curvilinear mesh node.
     */
    int d_node_idx;

    /*!
     * \brief A factory class to rebuild IBAnchorPointSpec objects from
     * SAMRAI::tbox::AbstractStream data streams.
     */
    class Factory : public IBTK::StreamableFactory
    {
    public:
        /*!
         * \brief Destructor.
         */
        ~Factory() = default;

        /*!
         * \brief Return the unique identifier used to specify the
         * IBTK::StreamableFactory object used by the IBTK::StreamableManager to
         * extract IBAnchorPointSpec objects from data streams.
         */
        int getStreamableClassID() const override;

        /*!
         * \brief Set the unique identifier used to specify the
         * IBTK::StreamableFactory object used by the IBTK::StreamableManager to
         * extract IBAnchorPointSpec objects from data streams.
         */
        void setStreamableClassID(int class_id) override;

        /*!
         * \brief Build an IBAnchorPointSpec object by unpacking data from the data
         * stream.
         */
        SAMRAI::tbox::Pointer<IBTK::Streamable> unpackStream(SAMRAI::tbox::AbstractStream& stream,
                                                             const SAMRAI::hier::IntVector<NDIM>& offset) override;

    private:
        /*!
         * \brief Default constructor.
         */
        Factory();

        /*!
         * \brief Copy constructor.
         *
         * \note This constructor is not implemented and should not be used.
         *
         * \param from The value to copy to this object.
         */
        Factory(const Factory& from) = delete;

        /*!
         * \brief Assignment operator.
         *
         * \note This operator is not implemented and should not be used.
         *
         * \param that The value to assign to this object.
         *
         * \return A reference to this object.
         */
        Factory& operator=(const Factory& that) = delete;

        friend class IBAnchorPointSpec;
    };
    using IBAnchorPointSpecFactory = IBAnchorPointSpec::Factory;
};
} // namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibamr/private/IBAnchorPointSpec-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_IBAnchorPointSpec
