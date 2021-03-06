// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2020 by the IBAMR developers
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

#ifndef included_IBTK_StreamableFactory
#define included_IBTK_StreamableFactory

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

namespace IBTK
{
class Streamable;
} // namespace IBTK
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

namespace IBTK
{
/*!
 * \brief Class StreamableFactory is an abstract interface for classes that can
 * unpack particular concrete Streamable objects from
 * SAMRAI::tbox::AbstractStream data streams.
 *
 * \note Each concrete Streamable class must have a corresponding concrete
 * StreamableFactory class.  Classes that implement the Streamable interface are
 * able to pack themselves into a stream; the corresponding StreamableFactory
 * class is used to unpack that data and construct the corresponding Streamable
 * object.
 *
 * \note Class StreamableManager should be used for all communications and
 * storage operations.
 *
 * \see Streamable
 * \see StreamableManager
 */
class StreamableFactory : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default empty constructor.
     */
    StreamableFactory() = default;

    /*!
     * \brief Virtual destructor.
     */
    virtual ~StreamableFactory() = default;

    /*!
     * \brief Return the unique class identifier used to specify the
     * StreamableFactory object used by the StreamableManager to extract
     * Streamable objects from data streams.
     */
    virtual int getStreamableClassID() const = 0;

    /*!
     * \brief Set the unique identifier used to specify the StreamableFactory
     * object used by the StreamableManager to extract Streamable objects from
     * data streams.
     */
    virtual void setStreamableClassID(int class_id) = 0;

    /*!
     * \brief Build a Streamable object by unpacking data from the data stream.
     */
    virtual SAMRAI::tbox::Pointer<Streamable> unpackStream(SAMRAI::tbox::AbstractStream& stream,
                                                           const SAMRAI::hier::IntVector<NDIM>& offset) = 0;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StreamableFactory(const StreamableFactory& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StreamableFactory& operator=(const StreamableFactory& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_StreamableFactory
