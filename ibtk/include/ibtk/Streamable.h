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

#ifndef included_IBTK_Streamable
#define included_IBTK_Streamable

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/ibtk_utilities.h"

#include "tbox/DescribedClass.h"

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
 * \brief Class Streamable is an abstract interface for objects that can be
 * packed into SAMRAI::tbox::AbstractStream data streams.
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
 * \see StreamableFactory
 * \see StreamableManager
 */
class Streamable : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default empty constructor.
     */
    Streamable() = default;

    /*!
     * \brief Virtual destructor.
     */
    virtual ~Streamable() = default;

    /*!
     * \brief Return the unique class identifier used to specify the
     * StreamableFactory object used by the StreamableManager to extract
     * Streamable objects from data streams.
     */
    virtual int getStreamableClassID() const = 0;

    /*!
     * \brief Return an upper bound on the amount of space required to pack the
     * object to a buffer.
     */
    virtual size_t getDataStreamSize() const = 0;

    /*!
     * \brief Pack data into the output stream.
     */
    virtual void packStream(SAMRAI::tbox::AbstractStream& stream) = 0;

    /*!
     * \brief Indicate that the Streamable object has been shifted across a
     * periodic boundary.
     *
     * \note A default empty implementation is provided.
     */
    virtual void registerPeriodicShift(const SAMRAI::hier::IntVector<NDIM>& offset, const Vector& displacement);

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    Streamable(const Streamable& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    Streamable& operator=(const Streamable& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_Streamable
