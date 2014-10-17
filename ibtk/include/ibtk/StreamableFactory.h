// Filename: StreamableFactory.h
// Created on 14 Jun 2004 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_StreamableFactory
#define included_StreamableFactory

/////////////////////////////// INCLUDES /////////////////////////////////////

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
    StreamableFactory();

    /*!
     * \brief Virtual destructor.
     */
    virtual ~StreamableFactory();

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
    StreamableFactory(const StreamableFactory& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StreamableFactory& operator=(const StreamableFactory& that);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_StreamableFactory
