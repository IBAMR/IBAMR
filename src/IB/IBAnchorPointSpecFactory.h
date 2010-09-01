// Filename: IBAnchorPointSpecFactory.h
// Created on 18 Aug 2008 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#ifndef included_IBAnchorPointSpecFactory
#define included_IBAnchorPointSpecFactory

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
