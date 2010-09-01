// Filename: IBAnchorPointSpec.h
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

#ifndef included_IBAnchorPointSpec
#define included_IBAnchorPointSpec

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
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
     * The Lagrangian index of the anchored curvilinear mesh node.
     */
    int d_node_idx;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibamr/IBAnchorPointSpec.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAnchorPointSpec
