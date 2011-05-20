// Filename: IBSourceSpec.h
// Created on 28 Apr 2011 by Boyce Griffith
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

#ifndef included_IBSourceSpec
#define included_IBSourceSpec

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/Streamable.h>

// SAMRAI INCLUDES
#include <tbox/AbstractStream.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// FORWARD DECLARATION //////////////////////////

namespace IBAMR
{
class IBSourceSpecFactory;
}

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBSourceSpec encapsulates the data required to initialize
 * distributed internal sources and sinks.
 */
class IBSourceSpec
    : public IBTK::Streamable
{
public:
    friend class IBSourceSpecFactory;

    /*!
     * \brief Register this class and its factory class with the singleton
     * IBTK::StreamableManager object.  This method must be called before any
     * IBSourceSpec objects are created.
     *
     * \note This method is collective on all MPI processes.  This is done to
     * ensure that all processes employ the same class ID for the
     * IBSourceSpec class.
     */
    static void
    registerWithStreamableManager();

    /*!
     * \brief Returns a boolean indicating whether the class has been registered
     * with the singleton IBTK::StreamableManager object.
     */
    static bool
    getIsRegisteredWithStreamableManager();

    /*!
     * \brief Set the number of internal sources and sinks on the specified
     * level of the patch hierarchy.
     */
    static void
    setNumSources(
        const int ln,
        const unsigned int num_sources);

    /*!
     * \brief Get the number of internal sources and sinks on the specified
     * level of the patch hierarchy.
     */
    static unsigned int
    getNumSources(
        const int ln);

    /*!
     * \brief Set the names of the internal sources and sinks on the specified
     * level of the patch hierarchy.
     */
    static void
    setSourceNames(
        const int ln,
        const std::vector<std::string>& names);

    /*!
     * \brief Get the names of the internal sources and sinks on the specified
     * level of the patch hierarchy.
     */
    static const std::vector<std::string>&
    getSourceNames(
        const int ln);

    /*!
     * \brief Set the sizes of the internal sources and sinks on the specified
     * level of the patch hierarchy.
     */
    static void
    setSourceRadii(
        const int ln,
        const std::vector<double>& radii);

    /*!
     * \brief Get the sizes of the internal sources and sinks on the specified
     * level of the patch hierarchy.
     */
    static const std::vector<double>&
    getSourceRadii(
        const int ln);

    /*!
     * \brief Default constructor.
     */
    IBSourceSpec(
        const int master_idx=-1,
        const int source_idx=-1);

    /*!
     * \brief Destructor.
     */
    ~IBSourceSpec();

    /*!
     * \return A const reference to the master node index.
     */
    const int&
    getMasterNodeIndex() const;

    /*!
     * \return A non-const reference to the master node index.
     */
    int&
    getMasterNodeIndex();

    /*!
     * \return A const reference to the source index associated with the master
     * node.
     */
    const int&
    getSourceIndex() const;

    /*!
     * \return A non-const reference to the source index associated with the
     * master node.
     */
    int&
    getSourceIndex();

    /*!
     * \brief Return the unique identifier used to specify the
     * IBTK::StreamableFactory object used by the IBTK::StreamableManager to
     * extract Streamable objects from data streams.
     */
    int
    getStreamableClassID() const;

    /*!
     * \brief Return an upper bound on the amount of space required to pack the
     * object to a buffer.
     */
    size_t
    getDataStreamSize() const;

    /*!
     * \brief Pack data into the output stream.
     */
    void
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
    IBSourceSpec(
        const IBSourceSpec& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBSourceSpec&
    operator=(
        const IBSourceSpec& that);

    /*!
     * Indicates whether the factory has been registered with the
     * IBTK::StreamableManager.
     */
    static bool s_registered_factory;

    /*!
     * The class ID for this object type assigned by the
     * IBTK::StreamableManager.
     */
    static int s_class_id;

    /*!
     * The numbers of sources/sinks on each level of the patch hierarchy.
     */
    static std::vector<unsigned int> s_num_sources;

    /*!
     * The names of the sources and sinks.
     */
    static std::vector<std::vector<std::string> > s_source_names;

    /*!
     * The sizes of the sources and sinks.
     */
    static std::vector<std::vector<double> > s_source_radii;

    /*!
     * Data required to define the source.
     */
    int d_master_idx, d_source_idx;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

#include "IBSourceSpec.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBSourceSpec
