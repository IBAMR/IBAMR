// Filename: IndexSet.h
// Created on 30 Apr 2010 by Boyce Griffith
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

#ifndef included_IndexSet
#define included_IndexSet

/////////////////////////////// INCLUDES /////////////////////////////////////

// C++ STDLIB INCLUDES
#include <set>

// SAMRAI INCLUDES
#include <Index.h>
#include <IntVector.h>
#include <tbox/AbstractStream.h>
#include <tbox/Database.h>
#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

/////////////////////////////// FORWARD DECLARATIONS /////////////////////////

namespace IBTK
{
class LDataManager;
}// namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class IndexSet provides inter-processor communications and database
 * access functionality to a simple index set class which is intended to be used
 * with class SAMRAI::pdat::IndexData.
 *
 * \note This class meets the required specification for use with the templated
 * class SAMRAI::pdat::IndexData.
 */
class IndexSet
    : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     */
    IndexSet();

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     */
    IndexSet(
        const IndexSet& from);

    /*!
     * \brief Destructor.
     */
    ~IndexSet();

    /*!
     * \brief Assignment operator.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IndexSet&
    operator=(
        const IndexSet& that);

    /*!
     * \brief The index set wrapped by this class.
     */
    std::set<int> index_set;

    /*!
     * \brief Copy data from the source.
     *
     * \note The index of the destination object is src_index + src_offset.
     */
    void
    copySourceItem(
        const SAMRAI::hier::Index<NDIM>& src_index,
        const SAMRAI::hier::IntVector<NDIM>& src_offset,
        const IndexSet& src_item);

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
        SAMRAI::tbox::AbstractStream& stream) const;

    /*!
     * \brief Unpack data from the input stream.
     */
    void
    unpackStream(
        SAMRAI::tbox::AbstractStream& stream,
        const SAMRAI::hier::IntVector<NDIM>& offset);

    /*!
     * \brief Pack data into a database.
     */
    void
    putToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& database) const;

    /*!
     * \brief Unpack data from a database.
     */
    void
    getFromDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& database);
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibtk/IndexSet.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IndexSet
