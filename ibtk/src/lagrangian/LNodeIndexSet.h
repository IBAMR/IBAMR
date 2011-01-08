// Filename: LNodeIndexSet.h
// Created on 29 Feb 2004 by Boyce Griffith
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

#ifndef included_LNodeIndexSet
#define included_LNodeIndexSet

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/LNodeIndex.h>

// SAMRAI INCLUDES
#include <Index.h>
#include <IntVector.h>
#include <tbox/AbstractStream.h>
#include <tbox/Database.h>
#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// FORWARD DECLARATIONS /////////////////////////

namespace IBTK
{
class LDataManager;
}// namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LNodeIndexSet provides inter-processor communications and
 * database access functionality to a collection of LNodexIndex objects.
 *
 * \note This class meets the required specification for use with the templated
 * class SAMRAI::pdat::IndexData.
 */
class LNodeIndexSet
    : public SAMRAI::tbox::DescribedClass
{
public:
    friend class LDataManager;

    /*!
     * \brief The type of the collection.
     */
    typedef std::vector<SAMRAI::tbox::Pointer<LNodeIndex> > IndexSet;

    /*!
     * \brief The type of object, T, stored in the collection.
     */
    typedef IndexSet::value_type value_type;

    /*!
     * \brief SAMRAI::tbox::Pointer to T.
     */
    typedef IndexSet::pointer pointer;

    /*!
     * \brief Reference to T.
     */
    typedef IndexSet::reference reference;

    /*!
     * \brief Const reference to T.
     */
    typedef IndexSet::const_reference const_reference;

    /*!
     * \brief An unsigned integral type.
     */
    typedef IndexSet::size_type size_type;

    /*!
     * \brief A signed integral type.
     */
    typedef IndexSet::difference_type difference_type;

    /*!
     * \brief Iterator used to iterate through the set.
     */
    typedef IndexSet::iterator iterator;

    /*!
     * \brief Const iterator used to iterate through the collection.
     */
    typedef IndexSet::const_iterator const_iterator;

    /*!
     * \brief Default constructor.
     */
    LNodeIndexSet();

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     */
    LNodeIndexSet(
        const LNodeIndexSet& from);

    /*!
     * \brief Destructor.
     */
    ~LNodeIndexSet();

    /*!
     * \brief Assignment operator.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LNodeIndexSet&
    operator=(
        const LNodeIndexSet& that);

    /*!
     * \return A const_iterator pointing to the beginning of the set of indices.
     */
    const_iterator
    begin() const;

    /*!
     * \return An iterator pointing to the beginning of the set of indices.
     */
    iterator
    begin();

    /*!
     * \return A const_iterator pointing to the end of the set of indices.
     */
    const_iterator
    end() const;

    /*!
     * \return An iterator pointing to the end of the set of indices.
     */
    iterator
    end();

    /*!
     * \return The size of the set.
     */
    size_type
    size() const;

    /*!
     * \return Whether the set is empty.
     */
    bool
    empty() const;

    /*!
     * \brief Insert a new element at the end (of the set).
     */
    void
    push_back(
        const value_type& value);

    /*!
     * \return A const reference to the periodic offset.
     *
     * \note If the LNodeIndexSet lives in cell i, the index of the source
     * object is src_index = i - offset.
     */
    const SAMRAI::hier::IntVector<NDIM>&
    getPeriodicOffset() const;

    /*!
     * \brief Set the value of the periodic offset.
     *
     * \note If the LNodeIndexSet lives in cell i, the index of the source
     * object is src_index = i - offset.
     */
    void
    setPeriodicOffset(
        const SAMRAI::hier::IntVector<NDIM>& offset);

    /*!
     * \brief Copy data from the source.
     *
     * \note The index of the destination object is src_index + src_offset.
     */
    void
    copySourceItem(
        const SAMRAI::hier::Index<NDIM>& src_index,
        const SAMRAI::hier::IntVector<NDIM>& src_offset,
        const LNodeIndexSet& src_item);

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

private:
    /*!
     * \brief The collection of indices.
     */
    IndexSet d_set;

    /*!
     * \brief The periodic offset.
     */
    SAMRAI::hier::IntVector<NDIM> d_offset;

    /*!
     * \brief Misc. nested structs and classes.
     */
    struct LNodeIndexGetDataStreamSizeSum
        : std::binary_function<size_t,SAMRAI::tbox::Pointer<LNodeIndex>,size_t>
    {
        inline size_t
        operator()(
            size_t size_so_far,
            const SAMRAI::tbox::Pointer<LNodeIndex>& index) const
            {
                return size_so_far+index->getDataStreamSize();
            }
    };

    class LNodeIndexPackStream
        : public std::unary_function<SAMRAI::tbox::Pointer<LNodeIndex>,void>
    {
    public:
        inline
        LNodeIndexPackStream(
            SAMRAI::tbox::AbstractStream* const stream)
            : d_stream(stream)
            {
                return;
            }

        inline void
        operator()(
            const SAMRAI::tbox::Pointer<LNodeIndex>& index) const
            {
                index->packStream(*d_stream);
                return;
            }
    private:
        SAMRAI::tbox::AbstractStream* const d_stream;
    };

    class LNodeIndexUnpackStream
        : public std::unary_function<void,SAMRAI::tbox::Pointer<LNodeIndex> >
    {
    public:
        inline
        LNodeIndexUnpackStream(
            SAMRAI::tbox::AbstractStream* const stream,
            const SAMRAI::hier::IntVector<NDIM>& offset)
            : d_stream(stream),
              d_offset(offset)
            {
                return;
            }

        inline SAMRAI::tbox::Pointer<LNodeIndex>
        operator()() const
            {
                SAMRAI::tbox::Pointer<LNodeIndex> index_out = new LNodeIndex();
                index_out->unpackStream(*d_stream,d_offset);
                return index_out;
            }
    private:
        SAMRAI::tbox::AbstractStream* const d_stream;
        const SAMRAI::hier::IntVector<NDIM>& d_offset;
    };
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibtk/LNodeIndexSet.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LNodeIndexSet
