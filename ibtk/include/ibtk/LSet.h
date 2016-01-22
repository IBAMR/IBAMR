// Filename: LSet.h
// Created on 29 Feb 2004 by Boyce Griffith
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

#ifndef included_LSet
#define included_LSet

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <vector>

#include "IntVector.h"
#include "ibtk/LMarker.h"
#include "ibtk/LNode.h"
#include "ibtk/LNodeIndex.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Index;
} // namespace hier
namespace tbox
{
class AbstractStream;
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LSet provides inter-processor communications and database access
 * functionality to a collection of Lagrangian objects.
 */
template <class T>
class LSet : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief The container class.
     */
    typedef std::vector<SAMRAI::tbox::Pointer<T> > DataSet;

    /*!
     * \brief The type of object, T, stored in the collection.
     */
    typedef typename DataSet::value_type value_type;

    /*!
     * \brief Pointer to T.
     */
    typedef typename DataSet::pointer pointer;

    /*!
     * \brief Reference to T.
     */
    typedef typename DataSet::reference reference;

    /*!
     * \brief Const reference to T.
     */
    typedef typename DataSet::const_reference const_reference;

    /*!
     * \brief An unsigned integral type.
     */
    typedef typename DataSet::size_type size_type;

    /*!
     * \brief A signed integral type.
     */
    typedef typename DataSet::difference_type difference_type;

    /*!
     * \brief Iterator used to iterate through the set.
     */
    typedef typename DataSet::iterator iterator;

    /*!
     * \brief Const iterator used to iterate through the collection.
     */
    typedef typename DataSet::const_iterator const_iterator;

    /*!
     * \brief Default constructor.
     */
    LSet();

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     */
    LSet(const LSet& from);

    /*!
     * \brief Destructor.
     */
    ~LSet();

    /*!
     * \brief Assignment operator.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LSet& operator=(const LSet& that);

    /*!
     * \return A reference to the nth element of the set.
     */
    reference operator[](size_type n);

    /*!
     * \return A const reference to the nth element of the set.
     */
    const_reference operator[](size_type n) const;

    /*!
     * \return A const_iterator pointing to the beginning of the set of indices.
     */
    const_iterator begin() const;

    /*!
     * \return An iterator pointing to the beginning of the set of indices.
     */
    iterator begin();

    /*!
     * \return A const_iterator pointing to the end of the set of indices.
     */
    const_iterator end() const;

    /*!
     * \return An iterator pointing to the end of the set of indices.
     */
    iterator end();

    /*!
     * \return The size of the set.
     */
    size_type size() const;

    /*!
     * \return Whether the set is empty.
     */
    bool empty() const;

    /*!
     * \brief Insert a new element at the end (of the set).
     */
    void push_back(const value_type& value);

    /*!
     * \brief Inserts x before pos.
     */
    iterator insert(iterator pos, const typename LSet<T>::value_type& x);

    /*!
     * \brief Inserts the range [first,last) before pos.
     */
    template <class InputIterator>
    void insert(iterator pos, InputIterator first, InputIterator last);

    /*!
     * \brief Inserts n copies of x before pos.
     */
    void insert(iterator pos, size_type n, const typename LSet<T>::value_type& x);

    /*!
     * \brief Return a const reference to the set of data items.
     */
    const DataSet& getDataSet() const;

    /*!
     * \brief Return a non-const reference to the set of data items.
     */
    DataSet& getDataSet();

    /*!
     * \brief Reset the set of data items.
     */
    void setDataSet(const DataSet& set);

    /*!
     * \return A const reference to the periodic offset.
     *
     * \note If the LSet lives in cell i, the index of the source
     * object is src_index = i - offset.
     */
    const SAMRAI::hier::IntVector<NDIM>& getPeriodicOffset() const;

    /*!
     * \brief Set the value of the periodic offset.
     *
     * \note If the LSet lives in cell i, the index of the source
     * object is src_index = i - offset.
     */
    void setPeriodicOffset(const SAMRAI::hier::IntVector<NDIM>& offset);

    /*!
     * \brief Copy data from the source.
     *
     * \note The index of the destination object is src_index + src_offset.
     */
    void copySourceItem(const SAMRAI::hier::Index<NDIM>& src_index,
                        const SAMRAI::hier::IntVector<NDIM>& src_offset,
                        const LSet& src_item);

    /*!
     * \brief Return an upper bound on the amount of space required to pack the
     * object to a buffer.
     */
    size_t getDataStreamSize() const;

    /*!
     * \brief Pack data into the output stream.
     */
    void packStream(SAMRAI::tbox::AbstractStream& stream);

    /*!
     * \brief Unpack data from the input stream.
     */
    void unpackStream(SAMRAI::tbox::AbstractStream& stream, const SAMRAI::hier::IntVector<NDIM>& offset);

    /*!
     * \brief Pack data into a database.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> database);

    /*!
     * \brief Unpack data from a database.
     */
    void getFromDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> database);

private:
    /*!
     * \brief The collection of data items.
     */
    typename LSet<T>::DataSet d_set;

    /*!
     * \brief The periodic offset.
     */
    SAMRAI::hier::IntVector<NDIM> d_offset;
};
} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/LSet-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LSet
