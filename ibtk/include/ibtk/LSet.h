// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#ifndef included_IBTK_LSet
#define included_IBTK_LSet

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/LMarker.h"
#include "ibtk/LNode.h"
#include "ibtk/LNodeIndex.h"

#include "IntVector.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

#include <vector>

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
    using DataSet = std::vector<SAMRAI::tbox::Pointer<T> >;

    /*!
     * \brief The type of object, T, stored in the collection.
     */
    using value_type = typename DataSet::value_type;

    /*!
     * \brief Pointer to T.
     */
    using pointer = typename DataSet::pointer;

    /*!
     * \brief Reference to T.
     */
    using reference = typename DataSet::reference;

    /*!
     * \brief Const reference to T.
     */
    using const_reference = typename DataSet::const_reference;

    /*!
     * \brief An unsigned integral type.
     */
    using size_type = typename DataSet::size_type;

    /*!
     * \brief A signed integral type.
     */
    using difference_type = typename DataSet::difference_type;

    /*!
     * \brief Iterator used to iterate through the set.
     */
    using iterator = typename DataSet::iterator;

    /*!
     * \brief Const iterator used to iterate through the collection.
     */
    using const_iterator = typename DataSet::const_iterator;

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
    virtual ~LSet();

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

#endif //#ifndef included_IBTK_LSet
