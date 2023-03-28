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

#ifndef included_IBTK_LSet_inl_h
#define included_IBTK_LSet_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/FixedSizedStream.h"
#include "ibtk/LSet.h"

#include "tbox/Database.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

template <class T>
inline LSet<T>::LSet() : d_set(), d_offset(0)
{
    // intentionally blank
    return;
} // LSet

template <class T>
inline LSet<T>::LSet(const LSet<T>& from) : d_set(from.d_set), d_offset(from.d_offset)
{
    // intentionally blank
    return;
} // LSet

template <class T>
inline LSet<T>::~LSet()
{
    d_set.clear();
    return;
} // ~LSet

template <class T>
inline LSet<T>&
LSet<T>::operator=(const LSet<T>& that)
{
    if (this == &that) return *this; // check for self-assignment

    d_set = that.d_set;
    d_offset = that.d_offset;
    return *this;
} // operator=

template <class T>
inline typename LSet<T>::reference
LSet<T>::operator[](typename LSet<T>::size_type n)
{
    return d_set[n];
} // operator[]

template <class T>
inline typename LSet<T>::const_reference
LSet<T>::operator[](typename LSet<T>::size_type n) const
{
    return d_set[n];
} // operator[]

template <class T>
inline typename LSet<T>::const_iterator
LSet<T>::begin() const
{
    return d_set.begin();
} // begin

template <class T>
inline typename LSet<T>::iterator
LSet<T>::begin()
{
    return d_set.begin();
} // begin

template <class T>
inline typename LSet<T>::const_iterator
LSet<T>::end() const
{
    return d_set.end();
} // end

template <class T>
inline typename LSet<T>::iterator
LSet<T>::end()
{
    return d_set.end();
} // end

template <class T>
inline typename LSet<T>::size_type
LSet<T>::size() const
{
    return d_set.size();
} // size

template <class T>
inline bool
LSet<T>::empty() const
{
    return d_set.empty();
} // empty

template <class T>
inline void
LSet<T>::push_back(const LSet<T>::value_type& value)
{
    d_set.push_back(value);
    return;
} // push_back

template <class T>
inline typename LSet<T>::iterator
LSet<T>::insert(typename LSet<T>::iterator pos, const typename LSet<T>::value_type& x)
{
    return d_set.insert(pos, x);
} // insert

template <class T>
template <class InputIterator>
inline void
LSet<T>::insert(typename LSet<T>::iterator pos, InputIterator first, InputIterator last)
{
    d_set.insert(pos, first, last);
    return;
} // insert

template <class T>
inline void
LSet<T>::insert(typename LSet<T>::iterator pos, typename LSet<T>::size_type n, const typename LSet<T>::value_type& x)
{
    d_set.insert(pos, n, x);
    return;
} // insert

template <class T>
inline const typename LSet<T>::DataSet&
LSet<T>::getDataSet() const
{
    return d_set;
} // getDataSet

template <class T>
inline typename LSet<T>::DataSet&
LSet<T>::getDataSet()
{
    return d_set;
} // getDataSet

template <class T>
inline void
LSet<T>::setDataSet(const typename LSet<T>::DataSet& set)
{
    d_set = set;
    return;
} // setDataSet

template <class T>
inline const SAMRAI::hier::IntVector<NDIM>&
LSet<T>::getPeriodicOffset() const
{
    return d_offset;
} // getPeriodicOffset

template <class T>
inline void
LSet<T>::setPeriodicOffset(const SAMRAI::hier::IntVector<NDIM>& offset)
{
    d_offset = offset;
    return;
} // setPeriodicOffset

template <class T>
inline void
LSet<T>::copySourceItem(const SAMRAI::hier::Index<NDIM>& /*src_index*/,
                        const SAMRAI::hier::IntVector<NDIM>& src_offset,
                        const LSet<T>& src_item)
{
    d_set = src_item.d_set;
    d_offset = src_offset;
    return;
} // copySourceItem

template <class T>
inline size_t
LSet<T>::getDataStreamSize() const
{
    size_t size = SAMRAI::tbox::AbstractStream::sizeofInt();
    for (unsigned int k = 0; k < d_set.size(); ++k)
    {
        size += d_set[k]->getDataStreamSize();
    }
    return size;
} // getDataStreamSize

template <class T>
inline void
LSet<T>::packStream(SAMRAI::tbox::AbstractStream& stream)
{
    int num_idx = static_cast<int>(d_set.size());
    stream.pack(&num_idx, 1);
    for (unsigned int k = 0; k < d_set.size(); ++k)
    {
        d_set[k]->packStream(stream);
    }
    return;
} // packStream

template <class T>
inline void
LSet<T>::unpackStream(SAMRAI::tbox::AbstractStream& stream, const SAMRAI::hier::IntVector<NDIM>& offset)
{
    d_offset = offset;
    int num_idx;
    stream.unpack(&num_idx, 1);
    d_set.resize(num_idx);
    for (unsigned int k = 0; k < d_set.size(); ++k)
    {
        d_set[k] = new T(stream, offset);
    }
    typename LSet<T>::DataSet(d_set).swap(d_set); // trim-to-fit
    return;
} // unpackStream

template <class T>
inline void
LSet<T>::putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> database)
{
    const int data_sz = static_cast<int>(getDataStreamSize());
    FixedSizedStream stream(data_sz);
    packStream(stream);
    database->putInteger("data_sz", data_sz);
    database->putCharArray("data", static_cast<char*>(stream.getBufferStart()), data_sz);
    database->putIntegerArray("d_offset", d_offset, NDIM);
    return;
} // putToDatabase

template <class T>
inline void
LSet<T>::getFromDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> database)
{
    database->getIntegerArray("d_offset", d_offset, NDIM);
    const int data_sz = database->getInteger("data_sz");
    std::vector<char> data(data_sz);
    database->getCharArray("data", &data[0], data_sz);
    FixedSizedStream stream(&data[0], data_sz);
    unpackStream(stream, d_offset);
    return;
} // getFromDatabase

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LSet_inl_h
