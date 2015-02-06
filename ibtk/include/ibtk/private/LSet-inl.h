// Filename: LSet-inl.h
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

#ifndef included_LSet_inl_h
#define included_LSet_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/LSet.h"
#include "ibtk/ibtk_utilities.h"
#include "SAMRAI/tbox/Database.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

template <class T>
inline LSet<T>::LSet()
    : d_set(), d_offset(DIM, 0)
{
    // intentionally blank
    return;
} // LSet

template <class T>
inline LSet<T>::LSet(const LSet<T>& from)
    : d_set(from.d_set), d_offset(from.d_offset)
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
inline LSet<T>& LSet<T>::operator=(const LSet<T>& that)
{
    if (this == &that) return *this; // check for self-assignment

    d_set = that.d_set;
    d_offset = that.d_offset;
    return *this;
} // operator=

template <class T>
inline typename LSet<T>::reference LSet<T>::operator[](typename LSet<T>::size_type n)
{
    return d_set[n];
} // operator[]

template <class T>
inline typename LSet<T>::const_reference LSet<T>::operator[](typename LSet<T>::size_type n) const
{
    return d_set[n];
} // operator[]

template <class T>
inline typename LSet<T>::const_iterator LSet<T>::begin() const
{
    return d_set.begin();
} // begin

template <class T>
inline typename LSet<T>::iterator LSet<T>::begin()
{
    return d_set.begin();
} // begin

template <class T>
inline typename LSet<T>::const_iterator LSet<T>::end() const
{
    return d_set.end();
} // end

template <class T>
inline typename LSet<T>::iterator LSet<T>::end()
{
    return d_set.end();
} // end

template <class T>
inline typename LSet<T>::size_type LSet<T>::size() const
{
    return d_set.size();
} // size

template <class T>
inline bool LSet<T>::empty() const
{
    return d_set.empty();
} // empty

template <class T>
inline void LSet<T>::push_back(const LSet<T>::value_type& value)
{
    d_set.push_back(value);
    return;
} // push_back

template <class T>
inline typename LSet<T>::iterator LSet<T>::insert(typename LSet<T>::iterator pos, const typename LSet<T>::value_type& x)
{
    return d_set.insert(pos, x);
} // insert

template <class T>
template <class InputIterator>
inline void LSet<T>::insert(typename LSet<T>::iterator pos, InputIterator first, InputIterator last)
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
inline const typename LSet<T>::DataSet& LSet<T>::getDataSet() const
{
    return d_set;
} // getDataSet

template <class T>
inline typename LSet<T>::DataSet& LSet<T>::getDataSet()
{
    return d_set;
} // getDataSet

template <class T>
inline void LSet<T>::setDataSet(const typename LSet<T>::DataSet& set)
{
    d_set = set;
    return;
} // setDataSet

template <class T>
inline const SAMRAI::hier::IntVector& LSet<T>::getPeriodicOffset() const
{
    return d_offset;
} // getPeriodicOffset

template <class T>
inline void LSet<T>::setPeriodicOffset(const SAMRAI::hier::IntVector& offset)
{
    d_offset = offset;
    return;
} // setPeriodicOffset

template <class T>
inline void LSet<T>::copySourceItem(const SAMRAI::hier::Index& /*src_index*/,
                                    const SAMRAI::hier::IntVector& src_offset,
                                    const LSet<T>& src_item)
{
    d_set = src_item.d_set;
    d_offset = src_offset;
    return;
} // copySourceItem

template <class T>
inline size_t LSet<T>::getDataStreamSize() const
{
    size_t size = SAMRAI::tbox::MessageStream::getSizeof<int>();
    for (unsigned int k = 0; k < d_set.size(); ++k)
    {
        size += d_set[k]->getDataStreamSize();
    }
    return size;
} // getDataStreamSize

template <class T>
inline void LSet<T>::packStream(SAMRAI::tbox::MessageStream& stream)
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
inline void LSet<T>::unpackStream(SAMRAI::tbox::MessageStream& stream, const SAMRAI::hier::IntVector& offset)
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
inline void LSet<T>::putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> database)
{
    const int data_sz = static_cast<int>(getDataStreamSize());
    SAMRAI::tbox::MessageStream stream(data_sz, SAMRAI::tbox::MessageStream::Write);
    packStream(stream);
    database->putInteger("data_sz", data_sz);
    database->putCharArray("data", static_cast<char*>(stream.getBufferStart()), data_sz);
    database->putIntegerArray("d_offset", &d_offset[0], NDIM);
    return;
} // putToDatabase

template <class T>
inline void LSet<T>::getFromDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> database)
{
    database->getIntegerArray("d_offset", &d_offset[0], NDIM);
    const int data_sz = database->getInteger("data_sz");
    SAMRAI::tbox::MessageStream stream(data_sz, SAMRAI::tbox::MessageStream::Read);
    database->getCharArray("data", static_cast<char*>(stream.getBufferStart()), data_sz);
    unpackStream(stream, d_offset);
    return;
} // getFromDatabase

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LSet_inl_h
