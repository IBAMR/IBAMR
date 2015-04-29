// Filename: LSetDataIterator-inl.h
// Created on 11 Dec 2009 by Boyce Griffith
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

#ifndef included_LSetDataIterator_inl_h
#define included_LSetDataIterator_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/LSetDataIterator.h"
#include "ibtk/ibtk_utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

template <class T>
inline LSetDataIterator<T>::LSetDataIterator(const LSetData<T>& data)
    : d_box(DIM), d_index_it(data, /*begin*/ true), d_index_end(data, /*begin*/ false), d_node_set(NULL), d_node_it()
{
    // intentionally blank
    return;
}

template <class T>
inline LSetDataIterator<T>::LSetDataIterator(const LSetDataIterator<T>& that)
    : d_box(that.d_box), d_index_it(that.d_index_it), d_index_end(that.d_index_end), d_node_set(that.d_node_set),
      d_node_it(that.d_node_it)
{
    // intentionally blank
    return;
}

template <class T>
inline LSetDataIterator<T>::~LSetDataIterator()
{
    // intentionally blank
    return;
}

template <class T>
inline LSetDataIterator<T>& LSetDataIterator<T>::operator=(const LSetDataIterator<T>& that)
{
    if (this != &that)
    {
        d_box = that.d_box;
        d_index_it = that.d_index_it;
        d_index_end = that.d_index_end;
        d_node_set = that.d_node_set;
        d_node_it = that.d_node_it;
    }
    return *this;
}

template <class T>
inline bool LSetDataIterator<T>::operator==(const LSetDataIterator<T>& that)
{
    return ((!d_node_set && !that.d_node_set) || (d_box.isSpatiallyEqual(that.d_box) && d_index_it == that.d_index_it &&
                                                  d_node_set == that.d_node_set && d_node_it == that.d_node_it));
}

template <class T>
inline bool LSetDataIterator<T>::operator!=(const LSetDataIterator<T>& that)
{
    return !(*this == that);
}

template <class T>
inline LSetDataIterator<T>& LSetDataIterator<T>::operator++()
{
    if (!d_node_set) return *this;
    ++d_node_it;
    if (d_node_it != d_node_set->end()) return *this;
    d_node_set = NULL;
    d_index_it++;
    while ((d_index_it != d_index_end) && !d_box.contains(d_index_it.getIndex()))
    {
        ++d_index_it;
    }
    if (d_index_it != d_index_end)
    {
        d_node_set = &(*d_index_it);
        TBOX_ASSERT(!d_node_set->empty());
        d_node_it = d_node_set->begin();
    }
    return *this;
}

template <class T>
inline LSetDataIterator<T> LSetDataIterator<T>::operator++(int)
{
    LSetDataIterator<T> tmp(*this);
    ++(*this);
    return tmp;
}

template <class T>
inline typename LSet<T>::value_type& LSetDataIterator<T>::operator*() const
{
    return getDataItem();
}

template <class T>
inline typename LSet<T>::value_type& LSetDataIterator<T>::getDataItem() const
{
    TBOX_ASSERT((d_index_it != d_index_end) && d_node_set);
    return *d_node_it;
}

template <class T>
inline const SAMRAI::hier::Index& LSetDataIterator<T>::getCellIndex() const
{
    return d_index_it.getIndex();
}

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LSetDataIterator_inl_h
