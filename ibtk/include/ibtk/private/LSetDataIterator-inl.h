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

#ifndef included_IBTK_LSetDataIterator_inl_h
#define included_IBTK_LSetDataIterator_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/LSetDataIterator.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

template <class T>
inline LSetDataIterator<T>::LSetDataIterator() : d_box(), d_index_it(), d_node_set(), d_node_it()
{
    // intentionally blank
    return;
} // LSetDataIterator

template <class T>
inline LSetDataIterator<T>::LSetDataIterator(const LSetDataIterator<T>& that)
    : d_box(that.d_box), d_index_it(that.d_index_it), d_node_set(that.d_node_set), d_node_it(that.d_node_it)
{
    // intentionally blank
    return;
} // LSetDataIterator

template <class T>
inline LSetDataIterator<T>::~LSetDataIterator()
{
    // intentionally blank
    return;
} // ~LSetDataIterator

template <class T>
inline LSetDataIterator<T>&
LSetDataIterator<T>::operator=(const LSetDataIterator<T>& that)
{
    if (this != &that)
    {
        d_box = that.d_box;
        d_index_it = that.d_index_it;
        d_node_set = that.d_node_set;
        d_node_it = that.d_node_it;
    }
    return *this;
} // operator=

template <class T>
inline bool
LSetDataIterator<T>::operator==(const LSetDataIterator<T>& that)
{
    return ((!d_node_set && !that.d_node_set) || (d_box == that.d_box && d_index_it == that.d_index_it &&
                                                  d_node_set == that.d_node_set && d_node_it == that.d_node_it));
} // operator==

template <class T>
inline bool
LSetDataIterator<T>::operator!=(const LSetDataIterator<T>& that)
{
    return !(*this == that);
} // operator!=

template <class T>
inline LSetDataIterator<T>&
LSetDataIterator<T>::operator++()
{
    if (!d_node_set) return *this;
    ++d_node_it;
    if (d_node_it != d_node_set->end()) return *this;
    d_node_set = nullptr;
    d_index_it++;
    while (d_index_it && !d_box.contains(d_index_it.getIndex()))
    {
        d_index_it++;
    }
    if (d_index_it)
    {
        d_node_set = &(*d_index_it);
#if !defined(NDEBUG)
        TBOX_ASSERT(!d_node_set->empty());
#endif
        d_node_it = d_node_set->begin();
    }
    return *this;
} // operator++

template <class T>
inline LSetDataIterator<T>
LSetDataIterator<T>::operator++(int)
{
    LSetDataIterator<T> tmp(*this);
    ++(*this);
    return tmp;
} // operator++

template <class T>
inline typename LSet<T>::value_type&
LSetDataIterator<T>::operator*() const
{
    return getDataItem();
} // operator*

template <class T>
inline typename LSet<T>::value_type&
LSetDataIterator<T>::getDataItem() const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_index_it && d_node_set);
#endif
    return *d_node_it;
} // getDataItem

template <class T>
inline const SAMRAI::hier::Index<NDIM>&
LSetDataIterator<T>::getCellIndex() const
{
    return d_index_it.getIndex();
} // getCellIndex

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LSetDataIterator_inl_h
