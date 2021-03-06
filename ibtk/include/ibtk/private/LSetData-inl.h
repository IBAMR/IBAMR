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

#ifndef included_IBTK_LSetData_inl_h
#define included_IBTK_LSetData_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/LSetData.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

template <class T>
inline typename LSetData<T>::DataIterator
LSetData<T>::data_begin(const SAMRAI::hier::Box<NDIM>& box)
{
    typename LSetData<T>::DataIterator it;
    it.d_box = box * this->getGhostBox();
    it.d_index_it = SAMRAI::pdat::IndexIterator<NDIM, LSet<T>, SAMRAI::pdat::CellGeometry<NDIM> >(*this);
    if (it.d_index_it)
    {
        it.d_node_set = &(*it.d_index_it);
        while (it.d_index_it && !it.d_box.contains(it.d_index_it.getIndex()))
        {
            it.d_index_it++;
        }

        if (it.d_index_it && it.d_box.contains(it.d_index_it.getIndex()))
        {
            it.d_node_set = &(*it.d_index_it);
#if !defined(NDEBUG)
            TBOX_ASSERT(!it.d_node_set->empty());
#endif
            it.d_node_it = it.d_node_set->begin();
        }
        else
        {
            it.d_node_set = nullptr;
        }
    }
    else
    {
        it.d_node_set = nullptr;
    }
    return it;
} // data_begin

template <class T>
inline typename LSetData<T>::DataIterator
LSetData<T>::data_end()
{
    typename LSetData<T>::DataIterator it;
    it.d_node_set = nullptr;
    return it;
} // data_end

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LSetData_inl_h
