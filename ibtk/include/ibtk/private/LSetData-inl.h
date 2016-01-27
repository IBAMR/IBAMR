// Filename: LSetData-inl.h
// Created on 04 Jun 2007 by Boyce Griffith
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

#ifndef included_LSetData_inl_h
#define included_LSetData_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

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
            it.d_node_set = NULL;
        }
    }
    else
    {
        it.d_node_set = NULL;
    }
    return it;
} // data_begin

template <class T>
inline typename LSetData<T>::DataIterator
LSetData<T>::data_end()
{
    typename LSetData<T>::DataIterator it;
    it.d_node_set = NULL;
    return it;
} // data_end

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LSetData_inl_h
