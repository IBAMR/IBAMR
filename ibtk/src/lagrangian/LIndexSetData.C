// Filename: LIndexSetData.C
// Created on 13 May 2011 by Boyce Griffith
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

#include "LIndexSetData.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

template<class T>
LIndexSetData<T>::LIndexSetData(
    const Box<NDIM>& box,
    const IntVector<NDIM>& ghosts)
    : LSetData<T>(box,ghosts),
      d_interior_lag_indices(),
      d_ghost_lag_indices(),
      d_interior_global_petsc_indices(),
      d_ghost_global_petsc_indices(),
      d_interior_local_petsc_indices(),
      d_ghost_local_petsc_indices()
{
    // intentionally blank
    return;
}// LIndexSetData

template<class T>
LIndexSetData<T>::~LIndexSetData()
{
    // intentionally blank
    return;
}// ~LIndexSetData

template<class T>
void
LIndexSetData<T>::cacheLocalIndices()
{
    d_interior_lag_indices         .clear();
    d_ghost_lag_indices            .clear();
    d_interior_global_petsc_indices.clear();
    d_ghost_global_petsc_indices   .clear();
    d_interior_local_petsc_indices .clear();
    d_ghost_local_petsc_indices    .clear();
    const Box<NDIM>& patch_box = LSetData<T>::getBox();
    for (typename LSetData<T>::SetIterator it(*this); it; it++)
    {
        const CellIndex<NDIM>& i = it.getIndex();
        const LSet<T>& idx_set = *it;
        const bool patch_owns_idx_set = patch_box.contains(i);
        for (typename LSet<T>::const_iterator n = idx_set.begin(); n != idx_set.end(); ++n)
        {
            const typename LSet<T>::value_type& idx = *n;
            const int          lag_idx = idx.getLagrangianIndex();
            const int global_petsc_idx = idx.getGlobalPETScIndex();
            const int  local_petsc_idx = idx.getLocalPETScIndex();
            if (patch_owns_idx_set)
            {
                d_interior_lag_indices         .push_back(         lag_idx);
                d_interior_global_petsc_indices.push_back(global_petsc_idx);
                d_interior_local_petsc_indices .push_back( local_petsc_idx);
            }
            else
            {
                d_ghost_lag_indices         .push_back(         lag_idx);
                d_ghost_global_petsc_indices.push_back(global_petsc_idx);
                d_ghost_local_petsc_indices .push_back( local_petsc_idx);
            }
        }
    }
    return;
}// cacheLocalIndices

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <ibtk/LNodeSet.h>
template class IBTK::LIndexSetData<IBTK::LNode>;
template class Pointer<IBTK::LIndexSetData<IBTK::LNode> >;

#include <ibtk/LNodeIndexSet.h>
template class IBTK::LIndexSetData<IBTK::LNodeIndex>;
template class Pointer<IBTK::LIndexSetData<IBTK::LNodeIndex> >;

//////////////////////////////////////////////////////////////////////////////
