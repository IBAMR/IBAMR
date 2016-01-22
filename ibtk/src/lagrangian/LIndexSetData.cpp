// Filename: LIndexSetData.cpp
// Created on 13 May 2011 by Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <vector>

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellIndex.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "boost/array.hpp"
#include "ibtk/LIndexSetData.h"
#include "ibtk/LSet.h"
#include "ibtk/LSetData.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Pointer.h"

namespace IBTK
{
class LNode;
class LNodeIndex;
} // namespace IBTK

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

template <class T>
LIndexSetData<T>::LIndexSetData(const Box<NDIM>& box, const IntVector<NDIM>& ghosts)
    : LSetData<T>(box, ghosts),
      d_lag_indices(),
      d_interior_lag_indices(),
      d_ghost_lag_indices(),
      d_global_petsc_indices(),
      d_interior_global_petsc_indices(),
      d_ghost_global_petsc_indices(),
      d_local_petsc_indices(),
      d_interior_local_petsc_indices(),
      d_ghost_local_petsc_indices(),
      d_periodic_shifts(),
      d_interior_periodic_shifts(),
      d_ghost_periodic_shifts()
{
    // intentionally blank
    return;
} // LIndexSetData

template <class T>
LIndexSetData<T>::~LIndexSetData()
{
    // intentionally blank
    return;
} // ~LIndexSetData

template <class T>
void
LIndexSetData<T>::cacheLocalIndices(Pointer<Patch<NDIM> > patch, const IntVector<NDIM>& periodic_shift)
{
    d_lag_indices.clear();
    d_interior_lag_indices.clear();
    d_ghost_lag_indices.clear();
    d_global_petsc_indices.clear();
    d_interior_global_petsc_indices.clear();
    d_ghost_global_petsc_indices.clear();
    d_local_petsc_indices.clear();
    d_interior_local_petsc_indices.clear();
    d_ghost_local_petsc_indices.clear();
    d_periodic_shifts.clear();
    d_interior_periodic_shifts.clear();
    d_ghost_periodic_shifts.clear();

    const Box<NDIM>& patch_box = patch->getBox();
    const Index<NDIM>& ilower = patch_box.lower();
    const Index<NDIM>& iupper = patch_box.upper();

    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    boost::array<bool, NDIM> patch_touches_lower_periodic_bdry, patch_touches_upper_periodic_bdry;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        patch_touches_lower_periodic_bdry[axis] = pgeom->getTouchesPeriodicBoundary(axis, 0);
        patch_touches_upper_periodic_bdry[axis] = pgeom->getTouchesPeriodicBoundary(axis, 1);
    }

    for (typename LSetData<T>::SetIterator it(*this); it; it++)
    {
        const CellIndex<NDIM>& i = it.getIndex();
        boost::array<int, NDIM> offset;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (patch_touches_lower_periodic_bdry[d] && i(d) < ilower(d))
            {
                offset[d] = -periodic_shift(d); // X is ABOVE the top    of the patch --- need
                                                // to shift DOWN
            }
            else if (patch_touches_upper_periodic_bdry[d] && i(d) > iupper(d))
            {
                offset[d] = +periodic_shift(d); // X is BELOW the bottom of the patch --- need to shift UP
            }
            else
            {
                offset[d] = 0;
            }
        }
        const LSet<T>& idx_set = *it;
        const bool patch_owns_idx_set = patch_box.contains(i);
        for (typename LSet<T>::const_iterator n = idx_set.begin(); n != idx_set.end(); ++n)
        {
            const typename LSet<T>::value_type& idx = *n;
            const int lag_idx = idx->getLagrangianIndex();
            const int global_petsc_idx = idx->getGlobalPETScIndex();
            const int local_petsc_idx = idx->getLocalPETScIndex();
            d_lag_indices.push_back(lag_idx);
            d_global_petsc_indices.push_back(global_petsc_idx);
            d_local_petsc_indices.push_back(local_petsc_idx);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                d_periodic_shifts.push_back(static_cast<double>(offset[d]) * dx[d]);
            }
            if (patch_owns_idx_set)
            {
                d_interior_lag_indices.push_back(lag_idx);
                d_interior_global_petsc_indices.push_back(global_petsc_idx);
                d_interior_local_petsc_indices.push_back(local_petsc_idx);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    d_interior_periodic_shifts.push_back(static_cast<double>(offset[d]) * dx[d]);
                }
            }
            else
            {
                d_ghost_lag_indices.push_back(lag_idx);
                d_ghost_global_petsc_indices.push_back(global_petsc_idx);
                d_ghost_local_petsc_indices.push_back(local_petsc_idx);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    d_ghost_periodic_shifts.push_back(static_cast<double>(offset[d]) * dx[d]);
                }
            }
        }
    }
    return;
} // cacheLocalIndices

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

template class IBTK::LIndexSetData<IBTK::LNode>;
template class IBTK::LIndexSetData<IBTK::LNodeIndex>;

//////////////////////////////////////////////////////////////////////////////
