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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/LIndexSetData.h"
#include "ibtk/LNode.h"
#include "ibtk/LNodeIndex.h"
#include "ibtk/LSet.h"
#include "ibtk/LSetData.h"

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellIndex.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "tbox/Pointer.h"

#include <algorithm>
#include <array>
#include <utility>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

template <class T>
LIndexSetData<T>::LIndexSetData(Box<NDIM> box, IntVector<NDIM> ghosts) : LSetData<T>(std::move(box), std::move(ghosts))
{
    // intentionally blank
    return;
} // LIndexSetData

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
    const hier::Index<NDIM>& ilower = patch_box.lower();
    const hier::Index<NDIM>& iupper = patch_box.upper();

    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    std::array<bool, NDIM> patch_touches_lower_periodic_bdry, patch_touches_upper_periodic_bdry;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        patch_touches_lower_periodic_bdry[axis] = pgeom->getTouchesPeriodicBoundary(axis, 0);
        patch_touches_upper_periodic_bdry[axis] = pgeom->getTouchesPeriodicBoundary(axis, 1);
    }

    for (typename LSetData<T>::SetIterator it(*this); it; it++)
    {
        const CellIndex<NDIM>& i = it.getIndex();
        std::array<int, NDIM> offset;
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
        for (auto n = idx_set.begin(); n != idx_set.end(); ++n)
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
