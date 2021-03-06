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

#ifndef included_IBTK_LIndexSetData_inl_h
#define included_IBTK_LIndexSetData_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/LIndexSetData.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

template <class T>
inline const std::vector<int>&
LIndexSetData<T>::getLagrangianIndices() const
{
    return d_lag_indices;
} // getLagrangianIndices

template <class T>
inline const std::vector<int>&
LIndexSetData<T>::getInteriorLagrangianIndices() const
{
    return d_interior_lag_indices;
} // getInteriorLagrangianIndices

template <class T>
inline const std::vector<int>&
LIndexSetData<T>::getGhostLagrangianIndices() const
{
    return d_ghost_lag_indices;
} // getGhostLagrangianIndices

template <class T>
inline const std::vector<int>&
LIndexSetData<T>::getGlobalPETScIndices() const
{
    return d_global_petsc_indices;
} // getGlobalPETScIndices

template <class T>
inline const std::vector<int>&
LIndexSetData<T>::getInteriorGlobalPETScIndices() const
{
    return d_interior_global_petsc_indices;
} // getInteriorGlobalPETScIndices

template <class T>
inline const std::vector<int>&
LIndexSetData<T>::getGhostGlobalPETScIndices() const
{
    return d_ghost_global_petsc_indices;
} // getGhostGlobalPETScIndices

template <class T>
inline const std::vector<int>&
LIndexSetData<T>::getLocalPETScIndices() const
{
    return d_local_petsc_indices;
} // getLocalPETScIndices

template <class T>
inline const std::vector<int>&
LIndexSetData<T>::getInteriorLocalPETScIndices() const
{
    return d_interior_local_petsc_indices;
} // getInteriorLocalPETScIndices

template <class T>
inline const std::vector<int>&
LIndexSetData<T>::getGhostLocalPETScIndices() const
{
    return d_ghost_local_petsc_indices;
} // getGhostLocalPETScIndices

template <class T>
const std::vector<double>&
LIndexSetData<T>::getPeriodicShifts() const
{
    return d_periodic_shifts;
} // getPeriodicShifts

template <class T>
const std::vector<double>&
LIndexSetData<T>::getInteriorPeriodicShifts() const
{
    return d_interior_periodic_shifts;
} // getInteriorPeriodicShifts

template <class T>
const std::vector<double>&
LIndexSetData<T>::getGhostPeriodicShifts() const
{
    return d_ghost_periodic_shifts;
} // getGhostPeriodicShifts

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LIndexSetData_inl_h
