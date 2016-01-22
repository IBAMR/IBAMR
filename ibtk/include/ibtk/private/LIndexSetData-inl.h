// Filename: LIndexSetData-inl.h
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

#ifndef included_LIndexSetData_inl_h
#define included_LIndexSetData_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

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

#endif //#ifndef included_LIndexSetData_inl_h
