// Filename: LData-inl.h
// Created on 17 Apr 2004 by Boyce Griffith
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

#ifndef included_LData_inl_h
#define included_LData_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "boost/multi_array.hpp"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/LData.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline const std::string&
LData::getName() const
{
    return d_name;
} // getName

inline unsigned int
LData::getGlobalNodeCount() const
{
    return d_global_node_count;
} // getGlobalNodeCount

inline unsigned int
LData::getLocalNodeCount() const
{
    return d_local_node_count;
} // getLocalNodeCount

inline unsigned int
LData::getGhostNodeCount() const
{
    return d_ghost_node_count;
} // getGhostNodeCount

inline unsigned int
LData::getDepth() const
{
    return d_depth;
} // getDepth

inline Vec
LData::getVec()
{
    restoreArrays();
    return d_global_vec;
} // getVec

inline boost::multi_array_ref<double, 1>*
LData::getArray()
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_depth == 1);
#endif
    if (!d_array) getArrayCommon();
    return d_boost_array;
} // getArray

inline boost::multi_array_ref<double, 1>*
LData::getLocalFormArray()
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_depth == 1);
#endif
    if (!d_array) getArrayCommon();
    return d_boost_local_array;
} // getLocalFormArray

inline boost::multi_array_ref<double, 1>*
LData::getGhostedLocalFormArray()
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_depth == 1);
#endif
    if (!d_ghosted_local_array) getGhostedLocalFormArrayCommon();
    return d_boost_ghosted_local_array;
} // getGhostedLocalFormArray

inline boost::multi_array_ref<double, 2>*
LData::getVecArray()
{
    if (!d_array) getArrayCommon();
    return d_boost_vec_array;
} // getVecArray

inline boost::multi_array_ref<double, 2>*
LData::getLocalFormVecArray()
{
    if (!d_array) getArrayCommon();
    return d_boost_local_vec_array;
} // getLocalFormVecArray

inline boost::multi_array_ref<double, 2>*
LData::getGhostedLocalFormVecArray()
{
    if (!d_ghosted_local_array) getGhostedLocalFormArrayCommon();
    return d_boost_vec_ghosted_local_array;
} // getGhostedLocalFormVecArray

inline void
LData::restoreArrays()
{
    int ierr;
    if (d_ghosted_local_array)
    {
        ierr = VecRestoreArray(d_ghosted_local_vec, &d_ghosted_local_array);
        IBTK_CHKERRQ(ierr);
        d_ghosted_local_array = NULL;
        delete d_boost_ghosted_local_array;
        delete d_boost_vec_ghosted_local_array;
        d_boost_ghosted_local_array = NULL;
        d_boost_vec_ghosted_local_array = NULL;
    }
    if (d_ghosted_local_vec)
    {
        ierr = VecGhostRestoreLocalForm(d_global_vec, &d_ghosted_local_vec);
        IBTK_CHKERRQ(ierr);
        d_ghosted_local_vec = NULL;
    }
    if (d_array)
    {
        ierr = VecRestoreArray(d_global_vec, &d_array);
        IBTK_CHKERRQ(ierr);
        d_array = NULL;
        delete d_boost_array;
        delete d_boost_local_array;
        delete d_boost_vec_array;
        delete d_boost_local_vec_array;
        d_boost_array = NULL;
        d_boost_local_array = NULL;
        d_boost_vec_array = NULL;
        d_boost_local_vec_array = NULL;
    }
    return;
} // restoreArray

inline void
LData::beginGhostUpdate()
{
    const int ierr = VecGhostUpdateBegin(getVec(), INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);
    return;
} // beginGhostUpdate

inline void
LData::endGhostUpdate()
{
    const int ierr = VecGhostUpdateEnd(getVec(), INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);
    return;
} // endGhostUpdate

/////////////////////////////// PRIVATE //////////////////////////////////////

inline void
LData::getArrayCommon()
{
    if (!d_array)
    {
        int ierr = VecGetArray(d_global_vec, &d_array);
        IBTK_CHKERRQ(ierr);
        int ilower, iupper;
        ierr = VecGetOwnershipRange(d_global_vec, &ilower, &iupper);
        IBTK_CHKERRQ(ierr);
        typedef boost::multi_array_types::extent_range range;
        if (d_depth == 1)
        {
            typedef boost::multi_array<double, 1> array_type;
            array_type::extent_gen extents;
            d_boost_array = new boost::multi_array_ref<double, 1>(d_array, extents[range(ilower, iupper)]);
            d_boost_local_array = new boost::multi_array_ref<double, 1>(d_array, extents[iupper - ilower]);
        }
        typedef boost::multi_array<double, 2> array_type;
        array_type::extent_gen extents;
        d_boost_vec_array =
            new boost::multi_array_ref<double, 2>(d_array, extents[range(ilower / d_depth, iupper / d_depth)][d_depth]);
        d_boost_local_vec_array =
            new boost::multi_array_ref<double, 2>(d_array, extents[(iupper - ilower) / d_depth][d_depth]);
    }
    return;
} // getArrayCommon

inline void
LData::getGhostedLocalFormArrayCommon()
{
    if (!d_ghosted_local_vec)
    {
        int ierr = VecGhostGetLocalForm(d_global_vec, &d_ghosted_local_vec);
        IBTK_CHKERRQ(ierr);
    }
    if (!d_ghosted_local_array)
    {
        int ierr = VecGetArray(d_ghosted_local_vec, &d_ghosted_local_array);
        IBTK_CHKERRQ(ierr);
        int ilower, iupper;
        ierr = VecGetOwnershipRange(d_ghosted_local_vec, &ilower, &iupper);
        IBTK_CHKERRQ(ierr);
        if (d_depth == 1)
        {
            d_boost_ghosted_local_array =
                new boost::multi_array_ref<double, 1>(d_ghosted_local_array, boost::extents[iupper - ilower]);
        }
        d_boost_vec_ghosted_local_array = new boost::multi_array_ref<double, 2>(
            d_ghosted_local_array, boost::extents[(iupper - ilower) / d_depth][d_depth]);
    }
    return;
} // getGhostedLocalFormArrayCommon

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LData_inl_h
