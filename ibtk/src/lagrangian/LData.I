// Filename: LData.I
// Created on 17 Apr 2004 by Boyce Griffith
//
// Copyright (c) 2002-2013, Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>

// SAMRAI INCLUDES
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline const std::string&
LData::getName() const
{
    return d_name;
}// getName

inline unsigned int
LData::getGlobalNodeCount() const
{
    return d_global_node_count;
}// getGlobalNodeCount

inline unsigned int
LData::getLocalNodeCount() const
{
    return d_local_node_count;
}// getLocalNodeCount

inline unsigned int
LData::getGhostNodeCount() const
{
    return d_ghost_node_count;
}// getGhostNodeCount

inline unsigned int
LData::getDepth() const
{
    return d_depth;
}// getDepth

inline Vec
LData::getVec()
{
    restoreArrays();
    return d_global_vec;
}// getVec

inline blitz::Array<double,1>*
LData::getArray()
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_depth == 1);
#endif
    if (!d_array) getArrayCommon();
    return &d_blitz_array;
}// getArray

inline blitz::Array<double,1>*
LData::getLocalFormArray()
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_depth == 1);
#endif
    if (!d_array) getArrayCommon();
    return &d_blitz_local_array;
}// getLocalFormArray

inline blitz::Array<double,1>*
LData::getGhostedLocalFormArray()
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_depth == 1);
#endif
    if (!d_ghosted_local_array) getGhostedLocalFormArrayCommon();
    return &d_blitz_ghosted_local_array;
}// getGhostedLocalFormArray

inline blitz::Array<double,2>*
LData::getVecArray()
{
    if (!d_array) getArrayCommon();
    return &d_blitz_vec_array;
}// getVecArray

inline blitz::Array<double,2>*
LData::getLocalFormVecArray()
{
    if (!d_array) getArrayCommon();
    return &d_blitz_local_vec_array;
}// getLocalFormVecArray

inline blitz::Array<double,2>*
LData::getGhostedLocalFormVecArray()
{
    if (!d_ghosted_local_array) getGhostedLocalFormArrayCommon();
    return &d_blitz_vec_ghosted_local_array;
}// getGhostedLocalFormVecArray

inline void
LData::restoreArrays()
{
    if (d_array)
    {
        int ierr = VecRestoreArray(d_global_vec, &d_array);  IBTK_CHKERRQ(ierr);
        d_array = NULL;
        d_blitz_array.free();
        d_blitz_local_array.free();
        d_blitz_vec_array.free();
        d_blitz_local_vec_array.free();
    }
    if (d_ghosted_local_array)
    {
        int ierr = VecRestoreArray(d_ghosted_local_vec, &d_ghosted_local_array);  IBTK_CHKERRQ(ierr);
        d_ghosted_local_array = NULL;
        d_blitz_ghosted_local_array.free();
        d_blitz_vec_ghosted_local_array.free();
    }
    if (d_ghosted_local_vec)
    {
        int ierr = VecGhostRestoreLocalForm(d_global_vec, &d_ghosted_local_vec);  IBTK_CHKERRQ(ierr);
        d_ghosted_local_vec = NULL;
    }
    return;
}// restoreArray

inline void
LData::beginGhostUpdate()
{
    const int ierr = VecGhostUpdateBegin(getVec(),INSERT_VALUES,SCATTER_FORWARD);  IBTK_CHKERRQ(ierr);
    return;
}// beginGhostUpdate

inline void
LData::endGhostUpdate()
{
    const int ierr = VecGhostUpdateEnd(getVec(),INSERT_VALUES,SCATTER_FORWARD);  IBTK_CHKERRQ(ierr);
    return;
}// endGhostUpdate

/////////////////////////////// PRIVATE //////////////////////////////////////

inline void
LData::getArrayCommon()
{
    if (!d_array)
    {
        int ierr = VecGetArray(d_global_vec, &d_array);  IBTK_CHKERRQ(ierr);
        int ilower, iupper;
        ierr = VecGetOwnershipRange(d_global_vec, &ilower, &iupper);  IBTK_CHKERRQ(ierr);
        if (iupper-ilower == 0)
        {
            d_blitz_array.free();
            d_blitz_local_array.free();
            d_blitz_vec_array.free();
            d_blitz_local_vec_array.free();
        }
        else
        {
            blitz::TinyVector<int,1> shape(iupper-ilower);
            blitz::TinyVector<int,2> vec_shape((iupper-ilower)/d_depth,d_depth);
            blitz::GeneralArrayStorage<1> global_storage;
            global_storage.base() = ilower;
            blitz::GeneralArrayStorage<2> vec_global_storage;
            vec_global_storage.base() = ilower/d_depth , 0;
            if (d_depth == 1)
            {
                d_blitz_array.reference(blitz::Array<double,1>(d_array, shape, blitz::neverDeleteData, global_storage));
                d_blitz_local_array.reference(blitz::Array<double,1>(d_array, shape, blitz::neverDeleteData));
            }
            d_blitz_vec_array.reference(blitz::Array<double,2>(d_array, vec_shape, blitz::neverDeleteData, vec_global_storage));
            d_blitz_local_vec_array.reference(blitz::Array<double,2>(d_array, vec_shape, blitz::neverDeleteData));
        }
    }
    return;
}// getArrayCommon

inline void
LData::getGhostedLocalFormArrayCommon()
{
    if (!d_ghosted_local_vec)
    {
        int ierr = VecGhostGetLocalForm(d_global_vec, &d_ghosted_local_vec);  IBTK_CHKERRQ(ierr);
    }
    if (!d_ghosted_local_array)
    {
        int ierr = VecGetArray(d_ghosted_local_vec, &d_ghosted_local_array);  IBTK_CHKERRQ(ierr);
        int ilower, iupper;
        ierr = VecGetOwnershipRange(d_ghosted_local_vec, &ilower, &iupper);  IBTK_CHKERRQ(ierr);
        if (iupper-ilower == 0)
        {
            d_blitz_ghosted_local_array.free();
            d_blitz_vec_ghosted_local_array.free();
        }
        else
        {
            blitz::TinyVector<int,1> shape(iupper-ilower);
            blitz::TinyVector<int,2> vec_shape((iupper-ilower)/d_depth,d_depth);
            if (d_depth == 1)
            {
                d_blitz_ghosted_local_array.reference(blitz::Array<double,1>(d_ghosted_local_array, shape, blitz::neverDeleteData));
            }
            d_blitz_vec_ghosted_local_array.reference(blitz::Array<double,2>(d_ghosted_local_array, vec_shape, blitz::neverDeleteData));
        }
    }
    return;
}// getGhostedLocalFormArrayCommon

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
