// Filename: LNodeLevelData.I
// Created on 17 Apr 2004 by Boyce Griffith
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

inline double&
LNodeLevelData::operator()(
    const int i,
    const int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(0<=depth && depth<d_depth);
#endif
    return getLocalFormArray()[d_depth*i+depth];
}// operator()

inline int
LNodeLevelData::getGlobalNodeCount()
{
    int ret_val;
    int ierr = VecGetSize(getGlobalVec(),&ret_val);  IBTK_CHKERRQ(ierr);
    return ret_val/d_depth;
}// getGlobalNodeCount

inline int
LNodeLevelData::getLocalNodeCount()
{
    int ret_val;
    int ierr = VecGetLocalSize(getGlobalVec(),&ret_val);  IBTK_CHKERRQ(ierr);
    return ret_val/d_depth;
}// getLocalNodeCount

inline int
LNodeLevelData::getLocalGhostNodeCount()
{
    return d_nonlocal_petsc_indices.size()/d_depth;
}// getLocalGhostNodeCount

inline int
LNodeLevelData::getDepth() const
{
    return d_depth;
}// getDepth

inline void
LNodeLevelData::beginGhostUpdate()
{
    int ierr;
    Vec v = getGlobalVec();
    ierr = VecGhostUpdateBegin(v, INSERT_VALUES, SCATTER_FORWARD);  IBTK_CHKERRQ(ierr);
    return;
}// beginGhostUpdate

inline void
LNodeLevelData::endGhostUpdate()
{
    int ierr;
    Vec v = getGlobalVec();
    ierr = VecGhostUpdateEnd(v, INSERT_VALUES, SCATTER_FORWARD);  IBTK_CHKERRQ(ierr);
    return;
}// endGhostUpdate

/////////////////////////////// PRIVATE //////////////////////////////////////

inline Vec&
LNodeLevelData::getGlobalVec()
{
    if (d_in_local_form)
    {
        restoreLocalFormVec();
    }
    return d_global_vec;
}// getGlobalVec

inline Vec&
LNodeLevelData::getLocalFormVec()
{
    if (!d_in_local_form)
    {
        const int ierr = VecGhostGetLocalForm(d_global_vec, &d_local_vec);  IBTK_CHKERRQ(ierr);
        d_in_local_form = true;
    }
    else if (d_extracted_local_array)
    {
        restoreLocalFormArray();
    }
    return d_local_vec;
}// getLocalFormVec

inline void
LNodeLevelData::restoreLocalFormVec()
{
    if (d_in_local_form)
    {
        if (d_extracted_local_array)
        {
            restoreLocalFormArray();
        }
        const int ierr = VecGhostRestoreLocalForm(d_global_vec, &d_local_vec);  IBTK_CHKERRQ(ierr);
        d_in_local_form = false;
    }
    return;
}// restoreLocalFormVec

inline PetscScalar*
LNodeLevelData::getLocalFormArray()
{
    if (!d_extracted_local_array)
    {
        if (!d_in_local_form)
        {
            getLocalFormVec();
        }
        const int ierr = VecGetArray(d_local_vec, &d_local_vec_array);  IBTK_CHKERRQ(ierr);
        d_extracted_local_array = true;
    }
    return d_local_vec_array;
}// getLocalFormArray

inline void
LNodeLevelData::restoreLocalFormArray()
{
    if (!d_in_local_form || !d_extracted_local_array) return;
    const int ierr = VecRestoreArray(d_local_vec, &d_local_vec_array);  IBTK_CHKERRQ(ierr);
    d_extracted_local_array = false;
    return;
}// restoreLocalFormArray

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
