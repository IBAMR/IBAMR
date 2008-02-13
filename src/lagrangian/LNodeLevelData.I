// Filename: LNodeLevelData.I
// Last modified: <12.Feb.2008 21:21:11 griffith@box221.cims.nyu.edu>
// Created on 17 Apr 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// STOOLS INCLUDES
#include <stools/PETSC_SAMRAI_ERROR.h>

// SAMRAI INCLUDES
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
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
    int ierr = VecGetSize(getGlobalVec(),&ret_val);  PETSC_SAMRAI_ERROR(ierr);
    return ret_val/d_depth;
}// getGlobalNodeCount

inline int
LNodeLevelData::getLocalNodeCount()
{
    int ret_val;
    int ierr = VecGetLocalSize(getGlobalVec(),&ret_val);  PETSC_SAMRAI_ERROR(ierr);
    return ret_val/d_depth;
}// getLocalNodeCount

inline int
LNodeLevelData::getDepth() const
{
    return d_depth;
}// getDepth

inline void
LNodeLevelData::beginGhostUpdate()
{
    const int ierr = VecGhostUpdateBegin(
        getGlobalVec(), INSERT_VALUES, SCATTER_FORWARD);  PETSC_SAMRAI_ERROR(ierr);
    return;
}// beginGhostUpdate

inline void
LNodeLevelData::endGhostUpdate()
{
    const int ierr = VecGhostUpdateEnd(
        getGlobalVec(), INSERT_VALUES, SCATTER_FORWARD);  PETSC_SAMRAI_ERROR(ierr);
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
        const int ierr = VecGhostGetLocalForm(d_global_vec, &d_local_vec);  PETSC_SAMRAI_ERROR(ierr);
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
        const int ierr = VecGhostRestoreLocalForm(d_global_vec, &d_local_vec);  PETSC_SAMRAI_ERROR(ierr);
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
        const int ierr = VecGetArray(d_local_vec, &d_local_vec_array);  PETSC_SAMRAI_ERROR(ierr);
        d_extracted_local_array = true;
    }
    return d_local_vec_array;
}// getLocalFormArray

inline void
LNodeLevelData::restoreLocalFormArray()
{
    if (!d_in_local_form || !d_extracted_local_array) return;
    const int ierr = VecRestoreArray(d_local_vec, &d_local_vec_array);  PETSC_SAMRAI_ERROR(ierr);
    d_extracted_local_array = false;
    return;
}// restoreLocalFormArray

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
