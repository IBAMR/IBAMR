//
// LNodeLevelData.C
//
// Created on 17 Apr 2004
//         by Boyce E. Griffith (boyce@trasnaform.speakeasy.net).
//
// Last modified: <27.Jun.2005 01:49:48 boyce@mstu1.cims.nyu.edu>
//

#include "LNodeLevelData.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

// SAMRAI-tools INCLUDES
//
#include "LDataManager.h"
#include "PETSC_SAMRAI_ERROR.h"

// SAMRAI INCLUDES
//
#include "tbox/Utilities.h"

/////////////////////////////// INLINE ///////////////////////////////////////

#ifdef DEBUG_NO_INLINE
#include "LNodeLevelData.I"
#endif

/////////////////////////////// PUBLIC ///////////////////////////////////////

LNodeLevelData::~LNodeLevelData()
{
    if (d_in_local_form) restoreLocalFormVec();
    const int ierr = VecDestroy(d_global_vec);  PETSC_SAMRAI_ERROR(ierr);
    return;
}// ~LNodeLevelData

LNodeLevelData& LNodeLevelData::operator=(
    const LNodeLevelData& that)
{
    if (that.d_in_local_form)
    {
        TBOX_ERROR("LNodeLevelData operator=:\n" <<
                   "  operator=() source object must NOT be in local form." << endl);
    }
    if (this == &that) return(*this);  // check for self-assignment

    restoreLocalFormVec();
    const int ierr = VecDestroy(d_global_vec);  PETSC_SAMRAI_ERROR(ierr);
    
    assignThatToThis(that);
    return(*this);
}// operator=

void LNodeLevelData::putToDatabase(
    tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!db.isNull());
#endif
    const bool restore_local_vec = !d_in_local_form;
    const bool restore_local_array = !d_extracted_local_array;

    const int num_local_nodes = getLocalNodeCount();
    const double* const local_form_array = getLocalFormArray();
    const int n_nonlocal_indices = static_cast<int>(d_nonlocal_petsc_indices.size());
    
    db->putString("d_name", d_name);
    db->putInteger("d_depth", d_depth);
    db->putInteger("num_local_nodes", num_local_nodes);
    db->putInteger("n_nonlocal_indices", n_nonlocal_indices);
    if (n_nonlocal_indices > 0)
    {
        db->putIntegerArray("d_nonlocal_petsc_indices",
                            &d_nonlocal_petsc_indices[0],
                            n_nonlocal_indices);
    }
    if (d_depth*num_local_nodes + n_nonlocal_indices > 0)
    {
        db->putDoubleArray("vals", local_form_array,
                           d_depth*num_local_nodes + n_nonlocal_indices);
    }
    
    if (restore_local_array) restoreLocalFormArray();
    if (restore_local_vec  ) restoreLocalFormVec();
    
    return;
}// putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

LNodeLevelData::LNodeLevelData(
    const string& name,
    const int num_local_nodes,
    const int depth,
    const vector<int>& nonlocal_petsc_indices)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!name.empty());
    assert(num_local_nodes >= 0);
    assert(depth > 0);
#endif
    d_name = name;
    d_depth = depth;
    d_nonlocal_petsc_indices = nonlocal_petsc_indices;
    
    // Create the PETSc Vec which actually provides the storage for
    // the Lagrangian data.
    int ierr;

    if (d_depth == 1)
    {
        ierr = VecCreateGhost(PETSC_COMM_WORLD,
                              num_local_nodes, PETSC_DECIDE,
                              d_nonlocal_petsc_indices.size(),
                              &d_nonlocal_petsc_indices[0],
                              &d_global_vec);
        PETSC_SAMRAI_ERROR(ierr);
    }
    else
    {
        ierr = VecCreateGhostBlock(PETSC_COMM_WORLD, d_depth,
                                   d_depth*num_local_nodes, PETSC_DECIDE,
                                   d_nonlocal_petsc_indices.size(),
                                   &d_nonlocal_petsc_indices[0],
                                   &d_global_vec);
        PETSC_SAMRAI_ERROR(ierr);
    }

    ierr = VecSetBlockSize(d_global_vec, d_depth);  PETSC_SAMRAI_ERROR(ierr);
    
    d_in_local_form = false;

    d_local_vec_array = NULL;
    d_extracted_local_array = false;
    
    return;
}// LNodeLevelData

LNodeLevelData::LNodeLevelData(
    tbox::Pointer<tbox::Database> db)
{
    d_name = db->getString("d_name");
    d_depth = db->getInteger("d_depth");
    int num_local_nodes = db->getInteger("num_local_nodes");
    int n_nonlocal_indices = db->getInteger("n_nonlocal_indices");
    d_nonlocal_petsc_indices.resize(n_nonlocal_indices);
    if (n_nonlocal_indices > 0)
    {
        db->getIntegerArray("d_nonlocal_petsc_indices",
                            &d_nonlocal_petsc_indices[0], n_nonlocal_indices);
    }
    
    // Create the PETSc Vec which actually provides the storage for
    // the Lagrangian data.
    int ierr;
    
    if (d_depth == 1)
    {
        ierr = VecCreateGhost(PETSC_COMM_WORLD,
                              num_local_nodes, PETSC_DECIDE,
                              d_nonlocal_petsc_indices.size(),
                              &d_nonlocal_petsc_indices[0],
                              &d_global_vec);
        PETSC_SAMRAI_ERROR(ierr);
    }
    else
    {
        ierr = VecCreateGhostBlock(PETSC_COMM_WORLD, d_depth,
                                   d_depth*num_local_nodes, PETSC_DECIDE,
                                   d_nonlocal_petsc_indices.size(),
                                   &d_nonlocal_petsc_indices[0],
                                   &d_global_vec);
        PETSC_SAMRAI_ERROR(ierr);
    }
    
    ierr = VecSetBlockSize(d_global_vec, d_depth);  PETSC_SAMRAI_ERROR(ierr);
    
    d_in_local_form = false;

    d_local_vec_array = NULL;
    d_extracted_local_array = false;

    // Extract the values from the database.
    double* local_form_array = getLocalFormArray();
    if (d_depth*num_local_nodes + n_nonlocal_indices > 0)
    {
        db->getDoubleArray("vals", local_form_array,
                           d_depth*num_local_nodes + n_nonlocal_indices);
    }
    restoreLocalFormVec();
    return;
}// LNodeLevelData

LNodeLevelData::LNodeLevelData(
    const LNodeLevelData& from)
{
    if (from.d_in_local_form)
    {
        TBOX_ERROR("LNodeLevelData copy constructor:\n" <<
                   "  copy constructor source object must NOT be in local form." << endl);
    }
    
    d_in_local_form = false;
    d_local_vec_array = NULL;
    d_extracted_local_array = false;
    
    assignThatToThis(from);    
    return;
}// LNodeLevelData

void LNodeLevelData::resetData(
    Vec& new_global_vec,
    const vector<int>& new_nonlocal_petsc_indices)
{
    if (d_in_local_form)
    {
        TBOX_ERROR("LNodeLevelData::resetData()\n" <<
                   "  object must NOT be in local form." << endl);
    }
    d_global_vec = new_global_vec;
    d_nonlocal_petsc_indices = new_nonlocal_petsc_indices;
    return;
}// resetData
    
/////////////////////////////// PRIVATE //////////////////////////////////////

void LNodeLevelData::assignThatToThis(
    const LNodeLevelData& that)
{
    d_name = that.d_name;
    d_depth = that.d_depth;
    
    int ierr;

    ierr = VecDuplicate(that.d_global_vec, &d_global_vec);
    PETSC_SAMRAI_ERROR(ierr);
    
    ierr = VecCopy(that.d_global_vec, d_global_vec);
    PETSC_SAMRAI_ERROR(ierr);
    
    d_in_local_form = false;
    
    d_local_vec_array = NULL;
    d_extracted_local_array = false;
    return;
}// assignThatToThis

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#ifndef LACKS_EXPLICIT_TEMPLATE_INSTANTIATION

#include "tbox/Pointer.C"

//////////////////////////////////////////////////////////////////////
///
/// These declarations are required to use the LNodeLevelData class.
///
//////////////////////////////////////////////////////////////////////

template class tbox::Pointer<LNodeLevelData>;

#endif

//////////////////////////////////////////////////////////////////////////////
