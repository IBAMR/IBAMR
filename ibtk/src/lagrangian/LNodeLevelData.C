// Filename: LNodeLevelData.C
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

#include "LNodeLevelData.h"

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
#include <ibtk/LDataManager.h>
#include <ibtk/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LNodeLevelData::~LNodeLevelData()
{
    if (d_in_local_form) restoreLocalFormVec();
    const int ierr = VecDestroy(d_global_vec);  IBTK_CHKERRQ(ierr);
    return;
}// ~LNodeLevelData

LNodeLevelData&
LNodeLevelData::operator=(
    const LNodeLevelData& that)
{
    if (that.d_in_local_form)
    {
        TBOX_ERROR("LNodeLevelData operator=():\n" <<
                   "  operator=() source object must NOT be in local form." << std::endl);
    }
    if (this == &that) return *this;  // check for self-assignment

    int ierr;
    restoreLocalFormVec();
    ierr = VecDestroy(d_global_vec);                        IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(that.d_global_vec, &d_global_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecCopy(that.d_global_vec, d_global_vec);        IBTK_CHKERRQ(ierr);

    d_name = that.d_name;
    d_depth = that.d_depth;
    d_nonlocal_petsc_indices = that.d_nonlocal_petsc_indices;
    d_local_vec = NULL;
    d_in_local_form = false;
    d_local_vec_array = NULL;
    d_extracted_local_array = false;
    return *this;
}// operator=

void
LNodeLevelData::putToDatabase(
    Pointer<Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    const bool restore_local_vec = !d_in_local_form;
    const bool restore_local_array = !d_extracted_local_array;

    const int num_local_nodes = getLocalNodeCount();
    const double* const local_form_array = getLocalFormArray();
    const int n_nonlocal_indices = d_nonlocal_petsc_indices.size();

    db->putString("d_name", d_name);
    db->putInteger("d_depth", d_depth);
    db->putInteger("num_local_nodes", num_local_nodes);
    db->putInteger("n_nonlocal_indices", n_nonlocal_indices);
    if (n_nonlocal_indices > 0)
    {
        db->putIntegerArray("d_nonlocal_petsc_indices",
                            d_nonlocal_petsc_indices.empty() ? NULL : &d_nonlocal_petsc_indices[0],
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
    const std::string& name,
    const int num_local_nodes,
    const int depth,
    const std::vector<int>& nonlocal_petsc_indices)
    : d_name(name),
      d_depth(depth),
      d_nonlocal_petsc_indices(nonlocal_petsc_indices),
      d_global_vec(NULL),
      d_local_vec(NULL),
      d_in_local_form(false),
      d_local_vec_array(NULL),
      d_extracted_local_array(false)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(num_local_nodes >= 0);
    TBOX_ASSERT(depth > 0);
#endif
    // Create the PETSc Vec which actually provides the storage for the
    // Lagrangian data.
    int ierr;
    if (d_depth == 1)
    {
        ierr = VecCreateGhost(
            PETSC_COMM_WORLD,
            num_local_nodes, PETSC_DECIDE,
            d_nonlocal_petsc_indices.size(),
            d_nonlocal_petsc_indices.empty() ? NULL : &d_nonlocal_petsc_indices[0],
            &d_global_vec);  IBTK_CHKERRQ(ierr);
    }
    else
    {
        ierr = VecCreateGhostBlock(
            PETSC_COMM_WORLD, d_depth,
            d_depth*num_local_nodes, PETSC_DECIDE,
            d_nonlocal_petsc_indices.size(),
            d_nonlocal_petsc_indices.empty() ? NULL : &d_nonlocal_petsc_indices[0],
            &d_global_vec);  IBTK_CHKERRQ(ierr);
    }
    ierr = VecSetBlockSize(d_global_vec, d_depth);  IBTK_CHKERRQ(ierr);

    return;
}// LNodeLevelData

LNodeLevelData::LNodeLevelData(
    Pointer<Database> db)
    : d_name(db->getString("d_name")),
      d_depth(db->getInteger("d_depth")),
      d_nonlocal_petsc_indices(),
      d_global_vec(NULL),
      d_local_vec(NULL),
      d_in_local_form(false),
      d_local_vec_array(NULL),
      d_extracted_local_array(false)
{
    int num_local_nodes = db->getInteger("num_local_nodes");
    int n_nonlocal_indices = db->getInteger("n_nonlocal_indices");
    d_nonlocal_petsc_indices.resize(n_nonlocal_indices);
    if (n_nonlocal_indices > 0)
    {
        db->getIntegerArray("d_nonlocal_petsc_indices",
                            d_nonlocal_petsc_indices.empty() ? NULL : &d_nonlocal_petsc_indices[0],
                            n_nonlocal_indices);
    }

    // Create the PETSc Vec which actually provides the storage for the
    // Lagrangian data.
    int ierr;

    if (d_depth == 1)
    {
        ierr = VecCreateGhost(PETSC_COMM_WORLD,
                              num_local_nodes, PETSC_DECIDE,
                              d_nonlocal_petsc_indices.size(),
                              d_nonlocal_petsc_indices.empty() ? NULL : &d_nonlocal_petsc_indices[0],
                              &d_global_vec);
        IBTK_CHKERRQ(ierr);
    }
    else
    {
        ierr = VecCreateGhostBlock(PETSC_COMM_WORLD, d_depth,
                                   d_depth*num_local_nodes, PETSC_DECIDE,
                                   d_nonlocal_petsc_indices.size(),
                                   d_nonlocal_petsc_indices.empty() ? NULL : &d_nonlocal_petsc_indices[0],
                                   &d_global_vec);
        IBTK_CHKERRQ(ierr);
    }

    ierr = VecSetBlockSize(d_global_vec, d_depth);  IBTK_CHKERRQ(ierr);

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
    : d_name(from.d_name),
      d_depth(from.d_depth),
      d_nonlocal_petsc_indices(from.d_nonlocal_petsc_indices),
      d_global_vec(NULL),
      d_local_vec(NULL),
      d_in_local_form(false),
      d_local_vec_array(NULL),
      d_extracted_local_array(false)
{
    if (from.d_in_local_form)
    {
        TBOX_ERROR("LNodeLevelData copy constructor:\n" <<
                   "  copy constructor source object must NOT be in local form." << std::endl);
    }

    int ierr;
    ierr = VecDuplicate(from.d_global_vec, &d_global_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecCopy(from.d_global_vec, d_global_vec);        IBTK_CHKERRQ(ierr);

    return;
}// LNodeLevelData

void
LNodeLevelData::resetData(
    Vec& new_global_vec,
    const std::vector<int>& new_nonlocal_petsc_indices)
{
    if (d_in_local_form)
    {
        TBOX_ERROR("LNodeLevelData::resetData()\n" <<
                   "  object must NOT be in local form." << std::endl);
    }
    d_global_vec = new_global_vec;
    d_nonlocal_petsc_indices = new_nonlocal_petsc_indices;
    return;
}// resetData

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBTK::LNodeLevelData>;

//////////////////////////////////////////////////////////////////////////////
