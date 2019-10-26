// Filename: LData.cpp
// Created on 17 Apr 2004 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
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

#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/LData.h"
#include "ibtk/ibtk_macros.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include "petscsys.h"
#include "petscvec.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include "boost/multi_array.hpp"
IBTK_ENABLE_EXTRA_WARNINGS

#include <ostream>
#include <string>
#include <utility>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LData::LData(std::string name,
             const unsigned int num_local_nodes,
             const unsigned int depth,
             std::vector<int> nonlocal_petsc_indices)
    : d_name(std::move(name)), d_depth(depth), d_nonlocal_petsc_indices(std::move(nonlocal_petsc_indices))
{
    // Create the PETSc Vec that provides storage for the Lagrangian data.
    d_global_vec = setup_petsc_vector(num_local_nodes, d_depth, d_nonlocal_petsc_indices);

    int global_node_count = 0;
    const int ierr = VecGetSize(d_global_vec, &global_node_count);
    IBTK_CHKERRQ(ierr);
#if !defined(NDEBUG)
    TBOX_ASSERT(global_node_count >= 0);
#endif
    d_global_node_count = global_node_count;
    d_global_node_count /= d_depth;
    d_local_node_count = num_local_nodes;
    d_ghost_node_count = static_cast<int>(d_nonlocal_petsc_indices.size());
    return;
} // LData

LData::LData(std::string name, Vec vec, std::vector<int> nonlocal_petsc_indices, const bool manage_petsc_vec)
    : d_name(std::move(name)),
      d_nonlocal_petsc_indices(std::move(nonlocal_petsc_indices)),
      d_global_vec(vec),
      d_managing_petsc_vec(manage_petsc_vec)
{
    int ierr;
    int depth;
    ierr = VecGetBlockSize(d_global_vec, &depth);
    IBTK_CHKERRQ(ierr);
#if !defined(NDEBUG)
    TBOX_ASSERT(depth >= 0);
#endif
    d_depth = depth;
    int global_node_count;
    ierr = VecGetSize(d_global_vec, &global_node_count);
    IBTK_CHKERRQ(ierr);
#if !defined(NDEBUG)
    TBOX_ASSERT(global_node_count >= 0);
#endif
    d_global_node_count = global_node_count;
    d_global_node_count /= d_depth;
    int local_node_count;
    ierr = VecGetLocalSize(d_global_vec, &local_node_count);
    IBTK_CHKERRQ(ierr);
#if !defined(NDEBUG)
    TBOX_ASSERT(local_node_count >= 0);
#endif
    d_local_node_count = local_node_count;
    d_local_node_count /= d_depth;
    d_ghost_node_count = static_cast<int>(d_nonlocal_petsc_indices.size());
    return;
} // LData

LData::LData(Pointer<Database> db) : d_name(db->getString("d_name")), d_depth(db->getInteger("d_depth"))
{
    int num_local_nodes = db->getInteger("num_local_nodes");
    int num_ghost_nodes = db->getInteger("num_ghost_nodes");
    d_nonlocal_petsc_indices.resize(num_ghost_nodes);
    if (num_ghost_nodes > 0)
    {
        db->getIntegerArray("d_nonlocal_petsc_indices", d_nonlocal_petsc_indices.data(), num_ghost_nodes);
    }

    d_global_vec = setup_petsc_vector(num_local_nodes, d_depth, d_nonlocal_petsc_indices);

    int global_node_count;
    int ierr = VecGetSize(d_global_vec, &global_node_count);
    IBTK_CHKERRQ(ierr);
#if !defined(NDEBUG)
    TBOX_ASSERT(global_node_count >= 0);
#endif
    d_global_node_count = global_node_count;
    d_global_node_count /= d_depth;
    d_local_node_count = num_local_nodes;
    d_ghost_node_count = static_cast<int>(d_nonlocal_petsc_indices.size());

    // Extract the values from the database.
    double* ghosted_local_vec_array = getGhostedLocalFormVecArray()->data();
    if (num_local_nodes + num_ghost_nodes > 0)
    {
        db->getDoubleArray("vals", ghosted_local_vec_array, d_depth * (num_local_nodes + num_ghost_nodes));
    }
    restoreArrays();
    return;
} // LData

LData::~LData()
{
    restoreArrays();
    if (d_managing_petsc_vec)
    {
        const int ierr = VecDestroy(&d_global_vec);
        IBTK_CHKERRQ(ierr);
    }
    return;
} // ~LData

void
LData::resetData(Vec vec, const std::vector<int>& nonlocal_petsc_indices, const bool manage_petsc_vec)
{
    restoreArrays();
    int ierr;
    if (d_managing_petsc_vec)
    {
        ierr = VecDestroy(&d_global_vec);
        IBTK_CHKERRQ(ierr);
    }

    // Take ownership of new Vec
    d_global_vec = vec;
    d_managing_petsc_vec = manage_petsc_vec;

    int depth;
    ierr = VecGetBlockSize(d_global_vec, &depth);
    IBTK_CHKERRQ(ierr);
#if !defined(NDEBUG)
    TBOX_ASSERT(depth >= 0);
#endif
    d_depth = depth;
    int global_node_count;
    ierr = VecGetSize(d_global_vec, &global_node_count);
    IBTK_CHKERRQ(ierr);
#if !defined(NDEBUG)
    TBOX_ASSERT(global_node_count >= 0);
#endif
    d_global_node_count = global_node_count;
    d_global_node_count /= d_depth;
    int local_node_count;
    ierr = VecGetLocalSize(d_global_vec, &local_node_count);
    IBTK_CHKERRQ(ierr);
#if !defined(NDEBUG)
    TBOX_ASSERT(local_node_count >= 0);
#endif
    d_local_node_count = local_node_count;
    d_local_node_count /= d_depth;
    d_nonlocal_petsc_indices = nonlocal_petsc_indices;
    d_ghost_node_count = static_cast<int>(d_nonlocal_petsc_indices.size());
    return;
} // resetData

void
LData::putToDatabase(Pointer<Database> db)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(db);
#endif
    const int num_local_nodes = getLocalNodeCount();
    const int num_ghost_nodes = static_cast<int>(d_nonlocal_petsc_indices.size());
    db->putString("d_name", d_name);
    db->putInteger("d_depth", d_depth);
    db->putInteger("num_local_nodes", num_local_nodes);
    db->putInteger("num_ghost_nodes", num_ghost_nodes);
    if (num_ghost_nodes > 0)
    {
        db->putIntegerArray("d_nonlocal_petsc_indices", &d_nonlocal_petsc_indices[0], num_ghost_nodes);
    }
    const double* const ghosted_local_vec_array = getGhostedLocalFormVecArray()->data();
    if (num_local_nodes + num_ghost_nodes > 0)
    {
        db->putDoubleArray("vals", ghosted_local_vec_array, d_depth * (num_local_nodes + num_ghost_nodes));
    }
    restoreArrays();
    return;
} // putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
