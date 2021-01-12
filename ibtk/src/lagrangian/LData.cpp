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

#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/LData.h"

#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include "petscvec.h"

#include "ibtk/namespaces.h" // IWYU pragma: keep

IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/multi_array.hpp>
IBTK_ENABLE_EXTRA_WARNINGS

#include <algorithm>
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
