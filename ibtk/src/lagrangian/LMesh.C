// Filename: LMesh.C
// Created on 05 May 2011 by Boyce Griffith
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

#include "LMesh.h"

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
#include <ibtk/namespaces.h>

// C++ STDLIB INCLUDES
#include <algorithm>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LMesh::LMesh(
    const std::string& object_name)
    : d_object_name(object_name),
      d_nodes(),
      d_lag_idxs(),
      d_global_petsc_idxs(),
      d_local_petsc_idxs(),
      d_ghost_nodes(),
      d_ghost_lag_idxs(),
      d_ghost_global_petsc_idxs(),
      d_ghost_local_petsc_idxs()
{
    // intentionally blank
    return;
}// LMesh

LMesh::~LMesh()
{
    // intentionally blank
    return;
}// ~LMesh

void
LMesh::setNodes(
    const std::vector<LNode*>& nodes,
    const bool sorted)
{
    d_nodes = nodes;
    if (!sorted) std::sort(d_nodes.begin(), d_nodes.end(), LNodeIndexLocalPETScIndexComp());

    int node_counter = 0;
    const int num_nodes = d_nodes.size();
    d_lag_idxs.resize(num_nodes);
    d_global_petsc_idxs.resize(num_nodes);
    d_local_petsc_idxs.resize(num_nodes);
    if (num_nodes > 0)
    {
        for (std::vector<LNode*>::iterator it = d_nodes.begin(); it != d_nodes.end();
             ++it, ++node_counter)
        {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(node_counter == (*it)->getLocalPETScIndex());
#endif
            d_lag_idxs         [node_counter] = (*it)->getLagrangianIndex();
            d_global_petsc_idxs[node_counter] = (*it)->getGlobalPETScIndex();
            d_local_petsc_idxs [node_counter] = (*it)->getLocalPETScIndex();
        }
    }

    int ghost_node_counter = 0;
    const int num_ghost_nodes = d_ghost_nodes.size();
    if (num_ghost_nodes > 0)
    {
        for (std::vector<LNode*>::iterator it = d_ghost_nodes.begin(); it != d_ghost_nodes.end();
             ++it, ++node_counter, ++ghost_node_counter)
        {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(node_counter == (*it)->getLocalPETScIndex());
#endif
            d_ghost_lag_idxs         [ghost_node_counter] = (*it)->getLagrangianIndex();
            d_ghost_global_petsc_idxs[ghost_node_counter] = (*it)->getGlobalPETScIndex();
            d_ghost_local_petsc_idxs [ghost_node_counter] = (*it)->getLocalPETScIndex();
        }
    }
    return;
}// setNodes

void
LMesh::setGhostNodes(
    const std::vector<LNode*>& ghost_nodes,
    const bool sorted)
{
    d_ghost_nodes = ghost_nodes;
    if (!sorted) std::sort(d_ghost_nodes.begin(), d_ghost_nodes.end(), LNodeIndexLocalPETScIndexComp());

    int node_counter = d_nodes.size();
    int ghost_node_counter = 0;
    const int num_ghost_nodes = d_ghost_nodes.size();
    d_ghost_lag_idxs.resize(num_ghost_nodes);
    d_ghost_global_petsc_idxs.resize(num_ghost_nodes);
    d_ghost_local_petsc_idxs.resize(num_ghost_nodes);
    if (num_ghost_nodes > 0)
    {
        for (std::vector<LNode*>::iterator it = d_ghost_nodes.begin(); it != d_ghost_nodes.end();
             ++it, ++node_counter, ++ghost_node_counter)
        {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(node_counter == (*it)->getLocalPETScIndex());
#endif
            d_ghost_lag_idxs         [ghost_node_counter] = (*it)->getLagrangianIndex();
            d_ghost_global_petsc_idxs[ghost_node_counter] = (*it)->getGlobalPETScIndex();
            d_ghost_local_petsc_idxs [ghost_node_counter] = (*it)->getLocalPETScIndex();
        }
    }
    return;
}// setGhostNodes

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBTK::LMesh>;

//////////////////////////////////////////////////////////////////////////////
