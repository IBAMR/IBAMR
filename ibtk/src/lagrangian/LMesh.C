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
#include <ibtk/FixedSizedStream.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <tbox/RestartManager.h>

// C++ STDLIB INCLUDES
#include <algorithm>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LMesh::LMesh(
    const std::string& object_name,
    const bool register_for_restart)
    : d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_nodes(),
      d_lag_idxs(),
      d_global_petsc_idxs(),
      d_local_petsc_idxs(),
      d_ghost_nodes(),
      d_ghost_lag_idxs(),
      d_ghost_global_petsc_idxs(),
      d_ghost_local_petsc_idxs()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    return;
}// LMesh

LMesh::~LMesh()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }
    return;
}// ~LMesh

void
LMesh::setNodes(
    const std::vector<LNode>& nodes,
    const bool sorted)
{
    d_nodes = nodes;
    if (!sorted) std::sort(d_nodes.begin(), d_nodes.end());

    int node_counter = 0;
    const int num_nodes = d_nodes.size();
    d_lag_idxs.resize(num_nodes);
    d_global_petsc_idxs.resize(num_nodes);
    d_local_petsc_idxs.resize(num_nodes);
    if (num_nodes > 0)
    {
        for (std::vector<LNode>::iterator it = d_nodes.begin(); it != d_nodes.end();
             ++it, ++node_counter)
        {
            it->setLocalPETScIndex(node_counter);
            d_lag_idxs         [node_counter] = it->getLagrangianIndex();
            d_global_petsc_idxs[node_counter] = it->getGlobalPETScIndex();
            d_local_petsc_idxs [node_counter] = it->getLocalPETScIndex();
        }
    }

    int ghost_node_counter = 0;
    const int num_ghost_nodes = d_ghost_nodes.size();
    if (num_ghost_nodes > 0)
    {
        for (std::vector<LNode>::iterator it = d_ghost_nodes.begin(); it != d_ghost_nodes.end();
             ++it, ++node_counter, ++ghost_node_counter)
        {
            it->setLocalPETScIndex(node_counter);
            d_ghost_lag_idxs         [ghost_node_counter] = it->getLagrangianIndex();
            d_ghost_global_petsc_idxs[ghost_node_counter] = it->getGlobalPETScIndex();
            d_ghost_local_petsc_idxs [ghost_node_counter] = it->getLocalPETScIndex();
        }
    }
    return;
}// setNodes

void
LMesh::setGhostNodes(
    const std::vector<LNode>& ghost_nodes,
    const bool sorted)
{
    d_ghost_nodes = ghost_nodes;
    if (!sorted) std::sort(d_ghost_nodes.begin(), d_ghost_nodes.end());

    int node_counter = d_nodes.size();
    int ghost_node_counter = 0;
    const int num_ghost_nodes = d_ghost_nodes.size();
    d_ghost_lag_idxs.resize(num_ghost_nodes);
    d_ghost_global_petsc_idxs.resize(num_ghost_nodes);
    d_ghost_local_petsc_idxs.resize(num_ghost_nodes);
    if (num_ghost_nodes > 0)
    {
        for (std::vector<LNode>::iterator it = d_ghost_nodes.begin(); it != d_ghost_nodes.end();
             ++it, ++node_counter, ++ghost_node_counter)
        {
            it->setLocalPETScIndex(node_counter);
            d_ghost_lag_idxs         [ghost_node_counter] = it->getLagrangianIndex();
            d_ghost_global_petsc_idxs[ghost_node_counter] = it->getGlobalPETScIndex();
            d_ghost_local_petsc_idxs [ghost_node_counter] = it->getLocalPETScIndex();
        }
    }
    return;
}// setGhostNodes

void
LMesh::putToDatabase(
    Pointer<Database> db)
{
    const int num_nodes = d_nodes.size();
    db->putInteger("num_nodes", num_nodes);
    if (num_nodes > 0)
    {
        size_t data_sz = 0;
        for (std::vector<LNode>::const_iterator cit = d_nodes.begin(); cit != d_nodes.end(); ++cit)
        {
            data_sz += cit->getDataStreamSize();
        }
        FixedSizedStream stream(data_sz);
        for (std::vector<LNode>::iterator it = d_nodes.begin(); it != d_nodes.end(); ++it)
        {
            it->packStream(stream);
        }
        db->putInteger("data_sz", data_sz);
        db->putCharArray("data", static_cast<char*>(stream.getBufferStart()), data_sz);
    }

    const int num_ghost_nodes = d_ghost_nodes.size();
    db->putInteger("num_ghost_nodes", num_ghost_nodes);
    if (num_ghost_nodes > 0)
    {
        size_t ghost_data_sz = 0;
        for (std::vector<LNode>::const_iterator cit = d_ghost_nodes.begin(); cit != d_ghost_nodes.end(); ++cit)
        {
            ghost_data_sz += cit->getDataStreamSize();
        }
        FixedSizedStream ghost_stream(ghost_data_sz);
        for (std::vector<LNode>::iterator it = d_ghost_nodes.begin(); it != d_ghost_nodes.end(); ++it)
        {
            it->packStream(ghost_stream);
        }
        db->putInteger("ghost_data_sz", ghost_data_sz);
        db->putCharArray("ghost_data", static_cast<char*>(ghost_stream.getBufferStart()), ghost_data_sz);
    }
    return;
}// putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
LMesh::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR("Restart database corresponding to " << d_object_name << " not found in restart file.");
    }

    int node_counter = 0;
    const int num_nodes = db->getInteger("num_nodes");
    d_nodes.resize(num_nodes);
    d_lag_idxs.resize(num_nodes);
    d_global_petsc_idxs.resize(num_nodes);
    d_local_petsc_idxs.resize(num_nodes);
    if (num_nodes > 0)
    {
        const size_t data_sz = db->getInteger("data_sz");
        std::vector<char> data(data_sz);
        db->getCharArray("data", &data[0], data_sz);
        FixedSizedStream stream(&data[0], data_sz);
        for (std::vector<LNode>::iterator it = d_nodes.begin(); it != d_nodes.end();
             ++it, ++node_counter)
        {
            it->unpackStream(stream, IntVector<NDIM>(0));
            it->setLocalPETScIndex(node_counter);
            d_lag_idxs         [node_counter] = it->getLagrangianIndex();
            d_global_petsc_idxs[node_counter] = it->getGlobalPETScIndex();
            d_local_petsc_idxs [node_counter] = it->getLocalPETScIndex();
        }
    }

    int ghost_node_counter = 0;
    const int num_ghost_nodes = db->getInteger("num_ghost_nodes");
    d_ghost_nodes.resize(num_ghost_nodes);
    d_ghost_lag_idxs.resize(num_ghost_nodes);
    d_ghost_global_petsc_idxs.resize(num_ghost_nodes);
    d_ghost_local_petsc_idxs.resize(num_ghost_nodes);
    if (num_ghost_nodes > 0)
    {
        const size_t ghost_data_sz = db->getInteger("ghost_data_sz");
        std::vector<char> ghost_data(ghost_data_sz);
        db->getCharArray("ghost_data", &ghost_data[0], ghost_data_sz);
        FixedSizedStream ghost_stream(&ghost_data[0], ghost_data_sz);
        for (std::vector<LNode>::iterator it = d_ghost_nodes.begin(); it != d_ghost_nodes.end();
             ++it, ++node_counter, ++ghost_node_counter)
        {
            it->unpackStream(ghost_stream, IntVector<NDIM>(0));
            it->setLocalPETScIndex(node_counter);
            d_ghost_lag_idxs         [ghost_node_counter] = it->getLagrangianIndex();
            d_ghost_global_petsc_idxs[ghost_node_counter] = it->getGlobalPETScIndex();
            d_ghost_local_petsc_idxs [ghost_node_counter] = it->getLocalPETScIndex();
        }
    }
    return;
}// getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBTK::LMesh>;

//////////////////////////////////////////////////////////////////////////////
