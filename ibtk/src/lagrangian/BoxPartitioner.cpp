// Filename: BoxPartitioner.cpp
// Created on 6 Dec 2018 by David Wells
//
// Copyright (c) 2018, Boyce Griffith
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
#include <ibtk/BoxPartitioner.h>
#include <ibtk/PartitioningBox.h>

#include <libmesh/elem.h>
#include <libmesh/mesh_base.h>
#include <libmesh/node.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/point.h>

#include <tbox/SAMRAI_MPI.h>

#include <algorithm>
#include <vector>

#include <mpi.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

BoxPartitioner::BoxPartitioner(const PartitioningBoxes& bounding_boxes)
    : d_partitioning_boxes(bounding_boxes)
{
} // BoxPartitioner

BoxPartitioner::BoxPartitioner(const PartitioningBoxes& bounding_boxes,
                               const libMesh::System &position_system)
    : d_partitioning_boxes(bounding_boxes),
      d_position_system(&position_system)
{
} // BoxPartitioner

std::unique_ptr<libMesh::Partitioner> BoxPartitioner::clone() const
{
    return std::unique_ptr<libMesh::Partitioner>(new BoxPartitioner(d_partitioning_boxes,
                                                                    *d_position_system));
} // clone

void BoxPartitioner::_do_partition(libMesh::MeshBase& mesh, const unsigned int n)
{
    // We assume every cell is on every processor
    TBOX_ASSERT(mesh.is_replicated());
    // only implemented when we use SAMRAI's partitioning
    TBOX_ASSERT(n == static_cast<unsigned int>(SAMRAI::tbox::SAMRAI_MPI::getNodes()));

    // convert the libMesh type to an MPI type
    MPI_Datatype pid_integral_type = 0;
    switch (sizeof(libMesh::processor_id_type))
    {
    case 1:
        pid_integral_type = MPI_UNSIGNED_CHAR;
        break;
    case 2:
        pid_integral_type = MPI_UNSIGNED_SHORT;
        break;
    case 4:
        pid_integral_type = MPI_UNSIGNED;
        break;
    case 8:
        pid_integral_type = MPI_UNSIGNED_LONG;
        break;
    }

    const int current_rank = SAMRAI::tbox::SAMRAI_MPI::getRank();
    auto to_ibtk_point = [](const libMesh::Point& p) -> IBTK::Point
    {
        IBTK::Point point;
        for (unsigned int d = 0; d < NDIM; ++d)
            point[d] = p(d);
        return point;
    };

    // TODO this partitioning algorithm is not scalable since every processor
    // looks at every cell. This will ultimately be replaced by an approach
    // driven by SAMRAI, where SAMRAI propagates nodes and cells in a
    // particle-like fashion.
    //
    // Step 0: determine the current location of the Mesh nodes.
    const bool use_position_vector = d_position_system != nullptr;
    const unsigned int position_system_n = use_position_vector
        ? d_position_system->number() : 0;
    std::vector<double> position;
    if (use_position_vector)
    {
        TBOX_ASSERT(&d_position_system->get_mesh() == &mesh);
        libMesh::NumericVector<double> *position_solution
            = d_position_system->solution.get();
        position.resize(position_solution->size());
        position_solution->localize(position);
    }
    std::vector<libMesh::Point> node_positions(mesh.parallel_n_nodes());
    std::size_t node_n = 0;
    for (auto node_it = mesh.nodes_begin(); node_it != mesh.nodes_end(); ++node_it)
    {
        const libMesh::Node* const node = *node_it;
        TBOX_ASSERT(node->id() == node_n);
        libMesh::Point node_position;
        if (use_position_vector)
        {
            if (node->n_vars(position_system_n))
            {
                TBOX_ASSERT(node->n_vars(position_system_n) == NDIM);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    node_position(d)
                        = position[node->dof_number(position_system_n, d, 0)];
                }
            }
        }
        else
        {
            node_position = *node;
        }
        node_positions[node_n] = node_position;
        ++node_n;
    }
    TBOX_ASSERT(node_n == node_positions.size());

    // Step 1: determine which elements belong to which processor and
    // communicate the partitioning information across the network.
    //
    // We offset the processor ids by 1 so that we can check, by summation,
    // that each Elem and Node is given exactly one processor id.
    std::vector<libMesh::processor_id_type> elem_ids(mesh.parallel_n_elem());
    std::vector<libMesh::dof_id_type> local_elem_ids;
    std::size_t elem_n = 0;
    const auto end_elem = mesh.active_elements_end();
    for (auto elem = mesh.active_elements_begin(); elem != end_elem; ++elem)
    {
        libMesh::Point centroid;
        const unsigned int n_nodes = (*elem)->n_nodes();
        for (unsigned int node_n = 0; node_n < n_nodes; ++node_n)
        {
            centroid += node_positions[(*elem)->node_id(node_n)];
        }
        centroid *= 1.0/n_nodes;

        if (d_partitioning_boxes.contains(to_ibtk_point(centroid)))
        {
            local_elem_ids.push_back(elem_n);
            TBOX_ASSERT(elem_n < elem_ids.size());
            elem_ids[elem_n] = current_rank + 1;
        }
        ++elem_n;
    }

    std::vector<libMesh::processor_id_type> node_ids(mesh.parallel_n_nodes());
    std::vector<libMesh::dof_id_type> local_node_ids;
    node_n = 0;
    const auto end_node = mesh.active_nodes_end();
    for (auto node = mesh.active_nodes_begin(); node != end_node; ++node)
    {
        if (d_partitioning_boxes.contains(to_ibtk_point(node_positions[node_n])))
        {
            local_node_ids.push_back(node_n);
            TBOX_ASSERT(node_n < node_ids.size());
            node_ids[node_n] = current_rank + 1;
        }
        ++node_n;
    }

    // step 2: communicate the partitioning across all processors:
    int ierr = MPI_Allreduce(
        MPI_IN_PLACE, elem_ids.data(), elem_ids.size(), pid_integral_type,
        MPI_SUM, SAMRAI::tbox::SAMRAI_MPI::commWorld);
    TBOX_ASSERT(ierr == 0);
    ierr = MPI_Allreduce(
        MPI_IN_PLACE, node_ids.data(), node_ids.size(), pid_integral_type,
        MPI_SUM, SAMRAI::tbox::SAMRAI_MPI::commWorld);
    TBOX_ASSERT(ierr == 0);

    // step 3: verify that we partitioned each elem and node exactly once:
    for (const libMesh::dof_id_type id : local_elem_ids)
    {
        TBOX_ASSERT(id < elem_ids.size());
        TBOX_ASSERT(elem_ids[id] == current_rank + 1);
    }
    for (const libMesh::dof_id_type id : local_node_ids)
    {
        TBOX_ASSERT(id < node_ids.size());
        TBOX_ASSERT(node_ids[id] == current_rank + 1);
    }

    TBOX_ASSERT(std::find(elem_ids.begin(), elem_ids.end(), 0) == elem_ids.end());
    TBOX_ASSERT(std::find(node_ids.begin(), node_ids.end(), 0) == node_ids.end());

    // Step 4: label all elements and nodes with the correct processor id.
    elem_n = 0;
    for (auto elem = mesh.active_elements_begin(); elem != end_elem; ++elem)
    {
        TBOX_ASSERT(elem_n < elem_ids.size());
        (*elem)->processor_id() = elem_ids[elem_n] - 1;
        ++elem_n;
    }

    node_n = 0;
    for (auto node = mesh.active_nodes_begin(); node != end_node; ++node)
    {
        TBOX_ASSERT(node_n < node_ids.size());
        (*node)->processor_id() = node_ids[node_n] - 1;
        ++node_n;
    }
} // _do_partition

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////////////////////////////////////////////////////
