// Filename: SAMRAIPartitioner.h
// Created on 6 Nov 2018 by David Wells
//
// Copyright (c) 2018, Boyce Griffith, David Wells
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

#ifndef included_IBTK_ibtk_samraipartitioner
#define included_IBTK_ibtk_samraipartitioner

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/libmesh_utilities.h"
#include <ibtk/BoundingBox.h>

#include <libmesh/mesh_base.h>
#include <libmesh/partitioner.h>

#include <PatchHierarchy.h>

#include <mpi.h>

#include <algorithm>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
// TODO: move all code to a proper source file
class SAMRAIPartitioner : public libMesh::Partitioner
{
public:
    SAMRAIPartitioner(const PatchHierarchy<NDIM>& patch_hierarchy)
    {
        std::vector<IBTK::BoundingBox> boxes;
        const int finest_level = patch_hierarchy.getFinestLevelNumber();
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy.getPatchLevel(finest_level);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const Patch<NDIM>& patch = *level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geometry = patch.getPatchGeometry();
            boxes.emplace_back(*patch_geometry);
        }

        d_bounding_boxes = BoundingBoxes(boxes.begin(), boxes.end());
    }

    SAMRAIPartitioner(const BoundingBoxes& bounding_boxes) : d_bounding_boxes(bounding_boxes)
    {
    }

    virtual std::unique_ptr<libMesh::Partitioner> clone() const override
    {
        return std::unique_ptr<libMesh::Partitioner>(new SAMRAIPartitioner(d_bounding_boxes));
    }

    virtual void _do_partition(MeshBase& mesh, const unsigned int n) override
    {
        // We assume every cell is on every processor
        TBOX_ASSERT(mesh.is_replicated());
        // only implemented when we use SAMRAI's partitioning
        TBOX_ASSERT(n == static_cast<unsigned int>(SAMRAI_MPI::getNodes()));
        const int current_rank = SAMRAI_MPI::getRank();

        auto to_ibtk_point = [](const libMesh::Point& p) -> IBTK::Point {
            IBTK::Point point;
            for (unsigned int d = 0; d < NDIM; ++d) point[d] = p(d);
            return point;
        };

        // TODO this partitioning algorithm is not scalable since every processor
        // looks at every cell.
        std::vector<dof_id_type> local_elems;
        std::vector<dof_id_type> local_nodes;

        unsigned int elem_n = 0;
        unsigned int node_n = 0;
        for (auto el = mesh.active_elements_begin(); el != mesh.active_elements_end(); ++el)
        {
            auto& elem = *el;
            const IBTK::Point centroid = to_ibtk_point(elem->centroid());
            const unsigned int n_nodes = elem->n_nodes();
            for (unsigned int elem_node_n = 0; elem_node_n < n_nodes; ++elem_node_n)
            {
                const IBTK::Point vertex = to_ibtk_point(elem->point(elem_node_n));

                if (d_bounding_boxes.point_inside(vertex))
                {
                    elem->node_ptr(elem_node_n)->processor_id() = cast_int<processor_id_type>(current_rank);
                    local_nodes.push_back(node_n);
                }

                ++node_n;
            }

            if (d_bounding_boxes.point_inside(centroid))
            {
                elem->processor_id() = cast_int<processor_id_type>(current_rank);
                local_elems.push_back(elem_n);
            }
            ++elem_n;
        }

        // convert the libMesh type to an MPI type
        MPI_Datatype integral_type = 0;
        switch (sizeof(dof_id_type))
        {
        case 1:
            integral_type = MPI_UNSIGNED_CHAR;
            break;
        case 2:
            integral_type = MPI_UNSIGNED_SHORT;
            break;
        case 4:
            integral_type = MPI_UNSIGNED;
            break;
        case 8:
            integral_type = MPI_UNSIGNED_LONG;
            break;
        }

        // At this point we have partitioned the mesh (i.e., all Nodes and Elems
        // have been assigned a processor id). We now check our work.
        //
        // 1. Verify that we partitioned each Elem uniquely:
        {
            std::vector<processor_id_type> elem_ids(mesh.parallel_n_elem());
            for (const dof_id_type id : local_elems)
                // Ensure that we don't partition elements on rank 0 multiple times by
                // adding 1 to all ranks.
                elem_ids[id] = current_rank + 1;

            const int ierr = MPI_Allreduce(
                elem_ids.data(), elem_ids.data(), elem_ids.size(), integral_type, MPI_SUM, SAMRAI_MPI::commWorld);
            TBOX_ASSERT(ierr == 0);

            // Verify that no cells were assigned to two subdomains:
            for (const dof_id_type id : local_elems)
            {
                TBOX_ASSERT(elem_ids[id] == current_rank + 1);
            }
            // Verify that all cells were partitioned
            TBOX_ASSERT(std::find(elem_ids.begin(), elem_ids.end(), 0) == elem_ids.end());
        }

        // 2. Verify that we partitioned each Node uniquely:
        {
            std::vector<processor_id_type> node_ids(mesh.parallel_n_nodes());
            for (const dof_id_type id : local_nodes) node_ids[id] = current_rank + 1;

            const int ierr = MPI_Allreduce(
                node_ids.data(), node_ids.data(), node_ids.size(), integral_type, MPI_SUM, SAMRAI_MPI::commWorld);
            TBOX_ASSERT(ierr == 0);

            // Verify that no nodes were assigned to two subdomains:
            for (const dof_id_type id : local_nodes)
            {
                TBOX_ASSERT(node_ids[id] == current_rank + 1);
            }
            // Verify that all nodes were partitioned:
            TBOX_ASSERT(std::find(node_ids.begin(), node_ids.end(), 0) == node_ids.end());
        }
    }

protected:
    BoundingBoxes d_bounding_boxes;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_ibtk_samraipartitioner
