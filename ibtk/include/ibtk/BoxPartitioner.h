// Filename: BoxPartitioner.h
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
#ifndef included_IBTK_ibtk_boxpartitioner
#define included_IBTK_ibtk_boxpartitioner

/////////////////////////////// INCLUDES /////////////////////////////////////
#include <ibtk/PartitioningBox.h>

#include <libmesh/partitioner.h>
#include <libmesh/mesh_base.h>
#include <libmesh/system.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////
namespace IBTK
{
/*!
 * @brief A libMesh partitioner that partitions a mesh based on
 * PartitioningBox objects owned by each processor.
 *
 * IBAMR assumes that all Lagrangian structures exist on the finest level of
 * the fluid mesh. This partitioner uses the data stored by the SAMRAI
 * PatchHierarchy object (which is used, under most circumstances, to create
 * the PartitioningBoxes object) to compute partitioning boxes corresponding
 * to the locally owned fluid patches on the finest level: since these
 * partition the part of the domain over which the structure exists, they can
 * be used to partition the libMesh Mesh object across the MPI network. Put
 * another way: this Partitioner uses Eulerian grid data to partition the
 * structural meshes.
 */
class BoxPartitioner : public libMesh::Partitioner
{
public:
    /// Constructor from a PartitioningBoxes instance.
    BoxPartitioner(const PartitioningBoxes& partitioning_boxes);

    /*!
     * Constructor. This is like the other constructor that takes a
     * PartitioningBoxes object, but it permits the use of a background mesh
     * that is displaced by a vector finite element field.
     *
     * @param partitioning_boxes the boxes comprising the partition.
     *
     * @param position_system the libMesh::System object whose current
     * solution is the position of the Mesh which will subsequently be
     * partitioned.
     */
    BoxPartitioner(const PartitioningBoxes& partitioning_boxes,
                   const libMesh::System &position_system);

    virtual std::unique_ptr<libMesh::Partitioner> clone() const override;

    virtual void _do_partition(libMesh::MeshBase& mesh, const unsigned int n) override;

protected:
    /// The PartitioningBoxes object used to establish whether or not an Elem
    /// (via its centroid) or Node is owned by the current processor.
    PartitioningBoxes d_partitioning_boxes;

    /// Pointer, if relevant, to the libMesh mesh position system.
    const libMesh::System* const d_position_system = nullptr;
};
} // namespace IBTK
//////////////////////////////////////////////////////////////////////////////
#endif //#ifndef included_IBTK_ibtk_boxpartitioner
