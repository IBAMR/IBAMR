// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_boxpartitioner
#define included_IBTK_boxpartitioner

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#ifdef IBTK_HAVE_LIBMESH

#include <ibtk/PartitioningBox.h>

#include <libmesh/mesh_base.h>
#include <libmesh/partitioner.h>
#include <libmesh/system.h>

#include <memory>
#include <string>

namespace libMesh
{
class MeshBase;
class System;
} // namespace libMesh

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
    BoxPartitioner(const PartitioningBoxes& partitioning_boxes, const libMesh::System& position_system);

    /*!
     * \brief Enable or disable logging.
     */
    void setLoggingEnabled(bool enable_logging = true);

    /*!
     * \brief Determine whether logging is enabled or disabled.
     */
    bool getLoggingEnabled() const;

    /*!
     * Write the partitioning to a file in a simple point-based format: for
     * each patch several points are printed to the specified file in the
     * format
     *
     * x,y,z,r
     *
     * format.
     */
    void writePartitioning(const std::string& file_name) const;

    virtual std::unique_ptr<libMesh::Partitioner> clone() const override;

protected:
    /// The function used to actually do the partitioning.
    virtual void _do_partition(libMesh::MeshBase& mesh, const unsigned int n) override;

    /// Logging configuration.
    bool d_enable_logging = false;

    /// The PartitioningBoxes object used to establish whether or not an Elem
    /// (via its centroid) or Node is owned by the current processor.
    PartitioningBoxes d_partitioning_boxes;

    /// Pointer, if relevant, to the libMesh mesh position system.
    const libMesh::System* const d_position_system = nullptr;
};
} // namespace IBTK
//////////////////////////////////////////////////////////////////////////////
#endif //#ifdef IBTK_HAVE_LIBMESH
#endif //#ifndef included_IBTK_boxpartitioner
