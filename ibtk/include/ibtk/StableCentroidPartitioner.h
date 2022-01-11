// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2021 by the IBAMR developers
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

#ifndef included_IBTK_ibtk_stablecentroidpartitioner
#define included_IBTK_ibtk_stablecentroidpartitioner

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
 * the x coordinate of each centroid in a numerically stable way.
 *
 * This partitioner is intended only for internal use in IBAMR (e.g., to
 * generate platform-independent partitionings of meshes). It differs from
 * libMesh::CentroidPartitioner in a few significant ways:
 * <ol>
 *   <li>Centroids are sorted lexically instead of by a specified coordinate.</li>
 *   <li>To improve stability the least significant bits are cleared by casting
 *   the coordinates to and from single precision.</li>
 *   <li>To improve stability a stable sorting algorithm is used.</li>
 */
class StableCentroidPartitioner : public libMesh::Partitioner
{
public:
    /// Constructor
    StableCentroidPartitioner() = default;

    virtual std::unique_ptr<libMesh::Partitioner> clone() const override;

protected:
    /// The function used to actually do the partitioning.
    virtual void _do_partition(libMesh::MeshBase& mesh, const unsigned int n) override;
};
} // namespace IBTK
//////////////////////////////////////////////////////////////////////////////
#endif //#ifdef IBTK_HAVE_LIBMESH
#endif //#ifndef included_IBTK_ibtk_stablecentroidpartitioner
