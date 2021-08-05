// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2020 by the IBAMR developers
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
#include "ibtk/IBTK_MPI.h"
#include "ibtk/ibtk_utilities.h"
#include <ibtk/StableCentroidPartitioner.h>

#include <tbox/PIO.h>

#include <libmesh/elem.h>
#include <libmesh/id_types.h>
#include <libmesh/libmesh_config.h>
#include <libmesh/mesh_base.h>
#include <libmesh/point.h>

#include <algorithm>
#include <array>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

std::unique_ptr<Partitioner>
StableCentroidPartitioner::clone() const
{
    return std::unique_ptr<Partitioner>(new StableCentroidPartitioner());
} // clone

/////////////////////////////// PROTECTED ////////////////////////////////////

void
StableCentroidPartitioner::_do_partition(MeshBase& mesh, const unsigned int n)
{
    // We assume every cell is on every processor: this function is only in
    // libMesh 1.2.0 and newer
#if 1 < LIBMESH_MINOR_VERSION
    TBOX_ASSERT(mesh.is_replicated());
#endif
    // only implemented when we use SAMRAI's partitioning
    TBOX_ASSERT(n == static_cast<unsigned int>(IBTK_MPI::getNodes()));

    std::vector<std::pair<std::array<float, LIBMESH_DIM>, libMesh::Elem*> > centroids;
    auto el_end = mesh.elements_end();
    for (auto it = mesh.elements_begin(); it != el_end; ++it)
    {
        const libMesh::Point centroid = (*it)->centroid();

        std::array<float, LIBMESH_DIM> rounded_centroid = { 0.0f };
        for (unsigned int d = 0; d < LIBMESH_DIM; ++d)
        {
            rounded_centroid[d] = centroid(d);
            // also try to deal with values near zero
            if (std::abs(rounded_centroid[d]) < std::numeric_limits<float>::epsilon()) rounded_centroid[d] = 0.0f;
        }

        centroids.push_back(std::make_pair(rounded_centroid, *it));
    }
    std::stable_sort(
        centroids.begin(),
        centroids.end(),
        [](const std::pair<std::array<float, LIBMESH_DIM>, libMesh::Elem*>& a,
           const std::pair<std::array<float, LIBMESH_DIM>, libMesh::Elem*>& b)
        { return std::lexicographical_compare(a.first.begin(), a.first.end(), b.first.begin(), b.first.end()); });

    // proceed as libMesh would with CentroidPartitioner:
    const auto target_size = std::size_t(centroids.size() / n);
    for (std::size_t elem_n = 0; elem_n < centroids.size(); ++elem_n)
        centroids[elem_n].second->processor_id() = std::min<libMesh::processor_id_type>(elem_n / target_size, n - 1);
} // _do_partition

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////
} // namespace IBTK

/////////////////////////////////////////////////////////////////////////////
