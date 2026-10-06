// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2026 by the IBAMR developers
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
#include <ibtk/IBTK_MPI.h>
#include <ibtk/StableCentroidPartitioner.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/PIO.h>

#include <libmesh/elem.h>
#include <libmesh/id_types.h>
#include <libmesh/libmesh_config.h>
#include <libmesh/libmesh_version.h>
#include <libmesh/mesh_base.h>
#include <libmesh/point.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

std::unique_ptr<Partitioner>
StableCentroidPartitioner::clone() const
{
    return std::make_unique<StableCentroidPartitioner>();
} // clone

/////////////////////////////// PROTECTED ////////////////////////////////////

void
StableCentroidPartitioner::_do_partition(MeshBase& mesh, const unsigned int n)
{
    TBOX_ASSERT(mesh.is_replicated());
    // only implemented when we use SAMRAI's partitioning
    TBOX_ASSERT(n == static_cast<unsigned int>(IBTK_MPI::getNodes()));

    // Every process holds the whole mesh, so every process finds the same bounding box of the centroids without any
    // communication.
    std::vector<std::pair<std::array<float, LIBMESH_DIM>, libMesh::Elem*>> centroids;
    std::array<double, LIBMESH_DIM> x_min;
    std::array<double, LIBMESH_DIM> x_max;
    x_min.fill(std::numeric_limits<double>::max());
    x_max.fill(std::numeric_limits<double>::lowest());
    auto el_end = mesh.elements_end();
    for (auto it = mesh.elements_begin(); it != el_end; ++it)
    {
        const libMesh::Point centroid = (*it)->vertex_average();

        std::array<float, LIBMESH_DIM> rounded_centroid = { 0.0f };
        for (unsigned int d = 0; d < LIBMESH_DIM; ++d)
        {
            rounded_centroid[d] = centroid(d);
            x_min[d] = std::min(x_min[d], centroid(d));
            x_max[d] = std::max(x_max[d], centroid(d));
        }

        centroids.push_back(std::make_pair(rounded_centroid, *it));
    }

    // Set centroid coordinates that are zero up to rounding errors to zero so that they compare equal. The threshold
    // is relative to the largest extent of the centroids so that it does not depend on the units of the mesh.
    double extent = 0.0;
    for (unsigned int d = 0; d < LIBMESH_DIM; ++d)
    {
        extent = std::max(extent, x_max[d] - x_min[d]);
    }
    const auto zero_threshold = static_cast<float>(std::numeric_limits<float>::epsilon() * extent);
    for (auto& centroid : centroids)
    {
        for (float& x : centroid.first)
        {
            if (std::abs(x) < zero_threshold)
            {
                x = 0.0f;
            }
        }
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
