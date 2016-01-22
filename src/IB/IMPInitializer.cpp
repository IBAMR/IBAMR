// Filename: IMPInitializer.cpp
// Created on 16 Oct 2012 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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

#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <algorithm>
#include <map>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "boost/multi_array.hpp"
#include "ibamr/IMPInitializer.h"
#include "ibamr/MaterialPointSpec.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/IndexUtilities.h"
#include "ibtk/LData.h"
#include "ibtk/LIndexSetData.h"
#include "ibtk/LNode.h"
#include "ibtk/LNodeSet.h"
#include "ibtk/LNodeSetData.h"
#include "ibtk/LSiloDataWriter.h"
#include "ibtk/Streamable.h"
#include "ibtk/ibtk_utilities.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/elem.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_type.h"
#include "libmesh/id_types.h"
#include "libmesh/mesh_base.h"
#include "libmesh/point.h"
#include "libmesh/quadrature.h"
#include "libmesh/type_vector.h"
#include "libmesh/variant_filter_iterator.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

namespace IBTK
{
class LDataManager;
} // namespace IBTK

using namespace libMesh;

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const double MIN_POINTS = 1.0;
static const double POINT_FACTOR = 2.0;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IMPInitializer::IMPInitializer(const std::string& object_name,
                               Pointer<Database> input_db,
                               Pointer<PatchHierarchy<NDIM> > hierarchy,
                               Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
    : d_object_name(object_name),
      d_hierarchy(hierarchy),
      d_gridding_alg(gridding_alg),
      d_level_is_initialized(d_gridding_alg->getMaxLevels(), false),
      d_meshes(d_gridding_alg->getMaxLevels()),
      d_num_vertex(d_gridding_alg->getMaxLevels()),
      d_vertex_offset(d_gridding_alg->getMaxLevels()),
      d_vertex_posn(d_gridding_alg->getMaxLevels()),
      d_vertex_wgt(d_gridding_alg->getMaxLevels()),
      d_vertex_subdomain_id(d_gridding_alg->getMaxLevels()),
      d_silo_writer(NULL)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(input_db);
#endif

    // Register the specification objects with the StreamableManager class.
    MaterialPointSpec::registerWithStreamableManager();

    // Initialize object with data read from the input database.
    getFromInput(input_db);
    return;
} // IMPInitializer

IMPInitializer::~IMPInitializer()
{
    pout << d_object_name << ":  Deallocating initialization data.\n";
    return;
} // ~IMPInitializer

void
IMPInitializer::registerMesh(MeshBase* mesh, int level_number)
{
    const int max_levels = d_gridding_alg->getMaxLevels();
    if (level_number < 0) level_number = max_levels - 1;
    level_number = std::min(level_number, max_levels - 1);
    const unsigned int mesh_idx = static_cast<unsigned int>(d_meshes[level_number].size());
    d_meshes[level_number].push_back(mesh);

    // Compute the Cartesian grid spacing on the specified level of the mesh.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    const double* const dx0 = grid_geom->getDx();
    double dx[NDIM];
    std::copy(dx0, dx0 + NDIM, dx);
    for (int ln = 1; ln <= level_number; ++ln)
    {
        const IntVector<NDIM> ratio = d_gridding_alg->getRatioToCoarserLevel(ln);
        for (unsigned int d = 0; d < NDIM; ++d) dx[d] /= static_cast<double>(ratio(d));
    }
    const double dx_min = *std::min_element(dx, dx + NDIM);

    // Setup data structures for computing the positions of the material points
    // and weighting factors.
    const int dim = mesh->mesh_dimension();
    FEType fe_type(FIRST, LAGRANGE);
    AutoPtr<QBase> qrule = QBase::build(QGAUSS, dim, FIRST);
    AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));
    fe->attach_quadrature_rule(qrule.get());
    const std::vector<libMesh::Point>& q_point = fe->get_xyz();
    const std::vector<double>& JxW = fe->get_JxW();
    const MeshBase::const_element_iterator el_begin = mesh->active_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh->active_elements_end();

    // Count the number of material points.
    d_num_vertex[level_number].push_back(0);
    d_vertex_offset[level_number].push_back(0);
    if (mesh_idx > 0)
    {
        d_vertex_offset[level_number][mesh_idx] =
            d_vertex_offset[level_number][mesh_idx - 1] + d_num_vertex[level_number][mesh_idx - 1];
    }
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        const double hmax = elem->hmax();
        const int npts = std::max(MIN_POINTS, std::ceil(POINT_FACTOR * hmax / dx_min));
        const Order order = static_cast<Order>(std::min(2 * npts - 1, static_cast<int>(FORTYTHIRD)));
        if (order != qrule->get_order())
        {
            qrule = QBase::build(QGAUSS, dim, order);
            fe->attach_quadrature_rule(qrule.get());
        }
        fe->reinit(elem);
        d_num_vertex[level_number][mesh_idx] += qrule->n_points();
    }

    // Initialize the material points.
    d_vertex_posn[level_number].resize(d_vertex_posn[level_number].size() + 1);
    d_vertex_wgt[level_number].resize(d_vertex_wgt[level_number].size() + 1);
    d_vertex_subdomain_id[level_number].resize(d_vertex_subdomain_id[level_number].size() + 1);
    d_vertex_posn[level_number][mesh_idx].resize(d_num_vertex[level_number][mesh_idx]);
    d_vertex_wgt[level_number][mesh_idx].resize(d_num_vertex[level_number][mesh_idx]);
    d_vertex_subdomain_id[level_number][mesh_idx].resize(d_num_vertex[level_number][mesh_idx]);
    unsigned int k = 0;
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        const double hmax = elem->hmax();
        const int npts = std::max(MIN_POINTS, std::ceil(POINT_FACTOR * hmax / dx_min));
        const Order order = static_cast<Order>(std::min(2 * npts - 1, static_cast<int>(FORTYTHIRD)));
        if (order != qrule->get_order())
        {
            qrule = QBase::build(QGAUSS, dim, order);
            fe->attach_quadrature_rule(qrule.get());
        }
        fe->reinit(elem);
        for (unsigned int qp = 0; qp < qrule->n_points(); ++qp, ++k)
        {
            d_vertex_posn[level_number][mesh_idx][k] = q_point[qp];
            d_vertex_wgt[level_number][mesh_idx][k] = JxW[qp];
            d_vertex_subdomain_id[level_number][mesh_idx][k] = elem->subdomain_id();
        }
    }
    return;
} // registerMesh

void
IMPInitializer::registerLSiloDataWriter(Pointer<LSiloDataWriter> silo_writer)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(silo_writer);
#endif

    // Cache a pointer to the data writer.
    d_silo_writer = silo_writer;

    // Check to see if we are starting from a restart file.
    RestartManager* restart_manager = RestartManager::getManager();
    const bool is_from_restart = restart_manager->isFromRestart();

    // Initialize the Silo data writer only if we are not starting from a
    // restart file.
    if (!is_from_restart)
    {
        for (int ln = 0; ln < d_gridding_alg->getMaxLevels(); ++ln)
        {
            if (d_level_is_initialized[ln]) initializeLSiloDataWriter(ln);
        }
    }
    return;
} // registerLSiloDataWriter

bool
IMPInitializer::getLevelHasLagrangianData(const int level_number, const bool /*can_be_refined*/) const
{
    return !d_meshes[level_number].empty();
} // getLevelHasLagrangianData

unsigned int
IMPInitializer::computeGlobalNodeCountOnPatchLevel(const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                                   const int level_number,
                                                   const double /*init_data_time*/,
                                                   const bool /*can_be_refined*/,
                                                   const bool /*initial_time*/)
{
    return std::accumulate(d_num_vertex[level_number].begin(), d_num_vertex[level_number].end(), 0);
}

unsigned int
IMPInitializer::computeLocalNodeCountOnPatchLevel(const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                  const int level_number,
                                                  const double /*init_data_time*/,
                                                  const bool can_be_refined,
                                                  const bool /*initial_time*/)
{
    // Loop over all patches in the specified level of the patch level and count
    // the number of local vertices.
    int local_node_count = 0;
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());

        // Count the number of vertices whose initial locations will be within
        // the given patch.
        std::vector<std::pair<int, int> > patch_vertices;
        getPatchVertices(patch_vertices, patch, level_number, can_be_refined);
        local_node_count += patch_vertices.size();
    }
    return local_node_count;
} // computeLocalNodeCountOnPatchLevel

void
IMPInitializer::initializeStructureIndexingOnPatchLevel(
    std::map<int, std::string>& strct_id_to_strct_name_map,
    std::map<int, std::pair<int, int> >& strct_id_to_lag_idx_range_map,
    const int level_number,
    const double /*init_data_time*/,
    const bool /*can_be_refined*/,
    const bool /*initial_time*/,
    LDataManager* const /*l_data_manager*/)
{
    int offset = 0;
    for (unsigned int j = 0; j < d_meshes[level_number].size(); ++j)
    {
        std::ostringstream name_stream;
        name_stream << "mesh_" << j;
        strct_id_to_strct_name_map[j] = name_stream.str();
        strct_id_to_lag_idx_range_map[j] = std::make_pair(offset, offset + d_num_vertex[level_number][j]);
        offset += d_num_vertex[level_number][j];
    }
    return;
} // initializeStructureIndexingOnPatchLevel

unsigned int
IMPInitializer::initializeDataOnPatchLevel(const int lag_node_index_idx,
                                           const unsigned int global_index_offset,
                                           const unsigned int local_index_offset,
                                           Pointer<LData> X_data,
                                           Pointer<LData> U_data,
                                           const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                           const int level_number,
                                           const double /*init_data_time*/,
                                           const bool can_be_refined,
                                           const bool /*initial_time*/,
                                           LDataManager* const /*l_data_manager*/)
{
    // Determine the extents of the physical domain.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    const double* const grid_x_lower = grid_geom->getXLower();
    const double* const grid_x_upper = grid_geom->getXUpper();

    // Loop over all patches in the specified level of the patch level and
    // initialize the local vertices.
    boost::multi_array_ref<double, 2>& X_array = *X_data->getLocalFormVecArray();
    boost::multi_array_ref<double, 2>& U_array = *U_data->getLocalFormVecArray();
    int local_idx = -1;
    int local_node_count = 0;
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    const IntVector<NDIM>& ratio = level->getRatio();
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();

        Pointer<LNodeSetData> index_data = patch->getPatchData(lag_node_index_idx);

        // Initialize the vertices whose initial locations will be within the
        // given patch.
        std::vector<std::pair<int, int> > patch_vertices;
        getPatchVertices(patch_vertices, patch, level_number, can_be_refined);
        local_node_count += patch_vertices.size();
        for (std::vector<std::pair<int, int> >::const_iterator it = patch_vertices.begin(); it != patch_vertices.end();
             ++it)
        {
            const std::pair<int, int>& point_idx = (*it);
            const int lagrangian_idx = getCanonicalLagrangianIndex(point_idx, level_number) + global_index_offset;
            const int local_petsc_idx = ++local_idx + local_index_offset;
            const int global_petsc_idx = local_petsc_idx + global_index_offset;

            // Get the coordinates of the present vertex.
            const libMesh::Point& X = getVertexPosn(point_idx, level_number);
            const CellIndex<NDIM> idx = IndexUtilities::getCellIndex(&X(0), grid_geom, ratio);
            for (int d = 0; d < NDIM; ++d)
            {
                X_array[local_petsc_idx][d] = X(d);
                if (X(d) <= grid_x_lower[d])
                {
                    TBOX_ERROR(d_object_name << "::initializeDataOnPatchLevel():\n"
                                             << "  encountered node below lower physical boundary\n"
                                             << "  please ensure that all nodes are within the "
                                                "computational domain."
                                             << std::endl);
                }
                if (X(d) >= grid_x_upper[d])
                {
                    TBOX_ERROR(d_object_name << "::initializeDataOnPatchLevel():\n"
                                             << "  encountered node above upper physical boundary\n"
                                             << "  please ensure that all nodes are within the "
                                                "computational domain."
                                             << std::endl);
                }
            }

            // Create or retrieve a pointer to the LNodeSet associated with the
            // current Cartesian grid cell.
            if (!index_data->isElement(idx))
            {
                index_data->appendItemPointer(idx, new LNodeSet());
            }
            LNodeSet* const node_set = index_data->getItem(idx);
            static const IntVector<NDIM> periodic_offset(0);
            static const IBTK::Point periodic_displacement(IBTK::Point::Zero());
            Pointer<MaterialPointSpec> point_spec =
                new MaterialPointSpec(lagrangian_idx,
                                      d_vertex_wgt[level_number][point_idx.first][point_idx.second],
                                      d_vertex_subdomain_id[level_number][point_idx.first][point_idx.second]);
            std::vector<Pointer<Streamable> > node_data(1, point_spec);
            node_set->push_back(new LNode(
                lagrangian_idx, global_petsc_idx, local_petsc_idx, periodic_offset, periodic_displacement, node_data));

            // Initialize the velocity of the present vertex.
            std::fill(&U_array[local_petsc_idx][0], &U_array[local_petsc_idx][0] + NDIM, 0.0);
        }
    }
    X_data->restoreArrays();
    U_data->restoreArrays();

    d_level_is_initialized[level_number] = true;

    // If a Lagrangian Silo data writer is registered with the initializer,
    // setup the visualization data corresponding to the present level of the
    // locally refined grid.
    if (d_silo_writer) initializeLSiloDataWriter(level_number);
    return local_node_count;
} // initializeDataOnPatchLevel

void
IMPInitializer::tagCellsForInitialRefinement(const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                             const int level_number,
                                             const double /*error_data_time*/,
                                             const int tag_index)
{
    // Loop over all patches in the specified level of the patch level and tag
    // cells for refinement wherever there are vertices assigned to a finer
    // level of the Cartesian grid.
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    const Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
    const IntVector<NDIM>& ratio = level->getRatio();
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const Box<NDIM>& patch_box = patch->getBox();

        Pointer<CellData<NDIM, int> > tag_data = patch->getPatchData(tag_index);

        // Tag cells for refinement whenever there are vertices whose initial
        // locations will be within the index space of the given patch, but on
        // the finer levels of the AMR patch hierarchy.
        const int max_levels = d_gridding_alg->getMaxLevels();
        const bool can_be_refined = level_can_be_refined(level_number, max_levels);
        for (int ln = level_number + 1; ln < max_levels; ++ln)
        {
            std::vector<std::pair<int, int> > patch_vertices;
            getPatchVertices(patch_vertices, patch, ln, can_be_refined);
            for (std::vector<std::pair<int, int> >::const_iterator it = patch_vertices.begin();
                 it != patch_vertices.end();
                 ++it)
            {
                const std::pair<int, int>& point_idx = (*it);
                const libMesh::Point& X = getVertexPosn(point_idx, ln);
                const CellIndex<NDIM> i = IndexUtilities::getCellIndex(&X(0), grid_geom, ratio);
                if (patch_box.contains(i)) (*tag_data)(i) = 1;
            }
        }
    }
    return;
} // tagCellsForInitialRefinement

void
IMPInitializer::writeVertexFile(std::string filename, int mesh_no, int level_number)
{
    const int max_levels = d_gridding_alg->getMaxLevels();
    if (level_number < 0) level_number = max_levels - 1;
    level_number = std::min(level_number, max_levels - 1);

    const int mpi_rank = SAMRAI_MPI::getRank();
    if (mpi_rank == 0)
    {
        if (filename.find(".vertex") == std::string::npos)
        {
            filename += ".vertex";
        }
        std::ofstream vertex_file(filename.c_str(), std::fstream::out);

        const int num_vertices = d_num_vertex[level_number][mesh_no];
        vertex_file << num_vertices << "\n";
        for (int k = 0; k < num_vertices; ++k)
        {
            const libMesh::Point& X = d_vertex_posn[level_number][mesh_no][k];
            for (int d = 0; d < NDIM; ++d)
            {
                vertex_file << X(d) << "\t";
            }
            vertex_file << "\n";
        }

        vertex_file.close();
    }

    return;
} // writeVertexFile

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IMPInitializer::initializeLSiloDataWriter(const int level_number)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(level_number >= 0);
    TBOX_ASSERT(level_number < d_gridding_alg->getMaxLevels());
    TBOX_ASSERT(d_level_is_initialized[level_number]);
#endif
    // WARNING: For now, we just register the visualization data on MPI process
    // 0.  This will fail if the structure is too large to be stored in the
    // memory available to a single MPI process.
    if (SAMRAI_MPI::getRank() == 0)
    {
        for (unsigned int j = 0; j < d_num_vertex[level_number].size(); ++j)
        {
            if (d_num_vertex[level_number][j])
            {
                std::ostringstream name_stream;
                name_stream << "mesh_" << j;
                const std::string vertices_name = name_stream.str() + "_vertices";
                d_silo_writer->registerMarkerCloud(
                    vertices_name, d_num_vertex[level_number][j], d_vertex_offset[level_number][j], level_number);
            }
        }
    }
    return;
} // initializeLSiloDataWriter

void
IMPInitializer::getPatchVertices(std::vector<std::pair<int, int> >& patch_vertices,
                                 const Pointer<Patch<NDIM> > patch,
                                 const int level_number,
                                 const bool /*can_be_refined*/) const
{
    // Loop over all of the vertices to determine the indices of those vertices
    // within the present patch.
    //
    // NOTE: This is clearly not the best way to do this, but it will work for
    // now.
    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
    const double* const patch_x_lower = patch_geom->getXLower();
    const double* const patch_x_upper = patch_geom->getXUpper();

    // Count the number of patch vertices.
    unsigned int num_patch_vertices = 0;
    for (unsigned int j = 0; j < d_num_vertex[level_number].size(); ++j)
    {
        for (int k = 0; k < d_num_vertex[level_number][j]; ++k)
        {
            const libMesh::Point& X = d_vertex_posn[level_number][j][k];
            bool patch_owns_node = true;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                patch_owns_node = patch_owns_node && (patch_x_lower[d] <= X(d)) && (X(d) < patch_x_upper[d]);
            }
            if (patch_owns_node) ++num_patch_vertices;
        }
    }

    // Quick return if there are no local patch vertices.
    if (!num_patch_vertices) return;

    // Get the index data for the local patch vertices.
    patch_vertices.reserve(patch_vertices.size() + num_patch_vertices);
    for (unsigned int j = 0; j < d_num_vertex[level_number].size(); ++j)
    {
        for (int k = 0; k < d_num_vertex[level_number][j]; ++k)
        {
            const libMesh::Point& X = d_vertex_posn[level_number][j][k];
            bool patch_owns_node = true;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                patch_owns_node = patch_owns_node && (patch_x_lower[d] <= X(d)) && (X(d) < patch_x_upper[d]);
            }
            if (patch_owns_node) patch_vertices.push_back(std::make_pair(j, k));
        }
    }
    return;
} // getPatchVertices

int
IMPInitializer::getCanonicalLagrangianIndex(const std::pair<int, int>& point_index, const int level_number) const
{
    return d_vertex_offset[level_number][point_index.first] + point_index.second;
} // getCanonicalLagrangianIndex

const libMesh::Point&
IMPInitializer::getVertexPosn(const std::pair<int, int>& point_index, const int level_number) const
{
    return d_vertex_posn[level_number][point_index.first][point_index.second];
} // getVertexPosn

void IMPInitializer::getFromInput(Pointer<Database> /*db*/)
{
    return;
} // getFromInput

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
