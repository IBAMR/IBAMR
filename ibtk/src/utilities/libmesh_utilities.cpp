// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2019 by the IBAMR developers
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

#include "ibtk/FECache.h"
#include "ibtk/QuadratureCache.h"
#include "ibtk/libmesh_utilities.h"

#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/mesh_base.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parallel.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/point.h"
#include "libmesh/system.h"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>
getQuadratureKey(const libMesh::QuadratureType quad_type,
                 libMesh::Order order,
                 const bool use_adaptive_quadrature,
                 const double point_density,
                 const libMesh::Elem* const elem,
                 const boost::multi_array<double, 2>& X_node,
                 const double dx_min)
{
    const libMesh::ElemType elem_type = elem->type();
#ifndef NDEBUG
    TBOX_ASSERT(elem->p_level() == 0); // higher levels are not implemented
#endif
    if (use_adaptive_quadrature)
    {
        const double hmax = get_max_edge_length(elem, X_node);
        int npts = int(std::ceil(point_density * hmax / dx_min));
        if (npts < 3)
        {
            if (elem->default_order() == libMesh::FIRST)
                npts = 2;
            else
                npts = 3;
        }
        switch (quad_type)
        {
        case libMesh::QGAUSS:
            order = static_cast<libMesh::Order>(std::min(2 * npts - 1, static_cast<int>(libMesh::FORTYTHIRD)));
            break;
        case libMesh::QGRID:
            order = static_cast<libMesh::Order>(npts);
            break;
        default:
            TBOX_ERROR("IBTK::getQuadratureKey():\n"
                       << "  adaptive quadrature rules are available only for quad_type = QGAUSS "
                          "or QGRID\n");
        }
    }

    return std::make_tuple(elem_type, quad_type, order);
}

void
write_elem_partitioning(const std::string& file_name, const libMesh::System& position_system)
{
    const int current_rank = position_system.comm().rank();
    const unsigned int position_system_n = position_system.number();
    const libMesh::NumericVector<double>& local_position = *position_system.solution.get();
    const libMesh::MeshBase& mesh = position_system.get_mesh();
    const unsigned int spacedim = mesh.spatial_dimension();
    // TODO: there is something wrong with the way we set up the ghost data in
    // the position vectors: not all locally owned nodes are, in fact, locally
    // available. Get around this by localizing first. Since all processes
    // write to the same file this isn't the worst bottleneck in this
    // function, anyway.
    std::vector<double> position(local_position.size());
    local_position.localize(position);
    std::stringstream current_processor_output;

    const auto end_elem = mesh.local_elements_end();
    for (auto elem = mesh.local_elements_begin(); elem != end_elem; ++elem)
    {
        const unsigned int n_nodes = (*elem)->n_nodes();
        libMesh::Point center;
        // TODO: this is a bit crude: if we use isoparametric elements (e.g.,
        // Tri6) then this is not a very accurate representation of the center
        // of the element. We should replace this with something more accurate.
        for (unsigned int node_n = 0; node_n < n_nodes; ++node_n)
        {
            const libMesh::Node& node = (*elem)->node_ref(node_n);
            TBOX_ASSERT(node.n_vars(position_system_n) == spacedim);
            for (unsigned int d = 0; d < spacedim; ++d)
            {
                center(d) += position[node.dof_number(position_system_n, d, 0)];
            }
        }
        center *= 1.0 / n_nodes;

        for (unsigned int d = 0; d < spacedim; ++d)
        {
            current_processor_output << center(d) << ',';
        }
        if (spacedim == 2)
        {
            current_processor_output << 0.0 << ',';
        }
        current_processor_output << current_rank << '\n';
    }

    // clear the file before we append to it
    if (current_rank == 0)
    {
        std::remove(file_name.c_str());
    }
    const int n_processes = position_system.comm().size();
    for (int rank = 0; rank < n_processes; ++rank)
    {
        if (rank == current_rank)
        {
            std::ofstream out(file_name, std::ios_base::app);
            if (rank == 0)
            {
                out << "x,y,z,r\n";
            }
            out << current_processor_output.rdbuf();
        }
        position_system.comm().barrier();
    }
}

void
write_node_partitioning(const std::string& file_name, const libMesh::System& position_system)
{
    const int current_rank = position_system.comm().rank();
    const unsigned int position_system_n = position_system.number();
    const libMesh::NumericVector<double>& local_position = *position_system.solution.get();
    const libMesh::MeshBase& mesh = position_system.get_mesh();
    const unsigned int spacedim = mesh.spatial_dimension();

    // TODO: there is something wrong with the way we set up the ghost data in
    // the position vectors: not all locally owned nodes are, in fact,
    // locally available. Get around this by localizing first. Since all
    // processes write to the same file this isn't the worst bottleneck in
    // this function, anyway.
    std::vector<double> position(local_position.size());
    local_position.localize(position);
    std::stringstream current_processor_output;

    const auto end_node = mesh.local_nodes_end();
    for (auto node_it = mesh.local_nodes_begin(); node_it != end_node; ++node_it)
    {
        const libMesh::Node* const node = *node_it;
        if (node->n_vars(position_system_n))
        {
            TBOX_ASSERT(node->n_vars(position_system_n) == spacedim);
            for (unsigned int d = 0; d < spacedim; ++d)
            {
                current_processor_output << position[node->dof_number(position_system_n, d, 0)] << ',';
            }
            if (spacedim == 2)
            {
                current_processor_output << 0.0 << ',';
            }
            current_processor_output << current_rank << '\n';
        }
    }

    // clear the file before we append to it
    if (current_rank == 0)
    {
        std::remove(file_name.c_str());
    }
    const int n_processes = position_system.comm().size();
    for (int rank = 0; rank < n_processes; ++rank)
    {
        if (rank == current_rank)
        {
            std::ofstream out(file_name, std::ios_base::app);
            if (rank == 0)
            {
                out << "x,y,z,r\n";
            }
            out << current_processor_output.rdbuf();
        }
        position_system.comm().barrier();
    }
}

std::vector<libMeshWrappers::BoundingBox>
get_local_active_element_bounding_boxes(const libMesh::MeshBase& mesh,
                                        const libMesh::System& X_system,
                                        const libMesh::QuadratureType quad_type,
                                        const libMesh::Order quad_order,
                                        const bool use_adaptive_quadrature,
                                        const double point_density,
                                        const double patch_dx_min)
{
    const unsigned int dim = mesh.mesh_dimension();
    TBOX_ASSERT(dim <= LIBMESH_DIM);
    const unsigned int X_sys_num = X_system.number();
    auto X_ghost_vec_ptr = X_system.current_local_solution->zero_clone();
    auto& X_ghost_vec = dynamic_cast<libMesh::PetscVector<double>&>(*X_ghost_vec_ptr);
    X_ghost_vec = *X_system.solution;

    std::vector<libMeshWrappers::BoundingBox> bboxes;

    std::vector<std::vector<libMesh::dof_id_type> > dof_indices(dim);
    boost::multi_array<double, 2> X_node;
    QuadratureCache quad_cache(dim);
    FECache fe_cache(dim, X_system.get_dof_map().variable_type(0), update_phi);
    using quad_key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;
    const auto el_begin = mesh.active_local_elements_begin();
    const auto el_end = mesh.active_local_elements_end();
    for (auto el_it = el_begin; el_it != el_end; ++el_it)
    {
        // 0. Set up bounding box
        bboxes.emplace_back();
        auto& box = bboxes.back();
        libMesh::Point& lower_bound = box.first;
        libMesh::Point& upper_bound = box.second;
        for (unsigned int d = 0; d < LIBMESH_DIM; ++d)
        {
            lower_bound(d) = std::numeric_limits<double>::max();
            upper_bound(d) = -std::numeric_limits<double>::max();
        }

        // 1. extract node locations
        const libMesh::Elem* const elem = *el_it;
        const unsigned int n_nodes = elem->n_nodes();
        for (unsigned int d = 0; d < dim; ++d) dof_indices[d].clear();

        for (unsigned int k = 0; k < n_nodes; ++k)
        {
            const libMesh::Node* const node = elem->node_ptr(k);
            for (unsigned int d = 0; d < dim; ++d)
            {
                TBOX_ASSERT(node->n_dofs(X_sys_num, d) == 1);
                dof_indices[d].push_back(node->dof_number(X_sys_num, d, 0));
            }
        }
        get_values_for_interpolation(X_node, X_ghost_vec, dof_indices);

        // 2. compute mapped quadrature points in the deformed configuration
        const quad_key_type key =
            getQuadratureKey(quad_type, quad_order, use_adaptive_quadrature, point_density, elem, X_node, patch_dx_min);
        libMesh::QBase& quadrature = quad_cache[key];
        libMesh::FEBase& fe = fe_cache(key, elem);
        const std::vector<std::vector<double> >& phi_X = fe.get_phi();
        const std::vector<libMesh::Point>& q_points = quadrature.get_points();
        for (unsigned int qp = 0; qp < q_points.size(); ++qp)
        {
            libMesh::Point mapped_point;
            for (unsigned int k = 0; k < n_nodes; ++k)
            {
                for (unsigned int d = 0; d < dim; ++d)
                {
                    mapped_point(d) += phi_X[k][qp] * X_node[k][d];
                }
            }
            for (unsigned int d = 0; d < dim; ++d)
            {
                lower_bound(d) = std::min(lower_bound(d), mapped_point(d));
                upper_bound(d) = std::max(upper_bound(d), mapped_point(d));
            }

            // fill extra dimension with 0.0, which is libMesh's convention
            for (unsigned int d = dim; d < LIBMESH_DIM; ++d)
            {
                lower_bound(d) = 0.0;
                upper_bound(d) = 0.0;
            }
        }
    }

    return bboxes;
} // get_local_active_element_bounding_boxes

std::vector<libMeshWrappers::BoundingBox>
get_local_active_element_bounding_boxes(const libMesh::MeshBase& mesh, const libMesh::System& X_system)
{
    const unsigned int dim = mesh.mesh_dimension();
    const unsigned int X_sys_num = X_system.number();
    auto X_ghost_vec_ptr = X_system.current_local_solution->zero_clone();
    auto& X_ghost_vec = dynamic_cast<libMesh::PetscVector<double>&>(*X_ghost_vec_ptr);
    X_ghost_vec = *X_system.solution;

    std::vector<libMeshWrappers::BoundingBox> bboxes;
    bboxes.reserve(mesh.n_local_elem());

    std::vector<libMesh::dof_id_type> dof_indices;
    std::vector<double> X_node;
    const auto el_begin = mesh.active_local_elements_begin();
    const auto el_end = mesh.active_local_elements_end();
    for (auto el_it = el_begin; el_it != el_end; ++el_it)
    {
        const libMesh::Elem* const elem = *el_it;
        const unsigned int n_nodes = elem->n_nodes();
        bboxes.emplace_back();
        auto& box = bboxes.back();
        libMesh::Point& lower_bound = box.first;
        libMesh::Point& upper_bound = box.second;
        for (unsigned int d = 0; d < LIBMESH_DIM; ++d)
        {
            lower_bound(d) = std::numeric_limits<double>::max();
            upper_bound(d) = -std::numeric_limits<double>::max();
        }

        dof_indices.clear();
        for (unsigned int k = 0; k < n_nodes; ++k)
        {
            const libMesh::Node* const node = elem->node_ptr(k);
            for (unsigned int d = 0; d < dim; ++d)
            {
                TBOX_ASSERT(node->n_dofs(X_sys_num, d) == 1);
                dof_indices.push_back(node->dof_number(X_sys_num, d, 0));
            }
        }

        X_node.resize(dof_indices.size());
        X_ghost_vec.get(dof_indices, X_node.data());
        for (unsigned int k = 0; k < n_nodes; ++k)
        {
            for (unsigned int d = 0; d < dim; ++d)
            {
                const double& X = X_node[k * dim + d];
                lower_bound(d) = std::min(lower_bound(d), X);
                upper_bound(d) = std::max(upper_bound(d), X);
            }
        }

        // fill extra dimension with 0.0, which is libMesh's convention
        for (unsigned int d = dim; d < LIBMESH_DIM; ++d)
        {
            lower_bound(d) = 0.0;
            upper_bound(d) = 0.0;
        }
    }

    return bboxes;
} // get_local_active_element_bounding_boxes

std::vector<libMeshWrappers::BoundingBox>
get_global_active_element_bounding_boxes(const libMesh::MeshBase& mesh,
                                         const std::vector<libMeshWrappers::BoundingBox>& local_bboxes)
{
    const std::size_t n_elem = mesh.n_active_elem();
    std::vector<double> flattened_bboxes(2 * LIBMESH_DIM * n_elem);
    std::size_t elem_n = 0;
    const auto el_begin = mesh.active_local_elements_begin();
    const auto el_end = mesh.active_local_elements_end();
    for (auto el_it = el_begin; el_it != el_end; ++el_it)
    {
        const auto id = (*el_it)->id();
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            flattened_bboxes[2 * id * NDIM + d] = local_bboxes[elem_n].first(d);
            flattened_bboxes[(2 * id + 1) * NDIM + d] = local_bboxes[elem_n].second(d);
        }
        ++elem_n;
    }
    const int ierr = MPI_Allreduce(
        MPI_IN_PLACE, flattened_bboxes.data(), flattened_bboxes.size(), MPI_DOUBLE, MPI_SUM, mesh.comm().get());
    TBOX_ASSERT(ierr == 0);

    std::vector<libMeshWrappers::BoundingBox> global_bboxes(n_elem);
    for (unsigned int e = 0; e < n_elem; ++e)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            global_bboxes[e].first(d) = flattened_bboxes[2 * e * NDIM + d];
            global_bboxes[e].second(d) = flattened_bboxes[(2 * e + 1) * NDIM + d];
        }
    }
    return global_bboxes;
} // get_local_active_element_bounding_boxes

std::vector<libMeshWrappers::BoundingBox>
get_global_active_element_bounding_boxes(const libMesh::MeshBase& mesh, const libMesh::System& X_system)
{
    static_assert(NDIM <= LIBMESH_DIM,
                  "NDIM should be no more than LIBMESH_DIM for this function to "
                  "work correctly.");
    return get_global_active_element_bounding_boxes(mesh, get_local_active_element_bounding_boxes(mesh, X_system));
} // get_global_active_element_bounding_boxes
//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
