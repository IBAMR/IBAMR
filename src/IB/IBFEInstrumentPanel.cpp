// Filename: IBFEInstrumentPanel.cpp
// Created on 12 April 2018 by Charles Puelz
//
// Copyright (c) 2002-2018, Charles Puelz and Boyce Griffith
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

#include <array>
#include <algorithm>
#include <fstream>
#include <limits>
#include <map>
#include <cmath>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "BasePatchLevel.h"
#include "Box.h"
#include "BoxArray.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "Eigen/Geometry" // IWYU pragma: keep
#include "IBAMR_config.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "SideData.h"
#include "SideIndex.h"
#include "boost/multi_array.hpp"
#include "ibamr/IBFEInstrumentPanel.h"
#include "ibamr/IBFEMethod.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/FEDataManager.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LMesh.h"
#include "ibtk/LNode.h"
#include "ibtk/ibtk_utilities.h"
#include "libmesh/boundary_info.h"
#include "libmesh/dense_vector.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/equation_systems.h"
#include "libmesh/face_tri3.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_function.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/point.h"
#include "libmesh/quadrature_grid.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/string_to_enum.h"
#include "petscvec.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

using namespace libMesh;

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
double
linear_interp(const Vector& X,
              const Index<NDIM>& i_cell,
              const Vector& X_cell,
              const CellData<NDIM, double>& v,
              const Index<NDIM>& /*patch_lower*/,
              const Index<NDIM>& /*patch_upper*/,
              const double* const /*x_lower*/,
              const double* const /*x_upper*/,
              const double* const dx)
{
    std::array<bool, NDIM> is_lower;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        is_lower[d] = X[d] < X_cell[d];
    }
    double U = 0.0;
#if (NDIM == 3)
    for (int i_shift2 = (is_lower[2] ? -1 : 0); i_shift2 <= (is_lower[2] ? 0 : 1); ++i_shift2)
    {
#endif
        for (int i_shift1 = (is_lower[1] ? -1 : 0); i_shift1 <= (is_lower[1] ? 0 : 1); ++i_shift1)
        {
            for (int i_shift0 = (is_lower[0] ? -1 : 0); i_shift0 <= (is_lower[0] ? 0 : 1); ++i_shift0)
            {
                const Vector X_center(X_cell[0] + static_cast<double>(i_shift0) * dx[0],
                                      X_cell[1] + static_cast<double>(i_shift1) * dx[1]
#if (NDIM == 3)
                                      ,
                                      X_cell[2] + static_cast<double>(i_shift2) * dx[2]
#endif
                );
                const double wgt =
                    (((X[0] < X_center[0] ? X[0] - (X_center[0] - dx[0]) : (X_center[0] + dx[0]) - X[0]) / dx[0]) *
                     ((X[1] < X_center[1] ? X[1] - (X_center[1] - dx[1]) : (X_center[1] + dx[1]) - X[1]) / dx[1])
#if (NDIM == 3)
                     * ((X[2] < X_center[2] ? X[2] - (X_center[2] - dx[2]) : (X_center[2] + dx[2]) - X[2]) / dx[2])
#endif
                    );
                const Index<NDIM> i(i_shift0 + i_cell(0),
                                    i_shift1 + i_cell(1)
#if (NDIM == 3)
                                        ,
                                    i_shift2 + i_cell(2)
#endif
                );
                const CellIndex<NDIM> i_c(i);
                U += v(i_c) * wgt;
            }
        }
#if (NDIM == 3)
    }
#endif
    return U;
}

template <int N>
Eigen::Matrix<double, N, 1>
linear_interp(const Vector& X,
              const Index<NDIM>& i_cell,
              const Vector& X_cell,
              const CellData<NDIM, double>& v,
              const Index<NDIM>& /*patch_lower*/,
              const Index<NDIM>& /*patch_upper*/,
              const double* const /*x_lower*/,
              const double* const /*x_upper*/,
              const double* const dx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(v.getDepth() == N);
#endif
    std::array<bool, NDIM> is_lower;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        is_lower[d] = X[d] < X_cell[d];
    }
    Eigen::Matrix<double, N, 1> U(Eigen::Matrix<double, N, 1>::Zero());
#if (NDIM == 3)
    for (int i_shift2 = (is_lower[2] ? -1 : 0); i_shift2 <= (is_lower[2] ? 0 : 1); ++i_shift2)
    {
#endif
        for (int i_shift1 = (is_lower[1] ? -1 : 0); i_shift1 <= (is_lower[1] ? 0 : 1); ++i_shift1)
        {
            for (int i_shift0 = (is_lower[0] ? -1 : 0); i_shift0 <= (is_lower[0] ? 0 : 1); ++i_shift0)
            {
                const Vector X_center(X_cell[0] + static_cast<double>(i_shift0) * dx[0],
                                      X_cell[1] + static_cast<double>(i_shift1) * dx[1]
#if (NDIM == 3)
                                      ,
                                      X_cell[2] + static_cast<double>(i_shift2) * dx[2]
#endif
                );
                const double wgt =
                    (((X[0] < X_center[0] ? X[0] - (X_center[0] - dx[0]) : (X_center[0] + dx[0]) - X[0]) / dx[0]) *
                     ((X[1] < X_center[1] ? X[1] - (X_center[1] - dx[1]) : (X_center[1] + dx[1]) - X[1]) / dx[1])
#if (NDIM == 3)
                     * ((X[2] < X_center[2] ? X[2] - (X_center[2] - dx[2]) : (X_center[2] + dx[2]) - X[2]) / dx[2])
#endif
                    );
                const Index<NDIM> i(i_shift0 + i_cell(0),
                                    i_shift1 + i_cell(1)
#if (NDIM == 3)
                                        ,
                                    i_shift2 + i_cell(2)
#endif
                );
                const CellIndex<NDIM> i_c(i);
                for (int k = 0; k < N; ++k)
                {
                    U[k] += v(i_c, k) * wgt;
                }
            }
        }
#if (NDIM == 3)
    }
#endif
    return U;
}

Vector
linear_interp(const Vector& X,
              const Index<NDIM>& i_cell,
              const Vector& X_cell,
              const SideData<NDIM, double>& v,
              const Index<NDIM>& /*patch_lower*/,
              const Index<NDIM>& /*patch_upper*/,
              const double* const /*x_lower*/,
              const double* const /*x_upper*/,
              const double* const dx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(v.getDepth() == 1);
#endif
    Vector U(Vector::Zero());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        std::array<bool, NDIM> is_lower;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (d == axis)
            {
                is_lower[d] = false;
            }
            else
            {
                is_lower[d] = X[d] < X_cell[d];
            }
        }
#if (NDIM == 3)
        for (int i_shift2 = (is_lower[2] ? -1 : 0); i_shift2 <= (is_lower[2] ? 0 : 1); ++i_shift2)
        {
#endif
            for (int i_shift1 = (is_lower[1] ? -1 : 0); i_shift1 <= (is_lower[1] ? 0 : 1); ++i_shift1)
            {
                for (int i_shift0 = (is_lower[0] ? -1 : 0); i_shift0 <= (is_lower[0] ? 0 : 1); ++i_shift0)
                {
                    const Vector X_side(X_cell[0] + (static_cast<double>(i_shift0) + (axis == 0 ? -0.5 : 0.0)) * dx[0],
                                        X_cell[1] + (static_cast<double>(i_shift1) + (axis == 1 ? -0.5 : 0.0)) * dx[1]
#if (NDIM == 3)
                                        ,
                                        X_cell[2] + (static_cast<double>(i_shift2) + (axis == 2 ? -0.5 : 0.0)) * dx[2]
#endif
                    );
                    const double wgt =
                        (((X[0] < X_side[0] ? X[0] - (X_side[0] - dx[0]) : (X_side[0] + dx[0]) - X[0]) / dx[0]) *
                         ((X[1] < X_side[1] ? X[1] - (X_side[1] - dx[1]) : (X_side[1] + dx[1]) - X[1]) / dx[1])
#if (NDIM == 3)
                         * ((X[2] < X_side[2] ? X[2] - (X_side[2] - dx[2]) : (X_side[2] + dx[2]) - X[2]) / dx[2])
#endif
                        );
                    const Index<NDIM> i(i_shift0 + i_cell(0),
                                        i_shift1 + i_cell(1)
#if (NDIM == 3)
                                            ,
                                        i_shift2 + i_cell(2)
#endif
                    );
                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                    U[axis] += v(i_s) * wgt;
                }
            }
#if (NDIM == 3)
        }
#endif
    }
    return U;
}
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBFEInstrumentPanel::IBFEInstrumentPanel(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db, const int part)
    : d_part(part), d_plot_directory_name(NDIM == 2 ? "viz_inst2d" : "viz_inst3d")
{
    // get input data
    IBFEInstrumentPanel::getFromInput(input_db);

    // make plot directory
    Utilities::recursiveMkdir(d_plot_directory_name);

    // set up file streams
    if (SAMRAI_MPI::getRank() == 0)
    {
        d_mean_pressure_stream.open(d_plot_directory_name + "/mean_pressure.dat");
        d_flux_stream.open(d_plot_directory_name + "/flux.dat");
        d_mean_pressure_stream.precision(15);
        d_flux_stream.precision(15);
    }
}


void
IBFEInstrumentPanel::initializeHierarchyIndependentData(IBFEMethod* ib_method_ops)
{
    // get relevant things for corresponding part
    const FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager(d_part);
    const EquationSystems* equation_systems = fe_data_manager->getEquationSystems();
    const MeshBase* mesh = &equation_systems->get_mesh();
    const BoundaryInfo& boundary_info = *mesh->boundary_info;
    const libMesh::Parallel::Communicator& comm_in = mesh->comm();

    // get equation systems from the mesh we will need.
    const System& dX_system = equation_systems->get_system(IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
    const unsigned int dX_sys_num = dX_system.number();
    const System& U_system = equation_systems->get_system(IBFEMethod::VELOCITY_SYSTEM_NAME);
    const unsigned int U_sys_num = U_system.number();

    // some local variables
    std::vector<dof_id_type> nodes;
    std::vector<boundary_id_type> bcs;
    std::vector<std::vector<dof_id_type> > temp_node_dof_IDs;
    std::vector<std::set<dof_id_type> > temp_node_dof_ID_sets;
    std::vector<std::vector<libMesh::Point> > temp_nodes;
    std::vector<libMesh::Point> meter_centroids;
    // new API in 1.4.0
#if 1 <= LIBMESH_MAJOR_VERSION && 4 <= LIBMESH_MINOR_VERSION
    const std::vector<std::tuple<dof_id_type, boundary_id_type> > node_list = boundary_info.build_node_list();
    for (const std::tuple<dof_id_type, boundary_id_type>& pair : node_list)
    {
        nodes.push_back(std::get<0>(pair));
        bcs.push_back(std::get<1>(pair));
    }
#else
    boundary_info.build_node_list(nodes, bcs);
#endif

    // check to make sure there are node sets to work with
    if (nodes.size() == 0 || bcs.size() == 0 || (nodes.size() != bcs.size()))
    {
        TBOX_ERROR("IBFEInstrumentPanel::initializeHierarchyIndependentData : "
                   << "nodesets not set up correctly or don't exist in FE part with number " << d_part);
    }

    // resize members and local variables
    d_num_meters = d_nodeset_IDs_for_meters.size();
    d_meter_radii.resize(d_num_meters);
    d_quad_order.resize(d_num_meters);
    d_U_dof_idx.resize(d_num_meters);
    d_dX_dof_idx.resize(d_num_meters);
    d_node_dof_IDs.resize(d_num_meters);
    d_nodes.resize(d_num_meters);
    d_num_nodes.resize(d_num_meters);
    d_mean_pressure_values.resize(d_num_meters);
    d_flow_values.resize(d_num_meters);
    temp_node_dof_IDs.resize(d_num_meters);
    temp_node_dof_ID_sets.resize(d_num_meters);
    temp_nodes.resize(d_num_meters);
    meter_centroids.resize(d_num_meters);

    // populate temp vectors
    for (unsigned int ii = 0; ii < nodes.size(); ++ii)
    {
        for (int jj = 0; jj < d_nodeset_IDs_for_meters.size(); ++jj)
        {
            if (d_nodeset_IDs_for_meters[jj] == bcs[ii])
            {
                temp_node_dof_ID_sets[jj].insert(nodes[ii]);
            }
        }
    }

    // make sure nodes and node dof IDs have the same order
    // on all the processes.
    for (int jj = 0; jj < d_nodeset_IDs_for_meters.size(); ++jj)
    {
        for (const auto& node_id : temp_node_dof_ID_sets[jj])
        {
            temp_node_dof_IDs[jj].push_back(node_id);
            const Node* node = &mesh->node_ref(node_id);
            temp_nodes[jj].push_back(*node);
            meter_centroids[jj] += *node;
        }
    }

    // loop over meters and sort the nodes
    for (unsigned int jj = 0; jj < d_num_meters; ++jj)
    {
        // finish computing centroid
        meter_centroids[jj] /= static_cast<double>(temp_nodes[jj].size());
        // make sure nodes are sorted with a particular orientation
        libMesh::Point temp_point;
        dof_id_type temp_dof_id;
        double max_dist = std::numeric_limits<double>::max();
        d_nodes[jj].push_back(temp_nodes[jj][0]);
        d_node_dof_IDs[jj].push_back(temp_node_dof_IDs[jj][0]);
        max_dist = std::numeric_limits<double>::max();
        for (unsigned int kk = 1; kk < temp_nodes[jj].size(); ++kk)
        {
            for (unsigned int ll = 0; ll < temp_nodes[jj].size(); ++ll)
            {
                // here we find the closest node to the previous one
                // with some orientation
                libMesh::Point dist = temp_nodes[jj][ll] - d_nodes[jj][kk - 1];

                if (dist.norm() < max_dist)
                {
                    // make sure we haven't already added this node
                    bool added = false;
                    for (unsigned int ii = 1; ii < kk + 1; ++ii)
                    {
                        if (temp_nodes[jj][ll] == d_nodes[jj][ii - 1]) added = true;
                    }
                    if (!added)
                    {
                        temp_point = temp_nodes[jj][ll];
                        temp_dof_id = temp_node_dof_IDs[jj][ll];
                        max_dist = dist.norm();
                    }
                }
            }
            d_nodes[jj].push_back(temp_point);
            d_node_dof_IDs[jj].push_back(temp_dof_id);
            max_dist = std::numeric_limits<double>::max();
        }
    } // loop over meters

    // initialize meshes and number of nodes
    for (unsigned int jj = 0; jj < d_num_meters; ++jj)
    {
        d_meter_meshes.emplace_back(new SerialMesh(comm_in, NDIM));
        d_meter_mesh_names.push_back("meter_mesh_" + std::to_string(d_nodeset_IDs_for_meters[jj]));
        d_num_nodes[jj] = d_nodes[jj].size();
    }

    // build the meshes
    for (unsigned int ii = 0; ii < d_num_meters; ++ii)
    {
        d_meter_meshes[ii]->set_spatial_dimension(NDIM);
        d_meter_meshes[ii]->set_mesh_dimension(NDIM - 1);
        d_meter_meshes[ii]->reserve_nodes(d_num_nodes[ii] + 1);
        d_meter_meshes[ii]->reserve_elem(d_num_nodes[ii]);

        // add nodes
        for (unsigned int jj = 0; jj < d_num_nodes[ii]; ++jj)
        {
            d_meter_meshes[ii]->add_point(d_nodes[ii][jj], jj);
        }
        // add centroid
        d_meter_meshes[ii]->add_point(meter_centroids[ii], d_num_nodes[ii]);

        for (unsigned int jj = 0; jj < d_num_nodes[ii]; ++jj)
        {
            Elem* elem = new Tri3;
            elem->set_id(jj);
            // libMesh will delete elem
            elem = d_meter_meshes[ii]->add_elem(elem);
            elem->set_node(0) = d_meter_meshes[ii]->node_ptr(d_num_nodes[ii]);
            elem->set_node(1) = d_meter_meshes[ii]->node_ptr(jj);
            elem->set_node(2) = d_meter_meshes[ii]->node_ptr((jj + 1) % d_num_nodes[ii]);
        }
        d_meter_meshes[ii]->allow_renumbering(false);
        d_meter_meshes[ii]->prepare_for_use();
    } // loop over meters

    // initialize meter mesh equation systems, for both velocity and displacement
    for (unsigned int jj = 0; jj < d_num_meters; ++jj)
    {
        d_meter_systems.emplace_back(new EquationSystems(*d_meter_meshes[jj]));
        auto& velocity_sys = d_meter_systems[jj]->add_system<LinearImplicitSystem>(IBFEMethod::VELOCITY_SYSTEM_NAME);
        velocity_sys.add_variable("U_0", static_cast<Order>(1), LAGRANGE);
        velocity_sys.add_variable("U_1", static_cast<Order>(1), LAGRANGE);
        velocity_sys.add_variable("U_2", static_cast<Order>(1), LAGRANGE);
        auto& displacement_sys =
            d_meter_systems[jj]->add_system<LinearImplicitSystem>(IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
        displacement_sys.add_variable("dX_0", static_cast<Order>(1), LAGRANGE);
        displacement_sys.add_variable("dX_1", static_cast<Order>(1), LAGRANGE);
        displacement_sys.add_variable("dX_2", static_cast<Order>(1), LAGRANGE);

        // we use these serialized vectors to store the DOFs for the systems on all processes.
        velocity_sys.add_vector("serial solution", false, libMesh::SERIAL);
        displacement_sys.add_vector("serial solution", false, libMesh::SERIAL);

        d_meter_systems[jj]->init();
    }

    // store dof indices for the velocity and displacement systems that we will use later
    for (unsigned int jj = 0; jj < d_num_meters; ++jj)
    {
        for (unsigned int ii = 0; ii < d_num_nodes[jj]; ++ii)
        {
            const Node* node = &mesh->node_ref(d_node_dof_IDs[jj][ii]);
            std::vector<dof_id_type> dX_dof_index;
            std::vector<dof_id_type> U_dof_index;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dX_dof_index.push_back(node->dof_number(dX_sys_num, d, 0));
                U_dof_index.push_back(node->dof_number(U_sys_num, d, 0));
            }
            d_dX_dof_idx[jj].push_back(dX_dof_index);
            d_U_dof_idx[jj].push_back(U_dof_index);
        }
    }
    d_initialized = true;
}

void
IBFEInstrumentPanel::initializeHierarchyDependentData(IBFEMethod* ib_method_ops,
                                                      Pointer<PatchHierarchy<NDIM> > hierarchy)
{
    if (!d_initialized)
    {
        initializeHierarchyIndependentData(ib_method_ops);
    }
    if (d_num_meters == 0) return;

    // loop over meters, update system data, and get the maximum
    // radius for each meter in this FE part.
    // the radius of a meter is defined to be the largest distance
    // from the centroid to a node.
    for (unsigned int jj = 0; jj < d_num_meters; ++jj)
    {
        // update FE system data for meter_mesh
        initializeSystemDependentData(ib_method_ops, jj);
    }

    // get info about levels in AMR mesh
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    // Determine the finest grid spacing in the Cartesian grid hierarchy.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    const double* const domainXLower = grid_geom->getXLower();
    const double* const domainXUpper = grid_geom->getXUpper();
    const double* const dx_coarsest = grid_geom->getDx();
    TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
    const Box<NDIM> domain_box = grid_geom->getPhysicalDomain()[0];

    // get the finest spacing of fluid grid
    const IntVector<NDIM>& ratio_to_level_zero = hierarchy->getPatchLevel(finest_ln)->getRatio();
    std::array<double, NDIM> dx_finest;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        dx_finest[d] = dx_coarsest[d] / static_cast<double>(ratio_to_level_zero(d));
    }
    const double h_finest = *std::min_element(dx_finest.begin(), dx_finest.end());

    for (unsigned int jj = 0; jj < d_num_meters; ++jj)
    {
        // set the quadrature rule adaptively according to the fluid mesh size
        if (d_use_adaptive_quadrature)
            d_quad_order[jj] = static_cast<Order>(d_meter_radii[jj] / (0.25 * h_finest));
        else
            d_quad_order[jj] = d_input_quad_order;
        // print a warning
        if (d_quad_type == libMesh::QGRID && d_quad_order[jj] > libMesh::FORTYTHIRD)
        {
            TBOX_WARNING("IBFEInstrumentPanel::initializeHierarchyDependentData : "
                         << "QGrid quadrature order exceeds 43 for meter mesh in IBFE part " << d_part << "."
                         << " there may be undefined behavior in casting to this"
                         << " Order in older versions of libMesh.");
        }
    }

    // reset the quad point maps
    d_quad_point_map.clear();
    d_quad_point_map.resize(finest_ln + 1);

    // loop over levels and assign each quadrature point to one level
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        const IntVector<NDIM>& ratio = level->getRatio();
        const Box<NDIM> domain_box_level = Box<NDIM>::refine(domain_box, ratio);
        const Index<NDIM>& domain_box_level_lower = domain_box_level.lower();
        const Index<NDIM>& domain_box_level_upper = domain_box_level.upper();
        std::array<double, NDIM> dx;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dx[d] = dx_coarsest[d] / static_cast<double>(ratio(d));
        }

        Pointer<PatchLevel<NDIM> > finer_level =
            (ln < finest_ln ? hierarchy->getPatchLevel(ln + 1) : Pointer<BasePatchLevel<NDIM> >(nullptr));
        const IntVector<NDIM>& finer_ratio = (ln < finest_ln ? finer_level->getRatio() : IntVector<NDIM>(1));
        const Box<NDIM> finer_domain_box_level = Box<NDIM>::refine(domain_box, finer_ratio);
        const Index<NDIM>& finer_domain_box_level_lower = finer_domain_box_level.lower();
        const Index<NDIM>& finer_domain_box_level_upper = finer_domain_box_level.upper();
        std::array<double, NDIM> finer_dx;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            finer_dx[d] = dx_coarsest[d] / static_cast<double>(finer_ratio(d));
        }

        for (unsigned int jj = 0; jj < d_num_meters; ++jj)
        {
            const LinearImplicitSystem& displacement_sys =
                d_meter_systems[jj]->get_system<LinearImplicitSystem>(IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
            const NumericVector<double>& displacement_coords = displacement_sys.get_vector("serial solution");
            const DofMap& dof_map = displacement_sys.get_dof_map();
            FEType fe_type = displacement_sys.variable_type(0);

            // set up FE objects
            std::unique_ptr<FEBase> fe_elem(FEBase::build(NDIM - 1, fe_type));
            std::unique_ptr<QBase> qrule(QBase::build(d_quad_type, NDIM - 1, d_quad_order[jj]));
            fe_elem->attach_quadrature_rule(qrule.get());
            //  for evaluating the displacement system
            const std::vector<Real>& JxW = fe_elem->get_JxW();
            const std::vector<std::vector<Real> >& phi = fe_elem->get_phi();
            const std::vector<libMesh::Point>& qp_points = fe_elem->get_xyz();
            std::vector<dof_id_type> dof_indices;
            DenseMatrix<double> disp_coords;

            // loop over ALL elements in meter mesh, not just the local ones on this process!!
            MeshBase::const_element_iterator el = d_meter_meshes[jj]->active_elements_begin();
            const MeshBase::const_element_iterator end_el = d_meter_meshes[jj]->active_elements_end();
            for (; el != end_el; ++el)
            {
                const Elem* elem = *el;
                fe_elem->reinit(elem);

                // get dofs for displacement system and store in dense matrix
                disp_coords.resize(NDIM, phi.size());
                for (unsigned int d = 0; d < NDIM; ++d) // here d is the "variable number"
                {
                    dof_map.dof_indices(elem, dof_indices, d);
                    for (unsigned int nn = 0; nn < dof_indices.size(); ++nn)
                    {
                        disp_coords(d, nn) = displacement_coords(dof_indices[nn]);
                    }
                }

                // compute normal vector to element
                const libMesh::Point tau1 = *elem->node_ptr(1) - *elem->node_ptr(0);
                const libMesh::Point tau2 = *elem->node_ptr(2) - *elem->node_ptr(1);
                libMesh::Point normal_temp = tau1.cross(tau2).unit();
                Vector normal;
                for (unsigned int d = 0; d < NDIM; ++d) normal[d] = normal_temp(d);

                // loop over quadrature points, compute their physical locations
                // after displacement, and store their indices.
                for (unsigned int qp = 0; qp < qp_points.size(); ++qp)
                {
                    Vector qp_temp;
                    double disp_comp = 0.0;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        disp_comp = 0.0;
                        for (unsigned int nn = 0; nn < phi.size(); ++nn)
                        {
                            disp_comp += disp_coords(d, nn) * phi[nn][qp];
                        }
                        // calculating physical location of the quadrature point
                        qp_temp[d] = qp_points[qp](d) + disp_comp;
                    }

                    const Index<NDIM> i = IndexUtilities::getCellIndex(&qp_temp[0],
                                                                       domainXLower,
                                                                       domainXUpper,
                                                                       dx.data(),
                                                                       domain_box_level_lower,
                                                                       domain_box_level_upper);

                    const Index<NDIM> finer_i = IndexUtilities::getCellIndex(&qp_temp[0],
                                                                             domainXLower,
                                                                             domainXUpper,
                                                                             finer_dx.data(),
                                                                             finer_domain_box_level_lower,
                                                                             finer_domain_box_level_upper);

                    if (level->getBoxes().contains(i) &&
                        (ln == finest_ln || !finer_level->getBoxes().contains(finer_i)))
                    {
                        QuadPointStruct q;
                        q.meter_num = jj;
                        q.qp_xyz_current = qp_temp;
                        q.JxW = JxW[qp];
                        q.normal = normal;
                        d_quad_point_map[ln].insert(std::make_pair(i, q));
                    }
                }
            }
        }
    }
}

void
IBFEInstrumentPanel::readInstrumentData(const int U_data_idx,
                                        const int P_data_idx,
                                        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                        const double data_time)
{
    if (d_num_meters == 0) return;

    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    // Reset the instrument values.
    std::fill(d_flow_values.begin(), d_flow_values.end(), 0.0);
    std::fill(d_mean_pressure_values.begin(), d_mean_pressure_values.end(), 0.0);
    std::vector<double> A(d_num_meters, 0.0);

    // local counters for checking whether we have consistent
    // values for the number of quadrature points.
    int count_qp_1 = 0;
    int count_qp_2 = 0;

    // compute flow and mean pressure on mesh meters
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        count_qp_1 += d_quad_point_map[ln].size();

        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Index<NDIM>& patch_lower = patch_box.lower();
            const Index<NDIM>& patch_upper = patch_box.upper();

            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const x_lower = pgeom->getXLower();
            const double* const x_upper = pgeom->getXUpper();
            const double* const dx = pgeom->getDx();

            Pointer<CellData<NDIM, double> > U_cc_data = patch->getPatchData(U_data_idx);
            Pointer<SideData<NDIM, double> > U_sc_data = patch->getPatchData(U_data_idx);
            Pointer<CellData<NDIM, double> > P_cc_data = patch->getPatchData(P_data_idx);

            for (Box<NDIM>::Iterator b(patch_box); b; b++)
            {
                const Index<NDIM>& i = b();
                std::pair<QuadPointMap::const_iterator, QuadPointMap::const_iterator> qp_range =
                    d_quad_point_map[ln].equal_range(i);
                if (qp_range.first != qp_range.second)
                {
                    const Vector X_cell(x_lower[0] + dx[0] * (static_cast<double>(i(0) - patch_lower(0)) + 0.5),
                                        x_lower[1] + dx[1] * (static_cast<double>(i(1) - patch_lower(1)) + 0.5)
#if (NDIM == 3)
                                            ,
                                        x_lower[2] + dx[2] * (static_cast<double>(i(2) - patch_lower(2)) + 0.5)
#endif
                    );
                    if (U_cc_data)
                    {
                        for (auto it = qp_range.first; it != qp_range.second; ++it)
                        {
                            const int& meter_num = it->second.meter_num;
                            const double& JxW = it->second.JxW;
                            const Vector& X = it->second.qp_xyz_current;
                            const Vector& normal = it->second.normal;
                            const Vector U = linear_interp<NDIM>(
                                X, i, X_cell, *U_cc_data, patch_lower, patch_upper, x_lower, x_upper, dx);
                            d_flow_values[meter_num] += (U.dot(normal)) * JxW;
                        }
                    }
                    if (U_sc_data)
                    {
                        for (auto it = qp_range.first; it != qp_range.second; ++it)
                        {
                            const int& meter_num = it->second.meter_num;
                            const double& JxW = it->second.JxW;
                            const Vector& X = it->second.qp_xyz_current;
                            const Vector& normal = it->second.normal;
                            const Vector U =
                                linear_interp(X, i, X_cell, *U_sc_data, patch_lower, patch_upper, x_lower, x_upper, dx);
                            d_flow_values[meter_num] += (U.dot(normal)) * JxW;
                        }
                    }
                    if (P_cc_data)
                    {
                        for (auto it = qp_range.first; it != qp_range.second; ++it)
                        {
                            const int& meter_num = it->second.meter_num;
                            const double& JxW = it->second.JxW;
                            const Vector& X = it->second.qp_xyz_current;
                            double P =
                                linear_interp(X, i, X_cell, *P_cc_data, patch_lower, patch_upper, x_lower, x_upper, dx);
                            d_mean_pressure_values[meter_num] += P * JxW;
                            A[meter_num] += JxW;
                            count_qp_2 += 1;
                        }
                    }
                }
            }
        }
    }

    // check to make sure we don't double count quadrature points because
    // of overlapping patches or something else.
    const int count_qp_3 = SAMRAI_MPI::sumReduction(count_qp_2);
    if (count_qp_1 != count_qp_3)
    {
        TBOX_WARNING("IBFEInstrumentPanel::readInstrumentData :"
                     << " the total number of quadrature points in the meter meshes"
                     << " is not consistent with the number used in the"
                     << " calculations, for IBFE part " << d_part << "."
                     << " there may be overlapping patches in the AMR grid.");
    }

    // Synchronize the values across all processes.
    SAMRAI_MPI::sumReduction(&d_flow_values[0], d_num_meters);
    SAMRAI_MPI::sumReduction(&d_mean_pressure_values[0], d_num_meters);
    SAMRAI_MPI::sumReduction(&A[0], d_num_meters);

    // Normalize the mean pressure.
    for (unsigned int jj = 0; jj < d_num_meters; ++jj)
    {
        d_mean_pressure_values[jj] /= A[jj];
    }

    // we need to compute the flow correction by calculating the contribution
    // from the velocity of each meter mesh.
    for (unsigned int jj = 0; jj < d_num_meters; ++jj)
    {
        // get displacement and velocity systems for meter mesh
        const LinearImplicitSystem& velocity_sys =
            d_meter_systems[jj]->get_system<LinearImplicitSystem>(IBFEMethod::VELOCITY_SYSTEM_NAME);
        const NumericVector<double>& velocity_coords = velocity_sys.get_vector("serial solution");
        const DofMap& dof_map = velocity_sys.get_dof_map();
        FEType fe_type = velocity_sys.variable_type(0);

        // set up FE objects
        std::unique_ptr<FEBase> fe_elem(FEBase::build(NDIM - 1, fe_type));
        QGauss qrule(NDIM - 1, fe_type.default_quadrature_order());
        fe_elem->attach_quadrature_rule(&qrule);

        //  for evaluating the velocity system
        const std::vector<Real>& JxW = fe_elem->get_JxW();
        const std::vector<std::vector<Real> >& phi = fe_elem->get_phi();
        const std::vector<libMesh::Point>& qp_points = fe_elem->get_xyz();
        std::vector<dof_id_type> dof_indices;
        DenseMatrix<double> vel_coords;

        // loop over elements again to compute mass flux and mean pressure
        double flux_correction = 0.0;
        MeshBase::const_element_iterator el = d_meter_meshes[jj]->active_local_elements_begin();
        const MeshBase::const_element_iterator end_el = d_meter_meshes[jj]->active_local_elements_end();
        for (; el != end_el; ++el)
        {
            const Elem* elem = *el;
            fe_elem->reinit(elem);

            // get dofs for displacement system and store in dense matrix
            vel_coords.resize(NDIM, phi.size());
            for (unsigned int d = 0; d < NDIM; ++d) // here d is the "variable number"
            {
                dof_map.dof_indices(elem, dof_indices, d);
                for (unsigned int nn = 0; nn < dof_indices.size(); ++nn)
                {
                    vel_coords(d, nn) = velocity_coords(dof_indices[nn]);
                }
            }

            // compute normal vector to element
            const libMesh::Point tau1 = *elem->node_ptr(1) - *elem->node_ptr(0);
            const libMesh::Point tau2 = *elem->node_ptr(2) - *elem->node_ptr(1);
            const libMesh::Point normal = (tau1.cross(tau2)).unit();

            // loop over quadrature points
            double vel_comp;
            for (unsigned int qp = 0; qp < qp_points.size(); ++qp)
            {
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    vel_comp = 0.0;
                    for (unsigned int nn = 0; nn < phi.size(); ++nn)
                    {
                        vel_comp += vel_coords(d, nn) * phi[nn][qp];
                    }
                    flux_correction += vel_comp * normal(d) * JxW[qp];
                }
            }
        }

        const double total_correction = SAMRAI_MPI::sumReduction(flux_correction);
        d_flow_values[jj] -= total_correction;

    } // loop over meters

    // write data
    outputData(data_time);
}

void
IBFEInstrumentPanel::getFromInput(Pointer<Database> db)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(db);
#endif
    d_plot_directory_name = db->getString("meters_directory_name");
    d_instrument_dump_interval = db->getIntegerWithDefault("meters_dump_interval", 1);
    d_nodeset_IDs_for_meters = db->getIntegerArray("nodeset_IDs_for_meters");
    d_use_adaptive_quadrature = db->getBoolWithDefault("meters_adaptive_quadrature", false);
    d_quad_type = Utility::string_to_enum<QuadratureType>(db->getStringWithDefault("meters_quad_type", "QGAUSS"));
    d_input_quad_order = Utility::string_to_enum<Order>(db->getStringWithDefault("meters_quad_order", "FORTIETH"));
    if (d_use_adaptive_quadrature && d_quad_type != libMesh::QGRID)
    {
        TBOX_ERROR("IBFEInstrumentPanel::getFromInput :"
                   << " Adaptive quadrature for the meters"
                   << " is only supported with QuadratureType QGRID.");
    }
    return;
}

int
IBFEInstrumentPanel::getInstrumentDumpInterval() const
{
    return d_instrument_dump_interval;
}

std::string
IBFEInstrumentPanel::getPlotDirectoryName() const
{
    return d_plot_directory_name;
}

int
IBFEInstrumentPanel::getNumberOfMeterMeshes() const
{
    return d_num_meters;
}

MeshBase&
IBFEInstrumentPanel::getMeterMesh(const unsigned int jj) const
{
    return *d_meter_meshes[jj];
}

EquationSystems&
IBFEInstrumentPanel::getMeterMeshEquationSystems(const unsigned int jj) const
{
    return *d_meter_systems[jj];
}

const std::string&
IBFEInstrumentPanel::getMeterMeshName(const unsigned int jj) const
{
    return d_meter_mesh_names[jj];
}

const std::vector<libMesh::Point>&
IBFEInstrumentPanel::getMeterMeshNodes(const unsigned int jj) const
{
    return d_nodes[jj];
}

Order
IBFEInstrumentPanel::getMeterMeshQuadOrder(const unsigned int jj) const
{
    return d_quad_order[jj];
}

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBFEInstrumentPanel::initializeSystemDependentData(IBFEMethod* ib_method_ops, const int meter_mesh_number)
{
    // get the coordinate mapping system and velocity systems for the parent mesh
    const FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager(d_part);
    const EquationSystems* equation_systems = fe_data_manager->getEquationSystems();
    const System& dX_system = equation_systems->get_system(IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
    // \todo: find a better way to do this
    std::vector<double> dX_coords_parent;
    dX_system.update_global_solution(dX_coords_parent);
    const System& U_system = equation_systems->get_system(IBFEMethod::VELOCITY_SYSTEM_NAME);
    std::vector<double> U_coords_parent;
    U_system.update_global_solution(U_coords_parent);

    // get displacement and velocity systems for meter mesh
    auto& velocity_sys =
        d_meter_systems[meter_mesh_number]->get_system<LinearImplicitSystem>(IBFEMethod::VELOCITY_SYSTEM_NAME);
    const unsigned int velocity_sys_num = velocity_sys.number();
    NumericVector<double>& velocity_solution = *velocity_sys.solution;
    NumericVector<double>& velocity_coords = velocity_sys.get_vector("serial solution");

    auto& displacement_sys =
        d_meter_systems[meter_mesh_number]->get_system<LinearImplicitSystem>(IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
    const unsigned int displacement_sys_num = displacement_sys.number();
    NumericVector<double>& displacement_solution = *displacement_sys.solution;
    NumericVector<double>& displacement_coords = displacement_sys.get_vector("serial solution");

    // loop over the (perimeter) nodes in the meter mesh
    std::vector<double> mean_U_dofs;
    mean_U_dofs.resize(NDIM);
    std::vector<double> mean_dX_dofs;
    mean_dX_dofs.resize(NDIM);
    for (unsigned int ii = 0; ii < d_num_nodes[meter_mesh_number]; ++ii)
    {
        // get node on meter mesh
        const Node* node = &d_meter_meshes[meter_mesh_number]->node_ref(ii);

        // get corresponding dofs on parent mesh
        std::vector<double> U_dofs;
        U_dofs.resize(NDIM);
        std::vector<double> dX_dofs;
        dX_dofs.resize(NDIM);

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            U_dofs[d] = U_coords_parent[d_U_dof_idx[meter_mesh_number][ii][d]];
            dX_dofs[d] = dX_coords_parent[d_dX_dof_idx[meter_mesh_number][ii][d]];
            mean_U_dofs[d] += U_dofs[d] / static_cast<double>(d_num_nodes[meter_mesh_number]);
            mean_dX_dofs[d] += dX_dofs[d] / static_cast<double>(d_num_nodes[meter_mesh_number]);
        }

        // set dofs in meter mesh to correspond to the same values
        // as in the parent mesh
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            const int vel_dof_idx = node->dof_number(velocity_sys_num, d, 0);
            velocity_coords.set(vel_dof_idx, U_dofs[d]);
            const int disp_dof_idx = node->dof_number(displacement_sys_num, d, 0);
            displacement_coords.set(disp_dof_idx, dX_dofs[d]);
        }
    }

    // set dofs for the centroid node in the meter mesh
    const Node* centroid_node = &d_meter_meshes[meter_mesh_number]->node_ref(d_num_nodes[meter_mesh_number]);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        const int vel_dof_idx = centroid_node->dof_number(velocity_sys_num, d, 0);
        velocity_coords.set(vel_dof_idx, mean_U_dofs[d]);
        const int disp_dof_idx = centroid_node->dof_number(displacement_sys_num, d, 0);
        displacement_coords.set(disp_dof_idx, mean_dX_dofs[d]);
    }

    // also populate solution vector in the system
    MeshBase::const_node_iterator node_it = d_meter_meshes[meter_mesh_number]->local_nodes_begin();
    const MeshBase::const_node_iterator end_node_it = d_meter_meshes[meter_mesh_number]->local_nodes_end();
    for (; node_it != end_node_it; ++node_it)
    {
        const Node* node = *node_it;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            const int vel_dof_idx = node->dof_number(velocity_sys_num, d, 0);
            velocity_solution.set(vel_dof_idx, velocity_coords(vel_dof_idx));
            const int disp_dof_idx = node->dof_number(displacement_sys_num, d, 0);
            displacement_solution.set(disp_dof_idx, displacement_coords(disp_dof_idx));
        }
    }
    velocity_solution.close();
    displacement_solution.close();

    // compute the meter radius
    double max_meter_radius = 0.0;
    for (unsigned int ii = 0; ii < d_num_nodes[meter_mesh_number]; ++ii)
    {
        // get node on meter mesh
        const Node* node = &d_meter_meshes[meter_mesh_number]->node_ref(ii);
        // get the centroid
        const Node* centroid_node = &d_meter_meshes[meter_mesh_number]->node_ref(d_num_nodes[meter_mesh_number]);

        std::vector<double> node_disp(NDIM, 0.0);
        std::vector<double> centroid_disp(NDIM, 0.0);
        std::vector<numeric_index_type> node_disp_dof_idx(NDIM, 0);
        std::vector<numeric_index_type> centroid_disp_dof_idx(NDIM, 0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            node_disp_dof_idx[d] = node->dof_number(displacement_sys_num, d, 0);
            centroid_disp_dof_idx[d] = centroid_node->dof_number(displacement_sys_num, d, 0);
        }
        displacement_coords.get(node_disp_dof_idx, node_disp);
        displacement_coords.get(centroid_disp_dof_idx, centroid_disp);

        double radius_squared = 0.0;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            radius_squared += std::pow((*centroid_node)(d) + centroid_disp[d] - (*node)(d) + node_disp[d], 2.0);
        }
        max_meter_radius = std::max(std::pow(radius_squared, 0.5), max_meter_radius);
    }
    d_meter_radii[meter_mesh_number] = max_meter_radius;
}

double
IBFEInstrumentPanel::getMeterRadius(const int meter_mesh_number)
{
    // NOTE: this function should be called **after** updating the meter system data,
    // with initializeSystemDependentData.
    return d_meter_radii[meter_mesh_number];
}

void
IBFEInstrumentPanel::outputData(const double data_time)
{
    if (SAMRAI_MPI::getRank() == 0)
    {
        d_mean_pressure_stream << data_time;
        d_flux_stream << data_time;
        for (unsigned int jj = 0; jj < d_num_meters; ++jj)
        {
            d_mean_pressure_stream << " " << d_mean_pressure_values[jj];
            d_flux_stream << " " << d_flow_values[jj];
        }
        d_mean_pressure_stream << "\n";
        d_flux_stream << "\n";
    }
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
