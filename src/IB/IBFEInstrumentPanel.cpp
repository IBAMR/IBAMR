// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2022 by the IBAMR developers
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

#include "ibamr/IBFEInstrumentPanel.h"

#include "ibtk/IBTK_MPI.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/ibtk_utilities.h"

#include "libmesh/boundary_info.h"
#include "libmesh/face_tri3.h"
#include "libmesh/id_types.h"

#include <limits>
#include <set>

#include "ibamr/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline libMesh::Point
put_point_in_domain(const libMesh::Point& x, const double* const domain_x_lower, const double* const domain_x_upper)
{
    const double TOL = sqrt(std::numeric_limits<double>::epsilon());
    libMesh::Point x_corrected = x;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (x(d) <= domain_x_lower[d]) x_corrected(d) = domain_x_lower[d] + TOL;
        if (x(d) >= domain_x_upper[d]) x_corrected(d) = domain_x_upper[d] - TOL;
    }
    return x_corrected;
}
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBFEInstrumentPanel::IBFEInstrumentPanel(tbox::Pointer<tbox::Database> input_db, const int part)
    : d_part(part), d_plot_directory_name(NDIM == 2 ? "viz_inst2d" : "viz_inst3d")
{
    // get input data
    IBFEInstrumentPanel::getFromInput(input_db);

    // make plot directory
    Utilities::recursiveMkdir(d_plot_directory_name);
}

void
IBFEInstrumentPanel::initializeHierarchyIndependentData(IBFEMethod* const ib_method_ops)
{
    if (d_initialized) return;

    // Get major FE data structures.
    const FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager(d_part);
    const EquationSystems* equation_systems = fe_data_manager->getEquationSystems();
    const MeshBase& structure_mesh = equation_systems->get_mesh();
    const BoundaryInfo& boundary_info = structure_mesh.get_boundary_info();
    const libMesh::Parallel::Communicator& comm_in = structure_mesh.comm();

    const auto& x_mesh_system = equation_systems->get_system(ib_method_ops->getCurrentCoordinatesSystemName());
    const unsigned int x_mesh_sys_num = x_mesh_system.number();
    const auto& u_mesh_system = equation_systems->get_system(ib_method_ops->getVelocitySystemName());
    const unsigned int u_mesh_sys_num = u_mesh_system.number();

    // Get boundary information.
    std::vector<dof_id_type> nodes;
    std::vector<boundary_id_type> bcs;
    // new API in 1.4.0
#if LIBMESH_VERSION_LESS_THAN(1, 4, 0)
    boundary_info.build_node_list(nodes, bcs);
#else
    const std::vector<std::tuple<dof_id_type, boundary_id_type> > node_list = boundary_info.build_node_list();
    for (const std::tuple<dof_id_type, boundary_id_type>& pair : node_list)
    {
        nodes.push_back(std::get<0>(pair));
        bcs.push_back(std::get<1>(pair));
    }
#endif

    // Check to make sure there are node sets to work with.
    if (nodes.size() == 0 || bcs.size() == 0 || (nodes.size() != bcs.size()))
    {
        TBOX_ERROR("IBFEInstrumentPanel::initializeHierarchyIndependentData : "
                   << "nodesets not set up correctly or don't exist in FE part with number " << d_part);
    }

    // Resize members and local variables.
    d_num_meters = d_perimeter_nodeset_ids.size();
    d_meter_radii.resize(d_num_meters);
    d_quad_order.resize(d_num_meters);
    d_x_dof_idx.resize(d_num_meters);
    d_u_dof_idx.resize(d_num_meters);
    d_num_perim_nodes.resize(d_num_meters);
    d_flow_rate_values.resize(d_num_meters);
    d_mean_pressure_values.resize(d_num_meters);
    d_centroid_pressure_values.resize(d_num_meters);

    // Determine the number of points in each collection of perimeter nodes, keep track of the node indices for the
    // perimeter nodes, and evaluate the position of the centroid of each meter mesh.
    std::vector<libMesh::Point> meter_centroid(d_num_meters);
    std::vector<std::set<libMesh::dof_id_type> > structure_perimeter_node_ids(d_num_meters);
    for (unsigned int i = 0; i < structure_mesh.n_nodes(); ++i)
    {
        const Node* const node_ptr = structure_mesh.node_ptr(i);
        std::vector<short int> bdry_ids;
        boundary_info.boundary_ids(node_ptr, bdry_ids);
        for (unsigned int meter_idx = 0; meter_idx < d_num_meters; ++meter_idx)
        {
            const int nodeset_id = d_perimeter_nodeset_ids[meter_idx];
            if (find(bdry_ids.begin(), bdry_ids.end(), nodeset_id) != bdry_ids.end())
            {
                d_num_perim_nodes[meter_idx] += 1;
                structure_perimeter_node_ids[meter_idx].insert(node_ptr->id());
                meter_centroid[meter_idx] += *node_ptr;
            }
        }
    }
    for (unsigned int meter_idx = 0; meter_idx < d_num_meters; ++meter_idx)
    {
        meter_centroid[meter_idx] /= d_num_perim_nodes[meter_idx];
    }

    // Keep track of all elements associated with each perimeter node.
    std::vector<std::map<libMesh::dof_id_type, std::set<libMesh::dof_id_type> > > structure_elem_to_node_map(
        d_num_meters);
    std::vector<std::map<libMesh::dof_id_type, std::set<libMesh::dof_id_type> > > structure_node_to_elem_map(
        d_num_meters);
    for (unsigned int e = 0; e < structure_mesh.n_elem(); ++e)
    {
        const Elem* const elem_ptr = structure_mesh.elem_ptr(e);
        const unsigned int elem_id = elem_ptr->id();
        for (unsigned int i = 0; i < elem_ptr->n_nodes(); ++i)
        {
            const Node* const node_ptr = elem_ptr->node_ptr(i);
            const unsigned int node_id = node_ptr->id();
            std::vector<short int> bdry_ids;
            boundary_info.boundary_ids(node_ptr, bdry_ids);
            for (unsigned int meter_idx = 0; meter_idx < d_num_meters; ++meter_idx)
            {
                const int nodeset_id = d_perimeter_nodeset_ids[meter_idx];
                if (find(bdry_ids.begin(), bdry_ids.end(), nodeset_id) != bdry_ids.end())
                {
                    structure_elem_to_node_map[meter_idx][elem_id].insert(node_id);
                    structure_node_to_elem_map[meter_idx][node_id].insert(elem_id);
                }
            }
        }
    }

    // Create an ordered list of perimeter nodes.
    std::vector<std::vector<const Node*> > structure_nodes(d_num_meters);
    for (unsigned int meter_idx = 0; meter_idx < d_num_meters; ++meter_idx)
    {
        std::set<libMesh::dof_id_type> unassigned_node_ids = structure_perimeter_node_ids[meter_idx];
        auto first_node_id = *unassigned_node_ids.begin();
        auto first_node_ptr = structure_mesh.node_ptr(first_node_id);
        structure_nodes[meter_idx].push_back(first_node_ptr);
        unassigned_node_ids.erase(first_node_id);
        while (!unassigned_node_ids.empty())
        {
            const auto current_node_ptr = structure_nodes[meter_idx].back();
            const auto current_node_id = current_node_ptr->id();
            std::set<libMesh::dof_id_type> candidate_node_ids;

            // First check to see if there are unassigned nodes that share an element with the current node.
            for (const auto& candidate_elem_id : structure_node_to_elem_map[meter_idx][current_node_id])
            {
                for (const auto& node_id : structure_elem_to_node_map[meter_idx][candidate_elem_id])
                {
                    if (unassigned_node_ids.count(node_id) > 0)
                    {
                        candidate_node_ids.insert(node_id);
                    }
                }
            }

            // If there are no candidate node IDs identified by the mappings between elements and nodes, then we
            // simply search over *all* of the unassigned nodes.
            if (candidate_node_ids.empty())
            {
                candidate_node_ids = unassigned_node_ids;
            }
            TBOX_ASSERT(!candidate_node_ids.empty());

            if (candidate_node_ids.size() == 1)
            {
                // If we find only one node that shares an element with the current node, then that node must be the
                // next one in the perimeter mesh.
                auto next_node_id = *candidate_node_ids.begin();
                auto next_node_ptr = structure_mesh.node_ptr(next_node_id);
                structure_nodes[meter_idx].push_back(next_node_ptr);
                unassigned_node_ids.erase(next_node_id);
            }
            else
            {
                // There are multiple candidate nodes; choose the next node to be the closest one.
                //
                // If multiple nodes have the same distance from the current one to machine precision, pick the node
                // with the smallest ID.
                double min_distance = std::numeric_limits<double>::max();
                static const auto invalid_node_id = std::numeric_limits<libMesh::dof_id_type>::max();
                auto next_node_id = invalid_node_id;
                const Node* next_node_ptr = nullptr;
                for (const auto& candidate_node_id : candidate_node_ids)
                {
                    const auto candidate_node_ptr = structure_mesh.node_ptr(candidate_node_id);
                    const double distance = ((*candidate_node_ptr) - (*current_node_ptr)).norm();
                    if (((distance < min_distance) && !rel_equal_eps(distance, min_distance)) ||
                        (rel_equal_eps(distance, min_distance) && (candidate_node_id < next_node_id)))
                    {
                        next_node_id = candidate_node_id;
                        next_node_ptr = candidate_node_ptr;
                        min_distance = distance;
                    }
                }
                TBOX_ASSERT(next_node_ptr != nullptr);
                structure_nodes[meter_idx].push_back(next_node_ptr);
                unassigned_node_ids.erase(next_node_id);
            }
        }
        TBOX_ASSERT(structure_nodes[meter_idx].size() == d_num_perim_nodes[meter_idx]);
    }

    // Initialize the meter meshes and number of nodes
    for (unsigned int meter_idx = 0; meter_idx < d_num_meters; ++meter_idx)
    {
        d_num_perim_nodes[meter_idx] = structure_nodes[meter_idx].size();
        d_meter_meshes.emplace_back(new SerialMesh(comm_in, NDIM));
        d_meter_meshes[meter_idx]->set_spatial_dimension(NDIM);
        d_meter_meshes[meter_idx]->set_mesh_dimension(NDIM - 1);
        d_meter_meshes[meter_idx]->reserve_nodes(d_num_perim_nodes[meter_idx] + 1);
        d_meter_meshes[meter_idx]->reserve_elem(d_num_perim_nodes[meter_idx]);
        d_meter_mesh_names.push_back("meter_mesh_" + std::to_string(d_perimeter_nodeset_ids[meter_idx]));

        // Insert the perimeter nodes first, then add the centroid last.
        for (unsigned int i = 0; i < d_num_perim_nodes[meter_idx]; ++i)
        {
            d_meter_meshes[meter_idx]->add_point(*structure_nodes[meter_idx][i], i);
        }
        d_meter_meshes[meter_idx]->add_point(meter_centroid[meter_idx], d_num_perim_nodes[meter_idx]);

        // Create triangular elements.
        for (unsigned int i = 0; i < d_num_perim_nodes[meter_idx]; ++i)
        {
            Elem* elem = new Tri3;
            elem->set_id(i);
            // libMesh will delete elem
            elem = d_meter_meshes[meter_idx]->add_elem(elem);
            elem->set_node(0) = d_meter_meshes[meter_idx]->node_ptr(d_num_perim_nodes[meter_idx]);
            elem->set_node(1) = d_meter_meshes[meter_idx]->node_ptr(i);
            elem->set_node(2) = d_meter_meshes[meter_idx]->node_ptr((i + 1) % d_num_perim_nodes[meter_idx]);
        }
        d_meter_meshes[meter_idx]->allow_renumbering(false);
        d_meter_meshes[meter_idx]->prepare_for_use();

        d_meter_systems.emplace_back(new EquationSystems(*d_meter_meshes[meter_idx]));
        auto& velocity_sys = d_meter_systems[meter_idx]->add_system<System>(ib_method_ops->getVelocitySystemName());
        velocity_sys.add_variable("U_0", static_cast<Order>(1), LAGRANGE);
        velocity_sys.add_variable("U_1", static_cast<Order>(1), LAGRANGE);
        velocity_sys.add_variable("U_2", static_cast<Order>(1), LAGRANGE);
        velocity_sys.add_vector("serial solution", false, libMesh::SERIAL);

        d_meter_systems[meter_idx]->init();

        for (unsigned int i = 0; i < d_num_perim_nodes[meter_idx]; ++i)
        {
            const Node* const node = structure_nodes[meter_idx][i];
            std::array<dof_id_type, NDIM> x_dof_idx, u_dof_idx;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                x_dof_idx[d] = node->dof_number(x_mesh_sys_num, d, 0);
                u_dof_idx[d] = node->dof_number(u_mesh_sys_num, d, 0);
            }
            d_x_dof_idx[meter_idx].push_back(x_dof_idx);
            d_u_dof_idx[meter_idx].push_back(u_dof_idx);
        }
    }
    d_initialized = true;
    return;
}

void
IBFEInstrumentPanel::readInstrumentData(const int U_data_idx,
                                        const int P_data_idx,
                                        const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
                                        IBFEMethod* const ib_method_ops,
                                        const double data_time)
{
    if (d_num_meters == 0) return;

    if (!d_initialized) initializeHierarchyIndependentData(ib_method_ops);

    // Update the meter meshes to correspond to the current configuration of the FE mesh.
    for (unsigned int meter_idx = 0; meter_idx < d_num_meters; ++meter_idx)
    {
        resetMeterConfiguration(ib_method_ops, meter_idx);
    }

    // Evaluate the meter quadrature points and determine mappings from patch numbers to particular quadrature points.
    std::vector<std::map<int, std::vector<int> > > meter_idx_map;
    std::vector<std::map<int, std::vector<Vector> > > meter_x_map, meter_u_corr_map, meter_normal_map;
    std::vector<std::map<int, std::vector<double> > > meter_JxW_map;
    computeMeterQuadratureData(
        meter_idx_map, meter_x_map, meter_u_corr_map, meter_normal_map, meter_JxW_map, hierarchy, ib_method_ops);

    // Compute flow rate and pressure on mesh meter faces.
    std::fill(d_flow_rate_values.begin(), d_flow_rate_values.end(), 0.0);
    std::fill(d_mean_pressure_values.begin(), d_mean_pressure_values.end(), 0.0);
    std::fill(d_centroid_pressure_values.begin(), d_centroid_pressure_values.end(), 0.0);
    std::vector<double> A(d_num_meters, 0.0);
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            if (meter_x_map[ln].find(p()) == meter_x_map[ln].end()) continue;

            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const std::vector<int>& meter_idx = meter_idx_map[ln][p()];
            const std::vector<Vector>& meter_x = meter_x_map[ln][p()];
            const std::vector<Vector>& meter_u_corr = meter_u_corr_map[ln][p()];
            const std::vector<Vector>& meter_normal = meter_normal_map[ln][p()];
            const std::vector<double>& meter_JxW = meter_JxW_map[ln][p()];

            const auto n_qp = meter_idx.size();
            std::vector<double> u_qp(NDIM * n_qp, 0.0);
            std::vector<double> p_qp(n_qp, 0.0);
            std::vector<double> x_qp;
            x_qp.reserve(NDIM * n_qp);
            for (const auto& x : meter_x)
            {
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    x_qp.push_back(x[d]);
                }
            }

            Pointer<CellData<NDIM, double> > u_cc_data = patch->getPatchData(U_data_idx);
            Pointer<SideData<NDIM, double> > u_sc_data = patch->getPatchData(U_data_idx);
            if (u_cc_data)
            {
                LEInteractor::interpolate(
                    u_qp, NDIM, x_qp, NDIM, u_cc_data, patch, u_cc_data->getGhostBox(), d_u_interp_fcn);
            }
            else if (u_sc_data)
            {
                LEInteractor::interpolate(
                    u_qp, NDIM, x_qp, NDIM, u_sc_data, patch, u_sc_data->getGhostBox(), d_u_interp_fcn);
            }
            else
            {
                TBOX_ERROR("no velocity data!\n");
            }
            Pointer<CellData<NDIM, double> > p_cc_data = patch->getPatchData(P_data_idx);
            if (p_cc_data)
            {
                LEInteractor::interpolate(
                    p_qp, 1, x_qp, NDIM, p_cc_data, patch, p_cc_data->getGhostBox(), d_p_interp_fcn);
            }
            else
            {
                TBOX_ERROR("no pressure data!\n");
            }

            for (unsigned qp = 0; qp < n_qp; ++qp)
            {
                const Vector u(u_qp.data() + NDIM * qp);
                const Vector& u_corr = meter_u_corr[qp];
                const double p = p_qp[qp];
                const Vector& n = meter_normal[qp];
                const double JxW = meter_JxW[qp];
                d_flow_rate_values[meter_idx[qp]] += (u - u_corr).dot(n) * JxW;
                d_mean_pressure_values[meter_idx[qp]] += p * JxW;
                A[meter_idx[qp]] += JxW;
            }
        }
    }

    // Interpolate the pressure at the centroid of each meter mesh.
    //
    // TODO: Factor out common code for finding assignments of points to patches.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    const double* const domain_x_lower = grid_geom->getXLower();
    const double* const domain_x_upper = grid_geom->getXUpper();
    const double* const dx_coarsest = grid_geom->getDx();
    TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
    const Box<NDIM> domain_box = grid_geom->getPhysicalDomain()[0];
    for (unsigned int meter_idx = 0; meter_idx < d_num_meters; ++meter_idx)
    {
        const libMesh::Point x_centroid = put_point_in_domain(
            *d_meter_meshes[meter_idx]->node_ptr(d_num_perim_nodes[meter_idx]), domain_x_lower, domain_x_upper);

        // Find the level that contains the centroid.
        int x_centroid_global_ln = IBTK::invalid_level_number, x_centroid_local_ln = IBTK::invalid_level_number,
            x_centroid_local_patch_idx = IBTK::invalid_index;
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            Pointer<BoxTree<NDIM> > box_tree = level->getBoxTree();
            const IntVector<NDIM>& ratio = level->getRatio();
            const Box<NDIM> domain_box_level = Box<NDIM>::refine(domain_box, ratio);
            const hier::Index<NDIM>& domain_box_level_lower = domain_box_level.lower();
            const hier::Index<NDIM>& domain_box_level_upper = domain_box_level.upper();
            std::array<double, NDIM> dx;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dx[d] = dx_coarsest[d] / static_cast<double>(ratio(d));
            }
            const auto i = IndexUtilities::getCellIndex(&x_centroid(0),
                                                        domain_x_lower,
                                                        domain_x_upper,
                                                        dx.data(),
                                                        domain_box_level_lower,
                                                        domain_box_level_upper);
            Box<NDIM> cell_box(i, i);
            tbox::Array<int> local_patch_idxs;
            box_tree->findLocalOverlapIndices(local_patch_idxs, cell_box);
            if (local_patch_idxs.size() != 0)
            {
                TBOX_ASSERT(local_patch_idxs.size() == 1);
                x_centroid_local_patch_idx = local_patch_idxs[0];
                x_centroid_local_ln = ln;
            }
        }
        int bcast_root = -1;
        x_centroid_global_ln = IBTK_MPI::maxReduction(x_centroid_local_ln, &bcast_root);
        if (x_centroid_global_ln == x_centroid_local_ln)
        {
            // The centroid should be on a patch owned by exactly one processor
            TBOX_ASSERT(bcast_root == IBTK_MPI::getRank());
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(x_centroid_local_ln);
            if (x_centroid_local_patch_idx != IBTK::invalid_index)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(x_centroid_local_patch_idx);
                Pointer<CellData<NDIM, double> > p_cc_data = patch->getPatchData(P_data_idx);
                if (p_cc_data)
                {
                    LEInteractor::interpolate(&d_centroid_pressure_values[meter_idx],
                                              1,
                                              1,
                                              &x_centroid(0),
                                              NDIM,
                                              NDIM,
                                              p_cc_data,
                                              patch,
                                              p_cc_data->getGhostBox(),
                                              d_p_interp_fcn);
                }
                else
                {
                    TBOX_ERROR("no pressure data!\n");
                }
            }
        }
        d_centroid_pressure_values[meter_idx] = IBTK_MPI::bcast(d_centroid_pressure_values[meter_idx], bcast_root);
    }

    // Synchronize accumulated data.
    IBTK_MPI::sumReduction(d_flow_rate_values.data(), d_num_meters);
    IBTK_MPI::sumReduction(d_mean_pressure_values.data(), d_num_meters);
    IBTK_MPI::sumReduction(A.data(), d_num_meters);

    // Normalize the mean pressure.
    for (unsigned int meter_idx = 0; meter_idx < d_num_meters; ++meter_idx)
    {
        d_mean_pressure_values[meter_idx] /= A[meter_idx];
    }

    // Write data.
    outputData(data_time);
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
IBFEInstrumentPanel::getMeterMesh(const unsigned int meter_idx) const
{
    return *d_meter_meshes[meter_idx];
}

EquationSystems&
IBFEInstrumentPanel::getMeterMeshEquationSystems(const unsigned int meter_idx) const
{
    return *d_meter_systems[meter_idx];
}

const std::string&
IBFEInstrumentPanel::getMeterMeshName(const unsigned int meter_idx) const
{
    return d_meter_mesh_names[meter_idx];
}

QuadratureType
IBFEInstrumentPanel::getMeterMeshQuadType(const unsigned int /*meter_idx*/) const
{
    return d_quad_type;
}

Order
IBFEInstrumentPanel::getMeterMeshQuadOrder(const unsigned int meter_idx) const
{
    return d_quad_order[meter_idx];
}

const std::vector<double>&
IBFEInstrumentPanel::getMeterFlowRates() const
{
    return d_flow_rate_values;
}

const std::vector<double>&
IBFEInstrumentPanel::getMeterMeanPressures() const
{
    return d_mean_pressure_values;
}

const std::vector<double>&
IBFEInstrumentPanel::getMeterCentroidPressures() const
{
    return d_centroid_pressure_values;
}

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBFEInstrumentPanel::getFromInput(Pointer<Database> db)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(db);
#endif
    d_plot_directory_name = db->getString("meters_directory_name");
    d_instrument_dump_interval = db->getIntegerWithDefault("meters_dump_interval", 1);
    d_perimeter_nodeset_ids = db->getIntegerArray("nodeset_IDs_for_meters");
    d_use_adaptive_quadrature = db->getBoolWithDefault("meters_adaptive_quadrature", d_use_adaptive_quadrature);
    d_quad_type = Utility::string_to_enum<QuadratureType>(
        db->getStringWithDefault("meters_quad_type", Utility::enum_to_string<QuadratureType>(d_quad_type)));
    d_default_quad_order = Utility::string_to_enum<Order>(
        db->getStringWithDefault("meters_quad_order", Utility::enum_to_string<Order>(d_default_quad_order)));
    if (d_use_adaptive_quadrature && d_quad_type != libMesh::QGRID)
    {
        TBOX_ERROR("IBFEInstrumentPanel::getFromInput :"
                   << " Adaptive quadrature for the meters"
                   << " is only supported with QuadratureType QGRID.");
    }
    d_u_interp_fcn = db->getStringWithDefault("meters_velocity_interp_fcn", d_u_interp_fcn);
    d_p_interp_fcn = db->getStringWithDefault("meters_pressure_interp_fcn", d_p_interp_fcn);
    return;
}

void
IBFEInstrumentPanel::resetMeterConfiguration(IBFEMethod* const ib_method_ops, const int meter_idx)
{
    // Get systems that provide the current structural mesh coordinates and velocity.
    const FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager(d_part);
    const EquationSystems* equation_systems = fe_data_manager->getEquationSystems();
    const auto& x_mesh_system = equation_systems->get_system(ib_method_ops->getCurrentCoordinatesSystemName());
    const auto& x_mesh_vec = x_mesh_system.solution;
    const auto& u_mesh_system = equation_systems->get_system(ib_method_ops->getVelocitySystemName());
    const auto& u_mesh_vec = u_mesh_system.solution;

    // Get systems that provide the current meter mesh velocity.
    auto& u_meter_system = d_meter_systems[meter_idx]->get_system(ib_method_ops->getVelocitySystemName());
    auto& u_meter_vec = u_meter_system.solution;
    const int u_meter_sys_num = u_meter_system.number();

    const unsigned int n_dofs = u_meter_vec->size();
    TBOX_ASSERT(n_dofs == NDIM * (d_num_perim_nodes[meter_idx] + 1));
    std::vector<double> x_meter_dofs(n_dofs, 0.0), u_meter_dofs(n_dofs, 0.0);

    // Loop over the (perimeter) nodes in the meter and get the *local* data from the structural mesh.
    for (unsigned int i = 0; i < d_num_perim_nodes[meter_idx]; ++i)
    {
        // Look up the corresponding DOFs on the parent mesh.
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            const auto x_dof_idx = d_x_dof_idx[meter_idx][i][d];
            if (x_mesh_vec->first_local_index() <= x_dof_idx && x_dof_idx < x_mesh_vec->last_local_index())
            {
                x_meter_dofs[NDIM * i + d] = x_mesh_vec->el(x_dof_idx);
            }
            const auto u_dof_idx = d_u_dof_idx[meter_idx][i][d];
            if (u_mesh_vec->first_local_index() <= u_dof_idx && u_dof_idx < u_mesh_vec->last_local_index())
            {
                u_meter_dofs[NDIM * i + d] = u_mesh_vec->el(u_dof_idx);
            }
        }
    }

    // Broadcast the local data.
    IBTK_MPI::sumReduction(x_meter_dofs.data(), n_dofs);
    IBTK_MPI::sumReduction(u_meter_dofs.data(), n_dofs);

    // Compute the mean values and store them in the vector of DOFs.
    std::vector<double> x_mean(NDIM, 0.0), u_mean(NDIM, 0.0);
    for (unsigned int i = 0; i < d_num_perim_nodes[meter_idx]; ++i)
    {
        for (int d = 0; d < NDIM; ++d)
        {
            x_mean[d] += x_meter_dofs[NDIM * i + d] / static_cast<double>(d_num_perim_nodes[meter_idx]);
            u_mean[d] += u_meter_dofs[NDIM * i + d] / static_cast<double>(d_num_perim_nodes[meter_idx]);
        }
    }
    for (int d = 0; d < NDIM; ++d)
    {
        x_meter_dofs[NDIM * d_num_perim_nodes[meter_idx] + d] = x_mean[d];
        u_meter_dofs[NDIM * d_num_perim_nodes[meter_idx] + d] = u_mean[d];
    }

    // Set the values in the meter meshes.
    for (unsigned int i = 0; i <= d_num_perim_nodes[meter_idx]; ++i)
    {
        Node* node = &d_meter_meshes[meter_idx]->node_ref(i);
        for (int d = 0; d < NDIM; ++d)
        {
            (*node)(d) = x_meter_dofs[NDIM * i + d];
            const int u_dof_idx = node->dof_number(u_meter_sys_num, d, 0);
            u_meter_vec->set(u_dof_idx, u_meter_dofs[NDIM * i + d]);
        }
    }

    u_meter_vec->close();
    NumericVector<double>& u_meter_serial_vec = u_meter_system.get_vector("serial solution");
    u_meter_vec->localize(u_meter_serial_vec);
    u_meter_serial_vec.close();

    // Compute the meter radius.
    double max_meter_radius = 0.0;
    const libMesh::Point centroid_node = d_meter_meshes[meter_idx]->node_ref(d_num_perim_nodes[meter_idx]);
    for (unsigned int i = 0; i < d_num_perim_nodes[meter_idx]; ++i)
    {
        const libMesh::Point perim_node = d_meter_meshes[meter_idx]->node_ref(i);
        max_meter_radius = std::max(max_meter_radius, (perim_node - centroid_node).norm());
    }
    d_meter_radii[meter_idx] = max_meter_radius;
    return;
}

void
IBFEInstrumentPanel::computeMeterQuadratureData(std::vector<std::map<int, std::vector<int> > >& meter_idx_map,
                                                std::vector<std::map<int, std::vector<Vector> > >& meter_x_map,
                                                std::vector<std::map<int, std::vector<Vector> > >& meter_u_corr_map,
                                                std::vector<std::map<int, std::vector<Vector> > >& meter_normal_map,
                                                std::vector<std::map<int, std::vector<double> > >& meter_JxW_map,
                                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                                const IBFEMethod* const ib_method_ops)
{
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    // Ensure all data structures are the correct size.
    meter_idx_map.clear();
    meter_x_map.clear();
    meter_u_corr_map.clear();
    meter_normal_map.clear();
    meter_JxW_map.clear();

    meter_idx_map.resize(finest_ln + 1);
    meter_x_map.resize(finest_ln + 1);
    meter_u_corr_map.resize(finest_ln + 1);
    meter_normal_map.resize(finest_ln + 1);
    meter_JxW_map.resize(finest_ln + 1);

    // Determine the finest grid spacing in the Cartesian grid hierarchy.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    const double* const domain_x_lower = grid_geom->getXLower();
    const double* const domain_x_upper = grid_geom->getXUpper();
    const double* const dx_coarsest = grid_geom->getDx();
    TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
    const Box<NDIM> domain_box = grid_geom->getPhysicalDomain()[0];

    const IntVector<NDIM>& ratio_to_level_zero = hierarchy->getPatchLevel(finest_ln)->getRatio();
    std::array<double, NDIM> dx_finest;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        dx_finest[d] = dx_coarsest[d] / static_cast<double>(ratio_to_level_zero(d));
    }
    const double h_finest = *std::min_element(dx_finest.begin(), dx_finest.end());

    // Set quadrature rules for each meter mesh.
    for (unsigned int meter_idx = 0; meter_idx < d_num_meters; ++meter_idx)
    {
        if (d_use_adaptive_quadrature)
        {
            d_quad_order[meter_idx] = static_cast<Order>(d_meter_radii[meter_idx] / (0.25 * h_finest));
        }
        else
        {
            d_quad_order[meter_idx] = d_default_quad_order;
        }
        if (d_quad_type == libMesh::QGRID && d_quad_order[meter_idx] > libMesh::FORTYTHIRD)
        {
            TBOX_WARNING("IBFEInstrumentPanel::initializeHierarchyDependentData : "
                         << "QGrid quadrature order exceeds 43 for meter mesh in IBFE part " << d_part << "."
                         << " there may be undefined behavior in casting to this"
                         << " Order in older versions of libMesh.");
        }
    }

    // Loop over all levels and try to assign each quadrature point to a patch in each level.
    std::vector<std::vector<int> > meter_qp_global_ln(d_num_meters), meter_qp_local_ln(d_num_meters),
        meter_qp_local_patch_idx(d_num_meters);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        Pointer<BoxTree<NDIM> > box_tree = level->getBoxTree();
        const IntVector<NDIM>& ratio = level->getRatio();
        const Box<NDIM> domain_box_level = Box<NDIM>::refine(domain_box, ratio);
        const hier::Index<NDIM>& domain_box_level_lower = domain_box_level.lower();
        const hier::Index<NDIM>& domain_box_level_upper = domain_box_level.upper();
        std::array<double, NDIM> dx;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dx[d] = dx_coarsest[d] / static_cast<double>(ratio(d));
        }
        for (unsigned int meter_idx = 0; meter_idx < d_num_meters; ++meter_idx)
        {
            // Setup FE objects.
            auto& u_meter_system = d_meter_systems[meter_idx]->get_system(ib_method_ops->getVelocitySystemName());
            FEType fe_type = u_meter_system.variable_type(0);
            std::unique_ptr<FEBase> fe_elem(FEBase::build(NDIM - 1, fe_type));
            std::unique_ptr<QBase> qrule(QBase::build(d_quad_type, NDIM - 1, d_quad_order[meter_idx]));
            fe_elem->attach_quadrature_rule(qrule.get());
            const std::vector<libMesh::Point>& q_point = fe_elem->get_xyz();

            // Loop over ALL meter mesh elements.
            unsigned int meter_qp_idx = 0;
            MeshBase::const_element_iterator el = d_meter_meshes[meter_idx]->active_elements_begin();
            const MeshBase::const_element_iterator end_el = d_meter_meshes[meter_idx]->active_elements_end();
            for (; el != end_el; ++el)
            {
                const Elem* elem = *el;
                fe_elem->reinit(elem);
                const auto n_qp = q_point.size();

                meter_qp_global_ln[meter_idx].resize(
                    std::max(meter_qp_global_ln[meter_idx].size(), meter_qp_idx + n_qp), IBTK::invalid_level_number);
                meter_qp_local_ln[meter_idx].resize(std::max(meter_qp_local_ln[meter_idx].size(), meter_qp_idx + n_qp),
                                                    IBTK::invalid_level_number);
                meter_qp_local_patch_idx[meter_idx].resize(
                    std::max(meter_qp_local_patch_idx[meter_idx].size(), meter_qp_idx + n_qp), IBTK::invalid_index);

                // Loop over all quadrature points to determine if they should be "assigned" to the current level.
                for (unsigned int qp = 0; qp < n_qp; ++qp, ++meter_qp_idx)
                {
                    // Check to see if the quadrature point is in a local patch.
                    const libMesh::Point x = put_point_in_domain(q_point[qp], domain_x_lower, domain_x_upper);
                    const auto i = IndexUtilities::getCellIndex(&x(0),
                                                                domain_x_lower,
                                                                domain_x_upper,
                                                                dx.data(),
                                                                domain_box_level_lower,
                                                                domain_box_level_upper);
                    Box<NDIM> cell_box(i, i);
                    tbox::Array<int> local_patch_idxs;
                    box_tree->findLocalOverlapIndices(local_patch_idxs, cell_box);
                    if (local_patch_idxs.size() != 0)
                    {
                        TBOX_ASSERT(local_patch_idxs.size() == 1);
                        meter_qp_global_ln[meter_idx][meter_qp_idx] = ln;
                        meter_qp_local_ln[meter_idx][meter_qp_idx] = ln;
                        meter_qp_local_patch_idx[meter_idx][meter_qp_idx] = local_patch_idxs[0];
                    }
                }
            }
        }
    }

    // Determine level number assignments for all meter quadrature points and store data associated with each meter
    // quadrature point to the appropriate data structures.
    //
    // We want to assign each point to the finest available patch level.
    std::vector<std::vector<int> > meter_qp_assignment_count(d_num_meters);
    for (unsigned int meter_idx = 0; meter_idx < d_num_meters; ++meter_idx)
    {
        IBTK_MPI::maxReduction(meter_qp_global_ln[meter_idx].data(), meter_qp_global_ln[meter_idx].size());

        // Setup FE objects.
        auto& u_meter_system = d_meter_systems[meter_idx]->get_system(ib_method_ops->getVelocitySystemName());
        NumericVector<double>& u_meter_serial_vec = u_meter_system.get_vector("serial solution");
        const DofMap& u_dof_map = u_meter_system.get_dof_map();
        FEType fe_type = u_meter_system.variable_type(0);
        std::unique_ptr<FEBase> fe_elem(FEBase::build(NDIM - 1, fe_type));
        std::unique_ptr<QBase> qrule(QBase::build(d_quad_type, NDIM - 1, d_quad_order[meter_idx]));
        fe_elem->attach_quadrature_rule(qrule.get());

        const std::vector<libMesh::Point>& q_point = fe_elem->get_xyz();
        const std::vector<std::vector<double> >& phi = fe_elem->get_phi();
        const std::vector<Real>& JxW = fe_elem->get_JxW();
        std::vector<std::vector<dof_id_type> > u_dof_indices(NDIM);

        // Loop over ALL meter mesh elements.
        unsigned int meter_qp_idx = 0;
        MeshBase::const_element_iterator el = d_meter_meshes[meter_idx]->active_elements_begin();
        const MeshBase::const_element_iterator end_el = d_meter_meshes[meter_idx]->active_elements_end();
        for (; el != end_el; ++el)
        {
            const Elem* elem = *el;
            fe_elem->reinit(elem);
            const auto n_qp = q_point.size();
            for (int d = 0; d < NDIM; ++d)
            {
                u_dof_map.dof_indices(elem, u_dof_indices[d], d);
            }
            boost::multi_array<double, 2> u_node;
            get_values_for_interpolation(u_node, u_meter_serial_vec, u_dof_indices);

            // Keep track of how many patches each quadrature point is assigned to.
            //
            // (Each should be assigned to exactly one!)
            meter_qp_assignment_count[meter_idx].resize(
                std::max(meter_qp_assignment_count[meter_idx].size(), meter_qp_idx + n_qp), 0);

            // Compute the normal vector to the element.
            //
            // Because we specialize to P1 meter meshes, these are constant on each element.
            const libMesh::Point tau1 = *elem->node_ptr(1) - *elem->node_ptr(0);
            const libMesh::Point tau2 = *elem->node_ptr(2) - *elem->node_ptr(1);
            const libMesh::Point normal = tau1.cross(tau2).unit();

            // Loop over all quadrature points and set up data for the points that are assigned to a local patch.
            for (unsigned int qp = 0; qp < n_qp; ++qp, ++meter_qp_idx)
            {
                if ((meter_qp_local_ln[meter_idx][meter_qp_idx] == meter_qp_global_ln[meter_idx][meter_qp_idx]) &&
                    (meter_qp_local_patch_idx[meter_idx][meter_qp_idx] != IBTK::invalid_index))
                {
                    const int ln = meter_qp_local_ln[meter_idx][meter_qp_idx];
                    const int local_patch_idx = meter_qp_local_patch_idx[meter_idx][meter_qp_idx];

                    meter_idx_map[ln][local_patch_idx].push_back(meter_idx);
                    const libMesh::Point x = put_point_in_domain(q_point[qp], domain_x_lower, domain_x_upper);
                    meter_x_map[ln][local_patch_idx].push_back(Vector(&x(0)));
                    Vector u_corr;
                    interpolate(u_corr.data(), qp, u_node, phi);
                    meter_u_corr_map[ln][local_patch_idx].push_back(u_corr);
                    meter_normal_map[ln][local_patch_idx].push_back(Vector(&normal(0)));
                    meter_JxW_map[ln][local_patch_idx].push_back(JxW[qp]);
                    meter_qp_assignment_count[meter_idx][meter_qp_idx] += 1;
                }
            }
        }
    }

    // Confirm that each quadrature point is assigned to exactly one patch.
    for (unsigned int meter_idx = 0; meter_idx < d_num_meters; ++meter_idx)
    {
        const int num_meter_qp_points = meter_qp_assignment_count[meter_idx].size();
        TBOX_ASSERT(IBTK_MPI::minReduction(num_meter_qp_points) == IBTK_MPI::maxReduction(num_meter_qp_points));
        IBTK_MPI::sumReduction(meter_qp_assignment_count[meter_idx].data(), num_meter_qp_points);
        if ((*std::min_element(begin(meter_qp_assignment_count[meter_idx]),
                               end(meter_qp_assignment_count[meter_idx])) != 1) ||
            (*std::max_element(begin(meter_qp_assignment_count[meter_idx]),
                               end(meter_qp_assignment_count[meter_idx])) != 1))
        {
            // Print data about missing points.
            for (unsigned int meter_idx = 0; meter_idx < d_num_meters; ++meter_idx)
            {
                // Setup FE objects.
                auto& u_meter_system = d_meter_systems[meter_idx]->get_system(ib_method_ops->getVelocitySystemName());
                FEType fe_type = u_meter_system.variable_type(0);
                std::unique_ptr<FEBase> fe_elem(FEBase::build(NDIM - 1, fe_type));
                std::unique_ptr<QBase> qrule(QBase::build(d_quad_type, NDIM - 1, d_quad_order[meter_idx]));
                fe_elem->attach_quadrature_rule(qrule.get());
                const std::vector<libMesh::Point>& q_point = fe_elem->get_xyz();
                unsigned int meter_qp_idx = 0;
                MeshBase::const_element_iterator el = d_meter_meshes[meter_idx]->active_elements_begin();
                const MeshBase::const_element_iterator end_el = d_meter_meshes[meter_idx]->active_elements_end();
                for (; el != end_el; ++el)
                {
                    const Elem* elem = *el;
                    fe_elem->reinit(elem);
                    const auto n_qp = q_point.size();
                    for (unsigned int qp = 0; qp < n_qp; ++qp, ++meter_qp_idx)
                    {
                        if (meter_qp_assignment_count[meter_idx][meter_qp_idx] != 1)
                        {
                            pout << "misassigned meter qp with index=" << meter_qp_idx << " x=" << q_point[qp]
                                 << " assignment count=" << meter_qp_assignment_count[meter_idx][meter_qp_idx] << "\n";
                        }
                    }
                }
            }
            TBOX_ERROR("misassigned meter quadrature points detected\n");
        }
    }
    return;
}

void
IBFEInstrumentPanel::outputData(const double data_time)
{
    static const int mpi_root = 0;
    if (IBTK_MPI::getRank() == mpi_root)
    {
        d_flux_stream.open(d_plot_directory_name + "/flux.dat", std::ofstream::app);
        d_flux_stream.precision(15);
        d_flux_stream << data_time;

        d_mean_pressure_stream.open(d_plot_directory_name + "/mean_pressure.dat", std::ofstream::app);
        d_mean_pressure_stream.precision(15);
        d_mean_pressure_stream << data_time;

        d_centroid_pressure_stream.open(d_plot_directory_name + "/centroid_pressure.dat", std::ofstream::app);
        d_centroid_pressure_stream.precision(15);
        d_centroid_pressure_stream << data_time;

        for (unsigned int meter_idx = 0; meter_idx < d_num_meters; ++meter_idx)
        {
            d_flux_stream << " " << d_flow_rate_values[meter_idx];
            d_mean_pressure_stream << " " << d_mean_pressure_values[meter_idx];
            d_centroid_pressure_stream << " " << d_centroid_pressure_values[meter_idx];
        }

        d_flux_stream << "\n";
        d_flux_stream.close();

        d_mean_pressure_stream << "\n";
        d_mean_pressure_stream.close();

        d_centroid_pressure_stream << "\n";
        d_centroid_pressure_stream.close();
    }
}

double
IBFEInstrumentPanel::getMeterRadius(const int meter_idx)
{
    // NOTE: this function should be called **after** updating the meter system data,
    // with initializeSystemDependentData.
    return d_meter_radii[meter_idx];
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
