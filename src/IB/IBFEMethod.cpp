// Filename: IBFEMethod.cpp
// Created on 5 Oct 2011 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
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

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <ostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellIndex.h"
#include "GriddingAlgorithm.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchyDataOpsReal.h"
#include "Index.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "SideData.h"
#include "SideIndex.h"
#include "Variable.h"
#include "VariableDatabase.h"
#include "boost/multi_array.hpp"
#include "ibamr/IBFEMethod.h"
#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/FEDataInterpolation.h"
#include "ibtk/FEDataManager.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/libmesh_utilities.h"
#include "ibtk/BoxPartitioner.h"
#include "libmesh/boundary_info.h"
#include "libmesh/compare_types.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/edge.h"
#include "libmesh/elem.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_point_locator_type.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/enum_xdr_mode.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe_compute_data.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_type.h"
#include "libmesh/fem_context.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/point.h"
#include "libmesh/point_locator_base.h"
#include "libmesh/quadrature.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/system.h"
#include "libmesh/tensor_value.h"
#include "libmesh/type_tensor.h"
#include "libmesh/type_vector.h"
#include "libmesh/variant_filter_iterator.h"
#include "libmesh/vector_value.h"
#include "petscvec.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

namespace SAMRAI
{
namespace xfer
{
template <int DIM>
class RefineSchedule;
template <int DIM>
class CoarsenSchedule;
} // namespace xfer
} // namespace SAMRAI

using namespace libMesh;

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of IBFEMethod restart file data.
static const int IBFE_METHOD_VERSION = 1;

inline short int
get_dirichlet_bdry_ids(const std::vector<short int>& bdry_ids)
{
    short int dirichlet_bdry_ids = 0;
    for (const auto& bdry_id : bdry_ids)
    {
        if (bdry_id == FEDataManager::ZERO_DISPLACEMENT_X_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_Y_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_Z_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_XY_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_XZ_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_YZ_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_XYZ_BDRY_ID)
        {
            dirichlet_bdry_ids |= bdry_id;
        }
    }
    return dirichlet_bdry_ids;
}

inline bool
is_physical_bdry(const Elem* elem,
                 const unsigned short int side,
                 const BoundaryInfo& boundary_info,
                 const DofMap& dof_map)
{
    const std::vector<short int>& bdry_ids = boundary_info.boundary_ids(elem, side);
    bool at_physical_bdry = !elem->neighbor(side);
    for (const auto& bdry_id : bdry_ids)
    {
        if (dof_map.is_periodic_boundary(bdry_id)) at_physical_bdry = false;
    }
    return at_physical_bdry;
}

inline bool
is_dirichlet_bdry(const Elem* elem,
                  const unsigned short int side,
                  const BoundaryInfo& boundary_info,
                  const DofMap& dof_map)
{
    if (!is_physical_bdry(elem, side, boundary_info, dof_map)) return false;
    const std::vector<short int>& bdry_ids = boundary_info.boundary_ids(elem, side);
    return get_dirichlet_bdry_ids(bdry_ids) != 0;
}

inline bool
has_physical_bdry(const Elem* elem, const BoundaryInfo& boundary_info, const DofMap& dof_map)
{
    bool has_physical_bdry = false;
    for (unsigned short int side = 0; side < elem->n_sides() && !has_physical_bdry; ++side)
    {
        has_physical_bdry = has_physical_bdry || is_physical_bdry(elem, side, boundary_info, dof_map);
    }
    return has_physical_bdry;
}

inline void
get_x_and_FF(libMesh::VectorValue<double>& x,
             libMesh::TensorValue<double>& FF,
             const std::vector<double>& x_data,
             const std::vector<VectorValue<double> >& grad_x_data,
             const unsigned int dim = NDIM)
{
    x.zero();
    FF.zero();
    for (unsigned int i = 0; i < dim; ++i)
    {
        x(i) = x_data[i];
        for (unsigned int j = 0; j < dim; ++j)
        {
            FF(i, j) = grad_x_data[i](j);
        }
    }
    for (unsigned int i = dim; i < LIBMESH_DIM; ++i)
    {
        FF(i, i) = 1.0;
    }
    return;
}

static const Real PENALTY = 1.e10;

void
assemble_poisson(EquationSystems& es, const std::string& /*system_name*/)
{
    const MeshBase& mesh = es.get_mesh();
    const BoundaryInfo& boundary_info = *mesh.boundary_info;
    const unsigned int dim = mesh.mesh_dimension();
    auto& system = es.get_system<LinearImplicitSystem>(IBFEMethod::PHI_SYSTEM_NAME);
    const DofMap& dof_map = system.get_dof_map();
    FEType fe_type = dof_map.variable_type(0);

    std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));
    QGauss qrule(dim, FIFTH);
    fe->attach_quadrature_rule(&qrule);
    std::unique_ptr<FEBase> fe_face(FEBase::build(dim, fe_type));
    QGauss qface(dim - 1, FIFTH);
    fe_face->attach_quadrature_rule(&qface);

    const std::vector<Real>& JxW = fe->get_JxW();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

    const std::vector<std::vector<Real> >& phi_face = fe_face->get_phi();
    const std::vector<Real>& JxW_face = fe_face->get_JxW();

    DenseMatrix<Number> Ke;
    std::vector<dof_id_type> dof_indices;

    const double epsilon = es.parameters.get<Real>("Phi_epsilon");
    const double epsilon_inv = (std::abs(epsilon) > std::numeric_limits<double>::epsilon() ? 1.0 / epsilon : 0.0);

    MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    for (; el != end_el; ++el)
    {
        const Elem* elem = *el;
        dof_map.dof_indices(elem, dof_indices);
        fe->reinit(elem);
        auto Ke_size = static_cast<unsigned int>(dof_indices.size());
        Ke.resize(Ke_size, Ke_size);
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {
            for (unsigned int i = 0; i < phi.size(); i++)
            {
                for (unsigned int j = 0; j < phi.size(); j++)
                {
                    Ke(i, j) += (epsilon_inv * phi[i][qp] * phi[j][qp] + (dphi[i][qp] * dphi[j][qp])) * JxW[qp];
                }
            }
        }
        for (unsigned int side = 0; side < elem->n_sides(); side++)
        {
            if (is_physical_bdry(elem, side, boundary_info, dof_map))
            {
                fe_face->reinit(elem, side);
                for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                {
                    for (unsigned int i = 0; i < phi_face.size(); ++i)
                    {
                        for (unsigned int j = 0; j < phi_face.size(); ++j)
                        {
                            Ke(i, j) += PENALTY * phi_face[i][qp] * phi_face[j][qp] * JxW_face[qp];
                        }
                    }
                }
            }
        }
        dof_map.constrain_element_matrix(Ke, dof_indices);
        system.matrix->add_matrix(Ke, dof_indices);
    }
}

std::string
libmesh_restart_file_name(const std::string& restart_dump_dirname,
                          unsigned int time_step_number,
                          unsigned int part,
                          const std::string& extension)
{
    std::ostringstream file_name_prefix;
    file_name_prefix << restart_dump_dirname << "/libmesh_data_part_" << part << "." << std::setw(6)
                     << std::setfill('0') << std::right << time_step_number << "." << extension;
    return file_name_prefix.str();
}
} // namespace

const std::string IBFEMethod::COORDS_SYSTEM_NAME = "IB coordinates system";
const std::string IBFEMethod::COORD_MAPPING_SYSTEM_NAME = "IB coordinate mapping system";
const std::string IBFEMethod::FORCE_SYSTEM_NAME = "IB force system";
const std::string IBFEMethod::PHI_SYSTEM_NAME = "IB stress normalization system";
const std::string IBFEMethod::SOURCE_SYSTEM_NAME = "IB source system";
const std::string IBFEMethod::VELOCITY_SYSTEM_NAME = "IB velocity system";

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBFEMethod::IBFEMethod(const std::string& object_name,
                       Pointer<Database> input_db,
                       MeshBase* mesh,
                       int max_level_number,
                       bool register_for_restart,
                       const std::string& restart_read_dirname,
                       unsigned int restart_restore_number)
{
    commonConstructor(object_name,
                      input_db,
                      std::vector<MeshBase*>(1, mesh),
                      max_level_number,
                      register_for_restart,
                      restart_read_dirname,
                      restart_restore_number);
    return;
} // IBFEMethod

IBFEMethod::IBFEMethod(const std::string& object_name,
                       Pointer<Database> input_db,
                       const std::vector<MeshBase*>& meshes,
                       int max_level_number,
                       bool register_for_restart,
                       const std::string& restart_read_dirname,
                       unsigned int restart_restore_number)
    : d_num_parts(static_cast<int>(meshes.size()))
{
    commonConstructor(object_name,
                      input_db,
                      meshes,
                      max_level_number,
                      register_for_restart,
                      restart_read_dirname,
                      restart_restore_number);
    return;
} // IBFEMethod

IBFEMethod::~IBFEMethod()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }
    return;
} // ~IBFEMethod

FEDataManager*
IBFEMethod::getFEDataManager(const unsigned int part) const
{
    TBOX_ASSERT(d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_num_parts);
    return d_fe_data_managers[part];
} // getFEDataManager

void
IBFEMethod::registerStressNormalizationPart(unsigned int part)
{
    TBOX_ASSERT(d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_num_parts);
    if (d_is_stress_normalization_part[part]) return;
    d_has_stress_normalization_parts = true;
    d_is_stress_normalization_part[part] = true;
    System& Phi_system = d_equation_systems[part]->add_system<LinearImplicitSystem>(PHI_SYSTEM_NAME);
    d_equation_systems[part]->parameters.set<Real>("Phi_epsilon") = d_epsilon;
    Phi_system.attach_assemble_function(assemble_poisson);
    Phi_system.add_variable("Phi", d_fe_order[part], d_fe_family[part]);
    return;
} // registerStressNormalizationPart

void
IBFEMethod::registerInitialCoordinateMappingFunction(const CoordinateMappingFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_coordinate_mapping_fcn_data[part] = data;
    return;
} // registerInitialCoordinateMappingFunction

IBFEMethod::CoordinateMappingFcnData
IBFEMethod::getInitialCoordinateMappingFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_num_parts);
    return d_coordinate_mapping_fcn_data[part];
} // getInitialCoordinateMappingFunction

void
IBFEMethod::registerInitialVelocityFunction(const InitialVelocityFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_initial_velocity_fcn_data[part] = data;
    return;
} // registerInitialVelocityFunction

IBFEMethod::InitialVelocityFcnData
IBFEMethod::getInitialVelocityFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_num_parts);
    return d_initial_velocity_fcn_data[part];
} // getInitialVelocityFunction

void
IBFEMethod::registerPK1StressFunction(const PK1StressFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_PK1_stress_fcn_data[part].push_back(data);
    if (data.quad_type == INVALID_Q_RULE)
    {
        d_PK1_stress_fcn_data[part].back().quad_type = d_default_quad_type[part];
    }
    if (data.quad_order == INVALID_ORDER)
    {
        d_PK1_stress_fcn_data[part].back().quad_order = d_default_quad_order[part];
    }
    return;
} // registerPK1StressFunction

std::vector<IBFEMethod::PK1StressFcnData>
IBFEMethod::getPK1StressFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_num_parts);
    return d_PK1_stress_fcn_data[part];
} // getPK1StressFunction

void
IBFEMethod::registerLagBodyForceFunction(const LagBodyForceFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_lag_body_force_fcn_data[part] = data;
    return;
} // registerLagBodyForceFunction

IBFEMethod::LagBodyForceFcnData
IBFEMethod::getLagBodyForceFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_num_parts);
    return d_lag_body_force_fcn_data[part];
} // getLagBodyForceFunction

void
IBFEMethod::registerLagSurfacePressureFunction(const LagSurfacePressureFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_lag_surface_pressure_fcn_data[part] = data;
    return;
} // registerLagSurfacePressureFunction

IBFEMethod::LagSurfacePressureFcnData
IBFEMethod::getLagSurfacePressureFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_num_parts);
    return d_lag_surface_pressure_fcn_data[part];
} // getLagSurfacePressureFunction

void
IBFEMethod::registerLagSurfaceForceFunction(const LagSurfaceForceFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_lag_surface_force_fcn_data[part] = data;
    return;
} // registerLagSurfaceForceFunction

IBFEMethod::LagSurfaceForceFcnData
IBFEMethod::getLagSurfaceForceFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_num_parts);
    return d_lag_surface_force_fcn_data[part];
} // getLagSurfaceForceFunction

void
IBFEMethod::registerLagBodySourceFunction(const LagBodySourceFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    if (d_lag_body_source_part[part]) return;
    d_has_lag_body_source_parts = true;
    d_lag_body_source_part[part] = true;
    d_lag_body_source_fcn_data[part] = data;
    System& Q_system = d_equation_systems[part]->add_system<ExplicitSystem>(SOURCE_SYSTEM_NAME);
    Q_system.add_variable("Q", d_fe_order[part], d_fe_family[part]);
    return;
} // registerLagBodySourceFunction

IBFEMethod::LagBodySourceFcnData
IBFEMethod::getLagBodySourceFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_num_parts);
    return d_lag_body_source_fcn_data[part];
} // getLagBodySourceFunction

void
IBFEMethod::registerOverlappingVelocityReset(const unsigned int part1, const unsigned int part2)
{
    TBOX_ASSERT(part1 < d_num_parts);
    TBOX_ASSERT(part2 < d_num_parts);

    // Check to see if we have already registered this part.
    if (d_is_overlap_velocity_part[part1])
    {
        TBOX_ERROR("IBFEMethod::registerOverlappingVelocityReset(): already registered part " << part1 << "\n");
    }

    d_has_overlap_velocity_parts = true;
    d_is_overlap_velocity_part[part1] = true;
    d_overlap_velocity_master_part[part1] = part2;
    d_is_overlap_velocity_master_part[part2] = true;

    std::array<unsigned int, 2> part_idx;
    part_idx[0] = part1;
    part_idx[1] = part2;

    // Set up basic data structures.
    std::array<EquationSystems*, 2> es;
    std::array<MeshBase*, 2> mesh;
    std::array<const DofMap*, 2> U_dof_map, X_dof_map;
    std::array<std::unique_ptr<PointLocatorBase>, 2> ploc;
    std::array<numeric_index_type, 2> first_local_idx, last_local_idx;
    for (int k = 0; k < 2; ++k)
    {
        es[k] = d_fe_data_managers[part_idx[k]]->getEquationSystems();
        auto& U_system = es[k]->get_system<System>(VELOCITY_SYSTEM_NAME);
        U_dof_map[k] = &U_system.get_dof_map();
        NumericVector<double>& U_vec = *U_system.solution;
        auto& X_system = es[k]->get_system<System>(COORDS_SYSTEM_NAME);
        X_dof_map[k] = &X_system.get_dof_map();
        NumericVector<double>& X_vec = *X_system.solution;
        TBOX_ASSERT(U_vec.first_local_index() == X_vec.first_local_index());
        TBOX_ASSERT(U_vec.last_local_index() == X_vec.last_local_index());
        first_local_idx[k] = U_vec.first_local_index();
        last_local_idx[k] = U_vec.last_local_index();
        FEType fe_type = U_dof_map[k]->variable_type(0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(fe_type == U_dof_map[k]->variable_type(d));
            TBOX_ASSERT(fe_type == X_dof_map[k]->variable_type(d));
        }
        mesh[k] = &es[k]->get_mesh();
        ploc[k] = PointLocatorBase::build(TREE_ELEMENTS, *mesh[k]);
        ploc[k]->enable_out_of_mesh_mode();
        if (d_overlap_tolerance > 0.0) ploc[k]->set_close_to_point_tol(d_overlap_tolerance);
    }

    // Find the elements of part2 that contain nodes of part1.
    std::array<std::vector<std::vector<unsigned int> >, 2> U_dof_indices;
    for (int k = 0; k < 2; ++k)
    {
        U_dof_indices[k].resize(NDIM);
    }
    const unsigned int k_slave = 0;
    const unsigned int k_master = 1;
    std::map<dof_id_type, dof_id_type>& node_to_elem_map = d_overlap_velocity_part_node_to_elem_map[part_idx[k_slave]];
    MeshBase::node_iterator n = mesh[k_slave]->local_nodes_begin();
    const MeshBase::node_iterator end_n = mesh[k_slave]->local_nodes_end();
    for (; n != end_n; ++n)
    {
        const Node* const node = *n;
        const libMesh::Point& X = *node;
        const Elem* const other_elem = (*ploc[k_master])(X);
        if (other_elem)
        {
            node_to_elem_map[node->id()] = other_elem->id();
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                U_dof_map[k_master]->dof_indices(other_elem, U_dof_indices[k_master][d], d);
                for (const auto& idx : U_dof_indices[k_master][d])
                {
                    if (idx < first_local_idx[k_master] || idx >= last_local_idx[k_master])
                    {
                        d_overlap_velocity_part_ghost_idxs[part_idx[k_master]].insert(idx);
                    }
                }
            }
        }
    }
    return;
}

void
IBFEMethod::registerOverlappingForceConstraint(const unsigned int part1,
                                               const unsigned int part2,
                                               const double kappa,
                                               QBase* qrule1,
                                               QBase* qrule2)
{
    TBOX_ASSERT(part1 < d_num_parts);
    TBOX_ASSERT(part2 < d_num_parts);

    d_has_overlap_force_parts = true;
    d_is_overlap_force_part[part1] = true;
    d_is_overlap_force_part[part2] = true;

    // Double-check that this pair of parts hasn't already been registered.
    unsigned int n_pairs = d_overlap_force_part_idxs.size();
    for (unsigned int k = 0; k < n_pairs; ++k)
    {
        if ((d_overlap_force_part_idxs[k][0] == part1 && d_overlap_force_part_idxs[k][1] == part2) ||
            (d_overlap_force_part_idxs[k][0] == part2 && d_overlap_force_part_idxs[k][1] == part1))
            return;
    }

    // Set up basic data structures.
    n_pairs += 1;
    d_overlap_force_part_idxs.resize(n_pairs);
    std::array<unsigned int, 2>& part_idx = d_overlap_force_part_idxs.back();
    part_idx[0] = part1;
    part_idx[1] = part2;
    d_overlapping_elem_map.resize(n_pairs);
    std::array<std::map<libMesh::dof_id_type, std::map<unsigned int, libMesh::dof_id_type> >, 2>& elem_map =
        d_overlapping_elem_map.back();
    d_overlap_force_part_kappa.push_back(kappa);
    d_overlap_force_part_qrule.resize(n_pairs);
    std::array<QBase*, 2>& qrule = d_overlap_force_part_qrule.back();
    qrule[0] = qrule1;
    qrule[1] = qrule2;

    std::array<EquationSystems*, 2> es;
    std::array<MeshBase*, 2> mesh;
    std::array<const DofMap*, 2> F_dof_map, X_dof_map;
    std::array<std::unique_ptr<FEBase>, 2> fe;
    std::array<std::unique_ptr<PointLocatorBase>, 2> ploc;
    std::array<numeric_index_type, 2> first_local_idx, last_local_idx;
    for (int k = 0; k < 2; ++k)
    {
        es[k] = d_fe_data_managers[part_idx[k]]->getEquationSystems();
        auto& F_system = es[k]->get_system<System>(FORCE_SYSTEM_NAME);
        auto& X_system = es[k]->get_system<System>(COORDS_SYSTEM_NAME);
        F_dof_map[k] = &F_system.get_dof_map();
        X_dof_map[k] = &X_system.get_dof_map();
        NumericVector<double>& X_vec = *X_system.solution;
        first_local_idx[k] = X_vec.first_local_index();
        last_local_idx[k] = X_vec.last_local_index();
        FEType fe_type = F_dof_map[k]->variable_type(0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(fe_type == F_dof_map[k]->variable_type(d));
            TBOX_ASSERT(fe_type == X_dof_map[k]->variable_type(d));
        }
        mesh[k] = &es[k]->get_mesh();
        if (!qrule[k])
        {
            // \todo try to fix this when we update to C++11!
            qrule[k] = QBase::build(QGRID, mesh[k]->mesh_dimension(), FOURTH).release();
        }
        fe[k] = FEBase::build(mesh[k]->mesh_dimension(), fe_type);
        fe[k]->attach_quadrature_rule(qrule[k]);
        fe[k]->get_xyz();
        ploc[k] = PointLocatorBase::build(TREE_ELEMENTS, *mesh[k]);
        ploc[k]->enable_out_of_mesh_mode();
        if (d_overlap_tolerance > 0.0) ploc[k]->set_close_to_point_tol(d_overlap_tolerance);
    }

    // Find quadrature points that overlap with the other part.
    std::array<std::vector<std::vector<unsigned int> >, 2> X_dof_indices;
    for (int k = 0; k < 2; ++k)
    {
        X_dof_indices[k].resize(NDIM);
    }
    for (int k = 0; k < 2; ++k)
    {
        const int k_next = (k + 1) % 2;
        const std::vector<libMesh::Point>& q_point = fe[k]->get_xyz();
        MeshBase::const_element_iterator el = mesh[k]->active_local_elements_begin();
        const MeshBase::const_element_iterator end_el = mesh[k]->active_local_elements_end();
        for (; el != end_el; ++el)
        {
            const Elem* const elem = *el;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_dof_map[k]->dof_indices(elem, X_dof_indices[k][d], d);
                for (const auto& idx : X_dof_indices[k][d])
                {
                    if (idx < first_local_idx[k] || idx >= last_local_idx[k])
                    {
                        d_overlap_force_part_ghost_idxs[part_idx[k]].insert(idx);
                    }
                }
            }
            fe[k]->reinit(elem);
            for (unsigned int qp = 0; qp < qrule[k]->n_points(); qp++)
            {
                const Elem* const other_elem = (*ploc[k_next])(q_point[qp]);
                if (other_elem)
                {
                    elem_map[k][elem->id()][qp] = other_elem->id();
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        X_dof_map[k_next]->dof_indices(other_elem, X_dof_indices[k_next][d], d);
                        for (const auto& idx : X_dof_indices[k_next][d])
                        {
                            if (idx < first_local_idx[k_next] || idx >= last_local_idx[k_next])
                            {
                                d_overlap_force_part_ghost_idxs[part_idx[k_next]].insert(idx);
                            }
                        }
                    }
                }
            }
        }
    }
}

void
IBFEMethod::registerDirectForcingKinematics(Pointer<IBFEDirectForcingKinematics> data, unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_direct_forcing_kinematics_data[part] = data;
    return;
} // registerDirectForcingKinematics


const IntVector<NDIM>&
IBFEMethod::getMinimumGhostCellWidth() const
{
    return d_ghosts;
} // getMinimumGhostCellWidth

void
IBFEMethod::setupTagBuffer(Array<int>& tag_buffer, Pointer<GriddingAlgorithm<NDIM> > gridding_alg) const
{
    const int finest_hier_ln = gridding_alg->getMaxLevels() - 1;
    const int tsize = tag_buffer.size();
    tag_buffer.resizeArray(finest_hier_ln);
    for (int i = tsize; i < finest_hier_ln; ++i) tag_buffer[i] = 0;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        const int gcw = d_fe_data_managers[part]->getGhostCellWidth().max();
        const int tag_ln = d_fe_data_managers[part]->getLevelNumber() - 1;
        if (tag_ln >= 0 && tag_ln < finest_hier_ln)
        {
            tag_buffer[tag_ln] = std::max(tag_buffer[tag_ln], gcw);
        }
    }
    for (int ln = finest_hier_ln - 2; ln >= 0; --ln)
    {
        tag_buffer[ln] =
            std::max(tag_buffer[ln], tag_buffer[ln + 1] / gridding_alg->getRatioToCoarserLevel(ln + 1).max() + 1);
    }
    return;
} // setupTagBuffer

void
IBFEMethod::preprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    d_current_time = current_time;
    d_new_time = new_time;
    d_half_time = current_time + 0.5 * (new_time - current_time);

    // Extract the FE data.
    d_X_systems.resize(d_num_parts);
    d_X_current_vecs.resize(d_num_parts);
    d_X_new_vecs.resize(d_num_parts);
    d_X_half_vecs.resize(d_num_parts);
    d_X_IB_ghost_vecs.resize(d_num_parts);

    d_U_systems.resize(d_num_parts);
    d_U_current_vecs.resize(d_num_parts);
    d_U_new_vecs.resize(d_num_parts);
    d_U_half_vecs.resize(d_num_parts);

    d_F_systems.resize(d_num_parts);
    d_F_half_vecs.resize(d_num_parts);
    d_F_IB_ghost_vecs.resize(d_num_parts);

    d_Q_systems.resize(d_num_parts);
    d_Q_half_vecs.resize(d_num_parts);
    d_Q_IB_ghost_vecs.resize(d_num_parts);

    d_Phi_systems.resize(d_num_parts);
    d_Phi_half_vecs.resize(d_num_parts);

    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_X_systems[part] = &d_equation_systems[part]->get_system(COORDS_SYSTEM_NAME);
        d_X_current_vecs[part] = dynamic_cast<PetscVector<double>*>(d_X_systems[part]->current_local_solution.get());
        d_X_new_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_X_current_vecs[part]->clone().release()); // WARNING: must be manually deleted
        d_X_half_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_X_current_vecs[part]->clone().release()); // WARNING: must be manually deleted
        d_X_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_fe_data_managers[part]->buildGhostedCoordsVector(/*localize_data*/ false));

        d_U_systems[part] = &d_equation_systems[part]->get_system(VELOCITY_SYSTEM_NAME);
        d_U_current_vecs[part] = dynamic_cast<PetscVector<double>*>(d_U_systems[part]->current_local_solution.get());
        d_U_new_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_U_current_vecs[part]->clone().release()); // WARNING: must be manually deleted
        d_U_half_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_U_current_vecs[part]->clone().release()); // WARNING: must be manually deleted

        d_F_systems[part] = &d_equation_systems[part]->get_system(FORCE_SYSTEM_NAME);
        d_F_half_vecs[part] = dynamic_cast<PetscVector<double>*>(d_F_systems[part]->current_local_solution.get());
        d_F_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_fe_data_managers[part]->buildGhostedSolutionVector(FORCE_SYSTEM_NAME, /*localize_data*/ false));

        if (d_lag_body_source_part[part])
        {
            d_Q_systems[part] = &d_equation_systems[part]->get_system(SOURCE_SYSTEM_NAME);
            d_Q_half_vecs[part] = dynamic_cast<PetscVector<double>*>(d_Q_systems[part]->current_local_solution.get());
            d_Q_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
                d_fe_data_managers[part]->buildGhostedSolutionVector(SOURCE_SYSTEM_NAME, /*localize_data*/ false));
        }

        if (d_is_stress_normalization_part[part])
        {
            d_Phi_systems[part] = &d_equation_systems[part]->get_system(PHI_SYSTEM_NAME);
            d_Phi_half_vecs[part] =
                dynamic_cast<PetscVector<double>*>(d_Phi_systems[part]->current_local_solution.get());
        }

        // Initialize X^{n+1/2} and X^{n+1} to equal X^{n}, and initialize
        // U^{n+1/2} and U^{n+1} to equal U^{n}.
        copy_and_synch(*d_X_systems[part]->solution, *d_X_current_vecs[part]);
        *d_X_new_vecs[part] = *d_X_current_vecs[part];
        *d_X_half_vecs[part] = *d_X_current_vecs[part];

        copy_and_synch(*d_U_systems[part]->solution, *d_U_current_vecs[part]);
        *d_U_new_vecs[part] = *d_U_current_vecs[part];
        *d_U_half_vecs[part] = *d_U_current_vecs[part];

        copy_and_synch(*d_F_systems[part]->solution, *d_F_half_vecs[part]);

        if (d_lag_body_source_part[part])
        {
            copy_and_synch(*d_Q_systems[part]->solution, *d_Q_half_vecs[part]);
        }

        if (d_is_stress_normalization_part[part])
        {
            copy_and_synch(*d_Phi_systems[part]->solution, *d_Phi_half_vecs[part]);
        }

        if (d_direct_forcing_kinematics_data[part])
        {
            d_direct_forcing_kinematics_data[part]->preprocessIntegrateData(current_time, new_time, num_cycles);
        }
    }

    // Update the mask data.
    getVelocityHierarchyDataOps()->copyData(mask_new_idx, mask_current_idx);
    return;
} // preprocessIntegrateData

void
IBFEMethod::postprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    for (unsigned part = 0; part < d_num_parts; ++part)
    {
        // Reset time-dependent Lagrangian data.
        copy_and_synch(*d_X_new_vecs[part], *d_X_systems[part]->solution);
        *d_X_systems[part]->current_local_solution = *d_X_new_vecs[part];
        delete d_X_new_vecs[part];
        delete d_X_half_vecs[part];

        copy_and_synch(*d_U_new_vecs[part], *d_U_systems[part]->solution);
        *d_U_systems[part]->current_local_solution = *d_U_new_vecs[part];
        delete d_U_new_vecs[part];
        delete d_U_half_vecs[part];

        copy_and_synch(*d_F_half_vecs[part], *d_F_systems[part]->solution);
        *d_F_systems[part]->current_local_solution = *d_F_half_vecs[part];

        if (d_lag_body_source_part[part])
        {
            copy_and_synch(*d_Q_half_vecs[part], *d_Q_systems[part]->solution);
            *d_Q_systems[part]->current_local_solution = *d_Q_half_vecs[part];
        }

        if (d_is_stress_normalization_part[part])
        {
            copy_and_synch(*d_Phi_half_vecs[part], *d_Phi_systems[part]->solution);
            *d_Phi_systems[part]->current_local_solution = *d_Phi_half_vecs[part];
        }

        if (d_direct_forcing_kinematics_data[part])
        {
            d_direct_forcing_kinematics_data[part]->postprocessIntegrateData(current_time, new_time, num_cycles);
        }

        // Update the coordinate mapping dX = X - s.
        updateCoordinateMapping(part);
    }

    d_X_systems.clear();
    d_X_current_vecs.clear();
    d_X_new_vecs.clear();
    d_X_half_vecs.clear();
    d_X_IB_ghost_vecs.clear();

    d_U_systems.clear();
    d_U_current_vecs.clear();
    d_U_new_vecs.clear();
    d_U_half_vecs.clear();

    d_F_systems.clear();
    d_F_half_vecs.clear();
    d_F_IB_ghost_vecs.clear();

    d_Q_systems.clear();
    d_Q_half_vecs.clear();
    d_Q_IB_ghost_vecs.clear();

    d_Phi_systems.clear();
    d_Phi_half_vecs.clear();

    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time = std::numeric_limits<double>::quiet_NaN();
    d_half_time = std::numeric_limits<double>::quiet_NaN();
    return;
} // postprocessIntegrateData

void
IBFEMethod::interpolateVelocity(const int u_data_idx,
                                const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*u_synch_scheds*/,
                                const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                const double data_time)
{
    std::vector<NumericVector<double>*> U_vecs(d_num_parts), X_vecs(d_num_parts);
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        if (MathUtilities<double>::equalEps(data_time, d_current_time))
        {
            U_vecs[part] = d_U_current_vecs[part];
            X_vecs[part] = d_X_current_vecs[part];
        }
        else if (MathUtilities<double>::equalEps(data_time, d_half_time))
        {
            U_vecs[part] = d_U_half_vecs[part];
            X_vecs[part] = d_X_half_vecs[part];
        }
        else if (MathUtilities<double>::equalEps(data_time, d_new_time))
        {
            U_vecs[part] = d_U_new_vecs[part];
            X_vecs[part] = d_X_new_vecs[part];
        }
    }

    // Communicate ghost data.
    for (const auto& u_ghost_fill_sched : u_ghost_fill_scheds)
    {
        if (u_ghost_fill_sched) u_ghost_fill_sched->fillData(data_time);
    }
    std::vector<Pointer<RefineSchedule<NDIM> > > no_fill(u_ghost_fill_scheds.size(), nullptr);

    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        copy_and_synch(*X_vecs[part], *d_X_IB_ghost_vecs[part], /*close_v_in*/ false);
    }

    // Build the right-hand-sides to compute the interpolated data.
    std::vector<PetscVector<double>*> U_rhs_vecs(d_num_parts);
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        U_rhs_vecs[part] = static_cast<PetscVector<double>*>(U_vecs[part]->zero_clone().release());
        d_fe_data_managers[part]->interpWeighted(u_data_idx,
                                                 *U_rhs_vecs[part],
                                                 *d_X_IB_ghost_vecs[part],
                                                 VELOCITY_SYSTEM_NAME,
                                                 no_fill,
                                                 data_time,
                                                 /*close_F*/ false,
                                                 /*close_X*/ false);
        int ierr = VecAssemblyBegin(U_rhs_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
    }

    // Solve for the interpolated data.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        int ierr = VecAssemblyEnd(U_rhs_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        d_fe_data_managers[part]->computeL2Projection(*U_vecs[part],
                                                      *U_rhs_vecs[part],
                                                      VELOCITY_SYSTEM_NAME,
                                                      d_interp_spec[part].use_consistent_mass_matrix,
                                                      /*close_U*/ true,
                                                      /*close_F*/ false);
        delete U_rhs_vecs[part];
    }

    // Account for any velocity constraints.
    if (d_has_overlap_velocity_parts)
    {
        resetOverlapNodalValues(VELOCITY_SYSTEM_NAME, U_vecs);
    }
    return;
} // interpolateVelocity

void
IBFEMethod::forwardEulerStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    int ierr;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        ierr = VecWAXPY(d_X_new_vecs[part]->vec(), dt, d_U_current_vecs[part]->vec(), d_X_current_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(
            d_X_half_vecs[part]->vec(), 0.5, 0.5, 0.0, d_X_current_vecs[part]->vec(), d_X_new_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        d_X_new_vecs[part]->close();
        d_X_half_vecs[part]->close();

        if (d_direct_forcing_kinematics_data[part])
        {
            d_direct_forcing_kinematics_data[part]->forwardEulerStep(
                current_time, new_time, *d_X_current_vecs[part], *d_X_half_vecs[part], *d_X_new_vecs[part]);
        }
    }

    // Account for any velocity constraints.
    if (d_has_overlap_velocity_parts)
    {
        resetOverlapNodalValues(COORDS_SYSTEM_NAME, d_X_half_vecs);
        resetOverlapNodalValues(COORDS_SYSTEM_NAME, d_X_new_vecs);
    }
    return;
} // eulerStep

void
IBFEMethod::midpointStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    int ierr;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        ierr = VecWAXPY(d_X_new_vecs[part]->vec(), dt, d_U_half_vecs[part]->vec(), d_X_current_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(
            d_X_half_vecs[part]->vec(), 0.5, 0.5, 0.0, d_X_current_vecs[part]->vec(), d_X_new_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        d_X_new_vecs[part]->close();
        d_X_half_vecs[part]->close();
        if (d_direct_forcing_kinematics_data[part])
        {
            d_direct_forcing_kinematics_data[part]->midpointStep(
                current_time, new_time, *d_X_current_vecs[part], *d_X_half_vecs[part], *d_X_new_vecs[part]);
        }
    }

    // Account for any velocity constraints.
    if (d_has_overlap_velocity_parts)
    {
        resetOverlapNodalValues(COORDS_SYSTEM_NAME, d_X_half_vecs);
        resetOverlapNodalValues(COORDS_SYSTEM_NAME, d_X_new_vecs);
    }
    return;
} // midpointStep

void
IBFEMethod::trapezoidalStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    int ierr;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        ierr =
            VecWAXPY(d_X_new_vecs[part]->vec(), 0.5 * dt, d_U_current_vecs[part]->vec(), d_X_current_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPY(d_X_new_vecs[part]->vec(), 0.5 * dt, d_U_new_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(
            d_X_half_vecs[part]->vec(), 0.5, 0.5, 0.0, d_X_current_vecs[part]->vec(), d_X_new_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        d_X_new_vecs[part]->close();
        d_X_half_vecs[part]->close();
        if (d_direct_forcing_kinematics_data[part])
        {
            d_direct_forcing_kinematics_data[part]->trapezoidalStep(
                current_time, new_time, *d_X_current_vecs[part], *d_X_half_vecs[part], *d_X_new_vecs[part]);
        }
    }

    // Account for any velocity constraints.
    if (d_has_overlap_velocity_parts)
    {
        resetOverlapNodalValues(COORDS_SYSTEM_NAME, d_X_half_vecs);
        resetOverlapNodalValues(COORDS_SYSTEM_NAME, d_X_new_vecs);
    }
    return;
} // trapezoidalStep

void
IBFEMethod::computeLagrangianForce(const double data_time)
{
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
    for (unsigned part = 0; part < d_num_parts; ++part)
    {
        if (d_is_stress_normalization_part[part])
        {
            computeStressNormalization(*d_Phi_half_vecs[part], *d_X_half_vecs[part], data_time, part);
        }
        // Zero-out force vector
        d_F_half_vecs[part]->zero();
        computeInteriorForceDensity(*d_F_half_vecs[part], *d_X_half_vecs[part], d_Phi_half_vecs[part], data_time, part);
        if (d_direct_forcing_kinematics_data[part])
        {
            int ierr;
            PetscVector<double>* F_half_df = static_cast<PetscVector<double>*>(
                d_F_half_vecs[part]->zero_clone().release()); // WARNING: must be manually deleted
            d_direct_forcing_kinematics_data[part]->computeLagrangianForce(
                *F_half_df, *d_X_half_vecs[part], *d_U_half_vecs[part], data_time);
            ierr = VecAXPY(d_F_half_vecs[part]->vec(), 1.0, F_half_df->vec());
            IBTK_CHKERRQ(ierr);
            delete F_half_df;
        }
    }
    if (d_has_overlap_force_parts)
    {
        std::vector<std::unique_ptr<libMesh::PetscVector<double>>> X_half_ghost_vecs(d_X_half_vecs.size());
        for (unsigned int k = 0; k < d_X_half_vecs.size(); ++k)
        {
            if (d_is_overlap_force_part[k])
            {
                PetscVector<double>* X_half_vec = d_X_half_vecs[k];
                std::vector<numeric_index_type> ghost_idxs(d_overlap_force_part_ghost_idxs[k].begin(),
                                                           d_overlap_force_part_ghost_idxs[k].end());
                X_half_ghost_vecs[k] = std::unique_ptr<libMesh::PetscVector<double>>(
                    new libMesh::PetscVector<double>(
                    X_half_vec->comm(), X_half_vec->size(), X_half_vec->local_size(), ghost_idxs));
                copy_and_synch(*X_half_vec, *X_half_ghost_vecs[k], /*close_v_in*/ false);
            }
        }
        std::vector<libMesh::PetscVector<double> *> vec_pointers;
        for (std::unique_ptr<libMesh::PetscVector<double>> &vec : X_half_ghost_vecs)
        {
            vec_pointers.push_back(vec.get());
        }
        computeOverlapConstraintForceDensity(d_F_half_vecs, vec_pointers);
    }
    return;
} // computeLagrangianForce

void
IBFEMethod::spreadForce(const int f_data_idx,
                        RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                        const std::vector<Pointer<RefineSchedule<NDIM> > >& /*f_prolongation_scheds*/,
                        const double data_time)
{
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));

    // Communicate ghost data.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        PetscVector<double>* X_vec = d_X_half_vecs[part];
        PetscVector<double>* X_ghost_vec = d_X_IB_ghost_vecs[part];
        PetscVector<double>* F_vec = d_F_half_vecs[part];
        PetscVector<double>* F_ghost_vec = d_F_IB_ghost_vecs[part];
        copy_and_synch(*X_vec, *X_ghost_vec, /*close_v_in*/ false);
        copy_and_synch(*F_vec, *F_ghost_vec, /*close_v_in*/ false);
    }

    // Spread interior force density values.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        PetscVector<double>* X_ghost_vec = d_X_IB_ghost_vecs[part];
        PetscVector<double>* F_ghost_vec = d_F_IB_ghost_vecs[part];
        d_fe_data_managers[part]->spread(f_data_idx,
                                         *F_ghost_vec,
                                         *X_ghost_vec,
                                         FORCE_SYSTEM_NAME,
                                         f_phys_bdry_op,
                                         data_time,
                                         /*close_F*/ false,
                                         /*close_X*/ false);
    }

    // Handle any transmission conditions.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        PetscVector<double>* X_ghost_vec = d_X_IB_ghost_vecs[part];
        PetscVector<double>* F_ghost_vec = d_F_IB_ghost_vecs[part];
        if (d_split_normal_force || d_split_tangential_force)
        {
            if (d_use_jump_conditions && d_split_normal_force)
            {
                imposeJumpConditions(f_data_idx, *F_ghost_vec, *X_ghost_vec, data_time, part);
            }
            if (!d_use_jump_conditions || d_split_tangential_force)
            {
                spreadTransmissionForceDensity(f_data_idx, *X_ghost_vec, f_phys_bdry_op, data_time, part);
            }
        }
    }
    return;
} // spreadForce

bool
IBFEMethod::hasFluidSources() const
{
    return d_has_lag_body_source_parts;
}

void
IBFEMethod::computeLagrangianFluidSource(double data_time)
{
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        if (!d_lag_body_source_part[part]) continue;

        EquationSystems& equation_systems = *d_fe_data_managers[part]->getEquationSystems();
        const MeshBase& mesh = equation_systems.get_mesh();
        const unsigned int dim = mesh.mesh_dimension();

        // Extract the FE systems and DOF maps, and setup the FE object.
        auto& Q_system = equation_systems.get_system<ExplicitSystem>(SOURCE_SYSTEM_NAME);
        const DofMap& Q_dof_map = Q_system.get_dof_map();
        FEDataManager::SystemDofMapCache& Q_dof_map_cache =
            *d_fe_data_managers[part]->getDofMapCache(SOURCE_SYSTEM_NAME);
        FEType Q_fe_type = Q_dof_map.variable_type(0);
        std::vector<unsigned int> Q_dof_indices;
        System& X_system = equation_systems.get_system(COORDS_SYSTEM_NAME);
        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;

        FEDataInterpolation fe(dim, d_fe_data_managers[part]);
        std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, dim, FIFTH);
        fe.attachQuadratureRule(qrule.get());
        fe.evalQuadraturePoints();
        fe.evalQuadratureWeights();
        fe.registerSystem(Q_system);
        NumericVector<double>& X_vec = *d_X_half_vecs[part];
        const size_t X_sys_idx = fe.registerInterpolatedSystem(X_system, vars, vars, &X_vec);
        std::vector<size_t> Q_fcn_system_idxs;
        fe.setupInterpolatedSystemDataIndexes(
            Q_fcn_system_idxs, d_lag_body_source_fcn_data[part].system_data, &equation_systems);
        fe.init(/*use_IB_ghosted_vecs*/ false);

        const std::vector<libMesh::Point>& q_point = fe.getQuadraturePoints();
        const std::vector<double>& JxW = fe.getQuadratureWeights();
        const std::vector<std::vector<double> >& phi = fe.getPhi(Q_fe_type);

        const std::vector<std::vector<std::vector<double> > >& fe_interp_var_data = fe.getVarInterpolation();
        const std::vector<std::vector<std::vector<VectorValue<double> > > >& fe_interp_grad_var_data =
            fe.getGradVarInterpolation();

        std::vector<const std::vector<double>*> Q_var_data;
        std::vector<const std::vector<VectorValue<double> >*> Q_grad_var_data;

        // Setup global and elemental right-hand-side vectors.
        NumericVector<double>* Q_rhs_vec = Q_system.rhs;
        Q_rhs_vec->zero();
        DenseVector<double> Q_rhs_e;

        TensorValue<double> FF, FF_inv_trans;
        VectorValue<double> x;
        double Q;
        const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            Q_dof_map_cache.dof_indices(elem, Q_dof_indices, 0);
            Q_rhs_e.resize(static_cast<int>(Q_dof_indices.size()));
            fe.reinit(elem);
            fe.collectDataForInterpolation(elem);
            fe.interpolate(elem);
            const unsigned int n_qp = qrule->n_points();
            const size_t n_basis = phi.size();
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                const libMesh::Point& X = q_point[qp];
                const std::vector<double>& x_data = fe_interp_var_data[qp][X_sys_idx];
                const std::vector<VectorValue<double> >& grad_x_data = fe_interp_grad_var_data[qp][X_sys_idx];
                get_x_and_FF(x, FF, x_data, grad_x_data);

                fe.setInterpolatedDataPointers(Q_var_data, Q_grad_var_data, Q_fcn_system_idxs, elem, qp);
                d_lag_body_source_fcn_data[part].fcn(
                    Q, FF, x, X, elem, Q_var_data, Q_grad_var_data, data_time, d_lag_body_source_fcn_data[part].ctx);
                for (unsigned int k = 0; k < n_basis; ++k)
                {
                    Q_rhs_e(k) += Q * phi[k][qp] * JxW[qp];
                }
            }

            // Apply constraints (e.g., enforce periodic boundary conditions)
            // and add the elemental contributions to the global vector.
            Q_dof_map.constrain_element_vector(Q_rhs_e, Q_dof_indices);
            Q_rhs_vec->add_vector(Q_rhs_e, Q_dof_indices);
        }

        // Solve for Q.
        NumericVector<double>& Q_vec = *d_Q_half_vecs[part];
        d_fe_data_managers[part]->computeL2Projection(
            Q_vec, *Q_rhs_vec, SOURCE_SYSTEM_NAME, d_use_consistent_mass_matrix);
    }
}

void
IBFEMethod::spreadFluidSource(const int q_data_idx,
                              RobinPhysBdryPatchStrategy* q_phys_bdry_op,
                              const std::vector<Pointer<RefineSchedule<NDIM> > >& /*q_prolongation_scheds*/,
                              const double data_time)
{
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        if (!d_lag_body_source_part[part]) continue;
        PetscVector<double>* X_vec = d_X_half_vecs[part];
        PetscVector<double>* X_ghost_vec = d_X_IB_ghost_vecs[part];
        PetscVector<double>* Q_vec = d_Q_half_vecs[part];
        PetscVector<double>* Q_ghost_vec = d_Q_IB_ghost_vecs[part];
        copy_and_synch(*X_vec, *X_ghost_vec, /*close_v_in*/ false);
        copy_and_synch(*Q_vec, *Q_ghost_vec, /*close_v_in*/ false);
        d_fe_data_managers[part]->spread(q_data_idx,
                                         *Q_ghost_vec,
                                         *X_ghost_vec,
                                         SOURCE_SYSTEM_NAME,
                                         q_phys_bdry_op,
                                         data_time,
                                         /*close_Q*/ false,
                                         /*close_X*/ false);
    }
    return;
}

FEDataManager::InterpSpec
IBFEMethod::getDefaultInterpSpec() const
{
    return d_default_interp_spec;
}

FEDataManager::SpreadSpec
IBFEMethod::getDefaultSpreadSpec() const
{
    return d_default_spread_spec;
}

void
IBFEMethod::setInterpSpec(const FEDataManager::InterpSpec& interp_spec, const unsigned int part)
{
    TBOX_ASSERT(!d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_num_parts);
    d_interp_spec[part] = interp_spec;
    return;
}

void
IBFEMethod::setSpreadSpec(const FEDataManager::SpreadSpec& spread_spec, const unsigned int part)
{
    TBOX_ASSERT(!d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_num_parts);
    d_spread_spec[part] = spread_spec;
    return;
}

void
IBFEMethod::initializeFEEquationSystems()
{
    if (d_fe_equation_systems_initialized) return;

    const bool from_restart = RestartManager::getManager()->isFromRestart();

    // Create the FE data managers that manage mappings between the FE mesh
    // parts and the Cartesian grid.
    d_equation_systems.resize(d_num_parts);
    d_fe_data_managers.resize(d_num_parts, nullptr);
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        // Create FE data managers.
        std::ostringstream manager_stream;
        manager_stream << "IBFEMethod FEDataManager::" << part;
        const std::string& manager_name = manager_stream.str();
        d_fe_data_managers[part] = FEDataManager::getManager(manager_name, d_interp_spec[part], d_spread_spec[part]);
        if (d_workload_idx != IBTK::invalid_index)
        {
            d_fe_data_managers[part]->registerLoadBalancer(d_load_balancer, d_workload_idx);
        }

        d_ghosts = IntVector<NDIM>::max(d_ghosts, d_fe_data_managers[part]->getGhostCellWidth());

        // Create FE equation systems objects and corresponding variables.
        d_equation_systems[part] = std::unique_ptr<EquationSystems>(new EquationSystems(*d_meshes[part]));
        EquationSystems& equation_systems = *d_equation_systems[part];
        d_fe_data_managers[part]->setEquationSystems(&equation_systems, d_max_level_number - 1);
        d_fe_data_managers[part]->COORDINATES_SYSTEM_NAME = COORDS_SYSTEM_NAME;
        if (from_restart)
        {
            const std::string& file_name = libmesh_restart_file_name(
                d_libmesh_restart_read_dir, d_libmesh_restart_restore_number, part, d_libmesh_restart_file_extension);
            const XdrMODE xdr_mode = (d_libmesh_restart_file_extension == "xdr" ? DECODE : READ);
            const int read_mode =
                EquationSystems::READ_HEADER | EquationSystems::READ_DATA | EquationSystems::READ_ADDITIONAL_DATA;
            equation_systems.read(file_name, xdr_mode, read_mode, /*partition_agnostic*/ true);
        }
        else
        {
            auto& X_system = equation_systems.add_system<System>(COORDS_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream os;
                os << "X_" << d;
                X_system.add_variable(os.str(), d_fe_order[part], d_fe_family[part]);
            }

            auto& dX_system = equation_systems.add_system<System>(COORD_MAPPING_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream os;
                os << "dX_" << d;
                dX_system.add_variable(os.str(), d_fe_order[part], d_fe_family[part]);
            }

            auto& U_system = equation_systems.add_system<System>(VELOCITY_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream os;
                os << "U_" << d;
                U_system.add_variable(os.str(), d_fe_order[part], d_fe_family[part]);
            }

            auto& F_system = equation_systems.add_system<System>(FORCE_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream os;
                os << "F_" << d;
                F_system.add_variable(os.str(), d_fe_order[part], d_fe_family[part]);
            }
        }
    }
    d_fe_equation_systems_initialized = true;
    return;
}

void
IBFEMethod::initializeFEData()
{
    if (d_fe_data_initialized) return;

    initializeFEEquationSystems();
    doInitializeFEData(RestartManager::getManager()->isFromRestart());
    d_fe_data_initialized = true;
    return;
} // initializeFEData

void
IBFEMethod::reinitializeFEData()
{
    TBOX_ASSERT(d_fe_data_initialized);
    doInitializeFEData(true);
    return;
} // reinitializeFEData

void
IBFEMethod::doInitializeFEData(const bool use_present_data)
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        // Initialize FE equation systems.
        EquationSystems& equation_systems = *d_equation_systems[part];
        if (use_present_data)
        {
            equation_systems.reinit();
        }
        else
        {
            equation_systems.init();
            initializeCoordinates(part);
            initializeVelocity(part);
        }
        updateCoordinateMapping(part);

        // Assemble systems.
        auto& X_system = equation_systems.get_system<System>(COORDS_SYSTEM_NAME);
        auto& dX_system = equation_systems.get_system<System>(COORD_MAPPING_SYSTEM_NAME);
        auto& U_system = equation_systems.get_system<System>(VELOCITY_SYSTEM_NAME);
        auto& F_system = equation_systems.get_system<System>(FORCE_SYSTEM_NAME);

        X_system.assemble_before_solve = false;
        X_system.assemble();

        dX_system.assemble_before_solve = false;
        dX_system.assemble();

        U_system.assemble_before_solve = false;
        U_system.assemble();

        F_system.assemble_before_solve = false;
        F_system.assemble();

        if (d_is_stress_normalization_part[part])
        {
            auto& Phi_system = equation_systems.get_system<LinearImplicitSystem>(PHI_SYSTEM_NAME);
            Phi_system.assemble_before_solve = false;
            Phi_system.assemble();
        }

        if (d_direct_forcing_kinematics_data[part])
        {
            // TODO: none of the data computed by initializeKinematicsData
            // (e.g., the center of mass of the structure) should be different
            // after we repartition. The previous version of this function
            // used the same condition with restarts (i.e., the function
            // argument was !from_restart). We should double check this but I
            // believe this is correct.
            d_direct_forcing_kinematics_data[part]->initializeKinematicsData(!use_present_data);
        }

        // Set up boundary conditions.  Specifically, add appropriate boundary
        // IDs to the BoundaryInfo object associated with the mesh, and add DOF
        // constraints for the nodal forces and velocities.
        const MeshBase& mesh = equation_systems.get_mesh();

        DofMap& F_dof_map = F_system.get_dof_map();
        DofMap& U_dof_map = U_system.get_dof_map();
        const unsigned int F_sys_num = F_system.number();
        const unsigned int U_sys_num = U_system.number();
        MeshBase::const_element_iterator el_it = mesh.elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.elements_end();
        for (; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                const bool at_mesh_bdry = !elem->neighbor(side);
                if (!at_mesh_bdry) continue;

                static const short int dirichlet_bdry_id_set[3] = { FEDataManager::ZERO_DISPLACEMENT_X_BDRY_ID,
                                                                    FEDataManager::ZERO_DISPLACEMENT_Y_BDRY_ID,
                                                                    FEDataManager::ZERO_DISPLACEMENT_Z_BDRY_ID };
                const short int dirichlet_bdry_ids =
                    get_dirichlet_bdry_ids(mesh.boundary_info->boundary_ids(elem, side));
                if (!dirichlet_bdry_ids) continue;

                for (unsigned int n = 0; n < elem->n_nodes(); ++n)
                {
                    if (!elem->is_node_on_side(n, side)) continue;

                    Node* node = elem->get_node(n);
                    mesh.boundary_info->add_node(node, dirichlet_bdry_ids);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (!(dirichlet_bdry_ids & dirichlet_bdry_id_set[d])) continue;
                        if (node->n_dofs(F_sys_num))
                        {
                            const int F_dof_index = node->dof_number(F_sys_num, d, 0);
                            DofConstraintRow F_constraint_row;
                            F_constraint_row[F_dof_index] = 1.0;
                            F_dof_map.add_constraint_row(F_dof_index, F_constraint_row, 0.0, false);
                        }
                        if (node->n_dofs(U_sys_num))
                        {
                            const int U_dof_index = node->dof_number(U_sys_num, d, 0);
                            DofConstraintRow U_constraint_row;
                            U_constraint_row[U_dof_index] = 1.0;
                            U_dof_map.add_constraint_row(U_dof_index, U_constraint_row, 0.0, false);
                        }
                    }
                }
            }
        }
    }
    return;
} // doInitializeFEData

void
IBFEMethod::registerEulerianVariables()
{
    const IntVector<NDIM> ghosts = 1;
    mask_var = new SideVariable<NDIM, double>(d_object_name + "::mask");
    registerVariable(mask_current_idx,
                     mask_new_idx,
                     mask_scratch_idx,
                     mask_var,
                     ghosts,
                     "CONSERVATIVE_COARSEN",
                     "CONSERVATIVE_LINEAR_REFINE");
    return;
} // registerEulerianVariables

void
IBFEMethod::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                     Pointer<GriddingAlgorithm<NDIM> > gridding_alg,
                                     int /*u_data_idx*/,
                                     const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*u_synch_scheds*/,
                                     const std::vector<Pointer<RefineSchedule<NDIM> > >& /*u_ghost_fill_scheds*/,
                                     int /*integrator_step*/,
                                     double /*init_data_time*/,
                                     bool /*initial_time*/)
{
    // Cache pointers to the patch hierarchy and gridding algorithm.
    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;

    // Initialize the FE data manager.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->reinitElementMappings();
    }

    d_is_initialized = true;
    return;
} // initializePatchHierarchy

void
IBFEMethod::registerLoadBalancer(Pointer<LoadBalancer<NDIM> > load_balancer, int workload_data_idx)
{
    TBOX_ASSERT(load_balancer);
    d_load_balancer = load_balancer;
    d_workload_idx = workload_data_idx;
    // if there are parts then update them too
    for (IBTK::FEDataManager *data_manager : d_fe_data_managers)
    {
        TBOX_ASSERT(data_manager);
        data_manager->registerLoadBalancer(load_balancer, workload_data_idx);
    }
    return;
} // registerLoadBalancer

void
IBFEMethod::updateWorkloadEstimates(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/, int /*workload_data_idx*/)
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->updateWorkloadEstimates();
    }
    return;
} // updateWorkloadEstimates

void IBFEMethod::beginDataRedistribution(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                         Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    // intentionally blank
    return;
} // beginDataRedistribution

void IBFEMethod::endDataRedistribution(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                       Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    // if we are not initialized then there is nothing to do
    if (d_is_initialized)
    {
        // At this point in the code SAMRAI has already redistributed the
        // patches (usually by taking into account the number of IB points on
        // each patch). Here is the other half: we inform libMesh of the
        // updated partitioning so that libMesh Elems and Nodes are on the
        // same processor as the relevant SAMRAI patch.
        for (unsigned int part = 0; part < d_num_parts; ++part)
        {
            EquationSystems& equation_systems = *d_fe_data_managers[part]->getEquationSystems();
            MeshBase& mesh = equation_systems.get_mesh();
            BoxPartitioner partitioner(*hierarchy,
                                       equation_systems.get_system(COORDS_SYSTEM_NAME));
            partitioner.repartition(mesh);
        }

        reinitializeFEData();
        for (unsigned int part = 0; part < d_num_parts; ++part)
        {
            d_fe_data_managers[part]->reinitElementMappings();
        }
    }
    return;
} // endDataRedistribution

void
IBFEMethod::initializeLevelData(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                int level_number,
                                double init_data_time,
                                bool can_be_refined,
                                bool initial_time,
                                Pointer<BasePatchLevel<NDIM> > old_level,
                                bool allocate_data)
{
    const int finest_hier_level = hierarchy->getFinestLevelNumber();
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->setPatchHierarchy(hierarchy);
        d_fe_data_managers[part]->setPatchLevels(0, finest_hier_level);
        d_fe_data_managers[part]->initializeLevelData(
            hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);
        if (d_load_balancer && level_number == d_fe_data_managers[part]->getLevelNumber())
        {
            d_load_balancer->setWorkloadPatchDataIndex(d_workload_idx, level_number);
            d_fe_data_managers[part]->updateWorkloadEstimates(level_number, level_number);
        }
    }
    return;
} // initializeLevelData

void
IBFEMethod::resetHierarchyConfiguration(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                        int coarsest_level,
                                        int /*finest_level*/)
{
    const int finest_hier_level = hierarchy->getFinestLevelNumber();
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->setPatchHierarchy(hierarchy);
        d_fe_data_managers[part]->setPatchLevels(0, hierarchy->getFinestLevelNumber());
        d_fe_data_managers[part]->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_hier_level);
    }
    return;
} // resetHierarchyConfiguration

void
IBFEMethod::applyGradientDetector(Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
                                  int level_number,
                                  double error_data_time,
                                  int tag_index,
                                  bool initial_time,
                                  bool uses_richardson_extrapolation_too)
{
    Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->applyGradientDetector(
            hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    }
    return;
} // applyGradientDetector

void
IBFEMethod::putToDatabase(Pointer<Database> db)
{
    db->putInteger("IBFE_METHOD_VERSION", IBFE_METHOD_VERSION);
    db->putInteger("d_num_parts", d_num_parts);
    db->putIntegerArray("d_ghosts", d_ghosts, NDIM);
    db->putBool("d_split_normal_force", d_split_normal_force);
    db->putBool("d_split_tangential_force", d_split_tangential_force);
    db->putBool("d_use_jump_conditions", d_use_jump_conditions);
    db->putBool("d_use_consistent_mass_matrix", d_use_consistent_mass_matrix);
    return;
} // putToDatabase

void
IBFEMethod::writeFEDataToRestartFile(const std::string& restart_dump_dirname, unsigned int time_step_number)
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        const std::string& file_name =
            libmesh_restart_file_name(restart_dump_dirname, time_step_number, part, d_libmesh_restart_file_extension);
        const XdrMODE xdr_mode = (d_libmesh_restart_file_extension == "xdr" ? ENCODE : WRITE);
        const int write_mode = EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA;
        d_equation_systems[part]->write(file_name, xdr_mode, write_mode, /*partition_agnostic*/ true);
    }
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

void
IBFEMethod::computeStressNormalization(PetscVector<double>& Phi_vec,
                                       PetscVector<double>& X_vec,
                                       const double data_time,
                                       const unsigned int part)
{
    // Extract the mesh.
    EquationSystems& equation_systems = *d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems.get_mesh();
    const BoundaryInfo& boundary_info = *mesh.boundary_info;
    const unsigned int dim = mesh.mesh_dimension();

    // Setup extra data needed to compute stresses/forces.

    // Extract the FE systems and DOF maps, and setup the FE objects.
    auto& Phi_system = equation_systems.get_system<LinearImplicitSystem>(PHI_SYSTEM_NAME);
    const DofMap& Phi_dof_map = Phi_system.get_dof_map();
    FEDataManager::SystemDofMapCache& Phi_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(PHI_SYSTEM_NAME);
    std::vector<unsigned int> Phi_dof_indices;
    FEType Phi_fe_type = Phi_dof_map.variable_type(0);
    std::vector<int> Phi_vars(1, 0);

    System& X_system = equation_systems.get_system(COORDS_SYSTEM_NAME);
    std::vector<int> X_vars(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) X_vars[d] = d;

    FEDataInterpolation fe(dim, d_fe_data_managers[part]);
    std::unique_ptr<QBase> qrule_face = QBase::build(QGAUSS, dim - 1, FIFTH);
    fe.attachQuadratureRuleFace(qrule_face.get());
    fe.evalNormalsFace();
    fe.evalQuadraturePointsFace();
    fe.evalQuadratureWeightsFace();
    fe.registerSystem(Phi_system, Phi_vars, Phi_vars); // compute phi and dphi for the Phi system
    const size_t X_sys_idx = fe.registerInterpolatedSystem(X_system, X_vars, X_vars, &X_vec);
    const size_t num_PK1_fcns = d_PK1_stress_fcn_data[part].size();
    std::vector<std::vector<size_t> > PK1_fcn_system_idxs(num_PK1_fcns);
    for (unsigned int k = 0; k < num_PK1_fcns; ++k)
    {
        fe.setupInterpolatedSystemDataIndexes(
            PK1_fcn_system_idxs[k], d_PK1_stress_fcn_data[part][k].system_data, &equation_systems);
    }
    std::vector<size_t> surface_force_fcn_system_idxs;
    fe.setupInterpolatedSystemDataIndexes(
        surface_force_fcn_system_idxs, d_lag_surface_force_fcn_data[part].system_data, &equation_systems);
    std::vector<size_t> surface_pressure_fcn_system_idxs;
    fe.setupInterpolatedSystemDataIndexes(
        surface_pressure_fcn_system_idxs, d_lag_surface_pressure_fcn_data[part].system_data, &equation_systems);
    fe.init(/*use_IB_ghosted_vecs*/ false);

    const std::vector<libMesh::Point>& q_point_face = fe.getQuadraturePointsFace();
    const std::vector<double>& JxW_face = fe.getQuadratureWeightsFace();
    const std::vector<libMesh::Point>& normal_face = fe.getNormalsFace();
    const std::vector<std::vector<double> >& phi_face = fe.getPhiFace(Phi_fe_type);

    const std::vector<std::vector<std::vector<double> > >& fe_interp_var_data = fe.getVarInterpolation();
    const std::vector<std::vector<std::vector<VectorValue<double> > > >& fe_interp_grad_var_data =
        fe.getGradVarInterpolation();

    std::vector<std::vector<const std::vector<double>*> > PK1_var_data(num_PK1_fcns);
    std::vector<std::vector<const std::vector<VectorValue<double> >*> > PK1_grad_var_data(num_PK1_fcns);
    std::vector<const std::vector<double>*> surface_force_var_data, surface_pressure_var_data;
    std::vector<const std::vector<VectorValue<double> >*> surface_force_grad_var_data, surface_pressure_grad_var_data;

    // Setup global and elemental right-hand-side vectors.
    NumericVector<double>* Phi_rhs_vec = Phi_system.rhs;
    Phi_rhs_vec->zero();
    DenseVector<double> Phi_rhs_e;

    // Set up boundary conditions for Phi.
    TensorValue<double> PP, FF, FF_trans, FF_inv_trans;
    VectorValue<double> F, F_s, F_qp, n, x;
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;
        bool reinit_all_data = true;
        for (unsigned short int side = 0; side < elem->n_sides(); ++side)
        {
            // Skip non-physical boundaries.
            if (!is_physical_bdry(elem, side, boundary_info, Phi_dof_map)) continue;

            // Determine if we need to integrate surface forces along this
            // part of the physical boundary; if not, skip the present side.
            const bool at_dirichlet_bdry = is_dirichlet_bdry(elem, side, boundary_info, Phi_dof_map);
            if (at_dirichlet_bdry) continue;

            fe.reinit(elem, side);
            if (reinit_all_data)
            {
                Phi_dof_map_cache.dof_indices(elem, Phi_dof_indices);
                Phi_rhs_e.resize(static_cast<int>(Phi_dof_indices.size()));
                fe.collectDataForInterpolation(elem);
                reinit_all_data = false;
            }
            fe.interpolate(elem, side);
            const unsigned int n_qp = qrule_face->n_points();
            const size_t n_basis = phi_face.size();
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                // X:     reference coordinate
                // x:     current coordinate
                // FF:    deformation gradient associated with x = chi(X,t) (FF = dchi/dX)
                // J:     Jacobian determinant (J = det(FF))
                // N:     unit normal in the reference configuration
                // n:     unit normal in the current configuration
                // dA_da: reference surface area per current surface area (from Nanson's relation)
                const libMesh::Point& X = q_point_face[qp];
                const std::vector<double>& x_data = fe_interp_var_data[qp][X_sys_idx];
                const std::vector<VectorValue<double> >& grad_x_data = fe_interp_grad_var_data[qp][X_sys_idx];
                get_x_and_FF(x, FF, x_data, grad_x_data);
                const double J = std::abs(FF.det());
                FF_trans = FF.transpose();
                tensor_inverse_transpose(FF_inv_trans, FF, NDIM);
                const libMesh::VectorValue<double>& N = normal_face[qp];
                n = (FF_inv_trans * N).unit();
                const double dA_da = 1.0 / (J * (FF_inv_trans * N) * n);

                // Here we build up the boundary value for Phi.
                double Phi = 0.0;
                for (unsigned int k = 0; k < num_PK1_fcns; ++k)
                {
                    if (d_PK1_stress_fcn_data[part][k].fcn)
                    {
                        // Compute the value of the first Piola-Kirchhoff stress
                        // tensor at the quadrature point and add the corresponding
                        // traction force to the right-hand-side vector.
                        fe.setInterpolatedDataPointers(
                            PK1_var_data[k], PK1_grad_var_data[k], PK1_fcn_system_idxs[k], elem, qp);
                        d_PK1_stress_fcn_data[part][k].fcn(PP,
                                                           FF,
                                                           x,
                                                           X,
                                                           elem,
                                                           PK1_var_data[k],
                                                           PK1_grad_var_data[k],
                                                           data_time,
                                                           d_PK1_stress_fcn_data[part][k].ctx);
                        Phi += n * ((PP * FF_trans) * n) / J;
                    }
                }

                if (d_lag_surface_force_fcn_data[part].fcn)
                {
                    // Compute the value of the surface force at the
                    // quadrature point and add the corresponding force to
                    // the right-hand-side vector.
                    fe.setInterpolatedDataPointers(
                        surface_force_var_data, surface_force_grad_var_data, surface_force_fcn_system_idxs, elem, qp);
                    d_lag_surface_force_fcn_data[part].fcn(F_s,
                                                           n,
                                                           N,
                                                           FF,
                                                           x,
                                                           X,
                                                           elem,
                                                           side,
                                                           surface_force_var_data,
                                                           surface_force_grad_var_data,
                                                           data_time,
                                                           d_lag_surface_force_fcn_data[part].ctx);
                    Phi -= n * F_s * dA_da;
                }

                if (d_lag_surface_pressure_fcn_data[part].fcn)
                {
                    // Compute the value of the pressure at the quadrature
                    // point and add the corresponding force to the
                    // right-hand-side vector.
                    double P = 0.0;
                    fe.setInterpolatedDataPointers(surface_pressure_var_data,
                                                   surface_pressure_grad_var_data,
                                                   surface_pressure_fcn_system_idxs,
                                                   elem,
                                                   qp);
                    d_lag_surface_pressure_fcn_data[part].fcn(P,
                                                              n,
                                                              N,
                                                              FF,
                                                              x,
                                                              X,
                                                              elem,
                                                              side,
                                                              surface_pressure_var_data,
                                                              surface_pressure_grad_var_data,
                                                              data_time,
                                                              d_lag_surface_pressure_fcn_data[part].ctx);
                    Phi += P;
                }

                // Add the boundary forces to the right-hand-side vector.
                for (unsigned int i = 0; i < n_basis; ++i)
                {
                    Phi_rhs_e(i) += PENALTY * Phi * phi_face[i][qp] * JxW_face[qp];
                }
            }

            // Apply constraints (e.g., enforce periodic boundary conditions)
            // and add the elemental contributions to the global vector.
            Phi_dof_map.constrain_element_vector(Phi_rhs_e, Phi_dof_indices);
            Phi_rhs_vec->add_vector(Phi_rhs_e, Phi_dof_indices);
        }
    }

    // Solve for Phi.
    Phi_rhs_vec->close();
    Phi_system.solve();
    Phi_dof_map.enforce_constraints_exactly(Phi_system, &Phi_vec);
    copy_and_synch(*Phi_system.solution, Phi_vec, /*close_v_in*/ false);
    return;
}

void
IBFEMethod::computeInteriorForceDensity(PetscVector<double>& G_vec,
                                        PetscVector<double>& X_vec,
                                        PetscVector<double>* Phi_vec,
                                        const double data_time,
                                        const unsigned int part)
{
    // Extract the mesh.
    EquationSystems& equation_systems = *d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems.get_mesh();
    const BoundaryInfo& boundary_info = *mesh.boundary_info;
    const unsigned int dim = mesh.mesh_dimension();

    // Setup global and elemental right-hand-side vectors.
    std::unique_ptr<NumericVector<double> > G_rhs_vec = G_vec.zero_clone();
    DenseVector<double> G_rhs_e[NDIM];

    // First handle the stress contributions.  These are handled separately because
    // each stress function may use a different quadrature rule.
    const size_t num_PK1_fcns = d_PK1_stress_fcn_data[part].size();
    for (unsigned int k = 0; k < num_PK1_fcns; ++k)
    {
        if (!d_PK1_stress_fcn_data[part][k].fcn) continue;

        // Extract the FE systems and DOF maps, and setup the FE object.
        System& G_system = equation_systems.get_system(FORCE_SYSTEM_NAME);
        const DofMap& G_dof_map = G_system.get_dof_map();
        FEDataManager::SystemDofMapCache& G_dof_map_cache =
            *d_fe_data_managers[part]->getDofMapCache(FORCE_SYSTEM_NAME);
        FEType G_fe_type = G_dof_map.variable_type(0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(G_dof_map.variable_type(d) == G_fe_type);
        }
        std::vector<std::vector<unsigned int> > G_dof_indices(NDIM);
        System& X_system = equation_systems.get_system(COORDS_SYSTEM_NAME);
        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;

        FEDataInterpolation fe(dim, d_fe_data_managers[part]);
        std::unique_ptr<QBase> qrule =
            QBase::build(d_PK1_stress_fcn_data[part][k].quad_type, dim, d_PK1_stress_fcn_data[part][k].quad_order);
        std::unique_ptr<QBase> qrule_face =
            QBase::build(d_PK1_stress_fcn_data[part][k].quad_type, dim - 1, d_PK1_stress_fcn_data[part][k].quad_order);
        fe.attachQuadratureRule(qrule.get());
        fe.attachQuadratureRuleFace(qrule_face.get());
        fe.evalNormalsFace();
        fe.evalQuadraturePoints();
        fe.evalQuadraturePointsFace();
        fe.evalQuadratureWeights();
        fe.evalQuadratureWeightsFace();
        fe.registerSystem(G_system, std::vector<int>(), vars); // compute dphi for the force system
        const size_t X_sys_idx = fe.registerInterpolatedSystem(X_system, vars, vars, &X_vec);
        std::vector<size_t> PK1_fcn_system_idxs;
        fe.setupInterpolatedSystemDataIndexes(
            PK1_fcn_system_idxs, d_PK1_stress_fcn_data[part][k].system_data, &equation_systems);
        fe.init(/*use_IB_ghosted_vecs*/ false);

        const std::vector<libMesh::Point>& q_point = fe.getQuadraturePoints();
        const std::vector<double>& JxW = fe.getQuadratureWeights();
        const std::vector<std::vector<VectorValue<double> > >& dphi = fe.getDphi(G_fe_type);

        const std::vector<libMesh::Point>& q_point_face = fe.getQuadraturePointsFace();
        const std::vector<double>& JxW_face = fe.getQuadratureWeightsFace();
        const std::vector<libMesh::Point>& normal_face = fe.getNormalsFace();
        const std::vector<std::vector<double> >& phi_face = fe.getPhiFace(G_fe_type);

        const std::vector<std::vector<std::vector<double> > >& fe_interp_var_data = fe.getVarInterpolation();
        const std::vector<std::vector<std::vector<VectorValue<double> > > >& fe_interp_grad_var_data =
            fe.getGradVarInterpolation();

        std::vector<const std::vector<double>*> PK1_var_data;
        std::vector<const std::vector<VectorValue<double> >*> PK1_grad_var_data;

        // Loop over the elements to compute the right-hand side vector.  This
        // is computed via
        //
        //    rhs_k = -int{PP(s,t) grad phi_k(s)}ds + int{PP(s,t) N(s,t) phi_k(s)}dA(s)
        //
        // This right-hand side vector is used to solve for the nodal values of
        // the interior elastic force density.
        TensorValue<double> PP, FF, FF_inv_trans;
        VectorValue<double> F, F_qp, n, x;
        const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                G_dof_map_cache.dof_indices(elem, G_dof_indices[d], d);
                G_rhs_e[d].resize(static_cast<int>(G_dof_indices[d].size()));
            }
            fe.reinit(elem);
            fe.collectDataForInterpolation(elem);
            fe.interpolate(elem);
            const unsigned int n_qp = qrule->n_points();
            const size_t n_basis = dphi.size();
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                const libMesh::Point& X = q_point[qp];
                const std::vector<double>& x_data = fe_interp_var_data[qp][X_sys_idx];
                const std::vector<VectorValue<double> >& grad_x_data = fe_interp_grad_var_data[qp][X_sys_idx];
                get_x_and_FF(x, FF, x_data, grad_x_data);

                // Compute the value of the first Piola-Kirchhoff stress tensor
                // at the quadrature point and add the corresponding forces to
                // the right-hand-side vector.
                fe.setInterpolatedDataPointers(PK1_var_data, PK1_grad_var_data, PK1_fcn_system_idxs, elem, qp);
                d_PK1_stress_fcn_data[part][k].fcn(
                    PP, FF, x, X, elem, PK1_var_data, PK1_grad_var_data, data_time, d_PK1_stress_fcn_data[part][k].ctx);
                for (unsigned int k = 0; k < n_basis; ++k)
                {
                    F_qp = -PP * dphi[k][qp] * JxW[qp];
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        G_rhs_e[i](k) += F_qp(i);
                    }
                }
            }

            // Loop over the element boundaries.
            for (unsigned short int side = 0; side < elem->n_sides(); ++side)
            {
                // Skip non-physical boundaries.
                if (!is_physical_bdry(elem, side, boundary_info, G_dof_map)) continue;

                // Determine if we need to integrate surface forces along this
                // part of the physical boundary; if not, skip the present side.
                const bool at_dirichlet_bdry = is_dirichlet_bdry(elem, side, boundary_info, G_dof_map);
                const bool integrate_normal_force =
                    (d_split_normal_force && !at_dirichlet_bdry) || (!d_split_normal_force && at_dirichlet_bdry);
                const bool integrate_tangential_force = (d_split_tangential_force && !at_dirichlet_bdry) ||
                                                        (!d_split_tangential_force && at_dirichlet_bdry);
                if (!integrate_normal_force && !integrate_tangential_force) continue;

                fe.reinit(elem, side);
                fe.interpolate(elem, side);
                const unsigned int n_qp = qrule_face->n_points();
                const size_t n_basis = phi_face.size();
                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    const libMesh::Point& X = q_point_face[qp];
                    const std::vector<double>& x_data = fe_interp_var_data[qp][X_sys_idx];
                    const std::vector<VectorValue<double> >& grad_x_data = fe_interp_grad_var_data[qp][X_sys_idx];
                    get_x_and_FF(x, FF, x_data, grad_x_data);
                    tensor_inverse_transpose(FF_inv_trans, FF, NDIM);

                    F.zero();

                    // Compute the value of the first Piola-Kirchhoff stress
                    // tensor at the quadrature point and add the corresponding
                    // traction force to the right-hand-side vector.
                    if (d_PK1_stress_fcn_data[part][k].fcn)
                    {
                        fe.setInterpolatedDataPointers(PK1_var_data, PK1_grad_var_data, PK1_fcn_system_idxs, elem, qp);
                        d_PK1_stress_fcn_data[part][k].fcn(PP,
                                                           FF,
                                                           x,
                                                           X,
                                                           elem,
                                                           PK1_var_data,
                                                           PK1_grad_var_data,
                                                           data_time,
                                                           d_PK1_stress_fcn_data[part][k].ctx);
                        F += PP * normal_face[qp];
                    }

                    n = (FF_inv_trans * normal_face[qp]).unit();

                    if (!integrate_normal_force)
                    {
                        F -= (F * n) * n; // remove the normal component.
                    }

                    if (!integrate_tangential_force)
                    {
                        F -= (F - (F * n) * n); // remove the tangential component.
                    }

                    // Add the boundary forces to the right-hand-side vector.
                    for (unsigned int k = 0; k < n_basis; ++k)
                    {
                        F_qp = F * phi_face[k][qp] * JxW_face[qp];
                        for (unsigned int i = 0; i < NDIM; ++i)
                        {
                            G_rhs_e[i](k) += F_qp(i);
                        }
                    }
                }
            }

            // Apply constraints (e.g., enforce periodic boundary conditions)
            // and add the elemental contributions to the global vector.
            for (unsigned int i = 0; i < NDIM; ++i)
            {
                G_dof_map.constrain_element_vector(G_rhs_e[i], G_dof_indices[i]);
                G_rhs_vec->add_vector(G_rhs_e[i], G_dof_indices[i]);
            }
        }
    }

    // Now account for any additional force contributions.

    // Extract the FE systems and DOF maps, and setup the FE objects.
    System& G_system = equation_systems.get_system(FORCE_SYSTEM_NAME);
    const DofMap& G_dof_map = G_system.get_dof_map();
    FEDataManager::SystemDofMapCache& G_dof_map_cache = *d_fe_data_managers[part]->getDofMapCache(FORCE_SYSTEM_NAME);
    FEType G_fe_type = G_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(G_dof_map.variable_type(d) == G_fe_type);
    }
    std::vector<std::vector<unsigned int> > G_dof_indices(NDIM);
    System& X_system = equation_systems.get_system(COORDS_SYSTEM_NAME);
    System* Phi_system = Phi_vec ? &equation_systems.get_system(PHI_SYSTEM_NAME) : nullptr;
    std::vector<int> vars(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
    std::vector<int> Phi_vars(1, 0);
    std::vector<int> no_vars;

    FEDataInterpolation fe(dim, d_fe_data_managers[part]);
    std::unique_ptr<QBase> qrule = QBase::build(d_default_quad_type[part], dim, d_default_quad_order[part]);
    std::unique_ptr<QBase> qrule_face = QBase::build(d_default_quad_type[part], dim - 1, d_default_quad_order[part]);
    fe.attachQuadratureRule(qrule.get());
    fe.attachQuadratureRuleFace(qrule_face.get());
    fe.evalNormalsFace();
    fe.evalQuadraturePoints();
    fe.evalQuadraturePointsFace();
    fe.evalQuadratureWeights();
    fe.evalQuadratureWeightsFace();
    fe.registerSystem(G_system, vars, vars); // compute phi and dphi for the force system
    const size_t X_sys_idx = fe.registerInterpolatedSystem(X_system, vars, vars, &X_vec);
    const size_t Phi_sys_idx = Phi_vec ? fe.registerInterpolatedSystem(*Phi_system, Phi_vars, no_vars, Phi_vec) :
                                         std::numeric_limits<size_t>::max();
    std::vector<size_t> body_force_fcn_system_idxs;
    fe.setupInterpolatedSystemDataIndexes(
        body_force_fcn_system_idxs, d_lag_body_force_fcn_data[part].system_data, &equation_systems);
    std::vector<size_t> surface_force_fcn_system_idxs;
    fe.setupInterpolatedSystemDataIndexes(
        surface_force_fcn_system_idxs, d_lag_surface_force_fcn_data[part].system_data, &equation_systems);
    std::vector<size_t> surface_pressure_fcn_system_idxs;
    fe.setupInterpolatedSystemDataIndexes(
        surface_pressure_fcn_system_idxs, d_lag_surface_pressure_fcn_data[part].system_data, &equation_systems);
    fe.init(/*use_IB_ghosted_vecs*/ false);

    const std::vector<libMesh::Point>& q_point = fe.getQuadraturePoints();
    const std::vector<double>& JxW = fe.getQuadratureWeights();
    const std::vector<std::vector<double> >& phi = fe.getPhi(G_fe_type);
    const std::vector<std::vector<VectorValue<double> > >& dphi = fe.getDphi(G_fe_type);

    const std::vector<libMesh::Point>& q_point_face = fe.getQuadraturePointsFace();
    const std::vector<double>& JxW_face = fe.getQuadratureWeightsFace();
    const std::vector<libMesh::Point>& normal_face = fe.getNormalsFace();
    const std::vector<std::vector<double> >& phi_face = fe.getPhiFace(G_fe_type);

    const std::vector<std::vector<std::vector<double> > >& fe_interp_var_data = fe.getVarInterpolation();
    const std::vector<std::vector<std::vector<VectorValue<double> > > >& fe_interp_grad_var_data =
        fe.getGradVarInterpolation();

    std::vector<const std::vector<double>*> body_force_var_data, surface_force_var_data, surface_pressure_var_data;
    std::vector<const std::vector<VectorValue<double> >*> body_force_grad_var_data, surface_force_grad_var_data,
        surface_pressure_grad_var_data;

    // Loop over the elements to compute the right-hand side vector.
    TensorValue<double> PP, FF, FF_inv_trans;
    VectorValue<double> F, F_b, F_s, F_qp, n, x;
    boost::multi_array<double, 2> X_node;
    boost::multi_array<double, 1> Phi_node;
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            G_dof_map_cache.dof_indices(elem, G_dof_indices[d], d);
            G_rhs_e[d].resize(static_cast<int>(G_dof_indices[d].size()));
        }
        fe.reinit(elem);
        fe.collectDataForInterpolation(elem);
        fe.interpolate(elem);
        const unsigned int n_qp = qrule->n_points();
        const size_t n_basis = phi.size();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            const libMesh::Point& X = q_point[qp];
            const std::vector<double>& x_data = fe_interp_var_data[qp][X_sys_idx];
            const std::vector<VectorValue<double> >& grad_x_data = fe_interp_grad_var_data[qp][X_sys_idx];
            get_x_and_FF(x, FF, x_data, grad_x_data);
            const double J = std::abs(FF.det());
            tensor_inverse_transpose(FF_inv_trans, FF, NDIM);
            const double Phi =
                Phi_vec ? fe_interp_var_data[qp][Phi_sys_idx][0] : std::numeric_limits<double>::quiet_NaN();

            if (Phi_vec)
            {
                // Compute the value of the first Piola-Kirchhoff stress tensor
                // at the quadrature point and add the corresponding forces to
                // the right-hand-side vector.
                PP = -J * Phi * FF_inv_trans;
                for (unsigned int k = 0; k < n_basis; ++k)
                {
                    F_qp = -PP * dphi[k][qp] * JxW[qp];
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        G_rhs_e[i](k) += F_qp(i);
                    }
                }
            }

            if (d_lag_body_force_fcn_data[part].fcn)
            {
                // Compute the value of the body force at the quadrature
                // point and add the corresponding forces to the
                // right-hand-side vector.
                fe.setInterpolatedDataPointers(
                    body_force_var_data, body_force_grad_var_data, body_force_fcn_system_idxs, elem, qp);
                d_lag_body_force_fcn_data[part].fcn(F_b,
                                                    FF,
                                                    x,
                                                    X,
                                                    elem,
                                                    body_force_var_data,
                                                    body_force_grad_var_data,
                                                    data_time,
                                                    d_lag_body_force_fcn_data[part].ctx);
                for (unsigned int k = 0; k < n_basis; ++k)
                {
                    F_qp = F_b * phi[k][qp] * JxW[qp];
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        G_rhs_e[i](k) += F_qp(i);
                    }
                }
            }
        }

        // Loop over the element boundaries.
        for (unsigned short int side = 0; side < elem->n_sides(); ++side)
        {
            // Skip non-physical boundaries.
            if (!is_physical_bdry(elem, side, boundary_info, G_dof_map)) continue;

            // Determine if we need to compute surface forces along this
            // part of the physical boundary; if not, skip the present side.
            const bool at_dirichlet_bdry = is_dirichlet_bdry(elem, side, boundary_info, G_dof_map);
            const bool integrate_normal_force = !d_split_normal_force && !at_dirichlet_bdry;
            const bool integrate_tangential_force = !d_split_tangential_force && !at_dirichlet_bdry;
            if (!integrate_normal_force && !integrate_tangential_force) continue;

            fe.reinit(elem, side);
            fe.interpolate(elem, side);
            const unsigned int n_qp = qrule_face->n_points();
            const size_t n_basis = phi_face.size();
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                const libMesh::Point& X = q_point_face[qp];
                const std::vector<double>& x_data = fe_interp_var_data[qp][X_sys_idx];
                const std::vector<VectorValue<double> >& grad_x_data = fe_interp_grad_var_data[qp][X_sys_idx];
                get_x_and_FF(x, FF, x_data, grad_x_data);
                const double J = std::abs(FF.det());
                tensor_inverse_transpose(FF_inv_trans, FF, NDIM);
                const libMesh::VectorValue<double>& N = normal_face[qp];
                n = (FF_inv_trans * N).unit();

                F.zero();

                if (d_lag_surface_pressure_fcn_data[part].fcn)
                {
                    // Compute the value of the pressure at the quadrature
                    // point and add the corresponding force to the
                    // right-hand-side vector.
                    double P = 0;
                    fe.setInterpolatedDataPointers(surface_pressure_var_data,
                                                   surface_pressure_grad_var_data,
                                                   surface_pressure_fcn_system_idxs,
                                                   elem,
                                                   qp);
                    d_lag_surface_pressure_fcn_data[part].fcn(P,
                                                              n,
                                                              N,
                                                              FF,
                                                              x,
                                                              X,
                                                              elem,
                                                              side,
                                                              surface_pressure_var_data,
                                                              surface_pressure_grad_var_data,
                                                              data_time,
                                                              d_lag_surface_pressure_fcn_data[part].ctx);
                    F -= P * J * FF_inv_trans * normal_face[qp];
                }

                if (d_lag_surface_force_fcn_data[part].fcn)
                {
                    // Compute the value of the surface force at the
                    // quadrature point and add the corresponding force to
                    // the right-hand-side vector.
                    fe.setInterpolatedDataPointers(
                        surface_force_var_data, surface_force_grad_var_data, surface_force_fcn_system_idxs, elem, qp);
                    d_lag_surface_force_fcn_data[part].fcn(F_s,
                                                           n,
                                                           N,
                                                           FF,
                                                           x,
                                                           X,
                                                           elem,
                                                           side,
                                                           surface_force_var_data,
                                                           surface_force_grad_var_data,
                                                           data_time,
                                                           d_lag_surface_force_fcn_data[part].ctx);
                    F += F_s;
                }

                // Remote the normal component of the boundary force when needed.
                if (!integrate_normal_force) F -= (F * n) * n;

                // Remote the tangential component of the boundary force when needed.
                if (!integrate_tangential_force) F -= (F - (F * n) * n);

                // Add the boundary forces to the right-hand-side vector.
                for (unsigned int k = 0; k < n_basis; ++k)
                {
                    F_qp = F * phi_face[k][qp] * JxW_face[qp];
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        G_rhs_e[i](k) += F_qp(i);
                    }
                }
            }
        }

        // Apply constraints (e.g., enforce periodic boundary conditions)
        // and add the elemental contributions to the global vector.
        for (unsigned int i = 0; i < NDIM; ++i)
        {
            G_dof_map.constrain_element_vector(G_rhs_e[i], G_dof_indices[i]);
            G_rhs_vec->add_vector(G_rhs_e[i], G_dof_indices[i]);
        }
    }

    // Solve for G.
    d_fe_data_managers[part]->computeL2Projection(G_vec, *G_rhs_vec, FORCE_SYSTEM_NAME, d_use_consistent_mass_matrix);
    return;
} // computeInteriorForceDensity

void
IBFEMethod::resetOverlapNodalValues(const std::string& system_name, const std::vector<NumericVector<double>*>& F_vecs)
{
    std::vector<std::unique_ptr<libMesh::PetscVector<double>>> F_ghost_vecs(d_num_parts);
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        if (d_is_overlap_velocity_master_part[part])
        {
            std::vector<numeric_index_type> ghost_idxs(d_overlap_velocity_part_ghost_idxs[part].begin(),
                                                       d_overlap_velocity_part_ghost_idxs[part].end());
            F_ghost_vecs[part] = std::unique_ptr<libMesh::PetscVector<double>>
                (new libMesh::PetscVector<double>(
                F_vecs[part]->comm(), F_vecs[part]->size(), F_vecs[part]->local_size(), ghost_idxs));
            copy_and_synch(*F_vecs[part], *F_ghost_vecs[part], /*close_v_in*/ false);
        }
    }
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        if (d_is_overlap_velocity_part[part])
        {
            const unsigned int master_part = d_overlap_velocity_master_part[part];
            resetOverlapNodalValues(part, system_name, F_vecs[part], F_ghost_vecs[master_part].get());
        }
    }
    return;
} // resetOverlapNodalValues

void
IBFEMethod::resetOverlapNodalValues(const std::string& system_name, const std::vector<PetscVector<double>*>& F_vecs)
{
    resetOverlapNodalValues(system_name, std::vector<NumericVector<double>*>(F_vecs.begin(), F_vecs.end()));
    return;
}

void
IBFEMethod::resetOverlapNodalValues(const unsigned int part,
                                    const std::string& system_name,
                                    libMesh::NumericVector<double>* F_vec,
                                    libMesh::NumericVector<double>* F_master_vec)
{
    TBOX_ASSERT(part < d_num_parts);
    TBOX_ASSERT(d_is_overlap_velocity_part[part]);

    std::array<unsigned int, 2> part_idx;
    part_idx[0] = part;
    part_idx[1] = d_overlap_velocity_master_part[part];

    // Set up basic data structures.
    std::array<EquationSystems*, 2> es;
    std::array<MeshBase*, 2> mesh;
    std::array<const DofMap*, 2> F_dof_map;
    std::array<FEDataManager::SystemDofMapCache*, 2> F_dof_map_cache;
    FEType fe_type;
    std::array<unsigned int, 2> F_sys_num;
    for (int k = 0; k < 2; ++k)
    {
        es[k] = d_fe_data_managers[part_idx[k]]->getEquationSystems();
        auto& F_system = es[k]->get_system<System>(system_name);
        F_dof_map[k] = &F_system.get_dof_map();
        fe_type = F_dof_map[k]->variable_type(0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(fe_type == F_dof_map[k]->variable_type(d));
        }
        F_dof_map_cache[k] = d_fe_data_managers[part_idx[k]]->getDofMapCache(system_name);
        F_sys_num[k] = F_system.number();
        mesh[k] = &es[k]->get_mesh();
    }

    // Find the elements of part2 that contain nodes of part1.
    std::array<std::vector<std::vector<unsigned int> >, 2> F_dof_indices;
    for (int k = 0; k < 2; ++k)
    {
        F_dof_indices[k].resize(NDIM);
    }
    VectorValue<double> F;
    std::array<boost::multi_array<double, 2>, 2> F_node;
    const unsigned int k_slave = 0;
    const unsigned int k_master = 1;
    std::map<dof_id_type, dof_id_type>& node_to_elem_map = d_overlap_velocity_part_node_to_elem_map[part_idx[k_slave]];
    for (const auto& node_elem_pair : node_to_elem_map)
    {
        const Node* const node = mesh[k_slave]->node_ptr(node_elem_pair.first);
        const libMesh::Point& X = *node;
        const Elem* const other_elem = mesh[k_master]->elem_ptr(node_elem_pair.second);

        // Read the velocity from the master mesh.
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            F_dof_map_cache[k_master]->dof_indices(other_elem, F_dof_indices[k_master][d], d);
        }
        const libMesh::Point xi = FEInterface::inverse_map(other_elem->dim(), fe_type, other_elem, X);
        FEComputeData fe_data(*es[k_master], xi);
        FEInterface::compute_data(other_elem->dim(), fe_type, other_elem, fe_data);
        get_values_for_interpolation(F_node[k_master], *F_master_vec, F_dof_indices[k_master]);
        F.zero();
        for (unsigned int l = 0; l < fe_data.shape.size(); ++l)
        {
            const double& p = fe_data.shape[l];
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                F(d) += F_node[k_master][l][d] * p;
            }
        }
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            const int dof_index = node->dof_number(F_sys_num[k_slave], d, 0);
            F_vec->set(dof_index, F(d));
        }
    }
    F_vec->close();
    return;
}

void
IBFEMethod::computeOverlapConstraintForceDensity(std::vector<PetscVector<double>*>& F_vec,
                                                 std::vector<PetscVector<double>*>& X_vec)
{
    if (!d_has_overlap_force_parts) return;

    std::vector<NumericVector<double>*> F_rhs_vec(d_num_parts);
    for (unsigned int k = 0; k < d_num_parts; ++k)
    {
        if (d_is_overlap_force_part[k])
        {
            F_rhs_vec[k] = F_vec[k]->zero_clone().release(); // \todo fix this when we update to C++11
        }
    }
    d_overlap_force_part_max_displacement.resize(d_num_parts);
    std::fill(d_overlap_force_part_max_displacement.begin(),
              d_overlap_force_part_max_displacement.end(),
              std::vector<double>(d_num_parts, 0.0));

    // Add additional forces to constrain overlapping parts.
    unsigned int n_pairs = d_overlap_force_part_idxs.size();
    for (unsigned int k = 0; k < n_pairs; ++k)
    {
        // Setup basic data structures.
        std::array<unsigned int, 2>& part_idx = d_overlap_force_part_idxs[k];
        std::array<std::map<libMesh::dof_id_type, std::map<unsigned int, libMesh::dof_id_type> >, 2>& elem_map =
            d_overlapping_elem_map[k];
        const double kappa = d_overlap_force_part_kappa[k];
        std::array<QBase*, 2>& qrule = d_overlap_force_part_qrule[k];

        std::array<EquationSystems*, 2> es;
        std::array<MeshBase*, 2> mesh;
        std::array<System*, 2> F_system, X_system;
        std::array<const DofMap*, 2> F_dof_map, X_dof_map;
        std::array<FEDataManager::SystemDofMapCache*, 2> F_dof_map_cache, X_dof_map_cache;
        FEType fe_type;
        std::array<std::unique_ptr<FEBase>, 2> fe;
        for (int k = 0; k < 2; ++k)
        {
            es[k] = d_fe_data_managers[part_idx[k]]->getEquationSystems();
            F_system[k] = &es[k]->get_system<System>(FORCE_SYSTEM_NAME);
            X_system[k] = &es[k]->get_system<System>(COORDS_SYSTEM_NAME);
            F_dof_map[k] = &F_system[k]->get_dof_map();
            X_dof_map[k] = &X_system[k]->get_dof_map();
            FEType fe_type = F_dof_map[k]->variable_type(0);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                TBOX_ASSERT(fe_type == F_dof_map[k]->variable_type(d));
                TBOX_ASSERT(fe_type == X_dof_map[k]->variable_type(d));
            }
            F_dof_map_cache[k] = d_fe_data_managers[part_idx[k]]->getDofMapCache(FORCE_SYSTEM_NAME);
            X_dof_map_cache[k] = d_fe_data_managers[part_idx[k]]->getDofMapCache(COORDS_SYSTEM_NAME);
            mesh[k] = &es[k]->get_mesh();
            fe[k] = FEBase::build(mesh[k]->mesh_dimension(), fe_type);
            fe[k]->attach_quadrature_rule(qrule[k]);
            fe[k]->get_xyz();
            fe[k]->get_JxW();
            fe[k]->get_phi();
        }

        // Loop over the elements that overlap elements of the other part.
        std::array<std::vector<std::vector<unsigned int> >, 2> F_dof_indices, X_dof_indices;
        for (int k = 0; k < 2; ++k)
        {
            F_dof_indices[k].resize(NDIM);
            X_dof_indices[k].resize(NDIM);
        }
        std::array<std::array<DenseVector<double>, NDIM>, 2> F_rhs_e;
        std::array<boost::multi_array<double, 2>, 2> x_node;
        VectorValue<double> F, F_qp, other_x, x;
        for (int k = 0; k < 2; ++k)
        {
            const int k_next = (k + 1) % 2;
            const std::vector<libMesh::Point>& q_point = fe[k]->get_xyz();
            const std::vector<double>& JxW = fe[k]->get_JxW();
            const std::vector<std::vector<double> >& phi = fe[k]->get_phi();
            for (const auto& elem_pair : elem_map[k])
            {
                const Elem* const elem = mesh[k]->elem_ptr(elem_pair.first);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    F_dof_map_cache[k]->dof_indices(elem, F_dof_indices[k][d], d);
                    X_dof_map_cache[k]->dof_indices(elem, X_dof_indices[k][d], d);
                    F_rhs_e[k][d].resize(static_cast<int>(F_dof_indices[k][d].size()));
                }
                fe[k]->reinit(elem);
                get_values_for_interpolation(x_node[k], *X_vec[part_idx[k]], X_dof_indices[k]);
                const size_t n_basis = phi.size();
                for (auto qp_it = elem_pair.second.begin(); qp_it != elem_pair.second.end(); ++qp_it)
                {
                    const unsigned int qp = qp_it->first;
                    const libMesh::Point& X = q_point[qp];
                    interpolate(x, qp, x_node[k], phi);
                    const Elem* const other_elem = mesh[k_next]->elem_ptr(qp_it->second);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        X_dof_map_cache[k_next]->dof_indices(other_elem, X_dof_indices[k_next][d], d);
                    }
                    const libMesh::Point xi = FEInterface::inverse_map(other_elem->dim(), fe_type, other_elem, X);
                    FEComputeData fe_data(*es[k_next], xi);
                    FEInterface::compute_data(other_elem->dim(), fe_type, other_elem, fe_data);
                    get_values_for_interpolation(x_node[k_next], *X_vec[part_idx[k_next]], X_dof_indices[k_next]);
                    other_x.zero();
                    for (unsigned int l = 0; l < fe_data.shape.size(); ++l)
                    {
                        const double& p = fe_data.shape[l];
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            other_x(d) += x_node[k_next][l][d] * p;
                        }
                    }
                    VectorValue<double> disp = (other_x - x);
                    VectorValue<double> F = kappa * disp;
                    d_overlap_force_part_max_displacement[part_idx[k]][part_idx[k_next]] =
                        std::max(d_overlap_force_part_max_displacement[part_idx[k]][part_idx[k_next]], disp.norm());
                    for (unsigned int l = 0; l < n_basis; ++l)
                    {
                        F_qp = F * phi[l][qp] * JxW[qp];
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            F_rhs_e[k][d](l) += F_qp(d);
                        }
                    }
                }

                // Apply constraints (e.g., enforce periodic boundary
                // conditions) and add the elemental contributions to the global
                // vector.
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    F_dof_map[k]->constrain_element_vector(F_rhs_e[k][d], F_dof_indices[k][d]);
                    F_rhs_vec[k]->add_vector(F_rhs_e[k][d], F_dof_indices[k][d]);
                }
            }
        }

        // Report the maximum displacement between this pair of parts.
        SAMRAI_MPI::maxReduction(&d_overlap_force_part_max_displacement[part_idx[0]][part_idx[1]], 1);
        plog << "max displacement from part " << part_idx[0] << " to part " << part_idx[1] << ": "
             << d_overlap_force_part_max_displacement[part_idx[0]][part_idx[1]] << "\n";
        SAMRAI_MPI::maxReduction(&d_overlap_force_part_max_displacement[part_idx[1]][part_idx[0]], 1);
        plog << "max displacement from part " << part_idx[1] << " to part " << part_idx[0] << ": "
             << d_overlap_force_part_max_displacement[part_idx[1]][part_idx[0]] << "\n";
    }

    // Solve for the constraint force densities.
    for (unsigned int k = 0; k < d_num_parts; ++k)
    {
        if (d_is_overlap_force_part[k])
        {
            std::unique_ptr<NumericVector<double> > F_tmp_vec = F_vec[k]->zero_clone();
            d_fe_data_managers[k]->computeL2Projection(
                *F_tmp_vec, *F_rhs_vec[k], FORCE_SYSTEM_NAME, d_use_consistent_mass_matrix);
            F_vec[k]->add(*F_tmp_vec);
        }
    }

    // Cleanup.
    for (unsigned int k = 0; k < d_num_parts; ++k)
    {
        if (d_is_overlap_force_part[k])
        {
            delete F_rhs_vec[k];
        }
    }
    return;
} // computeOverlapConstraintForceDensity

void
IBFEMethod::spreadTransmissionForceDensity(const int f_data_idx,
                                           PetscVector<double>& X_ghost_vec,
                                           RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                                           const double data_time,
                                           const unsigned int part)
{
    if (!d_split_normal_force && !d_split_tangential_force) return;

    // Check to see if we need to integrate the surface forces.
    const bool integrate_normal_force =
        d_split_normal_force && !d_use_jump_conditions && !d_is_stress_normalization_part[part];
    const bool integrate_tangential_force = d_split_tangential_force;
    if (!integrate_normal_force && !integrate_tangential_force) return;

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Make a copy of the Eulerian data.
    Pointer<hier::Variable<NDIM> > f_var;
    var_db->mapIndexToVariable(f_data_idx, f_var);
    const int f_copy_data_idx = var_db->registerClonedPatchDataIndex(f_var, f_data_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(f_copy_data_idx);
    }
    Pointer<HierarchyDataOpsReal<NDIM, double> > f_data_ops =
        HierarchyDataOpsManager<NDIM>::getManager()->getOperationsDouble(f_var, d_hierarchy, true);
    f_data_ops->swapData(f_copy_data_idx, f_data_idx);
    f_data_ops->setToScalar(f_data_idx, 0.0, /*interior_only*/ false);

    // Extract the mesh.
    EquationSystems& equation_systems = *d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems.get_mesh();
    const BoundaryInfo& boundary_info = *mesh.boundary_info;
    const unsigned int dim = mesh.mesh_dimension();

    // Extract the FE systems and DOF maps, and setup the FE object.
    System& G_system = equation_systems.get_system(FORCE_SYSTEM_NAME);
    const DofMap& G_dof_map = G_system.get_dof_map();
    FEType G_fe_type = G_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(G_dof_map.variable_type(d) == G_fe_type);
    }
    std::vector<std::vector<unsigned int> > G_dof_indices(NDIM);
    System& X_system = equation_systems.get_system(COORDS_SYSTEM_NAME);
    std::vector<int> vars(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;

    FEDataInterpolation fe(dim, d_fe_data_managers[part]);
    std::unique_ptr<QBase> default_qrule_face =
        QBase::build(d_default_quad_type[part], dim - 1, d_default_quad_order[part]);
    std::unique_ptr<QBase> qrule_face;
    fe.attachQuadratureRuleFace(default_qrule_face.get());
    fe.evalNormalsFace();
    fe.evalQuadraturePointsFace();
    fe.evalQuadratureWeightsFace();
    const size_t X_sys_idx = fe.registerInterpolatedSystem(X_system, vars, vars, &X_ghost_vec);
    const size_t num_PK1_fcns = d_PK1_stress_fcn_data[part].size();
    std::vector<std::vector<size_t> > PK1_fcn_system_idxs(num_PK1_fcns);
    for (unsigned int k = 0; k < num_PK1_fcns; ++k)
    {
        fe.setupInterpolatedSystemDataIndexes(
            PK1_fcn_system_idxs[k], d_PK1_stress_fcn_data[part][k].system_data, &equation_systems);
    }
    std::vector<size_t> surface_force_fcn_system_idxs;
    fe.setupInterpolatedSystemDataIndexes(
        surface_force_fcn_system_idxs, d_lag_surface_force_fcn_data[part].system_data, &equation_systems);
    std::vector<size_t> surface_pressure_fcn_system_idxs;
    fe.setupInterpolatedSystemDataIndexes(
        surface_pressure_fcn_system_idxs, d_lag_surface_pressure_fcn_data[part].system_data, &equation_systems);
    fe.init(/*use_IB_ghosted_vecs*/ true);

    const std::vector<libMesh::Point>& q_point_face = fe.getQuadraturePointsFace();
    const std::vector<double>& JxW_face = fe.getQuadratureWeightsFace();
    const std::vector<libMesh::Point>& normal_face = fe.getNormalsFace();

    const std::vector<std::vector<std::vector<double> > >& fe_interp_var_data = fe.getVarInterpolation();
    const std::vector<std::vector<std::vector<VectorValue<double> > > >& fe_interp_grad_var_data =
        fe.getGradVarInterpolation();

    std::vector<std::vector<const std::vector<double>*> > PK1_var_data(num_PK1_fcns);
    std::vector<std::vector<const std::vector<VectorValue<double> >*> > PK1_grad_var_data(num_PK1_fcns);
    std::vector<const std::vector<double>*> surface_force_var_data, surface_pressure_var_data;
    std::vector<const std::vector<VectorValue<double> >*> surface_force_grad_var_data, surface_pressure_grad_var_data;

    // Loop over the patches to spread the transmission elastic force density
    // onto the grid.
    const std::vector<std::vector<Elem*> >& active_patch_element_map =
        d_fe_data_managers[part]->getActivePatchElementMap();
    const int level_num = d_fe_data_managers[part]->getLevelNumber();
    TensorValue<double> PP, FF, FF_inv_trans;
    VectorValue<double> F, F_s, n, x;
    double P;
    std::vector<double> T_bdry, x_bdry;
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const std::vector<Elem*>& patch_elems = active_patch_element_map[local_patch_num];
        const size_t num_active_patch_elems = patch_elems.size();
        if (num_active_patch_elems == 0) continue;

        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const patch_dx = patch_geom->getDx();
        const double patch_dx_min = *std::min_element(patch_dx, patch_dx + NDIM);

        // Loop over the elements and compute the values to be spread and the
        // positions of the quadrature points.
        T_bdry.clear();
        x_bdry.clear();
        int qp_offset = 0;
        for (size_t e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems[e_idx];
            const bool touches_physical_bdry = has_physical_bdry(elem, boundary_info, G_dof_map);
            if (!touches_physical_bdry) continue;

            fe.reinit(elem);
            fe.collectDataForInterpolation(elem);
            const boost::multi_array<double, 2>& X_node = fe.getElemData(elem, X_sys_idx);

            // Loop over the element boundaries.
            for (unsigned short int side = 0; side < elem->n_sides(); ++side)
            {
                // Skip non-physical boundaries.
                if (!is_physical_bdry(elem, side, boundary_info, G_dof_map)) continue;

                // Skip Dirichlet boundaries.
                if (is_dirichlet_bdry(elem, side, boundary_info, G_dof_map)) continue;

                // Construct a side element.
                std::unique_ptr<Elem> side_elem = elem->build_side(side, /*proxy*/ false);
                const bool qrule_changed = d_fe_data_managers[part]->updateSpreadQuadratureRule(
                    qrule_face, d_spread_spec[part], side_elem.get(), X_node, patch_dx_min);
                if (qrule_changed)
                {
                    fe.attachQuadratureRuleFace(qrule_face.get());
                }
                fe.reinit(elem, side);
                fe.interpolate(elem, side);
                const unsigned int n_qp = qrule_face->n_points();
                T_bdry.resize(T_bdry.size() + NDIM * n_qp);
                x_bdry.resize(x_bdry.size() + NDIM * n_qp);
                for (unsigned int qp = 0; qp < n_qp; ++qp, ++qp_offset)
                {
                    const libMesh::Point& X = q_point_face[qp];
                    const std::vector<double>& x_data = fe_interp_var_data[qp][X_sys_idx];
                    const std::vector<VectorValue<double> >& grad_x_data = fe_interp_grad_var_data[qp][X_sys_idx];
                    get_x_and_FF(x, FF, x_data, grad_x_data);
                    const double J = std::abs(FF.det());
                    tensor_inverse_transpose(FF_inv_trans, FF, NDIM);
                    const libMesh::VectorValue<double>& N = normal_face[qp];
                    n = (FF_inv_trans * N).unit();

                    F.zero();

                    for (unsigned int k = 0; k < num_PK1_fcns; ++k)
                    {
                        if (d_PK1_stress_fcn_data[part][k].fcn)
                        {
                            // Compute the value of the first Piola-Kirchhoff stress
                            // tensor at the quadrature point and compute the
                            // corresponding force.
                            fe.setInterpolatedDataPointers(
                                PK1_var_data[k], PK1_grad_var_data[k], PK1_fcn_system_idxs[k], elem, qp);
                            d_PK1_stress_fcn_data[part][k].fcn(PP,
                                                               FF,
                                                               x,
                                                               X,
                                                               elem,
                                                               PK1_var_data[k],
                                                               PK1_grad_var_data[k],
                                                               data_time,
                                                               d_PK1_stress_fcn_data[part][k].ctx);
                            F -= PP * normal_face[qp] * JxW_face[qp];
                        }
                    }

                    if (d_lag_surface_pressure_fcn_data[part].fcn)
                    {
                        // Compute the value of the pressure at the quadrature
                        // point and compute the corresponding force.
                        fe.setInterpolatedDataPointers(surface_pressure_var_data,
                                                       surface_pressure_grad_var_data,
                                                       surface_pressure_fcn_system_idxs,
                                                       elem,
                                                       qp);
                        d_lag_surface_pressure_fcn_data[part].fcn(P,
                                                                  n,
                                                                  N,
                                                                  FF,
                                                                  x,
                                                                  X,
                                                                  elem,
                                                                  side,
                                                                  surface_pressure_var_data,
                                                                  surface_pressure_grad_var_data,
                                                                  data_time,
                                                                  d_lag_surface_pressure_fcn_data[part].ctx);
                        F -= P * J * FF_inv_trans * normal_face[qp] * JxW_face[qp];
                    }

                    if (d_lag_surface_force_fcn_data[part].fcn)
                    {
                        // Compute the value of the surface force at the
                        // quadrature point and compute the corresponding force.
                        fe.setInterpolatedDataPointers(surface_force_var_data,
                                                       surface_force_grad_var_data,
                                                       surface_force_fcn_system_idxs,
                                                       elem,
                                                       qp);
                        d_lag_surface_force_fcn_data[part].fcn(F_s,
                                                               n,
                                                               N,
                                                               FF,
                                                               x,
                                                               X,
                                                               elem,
                                                               side,
                                                               surface_force_var_data,
                                                               surface_force_grad_var_data,
                                                               data_time,
                                                               d_lag_surface_force_fcn_data[part].ctx);
                        F += F_s * JxW_face[qp];
                    }

                    // Remote the normal component of the boundary force when needed.
                    if (!integrate_normal_force) F -= (F * n) * n;

                    // Remote the tangential component of the boundary force when needed.
                    if (!integrate_tangential_force) F -= (F - (F * n) * n);

                    const int idx = NDIM * qp_offset;
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        T_bdry[idx + i] = F(i);
                    }
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        x_bdry[idx + i] = x(i);
                    }
                }
            }
        }

        if (qp_offset == 0) continue;

        // Spread the boundary forces to the grid.
        const std::string& spread_kernel_fcn = d_spread_spec[part].kernel_fcn;
        const hier::IntVector<NDIM>& ghost_width = d_fe_data_managers[part]->getGhostCellWidth();
        const Box<NDIM> spread_box = Box<NDIM>::grow(patch->getBox(), ghost_width);
        Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(f_data_idx);
        LEInteractor::spread(f_data, T_bdry, NDIM, x_bdry, NDIM, patch, spread_box, spread_kernel_fcn);
        if (f_phys_bdry_op)
        {
            f_phys_bdry_op->setPatchDataIndex(f_data_idx);
            f_phys_bdry_op->accumulateFromPhysicalBoundaryData(*patch, data_time, f_data->getGhostCellWidth());
        }
    }

    // Accumulate data.
    f_data_ops->swapData(f_copy_data_idx, f_data_idx);
    f_data_ops->add(f_data_idx, f_data_idx, f_copy_data_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(f_copy_data_idx);
    }
    var_db->removePatchDataIndex(f_copy_data_idx);
    return;
} // spreadTransmissionForceDensity

void
IBFEMethod::imposeJumpConditions(const int f_data_idx,
                                 PetscVector<double>& G_ghost_vec,
                                 PetscVector<double>& X_ghost_vec,
                                 const double data_time,
                                 const unsigned int part)
{
    if (!d_split_normal_force && !d_split_tangential_force) return;

    // Check to see if we need to integrate the normal surface force.
    const bool integrate_normal_force =
        d_split_normal_force && d_use_jump_conditions && !d_is_stress_normalization_part[part];
    if (!integrate_normal_force) return;

    // Extract the mesh.
    EquationSystems& equation_systems = *d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems.get_mesh();
    const BoundaryInfo& boundary_info = *mesh.boundary_info;
    const unsigned int dim = mesh.mesh_dimension();
    TBOX_ASSERT(dim == NDIM);

    // Extract the FE systems and DOF maps, and setup the FE object.
    System& G_system = equation_systems.get_system(FORCE_SYSTEM_NAME);
    DofMap& G_dof_map = G_system.get_dof_map();
    FEType G_fe_type = G_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(G_dof_map.variable_type(d) == G_fe_type);
    }
    std::vector<std::vector<unsigned int> > G_dof_indices(NDIM);
    System& X_system = equation_systems.get_system(COORDS_SYSTEM_NAME);
    DofMap& X_dof_map = X_system.get_dof_map();
    std::vector<int> vars(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
    std::vector<int> no_vars;

    FEDataInterpolation fe(dim, d_fe_data_managers[part]);
    std::unique_ptr<QBase> qrule_face = QBase::build(d_default_quad_type[part], dim - 1, d_default_quad_order[part]);
    fe.attachQuadratureRuleFace(qrule_face.get());
    fe.evalQuadraturePointsFace();
    fe.evalNormalsFace();
    const size_t G_sys_idx = fe.registerInterpolatedSystem(G_system, vars, no_vars, &G_ghost_vec);
    const size_t X_sys_idx = fe.registerInterpolatedSystem(X_system, vars, vars, &X_ghost_vec);
    const size_t num_PK1_fcns = d_PK1_stress_fcn_data[part].size();
    std::vector<std::vector<size_t> > PK1_fcn_system_idxs(num_PK1_fcns);
    for (unsigned int k = 0; k < num_PK1_fcns; ++k)
    {
        fe.setupInterpolatedSystemDataIndexes(
            PK1_fcn_system_idxs[k], d_PK1_stress_fcn_data[part][k].system_data, &equation_systems);
    }
    std::vector<size_t> surface_force_fcn_system_idxs;
    fe.setupInterpolatedSystemDataIndexes(
        surface_force_fcn_system_idxs, d_lag_surface_force_fcn_data[part].system_data, &equation_systems);
    std::vector<size_t> surface_pressure_fcn_system_idxs;
    fe.setupInterpolatedSystemDataIndexes(
        surface_pressure_fcn_system_idxs, d_lag_surface_pressure_fcn_data[part].system_data, &equation_systems);
    fe.init(/*use_IB_ghosted_vecs*/ true);

    const std::vector<libMesh::Point>& q_point_face = fe.getQuadraturePointsFace();
    const std::vector<libMesh::Point>& normal_face = fe.getNormalsFace();

    const std::vector<std::vector<std::vector<double> > >& fe_interp_var_data = fe.getVarInterpolation();
    const std::vector<std::vector<std::vector<VectorValue<double> > > >& fe_interp_grad_var_data =
        fe.getGradVarInterpolation();

    std::vector<std::vector<const std::vector<double>*> > PK1_var_data(num_PK1_fcns);
    std::vector<std::vector<const std::vector<VectorValue<double> >*> > PK1_grad_var_data(num_PK1_fcns);
    std::vector<const std::vector<double>*> surface_force_var_data, surface_pressure_var_data;
    std::vector<const std::vector<VectorValue<double> >*> surface_force_grad_var_data, surface_pressure_grad_var_data;

    // Loop over the patches to impose jump conditions on the Eulerian grid that
    // are determined from the interior and transmission elastic force
    // densities.
    const std::vector<std::vector<Elem*> >& active_patch_element_map =
        d_fe_data_managers[part]->getActivePatchElementMap();
    const int level_num = d_fe_data_managers[part]->getLevelNumber();
    TensorValue<double> PP, FF, FF_inv_trans;
    VectorValue<double> G, F, F_s, n, x;
    std::vector<libMesh::Point> X_node_cache, x_node_cache;
    IBTK::Point x_min, x_max;
    std::vector<std::vector<unsigned int> > side_dof_indices(NDIM);
    std::vector<libMesh::Point> intersection_ref_coords;
    std::vector<SideIndex<NDIM> > intersection_indices;
    std::vector<std::pair<double, libMesh::Point> > intersections;
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    const IntVector<NDIM>& ratio = level->getRatio();
    const Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const std::vector<Elem*>& patch_elems = active_patch_element_map[local_patch_num];
        const size_t num_active_patch_elems = patch_elems.size();
        if (num_active_patch_elems == 0) continue;

        const Pointer<Patch<NDIM> > patch = level->getPatch(p());
        Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(f_data_idx);
        const Box<NDIM>& patch_box = patch->getBox();
        const CellIndex<NDIM>& patch_lower = patch_box.lower();
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const x_lower = patch_geom->getXLower();
        const double* const dx = patch_geom->getDx();

        SideData<NDIM, int> num_intersections(patch_box, 1, IntVector<NDIM>(0));
        num_intersections.fillAll(0);

        // Loop over the elements.
        for (size_t e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems[e_idx];
            const bool touches_physical_bdry = has_physical_bdry(elem, boundary_info, G_dof_map);
            if (!touches_physical_bdry) continue;

            fe.reinit(elem);
            fe.collectDataForInterpolation(elem);

            // Loop over the element boundaries.
            for (unsigned short int side = 0; side < elem->n_sides(); ++side)
            {
                // Skip non-physical boundaries.
                if (!is_physical_bdry(elem, side, boundary_info, G_dof_map)) continue;

                // Skip Dirichlet boundaries.
                if (is_dirichlet_bdry(elem, side, boundary_info, G_dof_map)) continue;

                // Construct a side element.
                std::unique_ptr<Elem> side_elem = elem->build_side(side, /*proxy*/ false);
                const unsigned int n_node_side = side_elem->n_nodes();
                for (int d = 0; d < NDIM; ++d)
                {
                    X_dof_map.dof_indices(side_elem.get(), side_dof_indices[d], d);
                }

                // Cache the nodal and physical coordinates of the side element,
                // determine the bounding box of the current configuration of
                // the side element, and set the nodal coordinates to correspond
                // to the physical coordinates.
                X_node_cache.resize(n_node_side);
                x_node_cache.resize(n_node_side);
                x_min = IBTK::Point::Constant(std::numeric_limits<double>::max());
                x_max = IBTK::Point::Constant(-std::numeric_limits<double>::max());
                for (unsigned int k = 0; k < n_node_side; ++k)
                {
                    X_node_cache[k] = side_elem->point(k);
                    libMesh::Point& x = x_node_cache[k];
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        x(d) = X_ghost_vec(side_dof_indices[d][k]);
                        x_min[d] = std::min(x_min[d], x(d));
                        x_max[d] = std::max(x_max[d], x(d));
                    }
                    side_elem->point(k) = x;
                }
                Box<NDIM> box(IndexUtilities::getCellIndex(&x_min[0], grid_geom, ratio),
                              IndexUtilities::getCellIndex(&x_max[0], grid_geom, ratio));
                box.grow(IntVector<NDIM>(1));
                box = box * patch_box;

                // Loop over coordinate directions and look for intersections
                // with the background fluid grid.
                intersection_ref_coords.clear();
                intersection_indices.clear();
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    // Setup a unit vector pointing in the coordinate direction
                    // of interest.
                    VectorValue<double> q;
                    q(axis) = 1.0;

                    // Loop over the relevant range of indices.
                    Box<NDIM> axis_box = box;
                    axis_box.lower(axis) = 0;
                    axis_box.upper(axis) = 0;
                    for (BoxIterator<NDIM> b(axis_box); b; b++)
                    {
                        const Index<NDIM>& i_c = b();
                        libMesh::Point r;
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            r(d) =
                                (d == axis ? 0.0 :
                                             x_lower[d] + dx[d] * (static_cast<double>(i_c(d) - patch_lower[d]) + 0.5));
                        }
#if (NDIM == 2)
                        intersect_line_with_edge(intersections, static_cast<Edge*>(side_elem.get()), r, q);
#endif
#if (NDIM == 3)
                        intersect_line_with_face(intersections, static_cast<Face*>(side_elem.get()), r, q);
#endif
                        for (const auto& intersection : intersections)
                        {
                            libMesh::Point x = r + intersection.first * q;
                            SideIndex<NDIM> i_s(i_c, axis, 0);
                            i_s(axis) = std::floor((x(axis) - x_lower[axis]) / dx[axis] + 0.5) + patch_lower[axis];
                            intersection_ref_coords.push_back(intersection.second);
                            intersection_indices.push_back(i_s);
                            num_intersections(i_s) += 1;
                        }
                    }
                }

                // Restore the element coordinates.
                for (unsigned int k = 0; k < n_node_side; ++k)
                {
                    side_elem->point(k) = X_node_cache[k];
                }

                // If there are no intersection points, then continue to the
                // next side.
                if (intersection_ref_coords.empty()) continue;

                // Evaluate the jump conditions and apply them to the Eulerian
                // grid.
                const bool impose_dp_dn_jumps = false;
                static const double TOL = std::sqrt(std::numeric_limits<double>::epsilon());
                fe.reinit(elem, side, TOL, &intersection_ref_coords);
                fe.interpolate(elem, side);
                const size_t n_qp = intersection_ref_coords.size();
                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    const SideIndex<NDIM>& i_s = intersection_indices[qp];
                    const unsigned int axis = i_s.getAxis();
                    const libMesh::Point& X = q_point_face[qp];
                    const std::vector<double>& x_data = fe_interp_var_data[qp][X_sys_idx];
                    const std::vector<VectorValue<double> >& grad_x_data = fe_interp_grad_var_data[qp][X_sys_idx];
                    get_x_and_FF(x, FF, x_data, grad_x_data);
                    const double J = std::abs(FF.det());
                    tensor_inverse_transpose(FF_inv_trans, FF, NDIM);
                    const libMesh::VectorValue<double>& N = normal_face[qp];
                    n = (FF_inv_trans * N).unit();
                    const double dA_da = 1.0 / (J * (FF_inv_trans * N) * n);
                    const std::vector<double>& G_data = fe_interp_var_data[qp][G_sys_idx];
                    std::copy(G_data.begin(), G_data.end(), &G(0));
#if !defined(NDEBUG)
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (d == axis)
                        {
                            const double x_lower_bound = x_lower[d] +
                                                         (static_cast<double>(i_s(d) - patch_lower[d]) - 0.5) * dx[d] -
                                                         std::sqrt(std::numeric_limits<double>::epsilon());
                            const double x_upper_bound = x_lower[d] +
                                                         (static_cast<double>(i_s(d) - patch_lower[d]) + 0.5) * dx[d] +
                                                         std::sqrt(std::numeric_limits<double>::epsilon());
                            TBOX_ASSERT(x_lower_bound <= x(d) && x(d) <= x_upper_bound);
                        }
                        else
                        {
                            const double x_intersection =
                                x_lower[d] + (static_cast<double>(i_s(d) - patch_lower[d]) + 0.5) * dx[d];
                            const double x_interp = x(d);
                            const double rel_diff =
                                std::abs(x_intersection - x_interp) /
                                std::max(1.0, std::max(std::abs(x_intersection), std::abs(x_interp)));
                            TBOX_ASSERT(rel_diff <= std::sqrt(std::numeric_limits<double>::epsilon()));
                        }
                    }
#endif
                    F.zero();

                    for (unsigned int k = 0; k < num_PK1_fcns; ++k)
                    {
                        if (d_PK1_stress_fcn_data[part][k].fcn)
                        {
                            // Compute the value of the first Piola-Kirchhoff
                            // stress tensor at the quadrature point and compute
                            // the corresponding force.
                            fe.setInterpolatedDataPointers(
                                PK1_var_data[k], PK1_grad_var_data[k], PK1_fcn_system_idxs[k], elem, qp);
                            d_PK1_stress_fcn_data[part][k].fcn(PP,
                                                               FF,
                                                               x,
                                                               X,
                                                               elem,
                                                               PK1_var_data[k],
                                                               PK1_grad_var_data[k],
                                                               data_time,
                                                               d_PK1_stress_fcn_data[part][k].ctx);
                            F -= PP * normal_face[qp];
                        }
                    }

                    if (d_lag_surface_pressure_fcn_data[part].fcn)
                    {
                        // Compute the value of the pressure at the quadrature
                        // point and compute the corresponding force.
                        double P = 0.0;
                        fe.setInterpolatedDataPointers(surface_pressure_var_data,
                                                       surface_pressure_grad_var_data,
                                                       surface_pressure_fcn_system_idxs,
                                                       elem,
                                                       qp);
                        d_lag_surface_pressure_fcn_data[part].fcn(P,
                                                                  n,
                                                                  N,
                                                                  FF,
                                                                  x,
                                                                  X,
                                                                  elem,
                                                                  side,
                                                                  surface_pressure_var_data,
                                                                  surface_pressure_grad_var_data,
                                                                  data_time,
                                                                  d_lag_surface_pressure_fcn_data[part].ctx);
                        F -= P * J * FF_inv_trans * normal_face[qp];
                    }

                    if (d_lag_surface_force_fcn_data[part].fcn)
                    {
                        // Compute the value of the surface force at the
                        // quadrature point and compute the corresponding force.
                        fe.setInterpolatedDataPointers(surface_force_var_data,
                                                       surface_force_grad_var_data,
                                                       surface_force_fcn_system_idxs,
                                                       elem,
                                                       qp);
                        d_lag_surface_force_fcn_data[part].fcn(F_s,
                                                               n,
                                                               N,
                                                               FF,
                                                               x,
                                                               X,
                                                               elem,
                                                               side,
                                                               surface_force_var_data,
                                                               surface_force_grad_var_data,
                                                               data_time,
                                                               d_lag_surface_force_fcn_data[part].ctx);
                        F += F_s;
                    }

                    F *= dA_da;

                    // Determine the value of the interior force density at the
                    // boundary, and convert it to force per unit volume in the
                    // current configuration.  This value determines the
                    // discontinuity in the normal derivative of the pressure at
                    // the fluid-structure interface.
                    if (impose_dp_dn_jumps)
                    {
                        G /= J;
                    }
                    else
                    {
                        G.zero();
                    }

                    // Impose the jump conditions.
                    const double x_cell_bdry =
                        x_lower[axis] + static_cast<double>(i_s(axis) - patch_lower[axis]) * dx[axis];
                    const double h = x_cell_bdry + (x(axis) > x_cell_bdry ? +0.5 : -0.5) * dx[axis] - x(axis);
                    const double C_p = F * n - h * G(axis);
                    (*f_data)(i_s) += (n(axis) > 0.0 ? +1.0 : -1.0) * (C_p / dx[axis]);
                }
            }
        }
    }
    return;
} // imposeJumpConditions

void
IBFEMethod::initializeCoordinates(const unsigned int part)
{
    EquationSystems& equation_systems = *d_fe_data_managers[part]->getEquationSystems();
    MeshBase& mesh = equation_systems.get_mesh();
    System& X_system = equation_systems.get_system(COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    const bool identity_mapping = !d_coordinate_mapping_fcn_data[part].fcn;
    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
    {
        Node* n = *it;
        if (n->n_vars(X_sys_num))
        {
            TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM);
            const libMesh::Point& X = *n;
            libMesh::Point x = X;
            if (!identity_mapping)
            {
                d_coordinate_mapping_fcn_data[part].fcn(x, X, d_coordinate_mapping_fcn_data[part].ctx);
            }
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const int dof_index = n->dof_number(X_sys_num, d, 0);
                X_coords.set(dof_index, x(d));
            }
        }
    }
    X_coords.close();
    X_system.get_dof_map().enforce_constraints_exactly(X_system, &X_coords);
    copy_and_synch(X_coords, *X_system.current_local_solution, /*close_v_in*/ false);
    return;
} // initializeCoordinates

void
IBFEMethod::updateCoordinateMapping(const unsigned int part)
{
    EquationSystems& equation_systems = *d_fe_data_managers[part]->getEquationSystems();
    MeshBase& mesh = equation_systems.get_mesh();
    System& X_system = equation_systems.get_system(COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    System& dX_system = equation_systems.get_system(COORD_MAPPING_SYSTEM_NAME);
    const unsigned int dX_sys_num = dX_system.number();
    NumericVector<double>& dX_coords = *dX_system.solution;
    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
    {
        Node* n = *it;
        if (n->n_vars(X_sys_num))
        {
            TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM);
            TBOX_ASSERT(n->n_vars(dX_sys_num) == NDIM);
            const libMesh::Point& X = *n;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const int X_dof_index = n->dof_number(X_sys_num, d, 0);
                const int dX_dof_index = n->dof_number(dX_sys_num, d, 0);
                dX_coords.set(dX_dof_index, X_coords(X_dof_index) - X(d));
            }
        }
    }
    copy_and_synch(dX_coords, *dX_system.current_local_solution);
    return;
} // updateCoordinateMapping

void
IBFEMethod::initializeVelocity(const unsigned int part)
{
    EquationSystems& equation_systems = *d_fe_data_managers[part]->getEquationSystems();
    MeshBase& mesh = equation_systems.get_mesh();
    System& U_system = equation_systems.get_system(VELOCITY_SYSTEM_NAME);
    const unsigned int U_sys_num = U_system.number();
    NumericVector<double>& U_vec = *U_system.solution;
    VectorValue<double> U;
    if (!d_initial_velocity_fcn_data[part].fcn)
    {
        U.zero();
    }
    else
    {
        for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
        {
            Node* n = *it;
            if (n->n_vars(U_sys_num))
            {
                TBOX_ASSERT(n->n_vars(U_sys_num) == NDIM);
                const libMesh::Point& X = *n;
                d_initial_velocity_fcn_data[part].fcn(U, X, d_initial_velocity_fcn_data[part].ctx);
            }
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const int dof_index = n->dof_number(U_sys_num, d, 0);
                U_vec.set(dof_index, U(d));
            }
        }
    }
    U_vec.close();
    U_system.get_dof_map().enforce_constraints_exactly(U_system, &U_vec);
    copy_and_synch(U_vec, *U_system.current_local_solution, /*close_v_in*/ false);
    return;
} // initializeVelocity

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBFEMethod::commonConstructor(const std::string& object_name,
                              Pointer<Database> input_db,
                              const std::vector<libMesh::MeshBase*>& meshes,
                              int max_level_number,
                              bool register_for_restart,
                              const std::string& restart_read_dirname,
                              unsigned int restart_restore_number)
{
    // Set the object name and register it with the restart manager.
    d_object_name = object_name;
    d_registered_for_restart = false;
    if (register_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
        d_registered_for_restart = true;
    }
    d_libmesh_restart_read_dir = restart_read_dirname;
    d_libmesh_restart_restore_number = restart_restore_number;

    // Store the mesh pointers.
    d_meshes = meshes;
    d_max_level_number = max_level_number;

    // Set some default values.
    const bool use_adaptive_quadrature = true;
    const int point_density = 2.0;
    const bool use_nodal_quadrature = false;
    const bool interp_use_consistent_mass_matrix = true;
    d_default_interp_spec = FEDataManager::InterpSpec("IB_4",
                                                      QGAUSS,
                                                      INVALID_ORDER,
                                                      use_adaptive_quadrature,
                                                      point_density,
                                                      interp_use_consistent_mass_matrix,
                                                      use_nodal_quadrature);
    d_default_spread_spec = FEDataManager::SpreadSpec(
        "IB_4", QGAUSS, INVALID_ORDER, use_adaptive_quadrature, point_density, use_nodal_quadrature);

    d_fe_family.resize(d_num_parts, INVALID_FE);
    d_fe_order.resize(d_num_parts, INVALID_ORDER);
    d_default_quad_type.resize(d_num_parts, INVALID_Q_RULE);
    d_default_quad_order.resize(d_num_parts, INVALID_ORDER);

    // Indicate that all of the parts do NOT use stress normalization by default
    // and set some default values.
    d_is_stress_normalization_part.resize(d_num_parts, false);

    // Indicate that there are no overlapping parts by default.
    d_is_overlap_velocity_part.resize(d_num_parts, false);
    d_is_overlap_velocity_master_part.resize(d_num_parts, false);
    d_overlap_velocity_master_part.resize(d_num_parts, -1);
    d_overlap_velocity_part_node_to_elem_map.resize(d_num_parts);
    d_overlap_velocity_part_ghost_idxs.resize(d_num_parts);

    d_is_overlap_force_part.resize(d_num_parts, false);
    d_overlap_force_part_ghost_idxs.resize(d_num_parts);

    // Initialize function data to NULL.
    d_coordinate_mapping_fcn_data.resize(d_num_parts);
    d_initial_velocity_fcn_data.resize(d_num_parts);
    d_PK1_stress_fcn_data.resize(d_num_parts);
    d_lag_body_force_fcn_data.resize(d_num_parts);
    d_lag_surface_pressure_fcn_data.resize(d_num_parts);
    d_lag_surface_force_fcn_data.resize(d_num_parts);
    d_lag_body_source_part.resize(d_num_parts, false);
    d_lag_body_source_fcn_data.resize(d_num_parts);
    d_direct_forcing_kinematics_data.resize(d_num_parts, Pointer<IBFEDirectForcingKinematics>(nullptr));

    // Determine whether we should use first-order or second-order shape
    // functions for each part of the structure.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        const MeshBase& mesh = *meshes[part];
        bool mesh_has_first_order_elems = false;
        bool mesh_has_second_order_elems = false;
        MeshBase::const_element_iterator el_it = mesh.elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.elements_end();
        for (; el_it != el_end; ++el_it)
        {
            const Elem* const elem = *el_it;
            mesh_has_first_order_elems = mesh_has_first_order_elems || elem->default_order() == FIRST;
            mesh_has_second_order_elems = mesh_has_second_order_elems || elem->default_order() == SECOND;
        }
        mesh_has_first_order_elems = SAMRAI_MPI::maxReduction(mesh_has_first_order_elems);
        mesh_has_second_order_elems = SAMRAI_MPI::maxReduction(mesh_has_second_order_elems);
        if ((mesh_has_first_order_elems && mesh_has_second_order_elems) ||
            (!mesh_has_first_order_elems && !mesh_has_second_order_elems))
        {
            TBOX_ERROR(d_object_name
                       << "::IBFEMethod():\n"
                       << "  each FE mesh part must contain only FIRST order elements or only SECOND order elements"
                       << std::endl);
        }
        d_fe_family[part] = LAGRANGE;
        d_default_quad_type[part] = QGAUSS;
        if (mesh_has_first_order_elems)
        {
            d_fe_order[part] = FIRST;
            d_default_quad_order[part] = THIRD;
        }
        if (mesh_has_second_order_elems)
        {
            d_fe_order[part] = SECOND;
            d_default_quad_order[part] = FIFTH;
        }

        // Report configuration.
        pout << "\n";
        pout << d_object_name << ": mesh part " << part << " is using "
             << Utility::enum_to_string<Order>(d_fe_order[part]) << " order "
             << Utility::enum_to_string<FEFamily>(d_fe_family[part]) << " finite elements.\n";
        pout << "\n";
    }

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    // Set up the interaction spec objects.
    d_interp_spec.resize(d_num_parts, d_default_interp_spec);
    d_spread_spec.resize(d_num_parts, d_default_spread_spec);

    return;
} // commonConstructor

void
IBFEMethod::getFromInput(Pointer<Database> db, bool /*is_from_restart*/)
{
    // Interpolation settings.
    if (db->isString("interp_delta_fcn"))
        d_default_interp_spec.kernel_fcn = db->getString("interp_delta_fcn");
    else if (db->isString("IB_delta_fcn"))
        d_default_interp_spec.kernel_fcn = db->getString("IB_delta_fcn");
    else if (db->isString("interp_kernel_fcn"))
        d_default_interp_spec.kernel_fcn = db->getString("interp_kernel_fcn");
    else if (db->isString("IB_kernel_fcn"))
        d_default_interp_spec.kernel_fcn = db->getString("IB_kernel_fcn");

    if (db->isString("interp_quad_type"))
        d_default_interp_spec.quad_type = Utility::string_to_enum<QuadratureType>(db->getString("interp_quad_type"));
    else if (db->isString("IB_quad_type"))
        d_default_interp_spec.quad_type = Utility::string_to_enum<QuadratureType>(db->getString("IB_quad_type"));

    if (db->isString("interp_quad_order"))
        d_default_interp_spec.quad_order = Utility::string_to_enum<Order>(db->getString("interp_quad_order"));
    else if (db->isString("IB_quad_order"))
        d_default_interp_spec.quad_order = Utility::string_to_enum<Order>(db->getString("IB_quad_order"));

    if (db->isBool("interp_use_adaptive_quadrature"))
        d_default_interp_spec.use_adaptive_quadrature = db->getBool("interp_use_adaptive_quadrature");
    else if (db->isBool("IB_use_adaptive_quadrature"))
        d_default_interp_spec.use_adaptive_quadrature = db->getBool("IB_use_adaptive_quadrature");

    if (db->isDouble("interp_point_density"))
        d_default_interp_spec.point_density = db->getDouble("interp_point_density");
    else if (db->isDouble("IB_point_density"))
        d_default_interp_spec.point_density = db->getDouble("IB_point_density");

    if (db->isBool("interp_use_consistent_mass_matrix"))
        d_default_interp_spec.use_consistent_mass_matrix = db->getBool("interp_use_consistent_mass_matrix");
    else if (db->isBool("IB_use_consistent_mass_matrix"))
        d_default_interp_spec.use_consistent_mass_matrix = db->getBool("IB_use_consistent_mass_matrix");

    if (db->isBool("interp_use_nodal_quadrature"))
        d_default_interp_spec.use_nodal_quadrature = db->getBool("interp_use_nodal_quadrature");
    else if (db->isBool("IB_use_nodal_quadrature"))
        d_default_interp_spec.use_nodal_quadrature = db->getBool("IB_use_nodal_quadrature");

    // Spreading settings.
    if (db->isString("spread_delta_fcn"))
        d_default_spread_spec.kernel_fcn = db->getString("spread_delta_fcn");
    else if (db->isString("IB_delta_fcn"))
        d_default_spread_spec.kernel_fcn = db->getString("IB_delta_fcn");
    else if (db->isString("spread_kernel_fcn"))
        d_default_spread_spec.kernel_fcn = db->getString("spread_kernel_fcn");
    else if (db->isString("IB_kernel_fcn"))
        d_default_spread_spec.kernel_fcn = db->getString("IB_kernel_fcn");

    if (db->isString("spread_quad_type"))
        d_default_spread_spec.quad_type = Utility::string_to_enum<QuadratureType>(db->getString("spread_quad_type"));
    else if (db->isString("IB_quad_type"))
        d_default_spread_spec.quad_type = Utility::string_to_enum<QuadratureType>(db->getString("IB_quad_type"));

    if (db->isString("spread_quad_order"))
        d_default_spread_spec.quad_order = Utility::string_to_enum<Order>(db->getString("spread_quad_order"));
    else if (db->isString("IB_quad_order"))
        d_default_spread_spec.quad_order = Utility::string_to_enum<Order>(db->getString("IB_quad_order"));

    if (db->isBool("spread_use_adaptive_quadrature"))
        d_default_spread_spec.use_adaptive_quadrature = db->getBool("spread_use_adaptive_quadrature");
    else if (db->isBool("IB_use_adaptive_quadrature"))
        d_default_spread_spec.use_adaptive_quadrature = db->getBool("IB_use_adaptive_quadrature");

    if (db->isDouble("spread_point_density"))
        d_default_spread_spec.point_density = db->getDouble("spread_point_density");
    else if (db->isDouble("IB_point_density"))
        d_default_spread_spec.point_density = db->getDouble("IB_point_density");

    if (db->isBool("spread_use_nodal_quadrature"))
        d_default_spread_spec.use_nodal_quadrature = db->getBool("spread_use_nodal_quadrature");
    else if (db->isBool("IB_use_nodal_quadrature"))
        d_default_spread_spec.use_nodal_quadrature = db->getBool("IB_use_nodal_quadrature");

    // Force computation settings.
    if (db->isBool("split_normal_force"))
        d_split_normal_force = db->getBool("split_normal_force");
    else if (db->isBool("split_forces"))
        d_split_normal_force = db->getBool("split_forces");
    if (db->isBool("split_tangential_force"))
        d_split_tangential_force = db->getBool("split_tangential_force");
    else if (db->isBool("split_forces"))
        d_split_tangential_force = db->getBool("split_forces");
    if (db->isBool("use_jump_conditions")) d_use_jump_conditions = db->getBool("use_jump_conditions");
    if (db->isBool("use_consistent_mass_matrix"))
        d_use_consistent_mass_matrix = db->getBool("use_consistent_mass_matrix");

    // Restart settings.
    if (db->isString("libmesh_restart_file_extension"))
    {
        d_libmesh_restart_file_extension = db->getString("libmesh_restart_file_extension");
    }
    else
    {
        d_libmesh_restart_file_extension = "xdr";
    }

    // Other settings.
    if (db->isDouble("overlap_tolerance")) d_overlap_tolerance = db->getDouble("overlap_tolerance");
    if (db->isInteger("min_ghost_cell_width"))
    {
        d_ghosts = db->getInteger("min_ghost_cell_width");
    }
    else if (db->isDouble("min_ghost_cell_width"))
    {
        d_ghosts = static_cast<int>(std::ceil(db->getDouble("min_ghost_cell_width")));
    }
    if (db->keyExists("do_log"))
        d_do_log = db->getBool("do_log");
    else if (db->keyExists("enable_logging"))
        d_do_log = db->getBool("enable_logging");

    if (db->isDouble("epsilon")) d_epsilon = db->getDouble("epsilon");
    return;
} // getFromInput

void
IBFEMethod::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to " << d_object_name
                                 << " not found in restart file." << std::endl);
    }
    int ver = db->getInteger("IBFE_METHOD_VERSION");
    if (ver != IBFE_METHOD_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    db->getIntegerArray("d_ghosts", d_ghosts, NDIM);
    d_split_normal_force = db->getBool("d_split_normal_force");
    d_split_tangential_force = db->getBool("d_split_tangential_force");
    d_use_jump_conditions = db->getBool("d_use_jump_conditions");
    d_use_consistent_mass_matrix = db->getBool("d_use_consistent_mass_matrix");
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
