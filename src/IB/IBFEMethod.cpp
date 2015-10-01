// Filename: IBFEMethod.cpp
// Created on 5 Oct 2011 by Boyce Griffith
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

#include <stdbool.h>
#include <stddef.h>
#include <algorithm>
#include <cmath>
#include <limits>
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
#include "ibtk/FEDataManager.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/libmesh_utilities.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/boundary_info.h"
#include "libmesh/compare_types.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/edge.h"
#include "libmesh/elem.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe_type.h"
#include "libmesh/fem_context.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/point.h"
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

inline short int get_dirichlet_bdry_ids(const std::vector<short int>& bdry_ids)
{
    short int dirichlet_bdry_ids = 0;
    for (std::vector<short int>::const_iterator cit = bdry_ids.begin(); cit != bdry_ids.end(); ++cit)
    {
        const short int bdry_id = *cit;
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

inline bool is_physical_bdry(const Elem* elem,
                             const unsigned short int side,
                             const BoundaryInfo& boundary_info,
                             const DofMap& dof_map)
{
    const std::vector<short int>& bdry_ids = boundary_info.boundary_ids(elem, side);
    bool at_physical_bdry = !elem->neighbor(side);
    for (std::vector<short int>::const_iterator cit = bdry_ids.begin(); cit != bdry_ids.end(); ++cit)
    {
        if (dof_map.is_periodic_boundary(*cit)) at_physical_bdry = false;
    }
    return at_physical_bdry;
}

inline bool is_dirichlet_bdry(const Elem* elem,
                              const unsigned short int side,
                              const BoundaryInfo& boundary_info,
                              const DofMap& dof_map)
{
    if (!is_physical_bdry(elem, side, boundary_info, dof_map)) return false;
    const std::vector<short int>& bdry_ids = boundary_info.boundary_ids(elem, side);
    return get_dirichlet_bdry_ids(bdry_ids) != 0;
}

inline bool has_physical_bdry(const Elem* elem, const BoundaryInfo& boundary_info, const DofMap& dof_map)
{
    bool has_physical_bdry = false;
    for (unsigned short int side = 0; side < elem->n_sides() && !has_physical_bdry; ++side)
    {
        has_physical_bdry = has_physical_bdry || is_physical_bdry(elem, side, boundary_info, dof_map);
    }
    return has_physical_bdry;
}

inline void get_x_and_FF(libMesh::VectorValue<double>& x,
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

void assemble_poisson(EquationSystems& es, const std::string& /*system_name*/)
{
    const MeshBase& mesh = es.get_mesh();
    const BoundaryInfo& boundary_info = *mesh.boundary_info;
    const unsigned int dim = mesh.mesh_dimension();
    LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>(IBFEMethod::PHI_SYSTEM_NAME);
    const DofMap& dof_map = system.get_dof_map();
    FEType fe_type = dof_map.variable_type(0);

    AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));
    QGauss qrule(dim, FIFTH);
    fe->attach_quadrature_rule(&qrule);
    AutoPtr<FEBase> fe_face(FEBase::build(dim, fe_type));
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
        unsigned int Ke_size = static_cast<unsigned int>(dof_indices.size());
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
}

const std::string IBFEMethod::BODY_VELOCITY_SYSTEM_NAME = "IB body velocity system";
const std::string IBFEMethod::COORDS_SYSTEM_NAME = "IB coordinates system";
const std::string IBFEMethod::COORD_MAPPING_SYSTEM_NAME = "IB coordinate mapping system";
const std::string IBFEMethod::FORCE_SYSTEM_NAME = "IB force system";
const std::string IBFEMethod::PHI_SYSTEM_NAME = "IB stress normalization system";
const std::string IBFEMethod::VELOCITY_SYSTEM_NAME = "IB velocity system";

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBFEMethod::IBFEMethod(const std::string& object_name,
                       Pointer<Database> input_db,
                       Mesh* mesh,
                       int max_level_number,
                       bool register_for_restart)
    : d_num_parts(1)
{
    commonConstructor(object_name, input_db, std::vector<Mesh*>(1, mesh), max_level_number, register_for_restart);
    return;
} // IBFEMethod

IBFEMethod::IBFEMethod(const std::string& object_name,
                       Pointer<Database> input_db,
                       const std::vector<Mesh*>& meshes,
                       int max_level_number,
                       bool register_for_restart)
    : d_num_parts(static_cast<int>(meshes.size()))
{
    commonConstructor(object_name, input_db, meshes, max_level_number, register_for_restart);
    return;
} // IBFEMethod

IBFEMethod::~IBFEMethod()
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        delete d_equation_systems[part];
    }
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }
    return;
} // ~IBFEMethod

FEDataManager* IBFEMethod::getFEDataManager(const unsigned int part) const
{
    TBOX_ASSERT(part < d_num_parts);
    return d_fe_data_managers[part];
} // getFEDataManager

void IBFEMethod::registerStressNormalizationPart(unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    if (d_stress_normalization_part[part]) return;
    d_has_stress_normalization_parts = true;
    d_stress_normalization_part[part] = true;
    System& Phi_system = d_equation_systems[part]->add_system<LinearImplicitSystem>(PHI_SYSTEM_NAME);
    d_equation_systems[part]->parameters.set<Real>("Phi_epsilon") = d_epsilon;
    Phi_system.attach_assemble_function(assemble_poisson);
    Phi_system.add_variable("Phi", d_fe_order, d_fe_family);
    return;
} // registerStressNormalizationPart

void IBFEMethod::registerConstrainedPart(unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    if (d_constrained_part[part]) return;
    d_has_constrained_parts = true;
    d_constrained_part[part] = true;
    System& U_b_system = d_equation_systems[part]->add_system<System>(BODY_VELOCITY_SYSTEM_NAME);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        std::ostringstream os;
        os << "U_b_" << d;
        U_b_system.add_variable(os.str(), d_fe_order, d_fe_family);
    }
    return;
} // registerConstrainedPart

void IBFEMethod::registerConstrainedVelocityFunction(ConstrainedVelocityFcnPtr fcn, void* ctx, unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    registerConstrainedVelocityFunction(ConstrainedVelocityFcnData(fcn, ctx), part);
    return;
} // registerConstrainedVelocityFunction

void IBFEMethod::registerConstrainedVelocityFunction(const ConstrainedVelocityFcnData& data, unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    registerConstrainedPart(part);
    d_constrained_velocity_fcn_data[part] = data;
    return;
} // registerConstrainedVelocityFunction

void IBFEMethod::registerInitialCoordinateMappingFunction(CoordinateMappingFcnPtr fcn,
                                                          void* ctx,
                                                          const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    registerInitialCoordinateMappingFunction(CoordinateMappingFcnData(fcn, ctx), part);
    return;
} // registerInitialCoordinateMappingFunction

void IBFEMethod::registerInitialCoordinateMappingFunction(const CoordinateMappingFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_coordinate_mapping_fcn_data[part] = data;
    return;
} // registerInitialCoordinateMappingFunction

void IBFEMethod::registerPK1StressFunction(PK1StressFcnPtr fcn,
                                           const std::vector<unsigned int>& systems,
                                           void* ctx,
                                           QuadratureType quad_type,
                                           Order quad_order,
                                           const unsigned int part)
{
    registerPK1StressFunction(PK1StressFcnData(fcn, systems, ctx, quad_type, quad_order), part);
    return;
} // registerPK1StressFunction

void IBFEMethod::registerPK1StressFunction(const PK1StressFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_PK1_stress_fcn_data[part].push_back(data);
    if (data.quad_type == INVALID_Q_RULE)
    {
        d_PK1_stress_fcn_data[part].back().quad_type = d_quad_type;
    }
    if (data.quad_order == INVALID_ORDER)
    {
        d_PK1_stress_fcn_data[part].back().quad_order = d_quad_order;
    }
    d_fcn_systems[part].insert(data.systems.begin(), data.systems.end());
    d_body_fcn_systems[part].insert(data.systems.begin(), data.systems.end());
    return;
} // registerPK1StressFunction

void IBFEMethod::registerLagBodyForceFunction(LagBodyForceFcnPtr fcn,
                                              const std::vector<unsigned int>& systems,
                                              void* ctx,
                                              const unsigned int part)
{
    registerLagBodyForceFunction(LagBodyForceFcnData(fcn, systems, ctx), part);
    return;
} // registerLagBodyForceFunction

void IBFEMethod::registerLagBodyForceFunction(const LagBodyForceFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_lag_body_force_fcn_data[part] = data;
    d_fcn_systems[part].insert(data.systems.begin(), data.systems.end());
    d_body_fcn_systems[part].insert(data.systems.begin(), data.systems.end());
    return;
} // registerLagBodyForceFunction

void IBFEMethod::registerLagSurfacePressureFunction(LagSurfacePressureFcnPtr fcn,
                                                    const std::vector<unsigned int>& systems,
                                                    void* ctx,
                                                    const unsigned int part)
{
    registerLagSurfacePressureFunction(LagSurfacePressureFcnData(fcn, systems, ctx), part);
    return;
} // registerLagSurfacePressureFunction

void IBFEMethod::registerLagSurfacePressureFunction(const LagSurfacePressureFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_lag_surface_pressure_fcn_data[part] = data;
    d_fcn_systems[part].insert(data.systems.begin(), data.systems.end());
    d_surface_fcn_systems[part].insert(data.systems.begin(), data.systems.end());
    return;
} // registerLagSurfacePressureFunction

void IBFEMethod::registerLagSurfaceForceFunction(LagSurfaceForceFcnPtr fcn,
                                                 const std::vector<unsigned int>& systems,
                                                 void* ctx,
                                                 const unsigned int part)
{
    registerLagSurfaceForceFunction(LagSurfaceForceFcnData(fcn, systems, ctx), part);
    return;
} // registerLagSurfaceForceFunction

void IBFEMethod::registerLagSurfaceForceFunction(const LagSurfaceForceFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_lag_surface_force_fcn_data[part] = data;
    d_fcn_systems[part].insert(data.systems.begin(), data.systems.end());
    d_surface_fcn_systems[part].insert(data.systems.begin(), data.systems.end());
    return;
} // registerLagSurfaceForceFunction

const IntVector<NDIM>& IBFEMethod::getMinimumGhostCellWidth() const
{
    return d_ghosts;
} // getMinimumGhostCellWidth

void IBFEMethod::setupTagBuffer(Array<int>& tag_buffer, Pointer<GriddingAlgorithm<NDIM> > gridding_alg) const
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

void IBFEMethod::preprocessIntegrateData(double current_time, double new_time, int /*num_cycles*/)
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

    d_Phi_systems.resize(d_num_parts);
    d_Phi_half_vecs.resize(d_num_parts);

    d_U_b_systems.resize(d_num_parts);
    d_U_b_current_vecs.resize(d_num_parts);
    d_U_b_new_vecs.resize(d_num_parts);
    d_U_b_half_vecs.resize(d_num_parts);

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

        if (d_stress_normalization_part[part])
        {
            d_Phi_systems[part] = &d_equation_systems[part]->get_system(PHI_SYSTEM_NAME);
            d_Phi_half_vecs[part] =
                dynamic_cast<PetscVector<double>*>(d_Phi_systems[part]->current_local_solution.get());
        }

        if (d_constrained_part[part])
        {
            d_U_b_systems[part] = &d_equation_systems[part]->get_system(BODY_VELOCITY_SYSTEM_NAME);
            d_U_b_current_vecs[part] =
                dynamic_cast<PetscVector<double>*>(d_U_b_systems[part]->current_local_solution.get());
            d_U_b_new_vecs[part] = dynamic_cast<PetscVector<double>*>(
                d_U_b_current_vecs[part]->clone().release()); // WARNING: must be manually deleted
            d_U_b_half_vecs[part] = dynamic_cast<PetscVector<double>*>(
                d_U_b_current_vecs[part]->clone().release()); // WARNING: must be manually deleted
        }

        // Initialize X^{n+1/2} and X^{n+1} to equal X^{n}, and initialize
        // U^{n+1/2} and U^{n+1} to equal U^{n}.
        d_X_systems[part]->solution->close();
        d_X_systems[part]->solution->localize(*d_X_current_vecs[part]);
        d_X_systems[part]->solution->localize(*d_X_new_vecs[part]);
        d_X_systems[part]->solution->localize(*d_X_half_vecs[part]);

        d_U_systems[part]->solution->close();
        d_U_systems[part]->solution->localize(*d_U_current_vecs[part]);
        d_U_systems[part]->solution->localize(*d_U_new_vecs[part]);
        d_U_systems[part]->solution->localize(*d_U_half_vecs[part]);

        d_F_systems[part]->solution->close();
        d_F_systems[part]->solution->localize(*d_F_half_vecs[part]);

        if (d_stress_normalization_part[part])
        {
            d_Phi_systems[part]->solution->close();
            d_Phi_systems[part]->solution->localize(*d_Phi_half_vecs[part]);
        }

        if (d_constrained_part[part])
        {
            d_U_b_systems[part]->solution->close();
            d_U_b_systems[part]->solution->localize(*d_U_b_current_vecs[part]);
            d_U_b_systems[part]->solution->localize(*d_U_b_new_vecs[part]);
            d_U_b_systems[part]->solution->localize(*d_U_b_half_vecs[part]);
        }
    }
    return;
} // preprocessIntegrateData

void IBFEMethod::postprocessIntegrateData(double /*current_time*/, double /*new_time*/, int /*num_cycles*/)
{
    for (unsigned part = 0; part < d_num_parts; ++part)
    {
        // Reset time-dependent Lagrangian data.
        d_X_new_vecs[part]->close();
        *d_X_systems[part]->solution = *d_X_new_vecs[part];
        d_X_systems[part]->solution->close();
        d_X_systems[part]->solution->localize(*d_X_systems[part]->current_local_solution);
        delete d_X_new_vecs[part];
        delete d_X_half_vecs[part];

        d_U_new_vecs[part]->close();
        *d_U_systems[part]->solution = *d_U_new_vecs[part];
        d_U_systems[part]->solution->close();
        d_U_systems[part]->solution->localize(*d_U_systems[part]->current_local_solution);
        delete d_U_new_vecs[part];
        delete d_U_half_vecs[part];

        d_F_half_vecs[part]->close();
        *d_F_systems[part]->solution = *d_F_half_vecs[part];
        d_F_systems[part]->solution->close();
        d_F_systems[part]->solution->localize(*d_F_systems[part]->current_local_solution);

        if (d_stress_normalization_part[part])
        {
            d_Phi_half_vecs[part]->close();
            *d_Phi_systems[part]->solution = *d_Phi_half_vecs[part];
            d_Phi_systems[part]->solution->close();
            d_Phi_systems[part]->solution->localize(*d_Phi_systems[part]->current_local_solution);
        }

        if (d_constrained_part[part])
        {
            d_U_b_new_vecs[part]->close();
            *d_U_b_systems[part]->solution = *d_U_b_new_vecs[part];
            d_U_b_systems[part]->solution->close();
            d_U_b_systems[part]->solution->localize(*d_U_b_systems[part]->current_local_solution);
            delete d_U_b_new_vecs[part];
            delete d_U_b_half_vecs[part];
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

    d_Phi_systems.clear();
    d_Phi_half_vecs.clear();

    d_U_b_systems.clear();
    d_U_b_current_vecs.clear();
    d_U_b_new_vecs.clear();
    d_U_b_half_vecs.clear();

    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time = std::numeric_limits<double>::quiet_NaN();
    d_half_time = std::numeric_limits<double>::quiet_NaN();
    return;
} // postprocessIntegrateData

void IBFEMethod::interpolateVelocity(const int u_data_idx,
                                     const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*u_synch_scheds*/,
                                     const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                     const double data_time)
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        NumericVector<double>* X_vec = NULL;
        NumericVector<double>* X_ghost_vec = d_X_IB_ghost_vecs[part];
        NumericVector<double>* U_vec = NULL;
        NumericVector<double>* U_b_vec = NULL;
        if (MathUtilities<double>::equalEps(data_time, d_current_time))
        {
            X_vec = d_X_current_vecs[part];
            U_vec = d_U_current_vecs[part];
            if (d_constrained_part[part]) U_b_vec = d_U_b_current_vecs[part];
        }
        else if (MathUtilities<double>::equalEps(data_time, d_half_time))
        {
            X_vec = d_X_half_vecs[part];
            U_vec = d_U_half_vecs[part];
            if (d_constrained_part[part]) U_b_vec = d_U_b_half_vecs[part];
        }
        else if (MathUtilities<double>::equalEps(data_time, d_new_time))
        {
            X_vec = d_X_new_vecs[part];
            U_vec = d_U_new_vecs[part];
            if (d_constrained_part[part]) U_b_vec = d_U_b_new_vecs[part];
        }
        X_vec->localize(*X_ghost_vec);
        if (d_use_IB_interp_operator)
        {
            d_fe_data_managers[part]->interp(u_data_idx, *U_vec, *X_ghost_vec, VELOCITY_SYSTEM_NAME,
                                             u_ghost_fill_scheds, data_time);
        }
        else
        {
            d_fe_data_managers[part]->restrictData(u_data_idx, *U_vec, *X_ghost_vec, VELOCITY_SYSTEM_NAME);
        }
        if (d_constrained_part[part] && d_constrained_velocity_fcn_data[part].fcn)
        {
            EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
            d_constrained_velocity_fcn_data[part].fcn(*U_b_vec, *U_vec, *X_vec, equation_systems, data_time,
                                                      d_constrained_velocity_fcn_data[part].ctx);
        }
    }
    return;
} // interpolateVelocity

void IBFEMethod::eulerStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    int ierr;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        PetscVector<double>* U_current_vec = NULL;
        if (d_constrained_part[part])
        {
            U_current_vec = d_U_b_current_vecs[part];
        }
        else
        {
            U_current_vec = d_U_current_vecs[part];
        }
        ierr = VecWAXPY(d_X_new_vecs[part]->vec(), dt, U_current_vec->vec(), d_X_current_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(d_X_half_vecs[part]->vec(), 0.5, 0.5, 0.0, d_X_current_vecs[part]->vec(),
                           d_X_new_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        d_X_new_vecs[part]->close();
        d_X_half_vecs[part]->close();
    }
    return;
} // eulerStep

void IBFEMethod::midpointStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    int ierr;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        PetscVector<double>* U_half_vec = NULL;
        if (d_constrained_part[part])
        {
            U_half_vec = d_U_b_half_vecs[part];
        }
        else
        {
            U_half_vec = d_U_half_vecs[part];
        }
        ierr = VecWAXPY(d_X_new_vecs[part]->vec(), dt, U_half_vec->vec(), d_X_current_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(d_X_half_vecs[part]->vec(), 0.5, 0.5, 0.0, d_X_current_vecs[part]->vec(),
                           d_X_new_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        d_X_new_vecs[part]->close();
        d_X_half_vecs[part]->close();
    }
    return;
} // midpointStep

void IBFEMethod::trapezoidalStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    int ierr;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        PetscVector<double>* U_current_vec = NULL;
        PetscVector<double>* U_new_vec = NULL;
        if (d_constrained_part[part])
        {
            U_current_vec = d_U_b_current_vecs[part];
            U_new_vec = d_U_b_half_vecs[part];
        }
        else
        {
            U_current_vec = d_U_current_vecs[part];
            U_new_vec = d_U_new_vecs[part];
        }
        ierr = VecWAXPY(d_X_new_vecs[part]->vec(), 0.5 * dt, U_current_vec->vec(), d_X_current_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPY(d_X_new_vecs[part]->vec(), 0.5 * dt, U_new_vec->vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(d_X_half_vecs[part]->vec(), 0.5, 0.5, 0.0, d_X_current_vecs[part]->vec(),
                           d_X_new_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        d_X_new_vecs[part]->close();
        d_X_half_vecs[part]->close();
    }
    return;
} // trapezoidalStep

void IBFEMethod::computeLagrangianForce(const double data_time)
{
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
    for (unsigned part = 0; part < d_num_parts; ++part)
    {
        if (d_constrained_part[part])
        {
            computeConstraintForceDensity(*d_F_half_vecs[part], *d_X_half_vecs[part], *d_U_half_vecs[part],
                                          *d_U_b_half_vecs[part], data_time, part);
        }
        else
        {
            if (d_stress_normalization_part[part])
            {
                computeStressNormalization(*d_Phi_half_vecs[part], *d_X_half_vecs[part], data_time, part);
            }
            computeInteriorForceDensity(*d_F_half_vecs[part], *d_X_half_vecs[part], d_Phi_half_vecs[part], data_time,
                                        part);
        }
    }
    return;
} // computeLagrangianForce

void IBFEMethod::spreadForce(const int f_data_idx,
                             RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                             const std::vector<Pointer<RefineSchedule<NDIM> > >& /*f_prolongation_scheds*/,
                             const double data_time)
{
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        PetscVector<double>* X_vec = d_X_half_vecs[part];
        PetscVector<double>* X_ghost_vec = d_X_IB_ghost_vecs[part];
        PetscVector<double>* F_vec = d_F_half_vecs[part];
        PetscVector<double>* F_ghost_vec = d_F_IB_ghost_vecs[part];
        X_vec->localize(*X_ghost_vec);
        F_vec->localize(*F_ghost_vec);
        if (d_use_IB_spread_operator)
        {
            d_fe_data_managers[part]->spread(f_data_idx, *F_ghost_vec, *X_ghost_vec, FORCE_SYSTEM_NAME, f_phys_bdry_op,
                                             data_time);
        }
        else
        {
            d_fe_data_managers[part]->prolongData(f_data_idx, *F_ghost_vec, *X_ghost_vec, FORCE_SYSTEM_NAME,
                                                  /*is_density*/ true,
                                                  /*accumulate_on_grid*/ true);
        }
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

void IBFEMethod::initializeFEData()
{
    if (d_fe_data_initialized) return;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        // Initialize FE equation systems.
        EquationSystems* equation_systems = d_equation_systems[part];
        equation_systems->init();
        initializeCoordinates(part);
        updateCoordinateMapping(part);

        // Assemble systems.
        System& X_system = equation_systems->get_system<System>(COORDS_SYSTEM_NAME);
        System& dX_system = equation_systems->get_system<System>(COORD_MAPPING_SYSTEM_NAME);
        System& U_system = equation_systems->get_system<System>(VELOCITY_SYSTEM_NAME);
        System& F_system = equation_systems->get_system<System>(FORCE_SYSTEM_NAME);

        X_system.assemble_before_solve = false;
        X_system.assemble();

        dX_system.assemble_before_solve = false;
        dX_system.assemble();

        U_system.assemble_before_solve = false;
        U_system.assemble();

        F_system.assemble_before_solve = false;
        F_system.assemble();

        if (d_stress_normalization_part[part])
        {
            LinearImplicitSystem& Phi_system = equation_systems->get_system<LinearImplicitSystem>(PHI_SYSTEM_NAME);
            Phi_system.assemble_before_solve = false;
            Phi_system.assemble();
        }

        if (d_constrained_part[part])
        {
            System& U_b_system = equation_systems->get_system<System>(BODY_VELOCITY_SYSTEM_NAME);
            U_b_system.assemble_before_solve = false;
            U_b_system.assemble();
        }

        // Set up boundary conditions.  Specifically, add appropriate boundary
        // IDs to the BoundaryInfo object associated with the mesh, and add DOF
        // constraints for the nodal forces and velocities.
        const MeshBase& mesh = equation_systems->get_mesh();
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
    d_fe_data_initialized = true;
    return;
} // initializeFEData

void IBFEMethod::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy,
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

void IBFEMethod::registerLoadBalancer(Pointer<LoadBalancer<NDIM> > load_balancer, int workload_data_idx)
{
    TBOX_ASSERT(load_balancer);
    d_load_balancer = load_balancer;
    d_workload_idx = workload_data_idx;

    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->registerLoadBalancer(load_balancer, workload_data_idx);
    }
    return;
} // registerLoadBalancer

void IBFEMethod::updateWorkloadEstimates(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/, int /*workload_data_idx*/)
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

void IBFEMethod::endDataRedistribution(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                       Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    if (d_is_initialized)
    {
        for (unsigned int part = 0; part < d_num_parts; ++part)
        {
            d_fe_data_managers[part]->reinitElementMappings();
        }
    }
    return;
} // endDataRedistribution

void IBFEMethod::initializeLevelData(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
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
        d_fe_data_managers[part]->initializeLevelData(hierarchy, level_number, init_data_time, can_be_refined,
                                                      initial_time, old_level, allocate_data);
        if (d_load_balancer && level_number == d_fe_data_managers[part]->getLevelNumber())
        {
            d_load_balancer->setWorkloadPatchDataIndex(d_workload_idx, level_number);
            d_fe_data_managers[part]->updateWorkloadEstimates(level_number, level_number);
        }
    }
    return;
} // initializeLevelData

void IBFEMethod::resetHierarchyConfiguration(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
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

void IBFEMethod::applyGradientDetector(Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
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
        d_fe_data_managers[part]->applyGradientDetector(hierarchy, level_number, error_data_time, tag_index,
                                                        initial_time, uses_richardson_extrapolation_too);
    }
    return;
} // applyGradientDetector

void IBFEMethod::putToDatabase(Pointer<Database> db)
{
    db->putInteger("IBFE_METHOD_VERSION", IBFE_METHOD_VERSION);
    db->putIntegerArray("d_ghosts", d_ghosts, NDIM);
    db->putBool("d_split_normal_force", d_split_normal_force);
    db->putBool("d_split_tangential_force", d_split_tangential_force);
    db->putBool("d_use_jump_conditions", d_use_jump_conditions);
    db->putString("d_fe_family", Utility::enum_to_string<FEFamily>(d_fe_family));
    db->putString("d_fe_order", Utility::enum_to_string<Order>(d_fe_order));
    db->putString("d_quad_type", Utility::enum_to_string<QuadratureType>(d_quad_type));
    db->putString("d_quad_order", Utility::enum_to_string<Order>(d_quad_order));
    db->putBool("d_use_consistent_mass_matrix", d_use_consistent_mass_matrix);
    return;
} // putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

void IBFEMethod::computeConstraintForceDensity(PetscVector<double>& F_vec,
                                               PetscVector<double>& /*X_vec*/,
                                               PetscVector<double>& U_vec,
                                               PetscVector<double>& U_b_vec,
                                               const double /*data_time*/,
                                               const unsigned int part)
{
    if (!d_constrained_part[part]) return;

    const double dt = d_new_time - d_current_time;
    const double rho = getINSHierarchyIntegrator()->getStokesSpecifications()->getRho();
    int ierr = VecAXPBYPCZ(F_vec.vec(), d_constraint_omega * rho / dt, -d_constraint_omega * rho / dt, 0.0,
                           U_b_vec.vec(), U_vec.vec());
    IBTK_CHKERRQ(ierr);
    F_vec.close();
    return;
} // computeConstraintForceDensity

void IBFEMethod::setupStressAndForceSystemData(std::vector<std::vector<NumericVector<double>*> >* PK1_stress_fcn_data,
                                               std::vector<NumericVector<double>*>* lag_body_force_fcn_data,
                                               std::vector<NumericVector<double>*>* lag_surface_force_fcn_data,
                                               std::vector<NumericVector<double>*>* lag_surface_pressure_fcn_data,
                                               const unsigned int part)
{
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    for (std::set<unsigned int>::const_iterator cit = d_fcn_systems[part].begin(); cit != d_fcn_systems[part].end();
         ++cit)
    {
        System& system = equation_systems->get_system(*cit);
        system.update();
    }

    if (PK1_stress_fcn_data)
    {
        const size_t num_PK1_stress_fcns = d_PK1_stress_fcn_data[part].size();
        PK1_stress_fcn_data->resize(num_PK1_stress_fcns);
        for (unsigned int k = 0; k < num_PK1_stress_fcns; ++k)
        {
            (*PK1_stress_fcn_data)[k].clear();
            std::vector<unsigned int>& PK1_stress_fcn_systems = d_PK1_stress_fcn_data[part][k].systems;
            for (std::vector<unsigned int>::const_iterator cit = PK1_stress_fcn_systems.begin();
                 cit != PK1_stress_fcn_systems.end(); ++cit)
            {
                System& system = equation_systems->get_system(*cit);
                (*PK1_stress_fcn_data)[k].push_back(system.current_local_solution.get());
            }
        }
    }

    if (lag_body_force_fcn_data)
    {
        lag_body_force_fcn_data->clear();
        std::vector<unsigned int>& lag_body_force_fcn_systems = d_lag_body_force_fcn_data[part].systems;
        for (std::vector<unsigned int>::const_iterator cit = lag_body_force_fcn_systems.begin();
             cit != lag_body_force_fcn_systems.end(); ++cit)
        {
            System& system = equation_systems->get_system(*cit);
            lag_body_force_fcn_data->push_back(system.current_local_solution.get());
        }
    }

    if (lag_surface_force_fcn_data)
    {
        lag_surface_force_fcn_data->clear();
        std::vector<unsigned int>& lag_surface_force_fcn_systems = d_lag_surface_force_fcn_data[part].systems;
        for (std::vector<unsigned int>::const_iterator cit = lag_surface_force_fcn_systems.begin();
             cit != lag_surface_force_fcn_systems.end(); ++cit)
        {
            System& system = equation_systems->get_system(*cit);
            lag_surface_force_fcn_data->push_back(system.current_local_solution.get());
        }
    }

    if (lag_surface_pressure_fcn_data)
    {
        lag_surface_pressure_fcn_data->clear();
        std::vector<unsigned int>& lag_surface_pressure_fcn_systems = d_lag_surface_pressure_fcn_data[part].systems;
        for (std::vector<unsigned int>::const_iterator cit = lag_surface_pressure_fcn_systems.begin();
             cit != lag_surface_pressure_fcn_systems.end(); ++cit)
        {
            System& system = equation_systems->get_system(*cit);
            lag_surface_pressure_fcn_data->push_back(system.current_local_solution.get());
        }
    }
    return;
}

class FEFunctionInterpolation
{
public:
    FEFunctionInterpolation(const unsigned int dim = NDIM)
        : d_dim(dim), d_initialized(false), d_eval_q_point(false), d_eval_JxW(false), d_eval_q_point_face(false),
          d_eval_JxW_face(false), d_eval_normal_face(false), d_qrule(NULL), d_qrule_face(NULL), d_q_point(NULL),
          d_q_point_face(NULL), d_JxW(NULL), d_JxW_face(NULL), d_normal_face(NULL), d_current_elem(NULL)
    {
        return;
    }

    ~FEFunctionInterpolation()
    {
        return;
    }

    void attachQuadratureRule(QBase* qrule)
    {
        d_qrule = qrule;
        for (size_t k = 0; k < d_fe.size(); ++k)
        {
            Pointer<FEBase> fe = d_fe[k];
            if (fe)
            {
                fe->attach_quadrature_rule(d_qrule);
            }
        }
        return;
    }

    void attachQuadratureRuleFace(QBase* qrule_face)
    {
        d_qrule_face = qrule_face;
        for (size_t k = 0; k < d_fe_face.size(); ++k)
        {
            Pointer<FEBase> fe_face = d_fe_face[k];
            if (fe_face)
            {
                fe_face->attach_quadrature_rule(d_qrule_face);
            }
        }
        return;
    }

    void evalQuadraturePoints()
    {
        TBOX_ASSERT(!d_initialized);
        d_eval_q_point = true;
        return;
    }

    void evalQuadratureWeights()
    {
        TBOX_ASSERT(!d_initialized);
        d_eval_JxW = true;
        return;
    }

    void evalQuadraturePointsFace()
    {
        TBOX_ASSERT(!d_initialized);
        d_eval_q_point_face = true;
        return;
    }

    void evalQuadratureWeightsFace()
    {
        TBOX_ASSERT(!d_initialized);
        d_eval_JxW_face = true;
        return;
    }

    void evalNormalsFace()
    {
        TBOX_ASSERT(!d_initialized);
        d_eval_normal_face = true;
        return;
    }

    const QBase* getQrule() const
    {
        return d_qrule;
    }

    const QBase* getQruleFace() const
    {
        return d_qrule_face;
    }

    const std::vector<libMesh::Point>& getQuadraturePoints() const
    {
        TBOX_ASSERT(d_initialized);
        TBOX_ASSERT(d_eval_q_point);
        return *d_q_point;
    }

    const std::vector<double>& getQuadratureWeights() const
    {
        TBOX_ASSERT(d_initialized);
        TBOX_ASSERT(d_eval_JxW);
        return *d_JxW;
    }

    const std::vector<libMesh::Point>& getQuadraturePointsFace() const
    {
        TBOX_ASSERT(d_initialized);
        TBOX_ASSERT(d_eval_q_point_face);
        return *d_q_point_face;
    }

    const std::vector<double>& getQuadratureWeightsFace() const
    {
        TBOX_ASSERT(d_initialized);
        TBOX_ASSERT(d_eval_JxW_face);
        return *d_JxW_face;
    }

    const std::vector<libMesh::Point>& getNormalsFace() const
    {
        TBOX_ASSERT(d_initialized);
        TBOX_ASSERT(d_eval_normal_face);
        return *d_normal_face;
    }

    const std::vector<std::vector<double> >& getPhi(const FEType& fe_type) const
    {
        TBOX_ASSERT(d_initialized);
        const size_t fe_type_idx = getFETypeIndex(fe_type);
        TBOX_ASSERT(fe_type_idx < d_fe_types.size());
        TBOX_ASSERT(d_eval_phi[fe_type_idx]);
        return *d_phi[fe_type_idx];
    }

    const std::vector<std::vector<VectorValue<double> > >& getDphi(const FEType& fe_type) const
    {
        TBOX_ASSERT(d_initialized);
        const size_t fe_type_idx = getFETypeIndex(fe_type);
        TBOX_ASSERT(fe_type_idx < d_fe_types.size());
        TBOX_ASSERT(d_eval_dphi[fe_type_idx]);
        return *d_dphi[fe_type_idx];
    }

    const std::vector<std::vector<double> >& getPhiFace(const FEType& fe_type) const
    {
        TBOX_ASSERT(d_initialized);
        const size_t fe_type_idx = getFETypeIndex(fe_type);
        TBOX_ASSERT(fe_type_idx < d_fe_types.size());
        TBOX_ASSERT(d_eval_phi[fe_type_idx]);
        return *d_phi_face[fe_type_idx];
    }

    const std::vector<std::vector<VectorValue<double> > >& getDphiFace(const FEType& fe_type) const
    {
        TBOX_ASSERT(d_initialized);
        const size_t fe_type_idx = getFETypeIndex(fe_type);
        TBOX_ASSERT(fe_type_idx < d_fe_types.size());
        TBOX_ASSERT(d_eval_dphi[fe_type_idx]);
        return *d_dphi_face[fe_type_idx];
    }

    // Configure the class to evaluate the requested shape functions / derivatives for the given system.
    //
    // NOTE: To interpolate variables associated with the system, use registerInterpolatedSystem().
    void registerSystem(const System& system,
                        const std::vector<int>& phi_vars = std::vector<int>(1, 0),
                        const std::vector<int>& dphi_vars = std::vector<int>(1, 0))
    {
        TBOX_ASSERT(!d_initialized && (!phi_vars.empty() || !dphi_vars.empty()));
        const unsigned int sys_num = system.number();
        for (std::vector<const System *>::iterator it = d_noninterp_systems.begin(), it_end = d_noninterp_systems.end();
             it != it_end; ++it)
        {
            if ((*it)->number() == sys_num)
            {
                // This system has already been registered.  If so, just merge the collection of variables.
                size_t system_idx = std::distance(d_noninterp_systems.begin(), it);

                std::vector<int>& orig_all_vars = d_noninterp_system_all_vars[system_idx];
                std::set<int> all_vars_set(orig_all_vars.begin(), orig_all_vars.end());
                all_vars_set.insert(phi_vars.begin(), phi_vars.end());
                all_vars_set.insert(dphi_vars.begin(), dphi_vars.end());
                d_noninterp_system_all_vars[system_idx].assign(all_vars_set.begin(), all_vars_set.end());

                std::vector<int>& orig_phi_vars = d_noninterp_system_phi_vars[system_idx];
                std::set<int> phi_vars_set(orig_phi_vars.begin(), orig_phi_vars.end());
                phi_vars_set.insert(phi_vars.begin(), phi_vars.end());
                d_noninterp_system_phi_vars[system_idx].assign(phi_vars_set.begin(), phi_vars_set.end());

                std::vector<int>& orig_dphi_vars = d_noninterp_system_dphi_vars[system_idx];
                std::set<int> dphi_vars_set(orig_dphi_vars.begin(), orig_dphi_vars.end());
                dphi_vars_set.insert(dphi_vars.begin(), dphi_vars.end());
                d_noninterp_system_dphi_vars[system_idx].assign(dphi_vars_set.begin(), dphi_vars_set.end());
                return;
            }
        }

        // We have not previously registered the system, so we need to register it here.
        d_noninterp_systems.push_back(&system);

        std::set<int> all_vars_set;
        all_vars_set.insert(phi_vars.begin(), phi_vars.end());
        all_vars_set.insert(dphi_vars.begin(), dphi_vars.end());
        d_noninterp_system_all_vars.push_back(std::vector<int>(all_vars_set.begin(), all_vars_set.end()));

        std::set<int> phi_vars_set(phi_vars.begin(), phi_vars.end());
        d_noninterp_system_phi_vars.push_back(std::vector<int>(phi_vars_set.begin(), phi_vars_set.end()));

        std::set<int> dphi_vars_set(dphi_vars.begin(), dphi_vars.end());
        d_noninterp_system_dphi_vars.push_back(std::vector<int>(dphi_vars_set.begin(), dphi_vars_set.end()));
        return;
    }

    // Configure the class to interpolate the requested variables / gradients for the given system.
    //
    // If the system data vector is NULL, then this class will use system.current_local_solution.
    //
    // NOTE: The same system can be registered multiple times with different sets of variables/gradients and system data
    // vectors.
    //
    // \return Returns an index that is used to access the interpolated data.
    size_t registerInterpolatedSystem(const System& system,
                                      const std::vector<int>& vars = std::vector<int>(1, 0),
                                      const std::vector<int>& grad_vars = std::vector<int>(),
                                      NumericVector<double>* system_data = NULL)
    {
        TBOX_ASSERT(!d_initialized && (!vars.empty() || !grad_vars.empty()));
        const unsigned int sys_num = system.number();
        for (std::vector<const System *>::iterator it = d_systems.begin(), it_end = d_systems.end(); it != it_end; ++it)
        {
            if ((*it)->number() == sys_num)
            {
                // This system has already been registered.  Check to see if the same collection of variables (etc.)
                // were used in the previous registration action; if so, do not re-register the system.
                //
                // NOTE: The same system may be registered multiple times with different collections of variables,
                // gradient variables, system data, etc.  We want to make sure we check all of the registered systems to
                // see if any of them match the function arguments.
                const size_t system_idx = std::distance(d_systems.begin(), it);
                bool same_data = true;
                same_data = same_data && (vars == d_system_vars[system_idx]);
                same_data = same_data && (grad_vars == d_system_grad_vars[system_idx]);
                same_data = same_data && (system_data == d_system_data[system_idx]);
                if (same_data) return system_idx;
            }
        }

        // Either we have not previously registered this system, or we are registering it here with a different
        // collection of variables/data.  In either case, we need to register it here.
        const size_t system_idx = d_systems.size();
        d_systems.push_back(&system);
        std::set<int> all_vars_set;
        all_vars_set.insert(vars.begin(), vars.end());
        all_vars_set.insert(grad_vars.begin(), grad_vars.end());
        std::vector<int> all_vars(all_vars_set.begin(), all_vars_set.end());
        d_system_all_vars.push_back(all_vars);
        d_system_vars.push_back(vars);
        std::vector<size_t> var_idx(vars.size());
        for (size_t k = 0; k < vars.size(); ++k)
        {
            var_idx[k] = std::distance(all_vars.begin(), std::find(all_vars.begin(), all_vars.end(), vars[k]));
        }
        d_system_var_idx.push_back(var_idx);
        d_system_grad_vars.push_back(grad_vars);
        std::vector<size_t> grad_var_idx(grad_vars.size());
        for (size_t k = 0; k < grad_vars.size(); ++k)
        {
            grad_var_idx[k] =
                std::distance(all_vars.begin(), std::find(all_vars.begin(), all_vars.end(), grad_vars[k]));
        }
        d_system_grad_var_idx.push_back(grad_var_idx);
        d_system_data.push_back(system_data);
        d_system_var_data.push_back(std::vector<std::vector<double> >(vars.size()));
        d_system_grad_var_data.push_back(std::vector<std::vector<VectorValue<double> > >(grad_vars.size()));
        return system_idx;
    }

    const std::vector<std::vector<double> >& getVarInterpolation(const size_t system_idx)
    {
        TBOX_ASSERT(system_idx < d_systems.size());
        return d_system_var_data[system_idx];
    }

    const std::vector<std::vector<VectorValue<double> > >& getGradVarInterpolation(const size_t system_idx)
    {
        TBOX_ASSERT(system_idx < d_systems.size());
        return d_system_grad_var_data[system_idx];
    }

    void init()
    {
        TBOX_ASSERT(!d_initialized);

        // Collect the distinct FETypes to be used.
        std::set<FEType> fe_type_set;
        const size_t num_systems = d_systems.size();
        for (size_t system_idx = 0; system_idx < num_systems; ++system_idx)
        {
            const System& system = *d_systems[system_idx];
            const DofMap& system_dof_map = system.get_dof_map();
            const std::vector<int>& all_vars = d_system_all_vars[system_idx];
            for (unsigned int k = 0; k < all_vars.size(); ++k)
            {
                fe_type_set.insert(system_dof_map.variable_type(all_vars[k]));
            }
        }
        const size_t num_noninterp_systems = d_noninterp_systems.size();
        for (size_t system_idx = 0; system_idx < num_noninterp_systems; ++system_idx)
        {
            const System& system = *d_noninterp_systems[system_idx];
            const DofMap& system_dof_map = system.get_dof_map();
            const std::vector<int>& all_vars = d_noninterp_system_all_vars[system_idx];
            for (unsigned int k = 0; k < all_vars.size(); ++k)
            {
                fe_type_set.insert(system_dof_map.variable_type(all_vars[k]));
            }
        }
        d_fe_types.assign(fe_type_set.begin(), fe_type_set.end());

        // Determine which FETypes are used for which variables, and which shape functions and shape function
        // derivatives need to be computed.
        const size_t num_fe_types = d_fe_types.size();
        d_eval_phi.resize(num_fe_types, false);
        d_eval_dphi.resize(num_fe_types, false);
        d_system_var_fe_type_idx.resize(num_systems);
        d_system_grad_var_fe_type_idx.resize(num_systems);
        for (size_t system_idx = 0; system_idx < num_systems; ++system_idx)
        {
            const System& system = *d_systems[system_idx];
            const DofMap& system_dof_map = system.get_dof_map();

            const std::vector<int>& vars = d_system_vars[system_idx];
            d_system_var_fe_type_idx[system_idx].resize(vars.size(), -1);
            for (unsigned int k = 0; k < vars.size(); ++k)
            {
                const FEType& fe_type = system_dof_map.variable_type(vars[k]);
                const size_t fe_type_idx = getFETypeIndex(fe_type);
                d_eval_phi[fe_type_idx] = true;
                d_system_var_fe_type_idx[system_idx][k] = fe_type_idx;
            }

            const std::vector<int>& grad_vars = d_system_grad_vars[system_idx];
            d_system_grad_var_fe_type_idx[system_idx].resize(grad_vars.size(), -1);
            for (unsigned int k = 0; k < grad_vars.size(); ++k)
            {
                const FEType& fe_type = system_dof_map.variable_type(grad_vars[k]);
                const size_t fe_type_idx = getFETypeIndex(fe_type);
                d_eval_dphi[fe_type_idx] = true;
                d_system_grad_var_fe_type_idx[system_idx][k] = fe_type_idx;
            }
        }
        for (size_t system_idx = 0; system_idx < num_noninterp_systems; ++system_idx)
        {
            const System& system = *d_noninterp_systems[system_idx];
            const DofMap& system_dof_map = system.get_dof_map();

            const std::vector<int>& phi_vars = d_noninterp_system_phi_vars[system_idx];
            for (unsigned int k = 0; k < phi_vars.size(); ++k)
            {
                const FEType& fe_type = system_dof_map.variable_type(phi_vars[k]);
                const size_t fe_type_idx = getFETypeIndex(fe_type);
                d_eval_phi[fe_type_idx] = true;
            }

            const std::vector<int>& dphi_vars = d_noninterp_system_dphi_vars[system_idx];
            for (unsigned int k = 0; k < dphi_vars.size(); ++k)
            {
                const FEType& fe_type = system_dof_map.variable_type(dphi_vars[k]);
                const size_t fe_type_idx = getFETypeIndex(fe_type);
                d_eval_dphi[fe_type_idx] = true;
            }
        }

        // Set up FE objects and request access to shape functions / gradients, as needed.
        d_fe.resize(num_fe_types);
        d_fe_face.resize(num_fe_types);
        d_phi.resize(num_fe_types, NULL);
        d_dphi.resize(num_fe_types, NULL);
        d_phi_face.resize(num_fe_types, NULL);
        d_dphi_face.resize(num_fe_types, NULL);
        for (unsigned int fe_type_idx = 0; fe_type_idx < num_fe_types; ++fe_type_idx)
        {
            Pointer<FEBase>& fe = d_fe[fe_type_idx];
            if (!fe)
            {
                const FEType& fe_type = d_fe_types[fe_type_idx];
                fe = FEBase::build(d_dim, fe_type).release();
                if (d_qrule) fe->attach_quadrature_rule(d_qrule);
                if (d_eval_q_point && !d_q_point) d_q_point = &fe->get_xyz();
                if (d_eval_JxW && !d_JxW) d_JxW = &fe->get_JxW();
                if (d_eval_phi[fe_type_idx]) d_phi[fe_type_idx] = &fe->get_phi();
                if (d_eval_dphi[fe_type_idx]) d_dphi[fe_type_idx] = &fe->get_dphi();
            }

            Pointer<FEBase>& fe_face = d_fe_face[fe_type_idx];
            if (!fe_face)
            {
                const FEType& fe_type = d_fe_types[fe_type_idx];
                fe_face = FEBase::build(d_dim, fe_type).release();
                if (d_qrule_face) fe_face->attach_quadrature_rule(d_qrule_face);
                if (d_eval_q_point_face && !d_q_point_face) d_q_point_face = &fe_face->get_xyz();
                if (d_eval_JxW_face && !d_JxW_face) d_JxW_face = &fe_face->get_JxW();
                if (d_eval_normal_face && !d_normal_face) d_normal_face = &fe_face->get_normals();
                if (d_eval_phi[fe_type_idx]) d_phi_face[fe_type_idx] = &fe_face->get_phi();
                if (d_eval_dphi[fe_type_idx]) d_dphi_face[fe_type_idx] = &fe_face->get_dphi();
            }
        }

        // Indicate that we have initialized the class.
        d_initialized = true;
        return;
    }

    void reinit(const Elem* elem)
    {
        TBOX_ASSERT(d_initialized);
        d_current_elem = elem;
        const size_t num_fe_types = d_fe_types.size();
        for (unsigned int fe_type_idx = 0; fe_type_idx < num_fe_types; ++fe_type_idx)
        {
            Pointer<FEBase>& fe = d_fe[fe_type_idx];
            fe->reinit(elem);
        }
        return;
    }

    void reinit(const Elem* const elem, const unsigned int side)
    {
        TBOX_ASSERT(d_initialized);
        d_current_elem = elem;
        d_current_side = side;
        const size_t num_fe_types = d_fe_types.size();
        for (unsigned int fe_type_idx = 0; fe_type_idx < num_fe_types; ++fe_type_idx)
        {
            Pointer<FEBase>& fe_face = d_fe_face[fe_type_idx];
            fe_face->reinit(elem, side);
        }
        return;
    }

    void collectDataForInterpolation(const Elem* const elem)
    {
        TBOX_ASSERT(d_initialized);
        TBOX_ASSERT(elem == d_current_elem);

        // Collect local DOF data for the element.
        const size_t num_systems = d_systems.size();
        d_system_elem_data.resize(num_systems);
        for (size_t system_idx = 0; system_idx < num_systems; ++system_idx)
        {
            const System& system = *d_systems[system_idx];
            const DofMap& system_dof_map = system.get_dof_map();
            const std::vector<int>& all_vars = d_system_all_vars[system_idx];
            const size_t num_vars = all_vars.size();

            // Get the DOF mappings and local data for all variables.
            std::vector<std::vector<unsigned int> > dof_indices(num_vars);
            NumericVector<double>* system_data = d_system_data[system_idx];
            if (!system_data) system_data = system.current_local_solution.get();
            for (size_t k = 0; k < num_vars; ++k)
            {
                system_dof_map.dof_indices(d_current_elem, dof_indices[k], all_vars[k]);
            }
            boost::multi_array<double, 2>& elem_data = d_system_elem_data[system_idx];
            get_values_for_interpolation(elem_data, *system_data, dof_indices);
        }
        return;
    }

    // NOTE: The data to be interpolated are set by calling reinit().
    void interpolate(const Elem* const elem)
    {
        TBOX_ASSERT(d_initialized);
        TBOX_ASSERT(d_qrule);
        TBOX_ASSERT(elem == d_current_elem);
        interpolateCommon(d_system_var_data, d_system_grad_var_data, d_qrule, d_phi, d_dphi);
        return;
    }

    // NOTE: The data to be interpolated are set by calling reinit().
    void interpolate(const Elem* const elem, const int side)
    {
        TBOX_ASSERT(d_initialized);
        TBOX_ASSERT(d_qrule_face);
        TBOX_ASSERT(elem == d_current_elem);
        TBOX_ASSERT(side == d_current_side);
        interpolateCommon(d_system_var_data, d_system_grad_var_data, d_qrule_face, d_phi_face, d_dphi_face);
        return;
    }

private:
    FEFunctionInterpolation(const FEFunctionInterpolation&);
    FEFunctionInterpolation& operator=(const FEFunctionInterpolation&);

    size_t getFETypeIndex(const FEType& fe_type) const
    {
        return std::distance(d_fe_types.begin(), std::find(d_fe_types.begin(), d_fe_types.end(), fe_type));
    }

    void interpolateCommon(std::vector<std::vector<std::vector<double> > >& system_var_data,
                           std::vector<std::vector<std::vector<VectorValue<double> > > >& system_grad_var_data,
                           const QBase* const qrule,
                           const std::vector<const std::vector<std::vector<double> >*>& phi_data,
                           const std::vector<const std::vector<std::vector<VectorValue<double> > >*>& dphi_data)
    {
        // Determine the number of quadrature points where the interpolation will be evaluated.
        const unsigned int n_qp = qrule->n_points();

        // Interpolate data for each system.
        for (size_t system_idx = 0; system_idx < d_systems.size(); ++system_idx)
        {
            const boost::multi_array<double, 2>& elem_data = d_system_elem_data[system_idx];

            // Interpolate regular variables.
            {
                const std::vector<size_t>& var_idxs = d_system_var_idx[system_idx];
                const unsigned int num_vars = static_cast<unsigned int>(var_idxs.size());
                std::vector<std::vector<double> >& var_data = system_var_data[system_idx];
                var_data.resize(n_qp);
                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    var_data[qp].resize(num_vars);
                }
                const std::vector<size_t>& var_fe_type_idxs = d_system_var_fe_type_idx[system_idx];
                for (unsigned int k = 0; k < num_vars; ++k)
                {
                    const size_t var_idx = var_idxs[k];
                    const size_t fe_type_idx = var_fe_type_idxs[k];
                    const std::vector<std::vector<double> >& phi = *phi_data[fe_type_idx];
                    for (unsigned int qp = 0; qp < n_qp; ++qp)
                    {
                        var_data[qp][k] = 0.0;
                        for (unsigned int i = 0; i < phi.size(); ++i)
                        {
                            var_data[qp][k] += elem_data[i][var_idx] * phi[i][qp];
                        }
                    }
                }
            }

            // Interpolate gradient variables.
            {
                const std::vector<size_t>& grad_var_idxs = d_system_grad_var_idx[system_idx];
                const unsigned int num_grad_vars = static_cast<unsigned int>(grad_var_idxs.size());
                std::vector<std::vector<VectorValue<double> > >& grad_var_data = d_system_grad_var_data[system_idx];
                grad_var_data.resize(n_qp);
                for (unsigned int qp = 0; qp < n_qp; ++qp)
                {
                    grad_var_data[qp].resize(num_grad_vars);
                }
                const std::vector<size_t>& grad_var_fe_type_idxs = d_system_grad_var_fe_type_idx[system_idx];
                for (unsigned int k = 0; k < num_grad_vars; ++k)
                {
                    const size_t var_idx = grad_var_idxs[k];
                    const size_t fe_type_idx = grad_var_fe_type_idxs[k];
                    const std::vector<std::vector<VectorValue<double> > >& dphi = *dphi_data[fe_type_idx];
                    for (unsigned int qp = 0; qp < n_qp; ++qp)
                    {
                        grad_var_data[qp][k].zero();
                        for (unsigned int i = 0; i < dphi.size(); ++i)
                        {
                            grad_var_data[qp][k] += elem_data[i][var_idx] * dphi[i][qp];
                        }
                    }
                }
            }
        }
        return;
    }

    const unsigned int d_dim;
    bool d_initialized;
    bool d_eval_q_point, d_eval_JxW, d_eval_q_point_face, d_eval_JxW_face, d_eval_normal_face;
    QBase *d_qrule, *d_qrule_face;
    const std::vector<libMesh::Point> *d_q_point, *d_q_point_face;
    const std::vector<double> *d_JxW, *d_JxW_face;
    const std::vector<libMesh::Point>* d_normal_face;

    // Data associated with systems.
    std::vector<const System*> d_systems;
    std::vector<std::vector<int> > d_system_all_vars, d_system_vars, d_system_grad_vars;
    std::vector<std::vector<size_t> > d_system_var_idx, d_system_grad_var_idx;
    std::vector<NumericVector<double>*> d_system_data;
    std::vector<std::vector<size_t> > d_system_var_fe_type_idx, d_system_grad_var_fe_type_idx;
    std::vector<std::vector<std::vector<double> > > d_system_var_data;
    std::vector<std::vector<std::vector<VectorValue<double> > > > d_system_grad_var_data;
    std::vector<const System*> d_noninterp_systems;
    std::vector<std::vector<int> > d_noninterp_system_all_vars, d_noninterp_system_phi_vars,
        d_noninterp_system_dphi_vars;

    // Data associated with FETypes.
    std::vector<FEType> d_fe_types;
    std::vector<Pointer<FEBase> > d_fe, d_fe_face;
    std::vector<bool> d_eval_phi, d_eval_dphi;
    std::vector<const std::vector<std::vector<double> > *> d_phi, d_phi_face;
    std::vector<const std::vector<std::vector<VectorValue<double> > > *> d_dphi, d_dphi_face;

    // Data associated with the current element.
    const Elem* d_current_elem;
    unsigned int d_current_side;
    std::vector<boost::multi_array<double, 2> > d_system_elem_data;
};

void IBFEMethod::computeStressNormalization(PetscVector<double>& Phi_vec,
                                            PetscVector<double>& x_vec,
                                            const double data_time,
                                            const unsigned int part)
{
    if (d_constrained_part[part]) return;

    // Extract the mesh.
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems->get_mesh();
    const BoundaryInfo& boundary_info = *mesh.boundary_info;
    const unsigned int dim = mesh.mesh_dimension();

    // Setup extra data needed to compute stresses/forces.
    const size_t num_PK1_stress_fcns = d_PK1_stress_fcn_data[part].size();
    std::vector<std::vector<DenseVector<double> > > PK1_stress_fcn_data(num_PK1_stress_fcns);
    std::vector<DenseVector<double> > lag_surface_force_fcn_data, lag_surface_pressure_fcn_data;

    // Extract the FE systems and DOF maps, and setup the FE objects.
    LinearImplicitSystem& Phi_system = equation_systems->get_system<LinearImplicitSystem>(PHI_SYSTEM_NAME);
    const DofMap& Phi_dof_map = Phi_system.get_dof_map();
    std::vector<unsigned int> Phi_dof_indices;
    FEType Phi_fe_type = Phi_dof_map.variable_type(0);
    std::vector<int> Phi_vars(1, 0);

    System& x_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    std::vector<int> x_vars(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) x_vars[d] = d;

    FEFunctionInterpolation fe(dim);
    AutoPtr<QBase> qrule_face = QBase::build(QGAUSS, dim - 1, FIFTH);
    fe.attachQuadratureRuleFace(qrule_face.get());
    fe.evalNormalsFace();
    fe.evalQuadraturePointsFace();
    fe.evalQuadratureWeightsFace();
    fe.registerSystem(Phi_system, Phi_vars, Phi_vars); // compute phi and dphi for the Phi system
    const size_t x_sys_idx = fe.registerInterpolatedSystem(x_system, x_vars, x_vars, &x_vec);
    fe.init();

    const std::vector<libMesh::Point>& q_point_face = fe.getQuadraturePointsFace();
    const std::vector<double>& JxW_face = fe.getQuadratureWeightsFace();
    const std::vector<libMesh::Point>& normal_face = fe.getNormalsFace();
    const std::vector<std::vector<double> >& phi_face = fe.getPhiFace(Phi_fe_type);

    const std::vector<std::vector<double> >& x_data = fe.getVarInterpolation(x_sys_idx);
    const std::vector<std::vector<VectorValue<double> > >& grad_x_data = fe.getGradVarInterpolation(x_sys_idx);

    // Setup global and elemental right-hand-side vectors.
    NumericVector<double>* Phi_rhs_vec = Phi_system.rhs;
    Phi_rhs_vec->zero();
    Phi_rhs_vec->close();
    DenseVector<double> Phi_rhs_e;

    // Set up boundary conditions for Phi.
    TensorValue<double> PP, FF, FF_trans, FF_inv_trans;
    VectorValue<double> F, F_s, F_qp, n, x;
    double P;
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
                Phi_dof_map.dof_indices(elem, Phi_dof_indices);
                Phi_rhs_e.resize(static_cast<int>(Phi_dof_indices.size()));
                fe.collectDataForInterpolation(elem);
                reinit_all_data = false;
            }
            fe.interpolate(elem, side);
            const unsigned int n_qp = qrule_face->n_points();
            const size_t n_basis = phi_face.size();
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                const libMesh::Point& X = q_point_face[qp];
                get_x_and_FF(x, FF, x_data[qp], grad_x_data[qp]);
                const double J = std::abs(FF.det());
                FF_trans = FF.transpose();
                tensor_inverse_transpose(FF_inv_trans, FF, NDIM);
                n = (FF_inv_trans * normal_face[qp]).unit();
                const double dA_da = 1.0 / (J * (FF_inv_trans * normal_face[qp]) * n);

                double Phi = 0.0;
                for (unsigned int k = 0; k < num_PK1_stress_fcns; ++k)
                {
                    if (!d_PK1_stress_fcn_data[part][k].fcn) continue;

                    // Compute the value of the first Piola-Kirchhoff stress
                    // tensor at the quadrature point and add the corresponding
                    // traction force to the right-hand-side vector.
                    if (d_PK1_stress_fcn_data[part][k].fcn)
                    {
                        d_PK1_stress_fcn_data[part][k].fcn(PP, FF, x, X, elem, PK1_stress_fcn_data[k], data_time,
                                                           d_PK1_stress_fcn_data[part][k].ctx);
                        Phi += n * ((PP * FF_trans) * n) / J;
                    }
                }

                if (d_lag_surface_force_fcn_data[part].fcn)
                {
                    // Compute the value of the surface force at the
                    // quadrature point and add the corresponding force to
                    // the right-hand-side vector.
                    d_lag_surface_force_fcn_data[part].fcn(F_s, FF, x, X, elem, side, lag_surface_force_fcn_data,
                                                           data_time, d_lag_surface_force_fcn_data[part].ctx);
                    Phi -= n * F_s * dA_da;
                }

                if (d_lag_surface_pressure_fcn_data[part].fcn)
                {
                    // Compute the value of the pressure at the quadrature
                    // point and add the corresponding force to the
                    // right-hand-side vector.
                    d_lag_surface_pressure_fcn_data[part].fcn(P, FF, x, X, elem, side, lag_surface_pressure_fcn_data,
                                                              data_time, d_lag_surface_pressure_fcn_data[part].ctx);
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
    Phi_system.solution->close();
    Phi_system.solution->localize(Phi_vec);
    Phi_dof_map.enforce_constraints_exactly(Phi_system, &Phi_vec);
    return;
}

void IBFEMethod::computeInteriorForceDensity(PetscVector<double>& G_vec,
                                             PetscVector<double>& x_vec,
                                             PetscVector<double>* Phi_vec,
                                             const double data_time,
                                             const unsigned int part)
{
    if (d_constrained_part[part]) return;

    // Extract the mesh.
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems->get_mesh();
    const BoundaryInfo& boundary_info = *mesh.boundary_info;
    const unsigned int dim = mesh.mesh_dimension();

    // Setup extra data needed to compute stresses/forces.
    const size_t num_PK1_stress_fcns = d_PK1_stress_fcn_data[part].size();
    std::vector<std::vector<DenseVector<double> > > PK1_stress_fcn_data(num_PK1_stress_fcns);
    std::vector<DenseVector<double> > lag_body_force_fcn_data, lag_surface_force_fcn_data,
        lag_surface_pressure_fcn_data;

    // Setup global and elemental right-hand-side vectors.
    AutoPtr<NumericVector<double> > G_rhs_vec = G_vec.zero_clone();
    DenseVector<double> G_rhs_e[NDIM];

    // First handle the stress contributions.  These are handled separately because
    // each stress function may use a different quadrature rule.
    for (unsigned int k = 0; k < num_PK1_stress_fcns; ++k)
    {
        if (!d_PK1_stress_fcn_data[part][k].fcn) continue;

        // Extract the FE systems and DOF maps, and setup the FE object.
        System& G_system = equation_systems->get_system(FORCE_SYSTEM_NAME);
        const DofMap& G_dof_map = G_system.get_dof_map();
        FEType G_fe_type = G_dof_map.variable_type(0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(G_dof_map.variable_type(d) == G_fe_type);
        }
        std::vector<std::vector<unsigned int> > G_dof_indices(NDIM);
        System& x_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;

        FEFunctionInterpolation fe(dim);
        AutoPtr<QBase> qrule =
            QBase::build(d_PK1_stress_fcn_data[part][k].quad_type, dim, d_PK1_stress_fcn_data[part][k].quad_order);
        AutoPtr<QBase> qrule_face =
            QBase::build(d_PK1_stress_fcn_data[part][k].quad_type, dim - 1, d_PK1_stress_fcn_data[part][k].quad_order);
        fe.attachQuadratureRule(qrule.get());
        fe.attachQuadratureRuleFace(qrule_face.get());
        fe.evalNormalsFace();
        fe.evalQuadraturePoints();
        fe.evalQuadraturePointsFace();
        fe.evalQuadratureWeights();
        fe.evalQuadratureWeightsFace();
        fe.registerSystem(G_system, std::vector<int>(), vars); // compute dphi for the force system
        const size_t x_sys_idx = fe.registerInterpolatedSystem(x_system, vars, vars, &x_vec);
        fe.init();

        const std::vector<libMesh::Point>& q_point = fe.getQuadraturePoints();
        const std::vector<double>& JxW = fe.getQuadratureWeights();
        const std::vector<std::vector<VectorValue<double> > >& dphi = fe.getDphi(G_fe_type);

        const std::vector<libMesh::Point>& q_point_face = fe.getQuadraturePointsFace();
        const std::vector<double>& JxW_face = fe.getQuadratureWeightsFace();
        const std::vector<libMesh::Point>& normal_face = fe.getNormalsFace();
        const std::vector<std::vector<double> >& phi_face = fe.getPhiFace(G_fe_type);

        const std::vector<std::vector<double> >& x_data = fe.getVarInterpolation(x_sys_idx);
        const std::vector<std::vector<VectorValue<double> > >& grad_x_data = fe.getGradVarInterpolation(x_sys_idx);

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
                G_dof_map.dof_indices(elem, G_dof_indices[d], d);
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
                get_x_and_FF(x, FF, x_data[qp], grad_x_data[qp]);

                // Compute the value of the first Piola-Kirchhoff stress tensor
                // at the quadrature point and add the corresponding forces to
                // the right-hand-side vector.
                d_PK1_stress_fcn_data[part][k].fcn(PP, FF, x, X, elem, PK1_stress_fcn_data[k], data_time,
                                                   d_PK1_stress_fcn_data[part][k].ctx);
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
                    get_x_and_FF(x, FF, x_data[qp], grad_x_data[qp]);
                    tensor_inverse_transpose(FF_inv_trans, FF, NDIM);

                    F.zero();

                    // Compute the value of the first Piola-Kirchhoff stress
                    // tensor at the quadrature point and add the corresponding
                    // traction force to the right-hand-side vector.
                    if (d_PK1_stress_fcn_data[part][k].fcn)
                    {
                        d_PK1_stress_fcn_data[part][k].fcn(PP, FF, x, X, elem, PK1_stress_fcn_data[k], data_time,
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
    System& G_system = equation_systems->get_system(FORCE_SYSTEM_NAME);
    const DofMap& G_dof_map = G_system.get_dof_map();
    FEType G_fe_type = G_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(G_dof_map.variable_type(d) == G_fe_type);
    }
    std::vector<std::vector<unsigned int> > G_dof_indices(NDIM);
    System& x_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    System* Phi_system = Phi_vec ? &equation_systems->get_system(PHI_SYSTEM_NAME) : NULL;
    std::vector<int> vars(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
    std::vector<int> Phi_vars(1, 0);
    std::vector<int> no_vars;

    FEFunctionInterpolation fe(dim);
    AutoPtr<QBase> qrule = QBase::build(d_quad_type, dim, d_quad_order);
    AutoPtr<QBase> qrule_face = QBase::build(d_quad_type, dim - 1, d_quad_order);
    fe.attachQuadratureRule(qrule.get());
    fe.attachQuadratureRuleFace(qrule_face.get());
    fe.evalNormalsFace();
    fe.evalQuadraturePoints();
    fe.evalQuadraturePointsFace();
    fe.evalQuadratureWeights();
    fe.evalQuadratureWeightsFace();
    fe.registerSystem(G_system, vars, vars); // compute phi and dphi for the force system
    const size_t x_sys_idx = fe.registerInterpolatedSystem(x_system, vars, vars, &x_vec);
    const size_t Phi_sys_idx =
        Phi_vec ? fe.registerInterpolatedSystem(*Phi_system, Phi_vars, no_vars, Phi_vec) : SIZE_T_MAX;
    fe.init();

    const std::vector<libMesh::Point>& q_point = fe.getQuadraturePoints();
    const std::vector<double>& JxW = fe.getQuadratureWeights();
    const std::vector<std::vector<double> >& phi = fe.getPhi(G_fe_type);
    const std::vector<std::vector<VectorValue<double> > >& dphi = fe.getDphi(G_fe_type);

    const std::vector<libMesh::Point>& q_point_face = fe.getQuadraturePointsFace();
    const std::vector<double>& JxW_face = fe.getQuadratureWeightsFace();
    const std::vector<libMesh::Point>& normal_face = fe.getNormalsFace();
    const std::vector<std::vector<double> >& phi_face = fe.getPhiFace(G_fe_type);

    const std::vector<std::vector<double> >& x_data = fe.getVarInterpolation(x_sys_idx);
    const std::vector<std::vector<VectorValue<double> > >& grad_x_data = fe.getGradVarInterpolation(x_sys_idx);
    const std::vector<std::vector<double> >* const phi_data = Phi_vec ? &fe.getVarInterpolation(Phi_sys_idx) : NULL;

    // Loop over the elements to compute the right-hand side vector.
    TensorValue<double> PP, FF, FF_inv_trans;
    VectorValue<double> F, F_b, F_s, F_qp, n, x;
    double J, P, Phi;
    boost::multi_array<double, 2> X_node;
    boost::multi_array<double, 1> Phi_node;
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            G_dof_map.dof_indices(elem, G_dof_indices[d], d);
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
            get_x_and_FF(x, FF, x_data[qp], grad_x_data[qp]);
            J = std::abs(FF.det());
            tensor_inverse_transpose(FF_inv_trans, FF, NDIM);
            Phi = Phi_vec ? (*phi_data)[qp][0] : std::numeric_limits<double>::quiet_NaN();

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
                d_lag_body_force_fcn_data[part].fcn(F_b, FF, x, X, elem, lag_body_force_fcn_data, data_time,
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
                get_x_and_FF(x, FF, x_data[qp], grad_x_data[qp]);
                J = std::abs(FF.det());
                tensor_inverse_transpose(FF_inv_trans, FF, NDIM);

                F.zero();

                if (d_lag_surface_pressure_fcn_data[part].fcn)
                {
                    // Compute the value of the pressure at the quadrature
                    // point and add the corresponding force to the
                    // right-hand-side vector.
                    d_lag_surface_pressure_fcn_data[part].fcn(P, FF, x, X, elem, side, lag_surface_pressure_fcn_data,
                                                              data_time, d_lag_surface_pressure_fcn_data[part].ctx);
                    F -= P * J * FF_inv_trans * normal_face[qp];
                }

                if (d_lag_surface_force_fcn_data[part].fcn)
                {
                    // Compute the value of the surface force at the
                    // quadrature point and add the corresponding force to
                    // the right-hand-side vector.
                    d_lag_surface_force_fcn_data[part].fcn(F_s, FF, x, X, elem, side, lag_surface_force_fcn_data,
                                                           data_time, d_lag_surface_force_fcn_data[part].ctx);
                    F += F_s;
                }

                // Compute the outward unit normal vector in the current configuration.
                n = (FF_inv_trans * normal_face[qp]).unit();

                // Remote the normal component of the boundary force when needed.
                if (!integrate_normal_force)
                {
                    F -= (F * n) * n;
                }

                // Remote the tangential component of the boundary force when needed.
                if (!integrate_tangential_force)
                {
                    F -= (F - (F * n) * n);
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

    // Solve for G.
    d_fe_data_managers[part]->computeL2Projection(G_vec, *G_rhs_vec, FORCE_SYSTEM_NAME, d_use_consistent_mass_matrix);
    return;
} // computeInteriorForceDensity

void IBFEMethod::spreadTransmissionForceDensity(const int f_data_idx,
                                                PetscVector<double>& X_ghost_vec,
                                                RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                                                const double data_time,
                                                const unsigned int part)
{
    if (d_constrained_part[part] || (!d_split_normal_force && !d_split_tangential_force)) return;

    // Check to see if we need to integrate the surface forces.
    const bool integrate_normal_force =
        d_split_normal_force && !d_use_jump_conditions && !d_stress_normalization_part[part];
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
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems->get_mesh();
    const BoundaryInfo& boundary_info = *mesh.boundary_info;
    const unsigned int dim = mesh.mesh_dimension();
    AutoPtr<QBase> qrule_face;

    // Extract the FE systems and DOF maps, and setup the FE object.
    System& system = equation_systems->get_system(FORCE_SYSTEM_NAME);
    const DofMap& dof_map = system.get_dof_map();
    std::vector<std::vector<unsigned int> > dof_indices(NDIM);
    std::vector<std::vector<unsigned int> > side_dof_indices(NDIM);
    FEType fe_type = dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(dof_map.variable_type(d) == fe_type);
    }
    AutoPtr<FEBase> fe_face(FEBase::build(dim, fe_type));
    const std::vector<libMesh::Point>& q_point_face = fe_face->get_xyz();
    const std::vector<double>& JxW_face = fe_face->get_JxW();
    const std::vector<libMesh::Point>& normal_face = fe_face->get_normals();
    const std::vector<std::vector<double> >& phi_face = fe_face->get_phi();
    const std::vector<std::vector<VectorValue<double> > >& dphi_face = fe_face->get_dphi();

    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(X_dof_map.variable_type(d) == fe_type);
    }

    // Setup extra data needed to compute stresses/forces.
    const size_t num_PK1_stress_fcns = d_PK1_stress_fcn_data[part].size();
    std::vector<std::vector<DenseVector<double> > > PK1_stress_fcn_data(num_PK1_stress_fcns);
    std::vector<DenseVector<double> > lag_body_force_fcn_data, lag_surface_force_fcn_data,
        lag_surface_pressure_fcn_data;

    // Extract the underlying solution data.
    PetscVector<double>* X_petsc_vec = dynamic_cast<PetscVector<double>*>(&X_ghost_vec);
    Vec X_global_vec = X_petsc_vec->vec();
    Vec X_local_vec;
    VecGhostGetLocalForm(X_global_vec, &X_local_vec);
    double* X_local_soln;
    VecGetArray(X_local_vec, &X_local_soln);

    // Loop over the patches to spread the transmission elastic force density
    // onto the grid.
    const std::vector<std::vector<Elem*> >& active_patch_element_map =
        d_fe_data_managers[part]->getActivePatchElementMap();
    const int level_num = d_fe_data_managers[part]->getLevelNumber();
    TensorValue<double> PP, FF, FF_inv_trans;
    VectorValue<double> F, F_s, n;
    libMesh::Point X_qp;
    double P;
    boost::multi_array<double, 2> X_node;
    std::vector<double> T_bdry, X_bdry;
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
        X_bdry.clear();
        int qp_offset = 0;
        for (size_t e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems[e_idx];
            const bool touches_physical_bdry = has_physical_bdry(elem, boundary_info, dof_map);
            if (!touches_physical_bdry) continue;

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dof_map.dof_indices(elem, dof_indices[d], d);
            }
            get_values_for_interpolation(X_node, *X_petsc_vec, X_local_soln, dof_indices);

            // Loop over the element boundaries.
            for (unsigned short int side = 0; side < elem->n_sides(); ++side)
            {
                // Skip non-physical boundaries.
                if (!is_physical_bdry(elem, side, boundary_info, dof_map)) continue;

                // Skip Dirichlet boundaries.
                if (is_dirichlet_bdry(elem, side, boundary_info, dof_map)) continue;

                // Construct a side element.
                AutoPtr<Elem> side_elem = elem->build_side(side, /*proxy*/ false);
                const bool qrule_needs_reinit = d_fe_data_managers[part]->updateSpreadQuadratureRule(
                    qrule_face, d_spread_spec, side_elem.get(), X_node, patch_dx_min);
                if (qrule_needs_reinit)
                {
                    fe_face->attach_quadrature_rule(qrule_face.get());
                }
                fe_face->reinit(elem, side);
                const unsigned int n_qp = qrule_face->n_points();
                T_bdry.resize(T_bdry.size() + NDIM * n_qp);
                X_bdry.resize(X_bdry.size() + NDIM * n_qp);
                for (unsigned int qp = 0; qp < n_qp; ++qp, ++qp_offset)
                {
                    const libMesh::Point& s_qp = q_point_face[qp];
                    interpolate(X_qp, qp, X_node, phi_face);
                    jacobian(FF, qp, X_node, dphi_face);
                    const double J = std::abs(FF.det());
                    tensor_inverse_transpose(FF_inv_trans, FF, NDIM);
                    F.zero();

                    for (unsigned int k = 0; k < num_PK1_stress_fcns; ++k)
                    {
                        if (d_PK1_stress_fcn_data[part][k].fcn)
                        {
                            // Compute the value of the first Piola-Kirchhoff stress
                            // tensor at the quadrature point and compute the
                            // corresponding force.
                            d_PK1_stress_fcn_data[part][k].fcn(PP, FF, X_qp, s_qp, elem, PK1_stress_fcn_data[k],
                                                               data_time, d_PK1_stress_fcn_data[part][k].ctx);
                            F -= PP * normal_face[qp] * JxW_face[qp];
                        }
                    }

                    if (d_lag_surface_pressure_fcn_data[part].fcn)
                    {
                        // Compute the value of the pressure at the quadrature
                        // point and compute the corresponding force.
                        d_lag_surface_pressure_fcn_data[part].fcn(P, FF, X_qp, s_qp, elem, side,
                                                                  lag_surface_pressure_fcn_data, data_time,
                                                                  d_lag_surface_pressure_fcn_data[part].ctx);
                        F -= P * J * FF_inv_trans * normal_face[qp] * JxW_face[qp];
                    }

                    if (d_lag_surface_force_fcn_data[part].fcn)
                    {
                        // Compute the value of the surface force at the
                        // quadrature point and compute the corresponding force.
                        d_lag_surface_force_fcn_data[part].fcn(F_s, FF, X_qp, s_qp, elem, side,
                                                               lag_surface_force_fcn_data, data_time,
                                                               d_lag_surface_force_fcn_data[part].ctx);
                        F += F_s * JxW_face[qp];
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

                    const int idx = NDIM * qp_offset;
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        T_bdry[idx + i] = F(i);
                    }
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        X_bdry[idx + i] = X_qp(i);
                    }
                }
            }
        }

        if (qp_offset == 0) continue;

        // Spread the boundary forces to the grid.
        const std::string& spread_kernel_fcn = d_spread_spec.kernel_fcn;
        const hier::IntVector<NDIM>& ghost_width = d_fe_data_managers[part]->getGhostCellWidth();
        const Box<NDIM> spread_box = Box<NDIM>::grow(patch->getBox(), ghost_width);
        Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(f_data_idx);
        LEInteractor::spread(f_data, T_bdry, NDIM, X_bdry, NDIM, patch, spread_box, spread_kernel_fcn);
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

    VecRestoreArray(X_local_vec, &X_local_soln);
    VecGhostRestoreLocalForm(X_global_vec, &X_local_vec);
    return;
} // spreadTransmissionForceDensity

void IBFEMethod::imposeJumpConditions(const int f_data_idx,
                                      PetscVector<double>& F_ghost_vec,
                                      PetscVector<double>& X_ghost_vec,
                                      const double data_time,
                                      const unsigned int part)
{
    if (d_constrained_part[part] || (!d_split_normal_force && !d_split_tangential_force)) return;

    // Check to see if we need to integrate the normal surface force.
    const bool integrate_normal_force =
        d_split_normal_force && d_use_jump_conditions && !d_stress_normalization_part[part];
    if (!integrate_normal_force) return;

    // Extract the mesh.
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems->get_mesh();
    const BoundaryInfo& boundary_info = *mesh.boundary_info;
    const unsigned int dim = mesh.mesh_dimension();
    TBOX_ASSERT(dim == NDIM);

    // Extract the FE systems and DOF maps, and setup the FE object.
    System& system = equation_systems->get_system(FORCE_SYSTEM_NAME);
    const DofMap& dof_map = system.get_dof_map();
    std::vector<std::vector<unsigned int> > dof_indices(NDIM);
    std::vector<std::vector<unsigned int> > side_dof_indices(NDIM);
    FEType fe_type = dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(dof_map.variable_type(d) == fe_type);
    }
    AutoPtr<FEBase> fe_face(FEBase::build(dim, fe_type));
    const std::vector<libMesh::Point>& q_point_face = fe_face->get_xyz();
    const std::vector<libMesh::Point>& normal_face = fe_face->get_normals();
    const std::vector<std::vector<double> >& phi_face = fe_face->get_phi();
    const std::vector<std::vector<VectorValue<double> > >& dphi_face = fe_face->get_dphi();

    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(X_dof_map.variable_type(d) == fe_type);
    }

    // Setup extra data needed to compute stresses/forces.
    const size_t num_PK1_stress_fcns = d_PK1_stress_fcn_data[part].size();
    std::vector<std::vector<DenseVector<double> > > PK1_stress_fcn_data(num_PK1_stress_fcns);
    std::vector<DenseVector<double> > lag_body_force_fcn_data, lag_surface_pressure_fcn_data,
        lag_surface_force_fcn_data;

    // Extract the underlying solution data.
    PetscVector<double>* F_petsc_vec = dynamic_cast<PetscVector<double>*>(&F_ghost_vec);
    Vec F_global_vec = F_petsc_vec->vec();
    Vec F_local_vec;
    VecGhostGetLocalForm(F_global_vec, &F_local_vec);
    double* F_local_soln;
    VecGetArray(F_local_vec, &F_local_soln);

    PetscVector<double>* X_petsc_vec = dynamic_cast<PetscVector<double>*>(&X_ghost_vec);
    Vec X_global_vec = X_petsc_vec->vec();
    Vec X_local_vec;
    VecGhostGetLocalForm(X_global_vec, &X_local_vec);
    double* X_local_soln;
    VecGetArray(X_local_vec, &X_local_soln);

    // Loop over the patches to impose jump conditions on the Eulerian grid that
    // are determined from the interior and transmission elastic force
    // densities.
    const std::vector<std::vector<Elem*> >& active_patch_element_map =
        d_fe_data_managers[part]->getActivePatchElementMap();
    const int level_num = d_fe_data_managers[part]->getLevelNumber();
    TensorValue<double> PP, FF, FF_inv_trans;
    VectorValue<double> F, F_s, F_qp, n;
    libMesh::Point X_qp;
    double P;
    boost::multi_array<double, 2> F_node, X_node;
    std::vector<libMesh::Point> s_node_cache, X_node_cache;
    IBTK::Point X_min, X_max;
    std::vector<libMesh::Point> intersection_ref_coords;
    std::vector<SideIndex<NDIM> > intersection_indices;
    std::vector<std::pair<double, libMesh::Point> > intersections;
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
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
        const CellIndex<NDIM>& patch_upper = patch_box.upper();
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const x_lower = patch_geom->getXLower();
        const double* const x_upper = patch_geom->getXUpper();
        const double* const dx = patch_geom->getDx();

        SideData<NDIM, int> num_intersections(patch_box, 1, IntVector<NDIM>(0));
        num_intersections.fillAll(0);

        // Loop over the elements.
        for (size_t e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
        {
            Elem* const elem = patch_elems[e_idx];
            const bool touches_physical_bdry = has_physical_bdry(elem, boundary_info, dof_map);
            if (!touches_physical_bdry) continue;

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dof_map.dof_indices(elem, dof_indices[d], d);
            }

            // Loop over the element boundaries.
            for (unsigned short int side = 0; side < elem->n_sides(); ++side)
            {
                // Skip non-physical boundaries.
                if (!is_physical_bdry(elem, side, boundary_info, dof_map)) continue;

                // Skip Dirichlet boundaries.
                if (is_dirichlet_bdry(elem, side, boundary_info, dof_map)) continue;

                // Construct a side element.
                AutoPtr<Elem> side_elem = elem->build_side(side, /*proxy*/ false);
                const unsigned int n_node_side = side_elem->n_nodes();
                for (int d = 0; d < NDIM; ++d)
                {
                    dof_map.dof_indices(side_elem.get(), side_dof_indices[d], d);
                }

                // Cache the nodal and physical coordinates of the side element,
                // determine the bounding box of the current configuration of
                // the side element, and set the nodal coordinates to correspond
                // to the physical coordinates.
                s_node_cache.resize(n_node_side);
                X_node_cache.resize(n_node_side);
                X_min = IBTK::Point::Constant(+0.5 * std::numeric_limits<double>::max());
                X_max = IBTK::Point::Constant(-0.5 * std::numeric_limits<double>::max());
                for (unsigned int k = 0; k < n_node_side; ++k)
                {
                    s_node_cache[k] = side_elem->point(k);
                    libMesh::Point& X = X_node_cache[k];
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        X(d) = X_ghost_vec(side_dof_indices[d][k]);
                        X_min[d] = std::min(X_min[d], X(d));
                        X_max[d] = std::max(X_max[d], X(d));
                    }
                    side_elem->point(k) = X;
                }
                Box<NDIM> box(IndexUtilities::getCellIndex(&X_min[0], x_lower, x_upper, dx, patch_lower, patch_upper),
                              IndexUtilities::getCellIndex(&X_max[0], x_lower, x_upper, dx, patch_lower, patch_upper));
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
                            r(d) = (d == axis ? 0.0 : x_lower[d] +
                                                          dx[d] * (static_cast<double>(i_c(d) - patch_lower[d]) + 0.5));
                        }
#if (NDIM == 2)
                        intersect_line_with_edge(intersections, static_cast<Edge*>(side_elem.get()), r, q);
#endif
#if (NDIM == 3)
                        intersect_line_with_face(intersections, static_cast<Face*>(side_elem.get()), r, q);
#endif
                        for (unsigned int k = 0; k < intersections.size(); ++k)
                        {
                            libMesh::Point X = r + intersections[k].first * q;
                            SideIndex<NDIM> i_s(i_c, axis, 0);
                            i_s(axis) = std::floor((X(axis) - x_lower[axis]) / dx[axis] + 0.5) + patch_lower[axis];
                            intersection_ref_coords.push_back(intersections[k].second);
                            intersection_indices.push_back(i_s);
                            num_intersections(i_s) += 1;
                        }
                    }
                }

                // Restore the element coordinates.
                for (unsigned int k = 0; k < n_node_side; ++k)
                {
                    side_elem->point(k) = s_node_cache[k];
                }

                // If there are no intersection points, then continue to the
                // next side.
                if (intersection_ref_coords.empty()) continue;

                // Evaluate the jump conditions and apply them to the Eulerian
                // grid.
                const bool impose_dp_dn_jumps = !d_use_IB_spread_operator;
                static const double TOL = sqrt(std::numeric_limits<double>::epsilon());
                fe_face->reinit(elem, side, TOL, &intersection_ref_coords);
                if (impose_dp_dn_jumps) get_values_for_interpolation(F_node, *F_petsc_vec, F_local_soln, dof_indices);
                get_values_for_interpolation(X_node, *X_petsc_vec, X_local_soln, dof_indices);
                for (unsigned int qp = 0; qp < intersection_ref_coords.size(); ++qp)
                {
                    const SideIndex<NDIM>& i_s = intersection_indices[qp];
                    const unsigned int axis = i_s.getAxis();
                    interpolate(X_qp, qp, X_node, phi_face);
#if !defined(NDEBUG)
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (d == axis)
                        {
                            const double X_lower_bound = x_lower[d] +
                                                         (static_cast<double>(i_s(d) - patch_lower[d]) - 0.5) * dx[d] -
                                                         sqrt(std::numeric_limits<double>::epsilon());
                            const double X_upper_bound = x_lower[d] +
                                                         (static_cast<double>(i_s(d) - patch_lower[d]) + 0.5) * dx[d] +
                                                         sqrt(std::numeric_limits<double>::epsilon());
                            TBOX_ASSERT(X_lower_bound <= X_qp(d) && X_qp(d) <= X_upper_bound);
                        }
                        else
                        {
                            const double X_intersection =
                                x_lower[d] + (static_cast<double>(i_s(d) - patch_lower[d]) + 0.5) * dx[d];
                            const double X_interp = X_qp(d);
                            const double rel_diff =
                                std::abs(X_intersection - X_interp) /
                                std::max(1.0, std::max(std::abs(X_intersection), std::abs(X_interp)));
                            TBOX_ASSERT(rel_diff <= sqrt(std::numeric_limits<double>::epsilon()));
                        }
                    }
#endif
                    const libMesh::Point& s_qp = q_point_face[qp];
                    interpolate(X_qp, qp, X_node, phi_face);
                    jacobian(FF, qp, X_node, dphi_face);
                    const double J = std::abs(FF.det());
                    tensor_inverse_transpose(FF_inv_trans, FF, NDIM);
                    F.zero();

                    for (unsigned int k = 0; k < num_PK1_stress_fcns; ++k)
                    {
                        if (d_PK1_stress_fcn_data[part][k].fcn)
                        {
                            // Compute the value of the first Piola-Kirchhoff
                            // stress tensor at the quadrature point and compute
                            // the corresponding force.
                            d_PK1_stress_fcn_data[part][k].fcn(PP, FF, X_qp, s_qp, elem, PK1_stress_fcn_data[k],
                                                               data_time, d_PK1_stress_fcn_data[part][k].ctx);
                            F -= PP * normal_face[qp];
                        }
                    }

                    if (d_lag_surface_pressure_fcn_data[part].fcn)
                    {
                        // Compute the value of the pressure at the quadrature
                        // point and compute the corresponding force.
                        d_lag_surface_pressure_fcn_data[part].fcn(P, FF, X_qp, s_qp, elem, side,
                                                                  lag_surface_pressure_fcn_data, data_time,
                                                                  d_lag_surface_pressure_fcn_data[part].ctx);
                        F -= P * J * FF_inv_trans * normal_face[qp];
                    }

                    if (d_lag_surface_force_fcn_data[part].fcn)
                    {
                        // Compute the value of the surface force at the
                        // quadrature point and compute the corresponding force.
                        d_lag_surface_force_fcn_data[part].fcn(F_s, FF, X_qp, s_qp, elem, side,
                                                               lag_surface_force_fcn_data, data_time,
                                                               d_lag_surface_force_fcn_data[part].ctx);
                        F += F_s;
                    }

                    // Use Nanson's formula (n da = J FF^{-T} N dA) to convert
                    // force per unit area in the reference configuration into
                    // force per unit area in the current configuration.  This
                    // value determines the discontinuity in the pressure at the
                    // fluid-structure interface.
                    n = (FF_inv_trans * normal_face[qp]).unit();
                    const double dA_da = 1.0 / (J * (FF_inv_trans * normal_face[qp]) * n);
                    F *= dA_da;

                    // Determine the value of the interior force density at the
                    // boundary, and convert it to force per unit volume in the
                    // current configuration.  This value determines the
                    // discontinuity in the normal derivative of the pressure at
                    // the fluid-structure interface.
                    if (!impose_dp_dn_jumps)
                    {
                        F_qp.zero();
                    }
                    else
                    {
                        interpolate(F_qp, qp, F_node, phi_face);
                        F_qp /= J;
                    }

                    // Impose the jump conditions.
                    const double X = X_qp(axis);
                    const double x_cell_bdry =
                        x_lower[axis] + static_cast<double>(i_s(axis) - patch_lower[axis]) * dx[axis];
                    const double h = x_cell_bdry + (X > x_cell_bdry ? +0.5 : -0.5) * dx[axis] - X;
                    const double C_p = F * n - h * F_qp(axis);
                    (*f_data)(i_s) += (n(axis) > 0.0 ? +1.0 : -1.0) * (C_p / dx[axis]);
                }
            }
        }
    }

    VecRestoreArray(F_local_vec, &F_local_soln);
    VecGhostRestoreLocalForm(F_global_vec, &F_local_vec);

    VecRestoreArray(X_local_vec, &X_local_soln);
    VecGhostRestoreLocalForm(X_global_vec, &X_local_vec);
    return;
} // imposeJumpConditions

void IBFEMethod::initializeCoordinates(const unsigned int part)
{
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    MeshBase& mesh = equation_systems->get_mesh();
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    const bool identity_mapping = !d_coordinate_mapping_fcn_data[part].fcn;
    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
    {
        Node* n = *it;
        if (n->n_vars(X_sys_num))
        {
            TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM);
            const libMesh::Point& s = *n;
            libMesh::Point X = s;
            if (!identity_mapping)
            {
                d_coordinate_mapping_fcn_data[part].fcn(X, s, d_coordinate_mapping_fcn_data[part].ctx);
            }
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const int dof_index = n->dof_number(X_sys_num, d, 0);
                X_coords.set(dof_index, X(d));
            }
        }
    }
    X_coords.close();
    X_system.get_dof_map().enforce_constraints_exactly(X_system, &X_coords);
    return;
} // initializeCoordinates

void IBFEMethod::updateCoordinateMapping(const unsigned int part)
{
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    MeshBase& mesh = equation_systems->get_mesh();
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    System& dX_system = equation_systems->get_system(COORD_MAPPING_SYSTEM_NAME);
    const unsigned int dX_sys_num = dX_system.number();
    NumericVector<double>& dX_coords = *dX_system.solution;
    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
    {
        Node* n = *it;
        if (n->n_vars(X_sys_num))
        {
            TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM);
            TBOX_ASSERT(n->n_vars(dX_sys_num) == NDIM);
            const libMesh::Point& s = *n;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const int X_dof_index = n->dof_number(X_sys_num, d, 0);
                const int dX_dof_index = n->dof_number(dX_sys_num, d, 0);
                dX_coords.set(dX_dof_index, X_coords(X_dof_index) - s(d));
            }
        }
    }
    dX_coords.close();
    return;
} // updateCoordinateMapping

/////////////////////////////// PRIVATE //////////////////////////////////////

void IBFEMethod::commonConstructor(const std::string& object_name,
                                   Pointer<Database> input_db,
                                   const std::vector<libMesh::Mesh*>& meshes,
                                   int max_level_number,
                                   bool register_for_restart)
{
    // Set the object name and register it with the restart manager.
    d_object_name = object_name;
    d_registered_for_restart = false;
    if (register_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
        d_registered_for_restart = true;
    }

    // Set some default values.
    const bool use_adaptive_quadrature = true;
    const int point_density = 2.0;
    const bool interp_use_consistent_mass_matrix = true;
    d_use_IB_interp_operator = true;
    d_interp_spec = FEDataManager::InterpSpec("IB_4", QGAUSS, INVALID_ORDER, use_adaptive_quadrature, point_density,
                                              interp_use_consistent_mass_matrix);
    d_use_IB_spread_operator = true;
    d_spread_spec = FEDataManager::SpreadSpec("IB_4", QGAUSS, INVALID_ORDER, use_adaptive_quadrature, point_density);
    d_ghosts = 0;
    d_split_normal_force = false;
    d_split_tangential_force = false;
    d_use_jump_conditions = false;
    d_fe_family = LAGRANGE;
    d_fe_order = INVALID_ORDER;
    d_quad_type = QGAUSS;
    d_quad_order = INVALID_ORDER;
    d_use_consistent_mass_matrix = true;
    d_do_log = false;

    // Indicate that all of the parts do NOT use stress normalization by default
    // and set some default values.
    d_epsilon = 0.0;
    d_has_stress_normalization_parts = false;
    d_stress_normalization_part.resize(d_num_parts, false);

    // Indicate that all of the parts are unconstrained by default and set some
    // default values.
    d_has_constrained_parts = false;
    d_constrained_part.resize(d_num_parts, false);
    d_constrained_velocity_fcn_data.resize(d_num_parts);
    d_constraint_omega = 2.0;

    // Initialize function data to NULL.
    d_coordinate_mapping_fcn_data.resize(d_num_parts);
    d_PK1_stress_fcn_data.resize(d_num_parts);
    d_lag_body_force_fcn_data.resize(d_num_parts);
    d_lag_surface_pressure_fcn_data.resize(d_num_parts);
    d_lag_surface_force_fcn_data.resize(d_num_parts);
    d_fcn_systems.resize(d_num_parts);
    d_body_fcn_systems.resize(d_num_parts);
    d_surface_fcn_systems.resize(d_num_parts);

    // Determine whether we should use first-order or second-order shape
    // functions for each part of the structure.
    bool mesh_has_first_order_elems = false;
    bool mesh_has_second_order_elems = false;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        const MeshBase& mesh = *meshes[part];
        MeshBase::const_element_iterator el_it = mesh.elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.elements_end();
        for (; el_it != el_end; ++el_it)
        {
            const Elem* const elem = *el_it;
            mesh_has_first_order_elems = mesh_has_first_order_elems || elem->default_order() == FIRST;
            mesh_has_second_order_elems = mesh_has_second_order_elems || elem->default_order() == SECOND;
        }
    }
    mesh_has_first_order_elems = SAMRAI_MPI::maxReduction(mesh_has_first_order_elems);
    mesh_has_second_order_elems = SAMRAI_MPI::maxReduction(mesh_has_second_order_elems);
    if ((mesh_has_first_order_elems && mesh_has_second_order_elems) ||
        (!mesh_has_first_order_elems && !mesh_has_second_order_elems))
    {
        TBOX_ERROR(d_object_name << "::IBFEMethod():\n"
                                 << "  all parts of FE mesh must contain only FIRST order elements "
                                    "or only SECOND order elements"
                                 << std::endl);
    }
    if (mesh_has_first_order_elems)
    {
        d_fe_order = FIRST;
        d_quad_order = THIRD;
    }
    if (mesh_has_second_order_elems)
    {
        d_fe_order = SECOND;
        d_quad_order = FIFTH;
    }

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    // Report configuration.
    pout << "\n";
    pout << d_object_name << ": using " << Utility::enum_to_string<Order>(d_fe_order) << " order "
         << Utility::enum_to_string<FEFamily>(d_fe_family) << " finite elements.\n";
    pout << "\n";

    // Check the choices for the kernel function.
    if (d_interp_spec.kernel_fcn != d_spread_spec.kernel_fcn)
    {
        pout << "WARNING: different kernel functions are being used for velocity "
                "interpolation and "
                "force spreading.\n"
             << "         recommended usage is to employ the same kernel functions for both "
                "interpolation and spreading.\n";
    }

    // Create the FE data managers that manage mappings between the FE mesh
    // parts and the Cartesian grid.
    d_meshes = meshes;
    d_equation_systems.resize(d_num_parts, NULL);
    d_fe_data_managers.resize(d_num_parts, NULL);
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        // Create FE data managers.
        std::ostringstream manager_stream;
        manager_stream << "IBFEMethod FEDataManager::" << part;
        const std::string& manager_name = manager_stream.str();
        d_fe_data_managers[part] = FEDataManager::getManager(manager_name, d_interp_spec, d_spread_spec);
        d_ghosts = IntVector<NDIM>::max(d_ghosts, d_fe_data_managers[part]->getGhostCellWidth());

        // Create FE equation systems objects and corresponding variables.
        d_equation_systems[part] = new EquationSystems(*d_meshes[part]);
        EquationSystems* equation_systems = d_equation_systems[part];
        d_fe_data_managers[part]->setEquationSystems(equation_systems, max_level_number - 1);

        d_fe_data_managers[part]->COORDINATES_SYSTEM_NAME = COORDS_SYSTEM_NAME;
        System& X_system = equation_systems->add_system<System>(COORDS_SYSTEM_NAME);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "X_" << d;
            X_system.add_variable(os.str(), d_fe_order, d_fe_family);
        }

        System& dX_system = equation_systems->add_system<System>(COORD_MAPPING_SYSTEM_NAME);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "dX_" << d;
            dX_system.add_variable(os.str(), d_fe_order, d_fe_family);
        }

        System& U_system = equation_systems->add_system<System>(VELOCITY_SYSTEM_NAME);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "U_" << d;
            U_system.add_variable(os.str(), d_fe_order, d_fe_family);
        }

        System& F_system = equation_systems->add_system<System>(FORCE_SYSTEM_NAME);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "F_" << d;
            F_system.add_variable(os.str(), d_fe_order, d_fe_family);
        }
    }

    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time = std::numeric_limits<double>::quiet_NaN();
    d_half_time = std::numeric_limits<double>::quiet_NaN();

    // Keep track of the initialization state.
    d_fe_data_initialized = false;
    d_is_initialized = false;
    return;
} // commonConstructor

void IBFEMethod::getFromInput(Pointer<Database> db, bool /*is_from_restart*/)
{
    // Interpolation settings.
    if (db->isBool("use_IB_interp_operator"))
        d_use_IB_interp_operator = db->getBool("use_IB_interp_operator");
    else if (db->isBool("use_IB_interaction_operators"))
        d_use_IB_interp_operator = db->getBool("use_IB_interaction_operators");

    if (db->isString("interp_delta_fcn"))
        d_interp_spec.kernel_fcn = db->getString("interp_delta_fcn");
    else if (db->isString("IB_delta_fcn"))
        d_interp_spec.kernel_fcn = db->getString("IB_delta_fcn");
    else if (db->isString("interp_kernel_fcn"))
        d_interp_spec.kernel_fcn = db->getString("interp_kernel_fcn");
    else if (db->isString("IB_kernel_fcn"))
        d_interp_spec.kernel_fcn = db->getString("IB_kernel_fcn");

    if (db->isString("interp_quad_type"))
        d_interp_spec.quad_type = Utility::string_to_enum<QuadratureType>(db->getString("interp_quad_type"));
    else if (db->isString("IB_quad_type"))
        d_interp_spec.quad_type = Utility::string_to_enum<QuadratureType>(db->getString("IB_quad_type"));

    if (db->isString("interp_quad_order"))
        d_interp_spec.quad_order = Utility::string_to_enum<Order>(db->getString("interp_quad_order"));
    else if (db->isString("IB_quad_order"))
        d_interp_spec.quad_order = Utility::string_to_enum<Order>(db->getString("IB_quad_order"));

    if (db->isBool("interp_use_adaptive_quadrature"))
        d_interp_spec.use_adaptive_quadrature = db->getBool("interp_use_adaptive_quadrature");
    else if (db->isBool("IB_use_adaptive_quadrature"))
        d_interp_spec.use_adaptive_quadrature = db->getBool("IB_use_adaptive_quadrature");

    if (db->isDouble("interp_point_density"))
        d_interp_spec.point_density = db->getDouble("interp_point_density");
    else if (db->isDouble("IB_point_density"))
        d_interp_spec.point_density = db->getDouble("IB_point_density");

    if (db->isBool("interp_use_consistent_mass_matrix"))
        d_interp_spec.use_consistent_mass_matrix = db->getBool("interp_use_consistent_mass_matrix");
    else if (db->isBool("IB_use_consistent_mass_matrix"))
        d_interp_spec.use_consistent_mass_matrix = db->getBool("IB_use_consistent_mass_matrix");

    // Spreading settings.
    if (db->isBool("use_IB_spread_operator"))
        d_use_IB_spread_operator = db->getBool("use_IB_spread_operator");
    else if (db->isBool("use_IB_interaction_operators"))
        d_use_IB_spread_operator = db->getBool("use_IB_interaction_operators");

    if (db->isString("spread_delta_fcn"))
        d_spread_spec.kernel_fcn = db->getString("spread_delta_fcn");
    else if (db->isString("IB_delta_fcn"))
        d_spread_spec.kernel_fcn = db->getString("IB_delta_fcn");
    else if (db->isString("spread_kernel_fcn"))
        d_spread_spec.kernel_fcn = db->getString("spread_kernel_fcn");
    else if (db->isString("IB_kernel_fcn"))
        d_spread_spec.kernel_fcn = db->getString("IB_kernel_fcn");

    if (db->isString("spread_quad_type"))
        d_spread_spec.quad_type = Utility::string_to_enum<QuadratureType>(db->getString("spread_quad_type"));
    else if (db->isString("IB_quad_type"))
        d_spread_spec.quad_type = Utility::string_to_enum<QuadratureType>(db->getString("IB_quad_type"));

    if (db->isString("spread_quad_order"))
        d_spread_spec.quad_order = Utility::string_to_enum<Order>(db->getString("spread_quad_order"));
    else if (db->isString("IB_quad_order"))
        d_spread_spec.quad_order = Utility::string_to_enum<Order>(db->getString("IB_quad_order"));

    if (db->isBool("spread_use_adaptive_quadrature"))
        d_spread_spec.use_adaptive_quadrature = db->getBool("spread_use_adaptive_quadrature");
    else if (db->isBool("IB_use_adaptive_quadrature"))
        d_spread_spec.use_adaptive_quadrature = db->getBool("IB_use_adaptive_quadrature");

    if (db->isDouble("spread_point_density"))
        d_spread_spec.point_density = db->getDouble("spread_point_density");
    else if (db->isDouble("IB_point_density"))
        d_spread_spec.point_density = db->getDouble("IB_point_density");

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
    if (db->isString("quad_type")) d_quad_type = Utility::string_to_enum<QuadratureType>(db->getString("quad_type"));
    if (db->isString("quad_order")) d_quad_order = Utility::string_to_enum<Order>(db->getString("quad_order"));
    if (db->isBool("use_consistent_mass_matrix"))
        d_use_consistent_mass_matrix = db->getBool("use_consistent_mass_matrix");

    // Other settings.
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

    if (db->isDouble("constraint_omega")) d_constraint_omega = db->getDouble("constraint_omega");
    return;
} // getFromInput

void IBFEMethod::getFromRestart()
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
    d_fe_family = Utility::string_to_enum<FEFamily>(db->getString("d_fe_family"));
    d_fe_order = Utility::string_to_enum<Order>(db->getString("d_fe_order"));
    d_quad_type = Utility::string_to_enum<QuadratureType>(db->getString("d_quad_type"));
    d_quad_order = Utility::string_to_enum<Order>(db->getString("d_quad_order"));
    d_use_consistent_mass_matrix = db->getBool("d_use_consistent_mass_matrix");
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
