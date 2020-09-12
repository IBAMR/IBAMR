// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
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

#include "ibamr/FEMechanicsBase.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep

#include "ibtk/FEDataInterpolation.h"
#include "ibtk/FEDataManager.h"
#include "ibtk/FEProjector.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/LibMeshSystemVectors.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/libmesh_utilities.h"

#include "tbox/RestartManager.h"

#include "libmesh/boundary_info.h"
#include "libmesh/compare_types.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/edge.h"
#include "libmesh/elem.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_parallel_type.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/enum_xdr_mode.h"
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"
#include "libmesh/face.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_type.h"
#include "libmesh/fem_context.h"
#include "libmesh/id_types.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh_config.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parameters.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/point.h"
#include "libmesh/quadrature.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/system.h"
#include "libmesh/tensor_value.h"
#include "libmesh/type_tensor.h"
#include "libmesh/type_vector.h"
#include "libmesh/variant_filter_iterator.h"
#include "libmesh/vector_value.h"

#include <iterator>
#include <utility>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of FEMechanicsBase restart file data.
const int FE_MECHANICS_BASE_VERSION = 1;

inline boundary_id_type
get_dirichlet_bdry_ids(const std::vector<boundary_id_type>& bdry_ids)
{
    boundary_id_type dirichlet_bdry_ids = 0;
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
    std::vector<boundary_id_type> bdry_ids;
    boundary_info.boundary_ids(elem, side, bdry_ids);
    bool at_physical_bdry = !elem->neighbor_ptr(side);
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
    std::vector<boundary_id_type> bdry_ids;
    boundary_info.boundary_ids(elem, side, bdry_ids);
    return get_dirichlet_bdry_ids(bdry_ids) != 0;
}

inline void
get_FF(libMesh::TensorValue<double>& FF,
       const std::vector<VectorValue<double> >& grad_x_data,
       const unsigned int dim = NDIM)
{
    FF.zero();
    for (unsigned int i = 0; i < dim; ++i)
    {
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
}
} // namespace

const std::string FEMechanicsBase::COORDS_SYSTEM_NAME = "IB coordinates system";
const std::string FEMechanicsBase::COORD_MAPPING_SYSTEM_NAME = "IB coordinate mapping system";
const std::string FEMechanicsBase::FORCE_SYSTEM_NAME = "IB force system";
const std::string FEMechanicsBase::PRESSURE_SYSTEM_NAME = "IB pressure system";
const std::string FEMechanicsBase::VELOCITY_SYSTEM_NAME = "IB velocity system";

/////////////////////////////// PUBLIC ///////////////////////////////////////

FEMechanicsBase::FEMechanicsBase(const std::string& object_name,
                                 const Pointer<Database>& input_db,
                                 MeshBase* mesh,
                                 bool register_for_restart,
                                 const std::string& restart_read_dirname,
                                 unsigned int restart_restore_number)
{
    commonConstructor(object_name,
                      input_db,
                      std::vector<MeshBase*>(1, mesh),
                      register_for_restart,
                      restart_read_dirname,
                      restart_restore_number);
}

FEMechanicsBase::FEMechanicsBase(const std::string& object_name,
                                 const Pointer<Database>& input_db,
                                 const std::vector<MeshBase*>& meshes,
                                 bool register_for_restart,
                                 const std::string& restart_read_dirname,
                                 unsigned int restart_restore_number)
{
    commonConstructor(
        object_name, input_db, meshes, register_for_restart, restart_read_dirname, restart_restore_number);
}

FEMechanicsBase::~FEMechanicsBase()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }
}

std::shared_ptr<FEData>
FEMechanicsBase::getFEData(const unsigned int part) const
{
    TBOX_ASSERT(d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_meshes.size());
    return d_fe_data[part];
}

void
FEMechanicsBase::registerInitialCoordinateMappingFunction(const CoordinateMappingFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_meshes.size());
    d_coordinate_mapping_fcn_data[part] = data;
}

FEMechanicsBase::CoordinateMappingFcnData
FEMechanicsBase::getInitialCoordinateMappingFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_meshes.size());
    return d_coordinate_mapping_fcn_data[part];
}

void
FEMechanicsBase::registerInitialVelocityFunction(const InitialVelocityFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_meshes.size());
    d_initial_velocity_fcn_data[part] = data;
}

FEMechanicsBase::InitialVelocityFcnData
FEMechanicsBase::getInitialVelocityFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_meshes.size());
    return d_initial_velocity_fcn_data[part];
}

void
FEMechanicsBase::registerPK1StressFunction(const PK1StressFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_meshes.size());
    d_PK1_stress_fcn_data[part].push_back(data);
    if (data.quad_type == INVALID_Q_RULE)
    {
        d_PK1_stress_fcn_data[part].back().quad_type = d_default_quad_type_stress[part];
    }
    if (data.quad_order == INVALID_ORDER)
    {
        d_PK1_stress_fcn_data[part].back().quad_order = d_default_quad_order_stress[part];
    }
}

std::vector<FEMechanicsBase::PK1StressFcnData>
FEMechanicsBase::getPK1StressFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_meshes.size());
    return d_PK1_stress_fcn_data[part];
}

void
FEMechanicsBase::registerLagBodyForceFunction(const LagBodyForceFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_meshes.size());
    d_lag_body_force_fcn_data[part] = data;
}

FEMechanicsBase::LagBodyForceFcnData
FEMechanicsBase::getLagBodyForceFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_meshes.size());
    return d_lag_body_force_fcn_data[part];
}

void
FEMechanicsBase::registerLagSurfacePressureFunction(const LagSurfacePressureFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_meshes.size());
    d_lag_surface_pressure_fcn_data[part] = data;
}

FEMechanicsBase::LagSurfacePressureFcnData
FEMechanicsBase::getLagSurfacePressureFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_meshes.size());
    return d_lag_surface_pressure_fcn_data[part];
}

void
FEMechanicsBase::registerLagSurfaceForceFunction(const LagSurfaceForceFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_meshes.size());
    d_lag_surface_force_fcn_data[part] = data;
}

FEMechanicsBase::LagSurfaceForceFcnData
FEMechanicsBase::getLagSurfaceForceFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_meshes.size());
    return d_lag_surface_force_fcn_data[part];
}

void
FEMechanicsBase::registerStaticPressurePart(PressureProjectionType projection_type,
                                            VolumetricEnergyDerivativeFcn U_prime_fcn,
                                            unsigned int part)
{
    TBOX_ASSERT(d_fe_equation_systems_initialized);
    TBOX_ASSERT(part < d_meshes.size());
    if (d_static_pressure_part[part]) return;
    d_has_static_pressure_parts = true;
    d_static_pressure_part[part] = true;
    auto& P_system = d_equation_systems[part]->add_system<ExplicitSystem>(PRESSURE_SYSTEM_NAME);
    // This system has a single variable so we don't need to also specify diagonal coupling
    P_system.add_variable("P", d_fe_order_pressure[part], d_fe_family_pressure[part]);
    // Setup cached system vectors at restart.
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    IBTK::setup_system_vectors(d_equation_systems[part].get(),
                               { PRESSURE_SYSTEM_NAME },
                               { "current", "half", "new", "tmp", "RHS Vector" },
                               from_restart);
    // Keep track of method parameters.
    d_static_pressure_proj_type[part] = projection_type;
    d_static_pressure_vol_energy_deriv_fcn[part] = std::move(U_prime_fcn);
}

void
FEMechanicsBase::preprocessIntegrateData(double current_time, double new_time, int /*num_cycles*/)
{
    // Keep track of the current time step interval.
    d_current_time = current_time;
    d_new_time = new_time;
    d_half_time = current_time + 0.5 * (new_time - current_time);
}

void
FEMechanicsBase::postprocessIntegrateData(double /*current_time*/, double /*new_time*/, int /*num_cycles*/)
{
    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time = std::numeric_limits<double>::quiet_NaN();
    d_half_time = std::numeric_limits<double>::quiet_NaN();

    // Update the coordinate mapping dX = X - s.
    for (unsigned part = 0; part < d_meshes.size(); ++part)
    {
        updateCoordinateMapping(part);
    }
}

void
FEMechanicsBase::initializeFEEquationSystems()
{
    if (d_fe_equation_systems_initialized) return;

    // Set up the coupling matrix that will be used by each system.
    d_diagonal_system_coupling.resize(NDIM);
    for (unsigned int i = 0; i < NDIM; ++i)
        for (unsigned int j = 0; j < NDIM; ++j) d_diagonal_system_coupling(i, j) = i == j ? 1 : 0;

    doInitializeFEEquationSystems();
    d_fe_equation_systems_initialized = true;
}

void
FEMechanicsBase::initializeFEData()
{
    if (d_fe_data_initialized) return;

    initializeFEEquationSystems();
    doInitializeFEData(RestartManager::getManager()->isFromRestart());
    d_fe_data_initialized = true;
}

void
FEMechanicsBase::reinitializeFEData()
{
    TBOX_ASSERT(d_fe_data_initialized);
    doInitializeFEData(true);
}

void
FEMechanicsBase::putToDatabase(Pointer<Database> db)
{
    db->putInteger("FE_MECHANICS_BASE_VERSION", FE_MECHANICS_BASE_VERSION);
    db->putBool("d_use_consistent_mass_matrix", d_use_consistent_mass_matrix);
    db->putString("d_libmesh_partitioner_type", enum_to_string<LibmeshPartitionerType>(d_libmesh_partitioner_type));
}

void
FEMechanicsBase::writeFEDataToRestartFile(const std::string& restart_dump_dirname, unsigned int time_step_number)
{
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        const std::string& file_name =
            libmesh_restart_file_name(restart_dump_dirname, time_step_number, part, d_libmesh_restart_file_extension);
        const XdrMODE xdr_mode = (d_libmesh_restart_file_extension == "xdr" ? ENCODE : WRITE);
        const int write_mode = EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA;
        d_equation_systems[part]->write(file_name,
                                        xdr_mode,
                                        write_mode,
                                        /*partition_agnostic*/ true);
    }
}

/////////////////////////////// PROTECTED ////////////////////////////////////

void
FEMechanicsBase::computeStaticPressure(PetscVector<double>& P_vec,
                                       PetscVector<double>& X_vec,
                                       const double /*data_time*/,
                                       const unsigned int part)
{
    if (!d_static_pressure_part[part]) return;
    const PressureProjectionType& proj_type = d_static_pressure_proj_type[part];
    const VolumetricEnergyDerivativeFcn& U_prime_fcn = d_static_pressure_vol_energy_deriv_fcn[part];

    // Extract the mesh.
    EquationSystems& equation_systems = *d_equation_systems[part];
    const MeshBase& mesh = equation_systems.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    // Setup extra data needed to compute stresses/forces.

    // Extract the FE systems and DOF maps, and setup the FE objects.
    auto& P_system = equation_systems.get_system<ExplicitSystem>(PRESSURE_SYSTEM_NAME);
    const DofMap& P_dof_map = P_system.get_dof_map();
    FEDataManager::SystemDofMapCache& P_dof_map_cache = *d_fe_data[part]->getDofMapCache(PRESSURE_SYSTEM_NAME);
    FEType P_fe_type = P_dof_map.variable_type(0);
    std::vector<int> P_vars = { 0 };
    std::vector<int> no_vars = {};
    auto& X_system = equation_systems.get_system<ExplicitSystem>(COORDS_SYSTEM_NAME);
    std::vector<int> X_vars(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) X_vars[d] = d;

    FEDataInterpolation fe(dim, d_fe_data[part]);
    std::unique_ptr<QBase> qrule =
        QBase::build(d_default_quad_type_pressure[part], dim, d_default_quad_order_pressure[part]);
    fe.attachQuadratureRule(qrule.get());
    fe.evalQuadraturePoints();
    fe.evalQuadratureWeights();
    fe.registerSystem(P_system, P_vars, no_vars);
    const size_t X_sys_idx = fe.registerInterpolatedSystem(X_system, no_vars, X_vars, &X_vec);
    fe.init();

    const std::vector<libMesh::Point>& q_point = fe.getQuadraturePoints();
    const std::vector<double>& JxW = fe.getQuadratureWeights();
    const std::vector<std::vector<double> >& phi = fe.getPhi(P_fe_type);

    const std::vector<std::vector<std::vector<double> > >& fe_interp_var_data = fe.getVarInterpolation();
    const std::vector<std::vector<std::vector<VectorValue<double> > > >& fe_interp_grad_var_data =
        fe.getGradVarInterpolation();

    // Setup global and elemental right-hand-side vectors.
    auto* P_rhs_vec = static_cast<PetscVector<double>*>(P_system.rhs);
    P_rhs_vec->zero();
    DenseVector<double> P_rhs_e;

    TensorValue<double> FF;
    std::vector<libMesh::dof_id_type> dof_id_scratch;
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;
        const auto& P_dof_indices = P_dof_map_cache.dof_indices(elem);
        P_rhs_e.resize(static_cast<int>(P_dof_indices.size()));
        fe.reinit(elem);
        fe.collectDataForInterpolation(elem);
        fe.interpolate(elem);
        const unsigned int n_qp = qrule->n_points();
        const size_t n_basis = phi.size();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            const std::vector<VectorValue<double> >& grad_x_data = fe_interp_grad_var_data[qp][X_sys_idx];
            get_FF(FF, grad_x_data);
            double J = FF.det();
            const double P = (U_prime_fcn ? U_prime_fcn(J) : -d_static_pressure_kappa * std::log(J));
            for (unsigned int k = 0; k < n_basis; ++k)
            {
                P_rhs_e(k) += P * phi[k][qp] * JxW[qp];
            }
        }

        // Apply constraints (e.g., enforce periodic boundary conditions)
        // and add the elemental contributions to the global vector.
        copy_dof_ids_to_vector(/*var_num*/ 0, P_dof_indices, dof_id_scratch);
        P_dof_map.constrain_element_vector(P_rhs_e, dof_id_scratch);
        P_rhs_vec->add_vector(P_rhs_e, dof_id_scratch);
    }

    // Solve for P.
    switch (proj_type)
    {
    case CONSISTENT_PROJECTION:
        d_fe_projectors[part]->computeL2Projection(
            P_vec, *P_rhs_vec, PRESSURE_SYSTEM_NAME, /*use_consistent_mass_matrix*/ true);
        break;
    case LUMPED_PROJECTION:
        d_fe_projectors[part]->computeL2Projection(
            P_vec, *P_rhs_vec, PRESSURE_SYSTEM_NAME, /*use_consistent_mass_matrix*/ false);
        break;
    case STABILIZED_PROJECTION:
        d_fe_projectors[part]->computeStabilizedL2Projection(
            P_vec, *P_rhs_vec, PRESSURE_SYSTEM_NAME, d_static_pressure_stab_param);
        break;
    default:
        TBOX_ERROR("unsupported pressure projection type\n");
    }
}

void
FEMechanicsBase::assembleInteriorForceDensityRHS(PetscVector<double>& F_rhs_vec,
                                                 PetscVector<double>& X_vec,
                                                 PetscVector<double>* P_vec,
                                                 const double data_time,
                                                 const unsigned int part)
{
    const bool using_pressure = P_vec != nullptr;

    // Extract the mesh.
    EquationSystems& equation_systems = *d_equation_systems[part];
    const MeshBase& mesh = equation_systems.get_mesh();
    const BoundaryInfo& boundary_info = *mesh.boundary_info;
    const unsigned int dim = mesh.mesh_dimension();

    // Setup global and elemental right-hand-side vectors.
    auto& F_system = equation_systems.get_system<ExplicitSystem>(FORCE_SYSTEM_NAME);
    // During assembly we sum into ghost regions - this only makes sense if we
    // have a ghosted vector.
    int ierr;
    TBOX_ASSERT(F_rhs_vec.type() == GHOSTED);
    Vec F_rhs_vec_local;
    ierr = VecGhostGetLocalForm(F_rhs_vec.vec(), &F_rhs_vec_local);
    IBTK_CHKERRQ(ierr);
    double* F_rhs_local_soln = nullptr;
    ierr = VecGetArray(F_rhs_vec_local, &F_rhs_local_soln);
    IBTK_CHKERRQ(ierr);
    std::array<DenseVector<double>, NDIM> F_rhs_e;
    std::vector<libMesh::dof_id_type> dof_id_scratch;

    // First handle the stress contributions.  These are handled separately
    // because each stress function may use a different quadrature rule.
    const size_t num_PK1_fcns = d_PK1_stress_fcn_data[part].size();
    for (unsigned int k = 0; k < num_PK1_fcns; ++k)
    {
        if (!d_PK1_stress_fcn_data[part][k].fcn) continue;

        // Extract the FE systems and DOF maps, and setup the FE object.
        const DofMap& F_dof_map = F_system.get_dof_map();
        FEDataManager::SystemDofMapCache& F_dof_map_cache = *d_fe_data[part]->getDofMapCache(FORCE_SYSTEM_NAME);
        FEType F_fe_type = F_dof_map.variable_type(0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(F_dof_map.variable_type(d) == F_fe_type);
        }
        auto& X_system = equation_systems.get_system<ExplicitSystem>(COORDS_SYSTEM_NAME);
        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;

        FEDataInterpolation fe(dim, d_fe_data[part]);
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
        fe.registerSystem(F_system, std::vector<int>(), vars); // compute dphi for the force system
        const size_t X_sys_idx = fe.registerInterpolatedSystem(X_system, vars, vars, &X_vec);
        std::vector<size_t> PK1_fcn_system_idxs;
        fe.setupInterpolatedSystemDataIndexes(
            PK1_fcn_system_idxs, d_PK1_stress_fcn_data[part][k].system_data, &equation_systems);
        fe.init();

        const std::vector<libMesh::Point>& q_point = fe.getQuadraturePoints();
        const std::vector<double>& JxW = fe.getQuadratureWeights();
        const std::vector<std::vector<VectorValue<double> > >& dphi = fe.getDphi(F_fe_type);

        const std::vector<libMesh::Point>& q_point_face = fe.getQuadraturePointsFace();
        const std::vector<double>& JxW_face = fe.getQuadratureWeightsFace();
        const std::vector<libMesh::Point>& normal_face = fe.getNormalsFace();
        const std::vector<std::vector<double> >& phi_face = fe.getPhiFace(F_fe_type);

        const std::vector<std::vector<std::vector<double> > >& fe_interp_var_data = fe.getVarInterpolation();
        const std::vector<std::vector<std::vector<VectorValue<double> > > >& fe_interp_grad_var_data =
            fe.getGradVarInterpolation();

        std::vector<const std::vector<double>*> PK1_var_data;
        std::vector<const std::vector<VectorValue<double> >*> PK1_grad_var_data;

        // Loop over the elements to compute the right-hand side vector.  This
        // is computed via
        //
        //    rhs_k = -int{PP(s,t) grad phi_k(s)}ds + int{PP(s,t) N(s,t)
        //    phi_k(s)}dA(s)
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
            const auto& F_dof_indices = F_dof_map_cache.dof_indices(elem);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                F_rhs_e[d].resize(static_cast<int>(F_dof_indices[d].size()));
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
                for (unsigned int basis_n = 0; basis_n < n_basis; ++basis_n)
                {
                    F_qp = -PP * dphi[basis_n][qp] * JxW[qp];
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        F_rhs_e[i](basis_n) += F_qp(i);
                    }
                }
            }

            // Loop over the element boundaries.
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                // Skip non-physical boundaries.
                if (!is_physical_bdry(elem, side, boundary_info, F_dof_map)) continue;

                // Determine if we need to integrate surface forces along this
                // part of the physical boundary; if not, skip the present side.
                const bool at_dirichlet_bdry = is_dirichlet_bdry(elem, side, boundary_info, F_dof_map);
                const bool integrate_normal_stress = (d_include_normal_stress_in_weak_form && !at_dirichlet_bdry) ||
                                                     (!d_include_normal_stress_in_weak_form && at_dirichlet_bdry);
                const bool integrate_tangential_stress =
                    (d_include_tangential_stress_in_weak_form && !at_dirichlet_bdry) ||
                    (!d_include_tangential_stress_in_weak_form && at_dirichlet_bdry);
                if (!integrate_normal_stress && !integrate_tangential_stress) continue;

                fe.reinit(elem, side);
                fe.interpolate(elem, side);
                const unsigned int n_qp_face = qrule_face->n_points();
                const size_t n_basis_face = phi_face.size();
                for (unsigned int qp = 0; qp < n_qp_face; ++qp)
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

                    if (!integrate_normal_stress)
                    {
                        F -= (F * n) * n; // remove the normal component.
                    }

                    if (!integrate_tangential_stress)
                    {
                        F -= (F - (F * n) * n); // remove the tangential component.
                    }

                    // Add the boundary forces to the right-hand-side vector.
                    for (unsigned int basis_face_n = 0; basis_face_n < n_basis_face; ++basis_face_n)
                    {
                        F_qp = F * phi_face[basis_face_n][qp] * JxW_face[qp];
                        for (unsigned int i = 0; i < NDIM; ++i)
                        {
                            F_rhs_e[i](basis_face_n) += F_qp(i);
                        }
                    }
                }
            }

            // Apply constraints (e.g., enforce periodic boundary conditions)
            // and add the elemental contributions to the global vector.
            for (unsigned int var_n = 0; var_n < NDIM; ++var_n)
            {
                copy_dof_ids_to_vector(var_n, F_dof_indices, dof_id_scratch);
                F_dof_map.constrain_element_vector(F_rhs_e[var_n], dof_id_scratch);
                for (unsigned int j = 0; j < dof_id_scratch.size(); ++j)
                {
                    F_rhs_local_soln[F_rhs_vec.map_global_to_local_index(dof_id_scratch[j])] += F_rhs_e[var_n](j);
                }
            }
        }
    }

    // Now account for any additional force contributions.

    // Extract the FE systems and DOF maps, and setup the FE objects.
    const DofMap& F_dof_map = F_system.get_dof_map();
    FEDataManager::SystemDofMapCache& F_dof_map_cache = *d_fe_data[part]->getDofMapCache(FORCE_SYSTEM_NAME);
    FEType F_fe_type = F_dof_map.variable_type(0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(F_dof_map.variable_type(d) == F_fe_type);
    }
    auto& X_system = equation_systems.get_system<ExplicitSystem>(COORDS_SYSTEM_NAME);
    System* P_system = using_pressure ? &equation_systems.get_system<ExplicitSystem>(PRESSURE_SYSTEM_NAME) : nullptr;
    std::vector<int> vars(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
    std::vector<int> P_vars(1, 0);
    std::vector<int> no_vars;

    FEDataInterpolation fe(dim, d_fe_data[part]);
    std::unique_ptr<QBase> qrule = QBase::build(d_default_quad_type_force[part], dim, d_default_quad_order_force[part]);
    std::unique_ptr<QBase> qrule_face =
        QBase::build(d_default_quad_type_force[part], dim - 1, d_default_quad_order_force[part]);
    fe.attachQuadratureRule(qrule.get());
    fe.attachQuadratureRuleFace(qrule_face.get());
    fe.evalNormalsFace();
    fe.evalQuadraturePoints();
    fe.evalQuadraturePointsFace();
    fe.evalQuadratureWeights();
    fe.evalQuadratureWeightsFace();
    fe.registerSystem(F_system, vars, vars); // compute phi and dphi for the force system
    const size_t X_sys_idx = fe.registerInterpolatedSystem(X_system, vars, vars, &X_vec);
    const size_t P_sys_idx = using_pressure ? fe.registerInterpolatedSystem(*P_system, P_vars, no_vars, P_vec) :
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
    fe.init();

    const std::vector<libMesh::Point>& q_point = fe.getQuadraturePoints();
    const std::vector<double>& JxW = fe.getQuadratureWeights();
    const std::vector<std::vector<double> >& phi = fe.getPhi(F_fe_type);
    const std::vector<std::vector<VectorValue<double> > >& dphi = fe.getDphi(F_fe_type);

    const std::vector<libMesh::Point>& q_point_face = fe.getQuadraturePointsFace();
    const std::vector<double>& JxW_face = fe.getQuadratureWeightsFace();
    const std::vector<libMesh::Point>& normal_face = fe.getNormalsFace();
    const std::vector<std::vector<double> >& phi_face = fe.getPhiFace(F_fe_type);

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
    boost::multi_array<double, 1> P_node;
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;
        const auto& F_dof_indices = F_dof_map_cache.dof_indices(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            F_rhs_e[d].resize(static_cast<int>(F_dof_indices[d].size()));
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

            if (using_pressure)
            {
                const double P = fe_interp_var_data[qp][P_sys_idx][0];

                // Compute the value of the first Piola-Kirchhoff stress tensor
                // at the quadrature point and add the corresponding forces to
                // the right-hand-side vector.
                PP = -J * P * FF_inv_trans;
                for (unsigned int k = 0; k < n_basis; ++k)
                {
                    F_qp = -PP * dphi[k][qp] * JxW[qp];
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        F_rhs_e[i](k) += F_qp(i);
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
                        F_rhs_e[i](k) += F_qp(i);
                    }
                }
            }
        }

        // Loop over the element boundaries.
        for (unsigned int side = 0; side < elem->n_sides(); ++side)
        {
            // Skip non-physical boundaries.
            if (!is_physical_bdry(elem, side, boundary_info, F_dof_map)) continue;

            // Determine if we need to compute surface forces along this
            // part of the physical boundary; if not, skip the present side.
            const bool at_dirichlet_bdry = is_dirichlet_bdry(elem, side, boundary_info, F_dof_map);
            const bool integrate_normal_force = d_include_normal_surface_forces_in_weak_form && !at_dirichlet_bdry;
            const bool integrate_tangential_force =
                d_include_tangential_surface_forces_in_weak_form && !at_dirichlet_bdry;
            if (!integrate_normal_force && !integrate_tangential_force) continue;

            fe.reinit(elem, side);
            fe.interpolate(elem, side);
            const unsigned int n_qp_face = qrule_face->n_points();
            const size_t n_basis_face = phi_face.size();
            for (unsigned int qp = 0; qp < n_qp_face; ++qp)
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
                for (unsigned int k = 0; k < n_basis_face; ++k)
                {
                    F_qp = F * phi_face[k][qp] * JxW_face[qp];
                    for (unsigned int i = 0; i < NDIM; ++i)
                    {
                        F_rhs_e[i](k) += F_qp(i);
                    }
                }
            }
        }

        // Apply constraints (e.g., enforce periodic boundary conditions)
        // and add the elemental contributions to the global vector.
        for (unsigned int var_n = 0; var_n < NDIM; ++var_n)
        {
            copy_dof_ids_to_vector(var_n, F_dof_indices, dof_id_scratch);
            F_dof_map.constrain_element_vector(F_rhs_e[var_n], dof_id_scratch);
            for (unsigned int j = 0; j < dof_id_scratch.size(); ++j)
            {
                F_rhs_local_soln[F_rhs_vec.map_global_to_local_index(dof_id_scratch[j])] += F_rhs_e[var_n](j);
            }
        }
    }

    ierr = VecRestoreArray(F_rhs_vec_local, &F_rhs_local_soln);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostRestoreLocalForm(F_rhs_vec.vec(), &F_rhs_vec_local);
    IBTK_CHKERRQ(ierr);
}

void
FEMechanicsBase::initializeCoordinates(const unsigned int part)
{
    EquationSystems& equation_systems = *d_fe_data[part]->getEquationSystems();
    MeshBase& mesh = equation_systems.get_mesh();
    auto& X_system = equation_systems.get_system<ExplicitSystem>(COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    const bool identity_mapping = !d_coordinate_mapping_fcn_data[part].fcn;
    auto it = mesh.local_nodes_begin();
    const auto end_it = mesh.local_nodes_end();
    for (; it != end_it; ++it)
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
    copy_and_synch(X_coords,
                   *X_system.current_local_solution,
                   /*close_v_in*/ false);
}

void
FEMechanicsBase::updateCoordinateMapping(const unsigned int part)
{
    EquationSystems& equation_systems = *d_fe_data[part]->getEquationSystems();
    MeshBase& mesh = equation_systems.get_mesh();
    auto& X_system = equation_systems.get_system<ExplicitSystem>(COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    auto& dX_system = equation_systems.get_system<ExplicitSystem>(COORD_MAPPING_SYSTEM_NAME);
    const unsigned int dX_sys_num = dX_system.number();
    NumericVector<double>& dX_coords = *dX_system.solution;
    auto it = mesh.local_nodes_begin();
    const auto end_it = mesh.local_nodes_end();
    for (; it != end_it; ++it)
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
}

void
FEMechanicsBase::initializeVelocity(const unsigned int part)
{
    EquationSystems& equation_systems = *d_fe_data[part]->getEquationSystems();
    MeshBase& mesh = equation_systems.get_mesh();
    auto& U_system = equation_systems.get_system<ExplicitSystem>(VELOCITY_SYSTEM_NAME);
    const unsigned int U_sys_num = U_system.number();
    NumericVector<double>& U_vec = *U_system.solution;
    VectorValue<double> U;
    if (!d_initial_velocity_fcn_data[part].fcn)
    {
        U.zero();
    }
    else
    {
        auto it = mesh.local_nodes_begin();
        const auto end_it = mesh.local_nodes_end();
        for (; it != end_it; ++it)
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
}

void
FEMechanicsBase::setup_system_vectors(EquationSystems* equation_systems,
                                      const std::vector<std::string>& system_names,
                                      const std::vector<std::string>& vector_names)
{
    IBAMR_DEPRECATED_MEMBER_FUNCTION1("FEMechanicsBase", "setup_system_vectors");
    IBTK::setup_system_vectors(
        equation_systems, system_names, vector_names, RestartManager::getManager()->isFromRestart());
}

void
FEMechanicsBase::setup_system_vector(System& system, const std::string& vector_name)
{
    IBAMR_DEPRECATED_MEMBER_FUNCTION1("FEMechanicsBase", "setup_system_vector");
    IBTK::setup_system_vector(system, vector_name, RestartManager::getManager()->isFromRestart());
}

std::string
FEMechanicsBase::libmesh_restart_file_name(const std::string& restart_dump_dirname,
                                           unsigned int time_step_number,
                                           unsigned int part,
                                           const std::string& extension)
{
    std::ostringstream file_name_prefix;
    file_name_prefix << restart_dump_dirname << "/libmesh_data_part_" << part << "." << std::setw(6)
                     << std::setfill('0') << std::right << time_step_number << "." << extension;
    return file_name_prefix.str();
}

/////////////////////////////// PRIVATE //////////////////////////////////////

void
FEMechanicsBase::commonConstructor(const std::string& object_name,
                                   const Pointer<Database>& input_db,
                                   const std::vector<libMesh::MeshBase*>& meshes,
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
    const auto n_parts = d_meshes.size();

    // Set some default values.
    d_fe_order_position.resize(n_parts, INVALID_ORDER);
    d_fe_order_force.resize(n_parts, INVALID_ORDER);
    d_fe_order_pressure.resize(n_parts, INVALID_ORDER);
    d_fe_family_position.resize(n_parts, INVALID_FE);
    d_fe_family_force.resize(n_parts, INVALID_FE);
    d_fe_family_pressure.resize(n_parts, INVALID_FE);
    d_default_quad_type_stress.resize(n_parts, INVALID_Q_RULE);
    d_default_quad_type_force.resize(n_parts, INVALID_Q_RULE);
    d_default_quad_type_pressure.resize(n_parts, INVALID_Q_RULE);
    d_default_quad_order_stress.resize(n_parts, INVALID_ORDER);
    d_default_quad_order_force.resize(n_parts, INVALID_ORDER);
    d_default_quad_order_pressure.resize(n_parts, INVALID_ORDER);

    // Initialize function data to NULL.
    d_coordinate_mapping_fcn_data.resize(n_parts);
    d_initial_velocity_fcn_data.resize(n_parts);
    d_PK1_stress_fcn_data.resize(n_parts);
    d_lag_body_force_fcn_data.resize(n_parts);
    d_lag_surface_pressure_fcn_data.resize(n_parts);
    d_lag_surface_force_fcn_data.resize(n_parts);

    // Indicate that all of the parts do NOT use static pressures by default.
    d_static_pressure_part.resize(n_parts, false);
    d_static_pressure_proj_type.resize(n_parts, UNKNOWN_PRESSURE_TYPE);
    d_static_pressure_vol_energy_deriv_fcn.resize(n_parts, nullptr);

    // Determine whether we should use first-order or second-order shape
    // functions for each part of the structure.
    for (unsigned int part = 0; part < n_parts; ++part)
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
        mesh_has_first_order_elems = IBTK_MPI::maxReduction(mesh_has_first_order_elems ? 1 : 0);
        mesh_has_second_order_elems = IBTK_MPI::maxReduction(mesh_has_second_order_elems ? 1 : 0);
        if ((mesh_has_first_order_elems && mesh_has_second_order_elems) ||
            (!mesh_has_first_order_elems && !mesh_has_second_order_elems))
        {
            TBOX_ERROR(d_object_name << "::FEMechanicsBase():\n"
                                     << "  each FE mesh part must contain only FIRST "
                                        "order elements or only SECOND order elements"
                                     << std::endl);
        }
        d_fe_family_position[part] = LAGRANGE;
        d_fe_family_force[part] = LAGRANGE;
        d_fe_family_pressure[part] = LAGRANGE;
        d_default_quad_type_stress[part] = QGAUSS;
        d_default_quad_type_force[part] = QGAUSS;
        d_default_quad_type_pressure[part] = QGAUSS;
        if (mesh_has_first_order_elems)
        {
            d_fe_order_position[part] = FIRST;
            d_fe_order_force[part] = FIRST;
            d_fe_order_pressure[part] = FIRST;
            d_default_quad_order_stress[part] = THIRD;
            d_default_quad_order_force[part] = THIRD;
            d_default_quad_order_pressure[part] = THIRD;
        }
        if (mesh_has_second_order_elems)
        {
            d_fe_order_position[part] = SECOND;
            d_fe_order_force[part] = SECOND;
            d_fe_order_pressure[part] = SECOND;
            d_default_quad_order_stress[part] = FIFTH;
            d_default_quad_order_force[part] = FIFTH;
            d_default_quad_order_pressure[part] = FIFTH;
        }
    }

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    // Report configuration for each part.
    for (unsigned int part = 0; part < n_parts; ++part)
    {
        pout << "\n";
        pout << d_object_name << ": mesh part " << part << " is using ";
        if ((d_fe_family_position[part] == d_fe_family_force[part]) &&
            (d_fe_family_force[part] == d_fe_family_pressure[part]) &&
            (d_fe_order_position[part] == d_fe_order_force[part]) &&
            (d_fe_order_force[part] == d_fe_order_pressure[part]))
        {
            pout << Utility::enum_to_string<Order>(d_fe_order_position[part]) << " order "
                 << Utility::enum_to_string<FEFamily>(d_fe_family_position[part]) << " finite elements.\n";
        }
        else
        {
            pout << "  " << Utility::enum_to_string<Order>(d_fe_order_position[part]) << " order "
                 << Utility::enum_to_string<FEFamily>(d_fe_family_position[part]) << " elements for positions\n"
                 << "  " << Utility::enum_to_string<Order>(d_fe_order_force[part]) << " order "
                 << Utility::enum_to_string<FEFamily>(d_fe_family_pressure[part]) << " elements for forces\n"
                 << "  " << Utility::enum_to_string<Order>(d_fe_order_pressure[part]) << " order "
                 << Utility::enum_to_string<FEFamily>(d_fe_family_pressure[part]) << " elements for pressures\n";
        }
        pout << "\n";
    }
}

void
FEMechanicsBase::getFromInput(const Pointer<Database>& db, bool /*is_from_restart*/)
{
    // libMesh parallelization settings.
    d_libmesh_partitioner_type =
        string_to_enum<LibmeshPartitionerType>(db->getStringWithDefault("libmesh_partitioner_type", "LIBMESH_DEFAULT"));

    // Force computation settings.
    if (db->isBool("use_consistent_mass_matrix"))
        d_use_consistent_mass_matrix = db->getBool("use_consistent_mass_matrix");

    // Pressure settings.
    if (db->isDouble("static_pressure_kappa")) d_static_pressure_kappa = db->getDouble("static_pressure_kappa");
    if (db->isDouble("static_pressure_stab_param"))
        d_static_pressure_stab_param = db->getDouble("static_pressure_stab_param");

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
    if (db->keyExists("do_log"))
        d_do_log = db->getBool("do_log");
    else if (db->keyExists("enable_logging"))
        d_do_log = db->getBool("enable_logging");
}

void
FEMechanicsBase::getFromRestart()
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
    int ver = db->getInteger("FE_MECHANICS_BASE_VERSION");
    if (ver != FE_MECHANICS_BASE_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    d_use_consistent_mass_matrix = db->getBool("d_use_consistent_mass_matrix");
    d_libmesh_partitioner_type = string_to_enum<LibmeshPartitionerType>(db->getString("d_libmesh_partitioner_type"));
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
