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

#include <IBTK_config.h>

#include "ibamr/FEMechanicsExplicitIntegrator.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep

#include "ibtk/FEDataInterpolation.h"
#include "ibtk/FEDataManager.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/LibMeshSystemVectors.h"
#include "ibtk/QuadratureCache.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/libmesh_utilities.h"

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

#include "petscvec.h"

using namespace libMesh;

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of FEMechanicsExplicitIntegrator restart file data.
static const int FE_MECHANICS_EXPLICIT_INTEGRATOR_VERSION = 0;

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

const std::string FEMechanicsExplicitIntegrator::COORDS_SYSTEM_NAME = "current coordinates system";
const std::string FEMechanicsExplicitIntegrator::COORD_MAPPING_SYSTEM_NAME = "displacement system";
const std::string FEMechanicsExplicitIntegrator::FORCE_SYSTEM_NAME = "force system";
const std::string FEMechanicsExplicitIntegrator::PRESSURE_SYSTEM_NAME = "pressure system";
const std::string FEMechanicsExplicitIntegrator::VELOCITY_SYSTEM_NAME = "velocity system";

/////////////////////////////// PUBLIC ///////////////////////////////////////

FEMechanicsExplicitIntegrator::FEMechanicsExplicitIntegrator(const std::string& object_name,
                                                             Pointer<Database> input_db,
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

FEMechanicsExplicitIntegrator::FEMechanicsExplicitIntegrator(const std::string& object_name,
                                                             Pointer<Database> input_db,
                                                             const std::vector<MeshBase*>& meshes,
                                                             bool register_for_restart,
                                                             const std::string& restart_read_dirname,
                                                             unsigned int restart_restore_number)
    : d_num_parts(static_cast<int>(meshes.size()))
{
    commonConstructor(
        object_name, input_db, meshes, register_for_restart, restart_read_dirname, restart_restore_number);
}

FEMechanicsExplicitIntegrator::~FEMechanicsExplicitIntegrator()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }
}

FEData*
FEMechanicsExplicitIntegrator::getFEData(unsigned int part) const
{
    TBOX_ASSERT(part < d_num_parts);
    return d_fe_data[part].get();
}

void
FEMechanicsExplicitIntegrator::registerInitialCoordinateMappingFunction(const CoordinateMappingFcnData& data,
                                                                        const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_coordinate_mapping_fcn_data[part] = data;
}

FEMechanicsExplicitIntegrator::CoordinateMappingFcnData
FEMechanicsExplicitIntegrator::getInitialCoordinateMappingFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_num_parts);
    return d_coordinate_mapping_fcn_data[part];
}

void
FEMechanicsExplicitIntegrator::registerInitialVelocityFunction(const InitialVelocityFcnData& data,
                                                               const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_initial_velocity_fcn_data[part] = data;
}

FEMechanicsExplicitIntegrator::InitialVelocityFcnData
FEMechanicsExplicitIntegrator::getInitialVelocityFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_num_parts);
    return d_initial_velocity_fcn_data[part];
}

void
FEMechanicsExplicitIntegrator::registerPK1StressFunction(const PK1StressFcnData& data, const unsigned int part)
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
}

std::vector<FEMechanicsExplicitIntegrator::PK1StressFcnData>
FEMechanicsExplicitIntegrator::getPK1StressFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_num_parts);
    return d_PK1_stress_fcn_data[part];
}

void
FEMechanicsExplicitIntegrator::registerLagBodyForceFunction(const LagBodyForceFcnData& data, const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_lag_body_force_fcn_data[part] = data;
}

FEMechanicsExplicitIntegrator::LagBodyForceFcnData
FEMechanicsExplicitIntegrator::getLagBodyForceFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_num_parts);
    return d_lag_body_force_fcn_data[part];
}

void
FEMechanicsExplicitIntegrator::registerLagSurfacePressureFunction(const LagSurfacePressureFcnData& data,
                                                                  const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_lag_surface_pressure_fcn_data[part] = data;
}

FEMechanicsExplicitIntegrator::LagSurfacePressureFcnData
FEMechanicsExplicitIntegrator::getLagSurfacePressureFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_num_parts);
    return d_lag_surface_pressure_fcn_data[part];
}

void
FEMechanicsExplicitIntegrator::registerLagSurfaceForceFunction(const LagSurfaceForceFcnData& data,
                                                               const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_lag_surface_force_fcn_data[part] = data;
}

FEMechanicsExplicitIntegrator::LagSurfaceForceFcnData
FEMechanicsExplicitIntegrator::getLagSurfaceForceFunction(unsigned int part) const
{
    TBOX_ASSERT(part < d_num_parts);
    return d_lag_surface_force_fcn_data[part];
}

void
FEMechanicsExplicitIntegrator::preprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    d_started_time_integration = true;

    d_current_time = current_time;
    d_new_time = new_time;
    d_half_time = current_time + 0.5 * (new_time - current_time);

    // Initialize vector data.
    d_X_vecs->copy("solution", { "current", "new" });
    d_U_vecs->copy("solution", { "current", "new" });
    d_F_vecs->copy("solution", { "current", "new" });
    d_P_vecs->copy("solution", { "current", "new" });
}

void
FEMechanicsExplicitIntegrator::postprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    batch_vec_ghost_update({ d_X_vecs->get("new"), d_U_vecs->get("new"), d_F_vecs->get("new"), d_P_vecs->get("new") },
                           INSERT_VALUES,
                           SCATTER_FORWARD);

    d_X_vecs->copy("new", { "solution", "current" });
    d_U_vecs->copy("new", { "solution", "current" });
    d_F_vecs->copy("new", { "solution", "current" });
    d_P_vecs->copy("new", { "solution", "current" });

    // Update the coordinate mapping dX = X - X0.
    for (unsigned part = 0; part < d_num_parts; ++part)
    {
        updateCoordinateMapping(part);
    }

    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time = std::numeric_limits<double>::quiet_NaN();
    d_half_time = std::numeric_limits<double>::quiet_NaN();
}

void
FEMechanicsExplicitIntegrator::forwardEulerStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        computeLagrangianForce(current_time);

        int ierr;
        ierr = VecAXPBY(
            d_U_vecs->get("current", part).vec(), 1.0 / d_rho0[part], 0.0, d_F_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecWAXPY(d_X_vecs->get("new", part).vec(),
                        dt,
                        d_U_vecs->get("current", part).vec(),
                        d_X_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(d_X_vecs->get("half", part).vec(),
                           0.5,
                           0.5,
                           0.0,
                           d_X_vecs->get("current", part).vec(),
                           d_X_vecs->get("new", part).vec());
        IBTK_CHKERRQ(ierr);
    }
}

void
FEMechanicsExplicitIntegrator::computeLagrangianForce(const double data_time)
{
    const std::string data_time_str = std::move(get_data_time_str(data_time, d_current_time, d_new_time));
    batch_vec_ghost_update(d_X_vecs->get(data_time_str), INSERT_VALUES, SCATTER_FORWARD);
    d_F_vecs->zero("RHS Vector");
    d_F_vecs->zero("tmp");
    for (unsigned part = 0; part < d_num_parts; ++part)
    {
        assembleInteriorForceDensityRHS(d_F_vecs->get("RHS Vector", part),
                                        d_X_vecs->get(data_time_str, part),
                                        d_U_vecs->get(data_time_str, part),
                                        d_P_vecs->get(data_time_str, part),
                                        data_time,
                                        part);
    }
    batch_vec_ghost_update(d_F_vecs->get("RHS Vector"), ADD_VALUES, SCATTER_REVERSE);

    // RHSs don't need ghost entries for the solve so we can skip the other scatter
    for (unsigned part = 0; part < d_num_parts; ++part)
    {
        d_fe_projector[part]->computeL2Projection(d_F_vecs->get("solution", part),
                                                  d_F_vecs->get("RHS Vector", part),
                                                  FORCE_SYSTEM_NAME,
                                                  d_use_consistent_mass_matrix,
                                                  /*close_U*/ false,
                                                  /*close_F*/ false);
        d_F_vecs->copy("solution", { data_time_str });
    }
}

void
FEMechanicsExplicitIntegrator::initializeFEEquationSystems()
{
    if (d_fe_equation_systems_initialized) return;

    const bool from_restart = RestartManager::getManager()->isFromRestart();

    // Set up the coupling matrix which will be used by each system.
    d_diagonal_system_coupling.resize(NDIM);
    for (unsigned int i = 0; i < NDIM; ++i)
        for (unsigned int j = 0; j < NDIM; ++j) d_diagonal_system_coupling(i, j) = i == j ? 1 : 0;

    // Create the EquationSystems objects
    d_equation_systems.resize(d_num_parts);
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        // Create FE equation systems objects and corresponding variables.
        d_equation_systems[part] = std::unique_ptr<EquationSystems>(new EquationSystems(*d_meshes[part]));
        EquationSystems& equation_systems = *d_equation_systems[part];
        if (from_restart)
        {
            const std::string& file_name = libmesh_restart_file_name(
                d_libmesh_restart_read_dir, d_libmesh_restart_restore_number, part, d_libmesh_restart_file_extension);
            const XdrMODE xdr_mode = (d_libmesh_restart_file_extension == "xdr" ? DECODE : READ);
            const int read_mode =
                EquationSystems::READ_HEADER | EquationSystems::READ_DATA | EquationSystems::READ_ADDITIONAL_DATA;
            equation_systems.read(file_name,
                                  xdr_mode,
                                  read_mode,
                                  /*partition_agnostic*/ true);
        }
        else
        {
            auto& X_system = equation_systems.add_system<ExplicitSystem>(COORDS_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_system.add_variable("X_" + std::to_string(d), d_fe_order[part], d_fe_family[part]);
            }
            X_system.get_dof_map()._dof_coupling = &d_diagonal_system_coupling;

            auto& dX_system = equation_systems.add_system<ExplicitSystem>(COORD_MAPPING_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dX_system.add_variable("dX_" + std::to_string(d), d_fe_order[part], d_fe_family[part]);
            }
            dX_system.get_dof_map()._dof_coupling = &d_diagonal_system_coupling;

            auto& U_system = equation_systems.add_system<ExplicitSystem>(VELOCITY_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                U_system.add_variable("U_" + std::to_string(d), d_fe_order[part], d_fe_family[part]);
            }
            U_system.get_dof_map()._dof_coupling = &d_diagonal_system_coupling;

            auto& F_system = equation_systems.add_system<ExplicitSystem>(FORCE_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                F_system.add_variable("F_" + std::to_string(d), d_fe_order[part], d_fe_family[part]);
            }
            F_system.get_dof_map()._dof_coupling = &d_diagonal_system_coupling;

            auto& P_system = equation_systems.add_system<ExplicitSystem>(PRESSURE_SYSTEM_NAME);
            P_system.add_variable("P", d_fe_order[part], d_fe_family[part]);
        }

        // Setup cached system vectors.
        //
        // NOTE: libMesh does not appear to preserve the type of the vector
        // after restart (GHOSTED vectors are now PARALLEL), and so we
        // manually reset these vectors here.
        auto insert_parallel_into_ghosted = [](const PetscVector<Number>& parallel_vector,
                                               PetscVector<Number>& ghosted_vector) {
            TBOX_ASSERT(parallel_vector.size() == ghosted_vector.size());
            TBOX_ASSERT(parallel_vector.local_size() == ghosted_vector.local_size());
            ghosted_vector = parallel_vector;
            ghosted_vector.close();
        };

        auto convert_parallel_to_ghosted = [&](const std::string& vector_name, System& system) {
            std::unique_ptr<NumericVector<double> > clone_vector;
            if (from_restart)
            {
                NumericVector<double>* current = system.request_vector(vector_name);
                if (current != nullptr)
                {
                    clone_vector = current->clone();
                }
            }
            system.remove_vector(vector_name);
            system.add_vector(vector_name, /*projections*/ true, /*type*/ GHOSTED);

            if (clone_vector != nullptr)
            {
                const auto& parallel_vector = dynamic_cast<const PetscVector<Number>&>(*clone_vector);
                auto& ghosted_vector = dynamic_cast<PetscVector<Number>&>(system.get_vector(vector_name));
                insert_parallel_into_ghosted(parallel_vector, ghosted_vector);
            }
        };

        const std::array<std::string, 2> system_names{ { COORDS_SYSTEM_NAME, VELOCITY_SYSTEM_NAME } };
        const std::array<std::string, 2> vector_names{ { "new", "half" } };
        for (const std::string& system_name : system_names)
        {
            auto& system = equation_systems.get_system(system_name);
            for (const std::string& vector_name : vector_names)
            {
                convert_parallel_to_ghosted(vector_name, system);
            }
        }

        {
            auto& system = equation_systems.get_system<ExplicitSystem>(FORCE_SYSTEM_NAME);
            const std::array<std::string, 3> vector_names{ { "half", "tmp", "RHS Vector" } };
            for (const std::string& vector_name : vector_names)
            {
                convert_parallel_to_ghosted(vector_name, system);
            }
            // libMesh also stores a convenience pointer to the RHS that we need to reset:
            system.rhs = &system.get_vector("RHS Vector");
        }

        {
            auto& system = equation_systems.get_system<ExplicitSystem>(PRESSURE_SYSTEM_NAME);
            const std::array<std::string, 3> vector_names{ { "half", "tmp", "RHS Vector" } };
            for (const std::string& vector_name : vector_names)
            {
                convert_parallel_to_ghosted(vector_name, system);
            }
            // libMesh also stores a convenience pointer to the RHS that we need to reset:
            system.rhs = &system.get_vector("RHS Vector");
        }
    }
    d_fe_equation_systems_initialized = true;
}

void
FEMechanicsExplicitIntegrator::initializeFEData()
{
    if (d_fe_data_initialized) return;

    initializeFEEquationSystems();
    doInitializeFEData(RestartManager::getManager()->isFromRestart());
    d_fe_data_initialized = true;
}

void
FEMechanicsExplicitIntegrator::reinitializeFEData()
{
    TBOX_ASSERT(d_fe_data_initialized);
    doInitializeFEData(true);
}

void
FEMechanicsExplicitIntegrator::doInitializeFEData(const bool use_present_data)
{
    // The choice of FEDataManager set is important here since it determines the
    // IB ghost regions.
    std::vector<EquationSystems*> equation_systems;
    std::transform(d_equation_systems.begin(),
                   d_equation_systems.end(),
                   std::back_inserter(equation_systems),
                   [](const std::unique_ptr<EquationSystems>& es) { return es.get(); });

    d_X_vecs.reset(new LibMeshSystemVectors(equation_systems, COORDS_SYSTEM_NAME));
    d_U_vecs.reset(new LibMeshSystemVectors(equation_systems, VELOCITY_SYSTEM_NAME));
    d_F_vecs.reset(new LibMeshSystemVectors(equation_systems, FORCE_SYSTEM_NAME));
    d_P_vecs.reset(new LibMeshSystemVectors(equation_systems, PRESSURE_SYSTEM_NAME));

    d_fe_data.resize(d_num_parts);
    d_fe_projector.resize(d_num_parts);
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
        auto& P_system = equation_systems.get_system<System>(PRESSURE_SYSTEM_NAME);

        X_system.assemble_before_solve = false;
        X_system.assemble();

        dX_system.assemble_before_solve = false;
        dX_system.assemble();

        U_system.assemble_before_solve = false;
        U_system.assemble();

        F_system.assemble_before_solve = false;
        F_system.assemble();

        P_system.assemble_before_solve = false;
        P_system.assemble();

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
                const bool at_mesh_bdry = !elem->neighbor_ptr(side);
                if (!at_mesh_bdry) continue;

                static const std::array<boundary_id_type,3> dirichlet_bdry_id_set = { FEDataManager::ZERO_DISPLACEMENT_X_BDRY_ID,
                                                                           FEDataManager::ZERO_DISPLACEMENT_Y_BDRY_ID,
                                                                           FEDataManager::ZERO_DISPLACEMENT_Z_BDRY_ID };
                std::vector<boundary_id_type> bdry_ids;
                mesh.boundary_info->boundary_ids(elem, side, bdry_ids);
                const boundary_id_type dirichlet_bdry_ids = get_dirichlet_bdry_ids(bdry_ids);
                if (!dirichlet_bdry_ids) continue;

                for (unsigned int n = 0; n < elem->n_nodes(); ++n)
                {
                    if (!elem->is_node_on_side(n, side)) continue;

                    const Node* const node = elem->node_ptr(n);
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

        d_fe_data[part] = std::make_shared<FEData>(d_object_name + "::FE data", d_registered_for_restart);
        d_fe_data[part]->setEquationSystems(d_equation_systems[part].get(), /*level_number*/ 0);
        d_fe_projector[part] = std::make_shared<FEProjector>(d_fe_data[part], /*enable_logging*/ true);
    }
}

void
FEMechanicsExplicitIntegrator::putToDatabase(Pointer<Database> db)
{
    db->putInteger("FE_MECHANICS_EXPLICIT_INTEGRATOR_VERSION", FE_MECHANICS_EXPLICIT_INTEGRATOR_VERSION);
    db->putInteger("d_num_parts", d_num_parts);
    db->putBool("d_use_consistent_mass_matrix", d_use_consistent_mass_matrix);
    db->putString("d_libmesh_partitioner_type", enum_to_string<LibmeshPartitionerType>(d_libmesh_partitioner_type));
}

void
FEMechanicsExplicitIntegrator::writeFEDataToRestartFile(const std::string& restart_dump_dirname,
                                                        unsigned int time_step_number)
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
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
FEMechanicsExplicitIntegrator::assembleInteriorForceDensityRHS(PetscVector<double>& F_rhs_vec,
                                                               PetscVector<double>& X_vec,
                                                               PetscVector<double>& U_vec,
                                                               PetscVector<double>& P_vec,
                                                               const double data_time,
                                                               const unsigned int part)
{
    // Extract the mesh.
    EquationSystems& equation_systems = *d_fe_data[part]->getEquationSystems();
    const MeshBase& mesh = equation_systems.get_mesh();
    const BoundaryInfo& boundary_info = *mesh.boundary_info;
    const unsigned int dim = mesh.mesh_dimension();

    // Setup global and elemental right-hand-side vectors.
    auto& F_system = equation_systems.get_system<ExplicitSystem>(FORCE_SYSTEM_NAME);

    // During assembly we sum into ghost regions - this only makes sense if we
    // have a ghosted vector.
    TBOX_ASSERT(F_rhs_vec.type() == GHOSTED);
    Vec F_rhs_vec_local;
    int ierr = VecGhostGetLocalForm(F_rhs_vec.vec(), &F_rhs_vec_local);
    IBTK_CHKERRQ(ierr);
    double* F_rhs_local_soln = nullptr;
    ierr = VecGetArray(F_rhs_vec_local, &F_rhs_local_soln);
    IBTK_CHKERRQ(ierr);
    std::array<DenseVector<double>,NDIM> F_rhs_e;
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
        fe.registerSystem(F_system, std::vector<int>(),
                          vars); // compute dphi for the force system
        const size_t X_sys_idx = fe.registerInterpolatedSystem(X_system, vars, vars, &X_vec);
        std::vector<size_t> PK1_fcn_system_idxs;
        fe.setupInterpolatedSystemDataIndexes(
            PK1_fcn_system_idxs, d_PK1_stress_fcn_data[part][k].system_data, &equation_systems);
        fe.init(/*use_IB_ghosted_vecs*/ false);

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

            // Apply constraints (e.g., enforce periodic boundary conditions)
            // and add the elemental contributions to the global vector.
            for (unsigned int i = 0; i < NDIM; ++i)
            {
                dof_id_scratch = F_dof_indices[i];
                F_dof_map.constrain_element_vector(F_rhs_e[i], dof_id_scratch);
                for (unsigned int j = 0; j < dof_id_scratch.size(); ++j)
                {
                    F_rhs_local_soln[F_rhs_vec.map_global_to_local_index(dof_id_scratch[j])] += F_rhs_e[i](j);
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
    std::vector<int> vars(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
    std::vector<int> no_vars;

    FEDataInterpolation fe(dim, d_fe_data[part]);
    std::unique_ptr<QBase> qrule = QBase::build(d_default_quad_type[part], dim, d_default_quad_order[part]);
    std::unique_ptr<QBase> qrule_face = QBase::build(d_default_quad_type[part], dim - 1, d_default_quad_order[part]);
    fe.attachQuadratureRule(qrule.get());
    fe.attachQuadratureRuleFace(qrule_face.get());
    fe.evalNormalsFace();
    fe.evalQuadraturePoints();
    fe.evalQuadraturePointsFace();
    fe.evalQuadratureWeights();
    fe.evalQuadratureWeightsFace();
    fe.registerSystem(F_system, vars,
                      vars); // compute phi and dphi for the force system
    const size_t X_sys_idx = fe.registerInterpolatedSystem(X_system, vars, vars, &X_vec);
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
    boost::multi_array<double, 1> Phi_node;
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
        for (unsigned short int side = 0; side < elem->n_sides(); ++side)
        {
            // Skip non-physical boundaries.
            if (!is_physical_bdry(elem, side, boundary_info, F_dof_map)) continue;

            // Determine if we need to compute surface forces along this
            // part of the physical boundary; if not, skip the present side.
            const bool at_dirichlet_bdry = is_dirichlet_bdry(elem, side, boundary_info, F_dof_map);
            if (at_dirichlet_bdry) continue;

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
        for (unsigned int i = 0; i < NDIM; ++i)
        {
            dof_id_scratch = F_dof_indices[i];
            F_dof_map.constrain_element_vector(F_rhs_e[i], dof_id_scratch);
            for (unsigned int j = 0; j < dof_id_scratch.size(); ++j)
            {
                F_rhs_local_soln[F_rhs_vec.map_global_to_local_index(dof_id_scratch[j])] += F_rhs_e[i](j);
            }
        }
    }

    ierr = VecRestoreArray(F_rhs_vec_local, &F_rhs_local_soln);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostRestoreLocalForm(F_rhs_vec.vec(), &F_rhs_vec_local);
    IBTK_CHKERRQ(ierr);
}

void
FEMechanicsExplicitIntegrator::initializeCoordinates(const unsigned int part)
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
FEMechanicsExplicitIntegrator::updateCoordinateMapping(const unsigned int part)
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
FEMechanicsExplicitIntegrator::initializeVelocity(const unsigned int part)
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

/////////////////////////////// PRIVATE //////////////////////////////////////

void
FEMechanicsExplicitIntegrator::commonConstructor(const std::string& object_name,
                                                 Pointer<Database> input_db,
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

    // Set some default values.
    d_rho0.resize(d_num_parts, 1.0);
    d_fe_family.resize(d_num_parts, INVALID_FE);
    d_fe_order.resize(d_num_parts, INVALID_ORDER);
    d_default_quad_type.resize(d_num_parts, INVALID_Q_RULE);
    d_default_quad_order.resize(d_num_parts, INVALID_ORDER);

    // Initialize function data to NULL.
    d_coordinate_mapping_fcn_data.resize(d_num_parts);
    d_initial_velocity_fcn_data.resize(d_num_parts);
    d_PK1_stress_fcn_data.resize(d_num_parts);
    d_lag_body_force_fcn_data.resize(d_num_parts);
    d_lag_surface_pressure_fcn_data.resize(d_num_parts);
    d_lag_surface_force_fcn_data.resize(d_num_parts);

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
            TBOX_ERROR(d_object_name << "::FEMechanicsExplicitIntegrator():\n"
                                     << "  each FE mesh part must contain only FIRST "
                                        "order elements or only SECOND order elements"
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
}

void
FEMechanicsExplicitIntegrator::getFromInput(Pointer<Database> db, bool /*is_from_restart*/)
{
    // Force computation settings.
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
    if (db->keyExists("do_log"))
        d_do_log = db->getBool("do_log");
    else if (db->keyExists("enable_logging"))
        d_do_log = db->getBool("enable_logging");

    d_libmesh_partitioner_type =
        string_to_enum<LibmeshPartitionerType>(db->getStringWithDefault("libmesh_partitioner_type", "LIBMESH_DEFAULT"));
}

void
FEMechanicsExplicitIntegrator::getFromRestart()
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
    int ver = db->getInteger("FE_MECHANICS_EXPLICIT_INTEGRATOR_VERSION");
    if (ver != FE_MECHANICS_EXPLICIT_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    d_use_consistent_mass_matrix = db->getBool("d_use_consistent_mass_matrix");
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
