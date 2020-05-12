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

#include "ibamr/FEMechanicsBase.h"
#include "ibamr/FEMechanicsExplicitIntegrator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep

#include "ibtk/FEProjector.h"

#include "libmesh/boundary_info.h"
#include "libmesh/enum_xdr_mode.h"
#include "libmesh/equation_systems.h"
#include "libmesh/quadrature_gauss.h"

#include <algorithm>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of FEMechanicsExplicitIntegrator restart file data.
const int EXPLICIT_FE_MECHANICS_INTEGRATOR_VERSION = 0;

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
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

FEMechanicsExplicitIntegrator::FEMechanicsExplicitIntegrator(const std::string& object_name,
                                                             const Pointer<Database>& input_db,
                                                             MeshBase* mesh,
                                                             bool register_for_restart,
                                                             const std::string& restart_read_dirname,
                                                             unsigned int restart_restore_number)
    : FEMechanicsBase(object_name, input_db, mesh, register_for_restart, restart_read_dirname, restart_restore_number)
{
    commonConstructor(object_name, input_db, { mesh }, restart_read_dirname, restart_restore_number);
}

FEMechanicsExplicitIntegrator::FEMechanicsExplicitIntegrator(const std::string& object_name,
                                                             const Pointer<Database>& input_db,
                                                             const std::vector<MeshBase*>& meshes,
                                                             bool register_for_restart,
                                                             const std::string& restart_read_dirname,
                                                             unsigned int restart_restore_number)
    : FEMechanicsBase(object_name, input_db, meshes, register_for_restart, restart_read_dirname, restart_restore_number)
{
    commonConstructor(object_name, input_db, meshes, restart_read_dirname, restart_restore_number);
}

void
FEMechanicsExplicitIntegrator::preprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    FEMechanicsBase::preprocessIntegrateData(current_time, new_time, num_cycles);

    // Initialize variables.
    d_X_vecs->copy("solution", { "current", "new", "half" });
    d_U_vecs->copy("solution", { "current", "new", "half" });
    d_F_vecs->copy("solution", { "current", "new", "half" });
    if (d_P_vecs) d_P_vecs->copy("solution", { "current", "new", "half" });
}

void
FEMechanicsExplicitIntegrator::postprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    std::vector<std::vector<PetscVector<double>*> > vecs{ d_X_vecs->get("new"),
                                                          d_U_vecs->get("new"),
                                                          d_F_vecs->get("new") };
    if (d_P_vecs) vecs.push_back(d_P_vecs->get("new"));
    batch_vec_ghost_update(vecs, INSERT_VALUES, SCATTER_FORWARD);

    d_X_vecs->copy("new", { "solution", "current" });
    d_U_vecs->copy("new", { "solution", "current" });
    d_F_vecs->copy("new", { "solution", "current" });
    if (d_P_vecs) d_P_vecs->copy("new", { "solution", "current" });

    FEMechanicsBase::postprocessIntegrateData(current_time, new_time, num_cycles);
}

void
FEMechanicsExplicitIntegrator::forwardEulerStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        computeLagrangianForce(current_time);
        int ierr;
        ierr = VecWAXPY(d_U_vecs->get("new", part).vec(),
                        dt / d_rhos[part],
                        d_F_vecs->get("current", part).vec(),
                        d_U_vecs->get("current", part).vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(d_U_vecs->get("half", part).vec(),
                           0.5,
                           0.5,
                           0.0,
                           d_U_vecs->get("current", part).vec(),
                           d_U_vecs->get("new", part).vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecWAXPY(d_X_vecs->get("new", part).vec(),
                        dt,
                        d_U_vecs->get("new", part).vec(),
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
    d_F_vecs->copy("current", { "new", "half" });
}

void
FEMechanicsExplicitIntegrator::computeLagrangianForce(const double data_time)
{
    const std::string data_time_str = get_data_time_str(data_time, d_current_time, d_new_time);
    batch_vec_ghost_update(d_X_vecs->get(data_time_str), INSERT_VALUES, SCATTER_FORWARD);
    d_F_vecs->zero("RHS Vector");
    d_F_vecs->zero("tmp");
    for (unsigned part = 0; part < d_meshes.size(); ++part)
    {
        assembleInteriorForceDensityRHS(d_F_vecs->get("RHS Vector", part),
                                        d_X_vecs->get(data_time_str, part),
                                        d_P_vecs ? &d_P_vecs->get(data_time_str, part) : nullptr,
                                        data_time,
                                        part);
    }
    batch_vec_ghost_update(d_F_vecs->get("RHS Vector"), ADD_VALUES, SCATTER_REVERSE);
    for (unsigned part = 0; part < d_meshes.size(); ++part)
    {
        d_fe_projectors[part]->computeL2Projection(d_F_vecs->get("solution", part),
                                                   d_F_vecs->get("RHS Vector", part),
                                                   FORCE_SYSTEM_NAME,
                                                   d_use_consistent_mass_matrix,
                                                   /*close_U*/ false,
                                                   /*close_F*/ false);
    }
    d_F_vecs->copy("solution", { data_time_str });
}

void
FEMechanicsExplicitIntegrator::putToDatabase(Pointer<Database> db)
{
    FEMechanicsBase::putToDatabase(db);
    db->putInteger("EXPLICIT_FE_MECHANICS_INTEGRATOR_VERSION", EXPLICIT_FE_MECHANICS_INTEGRATOR_VERSION);
}

/////////////////////////////// PROTECTED ////////////////////////////////////

void
FEMechanicsExplicitIntegrator::doInitializeFEEquationSystems()
{
    const bool from_restart = RestartManager::getManager()->isFromRestart();

    // Create the FE data managers that manage mappings between the FE mesh
    // parts and the Cartesian grid.
    d_equation_systems.resize(d_meshes.size());
    d_fe_data.resize(d_meshes.size());
    d_fe_projectors.resize(d_meshes.size());
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
    {
        d_fe_data[part] =
            std::make_shared<FEData>(d_object_name + "::FEData::" + std::to_string(part), d_registered_for_restart);
        d_fe_projectors[part].reset(new FEProjector(d_fe_data[part]));

        // Create FE equation systems objects and corresponding variables.
        d_equation_systems[part].reset(new EquationSystems(*d_meshes[part]));
        EquationSystems& equation_systems = *d_equation_systems[part];
        d_fe_data[part]->setEquationSystems(d_equation_systems[part].get(),
                                            /*level_number*/ 0); // TODO: remove level_number as a required parameter
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
                X_system.add_variable("X_" + std::to_string(d), d_fe_order_position[part], d_fe_family_position[part]);
            }
            X_system.get_dof_map()._dof_coupling = &d_diagonal_system_coupling;

            auto& dX_system = equation_systems.add_system<ExplicitSystem>(COORD_MAPPING_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dX_system.add_variable(
                    "dX_" + std::to_string(d), d_fe_order_position[part], d_fe_family_position[part]);
            }
            dX_system.get_dof_map()._dof_coupling = &d_diagonal_system_coupling;

            auto& U_system = equation_systems.add_system<ExplicitSystem>(VELOCITY_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                U_system.add_variable("U_" + std::to_string(d), d_fe_order_position[part], d_fe_family_position[part]);
            }
            U_system.get_dof_map()._dof_coupling = &d_diagonal_system_coupling;

            auto& F_system = equation_systems.add_system<ExplicitSystem>(FORCE_SYSTEM_NAME);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                F_system.add_variable("F_" + std::to_string(d), d_fe_order_force[part], d_fe_family_force[part]);
            }
            F_system.get_dof_map()._dof_coupling = &d_diagonal_system_coupling;
        }

        setup_system_vectors(
            &equation_systems, { COORDS_SYSTEM_NAME, VELOCITY_SYSTEM_NAME }, { "current", "half", "new" });
        setup_system_vectors(
            &equation_systems, { FORCE_SYSTEM_NAME }, { "current", "half", "new", "tmp", "RHS Vector" });
    }
}

void
FEMechanicsExplicitIntegrator::doInitializeFEData(const bool use_present_data)
{
    // The choice of FEDataManager set is important here since it determines the
    // IB ghost regions.
    std::vector<EquationSystems*> equation_systems;
    for (const auto& es : d_equation_systems) equation_systems.push_back(es.get());
    d_X_vecs.reset(new LibMeshSystemVectors(equation_systems, COORDS_SYSTEM_NAME));
    d_U_vecs.reset(new LibMeshSystemVectors(equation_systems, VELOCITY_SYSTEM_NAME));
    d_F_vecs.reset(new LibMeshSystemVectors(equation_systems, FORCE_SYSTEM_NAME));
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
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

                static const std::array<boundary_id_type, 3> dirichlet_bdry_id_set{
                    FEDataManager::ZERO_DISPLACEMENT_X_BDRY_ID,
                    FEDataManager::ZERO_DISPLACEMENT_Y_BDRY_ID,
                    FEDataManager::ZERO_DISPLACEMENT_Z_BDRY_ID
                };
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
    }
}

/////////////////////////////// PRIVATE //////////////////////////////////////

void
FEMechanicsExplicitIntegrator::commonConstructor(const std::string& object_name,
                                                 const Pointer<Database>& input_db,
                                                 const std::vector<libMesh::MeshBase*>& meshes,
                                                 const std::string& restart_read_dirname,
                                                 unsigned int restart_restore_number)
{
    // Set the object name and register it with the restart manager.
    d_object_name = object_name;
    d_registered_for_restart = false;
    d_libmesh_restart_read_dir = restart_read_dirname;
    d_libmesh_restart_restore_number = restart_restore_number;

    // Store the mesh pointers.
    d_meshes = meshes;

    // Initialize function data to NULL.
    d_coordinate_mapping_fcn_data.resize(d_meshes.size());
    d_initial_velocity_fcn_data.resize(d_meshes.size());
    d_PK1_stress_fcn_data.resize(d_meshes.size());
    d_lag_body_force_fcn_data.resize(d_meshes.size());
    d_lag_surface_pressure_fcn_data.resize(d_meshes.size());
    d_lag_surface_force_fcn_data.resize(d_meshes.size());

    // Set some default values.
    d_rhos.resize(d_meshes.size(), 1.0);
    d_libmesh_restart_file_extension = "xdr";
    d_libmesh_partitioner_type = LIBMESH_DEFAULT;

    // Determine whether we should use first-order or second-order shape
    // functions for each part of the structure.
    for (unsigned int part = 0; part < d_meshes.size(); ++part)
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
        d_fe_family_position[part] = LAGRANGE;
        d_fe_family_force[part] = LAGRANGE;
        d_fe_family_pressure[part] = LAGRANGE;
        d_default_quad_type_stress[part] = QGAUSS;
        d_default_quad_type_force[part] = QGAUSS;
        if (mesh_has_first_order_elems)
        {
            d_fe_order_position[part] = FIRST;
            d_fe_order_force[part] = FIRST;
            d_fe_order_pressure[part] = FIRST;
            d_default_quad_order_stress[part] = THIRD;
            d_default_quad_order_stress[part] = THIRD;
        }
        else if (mesh_has_second_order_elems)
        {
            d_fe_order_position[part] = SECOND;
            d_fe_order_force[part] = SECOND;
            d_fe_order_pressure[part] = SECOND;
            d_default_quad_order_stress[part] = FIFTH;
            d_default_quad_order_stress[part] = FIFTH;
        }
    }

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    // Determine how to compute weak forms.
    d_include_normal_stress_in_weak_form = false;
    d_include_tangential_stress_in_weak_form = false;
    d_include_normal_surface_forces_in_weak_form = true;
    d_include_tangential_surface_forces_in_weak_form = true;
}

void
FEMechanicsExplicitIntegrator::getFromInput(const Pointer<Database>& db, bool /*is_from_restart*/)
{
    // Problem parameters.
    if (db->isDouble("mass_density"))
        std::fill(d_rhos.begin(), d_rhos.end(), db->getDouble("mass_density"));
    else if (db->keyExists("mass_density"))
    {
        TBOX_ASSERT(db->getArraySize("mass_density") == d_rhos.size());
        db->getDoubleArray("mass_density", d_rhos.data(), db->getArraySize("mass_density"));
    }

    // libMesh settings.
    if (db->isString("libmesh_restart_file_extension"))
        d_libmesh_restart_file_extension = db->getString("libmesh_restart_file_extension");
    if (db->isString("libmesh_partitioner_type"))
        d_libmesh_partitioner_type = string_to_enum<LibmeshPartitionerType>(db->getString("libmesh_partitioner_type"));

    // Other settings.
    if (db->keyExists("do_log"))
        d_do_log = db->getBool("do_log");
    else if (db->keyExists("enable_logging"))
        d_do_log = db->getBool("enable_logging");
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
    int ver = db->getInteger("EXPLICIT_FE_MECHANICS_INTEGRATOR_VERSION");
    if (ver != EXPLICIT_FE_MECHANICS_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
