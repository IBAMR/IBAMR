// Filename: ExplicitFEMechanicsSolver.C
// Created on 12 Mar 2012 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#include "ExplicitFEMechanicsSolver.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/libmesh_utilities.h>

// LIBMESH INCLUDES
#include <boundary_info.h>
#include <dense_matrix.h>
#include <dense_vector.h>
#include <dof_map.h>
#include <fe_interface.h>
#include <mesh.h>
#include <petsc_linear_solver.h>
#include <petsc_matrix.h>
#include <petsc_vector.h>
#include <string_to_enum.h>
using namespace libMesh;

// SAMRAI INCLUDES
#include <tbox/RestartManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of ExplicitFEMechanicsSolver restart file data.
static const int EXPLICIT_FE_MECHANICS_SOLVER_VERSION = 1;
}

const std::string ExplicitFEMechanicsSolver::       COORDS_SYSTEM_NAME = "IB coordinates system";
const std::string ExplicitFEMechanicsSolver::COORD_MAPPING_SYSTEM_NAME = "IB coordinate mapping system";
const std::string ExplicitFEMechanicsSolver::        FORCE_SYSTEM_NAME = "IB force system";
const std::string ExplicitFEMechanicsSolver::     VELOCITY_SYSTEM_NAME = "IB velocity system";
const std::string ExplicitFEMechanicsSolver::    F_DIL_BAR_SYSTEM_NAME = "IB F_dil_bar system";

const short int ExplicitFEMechanicsSolver::DIRICHLET_BDRY_ID;

/////////////////////////////// PUBLIC ///////////////////////////////////////

ExplicitFEMechanicsSolver::ExplicitFEMechanicsSolver(
    const std::string& object_name,
    Pointer<Database> input_db,
    Mesh* mesh,
    bool register_for_restart)
    : d_num_parts(1)
{
    commonConstructor(object_name, input_db, std::vector<Mesh*>(1,mesh), register_for_restart);
    return;
}// ExplicitFEMechanicsSolver

ExplicitFEMechanicsSolver::ExplicitFEMechanicsSolver(
    const std::string& object_name,
    Pointer<Database> input_db,
    const std::vector<Mesh*>& meshes,
    bool register_for_restart)
    : d_num_parts(meshes.size())
{
    commonConstructor(object_name, input_db, meshes, register_for_restart);
    return;
}// ExplicitFEMechanicsSolver

ExplicitFEMechanicsSolver::~ExplicitFEMechanicsSolver()
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        delete d_equation_systems[part];
        std::map<std::string,libMesh::LinearSolver<double>*>& L2_proj_solver = d_L2_proj_solver[part];
        for (std::map<std::string,LinearSolver<double>*>::iterator it = L2_proj_solver.begin();
             it != L2_proj_solver.end(); ++it)
        {
            delete it->second;
        }
        std::map<std::string,libMesh::SparseMatrix<double>*>& L2_proj_matrix = d_L2_proj_matrix[part];
        for (std::map<std::string,SparseMatrix<double>*>::iterator it = L2_proj_matrix.begin();
             it != L2_proj_matrix.end(); ++it)
        {
            delete it->second;
        }
        std::map<std::string,libMesh::NumericVector<double>*>& L2_proj_matrix_diag = d_L2_proj_matrix_diag[part];
        for (std::map<std::string,NumericVector<double>*>::iterator it = L2_proj_matrix_diag.begin();
             it != L2_proj_matrix_diag.end(); ++it)
        {
            delete it->second;
        }
    }
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }
    return;
}// ~ExplicitFEMechanicsSolver

EquationSystems*
ExplicitFEMechanicsSolver::getEquationSystems(
    const unsigned int part) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(part < d_num_parts);
#endif
    return d_equation_systems[part];
}// getEquationSystems

void
ExplicitFEMechanicsSolver::registerInitialCoordinateMappingFunction(
    CoordinateMappingFcnPtr coordinate_mapping_fcn,
    void* coordinate_mapping_fcn_ctx,
    const unsigned int part)
{
    d_coordinate_mapping_fcns    [part] = coordinate_mapping_fcn;
    d_coordinate_mapping_fcn_ctxs[part] = coordinate_mapping_fcn_ctx;
    return;
}// registerInitialCoordinateMappingFunction

void
ExplicitFEMechanicsSolver::registerPK1StressTensorFunction(
    PK1StressFcnPtr PK1_stress_fcn,
    std::vector<unsigned int> PK1_stress_fcn_systems,
    void* PK1_stress_fcn_ctx,
    const unsigned int part)
{
    d_PK1_stress_fcns       [part] = PK1_stress_fcn;
    d_PK1_stress_fcn_systems[part] = PK1_stress_fcn_systems;
    d_PK1_stress_fcn_ctxs   [part] = PK1_stress_fcn_ctx;
    return;
}// registerPK1StressTensorFunction

void
ExplicitFEMechanicsSolver::registerLagBodyForceFunction(
    LagBodyForceFcnPtr lag_body_force_fcn,
    std::vector<unsigned int> lag_body_force_fcn_systems,
    void* lag_body_force_fcn_ctx,
    const unsigned int part)
{
    d_lag_body_force_fcns       [part] = lag_body_force_fcn;
    d_lag_body_force_fcn_systems[part] = lag_body_force_fcn_systems;
    d_lag_body_force_fcn_ctxs   [part] = lag_body_force_fcn_ctx;
    return;
}// registerLagBodyForceFunction

void
ExplicitFEMechanicsSolver::registerLagPressureFunction(
    LagPressureFcnPtr lag_pressure_fcn,
    std::vector<unsigned int> lag_pressure_fcn_systems,
    void* lag_pressure_fcn_ctx,
    const unsigned int part)
{
    d_lag_pressure_fcns       [part] = lag_pressure_fcn;
    d_lag_pressure_fcn_systems[part] = lag_pressure_fcn_systems;
    d_lag_pressure_fcn_ctxs   [part] = lag_pressure_fcn_ctx;
    return;
}// registerLagPressureFunction

void
ExplicitFEMechanicsSolver::registerLagSurfaceForceFunction(
    LagSurfaceForceFcnPtr lag_surface_force_fcn,
    std::vector<unsigned int> lag_surface_force_fcn_systems,
    void* lag_surface_force_fcn_ctx,
    const unsigned int part)
{
    d_lag_surface_force_fcns       [part] = lag_surface_force_fcn;
    d_lag_surface_force_fcn_systems[part] = lag_surface_force_fcn_systems;
    d_lag_surface_force_fcn_ctxs   [part] = lag_surface_force_fcn_ctx;
    return;
}// registerLagSurfaceForceFunction

void
ExplicitFEMechanicsSolver::preprocessIntegrateData(
    double current_time,
    double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    d_half_time = current_time+0.5*(new_time-current_time);

    // Ensure there is enough space for the solver components.
    d_L2_proj_solver     .resize(d_num_parts);
    d_L2_proj_matrix     .resize(d_num_parts);
    d_L2_proj_matrix_diag.resize(d_num_parts);
    d_L2_proj_quad_type  .resize(d_num_parts);
    d_L2_proj_quad_order .resize(d_num_parts);

    // Extract the FE data.
    d_X_systems          .resize(d_num_parts);
    d_X_current_vecs     .resize(d_num_parts);
    d_X_new_vecs         .resize(d_num_parts);
    d_X_half_vecs        .resize(d_num_parts);
    d_U_systems          .resize(d_num_parts);
    d_U_current_vecs     .resize(d_num_parts);
    d_U_new_vecs         .resize(d_num_parts);
    d_U_half_vecs        .resize(d_num_parts);
    d_F_systems          .resize(d_num_parts);
    d_F_half_vecs        .resize(d_num_parts);
    d_F_dil_bar_systems  .resize(d_num_parts);
    d_F_dil_bar_half_vecs.resize(d_num_parts);
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_X_systems          [part] = &d_equation_systems[part]->get_system(  COORDS_SYSTEM_NAME);
        d_X_current_vecs     [part] = dynamic_cast<PetscVector<double>*>(d_X_systems     [part]->solution.get());
        d_X_new_vecs         [part] = dynamic_cast<PetscVector<double>*>(d_X_current_vecs[part]->clone().release());  // WARNING: must be manually deleted
        d_X_half_vecs        [part] = dynamic_cast<PetscVector<double>*>(d_X_systems     [part]->current_local_solution.get());
        d_U_systems          [part] = &d_equation_systems[part]->get_system(VELOCITY_SYSTEM_NAME);
        d_U_current_vecs     [part] = dynamic_cast<PetscVector<double>*>(d_U_systems     [part]->solution.get());
        d_U_new_vecs         [part] = dynamic_cast<PetscVector<double>*>(d_U_current_vecs[part]->clone().release());  // WARNING: must be manually deleted
        d_U_half_vecs        [part] = dynamic_cast<PetscVector<double>*>(d_U_systems     [part]->current_local_solution.get());
        d_F_systems          [part] = &d_equation_systems[part]->get_system(   FORCE_SYSTEM_NAME);
        d_F_half_vecs        [part] = dynamic_cast<PetscVector<double>*>(d_F_systems     [part]->solution.get());
        d_F_dil_bar_systems  [part] = NULL;
        d_F_dil_bar_half_vecs[part] = NULL;
        if (d_use_Fbar_projection)
        {
            d_F_dil_bar_systems  [part] = &d_equation_systems[part]->get_system(F_DIL_BAR_SYSTEM_NAME);
            d_F_dil_bar_half_vecs[part] = dynamic_cast<PetscVector<double>*>(d_F_dil_bar_systems[part]->current_local_solution.get());
        }

        // Initialize X^{n+1/2} and X^{n+1} to equal X^{n}, and initialize
        // U^{n+1/2} and U^{n+1} to equal U^{n}.
        d_X_current_vecs[part]->localize(*d_X_half_vecs[part]);
        d_X_half_vecs[part]->close();
        d_X_current_vecs[part]->localize(*d_X_new_vecs[part]);
        d_X_new_vecs[part]->close();
        d_U_current_vecs[part]->localize(*d_U_half_vecs[part]);
        d_U_half_vecs[part]->close();
        d_U_current_vecs[part]->localize(*d_U_new_vecs[part]);
        d_U_new_vecs[part]->close();
    }
    return;
}// preprocessIntegrateData

void
ExplicitFEMechanicsSolver::integrateData(
    const double current_time,
    const double new_time)
{
    int ierr;
    const double half_time = current_time+0.5*(new_time-current_time);
    const double dt = new_time - current_time;

    // Advance X to time t^{n+1/2}.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        ierr = VecWAXPY(d_X_half_vecs[part]->vec(), 0.5*dt, d_U_current_vecs[part]->vec(), d_X_current_vecs[part]->vec()); IBTK_CHKERRQ(ierr);
    }

    // Compute F at time t^{n+1/2}.
    for (unsigned part = 0; part < d_num_parts; ++part)
    {
        if (d_use_Fbar_projection)
        {
            computeProjectedDilatationalStrain(*d_F_dil_bar_half_vecs[part], *d_X_half_vecs[part], part);
        }
        computeInteriorForceDensity(*d_F_half_vecs[part], *d_X_half_vecs[part], d_F_dil_bar_half_vecs[part], half_time, part);
    }

    // Advance U to time t^{n+1} and compute U at time t^{n+1/2}.
    for (unsigned part = 0; part < d_num_parts; ++part)
    {
        ierr = VecWAXPY(d_U_new_vecs[part]->vec(), dt/d_rho0, d_F_half_vecs[part]->vec(), d_U_current_vecs[part]->vec()); IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(d_U_half_vecs[part]->vec(), 0.5, 0.5, 0.0, d_U_current_vecs[part]->vec(), d_U_new_vecs[part]->vec());  IBTK_CHKERRQ(ierr);
    }

    // Advance X to time t^{n+1}.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        ierr = VecWAXPY(d_X_new_vecs[part]->vec(), dt, d_U_half_vecs[part]->vec(), d_X_current_vecs[part]->vec()); IBTK_CHKERRQ(ierr);
    }
    return;
}// integrateData

void
ExplicitFEMechanicsSolver::postprocessIntegrateData(
    double /*current_time*/,
    double /*new_time*/)
{
    for (unsigned part = 0; part < d_num_parts; ++part)
    {
        // Reset time-dependent Lagrangian data.
        (*d_X_current_vecs[part]) = (*d_X_new_vecs[part]);
        (*d_U_current_vecs[part]) = (*d_U_new_vecs[part]);
        if (d_use_Fbar_projection)
        {
            (*d_F_dil_bar_systems[part]->solution) = (*d_F_dil_bar_systems[part]->current_local_solution);
        }

        d_X_systems[part]->solution->localize(*d_X_systems[part]->current_local_solution);
        d_U_systems[part]->solution->localize(*d_U_systems[part]->current_local_solution);
        d_F_systems[part]->solution->localize(*d_F_systems[part]->current_local_solution);
        if (d_use_Fbar_projection)
        {
            d_F_dil_bar_systems[part]->solution->localize(*d_F_dil_bar_systems[part]->current_local_solution);
        }

        // Update the coordinate mapping dX = X - s.
        updateCoordinateMapping(part);

        // Deallocate Lagrangian scratch data.
        delete d_X_new_vecs[part];
        delete d_U_new_vecs[part];
    }
    d_X_systems          .clear();
    d_X_current_vecs     .clear();
    d_X_new_vecs         .clear();
    d_X_half_vecs        .clear();
    d_U_systems          .clear();
    d_U_current_vecs     .clear();
    d_U_new_vecs         .clear();
    d_U_half_vecs        .clear();
    d_F_systems          .clear();
    d_F_half_vecs        .clear();
    d_F_dil_bar_systems  .clear();
    d_F_dil_bar_half_vecs.clear();

    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time     = std::numeric_limits<double>::quiet_NaN();
    d_half_time    = std::numeric_limits<double>::quiet_NaN();
    return;
}// postprocessIntegrateData

void
ExplicitFEMechanicsSolver::initializeFEData()
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        // Initialize FE equation systems.
        EquationSystems* equation_systems = d_equation_systems[part];
        equation_systems->init();
        initializeCoordinates(part);
        updateCoordinateMapping(part);

        // Assemble systems.
        System& X_system = equation_systems->get_system<System>(COORDS_SYSTEM_NAME);
        X_system.assemble_before_solve = false;
        X_system.assemble();

        System& U_system = equation_systems->get_system<System>(VELOCITY_SYSTEM_NAME);
        U_system.assemble_before_solve = false;
        U_system.assemble();

        System& F_system = equation_systems->get_system<System>(FORCE_SYSTEM_NAME);
        F_system.assemble_before_solve = false;
        F_system.assemble();

        System& X_mapping_system = equation_systems->get_system<System>(COORD_MAPPING_SYSTEM_NAME);
        X_mapping_system.assemble_before_solve = false;
        X_mapping_system.assemble();

        // Set up boundary conditions.  Specifically, add appropriate boundary
        // IDs to the BoundaryInfo object associated with the mesh, and add DOF
        // constraints for the nodal forces and velocities.
        const MeshBase& mesh = equation_systems->get_mesh();
        const unsigned int dim = mesh.mesh_dimension();
        DofMap& F_dof_map = F_system.get_dof_map();
        DofMap& U_dof_map = U_system.get_dof_map();
        const unsigned int F_sys_num = F_system.number();
        const unsigned int U_sys_num = U_system.number();
        MeshBase::const_element_iterator       el_it  = mesh.elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.elements_end();
        for ( ; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                const bool at_mesh_bdry = elem->neighbor(side) == NULL;
                if (!at_mesh_bdry) continue;

                const std::vector<short int>& bdry_ids = mesh.boundary_info->boundary_ids(elem, side);
                const bool at_dirichlet_bdry = std::find(bdry_ids.begin(), bdry_ids.end(), DIRICHLET_BDRY_ID) != bdry_ids.end();
                if (!at_dirichlet_bdry) continue;

                for (unsigned int n = 0; n < elem->n_nodes(); ++n)
                {
                    if (!elem->is_node_on_side(n, side)) continue;

                    Node* node = elem->get_node(n);
                    mesh.boundary_info->add_node(node, DIRICHLET_BDRY_ID);
                    if (node->n_dofs(F_sys_num) > 0)
                    {
                        for (unsigned int d = 0; d < dim; ++d)
                        {
                            const int F_dof_index = node->dof_number(F_sys_num,d,0);
                            DofConstraintRow F_constraint_row;
                            F_constraint_row[F_dof_index] = 1.0;
                            F_dof_map.add_constraint_row(F_dof_index, F_constraint_row, 0.0, false);
                        }
                    }
                    if (node->n_dofs(U_sys_num) > 0)
                    {
                        for (unsigned int d = 0; d < dim; ++d)
                        {
                            const int U_dof_index = node->dof_number(U_sys_num,d,0);
                            DofConstraintRow U_constraint_row;
                            U_constraint_row[U_dof_index] = 1.0;
                            U_dof_map.add_constraint_row(U_dof_index, U_constraint_row, 0.0, false);
                        }
                    }
                }
            }
        }
    }

    d_is_initialized = true;
    return;
}// initializeFEData

void
ExplicitFEMechanicsSolver::putToDatabase(
    Pointer<Database> db)
{
    db->putInteger("EXPLICIT_FE_MECHANICS_SOLVER_VERSION", EXPLICIT_FE_MECHANICS_SOLVER_VERSION);
    db->putDouble("d_rho0", d_rho0);
    db->putBool("d_use_consistent_mass_matrix", d_use_consistent_mass_matrix);
    db->putBool("d_use_Fbar_projection", d_use_Fbar_projection);
    db->putString("d_fe_family", Utility::enum_to_string<FEFamily>(d_fe_family));
    db->putString("d_fe_order", Utility::enum_to_string<Order>(d_fe_order));
    db->putString("d_F_dil_bar_fe_family", Utility::enum_to_string<FEFamily>(d_F_dil_bar_fe_family));
    db->putString("d_F_dil_bar_fe_order", Utility::enum_to_string<Order>(d_F_dil_bar_fe_order));
    db->putString("d_quad_type", Utility::enum_to_string<QuadratureType>(d_quad_type));
    db->putString("d_quad_order", Utility::enum_to_string<Order>(d_quad_order));
    return;
}// putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

void
ExplicitFEMechanicsSolver::computeProjectedDilatationalStrain(
    NumericVector<double>& F_dil_bar_vec,
    NumericVector<double>& X_vec,
    const unsigned int part)
{
    // Extract the mesh.
    EquationSystems* equation_systems = d_equation_systems[part];
    const MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    AutoPtr<QBase> qrule = QBase::build(d_quad_type, dim, d_quad_order);

    // Extract the FE systems and DOF maps, and setup the FE objects.
    System& system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const DofMap& dof_map = system.get_dof_map();
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned int d = 0; d < dim; ++d) TBOX_ASSERT(dof_map.variable_type(d) == dof_map.variable_type(0));
#endif
    blitz::Array<std::vector<unsigned int>,1> dof_indices(dim);
    for (unsigned int d = 0; d < dim; ++d) dof_indices(d).reserve(27);
    AutoPtr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
    fe->attach_quadrature_rule(qrule.get());
    const std::vector<double>& JxW = fe->get_JxW();
    const std::vector<std::vector<VectorValue<double> > >& dphi = fe->get_dphi();

    System& F_dil_bar_system = equation_systems->get_system(F_DIL_BAR_SYSTEM_NAME);
    const DofMap& F_dil_bar_dof_map = F_dil_bar_system.get_dof_map();
    std::vector<unsigned int> F_dil_bar_dof_indices;
    F_dil_bar_dof_indices.reserve(27);
    AutoPtr<FEBase> F_dil_bar_fe(FEBase::build(dim, F_dil_bar_dof_map.variable_type(0)));
    F_dil_bar_fe->attach_quadrature_rule(qrule.get());
    const std::vector<std::vector<double> >& F_dil_bar_phi = F_dil_bar_fe->get_phi();

    // Setup global and elemental right-hand-side vectors.
    AutoPtr<NumericVector<double> > F_dil_bar_rhs_vec = F_dil_bar_vec.zero_clone();
    DenseVector<double> F_dil_bar_rhs_e;

    // F_dil_bar is the projection of F_dil = J^(1/d) II onto a (lower-order) FE
    // space.
    TensorValue<double> FF;
    blitz::Array<double,2> X_node;
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end   = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;

        fe->reinit(elem);
        for (unsigned int d = 0; d < dim; ++d)
        {
            dof_map.dof_indices(elem, dof_indices(d), d);
        }

        F_dil_bar_fe->reinit(elem);
        F_dil_bar_dof_map.dof_indices(elem, F_dil_bar_dof_indices);
        if (F_dil_bar_rhs_e.size() != F_dil_bar_dof_indices.size())
        {
            F_dil_bar_rhs_e.resize(F_dil_bar_dof_indices.size());  // NOTE: DenseVector::resize() automatically zeroes the vector contents.
        }
        else
        {
            F_dil_bar_rhs_e.zero();
        }

        const unsigned int n_qp = qrule->n_points();
        const unsigned int n_basis = F_dil_bar_dof_indices.size();

        get_values_for_interpolation(X_node, X_vec, dof_indices);
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            jacobian(FF,qp,X_node,dphi);
            const double J = FF.det();
            for (unsigned int k = 0; k < n_basis; ++k)
            {
                F_dil_bar_rhs_e(k) += F_dil_bar_phi[k][qp]*JxW[qp]*pow(J,1.0/static_cast<double>(dim));
            }
        }

        // Apply constraints (e.g., enforce periodic boundary conditions) and
        // add the elemental contributions to the global vector.
        F_dil_bar_dof_map.constrain_element_vector(F_dil_bar_rhs_e, F_dil_bar_dof_indices);
        F_dil_bar_rhs_vec->add_vector(F_dil_bar_rhs_e, F_dil_bar_dof_indices);
    }

    // Solve for F_dil_bar.
    F_dil_bar_rhs_vec->close();
    static const bool use_consistent_mass_matrix = true;
    computeL2Projection(F_dil_bar_vec, *F_dil_bar_rhs_vec, F_DIL_BAR_SYSTEM_NAME, part, use_consistent_mass_matrix);
    return;
}// computeProjectedDilatationalStrain

void
ExplicitFEMechanicsSolver::computeInteriorForceDensity(
    NumericVector<double>& G_vec,
    NumericVector<double>& X_vec,
    NumericVector<double>* F_dil_bar_vec,
    const double time,
    const unsigned int part)
{
    // Extract the mesh.
    EquationSystems* equation_systems = d_equation_systems[part];
    const MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    AutoPtr<QBase> qrule = QBase::build(d_quad_type, dim, d_quad_order);
    AutoPtr<QBase> qrule_face = QBase::build(d_quad_type, dim-1, d_quad_order);

    // Extract the FE systems and DOF maps, and setup the FE objects.
    System& system = equation_systems->get_system(FORCE_SYSTEM_NAME);
    const DofMap& dof_map = system.get_dof_map();
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned int d = 0; d < dim; ++d) TBOX_ASSERT(dof_map.variable_type(d) == dof_map.variable_type(0));
#endif
    blitz::Array<std::vector<unsigned int>,1> dof_indices(dim);
    for (unsigned int d = 0; d < dim; ++d) dof_indices(d).reserve(27);
    AutoPtr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
    fe->attach_quadrature_rule(qrule.get());
    const std::vector<Point>& q_point = fe->get_xyz();
    const std::vector<double>& JxW = fe->get_JxW();
    const std::vector<std::vector<double> >& phi = fe->get_phi();
    const std::vector<std::vector<VectorValue<double> > >& dphi = fe->get_dphi();
    AutoPtr<FEBase> fe_face(FEBase::build(dim, dof_map.variable_type(0)));
    fe_face->attach_quadrature_rule(qrule_face.get());
    const std::vector<Point>& q_point_face = fe_face->get_xyz();
    const std::vector<Point>& normal_face = fe_face->get_normals();
    const std::vector<std::vector<double> >& phi_face = fe_face->get_phi();
    const std::vector<std::vector<VectorValue<double> > >& dphi_face = fe_face->get_dphi();

#ifdef DEBUG_CHECK_ASSERTIONS
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    for (unsigned int d = 0; d < dim; ++d) TBOX_ASSERT(dof_map.variable_type(d) == X_dof_map.variable_type(d));
#endif

    System* F_dil_bar_system = NULL;
    const DofMap* F_dil_bar_dof_map = NULL;
    std::vector<unsigned int> F_dil_bar_dof_indices;
    F_dil_bar_dof_indices.reserve(27);
    AutoPtr<FEBase> F_dil_bar_fe;
    const std::vector<std::vector<double> >* F_dil_bar_phi = NULL;
    AutoPtr<FEBase> F_dil_bar_fe_face;
    const std::vector<std::vector<double> >* F_dil_bar_phi_face = NULL;
    if (F_dil_bar_vec != NULL)
    {
        F_dil_bar_system = &equation_systems->get_system(F_DIL_BAR_SYSTEM_NAME);
        F_dil_bar_dof_map = &F_dil_bar_system->get_dof_map();
        F_dil_bar_fe = FEBase::build(dim, F_dil_bar_dof_map->variable_type(0));
        F_dil_bar_fe->attach_quadrature_rule(qrule.get());
        F_dil_bar_phi = &F_dil_bar_fe->get_phi();
        F_dil_bar_fe_face = FEBase::build(dim, F_dil_bar_dof_map->variable_type(0));
        F_dil_bar_fe_face->attach_quadrature_rule(qrule_face.get());
        F_dil_bar_phi_face = &F_dil_bar_fe_face->get_phi();
    }

    // Setup extra data needed to compute stresses/forces.
    std::vector<NumericVector<double>*> PK1_stress_fcn_data;
    for (std::vector<unsigned int>::const_iterator cit = d_PK1_stress_fcn_systems[part].begin();
         cit != d_PK1_stress_fcn_systems[part].end(); ++cit)
    {
        System& system = equation_systems->get_system(*cit);
        system.update();
        PK1_stress_fcn_data.push_back(system.current_local_solution.get());
    }

    std::vector<NumericVector<double>*> lag_body_force_fcn_data;
    for (std::vector<unsigned int>::const_iterator cit = d_lag_body_force_fcn_systems[part].begin();
         cit != d_lag_body_force_fcn_systems[part].end(); ++cit)
    {
        System& system = equation_systems->get_system(*cit);
        system.update();
        lag_body_force_fcn_data.push_back(system.current_local_solution.get());
    }

    std::vector<NumericVector<double>*> lag_pressure_fcn_data;
    for (std::vector<unsigned int>::const_iterator cit = d_lag_pressure_fcn_systems[part].begin();
         cit != d_lag_pressure_fcn_systems[part].end(); ++cit)
    {
        System& system = equation_systems->get_system(*cit);
        system.update();
        lag_pressure_fcn_data.push_back(system.current_local_solution.get());
    }

    std::vector<NumericVector<double>*> lag_surface_force_fcn_data;
    for (std::vector<unsigned int>::const_iterator cit = d_lag_surface_force_fcn_systems[part].begin();
         cit != d_lag_surface_force_fcn_systems[part].end(); ++cit)
    {
        System& system = equation_systems->get_system(*cit);
        system.update();
        lag_surface_force_fcn_data.push_back(system.current_local_solution.get());
    }

    // Setup global and elemental right-hand-side vectors.
    AutoPtr<NumericVector<double> > G_rhs_vec = G_vec.zero_clone();
    DenseVector<double> G_rhs_e[LIBMESH_DIM];

    // Loop over the elements to compute the right-hand side vector.  This is
    // computed via
    //
    //    rhs_k = -int{PP(s,t) grad phi_k(s)}ds + int{PP(s,t) N(s,t) phi_k(s)}dA(s)
    //
    // This right-hand side vector is used to solve for the nodal values of the
    // interior elastic force density.
    TensorValue<double> PP, FF, FF_inv_trans, FF_bar;
    VectorValue<double> F, F_b, F_s, F_qp, n;
    Point X_qp;
    double P;
    blitz::Array<double,2> X_node;
    blitz::Array<double,1> F_dil_bar_node;
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end   = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;

        fe->reinit(elem);
        for (unsigned int d = 0; d < dim; ++d)
        {
            dof_map.dof_indices(elem, dof_indices(d), d);
            if (G_rhs_e[d].size() != dof_indices(d).size())
            {
                G_rhs_e[d].resize(dof_indices(d).size());  // NOTE: DenseVector::resize() automatically zeroes the vector contents.
            }
            else
            {
                G_rhs_e[d].zero();
            }
        }

        if (F_dil_bar_vec != NULL)
        {
            F_dil_bar_fe->reinit(elem);
            F_dil_bar_dof_map->dof_indices(elem, F_dil_bar_dof_indices);
        }

        const unsigned int n_qp = qrule->n_points();
        const unsigned int n_basis = dof_indices(0).size();

        get_values_for_interpolation(X_node, X_vec, dof_indices);
        if (F_dil_bar_vec != NULL) get_values_for_interpolation(F_dil_bar_node, *F_dil_bar_vec, F_dil_bar_dof_indices);
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            const Point& s_qp = q_point[qp];
            interpolate(X_qp,qp,X_node,phi);
            jacobian(FF,qp,X_node,dphi);
            if (F_dil_bar_vec != NULL)
            {
                jacobian(FF_bar,qp,X_node,dphi,F_dil_bar_node,*F_dil_bar_phi);
            }
            else
            {
                FF_bar = FF;
            }

            if (d_PK1_stress_fcns[part] != NULL)
            {
                // Compute the value of the first Piola-Kirchhoff stress tensor
                // at the quadrature point and add the corresponding forces to
                // the right-hand-side vector.
                d_PK1_stress_fcns[part](PP,FF_bar,X_qp,s_qp,elem,X_vec,PK1_stress_fcn_data,time,d_PK1_stress_fcn_ctxs[part]);
                for (unsigned int k = 0; k < n_basis; ++k)
                {
                    F_qp = -PP*dphi[k][qp]*JxW[qp];
                    for (unsigned int i = 0; i < dim; ++i)
                    {
                        G_rhs_e[i](k) += F_qp(i);
                    }
                }
            }

            if (d_lag_body_force_fcns[part] != NULL)
            {
                // Compute the value of the body force at the quadrature point
                // and add the corresponding forces to the right-hand-side
                // vector.
                d_lag_body_force_fcns[part](F_b,FF,X_qp,s_qp,elem,X_vec,lag_body_force_fcn_data,time,d_lag_body_force_fcn_ctxs[part]);
                for (unsigned int k = 0; k < n_basis; ++k)
                {
                    F_qp = phi[k][qp]*JxW[qp]*F_b;
                    for (unsigned int i = 0; i < dim; ++i)
                    {
                        G_rhs_e[i](k) += F_qp(i);
                    }
                }
            }
        }

        // Loop over the element boundaries.
        for (unsigned short int side = 0; side < elem->n_sides(); ++side)
        {
            // Determine whether we are at a physical boundary and, if so,
            // whether it is a Dirichlet boundary.
            bool at_physical_bdry = elem->neighbor(side) == NULL;
            bool at_dirichlet_bdry = false;
            const std::vector<short int>& bdry_ids = mesh.boundary_info->boundary_ids(elem, side);
            for (std::vector<short int>::const_iterator cit = bdry_ids.begin(); cit != bdry_ids.end(); ++cit)
            {
                const short int bdry_id = *cit;
                at_physical_bdry  = at_physical_bdry  && !dof_map.is_periodic_boundary(bdry_id);
                at_dirichlet_bdry = at_dirichlet_bdry || (bdry_id == DIRICHLET_BDRY_ID);
            }

            // Skip non-physical boundaries.
            if (!at_physical_bdry) continue;

            // Determine whether we need to compute surface forces along this
            // part of the physical boundary; if not, skip the present side.
            const bool compute_pressure           = d_lag_pressure_fcns     [part] != NULL && !at_dirichlet_bdry;
            const bool compute_surface_force      = d_lag_surface_force_fcns[part] != NULL && !at_dirichlet_bdry;
            if (!(compute_pressure || compute_surface_force)) continue;

            fe_face->reinit(elem, side);
            if (F_dil_bar_vec != NULL) F_dil_bar_fe_face->reinit(elem, side);

            const unsigned int n_qp = qrule_face->n_points();

            get_values_for_interpolation(X_node, X_vec, dof_indices);
            if (F_dil_bar_vec != NULL) get_values_for_interpolation(F_dil_bar_node, *F_dil_bar_vec, F_dil_bar_dof_indices);
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                const Point& s_qp = q_point_face[qp];
                interpolate(X_qp,qp,X_node,phi_face);
                jacobian(FF,qp,X_node,dphi_face);
                const double J = std::abs(FF.det());
                tensor_inverse_transpose(FF_inv_trans,FF,dim);
                F.zero();
                if (compute_pressure && d_lag_pressure_fcns[part] != NULL)
                {
                    // Compute the value of the pressure at the quadrature point
                    // and add the corresponding force to the right-hand-side
                    // vector.
                    d_lag_pressure_fcns[part](P,FF,X_qp,s_qp,elem,side,X_vec,lag_pressure_fcn_data,time,d_lag_pressure_fcn_ctxs[part]);
                    F -= P*J*FF_inv_trans*normal_face[qp];
                }
                if (compute_surface_force)
                {
                    // Compute the value of the surface force at the quadrature
                    // point and add the corresponding force to the
                    // right-hand-side vector.
                    d_lag_surface_force_fcns[part](F_s,FF,X_qp,s_qp,elem,side,X_vec,lag_surface_force_fcn_data,time,d_lag_surface_force_fcn_ctxs[part]);
                    F += F_s;
                }
            }
        }

        // Apply constraints (e.g., enforce periodic boundary conditions) and
        // add the elemental contributions to the global vector.
        for (unsigned int i = 0; i < dim; ++i)
        {
            dof_map.constrain_element_vector(G_rhs_e[i], dof_indices(i));
            G_rhs_vec->add_vector(G_rhs_e[i], dof_indices(i));
        }
    }

    // Solve for G.
    G_rhs_vec->close();
    computeL2Projection(G_vec, *G_rhs_vec, FORCE_SYSTEM_NAME, part, d_use_consistent_mass_matrix);
    return;
}// computeInteriorForceDensity

void
ExplicitFEMechanicsSolver::initializeCoordinates(
    const unsigned int part)
{
    EquationSystems* equation_systems = d_equation_systems[part];
    MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    const bool identity_mapping = d_coordinate_mapping_fcns[part] == NULL;
    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
    {
        Node* n = *it;
        if (n->n_vars(X_sys_num) > 0)
        {
            libmesh_assert(n->n_vars(X_sys_num) == dim);
            const Point& s = *n;
            Point X = s;
            if (!identity_mapping)
            {
                d_coordinate_mapping_fcns[part](X, s, d_coordinate_mapping_fcn_ctxs[part]);
            }
            for (unsigned int d = 0; d < dim; ++d)
            {
                const int dof_index = n->dof_number(X_sys_num,d,0);
                X_coords.set(dof_index,X(d));
            }
        }
    }
    X_system.get_dof_map().enforce_constraints_exactly(X_system, &X_coords);
    X_coords.close();
    return;
}// initializeCoordinates

void
ExplicitFEMechanicsSolver::updateCoordinateMapping(
    const unsigned int part)
{
    EquationSystems* equation_systems = d_equation_systems[part];
    MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    System& X_mapping_system = equation_systems->get_system(COORD_MAPPING_SYSTEM_NAME);
    const unsigned int X_mapping_sys_num = X_mapping_system.number();
    NumericVector<double>& dX_coords = *X_mapping_system.solution;
    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
    {
        Node* n = *it;
        if (n->n_vars(X_sys_num) > 0)
        {
            libmesh_assert(n->n_vars(X_sys_num) == dim);
            libmesh_assert(n->n_vars(X_mapping_sys_num) == dim);
            const Point& s = *n;
            for (unsigned int d = 0; d < dim; ++d)
            {
                const int X_dof_index = n->dof_number(X_sys_num,d,0);
                const int dX_dof_index = n->dof_number(X_mapping_sys_num,d,0);
                dX_coords.set(dX_dof_index,X_coords(X_dof_index)-s(d));
            }
        }
    }
    dX_coords.close();
    return;
}// updateCoordinateMapping

std::pair<LinearSolver<double>*,SparseMatrix<double>*>
ExplicitFEMechanicsSolver::buildL2ProjectionSolver(
    const std::string& system_name,
    const unsigned int part,
    const QuadratureType quad_type,
    const Order quad_order)
{
    std::map<std::string,libMesh::LinearSolver<double>*>& L2_proj_solver     = d_L2_proj_solver    [part];
    std::map<std::string,libMesh::SparseMatrix<double>*>& L2_proj_matrix     = d_L2_proj_matrix    [part];
    std::map<std::string,libMeshEnums::QuadratureType>&   L2_proj_quad_type  = d_L2_proj_quad_type [part];
    std::map<std::string,libMeshEnums::Order>&            L2_proj_quad_order = d_L2_proj_quad_order[part];
    if ((L2_proj_solver.count(system_name) == 0 || L2_proj_matrix.count(system_name) == 0) ||
        (L2_proj_quad_type[system_name] != quad_type) || (L2_proj_quad_order[system_name] != quad_order))
    {
        EquationSystems* equation_systems = d_equation_systems[part];
        const MeshBase& mesh = equation_systems->get_mesh();
        const unsigned int dim = mesh.mesh_dimension();
        AutoPtr<QBase> qrule = QBase::build(quad_type, dim, quad_order);

        System& system = equation_systems->get_system(system_name);
        const int sys_num = system.number();
        DofMap& dof_map = system.get_dof_map();
        // Sparsity patterns are not automatically computed for all system
        // types.  If one has not been computed for this system, we compute it
        // now.
        if (dof_map.get_n_nz().size() != dof_map.n_local_dofs())
        {
            dof_map.compute_sparsity(mesh);
        }
        std::vector<unsigned int> dof_indices;
        AutoPtr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
        fe->attach_quadrature_rule(qrule.get());
        const std::vector<double>& JxW = fe->get_JxW();
        const std::vector<std::vector<double> >& phi = fe->get_phi();

        LinearSolver<double>* solver = LinearSolver<double>::build().release();
        solver->init();

        SparseMatrix<double>* M_mat = SparseMatrix<double>::build().release();
        M_mat->attach_dof_map(dof_map);
        M_mat->init();

        DenseMatrix<double> M_e;

        // Loop over the mesh to construct the system matrix.
        const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end   = mesh.active_local_elements_end();
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            const Elem* const elem = *el_it;
            fe->reinit(elem);
            for (unsigned int var_num = 0; var_num < dof_map.n_variables(); ++var_num)
            {
                dof_map.dof_indices(elem, dof_indices, var_num);
                M_e.resize(dof_indices.size(), dof_indices.size());
                const unsigned int n_basis = dof_indices.size();
                const unsigned int n_qp = qrule->n_points();
                for (unsigned int i = 0; i < n_basis; ++i)
                {
                    for (unsigned int j = 0; j < n_basis; ++j)
                    {
                        for (unsigned int qp = 0; qp < n_qp; ++qp)
                        {
                            M_e(i,j) += (phi[i][qp]*phi[j][qp])*JxW[qp];
                        }
                    }
                }
                dof_map.constrain_element_matrix(M_e, dof_indices);
                M_mat->add_matrix(M_e, dof_indices);
            }
        }

        // Flush assemble the matrix.
        M_mat->close();

        // Reset values at Dirichlet boundaries.
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                if (elem->neighbor(side) != NULL) continue;
                const std::vector<short int>& bdry_ids = mesh.boundary_info->boundary_ids(elem, side);
                const bool at_dirichlet_bdry = std::find(bdry_ids.begin(), bdry_ids.end(), DIRICHLET_BDRY_ID) != bdry_ids.end();
                if (!at_dirichlet_bdry) continue;
                fe->reinit(elem);
                for (unsigned int n = 0; n < elem->n_nodes(); ++n)
                {
                    if (elem->is_node_on_side(n, side))
                    {
                        Node* node = elem->get_node(n);
                        for (unsigned int var_num = 0; var_num < dof_map.n_variables(); ++var_num)
                        {
                            dof_map.dof_indices(elem, dof_indices, var_num);
                            const unsigned int n_comp = node->n_comp(sys_num, var_num);
                            for (unsigned int comp = 0; comp < n_comp; ++comp)
                            {
                                const unsigned int node_dof_index = node->dof_number(sys_num, var_num, comp);
                                if (dof_map.is_constrained_dof(node_dof_index))
                                {
                                    for (std::vector<unsigned int>::const_iterator cit = dof_indices.begin(); cit != dof_indices.end(); ++cit)
                                    {
                                        const unsigned int k = *cit;
                                        M_mat->set(node_dof_index, k, (node_dof_index == k ? 1.0 : 0.0));
                                        M_mat->set(k, node_dof_index, (node_dof_index == k ? 1.0 : 0.0));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Assemble the matrix.
        M_mat->close();

        // Setup the solver.
        solver->reuse_preconditioner(true);

        // Store the solver, mass matrix, and configuration options.
        L2_proj_solver[system_name] = solver;
        L2_proj_matrix[system_name] = M_mat;
        L2_proj_quad_type[system_name] = quad_type;
        L2_proj_quad_order[system_name] = quad_order;
    }
    return std::make_pair(L2_proj_solver[system_name], L2_proj_matrix[system_name]);
}// buildL2ProjectionSolver

NumericVector<double>*
ExplicitFEMechanicsSolver::buildDiagonalL2MassMatrix(
    const std::string& system_name,
    const unsigned int part)
{
    std::map<std::string,libMesh::NumericVector<double>*>& L2_proj_matrix_diag = d_L2_proj_matrix_diag[part];
    if (L2_proj_matrix_diag.count(system_name) == 0)
    {
        EquationSystems* equation_systems = d_equation_systems[part];
        const MeshBase& mesh = equation_systems->get_mesh();
        const unsigned int dim = mesh.mesh_dimension();
        AutoPtr<QBase> qrule_trap    = QBase::build(QTRAP   , dim, FIRST);
        AutoPtr<QBase> qrule_simpson = QBase::build(QSIMPSON, dim, THIRD);

        System& system = equation_systems->get_system(system_name);
        const int sys_num = system.number();
        DofMap& dof_map = system.get_dof_map();
        std::vector<unsigned int> dof_indices;

        AutoPtr<FEBase> fe_trap(   FEBase::build(dim, dof_map.variable_type(0)));
        AutoPtr<FEBase> fe_simpson(FEBase::build(dim, dof_map.variable_type(0)));
        fe_trap   ->attach_quadrature_rule(qrule_trap   .get());
        fe_simpson->attach_quadrature_rule(qrule_simpson.get());

        NumericVector<double>* M_vec = system.solution->zero_clone().release();
        DenseVector<double> M_diag_e;

        // Loop over the mesh to construct the (diagonal) system matrix.
        //
        // We construct diagonal elemental mass matrices by using low-order
        // nodal quadrature rules.
        const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end   = mesh.active_local_elements_end();
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            const Elem* const elem = *el_it;
            QBase* qrule = NULL;
            FEBase* fe = NULL;
            if (elem->default_order() == FIRST)
            {
                qrule = qrule_trap.get();
                fe = fe_trap.get();
            }
            else
            {
                ElemType elem_type = elem->type();
                if (elem_type == EDGE3   ||
                    elem_type == TRI6    ||
                    elem_type == QUAD9   ||
                    elem_type == TET10   ||
                    elem_type == PRISM18 ||
                    elem_type == HEX27)
                {
                    qrule = qrule_simpson.get();
                    fe = fe_simpson.get();
                }
                else
                {
                    TBOX_ERROR("FEDataManager::buildDiagonalL2MassMatrix():\n"
                               << "  unsupported element type: " << Utility::enum_to_string<ElemType>(elem_type) << "\n");
                }
            }
            const std::vector<double>& JxW = fe->get_JxW();
            const std::vector<std::vector<double> >& phi = fe->get_phi();
            fe->reinit(elem);
            for (unsigned int var_num = 0; var_num < dof_map.n_variables(); ++var_num)
            {
                dof_map.dof_indices(elem, dof_indices, var_num);
                M_diag_e.resize(dof_indices.size());
                const unsigned int n_basis = dof_indices.size();
                const unsigned int n_qp = qrule->n_points();
                for (unsigned int i = 0; i < n_basis; ++i)
                {
                    for (unsigned int j = 0; j < n_basis; ++j)
                    {
                        for (unsigned int qp = 0; qp < n_qp; ++qp)
                        {
                            const double integrand = (phi[i][qp]*phi[j][qp])*JxW[qp];
                            if (i == j) M_diag_e(i) += integrand;
#ifdef DEBUG_CHECK_ASSERTIONS
                            else TBOX_ASSERT(std::abs(integrand) < std::numeric_limits<double>::epsilon());
#endif
                        }
                    }
                }
                dof_map.constrain_element_vector(M_diag_e, dof_indices);
                M_vec->add_vector(M_diag_e, dof_indices);
            }
        }

        // Flush assemble the vector representation of the diagonal matrix.
        M_vec->close();

        // Reset values at Dirichlet boundaries.
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                if (elem->neighbor(side) != NULL) continue;
                const std::vector<short int>& bdry_ids = mesh.boundary_info->boundary_ids(elem, side);
                const bool at_dirichlet_bdry = std::find(bdry_ids.begin(), bdry_ids.end(), DIRICHLET_BDRY_ID) != bdry_ids.end();
                if (!at_dirichlet_bdry) continue;
                for (unsigned int n = 0; n < elem->n_nodes(); ++n)
                {
                    if (elem->is_node_on_side(n, side))
                    {
                        Node* node = elem->get_node(n);
                        for (unsigned int var_num = 0; var_num < dof_map.n_variables(); ++var_num)
                        {
                            const unsigned int n_comp = node->n_comp(sys_num, var_num);
                            for (unsigned int comp = 0; comp < n_comp; ++comp)
                            {
                                const unsigned int node_dof_index = node->dof_number(sys_num, var_num, comp);
                                if (dof_map.is_constrained_dof(node_dof_index)) M_vec->set(node_dof_index, 1.0);
                            }
                        }
                    }
                }
            }
        }

        // Assemble the vector representation of the diagonal matrix.
        M_vec->close();

        // Store the diagonal mass matrix.
        L2_proj_matrix_diag[system_name] = M_vec;
    }
    return L2_proj_matrix_diag[system_name];
}// buildDiagonalL2MassMatrix

bool
ExplicitFEMechanicsSolver::computeL2Projection(
    NumericVector<double>& U_vec,
    NumericVector<double>& F_vec,
    const std::string& system_name,
    const unsigned int part,
    const bool consistent_mass_matrix,
    const QuadratureType quad_type,
    const Order quad_order,
    const double tol,
    const unsigned int max_its)
{
    int ierr;
    bool converged = false;

    F_vec.close();
    EquationSystems* equation_systems = d_equation_systems[part];
    System& system = equation_systems->get_system(system_name);
    const DofMap& dof_map = system.get_dof_map();
    dof_map.enforce_constraints_exactly(system, &F_vec);
    if (consistent_mass_matrix)
    {
        std::pair<libMesh::LinearSolver<double>*,SparseMatrix<double>*> proj_solver_components = buildL2ProjectionSolver(system_name, part, quad_type, quad_order);
        PetscLinearSolver<double>* solver = dynamic_cast<PetscLinearSolver<double>*>(proj_solver_components.first);
        PetscMatrix<double>* M_mat = dynamic_cast<PetscMatrix<double>*>(proj_solver_components.second);
        PetscBool rtol_set;
        double runtime_rtol;
        ierr = PetscOptionsGetReal("","-ksp_rtol",&runtime_rtol,&rtol_set); IBTK_CHKERRQ(ierr);
        PetscBool max_it_set;
        int runtime_max_it;
        ierr = PetscOptionsGetInt("","-ksp_max_it",&runtime_max_it,&max_it_set); IBTK_CHKERRQ(ierr);
        ierr = KSPSetFromOptions(solver->ksp()); IBTK_CHKERRQ(ierr);
        solver->solve(*M_mat, *M_mat, U_vec, F_vec, rtol_set ? runtime_rtol : tol, max_it_set ? runtime_max_it : max_its);
        KSPConvergedReason reason;
        ierr = KSPGetConvergedReason(solver->ksp(), &reason); IBTK_CHKERRQ(ierr);
        converged = reason > 0;
    }
    else
    {
        PetscVector<double>* M_diag_vec = dynamic_cast<PetscVector<double>*>(buildDiagonalL2MassMatrix(system_name));
        Vec M_diag_petsc_vec = M_diag_vec->vec();
        Vec U_petsc_vec = dynamic_cast<PetscVector<double>*>(&U_vec)->vec();
        Vec F_petsc_vec = dynamic_cast<PetscVector<double>*>(&F_vec)->vec();
        ierr = VecPointwiseDivide(U_petsc_vec, F_petsc_vec, M_diag_petsc_vec); IBTK_CHKERRQ(ierr);
        converged = true;
    }
    dof_map.enforce_constraints_exactly(system, &U_vec);
    return converged;
}// computeL2Projection

/////////////////////////////// PRIVATE //////////////////////////////////////

void
ExplicitFEMechanicsSolver::commonConstructor(
    const std::string& object_name,
    Pointer<Database> input_db,
    const std::vector<libMesh::Mesh*>& meshes,
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
    d_use_consistent_mass_matrix = true;
    d_use_Fbar_projection = false;
    d_fe_family = LAGRANGE;
    d_fe_order = INVALID_ORDER;
    d_F_dil_bar_fe_family = MONOMIAL;
    d_F_dil_bar_fe_order = CONSTANT;
    d_quad_type = QGAUSS;
    d_quad_order = FIFTH;
    d_do_log = false;

    // Initialize function pointers to NULL.
    d_coordinate_mapping_fcns.resize(d_num_parts,NULL);
    d_coordinate_mapping_fcn_ctxs.resize(d_num_parts,NULL);
    d_PK1_stress_fcns.resize(d_num_parts,NULL);
    d_PK1_stress_fcn_systems.resize(d_num_parts);
    d_PK1_stress_fcn_ctxs.resize(d_num_parts,NULL);
    d_lag_body_force_fcns.resize(d_num_parts,NULL);
    d_lag_body_force_fcn_systems.resize(d_num_parts);
    d_lag_body_force_fcn_ctxs.resize(d_num_parts,NULL);
    d_lag_pressure_fcns.resize(d_num_parts,NULL);
    d_lag_pressure_fcn_systems.resize(d_num_parts);
    d_lag_pressure_fcn_ctxs.resize(d_num_parts,NULL);
    d_lag_surface_force_fcns.resize(d_num_parts,NULL);
    d_lag_surface_force_fcn_systems.resize(d_num_parts);
    d_lag_surface_force_fcn_ctxs.resize(d_num_parts,NULL);

    // Determine whether we should use first-order or second-order shape
    // functions for each part of the structure.
    bool mesh_has_first_order_elems  = false;
    bool mesh_has_second_order_elems = false;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        const MeshBase& mesh = *meshes[part];
        MeshBase::const_element_iterator       el_it  = mesh.elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.elements_end();
        for ( ; el_it != el_end; ++el_it)
        {
            const Elem* const elem = *el_it;
            mesh_has_first_order_elems  = mesh_has_first_order_elems  || elem->default_order() == FIRST ;
            mesh_has_second_order_elems = mesh_has_second_order_elems || elem->default_order() == SECOND;
        }
    }
    mesh_has_first_order_elems  = SAMRAI_MPI::maxReduction(mesh_has_first_order_elems );
    mesh_has_second_order_elems = SAMRAI_MPI::maxReduction(mesh_has_second_order_elems);
    if (( mesh_has_first_order_elems &&  mesh_has_second_order_elems) ||
        (!mesh_has_first_order_elems && !mesh_has_second_order_elems))
    {
        TBOX_ERROR(d_object_name << "::ExplicitFEMechanicsSolver():\n"
                   << "  all parts of FE mesh must contain only FIRST order elements or only SECOND order elements" << std::endl);
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
    pout << d_object_name << ": using " << Utility::enum_to_string<Order>(d_fe_order) << " order " << Utility::enum_to_string<FEFamily>(d_fe_family) << " finite elements.\n\n";

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, from_restart);

    // Setup EquationSystems objects for each part and setup Systems.
    d_meshes = meshes;
    d_equation_systems.resize(d_num_parts, NULL);
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        // Create FE equation systems object.
        d_equation_systems[part] = new EquationSystems(*d_meshes[part]);
        EquationSystems* equation_systems = d_equation_systems[part];
        const MeshBase& mesh = equation_systems->get_mesh();
        const unsigned int dim = mesh.mesh_dimension();

        // Create FE systems and corresponding variables.
        System& X_system = equation_systems->add_system<System>(COORDS_SYSTEM_NAME);
        for (unsigned int d = 0; d < dim; ++d)
        {
            std::ostringstream os;
            os << "X_" << d;
            X_system.add_variable(os.str(), d_fe_order, d_fe_family);
        }

        System& X_mapping_system = equation_systems->add_system<System>(COORD_MAPPING_SYSTEM_NAME);
        for (unsigned int d = 0; d < dim; ++d)
        {
            std::ostringstream os;
            os << "dX_" << d;
            X_mapping_system.add_variable(os.str(), d_fe_order, d_fe_family);
        }

        System& F_system = equation_systems->add_system<System>(FORCE_SYSTEM_NAME);
        for (unsigned int d = 0; d < dim; ++d)
        {
            std::ostringstream os;
            os << "F_" << d;
            F_system.add_variable(os.str(), d_fe_order, d_fe_family);
        }

        System& U_system = equation_systems->add_system<System>(VELOCITY_SYSTEM_NAME);
        for (unsigned int d = 0; d < dim; ++d)
        {
            std::ostringstream os;
            os << "U_" << d;
            U_system.add_variable(os.str(), d_fe_order, d_fe_family);
        }

        if (d_use_Fbar_projection)
        {
            System& F_dil_bar_system = equation_systems->add_system<System>(F_DIL_BAR_SYSTEM_NAME);
            F_dil_bar_system.add_variable("F_dil_bar", d_F_dil_bar_fe_order, d_F_dil_bar_fe_family);
        }
    }

    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time     = std::numeric_limits<double>::quiet_NaN();
    d_half_time    = std::numeric_limits<double>::quiet_NaN();

    d_is_initialized = false;
    return;
}// commonConstructor

void
ExplicitFEMechanicsSolver::getFromInput(
    Pointer<Database> db,
    bool is_from_restart)
{
    if (!is_from_restart)
    {
        d_rho0 = db->getDouble("rho0");
        if (db->isBool("use_consistent_mass_matrix")) d_use_consistent_mass_matrix = db->getBool("use_consistent_mass_matrix");
        if (db->isBool("use_Fbar_projection")) d_use_Fbar_projection = db->getBool("use_Fbar_projection");
        if (db->isString("F_dil_bar_fe_family")) d_F_dil_bar_fe_family = Utility::string_to_enum<FEFamily>(db->getString("F_dil_bar_fe_family"));
        if (db->isString("F_dil_bar_fe_order")) d_F_dil_bar_fe_order = Utility::string_to_enum<Order>(db->getString("F_dil_bar_fe_order"));
        if (db->isString("quad_type")) d_quad_type = Utility::string_to_enum<QuadratureType>(db->getString("quad_type"));
        if (db->isString("quad_order")) d_quad_order = Utility::string_to_enum<Order>(db->getString("quad_order"));
    }
    if      (db->keyExists("do_log"        )) d_do_log = db->getBool("do_log"        );
    else if (db->keyExists("enable_logging")) d_do_log = db->getBool("enable_logging");
    return;
}// getFromInput

void
ExplicitFEMechanicsSolver::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to "
                   << d_object_name << " not found in restart file." << std::endl);
    }
    int ver = db->getInteger("EXPLICIT_FE_MECHANICS_SOLVER_VERSION");
    if (ver != EXPLICIT_FE_MECHANICS_SOLVER_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    d_rho0 = db->getDouble("d_rho0");
    d_use_consistent_mass_matrix = db->getBool("d_use_consistent_mass_matrix");
    d_use_Fbar_projection = db->getBool("d_use_Fbar_projection");
    d_fe_family = Utility::string_to_enum<FEFamily>(db->getString("d_fe_family"));
    d_fe_order = Utility::string_to_enum<Order>(db->getString("d_fe_order"));
    d_F_dil_bar_fe_family = Utility::string_to_enum<FEFamily>(db->getString("d_F_dil_bar_fe_family"));
    d_F_dil_bar_fe_order = Utility::string_to_enum<Order>(db->getString("d_F_dil_bar_fe_order"));
    d_quad_type = Utility::string_to_enum<QuadratureType>(db->getString("d_quad_type"));
    d_quad_order = Utility::string_to_enum<Order>(db->getString("d_quad_order"));
    return;
}// getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
