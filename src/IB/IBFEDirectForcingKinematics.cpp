// Filename: IBFEDirectForcingKinematics.cpp
// Created on 24 May 2018 by Amneet Bhalla
//
// Copyright (c) 2002-2018, Amneet Bhalla
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

#include "ibamr/IBFEDirectForcingKinematics.h"
#include "ibamr/IBFEMethod.h"
#include "ibamr/namespaces.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "libmesh/equation_systems.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
#if (NDIM == 3)
inline Eigen::Vector3d
solve_3x3_system(const Eigen::Vector3d& L, const Eigen::Matrix3d& T)
{
    const double a1 = T(0, 0);
    const double a2 = T(0, 1);
    const double a3 = T(0, 2);
    const double b1 = T(1, 0);
    const double b2 = T(1, 1);
    const double b3 = T(1, 2);
    const double c1 = T(2, 0);
    const double c2 = T(2, 1);
    const double c3 = T(2, 2);

    const double d1 = L[0];
    const double d2 = L[1];
    const double d3 = L[2];

    const double dnr = (a3 * b2 * c1 - a2 * b3 * c1 - a3 * b1 * c2 + a1 * b3 * c2 + a2 * b1 * c3 - a1 * b2 * c3);

    Eigen::Vector3d W;

    W[0] = (b3 * c2 * d1 - b2 * c3 * d1 - a3 * c2 * d2 + a2 * c3 * d2 + a3 * b2 * d3 - a2 * b3 * d3) / dnr;
    W[1] = -(b3 * c1 * d1 - b1 * c3 * d1 - a3 * c1 * d2 + a1 * c3 * d2 + a3 * b1 * d3 - a1 * b3 * d3) / dnr;
    W[2] = (b2 * c1 * d1 - b1 * c2 * d1 - a2 * c1 * d2 + a1 * c2 * d2 + a2 * b1 * d3 - a1 * b2 * d3) / dnr;

    return W;
} // solve_3x3_system
#endif

void
set_rotation_matrix(const Eigen::Vector3d& rot_vel,
                    const Eigen::Quaterniond& q_old,
                    Eigen::Quaterniond& q_new,
                    Eigen::Matrix3d& rot_mat,
                    const double dt)
{
    const double norm = rot_vel.norm();
    if (!MathUtilities<double>::equalEps(norm, 0.0))
    {
        Eigen::Vector3d rot_axis = rot_vel / norm;
        Eigen::Quaterniond q(Eigen::AngleAxisd(norm * dt, rot_axis));
        q_new = (q.normalized() * q_old).normalized();
    }
    else
    {
        q_new = q_old;
    }

    rot_mat = q_new.toRotationMatrix();
    return;

} // set_rotation_matrix

} // namespace
/////////////////////////////// PUBLIC ///////////////////////////////////////

IBFEDirectForcingKinematics::IBFEDirectForcingKinematics(const std::string& object_name,
                                                         Pointer<Database> input_db,
                                                         Pointer<IBFEMethod> ibfe_method_ops,
                                                         int part,
                                                         bool register_for_restart)
{
    // Set the object name and register it with the restart manager.
    d_object_name = object_name;
    d_ibfe_method_ops = ibfe_method_ops;
    d_part = part;

    d_registered_for_restart = false;
    if (register_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
        d_registered_for_restart = true;
    }

    // Set some default values.
    d_quaternion_current = Eigen::Quaterniond::Identity();
    d_quaternion_half = Eigen::Quaterniond::Identity();
    d_quaternion_new = Eigen::Quaterniond::Identity();
    d_trans_vel_current = Eigen::Vector3d::Zero();
    d_trans_vel_half = Eigen::Vector3d::Zero();
    d_trans_vel_new = Eigen::Vector3d::Zero();
    d_rot_vel_current = Eigen::Vector3d::Zero();
    d_rot_vel_half = Eigen::Vector3d::Zero();
    d_rot_vel_new = Eigen::Vector3d::Zero();
    d_center_of_mass_initial = Eigen::Vector3d::Zero();
    d_center_of_mass_current = Eigen::Vector3d::Zero();
    d_center_of_mass_half = Eigen::Vector3d::Zero();
    d_center_of_mass_new = Eigen::Vector3d::Zero();
    d_inertia_tensor_initial = Eigen::Matrix3d::Zero();

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    return;
} // IBFEDirectForcingKinematics

IBFEDirectForcingKinematics::~IBFEDirectForcingKinematics()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }
    return;
} // ~IBFEDirectForcingKinematics

void
IBFEDirectForcingKinematics::registerKinematicsFunction(KinematicsFcnPtr comvelfcn, void* ctx)
{
    registerKinematicsFunction(KinematicsFcnData(comvelfcn, ctx));
    return;
} // registerKinematicsFunction

void
IBFEDirectForcingKinematics::registerKinematicsFunction(const KinematicsFcnData& data)
{
    d_kinematics = data;
    return;
} // registerKinematicsFunction

void
IBFEDirectForcingKinematics::setSolveRigidBodyVelocity(const IBTK::FreeRigidDOFVector& solve_rigid_dofs)
{
    d_solve_rigid_vel = solve_rigid_dofs;
    return;
} // setSolveRigidBodyVelocity

void
IBFEDirectForcingKinematics::initializeKinematicsData(bool initial_time)
{
    if (initial_time)
    {
        computeCOMOfStructure(d_center_of_mass_initial);
        d_center_of_mass_current = d_center_of_mass_initial;

        computeMOIOfStructure(d_inertia_tensor_initial, d_center_of_mass_initial);
    }

    return;
} // initializeKinematicsData

void
IBFEDirectForcingKinematics::preprocessIntegrateData(double current_time, double new_time, int /*num_cycles*/)
{
    const double dt = new_time - current_time;
    d_current_time = current_time;
    d_new_time = new_time;
    d_half_time = current_time + 0.5 * dt;

    // Obtain the COM velocities.
    // Note: Only imposed DOFs will be used from them.
    Eigen::Vector3d trans_vel_current, trans_vel_half, trans_vel_new;
    Eigen::Vector3d rot_vel_current, rot_vel_half, rot_vel_new;

    if (d_kinematics.comvelfcn)
    {
        d_kinematics.comvelfcn(d_current_time, trans_vel_current, rot_vel_current, d_kinematics.ctx);
        d_kinematics.comvelfcn(d_half_time, trans_vel_half, rot_vel_half, d_kinematics.ctx);
        d_kinematics.comvelfcn(d_new_time, trans_vel_new, rot_vel_new, d_kinematics.ctx);
    }

    // Override imposed DOFs in internal vectors.
    for (int k = 0; k < NDIM; ++k)
    {
        if (!d_solve_rigid_vel[k])
        {
            d_trans_vel_current[k] = trans_vel_current[k];
            d_trans_vel_half[k] = trans_vel_half[k];
            d_trans_vel_new[k] = trans_vel_new[k];
        }
    }
#if (NDIM == 2)
    if (!d_solve_rigid_vel[2])
    {
        d_rot_vel_current[2] = rot_vel_current[2];
        d_rot_vel_half[2] = rot_vel_half[2];
        d_rot_vel_new[2] = rot_vel_new[2];
    }
#elif (NDIM == 3)
    for (int k = 0; k < NDIM; ++k)
    {
        if (!d_solve_rigid_vel[k + NDIM])
        {
            d_rot_vel_current[k] = rot_vel_current[k];
            d_rot_vel_half[k] = rot_vel_half[k];
            d_rot_vel_new[k] = rot_vel_new[k];
        }
    }
#endif

    return;
} // preprocessIntegrateData

void
IBFEDirectForcingKinematics::forwardEulerStep(double current_time,
                                              double new_time,
                                              PetscVector<double>& /*X_current_petsc*/,
                                              PetscVector<double>& X_half_petsc,
                                              PetscVector<double>& X_new_petsc)
{
    const double dt = new_time - current_time;

    // Fill the rotation matrix of structure for half and new times.
    Eigen::Matrix3d rotation_mat_half = Eigen::Matrix3d::Identity(3, 3);
    Eigen::Matrix3d rotation_mat_new = Eigen::Matrix3d::Identity(3, 3);
    set_rotation_matrix(d_rot_vel_current, d_quaternion_current, d_quaternion_half, rotation_mat_half, 0.5 * dt);
    set_rotation_matrix(d_rot_vel_current, d_quaternion_current, d_quaternion_new, rotation_mat_new, dt);

    // Rotate the body with current rotational velocity about center of mass
    // and translate the body to predicted position.
    Eigen::Vector3d dr = Eigen::Vector3d::Zero();
    Eigen::Vector3d Rxdr_half = Eigen::Vector3d::Zero();
    Eigen::Vector3d Rxdr_new = Eigen::Vector3d::Zero();

    EquationSystems* equation_systems = d_ibfe_method_ops->getFEDataManager(d_part)->getEquationSystems();
    System& X_system = equation_systems->get_system(IBFEMethod::COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();

    MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int total_local_nodes = mesh.n_nodes_on_proc(SAMRAI_MPI::getRank());
    std::vector<std::vector<numeric_index_type> > nodal_X_indices(NDIM);
    std::vector<std::vector<double> > nodal_X0_values(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        nodal_X_indices[d].reserve(total_local_nodes);
        nodal_X0_values[d].reserve(total_local_nodes);
    }

    IBFEMethod::CoordinateMappingFcnData mapping = d_ibfe_method_ops->getInitialCoordinateMappingFunction(d_part);
    const bool identity_mapping = !(mapping.fcn);
    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
    {
        const Node* const n = *it;
        const libMesh::Point& s = *n;
        libMesh::Point X0 = s;
        if (!identity_mapping)
        {
            mapping.fcn(X0, s, mapping.ctx);
        }
        if (n->n_vars(X_sys_num))
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM);
#endif
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                nodal_X_indices[d].push_back(n->dof_number(X_sys_num, d, 0));
                nodal_X0_values[d].push_back(X0(d));
            }
        }
    }

    for (unsigned int k = 0; k < total_local_nodes; ++k)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dr[d] = nodal_X0_values[d][k] - d_center_of_mass_initial[d];
        }

        // Rotate dr vector using the rotation matrix.
        Rxdr_half = rotation_mat_half * dr;
        Rxdr_new = rotation_mat_new * dr;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_half_petsc.set(nodal_X_indices[d][k],
                             d_center_of_mass_current[d] + Rxdr_half[d] + 0.5 * dt * d_trans_vel_current[d]);
            X_new_petsc.set(nodal_X_indices[d][k],
                            d_center_of_mass_current[d] + Rxdr_new[d] + dt * d_trans_vel_current[d]);
        }
    }
    X_half_petsc.close();
    X_new_petsc.close();

    // Compute the COM at mid-step and new time.
    d_center_of_mass_half = d_center_of_mass_current + 0.5 * dt * d_trans_vel_current;
    d_center_of_mass_new = d_center_of_mass_current + dt * d_trans_vel_current;

    return;
} // forwardEulerStep

void
IBFEDirectForcingKinematics::midpointStep(double current_time,
                                          double new_time,
                                          PetscVector<double>& /*X_current_petsc*/,
                                          PetscVector<double>& X_half_petsc,
                                          PetscVector<double>& X_new_petsc)
{
    const double dt = new_time - current_time;

    // Fill the rotation matrix of structure for new time.
    Eigen::Matrix3d rotation_mat = Eigen::Matrix3d::Identity(3, 3);
    set_rotation_matrix(d_rot_vel_half, d_quaternion_current, d_quaternion_new, rotation_mat, dt);

    // Rotate the body with current rotational velocity about center of mass
    // and translate the body to predicted position.
    Eigen::Vector3d dr = Eigen::Vector3d::Zero();
    Eigen::Vector3d Rxdr = Eigen::Vector3d::Zero();

    EquationSystems* equation_systems = d_ibfe_method_ops->getFEDataManager(d_part)->getEquationSystems();
    System& X_system = equation_systems->get_system(IBFEMethod::COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();

    MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int total_local_nodes = mesh.n_nodes_on_proc(SAMRAI_MPI::getRank());
    std::vector<std::vector<numeric_index_type> > nodal_X_indices(NDIM);
    std::vector<std::vector<double> > nodal_X0_values(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        nodal_X_indices[d].reserve(total_local_nodes);
        nodal_X0_values[d].reserve(total_local_nodes);
    }

    IBFEMethod::CoordinateMappingFcnData mapping = d_ibfe_method_ops->getInitialCoordinateMappingFunction(d_part);
    const bool identity_mapping = !(mapping.fcn);
    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
    {
        const Node* const n = *it;
        const libMesh::Point& s = *n;
        libMesh::Point X0 = s;
        if (!identity_mapping)
        {
            mapping.fcn(X0, s, mapping.ctx);
        }
        if (n->n_vars(X_sys_num))
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM);
#endif
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                nodal_X_indices[d].push_back(n->dof_number(X_sys_num, d, 0));
                nodal_X0_values[d].push_back(X0(d));
            }
        }
    }

    for (unsigned int k = 0; k < total_local_nodes; ++k)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dr[d] = nodal_X0_values[d][k] - d_center_of_mass_initial[d];
        }

        // Rotate dr vector using the rotation matrix.
        Rxdr = rotation_mat * dr;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_half_petsc.set(nodal_X_indices[d][k], d_center_of_mass_current[d] + Rxdr[d] + dt * d_trans_vel_half[d]);
        }
    }
    X_half_petsc.close();
    X_new_petsc = X_half_petsc;
    X_new_petsc.close();

    // Move and rotate the structure.
    d_center_of_mass_half = d_center_of_mass_current + dt * d_trans_vel_half;
    d_center_of_mass_new = d_center_of_mass_half;

    return;
} // midpointStep

void
IBFEDirectForcingKinematics::trapezoidalStep(double /*current_time*/,
                                             double /*new_time*/,
                                             PetscVector<double>& /*X_current_petsc*/,
                                             PetscVector<double>& /*X_half_petsc*/,
                                             PetscVector<double>& /*X_new_petsc*/)
{
    TBOX_ERROR(d_object_name
               << "::IBFEDirectForcingKinematics does not support trapezoidal time-stepping rule for position update."
               << " Only midpoint rule is supported for position update." << std::endl);
    return;
} // trapezoidalStep

void
IBFEDirectForcingKinematics::computeLagrangianForce(PetscVector<double>& F_petsc,
                                                    PetscVector<double>& X_petsc,
                                                    PetscVector<double>& U_petsc,
                                                    const double data_time)
{
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));

    bool all_imposed_dofs = true;
    for (int k = 0; k < s_max_free_dofs; ++k)
    {
        all_imposed_dofs = all_imposed_dofs && !d_solve_rigid_vel[k];
    }

    if (all_imposed_dofs)
    {
        computeImposedLagrangianForceDensity(F_petsc, X_petsc, U_petsc, data_time);
    }
    else
    {
        computeMixedLagrangianForceDensity(F_petsc, X_petsc, U_petsc, data_time);
    }

    return;
} // computeLagrangianForce

void
IBFEDirectForcingKinematics::postprocessIntegrateData(double /*current_time*/, double /*new_time*/, int /*num_cycles*/)
{
    d_trans_vel_current = d_trans_vel_new;
    d_rot_vel_current = d_rot_vel_new;
    d_center_of_mass_current = d_center_of_mass_new;
    d_quaternion_current = d_quaternion_new;

    return;
} // postprocessIntegrateData

void
IBFEDirectForcingKinematics::putToDatabase(Pointer<Database> db)
{
    std::ostringstream U, W, C0, C, Q;
    U << "U";
    W << "W";
    C0 << "C0";
    C << "C";
    Q << "Q";

    double Q_coeffs[4] = {
        d_quaternion_current.w(), d_quaternion_current.x(), d_quaternion_current.y(), d_quaternion_current.z()
    };

    db->putDoubleArray(U.str(), &d_trans_vel_current[0], 3);
    db->putDoubleArray(W.str(), &d_rot_vel_current[0], 3);
    db->putDoubleArray(C0.str(), &d_center_of_mass_initial[0], 3);
    db->putDoubleArray(C.str(), &d_center_of_mass_current[0], 3);
    db->putDoubleArray(Q.str(), &Q_coeffs[0], 4);

    return;
} // putToDatabase

/////////////////////////////// PRIVATE ////////////////////////////////////

void
IBFEDirectForcingKinematics::getFromInput(Pointer<Database> input_db, bool is_from_restart)
{
    // Get some input values.
    d_rho = input_db->getDouble("rho_s");

    if (!is_from_restart && input_db->keyExists("init_quaternion"))
    {
        double Q_coeffs[4];
        input_db->getDoubleArray("init_quaternion", Q_coeffs, 4);

        d_quaternion_current.w() = Q_coeffs[0];
        d_quaternion_current.x() = Q_coeffs[1];
        d_quaternion_current.y() = Q_coeffs[2];
        d_quaternion_current.z() = Q_coeffs[3];

        d_quaternion_half = d_quaternion_current;
        d_quaternion_new = d_quaternion_current;

        d_quaternion_current.normalized();
        d_quaternion_half.normalized();
        d_quaternion_new.normalized();
    }

    return;
} // getFromInput

void
IBFEDirectForcingKinematics::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR("IBFEDirectForcingKinematics::getFromRestart(): Restart database corresponding to "
                   << d_object_name << " not found in restart file." << std::endl);
    }

    std::ostringstream U, W, C0, C, Q;
    U << "U";
    W << "W";
    C0 << "C0";
    C << "C";
    Q << "Q";

    double Q_coeffs[4];
    db->getDoubleArray(U.str(), &d_trans_vel_current[0], 3);
    db->getDoubleArray(W.str(), &d_rot_vel_current[0], 3);
    db->getDoubleArray(C0.str(), &d_center_of_mass_initial[0], 3);
    db->getDoubleArray(C.str(), &d_center_of_mass_current[0], 3);
    db->getDoubleArray(Q.str(), &Q_coeffs[0], 4);

    d_quaternion_current.w() = Q_coeffs[0];
    d_quaternion_current.x() = Q_coeffs[1];
    d_quaternion_current.y() = Q_coeffs[2];
    d_quaternion_current.z() = Q_coeffs[3];
    d_quaternion_current.normalized();

    return;
} // getFromRestart

void
IBFEDirectForcingKinematics::computeCOMOfStructure(Eigen::Vector3d& X0)
{
    // Get the FE data.
    EquationSystems* equation_systems = d_ibfe_method_ops->getFEDataManager(d_part)->getEquationSystems();
    MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    // Extract the FE system and DOF map, and setup the FE object.
    System& X_system = equation_systems->get_system(IBFEMethod::COORDS_SYSTEM_NAME);
    X_system.solution->localize(*X_system.current_local_solution);
    DofMap& X_dof_map = X_system.get_dof_map();
    std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
    FEType fe_type = X_dof_map.variable_type(0);
    UniquePtr<QBase> qrule = fe_type.default_quadrature_rule(dim);

    UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
    fe->attach_quadrature_rule(qrule.get());
    const std::vector<double>& JxW = fe->get_JxW();
    const std::vector<std::vector<double> >& phi = fe->get_phi();

    // Extract the nodal coordinates.
    PetscVector<double>& X_petsc = dynamic_cast<PetscVector<double>&>(*X_system.current_local_solution.get());
    X_petsc.close();
    Vec X_global_vec = X_petsc.vec();
    Vec X_local_ghost_vec;
    VecGhostGetLocalForm(X_global_vec, &X_local_ghost_vec);
    double* X_local_ghost_soln;
    VecGetArray(X_local_ghost_vec, &X_local_ghost_soln);

    // Loop over the local elements to compute the local integrals.
    boost::multi_array<double, 2> X_node;
    double X_qp[NDIM];
    double vol_part = 0.0;
    X0.setZero();
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_dof_map.dof_indices(elem, X_dof_indices[d], d);
        }
        get_values_for_interpolation(X_node, X_petsc, X_local_ghost_soln, X_dof_indices);

        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(X_qp, qp, X_node, phi);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X0[d] += X_qp[d] * JxW[qp];
            }
            vol_part += JxW[qp];
        }
    }
    SAMRAI_MPI::sumReduction(&X0[0], NDIM);
    vol_part = SAMRAI_MPI::sumReduction(vol_part);
    X0 /= vol_part;

    VecRestoreArray(X_local_ghost_vec, &X_local_ghost_soln);
    VecGhostRestoreLocalForm(X_global_vec, &X_local_ghost_vec);

    return;
} // computeCOMOfStructure

void
IBFEDirectForcingKinematics::computeMOIOfStructure(Eigen::Matrix3d& I, const Eigen::Vector3d& X0)
{
    // Get the FE data.
    EquationSystems* equation_systems = d_ibfe_method_ops->getFEDataManager(d_part)->getEquationSystems();
    MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    // Extract the FE system and DOF map, and setup the FE object.
    System& X_system = equation_systems->get_system(IBFEMethod::COORDS_SYSTEM_NAME);
    X_system.solution->localize(*X_system.current_local_solution);
    DofMap& X_dof_map = X_system.get_dof_map();
    std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
    FEType fe_type = X_dof_map.variable_type(0);
    UniquePtr<QBase> qrule = fe_type.default_quadrature_rule(dim);

    UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
    fe->attach_quadrature_rule(qrule.get());
    const std::vector<double>& JxW = fe->get_JxW();
    const std::vector<std::vector<double> >& phi = fe->get_phi();

    // Extract the nodal coordinates.
    int ierr;
    PetscVector<double>& X_petsc = dynamic_cast<PetscVector<double>&>(*X_system.current_local_solution.get());
    X_petsc.close();
    Vec X_global_vec = X_petsc.vec();
    Vec X_local_ghost_vec;
    ierr = VecGhostGetLocalForm(X_global_vec, &X_local_ghost_vec);
    IBTK_CHKERRQ(ierr);
    double* X_local_ghost_soln;
    ierr = VecGetArray(X_local_ghost_vec, &X_local_ghost_soln);
    IBTK_CHKERRQ(ierr);

    // Loop over the local elements to compute the local integrals.
    boost::multi_array<double, 2> X_node;
    double X_qp[NDIM];
    I.setZero();
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_dof_map.dof_indices(elem, X_dof_indices[d], d);
        }
        get_values_for_interpolation(X_node, X_petsc, X_local_ghost_soln, X_dof_indices);

        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(X_qp, qp, X_node, phi);
#if (NDIM == 2)
            I(0, 0) += std::pow(X_qp[1] - X0[1], 2) * JxW[qp];
            I(0, 1) += -(X_qp[0] - X0[0]) * (X_qp[1] - X0[1]) * JxW[qp];
            I(1, 1) += std::pow(X_qp[0] - X0[0], 2) * JxW[qp];
            I(2, 2) += (std::pow(X_qp[0] - X0[0], 2) + std::pow(X_qp[1] - X0[1], 2)) * JxW[qp];
#endif

#if (NDIM == 3)
            I(0, 0) += (std::pow(X_qp[1] - X0[1], 2) + std::pow(X_qp[2] - X0[2], 2)) * JxW[qp];
            I(0, 1) += -(X_qp[0] - X0[0]) * (X_qp[1] - X0[1]) * JxW[qp];
            I(0, 2) += -(X_qp[0] - X0[0]) * (X_qp[2] - X0[2]) * JxW[qp];
            I(1, 1) += (std::pow(X_qp[0] - X0[0], 2) + std::pow(X_qp[2] - X0[2], 2)) * JxW[qp];
            I(1, 2) += -(X_qp[1] - X0[1]) * (X_qp[2] - X0[2]) * JxW[qp];
            I(2, 2) += (std::pow(X_qp[0] - X0[0], 2) + std::pow(X_qp[1] - X0[1], 2)) * JxW[qp];

#endif
        }
    }
    SAMRAI_MPI::sumReduction(&I(0, 0), 9);

    // Fill-in the symmetric part of inertia tensor.
    I(1, 0) = I(0, 1);
    I(2, 0) = I(0, 2);
    I(2, 1) = I(1, 2);

    ierr = VecRestoreArray(X_local_ghost_vec, &X_local_ghost_soln);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostRestoreLocalForm(X_global_vec, &X_local_ghost_vec);
    IBTK_CHKERRQ(ierr);

    return;
} // computeMOIOfStructure

void
IBFEDirectForcingKinematics::computeImposedLagrangianForceDensity(PetscVector<double>& F_petsc,
                                                                  PetscVector<double>& X_petsc,
                                                                  PetscVector<double>& U_petsc,
                                                                  const double /*data_time*/)
{
    const Eigen::Vector3d& X_com = d_center_of_mass_half;
    EquationSystems* equation_systems = d_ibfe_method_ops->getFEDataManager(d_part)->getEquationSystems();
    System& X_system = equation_systems->get_system(IBFEMethod::COORDS_SYSTEM_NAME);
    System& U_system = equation_systems->get_system(IBFEMethod::VELOCITY_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    const unsigned int U_sys_num = U_system.number();

    MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int total_local_nodes = mesh.n_nodes_on_proc(SAMRAI_MPI::getRank());
    std::vector<std::vector<numeric_index_type> > nodal_indices(NDIM);
    std::vector<std::vector<double> > nodal_X_values(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        nodal_indices[d].reserve(total_local_nodes);
        nodal_X_values[d].resize(total_local_nodes);
    }

    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
    {
        const Node* const n = *it;
        if (n->n_vars(X_sys_num) && n->n_vars(U_sys_num))
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM && n->n_vars(U_sys_num) == NDIM);
#endif
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                nodal_indices[d].push_back(n->dof_number(X_sys_num, d, 0));
            }
        }
    }

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        X_petsc.get(nodal_indices[d], &nodal_X_values[d][0]);
    }

    // Set the cross-product matrix
    Eigen::Matrix3d W(Eigen::Matrix3d::Zero());
#if (NDIM == 2)
    W(0, 1) = -d_rot_vel_half[2]; // -U[2];
    W(1, 0) = d_rot_vel_half[2];  // U[2];
#elif (NDIM == 3)
    W(0, 1) = -d_rot_vel_half[2]; //-U[5];
    W(1, 0) = d_rot_vel_half[2];  // U[5];
    W(0, 2) = d_rot_vel_half[1];  // U[4];
    W(2, 0) = -d_rot_vel_half[1]; // -U[4];
    W(1, 2) = -d_rot_vel_half[0]; //-U[3];
    W(2, 1) = d_rot_vel_half[0];  // U[3];
#endif
    Eigen::Vector3d R(Eigen::Vector3d::Zero());
    Eigen::Vector3d WxR(Eigen::Vector3d::Zero());
    for (unsigned int k = 0; k < total_local_nodes; ++k)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            R[d] = nodal_X_values[d][k] - X_com[d];
        }

        WxR = W * R;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            F_petsc.set(nodal_indices[d][k], d_trans_vel_half[d] + WxR[d]);
        }
    }
    F_petsc.close();

    // Compute the direct forcing as F_df = rho/dt * (U_b - U)
    const double dt = d_new_time - d_current_time;
    int ierr;
    ierr = VecAXPY(F_petsc.vec(), -1.0, U_petsc.vec());
    IBTK_CHKERRQ(ierr);
    PetscScalar alpha = d_rho / dt;
    ierr = VecScale(F_petsc.vec(), alpha);
    IBTK_CHKERRQ(ierr);

    return;
} // computeImposedLagrangianForceDensity

void
IBFEDirectForcingKinematics::computeMixedLagrangianForceDensity(PetscVector<double>& F_petsc,
                                                                PetscVector<double>& X_petsc,
                                                                PetscVector<double>& U_petsc,
                                                                const double data_time)
{
    const Eigen::Vector3d& X_com = d_center_of_mass_half;

    // Get FE objects.
    EquationSystems* equation_systems = d_ibfe_method_ops->getFEDataManager(d_part)->getEquationSystems();
    MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    // Extract the FE system and DOF map, and setup the FE object.
    System& X_system = equation_systems->get_system(IBFEMethod::COORDS_SYSTEM_NAME);
    System& U_system = equation_systems->get_system(IBFEMethod::VELOCITY_SYSTEM_NAME);
    DofMap& X_dof_map = X_system.get_dof_map();
    DofMap& U_dof_map = U_system.get_dof_map();
    std::vector<std::vector<unsigned int> > U_dof_indices(NDIM);
    std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
    FEType fe_type = U_dof_map.variable_type(0);
    UniquePtr<QBase> qrule = fe_type.default_quadrature_rule(dim);
    UniquePtr<FEBase> fe_autoptr(FEBase::build(dim, fe_type));
    FEBase* fe = fe_autoptr.get();
    fe->attach_quadrature_rule(qrule.get());
    const std::vector<double>& JxW = fe->get_JxW();
    const std::vector<std::vector<double> >& phi = fe->get_phi();

    int ierr;
    X_petsc.close();
    Vec X_global_vec = X_petsc.vec();
    Vec X_local_ghost_vec;
    ierr = VecGhostGetLocalForm(X_global_vec, &X_local_ghost_vec);
    IBTK_CHKERRQ(ierr);
    double* X_local_ghost_soln;
    ierr = VecGetArray(X_local_ghost_vec, &X_local_ghost_soln);
    IBTK_CHKERRQ(ierr);

    U_petsc.close();
    Vec U_global_vec = U_petsc.vec();
    Vec U_local_ghost_vec;
    ierr = VecGhostGetLocalForm(U_global_vec, &U_local_ghost_vec);
    IBTK_CHKERRQ(ierr);
    double* U_local_ghost_soln;
    ierr = VecGetArray(U_local_ghost_vec, &U_local_ghost_soln);
    IBTK_CHKERRQ(ierr);

    Eigen::Vector3d L = Eigen::Vector3d::Zero();
    Eigen::Vector3d F = Eigen::Vector3d::Zero();
    boost::multi_array<double, 2> X_node, U_node;
    double X_qp[NDIM], U_qp[NDIM];
    double vol_mesh = 0.0;
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_dof_map.dof_indices(elem, X_dof_indices[d], d);
            U_dof_map.dof_indices(elem, U_dof_indices[d], d);
        }
        get_values_for_interpolation(U_node, U_petsc, U_local_ghost_soln, U_dof_indices);
        get_values_for_interpolation(X_node, X_petsc, X_local_ghost_soln, X_dof_indices);

        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(X_qp, qp, X_node, phi);
            interpolate(U_qp, qp, U_node, phi);

            vol_mesh += JxW[qp];
            for (int d = 0; d < NDIM; ++d)
            {
                F[d] += U_qp[d] * JxW[qp];
            }

#if (NDIM == 2)
            L[2] += (U_qp[1] * (X_qp[0] - X_com[0]) - U_qp[0] * (X_qp[1] - X_com[1])) * JxW[qp];
#elif (NDIM == 3)
            L[0] += (U_qp[2] * (X_qp[1] - X_com[1]) - U_qp[1] * (X_qp[2] - X_com[2])) * JxW[qp];
            L[1] += (U_qp[0] * (X_qp[2] - X_com[2]) - U_qp[2] * (X_qp[0] - X_com[0])) * JxW[qp];
            L[2] += (U_qp[1] * (X_qp[0] - X_com[0]) - U_qp[0] * (X_qp[1] - X_com[1])) * JxW[qp];
#endif
        }
    }
    SAMRAI_MPI::sumReduction(&F[0], 3);
    SAMRAI_MPI::sumReduction(&L[0], 3);
    vol_mesh = SAMRAI_MPI::sumReduction(vol_mesh);

    ierr = VecRestoreArray(X_local_ghost_vec, &X_local_ghost_soln);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostRestoreLocalForm(X_global_vec, &X_local_ghost_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecRestoreArray(U_local_ghost_vec, &U_local_ghost_soln);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostRestoreLocalForm(U_global_vec, &U_local_ghost_vec);
    IBTK_CHKERRQ(ierr);

    // Linear momentum balance.
    Eigen::Vector3d U_new = Eigen::Vector3d::Zero();
    U_new = F / vol_mesh;

    // Angular momentum balance.
    // NOTE: Here we rotate the angular momentum vector back to initial time
    // configuration and solve the system there. The angular velocity is then
    // rotated back to the new configuration.
    Eigen::Vector3d W_new = Eigen::Vector3d::Zero();
    Eigen::Matrix3d R = d_quaternion_half.toRotationMatrix();
    Eigen::Vector3d L0 = (R.transpose()) * L;
#if (NDIM == 2)
    W_new[2] = L0[2] / d_inertia_tensor_initial(2, 2);
#elif (NDIM == 3)
    W_new = solve_3x3_system(L0, d_inertia_tensor_initial);
#endif
    W_new = R * W_new;

    // Override free DOFs in new vector.
    for (int k = 0; k < NDIM; ++k)
    {
        if (d_solve_rigid_vel[k])
        {
            d_trans_vel_new[k] = U_new[k];
        }
    }
#if (NDIM == 2)
    if (d_solve_rigid_vel[2])
    {
        d_rot_vel_new[2] = W_new[2];
    }
#elif (NDIM == 3)
    for (int k = 0; k < NDIM; ++k)
    {
        if (d_solve_rigid_vel[k + NDIM])
        {
            d_rot_vel_new[k] = W_new[k];
        }
    }
#endif

    d_trans_vel_half = 0.5 * (d_trans_vel_new + d_trans_vel_current);
    d_rot_vel_half = 0.5 * (d_rot_vel_new + d_rot_vel_current);

    computeImposedLagrangianForceDensity(F_petsc, X_petsc, U_petsc, data_time);

    return;
} // computeMixedLagrangianForceDensity

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace

//////////////////////////////////////////////////////////////////////////////
