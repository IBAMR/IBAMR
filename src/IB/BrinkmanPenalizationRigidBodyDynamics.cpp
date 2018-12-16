// Filename BrinkmanPenalizationRigidBodyDynamics.cpp
// Created on Dec 05, 2018 by Amneet Bhalla
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/BrinkmanPenalizationRigidBodyDynamics.h"
#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/IBHydrodynamicSurfaceForceEvaluator.h"
#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"
#include "ibamr/namespaces.h"
#include "tbox/RestartManager.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////
namespace
{
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
        q.normalize();
        q_new = (q * q_old).normalized();
    }
    else
    {
        q_new = q_old;
    }

    rot_mat = q_new.toRotationMatrix();
    return;

} // set_rotation_matrix

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

inline void
get_physical_coordinate(Eigen::Vector3d& side_coord, Pointer<Patch<NDIM> > patch, const SideIndex<NDIM>& side_idx)
{
    const int axis = side_idx.getAxis();
    Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
    const double* patch_X_lower = patch_geom->getXLower();
    const Box<NDIM>& patch_box = patch->getBox();
    const Index<NDIM>& patch_lower_idx = patch_box.lower();
    const double* const patch_dx = patch_geom->getDx();

    for (int d = 0; d < NDIM; ++d)
    {
        if (d == axis)
        {
            side_coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(side_idx(d) - patch_lower_idx(d)));
        }
        else
        {
            side_coord[d] =
                patch_X_lower[d] + patch_dx[d] * (static_cast<double>(side_idx(d) - patch_lower_idx(d)) + 0.5);
        }
    }
    return;
} // get_physical_coordinate

} // namespace

/////////////////////////////// PUBLIC //////////////////////////////////////
BrinkmanPenalizationRigidBodyDynamics::BrinkmanPenalizationRigidBodyDynamics(
    std::string object_name,
    Pointer<CellVariable<NDIM, double> > ls_solid_var,
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
    Pointer<INSVCStaggeredHierarchyIntegrator> fluid_solver,
    Pointer<Database> input_db,
    bool register_for_restart)
    : BrinkmanPenalizationStrategy(object_name, register_for_restart)
{
    d_ls_solid_var = ls_solid_var;
    d_adv_diff_solver = adv_diff_solver;
    d_fluid_solver = fluid_solver;

    d_hydro_force_eval = new IBHydrodynamicSurfaceForceEvaluator(
        d_object_name + "::hydro_force", ls_solid_var, adv_diff_solver, fluid_solver);
    d_hydro_force_eval->setSurfaceContourLevel(0.0);
    d_hydro_force_eval->writeToFile(false);

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    return;
} // BrinkmanPenalizationRigidBodyDynamics

BrinkmanPenalizationRigidBodyDynamics::~BrinkmanPenalizationRigidBodyDynamics()
{
    // intentionally left blank.
    return;
} // ~BrinkmanPenalizationStrategy

void
BrinkmanPenalizationRigidBodyDynamics::setInitialConditions(const Eigen::Vector3d& X_com,
                                                            const Eigen::Vector3d& U_com,
                                                            const Eigen::Vector3d& W_com,
                                                            const double& mass,
                                                            const Eigen::Matrix3d& J_com,
                                                            const Eigen::Quaterniond& quaternion)
{
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) return;

    d_center_of_mass_initial = X_com;
    d_center_of_mass_current = X_com;
    d_inertia_tensor_initial = J_com;
    d_mass = mass;

    d_quaternion_current = quaternion.normalized();
    d_quaternion_new = quaternion.normalized();

    d_trans_vel_current = U_com;
    d_rot_vel_current = W_com;
    d_trans_vel_new = U_com;
    d_rot_vel_new = W_com;

    return;
} // setInitialCOMPosnVelocity

void
BrinkmanPenalizationRigidBodyDynamics::registerKinematicsFunction(KinematicsFcnPtr comvelfcn, void* ctx)
{
    registerKinematicsFunction(KinematicsFcnData(comvelfcn, ctx));
    return;
} // registerKinematicsFunction

void
BrinkmanPenalizationRigidBodyDynamics::registerKinematicsFunction(const KinematicsFcnData& data)
{
    d_kinematics_fcn_data = data;
    return;
} // registerKinematicsFunction

void
BrinkmanPenalizationRigidBodyDynamics::setSolveRigidBodyVelocity(const IBTK::FreeRigidDOFVector& solve_rigid_dofs)
{
    d_solve_rigid_vel = solve_rigid_dofs;
    return;
} // setSolveRigidBodyVelocity

void
BrinkmanPenalizationRigidBodyDynamics::registerExternalForceTorqueFunction(ExternalForceTorqueFcnPtr forcetorquefcn,
                                                                           void* ctx)
{
    registerExternalForceTorqueFunction(ExternalForceTorqueFcnData(forcetorquefcn, ctx));
    return;
} // registerExternalForceTorqueFunction

void
BrinkmanPenalizationRigidBodyDynamics::registerExternalForceTorqueFunction(const ExternalForceTorqueFcnData& data)
{
    d_ext_force_torque_fcn_data = data;
    return;
} // registerExternalForceTorqueFunction

void
BrinkmanPenalizationRigidBodyDynamics::preprocessComputeBrinkmanPenalization(double current_time,
                                                                             double new_time,
                                                                             int num_cycles)
{
    d_trans_vel_new = d_trans_vel_current;
    d_rot_vel_new = d_rot_vel_current;
    d_center_of_mass_new = d_center_of_mass_current;
    d_quaternion_new = d_quaternion_current;

    return;
} // preprocessComputeBrinkmanPenalization

void
BrinkmanPenalizationRigidBodyDynamics::computeBrinkmanVelocity(int u_idx, double time, int cycle_num)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(time, d_new_time));
#endif

    const double hydro_compute_time = cycle_num > 0 ? time : d_current_time;
    d_hydro_force_eval->computeHydrodynamicForceTorque(d_hydro_force_pressure,
                                                       d_hydro_force_viscous,
                                                       d_hydro_torque_pressure,
                                                       d_hydro_torque_viscous,
                                                       d_center_of_mass_new,
                                                       hydro_compute_time,
                                                       d_current_time,
                                                       d_new_time);

    if (d_ext_force_torque_fcn_data.forcetorquefcn)
    {
        d_ext_force_torque_fcn_data.forcetorquefcn(time, d_ext_force, d_ext_torque, d_ext_force_torque_fcn_data.ctx);
    }
    else
    {
        d_ext_force.setZero();
        d_ext_torque.setZero();
    }

    // Get imposed motion.
    Eigen::Vector3d U_imposed = Eigen::Vector3d::Zero(), W_imposed = Eigen::Vector3d::Zero();
    if (d_kinematics_fcn_data.comvelfcn)
    {
        d_kinematics_fcn_data.comvelfcn(time, U_imposed, W_imposed, d_kinematics_fcn_data.ctx);
    }

    // Integrate Newton's second law of motion to find updated COM velocity and position.
    // a) Translational motion
    const double dt = d_new_time - d_current_time;
    Eigen::Vector3d F_rigid;
    F_rigid = d_hydro_force_pressure + d_hydro_force_viscous + d_ext_force;
    d_trans_vel_new = d_trans_vel_current + (dt * F_rigid) / d_mass;
    for (unsigned s = 0; s < NDIM; ++s)
    {
        if (!d_solve_rigid_vel(s))
        {
            d_trans_vel_new(s) = U_imposed(s);
        }
    }
    d_center_of_mass_new = d_center_of_mass_current + dt * d_trans_vel_new;

    // b) Rotational motion
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity(3, 3);
    set_rotation_matrix(d_rot_vel_new, d_quaternion_current, d_quaternion_new, R, dt);
    Eigen::Vector3d T_rigid = d_hydro_torque_pressure + d_hydro_torque_viscous + d_ext_torque;
    d_rot_vel_new.setZero();
#if (NDIM == 2)
    if (d_solve_rigid_vel(2))
    {
        d_rot_vel_new(2) = d_rot_vel_current(2) + (dt * T_rigid(2)) / d_inertia_tensor_initial(2, 2);
    }
#elif (NDIM == 3)
    Eigen::Vector3d T0_rigid = solve_3x3_system((R.transpose()) * T_rigid, d_inertia_tensor_initial);
    d_rot_vel_new = R * (d_rot_vel_current + dt * T0_rigid);
    for (unsigned s = NDIM; s < s_max_free_dofs; ++s)
    {
        if (!d_solve_rigid_vel(s))
        {
            d_rot_vel_new(s) = W_imposed(s);
        }
    }
#endif

    // Get the interpolated density variable
    const int rho_ins_idx = d_fluid_solver->getLinearOperatorRhoPatchDataIndex();

    // Ghost fill the level set values.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int ls_solid_idx = var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getNewContext());
    const int ls_scratch_idx =
        var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getScratchContext());

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_adv_diff_solver->getPatchHierarchy();
    int finest_ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM> > finest_level = patch_hierarchy->getPatchLevel(finest_ln);
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_ops(patch_hierarchy);
    hier_cc_ops.copyData(ls_scratch_idx, ls_solid_idx);

    Pointer<RefineAlgorithm<NDIM> > ghost_fill_alg = new RefineAlgorithm<NDIM>();
    Pointer<RefineOperator<NDIM> > refine_op = NULL;
    ghost_fill_alg->registerRefine(ls_scratch_idx, ls_scratch_idx, ls_scratch_idx, refine_op);
    ghost_fill_alg->createSchedule(finest_level)->fillData(time);

    // Set the rigid body velocity in u_idx
    Eigen::Vector3d r = Eigen::Vector3d::Zero();
    Eigen::Vector3d dr = Eigen::Vector3d::Zero();
    Eigen::Vector3d Wxdr = Eigen::Vector3d::Zero();
    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* patch_dx = patch_geom->getDx();
        const double alpha = 2.0 * patch_dx[0];

        Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(ls_scratch_idx);
        Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_idx);
        Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_ins_idx);
        TBOX_ASSERT((ls_solid_data->getGhostCellWidth()).min() >= 1);

        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
            {
                SideIndex<NDIM> s_i(it(), axis, SideIndex<NDIM>::Lower);

                const double phi_lower = (*ls_solid_data)(s_i.toCell(0));
                const double phi_upper = (*ls_solid_data)(s_i.toCell(1));
                const double phi = 0.5 * (phi_lower + phi_upper);
                double Hphi;
                if (phi < -alpha)
                {
                    Hphi = 0.0;
                }
                else if (std::abs(phi) <= alpha)
                {
                    Hphi = 0.5 + 0.5 * phi / alpha + 1.0 / (2.0 * M_PI) * std::sin(M_PI * phi / alpha);
                }
                else
                {
                    Hphi = 1.0;
                }
                if (phi <= alpha)
                {
                    get_physical_coordinate(r, patch, s_i);
                    dr = r - d_center_of_mass_new;
                    Wxdr = d_rot_vel_new.cross(dr);
                    const double penalty = (*rho_data)(s_i) / dt;
                    (*u_data)(s_i) = d_trans_vel_new(axis) + Wxdr(axis);
                    (*u_data)(s_i) *= (1.0 - Hphi) * penalty;
                }
            }
        }
    }

    return;
} // computeBrinkmanVelocity

void
BrinkmanPenalizationRigidBodyDynamics::demarcateBrinkmanZone(int u_idx, double time, int /*cycle_num*/)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(time, d_new_time));
#endif

    const double dt = d_new_time - d_current_time;

    // Get the interpolated density variable
    const int rho_ins_idx = d_fluid_solver->getLinearOperatorRhoPatchDataIndex();

    // Get the solid level set data. We have already copied new data into scratch data and have also
    // filled its ghost cells on the finest level.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int ls_scratch_idx =
        var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getScratchContext());
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_adv_diff_solver->getPatchHierarchy();
    int finest_ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM> > finest_level = patch_hierarchy->getPatchLevel(finest_ln);
    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* patch_dx = patch_geom->getDx();
        const double alpha = 2.0 * patch_dx[0];

        Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(ls_scratch_idx);
        Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_idx);
        Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_ins_idx);

        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
            {
                SideIndex<NDIM> s_i(it(), axis, SideIndex<NDIM>::Lower);

                const double phi_lower = (*ls_solid_data)(s_i.toCell(0));
                const double phi_upper = (*ls_solid_data)(s_i.toCell(1));
                const double phi = 0.5 * (phi_lower + phi_upper);
                double Hphi;
                if (phi < -alpha)
                {
                    Hphi = 0.0;
                }
                else if (std::abs(phi) <= alpha)
                {
                    Hphi = 0.5 + 0.5 * phi / alpha + 1.0 / (2.0 * M_PI) * std::sin(M_PI * phi / alpha);
                }
                else
                {
                    Hphi = 1.0;
                }

                if (phi <= alpha)
                {
                    const double penalty = (*rho_data)(s_i) / dt;
                    (*u_data)(s_i) = (1.0 - Hphi) * penalty;
                }
            }
        }
    }

    return;
} // demarcateBrinkmanZone

void
BrinkmanPenalizationRigidBodyDynamics::postprocessComputeBrinkmanPenalization(double current_time,
                                                                              double new_time,
                                                                              int num_cycles)
{
    BrinkmanPenalizationStrategy::postprocessComputeBrinkmanPenalization(current_time, new_time, num_cycles);

    d_trans_vel_current = d_trans_vel_new;
    d_rot_vel_current = d_rot_vel_new;
    d_center_of_mass_current = d_center_of_mass_new;
    d_quaternion_current = d_quaternion_new;

    return;
} // postprocessComputeBrinkmanPenalization

void
BrinkmanPenalizationRigidBodyDynamics::putToDatabase(Pointer<Database> db)
{
    std::ostringstream U, W, C0, C, J0, Q, M;
    U << "U";
    W << "W";
    C0 << "C0";
    C << "C";
    J0 << "J0";
    Q << "Q";
    M << "M";

    double Q_coeffs[4] = {
        d_quaternion_current.w(), d_quaternion_current.x(), d_quaternion_current.y(), d_quaternion_current.z()
    };

    db->putDoubleArray(U.str(), &d_trans_vel_current[0], 3);
    db->putDoubleArray(W.str(), &d_rot_vel_current[0], 3);
    db->putDoubleArray(C0.str(), &d_center_of_mass_initial[0], 3);
    db->putDoubleArray(C.str(), &d_center_of_mass_current[0], 3);
    db->putDoubleArray(J0.str(), &d_inertia_tensor_initial(0, 0), 9);
    db->putDoubleArray(Q.str(), &Q_coeffs[0], 4);
    db->putDouble(M.str(), d_mass);

    return;
} // postprocessComputeBrinkmanVelocity

/////////////////////////////// PRIVATE //////////////////////////////////////

void
BrinkmanPenalizationRigidBodyDynamics::getFromInput(Pointer<Database> input_db, bool is_from_restart)
{
    if (!is_from_restart && input_db->keyExists("init_quaternion"))
    {
        double Q_coeffs[4];
        input_db->getDoubleArray("init_quaternion", Q_coeffs, 4);

        d_quaternion_current.w() = Q_coeffs[0];
        d_quaternion_current.x() = Q_coeffs[1];
        d_quaternion_current.y() = Q_coeffs[2];
        d_quaternion_current.z() = Q_coeffs[3];

        d_quaternion_current.normalize();
        d_quaternion_new = d_quaternion_current;
    }

    if (input_db->keyExists("chi"))
    {
        d_chi = input_db->getDouble("chi");
    }
    return;
} // getFromInput

void
BrinkmanPenalizationRigidBodyDynamics::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR("BrinkmanPenalizationRigidBodyDynamics::getFromRestart(): Restart database corresponding to "
                   << d_object_name << " not found in restart file." << std::endl);
    }

    std::ostringstream U, W, C0, C, J0, Q, M;
    U << "U";
    W << "W";
    C0 << "C0";
    C << "C";
    J0 << "J0";
    Q << "Q";
    M << "M";

    double Q_coeffs[4];
    db->getDoubleArray(U.str(), &d_trans_vel_current[0], 3);
    db->getDoubleArray(W.str(), &d_rot_vel_current[0], 3);
    db->getDoubleArray(C0.str(), &d_center_of_mass_initial[0], 3);
    db->getDoubleArray(C.str(), &d_center_of_mass_current[0], 3);
    db->getDoubleArray(J0.str(), &d_inertia_tensor_initial(0, 0), 9);
    db->getDoubleArray(Q.str(), &Q_coeffs[0], 4);
    d_mass = db->getDouble(M.str());

    d_quaternion_current.w() = Q_coeffs[0];
    d_quaternion_current.x() = Q_coeffs[1];
    d_quaternion_current.y() = Q_coeffs[2];
    d_quaternion_current.z() = Q_coeffs[3];
    d_quaternion_current.normalize();

    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
