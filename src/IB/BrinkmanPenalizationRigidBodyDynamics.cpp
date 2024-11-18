// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2024 by the IBAMR developers
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

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/BrinkmanPenalizationRigidBodyDynamics.h"
#include "ibamr/BrinkmanPenalizationStrategy.h"
#include "ibamr/IBHydrodynamicSurfaceForceEvaluator.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"

#include "ibtk/IndexUtilities.h"
#include "ibtk/ibtk_utilities.h"

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellVariable.h"
#include "HierarchyCellDataOpsReal.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefineSchedule.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#include "Eigen/src/Geometry/Quaternion.h"
#include "Eigen/src/LU/FullPivLU.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

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
    if (!IBTK::abs_equal_eps(norm, 0.0))
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
    Eigen::FullPivLU<Eigen::Matrix3d> lu(T);
    return lu.solve(L);
} // solve_3x3_system
} // namespace

/////////////////////////////// PUBLIC //////////////////////////////////////
BrinkmanPenalizationRigidBodyDynamics::BrinkmanPenalizationRigidBodyDynamics(
    std::string object_name,
    Pointer<CellVariable<NDIM, double> > ls_solid_var,
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
    Pointer<INSVCStaggeredHierarchyIntegrator> fluid_solver,
    Pointer<Database> input_db,
    bool register_for_restart)
    : BrinkmanPenalizationStrategy(std::move(object_name), register_for_restart),
      d_adv_diff_solver(adv_diff_solver),
      d_fluid_solver(fluid_solver),
      d_ls_solid_var(ls_solid_var),
      d_hydro_force_eval(new IBHydrodynamicSurfaceForceEvaluator(d_object_name + "::hydro_force",
                                                                 d_ls_solid_var,
                                                                 d_adv_diff_solver,
                                                                 d_fluid_solver))
{
    d_hydro_force_eval->writeToFile(false);

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    d_hydro_force_eval->setSurfaceContourLevel(d_contour_level);

    return;
} // BrinkmanPenalizationRigidBodyDynamics

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
BrinkmanPenalizationRigidBodyDynamics::preprocessComputeBrinkmanPenalization(double /*current_time*/,
                                                                             double /*new_time*/,
                                                                             int /*num_cycles*/)
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
    TBOX_ASSERT(IBTK::rel_equal_eps(time, d_new_time));
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
        d_ext_force_torque_fcn_data.forcetorquefcn(
            time, cycle_num, d_ext_force, d_ext_torque, d_ext_force_torque_fcn_data.ctx);
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
        d_kinematics_fcn_data.comvelfcn(time, cycle_num, U_imposed, W_imposed, d_kinematics_fcn_data.ctx);
    }

    // Integrate Newton's second law of motion to find updated COM velocity and
    // position.
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
    const Eigen::Matrix3d R_current = d_quaternion_current.toRotationMatrix();
    Eigen::Matrix3d R_new = Eigen::Matrix3d::Identity(3, 3);
    set_rotation_matrix(d_rot_vel_new, d_quaternion_current, d_quaternion_new, R_new, dt);

    Eigen::Vector3d T_rigid = d_hydro_torque_pressure + d_hydro_torque_viscous + d_ext_torque;
    d_rot_vel_new.setZero();
    if (NDIM == 2)
    {
        d_rot_vel_new(2) = d_rot_vel_current(2) + (dt * T_rigid(2)) / d_inertia_tensor_initial(2, 2);
    }
    else if (NDIM == 3)
    {
        const Eigen::Vector3d IW_new =
            (T_rigid * dt + R_current * d_inertia_tensor_initial * R_current.transpose() * d_rot_vel_current);
        d_rot_vel_new = R_new * solve_3x3_system(R_new.transpose() * IW_new, d_inertia_tensor_initial);
    }
    for (unsigned s = NDIM; s < s_max_free_dofs; ++s)
    {
        if (!d_solve_rigid_vel(s))
        {
            if (NDIM == 2)
            {
                d_rot_vel_new(s) = W_imposed(s);
            }
            else if (NDIM == 3)
            {
                d_rot_vel_new(s - NDIM) = W_imposed(s - NDIM);
            }
        }
    }

    // Ghost fill the level set values.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int ls_solid_idx = var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getNewContext());
    const int ls_scratch_idx =
        var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getScratchContext());

    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;

    std::vector<InterpolationTransactionComponent> transaction_comps;
    transaction_comps.emplace_back(
        InterpolationTransactionComponent(ls_scratch_idx,
                                          ls_solid_idx,
                                          "CONSERVATIVE_LINEAR_REFINE",
                                          false,
                                          "CONSERVATIVE_COARSEN",
                                          "LINEAR",
                                          false,
                                          d_adv_diff_solver->getPhysicalBcCoefs(d_ls_solid_var)));
    Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_adv_diff_solver->getPatchHierarchy();

    if (!d_split_penalty)
    {
        hier_bdry_fill->initializeOperatorState(transaction_comps, patch_hierarchy);
        hier_bdry_fill->fillData(time);

        computeBrinkmanVelocityWithoutSplitting(u_idx, time, cycle_num);
    }
    else
    {
        const int u_new_idx = var_db->mapVariableAndContextToIndex(d_fluid_solver->getVelocityVariable(),
                                                                   d_fluid_solver->getNewContext());
        const int u_scratch_idx = var_db->mapVariableAndContextToIndex(d_fluid_solver->getVelocityVariable(),
                                                                       d_fluid_solver->getScratchContext());

        transaction_comps.emplace_back(
            InterpolationTransactionComponent(u_scratch_idx,
                                              u_new_idx,
                                              "NONE",
                                              false,
                                              "CUBIC_COARSEN",
                                              "LINEAR",
                                              false,
                                              d_fluid_solver->getVelocityBoundaryConditions()));

        hier_bdry_fill->initializeOperatorState(transaction_comps, patch_hierarchy);
        hier_bdry_fill->fillData(time);

        computeBrinkmanVelocityWithSplitting(u_idx, time, cycle_num);
    }

    return;
} // computeBrinkmanVelocity

void
BrinkmanPenalizationRigidBodyDynamics::demarcateBrinkmanZone(int u_idx, double time, int cycle_num)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(IBTK::rel_equal_eps(time, d_new_time));
#else
    NULL_USE(time);
#endif

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int ls_solid_idx = var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getNewContext());
    const int ls_scratch_idx =
        var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getScratchContext());

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_adv_diff_solver->getPatchHierarchy();

    // Ghost fill the level set values.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> phi_transaction_comps(1);
    phi_transaction_comps[0] = InterpolationTransactionComponent(ls_scratch_idx,
                                                                 ls_solid_idx,
                                                                 "CONSERVATIVE_LINEAR_REFINE",
                                                                 false,
                                                                 "CONSERVATIVE_COARSEN",
                                                                 "LINEAR",
                                                                 false,
                                                                 d_adv_diff_solver->getPhysicalBcCoefs(d_ls_solid_var));
    Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_bdry_fill->initializeOperatorState(phi_transaction_comps, patch_hierarchy);
    hier_bdry_fill->fillData(time);

    if (!d_split_penalty)
    {
        demarcateBrinkmanZoneWithoutSplitting(u_idx, time, cycle_num);
    }
    else
    {
        demarcateBrinkmanZoneWithSplitting(u_idx, time, cycle_num);
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
    std::string U = "U", W = "W", C0 = "C0", C = "C", J0 = "J0", Q = "Q", M = "M";

    double Q_coeffs[4] = {
        d_quaternion_current.w(), d_quaternion_current.x(), d_quaternion_current.y(), d_quaternion_current.z()
    };

    db->putDoubleArray(U, &d_trans_vel_current[0], 3);
    db->putDoubleArray(W, &d_rot_vel_current[0], 3);
    db->putDoubleArray(C0, &d_center_of_mass_initial[0], 3);
    db->putDoubleArray(C, &d_center_of_mass_current[0], 3);
    db->putDoubleArray(J0, &d_inertia_tensor_initial(0, 0), 9);
    db->putDoubleArray(Q, &Q_coeffs[0], 4);
    db->putDouble(M, d_mass);

    db->putDouble("num_interface_cells", d_num_interface_cells);
    db->putDouble("contour_level", d_contour_level);
    db->putDouble("penalty_factor", d_penalty_factor);
    db->putDouble("penalty_factor_tangential", d_penalty_factor_tangential);
    db->putDouble("penalty_factor_normal", d_penalty_factor_normal);
    db->putBool("use_rho_scale", d_use_rho_scale);
    db->putBool("use_mu_scale", d_use_mu_scale);
    db->putBool("split_penalty", d_split_penalty);

    return;
} // putToDatabase

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
    if (!is_from_restart)
    {
        if (input_db->keyExists("penalty_factor"))
        {
            d_penalty_factor = input_db->getDouble("penalty_factor");
        }

        if (input_db->keyExists("use_rho_scale"))
        {
            d_use_rho_scale = input_db->getBool("use_rho_scale");
        }

        if (input_db->keyExists("use_mu_scale"))
        {
            d_use_mu_scale = input_db->getBool("use_mu_scale");
        }

        if (input_db->keyExists("contour_level"))
        {
            d_contour_level = input_db->getDouble("contour_level");
        }

        if (input_db->keyExists("num_interface_cells"))
        {
            d_num_interface_cells = input_db->getDouble("num_interface_cells");
        }

        if (input_db->keyExists("split_penalty"))
        {
            d_split_penalty = input_db->getBool("split_penalty");
        }

        if (input_db->keyExists("penalty_factor_tangential"))
        {
            d_penalty_factor_tangential = input_db->getDouble("penalty_factor_tangential");
        }

        if (input_db->keyExists("penalty_factor_normal"))
        {
            d_penalty_factor_normal = input_db->getDouble("penalty_factor_normal");
        }
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
        TBOX_ERROR(
            "BrinkmanPenalizationRigidBodyDynamics::getFromRestart(): "
            "Restart database corresponding to "
            << d_object_name << " not found in restart file." << std::endl);
    }

    std::string U = "U", W = "W", C0 = "C0", C = "C", J0 = "J0", Q = "Q", M = "M";

    double Q_coeffs[4];
    db->getDoubleArray(U, &d_trans_vel_current[0], 3);
    db->getDoubleArray(W, &d_rot_vel_current[0], 3);
    db->getDoubleArray(C0, &d_center_of_mass_initial[0], 3);
    db->getDoubleArray(C, &d_center_of_mass_current[0], 3);
    db->getDoubleArray(J0, &d_inertia_tensor_initial(0, 0), 9);
    db->getDoubleArray(Q, &Q_coeffs[0], 4);
    d_mass = db->getDouble(M);

    d_quaternion_current.w() = Q_coeffs[0];
    d_quaternion_current.x() = Q_coeffs[1];
    d_quaternion_current.y() = Q_coeffs[2];
    d_quaternion_current.z() = Q_coeffs[3];
    d_quaternion_current.normalize();

    d_num_interface_cells = db->getDouble("num_interface_cells");
    d_contour_level = db->getDouble("contour_level");
    d_penalty_factor = db->getDouble("penalty_factor");
    d_penalty_factor_tangential = db->getDouble("penalty_factor_tangential");
    d_penalty_factor_normal = db->getDouble("penalty_factor_normal");
    d_use_rho_scale = db->getBool("use_rho_scale");
    d_use_mu_scale = db->getBool("use_mu_scale");
    d_split_penalty = db->getBool("split_penalty");

    return;
} // getFromRestart

void
BrinkmanPenalizationRigidBodyDynamics::computeBrinkmanVelocityWithoutSplitting(int u_in_idx,
                                                                               double /*time*/,
                                                                               int /*cycle_num*/)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    const int ls_scratch_idx =
        var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getScratchContext());

    // Get the interpolated density variable.
    const int rho_ins_idx = d_fluid_solver->getLinearOperatorRhoPatchDataIndex();

    // Get the cell-centered viscosity patch data index. Returns mu_scratch with ghost cells filled.
    const int mu_ins_idx = d_fluid_solver->getLinearOperatorMuPatchDataIndex();

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_adv_diff_solver->getPatchHierarchy();
    int finest_ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM> > finest_level = patch_hierarchy->getPatchLevel(finest_ln);
    const double dt = d_new_time - d_current_time;

    // Set the rigid body velocity in u_in_idx
    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* patch_dx = patch_geom->getDx();
        const double h_min = *(std::min_element(patch_dx, patch_dx + NDIM));
        double vol_cell = 1.0;
        for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
        const double alpha = d_num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

        Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(ls_scratch_idx);
        Pointer<SideData<NDIM, double> > u_in_data = patch->getPatchData(u_in_idx);
        Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_ins_idx);
        Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_ins_idx);
        TBOX_ASSERT((ls_solid_data->getGhostCellWidth()).min() >= 1);

        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
            {
                SideIndex<NDIM> s_i(it(), axis, SideIndex<NDIM>::Lower);

                const double phi_lower = (*ls_solid_data)(s_i.toCell(0));
                const double phi_upper = (*ls_solid_data)(s_i.toCell(1));
                const double phi = 0.5 * (phi_lower + phi_upper);
                const double Hphi = IBTK::smooth_heaviside(phi, alpha);
                if (phi <= alpha)
                {
                    const Eigen::Vector3d r = IBTK::IndexUtilities::getSideCenter<Eigen::Vector3d>(*patch, s_i);
                    const Eigen::Vector3d dr = r - d_center_of_mass_new;
                    const Eigen::Vector3d Wxdr = d_rot_vel_new.cross(dr);

                    double penalty_rho_scale = 0.0, penalty_mu_scale = 0.0;
                    if (d_use_rho_scale)
                    {
                        penalty_rho_scale = (*rho_data)(s_i) / dt;
                    }
                    if (d_use_mu_scale)
                    {
                        const double mu_lower = (*mu_data)(s_i.toCell(0));
                        const double mu_upper = (*mu_data)(s_i.toCell(1));
                        const double mu = 0.5 * (mu_lower + mu_upper);
                        penalty_mu_scale = mu / (h_min * h_min);
                    }
                    const double penalty = d_penalty_factor * (penalty_rho_scale + penalty_mu_scale);

                    (*u_in_data)(s_i) = d_trans_vel_new(axis) + Wxdr(axis);
                    (*u_in_data)(s_i) *= (1.0 - Hphi) * penalty;
                }
            }
        }
    }

    return;
} // computeBrinkmanVelocityWithoutSplitting

void
BrinkmanPenalizationRigidBodyDynamics::computeBrinkmanVelocityWithSplitting(int u_in_idx,
                                                                            double /*time*/,
                                                                            int /*cycle_num*/)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    const int ls_scratch_idx =
        var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getScratchContext());
    const int u_scratch_idx = var_db->mapVariableAndContextToIndex(d_fluid_solver->getVelocityVariable(),
                                                                   d_fluid_solver->getScratchContext());

    // Get the interpolated density variable.
    const int rho_ins_idx = d_fluid_solver->getLinearOperatorRhoPatchDataIndex();

    // Get the cell-centered viscosity patch data index. Returns mu_scratch with ghost cells filled.
    const int mu_ins_idx = d_fluid_solver->getLinearOperatorMuPatchDataIndex();

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_adv_diff_solver->getPatchHierarchy();
    int finest_ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM> > finest_level = patch_hierarchy->getPatchLevel(finest_ln);
    const double dt = d_new_time - d_current_time;

    // Set the rigid body velocity in u_in_idx
    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* patch_dx = patch_geom->getDx();
        const double h_min = *(std::min_element(patch_dx, patch_dx + NDIM));
        double vol_cell = 1.0;
        for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
        const double alpha = d_num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

        Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(ls_scratch_idx);
        Pointer<SideData<NDIM, double> > u_new_data = patch->getPatchData(u_scratch_idx);
        Pointer<SideData<NDIM, double> > u_in_data = patch->getPatchData(u_in_idx);
        Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_ins_idx);
        Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_ins_idx);
        TBOX_ASSERT((ls_solid_data->getGhostCellWidth()).min() >= 1);

        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
            {
                SideIndex<NDIM> s_i(it(), axis, SideIndex<NDIM>::Lower);
                CellIndex<NDIM> c_l = s_i.toCell(0);
                CellIndex<NDIM> c_u = s_i.toCell(1);
                const double phi_lower = (*ls_solid_data)(c_l);
                const double phi_upper = (*ls_solid_data)(c_u);
                const double phi = 0.5 * (phi_lower + phi_upper);
                const double Hphi = IBTK::smooth_heaviside(phi, alpha);

                // Inside the structure do not split velocity into normal and
                // tangential components
                if (phi < -alpha)
                {
                    const Eigen::Vector3d r = IBTK::IndexUtilities::getSideCenter<Eigen::Vector3d>(*patch, s_i);
                    const Eigen::Vector3d dr = r - d_center_of_mass_new;
                    const Eigen::Vector3d Wxdr = d_rot_vel_new.cross(dr);

                    double penalty_rho_scale = 0.0, penalty_mu_scale = 0.0;
                    if (d_use_rho_scale)
                    {
                        penalty_rho_scale = (*rho_data)(s_i) / dt;
                    }
                    if (d_use_mu_scale)
                    {
                        const double mu_lower = (*mu_data)(s_i.toCell(0));
                        const double mu_upper = (*mu_data)(s_i.toCell(1));
                        const double mu = 0.5 * (mu_lower + mu_upper);
                        penalty_mu_scale = mu / (h_min * h_min);
                    }
                    const double penalty = d_penalty_factor * (penalty_rho_scale + penalty_mu_scale);

                    (*u_in_data)(s_i) = d_trans_vel_new(axis) + Wxdr(axis);
                    (*u_in_data)(s_i) *= (1.0 - Hphi) * penalty;
                }

                // Near the solid interface impose different penalty factors to normal
                // and tangential velocity components
                if (phi <= alpha && phi >= -alpha)
                {
                    const Eigen::Vector3d r = IBTK::IndexUtilities::getSideCenter<Eigen::Vector3d>(*patch, s_i);
                    const Eigen::Vector3d dr = r - d_center_of_mass_new;
                    const Eigen::Vector3d Wxdr = d_rot_vel_new.cross(dr);
                    const Eigen::Vector3d Ub = d_trans_vel_new + Wxdr;

                    // Compute normal at the side center from the solid level set. Interpolate
                    // tangential components of velocity at the side center. We treat tangential
                    // components of the velocity explicit in time
                    Eigen::Vector3d n = Eigen::Vector3d::Zero();
                    Eigen::Vector3d Ut = Eigen::Vector3d::Zero();
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (d == axis)
                        {
                            n(axis) = (phi_upper - phi_lower) / patch_dx[axis];
                            Ut(axis) = 0.0;
                        }
                        else
                        {
                            CellIndex<NDIM> offset(0);
                            offset(d) = 1;

                            const double h = patch_dx[d];
                            n(d) = ((*ls_solid_data)(c_l + offset) - (*ls_solid_data)(c_l)) +
                                   ((*ls_solid_data)(c_u + offset) - (*ls_solid_data)(c_u)) +
                                   ((*ls_solid_data)(c_l) - (*ls_solid_data)(c_l - offset)) +
                                   ((*ls_solid_data)(c_u) - (*ls_solid_data)(c_u - offset));
                            n(d) *= 1.0 / (4.0 * h);

                            Ut(d) = (*u_new_data)(SideIndex<NDIM>(c_u, d, SideIndex<NDIM>::Lower)) +
                                    (*u_new_data)(SideIndex<NDIM>(c_u, d, SideIndex<NDIM>::Upper)) +
                                    (*u_new_data)(SideIndex<NDIM>(c_l, d, SideIndex<NDIM>::Lower)) +
                                    (*u_new_data)(SideIndex<NDIM>(c_l, d, SideIndex<NDIM>::Upper));
                            Ut(d) *= 0.25;
                        }
                    }
                    const double n_norm = n.norm();
                    if (!IBTK::abs_equal_eps(n_norm, 0.0))
                    {
                        n /= n_norm;
                    }
                    else
                    {
                        n.setZero();
                    }
                    const double Ub_dot_n = Ub.dot(n);
                    const double Ut_dot_n = Ut.dot(n);

                    double penalty_rho_scale = 0.0, penalty_mu_scale = 0.0;
                    if (d_use_rho_scale)
                    {
                        penalty_rho_scale = (*rho_data)(s_i) / dt;
                    }
                    if (d_use_mu_scale)
                    {
                        const double mu_lower = (*mu_data)(c_l);
                        const double mu_upper = (*mu_data)(c_u);
                        const double mu = 0.5 * (mu_lower + mu_upper);
                        penalty_mu_scale = mu / (h_min * h_min);
                    }
                    const double k_n = d_penalty_factor_normal * (penalty_rho_scale + penalty_mu_scale);
                    const double k_t = d_penalty_factor_tangential * (penalty_rho_scale + penalty_mu_scale);

                    // The penalized term in the RHS of the Navier-Stokes equation is
                    //
                    //        f_brinkman =   k * (Ub - u)
                    //
                    //        Ub = velocity of the structure
                    //        u = Eulerian velocity
                    //
                    //  We break this force near the interface as
                    //
                    //       f_brinkman =  k_n * (Ub - u) . (n (x) n) + k_t (Ub - u) . (I - n (x) n)
                    //
                    //       (x) denotes outer product
                    //        I  denotes the identity tensor
                    //
                    //  Note the choice of sign of normal to the solid interface is irrelevant
                    //  in the form above thanks to the symmetric outer product.
                    //
                    // While solving INS equations, we treat u(axis)*n(axis) implicitly, and u(d)*n(d) explicitly in
                    // time.

                    // Get the normal velocity contribution at s_i
                    (*u_in_data)(s_i) = (1.0 - Hphi) * k_n * (Ub_dot_n - Ut_dot_n) * n(axis);

                    // Add in the tangential velocity contribution
                    (*u_in_data)(s_i) += (1.0 - Hphi) * k_t * (Ub(axis) - (Ub_dot_n - Ut_dot_n) * n(axis));
                }
            }
        }
    }

    return;
} // computeBrinkmanVelocityWithSplitting

void
BrinkmanPenalizationRigidBodyDynamics::demarcateBrinkmanZoneWithoutSplitting(int u_idx,
                                                                             double /*time*/,
                                                                             int /*cycle_num*/)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    const int ls_scratch_idx =
        var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getScratchContext());

    // Get the interpolated density variable.
    const int rho_ins_idx = d_fluid_solver->getLinearOperatorRhoPatchDataIndex();

    // Get the cell-centered viscosity patch data index. Returns mu_scratch with ghost cells filled.
    const int mu_ins_idx = d_fluid_solver->getLinearOperatorMuPatchDataIndex();

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_adv_diff_solver->getPatchHierarchy();
    int finest_ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM> > finest_level = patch_hierarchy->getPatchLevel(finest_ln);
    const double dt = d_new_time - d_current_time;

    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* patch_dx = patch_geom->getDx();
        const double h_min = *(std::min_element(patch_dx, patch_dx + NDIM));
        double vol_cell = 1.0;
        for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
        const double alpha = d_num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

        Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(ls_scratch_idx);
        Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_idx);
        Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_ins_idx);
        Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_ins_idx);

        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
            {
                SideIndex<NDIM> s_i(it(), axis, SideIndex<NDIM>::Lower);

                const double phi_lower = (*ls_solid_data)(s_i.toCell(0));
                const double phi_upper = (*ls_solid_data)(s_i.toCell(1));
                const double phi = 0.5 * (phi_lower + phi_upper);
                const double Hphi = IBTK::smooth_heaviside(phi, alpha);
                if (phi <= alpha)
                {
                    double penalty_rho_scale = 0.0, penalty_mu_scale = 0.0;
                    if (d_use_rho_scale)
                    {
                        penalty_rho_scale = (*rho_data)(s_i) / dt;
                    }
                    if (d_use_mu_scale)
                    {
                        const double mu_lower = (*mu_data)(s_i.toCell(0));
                        const double mu_upper = (*mu_data)(s_i.toCell(1));
                        const double mu = 0.5 * (mu_lower + mu_upper);
                        penalty_mu_scale = mu / (h_min * h_min);
                    }

                    const double penalty = d_penalty_factor * (penalty_rho_scale + penalty_mu_scale);
                    (*u_data)(s_i) = (1.0 - Hphi) * penalty;
                }
            }
        }
    }

    return;
} // demarcateBrinkmanZoneWithoutSplitting

void
BrinkmanPenalizationRigidBodyDynamics::demarcateBrinkmanZoneWithSplitting(int u_idx, double /*time*/, int /*cycle_num*/)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    const int ls_scratch_idx =
        var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getScratchContext());

    // Get the interpolated density variable.
    const int rho_ins_idx = d_fluid_solver->getLinearOperatorRhoPatchDataIndex();

    // Get the cell-centered viscosity patch data index. Returns mu_scratch with ghost cells filled.
    const int mu_ins_idx = d_fluid_solver->getLinearOperatorMuPatchDataIndex();

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_adv_diff_solver->getPatchHierarchy();
    int finest_ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM> > finest_level = patch_hierarchy->getPatchLevel(finest_ln);
    const double dt = d_new_time - d_current_time;

    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* patch_dx = patch_geom->getDx();
        const double h_min = *(std::min_element(patch_dx, patch_dx + NDIM));
        double vol_cell = 1.0;
        for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
        const double alpha = d_num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

        Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(ls_scratch_idx);
        Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_idx);
        Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_ins_idx);
        Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_ins_idx);

        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
            {
                SideIndex<NDIM> s_i(it(), axis, SideIndex<NDIM>::Lower);
                CellIndex<NDIM> c_l = s_i.toCell(0);
                CellIndex<NDIM> c_u = s_i.toCell(1);
                const double phi_lower = (*ls_solid_data)(c_l);
                const double phi_upper = (*ls_solid_data)(c_u);
                const double phi = 0.5 * (phi_lower + phi_upper);
                const double Hphi = IBTK::smooth_heaviside(phi, alpha);

                // Inside the structure do not split velocity into normal and
                // tangential components
                if (phi < -alpha)
                {
                    double penalty_rho_scale = 0.0, penalty_mu_scale = 0.0;
                    if (d_use_rho_scale)
                    {
                        penalty_rho_scale = (*rho_data)(s_i) / dt;
                    }
                    if (d_use_mu_scale)
                    {
                        const double mu_lower = (*mu_data)(s_i.toCell(0));
                        const double mu_upper = (*mu_data)(s_i.toCell(1));
                        const double mu = 0.5 * (mu_lower + mu_upper);
                        penalty_mu_scale = mu / (h_min * h_min);
                    }

                    const double penalty = d_penalty_factor * (penalty_rho_scale + penalty_mu_scale);
                    (*u_data)(s_i) = (1.0 - Hphi) * penalty;
                }

                // Near the solid interface impose different penalty factors to normal
                // and tangential velocity components
                if (phi <= alpha && phi >= -alpha)
                {
                    // Compute normal at the side center from the solid level set.
                    Eigen::Vector3d n = Eigen::Vector3d::Zero();
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (d == axis)
                        {
                            n(axis) = (phi_upper - phi_lower) / patch_dx[axis];
                        }
                        else
                        {
                            CellIndex<NDIM> offset(0);
                            offset(d) = 1;

                            const double h = patch_dx[d];
                            n(d) = ((*ls_solid_data)(c_l + offset) - (*ls_solid_data)(c_l)) +
                                   ((*ls_solid_data)(c_u + offset) - (*ls_solid_data)(c_u)) +
                                   ((*ls_solid_data)(c_l) - (*ls_solid_data)(c_l - offset)) +
                                   ((*ls_solid_data)(c_u) - (*ls_solid_data)(c_u - offset));
                            n(d) *= 1.0 / (4.0 * h);
                        }
                    }
                    const double n_norm = n.norm();
                    if (!IBTK::abs_equal_eps(n_norm, 0.0))
                    {
                        n /= n_norm;
                    }
                    else
                    {
                        n.setZero();
                    }

                    double penalty_rho_scale = 0.0, penalty_mu_scale = 0.0;
                    if (d_use_rho_scale)
                    {
                        penalty_rho_scale = (*rho_data)(s_i) / dt;
                    }
                    if (d_use_mu_scale)
                    {
                        const double mu_lower = (*mu_data)(c_l);
                        const double mu_upper = (*mu_data)(c_u);
                        const double mu = 0.5 * (mu_lower + mu_upper);
                        penalty_mu_scale = mu / (h_min * h_min);
                    }
                    const double k_n = d_penalty_factor_normal * (penalty_rho_scale + penalty_mu_scale);
                    const double k_t = d_penalty_factor_tangential * (penalty_rho_scale + penalty_mu_scale);

                    // The penalized term in the RHS of the Navier-Stokes equation is
                    //
                    //        f_brinkman =   k * (Ub - u)
                    //
                    //        Ub = velocity of the structure
                    //        u = Eulerian velocity
                    //
                    //  We break this force near the interface as
                    //
                    //       f_brinkman =  k_n * (Ub - u) . (n (x) n) + k_t (Ub - u) . (I - n (x) n)
                    //
                    //       (x) denotes outer product
                    //        I  denotes the identity tensor
                    //
                    //  Note the choice of sign of normal to the solid interface is irrelevant
                    //  in the form above thanks to the symmetric outer product.
                    //
                    //  While solving INS equations, we treat u(axis)*n(axis) implicitly, and u(d)*n(d) explicitly in
                    //  time.

                    // Get the normal velocity contribution at s_i
                    (*u_data)(s_i) = (1.0 - Hphi) * k_n * n(axis) * n(axis);

                    // Add in the tangential velocity contribution
                    (*u_data)(s_i) += (1.0 - Hphi) * k_t * (1.0 - n(axis) * n(axis));
                }
            }
        }
    }

    return;
} // demarcateBrinkmanZoneWithSplitting

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
