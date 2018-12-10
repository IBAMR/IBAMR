// Filename: BrinkmanPenalizationRigidBodyDynamics.h
// Created on 05 Dec 2018 by Amneet Bhalla
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

#ifndef included_BrinkmanPenalizationRigidBodyDynamics
#define included_BrinkmanPenalizationRigidBodyDynamics

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "ibamr/BrinkmanPenalizationStrategy.h"
#include "ibtk/ibtk_utilities.h"
#include "tbox/Pointer.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
class IBHydrodynamicSurfaceForceEvaluator;
class INSVCStaggeredHierarchyIntegrator;
class AdvDiffHierarchyIntegrator;
} // namespace IBAMR
namespace SAMRAI
{
namespace pdat
{
template <int DIM, class TYPE>
class CellVariable;
} // namespace pdat
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief BrinkmanPenalizationRigidBodyDynamics provides an implementation of Brinkman penalization
 * body force for a rigid body motion in the momentum equation.
 */
class BrinkmanPenalizationRigidBodyDynamics : public BrinkmanPenalizationStrategy
{
    ////////////////////////////// PUBLIC ////////////////////////////////////////

public:
    /*
     * \brief Constructor of the class.
     */
    BrinkmanPenalizationRigidBodyDynamics(std::string object_name,
                                          SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_solid_var,
                                          SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                                          SAMRAI::tbox::Pointer<IBAMR::INSVCStaggeredHierarchyIntegrator> fluid_solver,
                                          SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db = nullptr,
                                          bool register_for_restart = true);

    /*
     * \brief Destructor of the class.
     */
    ~BrinkmanPenalizationRigidBodyDynamics();

    /*!
     * \brief Set initial conditions for the rigid body.
     *
     * \note It is OK to set some or all of these quantities to zero, because the fully-prescribed
     * kinematics case does not need any information of mass and moment of inertia.
     */
    void setInitialConditions(const Eigen::Vector3d& X_com,
                              const Eigen::Vector3d& U_com = Eigen::Vector3d::Zero(),
                              const Eigen::Vector3d& W_com = Eigen::Vector3d::Zero(),
                              const double& mass = 0.0,
                              const Eigen::Matrix3d& J_com = Eigen::Matrix3d::Zero(),
                              const Eigen::Quaterniond& quaternion = Eigen::Quaterniond::Identity());

    /*!
     * \brief Typedef specifying interface for specifying rigid body velocities.
     */
    using KinematicsFcnPtr = void (*)(double data_time, Eigen::Vector3d& U_com, Eigen::Vector3d& W_com, void* ctx);

    /*
     * \brief Kinematics function data.
     */
    struct KinematicsFcnData
    {
        KinematicsFcnData(KinematicsFcnPtr comvelfcn = nullptr, void* ctx = nullptr) : comvelfcn(comvelfcn), ctx(ctx)
        {
            // intentionally left blank
        }
        KinematicsFcnPtr comvelfcn;
        void* ctx;
    };

    /*!
     * \brief Register user defined constrained velocity functions.
     */
    void registerKinematicsFunction(KinematicsFcnPtr comvelfcn, void* ctx = nullptr);

    /*!
     * \brief Register user defined velocity function data.
     */
    void registerKinematicsFunction(const KinematicsFcnData& data);

    /*!
     * \brief Set what rigid DOFs need to be solved for this
     * particular structure.
     */
    virtual void setSolveRigidBodyVelocity(const IBTK::FreeRigidDOFVector& solve_rigid_dofs);

    /*!
     * \brief Typedef specifying interface for specifying additional rigid body force and torque.
     */
    using ExternalForceTorqueFcnPtr = void (*)(double data_time, Eigen::Vector3d& F, Eigen::Vector3d& T, void* ctx);

    /*
     * \brief External force/torque function data.
     */
    struct ExternalForceTorqueFcnData
    {
        ExternalForceTorqueFcnData(ExternalForceTorqueFcnPtr forcetorquefcn = nullptr, void* ctx = nullptr)
            : forcetorquefcn(forcetorquefcn), ctx(ctx)
        {
            // intentionally left blank
        }
        ExternalForceTorqueFcnPtr forcetorquefcn;
        void* ctx;
    };

    /*!
     * \brief Register user defined external force/torque functions.
     */
    void registerExternalForceTorqueFunction(ExternalForceTorqueFcnPtr forcetorquefcn, void* ctx = nullptr);

    /*!
     * \brief Register user defined external force/torque function data.
     */
    void registerExternalForceTorqueFunction(const ExternalForceTorqueFcnData& data);

    /*!
     * \brief Preprocess routine before computing Brinkman penalization terms.
     *
     */
    void preprocessComputeBrinkmanPenalization(double current_time, double new_time, int num_cycles) override;

    /*!
     * \brief Compute the desired rigid body velocity in the Brinkman penalized zone.
     */
    void computeBrinkmanVelocity(int u_idx, double time, int cycle_num) override;

    /*!
     * \brief Demarcate the Brinkman zone with Brinkman penalty term.
     */
    void demarcateBrinkmanZone(int u_idx, double time, int cycle_num) override;

    /*!
     * \brief Postprocess routine after computing Brinkman penalization terms.
     */
    void postprocessComputeBrinkmanPenalization(double current_time, double new_time, int num_cycles) override;

    /*!
     * \brief Write out object state to the given database.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    /*!
     * \brief Get initial center of mass of the body.
     */
    const Eigen::Vector3d& getInitialCOMPosn() const
    {
        return d_center_of_mass_initial;
    } // getInitialCOMPosn

    /*!
     * \brief Get initial inertia tensor of the body.
     */
    const Eigen::Matrix3d& getInitialInteriaTensor() const
    {
        return d_inertia_tensor_initial;
    } // getInitialInteriaTensor

    /*!
     * Get COM position.
     */
    const Eigen::Vector3d& getCurrentCOMPosn() const
    {
        return d_center_of_mass_current;
    } // getCurrentCOMPosn

    const Eigen::Vector3d& getNewCOMPosn() const
    {
        return d_center_of_mass_new;
    } // getNewCOMPosn

    /*!
     * Get COM translational velocity.
     */
    const Eigen::Vector3d& getCurrentCOMTransVelocity() const
    {
        return d_trans_vel_current;
    } // getCurrentCOMTransVelocity

    const Eigen::Vector3d& getNewCOMTransVelocity() const
    {
        return d_trans_vel_new;
    } // getNewCOMTransVelocity

    /*!
     * Get COM rotational velocity.
     */
    const Eigen::Vector3d& getCurrentCOMRotVelocity() const
    {
        return d_rot_vel_current;
    } // getCurrentCOMRotVelocity

    const Eigen::Vector3d& getNewCOMRotVelocity() const
    {
        return d_rot_vel_new;
    } // getNewCOMRotVelocity

    /////////////////////////////// PROTECTED ////////////////////////////////////

protected:
    // Center of mass.
    Eigen::Vector3d d_center_of_mass_initial = Eigen::Vector3d::Zero(),
                    d_center_of_mass_current = Eigen::Vector3d::Zero(), d_center_of_mass_new = Eigen::Vector3d::Zero();

    // Quaternion of the body.
    Eigen::Quaterniond d_quaternion_current = Eigen::Quaterniond::Identity(),
                       d_quaternion_new = Eigen::Quaterniond::Identity();

    // Indicate which rigid degrees of freedom to solve.
    IBTK::FRDV d_solve_rigid_vel;

    // Rigid body velocity of the structure.
    Eigen::Vector3d d_trans_vel_current = Eigen::Vector3d::Zero(), d_trans_vel_new = Eigen::Vector3d::Zero();
    Eigen::Vector3d d_rot_vel_current = Eigen::Vector3d::Zero(), d_rot_vel_new = Eigen::Vector3d::Zero();

    // Mass and inertial of the body.
    double d_mass = 0.0;
    Eigen::Matrix3d d_inertia_tensor_initial = Eigen::Matrix3d::Zero();

    // Hydrodynamic force evaluator.
    SAMRAI::tbox::Pointer<IBAMR::IBHydrodynamicSurfaceForceEvaluator> d_hydro_force_eval;

    // Pointers to solvers.
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;
    SAMRAI::tbox::Pointer<IBAMR::INSVCStaggeredHierarchyIntegrator> d_fluid_solver;

    // Level set variable defining the solid.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_solid_var;

    // Forces and torques on the body.
    Eigen::Vector3d d_hydro_force_pressure, d_hydro_force_viscous, d_hydro_torque_pressure, d_hydro_torque_viscous,
        d_ext_force, d_ext_torque;

    // Routines to get prescribed kinematics and additional external forces and torques.
    KinematicsFcnData d_kinematics_fcn_data;
    ExternalForceTorqueFcnData d_ext_force_torque_fcn_data;

    /////////////////////////////// PROTECTED ////////////////////////////////////

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    BrinkmanPenalizationRigidBodyDynamics(const BrinkmanPenalizationRigidBodyDynamics& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    BrinkmanPenalizationRigidBodyDynamics& operator=(const BrinkmanPenalizationRigidBodyDynamics& that) = delete;

    /*!
     * \brief Get options from input database.
     */
    /*!
     * Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db, bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();
};

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_BrinkmanPenalizationRigidBodyDynamics
