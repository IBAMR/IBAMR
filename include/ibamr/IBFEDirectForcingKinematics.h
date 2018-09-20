// Filename: IBFEDirectForcingKinematics.h
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

#ifndef included_IBAMR_IBFEDirectForcingKinematics
#define included_IBAMR_IBFEDirectForcingKinematics

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <set>
#include <stdbool.h>
#include <stddef.h>
#include <string>
#include <vector>

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "ibtk/ibtk_utilities.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"

namespace SAMRAI
{
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI
namespace IBAMR
{
class IBFEMethod;
} // namespace IBAMR
namespace libMesh
{
template <typename T>
class PetscVector;
} // namespace libMesh

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBFEDirectForcingKinematics is a helper class that provides direct
 * forcing IBMethod functionality to the IBFEMethod class.
 */
class IBFEDirectForcingKinematics : public virtual SAMRAI::tbox::DescribedClass, public SAMRAI::tbox::Serializable
{
public:
    /*!
     * \brief Constructor.
     */
    IBFEDirectForcingKinematics(const std::string& object_name,
                                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                SAMRAI::tbox::Pointer<IBAMR::IBFEMethod> ibfe_method_ops,
                                int part,
                                bool register_for_restart = true);

    /*!
     * \brief Destructor.
     */
    ~IBFEDirectForcingKinematics();

    /*!
     * \brief Typedef specifying interface for specifying rigid body velocities.
     */
    typedef void (*KinematicsFcnPtr)(double data_time, Eigen::Vector3d& U_com, Eigen::Vector3d& W_com, void* ctx);

    /*
     * \brief Kinematics function data.
     */
    struct KinematicsFcnData
    {
        KinematicsFcnData(KinematicsFcnPtr comvelfcn = NULL, void* ctx = NULL) : comvelfcn(comvelfcn), ctx(ctx)
        {
            // intentionally left blank
        }
        KinematicsFcnPtr comvelfcn;
        void* ctx;
    };

    /*!
     * \brief Register user defined constrained velocity functions.
     */
    void registerKinematicsFunction(KinematicsFcnPtr comvelfcn, void* ctx = NULL);

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
     * Initialize kinematics data only initially.
     */
    virtual void initializeKinematicsData(bool initial_time = true);

    /*!
     * Preprocess kinematics before hierarchy integrates.
     */
    virtual void preprocessIntegrateData(double current_time, double new_time, int num_cycles);

    /*!
     * Advance the positions of the Lagrangian structure using the forward Euler
     * method.
     */
    virtual void forwardEulerStep(double current_time,
                                  double new_time,
                                  libMesh::PetscVector<double>& X_current_petsc,
                                  libMesh::PetscVector<double>& X_half_petsc,
                                  libMesh::PetscVector<double>& X_new_petsc);

    /*!
     * Advance the positions of the Lagrangian structure using the (explicit)
     * midpoint rule.
     */
    virtual void midpointStep(double current_time,
                              double new_time,
                              libMesh::PetscVector<double>& X_current_petsc,
                              libMesh::PetscVector<double>& X_half_petsc,
                              libMesh::PetscVector<double>& X_new_petsc);

    /*!
     * Advance the positions of the Lagrangian structure using the (explicit)
     * trapezoidal rule.
     */
    virtual void trapezoidalStep(double current_time,
                                 double new_time,
                                 libMesh::PetscVector<double>& X_current_petsc,
                                 libMesh::PetscVector<double>& X_half_petsc,
                                 libMesh::PetscVector<double>& X_new_petsc);

    /*!
     * Compute the Lagrangian force at the specified time within the current
     * time interval.
     */
    virtual void computeLagrangianForce(libMesh::PetscVector<double>& F_petsc,
                                        libMesh::PetscVector<double>& X_petsc,
                                        libMesh::PetscVector<double>& U_petsc,
                                        const double data_time);

    /*!
     * Postprocess kinematics after hierarchy integrates.
     */
    virtual void postprocessIntegrateData(double current_time, double new_time, int num_cycles);

    /*!
     * Get COM position.
     */
    const Eigen::Vector3d& getStructureCOM() const
    {
        return d_center_of_mass_current;
    } // getStructureCOM

    /*!
     * Get COM translational velocity.
     */
    const Eigen::Vector3d& getCOMTransVelocity() const
    {
        return d_trans_vel_current;
    } // getCOMTransVelocity

    /*!
     * Get COM rotational velocity.
     */
    const Eigen::Vector3d& getCOMRotVelocity() const
    {
        return d_rot_vel_current;
    } // getCOMRotVelocity

    /*!
     * \brief Write out object state to the given database.
     *
     * \note An empty default implementation is provided.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

protected:
    /*!
     * Compute center of mass of the structure.
     */
    void computeCOMOfStructure(Eigen::Vector3d& X0);

    /*!
     * Compute moment of inertia tensor for the structure.
     */
    void computeMOIOfStructure(Eigen::Matrix3d& I, const Eigen::Vector3d& X0);

    /*
     * The current time step interval.
     */
    double d_current_time, d_new_time, d_half_time;

    /*
     * IBFE method and the part number registered with IBFE method.
     */
    SAMRAI::tbox::Pointer<IBAMR::IBFEMethod> d_ibfe_method_ops;
    int d_part;

    /*
     * Book-keeping
     */
    std::string d_object_name;

    /*
     * Rigid body kinematics of the body.
     */
    KinematicsFcnData d_kinematics;

    /*
     * Density of the solid.
     */
    double d_rho;

    // Center of mass.
    Eigen::Vector3d d_center_of_mass_initial, d_center_of_mass_current, d_center_of_mass_half, d_center_of_mass_new;

    // Quaternion of the body.
    Eigen::Quaterniond d_quaternion_current, d_quaternion_half, d_quaternion_new;

    // Indicate which rigid degrees of freedom to solve.
    IBTK::FRDV d_solve_rigid_vel;

    // Rigid body velocity of the structure.
    Eigen::Vector3d d_trans_vel_current, d_trans_vel_half, d_trans_vel_new;
    Eigen::Vector3d d_rot_vel_current, d_rot_vel_half, d_rot_vel_new;
    Eigen::Matrix3d d_inertia_tensor_initial;

    /*
     * A boolean value indicating whether the class is registered with the
     * restart database.
     */
    bool d_registered_for_restart;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBFEDirectForcingKinematics();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBFEDirectForcingKinematics(const IBFEDirectForcingKinematics& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBFEDirectForcingKinematics& operator=(const IBFEDirectForcingKinematics& that);

    /*!
     * Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db, bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();

    /*!
     * Compute the Lagrangian force at the specified time within the current
     * time interval for structure that has all the imposed DOFs.
     */
    void computeImposedLagrangianForceDensity(libMesh::PetscVector<double>& F_vec,
                                              libMesh::PetscVector<double>& X_vec,
                                              libMesh::PetscVector<double>& U_vec,
                                              const double data_time);

    /*!
     * Compute the Lagrangian force at the specified time within the current
     * time interval for structure that has both kinds of DOFs (free and imposed).
     */
    void computeMixedLagrangianForceDensity(libMesh::PetscVector<double>& F_vec,
                                            libMesh::PetscVector<double>& X_vec,
                                            libMesh::PetscVector<double>& U_vec,
                                            const double data_time);
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_IBFEDirectForcingKinematics
