// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_IBFEDirectForcingKinematics
#define included_IBAMR_IBFEDirectForcingKinematics

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibtk/ibtk_utilities.h>

#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include "Eigen/Core"
#include "Eigen/Geometry"

IBTK_ENABLE_EXTRA_WARNINGS

#include <set>
#include <string>
#include <vector>

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
     * Since this class has Eigen object members, which have special alignment
     * requirements, we must explicitly override operator new to get the
     * correct aligment for the object as a whole.
     */
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /*!
     * \brief Constructor.
     */
    IBFEDirectForcingKinematics(std::string object_name,
                                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                SAMRAI::tbox::Pointer<IBAMR::IBFEMethod> ibfe_method_ops,
                                int part,
                                bool register_for_restart = true);

    /*!
     * \brief Destructor.
     */
    virtual ~IBFEDirectForcingKinematics();

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
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

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
    Eigen::Vector3d d_center_of_mass_initial = Eigen::Vector3d::Zero(),
                    d_center_of_mass_current = Eigen::Vector3d::Zero(), d_center_of_mass_half = Eigen::Vector3d::Zero(),
                    d_center_of_mass_new = Eigen::Vector3d::Zero();

    // Quaternion of the body.
    Eigen::Quaterniond d_quaternion_current = Eigen::Quaterniond::Identity(),
                       d_quaternion_half = Eigen::Quaterniond::Identity(),
                       d_quaternion_new = Eigen::Quaterniond::Identity();

    // Indicate which rigid degrees of freedom to solve.
    IBTK::FRDV d_solve_rigid_vel;

    // Rigid body velocity of the structure.
    Eigen::Vector3d d_trans_vel_current = Eigen::Vector3d::Zero(), d_trans_vel_half = Eigen::Vector3d::Zero(),
                    d_trans_vel_new = Eigen::Vector3d::Zero();
    Eigen::Vector3d d_rot_vel_current = Eigen::Vector3d::Zero(), d_rot_vel_half = Eigen::Vector3d::Zero(),
                    d_rot_vel_new = Eigen::Vector3d::Zero();
    Eigen::Matrix3d d_inertia_tensor_initial = Eigen::Matrix3d::Zero();

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
    IBFEDirectForcingKinematics() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBFEDirectForcingKinematics(const IBFEDirectForcingKinematics& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBFEDirectForcingKinematics& operator=(const IBFEDirectForcingKinematics& that) = delete;

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
