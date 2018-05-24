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
#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "ibtk/FEDataManager.h"
#include "ibtk/ibtk_utilities.h"
#include "tbox/Pointer.h"

namespace IBTK
{
class RobinPhysBdryPatchStrategy;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BasePatchHierarchy;
template <int DIM>
class BasePatchLevel;
} // namespace hier
namespace tbox
{
class Database;
template <class TYPE>
class Array;
} // namespace tbox
namespace xfer
{
template <int DIM>
class CoarsenSchedule;
template <int DIM>
class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI
namespace libMesh
{
class EquationSystems;
class Mesh;
class Point;
class System;
template <typename T>
class NumericVector;
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
                                libMesh::Mesh* mesh,
                                IBTK::FEDataManager* fe_data_manager,
                                int level_number,
                                double d_rho,
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
     * Return a pointer to the finite element data manager object.
     */
    IBTK::FEDataManager* getFEDataManager() const;

    /*!
     * Advance the positions of the Lagrangian structure using the forward Euler
     * method.
     */
    void forwardEulerStep(double current_time,
                          double new_time,
                          PetscVector<double>& X_current_vec,
                          PetscVector<double>& X_half_vec,
                          PetscVector<double>& X_new_vec);

    /*!
     * Advance the positions of the Lagrangian structure using the (explicit)
     * midpoint rule.
     */
    void midpointStep(double current_time,
                      double new_time,
                      PetscVector<double>& X_current_vec,
                      PetscVector<double>& X_half_vec,
                      PetscVector<double>& X_new_vec);

    /*!
     * Advance the positions of the Lagrangian structure using the (explicit)
     * trapezoidal rule.
     */
    void trapezoidalStep(double current_time,
                         double new_time,
                         PetscVector<double>& X_current_vec,
                         PetscVector<double>& X_half_vec,
                         PetscVector<double>& X_new_vec);

    /*!
     * Compute the Lagrangian force at the specified time within the current
     * time interval.
     */
    void computeLagrangianForce(PetscVector<double>& F_vec,
                                PetscVector<double>& X_vec,
                                PetscVector<double>& U_vec,
                                const double data_time);

    /*
     * Pointers to the patch hierarchy and gridding algorithm objects associated
     * with this object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > d_gridding_alg;

    /*
     * The current time step interval.
     */
    double d_current_time, d_new_time, d_half_time;

    /*
     * FE data associated with this object.
     */
    libMesh::Mesh* d_mesh;
    int d_level_number;
    libMesh::EquationSystems* d_equation_systems;
    IBTK::FEDataManager* d_fe_data_manager;
    SAMRAI::hier::IntVector<NDIM> d_ghosts;
    libMesh::System *d_X_systems, *d_U_systems, *d_F_systems;
    libMesh::PetscVector<double>*d_X_current_vec, *d_X_new_vec, *d_X_half_vec, *d_X_IB_ghost_vec;
    libMesh::PetscVector<double>*d_U_current_vec, *d_U_new_vecs, *d_U_half_vecs;
    libMesh::PetscVector<double>*d_F_half_vecs, *d_F_IB_ghost_vecs;

    std::string d_object_name;

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

    // Rigid body velocity of the structures.
    Eigen::Vector3d d_trans_vel_current, d_trans_vel_half, d_trans_vel_new;
    Eigen::Vector3d d_rot_vel_current, d_rot_vel_half, d_rot_vel_new;

    /*
     * A boolean value indicating whether the class is registered with the
     * restart database.
     */
    bool d_registered_for_restart;

    /*!
     * \brief Write out object state to the given database.
     *
     * \note An empty default implementation is provided.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

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
    void computeImposedLagrangianForceDensity(PetscVector<double>& F_vec,
                                              PetscVector<double>& X_vec,
                                              PetscVector<double>& U_vec,
                                              const double data_time);

    /*!
     * Compute the Lagrangian force at the specified time within the current
     * time interval for structure that has some free DOFs.
     */
    void computeFreeLagrangianForceDensity(PetscVector<double>& F_vec,
                                           PetscVector<double>& X_vec,
                                           PetscVector<double>& U_vec,
                                           const double data_time);
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_IBFEDirectForcingKinematics
