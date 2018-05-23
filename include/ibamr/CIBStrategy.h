// Filename: CIBStrategy.h
// Created on 8 Nov 2014 by Amneet Bhalla
//
// Copyright (c) 2002-2017, Amneet Bhalla and Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of its
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

#ifndef included_IBAMR_CIBStrategy
#define included_IBAMR_CIBStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <map>
#include <vector>

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "ibamr/ibamr_enums.h"
#include "petscmat.h"
#include "petscvec.h"
#include "tbox/DescribedClass.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
static const int s_max_free_dofs = NDIM * (NDIM + 1) / 2;
typedef Eigen::Matrix<double, s_max_free_dofs, 1> RigidDOFVector;
typedef Eigen::Matrix<int, s_max_free_dofs, 1> FreeRigidDOFVector;
typedef RigidDOFVector RDV;
typedef FreeRigidDOFVector FRDV;

/*!
 * \brief Class CIBStrategy is a lightweight abstract strategy class which
 * provides support for constraint based IB methods for rigid bodies.
 */
class CIBStrategy : public virtual SAMRAI::tbox::DescribedClass
{
    ////////////////////////////// PUBLIC ////////////////////////////////////////
public:
    /*!
     *  \brief Constructor of the class.
     */
    CIBStrategy(const unsigned int parts);

    /*!
     *  \brief Destructor of the class.
     */
    virtual ~CIBStrategy();

    /*!
     * \brief Prepare the implementation class for sprading constraint force.
     * In particular, set the constraint Lagrangian force in the internal data
     * structure of the class.
     *
     * \param L Vec containing the constraint force for all structures.
     *
     * \param data_time Time at which constraint force is to be spread.
     *
     * \param scale Scales the constraint force before spreading.
     */
    virtual void setConstraintForce(Vec L, const double data_time, const double scale = 1.0) = 0;

    /*!
     * \brief Get the constraint rigid body force at the specified time within
     * the current time interval. Generally, the implementation class maintains
     * and stores the constraint force. This routine is called by constraint solver
     * to update the contraint force after the (converged) solution is obtained.
     *
     * \param time Time (current_time or new_time) at which constraint force is
     * required.
     */
    virtual void getConstraintForce(Vec* L, const double data_time) = 0;

    /*!
     * \brief Get the free rigid velocities (DOFs) at the specified time within
     * the current time interval. This routine is called by constraint solver
     * to update the free rigid DOFs after the (converged) solution is obtained.
     *
     * \note A default implementation is provided that returns the vector of free
     * DOFs.
     *
     * \param time Time (current_time or new_time) at which constraint force is
     * required.
     */
    virtual void getFreeRigidVelocities(Vec* U, const double data_time);

    /*!
     * \brief Get net external force and torque at the specified time within
     * the current time interval. This routine is called by constraint solver
     * to form the appropriate RHS.
     *
     * \note A default implementation is provided that returns the vector of
     * net external force and torque for the corresponding free DOFs.
     *
     * \param time Time (current_time or new_time) at which external force and torque
     * is required.
     */
    virtual void getNetExternalForceTorque(Vec* F, const double data_time);

    /*!
     * \brief Subtract the mean of constraint force from the background Eulerian
     * grid. This is required for certain cases like periodic steady Stokes.
     *
     * \param L Vec containing the constraint force.
     *
     * \param f_data_idx Patch data index of Eulerian body force.
     *
     * \param scale Factor by which \p L is scaled.
     */
    virtual void subtractMeanConstraintForce(Vec L, int f_data_idx, const double scale = 1.0) = 0;

    /*!
     * \brief Prepare the implementation class for getting the interpolated fluid
     * velocity on the Lagrangian vector \p V.
     *
     * \param V Vector that should contain the interpolated velocity.
     *
     * \param data_time Time at which Eulerian velocity is to be interpolated.
     *
     * \note A default implementation is provided that does nothing.
     */
    virtual void setInterpolatedVelocityVector(Vec V, const double data_time);

    /*!
     * \brief Get the interpolated velocity from the Eulerian grid at the specified time.
     *
     * \param V Vector that should contain the interpolated velocity.
     *
     * \param data_time Time at which Eulerian velocity is to be interpolated.
     *
     * \param scale Scales the velocity vector after interpolating from the
     * Eulerian grid.
     */
    virtual void getInterpolatedVelocity(Vec V, const double data_time, const double scale = 1.0) = 0;

    /*!
     * \brief Compute regularization vector for the mobility problem.
     *
     * \param D Vector containing the regularization for the mobility
     * problem.
     *
     * \param L Vector from which regularization is to be computed.
     */
    virtual void computeMobilityRegularization(Vec D, Vec L, const double scale = 1.0) = 0;

    /*!
     * \brief Get number of rigid structures registered with this class.
     */
    unsigned int getNumberOfRigidStructures() const;

    /*!
     * \brief Get number of nodes associated with the particular structure.
     */
    virtual unsigned int getNumberOfNodes(const unsigned int part) const = 0;

    /*!
     * \brief Set what rigid DOFs need to be solved for this
     * particular structure.
     *
     * \param part The rigid body for which we are setting the free DOFs.
     */
    void setSolveRigidBodyVelocity(const unsigned int part, const FreeRigidDOFVector& solve_rigid_dofs);

    /*!
     * \brief Query what rigid DOFs need to be solved for.
     */
    const FreeRigidDOFVector& getSolveRigidBodyVelocity(const unsigned int part, int& num_free_dofs) const;

    /*!
     * \brief Set the rigid body velocity at the nodal/marker points
     * contained in the Vec \em V.
     *
     * \param part The rigid part for which velocity needs to be set.
     *
     * \param U RDV contains the rigid component of velocities. For
     * two-dimensions the vector contains the values \f$[u,v,\omega_z]\f$
     * and for three-dimensions the vector values are
     * \f$[u,v,w,\omega_x,\omega_y,\omega_z]\f$.
     */
    virtual void setRigidBodyVelocity(const unsigned int part, const RigidDOFVector& U, Vec V) = 0;

    /*!
     * \brief Set the rigid body velocity at the nodal/marker points
     * contained in the Vec \em V.
     *
     * \param part The rigid part for which velocity needs to be set.
     *
     * \param U Vec contains the rigid component of velocities. For
     * two-dimensions the vector contains the values \f$[u,v,\omega_z]\f$
     * and for three-dimensions the vector values are
     * \f$[u,v,w,\omega_x,\omega_y,\omega_z]\f$.
     */
    virtual void setRigidBodyVelocity(const unsigned int part, Vec U, Vec V);

    /*!
     * \brief Set the rigid body velocity at the nodal/marker points
     * contained in the Vec V.
     *
     * \param U Vec that contains the rigid component of velocities for the required
     * components. For two-dimensions each sub Vec contains the values \f$[u,v,\omega_z]\f$
     * and for three-dimensions the vector values are
     * \f$[u,v,w,\omega_x,\omega_y,\omega_z]\f$.
     *
     * \param only_free_dofs Boolean indicating if the rigid body velocity
     * is to be set only for free DOFS for all parts. The corresponding size
     * of U_sub would be \f$ U_{sub} \leq NDIM * (NDIM + 1) / 2 \f$.
     *
     * \param only_imposed_dofs Boolean indicating if the rigid body velocity
     * is to be set only for prescribed kinematics dofs for all parts. The
     * corresponding size of U_sub would be \f$ U_{sub} \leq NDIM * (NDIM + 1) / 2 \f$.
     *
     * \param all_dofs Boolean indicating if the rigid body velocity
     * is to be set for all parts. The corresponding size of U_sub would be
     * \f$ U_{sub} = NDIM * (NDIM + 1) / 2 \f$.
     *
     * \note User is responsible for setting correct number of subvecs in U
     * that corresponds to the particular combination of booleans.
     */
    virtual void setRigidBodyVelocity(Vec U,
                                      Vec V,
                                      const bool only_free_dofs,
                                      const bool only_imposed_dofs,
                                      const bool all_dofs = false);

    /*!
     * \brief Compute total force and torque on the structure.
     *
     * \param L The Lagrange multiplier vector.
     *
     * \param F Vector RDV storing the net generalized force.
     *
     * \param F RDV storing the net generalized force.
     */
    virtual void computeNetRigidGeneralizedForce(const unsigned int part, Vec L, RigidDOFVector& F) = 0;

    /*!
     * \brief Compute total force and torque on the structure.
     *
     * \param part The structure index.
     *
     * \param L The Lagrange multiplier vector.
     *
     * \param F Vec storing the net generalized force.
     */
    virtual void computeNetRigidGeneralizedForce(const unsigned int part, Vec L, Vec F);

    /*!
     * \brief Compute total force and torque on the structure.
     *
     * \param L The Lagrange multiplier vector.
     * \param F Vec storing the net generalized force.
     *
     * \param only_free_dofs Boolean indicating if the net generalized
     * force and torque is to be computed only for free dofs of all bodies.
     *
     * \param only_imposed_dofs Boolean indicating if the net generalized
     * force and torque is to be computed for imposed dofs of all bodies.
     *
     * \param all_dofs Boolean indicating if the net generalized
     * force and torque is to be computed for all dofs of all bodies.
     *
     * \note User is responsible for setting correct number of subvecs in F
     * that corresponds to the particular combination of booleans.
     */
    virtual void computeNetRigidGeneralizedForce(Vec L,
                                                 Vec F,
                                                 const bool only_free_dofs,
                                                 const bool only_imposed_dofs,
                                                 const bool all_dofs = false);

    /*!
     * \brief Get total torque and force on the structure at new_time within
     * the current time interval.
     *
     * \param part The rigid part.
     */
    const RigidDOFVector& getNetRigidGeneralizedForce(const unsigned int part);

    /*!
     * \brief Update the mapping of free DOFs for all structures if they are collected
     * in a global vector.
     */
    void updateFreeDOFsMapping();

    /*!
     * \brief Update the rigid body velocity obtained from the constraint Stokes
     * solver for free-moving case.
     */
    void updateNewRigidBodyVelocity(const unsigned int part, const RigidDOFVector& U);

    /*!
     * \brief Update the rigid body velocity obtained from the constraint Stokes
     * solver for free-moving case.
     */
    void updateNewRigidBodyVelocity(const unsigned int part, Vec U);

    /*!
     * \brief Update the rigid body velocity obtained from the constraint Stokes
     * solver for free-moving case.
     */
    void updateNewRigidBodyVelocity(Vec U,
                                    const bool only_free_dofs,
                                    const bool only_imposed_dofs,
                                    const bool all_dofs = false);

    /*!
     * \brief Copy data from distributed PETSc Vec for specified stucture indices
     * to an array defined on a single processor. A default empty implementation
     * is provided.
     *
     * \param b PETSc Vec to copy from. The Vec stores data for nodal/marker points.
     *
     * \param array Data pointer to copy to.
     *
     * \param struct_ids Vector of structure indices.
     *
     * \param data_depth Depth of the data stored at each Lagrangian node.
     *
     * \param array_rank Rank of the processor on which the array is located.
     *
     * \note The size of \em array is assummed to be sum of nodes
     * of all the structures given in \em struct_ids times the \em data_depth.
     */
    virtual void copyVecToArray(Vec b,
                                double* array,
                                const std::vector<unsigned>& struct_ids,
                                const int data_depth,
                                const int array_rank);

    /*!
     * \brief Copy data from distributed PETSc Vec for specified stucture indices
     * to an array defined on a single processor. A default implementation
     * is provided.
     *
     * \param b PETSc Vec to copy from. The Vec stores only free DOFs of \em all
     * the structures.
     *
     * \param array Data pointer to copy to. It is a linear array of maximum free DOFs
     * of the passed structure IDs.
     *
     * \param struct_ids Vector of structure indices.
     *
     * \param array_rank Rank of the processor on which the array is located.
     *
     * \note The size of \em array is assummed to be sum of maximum number of
     * free degrees of freedom of all the structures given in \em struct_ids. The
     * caller is responsible for allocating and destroying array memory outside of
     * this routine.
     */
    virtual void
    copyFreeDOFsVecToArray(Vec b, double* array, const std::vector<unsigned>& struct_ids, const int array_rank);

    /*!
     * \brief Copy data from array defined on a single processor for specified
     * stucture indices to distributed PETScVec. A default empty implementation
     * is provided.
     *
     * \param b Copy to PETSc Vec.
     *
     * \param array Copy from data pointer.
     *
     * \param struct_ids Vector of structure indices.
     *
     * \param data_depth Depth of the data stored at each Lagrangian node.
     *
     * \param array_rank Rank of the processor on which the array is located.
     *
     * \note The size of \em array is assummed to be sum of nodes
     * of all the structures given in \em struct_ids times the \em data_depth.
     * The caller is responsible for allocating and destroying array memory
     * outside of this routine.
     */
    virtual void copyArrayToVec(Vec b,
                                double* array,
                                const std::vector<unsigned>& struct_ids,
                                const int data_depth,
                                const int array_rank);

    /*!
     * \brief Copy data from array defined on a single processor for specified
     * stucture indices to distributed PETScVec. A default implementation
     * is provided.
     *
     * \param b Copy to PETSc Vec. The Vec stores only free DOFs of \em all
     * the structures.
     *
     * \param array Copy from data pointer. It is a linear array of maximum free DOFs
     * of the passed structure IDs.
     *
     * \param struct_ids Vector of structure indices.
     *
     * \param array_rank Rank of the processor on which the array is located.
     *
     * \note The size of \em array is assummed to be sum of maximum number of
     * free degrees of freedom of all the structures given in \em struct_ids. The
     * caller is responsible for allocating and destroying array memory outside of
     * this routine.
     */
    virtual void
    copyFreeDOFsArrayToVec(Vec b, double* array, const std::vector<unsigned>& struct_ids, const int array_rank);

    /*!
     * \brief Set the DOFs from PETSc Vec \p U to RigidDOFVector \p Ur.
     */
    static void vecToRDV(Vec U, RigidDOFVector& Ur);

    /*!
     * \brief Set the DOFs from RigidDOFVector \p Ur to PETSc Vec \p U.
     */
    static void rdvToVec(const RigidDOFVector& Ur, Vec& U);

    /*!
     * \brief Set the DOFs from Eigen::Vector3d \p U and \p W to RigidDOFVector \p UW.
     */
    static void eigenToRDV(const Eigen::Vector3d& U, const Eigen::Vector3d& W, RigidDOFVector& UW);

    /*!
     * \brief Set the DOFs from RigidDOFVector \p UW to Eigen::Vector3d \p U and \p W.
     */
    static void rdvToEigen(const RigidDOFVector& UW, Eigen::Vector3d& U, Eigen::Vector3d& W);

    /*!
     * \brief Get the rigid body translational velocity at the beginning of
     * the timestep.
     */
    void getCurrentRigidBodyVelocity(const unsigned int part, RigidDOFVector& U);

    /*!
     * \brief Get the rigid body translational velocity at the end of
     * the timestep.
     */
    void getNewRigidBodyVelocity(const unsigned int part, RigidDOFVector& U);

    /*!
     * \brief Get body center of mass at the current time step.
     */
    const Eigen::Vector3d& getCurrentBodyCenterOfMass(const unsigned int part);

    /*!
     * \brief Get body center of mass at half time step.
     */
    const Eigen::Vector3d& getNewBodyCenterOfMass(const unsigned int part);

    /*!
     * \brief Construct dense mobility matrix for the prototypical structures
     * identified by their indices.
     * \note A default empty implementation is provided
     * in this class. The derived class provides the actual implementation.
     *
     * \param mat_name Matrix handle.
     *
     * \param mat_type Mobility matrix type, e.g., RPY, EMPIRICAL, etc.
     *
     * \param mobility_mat Dense sequential mobility matrix. The matrix is
     * stored in column-major(FORTRAN) order.
     * \note Must be allocated prior to entering this routine.
     *
     * \param prototype_struct_ids Indices of the structures as registered
     * with \see IBAMR::IBStrategy class. A combined dense mobility matrix
     * will formed for multiple structures.
     *
     * \param grid_dx NDIM vector of grid spacing of structure level.
     *
     * \param domain_extents NDIM vector of domain length.
     *
     * \param initial_time Boolean to indicate if the mobility matrix is to be
     * generated for the initial position of the structures.
     *
     * \param rho Fluid density
     *
     * \param mu Fluid viscosity.
     *
     * \param scale Scale for improving the conditioning number of dense mobility
     * matrix. The matrix is scaled as:
     * \f$  = \alpha * mobility_mat + \beta * identity_mat. \f$
     *
     * \param managing_rank Rank of the processor managing this dense matrix.
     */
    virtual void constructMobilityMatrix(const std::string& mat_name,
                                         MobilityMatrixType mat_type,
                                         Mat& mobility_mat,
                                         const std::vector<unsigned>& prototype_struct_ids,
                                         const double* grid_dx,
                                         const double* domain_extents,
                                         const bool initial_time,
                                         double rho,
                                         double mu,
                                         const std::pair<double, double>& scale,
                                         double f_periodic_corr,
                                         const int managing_rank);

    /*!
     * \brief Construct a geometric matrix for the prototypical structures
     * identified by their indices. A geometric matrix maps center of mass rigid
     * body velocity to nodal velocities. Geometric matrix is generally used with
     * a dense mobility matrices to construct an associated body-mobility matrix
     * algebrically.
     * \note A default empty implementation is provided in this class. The
     * derived class provides the actual implementation.
     *
     * \param mat_name Matrix handle.
     *
     * \param geometric_mat Dense sequential geometric matrix. The matrix is
     * stored in column-major(FORTRAN) order.
     * \note Must be allocated prior to entering this routine.
     *
     * \param prototype_struct_ids Indices of the structures as registered
     * with \see IBAMR::IBStrategy class. A combined block-diagonal geometric
     * matrix will be formed for multiple structures.
     *
     * \param initial_time Boolean to indicate if the corresponding geometric matrix
     * is to be generated for the initial position of the structures.
     *
     * \param managing_rank Rank of the processor managing this dense matrix.
     */
    virtual void constructGeometricMatrix(const std::string& mat_name,
                                          Mat& geometric_mat,
                                          const std::vector<unsigned>& prototype_struct_ids,
                                          const bool initial_time,
                                          const int managing_rank);

    /*!
     * \brief Rotate vector using rotation matrix to/from the reference frame
     * of the structures (which is at the initial time of the simulation).
     *
     * \param array Raw data pointer containing the vector enteries.
     *
     * \param struct_ids Structure ID indices.
     *
     * \param use_transpose Use transpose of rotation matrix to rotate the vector.
     * \note Transpose of rotation matrix is its inverse and it takes the vector back
     * to its reference frame.
     *
     * \param managing_rank Rank of the processor managing the matrix.
     *
     * \param depth Depth of the data array components.
     *
     */
    virtual void rotateArray(double* array,
                             const std::vector<unsigned>& struct_ids,
                             const bool use_transpose,
                             const int managing_rank,
                             const int depth);

    /////////////////////////////// PROTECTED ////////////////////////////////////
protected:
    /*!
     * \brief Fill the rotation matrix.
     * \param q_old Previous applied quaternion.
     * \param q_new New quaternion to set.
     * \param rot_mat Matrix to set.
     * \param dt Time interval of rotation.
     */
    void setRotationMatrix(const std::vector<Eigen::Vector3d>& rot_vel,
                           const std::vector<Eigen::Quaterniond>& q_old,
                           std::vector<Eigen::Quaterniond>& q_new,
                           std::vector<Eigen::Matrix3d>& rot_mat,
                           const double dt);

    // Number of rigid parts.
    unsigned int d_num_rigid_parts;

    // Center of mass.
    std::vector<Eigen::Vector3d> d_center_of_mass_initial, d_center_of_mass_current, d_center_of_mass_half,
        d_center_of_mass_new;

    // Quaternion of the body.
    std::vector<Eigen::Quaterniond> d_quaternion_current, d_quaternion_half, d_quaternion_new;

    // Whether to solve for rigid body velocity.
    std::vector<FRDV> d_solve_rigid_vel;

    // PETSc wrappers for free rigid body velocities and external force/torque.
    Vec d_U, d_F;

    // Mapping of free DOFs in the global vector.
    std::vector<std::pair<int, int> > d_free_dofs_map;
    bool d_free_dofs_map_updated;

    // Rigid body velocity of the structures.
    std::vector<Eigen::Vector3d> d_trans_vel_current, d_trans_vel_half, d_trans_vel_new;
    std::vector<Eigen::Vector3d> d_rot_vel_current, d_rot_vel_half, d_rot_vel_new;

    // Net rigid generalized force.
    std::vector<RigidDOFVector> d_net_rigid_generalized_force;

    /////////////////////////////// PRIVATE //////////////////////////////////////
private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CIBStrategy(const CIBStrategy& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CIBStrategy& operator=(const CIBStrategy& that);

}; // CIBStrategy

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // ifndef included_IBAMR_CIBStrategy
