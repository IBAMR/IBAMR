// Filename: CIBStrategy.h
// Created on 8 Nov 2014 by Amneet Bhalla
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
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

#ifndef included_CIBStrategy
#define included_CIBStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <vector>

#include "tbox/DescribedClass.h"
#include "petscvec.h"
#include "Eigen/Core"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
	
typedef Eigen::Matrix<double, NDIM*(NDIM+1)/2, 1> RigidDOFVector;
typedef RigidDOFVector RDV;
	
/*!
 * \brief Class CIBStrategy is a lightweight abstract strategy class which 
 * provides support for constraint based IB methods for rigid bodies.
 */
class CIBStrategy : public virtual SAMRAI::tbox::DescribedClass
{
//////////////////////////////////////////////////////////////////////////////
public:
	/*!
	 * \brief Constructor of the class.
	 */
	CIBStrategy(
		const unsigned int parts);

	/*!
     * \brief Destructor of the class.
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
	virtual void
	setConstraintForce(
		Vec L,
		double data_time,
		double scale = 1.0) = 0;
	
	/*!
	 * \brief Get the constraint rigid body force at the specified time within
	 * the current time interval.
	 *
	 * \param time Time (current_time or new_time) at which constraint force is 
	 * required.
	 */
	virtual void
	getConstraintForce(
		Vec* L,
		const double time) = 0;
	
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
	virtual void
	subtractMeanConstraintForce(
		Vec L,
		int f_data_idx,
		const double scale = 1.0) = 0;
	
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
	virtual void
	setInterpolatedVelocityVector(
		Vec V,
		double data_time);
	
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
	virtual void
	getInterpolatedVelocity(
		Vec V,
		double data_time,
		double scale = 1.0) = 0;
	
	/*!
	 * \brief Get number of rigid structures registered with this class.
	 */
	unsigned int
	getNumberOfRigidStructures() const;
	
	/*!
	 * \brief Set if the rigid velocity is to be solved for this
	 * particular structure.
	 *
	 * \param part The rigid body for which boolean is to be set.
	 */
	void
	setSolveRigidBodyVelocity(
		const unsigned int part,
		const bool solve_rigid_vel);
	
	/*!
	 * \brief Query if the rigid velocity is to be solved for.
	 *
	 * \return true if rigid velocity is to be solved for.
	 */
	bool
	getSolveRigidBodyVelocity(
		const unsigned int part) const;
	
	/*!
	 * \brief Set the rigid body velocity at the nodal/marker points
	 * contained in the Vec V.
	 *
	 * \param U contains the rigid component of velocities. For
	 * two-dimensions the vector contains the values \f$[u,v,\omega_z]\f$
	 * and for three-dimensions the vector values are
	 * \f$[u,v,w,\omega_x,\omega_y,\omega_z]\f$.
	 */
	virtual void
	setRigidBodyVelocity(
		const unsigned int part,
		const RigidDOFVector& U,
		Vec V) = 0;
	
	virtual void
	setRigidBodyVelocity(
		const unsigned int part,
		Vec U,
		Vec V);
	
	/*!
	 * \brief Set the rigid body velocity at the nodal/marker points
	 * contained in the Vec V.
	 *
	 * \param U PetscMultiVec that contains the rigid component of velocities.
	 * For two-dimensions each sub Vec contains the values \f$[u,v,\omega_z]\f$
	 * and for three-dimensions the vector values are
	 * \f$[u,v,w,\omega_x,\omega_y,\omega_z]\f$.
	 *
	 * \param only_free_parts Boolean indicating if the rigid body velocity
	 * is to be set only for free moving bodies.
	 *
	 * \param only_imposed_parts Boolean indicating if the rigid body velocity
	 * is to be set only for prescribed kinematics bodies.
	 *
	 * \param all_parts Boolean indicating if the rigid body velocity
	 * is to be set for all bodies.
	 *
	 * \note User is responsible for setting correct number of subvecs in U
	 * that corresponds to the particular combination of booleans.
	 */
	virtual void
	setRigidBodyVelocity(
		Vec U,
		Vec V,
		const bool only_free_parts,
		const bool only_imposed_parts,
		const bool all_parts=false);
	
	/*!
	 * \brief Compute total force and torque on the structure.
	 *
	 * \param part The structure index.
	 *
	 * \param L The Lagrange multiplier vector.
	 *
	 * \param F RDV storing the net generalized force.
	 */
	virtual void
	computeNetRigidGeneralizedForce(
		const unsigned int part,
		Vec L,
		RigidDOFVector& F) = 0;
	
	/*!
	 * \brief Compute total force and torque on the structure.
	 *
	 * \param part The structure index.
	 *
	 * \param L The Lagrange multiplier vector.
	 *
	 * \param F Vec storing the net generalized force.
	 */
	virtual void
	computeNetRigidGeneralizedForce(
		const unsigned int part,
		Vec L,
		Vec F);
	
	/*!
	 * \brief Compute total force and torque on the structure.
	 *
	 * \param L The Lagrange multiplier vector.
	 * \param F PetscMultiVec storing the net generalized force.
	 *
	 * \param only_free_parts Boolean indicating if the net generalized
	 * force and torque is to be computed only for free moving bodies.
	 *
	 * \param only_imposed_parts Boolean indicating if the net generalized
	 * force and torque is to be computed only for prescribed kinematics bodies.
	 *
	 * \param all_parts Boolean indicating if the net generalized
	 * force and torque is to be computed for all bodies.
	 *
	 * \note User is responsible for setting correct number of subvecs in F
	 * that corresponds to the particular combination of booleans.
	 */
	virtual void
	computeNetRigidGeneralizedForce(
		Vec L,
		Vec F,
		const bool only_free_parts,
		const bool only_imposed_parts,
		const bool all_parts=false);
	
	/*!
	 * \brief Update the rigid body velocity obtained from the constraint Stokes
	 * solver for free-moving case.
	 */
	void
	updateNewRigidBodyVelocity(
		const unsigned int part,
		const RigidDOFVector& U);
	
	void
	updateNewRigidBodyVelocity(
		const unsigned int part,
		Vec U);
	
	/*!
	 * \brief Set the DOFs from PETSc Vec \p U to RigidDOFVector \p Ur.
	 */
	static void
	vecToRDV(
		Vec U,
		RigidDOFVector& Ur);
	
	/*!
	 * \brief Set the DOFs from RigidDOFVector \p Ur to PETSc Vec \p U.
	 */
	static void
	rdvToVec(
		const RigidDOFVector& Ur,
		Vec	U);
	
	/*!
	 * \brief Get the rigid body translational velocity at the beginning of
	 * the timestep.
	 */
	void
	getCurrentRigidBodyVelocity(
		const unsigned int part,
		RigidDOFVector& U);
	
	/*!
	 * \brief Get the rigid body translational velocity at the end of
	 * the timestep.
	 */
	void
	getNewRigidBodyVelocity(
		const unsigned int part,
		RigidDOFVector& U);
	
/////////////////////////////// PROTECTED ////////////////////////////////////
	
protected:
	/*!
	 * \brief Number of rigid parts.
	 */
	unsigned int d_num_rigid_parts;
	
	/*!
	 * Center of mass and moment of inertia.
	 */
	std::vector<Eigen::Vector3d> d_center_of_mass_current, d_center_of_mass_half;
	std::vector<Eigen::Matrix3d> d_moment_of_inertia_current, d_moment_of_inertia_half;
	
	/*!
	 * Whether to solve for rigid body velocity.
	 */
	std::vector<bool> d_solve_rigid_vel;
	
	/*!
	 * Rigid body velocity of the structures.
	 */
	std::vector<Eigen::Vector3d> d_trans_vel_current, d_trans_vel_half, d_trans_vel_new;
	std::vector<Eigen::Vector3d> d_rot_vel_current, d_rot_vel_half, d_rot_vel_new;
	
//////////////////////////////////////////////////////////////////////////////
private:
	/*!
	 * \brief Copy constructor.
	 *
	 * \note This constructor is not implemented and should not be used.
	 *
	 * \param from The value to copy to this object.
	 */
	CIBStrategy(
	    const CIBStrategy& from);
	
	/*!
	 * \brief Assignment operator.
	 *
	 * \note This operator is not implemented and should not be used.
	 *
	 * \param that The value to assign to this object.
	 *
	 * \return A reference to this object.
	 */
	CIBStrategy& operator=(
	    const CIBStrategy& that);
	
};// CIBStrategy
	
}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //ifndef included_CIBStrategy
