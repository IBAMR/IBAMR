// Filename: CIBStrategy.cpp
// Created on 9 Nov 2014 by Amneet Bhalla
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


/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/CIBStrategy.h"
#include "ibamr/namespaces.h"
#include "ibtk/PETScMultiVec.h"
#include "tbox/Utilities.h"
#include "tbox/SAMRAI_MPI.h"


/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{

/////////////////////////////// STATIC ///////////////////////////////////////
	
	
/////////////////////////////// PUBLIC ///////////////////////////////////////
	
CIBStrategy::CIBStrategy(
	const unsigned int parts)
	: d_num_rigid_parts(parts)
{
	// Resize some arrays
	d_center_of_mass_current.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());
	d_center_of_mass_half.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());
	d_moment_of_inertia_current.resize(d_num_rigid_parts, Eigen::Matrix3d::Zero());
	d_moment_of_inertia_half.resize(d_num_rigid_parts, Eigen::Matrix3d::Zero());
	d_trans_vel_current.resize(d_num_rigid_parts,Eigen::Vector3d::Zero());
	d_trans_vel_half.resize(d_num_rigid_parts,Eigen::Vector3d::Zero());
	d_trans_vel_new.resize(d_num_rigid_parts,Eigen::Vector3d::Zero());
	d_rot_vel_current.resize(d_num_rigid_parts,Eigen::Vector3d::Zero());
	d_rot_vel_half.resize(d_num_rigid_parts,Eigen::Vector3d::Zero());
	d_rot_vel_new.resize(d_num_rigid_parts,Eigen::Vector3d::Zero());
	d_solve_rigid_vel.resize(d_num_rigid_parts,false);
	
	return;
}// CIBStrategy
	
CIBStrategy::~CIBStrategy()
{
	//intentionally left blank
	return;
}

void
CIBStrategy::setInterpolatedVelocityVector(
	Vec /*V*/,
	double /*data_time*/)
{
	// intentionally left blank
	return;
}// setInterpolatedVelocityVector
	
int
CIBStrategy::getNumberOfRigidStructures() const
{
	return d_num_rigid_parts;
		
}// getNumberOfRigidStructures
	
void
CIBStrategy::setSolveRigidBodyVelocity(
	const unsigned int part,
	const bool solve_rigid_vel)
{
	d_solve_rigid_vel[part] = solve_rigid_vel;
	return;
	
}// setSolveRigidBodyVelocity
	
bool
CIBStrategy::getSolveRigidBodyVelocity(
	const unsigned int part) const
{
	return d_solve_rigid_vel[part];
		
}// getSolveRigidBodyVelocity
	
void
CIBStrategy::setRigidBodyVelocity(
	const unsigned int part,
	Vec U,
	Vec V)
{
	RigidDOFVector Ur;
	vecToRDV(U, Ur);
	setRigidBodyVelocity(part, Ur, V);
		
	return;
}// setRigidBodyVelocity
	
void
CIBStrategy::setRigidBodyVelocity(
	Vec U,
	Vec V,
	const bool only_free_parts,
	const bool only_imposed_parts,
	const bool all_parts)
{
	bool has_free_parts = false, has_imposed_parts = false;
	for (unsigned part = 0; part < d_num_rigid_parts; ++part)
	{
		has_free_parts    = has_free_parts    ||  d_solve_rigid_vel[part];
		has_imposed_parts = has_imposed_parts || !d_solve_rigid_vel[part];
	}
		
	if (has_free_parts && only_free_parts)
	{
		Vec* vU;
		VecMultiVecGetSubVecs(U, &vU);
		for (unsigned part = 0, k = 0; part < d_num_rigid_parts; ++part)
		{
			if (!d_solve_rigid_vel[part]) continue;
			setRigidBodyVelocity(part, vU[k], V);
			++k;
		}
	}
	else if (has_imposed_parts && only_imposed_parts)
	{
		Vec* vU;
		VecMultiVecGetSubVecs(U, &vU);
		for (unsigned part = 0, k = 0; part < d_num_rigid_parts; ++part)
		{
			if (d_solve_rigid_vel[part]) continue;
			setRigidBodyVelocity(part, vU[k], V);
			++k;
		}
	}
	else if (all_parts)
	{
		Vec* vU;
		VecMultiVecGetSubVecs(U, &vU);
		for (unsigned part = 0; part < d_num_rigid_parts; ++part)
		{
			setRigidBodyVelocity(part, vU[part], V);
		}
	}
		
	return;
}// setRigidBodyVelocity
	
void
CIBStrategy::computeNetRigidGeneralizedForce(
	const unsigned int part,
	Vec L,
	Vec F)
{
	RigidDOFVector f;
	computeNetRigidGeneralizedForce(part, L, f);
	rdvToVec(f, F);
		
	return;
}// computeNetRigidGeneralizedForce
	
void
CIBStrategy::computeNetRigidGeneralizedForce(
	Vec L,
    Vec F,
	const bool only_free_parts,
	const bool only_imposed_parts,
	const bool all_parts)
{
	bool has_free_parts = false, has_imposed_parts = false;
	for (unsigned part = 0; part < d_num_rigid_parts; ++part)
	{
		has_free_parts    = has_free_parts    ||  d_solve_rigid_vel[part];
		has_imposed_parts = has_imposed_parts || !d_solve_rigid_vel[part];
	}
		
	if (has_free_parts && only_free_parts)
	{
		Vec* vF;
		VecMultiVecGetSubVecs(F, &vF);
		for (unsigned part = 0, k = 0; part < d_num_rigid_parts; ++part)
		{
			if (!d_solve_rigid_vel[part]) continue;
			computeNetRigidGeneralizedForce(part, L, vF[k]);
			++k;
		}
	}
	else if (has_imposed_parts && only_imposed_parts)
	{
		Vec* vF;
		VecMultiVecGetSubVecs(F, &vF);
		for (unsigned part = 0, k = 0; part < d_num_rigid_parts; ++part)
		{
			if (d_solve_rigid_vel[part]) continue;
			computeNetRigidGeneralizedForce(part, L, vF[k]);
			++k;
		}
	}
	else if (all_parts)
	{
		Vec* vF;
		VecMultiVecGetSubVecs(F, &vF);
		for (unsigned part = 0; part < d_num_rigid_parts; ++part)
		{
			computeNetRigidGeneralizedForce(part, L, vF[part]);
		}
	}
		
	return;
}// computeNetRigidGeneralizedForce
	
void
CIBStrategy::updateNewRigidBodyVelocity(
	const unsigned int part,
    const RigidDOFVector& U)
{
	TBOX_ASSERT(part < d_num_rigid_parts);
	TBOX_ASSERT(d_solve_rigid_vel[part]);
		
	for (unsigned d = 0; d < NDIM; ++d)
	{
		d_trans_vel_new[part][d] = U[d];
	}
		
#if (NDIM == 2)
		d_rot_vel_new[part].setZero();
		d_rot_vel_new[part][NDIM] = U[NDIM];
#elif (NDIM == 3)
		for (unsigned d = 0; d < NDIM; ++d)
		{
			d_trans_vel_new[part][d] = U[NDIM+d];
		}
#endif
		
		d_trans_vel_half[part] = 0.5*(d_trans_vel_current[part] + d_trans_vel_new[part]);
		d_rot_vel_half[part]   = 0.5*(d_rot_vel_current[part]   + d_rot_vel_new[part]);
		
		return;
}// updateNewRigidBodyVelocity
	
void
CIBStrategy::updateNewRigidBodyVelocity(
	const unsigned int part,
	Vec U)
{
	RigidDOFVector Ur;
	vecToRDV(U, Ur);
	updateNewRigidBodyVelocity(part, Ur);
		
	return;
}// updateNewRigidBodyVelocity
	
void
CIBStrategy::vecToRDV(
	Vec U,
	RigidDOFVector& Ur)
{
	// Extract the underlying array.
	PetscScalar* a;
	PetscInt s;
	VecGetArray(U,&a);
	VecGetSize (U,&s);
		
	// Fill in the required vector.
	int rank = SAMRAI_MPI::getRank();
	if (!rank)
	{
		for (int i = 0; i < s; ++i)
			Ur[i] = a[i];
	}
	SAMRAI_MPI::sumReduction(&Ur[0],s);
	VecRestoreArray(U,&a);
		return;
}// vecToRDV
	
void
CIBStrategy::rdvToVec(
	const RigidDOFVector& Ur,
	Vec U)
{
	static const int s = NDIM*(NDIM+1)/2;
	std::vector<int> idx(s);
	for (int i = 0; i < s; ++i) idx[i] = i;
	VecSetValues(U, s, &idx[0], &Ur[0], INSERT_VALUES);
	VecAssemblyBegin(U);
	VecAssemblyEnd(U);
		
	return;
}// rdvToVec
	
void
CIBStrategy::getCurrentRigidBodyVelocity(
	const unsigned int part,
	RigidDOFVector& U)
{
	for (unsigned d = 0; d < NDIM; ++d)
	{
		U[d] = d_trans_vel_current[part][d];
	}
		
#if (NDIM == 2)
	U[NDIM] = d_rot_vel_current[part][NDIM];
#elif (NDIM == 3)
	for (unsigned d = 0; d < NDIM; ++d)
	{
		U[NDIM+d] = d_rot_vel_current[part][d];
	}
#endif
		
	return;
}// getCurrentRigidBodyVelocity
	
void
CIBStrategy::getNewRigidBodyVelocity(
	const unsigned int part,
	RigidDOFVector& U)
{
	for (unsigned d = 0; d < NDIM; ++d)
	{
		U[d] = d_trans_vel_new[part][d];
	}
		
#if (NDIM == 2)
	U[NDIM] = d_rot_vel_new[part][NDIM];
#elif (NDIM == 3)
	for (unsigned d = 0; d < NDIM; ++d)
	{
		U[NDIM+d] = d_rot_vel_new[part][d];
	}
#endif
		
	return;
}// getNewRigidBodyVelocity

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR