// Modified by Baky Aug 2015

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

CIBStrategy::CIBStrategy(const unsigned int parts) : d_num_rigid_parts(parts)
{
    // Resize some arrays
    d_center_of_mass_current.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());
    d_center_of_mass_half.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());
    d_quaternion_current.resize(d_num_rigid_parts, Eigen::Quaterniond::Identity());
    d_quaternion_half.resize(d_num_rigid_parts, Eigen::Quaterniond::Identity());
    d_trans_vel_current.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());
    d_trans_vel_half.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());
    d_trans_vel_new.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());
    d_rot_vel_current.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());
    d_rot_vel_half.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());
    d_rot_vel_new.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());
    d_solve_rigid_vel.resize(d_num_rigid_parts);
    d_net_rigid_generalized_force.resize(d_num_rigid_parts);

    return;
} // CIBStrategy

CIBStrategy::~CIBStrategy()
{
    // intentionally left blank
    return;
}

void CIBStrategy::setInterpolatedVelocityVector(Vec /*V*/, const double /*data_time*/)
{
    // intentionally left blank
    return;
} // setInterpolatedVelocityVector

unsigned int CIBStrategy::getNumberOfRigidStructures() const
{
    return d_num_rigid_parts;

} // getNumberOfRigidStructures

void CIBStrategy::setSolveRigidBodyVelocity(const unsigned int part, const FreeRigidDOFVector& solve_rigid_vel)
{
    d_solve_rigid_vel[part] = solve_rigid_vel;
    return;

} // setSolveRigidBodyVelocity

const FreeRigidDOFVector& CIBStrategy::getSolveRigidBodyVelocity(const unsigned int part, int& num_free_dofs) const
{
    num_free_dofs = 0;
    for (int i = 0; i < s_max_free_dofs; ++i)
    {
        if (d_solve_rigid_vel[part][i])
        {
            ++num_free_dofs;
        }
    }

    return d_solve_rigid_vel[part];
} // getSolveRigidBodyVelocity

void CIBStrategy::setRigidBodyVelocity(unsigned int part, const RigidDOFVector& U, Vec V)
{
    std::vector<bool> skip_comp;
    skip_comp.resize(d_num_rigid_parts);
    fill(skip_comp.begin(), skip_comp.end(), true);   
    skip_comp[part] = false;
    std::vector<RigidDOFVector> temp_U;
    temp_U.resize(d_num_rigid_parts);
    temp_U[part]= U;

    //collective call over all nodes
    setRigidBodyVelocity(temp_U, V, skip_comp);
}
void CIBStrategy::setRigidBodyVelocity(const unsigned int part, Vec U, Vec V)
{
    RigidDOFVector Ur;
    vecToRDV(U, Ur);
    setRigidBodyVelocity(part, Ur, V);

    return;
} // setRigidBodyVelocity

void CIBStrategy::setRigidBodyVelocity(Vec U,
                                       Vec V,
                                       const bool only_free_dofs,
                                       const bool only_imposed_dofs,
                                       const bool all_dofs)
{
    std::vector<bool> skip_comp;
    skip_comp.resize(d_num_rigid_parts);
    std::vector<RigidDOFVector> UDOF;
    UDOF.resize(d_num_rigid_parts);
    for (unsigned part = 0; part < d_num_rigid_parts; ++part)
    {   
	UDOF[part].setZero();
    }

    if (only_free_dofs || only_imposed_dofs)
    {
        for (unsigned part = 0, comp_part = 0; part < d_num_rigid_parts; ++part)
        {
            int num_free_dofs=0;
	    for (int k = 0; k < s_max_free_dofs; ++k) num_free_dofs +=d_solve_rigid_vel[part][k];
	    RigidDOFVector& U_part = UDOF[part];
	    
	    Vec U_sub;
	    VecMultiVecGetSubVec(U, comp_part, &U_sub);
	    
	    PetscInt s;
	    PetscScalar* a = NULL;
	    VecGetArray(U_sub, &a);
	    VecGetSize(U_sub, &s);
	    
	    std::vector<double> a_vec(s, 0.0);
	    if (a != NULL)
	    {
		std::copy(&a[0], &a[s], &a_vec[0]);
	    }
	    SAMRAI_MPI::sumReduction(&a_vec[0], s);
	    VecRestoreArray(U_sub, &a);
	    
	    for (int k = 0, p = 0; k < s_max_free_dofs; ++k)
	    {
		if (only_free_dofs && num_free_dofs)
		{
		    if (d_solve_rigid_vel[part][k])
		    {
			U_part[k] = a_vec[p];
			++p;
		    }
		    skip_comp[part]=false;
		} 
		else if (only_imposed_dofs && (num_free_dofs < s_max_free_dofs))
		{
		    if (!d_solve_rigid_vel[part][k])
		    {
			U_part[k] = a_vec[p];
			++p;
		    }
		    skip_comp[part]=false;
		}
		else
		{
		    skip_comp[part]=true;
		}
	    }
	    ++comp_part;
	}
    } 
    else if (all_dofs)
    {
        Vec* vU;
        VecMultiVecGetSubVecs(U, &vU);
        for (unsigned part = 0; part < d_num_rigid_parts; ++part)
        {
	    RigidDOFVector& U_part = UDOF[part];
	    vecToRDV(vU[part], U_part);
	}
	fill(skip_comp.begin(), skip_comp.end(), false);   
    }  
    else
    {
	TBOX_ERROR("CIBStrategy::setRigidBodyVelocity() incorrect parameters in the function call" << std::endl);
    }
   
    //collective call over all nodes
    setRigidBodyVelocity(UDOF, V, skip_comp);

    return;
} // setRigidBodyVelocity

void CIBStrategy::setRigidBodyDeformationVelocity(const unsigned int /*part*/, Vec /*y*/)
{
    // intentionally left blank.
    return;
}

void CIBStrategy::setRigidBodyDeformationVelocity(Vec W)
{
    // Vec* vW;
    // VecMultiVecGetSubVecs(W, &vW);

    for (unsigned part = 0; part < d_num_rigid_parts; ++part)
    {  
	setRigidBodyDeformationVelocity(part, W);
    }
    return;
}
void CIBStrategy::computeNetRigidGeneralizedForce(const unsigned int part, Vec L,  RigidDOFVector& F)
{
    std::vector<bool> skip_comp;
    skip_comp.resize(d_num_rigid_parts);
    fill(skip_comp.begin(), skip_comp.end(), true);   
    skip_comp[part] = false;

    //collective call over all nodes
    computeNetRigidGeneralizedForce(L, d_net_rigid_generalized_force, skip_comp);
    F=d_net_rigid_generalized_force[part];

    return;
} // computeNetRigidGeneralizedForce

void CIBStrategy::computeNetRigidGeneralizedForce(const unsigned int part, Vec L, Vec F)
{
    RigidDOFVector f;
    computeNetRigidGeneralizedForce(part, L, f);
    rdvToVec(f, F);

    return;
} // computeNetRigidGeneralizedForce

void CIBStrategy::computeNetRigidGeneralizedForce(Vec L,
                                                  Vec F,
                                                  const bool only_free_dofs,
                                                  const bool only_imposed_dofs,
                                                  const bool all_dofs)
{
    std::vector<bool> skip_comp;
    skip_comp.resize(d_num_rigid_parts);

    if (only_free_dofs || only_imposed_dofs)
    {
	for (unsigned part = 0; part < d_num_rigid_parts; ++part)
	{
	    int num_free_dofs=0;
	    for (int k = 0; k < s_max_free_dofs; ++k) num_free_dofs +=d_solve_rigid_vel[part][k];
	    
	    if (only_free_dofs && num_free_dofs) 
	    {
		skip_comp[part] = false;
		d_net_rigid_generalized_force[part].setZero();
	    }
	    else if (only_imposed_dofs && (num_free_dofs < s_max_free_dofs))
	    { 
		skip_comp[part] = false;
	    }
	    else
	    {
		skip_comp[part] = true;
	    }
	}
    }
    else if (all_dofs)
    {     
	fill(skip_comp.begin(), skip_comp.end(), false);   
    }
    else
    {
	TBOX_ERROR("CIBStrategy::computeNetRigidGeneralizedForce() incorrect parameters in the function call" << std::endl);
    }

    //collective call over all nodes
    computeNetRigidGeneralizedForce(L, d_net_rigid_generalized_force, skip_comp);

    //flll the vector
    if (only_free_dofs || only_imposed_dofs)
    {
	for (unsigned part = 0, comp_part = 0; part < d_num_rigid_parts; ++part)
	{ 
	    if (skip_comp[part]) continue;
	    
	    RigidDOFVector& F_part = d_net_rigid_generalized_force[part];
	    
	    Vec F_sub;
	    VecMultiVecGetSubVec(F, comp_part, &F_sub);
	    
	    PetscInt s;
	    PetscScalar* a = NULL;
	    VecGetArray(F_sub, &a);
	    VecGetSize(F_sub, &s);
	    
	    if (a != NULL)
	    {
		if (only_free_dofs)
		{
		    for (int k = 0, p = 0; k < s_max_free_dofs; ++k)
		    {
			if (d_solve_rigid_vel[part][k])
			{
			    a[p] = F_part[k];
			    ++p;
			}
		    }
				
		}else if (only_imposed_dofs)
		{
		    for (int k = 0, p = 0; k < s_max_free_dofs; ++k)
		    {
			if (!d_solve_rigid_vel[part][k])
			{
			    a[p] = F_part[k];
			    ++p;
			}
		    }
		}
	    }
	    VecRestoreArray(F_sub, &a);
	    ++comp_part;
	}
    } else if (all_dofs)
    {        
	Vec* vF;
	VecMultiVecGetSubVecs(F, &vF);
	for (unsigned part = 0; part < d_num_rigid_parts; ++part)
	{
	    rdvToVec(d_net_rigid_generalized_force[part], vF[part]);
	}
	    
    }

    return;
} // computeNetRigidGeneralizedForce

const RigidDOFVector& CIBStrategy::getNetRigidGeneralizedForce(const unsigned int part)
{
    return d_net_rigid_generalized_force[part];
} // getNetRigidGeneralizedForce

void CIBStrategy::updateNewRigidBodyVelocity(const unsigned int part, const RigidDOFVector& U)
{
    TBOX_ASSERT(part < d_num_rigid_parts);
    rdvToEigen(U, d_trans_vel_new[part], d_rot_vel_new[part]);

    d_trans_vel_half[part] = 0.5 * (d_trans_vel_current[part] + d_trans_vel_new[part]);
    d_rot_vel_half[part] = 0.5 * (d_rot_vel_current[part] + d_rot_vel_new[part]);

    return;
} // updateNewRigidBodyVelocity

void CIBStrategy::updateNewRigidBodyVelocity(const unsigned int part, Vec U)
{
    RigidDOFVector Ur;
    vecToRDV(U, Ur);
    updateNewRigidBodyVelocity(part, Ur);

    return;
} // updateNewRigidBodyVelocity

void CIBStrategy::copyVecToArray(Vec /*b*/,
                                 double* /*array*/,
                                 const std::vector<unsigned>& /*struct_ids*/,
                                 const int /*data_depth*/,
                                 const int /*array_rank*/)
{
    // intentionally left-blank
    return;
} // copyVecToArray

void CIBStrategy::copyArrayToVec(Vec /*b*/,
                                 double* /*array*/,
                                 const std::vector<unsigned>& /*struct_ids*/,
                                 const int /*data_depth*/,
                                 const int /*array_rank*/)
{
    // intentionally left-blank
    return;
} // copyArrayToVec

void CIBStrategy::copyAllArrayToVec(Vec /*b*/,
				    double* /*array*/,
				    const std::vector<unsigned>& /*rhs_struct_ids*/,
				    const int /*data_depth*/,
				    const int /*array_rank*/)
{
    // intentionally left-blank
    return;
} // copyArrayToVec

void CIBStrategy::copyAllVecToArray(Vec /*b*/,
				    double* /*array*/,
				    const std::vector<unsigned>& /*rhs_struct_ids*/,
				    const int /*data_depth*/,
				    const int /*array_rank*/)
{
    // intentionally left-blank
    return;
} // copyVecToArray

void CIBStrategy::vecToRDV(Vec U, RigidDOFVector& Ur)
{
    // Extract the underlying array.
    PetscScalar* a = NULL;
    PetscInt s;
    VecGetArray(U, &a);
    VecGetSize(U, &s);

#if !defined(NDEBUG)
    TBOX_ASSERT(s == s_max_free_dofs);
#endif

    // Fill in the required vector.
    if (a != NULL)
    {
        std::copy(&a[0], &a[s], &Ur[0]);
    }
    else
    {
        Ur.setZero();
    }

    SAMRAI_MPI::sumReduction(&Ur[0], s);
    VecRestoreArray(U, &a);
    return;
} // vecToRDV

void CIBStrategy::rdvToVec(const RigidDOFVector& Ur, Vec& U)
{
    if (U == NULL)
    {
        PetscInt n = 0, N = s_max_free_dofs;
        if (!SAMRAI_MPI::getRank()) n = N;
        VecCreateMPI(PETSC_COMM_WORLD, n, N, &U);
    }

    std::vector<int> idx(s_max_free_dofs);
    for (int i = 0; i < s_max_free_dofs; ++i) idx[i] = i;
    VecSetValues(U, s_max_free_dofs, &idx[0], &Ur[0], INSERT_VALUES);
    VecAssemblyBegin(U);
    VecAssemblyEnd(U);

    return;
} // rdvToVec

void CIBStrategy::eigenToRDV(const Eigen::Vector3d& U, const Eigen::Vector3d& W, RigidDOFVector& UW)
{
    for (unsigned d = 0; d < NDIM; ++d)
    {
        UW[d] = U[d];
    }

#if (NDIM == 2)
    UW[2] = W[2];
#elif(NDIM == 3)
    for (unsigned d = 0; d < NDIM; ++d)
    {
        UW[3 + d] = W[d];
    }
#endif

    return;
} // eigenToRDV

void CIBStrategy::rdvToEigen(const RigidDOFVector& UW, Eigen::Vector3d& U, Eigen::Vector3d& W)
{
    for (unsigned d = 0; d < NDIM; ++d)
    {
        U[d] = UW[d];
    }

#if (NDIM == 2)
    W.setZero();
    W[2] = UW[2];
#elif(NDIM == 3)
    for (unsigned d = 0; d < NDIM; ++d)
    {
        W[d] = UW[3 + d];
    }
#endif

    return;
} // rdvToEigen

void CIBStrategy::getCurrentRigidBodyVelocity(const unsigned int part, RigidDOFVector& U)
{

    eigenToRDV(d_trans_vel_current[part], d_rot_vel_current[part], U);

    return;
} // getCurrentRigidBodyVelocity

void CIBStrategy::getNewRigidBodyVelocity(const unsigned int part, RigidDOFVector& U)
{

    eigenToRDV(d_trans_vel_new[part], d_rot_vel_new[part], U);

    return;
} // getNewRigidBodyVelocity

void CIBStrategy::constructMobilityMatrix(std::map<std::string, double*>&  /*mat_map*/,
					  std::map<std::string, MobilityMatrixType>& /*mat_type_map*/,
					  std::map<std::string, std::vector<unsigned> >& /*mat_prototype_id_map*/,
					  std::map<std::string, unsigned int>& /*managed_mat_nodes_map*/,
					  std::map<std::string, std::pair<double, double> >& /*mat_scale_map*/,
					  std::map<std::string, int>& /*managed_mat_proc_map*/,
					  const double* /*grid_dx*/,
					  const double* /*domain_extents*/,
					  const bool /*initial_time*/,
					  double /*rho*/,
					  double /*mu*/,
					  double /*f_periodic_corr*/)

{
    // intentionally left blank.
    return;
} // constructMobilityMatrix

void CIBStrategy::constructKinematicMatrix(double* /*kinematic_mat*/,
					   const std::vector<unsigned>& /*prototype_struct_ids*/,
					   const bool /*initial_time*/,
					   const int /*managing_rank*/)
{
    // intentionally left blank.
    return;
} // constructBodyMobilityMatrix

void CIBStrategy::rotateArrayInitalBodyFrame(double* /*array*/,
					     const std::vector<unsigned>& /*struct_ids*/,
					     const bool /*isTranspose*/,
					     const int /*managing_rank*/,
					     const bool /*BodyComps*/)
{
    // intentionally left blank.
    return;
} // rotateArrayInitalBodyFrame

Eigen::Quaterniond* CIBStrategy::getBodyQuaternion(const unsigned int part, const bool halfStep)
{
    if (halfStep) return &d_quaternion_half[part];
    return &d_quaternion_current[part];
};

void CIBStrategy::setRigidBodyClonesParameters(const int num_structs_types, 
					       const std::vector<int>& structs_clones_num)
{
    d_num_structs_types = num_structs_types;
    d_structs_clones_num = structs_clones_num;
};

void CIBStrategy::getRigidBodyClonesParameters(int& num_structs_types, 
					       std::vector<int>& structs_clones_num)
{
    num_structs_types = d_num_structs_types;
    structs_clones_num = d_structs_clones_num;
};

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR
