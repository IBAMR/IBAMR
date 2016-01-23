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
#include "ibamr/ibamr_utilities.h"
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
    d_isFree_component.resize(d_num_rigid_parts);
    d_isImposed_component.resize(d_num_rigid_parts);

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

void CIBStrategy::setCloneFreeDOFs(const unsigned int part, const FreeRigidDOFVector& solve_rigid_vel)
{
    d_solve_rigid_vel[part] = solve_rigid_vel;
    int num_free_dofs = 0;
    getSolveRigidBodyVelocity(part, num_free_dofs);

    if (num_free_dofs)
    {
        d_isFree_component[part] = true;
        if (num_free_dofs < s_max_free_dofs) d_isImposed_component[part] = true;
    }
    else
    {
        d_isFree_component[part] = false;
        d_isImposed_component[part] = true;
    }

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

void CIBStrategy::setRigidBodyVelocity(unsigned int part, Vec U, Vec V)
{
    PetscInt size;
    VecGetSize(U, &size);
    TBOX_ASSERT(size == s_max_free_dofs);

    std::vector<bool> skip_comp;
    skip_comp.resize(d_num_rigid_parts);
    fill(skip_comp.begin(), skip_comp.end(), true);
    skip_comp[part] = false;

    // collective call over all nodes
    setRigidBodyVelocity(U, V, skip_comp);
}

void CIBStrategy::setRigidBodyVelocity(Vec U,
                                       Vec V,
                                       const bool only_free_dofs,
                                       const bool only_imposed_dofs,
                                       const bool all_dofs)
{
#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    clock_t end_t = 0, start_med = 0;
    if (SAMRAI_MPI::getRank() == 0) start_med = clock();
#endif
    std::vector<bool> skip_comp;
    skip_comp.resize(d_num_rigid_parts);

    if (only_free_dofs || only_imposed_dofs)
    {
        for (unsigned part = 0; part < d_num_rigid_parts; ++part)
        {
            int num_free_dofs = 0;
            for (int k = 0; k < s_max_free_dofs; ++k) num_free_dofs += d_solve_rigid_vel[part][k];

            if (only_free_dofs && num_free_dofs)
                skip_comp[part] = false;
            else if (only_imposed_dofs && (num_free_dofs < s_max_free_dofs))
                skip_comp[part] = false;
            else
                skip_comp[part] = true;
        }
    }
    else if (all_dofs)
    {
        fill(skip_comp.begin(), skip_comp.end(), false);
    }
    else
    {
        TBOX_ERROR("CIBStrategy::setRigidBodyVelocity() incorrect parameters in the function call" << std::endl);
    }
#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0)
    {
        end_t = clock();
        pout << std::setprecision(4)
             << "         setRigidBodyVelocity: setting parameters, CPU time taken for the time step is:"
             << double(end_t - start_med) / double(CLOCKS_PER_SEC) << std::endl;
        ;
    }
    if (SAMRAI_MPI::getRank() == 0) start_med = clock();
#endif

    // collective call over all nodes
    setRigidBodyVelocity(U, V, skip_comp);

#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0)
    {
        end_t = clock();
        pout << std::setprecision(4)
             << "         setRigidBodyVelocity: collective call, CPU time taken for the time step is:"
             << double(end_t - start_med) / double(CLOCKS_PER_SEC) << std::endl;
        ;
    }
#endif
    return;
} // setRigidBodyVelocity

void CIBStrategy::setRigidBodyDeformationVelocity(Vec W)
{

    // intentionally left blank.
    return;
}

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
            int num_free_dofs = 0;
            for (int k = 0; k < s_max_free_dofs; ++k) num_free_dofs += d_solve_rigid_vel[part][k];

            if (only_free_dofs && num_free_dofs)
            {
                skip_comp[part] = false;
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
        TBOX_ERROR("CIBStrategy::computeNetRigidGeneralizedForce() incorrect parameters in the function call"
                   << std::endl);
    }

    // collective call over all nodes
    computeNetRigidGeneralizedForce(L, F, skip_comp);

    return;
} // computeNetRigidGeneralizedForce

void CIBStrategy::updateNewRigidBodyVelocity(Vec U, const bool all_dofs)
{
    Vec U_all;
    VecScatter ctx;
    VecScatterCreateToAll(U, &ctx, &U_all);

    VecScatterBegin(ctx, U, U_all, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, U, U_all, INSERT_VALUES, SCATTER_FORWARD);

    PetscScalar* u_array = NULL;
    PetscInt comps;
    VecGetSize(U, &comps);
    VecGetArray(U_all, &u_array);

    int counter = 0;
    for (unsigned i = 0; i < d_num_rigid_parts; ++i)
    {
        if (d_isFree_component[i])
        {
            for (int d = 0; d < NDIM; ++d)
            {
                if (d_solve_rigid_vel[i][d]) d_trans_vel_new[i][d] = u_array[counter * s_max_free_dofs + d];
            }

#if (NDIM == 2)
            if (d_solve_rigid_vel[i][NDIM]) d_rot_vel_new[i][2] = u_array[counter * s_max_free_dofs + 3];
#elif(NDIM == 3)
            for (int d = 0; d < NDIM; ++d)
            {
                if (d_solve_rigid_vel[i][d + NDIM]) d_rot_vel_new[i][d] = u_array[counter * s_max_free_dofs + NDIM + d];
            }
#endif
            counter++;
        }
        d_trans_vel_half[i] = 0.5 * (d_trans_vel_current[i] + d_trans_vel_new[i]);
        d_rot_vel_half[i] = 0.5 * (d_rot_vel_current[i] + d_rot_vel_new[i]);
    }

    VecRestoreArray(U_all, &u_array);
    VecScatterDestroy(&ctx);
    VecDestroy(&U_all);

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
                                    const int /*data_depth*/)
{
    // intentionally left-blank
    return;
} // copyArrayToVec

void CIBStrategy::copyAllVecToArray(Vec /*b*/,
                                    double* /*array*/,
                                    const std::vector<unsigned>& /*rhs_struct_ids*/,
                                    const int /*data_depth*/)
{
    // intentionally left-blank
    return;
} // copyVecToArray

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

void CIBStrategy::constructMobilityMatrix(std::map<std::string, double*>& /*mat_map*/,
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

Eigen::Vector3d* CIBStrategy::getBodyCenterOfMass(const unsigned int part, const bool halfStep)
{
    if (halfStep) return &d_center_of_mass_half[part];
    return &d_center_of_mass_current[part];
};

Eigen::Quaterniond* CIBStrategy::getBodyQuaternion(const unsigned int part, const bool halfStep)
{
    if (halfStep) return &d_quaternion_half[part];
    return &d_quaternion_current[part];
};

void CIBStrategy::setRigidBodyClonesParameters(const int num_structs_types, const std::vector<int>& structs_clones_num)
{
    d_num_structs_types = num_structs_types;
    d_structs_clones_num = structs_clones_num;
};

void CIBStrategy::getRigidBodyClonesParameters(int& num_structs_types, std::vector<int>& structs_clones_num)
{
    num_structs_types = d_num_structs_types;
    structs_clones_num = d_structs_clones_num;
};

void CIBStrategy::setHomogeneousBc(bool homogeneous_bc)
{

    d_homogeneous_bc = homogeneous_bc;
    return;
} // setHomogeneousBc

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR
