// Filename: CIBStrategy.cpp
// Created on 9 Nov 2014 by Amneet Bhalla
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/CIBStrategy.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h"
#include "tbox/MathUtilities.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

CIBStrategy::CIBStrategy(const unsigned int parts) : d_num_rigid_parts(parts)
{
    // Resize some arrays
    d_center_of_mass_initial.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());
    d_center_of_mass_current.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());
    d_center_of_mass_half.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());
    d_center_of_mass_new.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());
    d_compute_center_of_mass_initial.resize(d_num_rigid_parts, true);

    d_quaternion_current.resize(d_num_rigid_parts, Eigen::Quaterniond::Identity());
    d_quaternion_half.resize(d_num_rigid_parts, Eigen::Quaterniond::Identity());
    d_quaternion_new.resize(d_num_rigid_parts, Eigen::Quaterniond::Identity());

    d_trans_vel_current.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());
    d_trans_vel_half.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());
    d_trans_vel_new.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());

    d_rot_vel_current.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());
    d_rot_vel_half.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());
    d_rot_vel_new.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());

    d_net_rigid_generalized_force.resize(d_num_rigid_parts);
    d_solve_rigid_vel.resize(d_num_rigid_parts);
    d_free_dofs_map.resize(d_num_rigid_parts, std::make_pair(-1, -1));
    d_free_dofs_map_updated = false;

    return;
} // CIBStrategy

CIBStrategy::~CIBStrategy()
{
    // intentionally left blank
    return;
} // ~CIBStrategy

void
CIBStrategy::getFreeRigidVelocities(Vec* U, const double /*data_time*/)
{
    *U = d_U;
    return;

} // getFreeRigidVelocities

void
CIBStrategy::getNetExternalForceTorque(Vec* F, const double /*data_time*/)
{
    *F = d_F;
    return;

} // getNetExternalForceTorque

void
CIBStrategy::setInterpolatedVelocityVector(Vec /*V*/, const double /*data_time*/)
{
    // intentionally left blank
    return;
} // setInterpolatedVelocityVector

unsigned int
CIBStrategy::getNumberOfRigidStructures() const
{
    return d_num_rigid_parts;

} // getNumberOfRigidStructures

void
CIBStrategy::setInitialCenterOfMass(const unsigned int part, const Eigen::Vector3d& XCOM_0)
{

#if !defined(NDEBUG)
    TBOX_ASSERT(part < d_num_rigid_parts);
#endif
    d_center_of_mass_initial[part] = XCOM_0;
    d_compute_center_of_mass_initial[part] = false;

}// setInitialCenterOfMass

void
CIBStrategy::setSolveRigidBodyVelocity(const unsigned int part, const FreeRigidDOFVector& solve_rigid_vel)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(part < d_num_rigid_parts);
#endif
    d_solve_rigid_vel[part] = solve_rigid_vel;
    d_free_dofs_map_updated = false;

    return;

} // setSolveRigidBodyVelocity

const FreeRigidDOFVector&
CIBStrategy::getSolveRigidBodyVelocity(const unsigned int part, int& num_free_dofs) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(part < d_num_rigid_parts);
#endif

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

void
CIBStrategy::updateFreeDOFsMapping()
{
    if (d_free_dofs_map_updated) return;

    int free_dofs_counter = 0;
    for (unsigned part = 0; part < d_num_rigid_parts; ++part)
    {
        int n_free_dofs = 0;
        getSolveRigidBodyVelocity(part, n_free_dofs);
        if (n_free_dofs)
        {
            d_free_dofs_map[part].first = free_dofs_counter;
            d_free_dofs_map[part].second = d_free_dofs_map[part].first + n_free_dofs;
            free_dofs_counter += n_free_dofs;
        }
        else
        {
            d_free_dofs_map[part].first = -1;
            d_free_dofs_map[part].second = -1;
        }
    }
    d_free_dofs_map_updated = true;

} // updateFreeDOFsMapping

void
CIBStrategy::setRigidBodyVelocity(unsigned int part, Vec U, Vec V)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(part < d_num_rigid_parts);
    PetscInt size;
    VecGetSize(U, &size);
    TBOX_ASSERT(size == s_max_free_dofs);
#endif

    RigidDOFVector Ur;
    vecToRDV(U, Ur);
    setRigidBodyVelocity(part, Ur, V);

    return;
} // setRigidBodyVelocity

void
CIBStrategy::setRigidBodyVelocity(Vec U,
                                  Vec V,
                                  const bool only_free_dofs,
                                  const bool only_imposed_dofs,
                                  const bool all_dofs)
{
    // Get rigid DOFs on all processors.
    Vec U_all;
    VecScatter ctx;
    VecScatterCreateToAll(U, &ctx, &U_all);
    VecScatterBegin(ctx, U, U_all, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, U, U_all, INSERT_VALUES, SCATTER_FORWARD);
    const PetscScalar* U_array;
    VecGetArrayRead(U_all, &U_array);

    if (only_free_dofs)
    {
        int part_free_dofs_begin = 0;
        for (unsigned part = 0; part < d_num_rigid_parts; ++part)
        {
            int num_free_dofs;
            const FreeRigidDOFVector& solve_dofs = getSolveRigidBodyVelocity(part, num_free_dofs);
            if (!num_free_dofs) continue;

            RigidDOFVector U_part;
            U_part.setZero();

            const PetscScalar* a = &U_array[part_free_dofs_begin];
            for (int k = 0, p = 0; k < s_max_free_dofs; ++k)
            {
                if (solve_dofs[k])
                {
                    U_part[k] = a[p];
                    ++p;
                }
            }
            setRigidBodyVelocity(part, U_part, V);
            part_free_dofs_begin += num_free_dofs;
        }
    }
    else if (only_imposed_dofs)
    {
        int part_imposed_dofs_begin = 0;
        for (unsigned part = 0; part < d_num_rigid_parts; ++part)
        {
            int num_free_dofs;
            const FreeRigidDOFVector& solve_dofs = getSolveRigidBodyVelocity(part, num_free_dofs);
            if (num_free_dofs == s_max_free_dofs) continue;

            RigidDOFVector U_part;
            U_part.setZero();

            const PetscScalar* a = &U_array[part_imposed_dofs_begin];
            for (int k = 0, p = 0; k < s_max_free_dofs; ++k)
            {
                if (!solve_dofs[k])
                {
                    U_part[k] = a[p];
                    ++p;
                }
            }
            setRigidBodyVelocity(part, U_part, V);
            part_imposed_dofs_begin += (s_max_free_dofs - num_free_dofs);
        }
    }
    else if (all_dofs)
    {
        int part_dofs_begin = 0;
        for (unsigned part = 0; part < d_num_rigid_parts; ++part)
        {
            RigidDOFVector U_part;
            U_part.setZero();

            const PetscScalar* a = &U_array[part_dofs_begin];
            for (int k = 0; k < s_max_free_dofs; ++k)
            {
                U_part[k] = a[k];
            }

            setRigidBodyVelocity(part, U_part, V);
            part_dofs_begin += s_max_free_dofs;
        }
    }

    VecRestoreArrayRead(U_all, &U_array);
    VecScatterDestroy(&ctx);
    VecDestroy(&U_all);

    return;

} // setRigidBodyVelocity

void
CIBStrategy::computeNetRigidGeneralizedForce(const unsigned int part, Vec L, Vec F)
{
    RigidDOFVector f;
    computeNetRigidGeneralizedForce(part, L, f);
    rdvToVec(f, F);

    return;
} // computeNetRigidGeneralizedForce

void
CIBStrategy::computeNetRigidGeneralizedForce(Vec L,
                                             Vec F,
                                             const bool only_free_dofs,
                                             const bool only_imposed_dofs,
                                             const bool all_dofs)
{
    // Here we use the fact that all vector enteries of F are on
    // a single processor, and we can set values directly in the array
    // rather than using the costly VecSetValues() followed by VecAssemblyBegin/End().
    PetscScalar* F_array = NULL;
    VecGetArray(F, &F_array);

    if (only_free_dofs)
    {
        int part_free_dofs_begin = 0;
        for (unsigned part = 0; part < d_num_rigid_parts; ++part)
        {
            int num_free_dofs;
            const FRDV& solve_dofs = getSolveRigidBodyVelocity(part, num_free_dofs);
            if (!num_free_dofs) continue;

            RigidDOFVector F_part;
            computeNetRigidGeneralizedForce(part, L, F_part);

            if (F_array != NULL)
            {
                PetscScalar* f = &F_array[part_free_dofs_begin];
                for (int k = 0, p = 0; k < s_max_free_dofs; ++k)
                {
                    if (solve_dofs[k])
                    {
                        f[p] = F_part[k];
                        ++p;
                    }
                }
            }
            part_free_dofs_begin += num_free_dofs;
        }
    }
    else if (only_imposed_dofs)
    {
        int part_imposed_dofs_begin = 0;
        for (unsigned part = 0; part < d_num_rigid_parts; ++part)
        {
            int num_free_dofs;
            const FRDV& solve_dofs = getSolveRigidBodyVelocity(part, num_free_dofs);
            if (num_free_dofs == s_max_free_dofs) continue;

            RigidDOFVector F_part;
            computeNetRigidGeneralizedForce(part, L, F_part);

            if (F_array != NULL)
            {
                PetscScalar* f = &F_array[part_imposed_dofs_begin];
                for (int k = 0, p = 0; k < s_max_free_dofs; ++k)
                {
                    if (!solve_dofs[k])
                    {
                        f[p] = F_part[k];
                        ++p;
                    }
                }
            }
            part_imposed_dofs_begin += (s_max_free_dofs - num_free_dofs);
        }
    }
    else if (all_dofs)
    {
        int part_dofs_begin = 0;
        for (unsigned part = 0; part < d_num_rigid_parts; ++part)
        {
            RigidDOFVector F_part;
            computeNetRigidGeneralizedForce(part, L, F_part);

            if (F_array != NULL)
            {
                PetscScalar* f = &F_array[part_dofs_begin];
                for (int k = 0; k < s_max_free_dofs; ++k)
                {
                    f[k] = F_part[k];
                }
            }
            part_dofs_begin += s_max_free_dofs;
        }
    }
    VecRestoreArray(F, &F_array);

    return;
} // computeNetRigidGeneralizedForce

const RigidDOFVector&
CIBStrategy::getNetRigidGeneralizedForce(const unsigned int part)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(part < d_num_rigid_parts);
#endif

    return d_net_rigid_generalized_force[part];
} // getNetRigidGeneralizedForce

void
CIBStrategy::updateNewRigidBodyVelocity(const unsigned int part, const RigidDOFVector& U)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(part < d_num_rigid_parts);
#endif

    rdvToEigen(U, d_trans_vel_new[part], d_rot_vel_new[part]);

    d_trans_vel_half[part] = 0.5 * (d_trans_vel_current[part] + d_trans_vel_new[part]);
    d_rot_vel_half[part] = 0.5 * (d_rot_vel_current[part] + d_rot_vel_new[part]);

    return;
} // updateNewRigidBodyVelocity

void
CIBStrategy::updateNewRigidBodyVelocity(const unsigned int part, Vec U)
{
    RigidDOFVector Ur;
    vecToRDV(U, Ur);
    updateNewRigidBodyVelocity(part, Ur);

    return;
} // updateNewRigidBodyVelocity

void
CIBStrategy::updateNewRigidBodyVelocity(Vec U,
                                        const bool only_free_dofs,
                                        const bool only_imposed_dofs,
                                        const bool all_dofs)
{
    // Get rigid DOFs on all processors.
    Vec U_all;
    VecScatter ctx;
    VecScatterCreateToAll(U, &ctx, &U_all);
    VecScatterBegin(ctx, U, U_all, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, U, U_all, INSERT_VALUES, SCATTER_FORWARD);
    const PetscScalar* U_array;
    VecGetArrayRead(U_all, &U_array);

    if (only_free_dofs)
    {
        int part_free_dofs_begin = 0;
        for (unsigned part = 0; part < d_num_rigid_parts; ++part)
        {
            int num_free_dofs;
            const FreeRigidDOFVector& solve_dofs = getSolveRigidBodyVelocity(part, num_free_dofs);
            if (!num_free_dofs) continue;

            RDV UW;
            eigenToRDV(d_trans_vel_new[part], d_rot_vel_new[part], UW);
            const PetscScalar* a = &U_array[part_free_dofs_begin];
            for (int k = 0, p = 0; k < s_max_free_dofs; ++k)
            {
                if (solve_dofs[k])
                {
                    UW[k] = a[p];
                    ++p;
                }
            }
            rdvToEigen(UW, d_trans_vel_new[part], d_rot_vel_new[part]);

            d_trans_vel_half[part] = 0.5 * (d_trans_vel_current[part] + d_trans_vel_new[part]);
            d_rot_vel_half[part] = 0.5 * (d_rot_vel_current[part] + d_rot_vel_new[part]);

            part_free_dofs_begin += num_free_dofs;
        }
    }
    else if (only_imposed_dofs)
    {
        int part_imposed_dofs_begin = 0;
        for (unsigned part = 0; part < d_num_rigid_parts; ++part)
        {
            int num_free_dofs;
            const FreeRigidDOFVector& solve_dofs = getSolveRigidBodyVelocity(part, num_free_dofs);
            if (num_free_dofs == s_max_free_dofs) continue;

            RDV UW;
            eigenToRDV(d_trans_vel_new[part], d_rot_vel_new[part], UW);
            const PetscScalar* a = &U_array[part_imposed_dofs_begin];
            for (int k = 0, p = 0; k < s_max_free_dofs; ++k)
            {
                if (!solve_dofs[k])
                {
                    UW[k] = a[p];
                    ++p;
                }
            }
            rdvToEigen(UW, d_trans_vel_new[part], d_rot_vel_new[part]);

            d_trans_vel_half[part] = 0.5 * (d_trans_vel_current[part] + d_trans_vel_new[part]);
            d_rot_vel_half[part] = 0.5 * (d_rot_vel_current[part] + d_rot_vel_new[part]);

            part_imposed_dofs_begin += (s_max_free_dofs - num_free_dofs);
        }
    }
    else if (all_dofs)
    {
        int part_dofs_begin = 0;
        for (unsigned part = 0; part < d_num_rigid_parts; ++part)
        {
            RDV UW;
            const PetscScalar* a = &U_array[part_dofs_begin];
            for (int k = 0; k < s_max_free_dofs; ++k)
            {
                UW[k] = a[k];
            }
            rdvToEigen(UW, d_trans_vel_new[part], d_rot_vel_new[part]);

            d_trans_vel_half[part] = 0.5 * (d_trans_vel_current[part] + d_trans_vel_new[part]);
            d_rot_vel_half[part] = 0.5 * (d_rot_vel_current[part] + d_rot_vel_new[part]);

            part_dofs_begin += s_max_free_dofs;
        }
    }

    VecRestoreArrayRead(U_all, &U_array);
    VecScatterDestroy(&ctx);
    VecDestroy(&U_all);

    return;
} // updateNewRigidBodyVelocity

void
CIBStrategy::copyVecToArray(Vec /*b*/,
                            double* /*array*/,
                            const std::vector<unsigned>& /*struct_ids*/,
                            const int /*data_depth*/,
                            const int /*array_rank*/)
{
    // intentionally left-blank
    return;
} // copyVecToArray

void
CIBStrategy::copyFreeDOFsVecToArray(Vec b, double* array, const std::vector<unsigned>& struct_ids, const int array_rank)
{
    if (struct_ids.empty()) return;
    const int num_structs = static_cast<int>(struct_ids.size());

    // Get the mapping of indices in the two vectors.
    updateFreeDOFsMapping();
    std::vector<int> map_vec, map_array;
    map_vec.reserve(num_structs * s_max_free_dofs);
    map_array.reserve(num_structs * s_max_free_dofs);
    int local_counter = 0;
    for (int part = 0; part < num_structs; ++part)
    {
        const unsigned struct_id = struct_ids[part];
        int struct_free_dofs = 0;
        const FreeRigidDOFVector& solve_dofs = getSolveRigidBodyVelocity(struct_id, struct_free_dofs);

        for (int k = 0; k < s_max_free_dofs; ++k, ++local_counter)
        {
            if (solve_dofs[k])
            {
                map_array.push_back(local_counter);
            }
        }

        if (!struct_free_dofs) continue;

        const int& idx_begin = d_free_dofs_map[struct_id].first;
        const int& idx_end = d_free_dofs_map[struct_id].second;
        for (int i = idx_begin; i < idx_end; ++i)
        {
            map_vec.push_back(i);
        }
    }
    int idx_size = static_cast<int>(map_vec.size());
#if !defined(NDEBUG)
    TBOX_ASSERT(idx_size == static_cast<int>(map_array.size()));
#endif

    // Wrap the raw data in a PETSc Vec.
    int rank = SAMRAI_MPI::getRank();
    PetscInt array_size = num_structs * s_max_free_dofs;
    PetscInt array_local_size = 0;
    if (rank == array_rank) array_local_size = array_size;
    Vec array_vec;
    VecCreateMPIWithArray(PETSC_COMM_WORLD, /*blocksize*/ 1, array_local_size, PETSC_DECIDE, array, &array_vec);
    VecSet(array_vec, 0.0);
    VecAssemblyBegin(array_vec);
    VecAssemblyEnd(array_vec);

    // Create index sets to define global index mapping of the two vectors.
    IS is_vec;
    IS is_array;
    ISCreateGeneral(PETSC_COMM_SELF, idx_size, &map_vec[0], PETSC_COPY_VALUES, &is_vec);
    ISCreateGeneral(PETSC_COMM_SELF, idx_size, &map_array[0], PETSC_COPY_VALUES, &is_array);

    // Scatter values
    VecScatter ctx;
    VecScatterCreate(b, is_vec, array_vec, is_array, &ctx);
    VecScatterBegin(ctx, b, array_vec, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, b, array_vec, INSERT_VALUES, SCATTER_FORWARD);

    // Cleanup temporary objects.
    VecScatterDestroy(&ctx);
    ISDestroy(&is_vec);
    ISDestroy(&is_array);
    VecDestroy(&array_vec);

    return;
} // copyFreeDOFsVecToArray

void
CIBStrategy::copyArrayToVec(Vec /*b*/,
                            double* /*array*/,
                            const std::vector<unsigned>& /*struct_ids*/,
                            const int /*data_depth*/,
                            const int /*array_rank*/)
{
    // intentionally left-blank
    return;
} // copyArrayToVec

void
CIBStrategy::copyFreeDOFsArrayToVec(Vec b, double* array, const std::vector<unsigned>& struct_ids, const int array_rank)
{
    if (struct_ids.empty()) return;
    const int num_structs = static_cast<int>(struct_ids.size());

    // Get the mapping of indices in the two vectors.
    updateFreeDOFsMapping();
    std::vector<int> map_vec, map_array;
    map_vec.reserve(num_structs * s_max_free_dofs);
    map_array.reserve(num_structs * s_max_free_dofs);
    int local_counter = 0;
    for (int part = 0; part < num_structs; ++part)
    {
        const unsigned struct_id = struct_ids[part];
        int struct_free_dofs = 0;
        const FreeRigidDOFVector& solve_dofs = getSolveRigidBodyVelocity(struct_id, struct_free_dofs);

        for (int k = 0; k < s_max_free_dofs; ++k, ++local_counter)
        {
            if (solve_dofs[k])
            {
                map_array.push_back(local_counter);
            }
        }

        if (!struct_free_dofs) continue;

        const int& idx_begin = d_free_dofs_map[struct_id].first;
        const int& idx_end = d_free_dofs_map[struct_id].second;
        for (int i = idx_begin; i < idx_end; ++i)
        {
            map_vec.push_back(i);
        }
    }
    int idx_size = static_cast<int>(map_vec.size());
#if !defined(NDEBUG)
    TBOX_ASSERT(idx_size == static_cast<int>(map_array.size()));
#endif

    // Wrap the raw data in a PETSc Vec.
    int rank = SAMRAI_MPI::getRank();
    PetscInt array_size = num_structs * s_max_free_dofs;
    PetscInt array_local_size = 0;
    if (rank == array_rank) array_local_size = array_size;
    Vec array_vec;
    VecCreateMPIWithArray(PETSC_COMM_WORLD, /*blocksize*/ 1, array_local_size, PETSC_DECIDE, array, &array_vec);

    // Create index sets to define global index mapping of the two vectors.
    IS is_vec;
    IS is_array;
    ISCreateGeneral(PETSC_COMM_SELF, idx_size, &map_vec[0], PETSC_COPY_VALUES, &is_vec);
    ISCreateGeneral(PETSC_COMM_SELF, idx_size, &map_array[0], PETSC_COPY_VALUES, &is_array);

    // Scatter values
    VecScatter ctx;
    VecScatterCreate(array_vec, is_array, b, is_vec, &ctx);
    VecScatterBegin(ctx, array_vec, b, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, array_vec, b, INSERT_VALUES, SCATTER_FORWARD);

    // Destroy temporary objects
    VecScatterDestroy(&ctx);
    ISDestroy(&is_vec);
    ISDestroy(&is_array);
    VecDestroy(&array_vec);

    return;
} // copyFreeDOFsArrayToVec

void
CIBStrategy::vecToRDV(Vec U, RigidDOFVector& Ur)
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

void
CIBStrategy::rdvToVec(const RigidDOFVector& Ur, Vec& U)
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

void
CIBStrategy::eigenToRDV(const Eigen::Vector3d& U, const Eigen::Vector3d& W, RigidDOFVector& UW)
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

void
CIBStrategy::rdvToEigen(const RigidDOFVector& UW, Eigen::Vector3d& U, Eigen::Vector3d& W)
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

void
CIBStrategy::getCurrentRigidBodyVelocity(const unsigned int part, RigidDOFVector& U)
{
    eigenToRDV(d_trans_vel_current[part], d_rot_vel_current[part], U);

    return;
} // getCurrentRigidBodyVelocity

void
CIBStrategy::getNewRigidBodyVelocity(const unsigned int part, RigidDOFVector& U)
{
    eigenToRDV(d_trans_vel_new[part], d_rot_vel_new[part], U);

    return;
} // getNewRigidBodyVelocity

const Eigen::Vector3d&
CIBStrategy::getCurrentBodyCenterOfMass(const unsigned int part)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(part < d_num_rigid_parts);
#endif
    return d_center_of_mass_current[part];

} // getCurrentBodyCenterOfMass

const Eigen::Vector3d&
CIBStrategy::getNewBodyCenterOfMass(const unsigned int part)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(part < d_num_rigid_parts);
#endif

    return d_center_of_mass_new[part];

} // getMidPointBodyCenterOfMass

void
CIBStrategy::constructMobilityMatrix(const std::string& /*mat_name*/,
                                     MobilityMatrixType /*mat_type*/,
                                     Mat& /*mobility_mat*/,
                                     const std::vector<unsigned>& /*prototype_struct_ids*/,
                                     const double* /*grid_dx*/,
                                     const double* /*domain_extents*/,
                                     const bool /*initial_time*/,
                                     double /*rho*/,
                                     double /*mu*/,
                                     const std::pair<double, double>& /*scale*/,
                                     double /*f_periodic_corr*/,
                                     const int /*managing_rank*/)
{
    // intentionally left blank.

    return;
} // constructMobilityMatrix

void
CIBStrategy::constructGeometricMatrix(const std::string& /*mat_name*/,
                                      Mat& /*geometric_mat*/,
                                      const std::vector<unsigned>& /*prototype_struct_ids*/,
                                      const bool /*initial_time*/,
                                      const int /*managing_rank*/)
{
    // intentionally left blank.

    return;
} // constructGeometricMatrix

void
CIBStrategy::rotateArray(double* /*array*/,
                         const std::vector<unsigned>& /*struct_ids*/,
                         const bool /*use_transpose*/,
                         const int /*managing_rank*/,
                         const int /*depth*/)
{
    // intentionally left blank.

    return;
} // rotateArray

/////////////////////////////// PROTECTED ////////////////////////////////////

void
CIBStrategy::setRotationMatrix(const std::vector<Eigen::Vector3d>& rot_vel,
                               const std::vector<Eigen::Quaterniond>& q_old,
                               std::vector<Eigen::Quaterniond>& q_new,
                               std::vector<Eigen::Matrix3d>& rot_mat,
                               const double dt)
{
    for (unsigned struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
    {
        const double norm = rot_vel[struct_no].norm();
        if (!MathUtilities<double>::equalEps(norm, 0.0))
        {
            Eigen::Vector3d rot_axis = rot_vel[struct_no] / norm;
            Eigen::Quaterniond q(Eigen::AngleAxisd(norm * dt, rot_axis));
            q_new[struct_no] = (q.normalized() * q_old[struct_no]).normalized();
        }
        else
        {
            q_new[struct_no] = q_old[struct_no];
        }

        rot_mat[struct_no] = q_new[struct_no].toRotationMatrix();
    }
} // setRotationMatrix

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR
