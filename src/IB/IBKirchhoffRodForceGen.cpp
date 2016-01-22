// Filename: IBKirchhoffRodForceGen.cpp
// Created on 22 Jun 2010 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <algorithm>
#include <ostream>
#include <vector>

#include "Eigen/Geometry"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "boost/array.hpp"
#include "ibamr/IBKirchhoffRodForceGen.h"
#include "ibamr/IBRodForceSpec.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LMesh.h"
#include "ibtk/LNode.h"
#include "ibtk/ibtk_utilities.h"
#include "petscmat.h"
#include "petscsys.h"
#include "petscvec.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#include "unsupported/Eigen/MatrixFunctions"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_compute_lagrangian_force_and_torque;
static Timer* t_initialize_level_data;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBKirchhoffRodForceGen::IBKirchhoffRodForceGen(Pointer<Database> input_db)
    : d_D_next_mats(),
      d_X_next_mats(),
      d_petsc_curr_node_idxs(),
      d_petsc_next_node_idxs(),
      d_material_params(),
      d_is_initialized()
{
    // Initialize object with data read from the input database.
    getFromInput(input_db);

    // Setup Timers.
    IBAMR_DO_ONCE(t_compute_lagrangian_force_and_torque = TimerManager::getManager()->getTimer(
                      "IBAMR::IBKirchhoffRodForceGen::computeLagrangianForceAndTorque()");
                  t_initialize_level_data =
                      TimerManager::getManager()->getTimer("IBAMR::IBKirchhoffRodForceGen::initializeLevelData()"););
    return;
} // IBKirchhoffRodForceGen

IBKirchhoffRodForceGen::~IBKirchhoffRodForceGen()
{
    int ierr;
    for (std::vector<Mat>::iterator it = d_D_next_mats.begin(); it != d_D_next_mats.end(); ++it)
    {
        if (*it)
        {
            ierr = MatDestroy(&*it);
            IBTK_CHKERRQ(ierr);
        }
    }
    for (std::vector<Mat>::iterator it = d_X_next_mats.begin(); it != d_X_next_mats.end(); ++it)
    {
        if (*it)
        {
            ierr = MatDestroy(&*it);
            IBTK_CHKERRQ(ierr);
        }
    }
    return;
} // ~IBKirchhoffRodForceGen

void
IBKirchhoffRodForceGen::initializeLevelData(const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                            const int level_number,
                                            const double /*init_data_time*/,
                                            const bool /*initial_time*/,
                                            LDataManager* const l_data_manager)
{
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return;

    IBAMR_TIMER_START(t_initialize_level_data);

#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
#endif
    int ierr;

    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Resize the vectors corresponding to data individually maintained for
    // separate levels of the patch hierarchy.
    const int level_num = level->getLevelNumber();
    const int new_size = std::max(level_num + 1, static_cast<int>(d_is_initialized.size()));

    d_D_next_mats.resize(new_size);
    d_X_next_mats.resize(new_size);
    d_petsc_curr_node_idxs.resize(new_size);
    d_petsc_next_node_idxs.resize(new_size);
    d_material_params.resize(new_size);
    d_is_initialized.resize(new_size, false);

    Mat& D_next_mat = d_D_next_mats[level_num];
    Mat& X_next_mat = d_X_next_mats[level_num];
    std::vector<int>& petsc_curr_node_idxs = d_petsc_curr_node_idxs[level_num];
    std::vector<int>& petsc_next_node_idxs = d_petsc_next_node_idxs[level_num];
    std::vector<boost::array<double, IBRodForceSpec::NUM_MATERIAL_PARAMS> >& material_params =
        d_material_params[level_num];

    if (D_next_mat)
    {
        ierr = MatDestroy(&D_next_mat);
        IBTK_CHKERRQ(ierr);
    }
    if (X_next_mat)
    {
        ierr = MatDestroy(&X_next_mat);
        IBTK_CHKERRQ(ierr);
    }
    petsc_curr_node_idxs.clear();
    petsc_next_node_idxs.clear();
    material_params.clear();

    // The LMesh object provides the set of local Lagrangian nodes.
    const Pointer<LMesh> mesh = l_data_manager->getLMesh(level_num);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

    // Determine the "next" node indices for all rods associated with the
    // present MPI process.
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
        const IBRodForceSpec* const force_spec = node_idx->getNodeDataItem<IBRodForceSpec>();
        if (force_spec)
        {
            const int& curr_idx = node_idx->getLagrangianIndex();
            const unsigned int num_rods = force_spec->getNumberOfRods();
#if !defined(NDEBUG)
            TBOX_ASSERT(curr_idx == force_spec->getMasterNodeIndex());
#endif
            const std::vector<int>& next_idxs = force_spec->getNextNodeIndices();
            const std::vector<boost::array<double, IBRodForceSpec::NUM_MATERIAL_PARAMS> >& params =
                force_spec->getMaterialParams();
#if !defined(NDEBUG)
            TBOX_ASSERT(num_rods == next_idxs.size());
#endif
            for (unsigned int k = 0; k < num_rods; ++k)
            {
                petsc_curr_node_idxs.push_back(curr_idx);
                petsc_next_node_idxs.push_back(next_idxs[k]);
                material_params.push_back(params[k]);
            }
        }
    }

    // Map the Lagrangian node indices to the PETSc indices corresponding to the
    // present data distribution.
    l_data_manager->mapLagrangianToPETSc(petsc_curr_node_idxs, level_num);
    l_data_manager->mapLagrangianToPETSc(petsc_next_node_idxs, level_num);

    // Determine the global node offset and the number of local nodes.
    const int global_node_offset = l_data_manager->getGlobalNodeOffset(level_num);
    const int num_local_nodes = l_data_manager->getNumberOfLocalNodes(level_num);

    // Determine the non-zero structure for the matrices.
    const int local_sz = static_cast<int>(petsc_curr_node_idxs.size());

    std::vector<int> next_d_nz(local_sz, 1), next_o_nz(local_sz, 0);
    for (int k = 0; k < local_sz; ++k)
    {
        const int& next_idx = petsc_next_node_idxs[k];
        if (next_idx >= global_node_offset && next_idx < global_node_offset + num_local_nodes)
        {
            ++next_d_nz[k]; // a "local"    next index
        }
        else
        {
            ++next_o_nz[k]; // a "nonlocal" next index
        }
    }

    // Create new MPI block AIJ matrices and set the values of the non-zero
    // entries.
    {
        ierr = MatCreateBAIJ(PETSC_COMM_WORLD,
                             3 * 3,
                             3 * 3 * local_sz,
                             3 * 3 * num_local_nodes,
                             PETSC_DETERMINE,
                             PETSC_DETERMINE,
                             0,
                             local_sz ? &next_d_nz[0] : NULL,
                             0,
                             local_sz ? &next_o_nz[0] : NULL,
                             &D_next_mat);
        IBTK_CHKERRQ(ierr);

        Eigen::Matrix<double, 3 * 3, 3 * 3, Eigen::RowMajor> curr_vals(
            Eigen::Matrix<double, 3 * 3, 3 * 3, Eigen::RowMajor>::Zero());
        Eigen::Matrix<double, 3 * 3, 3 * 3, Eigen::RowMajor> next_vals(
            Eigen::Matrix<double, 3 * 3, 3 * 3, Eigen::RowMajor>::Zero());
        for (unsigned int d = 0; d < 3 * 3; ++d)
        {
            curr_vals(d, d) = 0.0;
            next_vals(d, d) = +1.0;
        }

        int i_offset;

        ierr = MatGetOwnershipRange(D_next_mat, &i_offset, NULL);
        IBTK_CHKERRQ(ierr);
        i_offset /= 3 * 3;

        for (int k = 0; k < local_sz; ++k)
        {
            int i = i_offset + k;
            int j_curr = petsc_curr_node_idxs[k];
            int j_next = petsc_next_node_idxs[k];
            ierr = MatSetValuesBlocked(D_next_mat, 1, &i, 1, &j_curr, curr_vals.data(), INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetValuesBlocked(D_next_mat, 1, &i, 1, &j_next, next_vals.data(), INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
    }

    {
        ierr = MatCreateBAIJ(PETSC_COMM_WORLD,
                             NDIM,
                             NDIM * local_sz,
                             NDIM * num_local_nodes,
                             PETSC_DETERMINE,
                             PETSC_DETERMINE,
                             0,
                             local_sz ? &next_d_nz[0] : NULL,
                             0,
                             local_sz ? &next_o_nz[0] : NULL,
                             &X_next_mat);
        IBTK_CHKERRQ(ierr);

        Matrix curr_vals(Matrix::Zero());
        Matrix next_vals(Matrix::Zero());
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            curr_vals(d, d) = 0.0;
            next_vals(d, d) = +1.0;
        }

        int i_offset;

        ierr = MatGetOwnershipRange(X_next_mat, &i_offset, NULL);
        IBTK_CHKERRQ(ierr);
        i_offset /= NDIM;

        for (int k = 0; k < local_sz; ++k)
        {
            int i = i_offset + k;
            int j_curr = petsc_curr_node_idxs[k];
            int j_next = petsc_next_node_idxs[k];
            ierr = MatSetValuesBlocked(X_next_mat, 1, &i, 1, &j_curr, curr_vals.data(), INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetValuesBlocked(X_next_mat, 1, &i, 1, &j_next, next_vals.data(), INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
    }

    // Assemble the matrices.
    ierr = MatAssemblyBegin(D_next_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyBegin(X_next_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(D_next_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(X_next_mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);

    // Indicate that the level data has been initialized.
    d_is_initialized[level_num] = true;

    IBAMR_TIMER_STOP(t_initialize_level_data);
    return;
} // initializeLevelData

void
IBKirchhoffRodForceGen::computeLagrangianForceAndTorque(Pointer<LData> F_data,
                                                        Pointer<LData> N_data,
                                                        Pointer<LData> X_data,
                                                        Pointer<LData> D_data,
                                                        const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                                        const int level_number,
                                                        const double /*data_time*/,
                                                        LDataManager* const l_data_manager)
{
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return;

    IBAMR_TIMER_START(t_compute_lagrangian_force_and_torque);

#if !defined(NDEBUG)
    TBOX_ASSERT(level_number < static_cast<int>(d_is_initialized.size()));
    TBOX_ASSERT(d_is_initialized[level_number]);
#endif

    const int global_offset = l_data_manager->getGlobalNodeOffset(level_number);

    int ierr;

    // Create appropriately sized temporary vectors.
    int i_start, i_stop;

    ierr = MatGetOwnershipRange(d_D_next_mats[level_number], &i_start, &i_stop);
    IBTK_CHKERRQ(ierr);

    Vec D_vec = D_data->getVec();
    Vec D_next_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop - i_start, PETSC_DECIDE, &D_next_vec);
    IBTK_CHKERRQ(ierr);

    ierr = MatGetOwnershipRange(d_X_next_mats[level_number], &i_start, &i_stop);
    IBTK_CHKERRQ(ierr);

    Vec X_vec = X_data->getVec();
    Vec X_next_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop - i_start, PETSC_DECIDE, &X_next_vec);
    IBTK_CHKERRQ(ierr);

    // Compute the node displacements.
    ierr = MatMult(d_D_next_mats[level_number], D_vec, D_next_vec);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(d_X_next_mats[level_number], X_vec, X_next_vec);
    IBTK_CHKERRQ(ierr);

    // Compute the rod forces acting on the nodes of the Lagrangian mesh.
    double* D_vals;
    ierr = VecGetArray(D_vec, &D_vals);
    IBTK_CHKERRQ(ierr);

    double* D_next_vals;
    ierr = VecGetArray(D_next_vec, &D_next_vals);
    IBTK_CHKERRQ(ierr);

    double* X_vals;
    ierr = VecGetArray(X_vec, &X_vals);
    IBTK_CHKERRQ(ierr);

    double* X_next_vals;
    ierr = VecGetArray(X_next_vec, &X_next_vals);
    IBTK_CHKERRQ(ierr);

    std::vector<int>& petsc_curr_node_idxs = d_petsc_curr_node_idxs[level_number];
    std::vector<int>& petsc_next_node_idxs = d_petsc_next_node_idxs[level_number];
    const std::vector<boost::array<double, IBRodForceSpec::NUM_MATERIAL_PARAMS> >& material_params =
        d_material_params[level_number];

    const size_t local_sz = petsc_curr_node_idxs.size();
    std::vector<double> F_curr_node_vals(NDIM * local_sz, 0.0);
    std::vector<double> N_curr_node_vals(NDIM * local_sz, 0.0);
    std::vector<double> F_next_node_vals(NDIM * local_sz, 0.0);
    std::vector<double> N_next_node_vals(NDIM * local_sz, 0.0);

    for (unsigned int k = 0; k < local_sz; ++k)
    {
        // Compute the forces applied by the rod to the "current" and "next"
        // nodes.
        const int D1_offset = 0;
        Eigen::Map<const Vector3d> D1(&D_vals[(petsc_curr_node_idxs[k] - global_offset) * 3 * 3 + D1_offset]);
        Eigen::Map<const Vector3d> D1_next(&D_next_vals[k * 3 * 3 + D1_offset]);

        const int D2_offset = 3;
        Eigen::Map<const Vector3d> D2(&D_vals[(petsc_curr_node_idxs[k] - global_offset) * 3 * 3 + D2_offset]);
        Eigen::Map<const Vector3d> D2_next(&D_next_vals[k * 3 * 3 + D2_offset]);

        const int D3_offset = 6;
        Eigen::Map<const Vector3d> D3(&D_vals[(petsc_curr_node_idxs[k] - global_offset) * 3 * 3 + D3_offset]);
        Eigen::Map<const Vector3d> D3_next(&D_next_vals[k * 3 * 3 + D3_offset]);

        Eigen::Map<const Vector3d> X(&X_vals[(petsc_curr_node_idxs[k] - global_offset) * NDIM]);
        Eigen::Map<const Vector3d> X_next(&X_next_vals[k * NDIM]);

        boost::array<Eigen::Map<const Vector3d>*, 3> D = { { &D1, &D2, &D3 } };
        boost::array<Eigen::Map<const Vector3d>*, 3> D_next = { { &D1_next, &D2_next, &D3_next } };
        Matrix3d A(Matrix3d::Zero());
        for (int i = 0; i < 3; ++i)
        {
            A += (*D_next[i]) * (*D[i]).transpose();
        }
        Matrix3d sqrt_A = A.sqrt();
        Vector3d D1_half, D2_half, D3_half;
        boost::array<Vector3d*, 3> D_half = { { &D1_half, &D2_half, &D3_half } };
        for (int i = 0; i < 3; ++i)
        {
            *D_half[i] = sqrt_A * (*D[i]);
        }

        const double ds = material_params[k][0];
        const double a1 = material_params[k][1];
        const double a2 = material_params[k][2];
        const double a3 = material_params[k][3];
        const double b1 = material_params[k][4];
        const double b2 = material_params[k][5];
        const double b3 = material_params[k][6];
        const double kappa1 = material_params[k][7];
        const double kappa2 = material_params[k][8];
        const double tau = material_params[k][9];

        const Vector3d dX_ds((X_next - X) / ds);
        const double F1 = b1 * D1_half.dot(dX_ds);
        const double F2 = b2 * D2_half.dot(dX_ds);
        const double F3 = b3 * (D3_half.dot(dX_ds) - 1.0);
        const Vector3d F_half = F1 * D1_half + F2 * D2_half + F3 * D3_half;

        const Vector3d dD1_ds((D1_next - D1) / ds);
        const Vector3d dD2_ds((D2_next - D2) / ds);
        const Vector3d dD3_ds((D3_next - D3) / ds);
        const double N1 = a1 * (dD2_ds.dot(D3_half) - kappa1);
        const double N2 = a2 * (dD3_ds.dot(D1_half) - kappa2);
        const double N3 = a3 * (dD1_ds.dot(D2_half) - tau);
        const Vector3d N_half = N1 * D1_half + N2 * D2_half + N3 * D3_half;

        Eigen::Map<Vector3d> F_curr(&F_curr_node_vals[k * NDIM]);
        Eigen::Map<Vector3d> F_next(&F_next_node_vals[k * NDIM]);
        F_curr = F_half;
        F_next = -F_half;

        Eigen::Map<Vector3d> N_curr(&N_curr_node_vals[k * NDIM]);
        Eigen::Map<Vector3d> N_next(&N_next_node_vals[k * NDIM]);
        N_curr = N_half + 0.5 * (X_next - X).cross(F_half);
        N_next = -N_half + 0.5 * (X_next - X).cross(F_half);
    }

    ierr = VecRestoreArray(D_vec, &D_vals);
    IBTK_CHKERRQ(ierr);

    ierr = VecRestoreArray(D_next_vec, &D_next_vals);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&D_next_vec);
    IBTK_CHKERRQ(ierr);

    ierr = VecRestoreArray(X_vec, &X_vals);
    IBTK_CHKERRQ(ierr);

    ierr = VecRestoreArray(X_next_vec, &X_next_vals);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&X_next_vec);
    IBTK_CHKERRQ(ierr);

    Vec F_vec = F_data->getVec();
    Vec N_vec = N_data->getVec();
    if (local_sz > 0)
    {
        ierr = VecSetValuesBlocked(F_vec,
                                   static_cast<int>(petsc_curr_node_idxs.size()),
                                   &petsc_curr_node_idxs[0],
                                   &F_curr_node_vals[0],
                                   ADD_VALUES);
        IBTK_CHKERRQ(ierr);
        ierr = VecSetValuesBlocked(F_vec,
                                   static_cast<int>(petsc_next_node_idxs.size()),
                                   &petsc_next_node_idxs[0],
                                   &F_next_node_vals[0],
                                   ADD_VALUES);
        IBTK_CHKERRQ(ierr);
        ierr = VecSetValuesBlocked(N_vec,
                                   static_cast<int>(petsc_curr_node_idxs.size()),
                                   &petsc_curr_node_idxs[0],
                                   &N_curr_node_vals[0],
                                   ADD_VALUES);
        IBTK_CHKERRQ(ierr);
        ierr = VecSetValuesBlocked(N_vec,
                                   static_cast<int>(petsc_next_node_idxs.size()),
                                   &petsc_next_node_idxs[0],
                                   &N_next_node_vals[0],
                                   ADD_VALUES);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(F_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyBegin(N_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(N_vec);
    IBTK_CHKERRQ(ierr);

    IBAMR_TIMER_STOP(t_compute_lagrangian_force_and_torque);
    return;
} // computeLagrangianForceAndTorque

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBKirchhoffRodForceGen::getFromInput(Pointer<Database> db)
{
    if (db)
    {
        // intentionally blank
    }
    return;
} // getFromInput

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
