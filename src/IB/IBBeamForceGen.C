// Filename: IBBeamForceGen.C
// Created on 22 Mar 2007 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
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

#include "IBBeamForceGen.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/IBBeamForceSpec.h>
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/LNodeIndexSetData.h>
#include <ibtk/PETScVecOps.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

// BLITZ++ INCLUDES
#include <blitz/array.h>

// C++ STDLIB INCLUDES
#include <numeric>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_compute_lagrangian_force;
static Timer* t_compute_lagrangian_force_jacobian;
static Timer* t_compute_lagrangian_force_jacobian_nonzero_structure;
static Timer* t_initialize_level_data;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBBeamForceGen::IBBeamForceGen(
    Pointer<Database> input_db)
    : d_D_next_mats(),
      d_D_prev_mats(),
      d_petsc_mastr_node_idxs(),
      d_petsc_next_node_idxs(),
      d_petsc_prev_node_idxs(),
      d_bend_rigidities(),
      d_mesh_dependent_curvatures(),
      d_is_initialized()
{
    // Initialize object with data read from the input database.
    getFromInput(input_db);

    // Setup Timers.
    IBAMR_DO_ONCE(
        t_compute_lagrangian_force                            = TimerManager::getManager()->getTimer("IBAMR::IBBeamForceGen::computeLagrangianForce()");
        t_compute_lagrangian_force_jacobian                   = TimerManager::getManager()->getTimer("IBAMR::IBBeamForceGen::computeLagrangianForceJacobian()");
        t_compute_lagrangian_force_jacobian_nonzero_structure = TimerManager::getManager()->getTimer("IBAMR::IBBeamForceGen::computeLagrangianForceJacobianNonzeroStructure()");
        t_initialize_level_data                               = TimerManager::getManager()->getTimer("IBAMR::IBBeamForceGen::initializeLevelData()");
                  );
    return;
}// IBBeamForceGen

IBBeamForceGen::~IBBeamForceGen()
{
    int ierr;
    for (std::vector<Mat>::iterator it = d_D_next_mats.begin(); it != d_D_next_mats.end(); ++it)
    {
        if (*it)
        {
            ierr = MatDestroy(*it);  IBTK_CHKERRQ(ierr);
        }
    }
    for (std::vector<Mat>::iterator it = d_D_prev_mats.begin(); it != d_D_prev_mats.end(); ++it)
    {
        if (*it)
        {
            ierr = MatDestroy(*it);  IBTK_CHKERRQ(ierr);
        }
    }
    return;
}// ~IBBeamForceGen

void
IBBeamForceGen::initializeLevelData(
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool initial_time,
    LDataManager* const lag_manager)
{
    if (!lag_manager->levelContainsLagrangianData(level_number)) return;

    t_initialize_level_data->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
#endif
    (void) init_data_time;
    (void) initial_time;

    int ierr;

    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Resize the vectors corresponding to data individually maintained for
    // separate levels of the patch hierarchy.
    const int level_num = level->getLevelNumber();
    const int new_size = std::max(level_num+1, int(d_is_initialized.size()));

    d_D_next_mats.resize(new_size);
    d_D_prev_mats.resize(new_size);
    d_petsc_mastr_node_idxs.resize(new_size);
    d_petsc_next_node_idxs.resize(new_size);
    d_petsc_prev_node_idxs.resize(new_size);
    d_bend_rigidities.resize(new_size);
    d_mesh_dependent_curvatures.resize(new_size);
    d_is_initialized.resize(new_size, false);

    Mat& D_next_mat = d_D_next_mats[level_num];
    Mat& D_prev_mat = d_D_prev_mats[level_num];
    std::vector<int>& petsc_mastr_node_idxs = d_petsc_mastr_node_idxs[level_num];
    std::vector<int>& petsc_next_node_idxs = d_petsc_next_node_idxs[level_num];
    std::vector<int>& petsc_prev_node_idxs = d_petsc_prev_node_idxs[level_num];
    std::vector<double>& bend_rigidities = d_bend_rigidities[level_num];
    std::vector<std::vector<double> >& mesh_dependent_curvatures = d_mesh_dependent_curvatures[level_num];

    if (D_next_mat)
    {
        ierr = MatDestroy(D_next_mat);  IBTK_CHKERRQ(ierr);
    }
    if (D_prev_mat)
    {
        ierr = MatDestroy(D_prev_mat);  IBTK_CHKERRQ(ierr);
    }
    petsc_mastr_node_idxs.clear();
    petsc_next_node_idxs.clear();
    petsc_prev_node_idxs.clear();
    bend_rigidities.clear();
    mesh_dependent_curvatures.clear();

    // The patch data descriptor index for the LNodeIndexSetData.
    const int lag_node_index_idx = lag_manager->getLNodeIndexPatchDescriptorIndex();

    // Determine the "next" and "prev" node indices for all beams associated
    // with the present MPI process.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<LNodeIndexSetData> idx_data = patch->getPatchData(lag_node_index_idx);
        for (LNodeIndexSetData::DataIterator it = idx_data->data_begin(patch_box); it != idx_data->data_end(); ++it)
        {
            const LNodeIndex& node_idx = *it;
            const Pointer<IBBeamForceSpec> force_spec = node_idx.getNodeData<IBBeamForceSpec>();
            if (!force_spec.isNull())
            {
                const int& mastr_idx = node_idx.getLagrangianIndex();
                const unsigned num_beams = force_spec->getNumberOfBeams();
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(mastr_idx == force_spec->getMasterNodeIndex());
#endif
                const std::vector<std::pair<int,int> >& nghbrs = force_spec->getNeighborNodeIndices();
                const std::vector<double>& bend = force_spec->getBendingRigidities();
                const std::vector<std::vector<double> >& curv = force_spec->getMeshDependentCurvatures();
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(num_beams == nghbrs.size());
                TBOX_ASSERT(num_beams == bend.size());
                TBOX_ASSERT(num_beams == curv.size());
#endif
                for (unsigned k = 0; k < num_beams; ++k)
                {
                    petsc_mastr_node_idxs.push_back(mastr_idx);
                    petsc_next_node_idxs.push_back(nghbrs[k].first );
                    petsc_prev_node_idxs.push_back(nghbrs[k].second);
                    bend_rigidities.push_back(bend[k]);
                    mesh_dependent_curvatures.push_back(curv[k]);
                }
            }
        }
    }

    // Map the Lagrangian node indices to the PETSc indices corresponding to the
    // present data distribution.
    lag_manager->mapLagrangianToPETSc(petsc_mastr_node_idxs, level_num);
    lag_manager->mapLagrangianToPETSc(petsc_next_node_idxs, level_num);
    lag_manager->mapLagrangianToPETSc(petsc_prev_node_idxs, level_num);

    // Determine the global node offset and the number of local nodes.
    const int global_node_offset = lag_manager->getGlobalNodeOffset(level_num);
    const int num_local_nodes = lag_manager->getNumberOfLocalNodes(level_num);

    // Determine the non-zero structure for the matrices.
    const int local_sz = petsc_mastr_node_idxs.size();

    std::vector<int> next_d_nz(local_sz,1), next_o_nz(local_sz,0);
    for (int k = 0; k < local_sz; ++k)
    {
        const int& next_idx = petsc_next_node_idxs[k];
        if (next_idx >= global_node_offset &&
            next_idx <  global_node_offset+num_local_nodes)
        {
            ++next_d_nz[k]; // a "local"    next index
        }
        else
        {
            ++next_o_nz[k]; // a "nonlocal" next index
        }
    }

    std::vector<int> prev_d_nz(local_sz,1), prev_o_nz(local_sz,0);
    for (int k = 0; k < local_sz; ++k)
    {
        const int& prev_idx = petsc_prev_node_idxs[k];
        if (prev_idx >= global_node_offset &&
            prev_idx <  global_node_offset+num_local_nodes)
        {
            ++prev_d_nz[k]; // a "local"    prev index
        }
        else
        {
            ++prev_o_nz[k]; // a "nonlocal" prev index
        }
    }

    // Create new MPI block AIJ matrices and set the values of the non-zero
    // entries.
    ierr = MatCreateMPIBAIJ(PETSC_COMM_WORLD,
                            NDIM, NDIM*local_sz, NDIM*num_local_nodes,
                            PETSC_DETERMINE, PETSC_DETERMINE,
                            PETSC_DEFAULT, local_sz > 0 ? &next_d_nz[0] : PETSC_NULL,
                            PETSC_DEFAULT, local_sz > 0 ? &next_o_nz[0] : PETSC_NULL,
                            &D_next_mat);  IBTK_CHKERRQ(ierr);

    ierr = MatCreateMPIBAIJ(PETSC_COMM_WORLD,
                            NDIM, NDIM*local_sz, NDIM*num_local_nodes,
                            PETSC_DETERMINE, PETSC_DETERMINE,
                            PETSC_DEFAULT, local_sz > 0 ? &prev_d_nz[0] : PETSC_NULL,
                            PETSC_DEFAULT, local_sz > 0 ? &prev_o_nz[0] : PETSC_NULL,
                            &D_prev_mat);  IBTK_CHKERRQ(ierr);

    blitz::Array<double,2> mastr_vals(NDIM,NDIM);  mastr_vals = 0.0;
    blitz::Array<double,2> slave_vals(NDIM,NDIM);  slave_vals = 0.0;
    for (int d = 0; d < NDIM; ++d)
    {
        mastr_vals(d,d) = -1.0;
        slave_vals(d,d) = +1.0;
    }

    int i_offset;

    ierr = MatGetOwnershipRange(D_next_mat, &i_offset, PETSC_NULL);
    IBTK_CHKERRQ(ierr);
    i_offset /= NDIM;

    for (int k = 0; k < local_sz; ++k)
    {
        int i = i_offset + k;
        int j_mastr = petsc_mastr_node_idxs[k];
        int j_slave = petsc_next_node_idxs[k];
        ierr = MatSetValuesBlocked(D_next_mat,1,&i,1,&j_mastr,mastr_vals.data(),INSERT_VALUES);  IBTK_CHKERRQ(ierr);
        ierr = MatSetValuesBlocked(D_next_mat,1,&i,1,&j_slave,slave_vals.data(),INSERT_VALUES);  IBTK_CHKERRQ(ierr);
    }

    ierr = MatGetOwnershipRange(D_prev_mat, &i_offset, PETSC_NULL);
    IBTK_CHKERRQ(ierr);
    i_offset /= NDIM;

    for (int k = 0; k < local_sz; ++k)
    {
        int i = i_offset + k;
        int j_mastr = petsc_mastr_node_idxs[k];
        int j_slave = petsc_prev_node_idxs[k];
        ierr = MatSetValuesBlocked(D_prev_mat,1,&i,1,&j_mastr,mastr_vals.data(),INSERT_VALUES);  IBTK_CHKERRQ(ierr);
        ierr = MatSetValuesBlocked(D_prev_mat,1,&i,1,&j_slave,slave_vals.data(),INSERT_VALUES);  IBTK_CHKERRQ(ierr);
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(D_next_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyBegin(D_prev_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(D_next_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(D_prev_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);

    // Indicate that the level data has been initialized.
    d_is_initialized[level_num] = true;

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
IBBeamForceGen::computeLagrangianForce(
    Pointer<LData> F_data,
    Pointer<LData> X_data,
    Pointer<LData> U_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    if (!lag_manager->levelContainsLagrangianData(level_number)) return;

    t_compute_lagrangian_force->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(level_number < int(d_is_initialized.size()));
    TBOX_ASSERT(d_is_initialized[level_number]);
#endif

    int ierr;

    // Create an appropriately sized temporary vector to store the node
    // displacements.
    int i_start, i_stop;
    ierr = MatGetOwnershipRange(d_D_next_mats[level_number], &i_start, &i_stop);
    IBTK_CHKERRQ(ierr);

    Vec D_next_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop-i_start, PETSC_DECIDE, &D_next_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecSetBlockSize(D_next_vec, NDIM);                                          IBTK_CHKERRQ(ierr);

    Vec D_prev_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop-i_start, PETSC_DECIDE, &D_prev_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecSetBlockSize(D_prev_vec, NDIM);                                          IBTK_CHKERRQ(ierr);

    // Compute the node displacements.
    ierr = MatMult(d_D_next_mats[level_number], X_data->getVec(), D_next_vec);  IBTK_CHKERRQ(ierr);
    ierr = MatMult(d_D_prev_mats[level_number], X_data->getVec(), D_prev_vec);  IBTK_CHKERRQ(ierr);

    // Compute the beam forces acting on the nodes of the Lagrangian mesh.
    double* D_next_vals;
    ierr = VecGetArray(D_next_vec, &D_next_vals);  IBTK_CHKERRQ(ierr);

    double* D_prev_vals;
    ierr = VecGetArray(D_prev_vec, &D_prev_vals);  IBTK_CHKERRQ(ierr);

    std::vector<int>& petsc_mastr_node_idxs = d_petsc_mastr_node_idxs[level_number];
    std::vector<int>& petsc_next_node_idxs = d_petsc_next_node_idxs[level_number];
    std::vector<int>& petsc_prev_node_idxs = d_petsc_prev_node_idxs[level_number];
    std::vector<double>& bend_rigidities = d_bend_rigidities[level_number];
    std::vector<std::vector<double> >& mesh_dependent_curvatures = d_mesh_dependent_curvatures[level_number];

    const int local_sz = petsc_mastr_node_idxs.size();
    std::vector<double> F_mastr_node_vals(NDIM*local_sz,0.0);
    std::vector<double> F_nghbr_node_vals(NDIM*local_sz,0.0);

    for (int k = 0; k < local_sz; ++k)
    {
        // Compute the forces applied by the beam to the "master" and "slave"
        // nodes.
        double* const F_mastr_node = &F_mastr_node_vals[k*NDIM];
        double* const F_nghbr_node = &F_nghbr_node_vals[k*NDIM];
        const double* const D_next = &D_next_vals[k*NDIM];
        const double* const D_prev = &D_prev_vals[k*NDIM];
        const double& bend = bend_rigidities[k];
        const std::vector<double>& curv = mesh_dependent_curvatures[k];
        for (int d = 0; d < NDIM; ++d)
        {
            F_mastr_node[d] = +2.0*bend*(D_next[d]+D_prev[d]-curv[d]);
            F_nghbr_node[d] = -1.0*bend*(D_next[d]+D_prev[d]-curv[d]);
        }
    }

    ierr = VecRestoreArray(D_next_vec, &D_next_vals);  IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(D_next_vec);                    IBTK_CHKERRQ(ierr);

    ierr = VecRestoreArray(D_prev_vec, &D_prev_vals);  IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(D_prev_vec);                    IBTK_CHKERRQ(ierr);

    Vec F_vec = F_data->getVec();
#if 0
    ierr = VecSetValuesBlocked(F_vec,
                               petsc_mastr_node_idxs.size(),
                               !petsc_mastr_node_idxs.empty() ? &petsc_mastr_node_idxs[0] : PETSC_NULL,
                               !petsc_mastr_node_idxs.empty() ? &    F_mastr_node_vals[0] : PETSC_NULL,
                               ADD_VALUES);  IBTK_CHKERRQ(ierr);
    ierr = VecSetValuesBlocked(F_vec,
                               petsc_next_node_idxs.size(),
                               !petsc_next_node_idxs.empty() ? &petsc_next_node_idxs[0] : PETSC_NULL,
                               !petsc_next_node_idxs.empty() ? &   F_nghbr_node_vals[0] : PETSC_NULL,
                               ADD_VALUES);  IBTK_CHKERRQ(ierr);
    ierr = VecSetValuesBlocked(F_vec,
                               petsc_prev_node_idxs.size(),
                               !petsc_prev_node_idxs.empty() ? &petsc_prev_node_idxs[0] : PETSC_NULL,
                               !petsc_prev_node_idxs.empty() ? &   F_nghbr_node_vals[0] : PETSC_NULL,
                               ADD_VALUES);  IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F_vec);    IBTK_CHKERRQ(ierr);
#else
    ierr = PETScVecOps::VecSetValuesBlocked(F_vec,
                                            petsc_mastr_node_idxs.size(),
                                            !petsc_mastr_node_idxs.empty() ? &petsc_mastr_node_idxs[0] : PETSC_NULL,
                                            !petsc_mastr_node_idxs.empty() ? &    F_mastr_node_vals[0] : PETSC_NULL,
                                            ADD_VALUES);  IBTK_CHKERRQ(ierr);
    ierr = PETScVecOps::VecSetValuesBlocked(F_vec,
                                            petsc_next_node_idxs.size(),
                                            !petsc_next_node_idxs.empty() ? &petsc_next_node_idxs[0] : PETSC_NULL,
                                            !petsc_next_node_idxs.empty() ? &   F_nghbr_node_vals[0] : PETSC_NULL,
                                            ADD_VALUES);  IBTK_CHKERRQ(ierr);
    ierr = PETScVecOps::VecSetValuesBlocked(F_vec,
                                            petsc_prev_node_idxs.size(),
                                            !petsc_prev_node_idxs.empty() ? &petsc_prev_node_idxs[0] : PETSC_NULL,
                                            !petsc_prev_node_idxs.empty() ? &   F_nghbr_node_vals[0] : PETSC_NULL,
                                            ADD_VALUES);  IBTK_CHKERRQ(ierr);
    ierr = PETScVecOps::VecAssemblyBegin(F_vec);  IBTK_CHKERRQ(ierr);
    ierr = PETScVecOps::VecAssemblyEnd(F_vec);    IBTK_CHKERRQ(ierr);
#endif
    t_compute_lagrangian_force->stop();
    return;
}// computeLagrangianForce

void
IBBeamForceGen::computeLagrangianForceJacobianNonzeroStructure(
    std::vector<int>& d_nnz,
    std::vector<int>& o_nnz,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    if (!lag_manager->levelContainsLagrangianData(level_number)) return;

    t_compute_lagrangian_force_jacobian_nonzero_structure->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(level_number < int(d_is_initialized.size()));
    TBOX_ASSERT(d_is_initialized[level_number]);
#endif

    int ierr;

    // Look up the cached connectivity information.
    std::vector<int>& petsc_mastr_node_idxs = d_petsc_mastr_node_idxs[level_number];
    std::vector<int>& petsc_next_node_idxs = d_petsc_next_node_idxs[level_number];
    std::vector<int>& petsc_prev_node_idxs = d_petsc_prev_node_idxs[level_number];
    const int local_sz = petsc_mastr_node_idxs.size();

    // Determine the global node offset and the number of local nodes.
    const int global_node_offset = lag_manager->getGlobalNodeOffset(level_number);
    const int num_local_nodes = lag_manager->getNumberOfLocalNodes(level_number);

    // Determine the non-zero structure for the matrix used to store the
    // Jacobian of the force.
    Vec d_nnz_vec, o_nnz_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, num_local_nodes, PETSC_DETERMINE, &d_nnz_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, num_local_nodes, PETSC_DETERMINE, &o_nnz_vec);  IBTK_CHKERRQ(ierr);

    ierr = VecSet(d_nnz_vec, 1.0);  IBTK_CHKERRQ(ierr);
    ierr = VecSet(o_nnz_vec, 0.0);  IBTK_CHKERRQ(ierr);

    for (int k = 0; k < local_sz; ++k)
    {
        const int& mastr_idx = petsc_mastr_node_idxs[k];
        const int& next_idx = petsc_next_node_idxs[k];
        const int& prev_idx = petsc_prev_node_idxs[k];

        const bool next_is_local = (next_idx >= global_node_offset &&
                                    next_idx <  global_node_offset + num_local_nodes);
        const bool prev_is_local = (prev_idx >= global_node_offset &&
                                    prev_idx <  global_node_offset + num_local_nodes);
        if (next_is_local && prev_is_local)
        {
            static const int d_N = 3;
            const int d_idxs[d_N] = { mastr_idx , next_idx , prev_idx };
            const double d_vals[d_N] = { 2.0 , 2.0 , 2.0 };
            ierr = VecSetValues(d_nnz_vec, d_N, d_idxs, d_vals, ADD_VALUES);  IBTK_CHKERRQ(ierr);
        }
        else if (next_is_local)
        {
            static const int d_N = 2;
            const int d_idxs[d_N] = { mastr_idx , next_idx };
            const double d_vals[d_N] = { 1.0 , 1.0 };
            ierr = VecSetValues(d_nnz_vec, d_N, d_idxs, d_vals, ADD_VALUES);  IBTK_CHKERRQ(ierr);

            static const int o_N = 3;
            const int o_idxs[o_N] = { mastr_idx , next_idx , prev_idx };
            const double o_vals[o_N] = { 1.0 , 1.0 , 2.0 };
            ierr = VecSetValues(o_nnz_vec, o_N, o_idxs, o_vals, ADD_VALUES);  IBTK_CHKERRQ(ierr);
        }
        else if (prev_is_local)
        {
            static const int d_N = 2;
            const int d_idxs[d_N] = { mastr_idx , prev_idx };
            const double d_vals[d_N] = { 1.0 , 1.0 };
            ierr = VecSetValues(d_nnz_vec, d_N, d_idxs, d_vals, ADD_VALUES);  IBTK_CHKERRQ(ierr);

            static const int o_N = 3;
            const int o_idxs[o_N] = { mastr_idx , next_idx , prev_idx };
            const double o_vals[o_N] = { 1.0 , 2.0 , 1.0 };
            ierr = VecSetValues(o_nnz_vec, o_N, o_idxs, o_vals, ADD_VALUES);  IBTK_CHKERRQ(ierr);
        }
        else
        {
            // NOTE: Here, we are slightly pessimistically allocating space both
            // for the case that the previous and next nodes are on different
            // processors, and for the case that the previous and next nodes are
            // on the same processor.
            static const int d_N = 2;
            const int d_idxs[d_N] = { next_idx , prev_idx };
            const double d_vals[d_N] = { 2.0 , 2.0 };
            ierr = VecSetValues(d_nnz_vec, d_N, d_idxs, d_vals, ADD_VALUES);  IBTK_CHKERRQ(ierr);

            static const int o_N = 3;
            const int o_idxs[o_N] = { mastr_idx , next_idx , prev_idx };
            const double o_vals[o_N] = { 2.0 , 2.0 , 2.0 };
            ierr = VecSetValues(o_nnz_vec, o_N, o_idxs, o_vals, ADD_VALUES);  IBTK_CHKERRQ(ierr);
        }
    }

    ierr = VecAssemblyBegin(d_nnz_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyBegin(o_nnz_vec);  IBTK_CHKERRQ(ierr);

    ierr = VecAssemblyEnd(d_nnz_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(o_nnz_vec);  IBTK_CHKERRQ(ierr);

    double* d_nnz_vec_arr;
    double* o_nnz_vec_arr;

    ierr = VecGetArray(d_nnz_vec, &d_nnz_vec_arr);  IBTK_CHKERRQ(ierr);
    ierr = VecGetArray(o_nnz_vec, &o_nnz_vec_arr);  IBTK_CHKERRQ(ierr);

    for (int k = 0; k < num_local_nodes; ++k)
    {
        d_nnz[k] += int(d_nnz_vec_arr[k]);
        o_nnz[k] += int(o_nnz_vec_arr[k]);
    }

    ierr = VecRestoreArray(d_nnz_vec, &d_nnz_vec_arr);  IBTK_CHKERRQ(ierr);
    ierr = VecRestoreArray(o_nnz_vec, &o_nnz_vec_arr);  IBTK_CHKERRQ(ierr);

    ierr = VecDestroy(d_nnz_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(o_nnz_vec);  IBTK_CHKERRQ(ierr);

    t_compute_lagrangian_force_jacobian_nonzero_structure->stop();
    return;
}// computeLagrangianForceJacobianNonzeroStructure

void
IBBeamForceGen::computeLagrangianForceJacobian(
    Mat& J_mat,
    MatAssemblyType assembly_type,
    const double X_coef,
    Pointer<LData> X_data,
    const double U_coef,
    Pointer<LData> U_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    if (!lag_manager->levelContainsLagrangianData(level_number)) return;

    t_compute_lagrangian_force_jacobian->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(level_number < int(d_is_initialized.size()));
    TBOX_ASSERT(d_is_initialized[level_number]);
#endif

    int ierr;

    // Look up the cached connectivity information.
    std::vector<int>& petsc_mastr_node_idxs = d_petsc_mastr_node_idxs[level_number];
    std::vector<int>& petsc_next_node_idxs = d_petsc_next_node_idxs[level_number];
    std::vector<int>& petsc_prev_node_idxs = d_petsc_prev_node_idxs[level_number];
    const int local_sz = petsc_mastr_node_idxs.size();

    // Compute the matrix elements and add them to the Jacobian matrix.
    std::vector<double>& bend_rigidities = d_bend_rigidities[level_number];
    blitz::Array<double,2> dF_dX(NDIM,NDIM);  dF_dX = 0.0;
    for (int k = 0; k < local_sz; ++k)
    {
        const int& petsc_mastr_idx = petsc_mastr_node_idxs[k];
        const int& petsc_next_idx = petsc_next_node_idxs[k];
        const int& petsc_prev_idx = petsc_prev_node_idxs[k];
        const double& bend = bend_rigidities[k];

        for (int alpha = 0; alpha < NDIM; ++alpha)
        {
            dF_dX(alpha,alpha) = -1.0*bend*X_coef;
        }
        ierr = MatSetValuesBlocked(J_mat,1,&petsc_prev_idx,1,&petsc_prev_idx,dF_dX.data(),ADD_VALUES);  IBTK_CHKERRQ(ierr);
        ierr = MatSetValuesBlocked(J_mat,1,&petsc_prev_idx,1,&petsc_next_idx,dF_dX.data(),ADD_VALUES);  IBTK_CHKERRQ(ierr);
        ierr = MatSetValuesBlocked(J_mat,1,&petsc_next_idx,1,&petsc_prev_idx,dF_dX.data(),ADD_VALUES);  IBTK_CHKERRQ(ierr);
        ierr = MatSetValuesBlocked(J_mat,1,&petsc_next_idx,1,&petsc_next_idx,dF_dX.data(),ADD_VALUES);  IBTK_CHKERRQ(ierr);

        for (int alpha = 0; alpha < NDIM; ++alpha)
        {
            dF_dX(alpha,alpha) = +2.0*bend*X_coef;
        }
        ierr = MatSetValuesBlocked(J_mat,1,&petsc_prev_idx,1,&petsc_mastr_idx,dF_dX.data(),ADD_VALUES);  IBTK_CHKERRQ(ierr);
        ierr = MatSetValuesBlocked(J_mat,1,&petsc_next_idx,1,&petsc_mastr_idx,dF_dX.data(),ADD_VALUES);  IBTK_CHKERRQ(ierr);
        ierr = MatSetValuesBlocked(J_mat,1,&petsc_mastr_idx,1,&petsc_prev_idx,dF_dX.data(),ADD_VALUES);  IBTK_CHKERRQ(ierr);
        ierr = MatSetValuesBlocked(J_mat,1,&petsc_mastr_idx,1,&petsc_next_idx,dF_dX.data(),ADD_VALUES);  IBTK_CHKERRQ(ierr);

        for (int alpha = 0; alpha < NDIM; ++alpha)
        {
            dF_dX(alpha,alpha) = -4.0*bend*X_coef;
        }
        ierr = MatSetValuesBlocked(J_mat,1,&petsc_mastr_idx,1,&petsc_mastr_idx,dF_dX.data(),ADD_VALUES);  IBTK_CHKERRQ(ierr);
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(J_mat, assembly_type);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(J_mat, assembly_type);    IBTK_CHKERRQ(ierr);

    t_compute_lagrangian_force_jacobian->stop();
    return;
}// computeLagrangianForceJacobian

double
IBBeamForceGen::computeLagrangianEnergy(
    Pointer<LData> X_data,
    Pointer<LData> U_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    TBOX_WARNING("IBBeamForceGen::computeLagrangianEnergy():\n"
                 << "  unimplemented; returning 0.0." << std::endl);
    return 0.0;
}// computeLagrangianEnergy

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBBeamForceGen::getFromInput(
    Pointer<Database> db)
{
    if (!db.isNull())
    {
        // intentionally blank
    }
    return;
}// getFromInput

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBBeamForceGen>;

//////////////////////////////////////////////////////////////////////////////
