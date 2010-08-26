// Filename: IBSpringForceGen.C
// Created on 14 Jul 2004 by Boyce Griffith
//
// Copyright (c) 2002-2010 Boyce Griffith
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "IBSpringForceGen.h"

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
#include <ibamr/IBSpringForceSpec.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/LNodeIndexData.h>
#include <ibtk/PETScVecOps.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

// BLITZ++ INCLUDES
#include <blitz/array.h>
#include <blitz/tinyvec.h>

// C++ STDLIB INCLUDES
#include <algorithm>
#include <functional>
#include <numeric>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Pointer<Timer> t_compute_lagrangian_force;
static Pointer<Timer> t_compute_lagrangian_force_jacobian;
static Pointer<Timer> t_compute_lagrangian_force_jacobian_nonzero_structure;
static Pointer<Timer> t_compute_lagrangian_energy;
static Pointer<Timer> t_initialize_level_data;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBSpringForceGen::IBSpringForceGen(
    Pointer<Database> input_db)
    : d_D_mats(),
      d_lag_mastr_node_idxs(),
      d_lag_slave_node_idxs(),
      d_petsc_mastr_node_idxs(),
      d_petsc_slave_node_idxs(),
      d_force_fcn_idxs(),
      d_stiffnesses(),
      d_rest_lengths(),
      d_is_initialized(),
      d_force_fcn_map()
{
    // Initialize object with data read from the input database.
    getFromInput(input_db);

    // Setup the default force generation functions.
    registerSpringForceFunction(0, &IBAMR::default_linear_spring_force);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_compute_lagrangian_force                            = TimerManager::getManager()->getTimer("IBAMR::IBSpringForceGen::computeLagrangianForce()");
        t_compute_lagrangian_force_jacobian                   = TimerManager::getManager()->getTimer("IBAMR::IBSpringForceGen::computeLagrangianForceJacobian()");
        t_compute_lagrangian_force_jacobian_nonzero_structure = TimerManager::getManager()->getTimer("IBAMR::IBSpringForceGen::computeLagrangianForceJacobianNonzeroStructure()");
        t_initialize_level_data                               = TimerManager::getManager()->getTimer("IBAMR::IBSpringForceGen::initializeLevelData()");
        t_compute_lagrangian_energy                           = TimerManager::getManager()->getTimer("IBAMR::IBSpringForceGen::computeLagrangianEnergy()");
        timers_need_init = false;
    }
    return;
}// IBSpringForceGen

IBSpringForceGen::~IBSpringForceGen()
{
    int ierr;
    for (std::vector<Mat>::iterator it = d_D_mats.begin(); it != d_D_mats.end(); ++it)
    {
        if (*it)
        {
            ierr = MatDestroy(*it);  IBTK_CHKERRQ(ierr);
        }
    }
    return;
}// ~IBSpringForceGen

void
IBSpringForceGen::registerSpringForceFunction(
    const int force_fcn_index,
    void (*force_fcn)(double F[NDIM], const double D[NDIM], const double& stf, const double& rst, const int& lag_mastr_idx, const int& lag_slave_idx))
{
    d_force_fcn_map[force_fcn_index] = force_fcn;
    return;
}// registerSpringForceFunction

void
IBSpringForceGen::initializeLevelData(
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
    const int new_size = std::max(level_number+1, int(d_is_initialized.size()));

    d_D_mats.resize(new_size);
    d_lag_mastr_node_idxs.resize(new_size);
    d_lag_slave_node_idxs.resize(new_size);
    d_petsc_mastr_node_idxs.resize(new_size);
    d_petsc_slave_node_idxs.resize(new_size);
    d_force_fcn_idxs.resize(new_size);
    d_stiffnesses.resize(new_size);
    d_rest_lengths.resize(new_size);
    d_is_initialized.resize(new_size, false);

    Mat& D_mat = d_D_mats[level_number];
    std::vector<int>& lag_mastr_node_idxs = d_lag_mastr_node_idxs[level_number];
    std::vector<int>& lag_slave_node_idxs = d_lag_slave_node_idxs[level_number];
    std::vector<int>& petsc_mastr_node_idxs = d_petsc_mastr_node_idxs[level_number];
    std::vector<int>& petsc_slave_node_idxs = d_petsc_slave_node_idxs[level_number];
    std::vector<int>& force_fcn_idxs = d_force_fcn_idxs[level_number];
    std::vector<double>& stiffnesses = d_stiffnesses[level_number];
    std::vector<double>& rest_lengths = d_rest_lengths[level_number];

    if (D_mat)
    {
        ierr = MatDestroy(D_mat);  IBTK_CHKERRQ(ierr);
    }
    lag_mastr_node_idxs.clear();
    lag_slave_node_idxs.clear();
    petsc_mastr_node_idxs.clear();
    petsc_slave_node_idxs.clear();
    force_fcn_idxs.clear();
    stiffnesses.clear();
    rest_lengths.clear();

    // The patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = lag_manager->getLNodeIndexPatchDescriptorIndex();

    // Determine the "master" and "slave" node indices for all springs
    // associated with the present MPI process.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
        for (LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
             it != idx_data->lnode_index_end(); ++it)
        {
            const LNodeIndex& node_idx = *it;
            Pointer<IBSpringForceSpec> force_spec = node_idx.getStashData<IBSpringForceSpec>();
            if (!force_spec.isNull())
            {
                const int& mastr_idx = node_idx.getLagrangianIndex();
                const unsigned num_springs = force_spec->getNumberOfSprings();
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(mastr_idx == force_spec->getMasterNodeIndex());
#endif
                const std::vector<int>& slv = force_spec->getSlaveNodeIndices();
                const std::vector<int>& fcn = force_spec->getForceFunctionIndices();
                const std::vector<double>& stf = force_spec->getStiffnesses();
                const std::vector<double>& rst = force_spec->getRestingLengths();
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(num_springs == slv.size());
                TBOX_ASSERT(num_springs == fcn.size());
                TBOX_ASSERT(num_springs == stf.size());
                TBOX_ASSERT(num_springs == rst.size());
#endif
                if (num_springs > 0)
                {
                    lag_mastr_node_idxs.insert(lag_mastr_node_idxs.end(), num_springs, mastr_idx);
                    lag_slave_node_idxs.insert(lag_slave_node_idxs.end(), slv.begin(), slv.end());
                    force_fcn_idxs.insert(force_fcn_idxs.end(), fcn.begin(), fcn.end());
                    stiffnesses   .insert(stiffnesses   .end(), stf.begin(), stf.end());
                    rest_lengths  .insert(rest_lengths  .end(), rst.begin(), rst.end());
                }
            }
        }
    }

    // Map the Lagrangian master/slave node indices to the PETSc indices
    // corresponding to the present data distribution.
    petsc_mastr_node_idxs = lag_mastr_node_idxs;
    petsc_slave_node_idxs = lag_slave_node_idxs;
    lag_manager->mapLagrangianToPETSc(petsc_mastr_node_idxs, level_number);
    lag_manager->mapLagrangianToPETSc(petsc_slave_node_idxs, level_number);

    // Determine the global node offset and the number of local nodes.
    const int global_node_offset = lag_manager->getGlobalNodeOffset(level_number);
    const int num_local_nodes = lag_manager->getNumberOfLocalNodes(level_number);

    // Determine the non-zero structure for the matrix used to compute nodal
    // displacements.
    const int local_sz = petsc_mastr_node_idxs.size();
    std::vector<int> d_nnz(local_sz,1), o_nnz(local_sz,0);
    for (int k = 0; k < local_sz; ++k)
    {
        const int& slave_idx = petsc_slave_node_idxs[k];
        if (slave_idx >= global_node_offset &&
            slave_idx <  global_node_offset+num_local_nodes)
        {
            ++d_nnz[k]; // a "local"    slave index
        }
        else
        {
            ++o_nnz[k]; // a "nonlocal" slave index
        }
    }

    // Create a new MPI block AIJ matrix used to compute nodal displacements and
    // set the values of the non-zero entries.
    ierr = MatCreateMPIBAIJ(PETSC_COMM_WORLD,
                            NDIM, NDIM*local_sz, NDIM*num_local_nodes,
                            PETSC_DETERMINE, PETSC_DETERMINE,
                            PETSC_DEFAULT, local_sz > 0 ? &d_nnz[0] : PETSC_NULL,
                            PETSC_DEFAULT, local_sz > 0 ? &o_nnz[0] : PETSC_NULL,
                            &D_mat);  IBTK_CHKERRQ(ierr);

    blitz::Array<double,2> mastr_vals(NDIM,NDIM);  mastr_vals = 0.0;
    blitz::Array<double,2> slave_vals(NDIM,NDIM);  slave_vals = 0.0;
    for (int d = 0; d < NDIM; ++d)
    {
        mastr_vals(d,d) = -1.0;
        slave_vals(d,d) = +1.0;
    }

    int i_offset;
    ierr = MatGetOwnershipRange(D_mat, &i_offset, PETSC_NULL);
    IBTK_CHKERRQ(ierr);
    i_offset /= NDIM;

    for (int k = 0; k < local_sz; ++k)
    {
        int i = i_offset + k;
        int j_mastr = petsc_mastr_node_idxs[k];
        int j_slave = petsc_slave_node_idxs[k];
        ierr = MatSetValuesBlocked(D_mat,1,&i,1,&j_mastr,mastr_vals.data(),INSERT_VALUES);  IBTK_CHKERRQ(ierr);
        ierr = MatSetValuesBlocked(D_mat,1,&i,1,&j_slave,slave_vals.data(),INSERT_VALUES);  IBTK_CHKERRQ(ierr);
    }

    // Assemble the matrices.
    ierr = MatAssemblyBegin(D_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(D_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);

    // Indicate that the level data has been initialized.
    d_is_initialized[level_number] = true;

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
IBSpringForceGen::computeLagrangianForce(
    Pointer<LNodeLevelData> F_data,
    Pointer<LNodeLevelData> X_data,
    Pointer<LNodeLevelData> U_data,
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
    ierr = MatGetOwnershipRange(d_D_mats[level_number], &i_start, &i_stop);
    IBTK_CHKERRQ(ierr);

    Vec D_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop-i_start, PETSC_DECIDE, &D_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecSetBlockSize(D_vec, NDIM);                                          IBTK_CHKERRQ(ierr);

    // Compute the node displacements.
    ierr = MatMult(d_D_mats[level_number], X_data->getGlobalVec(), D_vec);
    IBTK_CHKERRQ(ierr);

    // Compute the spring forces acting on the nodes of the Lagrangian mesh.
    double* D_vals;
    ierr = VecGetArray(D_vec, &D_vals);  IBTK_CHKERRQ(ierr);

    std::vector<int>& lag_mastr_node_idxs = d_lag_mastr_node_idxs[level_number];
    std::vector<int>& lag_slave_node_idxs = d_lag_slave_node_idxs[level_number];
    std::vector<int>& petsc_mastr_node_idxs = d_petsc_mastr_node_idxs[level_number];
    std::vector<int>& petsc_slave_node_idxs = d_petsc_slave_node_idxs[level_number];
    std::vector<int>& force_fcn_idxs = d_force_fcn_idxs[level_number];
    std::vector<double>& stiffnesses = d_stiffnesses[level_number];
    std::vector<double>& rest_lengths = d_rest_lengths[level_number];

    const int local_sz = petsc_mastr_node_idxs.size();
    std::vector<double> F_mastr_node_vals(NDIM*local_sz,0.0);
    std::vector<double> F_slave_node_vals(NDIM*local_sz,0.0);

    for (int k = 0; k < local_sz; ++k)
    {
        // Compute the force applied by the spring to the "master" node.
        double* const F_mastr_node = &F_mastr_node_vals[k*NDIM];
        const double* const D = &D_vals[k*NDIM];
        const double& stf = stiffnesses[k];
        const double& rst = rest_lengths[k];
        const int& lag_mastr_idx = lag_mastr_node_idxs[k];
        const int& lag_slave_idx = lag_slave_node_idxs[k];
        const int& force_fcn_id = force_fcn_idxs[k];
        d_force_fcn_map[force_fcn_id](F_mastr_node,D,stf,rst,lag_mastr_idx,lag_slave_idx);

        // Compute the force applied by the spring to the corresponding "slave"
        // node.
        double* const F_slave_node = &F_slave_node_vals[k*NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            F_slave_node[d] = -F_mastr_node[d];
        }
    }

    ierr = VecRestoreArray(D_vec, &D_vals);  IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(D_vec);                IBTK_CHKERRQ(ierr);

    Vec F_vec = F_data->getGlobalVec();
#if 0
    ierr = VecSetValuesBlocked(F_vec,
                               petsc_mastr_node_idxs.size(),
                               !petsc_mastr_node_idxs.empty() ? &petsc_mastr_node_idxs[0] : PETSC_NULL,
                               !petsc_mastr_node_idxs.empty() ? &    F_mastr_node_vals[0] : PETSC_NULL,
                               ADD_VALUES);  IBTK_CHKERRQ(ierr);
    ierr = VecSetValuesBlocked(F_vec,
                               petsc_slave_node_idxs.size(),
                               !petsc_slave_node_idxs.empty() ? &petsc_slave_node_idxs[0] : PETSC_NULL,
                               !petsc_slave_node_idxs.empty() ? &    F_slave_node_vals[0] : PETSC_NULL,
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
                                            petsc_slave_node_idxs.size(),
                                            !petsc_slave_node_idxs.empty() ? &petsc_slave_node_idxs[0] : PETSC_NULL,
                                            !petsc_slave_node_idxs.empty() ? &    F_slave_node_vals[0] : PETSC_NULL,
                                            ADD_VALUES);  IBTK_CHKERRQ(ierr);
    ierr = PETScVecOps::VecAssemblyBegin(F_vec);  IBTK_CHKERRQ(ierr);
    ierr = PETScVecOps::VecAssemblyEnd(F_vec);    IBTK_CHKERRQ(ierr);
#endif
    t_compute_lagrangian_force->stop();
    return;
}// computeLagrangianForce

void
IBSpringForceGen::computeLagrangianForceJacobianNonzeroStructure(
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
    std::vector<int>& petsc_slave_node_idxs = d_petsc_slave_node_idxs[level_number];
    const int local_sz = petsc_mastr_node_idxs.size();

    // Determine the global node offset and the number of local nodes.
    const int global_node_offset = lag_manager->getGlobalNodeOffset(level_number);
    const int num_local_nodes = lag_manager->getNumberOfLocalNodes(level_number);

    // Determine the non-zero structure for the matrix used to store the
    // Jacobian of the force.
    //
    // NOTE: Each edge is *only* associated with a single node in the mesh.  We
    // must take this into account when determining the non-zero structure of
    // the matrix.
    Vec d_nnz_vec, o_nnz_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, num_local_nodes, PETSC_DETERMINE, &d_nnz_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, num_local_nodes, PETSC_DETERMINE, &o_nnz_vec);  IBTK_CHKERRQ(ierr);

    ierr = VecSet(d_nnz_vec, 1.0);  IBTK_CHKERRQ(ierr);
    ierr = VecSet(o_nnz_vec, 0.0);  IBTK_CHKERRQ(ierr);

    for (int k = 0; k < local_sz; ++k)
    {
        const int& mastr_idx = petsc_mastr_node_idxs[k];
        const int& slave_idx = petsc_slave_node_idxs[k];

        const bool slave_is_local = (slave_idx >= global_node_offset &&
                                     slave_idx <  global_node_offset + num_local_nodes);

        static const int N = 2;
        const int idxs[N] = { mastr_idx , slave_idx };
        const double vals[N] = { 1.0 , 1.0 };
        if (slave_is_local)
        {
            ierr = VecSetValues(d_nnz_vec, N, idxs, vals, ADD_VALUES);  IBTK_CHKERRQ(ierr);
        }
        else
        {
            ierr = VecSetValues(o_nnz_vec, N, idxs, vals, ADD_VALUES);  IBTK_CHKERRQ(ierr);
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
IBSpringForceGen::computeLagrangianForceJacobian(
    Mat& J_mat,
    MatAssemblyType assembly_type,
    const double X_coef,
    Pointer<LNodeLevelData> X_data,
    const double U_coef,
    Pointer<LNodeLevelData> U_data,
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

    // Create an appropriately sized temporary vector to store the node
    // displacements.
    int i_start, i_stop;
    ierr = MatGetOwnershipRange(d_D_mats[level_number], &i_start, &i_stop);
    IBTK_CHKERRQ(ierr);

    Vec D_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop-i_start, PETSC_DECIDE, &D_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecSetBlockSize(D_vec, NDIM);                                          IBTK_CHKERRQ(ierr);

    // Compute the node displacements.
    ierr = MatMult(d_D_mats[level_number], X_data->getGlobalVec(), D_vec);
    IBTK_CHKERRQ(ierr);

    // Compute the force Jacobians and insert them into the Jacobian matrix.
    double* D_vals;
    ierr = VecGetArray(D_vec, &D_vals);  IBTK_CHKERRQ(ierr);

    std::vector<int>& lag_mastr_node_idxs = d_lag_mastr_node_idxs[level_number];
    std::vector<int>& lag_slave_node_idxs = d_lag_slave_node_idxs[level_number];
    std::vector<int>& petsc_mastr_node_idxs = d_petsc_mastr_node_idxs[level_number];
    std::vector<int>& petsc_slave_node_idxs = d_petsc_slave_node_idxs[level_number];
    std::vector<int>& force_fcn_idxs = d_force_fcn_idxs[level_number];
    std::vector<double>& stiffnesses = d_stiffnesses[level_number];
    std::vector<double>& rest_lengths = d_rest_lengths[level_number];

    const int local_sz = petsc_mastr_node_idxs.size();
    for (int k = 0; k < local_sz; ++k)
    {
        // Compute the Jacobian of the force applied by the spring to the
        // "master" node with respect to the position of the "slave" node.
        blitz::Array<double,2> dF_dX(NDIM,NDIM);
        const double* const D = &D_vals[k*NDIM];
        const double& stf = stiffnesses[k];
        const double& rst = rest_lengths[k];
        const int& lag_mastr_idx = lag_mastr_node_idxs[k];
        const int& lag_slave_idx = lag_slave_node_idxs[k];
        const int& force_fcn_id = force_fcn_idxs[k];

        static const double eps = sqrt(std::numeric_limits<double>::epsilon());
        for (int j = 0; j < NDIM; ++j)
        {
            blitz::TinyVector<double,NDIM> D_eps(D);
            D_eps[j] += 1.0*eps;
            blitz::TinyVector<double,NDIM> F1;
            d_force_fcn_map[force_fcn_id](&F1[0],&D_eps[0],stf,rst,lag_mastr_idx,lag_slave_idx);
            D_eps[j] -= 2.0*eps;
            blitz::TinyVector<double,NDIM> F0;
            d_force_fcn_map[force_fcn_id](&F0[0],&D_eps[0],stf,rst,lag_mastr_idx,lag_slave_idx);
            for (int i = 0; i < NDIM; ++i)
            {
                dF_dX(i,j) = (F1[i]-F0[i])/(2.0*eps);
            }
        }

        // Scale the Jacobian entries appropriately.
        dF_dX *= X_coef;

        // Get the PETSc indices corresponding to the "master" and "slave"
        // nodes.
        const int& petsc_mastr_idx = petsc_mastr_node_idxs[k];
        const int& petsc_slave_idx = petsc_slave_node_idxs[k];

        // Accumulate the off-diagonal parts of the matrix.
        ierr = MatSetValuesBlocked(J_mat,1,&petsc_mastr_idx,1,&petsc_slave_idx,dF_dX.data(),ADD_VALUES);  IBTK_CHKERRQ(ierr);
        ierr = MatSetValuesBlocked(J_mat,1,&petsc_slave_idx,1,&petsc_mastr_idx,dF_dX.data(),ADD_VALUES);  IBTK_CHKERRQ(ierr);

        // Negate dF_dX to obtain the Jacobian of the force applied by the
        // spring to the "master" node with respect to the position of the
        // "master" node.
        dF_dX *= -1.0;

        // Accumulate the diagonal parts of the matrix.
        ierr = MatSetValuesBlocked(J_mat,1,&petsc_mastr_idx,1,&petsc_mastr_idx,dF_dX.data(),ADD_VALUES);  IBTK_CHKERRQ(ierr);
        ierr = MatSetValuesBlocked(J_mat,1,&petsc_slave_idx,1,&petsc_slave_idx,dF_dX.data(),ADD_VALUES);  IBTK_CHKERRQ(ierr);
    }

    ierr = VecRestoreArray(D_vec, &D_vals);  IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(D_vec);                IBTK_CHKERRQ(ierr);

    // Assemble the matrix.
    ierr = MatAssemblyBegin(J_mat, assembly_type);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(J_mat, assembly_type);    IBTK_CHKERRQ(ierr);

    t_compute_lagrangian_force_jacobian->stop();
    return;
}// computeLagrangianForceJacobian

double
IBSpringForceGen::computeLagrangianEnergy(
    Pointer<LNodeLevelData> X_data,
    Pointer<LNodeLevelData> U_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    TBOX_WARNING("IBSpringForceGen::computeLagrangianEnergy():\n"
                 << "  the implementation of this function is specialized to linear springs." << std::endl);

    if (!lag_manager->levelContainsLagrangianData(level_number)) return 0.0;

    t_compute_lagrangian_energy->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(level_number < int(d_is_initialized.size()));
    TBOX_ASSERT(d_is_initialized[level_number]);
#endif

    int ierr;

    // Create an appropriately sized temporary vector to store the node
    // displacements.
    int i_start, i_stop;
    ierr = MatGetOwnershipRange(d_D_mats[level_number], &i_start, &i_stop);
    IBTK_CHKERRQ(ierr);

    Vec D_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop-i_start, PETSC_DECIDE, &D_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecSetBlockSize(D_vec, NDIM);                                          IBTK_CHKERRQ(ierr);

    // Compute the node displacements.
    ierr = MatMult(d_D_mats[level_number], X_data->getGlobalVec(), D_vec);
    IBTK_CHKERRQ(ierr);

    // Compute the energy stored in a collection of linear springs with the
    // specified resting lengths.
    //
    // WARNING: This will not yield the correct result except for the standard
    // linear force function.
    double* D_vals;
    ierr = VecGetArray(D_vec, &D_vals);  IBTK_CHKERRQ(ierr);

    std::vector<double>& stiffnesses = d_stiffnesses[level_number];
    std::vector<double>& rest_lengths = d_rest_lengths[level_number];
    const int local_sz = stiffnesses.size();

    double energy = 0.0;
    for (int k = 0; k < local_sz; ++k)
    {
        const double* const D = &D_vals[k*NDIM];
        double R_sq = 0.0;
        for (int d = 0; d < NDIM; ++d)
        {
            R_sq += D[d]*D[d];
        }
        const double R = sqrt(R_sq);

        const double& kappa = stiffnesses[k];
        const double& R0 = rest_lengths[k];

        energy += 0.5*kappa*(R-R0)*(R-R0);
    }

    ierr = VecRestoreArray(D_vec, &D_vals);  IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(D_vec);                IBTK_CHKERRQ(ierr);

    t_compute_lagrangian_energy->stop();
    return SAMRAI_MPI::sumReduction(energy);
}// computeLagrangianEnergy

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBSpringForceGen::getFromInput(
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
template class Pointer<IBAMR::IBSpringForceGen>;

//////////////////////////////////////////////////////////////////////////////
