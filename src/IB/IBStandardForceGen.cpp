// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/IBBeamForceSpec.h"
#include "ibamr/IBSpringForceFunctions.h"
#include "ibamr/IBSpringForceSpec.h"
#include "ibamr/IBStandardForceGen.h"
#include "ibamr/IBTargetPointForceSpec.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep

#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LMesh.h"
#include "ibtk/LNode.h"
#include "ibtk/compiler_hints.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/private/LData-inl.h"
#include "ibtk/private/LDataManager-inl.h"
#include "ibtk/private/LMesh-inl.h"
#include "ibtk/private/LNode-inl.h"
#include "ibtk/private/LNodeIndex-inl.h"

#include "IntVector.h"
#include "PatchHierarchy.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include "petscmat.h"
#include "petscvec.h"
#include <petsclog.h>

IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/multi_array.hpp>
IBTK_ENABLE_EXTRA_WARNINGS

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
void
resetLocalPETScIndices(std::vector<int>& inds, const int global_node_offset, const int num_local_nodes)
{
#if defined(NDEBUG)
    NULL_USE(num_local_nodes);
#endif
    for (int& idx : inds)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(idx >= global_node_offset && idx < global_node_offset + num_local_nodes);
#endif
        idx -= global_node_offset;
    }
    return;
} // resetLocalPETScIndices

void
resetLocalOrNonlocalPETScIndices(std::vector<int>& inds,
                                 const int global_node_offset,
                                 const int num_local_nodes,
                                 const std::vector<int>& nonlocal_petsc_idxs)
{
    for (int& idx : inds)
    {
        if (idx >= global_node_offset && idx < global_node_offset + num_local_nodes)
        {
            // A local node.
            idx -= global_node_offset;
        }
        else
        {
            // A nonlocal node.
            //
            // First, lookup the slave node index in the set of ghost nodes.
            const std::vector<int>::const_iterator posn =
                std::lower_bound(nonlocal_petsc_idxs.begin(), nonlocal_petsc_idxs.end(), idx);
#if !defined(NDEBUG)
            TBOX_ASSERT(idx == *posn);
#endif
            // Second, set the local index via the offset of the ghost node
            // index within the set of ghost nodes.
            const int offset = static_cast<int>(std::distance(nonlocal_petsc_idxs.begin(), posn));
            idx = num_local_nodes + offset;
        }
    }
    return;
} // resetLocalOrNonlocalPETScIndices
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBStandardForceGen::IBStandardForceGen(Pointer<Database> input_db)
{
    // Setup the default force generation functions.
    registerSpringForceFunction(0, &default_spring_force, &default_spring_force_deriv);

    // Set up force generator from input.
    if (input_db)
    {
        if (input_db->keyExists("log_target_point_displacements"))
            d_log_target_point_displacements = input_db->getBool("log_target_point_displacements");
    }
    return;
} // IBStandardForceGen

void
IBStandardForceGen::registerSpringForceFunction(const int force_fcn_index,
                                                const SpringForceFcnPtr spring_force_fcn_ptr,
                                                const SpringForceDerivFcnPtr spring_force_deriv_fcn_ptr)
{
    d_spring_force_fcn_map[force_fcn_index] = spring_force_fcn_ptr;
    d_spring_force_deriv_fcn_map[force_fcn_index] = spring_force_deriv_fcn_ptr;
    return;
} // registerSpringForceFunction

void
IBStandardForceGen::initializeLevelData(const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                        const int level_number,
                                        const double init_data_time,
                                        const bool initial_time,
                                        LDataManager* const l_data_manager)
{
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return;

#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
#endif
    // Resize the vectors corresponding to data individually maintained for
    // separate levels of the patch hierarchy.
    const int new_size = std::max(level_number + 1, static_cast<int>(d_is_initialized.size()));

    d_spring_data.resize(new_size);
    d_beam_data.resize(new_size);
    d_target_point_data.resize(new_size);
    d_X_ghost_data.resize(new_size);
    d_F_ghost_data.resize(new_size);
    d_dX_data.resize(new_size);
    d_is_initialized.resize(new_size, false);

    // Keep track of all of the nonlocal PETSc indices required to compute the
    // forces.
    std::set<int> nonlocal_petsc_idx_set;

    // Setup the cached data.
    initializeSpringLevelData(
        nonlocal_petsc_idx_set, hierarchy, level_number, init_data_time, initial_time, l_data_manager);
    initializeBeamLevelData(
        nonlocal_petsc_idx_set, hierarchy, level_number, init_data_time, initial_time, l_data_manager);
    initializeTargetPointLevelData(
        nonlocal_petsc_idx_set, hierarchy, level_number, init_data_time, initial_time, l_data_manager);

    // Put the nonlocal PETSc indices into a vector.
    std::vector<int> nonlocal_petsc_idxs(nonlocal_petsc_idx_set.begin(), nonlocal_petsc_idx_set.end());

    // Put all cached PETSc node indices into local form.
    const int global_node_offset = l_data_manager->getGlobalNodeOffset(level_number);
    const int num_local_nodes = l_data_manager->getNumberOfLocalNodes(level_number);

    resetLocalPETScIndices(d_spring_data[level_number].petsc_mastr_node_idxs, global_node_offset, num_local_nodes);
    resetLocalPETScIndices(d_beam_data[level_number].petsc_mastr_node_idxs, global_node_offset, num_local_nodes);
    resetLocalPETScIndices(d_target_point_data[level_number].petsc_node_idxs, global_node_offset, num_local_nodes);

    resetLocalOrNonlocalPETScIndices(
        d_spring_data[level_number].petsc_slave_node_idxs, global_node_offset, num_local_nodes, nonlocal_petsc_idxs);
    resetLocalOrNonlocalPETScIndices(
        d_beam_data[level_number].petsc_next_node_idxs, global_node_offset, num_local_nodes, nonlocal_petsc_idxs);
    resetLocalOrNonlocalPETScIndices(
        d_beam_data[level_number].petsc_prev_node_idxs, global_node_offset, num_local_nodes, nonlocal_petsc_idxs);

    const std::string level_number_str = std::to_string(level_number);

    d_X_ghost_data[level_number] =
        new LData("IBStandardForceGen::X_ghost_" + level_number_str, num_local_nodes, NDIM, nonlocal_petsc_idxs);

    d_F_ghost_data[level_number] =
        new LData("IBStandardForceGen::F_ghost_" + level_number_str, num_local_nodes, NDIM, nonlocal_petsc_idxs);

    d_dX_data[level_number] = new LData("IBStandardForceGen::dX_" + level_number_str, num_local_nodes, NDIM);

    // Compute periodic displacements.
    boost::multi_array_ref<double, 2>& dX_array = *d_dX_data[level_number]->getLocalFormVecArray();
    const Pointer<LMesh> mesh = l_data_manager->getLMesh(level_number);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
    for (int k = 0; k < num_local_nodes; ++k)
    {
        const LNode* const node = local_nodes[k];
        const Vector& periodic_displacement = node->getPeriodicDisplacement();
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dX_array[k][d] = periodic_displacement[d];
        }
    }
    d_dX_data[level_number]->restoreArrays();

    // Transform all of the cached indices to correspond to a data depth of
    // NDIM.
    std::vector<std::vector<int>*> idx_vecs = { &d_spring_data[level_number].petsc_mastr_node_idxs,
                                                &d_spring_data[level_number].petsc_slave_node_idxs,
                                                &d_spring_data[level_number].petsc_global_mastr_node_idxs,
                                                &d_spring_data[level_number].petsc_global_slave_node_idxs,
                                                &d_beam_data[level_number].petsc_mastr_node_idxs,
                                                &d_beam_data[level_number].petsc_next_node_idxs,
                                                &d_beam_data[level_number].petsc_prev_node_idxs,
                                                &d_beam_data[level_number].petsc_global_mastr_node_idxs,
                                                &d_beam_data[level_number].petsc_global_next_node_idxs,
                                                &d_beam_data[level_number].petsc_global_prev_node_idxs,
                                                &d_target_point_data[level_number].petsc_node_idxs,
                                                &d_target_point_data[level_number].petsc_global_node_idxs };
    for (auto v_ptr : idx_vecs)
    {
        std::for_each(v_ptr->begin(), v_ptr->end(), [](int& i) { i *= NDIM; });
    }

    // Indicate that the level data has been initialized.
    d_is_initialized[level_number] = true;
    return;
} // initializeLevelData

void
IBStandardForceGen::computeLagrangianForce(Pointer<LData> F_data,
                                           Pointer<LData> X_data,
                                           Pointer<LData> U_data,
                                           const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                           const int level_number,
                                           const double data_time,
                                           LDataManager* const l_data_manager)
{
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return;

    int ierr;

    // Initialize ghost data.
    Pointer<LData> F_ghost_data = d_F_ghost_data[level_number];
    Vec F_ghost_local_form_vec;
    ierr = VecGhostGetLocalForm(F_ghost_data->getVec(), &F_ghost_local_form_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecSet(F_ghost_local_form_vec, 0.0);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostRestoreLocalForm(F_ghost_data->getVec(), &F_ghost_local_form_vec);
    IBTK_CHKERRQ(ierr);

    Pointer<LData> X_ghost_data = d_X_ghost_data[level_number];
    Pointer<LData> dX_data = d_dX_data[level_number];
    ierr = VecAXPBYPCZ(X_ghost_data->getVec(), 1.0, 1.0, 0.0, X_data->getVec(), dX_data->getVec());
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(X_ghost_data->getVec(), INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(X_ghost_data->getVec(), INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);

    // Compute the forces.
    computeLagrangianSpringForce(F_ghost_data, X_ghost_data, hierarchy, level_number, data_time, l_data_manager);
    computeLagrangianBeamForce(F_ghost_data, X_ghost_data, hierarchy, level_number, data_time, l_data_manager);
    computeLagrangianTargetPointForce(
        F_ghost_data, X_ghost_data, U_data, hierarchy, level_number, data_time, l_data_manager);

    // Add the locally computed forces to the Lagrangian force vector.
    //
    // WARNING: The following operations may yield nondeterministic results in
    // parallel environments (i.e., the order of summation may not be
    // consistent).
    ierr = VecGhostUpdateBegin(F_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);
    IBTK_CHKERRQ(ierr);
    ierr = VecAXPY(F_data->getVec(), 1.0, F_ghost_data->getVec());
    return;
} // computeLagrangianForce

void
IBStandardForceGen::computeLagrangianForceJacobianNonzeroStructure(std::vector<int>& d_nnz,
                                                                   std::vector<int>& o_nnz,
                                                                   const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                                                   const int level_number,
                                                                   LDataManager* const l_data_manager)
{
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return;

#if !defined(NDEBUG)
    TBOX_ASSERT(level_number < static_cast<int>(d_is_initialized.size()));
    TBOX_ASSERT(d_is_initialized[level_number]);
#endif

    int ierr;

    // Determine the global node offset and the number of local nodes.
    const int global_node_offset = l_data_manager->getGlobalNodeOffset(level_number);
    const int num_nodes = l_data_manager->getNumberOfNodes(level_number);
    const int num_local_nodes = l_data_manager->getNumberOfLocalNodes(level_number);

    // Determine the non-zero structure for the matrix used to store the
    // Jacobian of the force. Here, we are filling enteries for blocked rows.
    //
    // NOTE #1: Each spring and beam is *only* associated with a single node in
    // the mesh.  We must take this into account when determining the non-zero
    // structure of the matrix.
    //
    // NOTE #2: The following ensures only that sufficient space is allocated to
    // store the Jacobian matrix.  In general, this routine will request MORE
    // space than is ACTUALLY required.
    Vec d_nnz_vec, o_nnz_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, num_local_nodes, PETSC_DETERMINE, &d_nnz_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, num_local_nodes, PETSC_DETERMINE, &o_nnz_vec);
    IBTK_CHKERRQ(ierr);

    ierr = VecSet(d_nnz_vec, 1.0);
    IBTK_CHKERRQ(ierr);
    ierr = VecSet(o_nnz_vec, 0.0);
    IBTK_CHKERRQ(ierr);

    { // Spring forces.

        const std::vector<int>& petsc_global_mastr_node_idxs = d_spring_data[level_number].petsc_global_mastr_node_idxs;
        const std::vector<int>& petsc_global_slave_node_idxs = d_spring_data[level_number].petsc_global_slave_node_idxs;
        for (unsigned int k = 0; k < petsc_global_mastr_node_idxs.size(); ++k)
        {
            const int& mastr_idx = petsc_global_mastr_node_idxs[k];
            const int& slave_idx = petsc_global_slave_node_idxs[k];

            const bool slave_is_local =
                (slave_idx >= NDIM * global_node_offset && slave_idx < NDIM * (global_node_offset + num_local_nodes));

            static const int N = 2;
            const int idxs[N] = { mastr_idx / NDIM, slave_idx / NDIM };
            const double vals[N] = { 1.0, 1.0 };

            if (slave_is_local)
            {
                ierr = VecSetValues(d_nnz_vec, N, idxs, vals, ADD_VALUES);
                IBTK_CHKERRQ(ierr);
            }
            else
            {
                ierr = VecSetValues(o_nnz_vec, N, idxs, vals, ADD_VALUES);
                IBTK_CHKERRQ(ierr);
            }
        }
    }

    { // Beam forces.

        const std::vector<int>& petsc_global_mastr_node_idxs = d_beam_data[level_number].petsc_global_mastr_node_idxs;
        const std::vector<int>& petsc_global_next_node_idxs = d_beam_data[level_number].petsc_global_next_node_idxs;
        const std::vector<int>& petsc_global_prev_node_idxs = d_beam_data[level_number].petsc_global_prev_node_idxs;
        for (unsigned int k = 0; k < petsc_global_mastr_node_idxs.size(); ++k)
        {
            const int& mastr_idx = petsc_global_mastr_node_idxs[k];
            const int& next_idx = petsc_global_next_node_idxs[k];
            const int& prev_idx = petsc_global_prev_node_idxs[k];

            const bool next_is_local =
                (next_idx >= NDIM * global_node_offset && next_idx < NDIM * (global_node_offset + num_local_nodes));
            const bool prev_is_local =
                (prev_idx >= NDIM * global_node_offset && prev_idx < NDIM * (global_node_offset + num_local_nodes));

            if (next_is_local && prev_is_local)
            {
                static const int d_N = 3;
                const int d_idxs[d_N] = { mastr_idx / NDIM, next_idx / NDIM, prev_idx / NDIM };
                const double d_vals[d_N] = { 2.0, 2.0, 2.0 };
                ierr = VecSetValues(d_nnz_vec, d_N, d_idxs, d_vals, ADD_VALUES);
                IBTK_CHKERRQ(ierr);
            }
            else if (next_is_local && (!prev_is_local))
            {
                static const int d_N = 2;
                const int d_idxs[d_N] = { mastr_idx / NDIM, next_idx / NDIM };
                const double d_vals[d_N] = { 1.0, 1.0 };
                ierr = VecSetValues(d_nnz_vec, d_N, d_idxs, d_vals, ADD_VALUES);
                IBTK_CHKERRQ(ierr);

                static const int o_N = 3;
                const int o_idxs[o_N] = { mastr_idx / NDIM, next_idx / NDIM, prev_idx / NDIM };
                const double o_vals[o_N] = { 1.0, 1.0, 2.0 };
                ierr = VecSetValues(o_nnz_vec, o_N, o_idxs, o_vals, ADD_VALUES);
                IBTK_CHKERRQ(ierr);
            }
            else if ((!next_is_local) && prev_is_local)
            {
                static const int d_N = 2;
                const int d_idxs[d_N] = { mastr_idx / NDIM, prev_idx / NDIM };
                const double d_vals[d_N] = { 1.0, 1.0 };
                ierr = VecSetValues(d_nnz_vec, d_N, d_idxs, d_vals, ADD_VALUES);
                IBTK_CHKERRQ(ierr);

                static const int o_N = 3;
                const int o_idxs[o_N] = { mastr_idx / NDIM, next_idx / NDIM, prev_idx / NDIM };
                const double o_vals[o_N] = { 1.0, 2.0, 1.0 };
                ierr = VecSetValues(o_nnz_vec, o_N, o_idxs, o_vals, ADD_VALUES);
                IBTK_CHKERRQ(ierr);
            }
            else
            {
                // NOTE: Rather than trying to find out if the previous and next
                // nodes are assigned to the same processor, we instead allocate
                // space both for the case that the previous and next nodes are
                // on different processors, and for the case that the previous
                // and next nodes are on the same processor.
                static const int d_N = 2;
                const int d_idxs[d_N] = { next_idx / NDIM, prev_idx / NDIM };
                const double d_vals[d_N] = { 2.0, 2.0 };
                ierr = VecSetValues(d_nnz_vec, d_N, d_idxs, d_vals, ADD_VALUES);
                IBTK_CHKERRQ(ierr);

                static const int o_N = 3;
                const int o_idxs[o_N] = { mastr_idx / NDIM, next_idx / NDIM, prev_idx / NDIM };
                const double o_vals[o_N] = { 2.0, 2.0, 2.0 };
                ierr = VecSetValues(o_nnz_vec, o_N, o_idxs, o_vals, ADD_VALUES);
                IBTK_CHKERRQ(ierr);
            }
        }
    }

    ierr = VecAssemblyBegin(d_nnz_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyBegin(o_nnz_vec);
    IBTK_CHKERRQ(ierr);

    ierr = VecAssemblyEnd(d_nnz_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(o_nnz_vec);
    IBTK_CHKERRQ(ierr);

    double* d_nnz_vec_arr;
    ierr = VecGetArray(d_nnz_vec, &d_nnz_vec_arr);
    IBTK_CHKERRQ(ierr);

    double* o_nnz_vec_arr;
    ierr = VecGetArray(o_nnz_vec, &o_nnz_vec_arr);
    IBTK_CHKERRQ(ierr);

    d_nnz.resize(num_local_nodes);
    o_nnz.resize(num_local_nodes);
    for (int k = 0; k < num_local_nodes; ++k)
    {
        d_nnz[k] = std::min(static_cast<int>(d_nnz_vec_arr[k]), num_local_nodes);
        o_nnz[k] = std::min(static_cast<int>(o_nnz_vec_arr[k]), num_nodes - num_local_nodes);
    }

    ierr = VecRestoreArray(d_nnz_vec, &d_nnz_vec_arr);
    IBTK_CHKERRQ(ierr);
    ierr = VecRestoreArray(o_nnz_vec, &o_nnz_vec_arr);
    IBTK_CHKERRQ(ierr);

    ierr = VecDestroy(&d_nnz_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&o_nnz_vec);
    IBTK_CHKERRQ(ierr);
    return;
} // computeLagrangianForceJacobianNonzeroStructure

void
IBStandardForceGen::computeLagrangianForceJacobian(Mat& J_mat,
                                                   MatAssemblyType assembly_type,
                                                   const double X_coef,
                                                   Pointer<LData> X_data,
                                                   const double U_coef,
                                                   Pointer<LData> /*U_data*/,
                                                   const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                                   const int level_number,
                                                   const double /*data_time*/,
                                                   LDataManager* const l_data_manager)
{
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return;

#if !defined(NDEBUG)
    TBOX_ASSERT(level_number < static_cast<int>(d_is_initialized.size()));
    TBOX_ASSERT(d_is_initialized[level_number]);
#endif

    int ierr;
    Pointer<LData> X_ghost_data = d_X_ghost_data[level_number];
    Pointer<LData> dX_data = d_dX_data[level_number];
    ierr = VecAXPBYPCZ(X_ghost_data->getVec(), 1.0, 1.0, 0.0, X_data->getVec(), dX_data->getVec());
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(X_ghost_data->getVec(), INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(X_ghost_data->getVec(), INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);

    { // Spring forces.
        const std::vector<int>& lag_mastr_node_idxs = d_spring_data[level_number].lag_mastr_node_idxs;
        const std::vector<int>& lag_slave_node_idxs = d_spring_data[level_number].lag_slave_node_idxs;
        const std::vector<int>& petsc_mastr_node_idxs = d_spring_data[level_number].petsc_mastr_node_idxs;
        const std::vector<int>& petsc_slave_node_idxs = d_spring_data[level_number].petsc_slave_node_idxs;
        const std::vector<int>& petsc_global_mastr_node_idxs = d_spring_data[level_number].petsc_global_mastr_node_idxs;
        const std::vector<int>& petsc_global_slave_node_idxs = d_spring_data[level_number].petsc_global_slave_node_idxs;
        const std::vector<SpringForceFcnPtr>& force_fcns = d_spring_data[level_number].force_fcns;
        const std::vector<SpringForceDerivFcnPtr>& force_deriv_fcns = d_spring_data[level_number].force_deriv_fcns;
        const std::vector<const double*>& parameters = d_spring_data[level_number].parameters;
        const double* const X_node = X_ghost_data->getGhostedLocalFormVecArray()->data();
        MatrixNd dF_dX;
        Vector D;
        double dT_dR;
        for (unsigned int k = 0; k < petsc_mastr_node_idxs.size(); ++k)
        {
            // Compute the Jacobian of the force applied by the spring to the
            // "master" node with respect to the position of the "slave" node.
            const int& lag_mastr_idx = lag_mastr_node_idxs[k];
            const int& lag_slave_idx = lag_slave_node_idxs[k];
            int petsc_mastr_idx = petsc_mastr_node_idxs[k];
            int petsc_slave_idx = petsc_slave_node_idxs[k];
            int petsc_global_mastr_idx = petsc_global_mastr_node_idxs[k];
            int petsc_global_slave_idx = petsc_global_slave_node_idxs[k];
            const SpringForceFcnPtr force_fcn = force_fcns[k];
            const SpringForceDerivFcnPtr force_deriv_fcn = force_deriv_fcns[k];
            const double* const params = parameters[k];
            for (unsigned int i = 0; i < NDIM; ++i)
            {
                D(i) = X_node[petsc_slave_idx + i] - X_node[petsc_mastr_idx + i];
            }
            const double R = D.norm();
            const double T = force_fcn(R, params, lag_mastr_idx, lag_slave_idx);
            if (!force_deriv_fcn)
            {
                // Use finite differences to approximate dT/dR.
                const double eps = std::max(R, 1.0) * std::pow(std::numeric_limits<double>::epsilon(), 1.0 / 3.0);
                dT_dR = (force_fcn(R + eps, params, lag_mastr_idx, lag_slave_idx) -
                         force_fcn(R - eps, params, lag_mastr_idx, lag_slave_idx)) /
                        (2.0 * eps);
            }
            else
            {
                dT_dR = force_deriv_fcn(R, params, lag_mastr_idx, lag_slave_idx);
            }

            // F = T(R) D/R
            //
            // dF_k/dx_l = (dT/dR * dR/dx_l * D/R) + (T/R * dD/dx_l) + (T *D* (-1/R^2)* dR/dx_l)
            //
            // dR/dx_l = 1/R * D if x_l is "slave" and dR/dx_l = -1/R * D if x_l is "master"
            //
            // dD/dx_l = e_l outer e_k dD_k/dx_l = e_l outer e_l = I if x_l is "slave"
            // and e_l outer -e_l = -I if x_l is master.
            //
            // => dF_k/dx_l = (1/R^2 * dT/dR * D outer D) + (T/R * I) - (1/R^2 * T/R * D outer D) if x_l is "slave"
            // and -dF_k/dx_l if x_l is "master"

            for (unsigned int i = 0; i < NDIM; ++i)
            {
                for (unsigned int j = 0; j < NDIM; ++j)
                {
                    dF_dX(i, j) = X_coef * ((T / R) * ((i == j ? 1.0 : 0.0)) + (dT_dR - T / R) * D[i] * D[j] / (R * R));
                }
            }

            // Rows and cols for blocked matrix.
            petsc_global_mastr_idx /= NDIM;
            petsc_global_slave_idx /= NDIM;

            // Accumulate the off-diagonal parts of the matrix.
            ierr = MatSetValuesBlocked(
                J_mat, 1, &petsc_global_mastr_idx, 1, &petsc_global_slave_idx, dF_dX.data(), ADD_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetValuesBlocked(
                J_mat, 1, &petsc_global_slave_idx, 1, &petsc_global_mastr_idx, dF_dX.data(), ADD_VALUES);
            IBTK_CHKERRQ(ierr);

            // Negate dF_dX to obtain the Jacobian of the force applied by the
            // spring to the "master" node with respect to the position of the
            // "master" node.
            dF_dX *= -1.0;

            // Accumulate the diagonal parts of the matrix.
            ierr = MatSetValuesBlocked(
                J_mat, 1, &petsc_global_mastr_idx, 1, &petsc_global_mastr_idx, dF_dX.data(), ADD_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetValuesBlocked(
                J_mat, 1, &petsc_global_slave_idx, 1, &petsc_global_slave_idx, dF_dX.data(), ADD_VALUES);
            IBTK_CHKERRQ(ierr);
        }
    }

    { // Beam forces.
        const std::vector<int>& petsc_global_mastr_node_idxs = d_beam_data[level_number].petsc_global_mastr_node_idxs;
        const std::vector<int>& petsc_global_next_node_idxs = d_beam_data[level_number].petsc_global_next_node_idxs;
        const std::vector<int>& petsc_global_prev_node_idxs = d_beam_data[level_number].petsc_global_prev_node_idxs;
        const std::vector<const double*>& rigidities = d_beam_data[level_number].rigidities;
        MatrixNd dF_dX(MatrixNd::Zero());
        for (unsigned int k = 0; k < petsc_global_mastr_node_idxs.size(); ++k)
        {
            const int petsc_global_mastr_idx = petsc_global_mastr_node_idxs[k] / NDIM; // block indices
            const int petsc_global_next_idx = petsc_global_next_node_idxs[k] / NDIM;
            const int petsc_global_prev_idx = petsc_global_prev_node_idxs[k] / NDIM;
            const double& bend = *rigidities[k];

            for (unsigned int alpha = 0; alpha < NDIM; ++alpha)
            {
                dF_dX(alpha, alpha) = -1.0 * bend * X_coef;
            }
            ierr = MatSetValuesBlocked(
                J_mat, 1, &petsc_global_prev_idx, 1, &petsc_global_prev_idx, dF_dX.data(), ADD_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetValuesBlocked(
                J_mat, 1, &petsc_global_prev_idx, 1, &petsc_global_next_idx, dF_dX.data(), ADD_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetValuesBlocked(
                J_mat, 1, &petsc_global_next_idx, 1, &petsc_global_prev_idx, dF_dX.data(), ADD_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetValuesBlocked(
                J_mat, 1, &petsc_global_next_idx, 1, &petsc_global_next_idx, dF_dX.data(), ADD_VALUES);
            IBTK_CHKERRQ(ierr);

            for (unsigned int alpha = 0; alpha < NDIM; ++alpha)
            {
                dF_dX(alpha, alpha) = +2.0 * bend * X_coef;
            }
            ierr = MatSetValuesBlocked(
                J_mat, 1, &petsc_global_prev_idx, 1, &petsc_global_mastr_idx, dF_dX.data(), ADD_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetValuesBlocked(
                J_mat, 1, &petsc_global_next_idx, 1, &petsc_global_mastr_idx, dF_dX.data(), ADD_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetValuesBlocked(
                J_mat, 1, &petsc_global_mastr_idx, 1, &petsc_global_prev_idx, dF_dX.data(), ADD_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetValuesBlocked(
                J_mat, 1, &petsc_global_mastr_idx, 1, &petsc_global_next_idx, dF_dX.data(), ADD_VALUES);
            IBTK_CHKERRQ(ierr);

            for (unsigned int alpha = 0; alpha < NDIM; ++alpha)
            {
                dF_dX(alpha, alpha) = -4.0 * bend * X_coef;
            }
            ierr = MatSetValuesBlocked(
                J_mat, 1, &petsc_global_mastr_idx, 1, &petsc_global_mastr_idx, dF_dX.data(), ADD_VALUES);
            IBTK_CHKERRQ(ierr);
        }
    }

    { // Target point forces.
        const std::vector<int>& petsc_global_node_idxs = d_target_point_data[level_number].petsc_global_node_idxs;
        const std::vector<const double*>& kappa = d_target_point_data[level_number].kappa;
        const std::vector<const double*>& eta = d_target_point_data[level_number].eta;
        MatrixNd dF_dX(MatrixNd::Zero());
        for (unsigned int k = 0; k < petsc_global_node_idxs.size(); ++k)
        {
            const int petsc_global_node_idx = petsc_global_node_idxs[k] / NDIM; // block index
            const double& K = *kappa[k];
            const double& E = *eta[k];
            for (unsigned int alpha = 0; alpha < NDIM; ++alpha)
            {
                dF_dX(alpha, alpha) = -X_coef * K - U_coef * E;
            }
            ierr = MatSetValuesBlocked(
                J_mat, 1, &petsc_global_node_idx, 1, &petsc_global_node_idx, dF_dX.data(), ADD_VALUES);
            IBTK_CHKERRQ(ierr);
        }
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(J_mat, assembly_type);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(J_mat, assembly_type);
    IBTK_CHKERRQ(ierr);
    return;
} // computeLagrangianForceJacobian

double
IBStandardForceGen::computeLagrangianEnergy(Pointer<LData> /*X_data*/,
                                            Pointer<LData> /*U_data*/,
                                            const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                            const int level_number,
                                            const double /*data_time*/,
                                            LDataManager* const l_data_manager)
{
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return 0.0;

    // Compute the energy.
    TBOX_ERROR("not currently implemented\n");
    return std::numeric_limits<double>::quiet_NaN();
} // computeLagrangianEnergy

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBStandardForceGen::initializeSpringLevelData(std::set<int>& nonlocal_petsc_idx_set,
                                              const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                              const int level_number,
                                              const double /*init_data_time*/,
                                              const bool /*initial_time*/,
                                              LDataManager* const l_data_manager)
{
    std::vector<int>& lag_mastr_node_idxs = d_spring_data[level_number].lag_mastr_node_idxs;
    std::vector<int>& lag_slave_node_idxs = d_spring_data[level_number].lag_slave_node_idxs;
    std::vector<int>& petsc_mastr_node_idxs = d_spring_data[level_number].petsc_mastr_node_idxs;
    std::vector<int>& petsc_slave_node_idxs = d_spring_data[level_number].petsc_slave_node_idxs;
    std::vector<int>& petsc_global_mastr_node_idxs = d_spring_data[level_number].petsc_global_mastr_node_idxs;
    std::vector<int>& petsc_global_slave_node_idxs = d_spring_data[level_number].petsc_global_slave_node_idxs;
    std::vector<SpringForceFcnPtr>& force_fcns = d_spring_data[level_number].force_fcns;
    std::vector<SpringForceDerivFcnPtr>& force_deriv_fcns = d_spring_data[level_number].force_deriv_fcns;
    std::vector<const double*>& parameters = d_spring_data[level_number].parameters;

    // The LMesh object provides the set of local Lagrangian nodes.
    const Pointer<LMesh> mesh = l_data_manager->getLMesh(level_number);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
    const int num_local_nodes = static_cast<int>(local_nodes.size());

    // Determine how many springs are associated with the present MPI process.
    unsigned int total_num_springs = 0;
    for (const auto& node_idx : local_nodes)
    {
        const IBSpringForceSpec* const force_spec = node_idx->getNodeDataItem<IBSpringForceSpec>();
        if (force_spec) total_num_springs += force_spec->getNumberOfSprings();
    }

    // Resize arrays for storing cached values used to compute spring forces.
    lag_mastr_node_idxs.resize(total_num_springs);
    lag_slave_node_idxs.resize(total_num_springs);
    petsc_mastr_node_idxs.resize(total_num_springs);
    petsc_slave_node_idxs.resize(total_num_springs);
    petsc_global_mastr_node_idxs.resize(total_num_springs);
    petsc_global_slave_node_idxs.resize(total_num_springs);
    force_fcns.resize(total_num_springs);
    force_deriv_fcns.resize(total_num_springs);
    parameters.resize(total_num_springs);

    // Setup the data structures used to compute spring forces.
    int current_spring = 0;
    for (const auto& node_idx : local_nodes)
    {
        const IBSpringForceSpec* const force_spec = node_idx->getNodeDataItem<IBSpringForceSpec>();
        if (!force_spec) continue;

        const int lag_idx = node_idx->getLagrangianIndex();
#if !defined(NDEBUG)
        TBOX_ASSERT(lag_idx == force_spec->getMasterNodeIndex());
#endif
        const int petsc_idx = node_idx->getGlobalPETScIndex();
        const std::vector<int>& slv = force_spec->getSlaveNodeIndices();
        const std::vector<int>& fcn = force_spec->getForceFunctionIndices();
        const std::vector<std::vector<double> >& params = force_spec->getParameters();
        const unsigned int num_springs = force_spec->getNumberOfSprings();
#if !defined(NDEBUG)
        TBOX_ASSERT(num_springs == slv.size());
        TBOX_ASSERT(num_springs == fcn.size());
        TBOX_ASSERT(num_springs == params.size());
#endif
        for (unsigned int k = 0; k < num_springs; ++k)
        {
            lag_mastr_node_idxs[current_spring] = lag_idx;
            lag_slave_node_idxs[current_spring] = slv[k];
            petsc_mastr_node_idxs[current_spring] = petsc_idx;
            force_fcns[current_spring] = d_spring_force_fcn_map[fcn[k]];
            force_deriv_fcns[current_spring] = d_spring_force_deriv_fcn_map[fcn[k]];
            parameters[current_spring] = params.empty() ? nullptr : &params[k][0];
            ++current_spring;
        }
    }

    // Map the Lagrangian slave node indices to the PETSc indices corresponding
    // to the present data distribution.
    petsc_slave_node_idxs = lag_slave_node_idxs;
    l_data_manager->mapLagrangianToPETSc(petsc_slave_node_idxs, level_number);

    // Keep a copy of global PETSc indices.
    petsc_global_mastr_node_idxs = petsc_mastr_node_idxs;
    petsc_global_slave_node_idxs = petsc_slave_node_idxs;

    // Determine the ghost nodes required to compute spring forces.
    //
    // NOTE: Only slave nodes can be "off processor".  Master nodes are
    // guaranteed to be "on processor".
    const int global_node_offset = l_data_manager->getGlobalNodeOffset(level_number);
    for (int idx : petsc_slave_node_idxs)
    {
        if (UNLIKELY(idx < global_node_offset || idx >= global_node_offset + num_local_nodes))
        {
            nonlocal_petsc_idx_set.insert(idx);
        }
    }
    return;
} // initializeSpringLevelData

void
IBStandardForceGen::computeLagrangianSpringForce(Pointer<LData> F_data,
                                                 Pointer<LData> X_data,
                                                 const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                                 const int level_number,
                                                 const double /*data_time*/,
                                                 LDataManager* const /*l_data_manager*/)
{
    const int num_springs = static_cast<int>(d_spring_data[level_number].lag_mastr_node_idxs.size());
    const bool uses_springs = (num_springs > 0);
    const int* const lag_mastr_node_idxs = uses_springs ? &d_spring_data[level_number].lag_mastr_node_idxs[0] : nullptr;
    const int* const lag_slave_node_idxs = uses_springs ? &d_spring_data[level_number].lag_slave_node_idxs[0] : nullptr;
    const int* const petsc_mastr_node_idxs =
        uses_springs ? &d_spring_data[level_number].petsc_mastr_node_idxs[0] : nullptr;
    const int* const petsc_slave_node_idxs =
        uses_springs ? &d_spring_data[level_number].petsc_slave_node_idxs[0] : nullptr;
    const SpringForceFcnPtr* const force_fcns = uses_springs ? &d_spring_data[level_number].force_fcns[0] : nullptr;
    const double** const parameters = uses_springs ? &d_spring_data[level_number].parameters[0] : nullptr;
    double* const F_node = F_data->getLocalFormVecArray()->data();
    const double* const X_node = X_data->getGhostedLocalFormVecArray()->data();

    static const int BLOCKSIZE = 16; // this parameter needs to be tuned
    int k, kblock, kunroll, mastr_idx, slave_idx;
    double F[NDIM], D[NDIM], R, T_over_R;
    kblock = 0;
    for (; kblock < (num_springs - 1) / BLOCKSIZE;
         ++kblock) // ensure that the last block is NOT handled by this first loop
    {
        PREFETCH_READ_NTA_BLOCK(lag_mastr_node_idxs + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(lag_slave_node_idxs + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(petsc_mastr_node_idxs + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(petsc_slave_node_idxs + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(force_fcns + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(parameters + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        for (kunroll = 0; kunroll < BLOCKSIZE; ++kunroll)
        {
            k = kblock * BLOCKSIZE + kunroll;
            mastr_idx = petsc_mastr_node_idxs[k];
            slave_idx = petsc_slave_node_idxs[k];
#if !defined(NDEBUG)
            TBOX_ASSERT(mastr_idx != slave_idx);
#endif
            PREFETCH_READ_NTA_NDIM_BLOCK(F_node + petsc_mastr_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(F_node + petsc_slave_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(X_node + petsc_mastr_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(X_node + petsc_slave_node_idxs[k + 1]);
            PREFETCH_READ_NTA(parameters[k + 1]);
            D[0] = X_node[slave_idx + 0] - X_node[mastr_idx + 0];
            D[1] = X_node[slave_idx + 1] - X_node[mastr_idx + 1];
#if (NDIM == 3)
            D[2] = X_node[slave_idx + 2] - X_node[mastr_idx + 2];
#endif
#if (NDIM == 2)
            R = std::sqrt(D[0] * D[0] + D[1] * D[1]);
#endif
#if (NDIM == 3)
            R = std::sqrt(D[0] * D[0] + D[1] * D[1] + D[2] * D[2]);
#endif
            if (UNLIKELY(R < std::numeric_limits<double>::epsilon())) continue;
            T_over_R = (force_fcns[k])(R, parameters[k], lag_mastr_node_idxs[k], lag_slave_node_idxs[k]) / R;
            F[0] = T_over_R * D[0];
            F[1] = T_over_R * D[1];
#if (NDIM == 3)
            F[2] = T_over_R * D[2];
#endif
            F_node[mastr_idx + 0] += F[0];
            F_node[mastr_idx + 1] += F[1];
#if (NDIM == 3)
            F_node[mastr_idx + 2] += F[2];
#endif
            F_node[slave_idx + 0] -= F[0];
            F_node[slave_idx + 1] -= F[1];
#if (NDIM == 3)
            F_node[slave_idx + 2] -= F[2];
#endif
        }
    }
    for (k = kblock * BLOCKSIZE; k < num_springs; ++k)
    {
        mastr_idx = petsc_mastr_node_idxs[k];
        slave_idx = petsc_slave_node_idxs[k];
#if !defined(NDEBUG)
        TBOX_ASSERT(mastr_idx != slave_idx);
#endif
        D[0] = X_node[slave_idx + 0] - X_node[mastr_idx + 0];
        D[1] = X_node[slave_idx + 1] - X_node[mastr_idx + 1];
#if (NDIM == 3)
        D[2] = X_node[slave_idx + 2] - X_node[mastr_idx + 2];
#endif
#if (NDIM == 2)
        R = std::sqrt(D[0] * D[0] + D[1] * D[1]);
#endif
#if (NDIM == 3)
        R = std::sqrt(D[0] * D[0] + D[1] * D[1] + D[2] * D[2]);
#endif
        if (UNLIKELY(R < std::numeric_limits<double>::epsilon())) continue;
        T_over_R = (force_fcns[k])(R, parameters[k], lag_mastr_node_idxs[k], lag_slave_node_idxs[k]) / R;
        F[0] = T_over_R * D[0];
        F[1] = T_over_R * D[1];
#if (NDIM == 3)
        F[2] = T_over_R * D[2];
#endif
        F_node[mastr_idx + 0] += F[0];
        F_node[mastr_idx + 1] += F[1];
#if (NDIM == 3)
        F_node[mastr_idx + 2] += F[2];
#endif
        F_node[slave_idx + 0] -= F[0];
        F_node[slave_idx + 1] -= F[1];
#if (NDIM == 3)
        F_node[slave_idx + 2] -= F[2];
#endif
    }

    F_data->restoreArrays();
    X_data->restoreArrays();
    return;
} // computeLagrangianSpringForce

void
IBStandardForceGen::initializeBeamLevelData(std::set<int>& nonlocal_petsc_idx_set,
                                            const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                            const int level_number,
                                            const double /*init_data_time*/,
                                            const bool /*initial_time*/,
                                            LDataManager* const l_data_manager)
{
    std::vector<int>& petsc_mastr_node_idxs = d_beam_data[level_number].petsc_mastr_node_idxs;
    std::vector<int>& petsc_next_node_idxs = d_beam_data[level_number].petsc_next_node_idxs;
    std::vector<int>& petsc_prev_node_idxs = d_beam_data[level_number].petsc_prev_node_idxs;
    std::vector<int>& petsc_global_mastr_node_idxs = d_beam_data[level_number].petsc_global_mastr_node_idxs;
    std::vector<int>& petsc_global_next_node_idxs = d_beam_data[level_number].petsc_global_next_node_idxs;
    std::vector<int>& petsc_global_prev_node_idxs = d_beam_data[level_number].petsc_global_prev_node_idxs;

    std::vector<const double*>& rigidities = d_beam_data[level_number].rigidities;
    std::vector<const Vector*>& curvatures = d_beam_data[level_number].curvatures;

    // The LMesh object provides the set of local Lagrangian nodes.
    const Pointer<LMesh> mesh = l_data_manager->getLMesh(level_number);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

    // Determine how many beams are associated with the present MPI process.
    unsigned int total_num_beams = 0;
    for (const auto& node_idx : local_nodes)
    {
        const IBBeamForceSpec* const force_spec = node_idx->getNodeDataItem<IBBeamForceSpec>();
        if (force_spec) total_num_beams += force_spec->getNumberOfBeams();
    }
    petsc_mastr_node_idxs.resize(total_num_beams);
    petsc_next_node_idxs.resize(total_num_beams);
    petsc_prev_node_idxs.resize(total_num_beams);
    petsc_global_mastr_node_idxs.resize(total_num_beams);
    petsc_global_next_node_idxs.resize(total_num_beams);
    petsc_global_prev_node_idxs.resize(total_num_beams);
    rigidities.resize(total_num_beams);
    curvatures.resize(total_num_beams);

    // Setup the data structures used to compute beam forces.
    int current_beam = 0;
    for (const auto& node_idx : local_nodes)
    {
        const IBBeamForceSpec* const force_spec = node_idx->getNodeDataItem<IBBeamForceSpec>();
        if (!force_spec) continue;

#if !defined(NDEBUG)
        const int lag_idx = node_idx->getLagrangianIndex();
        TBOX_ASSERT(lag_idx == force_spec->getMasterNodeIndex());
#endif
        const int petsc_idx = node_idx->getGlobalPETScIndex();
        const std::vector<std::pair<int, int> >& nghbrs = force_spec->getNeighborNodeIndices();
        const std::vector<double>& bend = force_spec->getBendingRigidities();
        const std::vector<Vector>& curv = force_spec->getMeshDependentCurvatures();
        const unsigned int num_beams = force_spec->getNumberOfBeams();
#if !defined(NDEBUG)
        TBOX_ASSERT(num_beams == nghbrs.size());
        TBOX_ASSERT(num_beams == bend.size());
        TBOX_ASSERT(num_beams == curv.size());
#endif
        for (unsigned int k = 0; k < num_beams; ++k)
        {
            petsc_mastr_node_idxs[current_beam] = petsc_idx;
            petsc_next_node_idxs[current_beam] = nghbrs[k].first;
            petsc_prev_node_idxs[current_beam] = nghbrs[k].second;
            rigidities[current_beam] = &bend[k];
            curvatures[current_beam] = &curv[k];
            ++current_beam;
        }
    }

    // Map the Lagrangian neighbor node indices to the PETSc indices
    // corresponding to the present data distribution.
    l_data_manager->mapLagrangianToPETSc(petsc_next_node_idxs, level_number);
    l_data_manager->mapLagrangianToPETSc(petsc_prev_node_idxs, level_number);

    // Keep a copy of global PETSc indices.
    petsc_global_mastr_node_idxs = petsc_mastr_node_idxs;
    petsc_global_next_node_idxs = petsc_next_node_idxs;
    petsc_global_prev_node_idxs = petsc_prev_node_idxs;

    // Determine the global node offset and the number of local nodes.
    const int global_node_offset = l_data_manager->getGlobalNodeOffset(level_number);
    const int num_local_nodes = l_data_manager->getNumberOfLocalNodes(level_number);

    // Determine the ghost nodes required to compute beam forces.
    //
    // NOTE: Only neighbor nodes can be "off processor".  Master nodes are
    // guaranteed to be "on processor".
    for (int idx : petsc_next_node_idxs)
    {
        if (idx < global_node_offset || idx >= global_node_offset + num_local_nodes)
        {
            nonlocal_petsc_idx_set.insert(idx);
        }
    }
    for (int idx : petsc_prev_node_idxs)
    {
        if (idx < global_node_offset || idx >= global_node_offset + num_local_nodes)
        {
            nonlocal_petsc_idx_set.insert(idx);
        }
    }
    return;
} // initializeBeamLevelData

void
IBStandardForceGen::computeLagrangianBeamForce(Pointer<LData> F_data,
                                               Pointer<LData> X_data,
                                               const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                               const int level_number,
                                               const double /*data_time*/,
                                               LDataManager* const /*l_data_manager*/)
{
    const int num_beams = static_cast<int>(d_beam_data[level_number].petsc_mastr_node_idxs.size());
    const bool uses_beams = (num_beams > 0);
    const int* const petsc_mastr_node_idxs = uses_beams ? &d_beam_data[level_number].petsc_mastr_node_idxs[0] : nullptr;
    const int* const petsc_next_node_idxs = uses_beams ? &d_beam_data[level_number].petsc_next_node_idxs[0] : nullptr;
    const int* const petsc_prev_node_idxs = uses_beams ? &d_beam_data[level_number].petsc_prev_node_idxs[0] : nullptr;
    const double** const rigidities = uses_beams ? &d_beam_data[level_number].rigidities[0] : nullptr;
    const Vector** const curvatures = uses_beams ? &d_beam_data[level_number].curvatures[0] : nullptr;
    double* const F_node = F_data->getLocalFormVecArray()->data();
    const double* const X_node = X_data->getGhostedLocalFormVecArray()->data();

    static const int BLOCKSIZE = 16; // This parameter needs to be tuned.
    int k, kblock, kunroll, mastr_idx, next_idx, prev_idx;
    double K;
    const double* D2X0;
    double F[NDIM];
    kblock = 0;
    for (; kblock < (num_beams - 1) / BLOCKSIZE;
         ++kblock) // ensure that the last block is NOT handled by this first loop
    {
        PREFETCH_READ_NTA_BLOCK(petsc_mastr_node_idxs + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(petsc_next_node_idxs + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(petsc_prev_node_idxs + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(rigidities + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(curvatures + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        for (kunroll = 0; kunroll < BLOCKSIZE; ++kunroll)
        {
            k = kblock * BLOCKSIZE + kunroll;
            mastr_idx = petsc_mastr_node_idxs[k];
            next_idx = petsc_next_node_idxs[k];
            prev_idx = petsc_prev_node_idxs[k];
#if !defined(NDEBUG)
            TBOX_ASSERT(mastr_idx != next_idx);
            TBOX_ASSERT(mastr_idx != prev_idx);
#endif
            PREFETCH_READ_NTA_NDIM_BLOCK(F_node + petsc_mastr_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(F_node + petsc_next_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(F_node + petsc_prev_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(X_node + petsc_mastr_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(X_node + petsc_next_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(X_node + petsc_prev_node_idxs[k + 1]);
            PREFETCH_READ_NTA(rigidities[k + 1]);
            PREFETCH_READ_NTA(curvatures[k + 1]);
            K = *rigidities[k];
            D2X0 = curvatures[k]->data();
            F[0] = K * (X_node[next_idx + 0] + X_node[prev_idx + 0] - 2.0 * X_node[mastr_idx + 0] - D2X0[0]);
            F[1] = K * (X_node[next_idx + 1] + X_node[prev_idx + 1] - 2.0 * X_node[mastr_idx + 1] - D2X0[1]);
#if (NDIM == 3)
            F[2] = K * (X_node[next_idx + 2] + X_node[prev_idx + 2] - 2.0 * X_node[mastr_idx + 2] - D2X0[2]);
#endif
            F_node[mastr_idx + 0] += 2.0 * F[0];
            F_node[mastr_idx + 1] += 2.0 * F[1];
#if (NDIM == 3)
            F_node[mastr_idx + 2] += 2.0 * F[2];
#endif
            F_node[next_idx + 0] -= F[0];
            F_node[next_idx + 1] -= F[1];
#if (NDIM == 3)
            F_node[next_idx + 2] -= F[2];
#endif
            F_node[prev_idx + 0] -= F[0];
            F_node[prev_idx + 1] -= F[1];
#if (NDIM == 3)
            F_node[prev_idx + 2] -= F[2];
#endif
        }
    }
    for (k = kblock * BLOCKSIZE; k < num_beams; ++k)
    {
        mastr_idx = petsc_mastr_node_idxs[k];
        next_idx = petsc_next_node_idxs[k];
        prev_idx = petsc_prev_node_idxs[k];
#if !defined(NDEBUG)
        TBOX_ASSERT(mastr_idx != next_idx);
        TBOX_ASSERT(mastr_idx != prev_idx);
#endif
        K = *rigidities[k];
        D2X0 = curvatures[k]->data();
        F[0] = K * (X_node[next_idx + 0] + X_node[prev_idx + 0] - 2.0 * X_node[mastr_idx + 0] - D2X0[0]);
        F[1] = K * (X_node[next_idx + 1] + X_node[prev_idx + 1] - 2.0 * X_node[mastr_idx + 1] - D2X0[1]);
#if (NDIM == 3)
        F[2] = K * (X_node[next_idx + 2] + X_node[prev_idx + 2] - 2.0 * X_node[mastr_idx + 2] - D2X0[2]);
#endif
        F_node[mastr_idx + 0] += 2.0 * F[0];
        F_node[mastr_idx + 1] += 2.0 * F[1];
#if (NDIM == 3)
        F_node[mastr_idx + 2] += 2.0 * F[2];
#endif
        F_node[next_idx + 0] -= F[0];
        F_node[next_idx + 1] -= F[1];
#if (NDIM == 3)
        F_node[next_idx + 2] -= F[2];
#endif
        F_node[prev_idx + 0] -= F[0];
        F_node[prev_idx + 1] -= F[1];
#if (NDIM == 3)
        F_node[prev_idx + 2] -= F[2];
#endif
    }

    F_data->restoreArrays();
    X_data->restoreArrays();
    return;
} // computeLagrangianBeamForce

void
IBStandardForceGen::initializeTargetPointLevelData(std::set<int>& /*nonlocal_petsc_idx_set*/,
                                                   const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                                   const int level_number,
                                                   const double /*init_data_time*/,
                                                   const bool /*initial_time*/,
                                                   LDataManager* const l_data_manager)
{
    std::vector<int>& petsc_node_idxs = d_target_point_data[level_number].petsc_node_idxs;
    std::vector<int>& petsc_global_node_idxs = d_target_point_data[level_number].petsc_global_node_idxs;
    std::vector<const double*>& kappa = d_target_point_data[level_number].kappa;
    std::vector<const double*>& eta = d_target_point_data[level_number].eta;
    std::vector<const Point*>& X0 = d_target_point_data[level_number].X0;

    // The LMesh object provides the set of local Lagrangian nodes.
    const Pointer<LMesh> mesh = l_data_manager->getLMesh(level_number);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

    // Determine how many target points are associated with the present MPI
    // process.
    unsigned int total_num_target_points = 0;
    for (const auto& node_idx : local_nodes)
    {
        const IBTargetPointForceSpec* const force_spec = node_idx->getNodeDataItem<IBTargetPointForceSpec>();
        if (force_spec) total_num_target_points += 1;
    }

    // Resize arrays for storing cached values used to compute target point
    // forces.
    petsc_node_idxs.resize(total_num_target_points);
    petsc_global_node_idxs.resize(total_num_target_points);
    kappa.resize(total_num_target_points);
    eta.resize(total_num_target_points);
    X0.resize(total_num_target_points);

    // Setup the data structures used to compute target point forces.
    int current_target_point = 0;
    for (const auto& node_idx : local_nodes)
    {
        const IBTargetPointForceSpec* const force_spec = node_idx->getNodeDataItem<IBTargetPointForceSpec>();
        if (!force_spec) continue;
        petsc_global_node_idxs[current_target_point] = petsc_node_idxs[current_target_point] =
            node_idx->getGlobalPETScIndex();
        kappa[current_target_point] = &force_spec->getStiffness();
        eta[current_target_point] = &force_spec->getDamping();
        X0[current_target_point] = &force_spec->getTargetPointPosition();
        ++current_target_point;
    }

    return;
} // initializeTargetPointLevelData

void
IBStandardForceGen::computeLagrangianTargetPointForce(Pointer<LData> F_data,
                                                      Pointer<LData> X_data,
                                                      Pointer<LData> U_data,
                                                      const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                                      const int level_number,
                                                      const double /*data_time*/,
                                                      LDataManager* const /*l_data_manager*/)
{
    double max_displacement = 0.0;

    const int num_target_points = static_cast<int>(d_target_point_data[level_number].petsc_node_idxs.size());
    const bool uses_target_points = (num_target_points > 0);
    const int* const petsc_node_idxs =
        uses_target_points ? &d_target_point_data[level_number].petsc_node_idxs[0] : nullptr;
    const double** const kappa = uses_target_points ? &d_target_point_data[level_number].kappa[0] : nullptr;
    const double** const eta = uses_target_points ? &d_target_point_data[level_number].eta[0] : nullptr;
    const Point** const X0 = uses_target_points ? &d_target_point_data[level_number].X0[0] : nullptr;
    double* const F_node = F_data->getLocalFormVecArray()->data();
    const double* const X_node = X_data->getLocalFormVecArray()->data();
    const double* const U_node = U_data->getLocalFormVecArray()->data();

    static const int BLOCKSIZE = 16; // This parameter needs to be tuned.
    int k, kblock, kunroll, idx;
    double K, E, dX;
    const double* X_target;
    kblock = 0;
    for (; kblock < (num_target_points - 1) / BLOCKSIZE;
         ++kblock) // ensure that the last block is NOT handled by this first loop
    {
        PREFETCH_READ_NTA_BLOCK(petsc_node_idxs + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(kappa + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(eta + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(X0 + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        for (kunroll = 0; kunroll < BLOCKSIZE; ++kunroll)
        {
            k = kblock * BLOCKSIZE + kunroll;
            idx = petsc_node_idxs[k];
            PREFETCH_READ_NTA_NDIM_BLOCK(F_node + petsc_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(X_node + petsc_node_idxs[k + 1]);
            PREFETCH_READ_NTA(kappa[k + 1]);
            PREFETCH_READ_NTA(eta[k + 1]);
            PREFETCH_READ_NTA(X0[k + 1]);
            K = *kappa[k];
            E = *eta[k];
            X_target = X0[k]->data();
            F_node[idx + 0] += K * (X_target[0] - X_node[idx + 0]) - E * U_node[idx + 0];
            F_node[idx + 1] += K * (X_target[1] - X_node[idx + 1]) - E * U_node[idx + 1];
#if (NDIM == 3)
            F_node[idx + 2] += K * (X_target[2] - X_node[idx + 2]) - E * U_node[idx + 2];
#endif
            if (d_log_target_point_displacements)
            {
                dX = 0.0;
                dX += (X_target[0] - X_node[idx + 0]) * (X_target[0] - X_node[idx + 0]);
                dX += (X_target[1] - X_node[idx + 1]) * (X_target[1] - X_node[idx + 1]);
#if (NDIM == 3)
                dX += (X_target[2] - X_node[idx + 2]) * (X_target[2] - X_node[idx + 2]);
#endif
                dX = std::sqrt(dX);
                max_displacement = std::max(max_displacement, dX);
            }
        }
    }
    for (k = kblock * BLOCKSIZE; k < num_target_points; ++k)
    {
        idx = petsc_node_idxs[k];
        K = *kappa[k];
        E = *eta[k];
        X_target = X0[k]->data();
        F_node[idx + 0] += K * (X_target[0] - X_node[idx + 0]) - E * U_node[idx + 0];
        F_node[idx + 1] += K * (X_target[1] - X_node[idx + 1]) - E * U_node[idx + 1];
#if (NDIM == 3)
        F_node[idx + 2] += K * (X_target[2] - X_node[idx + 2]) - E * U_node[idx + 2];
#endif
        if (d_log_target_point_displacements)
        {
            dX = 0.0;
            dX += (X_target[0] - X_node[idx + 0]) * (X_target[0] - X_node[idx + 0]);
            dX += (X_target[1] - X_node[idx + 1]) * (X_target[1] - X_node[idx + 1]);
#if (NDIM == 3)
            dX += (X_target[2] - X_node[idx + 2]) * (X_target[2] - X_node[idx + 2]);
#endif
            dX = std::sqrt(dX);
            max_displacement = std::max(max_displacement, dX);
        }
    }

    if (d_log_target_point_displacements)
    {
        max_displacement = IBTK_MPI::maxReduction(max_displacement);
        plog << "IBStandardForceGen: maximum target point displacement: " << max_displacement << "\n";
    }

    F_data->restoreArrays();
    X_data->restoreArrays();
    U_data->restoreArrays();
    return;
} // computeLagrangianTargetPointForce

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
