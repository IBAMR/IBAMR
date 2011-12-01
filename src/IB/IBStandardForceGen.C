// Filename: IBStandardForceGen.C
// Created on 03 May 2005 by Boyce Griffith
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

#include "IBStandardForceGen.h"

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
#include <ibamr/IBSpringForceSpec.h>
#include <ibamr/IBTargetPointForceSpec.h>
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/compiler_hints.h>

// SAMRAI INCLUDES
#include <tbox/TimerManager.h>

// C++ STDLIB INCLUDES
#include <algorithm>

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
static Timer* t_compute_lagrangian_energy;
static Timer* t_initialize_level_data;

void
resetLocalPETScIndices(
    blitz::Array<int,1>& inds,
    const int global_node_offset,
    const int num_local_nodes)
{
#ifndef DEBUG_CHECK_ASSERTIONS
    NULL_USE(num_local_nodes);
#endif
    for (blitz::Array<int,1>::iterator it = inds.begin(); it !=inds.end(); ++it)
    {
        int& idx = *it;
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(idx >= global_node_offset && idx < global_node_offset+num_local_nodes);
#endif
        idx -= global_node_offset;
    }
    return;
}// resetLocalPETScIndices

void
resetLocalOrNonlocalPETScIndices(
    blitz::Array<int,1>& inds,
    const int global_node_offset,
    const int num_local_nodes,
    const std::vector<int>& nonlocal_petsc_idxs)
{
    for (blitz::Array<int,1>::iterator it = inds.begin(); it != inds.end(); ++it)
    {
        int& idx = *it;
        if (idx >= global_node_offset && idx < global_node_offset+num_local_nodes)
        {
            // A local node.
            idx -= global_node_offset;
        }
        else
        {
            // A nonlocal node.
            //
            // First, lookup the slave node index in the set of ghost nodes.
            const std::vector<int>::const_iterator posn = std::lower_bound(nonlocal_petsc_idxs.begin(), nonlocal_petsc_idxs.end(), idx);
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(idx == *posn);
#endif
            // Second, set the local index via the offset of the ghost node
            // index within the set of ghost nodes.
            idx = num_local_nodes + std::distance(nonlocal_petsc_idxs.begin(), posn);
        }
    }
    return;
}// resetLocalOrNonlocalPETScIndices
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBStandardForceGen::IBStandardForceGen(
    const bool constant_material_properties)
    : d_constant_material_properties(constant_material_properties)
{
    if (d_constant_material_properties)
    {
        pout << "IBStandardForceGen:  Using constant material properites.\n"
             << "NOTICE:  Material properties may only be updated at the regrid interval.\n\n";
    }
    else
    {
        pout << "IBStandardForceGen:  Using nonconstant material properites.\n"
             << "NOTICE:  Material properties may be updated at any time.\n"
             << "         Models with constant material properties may see increased performance by setting constant_material_properties = true.\n\n";
    }

    // Setup the default force generation functions.
    registerSpringForceFunction(0, &default_linear_spring_force);

    // Setup Timers.
    IBAMR_DO_ONCE(
        t_compute_lagrangian_force                            = TimerManager::getManager()->getTimer("IBAMR::IBStandardForceGen::computeLagrangianForce()");
        t_compute_lagrangian_force_jacobian                   = TimerManager::getManager()->getTimer("IBAMR::IBStandardForceGen::computeLagrangianForceJacobian()");
        t_compute_lagrangian_force_jacobian_nonzero_structure = TimerManager::getManager()->getTimer("IBAMR::IBStandardForceGen::computeLagrangianForceJacobianNonzeroStructure()");
        t_initialize_level_data                               = TimerManager::getManager()->getTimer("IBAMR::IBStandardForceGen::initializeLevelData()");
        t_compute_lagrangian_energy                           = TimerManager::getManager()->getTimer("IBAMR::IBStandardForceGen::computeLagrangianEnergy()");
                  );
    return;
}// IBStandardForceGen

IBStandardForceGen::~IBStandardForceGen()
{
    // intentionally blank
    return;
}// ~IBStandardForceGen

void
IBStandardForceGen::registerSpringForceFunction(
    const int force_fcn_index,
    const SpringForceFcnPtr spring_force_fcn_ptr)
{
    d_spring_force_fcn_map[force_fcn_index] = spring_force_fcn_ptr;
    return;
}// registerSpringForceFunction

void
IBStandardForceGen::initializeLevelData(
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool initial_time,
    LDataManager* const l_data_manager)
{
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return;

    IBAMR_TIMER_START(t_initialize_level_data);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
#endif
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Resize the vectors corresponding to data individually maintained for
    // separate levels of the patch hierarchy.
    const int new_size = std::max(level_number+1, static_cast<int>(d_is_initialized.size()));

    d_spring_data      .resize(new_size);
    d_beam_data        .resize(new_size);
    d_target_point_data.resize(new_size);
    d_X_ghost_data     .resize(new_size);
    d_F_ghost_data     .resize(new_size);
    d_is_initialized   .resize(new_size, false);

    // Keep track of all of the nonlocal PETSc indices required to compute the
    // forces.
    std::set<int> nonlocal_petsc_idx_set;

    // Setup the cached data.
    initializeSpringLevelData(     nonlocal_petsc_idx_set, hierarchy, level_number, init_data_time, initial_time, l_data_manager);
    initializeBeamLevelData(       nonlocal_petsc_idx_set, hierarchy, level_number, init_data_time, initial_time, l_data_manager);
    initializeTargetPointLevelData(nonlocal_petsc_idx_set, hierarchy, level_number, init_data_time, initial_time, l_data_manager);

    // Put the nonlocal PETSc indices into a vector.
    std::vector<int> nonlocal_petsc_idxs(nonlocal_petsc_idx_set.begin(),nonlocal_petsc_idx_set.end());

    // Put all cached PETSc node indices into local form.
    const int global_node_offset = l_data_manager->getGlobalNodeOffset(level_number);
    const int num_local_nodes = l_data_manager->getNumberOfLocalNodes(level_number);

    resetLocalPETScIndices(          d_spring_data      [level_number].petsc_mastr_node_idxs, global_node_offset, num_local_nodes);
    resetLocalPETScIndices(          d_beam_data        [level_number].petsc_mastr_node_idxs, global_node_offset, num_local_nodes);
    resetLocalPETScIndices(          d_target_point_data[level_number].petsc_node_idxs      , global_node_offset, num_local_nodes);

    resetLocalOrNonlocalPETScIndices(d_spring_data      [level_number].petsc_slave_node_idxs, global_node_offset, num_local_nodes, nonlocal_petsc_idxs);
    resetLocalOrNonlocalPETScIndices(d_beam_data        [level_number].petsc_next_node_idxs , global_node_offset, num_local_nodes, nonlocal_petsc_idxs);
    resetLocalOrNonlocalPETScIndices(d_beam_data        [level_number].petsc_prev_node_idxs , global_node_offset, num_local_nodes, nonlocal_petsc_idxs);

    std::ostringstream X_name_stream;
    X_name_stream << "IBStandardForceGen::X_ghost_" << level_number;
    d_X_ghost_data[level_number] = new LData(X_name_stream.str(), num_local_nodes, NDIM, nonlocal_petsc_idxs);

    std::ostringstream F_name_stream;
    F_name_stream << "IBStandardForceGen::F_ghost_" << level_number;
    d_F_ghost_data[level_number] = new LData(F_name_stream.str(), num_local_nodes, NDIM, nonlocal_petsc_idxs);

    // Transform all of the cached indices to correspond to a data depth of
    // NDIM.
    std::transform(d_spring_data[level_number].petsc_mastr_node_idxs.begin(), d_spring_data[level_number].petsc_mastr_node_idxs.end(),
                   d_spring_data[level_number].petsc_mastr_node_idxs.begin(), std::bind2nd(std::multiplies<int>(),NDIM));
    std::transform(d_spring_data[level_number].petsc_slave_node_idxs.begin(), d_spring_data[level_number].petsc_slave_node_idxs.end(),
                   d_spring_data[level_number].petsc_slave_node_idxs.begin(), std::bind2nd(std::multiplies<int>(),NDIM));
    std::transform(d_beam_data[level_number].petsc_mastr_node_idxs.begin(), d_beam_data[level_number].petsc_mastr_node_idxs.end(),
                   d_beam_data[level_number].petsc_mastr_node_idxs.begin(), std::bind2nd(std::multiplies<int>(),NDIM));
    std::transform(d_beam_data[level_number].petsc_next_node_idxs.begin(), d_beam_data[level_number].petsc_next_node_idxs.end(),
                   d_beam_data[level_number].petsc_next_node_idxs.begin(), std::bind2nd(std::multiplies<int>(),NDIM));
    std::transform(d_beam_data[level_number].petsc_prev_node_idxs.begin(), d_beam_data[level_number].petsc_prev_node_idxs.end(),
                   d_beam_data[level_number].petsc_prev_node_idxs.begin(), std::bind2nd(std::multiplies<int>(),NDIM));
    std::transform(d_target_point_data[level_number].petsc_node_idxs.begin(), d_target_point_data[level_number].petsc_node_idxs.end(),
                   d_target_point_data[level_number].petsc_node_idxs.begin(), std::bind2nd(std::multiplies<int>(),NDIM));

    // Indicate that the level data has been initialized.
    d_is_initialized[level_number] = true;

    IBAMR_TIMER_STOP(t_initialize_level_data);
    return;
}// initializeLevelData

void
IBStandardForceGen::computeLagrangianForce(
    Pointer<LData> F_data,
    Pointer<LData> X_data,
    Pointer<LData> U_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const l_data_manager)
{
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return;

    IBAMR_TIMER_START(t_compute_lagrangian_force);

    int ierr;

    // Initialize ghost data.
    Pointer<LData> F_ghost_data = d_F_ghost_data[level_number];
    Vec F_ghost_local_form_vec;
    ierr = VecGhostGetLocalForm(F_ghost_data->getVec(), &F_ghost_local_form_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecSet(F_ghost_local_form_vec, 0.0);  IBTK_CHKERRQ(ierr);
    ierr = VecGhostRestoreLocalForm(F_ghost_data->getVec(), &F_ghost_local_form_vec);  IBTK_CHKERRQ(ierr);

    Pointer<LData> X_ghost_data = d_X_ghost_data[level_number];
    ierr = VecCopy(X_data->getVec(), X_ghost_data->getVec());  IBTK_CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(X_ghost_data->getVec(), INSERT_VALUES, SCATTER_FORWARD);  IBTK_CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(  X_ghost_data->getVec(), INSERT_VALUES, SCATTER_FORWARD);  IBTK_CHKERRQ(ierr);

    // Compute the forces.
    computeLagrangianSpringForce(       F_ghost_data, X_ghost_data,         hierarchy, level_number, data_time, l_data_manager);
    computeLagrangianBeamForce(         F_ghost_data, X_ghost_data,         hierarchy, level_number, data_time, l_data_manager);
    computeLagrangianTargetPointForce(  F_ghost_data, X_ghost_data, U_data, hierarchy, level_number, data_time, l_data_manager);

    // Add the locally computed forces to the Lagrangian force vector.
    //
    // WARNING: The following operations may yield nondeterministic results in
    // parallel environments (i.e., the order of summation may not be
    // consistent).
    ierr = VecGhostUpdateBegin(F_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);  IBTK_CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(  F_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);  IBTK_CHKERRQ(ierr);
    ierr = VecAXPY(F_data->getVec(), 1.0, F_ghost_data->getVec());

    IBAMR_TIMER_STOP(t_compute_lagrangian_force);
    return;
}// computeLagrangianForce

void
IBStandardForceGen::computeLagrangianForceJacobianNonzeroStructure(
    std::vector<int>& /*d_nnz*/,
    std::vector<int>& /*o_nnz*/,
    const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    const int level_number,
    const double /*data_time*/,
    LDataManager* const l_data_manager)
{
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return;

    // Compute the nonzero structure of the Jacobian matrix.
    TBOX_ERROR("not currently implemented\n");
    return;
}// computeLagrangianForceJacobianNonzeroStructure

void
IBStandardForceGen::computeLagrangianForceJacobian(
    Mat& /*J_mat*/,
    MatAssemblyType /*assembly_type*/,
    const double /*X_coef*/,
    Pointer<LData> /*X_data*/,
    const double /*U_coef*/,
    Pointer<LData> /*U_data*/,
    const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    const int level_number,
    const double /*data_time*/,
    LDataManager* const l_data_manager)
{
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return;

    // Compute the Jacobian matrix.
    TBOX_ERROR("not currently implemented\n");
    return;
}// computeLagrangianForceJacobian

double
IBStandardForceGen::computeLagrangianEnergy(
    Pointer<LData> /*X_data*/,
    Pointer<LData> /*U_data*/,
    const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    const int level_number,
    const double /*data_time*/,
    LDataManager* const l_data_manager)
{
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return 0.0;

    // Compute the energy.
    TBOX_ERROR("not currently implemented\n");
    double ret_val = std::numeric_limits<double>::quiet_NaN();
    return ret_val;
}// computeLagrangianEnergy

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBStandardForceGen::initializeSpringLevelData(
    std::set<int>& nonlocal_petsc_idx_set,
    const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    const int level_number,
    const double /*init_data_time*/,
    const bool /*initial_time*/,
    LDataManager* const l_data_manager)
{
    blitz::Array<int,1>&            lag_mastr_node_idxs = d_spring_data[level_number].lag_mastr_node_idxs;
    blitz::Array<int,1>&            lag_slave_node_idxs = d_spring_data[level_number].lag_slave_node_idxs;
    blitz::Array<int,1>&          petsc_mastr_node_idxs = d_spring_data[level_number].petsc_mastr_node_idxs;
    blitz::Array<int,1>&          petsc_slave_node_idxs = d_spring_data[level_number].petsc_slave_node_idxs;
    blitz::Array<SpringForceFcnPtr,1>&       force_fcns = d_spring_data[level_number].force_fcns;
    blitz::Array<double,1>&                 stiffnesses = d_spring_data[level_number].stiffnesses;
    blitz::Array<double,1>&                rest_lengths = d_spring_data[level_number].rest_lengths;
    blitz::Array<const double*,1>&  dynamic_stiffnesses = d_spring_data[level_number].dynamic_stiffnesses;
    blitz::Array<const double*,1>& dynamic_rest_lengths = d_spring_data[level_number].dynamic_rest_lengths;

    // The LMesh object provides the set of local Lagrangian nodes.
    const Pointer<LMesh> mesh = l_data_manager->getLMesh(level_number);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
    const int num_local_nodes = local_nodes.size();

    // Quick return if local_nodes is empty.
    if (local_nodes.empty())
    {
        static const int num_springs = 0;
        lag_mastr_node_idxs  .resize(num_springs);
        lag_slave_node_idxs  .resize(num_springs);
        petsc_mastr_node_idxs.resize(num_springs);
        petsc_slave_node_idxs.resize(num_springs);
        force_fcns           .resize(num_springs);
        if (d_constant_material_properties)
        {
            stiffnesses         .resize(num_springs);
            rest_lengths        .resize(num_springs);
            dynamic_stiffnesses .resize(0);
            dynamic_rest_lengths.resize(0);
        }
        else
        {
            stiffnesses         .resize(0);
            rest_lengths        .resize(0);
            dynamic_stiffnesses .resize(num_springs);
            dynamic_rest_lengths.resize(num_springs);
        }
        return;
    }

    // Determine how many springs are associated with the present MPI process.
    int num_springs = 0;
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
        const IBSpringForceSpec* const force_spec = node_idx->getNodeDataItem<IBSpringForceSpec>();
        if (force_spec != NULL) num_springs += force_spec->getNumberOfSprings();
    }

    // Resize arrays for storing cached values used to compute spring forces.
    lag_mastr_node_idxs  .resize(num_springs);
    lag_slave_node_idxs  .resize(num_springs);
    petsc_mastr_node_idxs.resize(num_springs);
    petsc_slave_node_idxs.resize(num_springs);
    force_fcns           .resize(num_springs);
    if (d_constant_material_properties)
    {
        stiffnesses         .resize(num_springs);
        rest_lengths        .resize(num_springs);
        dynamic_stiffnesses .resize(0);
        dynamic_rest_lengths.resize(0);
    }
    else
    {
        stiffnesses         .resize(0);
        rest_lengths        .resize(0);
        dynamic_stiffnesses .resize(num_springs);
        dynamic_rest_lengths.resize(num_springs);
    }

    // Return early if there are no local springs.
    if (num_springs == 0) return;

    // Setup the data structures used to compute spring forces.
    int current_spring = 0;
    if (d_constant_material_properties)
    {
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const IBSpringForceSpec* const force_spec = node_idx->getNodeDataItem<IBSpringForceSpec>();
            if (force_spec == NULL) continue;
            const int lag_idx = node_idx->getLagrangianIndex();
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(lag_idx == force_spec->getMasterNodeIndex());
#endif
            const int petsc_idx = node_idx->getGlobalPETScIndex();
            const std::vector<int>& slv = force_spec->getSlaveNodeIndices();
            const std::vector<int>& fcn = force_spec->getForceFunctionIndices();
            const std::vector<double>& stf = force_spec->getStiffnesses();
            const std::vector<double>& rst = force_spec->getRestingLengths();
            const unsigned int num_springs = force_spec->getNumberOfSprings();
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(num_springs == slv.size());
            TBOX_ASSERT(num_springs == fcn.size());
            TBOX_ASSERT(num_springs == stf.size());
            TBOX_ASSERT(num_springs == rst.size());
#endif
            for (unsigned int k = 0; k < num_springs; ++k)
            {
                lag_mastr_node_idxs  (current_spring) = lag_idx;
                lag_slave_node_idxs  (current_spring) = slv[k];
                petsc_mastr_node_idxs(current_spring) = petsc_idx;
                force_fcns           (current_spring) = d_spring_force_fcn_map[fcn[k]];
                stiffnesses          (current_spring) = stf[k];
                rest_lengths         (current_spring) = rst[k];
                ++current_spring;
            }
        }
    }
    else
    {
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const IBSpringForceSpec* const force_spec = node_idx->getNodeDataItem<IBSpringForceSpec>();
            if (force_spec == NULL) continue;

            const int lag_idx = node_idx->getLagrangianIndex();
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(lag_idx == force_spec->getMasterNodeIndex());
#endif
            const int petsc_idx = node_idx->getGlobalPETScIndex();
            const std::vector<int>& slv = force_spec->getSlaveNodeIndices();
            const std::vector<int>& fcn = force_spec->getForceFunctionIndices();
            const std::vector<double>& stf = force_spec->getStiffnesses();
            const std::vector<double>& rst = force_spec->getRestingLengths();
            const unsigned int num_springs = force_spec->getNumberOfSprings();
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(num_springs == slv.size());
            TBOX_ASSERT(num_springs == fcn.size());
            TBOX_ASSERT(num_springs == stf.size());
            TBOX_ASSERT(num_springs == rst.size());
#endif
            for (unsigned int k = 0; k < num_springs; ++k)
            {
                lag_mastr_node_idxs  (current_spring) = lag_idx;
                lag_slave_node_idxs  (current_spring) = slv[k];
                petsc_mastr_node_idxs(current_spring) = petsc_idx;
                force_fcns           (current_spring) = d_spring_force_fcn_map[fcn[k]];
                dynamic_stiffnesses  (current_spring) = &stf[k];
                dynamic_rest_lengths (current_spring) = &rst[k];
                ++current_spring;
            }
        }
    }

    // Map the Lagrangian slave node indices to the PETSc indices corresponding
    // to the present data distribution.
    petsc_slave_node_idxs = lag_slave_node_idxs;
    l_data_manager->mapLagrangianToPETSc(petsc_slave_node_idxs, level_number);

    // Determine the ghost nodes required to compute spring forces.
    //
    // NOTE: Only slave nodes can be "off processor".  Master nodes are
    // guaranteed to be "on processor".
    const int global_node_offset = l_data_manager->getGlobalNodeOffset(level_number);
    for (int k = 0; k < petsc_slave_node_idxs.size(); ++k)
    {
        const int idx = petsc_slave_node_idxs(k);
        if (UNLIKELY(idx < global_node_offset || idx >= global_node_offset+num_local_nodes))
        {
            nonlocal_petsc_idx_set.insert(idx);
        }
    }
    return;
}// initializeSpringLevelData

void
IBStandardForceGen::computeLagrangianSpringForce(
    Pointer<LData> F_data,
    Pointer<LData> X_data,
    const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    const int level_number,
    const double /*data_time*/,
    LDataManager* const /*l_data_manager*/)
{
    const int num_springs = d_spring_data[level_number].lag_mastr_node_idxs.size();
    const int*               const restrict   lag_mastr_node_idxs = d_spring_data[level_number].lag_mastr_node_idxs  .data();
    const int*               const restrict   lag_slave_node_idxs = d_spring_data[level_number].lag_slave_node_idxs  .data();
    const int*               const restrict petsc_mastr_node_idxs = d_spring_data[level_number].petsc_mastr_node_idxs.data();
    const int*               const restrict petsc_slave_node_idxs = d_spring_data[level_number].petsc_slave_node_idxs.data();
    const SpringForceFcnPtr* const restrict            force_fcns = d_spring_data[level_number].force_fcns           .data();
    const double*            const restrict           stiffnesses = d_spring_data[level_number].stiffnesses          .data();
    const double*            const restrict          rest_lengths = d_spring_data[level_number].rest_lengths         .data();
    const double**           const restrict   dynamic_stiffnesses = d_spring_data[level_number].dynamic_stiffnesses  .data();
    const double**           const restrict  dynamic_rest_lengths = d_spring_data[level_number].dynamic_rest_lengths .data();
    double*                  const restrict                F_node = F_data->getGhostedLocalFormVecArray()->data();
    const double*            const restrict                X_node = X_data->getGhostedLocalFormVecArray()->data();

    static const int BLOCKSIZE = 16;  // This parameter needs to be tuned.
    int k, kblock, kunroll, mastr_idx, slave_idx;
    double F[NDIM], D[NDIM];
    kblock = 0;
    if (d_constant_material_properties)
    {
        for ( ; kblock < (num_springs-1)/BLOCKSIZE; ++kblock)  // ensure that the last block is NOT handled by this first loop
        {
            PREFETCH_READ_NTA_BLOCK(  lag_mastr_node_idxs+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK(  lag_slave_node_idxs+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK(petsc_mastr_node_idxs+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK(petsc_slave_node_idxs+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK(           force_fcns+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK(          stiffnesses+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK(         rest_lengths+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            for (kunroll = 0; kunroll < BLOCKSIZE; ++kunroll)
            {
                k = kblock*BLOCKSIZE+kunroll;
                mastr_idx = petsc_mastr_node_idxs[k];
                slave_idx = petsc_slave_node_idxs[k];
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(mastr_idx != slave_idx);
#endif
                PREFETCH_READ_NTA_NDIM_BLOCK(F_node+petsc_mastr_node_idxs[k+1]);
                PREFETCH_READ_NTA_NDIM_BLOCK(F_node+petsc_slave_node_idxs[k+1]);
                PREFETCH_READ_NTA_NDIM_BLOCK(X_node+petsc_mastr_node_idxs[k+1]);
                PREFETCH_READ_NTA_NDIM_BLOCK(X_node+petsc_slave_node_idxs[k+1]);
                D[0] = X_node[slave_idx+0] - X_node[mastr_idx+0];
                D[1] = X_node[slave_idx+1] - X_node[mastr_idx+1];
#if (NDIM == 3)
                D[2] = X_node[slave_idx+2] - X_node[mastr_idx+2];
#endif
                (force_fcns[k])(F,D,stiffnesses[k],rest_lengths[k],lag_mastr_node_idxs[k],lag_slave_node_idxs[k]);
                F_node[mastr_idx+0] += F[0];
                F_node[mastr_idx+1] += F[1];
#if (NDIM == 3)
                F_node[mastr_idx+2] += F[2];
#endif
                F_node[slave_idx+0] -= F[0];
                F_node[slave_idx+1] -= F[1];
#if (NDIM == 3)
                F_node[slave_idx+2] -= F[2];
#endif
            }
        }
        for (k = kblock*BLOCKSIZE; k < num_springs; ++k)
        {
            mastr_idx = petsc_mastr_node_idxs[k];
            slave_idx = petsc_slave_node_idxs[k];
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(mastr_idx != slave_idx);
#endif
            D[0] = X_node[slave_idx+0] - X_node[mastr_idx+0];
            D[1] = X_node[slave_idx+1] - X_node[mastr_idx+1];
#if (NDIM == 3)
            D[2] = X_node[slave_idx+2] - X_node[mastr_idx+2];
#endif
            (force_fcns[k])(F,D,stiffnesses[k],rest_lengths[k],lag_mastr_node_idxs[k],lag_slave_node_idxs[k]);
            F_node[mastr_idx+0] += F[0];
            F_node[mastr_idx+1] += F[1];
#if (NDIM == 3)
            F_node[mastr_idx+2] += F[2];
#endif
            F_node[slave_idx+0] -= F[0];
            F_node[slave_idx+1] -= F[1];
#if (NDIM == 3)
            F_node[slave_idx+2] -= F[2];
#endif
        }
    }
    else
    {
        for ( ; kblock < (num_springs-1)/BLOCKSIZE; ++kblock)  // ensure that the last block is NOT handled by this first loop
        {
            PREFETCH_READ_NTA_BLOCK(  lag_mastr_node_idxs+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK(  lag_slave_node_idxs+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK(petsc_mastr_node_idxs+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK(petsc_slave_node_idxs+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK(           force_fcns+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK(  dynamic_stiffnesses+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK( dynamic_rest_lengths+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            for (kunroll = 0; kunroll < BLOCKSIZE; ++kunroll)
            {
                k = kblock*BLOCKSIZE+kunroll;
                mastr_idx = petsc_mastr_node_idxs[k];
                slave_idx = petsc_slave_node_idxs[k];
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(mastr_idx != slave_idx);
#endif
                PREFETCH_READ_NTA_NDIM_BLOCK(F_node+petsc_mastr_node_idxs[k+1]);
                PREFETCH_READ_NTA_NDIM_BLOCK(F_node+petsc_slave_node_idxs[k+1]);
                PREFETCH_READ_NTA_NDIM_BLOCK(X_node+petsc_mastr_node_idxs[k+1]);
                PREFETCH_READ_NTA_NDIM_BLOCK(X_node+petsc_slave_node_idxs[k+1]);
                PREFETCH_READ_NTA(                    dynamic_stiffnesses[k+1]);
                PREFETCH_READ_NTA(                   dynamic_rest_lengths[k+1]);
                D[0] = X_node[slave_idx+0] - X_node[mastr_idx+0];
                D[1] = X_node[slave_idx+1] - X_node[mastr_idx+1];
#if (NDIM == 3)
                D[2] = X_node[slave_idx+2] - X_node[mastr_idx+2];
#endif
                (force_fcns[k])(F,D,*dynamic_stiffnesses[k],*dynamic_rest_lengths[k],lag_mastr_node_idxs[k],lag_slave_node_idxs[k]);
                F_node[mastr_idx+0] += F[0];
                F_node[mastr_idx+1] += F[1];
#if (NDIM == 3)
                F_node[mastr_idx+2] += F[2];
#endif
                F_node[slave_idx+0] -= F[0];
                F_node[slave_idx+1] -= F[1];
#if (NDIM == 3)
                F_node[slave_idx+2] -= F[2];
#endif
            }
        }
        for (k = kblock*BLOCKSIZE; k < num_springs; ++k)
        {
            mastr_idx = petsc_mastr_node_idxs[k];
            slave_idx = petsc_slave_node_idxs[k];
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(mastr_idx != slave_idx);
#endif
            D[0] = X_node[slave_idx+0] - X_node[mastr_idx+0];
            D[1] = X_node[slave_idx+1] - X_node[mastr_idx+1];
#if (NDIM == 3)
            D[2] = X_node[slave_idx+2] - X_node[mastr_idx+2];
#endif
            (force_fcns[k])(F,D,*dynamic_stiffnesses[k],*dynamic_rest_lengths[k],lag_mastr_node_idxs[k],lag_slave_node_idxs[k]);
            F_node[mastr_idx+0] += F[0];
            F_node[mastr_idx+1] += F[1];
#if (NDIM == 3)
            F_node[mastr_idx+2] += F[2];
#endif
            F_node[slave_idx+0] -= F[0];
            F_node[slave_idx+1] -= F[1];
#if (NDIM == 3)
            F_node[slave_idx+2] -= F[2];
#endif
        }
    }

    F_data->restoreArrays();
    X_data->restoreArrays();
    return;
}// computeLagrangianSpringForce

void
IBStandardForceGen::initializeBeamLevelData(
    std::set<int>& nonlocal_petsc_idx_set,
    const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    const int level_number,
    const double /*init_data_time*/,
    const bool /*initial_time*/,
    LDataManager* const l_data_manager)
{
    blitz::Array<int,1>&                                petsc_mastr_node_idxs = d_beam_data[level_number].petsc_mastr_node_idxs;
    blitz::Array<int,1>&                                 petsc_next_node_idxs = d_beam_data[level_number].petsc_next_node_idxs;
    blitz::Array<int,1>&                                 petsc_prev_node_idxs = d_beam_data[level_number].petsc_prev_node_idxs;
    blitz::Array<double,1>&                                        rigidities = d_beam_data[level_number].rigidities;
    blitz::Array<blitz::TinyVector<double,NDIM>,1>&                curvatures = d_beam_data[level_number].curvatures;
    blitz::Array<const double*,1>&                         dynamic_rigidities = d_beam_data[level_number].dynamic_rigidities;
    blitz::Array<const blitz::TinyVector<double,NDIM>*,1>& dynamic_curvatures = d_beam_data[level_number].dynamic_curvatures;

    // The LMesh object provides the set of local Lagrangian nodes.
    const Pointer<LMesh> mesh = l_data_manager->getLMesh(level_number);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

    // Quick return if local_nodes is empty.
    if (local_nodes.empty())
    {
        static const int num_beams = 0;
        petsc_mastr_node_idxs. resize(num_beams);
        petsc_next_node_idxs  .resize(num_beams);
        petsc_prev_node_idxs  .resize(num_beams);
        if (d_constant_material_properties)
        {
            rigidities        .resize(num_beams);
            curvatures        .resize(num_beams);
        }
        else
        {
            dynamic_rigidities.resize(num_beams);
            dynamic_curvatures.resize(num_beams);
        }
    }

    // Determine how many beams are associated with the present MPI process.
    int num_beams = 0;
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
        const IBBeamForceSpec* const force_spec = node_idx->getNodeDataItem<IBBeamForceSpec>();
        if (force_spec != NULL) num_beams += force_spec->getNumberOfBeams();
    }
    petsc_mastr_node_idxs. resize(num_beams);
    petsc_next_node_idxs  .resize(num_beams);
    petsc_prev_node_idxs  .resize(num_beams);
    if (d_constant_material_properties)
    {
        rigidities        .resize(num_beams);
        curvatures        .resize(num_beams);
    }
    else
    {
        dynamic_rigidities.resize(num_beams);
        dynamic_curvatures.resize(num_beams);
    }

    // Return early if there are no local beams.
    if (num_beams == 0) return;

    // Setup the data structures used to compute beam forces.
    int current_beam = 0;
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
        const IBBeamForceSpec* const force_spec = node_idx->getNodeDataItem<IBBeamForceSpec>();
        if (force_spec == NULL) continue;

#ifdef DEBUG_CHECK_ASSERTIONS
        const int lag_idx = node_idx->getLagrangianIndex();
        TBOX_ASSERT(lag_idx == force_spec->getMasterNodeIndex());
#endif
        const int petsc_idx = node_idx->getGlobalPETScIndex();
        const std::vector<std::pair<int,int> >& nghbrs = force_spec->getNeighborNodeIndices();
        const std::vector<double>& bend = force_spec->getBendingRigidities();
        const std::vector<blitz::TinyVector<double,NDIM> >& curv = force_spec->getMeshDependentCurvatures();
        const unsigned int num_beams = force_spec->getNumberOfBeams();
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(num_beams == nghbrs.size());
        TBOX_ASSERT(num_beams == bend.size());
        TBOX_ASSERT(num_beams == curv.size());
#endif
        for (unsigned int k = 0; k < num_beams; ++k)
        {
            petsc_mastr_node_idxs (current_beam) = petsc_idx;
            petsc_next_node_idxs  (current_beam) = nghbrs[k].first;
            petsc_prev_node_idxs  (current_beam) = nghbrs[k].second;
            if (d_constant_material_properties)
            {
                rigidities        (current_beam) = bend[k];
                curvatures        (current_beam) = curv[k];
            }
            else
            {
                dynamic_rigidities(current_beam) = &bend[k];
                dynamic_curvatures(current_beam) = &curv[k];
            }
            ++current_beam;
        }
    }

    // Map the Lagrangian neighbor node indices to the PETSc indices
    // corresponding to the present data distribution.
    l_data_manager->mapLagrangianToPETSc(petsc_next_node_idxs, level_number);
    l_data_manager->mapLagrangianToPETSc(petsc_prev_node_idxs, level_number);

    // Determine the global node offset and the number of local nodes.
    const int global_node_offset = l_data_manager->getGlobalNodeOffset(level_number);
    const int num_local_nodes = l_data_manager->getNumberOfLocalNodes(level_number);

    // Determine the ghost nodes required to compute beam forces.
    //
    // NOTE: Only neighbor nodes can be "off processor".  Master nodes are
    // guaranteed to be "on processor".
    for (blitz::Array<int,1>::const_iterator cit = petsc_next_node_idxs.begin();
         cit != petsc_next_node_idxs.end(); ++cit)
    {
        const int idx = *cit;
        if (idx < global_node_offset || idx >= global_node_offset+num_local_nodes)
        {
            nonlocal_petsc_idx_set.insert(idx);
        }
    }
    for (blitz::Array<int,1>::const_iterator cit = petsc_prev_node_idxs.begin();
         cit != petsc_prev_node_idxs.end(); ++cit)
    {
        const int idx = *cit;
        if (idx < global_node_offset || idx >= global_node_offset+num_local_nodes)
        {
            nonlocal_petsc_idx_set.insert(idx);
        }
    }
    return;
}// initializeBeamLevelData

void
IBStandardForceGen::computeLagrangianBeamForce(
    Pointer<LData> F_data,
    Pointer<LData> X_data,
    const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    const int level_number,
    const double /*data_time*/,
    LDataManager* const /*l_data_manager*/)
{
    const int num_beams = d_beam_data[level_number].petsc_mastr_node_idxs.size();
    const int*                             const restrict petsc_mastr_node_idxs = d_beam_data[level_number].petsc_mastr_node_idxs.data();
    const int*                             const restrict  petsc_next_node_idxs = d_beam_data[level_number].petsc_next_node_idxs .data();
    const int*                             const restrict  petsc_prev_node_idxs = d_beam_data[level_number].petsc_prev_node_idxs .data();
    const double*                          const restrict            rigidities = d_beam_data[level_number].rigidities           .data();
    const blitz::TinyVector<double,NDIM>*  const restrict            curvatures = d_beam_data[level_number].curvatures           .data();
    const double**                         const restrict    dynamic_rigidities = d_beam_data[level_number].dynamic_rigidities   .data();
    const blitz::TinyVector<double,NDIM>** const restrict    dynamic_curvatures = d_beam_data[level_number].dynamic_curvatures   .data();
    double*                                const restrict                F_node = F_data->getGhostedLocalFormVecArray()->data();
    const double*                          const restrict                X_node = X_data->getGhostedLocalFormVecArray()->data();

    static const int BLOCKSIZE = 16;  // This parameter needs to be tuned.
    int k, kblock, kunroll, mastr_idx, next_idx, prev_idx;
    double K;
    const double* restrict D2X0;
    double F[NDIM];
    kblock = 0;
    if (d_constant_material_properties)
    {
        for ( ; kblock < (num_beams-1)/BLOCKSIZE; ++kblock)  // ensure that the last block is NOT handled by this first loop
        {
            PREFETCH_READ_NTA_BLOCK(petsc_mastr_node_idxs+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK( petsc_next_node_idxs+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK( petsc_prev_node_idxs+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK(           rigidities+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK(           curvatures+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            for (kunroll = 0; kunroll < BLOCKSIZE; ++kunroll)
            {
                k = kblock*BLOCKSIZE+kunroll;
                mastr_idx = petsc_mastr_node_idxs[k];
                next_idx  = petsc_next_node_idxs [k];
                prev_idx  = petsc_prev_node_idxs [k];
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(mastr_idx != next_idx);
                TBOX_ASSERT(mastr_idx != prev_idx);
#endif
                PREFETCH_READ_NTA_NDIM_BLOCK(F_node+petsc_mastr_node_idxs[k+1]);
                PREFETCH_READ_NTA_NDIM_BLOCK(F_node+ petsc_next_node_idxs[k+1]);
                PREFETCH_READ_NTA_NDIM_BLOCK(F_node+ petsc_prev_node_idxs[k+1]);
                PREFETCH_READ_NTA_NDIM_BLOCK(X_node+petsc_mastr_node_idxs[k+1]);
                PREFETCH_READ_NTA_NDIM_BLOCK(X_node+ petsc_next_node_idxs[k+1]);
                PREFETCH_READ_NTA_NDIM_BLOCK(X_node+ petsc_prev_node_idxs[k+1]);
                K = rigidities[k];
                D2X0 = curvatures[k].data();
                F[0] = K*(X_node[next_idx+0]+X_node[prev_idx+0]-2.0*X_node[mastr_idx+0]-D2X0[0]);
                F[1] = K*(X_node[next_idx+1]+X_node[prev_idx+1]-2.0*X_node[mastr_idx+1]-D2X0[1]);
#if (NDIM == 3)
                F[2] = K*(X_node[next_idx+2]+X_node[prev_idx+2]-2.0*X_node[mastr_idx+2]-D2X0[2]);
#endif
                F_node[mastr_idx+0] += 2.0*F[0];
                F_node[mastr_idx+1] += 2.0*F[1];
#if (NDIM == 3)
                F_node[mastr_idx+2] += 2.0*F[2];
#endif
                F_node[next_idx +0] -=     F[0];
                F_node[next_idx +1] -=     F[1];
#if (NDIM == 3)
                F_node[next_idx +2] -=     F[2];
#endif
                F_node[prev_idx +0] -=     F[0];
                F_node[prev_idx +1] -=     F[1];
#if (NDIM == 3)
                F_node[prev_idx +2] -=     F[2];
#endif
            }
        }
        for (k = kblock*BLOCKSIZE; k < num_beams; ++k)
        {
            mastr_idx = petsc_mastr_node_idxs[k];
            next_idx  = petsc_next_node_idxs [k];
            prev_idx  = petsc_prev_node_idxs [k];
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(mastr_idx != next_idx);
            TBOX_ASSERT(mastr_idx != prev_idx);
#endif
            K = rigidities[k];
            D2X0 = curvatures[k].data();
            F[0] = K*(X_node[next_idx+0]+X_node[prev_idx+0]-2.0*X_node[mastr_idx+0]-D2X0[0]);
            F[1] = K*(X_node[next_idx+1]+X_node[prev_idx+1]-2.0*X_node[mastr_idx+1]-D2X0[1]);
#if (NDIM == 3)
            F[2] = K*(X_node[next_idx+2]+X_node[prev_idx+2]-2.0*X_node[mastr_idx+2]-D2X0[2]);
#endif
            F_node[mastr_idx+0] += 2.0*F[0];
            F_node[mastr_idx+1] += 2.0*F[1];
#if (NDIM == 3)
            F_node[mastr_idx+2] += 2.0*F[2];
#endif
            F_node[next_idx +0] -=     F[0];
            F_node[next_idx +1] -=     F[1];
#if (NDIM == 3)
            F_node[next_idx +2] -=     F[2];
#endif
            F_node[prev_idx +0] -=     F[0];
            F_node[prev_idx +1] -=     F[1];
#if (NDIM == 3)
            F_node[prev_idx +2] -=     F[2];
#endif
        }
    }
    else
    {
        for ( ; kblock < (num_beams-1)/BLOCKSIZE; ++kblock)  // ensure that the last block is NOT handled by this first loop
        {
            PREFETCH_READ_NTA_BLOCK(petsc_mastr_node_idxs+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK( petsc_next_node_idxs+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK( petsc_prev_node_idxs+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK(   dynamic_rigidities+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK(   dynamic_curvatures+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            for (kunroll = 0; kunroll < BLOCKSIZE; ++kunroll)
            {
                k = kblock*BLOCKSIZE+kunroll;
                mastr_idx = petsc_mastr_node_idxs[k];
                next_idx  = petsc_next_node_idxs [k];
                prev_idx  = petsc_prev_node_idxs [k];
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(mastr_idx != next_idx);
                TBOX_ASSERT(mastr_idx != prev_idx);
#endif
                PREFETCH_READ_NTA_NDIM_BLOCK(F_node+petsc_mastr_node_idxs[k+1]);
                PREFETCH_READ_NTA_NDIM_BLOCK(F_node+ petsc_next_node_idxs[k+1]);
                PREFETCH_READ_NTA_NDIM_BLOCK(F_node+ petsc_prev_node_idxs[k+1]);
                PREFETCH_READ_NTA_NDIM_BLOCK(X_node+petsc_mastr_node_idxs[k+1]);
                PREFETCH_READ_NTA_NDIM_BLOCK(X_node+ petsc_next_node_idxs[k+1]);
                PREFETCH_READ_NTA_NDIM_BLOCK(X_node+ petsc_prev_node_idxs[k+1]);
                PREFETCH_READ_NTA(                     dynamic_rigidities[k+1]);
                PREFETCH_READ_NTA(                     dynamic_curvatures[k+1]);
                K = *dynamic_rigidities[k];
                D2X0 = dynamic_curvatures[k]->data();
                F[0] = K*(X_node[next_idx+0]+X_node[prev_idx+0]-2.0*X_node[mastr_idx+0]-D2X0[0]);
                F[1] = K*(X_node[next_idx+1]+X_node[prev_idx+1]-2.0*X_node[mastr_idx+1]-D2X0[1]);
#if (NDIM == 3)
                F[2] = K*(X_node[next_idx+2]+X_node[prev_idx+2]-2.0*X_node[mastr_idx+2]-D2X0[2]);
#endif
                F_node[mastr_idx+0] += 2.0*F[0];
                F_node[mastr_idx+1] += 2.0*F[1];
#if (NDIM == 3)
                F_node[mastr_idx+2] += 2.0*F[2];
#endif
                F_node[next_idx +0] -=     F[0];
                F_node[next_idx +1] -=     F[1];
#if (NDIM == 3)
                F_node[next_idx +2] -=     F[2];
#endif
                F_node[prev_idx +0] -=     F[0];
                F_node[prev_idx +1] -=     F[1];
#if (NDIM == 3)
                F_node[prev_idx +2] -=     F[2];
#endif
            }
        }
        for (k = kblock*BLOCKSIZE; k < num_beams; ++k)
        {
            mastr_idx = petsc_mastr_node_idxs[k];
            next_idx  = petsc_next_node_idxs [k];
            prev_idx  = petsc_prev_node_idxs [k];
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(mastr_idx != next_idx);
            TBOX_ASSERT(mastr_idx != prev_idx);
#endif
            K = *dynamic_rigidities[k];
            D2X0 = dynamic_curvatures[k]->data();
            F[0] = K*(X_node[next_idx+0]+X_node[prev_idx+0]-2.0*X_node[mastr_idx+0]-D2X0[0]);
            F[1] = K*(X_node[next_idx+1]+X_node[prev_idx+1]-2.0*X_node[mastr_idx+1]-D2X0[1]);
#if (NDIM == 3)
            F[2] = K*(X_node[next_idx+2]+X_node[prev_idx+2]-2.0*X_node[mastr_idx+2]-D2X0[2]);
#endif
            F_node[mastr_idx+0] += 2.0*F[0];
            F_node[mastr_idx+1] += 2.0*F[1];
#if (NDIM == 3)
            F_node[mastr_idx+2] += 2.0*F[2];
#endif
            F_node[next_idx +0] -=     F[0];
            F_node[next_idx +1] -=     F[1];
#if (NDIM == 3)
            F_node[next_idx +2] -=     F[2];
#endif
            F_node[prev_idx +0] -=     F[0];
            F_node[prev_idx +1] -=     F[1];
#if (NDIM == 3)
            F_node[prev_idx +2] -=     F[2];
#endif
        }
    }

    F_data->restoreArrays();
    X_data->restoreArrays();
    return;
}// computeLagrangianBeamForce

void
IBStandardForceGen::initializeTargetPointLevelData(
    std::set<int>& /*nonlocal_petsc_idx_set*/,
    const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    const int level_number,
    const double /*init_data_time*/,
    const bool /*initial_time*/,
    LDataManager* const l_data_manager)
{
    blitz::Array<int,1>&                              petsc_node_idxs = d_target_point_data[level_number].petsc_node_idxs;
    blitz::Array<double,1>&                                     kappa = d_target_point_data[level_number].kappa;
    blitz::Array<double,1>&                                       eta = d_target_point_data[level_number].eta;
    blitz::Array<blitz::TinyVector<double,NDIM>,1>&                X0 = d_target_point_data[level_number].X0;
    blitz::Array<const double*,1>&                      dynamic_kappa = d_target_point_data[level_number].dynamic_kappa;
    blitz::Array<const double*,1>&                        dynamic_eta = d_target_point_data[level_number].dynamic_eta;
    blitz::Array<const blitz::TinyVector<double,NDIM>*,1>& dynamic_X0 = d_target_point_data[level_number].dynamic_X0;

    // The LMesh object provides the set of local Lagrangian nodes.
    const Pointer<LMesh> mesh = l_data_manager->getLMesh(level_number);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

    // Quick return if local_nodes is empty.
    if (local_nodes.empty())
    {
        static const int num_target_points = 0;
        petsc_node_idxs.resize(num_target_points);
        if (d_constant_material_properties)
        {
            kappa        .resize(num_target_points);
            eta          .resize(num_target_points);
            X0           .resize(num_target_points);
            dynamic_kappa.resize(0);
            dynamic_eta  .resize(0);
            dynamic_X0   .resize(0);
        }
        else
        {
            kappa        .resize(0);
            eta          .resize(0);
            X0           .resize(0);
            dynamic_kappa.resize(num_target_points);
            dynamic_eta  .resize(num_target_points);
            dynamic_X0   .resize(num_target_points);
        }
        return;
    }

    // Determine how many target points are associated with the present MPI
    // process.
    int num_target_points = 0;
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
        const IBTargetPointForceSpec* const force_spec = node_idx->getNodeDataItem<IBTargetPointForceSpec>();
        if (force_spec != NULL) num_target_points += 1;
    }

    // Resize arrays for storing cached values used to compute target point
    // forces.
    petsc_node_idxs.resize(num_target_points);
    if (d_constant_material_properties)
    {
        kappa        .resize(num_target_points);
        eta          .resize(num_target_points);
        X0           .resize(num_target_points);
        dynamic_kappa.resize(0);
        dynamic_eta  .resize(0);
        dynamic_X0   .resize(0);
    }
    else
    {
        kappa        .resize(0);
        eta          .resize(0);
        X0           .resize(0);
        dynamic_kappa.resize(num_target_points);
        dynamic_eta  .resize(num_target_points);
        dynamic_X0   .resize(num_target_points);
    }

    // Return early if there are no local target points.
    if (num_target_points == 0) return;

    // Setup the data structures used to compute target point forces.
    int current_target_point = 0;
    if (d_constant_material_properties)
    {
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const IBTargetPointForceSpec* const force_spec = node_idx->getNodeDataItem<IBTargetPointForceSpec>();
            if (force_spec == NULL) continue;
            petsc_node_idxs(current_target_point) = node_idx->getGlobalPETScIndex();
            kappa          (current_target_point) = force_spec->getStiffness();
            eta            (current_target_point) = force_spec->getDamping();
            X0             (current_target_point) = force_spec->getTargetPointPosition();
            ++current_target_point;
        }
    }
    else
    {
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const IBTargetPointForceSpec* const force_spec = node_idx->getNodeDataItem<IBTargetPointForceSpec>();
            if (force_spec == NULL) continue;
            petsc_node_idxs(current_target_point) = node_idx->getGlobalPETScIndex();
            dynamic_kappa  (current_target_point) = &force_spec->getStiffness();
            dynamic_eta    (current_target_point) = &force_spec->getDamping();
            dynamic_X0     (current_target_point) = &force_spec->getTargetPointPosition();
            ++current_target_point;
        }
    }
    return;
}// initializeTargetPointLevelData

void
IBStandardForceGen::computeLagrangianTargetPointForce(
    Pointer<LData> F_data,
    Pointer<LData> X_data,
    Pointer<LData> U_data,
    const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    const int level_number,
    const double /*data_time*/,
    LDataManager* const /*l_data_manager*/)
{
    const int num_target_points = d_target_point_data[level_number].petsc_node_idxs.size();
    const int*                             const restrict petsc_node_idxs = d_target_point_data[level_number].petsc_node_idxs.data();
    const double*                          const restrict           kappa = d_target_point_data[level_number].kappa          .data();
    const double*                          const restrict             eta = d_target_point_data[level_number].eta            .data();
    const blitz::TinyVector<double,NDIM>*  const restrict              X0 = d_target_point_data[level_number].X0             .data();
    const double**                         const restrict   dynamic_kappa = d_target_point_data[level_number].dynamic_kappa  .data();
    const double**                         const restrict     dynamic_eta = d_target_point_data[level_number].dynamic_eta    .data();
    const blitz::TinyVector<double,NDIM>** const restrict      dynamic_X0 = d_target_point_data[level_number].dynamic_X0     .data();
    double*                                const restrict          F_node = F_data->getGhostedLocalFormVecArray()->data();
    const double*                          const restrict          X_node = X_data->getGhostedLocalFormVecArray()->data();
    const double*                          const restrict          U_node = U_data->getGhostedLocalFormVecArray()->data();

    static const int BLOCKSIZE = 16;  // This parameter needs to be tuned.
    int k, kblock, kunroll, idx;
    double K, E;
    const double* restrict X_target;
    kblock = 0;
    if (d_constant_material_properties)
    {
        for ( ; kblock < (num_target_points-1)/BLOCKSIZE; ++kblock)  // ensure that the last block is NOT handled by this first loop
        {
            PREFETCH_READ_NTA_BLOCK(petsc_node_idxs+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK(          kappa+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK(            eta+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK(             X0+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            for (kunroll = 0; kunroll < BLOCKSIZE; ++kunroll)
            {
                k = kblock*BLOCKSIZE+kunroll;
                idx = petsc_node_idxs[k];
                PREFETCH_READ_NTA_NDIM_BLOCK(F_node+petsc_node_idxs[k+1]);
                PREFETCH_READ_NTA_NDIM_BLOCK(X_node+petsc_node_idxs[k+1]);
                K = kappa[k];
                E = eta[k];
                X_target = X0[k].data();
                F_node[idx+0] += K*(X_target[0] - X_node[idx+0]) - E*U_node[idx+0];
                F_node[idx+1] += K*(X_target[1] - X_node[idx+1]) - E*U_node[idx+1];
#if (NDIM == 3)
                F_node[idx+2] += K*(X_target[2] - X_node[idx+2]) - E*U_node[idx+2];
#endif
            }
        }
        for (k = kblock*BLOCKSIZE; k < num_target_points; ++k)
        {
            idx = petsc_node_idxs[k];
            K = kappa[k];
            E = eta[k];
            X_target = X0[k].data();
            F_node[idx+0] += K*(X_target[0] - X_node[idx+0]) - E*U_node[idx+0];
            F_node[idx+1] += K*(X_target[1] - X_node[idx+1]) - E*U_node[idx+1];
#if (NDIM == 3)
            F_node[idx+2] += K*(X_target[2] - X_node[idx+2]) - E*U_node[idx+2];
#endif
        }
    }
    else
    {
        for ( ; kblock < (num_target_points-1)/BLOCKSIZE; ++kblock)  // ensure that the last block is NOT handled by this first loop
        {
            PREFETCH_READ_NTA_BLOCK(petsc_node_idxs+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK(  dynamic_kappa+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK(    dynamic_eta+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            PREFETCH_READ_NTA_BLOCK(     dynamic_X0+BLOCKSIZE*(kblock+1), BLOCKSIZE);
            for (kunroll = 0; kunroll < BLOCKSIZE; ++kunroll)
            {
                k = kblock*BLOCKSIZE+kunroll;
                idx = petsc_node_idxs[k];
                PREFETCH_READ_NTA_NDIM_BLOCK(F_node+petsc_node_idxs[k+1]);
                PREFETCH_READ_NTA_NDIM_BLOCK(X_node+petsc_node_idxs[k+1]);
                PREFETCH_READ_NTA(                    dynamic_kappa[k+1]);
                PREFETCH_READ_NTA(                      dynamic_eta[k+1]);
                PREFETCH_READ_NTA(                       dynamic_X0[k+1]);
                K = *dynamic_kappa[k];
                E = *dynamic_eta[k];
                X_target = dynamic_X0[k]->data();
                F_node[idx+0] += K*(X_target[0] - X_node[idx+0]) - E*U_node[idx+0];
                F_node[idx+1] += K*(X_target[1] - X_node[idx+1]) - E*U_node[idx+1];
#if (NDIM == 3)
                F_node[idx+2] += K*(X_target[2] - X_node[idx+2]) - E*U_node[idx+2];
#endif
            }
        }
        for (k = kblock*BLOCKSIZE; k < num_target_points; ++k)
        {
            idx = petsc_node_idxs[k];
            K = *dynamic_kappa[k];
            E = *dynamic_eta[k];
            X_target = dynamic_X0[k]->data();
            F_node[idx+0] += K*(X_target[0] - X_node[idx+0]) - E*U_node[idx+0];
            F_node[idx+1] += K*(X_target[1] - X_node[idx+1]) - E*U_node[idx+1];
#if (NDIM == 3)
            F_node[idx+2] += K*(X_target[2] - X_node[idx+2]) - E*U_node[idx+2];
#endif
        }
    }

    F_data->restoreArrays();
    X_data->restoreArrays();
    U_data->restoreArrays();
    return;
}// computeLagrangianTargetPointForce

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
