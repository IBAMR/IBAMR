// Filename: IBPDForceGen.cpp
// Created on 08 Apr 2016 by Amneet Bhalla
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

#include <math.h>
#include <stddef.h>
#include <algorithm>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <ostream>
#include <set>
#include <utility>

#include "IntVector.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "boost/multi_array.hpp"
#include "ibamr/IBPDForceGen.h"
#include "ibamr/IBSpringForceSpec.h"
#include "ibamr/IBTargetPointForceSpec.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LMesh.h"
#include "ibtk/LNode.h"
#include "ibtk/compiler_hints.h"
#include "ibtk/ibtk_utilities.h"
#include "petscsys.h"
#include "petscvec.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

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
    for (std::vector<int>::iterator it = inds.begin(); it != inds.end(); ++it)
    {
        int& idx = *it;
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
    for (std::vector<int>::iterator it = inds.begin(); it != inds.end(); ++it)
    {
        int& idx = *it;
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

double
default_inf_fcn(double R0, double delta)
{
#if (NDIM == 2)
    static const double A = 2;
    static const double C = 60.0 / (7.0 * M_PI * A * A);
#elif(NDIM == 3)
    static const double A = 2;
    static const double C = 24.0 / (2.0 * M_PI * A * A * A);
#endif

    double W;

    double r = R0 / delta;
    if (r < 1)
    {
        W = C * (2.0 / 3.0 - r * r + 0.5 * r * r * r);
    }
    else if (r <= 2)
    {
        W = C * std::pow((2.0 - r), 3) / 6.0;
    }
    else
    {
        W = 0.0;
    }

    return W;

} // default_inf_fcn

double
default_vol_frac_fcn(double R0, double horizon, double delta)
{
    double vol_frac;
    if (R0 <= (horizon - delta))
    {
        vol_frac = 1.0;
    }
    else if (R0 <= (horizon + delta))
    {
        vol_frac = (horizon + delta - R0) / (2.0 * delta);
    }
    else
    {
        vol_frac = 0.0;
    }

    return vol_frac;

} // default_vol_frac_fcn

void
default_PK1_fcn(Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor>& PK1,
                const Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> >& FF,
                const Eigen::Map<const IBTK::Vector>& /*X0*/,
                int /*lag_idx*/)
{
    typedef Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> mat_type;
    static const mat_type II = mat_type::Identity();
    mat_type E = 0.5 * (FF.transpose() * FF - II);
    const double trE = E.trace();
    PK1 = IBPDForceGen::s_lame_0 * trE * FF + 2 * IBPDForceGen::s_lame_1 * FF * E;

    return;

} // default_PK1_fcn

Eigen::Vector4d
default_force_damage_fcn(const double /*horizon*/,
                         const double /*delta*/,
                         const double W,
                         const double vol_frac,
                         double* parameters,
                         const Eigen::Map<const IBTK::Vector>& X0_mastr,
                         const Eigen::Map<const IBTK::Vector>& X0_slave,
                         const Eigen::Map<const IBTK::Vector>& X_mastr,
                         const Eigen::Map<const IBTK::Vector>& X_slave,
                         const Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> >& FF_mastr,
                         const Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> >& FF_slave,
                         const Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> >& B_mastr,
                         const Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> >& B_slave,
                         Eigen::Map<IBTK::Vector>& F_mastr,
                         Eigen::Map<IBTK::Vector>& F_slave,
                         const int lag_mastr_node_idx,
                         const int lag_slave_node_idx)
{
    // Bond parameters
    // 0 --> Kappa, 1 --> R0, 2--> user defined
    const double& R0 = parameters[1];
    const double& vol_mastr = parameters[2];
    const double& vol_slave = parameters[3];
    double& fail = parameters[4];
    const double& critical_stretch = parameters[5];

    // Estimate failure.
    const double R = (X_slave - X_mastr).norm();
    const double stretch = (R - R0) / R0;
    if (!MathUtilities<double>::equalEps(fail, 0.0) && std::fabs(stretch) > critical_stretch)
    {
        fail = 0.0;
    }

    // PK1 stress tensor
    typedef IBTK::Vector vec_type;
    typedef Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> mat_type;
    mat_type PK1_mastr, PK1_slave;
    default_PK1_fcn(PK1_mastr, FF_mastr, X0_mastr, lag_mastr_node_idx);
    default_PK1_fcn(PK1_slave, FF_slave, X0_slave, lag_slave_node_idx);

    // Compute PD force.
    vec_type trac = W * (PK1_mastr * B_mastr + PK1_slave * B_slave) * (X0_slave - X0_mastr);
    F_mastr += fail * (vol_frac * vol_slave) * trac * (vol_frac * vol_slave);
    F_slave += -fail * (vol_frac * vol_mastr) * trac * (vol_frac * vol_mastr);

    // Compute damage.
    Eigen::Vector4d D;
    D(0) = vol_slave * vol_frac * fail;
    D(1) = vol_slave * vol_frac;
    D(2) = vol_mastr * vol_frac * fail;
    D(3) = vol_mastr * vol_frac;

    return D;

} // default_force_damage_fcn

void
default_target_point_force_fcn(const Eigen::Map<const IBTK::Vector>& X,
                               const Eigen::Map<const IBTK::Vector>& X_target,
                               const Eigen::Map<const IBTK::Vector>& U,
                               double K,
                               double E,
                               int /*lag_idx*/,
                               Eigen::Map<IBTK::Vector>& F)
{
    F += K * (X_target - X) - E * U;

    return;
} // default_target_point_force_fcn
}

double IBPDForceGen::s_lame_0 = 1.0;
double IBPDForceGen::s_lame_1 = 1.0;

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBPDForceGen::IBPDForceGen(Pointer<Database> input_db) : d_horizon(3.0), d_ds(1.0), d_use_mean_disp(false)
{
    // Get values from input database.
    if (input_db)
    {
        if (input_db->keyExists("horizon"))
        {
            d_horizon = input_db->getDouble("horizon");
        }

        if (input_db->keyExists("ds"))
        {
            d_ds = input_db->getDouble("ds");
        }

        if (input_db->keyExists("use_mean_displacement"))
        {
            d_use_mean_disp = input_db->getBool("use_mean_displacement");
        }
    }
    registerBondForceSpecificationFunction(
        0, &default_PK1_fcn, &default_force_damage_fcn, &default_inf_fcn, &default_vol_frac_fcn);
    registerTargetPointForceFunction(&default_target_point_force_fcn);

    return;
} // IBStandardForceGen

IBPDForceGen::~IBPDForceGen()
{
    // intentionally blank
    return;
} // ~IBPDForceGen

void
IBPDForceGen::registerBondForceSpecificationFunction(int force_fcn_index,
                                                     const BondPK1FcnPtr bond_PK1_fcn_ptr,
                                                     const BondForceDamageFcnPtr bond_force_damage_fcn_ptr,
                                                     const BondInfluenceFcnPtr bond_inf_fcn_ptr,
                                                     const BondVolFracFcnPtr bond_vol_frac_fcn_ptr)
{
    d_bond_PK1_fcn_map[force_fcn_index] = bond_PK1_fcn_ptr;
    d_bond_force_damage_fcn_map[force_fcn_index] = bond_force_damage_fcn_ptr;
    d_bond_inf_fcn_map[force_fcn_index] = bond_inf_fcn_ptr;
    d_bond_vol_frac_fcn_map[force_fcn_index] = bond_vol_frac_fcn_ptr;

    return;

} // registerBondForceSpecificationFunction

void
IBPDForceGen::registerTargetPointForceFunction(const TargetPointForceFcnPtr target_point_force_fcn_ptr)
{
    d_target_point_force_fcn = target_point_force_fcn_ptr;
    return;

} // registerTargetPointForceFunction

void
IBPDForceGen::initializeLevelData(const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                  const int level_number,
                                  const double init_data_time,
                                  const bool initial_time,
                                  LDataManager* const l_data_manager)
{
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return;
    plog << "IBPDForceGen::initializeLevelData(): initializing IBPDForceGen level data." << std::endl;

#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
#endif
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Resize the vectors corresponding to data individually maintained for
    // separate levels of the patch hierarchy.
    const int new_size = std::max(level_number + 1, static_cast<int>(d_is_initialized.size()));

    d_bond_data.resize(new_size);
    d_target_point_data.resize(new_size);
    d_X0_ghost_data.resize(new_size, Pointer<LData>(NULL));
    d_X_ghost_data.resize(new_size, Pointer<LData>(NULL));
    d_X_mean_ghost_data.resize(new_size, Pointer<LData>(NULL));
    d_F_ghost_data.resize(new_size, Pointer<LData>(NULL));
    d_N_ghost_data.resize(new_size, Pointer<LData>(NULL));
    d_dmg_ghost_data.resize(new_size, Pointer<LData>(NULL));
    d_B_ghost_data.resize(new_size, Pointer<LData>(NULL));
    d_FF_ghost_data.resize(new_size, Pointer<LData>(NULL));
    d_dX_data.resize(new_size, Pointer<LData>(NULL));
    d_is_initialized.resize(new_size, false);

    // Keep track of all of the nonlocal PETSc indices required to compute the
    // forces.
    std::set<int> nonlocal_petsc_idx_set;

    // Setup the cached data.
    initializeBondLevelData(
        nonlocal_petsc_idx_set, hierarchy, level_number, init_data_time, initial_time, l_data_manager);
    initializeTargetPointLevelData(
        nonlocal_petsc_idx_set, hierarchy, level_number, init_data_time, initial_time, l_data_manager);

    // Put the nonlocal PETSc indices into a vector.
    std::vector<int> nonlocal_petsc_idxs(nonlocal_petsc_idx_set.begin(), nonlocal_petsc_idx_set.end());

    // Put all cached PETSc node indices into local form.
    const int global_node_offset = l_data_manager->getGlobalNodeOffset(level_number);
    const int num_local_nodes = l_data_manager->getNumberOfLocalNodes(level_number);
    resetLocalPETScIndices(d_bond_data[level_number].petsc_mastr_node_idxs, global_node_offset, num_local_nodes);
    resetLocalOrNonlocalPETScIndices(
        d_bond_data[level_number].petsc_slave_node_idxs, global_node_offset, num_local_nodes, nonlocal_petsc_idxs);
    resetLocalPETScIndices(d_target_point_data[level_number].petsc_node_idxs, global_node_offset, num_local_nodes);

    std::ostringstream X0_name_stream;
    X0_name_stream << "IBPDForceGen::X0_ghost_" << level_number;
    d_X0_ghost_data[level_number] = new LData(X0_name_stream.str(), num_local_nodes, NDIM, nonlocal_petsc_idxs);

    std::ostringstream X_name_stream;
    X_name_stream << "IBPDForceGen::X_ghost_" << level_number;
    d_X_ghost_data[level_number] = new LData(X_name_stream.str(), num_local_nodes, NDIM, nonlocal_petsc_idxs);

    if (d_use_mean_disp)
    {
        std::ostringstream X_mean_name_stream;
        X_mean_name_stream << "IBPDForceGen::X_mean_ghost_" << level_number;
        d_X_mean_ghost_data[level_number] =
            new LData(X_mean_name_stream.str(), num_local_nodes, NDIM, nonlocal_petsc_idxs);

        std::ostringstream N_name_stream;
        N_name_stream << "IBPDForceGen::N_ghost_" << level_number;
        d_N_ghost_data[level_number] = new LData(N_name_stream.str(), num_local_nodes, 1, nonlocal_petsc_idxs);
    }

    std::ostringstream F_name_stream;
    F_name_stream << "IBPDForceGen::F_ghost_" << level_number;
    d_F_ghost_data[level_number] = new LData(F_name_stream.str(), num_local_nodes, NDIM, nonlocal_petsc_idxs);

    std::ostringstream dmg_name_stream;
    dmg_name_stream << "IBPDForceGen::dmg_ghost_" << level_number;
    d_dmg_ghost_data[level_number] = new LData(dmg_name_stream.str(), num_local_nodes, 2, nonlocal_petsc_idxs);

    std::ostringstream B_name_stream;
    B_name_stream << "IBPDForceGen::B_ghost_" << level_number;
    d_B_ghost_data[level_number] = new LData(B_name_stream.str(), num_local_nodes, NDIM * NDIM, nonlocal_petsc_idxs);

    std::ostringstream FF_name_stream;
    FF_name_stream << "IBPDForceGen::FF_ghost_" << level_number;
    d_FF_ghost_data[level_number] = new LData(FF_name_stream.str(), num_local_nodes, NDIM * NDIM, nonlocal_petsc_idxs);

    std::ostringstream dX_name_stream;
    dX_name_stream << "IBPDForceGen::dX_" << level_number;
    d_dX_data[level_number] = new LData(dX_name_stream.str(), num_local_nodes, NDIM);

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

    // Indicate that the level data has been initialized and we need to compute shape tensor.
    d_is_initialized[level_number] = true;
    d_compute_shape_tensor = true;

    return;
} // initializeLevelData

void
IBPDForceGen::computeLagrangianForceAndDamage(Pointer<LData> F_data,
                                              Pointer<LData> D_data,
                                              Pointer<LData> X_data,
                                              Pointer<LData> U_data,
                                              const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                              const int level_number,
                                              const double data_time,
                                              LDataManager* const l_data_manager)
{
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return;
    const Pointer<LMesh> mesh = l_data_manager->getLMesh(level_number);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
    int ierr;

    // Initialize various ghost data.
    Pointer<LData> F_ghost_data = d_F_ghost_data[level_number];
    Vec F_ghost_local_form_vec;
    ierr = VecGhostGetLocalForm(F_ghost_data->getVec(), &F_ghost_local_form_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecSet(F_ghost_local_form_vec, 0.0);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostRestoreLocalForm(F_ghost_data->getVec(), &F_ghost_local_form_vec);
    IBTK_CHKERRQ(ierr);

    Pointer<LData> D_ghost_data = d_dmg_ghost_data[level_number];
    Vec D_ghost_local_form_vec;
    ierr = VecGhostGetLocalForm(D_ghost_data->getVec(), &D_ghost_local_form_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecSet(D_ghost_local_form_vec, 0.0);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostRestoreLocalForm(D_ghost_data->getVec(), &D_ghost_local_form_vec);
    IBTK_CHKERRQ(ierr);

    Pointer<LData> X_ghost_data = d_X_ghost_data[level_number];
    Pointer<LData> dX_data = d_dX_data[level_number];
    ierr = VecAXPBYPCZ(X_ghost_data->getVec(), 1.0, 1.0, 0.0, X_data->getVec(), dX_data->getVec());
    TBOX_ASSERT(X_data);
    ierr = VecCopy(X_data->getVec(), X_ghost_data->getVec());
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(X_ghost_data->getVec(), INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(X_ghost_data->getVec(), INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);

    Pointer<LData> X_mean_ghost_data = d_X_mean_ghost_data[level_number];
    Pointer<LData> N_ghost_data = d_N_ghost_data[level_number];
    if (d_use_mean_disp)
    {
        Vec X_mean_ghost_local_form_vec;
        ierr = VecGhostGetLocalForm(X_mean_ghost_data->getVec(), &X_mean_ghost_local_form_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecSet(X_mean_ghost_local_form_vec, 0.0);
        IBTK_CHKERRQ(ierr);
        ierr = VecGhostRestoreLocalForm(X_mean_ghost_data->getVec(), &X_mean_ghost_local_form_vec);
        IBTK_CHKERRQ(ierr);

        Vec N_ghost_local_form_vec;
        ierr = VecGhostGetLocalForm(N_ghost_data->getVec(), &N_ghost_local_form_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecSet(N_ghost_local_form_vec, 0.0);
        IBTK_CHKERRQ(ierr);
        ierr = VecGhostRestoreLocalForm(N_ghost_data->getVec(), &N_ghost_local_form_vec);
        IBTK_CHKERRQ(ierr);
    }

    Pointer<LData> X0_ghost_data = d_X0_ghost_data[level_number];
    Pointer<LData> X0 = l_data_manager->getLData("X0", level_number);
    ierr = VecCopy(X0->getVec(), X0_ghost_data->getVec());
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(X0_ghost_data->getVec(), INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(X0_ghost_data->getVec(), INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);

    Pointer<LData> B_ghost_data = d_B_ghost_data[level_number];
    if (d_compute_shape_tensor)
    {
        Vec B_ghost_local_form_vec;
        ierr = VecGhostGetLocalForm(B_ghost_data->getVec(), &B_ghost_local_form_vec);
        IBTK_CHKERRQ(ierr);
        ierr = VecSet(B_ghost_local_form_vec, 0.0);
        IBTK_CHKERRQ(ierr);
        ierr = VecGhostRestoreLocalForm(B_ghost_data->getVec(), &B_ghost_local_form_vec);
        IBTK_CHKERRQ(ierr);
    }

    Pointer<LData> FF_ghost_data = d_FF_ghost_data[level_number];
    Vec FF_ghost_local_form_vec;
    ierr = VecGhostGetLocalForm(FF_ghost_data->getVec(), &FF_ghost_local_form_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecSet(FF_ghost_local_form_vec, 0.0);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostRestoreLocalForm(FF_ghost_data->getVec(), &FF_ghost_local_form_vec);
    IBTK_CHKERRQ(ierr);

    // Compute mean position
    if (d_use_mean_disp)
    {
        computeMeanPosition(
            X_mean_ghost_data, N_ghost_data, X_ghost_data, hierarchy, level_number, data_time, l_data_manager);
        ierr = VecGhostUpdateBegin(X_mean_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);
        IBTK_CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(X_mean_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);
        IBTK_CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(N_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);
        IBTK_CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(N_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);
        IBTK_CHKERRQ(ierr);

        boost::multi_array_ref<double, 1>& N_ghost_data_array = *N_ghost_data->getLocalFormArray();
        boost::multi_array_ref<double, 2>& X_mean_ghost_data_array = *X_mean_ghost_data->getLocalFormVecArray();
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int local_idx = node_idx->getLocalPETScIndex();
            double* X_mean = &X_mean_ghost_data_array[local_idx][0];
            double* N = &N_ghost_data_array[local_idx];

            for (int d = 0; d < NDIM; ++d)
            {
                X_mean[d] /= N[0];
            }
        }
        ierr = VecGhostUpdateBegin(X_mean_ghost_data->getVec(), INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(X_mean_ghost_data->getVec(), INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
        X_mean_ghost_data->restoreArrays();
        N_ghost_data->restoreArrays();
    }

    // Compute shape tensor.
    if (d_compute_shape_tensor)
    {
        computeShapeTensor(B_ghost_data, X0_ghost_data, hierarchy, level_number, data_time, l_data_manager);
        ierr = VecGhostUpdateBegin(B_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);
        IBTK_CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(B_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);
        IBTK_CHKERRQ(ierr);

        boost::multi_array_ref<double, 2>& B_ghost_data_array = *B_ghost_data->getLocalFormVecArray();
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int lag_idx = node_idx->getLagrangianIndex();
            const int local_idx = node_idx->getLocalPETScIndex();
            double* B = &B_ghost_data_array[local_idx][0];
            Eigen::Map<Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> > eig_B(B);

            // Scale the matrix.
            const double scale = eig_B.norm();
            eig_B *= (1.0 / scale);

            // Invert the scaled matrix.
            bool invertible;
            eig_B.computeInverseWithCheck(eig_B, invertible);
            if (!invertible)
            {
                TBOX_ERROR("IBPDForceGen::computeLagrangianForceAndDamage() : Matrix inversion failed.\n"
                           << " Lagrangian index = "
                           << lag_idx
                           << "\nScaled B tensor is \n"
                           << eig_B
                           << "\n"
                           << " Scaling factor  = "
                           << scale
                           << "\n");
            }

            // Scale back the inverse-matrix.
            eig_B *= (1.0 / scale);
        }

        ierr = VecGhostUpdateBegin(B_ghost_data->getVec(), INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(B_ghost_data->getVec(), INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
        B_ghost_data->restoreArrays();
        d_compute_shape_tensor = false;
    }

    // Compute the deformation gradient tensor.
    // FF = int_k {Y outer X}_k . Inv{B}
    if (d_use_mean_disp)
    {
        computeDeformationGradientTensor(
            FF_ghost_data, X_mean_ghost_data, X0_ghost_data, hierarchy, level_number, data_time, l_data_manager);
    }
    else
    {
        computeDeformationGradientTensor(
            FF_ghost_data, X_ghost_data, X0_ghost_data, hierarchy, level_number, data_time, l_data_manager);
    }
    ierr = VecGhostUpdateBegin(FF_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(FF_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);
    IBTK_CHKERRQ(ierr);

    {
        boost::multi_array_ref<double, 2>& FF_ghost_data_array = *FF_ghost_data->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& B_ghost_data_array = *B_ghost_data->getLocalFormVecArray();
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int local_idx = node_idx->getLocalPETScIndex();
            const double* B = &B_ghost_data_array[local_idx][0];
            double* FF = &FF_ghost_data_array[local_idx][0];
            Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> > eig_B(B);
            Eigen::Map<Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> > eig_FF(FF);
            eig_FF = eig_FF * eig_B;
        }
        ierr = VecGhostUpdateBegin(FF_ghost_data->getVec(), INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(FF_ghost_data->getVec(), INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
        FF_ghost_data->restoreArrays();
        B_ghost_data->restoreArrays();
    }

    // Compute the forces and damage functions.
    computeLagrangianBondForceAndDamage(F_ghost_data,
                                        D_ghost_data,
                                        X_ghost_data,
                                        X0_ghost_data,
                                        FF_ghost_data,
                                        B_ghost_data,
                                        hierarchy,
                                        level_number,
                                        data_time,
                                        l_data_manager);
    computeLagrangianTargetPointForce(
        F_ghost_data, X_ghost_data, U_data, hierarchy, level_number, data_time, l_data_manager);
    ierr = VecGhostUpdateBegin(F_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);
    IBTK_CHKERRQ(ierr);
    ierr = VecAXPY(F_data->getVec(), 1.0, F_ghost_data->getVec());
    ierr = VecGhostUpdateBegin(D_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);
    IBTK_CHKERRQ(ierr);

    // Compute the damage factor.
    boost::multi_array_ref<double, 1>& D_data_array = *D_data->getLocalFormArray();
    boost::multi_array_ref<double, 2>& D_ghost_data_array = *D_ghost_data->getLocalFormVecArray();
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
        const int local_idx = node_idx->getLocalPETScIndex();
        double* D_ghost = &D_ghost_data_array[local_idx][0];
        double* D = &D_data_array[local_idx];

        D[0] = 1.0 - D_ghost[0] / D_ghost[1];
    }
    D_data->restoreArrays();
    D_ghost_data->restoreArrays();

    return;
} // computeLagrangianForce

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBPDForceGen::computeMeanPosition(Pointer<LData> X_mean_data,
                                  Pointer<LData> N_data,
                                  Pointer<LData> X_data,
                                  Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                  int level_number,
                                  double /*data_time*/,
                                  IBTK::LDataManager* /*l_data_manager*/)
{
    const int num_bonds = static_cast<int>(d_bond_data[level_number].lag_mastr_node_idxs.size());
    if (num_bonds == 0) return;

    const int* const petsc_mastr_node_idxs = &d_bond_data[level_number].petsc_mastr_node_idxs[0];
    const int* const petsc_slave_node_idxs = &d_bond_data[level_number].petsc_slave_node_idxs[0];
    double** const parameters = &d_bond_data[level_number].parameters[0];
    double* const X_mean_node = X_mean_data->getGhostedLocalFormVecArray()->data();
    double* const N_node = N_data->getGhostedLocalFormArray()->data();
    const double* const X_node = X_data->getGhostedLocalFormVecArray()->data();

    static const int BLOCKSIZE = 16; // this parameter needs to be tuned
    int k, kblock, kunroll, X_mastr_idx, X_slave_idx, N_mastr_idx, N_slave_idx;
    kblock = 0;
    for (; kblock < (num_bonds - 1) / BLOCKSIZE;
         ++kblock) // ensure that the last block is NOT handled by this first loop
    {
        PREFETCH_READ_NTA_BLOCK(petsc_mastr_node_idxs + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(petsc_slave_node_idxs + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(parameters + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        for (kunroll = 0; kunroll < BLOCKSIZE; ++kunroll)
        {
            k = kblock * BLOCKSIZE + kunroll;
            X_mastr_idx = petsc_mastr_node_idxs[k] * NDIM;
            X_slave_idx = petsc_slave_node_idxs[k] * NDIM;
            N_mastr_idx = petsc_mastr_node_idxs[k];
            N_slave_idx = petsc_slave_node_idxs[k];
#if !defined(NDEBUG)
            TBOX_ASSERT(X_mastr_idx != X_slave_idx);
            TBOX_ASSERT(N_mastr_idx != N_slave_idx);
#endif
            PREFETCH_READ_NTA_NDIM_BLOCK(X_mean_node + NDIM * petsc_mastr_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(X_mean_node + NDIM * petsc_slave_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(X_node + NDIM * petsc_mastr_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(X_node + NDIM * petsc_slave_node_idxs[k + 1]);
            PREFETCH_READ_NTA(parameters[k + 1]);
            PREFETCH_READ_NTA(N_node + k + 1);

            const double* bond_params = parameters[k];
            const double& fail = bond_params[4];

            N_node[N_mastr_idx] += fail * 1.0;
            N_node[N_slave_idx] += fail * 1.0;
            for (int d = 0; d < NDIM; ++d)
            {
                X_mean_node[X_mastr_idx + d] += fail * X_node[X_slave_idx + d];
                X_mean_node[X_slave_idx + d] += fail * X_node[X_mastr_idx + d];
            }
        }
    }
    for (k = kblock * BLOCKSIZE; k < num_bonds; ++k)
    {
        X_mastr_idx = petsc_mastr_node_idxs[k] * NDIM;
        X_slave_idx = petsc_slave_node_idxs[k] * NDIM;
        N_mastr_idx = petsc_mastr_node_idxs[k];
        N_slave_idx = petsc_slave_node_idxs[k];
#if !defined(NDEBUG)
        TBOX_ASSERT(X_mastr_idx != X_slave_idx);
        TBOX_ASSERT(N_mastr_idx != N_slave_idx);
#endif

        const double* bond_params = parameters[k];
        const double& fail = bond_params[4];

        N_node[N_mastr_idx] += fail * 1.0;
        N_node[N_slave_idx] += fail * 1.0;
        for (int d = 0; d < NDIM; ++d)
        {
            X_mean_node[X_mastr_idx + d] += fail * X_node[X_slave_idx + d];
            X_mean_node[X_slave_idx + d] += fail * X_node[X_mastr_idx + d];
        }
    }

    X_data->restoreArrays();
    X_mean_data->restoreArrays();
    N_data->restoreArrays();

    return;
} // computeMeanPosition

void
IBPDForceGen::computeShapeTensor(Pointer<LData> B_data,
                                 Pointer<LData> X0_data,
                                 Pointer<PatchHierarchy<NDIM> > hierarchy,
                                 int level_number,
                                 double /*data_time*/,
                                 LDataManager* /*l_data_manager*/)
{
    const int num_bonds = static_cast<int>(d_bond_data[level_number].lag_mastr_node_idxs.size());
    if (num_bonds == 0) return;

    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    const IntVector<NDIM>& ratio = level->getRatio();
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
    const double* dx0 = grid_geom->getDx();
    double dx[NDIM];
    for (int d = 0; d < NDIM; ++d)
    {
        dx[d] = dx0[d] / ratio[d];
    }
    const double horizon = d_horizon * (*std::max_element(dx, dx + NDIM));
    const double delta = d_ds * (*std::min_element(dx, dx + NDIM));

    const int* const petsc_mastr_node_idxs = &d_bond_data[level_number].petsc_mastr_node_idxs[0];
    const int* const petsc_slave_node_idxs = &d_bond_data[level_number].petsc_slave_node_idxs[0];
    double** const parameters = &d_bond_data[level_number].parameters[0];
    const BondInfluenceFcnPtr* const force_inf_fcns = &d_bond_data[level_number].force_inf_fcns[0];
    const BondVolFracFcnPtr* const force_vol_frac_fcns = &d_bond_data[level_number].force_vol_frac_fcns[0];
    double* const B_node = B_data->getGhostedLocalFormVecArray()->data();
    const double* const X0_node = X0_data->getGhostedLocalFormVecArray()->data();

    static const int BLOCKSIZE = 16; // this parameter needs to be tuned
    int k, kblock, kunroll, X_mastr_idx, X_slave_idx, B_mastr_idx, B_slave_idx;
    double Q0[NDIM], R0, W, vol_frac;
    kblock = 0;
    for (; kblock < (num_bonds - 1) / BLOCKSIZE;
         ++kblock) // ensure that the last block is NOT handled by this first loop
    {
        PREFETCH_READ_NTA_BLOCK(petsc_mastr_node_idxs + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(petsc_slave_node_idxs + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(parameters + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        for (kunroll = 0; kunroll < BLOCKSIZE; ++kunroll)
        {
            k = kblock * BLOCKSIZE + kunroll;
            X_mastr_idx = petsc_mastr_node_idxs[k] * NDIM;
            X_slave_idx = petsc_slave_node_idxs[k] * NDIM;
            B_mastr_idx = petsc_mastr_node_idxs[k] * NDIM * NDIM;
            B_slave_idx = petsc_slave_node_idxs[k] * NDIM * NDIM;
#if !defined(NDEBUG)
            TBOX_ASSERT(X_mastr_idx != X_slave_idx);
            TBOX_ASSERT(B_mastr_idx != B_slave_idx);
#endif
            PREFETCH_READ_NTA_BLOCK(B_node + NDIM * NDIM * petsc_mastr_node_idxs[k + 1], NDIM * NDIM);
            PREFETCH_READ_NTA_BLOCK(B_node + NDIM * NDIM * petsc_slave_node_idxs[k + 1], NDIM * NDIM);
            PREFETCH_READ_NTA_NDIM_BLOCK(X0_node + NDIM * petsc_mastr_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(X0_node + NDIM * petsc_slave_node_idxs[k + 1]);
            PREFETCH_READ_NTA(parameters[k + 1]);

            Q0[0] = X0_node[X_slave_idx + 0] - X0_node[X_mastr_idx + 0];
            Q0[1] = X0_node[X_slave_idx + 1] - X0_node[X_mastr_idx + 1];
#if (NDIM == 3)
            Q0[2] = X0_node[X_slave_idx + 2] - X0_node[X_mastr_idx + 2];
#endif

            const double* bond_params = parameters[k];
            const double& vol_mastr = bond_params[2];
            const double& vol_slave = bond_params[3];
            const double& fail = bond_params[4];
            R0 = parameters[k][1];
            W = force_inf_fcns[k](R0, delta);
            vol_frac = force_vol_frac_fcns[k](R0, horizon, delta);

            Eigen::Map<IBTK::Vector> eig_Q0(Q0);
            Eigen::Map<Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> > eig_B_mastr(&B_node[B_mastr_idx]),
                eig_B_slave(&B_node[B_slave_idx]);
            eig_B_mastr.noalias() += fail * W * vol_frac * vol_slave * eig_Q0 * eig_Q0.transpose();
            eig_B_slave.noalias() += fail * W * vol_frac * vol_mastr * eig_Q0 * eig_Q0.transpose();
        }
    }
    for (k = kblock * BLOCKSIZE; k < num_bonds; ++k)
    {
        X_mastr_idx = petsc_mastr_node_idxs[k] * NDIM;
        X_slave_idx = petsc_slave_node_idxs[k] * NDIM;
        B_mastr_idx = petsc_mastr_node_idxs[k] * NDIM * NDIM;
        B_slave_idx = petsc_slave_node_idxs[k] * NDIM * NDIM;
#if !defined(NDEBUG)
        TBOX_ASSERT(X_mastr_idx != X_slave_idx);
        TBOX_ASSERT(B_mastr_idx != B_slave_idx);
#endif

        Q0[0] = X0_node[X_slave_idx + 0] - X0_node[X_mastr_idx + 0];
        Q0[1] = X0_node[X_slave_idx + 1] - X0_node[X_mastr_idx + 1];
#if (NDIM == 3)
        Q0[2] = X0_node[X_slave_idx + 2] - X0_node[X_mastr_idx + 2];
#endif

        const double* bond_params = parameters[k];
        const double& vol_mastr = bond_params[2];
        const double& vol_slave = bond_params[3];
        const double& fail = bond_params[4];
        R0 = parameters[k][1];
        W = force_inf_fcns[k](R0, delta);
        vol_frac = force_vol_frac_fcns[k](R0, horizon, delta);

        Eigen::Map<IBTK::Vector> eig_Q0(Q0);
        Eigen::Map<Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> > eig_B_mastr(&B_node[B_mastr_idx]),
            eig_B_slave(&B_node[B_slave_idx]);
        eig_B_mastr.noalias() += fail * W * vol_frac * vol_slave * eig_Q0 * eig_Q0.transpose();
        eig_B_slave.noalias() += fail * W * vol_frac * vol_mastr * eig_Q0 * eig_Q0.transpose();
    }

    X0_data->restoreArrays();
    B_data->restoreArrays();

    return;
} // computeShapeTensor

void
IBPDForceGen::computeDeformationGradientTensor(Pointer<LData> FF_data,
                                               Pointer<LData> X_data,
                                               Pointer<LData> X0_data,
                                               Pointer<PatchHierarchy<NDIM> > hierarchy,
                                               int level_number,
                                               double /*data_time*/,
                                               LDataManager* /*l_data_manager*/)
{
    const int num_bonds = static_cast<int>(d_bond_data[level_number].lag_mastr_node_idxs.size());
    if (num_bonds == 0) return;

    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    const IntVector<NDIM>& ratio = level->getRatio();
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
    const double* dx0 = grid_geom->getDx();
    double dx[NDIM];
    for (int d = 0; d < NDIM; ++d)
    {
        dx[d] = dx0[d] / ratio[d];
    }
    const double horizon = d_horizon * (*std::min_element(dx, dx + NDIM));
    const double delta = d_ds * (*std::min_element(dx, dx + NDIM));

    const int* const petsc_mastr_node_idxs = &d_bond_data[level_number].petsc_mastr_node_idxs[0];
    const int* const petsc_slave_node_idxs = &d_bond_data[level_number].petsc_slave_node_idxs[0];
    double** const parameters = &d_bond_data[level_number].parameters[0];
    const BondInfluenceFcnPtr* const force_inf_fcns = &d_bond_data[level_number].force_inf_fcns[0];
    const BondVolFracFcnPtr* const force_vol_frac_fcns = &d_bond_data[level_number].force_vol_frac_fcns[0];
    double* const FF_node = FF_data->getGhostedLocalFormVecArray()->data();
    const double* const X_node = X_data->getGhostedLocalFormVecArray()->data();
    const double* const X0_node = X0_data->getGhostedLocalFormVecArray()->data();

    static const int BLOCKSIZE = 16; // this parameter needs to be tuned
    int k, kblock, kunroll, X_mastr_idx, X_slave_idx, FF_mastr_idx, FF_slave_idx;
    double Q[NDIM], Q0[NDIM], R0, W, vol_frac;
    kblock = 0;
    for (; kblock < (num_bonds - 1) / BLOCKSIZE;
         ++kblock) // ensure that the last block is NOT handled by this first loop
    {
        PREFETCH_READ_NTA_BLOCK(petsc_mastr_node_idxs + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(petsc_slave_node_idxs + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(parameters + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        for (kunroll = 0; kunroll < BLOCKSIZE; ++kunroll)
        {
            k = kblock * BLOCKSIZE + kunroll;
            X_mastr_idx = petsc_mastr_node_idxs[k] * NDIM;
            X_slave_idx = petsc_slave_node_idxs[k] * NDIM;
            FF_mastr_idx = petsc_mastr_node_idxs[k] * NDIM * NDIM;
            FF_slave_idx = petsc_slave_node_idxs[k] * NDIM * NDIM;
#if !defined(NDEBUG)
            TBOX_ASSERT(X_mastr_idx != X_slave_idx);
            TBOX_ASSERT(FF_mastr_idx != FF_slave_idx);
#endif
            PREFETCH_READ_NTA_BLOCK(FF_node + NDIM * NDIM * petsc_mastr_node_idxs[k + 1], NDIM * NDIM);
            PREFETCH_READ_NTA_BLOCK(FF_node + NDIM * NDIM * petsc_slave_node_idxs[k + 1], NDIM * NDIM);
            PREFETCH_READ_NTA_NDIM_BLOCK(X_node + NDIM * petsc_mastr_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(X_node + NDIM * petsc_slave_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(X0_node + NDIM * petsc_mastr_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(X0_node + NDIM * petsc_slave_node_idxs[k + 1]);
            PREFETCH_READ_NTA(parameters[k + 1]);

            Q[0] = X_node[X_slave_idx + 0] - X_node[X_mastr_idx + 0];
            Q[1] = X_node[X_slave_idx + 1] - X_node[X_mastr_idx + 1];
#if (NDIM == 3)
            Q[2] = X_node[X_slave_idx + 2] - X_node[X_mastr_idx + 2];
#endif

            Q0[0] = X0_node[X_slave_idx + 0] - X0_node[X_mastr_idx + 0];
            Q0[1] = X0_node[X_slave_idx + 1] - X0_node[X_mastr_idx + 1];
#if (NDIM == 3)
            Q0[2] = X0_node[X_slave_idx + 2] - X0_node[X_mastr_idx + 2];
#endif

            const double* bond_params = parameters[k];
            const double& vol_mastr = bond_params[2];
            const double& vol_slave = bond_params[3];
            const double& fail = bond_params[4];
            R0 = parameters[k][1];
            W = force_inf_fcns[k](R0, delta);
            vol_frac = force_vol_frac_fcns[k](R0, horizon, delta);

            Eigen::Map<const IBTK::Vector> eig_Q(Q), eig_Q0(Q0);
            Eigen::Map<Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> > eig_FF_mastr(&FF_node[FF_mastr_idx]),
                eig_FF_slave(&FF_node[FF_slave_idx]);
            eig_FF_mastr.noalias() += fail * W * vol_frac * vol_slave * eig_Q * eig_Q0.transpose();
            eig_FF_slave.noalias() += fail * W * vol_frac * vol_mastr * eig_Q * eig_Q0.transpose();
        }
    }
    for (k = kblock * BLOCKSIZE; k < num_bonds; ++k)
    {
        X_mastr_idx = petsc_mastr_node_idxs[k] * NDIM;
        X_slave_idx = petsc_slave_node_idxs[k] * NDIM;
        FF_mastr_idx = petsc_mastr_node_idxs[k] * NDIM * NDIM;
        FF_slave_idx = petsc_slave_node_idxs[k] * NDIM * NDIM;
#if !defined(NDEBUG)
        TBOX_ASSERT(X_mastr_idx != X_slave_idx);
        TBOX_ASSERT(FF_mastr_idx != FF_slave_idx);
#endif

        Q[0] = X_node[X_slave_idx + 0] - X_node[X_mastr_idx + 0];
        Q[1] = X_node[X_slave_idx + 1] - X_node[X_mastr_idx + 1];
#if (NDIM == 3)
        Q[2] = X_node[X_slave_idx + 2] - X_node[X_mastr_idx + 2];
#endif

        Q0[0] = X0_node[X_slave_idx + 0] - X0_node[X_mastr_idx + 0];
        Q0[1] = X0_node[X_slave_idx + 1] - X0_node[X_mastr_idx + 1];
#if (NDIM == 3)
        Q0[2] = X0_node[X_slave_idx + 2] - X0_node[X_mastr_idx + 2];
#endif

        const double* bond_params = parameters[k];
        const double& vol_mastr = bond_params[2];
        const double& vol_slave = bond_params[3];
        const double& fail = bond_params[4];
        R0 = parameters[k][1];
        W = force_inf_fcns[k](R0, delta);
        vol_frac = force_vol_frac_fcns[k](R0, horizon, delta);

        Eigen::Map<const IBTK::Vector> eig_Q(Q), eig_Q0(Q0);
        Eigen::Map<Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> > eig_FF_mastr(&FF_node[FF_mastr_idx]),
            eig_FF_slave(&FF_node[FF_slave_idx]);
        eig_FF_mastr.noalias() += fail * W * vol_frac * vol_slave * eig_Q * eig_Q0.transpose();
        eig_FF_slave.noalias() += fail * W * vol_frac * vol_mastr * eig_Q * eig_Q0.transpose();
    }

    X_data->restoreArrays();
    X0_data->restoreArrays();
    FF_data->restoreArrays();

    return;
} // computeDeformationGradientTensor

void
IBPDForceGen::initializeBondLevelData(std::set<int>& nonlocal_petsc_idx_set,
                                      const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                      const int level_number,
                                      const double /*init_data_time*/,
                                      const bool /*initial_time*/,
                                      LDataManager* const l_data_manager)
{
    std::vector<int>& lag_mastr_node_idxs = d_bond_data[level_number].lag_mastr_node_idxs;
    std::vector<int>& lag_slave_node_idxs = d_bond_data[level_number].lag_slave_node_idxs;
    std::vector<int>& petsc_mastr_node_idxs = d_bond_data[level_number].petsc_mastr_node_idxs;
    std::vector<int>& petsc_slave_node_idxs = d_bond_data[level_number].petsc_slave_node_idxs;
    std::vector<int>& petsc_global_mastr_node_idxs = d_bond_data[level_number].petsc_global_mastr_node_idxs;
    std::vector<int>& petsc_global_slave_node_idxs = d_bond_data[level_number].petsc_global_slave_node_idxs;
    std::vector<BondForceDamageFcnPtr>& force_dmg_fcns = d_bond_data[level_number].force_dmg_fcns;
    std::vector<BondPK1FcnPtr>& force_PK1_fcns = d_bond_data[level_number].force_PK1_fcns;
    std::vector<BondInfluenceFcnPtr>& force_inf_fcns = d_bond_data[level_number].force_inf_fcns;
    std::vector<BondVolFracFcnPtr>& force_vol_frac_fcns = d_bond_data[level_number].force_vol_frac_fcns;
    std::vector<double*>& parameters = d_bond_data[level_number].parameters;

    // The LMesh object provides the set of local Lagrangian nodes.
    const Pointer<LMesh> mesh = l_data_manager->getLMesh(level_number);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
    const int num_local_nodes = static_cast<int>(local_nodes.size());

    // Determine how many bonds are associated with the present MPI process.
    unsigned int num_bonds = 0;
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
        const IBSpringForceSpec* const force_spec = node_idx->getNodeDataItem<IBSpringForceSpec>();
        if (force_spec) num_bonds += force_spec->getNumberOfSprings();
    }

    // Resize arrays for storing cached values used to compute bond forces.
    lag_mastr_node_idxs.resize(num_bonds);
    lag_slave_node_idxs.resize(num_bonds);
    petsc_mastr_node_idxs.resize(num_bonds);
    petsc_slave_node_idxs.resize(num_bonds);
    petsc_global_mastr_node_idxs.resize(num_bonds);
    petsc_global_slave_node_idxs.resize(num_bonds);
    force_dmg_fcns.resize(num_bonds);
    force_PK1_fcns.resize(num_bonds);
    force_inf_fcns.resize(num_bonds);
    force_vol_frac_fcns.resize(num_bonds);
    parameters.resize(num_bonds);

    // Setup the data structures used to compute bond forces.
    int current_bond = 0;
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
        IBSpringForceSpec* const force_spec = node_idx->getNodeDataItem<IBSpringForceSpec>();
        if (!force_spec) continue;

        const int lag_idx = node_idx->getLagrangianIndex();
#if !defined(NDEBUG)
        TBOX_ASSERT(lag_idx == force_spec->getMasterNodeIndex());
#endif
        const int petsc_idx = node_idx->getGlobalPETScIndex();
        const std::vector<int>& slv = force_spec->getSlaveNodeIndices();
        const std::vector<int>& fcn = force_spec->getForceFunctionIndices();
        std::vector<std::vector<double> >& params = force_spec->getParameters();
        const unsigned int n_mastr_bonds = force_spec->getNumberOfSprings();
#if !defined(NDEBUG)
        TBOX_ASSERT(n_mastr_bonds == slv.size());
        TBOX_ASSERT(n_mastr_bonds == params.size());
#endif
        for (unsigned int k = 0; k < n_mastr_bonds; ++k)
        {
            lag_mastr_node_idxs[current_bond] = lag_idx;
            lag_slave_node_idxs[current_bond] = slv[k];
            petsc_mastr_node_idxs[current_bond] = petsc_idx;
            force_dmg_fcns[current_bond] = d_bond_force_damage_fcn_map[fcn[k]];
            force_PK1_fcns[current_bond] = d_bond_PK1_fcn_map[fcn[k]];
            force_inf_fcns[current_bond] = d_bond_inf_fcn_map[fcn[k]];
            force_vol_frac_fcns[current_bond] = d_bond_vol_frac_fcn_map[fcn[k]];
            parameters[current_bond] = params.empty() ? NULL : &params[k][0];
            ++current_bond;
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
    for (unsigned int k = 0; k < petsc_slave_node_idxs.size(); ++k)
    {
        const int idx = petsc_slave_node_idxs[k];
        if (UNLIKELY(idx < global_node_offset || idx >= global_node_offset + num_local_nodes))
        {
            nonlocal_petsc_idx_set.insert(idx);
        }
    }
    return;
} // initializeBondLevelData

void
IBPDForceGen::computeLagrangianBondForceAndDamage(Pointer<LData> F_data,
                                                  Pointer<LData> D_data,
                                                  Pointer<LData> X_data,
                                                  Pointer<LData> X0_data,
                                                  Pointer<LData> FF_data,
                                                  Pointer<LData> B_data,
                                                  const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                  const int level_number,
                                                  const double /*data_time*/,
                                                  LDataManager* const /*l_data_manager*/)
{
    const int num_bonds = static_cast<int>(d_bond_data[level_number].lag_mastr_node_idxs.size());
    if (num_bonds == 0) return;

    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    const IntVector<NDIM>& ratio = level->getRatio();
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
    const double* dx0 = grid_geom->getDx();
    double dx[NDIM];
    for (int d = 0; d < NDIM; ++d)
    {
        dx[d] = dx0[d] / ratio[d];
    }
    const double horizon = d_horizon * (*std::max_element(dx, dx + NDIM));
    const double delta = 0.5 * d_ds * (*std::min_element(dx, dx + NDIM));

    const int* const lag_mastr_node_idxs = &d_bond_data[level_number].lag_mastr_node_idxs[0];
    const int* const lag_slave_node_idxs = &d_bond_data[level_number].lag_slave_node_idxs[0];
    const int* const petsc_mastr_node_idxs = &d_bond_data[level_number].petsc_mastr_node_idxs[0];
    const int* const petsc_slave_node_idxs = &d_bond_data[level_number].petsc_slave_node_idxs[0];
    double** const parameters = &d_bond_data[level_number].parameters[0];
    const BondForceDamageFcnPtr* const force_dmg_fcns = &d_bond_data[level_number].force_dmg_fcns[0];
    const BondInfluenceFcnPtr* const force_inf_fcns = &d_bond_data[level_number].force_inf_fcns[0];
    const BondVolFracFcnPtr* const force_vol_frac_fcns = &d_bond_data[level_number].force_vol_frac_fcns[0];
    double* const F_node = F_data->getGhostedLocalFormVecArray()->data();
    double* const D_node = D_data->getGhostedLocalFormVecArray()->data();
    const double* const X_node = X_data->getGhostedLocalFormVecArray()->data();
    const double* const X0_node = X0_data->getGhostedLocalFormVecArray()->data();
    const double* const FF_node = FF_data->getGhostedLocalFormVecArray()->data();
    const double* const B_node = B_data->getGhostedLocalFormVecArray()->data();

    static const int BLOCKSIZE = 16; // this parameter needs to be tuned
    int k, kblock, kunroll, X_F_mastr_idx, X_F_slave_idx, FF_B_mastr_idx, FF_B_slave_idx, dmg_mastr_idx, dmg_slave_idx;
    double Q[NDIM], R, R0, W, vol_frac;
    kblock = 0;
    for (; kblock < (num_bonds - 1) / BLOCKSIZE;
         ++kblock) // ensure that the last block is NOT handled by this first loop
    {
        PREFETCH_READ_NTA_BLOCK(lag_mastr_node_idxs + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(lag_slave_node_idxs + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(petsc_mastr_node_idxs + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(petsc_slave_node_idxs + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(parameters + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        for (kunroll = 0; kunroll < BLOCKSIZE; ++kunroll)
        {
            k = kblock * BLOCKSIZE + kunroll;
            X_F_mastr_idx = petsc_mastr_node_idxs[k] * NDIM;
            X_F_slave_idx = petsc_slave_node_idxs[k] * NDIM;
            dmg_mastr_idx = petsc_mastr_node_idxs[k] * 2;
            dmg_slave_idx = petsc_slave_node_idxs[k] * 2;
            FF_B_mastr_idx = petsc_mastr_node_idxs[k] * NDIM * NDIM;
            FF_B_slave_idx = petsc_slave_node_idxs[k] * NDIM * NDIM;
#if !defined(NDEBUG)
            TBOX_ASSERT(X_F_mastr_idx != X_F_slave_idx);
            TBOX_ASSERT(dmg_mastr_idx != dmg_slave_idx);
            TBOX_ASSERT(FF_B_mastr_idx != FF_B_slave_idx);
#endif
            PREFETCH_READ_NTA_NDIM_BLOCK(F_node + NDIM * petsc_mastr_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(F_node + NDIM * petsc_slave_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(X_node + NDIM * petsc_mastr_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(X_node + NDIM * petsc_slave_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(X0_node + NDIM * petsc_mastr_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(X0_node + NDIM * petsc_slave_node_idxs[k + 1]);
            PREFETCH_READ_NTA_BLOCK(FF_node + NDIM * NDIM * petsc_mastr_node_idxs[k + 1], NDIM * NDIM);
            PREFETCH_READ_NTA_BLOCK(FF_node + NDIM * NDIM * petsc_slave_node_idxs[k + 1], NDIM * NDIM);
            PREFETCH_READ_NTA_BLOCK(B_node + NDIM * NDIM * petsc_mastr_node_idxs[k + 1], NDIM * NDIM);
            PREFETCH_READ_NTA_BLOCK(B_node + NDIM * NDIM * petsc_slave_node_idxs[k + 1], NDIM * NDIM);
            PREFETCH_READ_NTA(parameters[k + 1]);

            Q[0] = X_node[X_F_slave_idx + 0] - X_node[X_F_mastr_idx + 0];
            Q[1] = X_node[X_F_slave_idx + 1] - X_node[X_F_mastr_idx + 1];
#if (NDIM == 3)
            Q[2] = X_node[X_F_slave_idx + 2] - X_node[X_F_mastr_idx + 2];
#endif

#if (NDIM == 2)
            R = sqrt(Q[0] * Q[0] + Q[1] * Q[1]);
#endif
#if (NDIM == 3)
            R = sqrt(Q[0] * Q[0] + Q[1] * Q[1] + Q[2] * Q[2]);
#endif

            if (UNLIKELY(R < std::numeric_limits<double>::epsilon())) continue;

            Eigen::Map<const IBTK::Vector> eig_X0_mastr(&X0_node[X_F_mastr_idx]), eig_X0_slave(&X0_node[X_F_slave_idx]);
            Eigen::Map<const IBTK::Vector> eig_X_mastr(&X_node[X_F_mastr_idx]), eig_X_slave(&X_node[X_F_slave_idx]);
            Eigen::Map<IBTK::Vector> eig_F_mastr(&F_node[X_F_mastr_idx]), eig_F_slave(&F_node[X_F_slave_idx]);
            Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> > eig_B_mastr(&B_node[FF_B_mastr_idx]),
                eig_B_slave(&B_node[FF_B_slave_idx]);
            Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> > eig_FF_mastr(
                &FF_node[FF_B_mastr_idx]),
                eig_FF_slave(&FF_node[FF_B_slave_idx]);

            R0 = parameters[k][1];
            W = force_inf_fcns[k](R0, delta);
            vol_frac = force_vol_frac_fcns[k](R0, horizon, delta);

            Eigen::Vector4d dmg = force_dmg_fcns[k](horizon,
                                                    delta,
                                                    W,
                                                    vol_frac,
                                                    parameters[k],
                                                    eig_X0_mastr,
                                                    eig_X0_slave,
                                                    eig_X_mastr,
                                                    eig_X_slave,
                                                    eig_FF_mastr,
                                                    eig_FF_slave,
                                                    eig_B_mastr,
                                                    eig_B_slave,
                                                    eig_F_mastr,
                                                    eig_F_slave,
                                                    lag_mastr_node_idxs[k],
                                                    lag_slave_node_idxs[k]);

            D_node[dmg_mastr_idx + 0] += dmg[0];
            D_node[dmg_mastr_idx + 1] += dmg[1];
            D_node[dmg_slave_idx + 0] += dmg[2];
            D_node[dmg_slave_idx + 1] += dmg[3];
        }
    }
    for (k = kblock * BLOCKSIZE; k < num_bonds; ++k)
    {
        X_F_mastr_idx = petsc_mastr_node_idxs[k] * NDIM;
        X_F_slave_idx = petsc_slave_node_idxs[k] * NDIM;
        dmg_mastr_idx = petsc_mastr_node_idxs[k] * 2;
        dmg_slave_idx = petsc_slave_node_idxs[k] * 2;
        FF_B_mastr_idx = petsc_mastr_node_idxs[k] * NDIM * NDIM;
        FF_B_slave_idx = petsc_slave_node_idxs[k] * NDIM * NDIM;
#if !defined(NDEBUG)
        TBOX_ASSERT(X_F_mastr_idx != X_F_slave_idx);
        TBOX_ASSERT(dmg_mastr_idx != dmg_slave_idx);
        TBOX_ASSERT(FF_B_mastr_idx != FF_B_slave_idx);
#endif
        Q[0] = X_node[X_F_slave_idx + 0] - X_node[X_F_mastr_idx + 0];
        Q[1] = X_node[X_F_slave_idx + 1] - X_node[X_F_mastr_idx + 1];
#if (NDIM == 3)
        Q[2] = X_node[X_F_slave_idx + 2] - X_node[X_F_mastr_idx + 2];
#endif

#if (NDIM == 2)
        R = sqrt(Q[0] * Q[0] + Q[1] * Q[1]);
#endif
#if (NDIM == 3)
        R = sqrt(Q[0] * Q[0] + Q[1] * Q[1] + Q[2] * Q[2]);
#endif

        if (UNLIKELY(R < std::numeric_limits<double>::epsilon())) continue;

        Eigen::Map<const IBTK::Vector> eig_X0_mastr(&X0_node[X_F_mastr_idx]), eig_X0_slave(&X0_node[X_F_slave_idx]);
        Eigen::Map<const IBTK::Vector> eig_X_mastr(&X0_node[X_F_mastr_idx]), eig_X_slave(&X0_node[X_F_slave_idx]);
        Eigen::Map<IBTK::Vector> eig_F_mastr(&F_node[X_F_mastr_idx]), eig_F_slave(&F_node[X_F_slave_idx]);
        Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> > eig_B_mastr(&B_node[FF_B_mastr_idx]),
            eig_B_slave(&B_node[FF_B_slave_idx]);
        Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> > eig_FF_mastr(&FF_node[FF_B_mastr_idx]),
            eig_FF_slave(&FF_node[FF_B_slave_idx]);

        R0 = parameters[k][1];
        W = force_inf_fcns[k](R0, delta);
        vol_frac = force_vol_frac_fcns[k](R0, horizon, delta);

        Eigen::Vector4d dmg = force_dmg_fcns[k](horizon,
                                                delta,
                                                W,
                                                vol_frac,
                                                parameters[k],
                                                eig_X0_mastr,
                                                eig_X0_slave,
                                                eig_X_mastr,
                                                eig_X_slave,
                                                eig_FF_mastr,
                                                eig_FF_slave,
                                                eig_B_mastr,
                                                eig_B_slave,
                                                eig_F_mastr,
                                                eig_F_slave,
                                                lag_mastr_node_idxs[k],
                                                lag_slave_node_idxs[k]);

        D_node[dmg_mastr_idx + 0] += dmg[0];
        D_node[dmg_mastr_idx + 1] += dmg[1];
        D_node[dmg_slave_idx + 0] += dmg[2];
        D_node[dmg_slave_idx + 1] += dmg[3];
    }

    F_data->restoreArrays();
    D_data->restoreArrays();
    B_data->restoreArrays();
    X_data->restoreArrays();
    X0_data->restoreArrays();

    return;
} // computeLagrangianBondForceAndDamage

void
IBPDForceGen::initializeTargetPointLevelData(std::set<int>& /*nonlocal_petsc_idx_set*/,
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
    unsigned int num_target_points = 0;
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
        const IBTargetPointForceSpec* const force_spec = node_idx->getNodeDataItem<IBTargetPointForceSpec>();
        if (force_spec) num_target_points += 1;
    }

    // Resize arrays for storing cached values used to compute target point
    // forces.
    petsc_node_idxs.resize(num_target_points);
    petsc_global_node_idxs.resize(num_target_points);
    kappa.resize(num_target_points);
    eta.resize(num_target_points);
    X0.resize(num_target_points);

    // Setup the data structures used to compute target point forces.
    int current_target_point = 0;
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
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
IBPDForceGen::computeLagrangianTargetPointForce(Pointer<LData> F_data,
                                                Pointer<LData> X_data,
                                                Pointer<LData> U_data,
                                                const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                                const int level_number,
                                                const double /*data_time*/,
                                                LDataManager* const /*l_data_manager*/)
{
    const int num_target_points = static_cast<int>(d_target_point_data[level_number].petsc_node_idxs.size());
    if (num_target_points == 0) return;
    const int* const petsc_node_idxs = &d_target_point_data[level_number].petsc_node_idxs[0];
    const int* const petsc_global_node_idxs = &d_target_point_data[level_number].petsc_global_node_idxs[0];
    const double** const kappa = &d_target_point_data[level_number].kappa[0];
    const double** const eta = &d_target_point_data[level_number].eta[0];
    const Point** const X0 = &d_target_point_data[level_number].X0[0];
    double* const F_node = F_data->getLocalFormVecArray()->data();
    const double* const X_node = X_data->getLocalFormVecArray()->data();
    const double* const U_node = U_data->getLocalFormVecArray()->data();

    static const int BLOCKSIZE = 16; // This parameter needs to be tuned.
    int k, kblock, kunroll, local_idx, lag_idx;
    double K, E;
    const double* X_target;
    kblock = 0;
    for (; kblock < (num_target_points - 1) / BLOCKSIZE;
         ++kblock) // ensure that the last block is NOT handled by this first loop
    {
        PREFETCH_READ_NTA_BLOCK(petsc_node_idxs + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(petsc_global_node_idxs + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(kappa + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(eta + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        PREFETCH_READ_NTA_BLOCK(X0 + BLOCKSIZE * (kblock + 1), BLOCKSIZE);
        for (kunroll = 0; kunroll < BLOCKSIZE; ++kunroll)
        {
            k = kblock * BLOCKSIZE + kunroll;
            local_idx = NDIM * petsc_node_idxs[k];
            lag_idx = petsc_global_node_idxs[k];
            PREFETCH_READ_NTA_NDIM_BLOCK(F_node + NDIM * petsc_node_idxs[k + 1]);
            PREFETCH_READ_NTA_NDIM_BLOCK(X_node + NDIM * petsc_node_idxs[k + 1]);
            PREFETCH_READ_NTA(kappa[k + 1]);
            PREFETCH_READ_NTA(eta[k + 1]);
            PREFETCH_READ_NTA(X0[k + 1]);
            K = *kappa[k];
            E = *eta[k];
            X_target = X0[k]->data();

            Eigen::Map<const IBTK::Vector> eig_X(&X_node[local_idx]);
            Eigen::Map<const IBTK::Vector> eig_X_target(X_target);
            Eigen::Map<const IBTK::Vector> eig_U(&U_node[local_idx]);
            Eigen::Map<IBTK::Vector> eig_F(&F_node[local_idx]);

            d_target_point_force_fcn(eig_X, eig_X_target, eig_U, K, E, lag_idx, eig_F);
        }
    }
    for (k = kblock * BLOCKSIZE; k < num_target_points; ++k)
    {
        local_idx = NDIM * petsc_node_idxs[k];
        lag_idx = petsc_global_node_idxs[k];
        K = *kappa[k];
        E = *eta[k];
        X_target = X0[k]->data();

        Eigen::Map<const IBTK::Vector> eig_X(&X_node[local_idx]);
        Eigen::Map<const IBTK::Vector> eig_X_target(X_target);
        Eigen::Map<const IBTK::Vector> eig_U(&U_node[local_idx]);
        Eigen::Map<IBTK::Vector> eig_F(&F_node[local_idx]);

        d_target_point_force_fcn(eig_X, eig_X_target, eig_U, K, E, lag_idx, eig_F);
    }

    F_data->restoreArrays();
    X_data->restoreArrays();
    U_data->restoreArrays();
    return;
} // computeLagrangianTargetPointForce

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
