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

#include "Eigen/Dense"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "boost/multi_array.hpp"
#include "ibamr/IBPDForceGen.h"
#include "ibamr/IBSpringForceSpec.h"
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

/*static const int interior_begin = 0;
static const int interior_end = 39999;
static const int bottom_begin = 40000;
static const int bottom_end = 40617;
static const int top_begin = 40618;
static const int top_end = 41235;
static const int left_begin = 41236;
static const int left_end = 41835;
static const int right_begin = 41836;
static const int right_end = 42435;

static const double dens = 1.0;
static const double DX  = 1.0/199.0;*/

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


template <typename mat_type, typename vec_type>
void
get_PK1(mat_type& PK1, const Eigen::Map<const mat_type>& FF, const vec_type& X, int lag_idx)
{
    const double L = 2.0;
    const double M = 1.0;

    static const mat_type II = mat_type::Identity();
    mat_type E = 0.5 * (FF.transpose() * FF - II);
    const double trE = E.trace();
    PK1 = L * trE * FF + 2 * M * FF * E;

    return;

} // get_PK1

template <typename mat_type>
void
get_trac(IBTK::Vector& trac, const Eigen::Map<const mat_type>& eig_FF, const IBTK::Vector& normal)
{
    mat_type PK1;
    IBTK::Vector X;
    get_PK1(PK1, eig_FF, X, -1);

    trac = PK1 * normal;

    return;
}

template <typename vec_type, typename mat_type>
Eigen::Vector4d
pd_pk1force_damage_function(const double horizon,
                            const double delta,
                            const double R,
                            double* parameters,
                            const Eigen::Map<const vec_type>& X0_mastr,
                            const Eigen::Map<const vec_type>& X0_slave,
                            const Eigen::Map<const mat_type>& FF_mastr,
                            const Eigen::Map<const mat_type>& FF_slave,
                            const Eigen::Map<const mat_type>& B_mastr,
                            const Eigen::Map<const mat_type>& B_slave,
                            Eigen::Map<vec_type>& F_mastr,
                            Eigen::Map<vec_type>& F_slave,
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

    // Volume correction
    double vol_frac = 1.0;
    /*if (R0  <= (horizon - delta))
    {
        vol_frac = 1.0;
    }
    else if(R0  <= (horizon + delta))
    {
        vol_frac = (horizon + delta - R0)/(2.0 * delta);
    }
    else
    {
        vol_frac = 0.0;
    }*/

    // PK1 stress tensor
    mat_type PK1_mastr, PK1_slave;
    get_PK1(PK1_mastr, FF_mastr, X0_mastr, lag_mastr_node_idx);
    get_PK1(PK1_slave, FF_slave, X0_slave, lag_slave_node_idx);

    // Compute PD force.
    vec_type trac = (PK1_mastr*B_mastr + PK1_slave*B_slave)*(X0_slave - X0_mastr);
    F_mastr +=   vol_slave * trac;
    F_slave += - vol_mastr * trac;
    //vec_type trac = (PK1_mastr + PK1_slave) * (X0_slave - X0_mastr);
    //F_mastr += std::pow(vol_slave, -2.0 / 3.0) * trac;
    //F_slave += -std::pow(vol_mastr, -2.0 / 3.0) * trac;

    // Compute damage.
    const double stretch = (R - R0) / R0;
    Eigen::Vector4d D;
    if (std::fabs(stretch) > critical_stretch)
    {
        fail = 1.0;
    }
    D(0) = vol_slave * vol_frac * fail;
    D(1) = vol_slave * vol_frac;
    D(2) = vol_mastr * vol_frac * fail;
    D(3) = vol_mastr * vol_frac;

    return D;

} // pd_pk1force_damage_function
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBPDForceGen::IBPDForceGen(Pointer<Database> input_db) : d_horizon(3.0), d_ds(1.0)
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
    }

    return;
} // IBStandardForceGen

IBPDForceGen::~IBPDForceGen()
{
    // intentionally blank
    return;
} // ~IBPDForceGen

void
IBPDForceGen::initializeLevelData(const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                  const int level_number,
                                  const double init_data_time,
                                  const bool initial_time,
                                  LDataManager* const l_data_manager)
{
    if (!l_data_manager->levelContainsLagrangianData(level_number)) return;

#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
#endif
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Resize the vectors corresponding to data individually maintained for
    // separate levels of the patch hierarchy.
    const int new_size = std::max(level_number + 1, static_cast<int>(d_is_initialized.size()));

    d_bond_data.resize(new_size);
    d_X0_ghost_data.resize(new_size);
    d_X_ghost_data.resize(new_size);
    d_F_ghost_data.resize(new_size);
    d_dmg_ghost_data.resize(new_size);
    d_B_ghost_data.resize(new_size);
    d_FF_ghost_data.resize(new_size);
    d_dX_data.resize(new_size);
    d_is_initialized.resize(new_size, false);

    // Keep track of all of the nonlocal PETSc indices required to compute the
    // forces.
    std::set<int> nonlocal_petsc_idx_set;

    // Setup the cached data.
    initializeBondLevelData(
        nonlocal_petsc_idx_set, hierarchy, level_number, init_data_time, initial_time, l_data_manager);

    // Put the nonlocal PETSc indices into a vector.
    std::vector<int> nonlocal_petsc_idxs(nonlocal_petsc_idx_set.begin(), nonlocal_petsc_idx_set.end());

    // Put all cached PETSc node indices into local form.
    const int global_node_offset = l_data_manager->getGlobalNodeOffset(level_number);
    const int num_local_nodes = l_data_manager->getNumberOfLocalNodes(level_number);
    resetLocalPETScIndices(d_bond_data[level_number].petsc_mastr_node_idxs, global_node_offset, num_local_nodes);
    resetLocalOrNonlocalPETScIndices(
        d_bond_data[level_number].petsc_slave_node_idxs, global_node_offset, num_local_nodes, nonlocal_petsc_idxs);

    std::ostringstream X0_name_stream;
    X0_name_stream << "IBPDForceGen::X0_ghost_" << level_number;
    d_X0_ghost_data[level_number] = new LData(X0_name_stream.str(), num_local_nodes, NDIM, nonlocal_petsc_idxs);

    std::ostringstream X_name_stream;
    X_name_stream << "IBPDForceGen::X_ghost_" << level_number;
    d_X_ghost_data[level_number] = new LData(X_name_stream.str(), num_local_nodes, NDIM, nonlocal_petsc_idxs);

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

    // Indicate that the level data has been initialized.
    d_is_initialized[level_number] = true;
    return;
} // initializeLevelData

void
IBPDForceGen::computeLagrangianForceAndDamage(Pointer<LData> F_data,
                                              Pointer<LData> D_data,
                                              Pointer<LData> X_data,
                                              Pointer<LData> /*U_data*/,
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

    Pointer<LData> X0_ghost_data = d_X0_ghost_data[level_number];
    Pointer<LData> X0 = l_data_manager->getLData("X0", level_number);
    ierr = VecCopy(X0->getVec(), X0_ghost_data->getVec());
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(X0_ghost_data->getVec(), INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(X0_ghost_data->getVec(), INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);

    Pointer<LData> B_ghost_data = d_B_ghost_data[level_number];
    Vec B_ghost_local_form_vec;
    ierr = VecGhostGetLocalForm(B_ghost_data->getVec(), &B_ghost_local_form_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecSet(B_ghost_local_form_vec, 0.0);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostRestoreLocalForm(B_ghost_data->getVec(), &B_ghost_local_form_vec);
    IBTK_CHKERRQ(ierr);

    Pointer<LData> FF_ghost_data = d_FF_ghost_data[level_number];
    Vec FF_ghost_local_form_vec;
    ierr = VecGhostGetLocalForm(FF_ghost_data->getVec(), &FF_ghost_local_form_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecSet(FF_ghost_local_form_vec, 0.0);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostRestoreLocalForm(FF_ghost_data->getVec(), &FF_ghost_local_form_vec);
    IBTK_CHKERRQ(ierr);

    // Compute shape tensor.
    computeShapeTensor(B_ghost_data, X0_ghost_data, hierarchy, level_number, data_time, l_data_manager);
    ierr = VecGhostUpdateBegin(B_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(B_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);
    IBTK_CHKERRQ(ierr);

    {
        boost::multi_array_ref<double, 2>& B_ghost_data_array = *B_ghost_data->getLocalFormVecArray();
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int lag_idx = node_idx->getLagrangianIndex();
            const int local_idx = node_idx->getLocalPETScIndex();
            double* B = &B_ghost_data_array[local_idx][0];
            Eigen::Map<Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> > eig_B(B);

            // Scale the matrix.
            const double B_00 = eig_B(0, 0);
            eig_B *= (1.0 / B_00);

            // Invert the scaled matrix.
            bool invertible;
            eig_B.computeInverseWithCheck(eig_B, invertible);
            if (!invertible)
            {
                pout << " Lagrangian index = " << lag_idx << "\nScaled B tensor is \n" << eig_B << "\n";
                TBOX_ERROR("IBPDForceGen::computeLagrangianForceAndDamage() : Matrix inverse failed.\n");
            }

            // Scale the matrix back.
            eig_B *= (1.0 / B_00);
        }
        ierr = VecGhostUpdateBegin(B_ghost_data->getVec(), INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(B_ghost_data->getVec(), INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
        B_ghost_data->restoreArrays();
    }

    // Compute the deformation gradient tensor.
    // FF = int_k {Y outer X}_k . Inv{B}
    computeDeformationGradientTensor(
        FF_ghost_data, X_ghost_data, X0_ghost_data, hierarchy, level_number, data_time, l_data_manager);
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
    ierr = VecGhostUpdateBegin(F_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);
    IBTK_CHKERRQ(ierr);
    ierr = VecAXPY(F_data->getVec(), 1.0, F_ghost_data->getVec());
    ierr = VecGhostUpdateBegin(D_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);
    IBTK_CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D_ghost_data->getVec(), ADD_VALUES, SCATTER_REVERSE);
    IBTK_CHKERRQ(ierr);

    /*{
        boost::multi_array_ref<double, 2>& FF_ghost_data_array = *FF_ghost_data->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& F_data_array = *F_data->getLocalFormVecArray();
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int local_idx = node_idx->getLocalPETScIndex();
            const int lag_idx = node_idx->getLagrangianIndex();
            const double* FF = &FF_ghost_data_array[local_idx][0];
            double* F = &F_data_array[local_idx][0];
            Eigen::Map<IBTK::Vector> eig_F(F);
            Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> > eig_FF(FF);

            if (lag_idx >= left_begin && lag_idx <= left_end)
            {
                IBTK::Vector normal, trac;
                normal << -1.0, 0.0;

                get_trac(trac, eig_FF, normal);
                trac /= DX;
                eig_F += trac;
            }
            else if (lag_idx >= right_begin && lag_idx <= right_end)
            {
                IBTK::Vector normal, trac;
                normal << 1.0, 0.0;

                get_trac(trac, eig_FF, normal);
                trac /= DX;
                eig_F += trac;
            }
            else if (lag_idx >= top_begin && lag_idx <= top_end)
            {
                IBTK::Vector normal, trac;
                normal << 0.0, 1.0;

                get_trac(trac, eig_FF, normal);
                trac /= DX;
                eig_F += trac;
            }
        }
        FF_ghost_data->restoreArrays();
        F_data->restoreArrays();
    }*/

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

    // Restore arrays.
    B_ghost_data->restoreArrays();
    FF_ghost_data->restoreArrays();
    D_data->restoreArrays();
    D_ghost_data->restoreArrays();

    return;
} // computeLagrangianForce

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

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
        std::vector<std::vector<double> >& params = force_spec->getParameters();
        const unsigned int num_bonds = force_spec->getNumberOfSprings();
#if !defined(NDEBUG)
        TBOX_ASSERT(num_bonds == slv.size());
        TBOX_ASSERT(num_bonds == params.size());
#endif
        for (unsigned int k = 0; k < num_bonds; ++k)
        {
            lag_mastr_node_idxs[current_bond] = lag_idx;
            lag_slave_node_idxs[current_bond] = slv[k];
            petsc_mastr_node_idxs[current_bond] = petsc_idx;
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
IBPDForceGen::computeShapeTensor(SAMRAI::tbox::Pointer<IBTK::LData> B_data,
                                 SAMRAI::tbox::Pointer<IBTK::LData> X0_data,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > /*hierarchy*/,
                                 int level_number,
                                 double /*data_time*/,
                                 IBTK::LDataManager* /*l_data_manager*/)
{
    const int num_bonds = static_cast<int>(d_bond_data[level_number].lag_mastr_node_idxs.size());
    if (num_bonds == 0) return;

    const int* const petsc_mastr_node_idxs = &d_bond_data[level_number].petsc_mastr_node_idxs[0];
    const int* const petsc_slave_node_idxs = &d_bond_data[level_number].petsc_slave_node_idxs[0];
    double** const parameters = &d_bond_data[level_number].parameters[0];
    double* const B_node = B_data->getGhostedLocalFormVecArray()->data();
    const double* const X0_node = X0_data->getGhostedLocalFormVecArray()->data();

    static const int BLOCKSIZE = 16; // this parameter needs to be tuned
    int k, kblock, kunroll, X_mastr_idx, X_slave_idx, B_mastr_idx, B_slave_idx;
    double Q0[NDIM];
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

            Eigen::Map<IBTK::Vector> eig_Q0(Q0);
            //eig_Q0 = eig_Q0.cwiseProduct(eig_Q0);
            Eigen::Map<Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> > eig_B_mastr(&B_node[B_mastr_idx]),
                eig_B_slave(&B_node[B_slave_idx]);
            //eig_B_mastr += fail * vol_slave * eig_Q0.asDiagonal();
            //eig_B_slave += fail * vol_mastr * eig_Q0.asDiagonal();
             eig_B_mastr.noalias() += fail * vol_slave * eig_Q0 * eig_Q0.transpose();
             eig_B_slave.noalias() += fail * vol_mastr * eig_Q0 * eig_Q0.transpose();
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

        Eigen::Map<IBTK::Vector> eig_Q0(Q0);
        //eig_Q0 = eig_Q0.cwiseProduct(eig_Q0);
        Eigen::Map<Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> > eig_B_mastr(&B_node[B_mastr_idx]),
            eig_B_slave(&B_node[B_slave_idx]);
        //eig_B_mastr += fail * vol_slave * eig_Q0.asDiagonal();
        //eig_B_slave += fail * vol_mastr * eig_Q0.asDiagonal();
         eig_B_mastr.noalias() += fail * vol_slave * eig_Q0 * eig_Q0.transpose();
         eig_B_slave.noalias() += fail * vol_mastr * eig_Q0 * eig_Q0.transpose();
    }

    X0_data->restoreArrays();
    B_data->restoreArrays();

    return;
} // computeShapeTensor

void
IBPDForceGen::computeDeformationGradientTensor(SAMRAI::tbox::Pointer<IBTK::LData> FF_data,
                                               SAMRAI::tbox::Pointer<IBTK::LData> X_data,
                                               SAMRAI::tbox::Pointer<IBTK::LData> X0_data,
                                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > /*hierarchy*/,
                                               int level_number,
                                               double /*data_time*/,
                                               IBTK::LDataManager* /*l_data_manager*/)
{
    const int num_bonds = static_cast<int>(d_bond_data[level_number].lag_mastr_node_idxs.size());
    if (num_bonds == 0) return;

    const int* const petsc_mastr_node_idxs = &d_bond_data[level_number].petsc_mastr_node_idxs[0];
    const int* const petsc_slave_node_idxs = &d_bond_data[level_number].petsc_slave_node_idxs[0];
    double** const parameters = &d_bond_data[level_number].parameters[0];
    double* const FF_node = FF_data->getGhostedLocalFormVecArray()->data();
    const double* const X_node = X_data->getGhostedLocalFormVecArray()->data();
    const double* const X0_node = X0_data->getGhostedLocalFormVecArray()->data();

    static const int BLOCKSIZE = 16; // this parameter needs to be tuned
    int k, kblock, kunroll, X_mastr_idx, X_slave_idx, FF_mastr_idx, FF_slave_idx;
    double Q[NDIM], Q0[NDIM];
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

            Eigen::Map<const IBTK::Vector> eig_Q(Q), eig_Q0(Q0);
            Eigen::Map<Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> > eig_FF_mastr(&FF_node[FF_mastr_idx]),
                eig_FF_slave(&FF_node[FF_slave_idx]);
            eig_FF_mastr.noalias() += fail * vol_slave * eig_Q * eig_Q0.transpose();
            eig_FF_slave.noalias() += fail * vol_mastr * eig_Q * eig_Q0.transpose();
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

        Eigen::Map<const IBTK::Vector> eig_Q(Q), eig_Q0(Q0);
        Eigen::Map<Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> > eig_FF_mastr(&FF_node[FF_mastr_idx]),
            eig_FF_slave(&FF_node[FF_slave_idx]);
        eig_FF_mastr.noalias() += fail * vol_slave * eig_Q * eig_Q0.transpose();
        eig_FF_slave.noalias() += fail * vol_mastr * eig_Q * eig_Q0.transpose();
    }

    X_data->restoreArrays();
    X0_data->restoreArrays();
    FF_data->restoreArrays();

    return;
} // computeDeformationGradientTensor

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
    double* const F_node = F_data->getGhostedLocalFormVecArray()->data();
    double* const D_node = D_data->getGhostedLocalFormVecArray()->data();
    const double* const X_node = X_data->getGhostedLocalFormVecArray()->data();
    const double* const X0_node = X0_data->getGhostedLocalFormVecArray()->data();
    const double* const FF_node = FF_data->getGhostedLocalFormVecArray()->data();
    const double* const B_node = B_data->getGhostedLocalFormVecArray()->data();

    static const int BLOCKSIZE = 16; // this parameter needs to be tuned
    int k, kblock, kunroll, X_F_mastr_idx, X_F_slave_idx, FF_B_mastr_idx, FF_B_slave_idx, dmg_mastr_idx, dmg_slave_idx;
    double Q[NDIM], R;
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
            Eigen::Map<IBTK::Vector> eig_F_mastr(&F_node[X_F_mastr_idx]), eig_F_slave(&F_node[X_F_slave_idx]);
            Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> > eig_B_mastr(&B_node[FF_B_mastr_idx]),
                eig_B_slave(&B_node[FF_B_slave_idx]);
            Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> > eig_FF_mastr(
                &FF_node[FF_B_mastr_idx]),
                eig_FF_slave(&FF_node[FF_B_slave_idx]);

            Eigen::Vector4d dmg = pd_pk1force_damage_function(horizon,
                                                              delta,
                                                              R,
                                                              parameters[k],
                                                              eig_X0_mastr,
                                                              eig_X0_slave,
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
        Eigen::Map<IBTK::Vector> eig_F_mastr(&F_node[X_F_mastr_idx]), eig_F_slave(&F_node[X_F_slave_idx]);
        Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> > eig_B_mastr(&B_node[FF_B_mastr_idx]),
            eig_B_slave(&B_node[FF_B_slave_idx]);
        Eigen::Map<const Eigen::Matrix<double, NDIM, NDIM, Eigen::RowMajor> > eig_FF_mastr(&FF_node[FF_B_mastr_idx]),
            eig_FF_slave(&FF_node[FF_B_slave_idx]);

        Eigen::Vector4d dmg = pd_pk1force_damage_function(horizon,
                                                          delta,
                                                          R,
                                                          parameters[k],
                                                          eig_X0_mastr,
                                                          eig_X0_slave,
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

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
