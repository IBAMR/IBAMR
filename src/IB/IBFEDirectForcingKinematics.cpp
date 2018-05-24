// Filename: IBFEDirectForcingKinematics.cpp
// Created on 24 May 2018 by Amneet Bhalla
//
// Copyright (c) 2002-2018, Amneet Bhalla
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

#include "ibamr/IBFEDirectForcingKinematics.h"
#include "ibamr/IBFEMethod.h"
#include "ibamr/namespaces.h"

namespace SAMRAI
{
namespace xfer
{
template <int DIM>
class RefineSchedule;
template <int DIM>
class CoarsenSchedule;
} // namespace xfer
} // namespace SAMRAI

using namespace libMesh;

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

IBFEDirectForcingKinematics::IBFEDirectForcingKinematics(const std::string& object_name,
                                                         Pointer<Database> input_db,
                                                         libMesh::Mesh* mesh,
                                                         IBTK::FEDataManager* fe_data_manager,
                                                         int level_number,
                                                         bool register_for_restart)
{
    // Set the object name and register it with the restart manager.
    d_object_name = object_name;
    d_registered_for_restart = false;
    if (register_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
        d_registered_for_restart = true;
    }

    d_mesh = mesh;
    d_level_number = level_number;
    d_fe_data_manager = fe_data_manager;

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    return;
} // IBFEDirectForcingKinematics

IBFEDirectForcingKinematics::~IBFEDirectForcingKinematics()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }
    return;
} // ~IBFEDirectForcingKinematics

FEDataManager*
IBFEDirectForcingKinematics::getFEDataManager() const
{
    return d_fe_data_manager;
} // getFEDataManager

void
IBFEMethod::forwardEulerStepdouble current_time, double new_time, PetscVector<double>& X_current_vec, PetscVector<double>& X_half_vec, PetscVector<double>& X_new_vec);
{
    const double dt = new_time - current_time;
    int ierr;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        ierr = VecWAXPY(d_X_new_vecs[part]->vec(), dt, d_U_current_vecs[part]->vec(), d_X_current_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(
            d_X_half_vecs[part]->vec(), 0.5, 0.5, 0.0, d_X_current_vecs[part]->vec(), d_X_new_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        d_X_new_vecs[part]->close();
        d_X_half_vecs[part]->close();
    }
    return;
} // eulerStep

void
IBFEDirectForcingKinematics::midpointStepdouble current_time, double new_time, PetscVector<double>& X_current_vec, PetscVector<double>& X_half_vec, PetscVector<double>& X_new_vec);
{
    const double dt = new_time - current_time;
    int ierr;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        ierr = VecWAXPY(d_X_new_vecs[part]->vec(), dt, d_U_half_vecs[part]->vec(), d_X_current_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(
            d_X_half_vecs[part]->vec(), 0.5, 0.5, 0.0, d_X_current_vecs[part]->vec(), d_X_new_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        d_X_new_vecs[part]->close();
        d_X_half_vecs[part]->close();
    }
    return;
} // midpointStep

void
IBFEDirectForcingKinematics::trapezoidalStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    int ierr;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        ierr =
            VecWAXPY(d_X_new_vecs[part]->vec(), 0.5 * dt, d_U_current_vecs[part]->vec(), d_X_current_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPY(d_X_new_vecs[part]->vec(), 0.5 * dt, d_U_new_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPBYPCZ(
            d_X_half_vecs[part]->vec(), 0.5, 0.5, 0.0, d_X_current_vecs[part]->vec(), d_X_new_vecs[part]->vec());
        IBTK_CHKERRQ(ierr);
        d_X_new_vecs[part]->close();
        d_X_half_vecs[part]->close();
    }
    return;
} // trapezoidalStep

void
IBFEDirectForcingKinematics::computeLagrangianForce(PetscVector<double>& F_vec,
                                                    PetscVector<double>& X_vec,
                                                    PetscVector<double>& U_vec,
                                                    const double data_time)
{
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
    bool all_imposed_dofs = false;

    for (int k = 0; k <= s_max_free_dofs; ++k)
    {
        all_imposed_dofs = (all_imposed_dofs || d_solve_rigid_vel[k]);
    }

    PetscVector<double>* Ub_vec = U_vec.clone().release();
    if (all_imposed_dofs)
    {
        computeImposedLagrangianForceDensity(F_vec, X_vec, U_vec, data_time);
    }
    else
    {
        computeFreeLagrangianForceDensity(F_vec, X_vec, U_vec, data_time);
    }

    return;
} // computeLagrangianForce

void
IBFEDirectForcingKinematics::putToDatabase(Pointer<Database> db)
{
    std::ostringstream U, W, C, Q;
    U << "U_" << d_object_name;
    W << "W_" << d_object_name;
    C << "C_" << d_object_name;
    Q << "Q_" << d_object_name;

    double Q_coeffs[4] = { d_quaternion_current[part].w(),
                           d_quaternion_current[part].x(),
                           d_quaternion_current[part].y(),
                           d_quaternion_current[part].z() };

    db->putDoubleArray(U.str(), &d_trans_vel_current[part][0], 3);
    db->putDoubleArray(W.str(), &d_rot_vel_current[part][0], 3);
    db->putDoubleArray(C.str(), &d_center_of_mass_current[part][0], 3);
    db->putDoubleArray(Q.str(), &Q_coeffs[0], 4);

    return;
} // putToDatabase

/////////////////////////////// PRIVATE ////////////////////////////////////

void
IBFEDirectForcingKinematics::computeImposedLagrangianForceDensity(PetscVector<double>& F_vec,
                                                                  PetscVector<double>& X_vec,
                                                                  PetscVector<double>& U_vec,
                                                                  const double data_time)
{
    const Eigen::Vector3d& X_com = d_center_of_mass_half;

    MeshBase& mesh = d_equation_systems->get_mesh();
    const unsigned int total_local_nodes = mesh.n_nodes_on_proc(SAMRAI_MPI::getRank());
    System& X_system = equation_systems->get_system<System>(IBFEMethod::COORDS_SYSTEM_NAME);
    System& U_system = equation_systems->get_system<System>(IBFEMethod::VELOCITY_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    const unsigned int U_sys_num = U_system.number();

    std::vector<std::vector<numeric_index_type> > nodal_indices(NDIM);
    std::vector<std::vector<double> > nodal_X_values(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        nodal_indices[d].reserve(total_local_nodes);
        nodal_values[d].reserve(total_local_nodes);
    }

    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
    {
        const Node* const n = *it;
        if (n->n_vars(X_sys_num) && n->n_vars(U_sys_num))
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM && n->n_vars(U_sys_num) == NDIM);
#endif
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                nodal_indices[d].push_back(n->dof_number(X_sys_num, d, 0));
            }
        }
    }

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        X_vec.get(nodal_indices[d], nodal_X_values[d]);
    }

    // Set the cross-product matrix
    Eigen::Matrix3d W(Eigen::Matrix3d::Zero());
#if (NDIM == 2)
    W(0, 1) = -d_rot_vel_half[2]; // -U[2];
    W(1, 0) = d_rot_vel_half[2];  // U[2];
#elif (NDIM == 3)
    W(0, 1) = -d_rot_vel_half[2]; //-U[5];
    W(1, 0) = d_rot_vel_half[2];  // U[5];
    W(0, 2) = d_rot_vel_half[1];  // U[4];
    W(2, 0) = -d_rot_vel_half[1]; // -U[4];
    W(1, 2) = -d_rot_vel_half[0]; //-U[3];
    W(2, 1) = d_rot_vel_half[0];  // U[3];
#endif
    Eigen::Vector3d R(Eigen::Vector3d::Zero());
    Eigen::Vector3d WxR(Eigen::Vector3d::Zero());
    for (unsigned int k = 0; k < total_local_nodes; ++k)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            R[d] = nodal_X_values[d][k] - X_com[d];
        }

        WxR = W * R;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            F_vec.set(nodal_indices[d][k], d_trans_vel_half[d] + WxR[d]);
        }
    }

    ierr = VecAXPY(F_vec.vec(), -1.0, U_vec.vec());
    ierr = VecScale(F_vec.vec(), d_rho / dt);
    F_vec.close();

    return;
} // ImposedLagrangianForceDensity

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace

//////////////////////////////////////////////////////////////////////////////
