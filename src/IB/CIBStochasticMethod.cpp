// Filename: CIBStochasticMethod.cpp
// Created on 26 Sept 2017 by Brennan Sprinkle
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith.
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

#include "ibamr/CIBStochasticMethod.h"
#include "ibamr/CIBMethod.h"
#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/MobilityFunctions.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/namespaces.h"
#include "ibtk/LSiloDataWriter.h"

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

CIBStochasticMethod::CIBStochasticMethod(const std::string& object_name,
                                         Pointer<Database> input_db,
                                         const int no_structures,
                                         bool register_for_restart)
    : CIBMethod(object_name, input_db, no_structures, register_for_restart)
{
    d_kT = 0.0;
    d_L_scale = 1.0;
    d_rf_delta = 1.0e-5;
    d_reject = false;
    d_num_reject = 0;
    return;

} // CIBMethod

CIBStochasticMethod::~CIBStochasticMethod()
{
    return;
} // ~CIBMethod

void
CIBStochasticMethod::setkT(double kT)
{
    d_kT = kT;
    return;
}

void
CIBStochasticMethod::setLScale(double scale)
{
    d_L_scale = scale;
    return;
}

void
CIBStochasticMethod::setRFdelta(double delta)
{
    d_rf_delta = delta;
    return;
}

double
CIBStochasticMethod::getkT() const
{
    return d_kT;
}

double
CIBStochasticMethod::getLScale() const
{
    return d_L_scale;
}

double
CIBStochasticMethod::getRFdelta() const
{
    return d_rf_delta;
}

void
CIBStochasticMethod::preprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    d_reject = false;
    CIBMethod::preprocessIntegrateData(current_time, new_time, num_cycles);
}

void
CIBStochasticMethod::midpointStep(double current_time, double new_time)
{
    if (!MathUtilities<double>::equalEps(current_time, d_ib_solver->getStartTime()))
    {
        checkLagUpdate();
    }
    if (d_reject)
    {
        CIBMethod::midpointStep(current_time, current_time);
        d_time_integrator_needs_regrid = false;
    }
    else
    {
        CIBMethod::midpointStep(current_time, new_time);
    }
}

void
CIBStochasticMethod::checkLagUpdate()
{
    int reject = 0;

    Eigen::Matrix3d rotation_mat_new = Eigen::Matrix3d::Zero();
    Eigen::Vector3d dr = Eigen::Vector3d::Zero();
    Eigen::Vector3d R_dr_new = Eigen::Vector3d::Zero();

    Eigen::Vector3d X_new = Eigen::Vector3d::Zero();

    // Get the grid extents.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    const double* const domain_x_lower = grid_geom->getXLower();
    const double* const domain_x_upper = grid_geom->getXUpper();
    double domain_length[NDIM];
    for (int d = 0; d < NDIM; ++d)
    {
        domain_length[d] = domain_x_upper[d] - domain_x_lower[d];
    }
    const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift();

    // Rotate the body with new rotational velocity about origin
    // and translate the body to newer position.
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        const boost::multi_array_ref<double, 2>& X0_array =
            *(d_l_data_manager->getLData("X0_unshifted", ln)->getLocalFormVecArray());
        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        // Get structures on this level.
        const std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(ln);
        const unsigned structs_on_this_ln = (unsigned)structIDs.size();
#if !defined(NDEBUG)
        TBOX_ASSERT(structs_on_this_ln == d_num_rigid_parts);
#endif

        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int lag_idx = node_idx->getLagrangianIndex();
            const int local_idx = node_idx->getLocalPETScIndex();
            const double* const X0 = &X0_array[local_idx][0];

            int struct_handle = 0;
            if (structs_on_this_ln > 1) struct_handle = getStructureHandle(lag_idx);

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dr[d] = X0[d] - d_center_of_mass_initial[struct_handle][d];
            }

            // Rotate dr vector using the rotation matrix.
            rotation_mat_new = d_quaternion_new[struct_handle].toRotationMatrix();
            R_dr_new = rotation_mat_new * dr;

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_new(d) = d_center_of_mass_new[struct_handle][d] + R_dr_new(d);
                if (!periodic_shift[d])
                {
                    if (X_new(d) < domain_x_lower[d] || domain_x_upper[d] < X_new(d))
                    {
                        reject++;
                        pout << "REJECTED BY WALL\n";
                    }
                }
            }
        }
        d_l_data_manager->getLData("X0_unshifted", ln)->restoreArrays();
    }

    reject = SAMRAI_MPI::sumReduction(reject);
    if (reject)
    {
        d_reject = true;
        d_num_reject++;
    }

    return;
} // checkLagUpdate

bool
CIBStochasticMethod::getReject()
{
    return d_reject;
}

int
CIBStochasticMethod::getNumReject()
{
    return d_num_reject;
}

void
CIBStochasticMethod::moveLagrangianData(double delta)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = delta;

    // Fill the rotation matrix of structures with rotation angle 0.5*(W^n)*dt.
    std::vector<Eigen::Matrix3d> rotation_mat(d_num_rigid_parts, Eigen::Matrix3d::Identity(3, 3));
    setRotationMatrix(d_rot_vel_half, d_quaternion_current, d_quaternion_half, rotation_mat, 0.5 * dt); // place to
                                                                                                        // check

    // Get the domain limits.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    const double* const domain_x_lower = grid_geom->getXLower();
    const double* const domain_x_upper = grid_geom->getXUpper();
    double domain_length[NDIM];
    for (int d = 0; d < NDIM; ++d)
    {
        domain_length[d] = domain_x_upper[d] - domain_x_lower[d];
    }
    const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift();

    // Rotate the body with current rotational velocity about origin
    // and translate the body to predicted position X^n+1/2.
    std::vector<Pointer<LData> >* X_half_data;
    bool* X_half_needs_ghost_fill;
    getPositionData(&X_half_data, &X_half_needs_ghost_fill, d_half_time);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        boost::multi_array_ref<double, 2>& X_half_array = *((*X_half_data)[ln]->getLocalFormVecArray());
        const boost::multi_array_ref<double, 2>& X0_array =
            *(d_l_data_manager->getLData("X0_unshifted", ln)->getLocalFormVecArray());
        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        // Get structures on this level.
        const std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(ln);
        const unsigned structs_on_this_ln = static_cast<unsigned>(structIDs.size());
#if !defined(NDEBUG)
        TBOX_ASSERT(structs_on_this_ln == d_num_rigid_parts);
#endif
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int lag_idx = node_idx->getLagrangianIndex();
            const int local_idx = node_idx->getLocalPETScIndex();
            double* const X_half = &X_half_array[local_idx][0];
            const double* const X0 = &X0_array[local_idx][0];
            Eigen::Vector3d dr = Eigen::Vector3d::Zero();
            Eigen::Vector3d R_dr = Eigen::Vector3d::Zero();

            int struct_handle = 0;
            if (structs_on_this_ln > 1) struct_handle = getStructureHandle(lag_idx);

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dr[d] = X0[d] - d_center_of_mass_initial[struct_handle][d];
            }

            // Rotate dr vector using the rotation matrix.
            R_dr = rotation_mat[struct_handle] * dr;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_half[d] = d_center_of_mass_current[struct_handle][d] + R_dr[d] +
                            0.5 * dt * d_trans_vel_half[struct_handle][d];

                if (periodic_shift[d])
                {
                    while (X_half[d] < domain_x_lower[d])
                    {
                        X_half[d] += domain_length[d];
                    }
                    while (X_half[d] >= domain_x_upper[d])
                    {
                        X_half[d] -= domain_length[d];
                    }
                }
            }
        }
        (*X_half_data)[ln]->restoreArrays();
        d_l_data_manager->getLData("X0_unshifted", ln)->restoreArrays();
    }
    *X_half_needs_ghost_fill = true;

    // Compute the COM at mid-step.
    for (unsigned struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            d_center_of_mass_half[struct_no][d] =
                d_center_of_mass_current[struct_no][d] + 0.5 * dt * d_trans_vel_half[struct_no][d];
            if (periodic_shift[d])
            {
                while (d_center_of_mass_half[struct_no][d] < domain_x_lower[d])
                {
                    d_center_of_mass_half[struct_no][d] += domain_length[d];
                }
                while (d_center_of_mass_half[struct_no][d] >= domain_x_upper[d])
                {
                    d_center_of_mass_half[struct_no][d] -= domain_length[d];
                }
            }
        }
    }
    return;
} // moveLagrangianData

void
CIBStochasticMethod::computeRFDforcesAndDisplacements()
{
    Vec F_rfd, U_rfd;
    double Vel_scale = (1.0 / d_kT);
    getNetExternalForceTorque(&F_rfd, d_current_time);
    VecDuplicate(F_rfd, &U_rfd);
    VecCopy(F_rfd, U_rfd);
    VecScale(U_rfd, Vel_scale);

    Vec U_all;
    VecScatter ctx;
    VecScatterCreateToAll(U_rfd, &ctx, &U_all);
    VecScatterBegin(ctx, U_rfd, U_all, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, U_rfd, U_all, INSERT_VALUES, SCATTER_FORWARD);
    const PetscScalar* U_array;
    VecGetArrayRead(U_all, &U_array);

    int part_free_dofs_begin = 0;
    for (unsigned part = 0; part < d_num_rigid_parts; ++part)
    {
        int num_free_dofs;
        const FreeRigidDOFVector& solve_dofs = getSolveRigidBodyVelocity(part, num_free_dofs);
        if (!num_free_dofs) continue;

        RDV UW;
        eigenToRDV(d_trans_vel_current[part], d_rot_vel_current[part], UW);
        const PetscScalar* a = &U_array[part_free_dofs_begin];
        for (int k = 0, p = 0; k < s_max_free_dofs; ++k)
        {
            if (solve_dofs[k])
            {
                UW[k] = a[p];
                ++p;
            }
        }
        rdvToEigen(UW, d_trans_vel_current[part], d_rot_vel_current[part]);
        d_trans_vel_current[part] *= (d_L_scale * d_L_scale);
        part_free_dofs_begin += num_free_dofs;
    }
    VecDestroy(&U_rfd);
    VecRestoreArrayRead(U_all, &U_array);
    VecScatterDestroy(&ctx);
    VecDestroy(&U_all);

    return;
} // computeRFDforcesAndDisplacements

void
CIBStochasticMethod::setHalfTimeVelocity(Vec U)
{
    // Get rigid DOFs on all processors.
    Vec U_all;
    VecScatter ctx;
    VecScatterCreateToAll(U, &ctx, &U_all);
    VecScatterBegin(ctx, U, U_all, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, U, U_all, INSERT_VALUES, SCATTER_FORWARD);
    const PetscScalar* U_array;
    VecGetArrayRead(U_all, &U_array);

    int part_free_dofs_begin = 0;
    for (unsigned part = 0; part < d_num_rigid_parts; ++part)
    {
        int num_free_dofs;
        const FreeRigidDOFVector& solve_dofs = getSolveRigidBodyVelocity(part, num_free_dofs);
        if (!num_free_dofs) continue;

        RDV UW;
        eigenToRDV(d_trans_vel_half[part], d_rot_vel_half[part], UW);
        const PetscScalar* a = &U_array[part_free_dofs_begin];
        for (int k = 0, p = 0; k < s_max_free_dofs; ++k)
        {
            if (solve_dofs[k])
            {
                UW[k] = a[p];
                ++p;
            }
        }
        rdvToEigen(UW, d_trans_vel_half[part], d_rot_vel_half[part]);

        part_free_dofs_begin += num_free_dofs;
    }

    VecRestoreArrayRead(U_all, &U_array);
    VecScatterDestroy(&ctx);
    VecDestroy(&U_all);

    return;
} // setHalfTimeVelocity (Brennan Sprinkle)

void
CIBStochasticMethod::resetRFDVelocity()
{
    d_trans_vel_half = d_trans_vel_current;
    d_rot_vel_half = d_rot_vel_current;
    return;
}

void
CIBStochasticMethod::getConfiguration(const unsigned int part, Eigen::Vector3d& C, Eigen::Quaterniond& Q)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(part < d_num_rigid_parts);
#endif
    C = d_center_of_mass_current[part];
    Q = d_quaternion_current[part];

    return;
}

void
CIBStochasticMethod::getInitialCOM(const unsigned int part, Eigen::Vector3d& C)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(part < d_num_rigid_parts);
#endif
    C = d_center_of_mass_initial[part];
    return;
}

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
