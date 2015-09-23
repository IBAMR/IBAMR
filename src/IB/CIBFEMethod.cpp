// Filename: CIBFEMethod.cpp
// Created on 14 Oct 2014 by Amneet Bhalla
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

#include "ibamr/CIBFEMethod.h"
#include "ibamr/namespaces.h"
#include "ibamr/IBHierarchyIntegrator.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/LSiloDataWriter.h"
#include "ibtk/PETScMultiVec.h"
#include "ibtk/HierarchyMathOps.h"
#include "VisItDataWriter.h"
#include "libmesh/equation_systems.h"
#include "libmesh/system.h"
#include "libmesh/petsc_matrix.h"

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of CIBFEMethod restart file data.
static const int CIBFE_METHOD_VERSION = 1;
}

const std::string CIBFEMethod::CONSTRAINT_VELOCITY_SYSTEM_NAME = "IB constrained velocity system";

/////////////////////////////// PUBLIC ///////////////////////////////////////

CIBFEMethod::CIBFEMethod(const std::string& object_name,
                         Pointer<Database> input_db,
                         Mesh* mesh,
                         int max_level_number,
                         bool register_for_restart)
    : IBFEMethod(object_name, input_db, mesh, max_level_number, register_for_restart), CIBStrategy(1)
{
    commonConstructor(input_db);
    return;
} // CIBFEMethod

CIBFEMethod::CIBFEMethod(const std::string& object_name,
                         Pointer<Database> input_db,
                         const std::vector<Mesh*>& meshes,
                         int max_level_number,
                         bool register_for_restart)
    : IBFEMethod(object_name, input_db, meshes, max_level_number, register_for_restart),
      CIBStrategy(static_cast<unsigned>(meshes.size()))
{
    commonConstructor(input_db);
    return;
} // CIBFEMethod

CIBFEMethod::~CIBFEMethod()
{
    // IBFEMethod class is responsible for unregistering
    // the object and deleting the equation systems.

    return;
} // ~CIBFEMethod

void CIBFEMethod::registerConstrainedVelocityFunction(ConstrainedNodalVelocityFcnPtr nodalvelfcn,
                                                      ConstrainedCOMVelocityFcnPtr comvelfcn,
                                                      void* ctx,
                                                      unsigned int part)
{
    TBOX_ASSERT(part < d_num_rigid_parts);
    registerConstrainedVelocityFunction(ConstrainedVelocityFcnsData(nodalvelfcn, comvelfcn, ctx), part);
    return;
} // registerConstrainedVelocityFunction

void CIBFEMethod::registerConstrainedVelocityFunction(const ConstrainedVelocityFcnsData& data, unsigned int part)
{
    TBOX_ASSERT(part < d_num_rigid_parts);
    d_constrained_velocity_fcns_data[part] = data;
    return;
} // registerConstrainedVelocityFunction

void CIBFEMethod::registerExternalForceTorqueFunction(ExternalForceTorqueFcnPtr forcetorquefcn,
                                                      void* ctx,
                                                      unsigned int part)
{
    registerExternalForceTorqueFunction(ExternalForceTorqueFcnData(forcetorquefcn, ctx), part);
    return;
} // registerExternalForceTorqueFunction

void CIBFEMethod::registerExternalForceTorqueFunction(const ExternalForceTorqueFcnData& data, unsigned int part)
{

#if !defined(NDEBUG)
    TBOX_ASSERT(part < d_num_rigid_parts);
#endif

    d_ext_force_torque_fcn_data[part] = data;

    return;
} // registerExternalForceTorqueFunction

void CIBFEMethod::preprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    // Create most of the FE data vecs in the base class.
    IBFEMethod::preprocessIntegrateData(current_time, new_time, num_cycles);

    // Create vecs for constraint force and velocity.
    d_F_current_vecs.resize(d_num_rigid_parts);
    d_F_new_vecs.resize(d_num_rigid_parts);
    d_U_constrained_systems.resize(d_num_rigid_parts);
    d_U_constrained_current_vecs.resize(d_num_rigid_parts);
    d_U_constrained_half_vecs.resize(d_num_rigid_parts);

    // PETSc wrappers.
    d_vL_current.resize(d_num_rigid_parts, PETSC_NULL);
    d_vL_new.resize(d_num_rigid_parts, PETSC_NULL);

    for (unsigned int part = 0; part < d_num_rigid_parts; ++part)
    {
        d_F_current_vecs[part] = dynamic_cast<PetscVector<double>*>(d_F_systems[part]->current_local_solution.get());
        d_F_new_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_F_current_vecs[part]->clone().release()); // WARNING: must be manually deleted

        d_vL_current[part] = d_F_current_vecs[part]->vec();
        d_vL_new[part] = d_F_new_vecs[part]->vec();

        d_U_constrained_systems[part] = &d_equation_systems[part]->get_system(CONSTRAINT_VELOCITY_SYSTEM_NAME);
        d_U_constrained_current_vecs[part] =
            dynamic_cast<PetscVector<double>*>(d_U_constrained_systems[part]->current_local_solution.get());
        d_U_constrained_half_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_U_constrained_current_vecs[part]->clone().release()); // WARNING: must be manually deleted

        // Initialize F^{n+1} to F^{n}, and U_k^{n+1/2} to equal U_k^{n}.
        d_F_current_vecs[part]->localize(*d_F_new_vecs[part]);
        d_U_constrained_current_vecs[part]->localize(*d_U_constrained_half_vecs[part]);
    }

    // Create data structures for Lagrange multiplier
    VecCreateMultiVec(PETSC_COMM_WORLD, d_num_rigid_parts, &d_vL_current[0], &d_mv_L_current);
    VecCreateMultiVec(PETSC_COMM_WORLD, d_num_rigid_parts, &d_vL_new[0], &d_mv_L_new);

    // Get data for free and prescribed bodies.
    int free_parts = 0;
    for (unsigned int part = 0; part < d_num_rigid_parts; ++part)
    {
        int num_free_dofs;
        const FRDV& solve_dofs = getSolveRigidBodyVelocity(part, num_free_dofs);
        const bool prescribed_velocity = num_free_dofs < s_max_free_dofs;
        if (prescribed_velocity)
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(d_constrained_velocity_fcns_data[part].comvelfcn);
#endif
            Eigen::Vector3d trans_vel_current, trans_vel_half, trans_vel_new, rot_vel_current, rot_vel_half,
                rot_vel_new;

            d_constrained_velocity_fcns_data[part].comvelfcn(d_current_time, trans_vel_current, rot_vel_current);
            d_constrained_velocity_fcns_data[part].comvelfcn(d_half_time, trans_vel_half, rot_vel_half);
            d_constrained_velocity_fcns_data[part].comvelfcn(d_new_time, trans_vel_new, rot_vel_new);

            // Update only prescribed velocities in the internal data structure.
            for (int d = 0; d < NDIM; ++d)
            {
                if (!solve_dofs[d])
                {
                    d_trans_vel_current[part][d] = trans_vel_current[d];
                    d_trans_vel_half[part][d] = trans_vel_half[d];
                    d_trans_vel_new[part][d] = trans_vel_new[d];
                }
            }
#if (NDIM == 2)
            if (!solve_dofs[2])
            {
                d_rot_vel_current[part][2] = rot_vel_current[2];
                d_rot_vel_half[part][2] = rot_vel_half[2];
                d_rot_vel_new[part][2] = rot_vel_new[2];
            }
#elif(NDIM == 3)
            for (int d = 0; d < NDIM; ++d)
            {
                if (!solve_dofs[3 + d])
                {
                    d_rot_vel_current[part][d] = rot_vel_current[d];
                    d_rot_vel_half[part][d] = rot_vel_half[d];
                    d_rot_vel_new[part][d] = rot_vel_new[d];
                }
            }
#endif
        }

        if (num_free_dofs)
        {
            d_U.push_back(static_cast<Vec>(NULL));
            d_F.push_back(static_cast<Vec>(NULL));

            Vec& U = d_U.back();
            Vec& F = d_F.back();

            PetscInt n = 0, N = num_free_dofs;
            if (!SAMRAI_MPI::getRank()) n = N;
            VecCreateMPI(PETSC_COMM_WORLD, n, N, &U);
            VecCreateMPI(PETSC_COMM_WORLD, n, N, &F);

            // Initialize U and F.
            // For U we use current timestep value as a guess.

            RDV Fr;
            Eigen::Vector3d F_ext, T_ext;
            if (d_ext_force_torque_fcn_data[part].forcetorquefcn)
            {
                d_ext_force_torque_fcn_data[part].forcetorquefcn(d_new_time, F_ext, T_ext);
            }
            else
            {
                F_ext.setZero();
                T_ext.setZero();
            }
            eigenToRDV(F_ext, T_ext, Fr);

            RDV Ur;
            eigenToRDV(d_trans_vel_current[part], d_rot_vel_current[part], Ur);

            for (int k = 0, p = 0; k < s_max_free_dofs; ++k)
            {
                if (solve_dofs[k])
                {
                    VecSetValue(U, p, Ur[k], INSERT_VALUES);
                    VecSetValue(F, p, Fr[k], INSERT_VALUES);
                    ++p;
                }
            }
            VecAssemblyBegin(U);
            VecAssemblyBegin(F);

            VecAssemblyEnd(U);
            VecAssemblyEnd(F);

            ++free_parts;
        }
    }
    VecCreateMultiVec(PETSC_COMM_WORLD, free_parts, &d_U[0], &d_mv_U);
    VecCreateMultiVec(PETSC_COMM_WORLD, free_parts, &d_F[0], &d_mv_F);

    return;
} // preprocessIntegrateData

void CIBFEMethod::postprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    // Clean the temporary FE data vecs.
    IBFEMethod::postprocessIntegrateData(current_time, new_time, num_cycles);

    // Destroy the free-parts PETSc Vecs.
    for (size_t k = 0; k < d_U.size(); ++k)
    {
        VecDestroy(&d_U[k]);
        VecDestroy(&d_F[k]);
    }
    VecDestroy(&d_mv_U);
    VecDestroy(&d_mv_F);

    for (unsigned part = 0; part < d_num_rigid_parts; ++part)
    {
        // Reset time-dependent Lagrangian data.
        *d_F_systems[part]->solution = *d_F_new_vecs[part];
        *d_U_constrained_systems[part]->solution = *d_U_constrained_half_vecs[part];

        // Deallocate Lagrangian scratch data.
        delete d_F_new_vecs[part];
        delete d_U_constrained_half_vecs[part];
    }

    d_F_current_vecs.clear();
    d_F_new_vecs.clear();
    d_vL_current.clear();
    d_vL_new.clear();
    VecDestroy(&d_mv_L_current);
    VecDestroy(&d_mv_L_new);

    d_U_constrained_systems.clear();
    d_U_constrained_current_vecs.clear();
    d_U_constrained_half_vecs.clear();

    // New center of mass translational and rotational velocity becomes
    // current velocity for the next time step.
    d_trans_vel_current = d_trans_vel_new;
    d_rot_vel_current = d_rot_vel_new;

    return;
} // postprocessIntegrateData

bool CIBFEMethod::setComputeVelL2Projection(const bool compute_L2_projection)
{
    bool cached = d_compute_L2_projection;
    d_compute_L2_projection = compute_L2_projection;

    return cached;
} // setVelocityL2Projection

void CIBFEMethod::interpolateVelocity(const int u_data_idx,
                                      const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*u_synch_scheds*/,
                                      const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                      const double data_time)
{
    if (d_lag_velvec_is_initialized)
    {
        for (unsigned int part = 0; part < d_num_rigid_parts; ++part)
        {
            NumericVector<double>* X_vec = d_X_half_vecs[part];
            NumericVector<double>* X_ghost_vec = d_X_IB_ghost_vecs[part];
            NumericVector<double>* U_vec = d_U_half_vecs[part];
            X_vec->localize(*X_ghost_vec);
            if (d_use_IB_interp_operator && d_compute_L2_projection)
            {
                d_fe_data_managers[part]->interp(u_data_idx, *U_vec, *X_ghost_vec, VELOCITY_SYSTEM_NAME,
                                                 u_ghost_fill_scheds, data_time);
            }
            else if (d_use_IB_interp_operator && !d_compute_L2_projection)
            {
                d_fe_data_managers[part]->interpWeighted(u_data_idx, *U_vec, *X_ghost_vec, VELOCITY_SYSTEM_NAME,
                                                         u_ghost_fill_scheds, data_time);
            }
        }
        d_lag_velvec_is_initialized = false;
    }
    return;
} // interpolateVelocity

void CIBFEMethod::eulerStep(const double current_time, const double new_time)
{   
    std::vector<Eigen::Matrix3d> moment_of_inertia_current, moment_of_inertia_half;
    moment_of_inertia_current.resize(d_num_rigid_parts, Eigen::Matrix3d::Zero());
    moment_of_inertia_half.resize(d_num_rigid_parts, Eigen::Matrix3d::Zero());
    
    const double dt = new_time - current_time;

    // Compute center of mass and moment of inertia of the body at t^n.
    for (unsigned part = 0; part < d_num_rigid_parts; ++part)
    {
        computeCOMandMOI(part, d_center_of_mass_current[part], moment_of_inertia_current[part],
                         d_X_current_vecs[part]);
    }

    // Fill the rotation matrix of structures with rotation angle 0.5*(W^n)*dt.
    std::vector<Eigen::Matrix3d> rotation_mat(d_num_rigid_parts, Eigen::Matrix3d::Zero());
    for (unsigned int part = 0; part < d_num_rigid_parts; ++part)
    {
        for (int i = 0; i < 3; ++i)
        {
            rotation_mat[part](i, i) = 1.0;
        }
        setRotationMatrix(d_rot_vel_current[part], rotation_mat[part], 0.5 * dt);
    }

    // Rotate the body with current rotational velocity about center of mass
    // and translate the body to predicted position X^n+1/2.
    Eigen::Vector3d dr = Eigen::Vector3d::Zero();
    Eigen::Vector3d Rxdr = Eigen::Vector3d::Zero();
    for (unsigned int part = 0; part < d_num_rigid_parts; ++part)
    {
        EquationSystems* equation_systems = d_equation_systems[part];
        MeshBase& mesh = equation_systems->get_mesh();
        const unsigned int total_local_nodes = mesh.n_nodes_on_proc(SAMRAI_MPI::getRank());
        System& X_system = *d_X_systems[part];
        const unsigned int X_sys_num = X_system.number();
        PetscVector<double>& X_current = *d_X_current_vecs[part];
        PetscVector<double>& X_half = *d_X_half_vecs[part];

        std::vector<std::vector<numeric_index_type> > nodal_X_indices(NDIM);
        std::vector<std::vector<double> > nodal_X_values(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            nodal_X_indices[d].reserve(total_local_nodes);
            nodal_X_values[d].reserve(total_local_nodes);
        }

        for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
        {
            const Node* const n = *it;
            if (n->n_vars(X_sys_num))
            {
#if !defined(NDEBUG)
                TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM);
#endif
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    nodal_X_indices[d].push_back(n->dof_number(X_sys_num, d, 0));
                }
            }
        }

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_current.get(nodal_X_indices[d], &nodal_X_values[d][0]);
        }

        for (unsigned int k = 0; k < total_local_nodes; ++k)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dr[d] = nodal_X_values[d][k] - d_center_of_mass_current[part][d];
            }

            // Rotate dr vector using the rotation matrix.
            Rxdr = rotation_mat[part] * dr;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_half.set(nodal_X_indices[d][k],
                           d_center_of_mass_current[part][d] + Rxdr[d] + 0.5 * dt * d_trans_vel_current[part][d]);
            }
        }
        X_half.close();
    }

    // Compute the COM and MOI at mid-point.
    for (unsigned part = 0; part < d_num_rigid_parts; ++part)
    {
        computeCOMandMOI(part, d_center_of_mass_half[part], moment_of_inertia_half[part], d_X_half_vecs[part]);
    }

    return;
} // eulerStep

void CIBFEMethod::midpointStep(const double current_time, const double new_time)
{
    const double dt = new_time - current_time;
    // Fill the rotation matrix of structures with rotation angle (W^n+1)*dt.
    std::vector<Eigen::Matrix3d> rotation_mat(d_num_rigid_parts, Eigen::Matrix3d::Zero());
    for (unsigned int part = 0; part < d_num_rigid_parts; ++part)
    {
        for (int i = 0; i < 3; ++i)
        {
            rotation_mat[part](i, i) = 1.0;
        }
        setRotationMatrix(d_rot_vel_half[part], rotation_mat[part], dt);
    }

    // Rotate the body with current rotational velocity about origin
    // and translate the body to predicted position X^n+1/2.
    Eigen::Vector3d dr = Eigen::Vector3d::Zero();
    Eigen::Vector3d Rxdr = Eigen::Vector3d::Zero();
    for (unsigned int part = 0; part < d_num_rigid_parts; ++part)
    {
        EquationSystems* equation_systems = d_equation_systems[part];
        MeshBase& mesh = equation_systems->get_mesh();
        const unsigned int total_local_nodes = mesh.n_nodes_on_proc(SAMRAI_MPI::getRank());
        System& X_system = *d_X_systems[part];
        const unsigned int X_sys_num = X_system.number();
        PetscVector<double>& X_current = *d_X_current_vecs[part];
        PetscVector<double>& X_new = *d_X_new_vecs[part];

        std::vector<std::vector<numeric_index_type> > nodal_X_indices(NDIM);
        std::vector<std::vector<double> > nodal_X_values(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            nodal_X_indices[d].reserve(total_local_nodes);
            nodal_X_values[d].reserve(total_local_nodes);
        }

        for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
        {
            const Node* const n = *it;
            if (n->n_vars(X_sys_num))
            {
#if !defined(NDEBUG)
                TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM);
#endif
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    nodal_X_indices[d].push_back(n->dof_number(X_sys_num, d, 0));
                }
            }
        }

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_current.get(nodal_X_indices[d], &nodal_X_values[d][0]);
        }

        for (unsigned int k = 0; k < total_local_nodes; ++k)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dr[d] = nodal_X_values[d][k] - d_center_of_mass_current[part][d];
            }

            // Rotate dr vector using the rotation matrix.
            Rxdr = rotation_mat[part] * dr;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_new.set(nodal_X_indices[d][k],
                          d_center_of_mass_current[part][d] + Rxdr[d] + dt * d_trans_vel_half[part][d]);
            }
        }
        X_new.close();
    }

    return;
} // midpointStep

void CIBFEMethod::trapezoidalStep(const double /*current_time*/, const double /*new_time*/)
{
    TBOX_ERROR("CIBFEMethod does not support trapezoidal time-stepping rule for position update."
               << " Only midpoint rule is supported for position update." << std::endl);

    return;
} // trapezoidalStep

void CIBFEMethod::computeLagrangianForce(double /*data_time*/)
{
    // intentionally left blank
    return;
} // computeLagrangianForce

void CIBFEMethod::spreadForce(
    int f_data_idx,
    RobinPhysBdryPatchStrategy* f_phys_bdry_op,
    const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds,
    double data_time)
{
    if (d_constraint_force_is_initialized)
    {
        IBFEMethod::spreadForce(f_data_idx, f_phys_bdry_op, f_prolongation_scheds, data_time);
        d_constraint_force_is_initialized = false;
    }
    return;

} // spreadForce

void CIBFEMethod::setConstraintForce(Vec L, const double data_time, const double scale)
{

#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
#endif

    // Unpack the Lambda vector.
    Vec* vL;
    VecMultiVecGetSubVecs(L, &vL);

#if !defined(NDEBUG)
    PetscInt size;
    VecMultiVecGetNumberOfSubVecs(L, &size);
    TBOX_ASSERT(size == (int)d_num_rigid_parts);
#endif

    for (unsigned part = 0; part < d_num_rigid_parts; ++part)
    {
        PetscVector<double>* F_vec = d_F_half_vecs[part];
        VecCopy(vL[part], F_vec->vec());
        VecScale(F_vec->vec(), scale);
    }
    d_constraint_force_is_initialized = true;

    return;
} // setConstraintForce

void CIBFEMethod::getConstraintForce(Vec* L, const double data_time)
{

    if (MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        *L = d_mv_L_current;
    }
    else if (MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        *L = d_mv_L_new;
    }
    else
    {
        plog << "Warning CIBFEMethod::getConstraintForce() : constraint force "
             << "enquired at some other time than current or new time.";
        return;
    }

    return;
} // getConstraintForce

void CIBFEMethod::getFreeRigidVelocities(Vec* U, const double data_time)
{

#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_current_time) ||
                MathUtilities<double>::equalEps(data_time, d_new_time));
#endif

    *U = d_mv_U;
    return;

} // getFreeRigidVelocities

void CIBFEMethod::getNetExternalForceTorque(Vec* F, const double data_time)
{

#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_current_time) ||
                MathUtilities<double>::equalEps(data_time, d_new_time));
#endif

    *F = d_mv_F;
    return;

} // getNetExternalForceTorque

void CIBFEMethod::subtractMeanConstraintForce(Vec L, int f_data_idx, const double scale)
{
    // Unpack the Lambda vector.
    Vec* vL;
    VecMultiVecGetSubVecs(L, &vL);

#if !defined(NDEBUG)
    PetscInt size;
    VecMultiVecGetNumberOfSubVecs(L, &size);
    TBOX_ASSERT(size == (int)d_num_rigid_parts);
#endif

    double F[NDIM] = { 0.0 };
    for (unsigned part = 0; part < d_num_rigid_parts; ++part)
    {
        // Copy L to internal data structure.
        PetscVector<double>& L_petsc = *d_F_half_vecs[part];
        Vec L_global_vec = L_petsc.vec();
        VecCopy(vL[part], L_global_vec);
        VecScale(L_global_vec, scale);
        L_petsc.close();
        Vec L_local_ghost_vec;
        VecGhostGetLocalForm(L_global_vec, &L_local_ghost_vec);
        double* L_local_ghost_soln;
        VecGetArray(L_local_ghost_vec, &L_local_ghost_soln);

        // Build quadrature rule.
        EquationSystems* equation_systems = d_equation_systems[part];
        MeshBase& mesh = equation_systems->get_mesh();
        const unsigned int dim = mesh.mesh_dimension();
        AutoPtr<QBase> qrule = QBase::build(d_quad_type, dim, d_quad_order);

        // Extract the FE system and DOF map, and setup the FE object.
        System& L_system = *d_F_systems[part];
        DofMap& L_dof_map = L_system.get_dof_map();
        std::vector<std::vector<unsigned int> > L_dof_indices(NDIM);
        FEType L_fe_type = L_dof_map.variable_type(0);
        AutoPtr<FEBase> L_fe_autoptr(FEBase::build(dim, L_fe_type));
        FEBase* L_fe = L_fe_autoptr.get();
        L_fe->attach_quadrature_rule(qrule.get());
        const std::vector<double>& JxW = L_fe->get_JxW();
        const std::vector<std::vector<double> >& phi = L_fe->get_phi();

        double L_qp[NDIM];
        boost::multi_array<double, 2> L_node;
        const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            const Elem* const elem = *el_it;
            L_fe->reinit(elem);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                L_dof_map.dof_indices(elem, L_dof_indices[d], d);
            }
            get_values_for_interpolation(L_node, L_petsc, L_local_ghost_soln, L_dof_indices);

            const unsigned int n_qp = qrule->n_points();
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {
                interpolate(L_qp, qp, L_node, phi);

                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    F[d] += L_qp[d] * JxW[qp];
                }
            }
        }
    }
    SAMRAI_MPI::sumReduction(F, NDIM);

    // Subtract the mean from Eulerian body force
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double vol_domain = getHierarchyMathOps()->getVolumeOfPhysicalDomain();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM, double> > p_data = patch->getPatchData(f_data_idx);
            const Box<NDIM>& patch_box = patch->getBox();
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> it(patch_box, axis); it; it++)
                {
                    (*p_data)(it()) -= F[axis] / vol_domain;
                }
            }
        }
    }

    return;
} // subtractMeanConstraintForce

void CIBFEMethod::setInterpolatedVelocityVector(Vec /*V*/, const double data_time)
{

#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
#endif
    d_lag_velvec_is_initialized = true;
    return;
} // setInterpolatedVelocityVector

void CIBFEMethod::getInterpolatedVelocity(Vec V, const double data_time, const double scale)
{

#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
#endif

    // Unpack the velocity vector.
    Vec* vV;
    VecMultiVecGetSubVecs(V, &vV);
    for (unsigned int part = 0; part < d_num_rigid_parts; ++part)
    {
        PetscVector<double>* U_vec = d_U_half_vecs[part];
        VecCopy(U_vec->vec(), vV[part]);
        VecScale(vV[part], scale);
    }
    return;

} // getInterpolatedVelocity

void CIBFEMethod::computeMobilityRegularization(Vec D, Vec L, const double scale)
{

    if (!d_compute_L2_projection)
    {
        Vec* vL, *vD;
        VecMultiVecGetSubVecs(L, &vL);
        VecMultiVecGetSubVecs(D, &vD);

        for (unsigned part = 0; part < d_num_rigid_parts; ++part)
        {
            std::pair<LinearSolver<double>*, SparseMatrix<double>*> filter =
                d_fe_data_managers[part]->buildL2ProjectionSolver(VELOCITY_SYSTEM_NAME);
            Mat mat = static_cast<PetscMatrix<double>*>(filter.second)->mat();
            MatMult(mat, vL[part], vD[part]);
        }
        VecScale(D, scale);
    }
    else
    {
        VecCopy(L, D);
        VecScale(D, scale);
    }

    return;

} // computeMobilityRegularization

unsigned int CIBFEMethod::getNumberOfNodes(const unsigned int part) const
{
    return d_num_nodes[part];

} // getNumberOfNodes

void CIBFEMethod::computeNetRigidGeneralizedForce(const unsigned int part, Vec L, RigidDOFVector& F)
{
    EquationSystems* equation_systems = d_equation_systems[part];
    MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    AutoPtr<QBase> qrule = QBase::build(d_quad_type, dim, d_quad_order);

    // Extract the FE system and DOF map, and setup the FE object.
    System& L_system = *d_F_systems[part];
    System& X_system = *d_X_systems[part];
    DofMap& L_dof_map = L_system.get_dof_map();
    DofMap& X_dof_map = X_system.get_dof_map();
    std::vector<std::vector<unsigned int> > L_dof_indices(NDIM);
    std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
    FEType L_fe_type = L_dof_map.variable_type(0);
    FEType X_fe_type = X_dof_map.variable_type(0);
    AutoPtr<FEBase> L_fe_autoptr(FEBase::build(dim, L_fe_type)), X_fe_autoptr(NULL);
    if (L_fe_type != X_fe_type)
    {
        X_fe_autoptr = AutoPtr<FEBase>(FEBase::build(dim, X_fe_type));
    }
    FEBase* L_fe = L_fe_autoptr.get();
    FEBase* X_fe = X_fe_autoptr.get() ? X_fe_autoptr.get() : L_fe_autoptr.get();
    const std::vector<double>& JxW_L = L_fe->get_JxW();
    const std::vector<std::vector<double> >& phi_L = L_fe->get_phi();
    const std::vector<std::vector<double> >& phi_X = X_fe->get_phi();

    const Eigen::Vector3d& X_com = d_center_of_mass_half[part];
    libMesh::PetscVector<double>& X_petsc = *d_X_half_vecs[part];
    if (!X_petsc.closed()) X_petsc.close();
    Vec X_global_vec = X_petsc.vec();
    Vec X_local_ghost_vec;
    VecGhostGetLocalForm(X_global_vec, &X_local_ghost_vec);
    double* X_local_ghost_soln;
    VecGetArray(X_local_ghost_vec, &X_local_ghost_soln);

    libMesh::PetscVector<double>& L_petsc = *d_F_half_vecs[part];
    Vec L_global_vec = L_petsc.vec();
    VecCopy(L, L_global_vec);
    L_petsc.close(); // Sync ghost values.
    Vec L_local_ghost_vec;
    VecGhostGetLocalForm(L_global_vec, &L_local_ghost_vec);
    double* L_local_ghost_soln;
    VecGetArray(L_local_ghost_vec, &L_local_ghost_soln);

    F.setZero();
    boost::multi_array<double, 2> X_node, L_node;
    double X_qp[NDIM], L_qp[NDIM];
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        L_fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_dof_map.dof_indices(elem, X_dof_indices[d], d);
            L_dof_map.dof_indices(elem, L_dof_indices[d], d);
        }
        get_values_for_interpolation(L_node, L_petsc, L_local_ghost_soln, L_dof_indices);
        get_values_for_interpolation(X_node, X_petsc, X_local_ghost_soln, X_dof_indices);

        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(X_qp, qp, X_node, phi_X);
            interpolate(L_qp, qp, L_node, phi_L);

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                F[d] += L_qp[d] * JxW_L[qp];
            }
#if (NDIM == 2)
            F[NDIM] += (L_qp[1] * (X_qp[0] - X_com[0]) - L_qp[0] * (X_qp[1] - X_com[1])) * JxW_L[qp];
#elif(NDIM == 3)
            F[NDIM] += (L_qp[2] * (X_qp[1] - X_com[1]) - L_qp[1] * (X_qp[2] - X_com[2])) * JxW_L[qp];
            F[NDIM + 1] += (L_qp[0] * (X_qp[2] - X_com[2]) - L_qp[2] * (X_qp[0] - X_com[0])) * JxW_L[qp];
            F[NDIM + 2] += (L_qp[1] * (X_qp[0] - X_com[0]) - L_qp[0] * (X_qp[1] - X_com[1])) * JxW_L[qp];
#endif
        }
    }
    SAMRAI_MPI::sumReduction(&F[0], NDIM * (NDIM + 1) / 2);
    return;

} // computeNetRigidGeneralizedForce

void CIBFEMethod::copyVecToArray(Vec b,
                                 double* array,
                                 const std::vector<unsigned>& struct_ids,
                                 const int data_depth,
                                 const int array_rank)
{
    if (struct_ids.empty()) return;

    // Wrap the raw data in a PETSc Vec
    PetscInt total_nodes = 0;
    for (unsigned k = 0; k < struct_ids.size(); ++k)
    {
        total_nodes += getNumberOfNodes(struct_ids[k]);
    }
    PetscInt size = total_nodes * data_depth;
    int rank = SAMRAI_MPI::getRank();
    PetscInt array_local_size = 0;
    if (rank == array_rank) array_local_size = size;
    Vec array_vec;
    VecCreateMPIWithArray(PETSC_COMM_WORLD, /*blocksize*/ 1, array_local_size, PETSC_DECIDE, array, &array_vec);

    // Get the components of PETScMultiVec.
    Vec* vb;
    VecMultiVecGetSubVecs(b, &vb);

    // Scatter values
    PetscInt offset = 0;
    for (unsigned k = 0; k < struct_ids.size(); ++k)
    {
        const unsigned struct_id = struct_ids[k];
        PetscInt nodes = getNumberOfNodes(struct_id);
        PetscInt size_vec = nodes * data_depth;

        IS is_vec;
        IS is_array;
        PetscInt step = 1;

        ISCreateStride(PETSC_COMM_SELF, size_vec, 0, step, &is_vec);
        ISCreateStride(PETSC_COMM_SELF, size_vec, offset, step, &is_array);

        VecScatter ctx;
        VecScatterCreate(vb[struct_id], is_vec, array_vec, is_array, &ctx);

        VecScatterBegin(ctx, vb[struct_id], array_vec, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(ctx, vb[struct_id], array_vec, INSERT_VALUES, SCATTER_FORWARD);

        VecScatterDestroy(&ctx);
        ISDestroy(&is_vec);
        ISDestroy(&is_array);

        offset += size_vec;
    }

    // Destroy the wrapped Vec
    VecDestroy(&array_vec);

    return;
} // copyVecToArray

void CIBFEMethod::copyArrayToVec(Vec b,
                                 double* array,
                                 const std::vector<unsigned>& struct_ids,
                                 const int data_depth,
                                 const int array_rank)
{
    if (struct_ids.empty()) return;

    // Wrap the raw data in a PETSc Vec
    PetscInt total_nodes = 0;
    for (unsigned k = 0; k < struct_ids.size(); ++k)
    {
        total_nodes += getNumberOfNodes(struct_ids[k]);
    }
    PetscInt size = total_nodes * data_depth;
    int rank = SAMRAI_MPI::getRank();
    PetscInt array_local_size = 0;
    if (rank == array_rank) array_local_size = size;
    Vec array_vec;
    VecCreateMPIWithArray(PETSC_COMM_WORLD, /*blocksize*/ 1, array_local_size, PETSC_DECIDE, array, &array_vec);

    // Get the components of PETScMultiVec.
    Vec* vb;
    VecMultiVecGetSubVecs(b, &vb);

    // Scatter values
    PetscInt offset = 0;
    for (unsigned k = 0; k < struct_ids.size(); ++k)
    {
        const unsigned struct_id = struct_ids[k];
        PetscInt nodes = getNumberOfNodes(struct_id);
        PetscInt size_vec = nodes * data_depth;

        IS is_vec;
        IS is_array;
        PetscInt step = 1;

        ISCreateStride(PETSC_COMM_SELF, size_vec, 0, step, &is_vec);
        ISCreateStride(PETSC_COMM_SELF, size_vec, offset, step, &is_array);

        VecScatter ctx;
        VecScatterCreate(array_vec, is_array, vb[struct_id], is_vec, &ctx);

        VecScatterBegin(ctx, array_vec, vb[struct_id], INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(ctx, array_vec, vb[struct_id], INSERT_VALUES, SCATTER_FORWARD);

        VecScatterDestroy(&ctx);
        ISDestroy(&is_vec);
        ISDestroy(&is_array);

        offset += size_vec;
    }

    // Destroy the wrapped Vec
    VecDestroy(&array_vec);

    return;
} // copyArrayToVec

void CIBFEMethod::setRigidBodyVelocity(const unsigned int part, const RigidDOFVector& U, Vec V)
{
    PetscVector<double>& U_k = *d_U_constrained_half_vecs[part];
    PetscVector<double>& X_half = *d_X_half_vecs[part];
    void* ctx = d_constrained_velocity_fcns_data[part].ctx;
    d_constrained_velocity_fcns_data[part].nodalvelfcn(U_k, U, X_half, d_center_of_mass_half[part],
                                                       d_equation_systems[part], d_new_time, ctx);
    if (!U_k.closed()) U_k.close();

    // We filter the rigid body velocity using the basis function of the
    // deformational field, in the case when L2-projection is not performed.
    Vec* vV;
    VecMultiVecGetSubVecs(V, &vV);
    if (!d_compute_L2_projection)
    {
        std::pair<LinearSolver<double>*, SparseMatrix<double>*> filter =
            d_fe_data_managers[part]->buildL2ProjectionSolver(VELOCITY_SYSTEM_NAME);
        Mat mat = static_cast<PetscMatrix<double>*>(filter.second)->mat();
        MatMult(mat, U_k.vec(), vV[part]);
    }
    else
    {
        VecCopy(U_k.vec(), vV[part]);
    }

    return;
} // setRigidBodyVelocity

void CIBFEMethod::initializeFEData()
{
    if (d_fe_data_initialized) return;
    IBFEMethod::initializeFEData();
    d_fe_data_initialized = false;

    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        // Get mesh info.
        EquationSystems* equation_systems = d_equation_systems[part];
        const MeshBase& mesh = equation_systems->get_mesh();
        d_num_nodes[part] = mesh.n_nodes();

        // Assemble additional systems.
        System& U_constraint_system = equation_systems->get_system<System>(CONSTRAINT_VELOCITY_SYSTEM_NAME);
        U_constraint_system.assemble_before_solve = false;
        U_constraint_system.assemble();
    }

    d_fe_data_initialized = true;
    return;
} // initializeFEData

void CIBFEMethod::registerEulerianVariables()
{
    // Register a cc variable for plotting nodal Lambda.
    const IntVector<NDIM> ib_ghosts = getMinimumGhostCellWidth();

    d_eul_lambda_var = new CellVariable<NDIM, double>(d_object_name + "::eul_lambda", NDIM);
    registerVariable(d_eul_lambda_idx, d_eul_lambda_var, ib_ghosts, d_ib_solver->getCurrentContext());

    return;
} // registerEulerianVariables

void CIBFEMethod::registerEulerianCommunicationAlgorithms()
{
    // intentionally left blank.
    return;
} // registerEulerianCommunicationAlgorithms

void CIBFEMethod::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                           Pointer<GriddingAlgorithm<NDIM> > gridding_alg,
                                           int u_data_idx,
                                           const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
                                           const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                           int integrator_step,
                                           double init_data_time,
                                           bool initial_time)
{
    IBFEMethod::initializePatchHierarchy(hierarchy, gridding_alg, u_data_idx, u_synch_scheds, u_ghost_fill_scheds,
                                         integrator_step, init_data_time, initial_time);
    d_is_initialized = false;

    // Zero-out Eulerian lambda.
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    if (initial_time)
    {
        // Initialize the S[lambda] variable.
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CellData<NDIM, double> > lambda_data = patch->getPatchData(d_eul_lambda_idx);
                lambda_data->fillAll(0.0);
            }
        }
    }

    // Register Eulerian lambda with visit.
    if (d_output_eul_lambda && d_visit_writer)
    {
        d_visit_writer->registerPlotQuantity("S_lambda", "VECTOR", d_eul_lambda_idx, 0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (d == 0) d_visit_writer->registerPlotQuantity("S_lambda_x", "SCALAR", d_eul_lambda_idx, d);
            if (d == 1) d_visit_writer->registerPlotQuantity("S_lambda_y", "SCALAR", d_eul_lambda_idx, d);
            if (d == 2) d_visit_writer->registerPlotQuantity("S_lambda_z", "SCALAR", d_eul_lambda_idx, d);
        }
    }

    d_is_initialized = true;
    return;
} // initializePatchHierarchy

void CIBFEMethod::registerPreProcessSolveFluidEquationsCallBackFcn(preprocessSolveFluidEqn_callbackfcn callback,
                                                                   void* ctx)
{
    d_prefluidsolve_callback_fcns.push_back(callback);
    d_prefluidsolve_callback_fcns_ctx.push_back(ctx);

    return;
} // registerPreProcessSolveFluidEquationsCallBackFcn

void CIBFEMethod::preprocessSolveFluidEquations(double current_time, double new_time, int cycle_num)
{
    // Call any registered pre-fluid solve callback functions.
    for (unsigned i = 0; i < d_prefluidsolve_callback_fcns.size(); ++i)
    {
        d_prefluidsolve_callback_fcns[i](current_time, new_time, cycle_num, d_prefluidsolve_callback_fcns_ctx[i]);
    }

    return;
} // preprocessSolveFluidEquations

void CIBFEMethod::registerVisItDataWriter(Pointer<VisItDataWriter<NDIM> > visit_writer)
{
    d_visit_writer = visit_writer;
    return;
} // registerVisItDataWriter

int CIBFEMethod::getStructuresLevelNumber()
{
    return d_hierarchy->getFinestLevelNumber();

} // getStructuresLevelNumber

Pointer<PatchHierarchy<NDIM> > CIBFEMethod::getPatchHierarchy()
{
    return d_hierarchy;
} // getPatchHierarchy

Pointer<IBTK::HierarchyMathOps> CIBFEMethod::getHierarchyMathOps()
{
    return IBStrategy::getHierarchyMathOps();
} // getHierarchyMathOps

void CIBFEMethod::putToDatabase(Pointer<Database> db)
{
    db->putInteger("CIBFE_METHOD_VERSION", CIBFE_METHOD_VERSION);
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        std::ostringstream U, W;
        U << "U_" << part;
        W << "W_" << part;
        db->putDoubleArray(U.str(), &d_trans_vel_current[part][0], 3);
        db->putDoubleArray(W.str(), &d_rot_vel_current[part][0], 3);
    }

    return;
} // putToDatabase

//////////////////////////////////////////// PRIVATE ////////////////////////

void CIBFEMethod::commonConstructor(Pointer<Database> input_db)
{
    // Resize some arrays.
    d_num_nodes.resize(d_num_rigid_parts);
    d_constrained_velocity_fcns_data.resize(d_num_rigid_parts);
    d_ext_force_torque_fcn_data.resize(d_num_rigid_parts);

    // Set some default values.
    d_compute_L2_projection = false;
    d_output_eul_lambda = false;

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    // Add additional variable corresponding to constraint velocity.
    for (unsigned int part = 0; part < d_num_rigid_parts; ++part)
    {
        EquationSystems* equation_systems = d_equation_systems[part];
        System& U_constraint_system = equation_systems->add_system<System>(CONSTRAINT_VELOCITY_SYSTEM_NAME);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "U_constraint_" << d;
            U_constraint_system.add_variable(os.str(), d_fe_order, d_fe_family);
        }
    }

    // Keep track of the initialization state.
    d_fe_data_initialized = false;
    d_is_initialized = false;
    d_constraint_force_is_initialized = false;
    d_lag_velvec_is_initialized = false;

    return;
} // commonConstructor

void CIBFEMethod::getFromInput(Pointer<Database> input_db, bool /*is_from_restart*/)
{
    // Get some input values.
    d_compute_L2_projection = input_db->getBoolWithDefault("compute_L2_projection", d_compute_L2_projection);
    d_output_eul_lambda = input_db->getBoolWithDefault("output_eul_lambda", d_output_eul_lambda);

    return;
} // getFromInput

void CIBFEMethod::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR("CIBFEMethod::getFromRestart(): Restart database corresponding to "
                   << d_object_name << " not found in restart file." << std::endl);
    }

    int ver = db->getInteger("CIBFE_METHOD_VERSION");
    if (ver != CIBFE_METHOD_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }

    for (unsigned int part = 0; part < d_num_rigid_parts; ++part)
    {
        std::ostringstream U, W;
        U << "U_" << part;
        W << "W_" << part;
        db->getDoubleArray(U.str(), &d_trans_vel_current[part][0], 3);
        db->getDoubleArray(W.str(), &d_rot_vel_current[part][0], 3);
    }

    return;
} // getFromRestart

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
