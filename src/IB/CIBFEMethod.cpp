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
#include "ibtk/ibtk_utilities.h"
#include "ibtk/LSiloDataWriter.h"
#include "ibtk/PETScMultiVec.h"
#include "libmesh/equation_systems.h"
#include "libmesh/system.h"

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of CIBFEMethod restart file data.
static const int CIBFE_METHOD_VERSION = 1;
}

const std::string CIBFEMethod::CONSTRAINT_FORCE_SYSTEM_NAME = "IB force system";
const std::string CIBFEMethod::CONSTRAINT_VELOCITY_SYSTEM_NAME =
	"IB constrained velocity system";

/////////////////////////////// PUBLIC ///////////////////////////////////////

CIBFEMethod::CIBFEMethod(
    const std::string& object_name,
    Pointer<Database> input_db,
    Mesh* mesh,
    int max_level_number,
    bool register_for_restart)
    : IBFEMethod(object_name,input_db,mesh,max_level_number,register_for_restart),
	  CIBStrategy(1)
{
    commonConstructor(input_db);
    return;
}// CIBFEMethod

CIBFEMethod::CIBFEMethod(
    const std::string& object_name,
    Pointer<Database> input_db,
    const std::vector<Mesh*>& meshes,
    int max_level_number,
    bool register_for_restart)
    : IBFEMethod(object_name, input_db, meshes, max_level_number, register_for_restart),
	  CIBStrategy(static_cast<unsigned>(meshes.size()))
{
    commonConstructor(input_db);
    return;
}// CIBFEMethod

CIBFEMethod::~CIBFEMethod()
{
    // IBFEMethod class is responsible for unregistering
	// the object and deleting the equation systems.
	return;
} // ~CIBFEMethod

void 
CIBFEMethod::registerConstrainedVelocityFunction(
    ConstrainedNodalVelocityFcnPtr nodalvelfcn,
    ConstrainedCOMVelocityFcnPtr   comvelfcn,
    void* ctx,
    unsigned int part)
{
    TBOX_ASSERT(part < d_num_rigid_parts);
    registerConstrainedVelocityFunction(ConstrainedVelocityFcnsData(nodalvelfcn, comvelfcn, ctx), part);
    return;
}// registerConstrainedVelocityFunction

void 
CIBFEMethod::registerConstrainedVelocityFunction(
    const ConstrainedVelocityFcnsData& data,
    unsigned int part)
{
    TBOX_ASSERT(part < d_num_rigid_parts);
    d_constrained_velocity_fcns_data[part] = data;
    return;
}// registerConstrainedVelocityFunction

void
CIBFEMethod::preprocessIntegrateData(
    double current_time, 
    double new_time, 
    int num_cycles)
{
	// Create most of the FE data vecs in the base class.
	IBFEMethod::preprocessIntegrateData(current_time, new_time, num_cycles);

	// Create vecs for constraint force and velocities.
    d_L_systems.resize(d_num_rigid_parts);
	d_L_current_vecs.resize(d_num_rigid_parts);
	d_L_new_vecs.resize(d_num_rigid_parts);
	d_L_half_vecs.resize(d_num_rigid_parts);
    d_L_IB_ghost_vecs.resize(d_num_rigid_parts);

    d_U_constrained_systems.resize(d_num_rigid_parts);
    d_U_constrained_current_vecs.resize(d_num_rigid_parts);
    d_U_constrained_new_vecs.resize(d_num_rigid_parts);
    d_U_constrained_half_vecs.resize(d_num_rigid_parts);
	
	// PETSc wrappers.
	d_vL_current.resize(d_num_rigid_parts, PETSC_NULL);
	d_vL_new.resize(d_num_rigid_parts, PETSC_NULL);
	
    for (unsigned int part = 0; part < d_num_rigid_parts; ++part)
    {

        d_L_systems[part] = &d_equation_systems[part]->get_system(CONSTRAINT_FORCE_SYSTEM_NAME);
		d_L_current_vecs[part] = dynamic_cast<PetscVector<double>*>(
		    d_U_systems[part]->current_local_solution.get());
		d_L_new_vecs[part] = dynamic_cast<PetscVector<double>*>(
			d_L_current_vecs[part]->clone().release()); // WARNING: must be manually deleted
		d_L_half_vecs[part] = dynamic_cast<PetscVector<double>*>(
			d_L_current_vecs[part]->clone().release()); // WARNING: must be manually deleted
        d_L_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_fe_data_managers[part]->buildGhostedSolutionVector(CONSTRAINT_FORCE_SYSTEM_NAME,
                                                                 /*localize_data*/ false));
		d_vL_current[part] = d_L_current_vecs[part]->vec();
        d_vL_new[part]     = d_L_new_vecs[part]->vec();
		
        d_U_constrained_systems[part] = &d_equation_systems[part]->get_system(CONSTRAINT_VELOCITY_SYSTEM_NAME);
        d_U_constrained_current_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_U_constrained_systems[part]->current_local_solution.get());
        d_U_constrained_new_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_U_constrained_current_vecs[part]->clone().release()); // WARNING: must be manually deleted
        d_U_constrained_half_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_U_constrained_current_vecs[part]->clone().release()); // WARNING: must be manually deleted


        // Initialize L^{n+1/2} and L^{n+1} to equal L^{n}, and initialize
        // U_k^{n+1/2} and U_k^{n+1} to equal U_k^{n}.
		d_L_current_vecs[part]->localize(*d_L_half_vecs[part]);
		d_L_current_vecs[part]->localize(*d_L_new_vecs[part]);
        d_U_constrained_current_vecs[part]->localize(*d_U_constrained_half_vecs[part]);
        d_U_constrained_current_vecs[part]->localize(*d_U_constrained_new_vecs[part]);
       
		if (!d_solve_rigid_vel[part])
        { 
            d_constrained_velocity_fcns_data[part].comvelfcn(d_current_time,d_trans_vel_current[part],
                d_rot_vel_current[part]);
            d_constrained_velocity_fcns_data[part].comvelfcn(d_half_time,d_trans_vel_half[part],
                d_rot_vel_half[part]);
            d_constrained_velocity_fcns_data[part].comvelfcn(d_new_time,d_trans_vel_new[part],
                d_rot_vel_new[part]);
        }	
    }
	
	VecCreateMultiVec(PETSC_COMM_WORLD, d_num_rigid_parts, &d_vL_current[0], &d_mv_L_current);
	VecCreateMultiVec(PETSC_COMM_WORLD, d_num_rigid_parts, &d_vL_new[0], &d_mv_L_new);
	
    return;
}// preprocessIntegrateData

void
CIBFEMethod::postprocessIntegrateData(
    double current_time,
    double new_time,
    int num_cycles)
{
	// Clean the temporary FE data vecs.
	IBFEMethod::postprocessIntegrateData(current_time, new_time, num_cycles);
    for (unsigned part = 0; part < d_num_rigid_parts; ++part)
    {
        // Reset time-dependent Lagrangian data.
        *d_L_systems[part]->solution             = *d_L_new_vecs[part];
        *d_U_constrained_systems[part]->solution = *d_U_constrained_new_vecs[part];
	
        // Deallocate Lagrangian scratch data.
		delete d_L_new_vecs[part];
		delete d_L_half_vecs[part];
        delete d_U_constrained_new_vecs[part];
        delete d_U_constrained_half_vecs[part];
    }

    d_L_systems.clear();
	d_L_current_vecs.clear();
	d_L_new_vecs.clear();
    d_L_half_vecs.clear();
    d_L_IB_ghost_vecs.clear();
	d_vL_current.clear();
	d_vL_new.clear();
	VecDestroy(&d_mv_L_current);
	VecDestroy(&d_mv_L_new);

    d_U_constrained_systems.clear();
    d_U_constrained_current_vecs.clear();
    d_U_constrained_new_vecs.clear();
    d_U_constrained_half_vecs.clear();
    
    // New center of mass translational and rotational velocity becomes
    // current velocity for the next time step.
    d_trans_vel_current = d_trans_vel_new;
    d_rot_vel_current   = d_rot_vel_new;
	
    return;
}// postprocessIntegrateData

void 
CIBFEMethod::interpolateVelocity(
    const int u_data_idx,
    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
    const std::vector<Pointer<RefineSchedule<NDIM> > >&  u_ghost_fill_scheds,
    const double data_time)
{
	if (d_lag_vel_is_initialized)
	{
		IBFEMethod::interpolateVelocity(u_data_idx, u_synch_scheds, u_ghost_fill_scheds, data_time);
		d_lag_vel_is_initialized = false;
	}
    return;
}// interpolateVelocity

void
CIBFEMethod::eulerStep(
    const double current_time,
    const double new_time)
{
    const double dt = new_time-current_time;

    // Compute center of mass and moment of inertia of the body at t^n.
    computeCOMandMOIOfStructures(d_center_of_mass_current,d_moment_of_inertia_current,d_X_current_vecs);

    // Fill the rotation matrix of structures with rotation angle 0.5*(W^n)*dt.
    std::vector<Eigen::Matrix3d> rotation_mat(d_num_rigid_parts, Eigen::Matrix3d::Zero());
    for (unsigned int part = 0; part < d_num_rigid_parts; ++part)
    {
        for (int i = 0; i < 3; ++i)
			rotation_mat[part](i,i) = 1.0;
    }  
    setRotationMatrix(d_rot_vel_current,rotation_mat,0.5*dt);

    // Rotate the body with current rotational velocity about origin
    // and translate the body to predicted position X^n+1/2.   
    Eigen::Vector3d dr    = Eigen::Vector3d::Zero();
    Eigen::Vector3d Rxdr  = Eigen::Vector3d::Zero();
    for (int part = 0; part < d_num_rigid_parts; ++part)
    {
		EquationSystems* equation_systems = d_equation_systems[part];
        MeshBase& mesh = equation_systems->get_mesh();
        const unsigned int total_local_nodes = mesh.n_nodes_on_proc(SAMRAI_MPI::getRank());
        System& X_system = *d_X_systems[part];
        const unsigned int X_sys_num = X_system.number();
	    PetscVector<double>& X_current = *d_X_current_vecs[part];
        PetscVector<double>& X_half    = *d_X_half_vecs[part];
	
	    std::vector<std::vector<unsigned int> > nodal_X_indices(NDIM);
	    std::vector<std::vector<double> > nodal_X_values(NDIM);
	    for (unsigned int d = 0; d < NDIM; ++d)
	    {
            nodal_X_indices[d].reserve(total_local_nodes);
            nodal_X_values[d].reserve(total_local_nodes);
        }
	
        for (MeshBase::node_iterator it = mesh.local_nodes_begin();
             it != mesh.local_nodes_end(); ++it)
        {
            const Node* const n = *it;
            if (n->n_vars(X_sys_num))
            {
                TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    nodal_X_indices[d].push_back(n->dof_number(X_sys_num, d, 0));  
                }
	        }
	    }
        
        for (unsigned int d = 0; d < NDIM; ++d)
        {
			X_current.get(nodal_X_indices[d],nodal_X_values[d]);
        }
        
        for (unsigned int k = 0; k < total_local_nodes; ++k)
	    {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dr[d] = nodal_X_values[d][k] - d_center_of_mass_current[part][d];
            }
	    
            // Rotate dr vector using the rotation matrix.
            Rxdr = rotation_mat[part]*dr; 
            for (unsigned int d = 0; d < NDIM; ++d) 
            { 
                X_half.set(nodal_X_indices[d][k], d_center_of_mass_current[part][d] + Rxdr[d] + 
						   0.5*dt*d_trans_vel_current[part][d]);
            }	    
        }
    }
    
    // Compute the COM and MOI at mid-point. 
    computeCOMandMOIOfStructures(d_center_of_mass_half, d_moment_of_inertia_half, d_X_half_vecs);    

    return;
}// eulerStep

void
CIBFEMethod::midpointStep(
    const double current_time,
    const double new_time)
{
    const double dt = new_time - current_time;
    // Fill the rotation matrix of structures with rotation angle (W^n+1)*dt.
    std::vector<Eigen::Matrix3d> rotation_mat(d_num_rigid_parts,Eigen::Matrix3d::Zero());
    for (unsigned int part = 0; part < d_num_rigid_parts; ++part)
    {
        for (int i = 0; i < 3; ++i)
			rotation_mat[part](i,i) = 1.0;
    }  
    setRotationMatrix(d_rot_vel_half,rotation_mat,dt);
    
    // Rotate the body with current rotational velocity about origin
    // and translate the body to predicted position X^n+1/2.   
    Eigen::Vector3d dr    = Eigen::Vector3d::Zero();
    Eigen::Vector3d Rxdr  = Eigen::Vector3d::Zero();
    for (int part = 0; part < d_num_rigid_parts; ++part)
    {
        EquationSystems* equation_systems = d_equation_systems[part];
        MeshBase& mesh = equation_systems->get_mesh();
        const unsigned int total_local_nodes = mesh.n_nodes_on_proc(SAMRAI_MPI::getRank());
        System& X_system = *d_X_systems[part];
        const unsigned int X_sys_num = X_system.number();
        PetscVector<double>& X_current = *d_X_current_vecs[part];
        PetscVector<double>& X_new     = *d_X_new_vecs[part];
	
	    std::vector<std::vector<unsigned int> > nodal_X_indices(NDIM);
        std::vector<std::vector<double> > nodal_X_values(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            nodal_X_indices[d].reserve(total_local_nodes);
            nodal_X_values[d].reserve(total_local_nodes);
        }
	
        for (MeshBase::node_iterator it = mesh.local_nodes_begin();
             it != mesh.local_nodes_end(); ++it)
        {
            const Node* const n = *it;
            if (n->n_vars(X_sys_num))
            {
                TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    nodal_X_indices[d].push_back(n->dof_number(X_sys_num, d, 0));  
                }
            }
        }
        
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_current.get(nodal_X_indices[d],nodal_X_values[d]);
        }
        
        for (unsigned int k = 0; k < total_local_nodes; ++k)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dr[d] = nodal_X_values[d][k] - d_center_of_mass_current[part][d];
            }
	    
            // Rotate dr vector using the rotation matrix.
            Rxdr = rotation_mat[part]*dr; 
            for (unsigned int d = 0; d < NDIM; ++d) 
            { 
                X_new.set(nodal_X_indices[d][k], d_center_of_mass_current[part][d] + Rxdr[d] + 
                          dt*d_trans_vel_half[part][d]);
            }	    
        }
    }

    return;
}// midpointStep

void
CIBFEMethod::trapezoidalStep(
    const double /*current_time*/,
    const double /*new_time*/)
{
    TBOX_ERROR("CIBFEMethod does not support trapezoidal time-stepping rule for position update."
               << " Only midpoint rule is supported for position update." << std::endl);
    
    return;    
}// trapezoidalStep

void 
CIBFEMethod::computeLagrangianForce(
    double data_time)
{
    //intentionally left blank
    return;
}// computeLagrangianForce

void
CIBFEMethod::spreadForce(
    int f_data_idx,
    RobinPhysBdryPatchStrategy* f_phys_bdry_op,
    const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >&
        f_prolongation_scheds,
    double data_time)
{
    if (d_constraint_force_is_initialized)
	{
		IBFEMethod::spreadForce(f_data_idx, f_phys_bdry_op, f_prolongation_scheds, data_time);
		d_constraint_force_is_initialized = false;
	}
    return;
    
}// spreadForce
	
void
CIBFEMethod::setConstraintForce(
	Vec L,
	double /*data_time*/,
	double scale)
{
	//Unpack the Lambda vector.
	Vec* vL;
	VecMultiVecGetSubVecs(L, &vL);
	
#if !defined(NDEBUG)
	PetscInt size;
	VecMultiVecGetNumberOfSubVecs(L, &size);
	TBOX_ASSERT(size == d_num_rigid_parts);
#endif
	
	for (unsigned part = 0; part < d_num_rigid_parts; ++part)
	{
		PetscVector<double>* F_vec = d_F_half_vecs[part];
		VecCopy(vL[part], F_vec->vec());
		VecScale(F_vec->vec(), scale);
	}
	d_constraint_force_is_initialized = true;
	
	return;
	
}// setConstraintForce

/*void
CIBFEMethod::spreadForce(
	int f_data_idx,
	Vec L,
    RobinPhysBdryPatchStrategy* f_phys_bdry_op,
	const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >&
		f_prolongation_scheds,
	double data_time,
	double scale)
{
	//Unpack the Lambda vector.
	Vec* vL;
	VecMultiVecGetSubVecs(L, &vL);
	for (unsigned int part = 0; part < d_num_parts; ++part)
	{
		spreadForce(f_data_idx, part, vL[part], f_phys_bdry_op, f_prolongation_scheds, data_time, scale);
	}
	return;
		
}// spreadForce
	
void
cRigidIBFEMethod::spreadForce(
	int f_data_idx,
	const unsigned int part,
    Vec L,
	RobinPhysBdryPatchStrategy* f_phys_bdry_op,
	const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >&
		f_prolongation_scheds,
	double data_time,
	double scale)
{

	PetscVector<double>* X_vec = d_X_half_vecs[part];
	PetscVector<double>* X_ghost_vec = d_X_IB_ghost_vecs[part];
	PetscVector<double>* L_vec = d_L_half_vecs[part];
	PetscVector<double>* L_ghost_vec = d_L_IB_ghost_vecs[part];
	VecCopy(L, L_vec->vec());
	VecScale(L_vec->vec(), scale);
	X_vec->localize(*X_ghost_vec);
	L_vec->localize(*L_ghost_vec);
	d_fe_data_managers[part]->spread(f_data_idx,
									*L_ghost_vec,
									*X_ghost_vec,
									CONSTRAINT_FORCE_SYSTEM_NAME,
									f_phys_bdry_op,
									data_time);
	
	return;
		
}// spreadForce
	
void
cRigidIBFEMethod::interpolateVelocity(
	int u_data_idx,
	Vec V,
	const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >&
		u_synch_scheds,
	const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >&
		u_ghost_fill_scheds,
	double data_time,
	double scale)
{
	//Unpack the velocity vector.
	Vec* vV;
	VecMultiVecGetSubVecs(V, &vV);
	for (unsigned int part = 0; part < d_num_parts; ++part)
	{
		interpolateVelocity(u_data_idx, part, vV[part], u_synch_scheds, u_ghost_fill_scheds, data_time, scale);
	}
	return;
	
}// interpolateVelocity

void
cRigidIBFEMethod::interpolateVelocity(
	int u_data_idx,
	const unsigned int part,
	Vec V,
	const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >&
		u_synch_scheds,
	const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >&
		u_ghost_fill_scheds,
	double data_time,
	double scale)
{
	PetscVector<double>* X_vec = d_X_half_vecs[part];
	PetscVector<double>* X_ghost_vec = d_X_IB_ghost_vecs[part];
	PetscVector<double>* U_vec = d_U_half_vecs[part];

	X_vec->localize(*X_ghost_vec);
    d_fe_data_managers[part]->interp(u_data_idx,
									 *U_vec,
									 *X_ghost_vec,
									 VELOCITY_SYSTEM_NAME,
									 u_ghost_fill_scheds,
									 data_time);
	VecCopy(U_vec->vec(), V);
	VecScale(V, scale);
	return;
	
}// interpolateVelocity
	
*/
void
CIBFEMethod::computeNetRigidGeneralizedForce(
	const unsigned int part,
    Vec L,
    RigidDOFVector& F)
{
	EquationSystems* equation_systems = d_equation_systems[part];
	MeshBase& mesh = equation_systems->get_mesh();
	const unsigned int dim = mesh.mesh_dimension();
	AutoPtr<QBase> qrule = QBase::build(d_quad_type, dim, d_quad_order);
	
	// Extract the FE system and DOF map, and setup the FE object.
	System& L_system = *d_L_systems[part];
	System& X_system = *d_X_systems[part];
	DofMap& L_dof_map = L_system.get_dof_map();
	DofMap& X_dof_map = X_system.get_dof_map();
	std::vector<std::vector<unsigned int> >L_dof_indices(NDIM);
	std::vector<std::vector<unsigned int> >X_dof_indices(NDIM);
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
	
	const Eigen::Vector3d& X_com          = d_center_of_mass_half[part];
	libMesh::PetscVector<double>& X_petsc = *d_X_half_vecs[part];
	if (!X_petsc.closed()) X_petsc.close();
	Vec X_global_vec = X_petsc.vec();
	Vec X_local_ghost_vec;
	VecGhostGetLocalForm(X_global_vec, &X_local_ghost_vec);
	double* X_local_ghost_soln;
	VecGetArray(X_local_ghost_vec, &X_local_ghost_soln);
	
	libMesh::PetscVector<double>& L_petsc = *d_L_half_vecs[part];
	Vec L_global_vec = L_petsc.vec();
	VecCopy(L, L_global_vec);
	L_petsc.close();                  // Sync ghost values.
	Vec L_local_ghost_vec;
	VecGhostGetLocalForm(L_global_vec, &L_local_ghost_vec);
	double* L_local_ghost_soln;
	VecGetArray(L_local_ghost_vec, &L_local_ghost_soln);
	
	F.setZero();
	boost::multi_array<double,2> X_node, L_node;
	double X_qp[NDIM], L_qp[NDIM];
	const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator el_end   = mesh.active_local_elements_end();
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
				F[d] += L_qp[d]*JxW_L[qp];
			}
#if (NDIM == 2)
			F[NDIM]   += (L_qp[1]*(X_qp[0]-X_com[0]) - L_qp[0]*(X_qp[1]-X_com[1]))*JxW_L[qp];
#elif (NDIM == 3)
			F[NDIM]   += (L_qp[2]*(X_qp[1]-X_com[1]) - L_qp[1]*(X_qp[2]-X_com[2]))*JxW_L[qp];
			F[NDIM+1] += (L_qp[0]*(X_qp[2]-X_com[2]) - L_qp[2]*(X_qp[0]-X_com[0]))*JxW_L[qp];
			F[NDIM+2] += (L_qp[1]*(X_qp[0]-X_com[0]) - L_qp[0]*(X_qp[1]-X_com[1]))*JxW_L[qp];
#endif
		}
	}
	SAMRAI_MPI::sumReduction(&F[0], NDIM*(NDIM+1)/2);
	return;
	
}// computeNetRigidGeneralizedForce
	
void
cRigidIBFEMethod::setRigidBodyVelocity(
	const unsigned int part,
	const RigidDOFVector& U,
	Vec V)
{
	PetscVector<double>& U_k    = *d_U_constrained_half_vecs[part];
	PetscVector<double>& X_half = *d_X_half_vecs[part];
	void* ctx = d_constrained_velocity_fcn_data[part].ctx;
	d_constrained_velocity_fcn_data[part].nodalvelfcn(U_k,
												      U,
													  X_half,
													  d_center_of_mass_half[part],
													  d_equation_systems[part],
													  d_new_time,
													  ctx);
	Vec* vV;
	VecMultiVecGetSubVecs(V, &vV);
	VecCopy(U_k.vec(), vV[part]);
	return;
	
}// setRigidBodyVelocity
	
void
cRigidIBFEMethod::getRigidBodyForce(
	Vec* L,
	const double time)
{
	
	if (MathUtilities<double>::equalEps(time, d_current_time))
	{
		*L = d_mv_L_current;
	}
	else if (MathUtilities<double>::equalEps(time, d_new_time))
	{
		*L = d_mv_L_new;
	}
	else
	{
		plog << "Warning cRigidIBFEMethod::getRigidBodyForce() : constraint force "
		     << "enquired at some other time than current or new time.";
		return;
	}
	
	return;
}// getRigidBodyForce
	
void
cRigidIBFEMethod::initializeFEData()
{
    if (d_fe_data_initialized) return;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        // Initialize FE equation systems.
        EquationSystems* equation_systems = d_equation_systems[part];
        equation_systems->init();
        initializeCoordinates(part);
        updateCoordinateMapping(part);
		
        // Assemble systems.
        System& X_system = equation_systems->get_system<System>(COORDS_SYSTEM_NAME);
        System& dX_system = equation_systems->get_system<System>(COORD_MAPPING_SYSTEM_NAME);
        System& U_system = equation_systems->get_system<System>(VELOCITY_SYSTEM_NAME);
        System& L_system = equation_systems->get_system<System>(CONSTRAINT_FORCE_SYSTEM_NAME);
        System& U_constraint_system = equation_systems->get_system<System>(CONSTRAINT_VELOCITY_SYSTEM_NAME);
	
        X_system.assemble_before_solve = false;
        X_system.assemble();

        dX_system.assemble_before_solve = false;
        dX_system.assemble();

        U_system.assemble_before_solve = false;
        U_system.assemble();

        L_system.assemble_before_solve = false;
        L_system.assemble();

        U_constraint_system.assemble_before_solve = false;
        U_constraint_system.assemble();
        
        // Set up boundary conditions and nodal information.  Specifically, 
	    // add appropriate boundary IDs to the BoundaryInfo object associated
	    // with the mesh, and add DOF constraints for the nodal forces and velocities.
        const MeshBase& mesh = equation_systems->get_mesh();
	    d_num_nodes[part]    = mesh.n_nodes();
        DofMap& U_dof_map = U_system.get_dof_map();
        const unsigned int U_sys_num = U_system.number();
        MeshBase::const_element_iterator el_it = mesh.elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.elements_end();
        for (; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                const bool at_mesh_bdry = !elem->neighbor(side);
                if (!at_mesh_bdry) continue;

                static const short int dirichlet_bdry_id_set[3] = {
                    FEDataManager::ZERO_DISPLACEMENT_X_BDRY_ID,
                    FEDataManager::ZERO_DISPLACEMENT_Y_BDRY_ID,
                    FEDataManager::ZERO_DISPLACEMENT_Z_BDRY_ID
                };
                const short int dirichlet_bdry_ids =
                    get_dirichlet_bdry_ids(mesh.boundary_info->boundary_ids(elem, side));
                if (!dirichlet_bdry_ids) continue;

                for (unsigned int n = 0; n < elem->n_nodes(); ++n)
                {
                    if (!elem->is_node_on_side(n, side)) continue;

                    Node* node = elem->get_node(n);
                    mesh.boundary_info->add_node(node, dirichlet_bdry_ids);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (!(dirichlet_bdry_ids & dirichlet_bdry_id_set[d])) continue;
                        if (node->n_dofs(U_sys_num))
                        {
                            const int U_dof_index = node->dof_number(U_sys_num, d, 0);
                            DofConstraintRow U_constraint_row;
                            U_constraint_row[U_dof_index] = 1.0;
                            U_dof_map.add_constraint_row(
                                U_dof_index, U_constraint_row, 0.0, false);
                        }
                    }
                }
            }
        }
    } 
    d_fe_data_initialized = true;
    return;
}// initializeFEData

void
cRigidIBFEMethod::registerEulerianVariables()
{
    // Register a cc variable for plotting nodal Lambda.
    const IntVector<NDIM> ib_ghosts = getMinimumGhostCellWidth(); 
    d_eul_lambda_var = new CellVariable<NDIM,double>(d_object_name + "::eul_lambda",NDIM);
    registerVariable(d_eul_lambda_idx, d_eul_lambda_var, ib_ghosts, d_ib_solver->getCurrentContext());
    
    return;  
}// registerEulerianVariables

void
cRigidIBFEMethod::registerEulerianCommunicationAlgorithms()
{
    //intentionally left blank.   
    return;
}// registerEulerianCommunicationAlgorithms

void 
cRigidIBFEMethod::initializePatchHierarchy(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg,
    int /*u_data_idx*/,
    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*u_synch_scheds*/,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& /*u_ghost_fill_scheds*/,
    int /*integrator_step*/,
    double /*init_data_time*/,
    bool initial_time)
{
    // Cache pointers to the patch hierarchy and gridding algorithm.
    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;

    // Initialize the FE data manager.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->reinitElementMappings();
    }
    
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
				Pointer<CellData<NDIM,double> > lambda_data = patch->getPatchData(d_eul_lambda_idx);
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
}// initializePatchHierarchy

void
CIBFEMethod::registerPreProcessSolveFluidEquationsCallBackFcn(
    preprocessSolveFluidEqn_callbackfcn callback, 
    void* ctx)
{
    d_prefluidsolve_callback_fcns.push_back(callback);
    d_prefluidsolve_callback_fcns_ctx.push_back(ctx);

    return;
}// registerPreProcessSolveFluidEquationsCallBackFcn

void
CIBFEMethod::preprocessSolveFluidEquations(
    double current_time, 
    double new_time, 
    int cycle_num)
{
    // Call any registered pre-fluid solve callback functions.
    for (unsigned i = 0; i < d_prefluidsolve_callback_fcns.size(); ++i) 
    {  
		d_prefluidsolve_callback_fcns[i](current_time, new_time, cycle_num,
										 d_prefluidsolve_callback_fcns_ctx[i]);
    }

    return;
}// preprocessSolveFluidEquations

void
cRigidIBFEMethod::registerVisItDataWriter(
    Pointer<VisItDataWriter<NDIM> > visit_writer)
{
    d_visit_writer = visit_writer;
    return;
}// registerVisItDataWriter

int
cRigidIBFEMethod::getStructuresLevelNumber()
{
    return d_hierarchy->getFinestLevelNumber();

}// getStructuresLevelNumber

Pointer<PatchHierarchy<NDIM> >
cRigidIBFEMethod::getPatchHierarchy()
{
    return d_hierarchy;
}// getPatchHierarchy

Pointer<::HierarchyMathOps>
cRigidIBFEMethod::getHierarchyMathOps()
{
    return getHierarchyMathOps();
}// getHierarchyMathOps 

void
cRigidIBFEMethod::putToDatabase(
    Pointer<Database> db)
{
    db->putInteger("cRIGID_IBFE_METHOD_VERSION", cRIGID_IBFE_METHOD_VERSION);
    db->putIntegerArray("d_ghosts", d_ghosts, NDIM);
    db->putString("d_fe_family", Utility::enum_to_string<FEFamily>(d_fe_family));
    db->putString("d_fe_order", Utility::enum_to_string<Order>(d_fe_order));
    db->putString("d_quad_type", Utility::enum_to_string<QuadratureType>(d_quad_type));
    db->putString("d_quad_order", Utility::enum_to_string<Order>(d_quad_order));
    db->putBool("d_use_consistent_mass_matrix", d_use_consistent_mass_matrix);
    
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
		std::ostringstream U, W;
		U << "U_" << part;
		W << "W_" << part;
        db->putDoubleArray(U.str(), &d_trans_vel_current[part][0],3);
		db->putDoubleArray(W.str(), &d_rot_vel_current  [part][0],3);
    }
    
    return;
}//putToDatabase

//////////////////////////////////////////// PRIVATE ////////////////////////

void CIBFEMethod::commonConstructor(
    Pointer<Database> input_db)
{
    // Resize some arrays.
    d_num_nodes.resize(d_num_parts);
	
    // Set some default values.
    d_output_eul_lambda    = false;

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);


    // Add additional variables corresponding to constraint velocity and
	// constraint force.
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        EquationSystems* equation_systems = d_equation_systems[part];
        System& U_constraint_system = equation_systems->add_system<System>(CONSTRAINT_VELOCITY_SYSTEM_NAME);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "U_constraint_" << d;
            U_constraint_system.add_variable(os.str(), d_fe_order, d_fe_family);
        }
        
        System& Lambda_system = equation_systems->add_system<System>(CONSTRAINT_FORCE_SYSTEM_NAME);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "Lambda_" << d;
            Lambda_system.add_variable(os.str(), d_fe_order, d_fe_family);
        }
    }
    
    // Keep track of the initialization state.
    d_fe_data_initialized = false;
    d_is_initialized = false;
	d_constraint_force_is_initialized = false;
	d_lag_vel_is_initialized = false;
	
    return;
}// commonConstructor

void
CIBFEMethod::getFromInput(
    Pointer<Database> input_db,
    bool /*is_from_restart*/)
{  
    // Outputing S[lambda] and Lagrange Multiplier
    d_output_eul_lambda  = input_db->getBoolWithDefault("output_eul_lambda", d_output_eul_lambda);

    return;
}// getFromInput

void
CIBFEMethod::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR("cRigidIBFEMethod::getFromRestart(): Restart database corresponding to "
                   << d_object_name << " not found in restart file." << std::endl);
    }
    
    int ver = db->getInteger("CIBFE_METHOD_VERSION");
    if (ver != CIBFE_METHOD_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version."
                                 << std::endl);
    }    
	
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {      
		std::ostringstream U, W;
		U << "U_" << part;
		W << "W_" << part;
        db->getDoubleArray(U.str(),&d_trans_vel_current[part][0],3);
		db->getDoubleArray(W.str(),&d_rot_vel_current  [part][0],3);
    }
  
    return;
}// getFromRestart

void
CIBFEMethod::computeCOMandMOIOfStructures(
    std::vector<Eigen::Vector3d>& center_of_mass,
    std::vector<Eigen::Matrix3d>& moment_of_inertia,
    std::vector<PetscVector<double>*> X)
{
    //Find center of mass of the structures.
    for (int part = 0; part < d_num_parts; ++part)
    {
        // Extract FE mesh
        EquationSystems* equation_systems = d_equation_systems[part];
        MeshBase& mesh = equation_systems->get_mesh();
	    const unsigned int dim = mesh.mesh_dimension();
        AutoPtr<QBase> qrule = QBase::build(d_quad_type, dim, d_quad_order);
	
        // Extract the FE system and DOF map, and setup the FE object.
        System& x_system = *d_X_systems[part];
        DofMap& X_dof_map = x_system.get_dof_map();
		std::vector<std::vector<unsigned int> >X_dof_indices(NDIM);
        FEType fe_type = X_dof_map.variable_type(0);
        AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));
        fe->attach_quadrature_rule(qrule.get());
        const std::vector<double>& JxW = fe->get_JxW();
        const std::vector<std::vector<double> >& phi = fe->get_phi();
	    const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
	
		// Extract the nodal coordinates.
		PetscVector<double>& X_petsc = *X[part];
		if (!X_petsc.closed()) X_petsc.close();
		Vec X_global_vec  = X_petsc.vec();
        Vec X_local_ghost_vec;
		VecGhostGetLocalForm(X_global_vec, &X_local_ghost_vec);
        double* X_local_ghost_soln;
        VecGetArray(X_local_ghost_vec, &X_local_ghost_soln);	
	
		// Loop over the local elements to compute the local integrals.
		libMesh::TensorValue<double> dX_ds;
        boost::multi_array<double,2> X_node;
		double X_qp[NDIM];
		double vol_part = 0.0;
		center_of_mass[part].setZero();
        const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            const Elem* const elem = *el_it;
            fe->reinit(elem);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_dof_map.dof_indices(elem, X_dof_indices[d], d);
			}
			get_values_for_interpolation(X_node, X_petsc, X_local_ghost_soln, X_dof_indices);
	     
            const unsigned int n_qp = qrule->n_points();
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {   
				interpolate(X_qp, qp, X_node, phi);
				jacobian(dX_ds, qp, X_node, dphi);
		        const double det_dX_ds = dX_ds.det();
                for (unsigned int d = 0; d < NDIM; ++d)
				{
					center_of_mass[part][d] += X_qp[d]*det_dX_ds*JxW[qp];
				}
		        vol_part += det_dX_ds*JxW[qp];
            }
        }
        SAMRAI_MPI::sumReduction(&center_of_mass[part][0],NDIM);
	    SAMRAI_MPI::sumReduction(vol_part);
	
        for (unsigned int d = 0; d < NDIM; ++d)
        {   
            center_of_mass[part][d] /= vol_part; 
		}
		VecRestoreArray(X_local_ghost_vec, &X_local_ghost_soln);
		VecGhostRestoreLocalForm(X_global_vec, &X_local_ghost_vec);
    }	
    
    // Find moment of inertia tensor of structures.
    for (int part = 0; part < d_num_parts; ++part)
    {
        // Extract FE mesh
        EquationSystems* equation_systems = d_equation_systems[part];
        MeshBase& mesh = equation_systems->get_mesh();
		const unsigned int dim = mesh.mesh_dimension();
        AutoPtr<QBase> qrule = QBase::build(d_quad_type, dim, d_quad_order);
	
        // Extract the FE system and DOF map, and setup the FE object.
        System& x_system = *d_X_systems[part];
        DofMap& X_dof_map = x_system.get_dof_map();
		std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
        FEType fe_type = X_dof_map.variable_type(0);
        AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));
        fe->attach_quadrature_rule(qrule.get());
        const std::vector<double>& JxW = fe->get_JxW();
        const std::vector<std::vector<double> >& phi = fe->get_phi();
		const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
	
		// Extract the nodal coordinates.
		PetscVector<double>& X_petsc = *X[part];
		if (!X_petsc.closed()) X_petsc.close();
		Vec X_global_vec      = X[part]->vec();
		Vec X_local_ghost_vec;
		VecGhostGetLocalForm(X_global_vec, &X_local_ghost_vec);
        double* X_local_ghost_soln;
        VecGetArray(X_local_ghost_vec, &X_local_ghost_soln);	
	
		// Loop over the local elements to compute the local integrals.
		libMesh::TensorValue<double> dX_ds;
        boost::multi_array<double, 2> X_node;
		double X_qp[NDIM];
		moment_of_inertia[part].setZero();
		const Eigen::Vector3d& X_com = center_of_mass[part];
        const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();        
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            const Elem* const elem = *el_it;
            fe->reinit(elem);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_dof_map.dof_indices(elem, X_dof_indices[d], d);
			}
			get_values_for_interpolation(X_node, X_petsc, X_local_ghost_soln, X_dof_indices);
           
			const unsigned int n_qp = qrule->n_points();
            for (unsigned int qp = 0; qp < n_qp; ++qp)
            {   
				interpolate(X_qp, qp, X_node, phi);
				jacobian(dX_ds, qp, X_node, dphi);
				const double det_dX_ds = dX_ds.det();
#if (NDIM == 2)
				moment_of_inertia[part](0,0) += std::pow(X_qp[1] - X_com[1],2)*det_dX_ds*JxW[qp];
				moment_of_inertia[part](0,1) += -(X_qp[0] - X_com[0])*(X_qp[1] - X_com[1])*det_dX_ds*JxW[qp];
				moment_of_inertia[part](1,1) += std::pow(X_qp[0] - X_com[0],2)*det_dX_ds*JxW[qp];
				moment_of_inertia[part](2,2) += (std::pow(X_qp[0] - X_com[0],2) + std::pow(X_qp[1] - X_com[1],2))*det_dX_ds*JxW[qp];
#endif
#if (NDIM == 3)
				moment_of_inertia[part](0,0) += (std::pow(X_qp[1] - X_com[1],2) + std::pow(X_qp[2] - X_com[2],2))*det_dX_ds*JxW[qp];
				moment_of_inertia[part](0,1) += -(X_qp[0] - X_com[0])*(X_qp[1] - X_com[1])*det_dX_ds*JxW[qp];
	            moment_of_inertia[part](0,2) += -(X_qp[0] -X_com[0])*(X_qp[2] -X_com[2])*det_dX_ds*JxW[qp];
	            moment_of_inertia[part](1,1) += (std::pow(X_qp[0] - X_com[0],2) + std::pow(X_qp[2] - X_com[2],2))*det_dX_ds*JxW[qp];
				moment_of_inertia[part](1,2) += (-(X_qp[1] -X_com[1])*(X_qp[2]-X_com[2]))*det_dX_ds*JxW[qp];
	            moment_of_inertia[part](2,2) += (std::pow(X_qp[0] - X_com[0],2) + std::pow(X_qp[1] - X_com[1],2))*det_dX_ds*JxW[qp];
#endif		
            }
        }
        SAMRAI_MPI::sumReduction(&moment_of_inertia[part](0,0),9);
	
		//Fill-in symmetric part of inertia tensor.
		moment_of_inertia[part](1,0) = moment_of_inertia[part](0,1);
		moment_of_inertia[part](2,0) = moment_of_inertia[part](0,2);
        moment_of_inertia[part](2,1) = moment_of_inertia[part](1,2); 
	
		VecRestoreArray(X_local_ghost_vec, &X_local_ghost_soln);
		VecGhostRestoreLocalForm(X_global_vec, &X_local_ghost_vec);
    }	
    return;
}// computeCOMandMOIOfStructures

void
CIBFEMethod::setRotationMatrix(
    const std::vector<Eigen::Vector3d>& rot_vel,
    std::vector<Eigen::Matrix3d>& rot_mat,
    const double dt)
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        Eigen::Vector3d  e  = rot_vel[part];
        Eigen::Matrix3d& R  = rot_mat[part];
        const double norm_e = sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);

		if (norm_e > std::numeric_limits<double>::epsilon())
		{
			const double theta = norm_e*dt;
			for (int i = 0; i < 3; ++i) e[i] /= norm_e;
			const double c_t = cos(theta);
			const double s_t = sin(theta);
		     
			R(0,0) = c_t + (1.0-c_t)*e[0]*e[0];      R(0,1) = (1.0-c_t)*e[0]*e[1] - s_t*e[2];
			R(0,2) = (1.0-c_t)*e[0]*e[2] + s_t*e[1]; R(1,0) = (1.0-c_t)*e[1]*e[0] + s_t*e[2];
			R(1,1) = c_t + (1.0-c_t)*e[1]*e[1];      R(1,2) = (1.0-c_t)*e[1]*e[2] - s_t*e[0];
			R(2,0) = (1.0-c_t)*e[2]*e[0] - s_t*e[1]; R(2,1) = (1.0-c_t)*e[2]*e[1] + s_t*e[0];
			R(2,2) = c_t + (1.0-c_t)*e[2]*e[2];
		}
    } 

    return;
}// setRotationMatrix

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
