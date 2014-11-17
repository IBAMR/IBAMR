// Filename: cRigidIBFEMethod.cpp
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

#include "ibamr/cRigidIBFEMethod.h"
#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/ibtk_utilities.h"
#include "ibtk/libmesh_utilities.h"
#include "ibtk/LData.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/LSiloDataWriter.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/PETScMultiVec.h"
#include "tbox/TimerManager.h"
#include "HierarchyNodeDataOpsReal.h"
#include "HierarchyEdgeDataOpsReal.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/boundary_info.h"
#include "libmesh/compare_types.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/edge.h"
#include "libmesh/elem.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe_type.h"
#include "libmesh/fem_context.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/point.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/system.h"
#include "libmesh/tensor_value.h"
#include "libmesh/type_tensor.h"
#include "libmesh/type_vector.h"
#include "libmesh/vector_value.h"

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of cRigidIBFEMethod restart file data.
static const int cRIGID_IBFE_METHOD_VERSION = 1;
}

const std::string cRigidIBFEMethod::COORDS_SYSTEM_NAME = "IB coordinates system";
const std::string cRigidIBFEMethod::COORD_MAPPING_SYSTEM_NAME = "IB coordinate mapping system";
const std::string cRigidIBFEMethod::CONSTRAINT_FORCE_SYSTEM_NAME = "IB force system";
const std::string cRigidIBFEMethod::VELOCITY_SYSTEM_NAME = "IB velocity system";
const std::string cRigidIBFEMethod::CONSTRAINT_VELOCITY_SYSTEM_NAME = "IB constrained velocity system";
const std::string cRigidIBFEMethod::REGULATOR_SYSTEM_NAME = "IB regulator system";

inline short int 
get_dirichlet_bdry_ids(
    const std::vector<short int>& bdry_ids)
{
    short int dirichlet_bdry_ids = 0;
    for (std::vector<short int>::const_iterator cit = bdry_ids.begin(); cit != bdry_ids.end();
         ++cit)
    {
        const short int bdry_id = *cit;
        if (bdry_id == FEDataManager::ZERO_DISPLACEMENT_X_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_Y_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_Z_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_XY_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_XZ_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_YZ_BDRY_ID ||
            bdry_id == FEDataManager::ZERO_DISPLACEMENT_XYZ_BDRY_ID)
        {
            dirichlet_bdry_ids |= bdry_id;
        }
    }
    return dirichlet_bdry_ids;
}


/////////////////////////////// PUBLIC ///////////////////////////////////////

cRigidIBFEMethod::cRigidIBFEMethod(
    const std::string& object_name,
    Pointer<Database> input_db,
    Mesh* mesh,
    int max_level_number,
    bool register_for_restart)
    : cRigidIBStrategy(1)
{
    commonConstructor(object_name,
                      input_db,
                      std::vector<Mesh*>(1, mesh),
                      max_level_number,
                      register_for_restart);
    return;
}// cRigidIBFEMethod

cRigidIBFEMethod::cRigidIBFEMethod(
    const std::string& object_name,
    Pointer<Database> input_db,
    const std::vector<Mesh*>& meshes,
    int max_level_number,
    bool register_for_restart)
    : cRigidIBStrategy(static_cast<unsigned>(meshes.size()))
{
    commonConstructor(object_name,
					  input_db,
					  meshes,
					  max_level_number,
					  register_for_restart);
    return;
}// cRigidIBFEMethod

cRigidIBFEMethod::~cRigidIBFEMethod()
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        delete d_equation_systems[part];
    }
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }
    return;
} // ~cRigidIBFEMethod

FEDataManager* cRigidIBFEMethod::getFEDataManager(
    const unsigned int part) const
{
    TBOX_ASSERT(part < d_num_parts);
    return d_fe_data_managers[part];
} // getFEDataManager

const IntVector<NDIM>& 
cRigidIBFEMethod::getMinimumGhostCellWidth() const
{
    return d_ghosts;
}// getMinimumGhostCellWidth

void 
cRigidIBFEMethod::registerConstrainedVelocityFunction(
    ConstrainedVelocityFcnPtr    nodalvelfcn,
    ConstrainedCOMVelocityFcnPtr comvelfcn,
    void* ctx,
    unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    registerConstrainedVelocityFunction(ConstrainedVelocityFcnData(nodalvelfcn, comvelfcn, ctx), part);
    return;
}// registerConstrainedVelocityFunction

void 
cRigidIBFEMethod::registerConstrainedVelocityFunction(
    const ConstrainedVelocityFcnData& data,
    unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_constrained_velocity_fcn_data[part] = data;
    return;
}// registerConstrainedVelocityFunction

void 
cRigidIBFEMethod::registerInitialCoordinateMappingFunction(
    CoordinateMappingFcnPtr fcn,
    void* ctx,
    const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    registerInitialCoordinateMappingFunction(CoordinateMappingFcnData(fcn, ctx), part);
    return;
} // registerInitialCoordinateMappingFunction

void 
cRigidIBFEMethod::registerInitialCoordinateMappingFunction(
    const CoordinateMappingFcnData& data,
    const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_coordinate_mapping_fcn_data[part] = data;
    return;
}// registerInitialCoordinateMappingFunction

void 
cRigidIBFEMethod::registerLagBodyForceFunction(
    LagBodyForceFcnPtr fcn,
    const std::vector<unsigned int>& systems,
    void* ctx,
    const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    registerLagBodyForceFunction(LagBodyForceFcnData(fcn, systems, ctx), part);
    return;
}// registerLagBodyForceFunction

void 
cRigidIBFEMethod::registerLagBodyForceFunction(
    const LagBodyForceFcnData& data,
    const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_lag_body_force_fcn_data[part] = data;
    d_fcn_systems[part].insert(data.systems.begin(), data.systems.end());
    d_body_fcn_systems[part].insert(data.systems.begin(), data.systems.end());
    return;
}// registerLagBodyForceFunction

void 
cRigidIBFEMethod::registerLagSurfaceForceFunction(
    LagSurfaceForceFcnPtr fcn,
    const std::vector<unsigned int>& systems,
    void* ctx,
    const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    registerLagSurfaceForceFunction(LagSurfaceForceFcnData(fcn, systems, ctx), part);
    return;
}// registerLagSurfaceForceFunction

void 
cRigidIBFEMethod::registerLagSurfaceForceFunction(
    const LagSurfaceForceFcnData& data,
    const unsigned int part)
{
    TBOX_ASSERT(part < d_num_parts);
    d_lag_surface_force_fcn_data[part] = data;
    d_fcn_systems[part].insert(data.systems.begin(), data.systems.end());
    d_surface_fcn_systems[part].insert(data.systems.begin(), data.systems.end());
    return;
} // registerLagSurfaceForceFunction

void 
cRigidIBFEMethod::setupTagBuffer(
    Array<int>& tag_buffer,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg) const
{
    const int finest_hier_ln = gridding_alg->getMaxLevels() - 1;
    const int tsize = tag_buffer.size();
    tag_buffer.resizeArray(finest_hier_ln);
    for (int i = tsize; i < finest_hier_ln; ++i) tag_buffer[i] = 0;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        const int gcw = d_fe_data_managers[part]->getGhostCellWidth().max();
        const int tag_ln = d_fe_data_managers[part]->getLevelNumber() - 1;
        if (tag_ln >= 0 && tag_ln < finest_hier_ln)
        {
            tag_buffer[tag_ln] = std::max(tag_buffer[tag_ln], gcw);
        }
    }
    for (int ln = finest_hier_ln - 2; ln >= 0; --ln)
    {
        tag_buffer[ln] = std::max(
            tag_buffer[ln],
            tag_buffer[ln + 1] / gridding_alg->getRatioToCoarserLevel(ln + 1).max() + 1);
    }
    return;
}// setupTagBuffer

void
cRigidIBFEMethod::preprocessIntegrateData(
    double current_time, 
    double new_time, 
    int /*num_cycles*/)
{
    d_current_time = current_time;
    d_new_time = new_time;
    d_half_time = current_time + 0.5 * (new_time - current_time);

    // Extract the FE data.
    d_X_systems.resize(d_num_parts);
    d_X_current_vecs.resize(d_num_parts);
    d_X_new_vecs.resize(d_num_parts);
    d_X_half_vecs.resize(d_num_parts);
    d_X_IB_ghost_vecs.resize(d_num_parts);

    d_U_systems.resize(d_num_parts);
    d_U_current_vecs.resize(d_num_parts);
    d_U_new_vecs.resize(d_num_parts);
    d_U_half_vecs.resize(d_num_parts);

    d_L_systems.resize(d_num_parts);
	d_L_current_vecs.resize(d_num_parts);
	d_L_new_vecs.resize(d_num_parts);
	d_L_half_vecs.resize(d_num_parts);
    d_L_IB_ghost_vecs.resize(d_num_parts);

    d_U_constrained_systems.resize(d_num_parts);
    d_U_constrained_current_vecs.resize(d_num_parts);
    d_U_constrained_new_vecs.resize(d_num_parts);
    d_U_constrained_half_vecs.resize(d_num_parts);
	
	// PETSc wrappers.
	d_vL_current.resize(d_num_parts, PETSC_NULL);
	d_vL_new.resize(d_num_parts, PETSC_NULL);
	
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_X_systems[part] = &d_equation_systems[part]->get_system(COORDS_SYSTEM_NAME);
        d_X_current_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_X_systems[part]->current_local_solution.get());
        d_X_new_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_X_current_vecs[part]->clone().release()); // WARNING: must be manually deleted
        d_X_half_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_X_current_vecs[part]->clone().release()); // WARNING: must be manually deleted
        d_X_IB_ghost_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_fe_data_managers[part]->buildGhostedCoordsVector(/*localize_data*/ false));

        d_U_systems[part] = &d_equation_systems[part]->get_system(VELOCITY_SYSTEM_NAME);
        d_U_current_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_U_systems[part]->current_local_solution.get());
        d_U_new_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_U_current_vecs[part]->clone().release()); // WARNING: must be manually deleted
        d_U_half_vecs[part] = dynamic_cast<PetscVector<double>*>(
            d_U_current_vecs[part]->clone().release()); // WARNING: must be manually deleted

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


        // Initialize X^{n+1/2} and X^{n+1} to equal X^{n}, and initialize
        // U^{n+1/2} and U^{n+1} to equal U^{n}.
        d_X_current_vecs[part]->localize(*d_X_half_vecs[part]);
        d_X_current_vecs[part]->localize(*d_X_new_vecs[part]);
        d_U_current_vecs[part]->localize(*d_U_half_vecs[part]);
        d_U_current_vecs[part]->localize(*d_U_new_vecs[part]);
		d_L_current_vecs[part]->localize(*d_L_half_vecs[part]);
		d_L_current_vecs[part]->localize(*d_L_new_vecs[part]);
        d_U_constrained_current_vecs[part]->localize(*d_U_constrained_half_vecs[part]);
        d_U_constrained_current_vecs[part]->localize(*d_U_constrained_new_vecs[part]);
       
		if (!d_solve_rigid_vel[part])
        { 
            d_constrained_velocity_fcn_data[part].comvelfcn(d_current_time,d_trans_vel_current[part],
                d_rot_vel_current[part]);
            d_constrained_velocity_fcn_data[part].comvelfcn(d_half_time,d_trans_vel_half[part],
                d_rot_vel_half[part]);
            d_constrained_velocity_fcn_data[part].comvelfcn(d_new_time,d_trans_vel_new[part],
                d_rot_vel_new[part]);
        }	
    }
	
	VecCreateMultiVec(PETSC_COMM_WORLD, d_num_parts, &d_vL_current[0], &d_mv_L_current);
	VecCreateMultiVec(PETSC_COMM_WORLD, d_num_parts, &d_vL_new[0], &d_mv_L_new);
	
    return;
}// preprocessIntegrateData

void cRigidIBFEMethod::postprocessIntegrateData(
    double /*current_time*/,
    double /*new_time*/,
    int /*num_cycles*/)
{
    for (unsigned part = 0; part < d_num_parts; ++part)
    {
        // Reset time-dependent Lagrangian data.
        *d_X_systems[part]->solution             = *d_X_new_vecs[part];
        *d_U_systems[part]->solution             = *d_U_new_vecs[part];
        *d_L_systems[part]->solution             = *d_L_new_vecs[part];
        *d_U_constrained_systems[part]->solution = *d_U_constrained_new_vecs[part];
	
        // Update the coordinate mapping dX = X - s.
        updateCoordinateMapping(part);

        // Deallocate Lagrangian scratch data.
        delete d_X_new_vecs[part];
        delete d_X_half_vecs[part];
        delete d_U_new_vecs[part];
        delete d_U_half_vecs[part];
		delete d_L_new_vecs[part];
		delete d_L_half_vecs[part];
        delete d_U_constrained_new_vecs[part];
        delete d_U_constrained_half_vecs[part];
   
    }

    d_X_systems.clear();
    d_X_current_vecs.clear();
    d_X_new_vecs.clear();
    d_X_half_vecs.clear();
    d_X_IB_ghost_vecs.clear();

    d_U_systems.clear();
    d_U_current_vecs.clear();
    d_U_new_vecs.clear();
    d_U_half_vecs.clear();

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
    
    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time = std::numeric_limits<double>::quiet_NaN();
    d_half_time = std::numeric_limits<double>::quiet_NaN();
    return;
}// postprocessIntegrateData

void 
cRigidIBFEMethod::interpolateVelocity(
    const int /*u_data_idx*/,
    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*u_synch_scheds*/,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& /*u_ghost_fill_scheds*/,
    const double /*data_time*/)
{
    //intentionally left blank
    return;
}// interpolateVelocity

void
cRigidIBFEMethod::eulerStep(
    const double current_time,
    const double new_time)
{
    const double dt = new_time-current_time;

    // Compute center of mass and moment of inertia of the body at t^n.
    computeCOMandMOIOfStructures(d_center_of_mass_current,d_moment_of_inertia_current,d_X_current_vecs);

    // Fill the rotation matrix of structures with rotation angle 0.5*(W^n)*dt.
    std::vector<Eigen::Matrix3d> rotation_mat(d_num_parts, Eigen::Matrix3d::Zero());
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        for (int i = 0; i < 3; ++i)
			rotation_mat[part](i,i) = 1.0;
    }  
    setRotationMatrix(d_rot_vel_current,rotation_mat,0.5*dt);

    // Rotate the body with current rotational velocity about origin
    // and translate the body to predicted position X^n+1/2.   
    Eigen::Vector3d dr    = Eigen::Vector3d::Zero();
    Eigen::Vector3d Rxdr  = Eigen::Vector3d::Zero();
    for (int part = 0; part < d_num_parts; ++part)
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
cRigidIBFEMethod::midpointStep(
    const double current_time,
    const double new_time)
{
    const double dt = new_time - current_time;
    // Fill the rotation matrix of structures with rotation angle (W^n+1)*dt.
    std::vector<Eigen::Matrix3d> rotation_mat(d_num_parts,Eigen::Matrix3d::Zero());
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        for (int i = 0; i < 3; ++i)
			rotation_mat[part](i,i) = 1.0;
    }  
    setRotationMatrix(d_rot_vel_half,rotation_mat,dt);
    
    // Rotate the body with current rotational velocity about origin
    // and translate the body to predicted position X^n+1/2.   
    Eigen::Vector3d dr    = Eigen::Vector3d::Zero();
    Eigen::Vector3d Rxdr  = Eigen::Vector3d::Zero();
    for (int part = 0; part < d_num_parts; ++part)
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
cRigidIBFEMethod::trapezoidalStep(
    const double /*current_time*/,
    const double /*new_time*/)
{
    TBOX_ERROR("RigidIBFE does not support trapezoidal time-stepping rule for position update."
               << " Only midpoint rule is supported for position update." << std::endl);
    
    return;    
}// trapezoidalStep

void 
cRigidIBFEMethod::computeLagrangianForce(
    double data_time)
{
    //intentionally left blank
    return;
}// computeLagrangianForce

void
cRigidIBFEMethod::spreadForce(
    int /*f_data_idx*/,
    RobinPhysBdryPatchStrategy* /*f_phys_bdry_op*/,
    const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >&
        /*f_prolongation_scheds*/,
    double /*data_time*/)
{
    //intentionally left blank
    return;  
    
}// spreadForce

void
cRigidIBFEMethod::spreadForce(
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
		/*f_prolongation_scheds*/,
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
	
void
cRigidIBFEMethod::computeNetRigidGeneralizedForce(
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
cRigidIBFEMethod::registerLoadBalancer(
    Pointer<LoadBalancer<NDIM> > load_balancer,
    int workload_data_idx)
{
    TBOX_ASSERT(load_balancer);
    d_load_balancer = load_balancer;
    d_workload_idx = workload_data_idx;

    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->registerLoadBalancer(load_balancer, workload_data_idx);
    }
    return;
}// registerLoadBalancer

void 
cRigidIBFEMethod::updateWorkloadEstimates(
    Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    int /*workload_data_idx*/)
{
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->updateWorkloadEstimates();
    }
    return;
}// updateWorkloadEstimates

void
cRigidIBFEMethod::beginDataRedistribution(
    Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    // intentionally left blank
    return;
}// beginDataRedistribution

void
cRigidIBFEMethod::endDataRedistribution(
    Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
    Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    if (d_is_initialized)
    {
        for (unsigned int part = 0; part < d_num_parts; ++part)
        {
            d_fe_data_managers[part]->reinitElementMappings();
        }
    }
    return;
}// endDataRedistribution

void 
cRigidIBFEMethod::initializeLevelData(
    Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    int level_number,
    double init_data_time,
    bool can_be_refined,
    bool initial_time,
    Pointer<BasePatchLevel<NDIM> > old_level,
   bool allocate_data)
{
    const int finest_hier_level = hierarchy->getFinestLevelNumber();
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->setPatchHierarchy(hierarchy);
        d_fe_data_managers[part]->setPatchLevels(0, finest_hier_level);
        d_fe_data_managers[part]->initializeLevelData(hierarchy,
                                                      level_number,
                                                      init_data_time,
                                                      can_be_refined,
                                                      initial_time,
                                                      old_level,
                                                      allocate_data);
        if (d_load_balancer && level_number == d_fe_data_managers[part]->getLevelNumber())
        {
            d_load_balancer->setWorkloadPatchDataIndex(d_workload_idx, level_number);
            d_fe_data_managers[part]->updateWorkloadEstimates(level_number, level_number);
        }
    }
    return;
}// initializeLevelData

void 
cRigidIBFEMethod::resetHierarchyConfiguration(
    Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    int coarsest_level,
    int /*finest_level*/)
{
    const int finest_hier_level = hierarchy->getFinestLevelNumber();
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->setPatchHierarchy(hierarchy);
        d_fe_data_managers[part]->setPatchLevels(0, hierarchy->getFinestLevelNumber());
        d_fe_data_managers[part]->resetHierarchyConfiguration(
            hierarchy, coarsest_level, finest_hier_level);
    }
    return;
}// resetHierarchyConfiguration

void 
cRigidIBFEMethod::applyGradientDetector(
    Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    int level_number,
    double error_data_time,
    int tag_index,
    bool initial_time,
    bool uses_richardson_extrapolation_too)
{
    Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        d_fe_data_managers[part]->applyGradientDetector(hierarchy,
                                                        level_number,
                                                        error_data_time,
                                                        tag_index,
                                                        initial_time,
                                                        uses_richardson_extrapolation_too);
    }
    return;
}// applyGradientDetector

void
cRigidIBFEMethod::registerPreProcessSolveFluidEquationsCallBackFcn(
    preprocessSolveFluidEqn_callbackfcn callback, 
    void* ctx)
{
    d_prefluidsolve_callback_fcns.push_back(callback);
    d_prefluidsolve_callback_fcns_ctx.push_back(ctx);

    return;
}// registerPreProcessSolveFluidEquationsCallBackFcn

void
cRigidIBFEMethod::preprocessSolveFluidEquations(
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

void cRigidIBFEMethod::commonConstructor(
    const std::string& object_name,
    Pointer<Database> input_db,
    const std::vector<libMesh::Mesh*>& meshes,
    int max_level_number,
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
    
    // Resize some arrays.
    d_num_nodes.resize(d_num_parts);
	

    // Set some default values.
    d_output_eul_lambda    = false;
    const bool use_adaptive_quadrature = true;
    const int point_density = 2.0;
    const bool interp_use_consistent_mass_matrix = true;
    d_interp_spec = FEDataManager::InterpSpec("IB_4",
                                              QGAUSS,
                                              INVALID_ORDER,
                                              use_adaptive_quadrature,
                                              point_density,
                                              interp_use_consistent_mass_matrix);
    d_spread_spec = FEDataManager::SpreadSpec(
        "IB_4", QGAUSS, INVALID_ORDER, use_adaptive_quadrature, point_density);
    d_ghosts = 0;
    d_fe_family = LAGRANGE;
    d_fe_order = INVALID_ORDER;
    d_quad_type = QGAUSS;
    d_quad_order = INVALID_ORDER;
    d_use_consistent_mass_matrix = true;
    d_do_log = false;

    // Initialize function data to NULL.
    d_coordinate_mapping_fcn_data.resize(d_num_parts);
    d_lag_body_force_fcn_data.resize(d_num_parts);
    d_lag_surface_force_fcn_data.resize(d_num_parts);
    d_fcn_systems.resize(d_num_parts);
    d_body_fcn_systems.resize(d_num_parts);
    d_surface_fcn_systems.resize(d_num_parts);

    // Determine whether we should use first-order or second-order shape
    // functions for each part of the structure.
    bool mesh_has_first_order_elems = false;
    bool mesh_has_second_order_elems = false;
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        const MeshBase& mesh = *meshes[part];
        MeshBase::const_element_iterator el_it = mesh.elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.elements_end();
        for (; el_it != el_end; ++el_it)
        {
            const Elem* const elem = *el_it;
            mesh_has_first_order_elems =
                mesh_has_first_order_elems || elem->default_order() == FIRST;
            mesh_has_second_order_elems =
                mesh_has_second_order_elems || elem->default_order() == SECOND;
        }
    }
    mesh_has_first_order_elems = SAMRAI_MPI::maxReduction(mesh_has_first_order_elems);
    mesh_has_second_order_elems = SAMRAI_MPI::maxReduction(mesh_has_second_order_elems);
    if ((mesh_has_first_order_elems && mesh_has_second_order_elems) ||
        (!mesh_has_first_order_elems && !mesh_has_second_order_elems))
    {
        TBOX_ERROR(d_object_name
                   << "::IBFEMethod():\n"
                   << "  all parts of FE mesh must contain only FIRST order elements "
                      "or only SECOND order elements" << std::endl);
    }
    if (mesh_has_first_order_elems)
    {
        d_fe_order = FIRST;
        d_quad_order = THIRD;
    }
    if (mesh_has_second_order_elems)
    {
        d_fe_order = SECOND;
        d_quad_order = FIFTH;
    }

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    // Report configuration.
    pout << "\n";
    pout << d_object_name << ": using " << Utility::enum_to_string<Order>(d_fe_order)
         << " order " << Utility::enum_to_string<FEFamily>(d_fe_family)
         << " finite elements.\n";
    pout << "\n";

    // Check the choices for the kernel function.
    if (d_interp_spec.kernel_fcn != d_spread_spec.kernel_fcn)
    {
        pout << "WARNING: different kernel functions are being used for velocity "
                "interpolation and "
                "force spreading.\n"
             << "         recommended usage is to employ the same kernel functions for both "
                "interpolation and spreading.\n";
    }

    // Create the FE data managers that manage mappings between the FE mesh
    // parts and the Cartesian grid.
    d_meshes = meshes;
    d_equation_systems.resize(d_num_parts, NULL);
    d_fe_data_managers.resize(d_num_parts, NULL);
    for (unsigned int part = 0; part < d_num_parts; ++part)
    {
        // Create FE data managers.
        std::ostringstream manager_stream;
        manager_stream << "cRigidIBFEMethod FEDataManager::" << part;
        const std::string& manager_name = manager_stream.str();
        d_fe_data_managers[part] =
            FEDataManager::getManager(manager_name, d_interp_spec, d_spread_spec);
        d_ghosts =
            IntVector<NDIM>::max(d_ghosts, d_fe_data_managers[part]->getGhostCellWidth());

        // Create FE equation systems objects and corresponding variables.
        d_equation_systems[part] = new EquationSystems(*d_meshes[part]);
        EquationSystems* equation_systems = d_equation_systems[part];
        d_fe_data_managers[part]->setEquationSystems(equation_systems, max_level_number - 1);

        d_fe_data_managers[part]->COORDINATES_SYSTEM_NAME = COORDS_SYSTEM_NAME;
        System& X_system = equation_systems->add_system<System>(COORDS_SYSTEM_NAME);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "X_" << d;
            X_system.add_variable(os.str(), d_fe_order, d_fe_family);
        }

        System& dX_system = equation_systems->add_system<System>(COORD_MAPPING_SYSTEM_NAME);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "dX_" << d;
            dX_system.add_variable(os.str(), d_fe_order, d_fe_family);
        }

        System& U_system = equation_systems->add_system<System>(VELOCITY_SYSTEM_NAME);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "U_" << d;
            U_system.add_variable(os.str(), d_fe_order, d_fe_family);
        }
        
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
        
        System& Delta_system = equation_systems->add_system<System>(REGULATOR_SYSTEM_NAME);       
        Delta_system.add_variable("Delta", d_fe_order, d_fe_family);           
    }

    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time = std::numeric_limits<double>::quiet_NaN();
    d_half_time = std::numeric_limits<double>::quiet_NaN();
    
    // Keep track of the initialization state.
    d_fe_data_initialized = false;
    d_is_initialized = false;
    return;
}// commonConstructor

void 
cRigidIBFEMethod::initializeCoordinates(
  const unsigned int part)
{
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    MeshBase& mesh = equation_systems->get_mesh();
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    const bool identity_mapping = !d_coordinate_mapping_fcn_data[part].fcn;
    for (MeshBase::node_iterator it = mesh.local_nodes_begin();
		 it != mesh.local_nodes_end(); ++it)
    {
        Node* n = *it;
        if (n->n_vars(X_sys_num))
        {
            TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM);
            const libMesh::Point& s = *n;
            libMesh::Point X = s;
            if (!identity_mapping)
            {
                d_coordinate_mapping_fcn_data[part].fcn(
                    X, s, d_coordinate_mapping_fcn_data[part].ctx);
            }
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const int dof_index = n->dof_number(X_sys_num, d, 0);
                X_coords.set(dof_index, X(d));
            }
        }
    }
    X_system.get_dof_map().enforce_constraints_exactly(X_system, &X_coords);
    X_coords.close();
    return;
}// initializeCoordinates

void
cRigidIBFEMethod::updateCoordinateMapping(
    const unsigned int part)
{
    EquationSystems* equation_systems = d_fe_data_managers[part]->getEquationSystems();
    MeshBase& mesh = equation_systems->get_mesh();
    System& X_system = equation_systems->get_system(COORDS_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    NumericVector<double>& X_coords = *X_system.solution;
    System& dX_system = equation_systems->get_system(COORD_MAPPING_SYSTEM_NAME);
    const unsigned int dX_sys_num = dX_system.number();
    NumericVector<double>& dX_coords = *dX_system.solution;
    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end();
         ++it)
    {
        Node* n = *it;
        if (n->n_vars(X_sys_num))
        {
            TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM);
            TBOX_ASSERT(n->n_vars(dX_sys_num) == NDIM);
            const libMesh::Point& s = *n;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const int X_dof_index = n->dof_number(X_sys_num, d, 0);
                const int dX_dof_index = n->dof_number(dX_sys_num, d, 0);
                dX_coords.set(dX_dof_index, X_coords(X_dof_index) - s(d));
            }
        }
    }
    dX_coords.close();
    return;
}// updateCoordinateMapping

void
cRigidIBFEMethod::getFromInput(
    Pointer<Database> input_db,
    bool /*is_from_restart*/)
{  
    // Outputing S[lambda] and Lagrange Multiplier
    d_output_eul_lambda    = input_db->getBoolWithDefault("output_eul_lambda", d_output_eul_lambda);
    
    // Interp setting
    if (input_db->isString("interp_delta_fcn"))
        d_interp_spec.kernel_fcn = input_db->getString("interp_delta_fcn");
    else if (input_db->isString("IB_delta_fcn"))
        d_interp_spec.kernel_fcn = input_db->getString("IB_delta_fcn");
    else if (input_db->isString("interp_kernel_fcn"))
        d_interp_spec.kernel_fcn = input_db->getString("interp_kernel_fcn");
    else if (input_db->isString("IB_kernel_fcn"))
        d_interp_spec.kernel_fcn = input_db->getString("IB_kernel_fcn");
    
    if (input_db->isString("interp_quad_type"))
        d_interp_spec.quad_type =
            Utility::string_to_enum<QuadratureType>(input_db->getString("interp_quad_type"));
    else if (input_db->isString("IB_quad_type"))
        d_interp_spec.quad_type =
            Utility::string_to_enum<QuadratureType>(input_db->getString("IB_quad_type"));

    if (input_db->isString("interp_quad_order"))
        d_interp_spec.quad_order =
            Utility::string_to_enum<Order>(input_db->getString("interp_quad_order"));
    else if (input_db->isString("IB_quad_order"))
        d_interp_spec.quad_order =
            Utility::string_to_enum<Order>(input_db->getString("IB_quad_order"));

    if (input_db->isBool("interp_use_adaptive_quadrature"))
        d_interp_spec.use_adaptive_quadrature = input_db->getBool("interp_use_adaptive_quadrature");
    else if (input_db->isBool("IB_use_adaptive_quadrature"))
        d_interp_spec.use_adaptive_quadrature = input_db->getBool("IB_use_adaptive_quadrature");

    if (input_db->isDouble("interp_point_density"))
        d_interp_spec.point_density = input_db->getDouble("interp_point_density");
    else if (input_db->isDouble("IB_point_density"))
        d_interp_spec.point_density = input_db->getDouble("IB_point_density");

    if (input_db->isBool("interp_use_consistent_mass_matrix"))
        d_interp_spec.use_consistent_mass_matrix =
            input_db->getBool("interp_use_consistent_mass_matrix");
    else if (input_db->isBool("IB_use_consistent_mass_matrix"))
        d_interp_spec.use_consistent_mass_matrix =
            input_db->getBool("IB_use_consistent_mass_matrix");
	    
    // Spread setting
    if (input_db->isString("spread_delta_fcn"))
        d_spread_spec.kernel_fcn = input_db->getString("spread_delta_fcn");
    else if (input_db->isString("IB_delta_fcn"))
        d_spread_spec.kernel_fcn = input_db->getString("IB_delta_fcn");
    else if (input_db->isString("spread_kernel_fcn"))
        d_spread_spec.kernel_fcn = input_db->getString("spread_kernel_fcn");
    else if (input_db->isString("IB_kernel_fcn"))
        d_spread_spec.kernel_fcn = input_db->getString("IB_kernel_fcn");

    if (input_db->isString("spread_quad_type"))
        d_spread_spec.quad_type =
            Utility::string_to_enum<QuadratureType>(input_db->getString("spread_quad_type"));
    else if (input_db->isString("IB_quad_type"))
        d_spread_spec.quad_type =
            Utility::string_to_enum<QuadratureType>(input_db->getString("IB_quad_type"));

    if (input_db->isString("spread_quad_order"))
        d_spread_spec.quad_order =
            Utility::string_to_enum<Order>(input_db->getString("spread_quad_order"));
    else if (input_db->isString("IB_quad_order"))
        d_spread_spec.quad_order =
            Utility::string_to_enum<Order>(input_db->getString("IB_quad_order"));

    if (input_db->isBool("spread_use_adaptive_quadrature"))
        d_spread_spec.use_adaptive_quadrature = input_db->getBool("spread_use_adaptive_quadrature");
    else if (input_db->isBool("IB_use_adaptive_quadrature"))
        d_spread_spec.use_adaptive_quadrature = input_db->getBool("IB_use_adaptive_quadrature");

    if (input_db->isDouble("spread_point_density"))
        d_spread_spec.point_density = input_db->getDouble("spread_point_density");
    else if (input_db->isDouble("IB_point_density"))
        d_spread_spec.point_density = input_db->getDouble("IB_point_density");
    
    // Force setting
    if (input_db->isString("quad_type"))
        d_quad_type = Utility::string_to_enum<QuadratureType>(input_db->getString("quad_type"));
    if (input_db->isString("quad_order"))
        d_quad_order = Utility::string_to_enum<Order>(input_db->getString("quad_order"));
    if (input_db->isBool("use_consistent_mass_matrix"))
        d_use_consistent_mass_matrix = input_db->getBool("use_consistent_mass_matrix"); 
    
    // Other settings
    if (input_db->isInteger("min_ghost_cell_width"))
    {
        d_ghosts = input_db->getInteger("min_ghost_cell_width");
    }
    else if (input_db->isDouble("min_ghost_cell_width"))
    {
        d_ghosts = static_cast<int>(std::ceil(input_db->getDouble("min_ghost_cell_width")));
    }
    if (input_db->keyExists("do_log"))
        d_do_log = input_db->getBool("do_log");
    else if (input_db->keyExists("enable_logging"))
        d_do_log = input_db->getBool("enable_logging"); 
    
    return;
}// getFromInput

void
cRigidIBFEMethod::getFromRestart()
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
    
    int ver = db->getInteger("cIBFE_METHOD_VERSION");
    if (ver != cRIGID_IBFE_METHOD_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version."
                                 << std::endl);
    }    
    
    db->getIntegerArray("d_ghosts", d_ghosts, NDIM);
    d_fe_family = Utility::string_to_enum<FEFamily>(db->getString("d_fe_family"));
    d_fe_order = Utility::string_to_enum<Order>(db->getString("d_fe_order"));
    d_quad_type = Utility::string_to_enum<QuadratureType>(db->getString("d_quad_type"));
    d_quad_order = Utility::string_to_enum<Order>(db->getString("d_quad_order"));
    d_use_consistent_mass_matrix = db->getBool("d_use_consistent_mass_matrix");
    
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
cRigidIBFEMethod::computeCOMandMOIOfStructures(
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
        System& system = *d_X_systems[part];
        DofMap& X_dof_map = system.get_dof_map();
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
        System& system = *d_X_systems[part];
        DofMap& X_dof_map = system.get_dof_map();
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
cRigidIBFEMethod::setRotationMatrix(
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
