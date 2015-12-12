// Filename: CIBMethod.cpp
// Created on 21 Apr 2015 by Amneet Bhalla and Bakytzhan Kallemov 
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

#include "IBAMR_config.h"
#include "LocationIndexRobinBcCoefs.h"
#include "ibamr/CIBMethod.h"
#include "ibamr/CIBStandardInitializer.h"
#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/namespaces.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/LSiloDataWriter.h"
#include "ibtk/PETScMultiVec.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/ibtk_utilities.h"

namespace IBAMR
{

extern "C" {

void getEmpiricalMobilityMatrix(const char* kernel_name,
                                const double mu,
                                const double rho,
                                const double dt,
                                const double dx,
                                const double* X,
                                const int n,
                                const bool reset_constants,
                                const double periodic_correction,
                                const double l_domain,
                                double* mm);

void getRPYMobilityMatrix(const char* kernel_name,
                          const double mu,
                          const double dx,
                          const double* X,
                          const int n,
                          const double periodic_correction,
                          double* mm);
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CIBMethod::CIBMethod(const std::string& object_name,
                     Pointer<Database> input_db,
                     const int no_structures,
                     bool register_for_restart)
    : IBMethod(object_name, input_db, register_for_restart), CIBStrategy(no_structures)
{
    // Set some default values
    d_eul_lambda_idx = -1;
    d_output_eul_lambda = false;
    d_lambda_dump_interval = 0;
    d_time_integrator_needs_regrid = false;

    // Resize some arrays.
    d_constrained_velocity_fcns_data.resize(d_num_rigid_parts);
    d_ext_force_torque_fcn_data.resize(d_num_rigid_parts);
    d_struct_lag_idx_range.resize(d_num_rigid_parts);
    d_lambda_filename.resize(d_num_rigid_parts);
    d_reg_filename.resize(d_num_rigid_parts);

    // Initialize object with data read from the input and restart databases.
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db);

    return;

} // cIBMethod

CIBMethod::~CIBMethod()
{

    if ((SAMRAI_MPI::getRank() == 0)&&(d_lambda_dump_interval))
    {
	d_netlambda_stream.close();
	if(d_output_all_lambdas) d_lambda_stream.close();
    }

    return;
} // ~cIBMethod

void CIBMethod::registerConstrainedVelocityFunction(ConstrainedNodalVelocityFcnPtr nodalvelfcn,
                                                    ConstrainedCOMVelocityFcnPtr comvelfcn,
                                                    void* ctx,
                                                    unsigned int part)
{

#if !defined(NDEBUG)
    TBOX_ASSERT(part < d_num_rigid_parts);
#endif
    registerConstrainedVelocityFunction(ConstrainedVelocityFcnsData(nodalvelfcn, comvelfcn, ctx), part);

    return;
} // registerConstrainedVelocityFunction

void CIBMethod::registerConstrainedVelocityFunction(const ConstrainedVelocityFcnsData& data, unsigned int part)
{

#if !defined(NDEBUG)
    TBOX_ASSERT(part < d_num_rigid_parts);
#endif
    d_constrained_velocity_fcns_data[part] = data;

    return;
} // registerConstrainedVelocityFunction

void CIBMethod::registerExternalForceTorqueFunction(ExternalForceTorqueFcnPtr forcetorquefcn,
                                                    void* ctx,
                                                    unsigned int part)
{
    registerExternalForceTorqueFunction(ExternalForceTorqueFcnData(forcetorquefcn, ctx), part);
    return;
} // registerExternalForceTorqueFunction

void CIBMethod::registerExternalForceTorqueFunction(const ExternalForceTorqueFcnData& data, unsigned int part)
{

#if !defined(NDEBUG)
    TBOX_ASSERT(part < d_num_rigid_parts);
#endif

    d_ext_force_torque_fcn_data[part] = data;

    return;
} // registerExternalForceTorqueFunction
void CIBMethod::registerRigidBodyVelocityDeformationFunction(VelocityDeformationFunctionPtr VelDefFun)
{
    d_VelDefFun=VelDefFun;

}//registerVelocityDeformationFunction

int CIBMethod::getStructuresLevelNumber() const
{
    return d_hierarchy->getFinestLevelNumber();

} // getStructuresLevelNumber

int CIBMethod::getStructureHandle(const int lag_idx) const
{
    for (unsigned struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
    {
        const std::pair<int, int>& lag_idx_range = d_struct_lag_idx_range[struct_no];
        if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second) return struct_no;
    }

    return -1;
} // getStructureHandle

void CIBMethod::registerPreProcessSolveFluidEquationsCallBackFcn(preprocessSolveFluidEqn_callbackfcn callback,
                                                                 void* ctx)
{
    d_prefluidsolve_callback_fcns.push_back(callback);
    d_prefluidsolve_callback_fcns_ctx.push_back(ctx);

    return;
} // registerPreProcessSolveFluidEquationsCallBackFcn

void CIBMethod::preprocessSolveFluidEquations(double current_time, double new_time, int cycle_num)
{
    IBMethod::preprocessSolveFluidEquations(current_time, new_time, cycle_num);

    // Call any registered pre-fluid solve callback functions.
    for (unsigned i = 0; i < d_prefluidsolve_callback_fcns.size(); ++i)
    {
        d_prefluidsolve_callback_fcns[i](current_time, new_time, cycle_num, d_prefluidsolve_callback_fcns_ctx[i]);
    }

    return;
} // preprocessSolveFluidEquations

void CIBMethod::registerEulerianVariables()
{
    IBMethod::registerEulerianVariables();

    const IntVector<NDIM> ib_ghosts = getMinimumGhostCellWidth();
    d_eul_lambda_var = new CellVariable<NDIM, double>(d_object_name + "::eul_lambda", NDIM);
    registerVariable(d_eul_lambda_idx, d_eul_lambda_var, ib_ghosts, d_ib_solver->getCurrentContext());

    return;
} // registerEulerianVariables

void CIBMethod::registerEulerianCommunicationAlgorithms()
{

    IBMethod::registerEulerianCommunicationAlgorithms();

    Pointer<RefineAlgorithm<NDIM> > refine_alg_lambda;
    Pointer<RefineOperator<NDIM> > refine_op;
    refine_alg_lambda = new RefineAlgorithm<NDIM>();
    refine_op = NULL;
    refine_alg_lambda->registerRefine(d_eul_lambda_idx, d_eul_lambda_idx, d_eul_lambda_idx, refine_op);
    registerGhostfillRefineAlgorithm(d_object_name + "::eul_lambda", refine_alg_lambda);

    return;
} // registerEulerianCommunicationAlgorithms

void CIBMethod::preprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
#ifdef TIME_REPORT
    clock_t end_t=0, start_med=0;
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0) start_med = clock();
#endif

    IBMethod::preprocessIntegrateData(current_time, new_time, num_cycles);

    // Get data for free and prescribed bodies.
    std::vector<PetscInt> indices;
    std::vector<PetscScalar> U_vec;
    std::vector<PetscScalar> F_vec;
    int counter=0;
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

	if (!SAMRAI_MPI::getRank())
	{
	    if (num_free_dofs)
	    {
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
		
		for (int k = 0; k < s_max_free_dofs; ++k)
		{
		    if (solve_dofs[k])
		    {
			U_vec.push_back(Ur[k]);
			F_vec.push_back(Fr[k]);
		    }else
		    {
			U_vec.push_back(0.0);
			F_vec.push_back(0.0);
		    }
		    indices.push_back(counter++);
		    
		}
	    }
	}
    }
    const int global_size = SAMRAI_MPI::sumReduction(counter);

    VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, global_size, &d_mv_U);
    VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, global_size, &d_mv_F);

    VecSetValues(d_mv_U, counter, &indices[0], &U_vec[0], INSERT_VALUES); 
    VecAssemblyBegin(d_mv_U);
    VecAssemblyEnd(d_mv_U);

    VecSetValues(d_mv_F, counter, &indices[0], &F_vec[0], INSERT_VALUES); 
    VecAssemblyBegin(d_mv_F);
    VecAssemblyEnd(d_mv_F);

    const double start_time = d_ib_solver->getStartTime();
    const bool initial_time = MathUtilities<double>::equalEps(current_time, start_time);
    if (initial_time)
    {
	Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
	const double* const domain_x_lower = grid_geom->getXLower();
	const double* const domain_x_upper = grid_geom->getXUpper();  
	double domain_length[NDIM];
	for (int d = 0; d < NDIM; ++d)
	{
	    domain_length[d] = domain_x_upper[d]-domain_x_lower[d];
	}
	//const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift();
	const int struct_ln = getStructuresLevelNumber();

	for (unsigned struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
	{
	    d_center_of_mass_current[struct_no] = getStandardInitializer()->getInitialCOMStructure(struct_ln,struct_no);

	    for (unsigned int d = 0; d < NDIM; ++d) 
	    { 
		while (d_center_of_mass_current[struct_no][d] < domain_x_lower[d]) d_center_of_mass_current[struct_no][d] += domain_length[d];
		while (d_center_of_mass_current[struct_no][d] >= domain_x_upper[d]) d_center_of_mass_current[struct_no][d] -= domain_length[d];
	    }
	}
	    
    }

    d_quaternion_half = d_quaternion_current; 

    d_time_integrator_needs_regrid = false;

#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0)
    {
	end_t = clock();
	pout<< std::setprecision(4)<<"    CIBMethod:preprocessIntegrateData, CPU time taken for the time step is:"<< double(end_t-start_med)/double(CLOCKS_PER_SEC)<<std::endl;;
    }
#endif

    return;
} // preprocessIntegrateData

void CIBMethod::postprocessIntegrateData(double current_time, double new_time, int num_cycles)
{    
#ifdef TIME_REPORT
    clock_t end_t=0, start_med=0;
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0) start_med = clock();
#endif
    // Compute net rigid generalized force for structures.
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    Pointer<LData> ptr_lagmultpr = d_l_data_manager->getLData("lambda", finest_ln);
    //Vec L_vec = ptr_lagmultpr->getVec();
    
    std::vector<bool> skip_comp;
    skip_comp.resize(d_num_rigid_parts);
    fill(skip_comp.begin(), skip_comp.end(), false);
    //computeNetRigidGeneralizedForce(L_vec, d_net_rigid_generalized_force, skip_comp);
    
    VecDestroy(&d_mv_U);
    VecDestroy(&d_mv_F);

    // Dump Lagrange multiplier data.
    if (d_lambda_dump_interval && ((d_ib_solver->getIntegratorStep() + 1) % d_lambda_dump_interval == 0))
    {
        Pointer<LData> ptr_lagmultpr = d_l_data_manager->getLData("lambda", finest_ln);
        Vec lambda_petsc_vec_parallel = ptr_lagmultpr->getVec();
        Vec lambda_lag_vec_parallel = NULL;
        Vec lambda_lag_vec_seq = NULL;

        VecDuplicate(lambda_petsc_vec_parallel, &lambda_lag_vec_parallel);
        d_l_data_manager->scatterPETScToLagrangian(lambda_petsc_vec_parallel, lambda_lag_vec_parallel, finest_ln);
        d_l_data_manager->scatterToZero(lambda_lag_vec_parallel, lambda_lag_vec_seq);

        if (SAMRAI_MPI::getRank() == 0)
        {
            PetscScalar* L;
            VecGetArray(lambda_lag_vec_seq, &L);
            int counter_L = 0;
            Eigen::Vector3d total_lambda = Eigen::Vector3d::Zero();
            Eigen::Vector3d total_torque = Eigen::Vector3d::Zero();
	    Eigen::Vector3d dr   = Eigen::Vector3d::Zero();
	    Eigen::Vector3d R_dr = Eigen::Vector3d::Zero();
	    Eigen::Vector3d dF   = Eigen::Vector3d::Zero();
	    Eigen::Matrix3d rotation_mat = Eigen::Matrix3d::Identity(3,3);
            for (unsigned int struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
            {
                const int no_ib_pts = getNumberOfNodes(struct_no);
		rotation_mat=d_quaternion_current[struct_no].toRotationMatrix();
                for (int i = 0; i < no_ib_pts; ++i)
                {
		    const IBTK::Point& X = getStandardInitializer()->getPrototypeVertexPosn(finest_ln,struct_no,i);
                    for (int d = 0; d < NDIM; ++d)
                    {
			if (d_output_all_lambdas) d_lambda_stream << L[counter_L] << "\t";
                        total_lambda[d] += L[counter_L];
			dF[d] = L[counter_L];
			dr[d] = X[d];
			counter_L++;
                    }
                    if (d_output_all_lambdas) d_lambda_stream << std::endl;
		    R_dr = rotation_mat*dr; 
		    total_torque += R_dr.cross(dF);

                }
		d_netlambda_stream <<struct_no <<"\t";

                for (int d = 0; d < NDIM; ++d) d_netlambda_stream << total_lambda[d] << "\t";
#if (NDIM == 2)
		d_netlambda_stream << total_torque[2];
#elif (NDIM == 3)
		for (int d = 0; d < NDIM; ++d) d_netlambda_stream << total_torque[d] << "\t";
#endif
                d_netlambda_stream << std::endl;
                total_lambda.setZero();
                total_torque.setZero();
            }
	    if (d_output_all_lambdas) d_lambda_stream << std::endl;
	    d_netlambda_stream << std::endl;
            VecRestoreArray(lambda_lag_vec_seq, &L);
        }
        VecDestroy(&lambda_lag_vec_parallel);
        VecDestroy(&lambda_lag_vec_seq);
    }

    if (d_output_eul_lambda)
    {
        std::cout << " --- CIBMethod::postprocessIntegrateData --- " << std::endl;
        // Prepare the LData to spread
        std::vector<Pointer<LData> > spread_lag_data(finest_ln + 1, Pointer<LData>(NULL)),
            position_lag_data(finest_ln + 1, Pointer<LData>(NULL));

        spread_lag_data[finest_ln] = d_l_data_manager->getLData("lambda", finest_ln);
        
        position_lag_data[finest_ln] = d_l_data_manager->getLData("X", finest_ln);

        // Initialize the S[lambda] variable to zero.
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CellData<NDIM, double> > lambda_data = patch->getPatchData(d_eul_lambda_idx);
                lambda_data->fillAll(0.0);
            }
        }

        d_l_data_manager->spread(d_eul_lambda_idx, spread_lag_data, position_lag_data, /*f_phys_bdry_op*/ NULL);
    }

    // New center of mass translational and rotational velocity becomes
    // current velocity for the next time step.
    d_trans_vel_current = d_trans_vel_new;
    d_rot_vel_current = d_rot_vel_new;

    // Do the base class cleanup here.
    IBMethod::postprocessIntegrateData(current_time, new_time, num_cycles);
#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0)
    {
	end_t = clock();
	pout<< std::setprecision(4)<<"    CIBMethod:postprocessIntegrateData, CPU time taken for the time step is:"<< double(end_t-start_med)/double(CLOCKS_PER_SEC)<<std::endl;;
    }
#endif

    return;
} // postprocessIntegrateData

void CIBMethod::initializeLevelData(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                    int level_number,
                                    double init_data_time,
                                    bool can_be_refined,
                                    bool initial_time,
                                    Pointer<BasePatchLevel<NDIM> > old_level,
                                    bool allocate_data)
{
#ifdef TIME_REPORT
    clock_t end_t=0, start_med=0;
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0) start_med = clock();
#endif
    IBMethod::initializeLevelData(hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level,
                                  allocate_data);

    // Allocate LData corresponding to the Lagrange multiplier.
    if (initial_time && d_l_data_manager->levelContainsLagrangianData(level_number))
    {
        // Set structure index info.
        std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(level_number);
        std::sort(structIDs.begin(), structIDs.end());
        const unsigned structs_on_this_ln = (unsigned)structIDs.size();

        for (unsigned struct_no = 0; struct_no < structs_on_this_ln; ++struct_no)
        {
            d_struct_lag_idx_range[struct_no] =
                d_l_data_manager->getLagrangianStructureIndexRange(structIDs[struct_no], level_number);
        }

        // Create Lagrange multiplier and regularization data.
        Pointer<IBTK::LData> lag_mul_data = d_l_data_manager->createLData("lambda", level_number, NDIM,
                                                                          /*manage_data*/ true);
        Pointer<IBTK::LData> regulator_data = d_l_data_manager->createLData("regulator", level_number, NDIM,
                                                                            /*manage_data*/ true);

        // Initialize the Lagrange multiplier to zero.
        // Specific value of lambda will be assigned from structure specific input file.
        VecSet(lag_mul_data->getVec(), 0.0);

        // Initialize the regulator data with default value of h^3.
        // Specific weights will be assigned from structure specific input file.
        Vec regulator_petsc_vec = regulator_data->getVec();
        VecSet(regulator_petsc_vec, 1.0);

        if (d_silo_writer)
        {
            d_silo_writer->registerVariableData("lambda", lag_mul_data, level_number);
        }
    }

#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0)
    {
	end_t = clock();
	pout<< std::setprecision(4)<<"    CIBMethod:initializeLevelData, CPU time taken for the time step is:"<< double(end_t-start_med)/double(CLOCKS_PER_SEC)<<std::endl;;
    }
#endif

    return;
} // initializeLevelData

void CIBMethod::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                         Pointer<GriddingAlgorithm<NDIM> > gridding_alg,
                                         int u_data_idx,
                                         const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
                                         const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                         int integrator_step,
                                         double init_data_time,
                                         bool initial_time)
{
#ifdef TIME_REPORT
    clock_t end_t=0, start_med=0;
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0) start_med = clock();
#endif
    // Initialize various Lagrangian data objects required by the conventional
    // IB method.
    IBMethod::initializePatchHierarchy(hierarchy, gridding_alg, u_data_idx, u_synch_scheds, u_ghost_fill_scheds,
                                       integrator_step, init_data_time, initial_time);

    // Lookup the range of hierarchy levels.
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

    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        if (d_silo_writer)
        {
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
                Pointer<LData> lag_mul_data = d_l_data_manager->getLData("lambda", ln);
                d_silo_writer->registerVariableData("lambda", lag_mul_data, ln);
            }
        }
    }

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

    // Set lambda and regularization weight from input file.
    if (initial_time)
    {
        setInitialLambda(finest_ln);
        setRegularizationWeight(finest_ln);

	const int struct_ln = getStructuresLevelNumber();
	for(unsigned i=0;i<d_num_rigid_parts;i++)
	{
	    Eigen::Quaterniond* initial_body_quatern = getStandardInitializer()->getStructureQuaternion(struct_ln, i);
	    d_quaternion_current[i] = (*initial_body_quatern);
	}
	
	// //Copy all quarteninons from StandardInitializer as current 
	// for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
	//     for(unsigned i=0;i<d_num_rigid_parts;i++)
	//     {
	// 	Eigen::Quaterniond* initial_body_quatern = getStandardInitializer()->getStructureQuaternion(ln, i);
	// 	(*initial_body_quatern) = Eigen::Quaterniond::Identity();
	//     }
    }

#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0)
    {
	end_t = clock();
	pout<< std::setprecision(4)<<"    CIBMethod:intializePatchHierarchy, CPU time taken for the time step is:"<< double(end_t-start_med)/double(CLOCKS_PER_SEC)<<std::endl;;
    }
#endif

    return;
} // initializePatchHierarchy

void CIBMethod::interpolateVelocity(const int u_data_idx,
                                    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
                                    const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                    const double data_time)
{

    std::cout << " CIBMethod::interpolateVelocity " << std::endl;

    if (d_lag_velvec_is_initialized)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
#endif
        std::vector<Pointer<LData> >* U_half_data, *X_half_data;
        bool* X_half_needs_ghost_fill;
        getVelocityData(&U_half_data, d_half_time);
        getPositionData(&X_half_data, &X_half_needs_ghost_fill, d_half_time);


	if( u_synch_scheds.empty() && u_ghost_fill_scheds.empty() ){ 
	    std::cout << "empty interpolation " << std::endl;
	    // If walls -> fill ghost cells before interpolation.
	    fillGhostCells(u_data_idx, data_time); 
	    d_l_data_manager->interp(u_data_idx, 
	    			   *U_half_data, 
	    			   *X_half_data, 
	    			   getCoarsenSchedules(d_object_name+"::u::CONSERVATIVE_COARSEN"),
	    			   getGhostfillRefineSchedules(d_object_name+"::u"),
	    			   data_time);
	}
	else{
	    d_l_data_manager->interp(u_data_idx, *U_half_data, *X_half_data, u_synch_scheds, u_ghost_fill_scheds, data_time);
	}

        d_lag_velvec_is_initialized = false;
    }

    return;
} // interpolateVelocity

void CIBMethod::spreadForce(
    int f_data_idx,
    IBTK::RobinPhysBdryPatchStrategy* f_phys_bdry_op,
    const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds,
    double data_time)
{
    if (d_constraint_force_is_initialized)
    {

#if !defined(NDEBUG)
        TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
#endif
	if( f_phys_bdry_op==NULL && f_prolongation_scheds.empty() ){ 
	    std::cout << "empty spreading " << (d_u_phys_bdry_op==NULL) << std::endl;
	    IBMethod::spreadForce(f_data_idx, 
	     			  d_u_phys_bdry_op,                               /* f_phys_bdry_op */ /* d_u_phys_bdry_op */
				  getProlongRefineSchedules(d_object_name+"::f"), /* f_prolongation_scheds */
				  data_time);
	}
	else{
	    IBMethod::spreadForce(f_data_idx, f_phys_bdry_op, f_prolongation_scheds, data_time);
	}
        d_constraint_force_is_initialized = false;
    }

    return;
} // spreadForce

void CIBMethod::eulerStep(const double current_time, const double new_time)
{
#ifdef TIME_REPORT
    clock_t end_t=0, start_med=0;
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0) start_med = clock();
#endif
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;
    //   int flag_regrid=0;

    // setup the quaternions of structures with rotation angle 0.5*(W^n)*dt.
     std::vector<Eigen::Matrix3d> rotation_mat(d_num_rigid_parts, Eigen::Matrix3d::Identity(3,3));
    for (unsigned struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
    {
	const double vecnorm = d_rot_vel_current[struct_no].norm();
	if (!MathUtilities<double>::equalEps(vecnorm,0))
	{
	    Eigen::Vector3d rot_vel_axis = d_rot_vel_current[struct_no]/vecnorm;
	    Eigen::Quaterniond q_rot(Eigen::AngleAxisd(vecnorm*0.5*dt, rot_vel_axis));
	    d_quaternion_half[struct_no] = (q_rot.normalized()*d_quaternion_current[struct_no]).normalized();
	}
	else
	{
	    d_quaternion_half[struct_no] = d_quaternion_current[struct_no];
	}
	
	rotation_mat[struct_no]=d_quaternion_half[struct_no].toRotationMatrix();
    }

    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    const double* const domain_x_lower = grid_geom->getXLower();
    const double* const domain_x_upper = grid_geom->getXUpper();  
    double domain_length[NDIM];
    for (int d = 0; d < NDIM; ++d)
    {
	domain_length[d] = domain_x_upper[d]-domain_x_lower[d];
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
            Eigen::Vector3d dr = Eigen::Vector3d::Zero();
            Eigen::Vector3d R_dr = Eigen::Vector3d::Zero();

            int struct_handle = 0;
            if (structs_on_this_ln > 1) struct_handle = getStructureHandle(lag_idx);

            const std::pair<int,int>& lag_idx_range = d_struct_lag_idx_range[struct_handle]; 
	    const IBTK::Point& X = getStandardInitializer()->getPrototypeVertexPosn(ln,struct_handle,lag_idx-lag_idx_range.first);
            for (unsigned int d = 0; d < NDIM; ++d)  dr[d] = X[d];
          
            // Rotate dr vector using the rotation matrix.
            R_dr = rotation_mat[struct_handle]*dr; 

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_half[d] = d_center_of_mass_current[struct_handle][d] + R_dr[d] +
                            0.5 * dt * d_trans_vel_current[struct_handle][d];
		if (periodic_shift[d])
		{
		    while (X_half[d] < domain_x_lower[d]) X_half[d] += domain_length[d];
		    while (X_half[d] >= domain_x_upper[d]) X_half[d] -= domain_length[d];
		}
            }
        }
        (*X_half_data)[ln]->restoreArrays();
    }
    *X_half_needs_ghost_fill = true;

    // Compute the COM at mid-step.
    for (unsigned struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
    {
	for (unsigned int d = 0; d < NDIM; ++d) 
	{ 
	    d_center_of_mass_half[struct_no][d] = d_center_of_mass_current[struct_no][d] + 0.5*dt*d_trans_vel_current[struct_no][d];
	}
    }
    d_X_new_needs_ghost_fill = true;
    d_X_half_needs_reinit = true;

#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0)
    {
	end_t = clock();
	pout<< std::setprecision(4)<<"    CIBMethod:eulerStep, CPU time taken for the time step is:"<< double(end_t-start_med)/double(CLOCKS_PER_SEC)<<std::endl;;
    }
#endif

    return;
} // eulerStep

void CIBMethod::midpointStep(const double current_time, const double new_time)
{
#ifdef TIME_REPORT
    clock_t end_t=0, start_med=0;
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0) start_med = clock();
#endif
    const double dt = new_time - current_time;
    int flag_regrid=0;

    // Fill the rotation matrix of structures with rotation angle (W^n+1)*dt.
    std::vector<Eigen::Matrix3d> rotation_mat(d_num_rigid_parts, Eigen::Matrix3d::Identity(3,3));
    for (unsigned struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
    {
	const double vecnorm = d_rot_vel_half[struct_no].norm();
	if (!MathUtilities<double>::equalEps(vecnorm,0))
	{
	    Eigen::Vector3d rot_vel_axis = d_rot_vel_half[struct_no]/vecnorm;
	    Eigen::Quaterniond q_rot(Eigen::AngleAxisd(vecnorm*dt, rot_vel_axis));
	    d_quaternion_current[struct_no] = (q_rot.normalized()*d_quaternion_current[struct_no]).normalized();
	}
	rotation_mat[struct_no]=d_quaternion_current[struct_no].toRotationMatrix();
    }

    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    const double* const domain_x_lower = grid_geom->getXLower();
    const double* const domain_x_upper = grid_geom->getXUpper();  
    double domain_length[NDIM];
    for (int d = 0; d < NDIM; ++d)
    {
	domain_length[d] = domain_x_upper[d]-domain_x_lower[d];
    }
    const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift();

    // Rotate the body with current rotational velocity about origin
    // and translate the body to newer position.
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        boost::multi_array_ref<double, 2>& X_new_array = *d_X_new_data[ln]->getLocalFormVecArray();
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
            double* const X_new = &X_new_array[local_idx][0];
            Eigen::Vector3d dr = Eigen::Vector3d::Zero();
            Eigen::Vector3d R_dr = Eigen::Vector3d::Zero();

            int struct_handle = 0;
            if (structs_on_this_ln > 1) struct_handle = getStructureHandle(lag_idx);
            const std::pair<int,int>& lag_idx_range = d_struct_lag_idx_range[struct_handle]; 
	    const IBTK::Point& X = getStandardInitializer()->getPrototypeVertexPosn(ln,struct_handle,lag_idx-lag_idx_range.first);
            for (unsigned int d = 0; d < NDIM; ++d)  dr[d] = X[d];
            
            // Rotate dr vector using the rotation matrix.
            R_dr = rotation_mat[struct_handle]*dr; 

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_new[d] =
                    d_center_of_mass_current[struct_handle][d] + R_dr[d] + dt * d_trans_vel_half[struct_handle][d];
		if (periodic_shift[d])
		{
		    while (X_new[d] < domain_x_lower[d]) X_new[d] += domain_length[d];
		    while (X_new[d] >= domain_x_upper[d]) X_new[d] -= domain_length[d];

		    const double X_check = (X_new[d]-dt * d_trans_vel_half[struct_handle][d]);
		    if ((X_check>domain_x_upper[d]) || (X_check < domain_x_lower[d])) flag_regrid = 1;
		}
            }
        }
        d_X_new_data[ln]->restoreArrays();
    }

    flag_regrid = SAMRAI_MPI::sumReduction(flag_regrid);
    if (flag_regrid) d_time_integrator_needs_regrid = true;

    for (unsigned struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
    {
    	for (unsigned int d = 0; d < NDIM; ++d) 
    	{ 
    	    d_center_of_mass_current[struct_no][d] += dt*d_trans_vel_half[struct_no][d]; 
    	    if (periodic_shift[d])
    	    {
    	    	while (d_center_of_mass_current[struct_no][d] < domain_x_lower[d]) d_center_of_mass_current[struct_no][d] += domain_length[d];
    	    	while (d_center_of_mass_current[struct_no][d] >= domain_x_upper[d]) d_center_of_mass_current[struct_no][d] -= domain_length[d];
    	    }
    	}
    }
    d_X_new_needs_ghost_fill = true;
    d_X_half_needs_reinit = true;

#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0)
    {
	end_t = clock();
	pout<< std::setprecision(4)<<"    CIBMethod:midpointStep, CPU time taken for the time step is:"<< double(end_t-start_med)/double(CLOCKS_PER_SEC)<<std::endl;;
    }
#endif

    return;
} // midpointStep

void CIBMethod::trapezoidalStep(const double /*current_time*/, const double /*new_time*/)
{
    TBOX_ERROR("CIBMethod does not support trapezoidal time-stepping rule for position update."
               << " Only mid-point rule is supported." << std::endl);

    return;
} // trapezoidalStep

void CIBMethod::registerVisItDataWriter(Pointer<VisItDataWriter<NDIM> > visit_writer)
{
    d_visit_writer = visit_writer;
    return;
} // registerVisItDataWriter

void CIBMethod::putToDatabase(Pointer<Database> db)
{
    IBMethod::putToDatabase(db);

    for (unsigned int struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
    {
        std::ostringstream U, W;
        U << "U_" << struct_no;
        W << "W_" << struct_no;
        db->putDoubleArray(U.str(), &d_trans_vel_current[struct_no][0], 3);
        db->putDoubleArray(W.str(), &d_rot_vel_current[struct_no][0], 3);
    }

    return;
} // putToDatabase

void CIBMethod::setConstraintForce(Vec L, const double data_time, const double scale)
{

#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
#endif

    const int struct_ln = getStructuresLevelNumber();

    std::vector<Pointer<LData> >* F_half_data;
    bool* F_half_needs_ghost_fill;
    getForceData(&F_half_data, &F_half_needs_ghost_fill, d_half_time);
    Vec F_half = (*F_half_data)[struct_ln]->getVec();
    VecCopy(L, F_half);
    VecScale(F_half, scale);
    *F_half_needs_ghost_fill = true;

    d_constraint_force_is_initialized = true;

    return;
} // setConstraintForce

void CIBMethod::getConstraintForce(Vec* L, const double data_time)
{

#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_current_time) ||
                MathUtilities<double>::equalEps(data_time, d_new_time));
#endif
    const int finest_ln = getStructuresLevelNumber();
    Pointer<LData> ptr_lagmultpr = d_l_data_manager->getLData("lambda", finest_ln);
    Vec lambda = ptr_lagmultpr->getVec();
    *L = lambda;

    return;
} // getConstraintForce

void CIBMethod::getFreeRigidVelocities(Vec* U, const double data_time)
{

#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_current_time) ||
                MathUtilities<double>::equalEps(data_time, d_new_time));
#endif

    *U = d_mv_U;
    return;

} // getFreeRigidVelocities

void CIBMethod::getNetExternalForceTorque(Vec* F, const double data_time)
{

#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_current_time) ||
                MathUtilities<double>::equalEps(data_time, d_new_time));
#endif

    *F = d_mv_F;
    return;

} // getNetExternalForceTorque

void CIBMethod::subtractMeanConstraintForce(Vec L, int f_data_idx, const double scale)
{
    // Temporarily scale the L Vec.
    VecScale(L, scale);

    // Get the underlying array
    PetscScalar* L_array;
    VecGetArray(L, &L_array);
    PetscInt local_size_L;
    VecGetLocalSize(L, &local_size_L);

    double F[NDIM] = { 0.0 };
    const int local_no_ib_pts = local_size_L / NDIM;

    for (int k = 0; k < local_no_ib_pts; ++k)
    {
        for (int d = 0; d < NDIM; ++d)
        {
            F[d] += L_array[k * NDIM + d];
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

    // Unscale the vector.
    VecScale(L, 1.0 / scale);

    return;
} // subtractMeanConstraintForce

void CIBMethod::setInterpolatedVelocityVector(Vec /*V*/, const double data_time)
{

#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
#endif
    d_lag_velvec_is_initialized = true;

    return;
} // setInterpolatedVelocityVector

void CIBMethod::getInterpolatedVelocity(Vec V, const double data_time, const double scale)
{

#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
#endif

    const int ln = getStructuresLevelNumber();
    std::vector<Pointer<LData> >* U_half_data;
    getVelocityData(&U_half_data, d_half_time);
    VecCopy((*U_half_data)[ln]->getVec(), V);
    VecScale(V, scale);

    return;
} // getInterpolatedVelocity

void CIBMethod::computeMobilityRegularization(Vec D, Vec L, const double scale)
{
    const int struct_ln = getStructuresLevelNumber();
    Pointer<LData> reg_data = d_l_data_manager->getLData("regulator", struct_ln);
    Vec W = reg_data->getVec();
    VecPointwiseMult(D, L, W);
    VecScale(D, scale);

    return;
} // computeMobilityRegularization

unsigned int CIBMethod::getNumberOfNodes(const unsigned int struct_no) const
{
    std::pair<int, int> lag_idx_range = d_struct_lag_idx_range[struct_no];
    return (lag_idx_range.second - lag_idx_range.first);

} // getNumberOfStructuresNodes

void CIBMethod::setRigidBodyVelocity(Vec Uvec, Vec V, const std::vector<bool>& skip_comp, const bool isHalfTimeStep)
{
    const unsigned num_procs = SAMRAI_MPI::getNodes();
    const unsigned rank      = SAMRAI_MPI::getRank();
    const int struct_ln = getStructuresLevelNumber();

    Vec V_lag_vec = PETSC_NULL;
    VecDuplicate(V, &V_lag_vec);
    d_l_data_manager->scatterPETScToLagrangian(V,V_lag_vec,struct_ln);

    std::map<unsigned,unsigned> struct_to_free_map;
    unsigned counter=0;
    for (unsigned struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
    {
    	if (!d_isFree_component[struct_no]) continue;
    	struct_to_free_map[struct_no]=counter++;
    }

    Vec U_all;
    VecScatter ctx;
    VecScatterCreateToAll(Uvec,&ctx,&U_all);

    VecScatterBegin(ctx, Uvec, U_all, INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd  (ctx, Uvec, U_all, INSERT_VALUES,SCATTER_FORWARD);

    PetscScalar* u_array = NULL;
    PetscInt comps;
    VecGetSize(Uvec, &comps);
    VecGetArray(U_all, &u_array);

    std::vector<PetscScalar> V_Array; 
    std::vector<PetscInt> indices; 
    counter=0;
    for (unsigned struct_no = rank; struct_no < d_num_rigid_parts; struct_no +=num_procs)
    {
	if (skip_comp[struct_no]) continue;

	RigidDOFVector U;
	U.setZero();
	
	//setting nodal velocity is not working yet
	// if (d_constrained_velocity_fcns_data[struct_id].nodalvelfcn)
	// {
	    
	//     void* ctx = d_constrained_velocity_fcns_data[struct_id].ctx;
	//     CIBMethod* cib_method_ptr = this;
	//     d_constrained_velocity_fcns_data[struct_id].nodalvelfcn(struct_id, V, U, d_X_half_data[struct_ln]->getVec(),
	// 						       d_center_of_mass_half[struct_id], d_new_time, ctx, 
	// 						       cib_method_ptr);
	// }

	Eigen::Vector3d dr;
	Eigen::Vector3d R_dr;
	if (d_isImposed_component[struct_no]) getNewRigidBodyVelocity(struct_no, U); 

	for (int idof=0; idof < s_max_free_dofs; ++idof) 
	{
	    if (d_solve_rigid_vel[struct_no][idof])  
		U[idof] = u_array[struct_to_free_map[struct_no]*s_max_free_dofs+idof];
	}
	const std::pair<int, int>& part_idx_range = d_struct_lag_idx_range[struct_no];
	Eigen::Matrix3d rotation_mat;
	
	if  (isHalfTimeStep)	
	    rotation_mat = d_quaternion_half[struct_no].toRotationMatrix();
	else 
	    rotation_mat = d_quaternion_current[struct_no].toRotationMatrix();

	const unsigned num_nodes = getNumberOfNodes(struct_no);

	for (unsigned inode = 0; inode < num_nodes; ++inode)
	{    
	    const IBTK::Point& X = getStandardInitializer()->getPrototypeVertexPosn(struct_ln,struct_no, inode);
	    for (unsigned int d = 0; d < NDIM; ++d)  
	    {
		dr[d] = X[d];
		indices.push_back((part_idx_range.first+inode)*NDIM+d);
		counter++;
	    }

	    R_dr = rotation_mat*dr; 
#if (NDIM == 2)
	    V_Array.push_back(U[0] - U[2] * R_dr[1]);
	    V_Array.push_back(U[1] + U[2] * R_dr[0]);
#elif(NDIM == 3)
	    V_Array.push_back(U[0] + U[4] * R_dr[2] - U[5] * R_dr[1]);
	    V_Array.push_back(U[1] + U[5] * R_dr[0] - U[3] * R_dr[2]);
	    V_Array.push_back(U[2] + U[3] * R_dr[1] - U[4] * R_dr[0]);
#endif
	}
    }
    
    VecSetValues(V_lag_vec, counter, &indices[0], &V_Array[0], INSERT_VALUES); 
 
    VecAssemblyBegin(V_lag_vec);
    VecAssemblyEnd(V_lag_vec);
    
    d_l_data_manager->scatterLagrangianToPETSc(V_lag_vec, V, struct_ln);

    VecRestoreArray(U_all, &u_array);
    VecScatterDestroy(&ctx);
    VecDestroy(&U_all);
    VecDestroy(&V_lag_vec); 

    return;
} // setRigidBodyVelocity

void CIBMethod::setRigidBodyDeformationVelocity(Vec W)
{
    const int struct_ln = getStructuresLevelNumber();

    if (!d_VelDefFun) return;
    CIBMethod* cib_method = this;
    (*d_VelDefFun)(W, d_X_half_data[struct_ln]->getVec(), d_center_of_mass_half, cib_method);
    return;
}//setRigidBodyDeformationVelocity

void CIBMethod::computeNetRigidGeneralizedForce(Vec L, Vec F, const std::vector<bool>& skip_comp, const bool isHalfTimeStep)
{
    const unsigned num_procs = SAMRAI_MPI::getNodes();
    const unsigned rank      = SAMRAI_MPI::getRank();
    const int struct_ln = getStructuresLevelNumber();

    Vec L_lag_vec     = PETSC_NULL;
    Vec L_lag_vec_seq = PETSC_NULL;   
    VecDuplicate(L, &L_lag_vec);
    d_l_data_manager->scatterPETScToLagrangian(L,L_lag_vec,struct_ln);
    d_l_data_manager->scatterToAll(L_lag_vec, L_lag_vec_seq);

    PetscScalar *L_array;
    VecGetArray(L_lag_vec_seq,&L_array);  

    std::map<unsigned,unsigned> struct_to_free_map;
    unsigned counter=0;
    for (unsigned struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
    {
	if (!d_isFree_component[struct_no]) continue;
	struct_to_free_map[struct_no]=counter++;
    }
    VecSet(F, 0.0);    
    PetscInt size;
    VecGetSize(F,&size);
    //   TBOX_ASSERT((unsigned)size == s_max_free_dofs*counter);

    std::vector<PetscScalar> F_Array; 
    std::vector<PetscInt> indices; 

    counter=0;
    for (unsigned struct_no = rank; struct_no < d_num_rigid_parts; struct_no +=num_procs)
    {
 	if (skip_comp[struct_no]) continue;
	
	Eigen::Matrix3d rotation_mat;
	
	if  (isHalfTimeStep)	
	    rotation_mat = d_quaternion_half[struct_no].toRotationMatrix();
	else 
	    rotation_mat = d_quaternion_current[struct_no].toRotationMatrix();

	Eigen::Vector3d dr;
	Eigen::Vector3d R_dr;
	double f_array[s_max_free_dofs];
	std::fill(f_array, f_array + s_max_free_dofs, 0.0);
	const unsigned num_nodes = getNumberOfNodes(struct_no);
        const std::pair<int, int>& lag_idx_range = d_struct_lag_idx_range[struct_no];

	for (unsigned inode = 0; inode < num_nodes; ++inode)
	{
	    const IBTK::Point& X = getStandardInitializer()->getPrototypeVertexPosn(struct_ln,struct_no, inode);
	    for (unsigned int d = 0; d < NDIM; ++d)  dr[d] = X[d];
	    R_dr = rotation_mat*dr; 

	    PetscScalar *P = L_array + (lag_idx_range.first+inode)*NDIM;
	    for (int d = 0; d < NDIM; ++d)
	    {
		if (d_solve_rigid_vel[struct_no][d]) f_array[d] += P[d];
	    }
#if (NDIM == 2)
	    if (d_solve_rigid_vel[struct_no][2]) f_array[2] += P[1] * R_dr[0] - P[0] * R_dr[1];
#elif(NDIM == 3)
	    if (d_solve_rigid_vel[struct_no][3]) f_array[3] += P[2] * R_dr[1] - P[1] * R_dr[2];
	    if (d_solve_rigid_vel[struct_no][4]) f_array[4] += P[0] * R_dr[2] - P[2] * R_dr[0];
	    if (d_solve_rigid_vel[struct_no][5]) f_array[5] += P[1] * R_dr[0] - P[0] * R_dr[1];
#endif
	}
	
        for (int d = 0; d < s_max_free_dofs; ++d)
	{
	    F_Array.push_back(f_array[d]);
	    indices.push_back(struct_to_free_map[struct_no]*s_max_free_dofs+d);
	    counter++;
	}
    }

    VecSetValues(F, counter, &indices[0], &F_Array[0], INSERT_VALUES); 
    VecAssemblyBegin(F);
    VecAssemblyEnd(F);

    VecRestoreArray(L_lag_vec_seq,&L_array); 
    VecDestroy(&L_lag_vec); 
    VecDestroy(&L_lag_vec_seq);

    return;
} // computeNetRigidGeneralizedForce


void CIBMethod::copyAllVecToArray(Vec b,
				  double* array,
				  const std::vector<unsigned> & all_rhs_struct_ids,
				  const int data_depth)
{
    // Map the Lagrangian indices into PETSc indices
    const int struct_ln = getStructuresLevelNumber();
    const unsigned num_structs= all_rhs_struct_ids.size();

    Vec b_lag_vec_parallel    = PETSC_NULL;
    Vec b_lag_vec_seq         = PETSC_NULL;    
    VecDuplicate(b, &b_lag_vec_parallel);

    d_l_data_manager->scatterPETScToLagrangian(b,b_lag_vec_parallel,struct_ln);
    d_l_data_manager->scatterToAll(b_lag_vec_parallel, b_lag_vec_seq);

    PetscScalar *petsc_array;
    VecGetArray(b_lag_vec_seq,&petsc_array);  

    int counter = 0;
    for (unsigned i = 0; i < num_structs; ++i) 
    {
        const std::pair<int, int>& lag_idx_range = d_struct_lag_idx_range[all_rhs_struct_ids[i]];

	for (int k = lag_idx_range.first; k < lag_idx_range.second; ++k) 	
	{
	    for (int idir=0; idir<data_depth; ++idir) 
	    {
		array[counter++]  = petsc_array[k*data_depth+idir]; 
	    }
	}
    }
    VecRestoreArray(b_lag_vec_seq,&petsc_array); 

    VecDestroy(&b_lag_vec_parallel); 
    VecDestroy(&b_lag_vec_seq);

    return;
} // copyAllVecToArray

void CIBMethod::copyAllArrayToVec(Vec b,
				  double* array,
				  const std::vector<unsigned>& all_rhs_struct_ids,
				  const int data_depth)
{
    const int struct_ln = getStructuresLevelNumber();
    const unsigned num_structs= all_rhs_struct_ids.size();

    Vec b_lag_vec = PETSC_NULL;
    VecDuplicate(b, &b_lag_vec);
    VecSet(b_lag_vec, 0.0);

    PetscInt local_size = 0;
    for (unsigned k = 0; k < num_structs; ++k)
    {
     	local_size += getNumberOfNodes(all_rhs_struct_ids[k]);
    }
    local_size *= data_depth;

    PetscInt *indices = new PetscInt[local_size]; 

    unsigned counter = 0;
    for (unsigned i = 0; i < num_structs; ++i) 
    {
        const std::pair<int, int>& lag_idx_range = d_struct_lag_idx_range[all_rhs_struct_ids[i]];
	const int struct_size = getNumberOfNodes(all_rhs_struct_ids[i])*data_depth;
    	
	for (int k = 0; k < struct_size; ++k) 	
    	{
	    indices[counter++] = lag_idx_range.first*data_depth + k;
	} 
    }

    VecSetValues(b_lag_vec, local_size, indices, array, INSERT_VALUES); 
 
    VecAssemblyBegin(b_lag_vec);
    VecAssemblyEnd(b_lag_vec);
    
    d_l_data_manager->scatterLagrangianToPETSc(b_lag_vec, b, struct_ln);
    
    VecDestroy(&b_lag_vec); 
    delete [] indices;

    return;
} // copyAllArrayToVec

void CIBMethod::constructMobilityMatrix(std::map<std::string, double*>& managed_mat_map,
					std::map<std::string, MobilityMatrixType>& managed_mat_type_map,
					std::map<std::string, std::vector<unsigned> >& managed_mat_prototype_id_map,
					std::map<std::string, unsigned int>& managed_mat_nodes_map,
					std::map<std::string, std::pair<double, double> >& managed_mat_scale_map,
					std::map<std::string, int>& managed_mat_proc_map,
					const double* grid_dx,
					const double* domain_extents,
					const bool initial_time,
					double rho,
					double mu,
					double f_periodic_corr)
{
    const double dt = d_new_time - d_current_time;
    const int struct_ln = getStructuresLevelNumber();
    const char* ib_kernel = d_l_data_manager->getDefaultInterpKernelFunction().c_str();
    const int rank = SAMRAI_MPI::getRank();
    std::vector<double*> W_vector;
    for (std::map<std::string, double*>::iterator it = managed_mat_map.begin(); it != managed_mat_map.end();
	 ++it)
    {
	const std::string& mat_name = it->first;
	const std::vector<unsigned>&  prototype_struct_ids = managed_mat_prototype_id_map[mat_name];
	const int managing_proc = managed_mat_proc_map[mat_name];
	const int size = managed_mat_nodes_map[mat_name] * NDIM;

	// Regularization for  the mobility matrix
	double* WArr =NULL;
	if (rank == managing_proc) WArr =  new double[size];
	Vec W = d_l_data_manager->getLData("regulator", struct_ln)->getVec();
	std::vector<unsigned> struct_ids;
	if (rank == managing_proc) struct_ids =  prototype_struct_ids;
	copyAllVecToArray(W, WArr, struct_ids, /*depth*/ NDIM);
	if (rank == managing_proc) W_vector.push_back(WArr);
    }

    unsigned counter=0;
    unsigned file_counter = 0;

    unsigned managed_mats = (unsigned)managed_mat_map.size();
    static std::vector<bool> read_files(managed_mats, false);

    for (std::map<std::string, double*>::iterator it = managed_mat_map.begin(); it != managed_mat_map.end();
	 ++it)
    {
	const std::string& mat_name = it->first;
	double* mobility_mat = it->second;
	MobilityMatrixType mat_type = managed_mat_type_map[mat_name];
	const std::vector<unsigned>& prototype_struct_ids = managed_mat_prototype_id_map[mat_name];
	const std::pair<double, double>& scale = managed_mat_scale_map[mat_name];
	const int managing_proc = managed_mat_proc_map[mat_name];
	// Get the position data
	if (rank == managing_proc)
	{
	    // Get the size of matrix
	    const int num_nodes = managed_mat_nodes_map[mat_name];
	    const int size = managed_mat_nodes_map[mat_name] * NDIM;
	    
	    if (mat_type == FILE && !read_files[file_counter])
	    {
		// Get matrix from file
		
		//Needs to be implemented

		read_files[file_counter] = true;
	    }
	    else
	    {
		double* XW =  new double[size];
		Eigen::Vector3d dr = Eigen::Vector3d::Zero();
		Eigen::Vector3d R_dr = Eigen::Vector3d::Zero();
		Eigen::Matrix3d rotation_mat = Eigen::Matrix3d::Identity(3,3);
		
		unsigned offset=0;
		for (unsigned i = 0; i < prototype_struct_ids.size(); ++i)
		{
		    unsigned num_struct_nodes = getNumberOfNodes(prototype_struct_ids[i]);
		    if (!initial_time) rotation_mat = d_quaternion_half[prototype_struct_ids[i]].toRotationMatrix();
		    
		    
		    for (unsigned k = 0; k < num_struct_nodes; ++k)
		    {
	

			IBTK::Point X;

			//check if the global dense mobility matrix is used
			if (mat_name=="GLOBAL_DENSE")
			{
			    X=getStandardInitializer()->getInitialVertexPosn(struct_ln,i,k);
			}else
			{
			    X = getStandardInitializer()->getPrototypeVertexPosn(struct_ln,prototype_struct_ids[i],k);
			}

			if (initial_time)
			{
			    for (int d = 0; d < NDIM; ++d)	XW[offset++] = X[d];
			    
			}else
			{
			    for (unsigned int d = 0; d < NDIM; ++d)  dr[d] = X[d];
			    R_dr = rotation_mat*dr; 
			    for (int d = 0; d < NDIM; ++d)	XW[offset++] = d_center_of_mass_half[prototype_struct_ids[i]][d]+R_dr[d];
			}
		    }
		}
	
	
		// Generate mobility matrix
		if (mat_type == RPY)
		{
		    getRPYMobilityMatrix(ib_kernel, mu, grid_dx[0], &XW[0], num_nodes, f_periodic_corr, mobility_mat);
		}
		else if (mat_type == EMPIRICAL)
		{
		    getEmpiricalMobilityMatrix(ib_kernel, mu, rho, dt, grid_dx[0], &XW[0], num_nodes, 0, f_periodic_corr,
					       domain_extents[0], mobility_mat);
		}
		else
		{
		    TBOX_ERROR("CIBMethod::generateMobilityMatrix(): Invalid type of a mobility matrix." << std::endl);
		}
		delete[] XW;
		

	    }
	    double *W = W_vector[counter];
	    for (int i = 0; i < size; ++i)
	    {
		for (int j = 0; j < size; ++j)
		{
		    mobility_mat[i * size + j] *= scale.first;
		    if (i == j)
		    {
			mobility_mat[i * size + j] += scale.second * W[i];
		    }
		}
	    }
	    delete[] W;
	    counter++;
	}//rank
	file_counter++;
    }//it

    return;
} // generateMobilityMatrix

void CIBMethod::constructKinematicMatrix(double* kinematic_mat,
					 const std::vector<unsigned>& prototype_struct_ids,
					 const bool initial_time,
					 const int managing_proc)
{
    const int struct_ln = getStructuresLevelNumber();
    const int rank = SAMRAI_MPI::getRank();

    // Get the size of matrix
    unsigned num_nodes = 0;
    for (unsigned i = 0; i < prototype_struct_ids.size(); ++i)
    {
        num_nodes += getNumberOfNodes(prototype_struct_ids[i]);
    }
    int size = num_nodes * NDIM;
    int bsize = s_max_free_dofs * prototype_struct_ids.size();


    const unsigned shft = s_max_free_dofs;
    // Get the position data
    if (rank == managing_proc)
    {
	std::fill(kinematic_mat, kinematic_mat + bsize*size, 0.0);
	
	Eigen::Vector3d dr = Eigen::Vector3d::Zero();
	Eigen::Vector3d R_dr = Eigen::Vector3d::Zero();
	Eigen::Matrix3d rotation_mat = Eigen::Matrix3d::Identity(3,3);
	unsigned q = 0;
	
	for (unsigned i = 0; i < prototype_struct_ids.size(); ++i)
	{
	    num_nodes = getNumberOfNodes(prototype_struct_ids[i]);
	    if (!initial_time) rotation_mat = d_quaternion_half[prototype_struct_ids[i]].toRotationMatrix();

	    for (unsigned k = 0; k < num_nodes; ++k)
	    {
		const IBTK::Point X = getStandardInitializer()->getPrototypeVertexPosn(struct_ln,prototype_struct_ids[i],k);
		for (unsigned int d = 0; d < NDIM; ++d)  dr[d] = X[d];

		if (!initial_time)
		    R_dr = rotation_mat*dr; 
		else
		    R_dr=dr;

		//column major representation for LAPACK
		         
		kinematic_mat[/*col*/(i*shft  )*size  + /*row*/ q*NDIM    ]      =  1.0; //(1,1)
		kinematic_mat[/*col*/(i*shft  )*size  + /*row*/ q*NDIM + 1]      =  0.0; //(2,1)
		
		kinematic_mat[/*col*/(i*shft+1)*size  + /*row*/ q*NDIM    ]      =  0.0; //(1,2)
		kinematic_mat[/*col*/(i*shft+1)*size  + /*row*/ q*NDIM + 1]      =  1.0; //(2,2)
#if (NDIM == 2)
		kinematic_mat[/*col*/(i*shft  )*size  + /*row*/ q*NDIM + 2]      = -R_dr[1]; //(3,1)
		kinematic_mat[/*col*/(i*shft+1)*size  + /*row*/ q*NDIM + 2]      =  R_dr[0]; //(3,2)
#elif(NDIM == 3)
		kinematic_mat[/*col*/(i*shft+2)*size  + /*row*/ q*NDIM    ]      =  0.0;//(1,3)
		kinematic_mat[/*col*/(i*shft+3)*size  + /*row*/ q*NDIM    ]      =  0.0;//(1,4)
		kinematic_mat[/*col*/(i*shft+4)*size  + /*row*/ q*NDIM    ]      =  R_dr[2];//(1,5)
		kinematic_mat[/*col*/(i*shft+5)*size  + /*row*/ q*NDIM    ]      = -R_dr[1];//(1,6)

		kinematic_mat[/*col*/(i*shft+2)*size  + /*row*/ q*NDIM + 1]      =  0.0;//(2,3)
		kinematic_mat[/*col*/(i*shft+3)*size  + /*row*/ q*NDIM + 1]      = -R_dr[2];//(2,4)
		kinematic_mat[/*col*/(i*shft+4)*size  + /*row*/ q*NDIM + 1]      =  0.0;//(2,5)
		kinematic_mat[/*col*/(i*shft+5)*size  + /*row*/ q*NDIM + 1]      =  R_dr[0];//(2,6)
		
		kinematic_mat[/*col*/(i*shft  )*size  + /*row*/ q*NDIM + 2]      =  0.0;//(3,1)
		kinematic_mat[/*col*/(i*shft+1)*size  + /*row*/ q*NDIM + 2]      =  0.0;//(3,2)
		kinematic_mat[/*col*/(i*shft+2)*size  + /*row*/ q*NDIM + 2]      =  1.0;//(3,3)
		kinematic_mat[/*col*/(i*shft+3)*size  + /*row*/ q*NDIM + 2]      =  R_dr[1];//(3,4)
		kinematic_mat[/*col*/(i*shft+4)*size  + /*row*/ q*NDIM + 2]      = -R_dr[0];//(3,5)
		kinematic_mat[/*col*/(i*shft+5)*size  + /*row*/ q*NDIM + 2]      =  0.0;//(3,6)
#endif
		q++;
	    }
	}
    }
    return;
} // constructKinematicMatrix


void CIBMethod::registerStandardInitializer(Pointer<IBAMR::CIBStandardInitializer> ib_initializer)
{
    d_ib_initializer = ib_initializer;
    d_ib_initializer-> getClonesParameters(d_num_structs_types, d_structs_clones_num);
};


Pointer<IBAMR::CIBStandardInitializer> CIBMethod::getStandardInitializer()
{
    return d_ib_initializer;
};

void CIBMethod::rotateArrayInitalBodyFrame(double* array, 
					   const std::vector<unsigned>& struct_ids,					
					   const bool isTranspose,
					   const int managing_proc,
					   const bool BodyMobility)
{
    if (struct_ids.empty()) return;
    const unsigned num_structs = (unsigned)struct_ids.size();
    
    if (SAMRAI_MPI::getRank()==managing_proc)
    {
	unsigned offset=0;
	for (unsigned k = 0; k < num_structs; ++k)
	{
	    unsigned structID=struct_ids[k];
	    unsigned number_of_nodes = getNumberOfNodes(structID);

	    Eigen::Matrix3d Rot;
	    if (isTranspose) 
		Rot=(d_quaternion_half[structID].toRotationMatrix()).transpose();
	    else 
		Rot=d_quaternion_half[structID].toRotationMatrix();

	    if  (BodyMobility)	 		
	    {
		Eigen::Vector3d F = Eigen::Vector3d::Zero();
		Eigen::Vector3d T = Eigen::Vector3d::Zero();
		Eigen::Vector3d R_F = Eigen::Vector3d::Zero();
		Eigen::Vector3d R_T = Eigen::Vector3d::Zero();

		for (unsigned int d = 0; d < NDIM; ++d)
		{
		    F[d] = array[offset+d];
		}
#if(NDIM == 3)
		for (unsigned d = 0; d < NDIM; ++d)
		{
		    T[d] = array[offset+d+3];;
		}
		R_T = Rot * T;
#endif
		
		R_F = Rot * F;

		for (unsigned int d = 0; d < NDIM; ++d)
		{
		    array[offset+d] = R_F[d];
		}

#if(NDIM == 3)
		for (unsigned d = 0; d < NDIM; ++d)
		{
		    array[offset+d+3]= R_T[d];
		}
#endif

		offset +=s_max_free_dofs;

	    }else
	    {
		Eigen::Vector3d dr = Eigen::Vector3d::Zero();
		Eigen::Vector3d R_dr = Eigen::Vector3d::Zero();

		for (unsigned node = 0; node < number_of_nodes; ++node)
		{
		    for (unsigned int d = 0; d < NDIM; ++d)
		    {
			dr[d] = array[offset+d];
		    }
		    // Rotate dr vector using the rotation matrix.
		    R_dr = Rot * dr;
		    for (unsigned int d = 0; d < NDIM; ++d)
		    {
			array[offset+d] = R_dr[d];
		    }
		    offset +=NDIM;
		}//node
	    }

	}//k
    }
    return;
} // rotateArrayInitalBodyFrame

bool CIBMethod::flagRegrid() const
{
    return d_time_integrator_needs_regrid;
};

/////////////////////////////// PRIVATE //////////////////////////////////////

void CIBMethod::getFromInput(Pointer<Database> input_db)
{
    d_output_eul_lambda = input_db->getBoolWithDefault("output_eul_lambda", d_output_eul_lambda);
    d_lambda_dump_interval = input_db->getIntegerWithDefault("lambda_dump_interval", d_lambda_dump_interval);
    if (d_lambda_dump_interval)
    {
        const bool from_restart = RestartManager::getManager()->isFromRestart();
        std::string dir_name = input_db->getStringWithDefault("lambda_dirname", "./Output");
        if (!from_restart) Utilities::recursiveMkdir(dir_name);

	d_output_all_lambdas = input_db->getBoolWithDefault("output_all_lambdas", false);
	
        if (SAMRAI_MPI::getRank() == 0)
        {
            std::string filename = dir_name + "/" + "lambdas.dat";
            std::string filename2 = dir_name + "/" + "NetLambda.dat";
            if (from_restart)
	    {
		d_netlambda_stream.open(filename2.c_str(), std::ofstream::app | std::ofstream::out);
		if(d_output_all_lambdas) d_lambda_stream.open(filename.c_str(), std::ofstream::app | std::ofstream::out);
	    }
            else
	    {
		d_netlambda_stream.open(filename2.c_str(), std::ofstream::out);
		if(d_output_all_lambdas) d_lambda_stream.open(filename.c_str(), std::ofstream::out);
	    }
	    if (d_output_all_lambdas)
	    {
		d_lambda_stream.precision(15);
		d_lambda_stream << std::scientific;
	    }
            d_netlambda_stream.precision(15);
            d_netlambda_stream << std::scientific;
        }
    }

    // if (input_db->keyExists("lambda_filenames"))
    // {
    //     tbox::Array<std::string> lambda_filenames = input_db->getStringArray("lambda_filenames");
    //     TBOX_ASSERT(lambda_filenames.size() == (int)d_num_rigid_parts);
    //     for (unsigned struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
    //     {
    //         d_lambda_filename[struct_no] = lambda_filenames[struct_no];
    //     }
    // }

    if (input_db->keyExists("weight_filenames"))
    {
        tbox::Array<std::string> weight_filenames = input_db->getStringArray("weight_filenames");
        TBOX_ASSERT(weight_filenames.size() == (int)d_num_rigid_parts);
        for (unsigned struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
        {
            d_reg_filename[struct_no] = weight_filenames[struct_no];
        }
    }

    return;
} // getFromInput

void CIBMethod::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR("CIBMethod::getFromRestart(): Restart database corresponding to "
                   << d_object_name << " not found in restart file." << std::endl);
    }

    for (unsigned int struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
    {
        std::ostringstream U, W;
        U << "U_" << struct_no;
        W << "W_" << struct_no;
        db->getDoubleArray(U.str(), &d_trans_vel_current[struct_no][0], 3);
        db->getDoubleArray(W.str(), &d_rot_vel_current[struct_no][0], 3);
    }

    return;
} // getFromRestart

void CIBMethod::computeInitialCOMOfStructures(std::vector<Eigen::Vector3d>& center_of_mass)
{
    const int struct_ln = getStructuresLevelNumber();

    for (unsigned int struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
    {
        center_of_mass[struct_no].setZero();
    }

    for (unsigned k = 0; k < d_num_rigid_parts; ++k)
    {
	int no_blobs = getNumberOfNodes(k);
	for (int i = 0; i < no_blobs; ++i)
	{
	    const IBTK::Point X = getStandardInitializer()->getInitialVertexPosn(struct_ln,k,i);
	    for (int d = 0; d < NDIM; ++d) center_of_mass[k][d] += X[d];
	}

	center_of_mass[k] /= no_blobs;
    }
    return;
} // calculateCOMandMOIOfStructures

void CIBMethod::computeInitialMOIOfStructures(std::vector<Eigen::Matrix3d>& moment_of_inertia)
{
    const int struct_ln = getStructuresLevelNumber();
    // Zero out the moment of inertia tensor.
    for (unsigned struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
    {
        moment_of_inertia[struct_no].setZero();
    }

    for (unsigned k = 0; k < d_num_rigid_parts; ++k)
    {
	int no_blobs = getNumberOfNodes(k);
	const Eigen::Vector3d X_com  =  Eigen::Vector3d::Zero();
	for (int i = 0; i < no_blobs; ++i)
	{
	    const IBTK::Point& X = getStandardInitializer()->getPrototypeVertexPosn(struct_ln,k,i);
#if (NDIM == 2)
            moment_of_inertia[k](0, 0) += std::pow(X[1] - X_com[1], 2);
            moment_of_inertia[k](0, 1) += -(X[0] - X_com[0]) * (X[1] - X_com[1]);
            moment_of_inertia[k](1, 1) += std::pow(X[0] - X_com[0], 2);
            moment_of_inertia[k](2, 2) += std::pow(X[0] - X_com[0], 2) + std::pow(X[1] - X_com[1], 2);
#endif

#if (NDIM == 3)
            moment_of_inertia[k](0, 0) += std::pow(X[1] - X_com[1], 2) + std::pow(X[2] - X_com[2], 2);
            moment_of_inertia[k](0, 1) += -(X[0] - X_com[0]) * (X[1] - X_com[1]);
            moment_of_inertia[k](0, 2) += -(X[0] - X_com[0]) * (X[2] - X_com[2]);
            moment_of_inertia[k](1, 1) += std::pow(X[0] - X_com[0], 2) + std::pow(X[2] - X_com[2], 2);
            moment_of_inertia[k](1, 2) += -(X[1] - X_com[1]) * (X[2] - X_com[2]);
            moment_of_inertia[k](2, 2) += std::pow(X[0] - X_com[0], 2) + std::pow(X[1] - X_com[1], 2);
#endif
        }
    }

    // Fill-in symmetric part of inertia tensor.
    for (unsigned int struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
    {
        moment_of_inertia[struct_no](1, 0) = moment_of_inertia[struct_no](0, 1);
        moment_of_inertia[struct_no](2, 0) = moment_of_inertia[struct_no](0, 2);
        moment_of_inertia[struct_no](2, 1) = moment_of_inertia[struct_no](1, 2);
    }

    return;
} // calculateCOMandMOIOfStructures

void CIBMethod::setRegularizationWeight(const int level_number)
{
    Pointer<LData> reg_data = d_l_data_manager->getLData("regulator", level_number);
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    const double* const dx = grid_geom->getDx();
#if (NDIM == 2)
    const double cell_volume = dx[0] * dx[1];
#elif(NDIM == 3)
    const double cell_volume = dx[0] * dx[1] * dx[2];
#endif

    boost::multi_array_ref<double, 2>& reg_data_array = *reg_data->getLocalFormVecArray();
    const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(level_number);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

    // Get structures on this level.
    const std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(level_number);
    const unsigned structs_on_this_ln = (unsigned)structIDs.size();
#if !defined(NDEBUG)
    TBOX_ASSERT(structs_on_this_ln == d_num_rigid_parts);
#endif

    for (unsigned struct_no = 0; struct_no < structs_on_this_ln; ++struct_no)
    {
        const std::pair<int, int>& lag_idx_range = d_struct_lag_idx_range[struct_no];
        if (d_reg_filename[struct_no].empty())
        {
            for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
            {
                const LNode* const node_idx = *cit;
                const int lag_idx = node_idx->getLagrangianIndex();
                if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
                {
                    const int local_idx = node_idx->getLocalPETScIndex();
                    double* const W = &reg_data_array[local_idx][0];
                    for (unsigned int d = 0; d < NDIM; ++d) W[d] = cell_volume;
                }
            }
            continue;
        }

        // Read weights from file and set it.
        std::ifstream reg_filestream(d_reg_filename[struct_no].c_str(), std::ifstream::in);
        if (!reg_filestream.is_open())
        {
            TBOX_ERROR("CIBMethod::setRegularizationWeight()"
                       << "could not open file" << d_reg_filename[struct_no] << std::endl);
        }

        std::string line_f;
        int lag_pts = -1;
        if (std::getline(reg_filestream, line_f))
        {
            std::istringstream iss(line_f);
            iss >> lag_pts;
            if (lag_pts != (lag_idx_range.second - lag_idx_range.first))
            {
                TBOX_ERROR("CIBMethod::setRegularizationWeight() Total no. of Lagrangian points in the weight file "
                           << d_reg_filename[struct_no] << " not equal to corresponding vertex file." << std::endl);
            }
        }
        else
        {
            TBOX_ERROR("CIBMethod::setRegularizationWeight() Error in the input regularization file "
                       << d_reg_filename[struct_no] << " at line number 0. Total number of Lagrangian  points required."
                       << std::endl);
        }

        std::vector<double> reg_weight(lag_pts);
        for (int k = 0; k < lag_pts; ++k)
        {
            if (std::getline(reg_filestream, line_f))
            {
                std::istringstream iss(line_f);
                iss >> reg_weight[k];
            }
            else
            {
                TBOX_ERROR("CIBMethod::setRegularizationWeight() Error in the input regularization file "
                           << d_reg_filename[struct_no] << " at line number " << k + 1 << std::endl);
            }
        }

        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int lag_idx = node_idx->getLagrangianIndex();
            if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
            {
                const int local_idx = node_idx->getLocalPETScIndex();
                double* const W = &reg_data_array[local_idx][0];
                const double& weight = reg_weight[lag_idx - lag_idx_range.first];

                // For zero weight we do not use any regularization
                if (!MathUtilities<double>::equalEps(weight, 0.0))
                {
                    for (unsigned int d = 0; d < NDIM; ++d)
                        W[d] = cell_volume / reg_weight[lag_idx - lag_idx_range.first];
                }
                else
                {
                    for (unsigned int d = 0; d < NDIM; ++d) W[d] = 0.0;
                }
            }
        }
    }
    reg_data->restoreArrays();

    return;
} // setRegularizationWeight

void CIBMethod::setInitialLambda(const int level_number)
{

    Pointer<IBTK::LData> lambda_data = d_l_data_manager->getLData("lambda", level_number);
    boost::multi_array_ref<double, 2>& lambda_data_array = *lambda_data->getLocalFormVecArray();
    const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(level_number);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

    // Get structures on this level.
    const std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(level_number);
    const unsigned structs_on_this_ln = (unsigned)structIDs.size();
#if !defined(NDEBUG)
    TBOX_ASSERT(structs_on_this_ln == d_num_rigid_parts);
#endif

    for (unsigned struct_no = 0; struct_no < structs_on_this_ln; ++struct_no)
    {
        const std::pair<int, int> lag_idx_range = d_struct_lag_idx_range[struct_no];
        if (d_lambda_filename[struct_no].empty()) continue;

        std::ifstream lambda_filestream(d_lambda_filename[struct_no].c_str(), std::ifstream::in);
        if (!lambda_filestream.is_open())
        {
            TBOX_ERROR("CIBMethod::setInitialLambda()"
                       << "could not open file" << d_lambda_filename[struct_no] << std::endl);
        }

        std::string line_f;
        int lag_pts = -1;
        if (std::getline(lambda_filestream, line_f))
        {
            std::istringstream iss(line_f);
            iss >> lag_pts;

            if (lag_pts != (lag_idx_range.second - lag_idx_range.first))
            {
                TBOX_ERROR("CIBMethod::setInitialLambda() Total no. of Lagrangian points in the lambda file "
                           << d_lambda_filename[struct_no] << " not equal to corresponding vertex file." << std::endl);
            }
        }
        else
        {
            TBOX_ERROR("CIBMethod::::setInitialLambda() Error in the input lambda file "
                       << d_lambda_filename[struct_no] << " at line number 0. Total number of Lag pts. required."
                       << std::endl);
        }

        std::vector<double> initial_lambda(lag_pts * NDIM);
        for (int k = 0; k < lag_pts; ++k)
        {
            if (std::getline(lambda_filestream, line_f))
            {
                std::istringstream iss(line_f);
                for (int d = 0; d < NDIM; ++d) iss >> initial_lambda[k * NDIM + d];
            }
            else
            {
                TBOX_ERROR("CIBMethod::setInitialLambda() Error in the input lambda file "
                           << d_lambda_filename[struct_no] << " at line number " << k + 1 << std::endl);
            }
        }

        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int lag_idx = node_idx->getLagrangianIndex();
            if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
            {
                const int local_idx = node_idx->getLocalPETScIndex();
                double* const L = &lambda_data_array[local_idx][0];
                for (unsigned int d = 0; d < NDIM; ++d)
                    L[d] = initial_lambda[(lag_idx - lag_idx_range.first) * NDIM + d];
            }
        }
    }
    lambda_data->restoreArrays();

    return;
} // setInitialLambda


void CIBMethod::setVelocityBC(std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> *u_bc_coefs){
  
    d_u_bc_coefs = u_bc_coefs;
    
    return;
} // setVelocityBC

void CIBMethod::setVelocityPhysBdryOp(IBTK::RobinPhysBdryPatchStrategy* u_phys_bdry_op){
  
    d_u_phys_bdry_op = u_phys_bdry_op;
    return;
}

void CIBMethod::fillGhostCells(int in, const double time){

    typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps;

    InterpolationTransactionComponent x_component(in, 
						  /*DATA_REFINE_TYPE*/ "NONE",           /* NONE */
						  /*USE_CF_INTERPOLATION*/ true,         /* true */
						  /*DATA_COARSEN_TYPE*/ "CUBIC_COARSEN", /* Cubic_coarsen */
						  /*BDRY_EXTRAP_TYPE*/ "LINEAR",         /* Linear */
						  /*CONSISTENT_TYPE_2_BDRY*/ false,      /* false */
						  d_u_bc_coefs[0],                       /* d_u_bc_coefs[0] */
						  /*d_fill_pattern*/ NULL);              /* Null */

    transaction_comps.push_back(x_component);
    
    // Setup 
    int d_coarsest_ln = 0;
    int d_finest_ln   = d_hierarchy->getFinestLevelNumber();   

    const bool homogeneous_bc = d_homogeneous_bc;
                                
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_hier_bdry_fill; 
    d_hier_bdry_fill = new IBTK::HierarchyGhostCellInterpolation(); 
    d_hier_bdry_fill->setHomogeneousBc(homogeneous_bc); 
    d_hier_bdry_fill->initializeOperatorState(transaction_comps, d_hierarchy, d_coarsest_ln, d_finest_ln);

    // Fill ghost cells of x
    d_hier_bdry_fill->fillData(time);

    //Deallocate data
    d_hier_bdry_fill->deallocateOperatorState();
    d_hier_bdry_fill.setNull();
    transaction_comps.clear();

    return;
} // fillGhostCells


} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
