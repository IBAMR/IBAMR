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
    // intentionally left blank.
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

    IBMethod::preprocessIntegrateData(current_time, new_time, num_cycles);

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
	
	computeInitialCOMOfStructures(d_center_of_mass_initial);
	d_center_of_mass_current=d_center_of_mass_initial;
	for (unsigned struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
	    for (unsigned int d = 0; d < NDIM; ++d) 
	    { 
		while (d_center_of_mass_current[struct_no][d] < domain_x_lower[d]) d_center_of_mass_current[struct_no][d] += domain_length[d];
		while (d_center_of_mass_current[struct_no][d] >= domain_x_upper[d]) d_center_of_mass_current[struct_no][d] -= domain_length[d];
	    }
	    
    	d_quaternion_half = d_quaternion_current; 
    }
    return;
} // preprocessIntegrateData

void CIBMethod::postprocessIntegrateData(double current_time, double new_time, int num_cycles)
{    
    // Compute net rigid generalized force for structures.
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    Pointer<LData> ptr_lagmultpr = d_l_data_manager->getLData("lambda", finest_ln);
    Vec L_vec = ptr_lagmultpr->getVec();
    for (unsigned int part = 0; part < d_num_rigid_parts; ++part)
    {
        computeNetRigidGeneralizedForce(part, L_vec, d_net_rigid_generalized_force[part]);
    }

    // Destroy the free-parts PETSc Vecs.
    for (size_t k = 0; k < d_U.size(); ++k)
    {
        VecDestroy(&d_U[k]);
        VecDestroy(&d_F[k]);
    }
    VecDestroy(&d_mv_U);
    VecDestroy(&d_mv_F);
    d_U.clear();
    d_F.clear();

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
            int counter_L = -1;
            Eigen::Vector3d total_lambda = Eigen::Vector3d::Zero();

            d_lambda_stream << new_time << std::endl << std::endl;
            for (unsigned int struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
            {
                const int no_ib_pts = getNumberOfNodes(struct_no);
                d_lambda_stream << "structure: " << struct_no << " ib_pts: " << no_ib_pts << std::endl;

                for (int i = 0; i < no_ib_pts; ++i)
                {
                    for (int d = 0; d < NDIM; ++d)
                    {
                        d_lambda_stream << L[++counter_L] << "\t";
                        total_lambda[d] += L[counter_L];
                    }
                    d_lambda_stream << std::endl;
                }
                d_lambda_stream << "Net resultant lambda for structure: " << struct_no << " ";

                for (int d = 0; d < NDIM; ++d) d_lambda_stream << total_lambda[d] << "\t";
                d_lambda_stream << std::endl;
                total_lambda.setZero();
            }
            VecRestoreArray(lambda_lag_vec_seq, &L);
        }
        VecDestroy(&lambda_lag_vec_parallel);
        VecDestroy(&lambda_lag_vec_seq);
    }

    if (d_output_eul_lambda)
    {
        // Prepare the LData to spread
        std::vector<Pointer<LData> > spread_lag_data(finest_ln + 1, Pointer<LData>(NULL)),
            position_lag_data(finest_ln + 1, Pointer<LData>(NULL));

        spread_lag_data[finest_ln] = d_l_data_manager->getLData("lambda", finest_ln);
        ;
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
	
	//Copy all quarteninons from StandardInitializer as current 
	for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
	    for(unsigned i=0;i<d_num_rigid_parts;i++)
	    {
		Eigen::Quaterniond* initial_body_quatern = getStandardInitializer()->getStructureQuaternion(ln, i);
		(*initial_body_quatern) = Eigen::Quaterniond::Identity();
	    }
    }
    return;
} // initializePatchHierarchy

void CIBMethod::interpolateVelocity(const int u_data_idx,
                                    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_synch_scheds,
                                    const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
                                    const double data_time)
{
    if (d_lag_velvec_is_initialized)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(MathUtilities<double>::equalEps(data_time, d_half_time));
#endif
        std::vector<Pointer<LData> >* U_half_data, *X_half_data;
        bool* X_half_needs_ghost_fill;
        getVelocityData(&U_half_data, d_half_time);
        getPositionData(&X_half_data, &X_half_needs_ghost_fill, d_half_time);
        d_l_data_manager->interp(u_data_idx, *U_half_data, *X_half_data, u_synch_scheds, u_ghost_fill_scheds,
                                 data_time);

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
        IBMethod::spreadForce(f_data_idx, f_phys_bdry_op, f_prolongation_scheds, data_time);
        d_constraint_force_is_initialized = false;
    }

    return;
} // spreadForce

void CIBMethod::eulerStep(const double current_time, const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;

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
	    const IBTK::Point X = getStandardInitializer()->getVertexPosn(ln,struct_handle,lag_idx-lag_idx_range.first);
            for (unsigned int d = 0; d < NDIM; ++d)  
                dr[d] = X[d] - d_center_of_mass_initial[struct_handle][d];
          
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
    return;
} // eulerStep

void CIBMethod::midpointStep(const double current_time, const double new_time)
{

    const double dt = new_time - current_time;

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
	    const IBTK::Point X = getStandardInitializer()->getVertexPosn(ln,struct_handle,lag_idx-lag_idx_range.first);
            for (unsigned int d = 0; d < NDIM; ++d)  
                dr[d] = X[d] - d_center_of_mass_initial[struct_handle][d];
            
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
		}
            }
        }
        d_X_new_data[ln]->restoreArrays();
    }
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

void CIBMethod::setRigidBodyVelocity(const unsigned int part, const RigidDOFVector& U, Vec V)
{
    void* ctx = d_constrained_velocity_fcns_data[part].ctx;
    const int struct_ln = getStructuresLevelNumber();

    if (d_constrained_velocity_fcns_data[part].nodalvelfcn)
    {
	CIBMethod* cib_method_ptr = this;
        d_constrained_velocity_fcns_data[part].nodalvelfcn(part, V, U, d_X_half_data[struct_ln]->getVec(),
                                                           d_center_of_mass_half[part], d_new_time, ctx, 
							   cib_method_ptr);
    }
    else
    {
	
	Eigen::Matrix3d rotation_mat=d_quaternion_half[part].toRotationMatrix();
    
        // Wrap the PETSc V into LData
        std::vector<int> nonlocal_indices;
        LData V_data("V", V, nonlocal_indices, false);

        boost::multi_array_ref<double, 2>& V_data_array = *V_data.getLocalFormVecArray();

        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(struct_ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
        const std::pair<int, int>& part_idx_range = d_struct_lag_idx_range[part];

        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int lag_idx = node_idx->getLagrangianIndex();
            Eigen::Vector3d dr;
            Eigen::Vector3d R_dr;


            if (part_idx_range.first <= lag_idx && lag_idx < part_idx_range.second)
            {
                const int local_idx = node_idx->getLocalPETScIndex();
                double* const V_node = &V_data_array[local_idx][0];

		const IBTK::Point X = getStandardInitializer()->getVertexPosn(struct_ln,part,lag_idx-part_idx_range.first);
		for (unsigned int d = 0; d < NDIM; ++d)  dr[d] = X[d] - d_center_of_mass_initial[part][d];
		
		R_dr = rotation_mat*dr; 

#if (NDIM == 2)
                V_node[0] = U[0] - U[2] * R_dr[1];
                V_node[1] = U[1] + U[2] * R_dr[0];
#elif(NDIM == 3)
                V_node[0] = U[0] + U[4] * R_dr[2] - U[5] * R_dr[1];
                V_node[1] = U[1] + U[5] * R_dr[0] - U[3] * R_dr[2];
                V_node[2] = U[2] + U[3] * R_dr[1] - U[4] * R_dr[0];
#endif
            }
        }

        // Restore underlying arrays
        V_data.restoreArrays();
    }

    return;
} // setRigidBodyVelocity

void CIBMethod::computeNetRigidGeneralizedForce(const unsigned int part, Vec L, RigidDOFVector& F)
{
    const int struct_ln = getStructuresLevelNumber();

    // Wrap the distributed PETSc Vec L into LData
    std::vector<int> nonlocal_indices;
    LData p_data("P", L, nonlocal_indices, false);
    const boost::multi_array_ref<double, 2>& p_data_array = *p_data.getLocalFormVecArray();

    F.setZero();
    const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(struct_ln);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
        const int lag_idx = node_idx->getLagrangianIndex();
        const int local_idx = node_idx->getLocalPETScIndex();
        const double* const P = &p_data_array[local_idx][0];
        const unsigned struct_id = getStructureHandle(lag_idx);
        if (struct_id != part) continue;

	Eigen::Matrix3d rotation_mat=d_quaternion_half[part].toRotationMatrix();
	Eigen::Vector3d dr;
	Eigen::Vector3d R_dr;
        const std::pair<int, int>& part_idx_range = d_struct_lag_idx_range[part];
	const IBTK::Point X = getStandardInitializer()->getVertexPosn(struct_ln,part,lag_idx-part_idx_range.first);
	for (unsigned int d = 0; d < NDIM; ++d)  dr[d] = X[d] - d_center_of_mass_initial[part][d];
	R_dr = rotation_mat*dr; 

#if (NDIM == 2)
        for (int d = 0; d < NDIM; ++d)
        {
            F[d] += P[d];
        }
        F[2] += P[1] * R_dr[0] - P[0] * R_dr[1];
#elif(NDIM == 3)
        for (int d = 0; d < NDIM; ++d)
        {
            F[d] += P[d];
        }
        F[3] += P[2] * R_dr[1] - P[1] * R_dr[2];
        F[4] += P[0] * R_dr[2] - P[2] * R_dr[0];
        F[5] += P[1] * R_dr[0] - P[0] * R_dr[1];
#endif
    }
    SAMRAI_MPI::sumReduction(&F[0], NDIM * (NDIM + 1) / 2);
    p_data.restoreArrays();

    return;
} // computeNetRigidGeneralizedForce

void CIBMethod::copyVecToArray(Vec b,
                               double* array,
                               const std::vector<unsigned int>& struct_ids,
                               const int data_depth,
                               const int array_rank)
{
    if (struct_ids.empty()) return;
    const unsigned num_structs = (unsigned)struct_ids.size();

    // Get the Lagrangian indices of the structures.
    std::vector<int> map;
    PetscInt total_nodes = 0;
    for (unsigned k = 0; k < num_structs; ++k)
    {
        total_nodes += getNumberOfNodes(struct_ids[k]);
    }
    map.reserve(total_nodes);
    for (unsigned k = 0; k < num_structs; ++k)
    {
        const std::pair<int, int>& lag_idx_range = d_struct_lag_idx_range[struct_ids[k]];
        const unsigned struct_nodes = getNumberOfNodes(struct_ids[k]);
        for (unsigned j = 0; j < struct_nodes; ++j)
        {
            map.push_back(lag_idx_range.first + j);
        }
    }

    // Map the Lagrangian indices into PETSc indices
    const int struct_ln = getStructuresLevelNumber();
    d_l_data_manager->mapLagrangianToPETSc(map, struct_ln);

    // Wrap the raw data in a PETSc Vec
    PetscInt size = total_nodes * data_depth;
    int rank = SAMRAI_MPI::getRank();
    PetscInt array_local_size = 0;
    if (rank == array_rank) array_local_size = size;
    Vec array_vec;
    VecCreateMPIWithArray(PETSC_COMM_WORLD, /*blocksize*/ 1, array_local_size, PETSC_DECIDE, array, &array_vec);

    // Create index sets to define global index mapping.
    std::vector<PetscInt> vec_indices, array_indices;
    vec_indices.reserve(size);
    array_indices.reserve(size);
    for (PetscInt j = 0; j < total_nodes; ++j)
    {
        PetscInt petsc_idx = map[j];
        for (int d = 0; d < data_depth; ++d)
        {
            array_indices.push_back(j * data_depth + d);
            vec_indices.push_back(petsc_idx * data_depth + d);
        }
    }
    IS is_vec;
    IS is_array;
    ISCreateGeneral(PETSC_COMM_SELF, size, &vec_indices[0], PETSC_COPY_VALUES, &is_vec);
    ISCreateGeneral(PETSC_COMM_SELF, size, &array_indices[0], PETSC_COPY_VALUES, &is_array);

    // Scatter values
    VecScatter ctx;
    VecScatterCreate(b, is_vec, array_vec, is_array, &ctx);
    VecScatterBegin(ctx, b, array_vec, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, b, array_vec, INSERT_VALUES, SCATTER_FORWARD);

    // Cleanup temporary objects.
    VecScatterDestroy(&ctx);
    ISDestroy(&is_vec);
    ISDestroy(&is_array);
    VecDestroy(&array_vec);

    return;
} // copyVecToArray

void CIBMethod::copyArrayToVec(Vec b,
                               double* array,
                               const std::vector<unsigned>& struct_ids,
                               const int data_depth,
                               const int array_rank)
{
    if (struct_ids.empty()) return;
    const unsigned num_structs = (unsigned)struct_ids.size();

    // Get the Lagrangian indices of the structures.
    std::vector<int> map;
    PetscInt total_nodes = 0;
    for (unsigned k = 0; k < num_structs; ++k)
    {
        total_nodes += getNumberOfNodes(struct_ids[k]);
    }
    map.reserve(total_nodes);
    for (unsigned k = 0; k < num_structs; ++k)
    {
        const std::pair<int, int>& lag_idx_range = d_struct_lag_idx_range[struct_ids[k]];
        const unsigned struct_nodes = getNumberOfNodes(struct_ids[k]);
        for (unsigned j = 0; j < struct_nodes; ++j)
        {
            map.push_back(lag_idx_range.first + j);
        }
    }

    // Map the Lagrangian indices into PETSc indices
    const int struct_ln = getStructuresLevelNumber();
    d_l_data_manager->mapLagrangianToPETSc(map, struct_ln);

    // Wrap the array in a PETSc Vec
    PetscInt size = total_nodes * data_depth;
    int rank = SAMRAI_MPI::getRank();
    PetscInt array_local_size = 0;
    if (rank == array_rank) array_local_size = size;
    Vec array_vec;
    VecCreateMPIWithArray(PETSC_COMM_WORLD, /*blocksize*/ 1, array_local_size, PETSC_DECIDE, array, &array_vec);

    // Create index sets to define global index mapping.
    std::vector<PetscInt> vec_indices, array_indices;
    vec_indices.reserve(size);
    array_indices.reserve(size);
    for (PetscInt j = 0; j < total_nodes; ++j)
    {
        PetscInt petsc_idx = map[j];
        for (int d = 0; d < data_depth; ++d)
        {
            array_indices.push_back(j * data_depth + d);
            vec_indices.push_back(petsc_idx * data_depth + d);
        }
    }
    IS is_vec;
    IS is_array;
    ISCreateGeneral(PETSC_COMM_SELF, size, &vec_indices[0], PETSC_COPY_VALUES, &is_vec);
    ISCreateGeneral(PETSC_COMM_SELF, size, &array_indices[0], PETSC_COPY_VALUES, &is_array);

    // Scatter values
    VecScatter ctx;
    VecScatterCreate(array_vec, is_array, b, is_vec, &ctx);
    VecScatterBegin(ctx, array_vec, b, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, array_vec, b, INSERT_VALUES, SCATTER_FORWARD);

    // Destroy temporary objects
    VecScatterDestroy(&ctx);
    ISDestroy(&is_vec);
    ISDestroy(&is_array);

    VecDestroy(&array_vec);
    VecDestroy(&array_vec);

    return;
} // copyArrayToVec

void CIBMethod::constructMobilityMatrix(const std::string& /*mat_name*/,
                                        MobilityMatrixType mat_type,
                                        double* mobility_mat,
                                        const std::vector<unsigned>& prototype_struct_ids,
                                        const double* grid_dx,
                                        const double* domain_extents,
                                        const bool initial_time,
                                        double rho,
                                        double mu,
                                        const std::pair<double, double>& scale,
                                        double f_periodic_corr,
                                        const int managing_proc)
{
    const double dt = d_new_time - d_current_time;
    const int struct_ln = getStructuresLevelNumber();
    const char* ib_kernel = d_l_data_manager->getDefaultInterpKernelFunction().c_str();
    const int rank = SAMRAI_MPI::getRank();

    // Get the size of matrix
    unsigned num_nodes = 0;
    for (unsigned i = 0; i < prototype_struct_ids.size(); ++i)
    {
        num_nodes += getNumberOfNodes(prototype_struct_ids[i]);
    }
    const int size = num_nodes * NDIM;



    // Get the position data
    double* XW = NULL;
    if (rank == managing_proc) XW = new double[size];
    if (initial_time)
    {
	if (rank == managing_proc)
	{
	    unsigned offset=0;
	    for (unsigned i = 0; i < prototype_struct_ids.size(); ++i)
	    {
		num_nodes = getNumberOfNodes(prototype_struct_ids[i]);
		for (unsigned k = 0; k < num_nodes; ++k)
		{
		    const IBTK::Point X = getStandardInitializer()->getVertexPosn(struct_ln,prototype_struct_ids[i],k);
		    for (int d = 0; d < NDIM; ++d)	XW[offset++] = X[d];
		}
	    }
	}
    }
    else
    {
	if (rank == managing_proc)
	{
	    unsigned offset=0;
	    Eigen::Vector3d dr;
	    Eigen::Vector3d R_dr;
	    
	    for (unsigned i = 0; i < prototype_struct_ids.size(); ++i)
	    {
		num_nodes = getNumberOfNodes(prototype_struct_ids[i]);
		Eigen::Matrix3d rotation_mat = d_quaternion_half[prototype_struct_ids[i]].toRotationMatrix();
		
		for (unsigned k = 0; k < num_nodes; ++k)
		{
		    const IBTK::Point X = getStandardInitializer()->getVertexPosn(struct_ln,prototype_struct_ids[i],k);
		    for (unsigned int d = 0; d < NDIM; ++d)  dr[d] = X[d] - d_center_of_mass_initial[prototype_struct_ids[i]][d];
		    R_dr = rotation_mat*dr; 
		    for (int d = 0; d < NDIM; ++d)	XW[offset++] = d_center_of_mass_half[prototype_struct_ids[i]][d]+R_dr[d];
		}
	    }
	}
    }
    
    
    // Generate mobility matrix
    if (rank == managing_proc)
    {
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
    }

    // Regularize the mobility matrix
    Vec W = d_l_data_manager->getLData("regulator", struct_ln)->getVec();
    copyVecToArray(W, XW, prototype_struct_ids, /*depth*/ NDIM, managing_proc);
    if (rank == managing_proc)
    {
        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j < size; ++j)
            {
                mobility_mat[i * size + j] *= scale.first;
                if (i == j)
                {
                    mobility_mat[i * size + j] += scale.second * XW[i];
                }
            }
        }
        delete[] XW;
    }

    return;
} // generateMobilityMatrix

void CIBMethod::registerStandardInitializer(Pointer<IBAMR::CIBStandardInitializer> ib_initializer)
{
//    registerLInitStrategy(ib_initializer);
    d_ib_initializer = ib_initializer;
};


Pointer<IBAMR::CIBStandardInitializer> CIBMethod::getStandardInitializer()
{
    return d_ib_initializer;
};

void CIBMethod::rotateArrayInitalBodyFrame(double* array, 
					   const std::vector<unsigned>& struct_ids,					
					   const bool isTranspose,
					   const int managing_proc)
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

	    Eigen::Vector3d dr = Eigen::Vector3d::Zero();
	    Eigen::Vector3d R_dr = Eigen::Vector3d::Zero();
	    Eigen::Matrix3d Rot;
	    if (isTranspose) 
		Rot=(d_quaternion_half[structID].toRotationMatrix()).transpose();
	    else 
		Rot=d_quaternion_half[structID].toRotationMatrix();
	 		
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
	}//k
    }
    return;
} // rotateArrayInitalBodyFrame

/////////////////////////////// PRIVATE //////////////////////////////////////

void CIBMethod::getFromInput(Pointer<Database> input_db)
{
    d_output_eul_lambda = input_db->getBoolWithDefault("output_eul_lambda", d_output_eul_lambda);
    d_lambda_dump_interval = input_db->getIntegerWithDefault("lambda_dump_interval", d_lambda_dump_interval);
    if (d_lambda_dump_interval)
    {
        const bool from_restart = RestartManager::getManager()->isFromRestart();
        std::string dir_name = input_db->getStringWithDefault("lambda_dirname", "./lambda");
        if (!from_restart) Utilities::recursiveMkdir(dir_name);

        if (SAMRAI_MPI::getRank() == 0)
        {
            std::string filename = dir_name + "/" + "lambda";
            if (from_restart)
                d_lambda_stream.open(filename.c_str(), std::ofstream::app | std::ofstream::out);
            else
                d_lambda_stream.open(filename.c_str(), std::ofstream::out);

            d_lambda_stream.precision(16);
            d_lambda_stream << std::scientific;
        }
    }

    if (input_db->keyExists("lambda_filenames"))
    {
        tbox::Array<std::string> lambda_filenames = input_db->getStringArray("lambda_filenames");
        TBOX_ASSERT(lambda_filenames.size() == (int)d_num_rigid_parts);
        for (unsigned struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
        {
            d_lambda_filename[struct_no] = lambda_filenames[struct_no];
        }
    }

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
	    const IBTK::Point X = getStandardInitializer()->getVertexPosn(struct_ln,k,i);
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
    std::vector<Eigen::Vector3d> center_of_mass_initial;
    center_of_mass_initial.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());
    computeInitialCOMOfStructures(center_of_mass_initial);

    for (unsigned k = 0; k < d_num_rigid_parts; ++k)
    {
	int no_blobs = getNumberOfNodes(k);
	const Eigen::Vector3d& X_com  =  center_of_mass_initial[k];
	for (int i = 0; i < no_blobs; ++i)
	{
	    const IBTK::Point X = getStandardInitializer()->getVertexPosn(struct_ln,k,i);
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

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
