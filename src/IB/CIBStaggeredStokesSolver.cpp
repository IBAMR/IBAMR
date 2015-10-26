// Filename: CIBStaggeredStokesSolver.cpp
// Created on 10 Nov 2014 by Amneet Bhalla
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
#include "ibamr/CIBSaddlePointSolver.h"
#include "ibamr/CIBStaggeredStokesSolver.h"
#include "ibamr/CIBStrategy.h"
#include "ibamr/INSStaggeredHierarchyIntegrator.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/PETScMultiVec.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "tbox/Database.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

CIBStaggeredStokesSolver::CIBStaggeredStokesSolver(const std::string& object_name,
                                                   Pointer<Database> input_db,
                                                   Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator,
                                                   Pointer<CIBStrategy> cib_strategy,
                                                   const std::string& default_options_prefix)
    : StaggeredStokesSolver(), d_cib_strategy(cib_strategy, false),
      d_num_rigid_parts(d_cib_strategy->getNumberOfRigidStructures())
{
    GeneralSolver::init(object_name, /*homogeneous bcs*/ false);

    d_sp_solver = NULL;
    d_wide_u_var = NULL;
    d_wide_f_var = NULL;
    d_wide_ctx = NULL;
    d_wide_u_idx = -1;
    d_wide_f_idx = -1;
    d_x_wide = NULL;
    d_b_wide = NULL;
    d_is_initialized = false;
    d_reinitializing_solver = false;

    // Create the saddle-point solver for solving constraint problem.
    d_sp_solver = new CIBSaddlePointSolver(object_name, input_db, navier_stokes_integrator, d_cib_strategy,
                                           default_options_prefix);

    // Create widened variables for IB operations.
    Pointer<IBStrategy> ib_method_ops = d_cib_strategy;
    const IntVector<NDIM> ghost_width = ib_method_ops->getMinimumGhostCellWidth();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_wide_u_var = new SideVariable<NDIM, double>(d_object_name + "::wide_u_var", 1);
    d_wide_f_var = new SideVariable<NDIM, double>(d_object_name + "::wide_f_var", 1);
    d_wide_ctx = var_db->getContext(object_name + "::wide_ctx");
    d_wide_u_idx = var_db->registerVariableAndContext(d_wide_u_var, d_wide_ctx, ghost_width);
    d_wide_f_idx = var_db->registerVariableAndContext(d_wide_f_var, d_wide_ctx, ghost_width);

    std::string dir_name = input_db->getStringWithDefault("free_velocity_dirname", "./Output");	
    Utilities::recursiveMkdir(dir_name);
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    
    if (!SAMRAI_MPI::getRank())
    {
	const std::string fname = dir_name +"/U_free.dat";
		
	if (from_restart)
	    d_U_body_out.open(fname.c_str(), std::ios::out | std::ios::app);
	else	
	    d_U_body_out.open(fname.c_str(), std::ios::out );
	d_U_body_out<<std::scientific;	
	d_U_body_out.precision(15);
    }
    d_U_dump_interval = input_db->getIntegerWithDefault("free_velocity_dump_interval", 1);
    return;
} // CIBStaggeredStokesSolver

CIBStaggeredStokesSolver::~CIBStaggeredStokesSolver()
{
    d_reinitializing_solver = false;
    deallocateSolverState();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    var_db->removePatchDataIndex(d_wide_u_idx);
    var_db->removePatchDataIndex(d_wide_f_idx);

    if (!SAMRAI_MPI::getRank()) d_U_body_out.close();
    return;
} // ~CIBStaggeredStokesSolver()

void CIBStaggeredStokesSolver::setSolutionTime(double solution_time)
{
    GeneralSolver::setSolutionTime(solution_time);
    d_sp_solver->setSolutionTime(solution_time);

    return;
} // setSolutionTime

void CIBStaggeredStokesSolver::setTimeInterval(double current_time, double new_time)
{
    GeneralSolver::setTimeInterval(current_time, new_time);
    d_sp_solver->setTimeInterval(current_time, new_time);

    return;
} // setTimeInterval

Pointer<CIBSaddlePointSolver> CIBStaggeredStokesSolver::getSaddlePointSolver() const
{
    return d_sp_solver;
} // getSaddlePointSolver

void CIBStaggeredStokesSolver::initializeSolverState(const SAMRAIVectorReal<NDIM, double>& x,
                                                     const SAMRAIVectorReal<NDIM, double>& b)
{
    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized)
    {
        d_reinitializing_solver = true;
        deallocateSolverState();
    }

    // Wrap Eulerian data into PETSc Vecs.
    Pointer<PatchHierarchy<NDIM> > hierarchy = x.getPatchHierarchy();
    const int coarsest_ln = x.getCoarsestLevelNumber();
    const int finest_ln = x.getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_wide_u_idx)) level->allocatePatchData(d_wide_u_idx);
        if (!level->checkAllocated(d_wide_f_idx)) level->allocatePatchData(d_wide_f_idx);
    }

    Pointer<CellVariable<NDIM, double> > x_p_cc_var = x.getComponentVariable(1);
    Pointer<CellVariable<NDIM, double> > b_p_cc_var = b.getComponentVariable(1);
    const int x_p_idx = x.getComponentDescriptorIndex(1);
    const int b_p_idx = b.getComponentDescriptorIndex(1);

    d_x_wide = new SAMRAIVectorReal<NDIM, double>(x.getName() + "_wide_x", hierarchy, coarsest_ln, finest_ln);
    d_b_wide = new SAMRAIVectorReal<NDIM, double>(b.getName() + "_wide_b", hierarchy, coarsest_ln, finest_ln);

    d_x_wide->addComponent(d_wide_u_var, d_wide_u_idx, x.getControlVolumeIndex(0));
    d_x_wide->addComponent(x_p_cc_var, x_p_idx, x.getControlVolumeIndex(1));
    d_b_wide->addComponent(d_wide_f_var, d_wide_f_idx, b.getControlVolumeIndex(0));
    d_b_wide->addComponent(b_p_cc_var, b_p_idx, b.getControlVolumeIndex(1));

    // Wrap SAMRAI vector into PETSc Vec
    Vec u_p = PETScSAMRAIVectorReal::createPETScVector(d_x_wide);
    Vec g_h = PETScSAMRAIVectorReal::createPETScVector(d_b_wide);

    // Get the Lagrange multiplier that maintains the rigidity constraint.
    // NOTE: The current time corresponds to the time at which solver is initialized
    // which maybe different from the current time of the timestep being integrated upon.
    Vec L;
    d_cib_strategy->getConstraintForce(&L, d_current_time);

    // Create a vector of the type imposed velocity at the material/nodal points for RHS.
    // NOTE: In the initialization stage we do not need an actual velocity vector, a reference
    // to L should suffice to know the required structure.
    Vec V = L;

    // Get the rigid body velocities that need to be solved for.
    // NOTE: The current time corresponds to the time at which solver is initialized
    // which maybe different from the current time of the timestep being integrated upon.
    Vec U;
    d_cib_strategy->getFreeRigidVelocities(&U, d_current_time);

    // Create a vector that contains net external force and torque on the body.
    // NOTE: In the initialization stage we do not need an actual force vector, a reference
    // to U should suffice to know the required structure.
    Vec F = U;

    // Create the composite vectors.
    Vec mv_x, mv_b;
    std::vector<Vec> vx(3), vb(3);
    vx[0] = u_p;
    vx[1] = L;
    vx[2] = U;
    vb[0] = g_h;
    vb[1] = V;
    vb[2] = F;
    VecCreateMultiVec(PETSC_COMM_WORLD, 3, &vx[0], &mv_x);
    VecCreateMultiVec(PETSC_COMM_WORLD, 3, &vb[0], &mv_b);

    // Initialize the saddle-point solver.
    d_sp_solver->initializeSolverState(mv_x, mv_b);

    // Destroy the temporay vectors.
    PETScSAMRAIVectorReal::destroyPETScVector(u_p);
    PETScSAMRAIVectorReal::destroyPETScVector(g_h);
    VecDestroy(&mv_x);
    VecDestroy(&mv_b);

    d_is_initialized = true;
    d_reinitializing_solver = false;

    return;
} // initializeSolverState

void CIBStaggeredStokesSolver::setVelocityPoissonSpecifications(const PoissonSpecifications& u_problem_coefs)
{
    d_sp_solver->setVelocityPoissonSpecifications(u_problem_coefs);
    return;
} // setVelocityPoissonSpecifications

void CIBStaggeredStokesSolver::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                                  RobinBcCoefStrategy<NDIM>* p_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(u_bc_coefs.size() == NDIM);
#endif

    d_sp_solver->setPhysicalBcCoefs(u_bc_coefs, p_bc_coef);

    return;
} // setPhysicalBcCoefs

void CIBStaggeredStokesSolver::setPhysicalBoundaryHelper(Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_helper);
#endif
    d_sp_solver->setPhysicalBoundaryHelper(bc_helper);
    return;
} // setPhysicalBoundaryHelper

bool CIBStaggeredStokesSolver::solveSystem(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& b)
{

#ifdef TIME_REPORT
    clock_t end_t=0, start_med=0;
    if (SAMRAI_MPI::getRank() == 0) start_med = clock();
#endif

    // Create packaged vectors for the Saddle point solver.
    d_x_wide->copyVector(Pointer<SAMRAIVectorReal<NDIM, double> >(&x, false));
    d_b_wide->copyVector(Pointer<SAMRAIVectorReal<NDIM, double> >(&b, false));

    // Wrap SAMRAI vector into PETSc Vec
    Vec u_p = PETScSAMRAIVectorReal::createPETScVector(d_x_wide);
    Vec g_h = PETScSAMRAIVectorReal::createPETScVector(d_b_wide);

    // Get the Lagrange multiplier that maintains the rigidity constraint.
    // NOTE: We need L at new time to solve for it.
    Vec L;
    d_cib_strategy->getConstraintForce(&L, d_new_time);

    // Get the free body velocities to solve for.
    // NOTE: We need U at new time in the solver.
    Vec U;
    d_cib_strategy->getFreeRigidVelocities(&U, d_new_time);


     // Set the imposed velocity for all bodies in the RHS.
    Vec V,W;
    VecDuplicate(L, &V);
    VecDuplicate(L, &W);
    VecSet(V, 0.0);
    VecSet(W, 0.0);

    d_cib_strategy->setRigidBodyDeformationVelocity(W);
    d_cib_strategy->setRigidBodyVelocity(U, V, false, true, false); //set only imposed parts and zero for free components
    VecAYPX(V,-1.0,W);
    const double interp_scale = d_sp_solver->getInterpScale();
    VecScale(V, interp_scale);
    VecDestroy(&W);
   
    // Get the net external force and torque on the bodies.
    Vec F;
    d_cib_strategy->getNetExternalForceTorque(&F, d_new_time);

    // Create multivector to pass it to the saddle point solver.
    std::vector<Vec> vx(3), vb(3);
    vx[0] = u_p;
    vx[1] = L;
    vx[2] = U;
    vb[0] = g_h;
    vb[1] = V;
    vb[2] = F;

    Vec mv_x, mv_b;
    VecCreateMultiVec(PETSC_COMM_WORLD, 3, &vx[0], &mv_x);
    VecCreateMultiVec(PETSC_COMM_WORLD, 3, &vb[0], &mv_b);

#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0)
    {
	end_t = clock();
	pout<< std::setprecision(4)<<"  StaggeredStokesSolver: generating rhs and initial values, CPU time taken for the time step is:"<< double(end_t-start_med)/double(CLOCKS_PER_SEC)<<std::endl;;
	if (SAMRAI_MPI::getRank() == 0) start_med = clock();
    }
#endif

    // Solve for velocity, pressure and Lagrange multipliers.
    // Notice that the state of d_U is maintained, and is passed
    // as an initial guess for the next solve. We do not need to do
    // anything special for its initial guess to the Krylov solver.
    bool converged = d_sp_solver->solveSystem(mv_x, mv_b);

#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0)
    {
	end_t = clock();
	pout<< std::setprecision(4)<<"  StaggeredStokesSolver: Total time for system solve, CPU time taken for the time step is:"<< double(end_t-start_med)/double(CLOCKS_PER_SEC)<<std::endl;;
	if (SAMRAI_MPI::getRank() == 0) start_med = clock();
    }
#endif

    // Extract solution.
    x.copyVector(d_x_wide);

    d_cib_strategy->updateNewRigidBodyVelocity(U);

    // for (unsigned part = 0, free_part = 0; part < d_num_rigid_parts; ++part)
    // {
    //     RigidDOFVector U_part;
    //     d_cib_strategy->getNewRigidBodyVelocity(part, U_part);

    //     // Update free DOFs.
    //     int num_free_dofs = 0;
    //     const FreeRigidDOFVector& solve_dofs = d_cib_strategy->getSolveRigidBodyVelocity(part, num_free_dofs);
    //     if (num_free_dofs)
    //     {
    //         Vec U_sub;
    //         VecMultiVecGetSubVec(U, free_part, &U_sub);

    //         PetscInt s;
    //         PetscScalar* a = NULL;
    //         VecGetArray(U_sub, &a);
    //         VecGetSize(U_sub, &s);

    //         std::vector<double> a_vec(s, 0.0);
    //         if (a != NULL)
    //         {
    //             std::copy(&a[0], &a[s], &a_vec[0]);
    //         }
    //         SAMRAI_MPI::sumReduction(&a_vec[0], s);
    //         VecRestoreArray(U_sub, &a);

    //         for (int k = 0, p = 0; k < s_max_free_dofs; ++k)
    //         {
    //             if (solve_dofs[k])
    //             {
    //                 U_part[k] = a_vec[p];
    //                 ++p;
    //             }
    //         }
    //         ++free_part;
    //     }
    //     d_cib_strategy->updateNewRigidBodyVelocity(part, U_part);
    // }

    double half_time = 0.5 * (d_new_time + d_current_time);
    pout << "\n"
         << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n"
         << "Interpolating velocity on structure at time  " << half_time << "....\n" << std::endl;

    RefineAlgorithm<NDIM> ghost_fill_alg;
    ghost_fill_alg.registerRefine(d_wide_u_idx, d_wide_u_idx, d_wide_u_idx, NULL);
    Pointer<PatchHierarchy<NDIM> > hierarchy = x.getPatchHierarchy();
    Pointer<RefineSchedule<NDIM> > ghost_fill_schd = ghost_fill_alg.createSchedule(hierarchy->getPatchLevel(0));
    ghost_fill_schd->fillData(half_time);
    d_cib_strategy->setInterpolatedVelocityVector(V, half_time);

// #if 0
//         Pointer<CIBFEMethod> ib_method_ops = d_cib_strategy;
//         bool cached_compute_L2_projection = ib_method_ops->setComputeVelL2Projection(true);
//         ib_method_ops->interpolateVelocity(d_wide_u_idx, std::vector<Pointer<CoarsenSchedule<NDIM> > >(),
//                                            std::vector<Pointer<RefineSchedule<NDIM> > >(), half_time);
//         ib_method_ops->setComputeVelL2Projection(cached_compute_L2_projection);
//         d_cib_strategy->getInterpolatedVelocity(V, half_time);
//         Vec* vV;
//         VecMultiVecGetSubVecs(V, &vV);
//         VecView(vV[0], PETSC_VIEWER_STDOUT_WORLD);
// #endif

// #if 0
//     Pointer<IBStrategy> ib_method_ops = d_cib_strategy;
//     ib_method_ops->interpolateVelocity(d_wide_u_idx, std::vector<Pointer<CoarsenSchedule<NDIM> > >(),
//                                        std::vector<Pointer<RefineSchedule<NDIM> > >(), half_time);
//     d_cib_strategy->getInterpolatedVelocity(V, half_time);
//     VecView(V, PETSC_VIEWER_STDOUT_WORLD);
// #endif

    {
	//free parts velocity output
	static unsigned output_counter=0;

	PetscInt free_comps;
	VecGetSize(vx[2], &free_comps);

	Vec U_all;
	VecScatter ctx;
	VecScatterCreateToAll(vx[2],&ctx,&U_all);
	VecScatterBegin(ctx, vx[2], U_all, INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd  (ctx, vx[2], U_all, INSERT_VALUES,SCATTER_FORWARD);
	
	PetscScalar* u_array = NULL;
	VecGetArray(U_all, &u_array);

	if (free_comps && d_U_dump_interval && (output_counter % d_U_dump_interval==0))
	{
	    if (SAMRAI_MPI::getRank()== 0)
	    {
		for (int part = 0; part < free_comps/s_max_free_dofs; ++part)
		{
		    d_U_body_out <<part<<"\t";	
		    for (int d = 0; d < s_max_free_dofs; ++d)
		    {  
			d_U_body_out << u_array[part*s_max_free_dofs +d] << "\t";
		    }
		    d_U_body_out <<std::endl;	
		}
		d_U_body_out <<std::endl;
	    }
	}
	output_counter++;

	VecRestoreArray(U_all, &u_array);
	VecScatterDestroy(&ctx);
	VecDestroy(&U_all);

	// Delete PETSc vectors.
	PETScSAMRAIVectorReal::destroyPETScVector(u_p);
	PETScSAMRAIVectorReal::destroyPETScVector(g_h);
	VecDestroy(&V);
	VecDestroy(&mv_x);
	VecDestroy(&mv_b);
    }
#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0)
    {
	end_t = clock();
	pout<< std::setprecision(4)<<"  StaggeredStokesSolver: Update velocities and output, CPU time taken for the time step is:"<< double(end_t-start_med)/double(CLOCKS_PER_SEC)<<std::endl;;
	pout<<std::endl;
    }
#endif

    return converged;

} // solveSystem

void CIBStaggeredStokesSolver::deallocateSolverState()
{
    // Deallocate the saddle-point solver if not re-initializing
    if (!d_reinitializing_solver)
    {
        d_sp_solver->deallocateSolverState();
    }

    // Deallocate widened patch data.
    Pointer<PatchHierarchy<NDIM> > hierarchy = d_x_wide->getPatchHierarchy();
    const int coarsest_ln = d_x_wide->getCoarsestLevelNumber();
    const int finest_ln = d_x_wide->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_wide_u_idx)) level->deallocatePatchData(d_wide_u_idx);
        if (level->checkAllocated(d_wide_f_idx)) level->deallocatePatchData(d_wide_f_idx);
    }

    // Free the vectors.
    d_x_wide.setNull();
    d_b_wide.setNull();

    return;
} // deallocateSolverState

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
