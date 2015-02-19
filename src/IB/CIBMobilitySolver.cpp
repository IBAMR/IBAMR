// Filename: CIBMobilitySolver.cpp
// Created on 19 Feb 2015 by Amneet Bhalla and Bakytzhan Kallemov
//
// Copyright (c) 2002-2015, Amneet Bhalla, Bakytzhan Kallemov,
// and Boyce Griffith
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
#include <limits>

#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h"
#include "ibtk/PETScMultiVec.h"
#include "tbox/TimerManager.h"
#include "petsc-private/petscimpl.h"
#include "CIBMobilitySolver.h"
#include "DenseMobilitySolver.h"

namespace IBAMR
{

/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_solve_system;
static Timer* t_initialize_solver_state;
static Timer* t_deallocate_solver_state;
}

CIBMobilitySolver::CIBMobilitySolver(  
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator,
    Pointer<CIBStrategy> cib_strategy):
    d_input_db(input_db),
    d_ins_integrator(navier_stokes_integrator),
    d_cib_strategy(cib_strategy),
    d_recompute_mobility_matrix(false),
    d_is_initialized(false)
{
	// Get from input.
	if (d_input_db) getFromInput(d_input_db);
	
    return;
}// CIBMobilitySolver

CIBMobilitySolver::~CIBMobilitySolver()
{
  return;
}// ~CIBMobilitySolver

void
CIBMobilitySolver::setSolutionTime(
    const double solution_time)
{
	d_solution_time = solution_time;
	return;
}// setSolutionTime

void
CIBMobilitySolver::setTimeInterval(
    double current_time,
    double new_time)
{
    d_current_time          = current_time;
    d_new_time              = new_time;
	
	return;
}// setTimeInterval
	
void
CIBMobilitySolver::initializeSolverState(
	Vec x,
	Vec b)
{
	if (d_is_initialized) return;
	
	if (d_mobility_solver_type == DIRECT)
	{
		initializeDirectSolver(x,b);
	}
	else if (d_mobility_solver_type == KRYLOV)
	{
		initializeKrylovSolver(x,b);
	}
	else
	{
		TBOX_ERROR("CIBMobilitySolver::initializeSolverState() Unknown mobility solver type"
				   << std::endl);
	}

    d_is_initialized = true;

    return;
}// initializeSolverState

void
CIBMobilitySolver::deallocateSolverState()
{
    if (!d_is_initialized) return;

	if (d_mobility_solver_type == DIRECT)
	{
		for (unsigned k = 0; k < d_dense_solvers.size(); ++k)
		{
			d_dense_solvers[k]->deallocateSolverState();
		}
	}
	else if (d_mobility_solver_type == ITERATIVE)
	{
		d_krylov_solver->deallocateSolverState();
	}
	else
	{
		TBOX_ERROR("CIBMobilitySolver::initializeSolverState() Unknown mobility solver type"
				   << std::endl);
	}
	
    d_is_initialized = false;
    return;
}// deallocateSolverState

void
CIBMobilitySolver::solveSystem(
    Vec x, 
    Vec b)
{ 
    int ierr;
    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    
    if (deallocate_after_solve) initializeSolverState(x,b);
    
    // Get components of x and b.
    Vec *vx, *vb;
    ierr = IBTK::VecMultiVecGetSubVecs(x,&vx); IBTK_CHKERRQ(ierr);
    ierr = IBTK::VecMultiVecGetSubVecs(b,&vb); IBTK_CHKERRQ(ierr);
 
    // Fill rhs in seq Vec on processor zero.
    generateRHS(vb[0]);
    
    //use default as Single Dense matrix
    int  Num_SubMatrices=1;
    int N_cycles=1;

    if(d_MobilitySolverType)
    {
	Num_SubMatrices = d_cib_method->getNumberOfStructures();
	N_cycles = Num_SubMatrices/SAMRAI_MPI::getNodes();
	if (Num_SubMatrices % SAMRAI_MPI::getNodes()) N_cycles++;
    }

    for (int icycle=0;icycle<N_cycles;icycle++)
    {

	const int struct_ID=SAMRAI_MPI::getRank()+icycle*SAMRAI_MPI::getNodes();
	if (struct_ID < Num_SubMatrices)
	{
	    //solve the system for rhs
	    LocalDenseMatrix[icycle]->computeSolution(d_all_b);
	}
    }

    // Scatter solution back to PETSc
    returnSolution(vx[0]);
    
    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();
    
    return;
}// solveSystem

////////////////////////////// PRIVATE ///////////////////////////////////////

void
CIBMobilitySolver::getFromInput(
	Pointer<Database> input_db)
{
	// Get the mobility solver type and subtype from input database.
	const std::string solver_type    = input_db->getString("mobility_solver_type");
	const std::string solver_subtype = input_db->getString("mobility_solver_subtype");
	
	if (solver_type == "DIRECT")
	{
		d_mobility_solver_type = DIRECT;
		if      (solver_subtype == "DENSE")          d_mobility_solver_subtype = DENSE;
		else if (solver_subtype == "BLOCK_DIAGONAL") d_mobility_solver_subtype = BLOCK_DIAGONAL;
		else
		{
			TBOX_ERROR("CIBMobilitySolver::getFromInput() Unknown direct mobility solver subtype = "
					   << solver_subtype << " provided." << std::endl);
		}
	}
	else if (solver_type == "KRYLOV")
	{
		d_mobility_solver_type = KRYLOV;
		if (solver_subtype == "KRYLOV_REGULAR")   d_mobility_solver_subtype = KRYLOV_REGULAR;
		else if (solver_subtype == "KRYLOV_FMM")  d_mobility_solver_subtype = KRYLOV_FMM;
		else
		{
			TBOX_ERROR("CIBMobilitySolver::getFromInput() Unknown Krylov mobility solver subtype = "
					   << solver_subtype << " provided." << std::endl);
		}
	}
	else
	{
		TBOX_ERROR("CIBMobilitySolver::getFromInput() Unknown mobility solver type = "
				   << solver_type << " provided." << std::endl);
	}

	// Get if the mobility matrix is to recomputed every timestep.
	d_recompute_mobility_matrix =  input_db->getBool("compute_mob_matrix_everystep");

	return;
}
	
void
CIBMobilitySolver::initializeDirect()
{
    if (d_is_initialized && !d_recompute_mobility_matrix) return; 
  
    //store number of structures
    d_num_rigid_parts =  d_cib_strategy->getNumberOfRigidStructures();


    //use default as Single Dense matrix
    int  Num_SubMatrices=1;
    int N_cycles=1;

    if(d_MobilitySolverType)
    {
	Num_SubMatrices = d_num_rigid_parts;
	N_cycles = Num_SubMatrices/SAMRAI_MPI::getNodes();
	if (Num_SubMatrices % SAMRAI_MPI::getNodes()) N_cycles++;
    }
    
    clock_t start=0, end=0;
    if (SAMRAI_MPI::getRank() == 0)
    {
	start = clock();
	pout << "MobilitySolver::";
	if (d_MobilitySolverType == DENSE) pout<<"Setting up global dense mobility matrix...";
	else pout<<"Setting up block-diagonal mobility matrix...";
    }

    LocalDenseMatrix= new DenseMobilityMatrix *[N_cycles];

    PetscScalar *X_array, *W_array;
    PetscInt s, s2;

    //get positions and weights from CIBStrategy
    d_nodes_positions.resize(d_num_rigid_parts, PETSC_NULL);
    d_nodes_weights.resize(d_num_rigid_parts, PETSC_NULL);
	
    for (unsigned int part = 0; part < d_num_rigid_parts; ++part)
    {
	d_nodes_positions[part]  = d_cib_strategy->getStructurePETScVecPositions(ipart);
        d_nodes_weights[part]    = d_cib_strategy->getStructurePETScVecWeights(ipart);
    }

    Vec            X_petsc_SEQ, W_petsc_SEQ;
    VecScatter     ctx, ctw;
    for (int istruct=0;istruct<d_num_rigid_parts;icycle++)
    {
	VecScatterCreateToAll(d_nodes_positions[istruct],&ctx,&X_petsc_SEQ);
	VecScatterBegin(ctx,d_nodes_positions[istruct],X_petsc_SEQ,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(ctx,d_nodes_positions[istruct],X_petsc_SEQ,INSERT_VALUES,SCATTER_FORWARD);
	VecGetArray(X_petsc_SEQ,&X_array);
	VecGetSize (X_petsc_SEQ,&s);
	
	VecScatterCreateToAll(d_nodes_weights[istruct],&ctw,&W_petsc_SEQ);
	VecScatterBegin(ctw,d_nodes_weights[istruct],W_petsc_SEQ,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(ctw,d_nodes_weights[istruct],W_petsc_SEQ,INSERT_VALUES,SCATTER_FORWARD);
	VecGetArray(W_petsc_SEQ,&W_array);
	VecGetSize (W_petsc_SEQ,&s2);

	const int icycle= istruct/SAMRAI_MPI::getNodes();
	const int struct_ID=SAMRAI_MPI::getRank()+icycle*SAMRAI_MPI::getNodes();
	if (struct_ID ==istruct)
	{
	    LocalDenseMatrix[icycle]= new DenseMobilityMatrix( d_input_db->getDatabase("DenseMobilityMatrix"), d_ins_integrator, struct_ID);
	    //generate Mobility Submatrices
	    LocalDenseMatrix[icycle]->initializeDenseMobilityMatrix(d_dt, X_Array, W_Array);
	}
	    
	VecRestoreArray(X_petsc_SEQ,&X_array);
	VecDestroy(&X_petsc_SEQ);
	VecScatterDestroy(&ctx);
	
	VecRestoreArray(W_petsc_SEQ,&W_array);
	VecDestroy(&W_petsc_SEQ);
	VecScatterDestroy(&ctw);
    }

    for (int icycle=0;icycle<N_cycles;icycle++)
    {	
	const int struct_ID=SAMRAI_MPI::getRank()+icycle*SAMRAI_MPI::getNodes();
	if (struct_ID < Num_SubMatrices)
	{
	    //generate Mobility Submatrices
	    LocalDenseMatrix[icycle]->generateMobilityMatrix();
	}
    }
    if (SAMRAI_MPI::getRank() == 0) 
    {
	pout << "Done."<< std::endl;
	pout << "Factorizing...";
    }
   
    for (int icycle=0;icycle<N_cycles;icycle++)
    {
	const int struct_ID=SAMRAI_MPI::getRank()+icycle*SAMRAI_MPI::getNodes();
	if (struct_ID < Num_SubMatrices)
	{
	    //invert mobility submatrices
	    LocalDenseMatrix[icycle]->generateFrictionMatrix();
	}
    }
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0)
    {
	end = clock();
	pout << "Done. Time taken: " << double(end-start)/double(CLOCKS_PER_SEC)<< std::endl;
	pout<<std::endl;
    }

    d_is_initialized = true;
    return;
}//initialize BlockDiagonalMatrix 

void
CIBMobilitySolver::generateRHS(
    Vec b)
{
    int ierr;

    PetscScalar* a;
    PetscInt s;
    Vec            V_SEQ;
    VecScatter     ctx;

    VecScatterCreateToAll(b,&ctx,&V_SEQ);
    VecScatterBegin(ctx,b,V_SEQ,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(ctx,b,V_SEQ,INSERT_VALUES,SCATTER_FORWARD);
    VecGetArray(V_SEQ,&a);
    VecGetSize (V_SEQ,&s);

    if (istruct == d_struct_ID )
    {
	for (i = 0; i < ; i++) X[i] = a[i];
    }
    VecRestoreArray(V_SEQ,&_a);
    VecScatterDestroy(&ctx);
    VecDestroy(&V_SEQ);


    int size = NDIM*d_no_all_blobs;
    if (d_MobilitySolverType || !SAMRAI_MPI::getRank())
    {

	if (!d_all_b) d_all_b = new double[size];
    
	PetscScalar *barray;
	ierr = VecGetArray(b_lag_vec_seq,&barray);     IBTK_CHKERRQ(ierr);
    
	for (int k = 0; k < size; ++k) 	d_all_b[k]  = barray[k]; 
	ierr = VecRestoreArray(b_lag_vec_seq,&barray); IBTK_CHKERRQ(ierr);
    }
    VecDestroy(&b_lag_vec_parallel); 
    VecDestroy(&b_lag_vec_seq);

    return;
}// generateRHS

  

void
CIBMobilitySolver::returnSolution(
    Vec x)
{

    //use default as Single Dense matrix
    int  Num_SubMatrices=1;
    int N_cycles=1;

    if(d_MobilitySolverType)
    {
	Num_SubMatrices = d_cib_method->getNumberOfStructures();
	N_cycles = Num_SubMatrices/SAMRAI_MPI::getNodes();
	if (Num_SubMatrices % SAMRAI_MPI::getNodes()) N_cycles++;
    }

    for (int icycle=0;icycle<N_cycles;icycle++)
    {
	const int struct_ID=SAMRAI_MPI::getRank()+icycle*SAMRAI_MPI::getNodes();
	if (struct_ID < Num_SubMatrices)
	{
	    //setting input vectors for positions, rhs and weights  
	    const std::pair<int,int>& range=(d_cib_method->getIBKinematics())[struct_ID]->getStructureParameters().getLagIdxRange();
	    const  int offset=range.first*NDIM;
	    
	    int no_blobs=d_no_all_blobs;
	    if(d_MobilitySolverType) no_blobs = d_cib_method->getNumberOfStructuresNodes(struct_ID);

	    const int size = no_blobs*NDIM;
	    PetscInt indices[size]; 
	    for (int ipart = 0; ipart < size; ++ipart) indices[ipart] = offset+ipart;
	    VecSetValues(x_lag_vec, size, indices, d_all_b+offset, INSERT_VALUES); 
	}
    } 
    
    VecAssemblyBegin(x_lag_vec);
    VecAssemblyEnd(x_lag_vec);
    
    l_data_manager->scatterLagrangianToPETSc(x_lag_vec, x, struct_level);
    
    VecDestroy(&x_lag_vec); 
    return;
}  
}

