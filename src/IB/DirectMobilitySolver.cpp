// Filename: DirectMobilitySolver.cpp
// Created on 20 Feb 2015 by Amneet Bhalla  and Bakytzhan Kallemov 
//
// Copyright (c) 2002-2015, Amneet Bhalla and Boyce Griffith
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

#include <math.h>
#include <algorithm>

#include "CartesianGridGeometry.h"
#include "PatchHierarchy.h"
#include "ibamr/CIBStrategy.h"
#include "ibamr/DirectMobilitySolver.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/ibamr_utilities.h"
#include "ibtk/PETScMultiVec.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/namespaces.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"

extern "C" {

// LAPACK function to do LU factorization
int dgetrf_(const int& n1, const int& n2, double* a, const int& lda, int* ipiv, int& info);

// LAPACK function to find soultion using the LU factorization
int dgetrs_(const char* trans,
            const int& n,
            const int& nrhs,
            const double* a,
            const int& lda,
            const int* ipiv,
            double* b,
            const int& ldb,
            int& info);

// LAPACK function to do Cholesky factorization
int dpotrf_(const char* uplo, const int& n, double* a, const int& lda, int& info);

// LAPACK function to find solution using Cholesky factorization
int dpotrs_(const char* uplo,
            const int& n,
            const int& nrhs,
            const double* a,
            const int& lda,
            double* b,
            const int& ldb,
            int& info);

// LAPACK function to do SVD factorization
void dsyevr_(const char* jobz,
             const char* range,
             const char* uplo,
             const int& n,
             double* a,
             const int& lda,
             const double& vl,
             const double& vu,
             const int& il,
             const int& iu,
             const double& abstol,
             int& m,
             double* w,
             double* z,
             const int& ldz,
             int* isuppz,
             double* work,
             const int& lwork,
             int* iwork,
             const int& liwork,
             int& info);
}

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

////////////////////////////// PUBLIC ////////////////////////////////////////

DirectMobilitySolver::DirectMobilitySolver(const std::string& object_name,
                                           Pointer<Database> input_db,
                                           Pointer<CIBStrategy> cib_strategy)
    : d_object_name(object_name), d_cib_strategy(cib_strategy)
{
    // Some default values
    d_is_initialized = false;
    d_recompute_mob_mat = false;
    d_f_periodic_corr = 0.0;

    // Get from input
    if (input_db) getFromInput(input_db);

    IBAMR_DO_ONCE(t_solve_system = TimerManager::getManager()->getTimer("IBAMR::DirectMobilitySolver::solveSystem()");
                  t_initialize_solver_state =
                      TimerManager::getManager()->getTimer("IBAMR::DirectMobilitySolver::initializeSolverState()");
                  t_deallocate_solver_state =
                      TimerManager::getManager()->getTimer("IBAMR::DirectMobilitySolver::deallocateSolverState()"););

    return;
} // DirectMobilitySolver

DirectMobilitySolver::~DirectMobilitySolver()
{
    for (std::map<std::string, double*>::iterator it = d_managed_mat_map.begin(); it != d_managed_mat_map.end(); ++it)
    {
        delete[] it->second;
        delete[] d_ipiv[it->first];
    }

    d_is_initialized = false;

    return;
} // ~DirectMobilitySolver

void DirectMobilitySolver::registerMobilityMat(const std::string& mat_name,
                                               const unsigned prototype_struct_id,
                                               MobilityMatrixType mat_type,
                                               MobilityMatrixInverseType inv_type,
                                               const int managing_proc,
                                               const std::string& filename,
                                               std::pair<double, double> scale)
{
    registerMobilityMat(mat_name, std::vector<unsigned int>(1, prototype_struct_id), mat_type, inv_type, managing_proc,
                        filename, scale);

    return;
} // registerMobilityMat

void DirectMobilitySolver::registerMobilityMat(const std::string& mat_name,
                                               const std::vector<unsigned>& prototype_struct_ids,
                                               MobilityMatrixType mat_type,
                                               MobilityMatrixInverseType inv_type,
                                               const int managing_proc,
                                               const std::string& filename,
                                               std::pair<double, double> scale)
{

#if !defined(NDEBUG)
    TBOX_ASSERT(!mat_name.empty());
    for (unsigned k = 0; k < prototype_struct_ids.size(); ++k)
    {
        TBOX_ASSERT(prototype_struct_ids[k] < d_cib_strategy->getNumberOfRigidStructures());
    }
    TBOX_ASSERT(d_managed_mat_map.find(mat_name) == d_managed_mat_map.end());
    TBOX_ASSERT(mat_type != UNKNOWN_MOBILITY_MATRIX_TYPE);
    TBOX_ASSERT(inv_type != UNKNOWN_MOBILITY_MATRIX_INVERSE_TYPE);
#endif

    unsigned int num_nodes = 0;
    for (unsigned k = 0; k < prototype_struct_ids.size(); ++k)
    {
        num_nodes += d_cib_strategy->getNumberOfNodes(prototype_struct_ids[k]);
    }

    // Fill-in various maps.
    d_managed_mat_prototype_id_map[mat_name] = prototype_struct_ids;
    d_managed_mat_proc_map[mat_name] = managing_proc;
    d_managed_mat_nodes_map[mat_name] = num_nodes;
    d_managed_mat_type_map[mat_name] = mat_type;
    d_managed_mat_inv_type_map[mat_name] = inv_type;
    d_managed_mat_filename_map[mat_name] = filename;
    d_managed_mat_scale_map[mat_name] = scale;
    d_managed_mat_map[mat_name] = NULL;
    d_ipiv[mat_name] = NULL;

    // Allocate the actual matrix
    const int size = num_nodes * NDIM;
    const int rank = SAMRAI_MPI::getRank();
    if (rank == managing_proc)
    {
        d_managed_mat_map[mat_name] = new double[size * size];
        if (inv_type == LAPACK_LU)
        {
            d_ipiv[mat_name] = new int[size];
        }
    }

    return;
} // registerMobilityMat

void DirectMobilitySolver::registerStructIDsWithMobilityMat(const std::string& mat_name,
                                                            const std::vector<std::vector<unsigned> >& struct_ids)
{

#if !defined(NDEBUG)
    TBOX_ASSERT(d_managed_mat_map.find(mat_name) != d_managed_mat_map.end());
    for (unsigned i = 0; i < struct_ids.size(); ++i)
    {
        TBOX_ASSERT(struct_ids[i].size() == d_managed_mat_prototype_id_map[mat_name].size());
        unsigned num_nodes = 0;
        for (unsigned j = 0; j < struct_ids[i].size(); ++j)
        {
            TBOX_ASSERT(struct_ids[i][j] < d_cib_strategy->getNumberOfRigidStructures());

            num_nodes += d_cib_strategy->getNumberOfNodes(struct_ids[i][j]);
        }
        TBOX_ASSERT(num_nodes == d_managed_mat_nodes_map[mat_name]);
    }
#endif

    d_managed_mat_actual_id_map[mat_name] = struct_ids;

    return;
} // registerStructIDsWithMobilityMat

void DirectMobilitySolver::getFrictionMat(const std::string& mat_name,
                                          double** fm,
                                          int* size,
                                          int* managing_proc,
                                          MobilityMatrixInverseType* inv_method)
{

#if !defined(NDEBUG)
    TBOX_ASSERT(d_is_initialized);
#endif

    *fm = d_managed_mat_map[mat_name];
    *size = d_managed_mat_nodes_map[mat_name] * NDIM;
    *managing_proc = d_managed_mat_proc_map[mat_name];
    *inv_method = d_managed_mat_inv_type_map[mat_name];

    return;
} // getFrictionMat

void DirectMobilitySolver::setStokesSpecifications(const StokesSpecifications& stokes_spec)
{
    d_rho = stokes_spec.getRho();
    d_mu = stokes_spec.getMu();

    return;
} // setStokesSpecifications

void DirectMobilitySolver::setSolutionTime(const double solution_time)
{
    d_solution_time = solution_time;

    return;
} // setSolutionTime

void DirectMobilitySolver::setTimeInterval(double current_time, double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;

    return;
} // setTimeInterval

    bool DirectMobilitySolver::solveSystem(Vec x, Vec b, const bool skip_nonfree_parts)
{
    unsigned managed_mats = (unsigned)d_managed_mat_map.size();
    if (!managed_mats) return true;

    IBAMR_TIMER_START(t_solve_system);

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x, b);

    //check recreation of mobility matrix 
    if (d_recompute_mob_mat)
    {
	static double check_time=-1.0e+15;
	if (d_current_time>check_time) createMobilityMatrix();
	check_time = d_current_time;
    }

    const int rank = SAMRAI_MPI::getRank();
    for (std::map<std::string, double*>::iterator it = d_managed_mat_map.begin(); it != d_managed_mat_map.end(); ++it)
    {
        const std::string& mat_name = it->first;
        std::vector<std::vector<unsigned> >& struct_ids = d_managed_mat_actual_id_map[mat_name];
        const int managing_proc = d_managed_mat_proc_map[mat_name];

        static const int data_depth = NDIM;
        const int size = d_managed_mat_nodes_map[mat_name] * data_depth;
        for (unsigned k = 0; k < struct_ids.size(); ++k)
        {
	    //skipping unspecified kinematics parts if solving only for free parts
	    if (skip_nonfree_parts) 
	    {
	    	bool skipflag=true;
	    	const unsigned local_num_structs = (unsigned)struct_ids[k].size();
	    	for (unsigned kk = 0; kk < local_num_structs; ++kk)
		{
		    int num_free_dofs = 0;
		    d_cib_strategy->getSolveRigidBodyVelocity(struct_ids[k][kk], num_free_dofs);
	    	    if (num_free_dofs)
	    	    {
	    		skipflag=false;
	    		break;
	    	    }
		}
	    	if (skipflag) continue;
	    }
	    
	    double* rhs = NULL;
            if (rank == managing_proc) rhs = new double[size];
            d_cib_strategy->copyVecToArray(b, rhs, struct_ids[k], data_depth, managing_proc);
	    
	    //Here rotate vector to initial frame
	    if (!d_recompute_mob_mat)
	    	d_cib_strategy->rotateArrayInitalBodyFrame(rhs, struct_ids[k], true, managing_proc);

	    //compute solution
	    if (rank == managing_proc) computeSolution(mat_name, rhs);

	    //Rotate solution vector back to structure frame
	    if (!d_recompute_mob_mat)
	    	d_cib_strategy->rotateArrayInitalBodyFrame(rhs, struct_ids[k], false, managing_proc);

            d_cib_strategy->copyArrayToVec(x, rhs, struct_ids[k], data_depth, managing_proc);
            delete[] rhs;
        }
    }
    PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x));

    IBAMR_TIMER_STOP(t_solve_system);

    return true;
} // solveSystem

void DirectMobilitySolver::initializeSolverState(Vec x, Vec /*b*/)
{
    if (d_is_initialized) return;
    unsigned managed_mats = (unsigned)d_managed_mat_map.size();
    if (!managed_mats) return;

    IBAMR_TIMER_START(t_initialize_solver_state);

    // Get grid-info
    Vec* vx;
    IBTK::VecMultiVecGetSubVecs(x, &vx);
    d_patch_hierarchy =	IBTK::PETScSAMRAIVectorReal::getSAMRAIVector(vx[0])->getPatchHierarchy();

    if (!d_recompute_mob_mat) createMobilityMatrix();

    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_solver_state);
    return;
} // initializeSolverState

void DirectMobilitySolver::deallocateSolverState()
{
    d_is_initialized = false;

    return;
} // deallocateSolverState

const std::vector<unsigned>& DirectMobilitySolver::getPrototypeStructIDs(const std::string& mat_name)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_managed_mat_prototype_id_map.find(mat_name) != d_managed_mat_prototype_id_map.end());
#endif
    return d_managed_mat_prototype_id_map[mat_name];

} // getPrototypeStructIDs

const std::vector<std::vector<unsigned> >& DirectMobilitySolver::getStructIDs(const std::string& mat_name)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_managed_mat_actual_id_map.find(mat_name) != d_managed_mat_actual_id_map.end());
#endif
    return d_managed_mat_actual_id_map[mat_name];

} // getStructIDs

///////////////////////////// PRIVATE ////////////////////////////////////////

void DirectMobilitySolver::getFromInput(Pointer<Database> input_db)
{
    Pointer<Database> comp_db;
    comp_db = input_db->isDatabase("LAPACK_SVD") ? input_db->getDatabase("LAPACK_SVD") : Pointer<Database>(NULL);
    if (comp_db)
    {
        d_svd_inv_tol = comp_db->getDouble("eigenvalue_replace_value");
        d_svd_inv_eps = comp_db->getDouble("min_eigenvalue_threshold");
    }

    // Other parameters
    d_f_periodic_corr = input_db->getDoubleWithDefault("f_periodic_correction", d_f_periodic_corr);
    d_recompute_mob_mat = input_db->getBool("recompute_mob_mat_perstep");

    return;
} // getFromInput

void DirectMobilitySolver::generateFrictionMatrix()
{
    int rank = SAMRAI_MPI::getRank();
    for (std::map<std::string, double*>::iterator it = d_managed_mat_map.begin(); it != d_managed_mat_map.end(); ++it)
    {
        const std::string& mat_name = it->first;
        if (rank != d_managed_mat_proc_map[mat_name]) continue;

        MobilityMatrixInverseType inv_type = d_managed_mat_inv_type_map[mat_name];
        double* managed_mat = d_managed_mat_map[mat_name];
        int mat_size = d_managed_mat_nodes_map[mat_name] * NDIM;

        if (inv_type == LAPACK_CHOLESKY)
        {
            int err = 0;
            dpotrf_((char*)"L", mat_size, managed_mat, mat_size, err);
            if (err)
            {
                TBOX_ERROR("DirectMobilityMatrix::generateFrictionMatrix() Matrix inversion "
                           << "failed for matrix handle " << mat_name << " with error code " << err
                           << " using LAPACK CHOLESKY." << std::endl);
            }
        }
        else if (inv_type == LAPACK_LU)
        {
            int err = 0;
            dgetrf_(mat_size, mat_size, managed_mat, mat_size, d_ipiv[mat_name], err);
            if (err)
            {
                TBOX_ERROR("DirectMobilityMatrix::generateFrictionMatrix() Matrix inversion "
                           << "failed for matrix handle " << mat_name << " with error code " << err
                           << " using LAPACK LU." << std::endl);
            }
        }
        else if (inv_type == LAPACK_SVD)
        {
            /* Locals */
            int il, iu, m, lwork, liwork, iwkopt, err = 0;
            double abstol, vl, vu, wkopt;

            /* Local arrays */
            std::vector<int> isuppz(2 * mat_size);
            std::vector<double> w(mat_size);
            std::vector<double> z(mat_size * mat_size);

            abstol = -1.0; // Negative abstol means using the default value
            lwork = -1;    // Query and allocate the optimal workspace
            liwork = -1;

            /* Initiate eigenvalue problem solve*/
            dsyevr_((char*)"V", (char*)"A", (char*)"L", mat_size, managed_mat, mat_size, vl, vu, il, iu, abstol, m,
                    &w[0], &z[0], mat_size, &isuppz[0], &wkopt, lwork, &iwkopt, liwork, err);
            if (err)
            {
                TBOX_ERROR("DirectMobilityMatrix::generateFrictionMatrix() Matrix inversion "
                           << "failed for matrix handle " << mat_name << " with error code " << err
                           << " using LAPACK SVD at first stage." << std::endl);
            }

            lwork = (int)wkopt;
            std::vector<double> work(lwork);
            liwork = iwkopt;
            std::vector<int> iwork(liwork);

            /* Finalize eigenvalue problem solve*/
            dsyevr_((char*)"V", (char*)"A", (char*)"L", mat_size, managed_mat, mat_size, vl, vu, il, iu, abstol, m,
                    &w[0], &z[0], mat_size, &isuppz[0], &work[0], lwork, &iwork[0], liwork, err);
            if (err)
            {
                TBOX_ERROR("DirectMobilityMatrix::generateFrictionMatrix() Matrix inversion "
                           << "failed for matrix handle " << mat_name << " with error code " << err
                           << " using LAPACK SVD at second stage." << std::endl);
            }

            // Make negative eigenvalues to be equal to min eigen value from
            // input option
            int counter = 0, counter_zero = 0;
            for (int i = 0; i < mat_size; ++i)
            {
                if (w[i] < d_svd_inv_eps)
                {
                    w[i] = d_svd_inv_tol;
                    counter++;
                }
            }
            for (int i = 0; i < mat_size; ++i)
            {
                if (MathUtilities<double>::equalEps(w[i], 0.0))
                {
                    counter_zero++;
                }

                for (int j = 0; j < mat_size; ++j)
                {
                    if (MathUtilities<double>::equalEps(w[j], 0.0))
                    {
                        managed_mat[j * mat_size + i] = 0.0;
                    }
                    else
                    {
                        if (i == j)
                            managed_mat[j * mat_size + i] = z[j * mat_size + i] / sqrt(w[j]);
                        else
                            managed_mat[j * mat_size + i] = z[j * mat_size + i] / sqrt(w[j]);
                    }
                }
            }

            plog << "DirectMobilityMatrix::generateFrictionMatrix(): " << counter
                 << " eigenvalues for dense mobility matrix with handle " << mat_name
                 << "have been changed. Number of zero eigenvalues placed are " << counter_zero << std::endl;
        }
        else
        {
            TBOX_ERROR("DirectMobilityMatrix::generateFrictionMatrix(): Unsupported dense "
                       << "matrix inversion method called" << std::endl);
        }
    }
    return;
} // generateFrictionMatrix

void DirectMobilitySolver::computeSolution(const std::string& mat_name, double* rhs)
{
    const int mat_size = d_managed_mat_nodes_map[mat_name] * NDIM;
    double* friction_mat = d_managed_mat_map[mat_name];
    MobilityMatrixInverseType inv_type = d_managed_mat_inv_type_map[mat_name];

    if (inv_type == LAPACK_CHOLESKY)
    {
        int err = 0;
        dpotrs_((char*)"L", mat_size, 1, friction_mat, mat_size, rhs, mat_size, err);
        if (err)
        {
            TBOX_ERROR("DirectMobilitySolver::computeSolution() Solution failed using "
                       << "LAPACK CHOLESKY with error code " << err << std::endl);
        }
    }
    else if (inv_type == LAPACK_LU)
    {
        int err = 0;
        dgetrs_((char*)"N", mat_size, 1, friction_mat, mat_size, d_ipiv[mat_name], rhs, mat_size, err);
        if (err)
        {
            TBOX_ERROR("DirectMobilitySolver::computeSolution() Solution failed using "
                       << "LAPACK LU with error code " << err << std::endl);
        }
    }
    else if (inv_type == LAPACK_SVD)
    {
        std::vector<double> temp(mat_size);
        for (int i = 0; i < mat_size; ++i)
        {
            temp[i] = 0.0;
            for (int j = 0; j < mat_size; ++j)
            {
                temp[i] += friction_mat[i * mat_size + j] * rhs[j];
            }
        }

        for (int i = 0; i < mat_size; ++i)
        {
            rhs[i] = 0.0;
            for (int j = 0; j < mat_size; ++j)
            {
                rhs[i] += friction_mat[j * mat_size + i] * temp[j];
            }
        }
    }

    return;
} // computeSolution
  
void DirectMobilitySolver::createMobilityMatrix()
{
    unsigned managed_mats = (unsigned)d_managed_mat_map.size();
    if (!managed_mats) return;

    int file_counter = 0;
    static std::vector<bool> read_files(managed_mats, false);
    bool initial_time = !d_recompute_mob_mat;

    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM> > struct_patch_level = d_patch_hierarchy->getPatchLevel(finest_ln);
    const IntVector<NDIM>& ratio = struct_patch_level->getRatio();
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_patch_hierarchy->getGridGeometry();
    const double* dx0 = grid_geom->getDx();
    const double* X_upper = grid_geom->getXUpper();
    const double* X_lower = grid_geom->getXLower();
    double domain_extents[NDIM], dx[NDIM];
    for (int d = 0; d < NDIM; ++d)
    {
	dx[d] = dx0[d] / ratio(d);
	domain_extents[d] = X_upper[d] - X_lower[d];
    }

    for (std::map<std::string, double*>::iterator it = d_managed_mat_map.begin(); it != d_managed_mat_map.end();
	 ++it)
    {
	const std::string& mat_name = it->first;
	double* mobility_mat = it->second;
	MobilityMatrixType mat_type = d_managed_mat_type_map[mat_name];
	const std::vector<unsigned>& struct_ids = d_managed_mat_prototype_id_map[mat_name];
	const std::pair<double, double>& scale = d_managed_mat_scale_map[mat_name];
	const int managing_proc = d_managed_mat_proc_map[mat_name];

	if (mat_type == FILE && !read_files[file_counter])
	{
	    // Get matrix from file

	    read_files[file_counter] = true;
	}
	else
	{
	    d_cib_strategy->constructMobilityMatrix(mat_name, mat_type, mobility_mat, struct_ids, dx,
						    domain_extents, initial_time, d_rho, d_mu, scale,
						    d_f_periodic_corr, managing_proc);
	}
	++file_counter;
    }
    generateFrictionMatrix();
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR
