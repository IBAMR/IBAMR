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

#include <algorithm>
#include <math.h>

#include "CartesianGridGeometry.h"
#include "PatchHierarchy.h"
#include "ibamr/CIBStrategy.h"
#include "ibamr/DirectMobilitySolver.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "petsc/private/petscimpl.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"

extern "C" {

// LAPACK function to do LU factorization.
int dgetrf_(const int& n1, const int& n2, double* a, const int& lda, int* ipiv, int& info);

// LAPACK function to find soultion using the LU factorization.
int dgetrs_(const char* trans,
            const int& n,
            const int& nrhs,
            const double* a,
            const int& lda,
            const int* ipiv,
            double* b,
            const int& ldb,
            int& info);

// LAPACK function to do Cholesky factorization.
int dpotrf_(const char* uplo, const int& n, double* a, const int& lda, int& info);

// LAPACK function to find solution using Cholesky factorization.
int dpotrs_(const char* uplo,
            const int& n,
            const int& nrhs,
            const double* a,
            const int& lda,
            double* b,
            const int& ldb,
            int& info);

// LAPACK function to do SVD factorization.
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
static Timer* t_solve_body_system;
static Timer* t_initialize_solver_state;
static Timer* t_deallocate_solver_state;
}

////////////////////////////// PUBLIC ////////////////////////////////////////

DirectMobilitySolver::DirectMobilitySolver(const std::string& object_name,
                                           Pointer<Database> input_db,
                                           Pointer<CIBStrategy> cib_strategy)

{
    d_object_name = object_name;
    d_cib_strategy = cib_strategy;

    // Some default values
    d_is_initialized = false;
    d_recompute_mob_mat = false;
    d_f_periodic_corr = 0.0;

    // Get from input
    if (input_db) getFromInput(input_db);

    IBAMR_DO_ONCE(t_solve_system = TimerManager::getManager()->getTimer("IBAMR::DirectMobilitySolver::solveSystem()");
                  t_solve_body_system =
                      TimerManager::getManager()->getTimer("IBAMR::DirectMobilitySolver::solveBodySystem()");
                  t_initialize_solver_state =
                      TimerManager::getManager()->getTimer("IBAMR::DirectMobilitySolver::initializeSolverState()");
                  t_deallocate_solver_state =
                      TimerManager::getManager()->getTimer("IBAMR::DirectMobilitySolver::deallocateSolverState()"););

    return;
} // DirectMobilitySolver

DirectMobilitySolver::~DirectMobilitySolver()
{
    for (std::map<std::string, std::pair<Mat, Mat> >::iterator it = d_petsc_mat_map.begin();
         it != d_petsc_mat_map.end();
         ++it)
    {
        const std::string& mat_name = it->first;
        Mat& mobility_mat = d_petsc_mat_map[mat_name].first;
        Mat& body_mobility_mat = d_petsc_mat_map[mat_name].second;
        MatDestroy(&mobility_mat);
        MatDestroy(&body_mobility_mat);
    }

    for (std::map<std::string, std::pair<double*, double*> >::iterator it = d_mat_map.begin(); it != d_mat_map.end();
         ++it)
    {
        delete[](it->second).first;
        delete[](it->second).second;
    }

    for (std::map<std::string, Mat>::iterator it = d_petsc_geometric_mat_map.begin();
         it != d_petsc_geometric_mat_map.end();
         ++it)
    {
        Mat& geometric_mat = d_petsc_geometric_mat_map[it->first];
        MatDestroy(&geometric_mat);
    }

    for (std::map<std::string, double*>::iterator it = d_geometric_mat_map.begin(); it != d_geometric_mat_map.end();
         ++it)
    {
        delete[] it->second;
    }

    for (std::map<std::string, std::pair<int*, int*> >::iterator it = d_ipiv_map.begin(); it != d_ipiv_map.end(); ++it)
    {
        delete[](it->second).first;
        delete[](it->second).second;
    }

    d_is_initialized = false;

    return;
} // ~DirectMobilitySolver

void
DirectMobilitySolver::registerMobilityMat(const std::string& mat_name,
                                          const unsigned prototype_struct_id,
                                          MobilityMatrixType mat_type,
                                          std::pair<MobilityMatrixInverseType, MobilityMatrixInverseType> inv_type,
                                          const int managing_proc,
                                          const std::string& filename,
                                          std::pair<double, double> scale)
{
    registerMobilityMat(mat_name,
                        std::vector<unsigned int>(1, prototype_struct_id),
                        mat_type,
                        inv_type,
                        managing_proc,
                        filename,
                        scale);

    return;
} // registerMobilityMat

void
DirectMobilitySolver::registerMobilityMat(const std::string& mat_name,
                                          const std::vector<unsigned>& prototype_struct_ids,
                                          MobilityMatrixType mat_type,
                                          std::pair<MobilityMatrixInverseType, MobilityMatrixInverseType> inv_type,
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
    TBOX_ASSERT(d_mat_map.find(mat_name) == d_mat_map.end());
    TBOX_ASSERT(mat_type != UNKNOWN_MOBILITY_MATRIX_TYPE);
    TBOX_ASSERT(inv_type.first != UNKNOWN_MOBILITY_MATRIX_INVERSE_TYPE);
    TBOX_ASSERT(inv_type.second != UNKNOWN_MOBILITY_MATRIX_INVERSE_TYPE);
#endif

    unsigned int num_nodes = 0;
    for (unsigned k = 0; k < prototype_struct_ids.size(); ++k)
    {
        num_nodes += d_cib_strategy->getNumberOfNodes(prototype_struct_ids[k]);
    }

    // Fill-in various maps.
    d_mat_prototype_id_map[mat_name] = prototype_struct_ids;
    d_mat_proc_map[mat_name] = managing_proc;
    d_mat_nodes_map[mat_name] = num_nodes;
    d_mat_parts_map[mat_name] = static_cast<unsigned>(prototype_struct_ids.size());
    d_mat_type_map[mat_name] = mat_type;
    d_mat_inv_type_map[mat_name] = inv_type;
    d_mat_filename_map[mat_name] = filename;
    d_mat_scale_map[mat_name] = scale;
    d_mat_map[mat_name] = std::make_pair<double*, double*>(NULL, NULL);
    d_geometric_mat_map[mat_name] = NULL;
    d_ipiv_map[mat_name] = std::make_pair<int*, int*>(NULL, NULL);
    d_petsc_mat_map[mat_name] = std::make_pair<Mat, Mat>(NULL, NULL);
    d_petsc_geometric_mat_map[mat_name] = NULL;

    // Allocate the actual matrices.
    const int mobility_mat_size = num_nodes * NDIM;
    const int body_mobility_mat_size = d_mat_parts_map[mat_name] * s_max_free_dofs;
    const int rank = SAMRAI_MPI::getRank();

    if (rank == managing_proc)
    {
        d_mat_map[mat_name].first = new double[mobility_mat_size * mobility_mat_size];
        MatCreateSeqDense(PETSC_COMM_SELF,
                          mobility_mat_size,
                          mobility_mat_size,
                          d_mat_map[mat_name].first,
                          &d_petsc_mat_map[mat_name].first);

        d_mat_map[mat_name].second = new double[body_mobility_mat_size * body_mobility_mat_size];
        MatCreateSeqDense(PETSC_COMM_SELF,
                          body_mobility_mat_size,
                          body_mobility_mat_size,
                          d_mat_map[mat_name].second,
                          &d_petsc_mat_map[mat_name].second);

        d_geometric_mat_map[mat_name] = new double[mobility_mat_size * body_mobility_mat_size];
        MatCreateSeqDense(PETSC_COMM_SELF,
                          mobility_mat_size,
                          body_mobility_mat_size,
                          d_geometric_mat_map[mat_name],
                          &d_petsc_geometric_mat_map[mat_name]);

        if (d_mat_inv_type_map[mat_name].first == LAPACK_LU)
        {
            d_ipiv_map[mat_name].first = new int[mobility_mat_size];
        }
        if (d_mat_inv_type_map[mat_name].second == LAPACK_LU)
        {
            d_ipiv_map[mat_name].second = new int[body_mobility_mat_size];
        }
    }

    return;
} // registerMobilityMat

void
DirectMobilitySolver::registerStructIDsWithMobilityMat(const std::string& mat_name,
                                                       const std::vector<std::vector<unsigned> >& struct_ids)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_mat_map.find(mat_name) != d_mat_map.end());
    for (unsigned i = 0; i < struct_ids.size(); ++i)
    {
        TBOX_ASSERT(struct_ids[i].size() == d_mat_prototype_id_map[mat_name].size());
        unsigned num_nodes = 0;
        for (unsigned j = 0; j < struct_ids[i].size(); ++j)
        {
            TBOX_ASSERT(struct_ids[i][j] < d_cib_strategy->getNumberOfRigidStructures());

            num_nodes += d_cib_strategy->getNumberOfNodes(struct_ids[i][j]);
        }
        TBOX_ASSERT(num_nodes == d_mat_nodes_map[mat_name]);
    }
#endif

    d_mat_actual_id_map[mat_name] = struct_ids;

    return;
} // registerStructIDsWithMobilityMat

void
DirectMobilitySolver::setStokesSpecifications(const StokesSpecifications& stokes_spec)
{
    d_rho = stokes_spec.getRho();
    d_mu = stokes_spec.getMu();

    return;
} // setStokesSpecifications

void
DirectMobilitySolver::setSolutionTime(const double solution_time)
{
    d_solution_time = solution_time;

    return;
} // setSolutionTime

void
DirectMobilitySolver::setTimeInterval(double current_time, double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;

    return;
} // setTimeInterval

bool
DirectMobilitySolver::solveSystem(Vec x, Vec b)
{
    IBAMR_TIMER_START(t_solve_system);

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x, b);

    const int rank = SAMRAI_MPI::getRank();
    static const int data_depth = NDIM;

    for (std::map<std::string, std::pair<Mat, Mat> >::iterator it = d_petsc_mat_map.begin();
         it != d_petsc_mat_map.end();
         ++it)
    {
        const std::string& mat_name = it->first;
        Mat& mat = d_petsc_mat_map[mat_name].first;
        const MobilityMatrixInverseType& inv_type = d_mat_inv_type_map[mat_name].first;
        const std::vector<std::vector<unsigned> >& struct_ids = d_mat_actual_id_map[mat_name];
        const int managing_proc = d_mat_proc_map[mat_name];
        const int mat_size = d_mat_nodes_map[mat_name] * data_depth;
        const int num_structs = static_cast<int>(struct_ids.size());

        for (int k = 0; k < num_structs; ++k)
        {
            double* rhs = NULL;
            if (rank == managing_proc) rhs = new double[mat_size];
            d_cib_strategy->copyVecToArray(b, rhs, struct_ids[k], data_depth, managing_proc);
            if (!d_recompute_mob_mat)
            {
                d_cib_strategy->rotateArray(rhs,
                                            struct_ids[k],
                                            /*use_transpose*/ true,
                                            managing_proc,
                                            data_depth);
            }
            if (rank == managing_proc) computeSolution(mat, inv_type, d_ipiv_map[mat_name].first, rhs);
            if (!d_recompute_mob_mat)
            {
                d_cib_strategy->rotateArray(rhs,
                                            struct_ids[k],
                                            /*use_transpose*/ false,
                                            managing_proc,
                                            data_depth);
            }
            d_cib_strategy->copyArrayToVec(x, rhs, struct_ids[k], data_depth, managing_proc);
            delete[] rhs;
        }
    }

    IBAMR_TIMER_STOP(t_solve_system);

    return true;
} // solveSystem

bool
DirectMobilitySolver::solveBodySystem(Vec x, Vec b)
{
    IBAMR_TIMER_START(t_solve_body_system);

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x, b);

    const int rank = SAMRAI_MPI::getRank();
    static const int data_depth = s_max_free_dofs;

    for (std::map<std::string, std::pair<Mat, Mat> >::iterator it = d_petsc_mat_map.begin();
         it != d_petsc_mat_map.end();
         ++it)
    {
        const std::string& mat_name = it->first;
        Mat& mat = d_petsc_mat_map[mat_name].second;
        const MobilityMatrixInverseType& inv_type = d_mat_inv_type_map[mat_name].second;
        const std::vector<std::vector<unsigned> >& struct_ids = d_mat_actual_id_map[mat_name];
        const int mat_size = d_mat_parts_map[mat_name] * data_depth;
        const int managing_proc = d_mat_proc_map[mat_name];
        const int num_structs = static_cast<int>(struct_ids.size());

        for (int k = 0; k < num_structs; ++k)
        {
            double* rhs = NULL;
            if (rank == managing_proc) rhs = new double[mat_size];
            d_cib_strategy->copyFreeDOFsVecToArray(b, rhs, struct_ids[k], managing_proc);
            if (!d_recompute_mob_mat)
            {
                d_cib_strategy->rotateArray(rhs,
                                            struct_ids[k],
                                            /*use_transpose*/ true,
                                            managing_proc,
                                            data_depth);
            }
            if (rank == managing_proc) computeSolution(mat, inv_type, d_ipiv_map[mat_name].second, rhs);
            if (!d_recompute_mob_mat)
            {
                d_cib_strategy->rotateArray(rhs,
                                            struct_ids[k],
                                            /*use_transpose*/ false,
                                            managing_proc,
                                            data_depth);
            }
            d_cib_strategy->copyFreeDOFsArrayToVec(x, rhs, struct_ids[k], managing_proc);
            delete[] rhs;
        }
    }

    IBAMR_TIMER_STOP(t_solve_body_system);

    return true;
} // solveBodySystem

void
DirectMobilitySolver::initializeSolverState(Vec x, Vec /*b*/)
{
    if (d_is_initialized) return;

    IBAMR_TIMER_START(t_initialize_solver_state);

    int rank = SAMRAI_MPI::getRank();
    unsigned managed_mats = static_cast<unsigned>(d_mat_map.size());

    static bool recreate_mobility_matrices = true;
    static std::vector<bool> read_files(managed_mats, false);
    bool initial_time = !d_recompute_mob_mat;

    if (recreate_mobility_matrices)
    {
        // Get grid-info
        Vec* vx;
        VecNestGetSubVecs(x, NULL, &vx);
        Pointer<SAMRAIVectorReal<NDIM, double> > vx0;
        IBTK::PETScSAMRAIVectorReal::getSAMRAIVectorRead(vx[0], &vx0);
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = vx0->getPatchHierarchy();
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        IBTK::PETScSAMRAIVectorReal::restoreSAMRAIVectorRead(vx[0], &vx0);
        Pointer<PatchLevel<NDIM> > struct_patch_level = patch_hierarchy->getPatchLevel(finest_ln);
        const IntVector<NDIM>& ratio = struct_patch_level->getRatio();
        Pointer<CartesianGridGeometry<NDIM> > grid_geom = patch_hierarchy->getGridGeometry();
        const double* dx0 = grid_geom->getDx();
        const double* X_upper = grid_geom->getXUpper();
        const double* X_lower = grid_geom->getXLower();
        double domain_extents[NDIM], dx[NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            dx[d] = dx0[d] / ratio(d);
            domain_extents[d] = X_upper[d] - X_lower[d];
        }

        int file_counter = 0;
        for (std::map<std::string, std::pair<Mat, Mat> >::iterator it = d_petsc_mat_map.begin();
             it != d_petsc_mat_map.end();
             ++it, ++file_counter)
        {
            const std::string& mat_name = it->first;
            Mat& mobility_mat = d_petsc_mat_map[mat_name].first;
            Mat& geometric_mat = d_petsc_geometric_mat_map[mat_name];
            const MobilityMatrixType& mat_type = d_mat_type_map[mat_name];
            const std::vector<unsigned>& struct_ids = d_mat_prototype_id_map[mat_name];
            const std::pair<double, double>& scale = d_mat_scale_map[mat_name];
            const int managing_proc = d_mat_proc_map[mat_name];

            if (mat_type == READ_FROM_FILE && !read_files[file_counter])
            {
                // Get the matrix from file.
                const std::string& filename = d_mat_filename_map[mat_name];
                if (rank == managing_proc)
                {
                    PetscViewer binary_viewer;
                    PetscViewerBinaryOpen(PETSC_COMM_SELF, filename.c_str(), FILE_MODE_READ, &binary_viewer);
                    MatLoad(mobility_mat, binary_viewer);
                    PetscViewerDestroy(&binary_viewer);
                }

                read_files[file_counter] = true;
            }
            else
            {
                d_cib_strategy->constructMobilityMatrix(mat_name,
                                                        mat_type,
                                                        mobility_mat,
                                                        struct_ids,
                                                        dx,
                                                        domain_extents,
                                                        initial_time,
                                                        d_rho,
                                                        d_mu,
                                                        scale,
                                                        d_f_periodic_corr,
                                                        managing_proc);
            }

            // Construct the geometric matrix that maps rigid body velocity to
            // nodal velocity.
            d_cib_strategy->constructGeometricMatrix(mat_name, geometric_mat, struct_ids, initial_time, managing_proc);
        }
        factorizeMobilityMatrix();
        constructBodyMobilityMatrix();
        factorizeBodyMobilityMatrix();
    }

    d_is_initialized = true;
    recreate_mobility_matrices = d_recompute_mob_mat;

    IBAMR_TIMER_STOP(t_initialize_solver_state);

    return;
} // initializeSolverState

void
DirectMobilitySolver::deallocateSolverState()
{
    d_is_initialized = false;

    return;
} // deallocateSolverState

const std::vector<unsigned>&
DirectMobilitySolver::getPrototypeStructIDs(const std::string& mat_name)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_mat_prototype_id_map.find(mat_name) != d_mat_prototype_id_map.end());
#endif
    return d_mat_prototype_id_map[mat_name];

} // getPrototypeStructIDs

const std::vector<std::vector<unsigned> >&
DirectMobilitySolver::getStructIDs(const std::string& mat_name)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_mat_actual_id_map.find(mat_name) != d_mat_actual_id_map.end());
#endif
    return d_mat_actual_id_map[mat_name];

} // getStructIDs

///////////////////////////// PRIVATE ////////////////////////////////////////

void
DirectMobilitySolver::getFromInput(Pointer<Database> input_db)
{
    Pointer<Database> comp_db;
    comp_db = input_db->isDatabase("LAPACK_SVD") ? input_db->getDatabase("LAPACK_SVD") : Pointer<Database>(NULL);
    if (comp_db)
    {
        d_svd_replace_value = comp_db->getDouble("eigenvalue_replace_value");
        d_svd_eps = comp_db->getDouble("min_eigenvalue_threshold");
    }

    // Other parameters
    d_f_periodic_corr = input_db->getDoubleWithDefault("f_periodic_correction", d_f_periodic_corr);
    d_recompute_mob_mat = input_db->getBoolWithDefault("recompute_mob_mat_perstep", d_recompute_mob_mat);

    return;
} // getFromInput

void
DirectMobilitySolver::factorizeMobilityMatrix()
{
    int rank = SAMRAI_MPI::getRank();
    for (std::map<std::string, std::pair<Mat, Mat> >::iterator it = d_petsc_mat_map.begin();
         it != d_petsc_mat_map.end();
         ++it)
    {
        const std::string& mat_name = it->first;
        if (rank != d_mat_proc_map[mat_name]) continue;

        Mat& mat = d_petsc_mat_map[mat_name].first;
        const MobilityMatrixInverseType& inv_type = d_mat_inv_type_map[mat_name].first;
        const int mat_size = d_mat_nodes_map[mat_name] * NDIM;
        double* mat_data = NULL;
        MatDenseGetArray(mat, &mat_data);
        factorizeDenseMatrix(mat_data, mat_size, inv_type, d_ipiv_map[mat_name].first, mat_name, "Mobility");
        MatDenseRestoreArray(mat, &mat_data);
    }
    return;

} // factorizeMobilityMatrix

void
DirectMobilitySolver::constructBodyMobilityMatrix()
{
    int rank = SAMRAI_MPI::getRank();
    for (std::map<std::string, std::pair<Mat, Mat> >::iterator it = d_petsc_mat_map.begin();
         it != d_petsc_mat_map.end();
         ++it)
    {
        const std::string& mat_name = it->first;
        if (rank != d_mat_proc_map[mat_name]) continue;

        const int row_size = d_mat_nodes_map[mat_name] * NDIM;
        const int col_size = d_mat_parts_map[mat_name] * s_max_free_dofs;
        const MobilityMatrixInverseType& mobility_inv_type = d_mat_inv_type_map[mat_name].first;

        Mat& mobility_mat = d_petsc_mat_map[mat_name].first;
        Mat& body_mob_mat = d_petsc_mat_map[mat_name].second;
        Mat& geometric_mat = d_petsc_geometric_mat_map[mat_name];

        // Allocate a temporary matrix that holds the Matrix-Matrix product.
        // Here we are multiplying inverse of mobility matrix with geometric matrix.
        double* product_mat_data = new double[row_size * col_size];
        Mat product_mat;
        MatCreateSeqDense(PETSC_COMM_SELF, row_size, col_size, product_mat_data, &product_mat);
        MatCopy(geometric_mat, product_mat, SAME_NONZERO_PATTERN);

        for (int col = 0; col < col_size; ++col)
        {
            double* col_data;
            MatDenseGetArray(product_mat, &col_data);
            computeSolution(mobility_mat, mobility_inv_type, d_ipiv_map[mat_name].first, &col_data[col * row_size]);
            MatDenseRestoreArray(product_mat, &col_data);
        }
        MatTransposeMatMult(geometric_mat, product_mat, MAT_REUSE_MATRIX, PETSC_DEFAULT, &body_mob_mat);

        MatDestroy(&product_mat);
        delete[] product_mat_data;
    }

    return;
} // generateBodyFrictionMatrix

void
DirectMobilitySolver::factorizeBodyMobilityMatrix()
{
    int rank = SAMRAI_MPI::getRank();
    for (std::map<std::string, std::pair<Mat, Mat> >::iterator it = d_petsc_mat_map.begin();
         it != d_petsc_mat_map.end();
         ++it)
    {
        const std::string& mat_name = it->first;
        if (rank != d_mat_proc_map[mat_name]) continue;

        Mat& mat = d_petsc_mat_map[mat_name].second;
        const MobilityMatrixInverseType& inv_type = d_mat_inv_type_map[mat_name].second;
        const int mat_size = d_mat_parts_map[mat_name] * s_max_free_dofs;

        double* mat_data = NULL;
        MatDenseGetArray(mat, &mat_data);
        factorizeDenseMatrix(mat_data, mat_size, inv_type, d_ipiv_map[mat_name].second, mat_name, "Body Mobility");
        MatDenseRestoreArray(mat, &mat_data);
    }
    return;

} // factorizeBodyMobilityMatrix

void
DirectMobilitySolver::factorizeDenseMatrix(double* mat_data,
                                           const int mat_size,
                                           const MobilityMatrixInverseType& inv_type,
                                           int* ipiv,
                                           const std::string& mat_name,
                                           const std::string& err_msg)
{
    int err = 0;
    if (inv_type == LAPACK_CHOLESKY)
    {
        dpotrf_((char*)"L", mat_size, mat_data, mat_size, err);
        if (err)
        {
            TBOX_ERROR("DirectMobilityMatrix::factorizeDenseMatrix(). " << err_msg << " matrix factorization "
                                                                        << " failed for matrix handle "
                                                                        << mat_name
                                                                        << " with error code "
                                                                        << err
                                                                        << " using LAPACK CHOLESKY."
                                                                        << std::endl);
        }
    }
    else if (inv_type == LAPACK_LU)
    {
        dgetrf_(mat_size, mat_size, mat_data, mat_size, ipiv, err);
        if (err)
        {
            TBOX_ERROR("DirectMobilityMatrix::factorizeDenseMatrix(). " << err_msg << " matrix factorization "
                                                                        << "failed for matrix handle "
                                                                        << mat_name
                                                                        << " with error code "
                                                                        << err
                                                                        << " using LAPACK LU."
                                                                        << std::endl);
        }
    }
    else if (inv_type == LAPACK_SVD)
    {
        // Locals.
        int il, iu, m, lwork, liwork, iwkopt;
        double abstol, vl, vu, wkopt;

        // Local arrays.
        std::vector<int> isuppz(2 * mat_size);
        std::vector<double> w(mat_size);
        std::vector<double> z(mat_size * mat_size);

        abstol = -1.0; // Negative abstol means using the default value.
        lwork = -1;    // Query and allocate the optimal workspace.
        liwork = -1;

        // Initiate eigenvalue problem solve.
        dsyevr_((char*)"V",
                (char*)"A",
                (char*)"L",
                mat_size,
                mat_data,
                mat_size,
                vl,
                vu,
                il,
                iu,
                abstol,
                m,
                &w[0],
                &z[0],
                mat_size,
                &isuppz[0],
                &wkopt,
                lwork,
                &iwkopt,
                liwork,
                err);
        if (err)
        {
            TBOX_ERROR("DirectMobilityMatrix::factorizeDenseMatrix(). " << err_msg << " matrix factorization "
                                                                        << "failed for matrix handle "
                                                                        << mat_name
                                                                        << " with error code "
                                                                        << err
                                                                        << " using LAPACK SVD at first stage."
                                                                        << std::endl);
        }

        lwork = static_cast<int>(wkopt);
        std::vector<double> work(lwork);
        liwork = iwkopt;
        std::vector<int> iwork(liwork);

        // Finalize eigenvalue problem solve.
        dsyevr_((char*)"V",
                (char*)"A",
                (char*)"L",
                mat_size,
                mat_data,
                mat_size,
                vl,
                vu,
                il,
                iu,
                abstol,
                m,
                &w[0],
                &z[0],
                mat_size,
                &isuppz[0],
                &work[0],
                lwork,
                &iwork[0],
                liwork,
                err);
        if (err)
        {
            TBOX_ERROR("DirectMobilityMatrix::factorizeDenseMatrix(). " << err_msg << " matrix factorization "
                                                                        << "failed for matrix handle "
                                                                        << mat_name
                                                                        << " with error code "
                                                                        << err
                                                                        << " using LAPACK SVD at second stage."
                                                                        << std::endl);
        }

        // Make negative eigenvalues to be equal to min eigen value from
        // input option
        int counter = 0, counter_zero = 0;
        for (int i = 0; i < mat_size; ++i)
        {
            if (w[i] < d_svd_eps)
            {
                w[i] = d_svd_replace_value;
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
                    mat_data[j * mat_size + i] = 0.0;
                }
                else
                {
                    mat_data[j * mat_size + i] = z[j * mat_size + i] / sqrt(w[j]);
                }
            }
        }

        plog << "DirectMobilityMatrix::factorizeDenseMatrix(): For " << err_msg << " matrix: " << counter
             << " eigenvalues for dense matrix with handle " << mat_name
             << "have been changed. Number of zero eigenvalues placed are " << counter_zero << std::endl;
    }
    else
    {
        TBOX_ERROR("DirectMobilityMatrix::factorizeDenseMatrix(): Unsupported dense "
                   << "matrix inversion method called for "
                   << err_msg
                   << std::endl);
    }

    return;
} // factorizeDenseMatrix

void
DirectMobilitySolver::computeSolution(Mat& mat, const MobilityMatrixInverseType& inv_type, int* ipiv, double* rhs)
{
    // Get pointer to matrix.
    int mat_size;
    double* mat_data = NULL;
    MatGetSize(mat, &mat_size, NULL);
    MatDenseGetArray(mat, &mat_data);

    int err = 0;
    if (inv_type == LAPACK_CHOLESKY)
    {
        dpotrs_((char*)"L", mat_size, 1, mat_data, mat_size, rhs, mat_size, err);
        if (err)
        {
            TBOX_ERROR("DirectMobilitySolver::computeSolution(). Solution failed using "
                       << "LAPACK CHOLESKY with error code "
                       << err
                       << std::endl);
        }
    }
    else if (inv_type == LAPACK_LU)
    {
        dgetrs_((char*)"N", mat_size, 1, mat_data, mat_size, ipiv, rhs, mat_size, err);

        if (err)
        {
            TBOX_ERROR("DirectMobilitySolver::computeSolution(). Solution failed using "
                       << "LAPACK LU with error code "
                       << err
                       << std::endl);
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
                temp[i] += mat_data[i * mat_size + j] * rhs[j];
            }
        }

        for (int i = 0; i < mat_size; ++i)
        {
            rhs[i] = 0.0;
            for (int j = 0; j < mat_size; ++j)
            {
                rhs[i] += mat_data[j * mat_size + i] * temp[j];
            }
        }
    }
    else
    {
        TBOX_ERROR("DirectMobilitySolver::computeSolution(). Inverse method not supported." << std::endl);
    }

    MatDenseRestoreArray(mat, &mat_data);

    return;
} // computeSolution

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR
