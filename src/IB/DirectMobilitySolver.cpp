// ---------------------------------------------------------------------
//
// Copyright (c) 2015 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/CIBStrategy.h"
#include "ibamr/DirectMobilitySolver.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/IBTK_MPI.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "ibtk/ibtk_utilities.h"

#include "CartesianGridGeometry.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

#include "petscmat.h"
#include "petscvec.h"
#include "petscviewer.h"
#include "petscviewertypes.h"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <map>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

extern "C"
{
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
} // namespace

////////////////////////////// PUBLIC ////////////////////////////////////////

DirectMobilitySolver::DirectMobilitySolver(std::string object_name,
                                           Pointer<Database> input_db,
                                           Pointer<CIBStrategy> cib_strategy)
    : d_object_name(std::move(object_name)), d_cib_strategy(cib_strategy)
{
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
    for (const auto& petsc_mat_pair : d_petsc_mat_map)
    {
        const std::string& mat_name = petsc_mat_pair.first;
        Mat& mobility_mat = d_petsc_mat_map[mat_name].first;
        Mat& body_mobility_mat = d_petsc_mat_map[mat_name].second;
        MatDestroy(&mobility_mat);
        MatDestroy(&body_mobility_mat);
    }

    for (const auto& mat_pair : d_petsc_geometric_mat_map)
    {
        Mat& geometric_mat = d_petsc_geometric_mat_map[mat_pair.first];
        MatDestroy(&geometric_mat);
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
    for (const auto& prototype_struct_id : prototype_struct_ids)
    {
        TBOX_ASSERT(prototype_struct_id < d_cib_strategy->getNumberOfRigidStructures());
    }
    TBOX_ASSERT(d_mat_map.find(mat_name) == d_mat_map.end());
    TBOX_ASSERT(mat_type != UNKNOWN_MOBILITY_MATRIX_TYPE);
    TBOX_ASSERT(inv_type.first != UNKNOWN_MOBILITY_MATRIX_INVERSE_TYPE);
    TBOX_ASSERT(inv_type.second != UNKNOWN_MOBILITY_MATRIX_INVERSE_TYPE);
#endif

    unsigned int num_nodes = 0;
    for (const auto& prototype_struct_id : prototype_struct_ids)
    {
        num_nodes += d_cib_strategy->getNumberOfNodes(prototype_struct_id);
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
    d_mat_map[mat_name] = { {}, {} };
    d_geometric_mat_map[mat_name] = {};
    d_ipiv_map[mat_name] = { {}, {} };
    d_petsc_mat_map[mat_name] = { nullptr, nullptr };
    d_petsc_geometric_mat_map[mat_name] = nullptr;

    // Allocate the actual matrices.
    const int mobility_mat_size = num_nodes * NDIM;
    const int body_mobility_mat_size = d_mat_parts_map[mat_name] * s_max_free_dofs;
    const int rank = IBTK_MPI::getRank();

    if (rank == managing_proc)
    {
        d_mat_map[mat_name].first.resize(mobility_mat_size * mobility_mat_size);
        MatCreateSeqDense(PETSC_COMM_SELF,
                          mobility_mat_size,
                          mobility_mat_size,
                          d_mat_map[mat_name].first.data(),
                          &d_petsc_mat_map[mat_name].first);

        d_mat_map[mat_name].second.resize(body_mobility_mat_size * body_mobility_mat_size);
        MatCreateSeqDense(PETSC_COMM_SELF,
                          body_mobility_mat_size,
                          body_mobility_mat_size,
                          d_mat_map[mat_name].second.data(),
                          &d_petsc_mat_map[mat_name].second);

        d_geometric_mat_map[mat_name].resize(mobility_mat_size * body_mobility_mat_size);
        MatCreateSeqDense(PETSC_COMM_SELF,
                          mobility_mat_size,
                          body_mobility_mat_size,
                          d_geometric_mat_map[mat_name].data(),
                          &d_petsc_geometric_mat_map[mat_name]);

        if (d_mat_inv_type_map[mat_name].first == LAPACK_LU)
        {
            d_ipiv_map[mat_name].first.resize(mobility_mat_size);
        }
        if (d_mat_inv_type_map[mat_name].second == LAPACK_LU)
        {
            d_ipiv_map[mat_name].second.resize(body_mobility_mat_size);
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
    for (const auto& struct_id : struct_ids)
    {
        TBOX_ASSERT(struct_id.size() == d_mat_prototype_id_map[mat_name].size());
        unsigned num_nodes = 0;
        for (const auto& id : struct_id)
        {
            TBOX_ASSERT(id < d_cib_strategy->getNumberOfRigidStructures());

            num_nodes += d_cib_strategy->getNumberOfNodes(id);
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

    const int rank = IBTK_MPI::getRank();
    static const int data_depth = NDIM;

    for (const auto& petsc_mat_pair : d_petsc_mat_map)
    {
        const std::string& mat_name = petsc_mat_pair.first;
        Mat& mat = d_petsc_mat_map[mat_name].first;
        const MobilityMatrixInverseType& inv_type = d_mat_inv_type_map[mat_name].first;
        const std::vector<std::vector<unsigned> >& struct_ids = d_mat_actual_id_map[mat_name];
        const int managing_proc = d_mat_proc_map[mat_name];
        const int mat_size = d_mat_nodes_map[mat_name] * data_depth;
        const int num_structs = static_cast<int>(struct_ids.size());

        for (int k = 0; k < num_structs; ++k)
        {
            std::vector<double> rhs;
            if (rank == managing_proc) rhs.resize(mat_size);
            d_cib_strategy->copyVecToArray(b, rhs.data(), struct_ids[k], data_depth, managing_proc);
            if (!d_recompute_mob_mat)
            {
                d_cib_strategy->rotateArray(rhs.data(),
                                            struct_ids[k],
                                            /*use_transpose*/ true,
                                            managing_proc,
                                            data_depth);
            }
            if (rank == managing_proc) computeSolution(mat, inv_type, d_ipiv_map[mat_name].first.data(), rhs.data());
            if (!d_recompute_mob_mat)
            {
                d_cib_strategy->rotateArray(rhs.data(),
                                            struct_ids[k],
                                            /*use_transpose*/ false,
                                            managing_proc,
                                            data_depth);
            }
            d_cib_strategy->copyArrayToVec(x, rhs.data(), struct_ids[k], data_depth, managing_proc);
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

    const int rank = IBTK_MPI::getRank();
    static const int data_depth = s_max_free_dofs;

    for (const auto& petsc_mat_pair : d_petsc_mat_map)
    {
        const std::string& mat_name = petsc_mat_pair.first;
        Mat& mat = d_petsc_mat_map[mat_name].second;
        const MobilityMatrixInverseType& inv_type = d_mat_inv_type_map[mat_name].second;
        const std::vector<std::vector<unsigned> >& struct_ids = d_mat_actual_id_map[mat_name];
        const int mat_size = d_mat_parts_map[mat_name] * data_depth;
        const int managing_proc = d_mat_proc_map[mat_name];
        const int num_structs = static_cast<int>(struct_ids.size());

        for (int k = 0; k < num_structs; ++k)
        {
            std::vector<double> rhs;
            if (rank == managing_proc) rhs.resize(mat_size);
            d_cib_strategy->copyFreeDOFsVecToArray(b, rhs.data(), struct_ids[k], managing_proc);
            if (!d_recompute_mob_mat)
            {
                d_cib_strategy->rotateArray(rhs.data(),
                                            struct_ids[k],
                                            /*use_transpose*/ true,
                                            managing_proc,
                                            data_depth);
            }
            if (rank == managing_proc) computeSolution(mat, inv_type, d_ipiv_map[mat_name].second.data(), rhs.data());
            if (!d_recompute_mob_mat)
            {
                d_cib_strategy->rotateArray(rhs.data(),
                                            struct_ids[k],
                                            /*use_transpose*/ false,
                                            managing_proc,
                                            data_depth);
            }
            d_cib_strategy->copyFreeDOFsArrayToVec(x, rhs.data(), struct_ids[k], managing_proc);
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

    int rank = IBTK_MPI::getRank();
    auto managed_mats = static_cast<unsigned>(d_mat_map.size());

    static bool recreate_mobility_matrices = true;
    static std::vector<bool> read_files(managed_mats, false);
    bool initial_time = !d_recompute_mob_mat;

    if (recreate_mobility_matrices)
    {
        // Get grid-info
        Vec* vx;
        VecNestGetSubVecs(x, nullptr, &vx);
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
        for (auto it = d_petsc_mat_map.begin(); it != d_petsc_mat_map.end(); ++it, ++file_counter)
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
    comp_db = input_db->isDatabase("LAPACK_SVD") ? input_db->getDatabase("LAPACK_SVD") : Pointer<Database>(nullptr);
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
    int rank = IBTK_MPI::getRank();
    for (const auto& petsc_mat_pair : d_petsc_mat_map)
    {
        const std::string& mat_name = petsc_mat_pair.first;
        if (rank != d_mat_proc_map[mat_name]) continue;

        Mat& mat = d_petsc_mat_map[mat_name].first;
        const MobilityMatrixInverseType& inv_type = d_mat_inv_type_map[mat_name].first;
        const int mat_size = d_mat_nodes_map[mat_name] * NDIM;
        double* mat_data = nullptr;
        MatDenseGetArray(mat, &mat_data);
        factorizeDenseMatrix(mat_data, mat_size, inv_type, d_ipiv_map[mat_name].first.data(), mat_name, "Mobility");
        MatDenseRestoreArray(mat, &mat_data);
    }
    return;

} // factorizeMobilityMatrix

void
DirectMobilitySolver::constructBodyMobilityMatrix()
{
    int rank = IBTK_MPI::getRank();
    for (const auto& petsc_mat_pair : d_petsc_mat_map)
    {
        const std::string& mat_name = petsc_mat_pair.first;
        if (rank != d_mat_proc_map[mat_name]) continue;

        const int row_size = d_mat_nodes_map[mat_name] * NDIM;
        const int col_size = d_mat_parts_map[mat_name] * s_max_free_dofs;
        const MobilityMatrixInverseType& mobility_inv_type = d_mat_inv_type_map[mat_name].first;

        Mat& mobility_mat = d_petsc_mat_map[mat_name].first;
        Mat& body_mob_mat = d_petsc_mat_map[mat_name].second;
        Mat& geometric_mat = d_petsc_geometric_mat_map[mat_name];

        // Allocate a temporary matrix that holds the Matrix-Matrix product.
        // Here we are multiplying inverse of mobility matrix with geometric matrix.
        std::vector<double> product_mat_data(row_size * col_size);
        Mat product_mat;
        MatCreateSeqDense(PETSC_COMM_SELF, row_size, col_size, product_mat_data.data(), &product_mat);
        MatCopy(geometric_mat, product_mat, SAME_NONZERO_PATTERN);

        for (int col = 0; col < col_size; ++col)
        {
            double* col_data;
            MatDenseGetArray(product_mat, &col_data);
            computeSolution(
                mobility_mat, mobility_inv_type, d_ipiv_map[mat_name].first.data(), &col_data[col * row_size]);
            MatDenseRestoreArray(product_mat, &col_data);
        }
        MatTransposeMatMult(geometric_mat, product_mat, MAT_REUSE_MATRIX, PETSC_DEFAULT, &body_mob_mat);

        MatDestroy(&product_mat);
    }

    return;
} // generateBodyFrictionMatrix

void
DirectMobilitySolver::factorizeBodyMobilityMatrix()
{
    int rank = IBTK_MPI::getRank();
    for (const auto& petsc_mat_pair : d_petsc_mat_map)
    {
        const std::string& mat_name = petsc_mat_pair.first;
        if (rank != d_mat_proc_map[mat_name]) continue;

        Mat& mat = d_petsc_mat_map[mat_name].second;
        const MobilityMatrixInverseType& inv_type = d_mat_inv_type_map[mat_name].second;
        const int mat_size = d_mat_parts_map[mat_name] * s_max_free_dofs;

        double* mat_data = nullptr;
        MatDenseGetArray(mat, &mat_data);
        factorizeDenseMatrix(
            mat_data, mat_size, inv_type, d_ipiv_map[mat_name].second.data(), mat_name, "Body Mobility");
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
                                                                        << " failed for matrix handle " << mat_name
                                                                        << " with error code " << err
                                                                        << " using LAPACK CHOLESKY." << std::endl);
        }
    }
    else if (inv_type == LAPACK_LU)
    {
        dgetrf_(mat_size, mat_size, mat_data, mat_size, ipiv, err);
        if (err)
        {
            TBOX_ERROR("DirectMobilityMatrix::factorizeDenseMatrix(). "
                       << err_msg << " matrix factorization "
                       << "failed for matrix handle " << mat_name << " with error code " << err << " using LAPACK LU."
                       << std::endl);
        }
    }
    else if (inv_type == LAPACK_SVD)
    {
        // Use the symmetric eigenvalue decomposition as a stand-in for the SVD.
        // In particular, since A = V D V^T, where V is the matrix of eigenvectors
        // and D is the matrix of eigenvalues, we can factorize A as
        //
        //   A = V sqrt(D) sqrt(D) V^T
        //     = V sqrt(D) (V sqrt(D))^T
        //
        // and instead store A <- V sqrt(D)
        using MatrixType = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
        // Older versions of Eigen don't make Eigen::Index publicly available
#if EIGEN_VERSION_AT_LEAST(3, 3, 0)
        using IndexType = Eigen::Index;
#else
        using IndexType = typename MatrixType::Index;
#endif
        Eigen::Map<MatrixType> mat_view(mat_data, IndexType(mat_size), IndexType(mat_size));
        // For compatibility with the old code, we copy the lower triangle into
        // the upper triangle, even if they aren't actually equal
        for (int i = 0; i < mat_size; ++i)
        {
            for (int j = i + 1; j < mat_size; ++j)
            {
                mat_view(i, j) = mat_view(j, i);
            }
        }
        Eigen::SelfAdjointEigenSolver<MatrixType> eigensolver(mat_view);
        Eigen::Matrix<double, Eigen::Dynamic, 1> eigenvalues = eigensolver.eigenvalues();
        const MatrixType eigenvectors = eigensolver.eigenvectors();
        // Make negative eigenvalues to be equal to min eigen value from
        // input option
        int counter = 0, counter_zero = 0;
        for (int i = 0; i < mat_size; ++i)
        {
            if (eigenvalues[i] < d_svd_eps)
            {
                eigenvalues[i] = d_svd_replace_value;
                counter++;
            }
        }
        for (int i = 0; i < mat_size; ++i)
        {
            if (MathUtilities<double>::equalEps(eigenvalues[i], 0.0))
            {
                counter_zero++;
            }

            for (int j = 0; j < mat_size; ++j)
            {
                if (MathUtilities<double>::equalEps(eigenvalues[j], 0.0))
                {
                    mat_view(i, j) = 0.0;
                }
                else
                {
                    mat_view(i, j) = eigenvectors(i, j) / std::sqrt(eigenvalues[j]);
                }
            }
        }

        plog << "DirectMobilityMatrix::factorizeDenseMatrix(): For " << err_msg << " matrix: " << counter
             << " eigenvalues for dense matrix with handle " << mat_name
             << " have been changed. Number of zero eigenvalues placed are " << counter_zero << std::endl;
    }
    else
    {
        TBOX_ERROR("DirectMobilityMatrix::factorizeDenseMatrix(): Unsupported dense "
                   << "matrix inversion method called for " << err_msg << std::endl);
    }

    return;
} // factorizeDenseMatrix

void
DirectMobilitySolver::computeSolution(Mat& mat, const MobilityMatrixInverseType& inv_type, int* ipiv, double* rhs)
{
    // Get pointer to matrix.
    int mat_size;
    double* mat_data = nullptr;
    MatGetSize(mat, &mat_size, nullptr);
    MatDenseGetArray(mat, &mat_data);

    int err = 0;
    if (inv_type == LAPACK_CHOLESKY)
    {
        dpotrs_((char*)"L", mat_size, 1, mat_data, mat_size, rhs, mat_size, err);
        if (err)
        {
            TBOX_ERROR("DirectMobilitySolver::computeSolution(). Solution failed using "
                       << "LAPACK CHOLESKY with error code " << err << std::endl);
        }
    }
    else if (inv_type == LAPACK_LU)
    {
        dgetrs_((char*)"N", mat_size, 1, mat_data, mat_size, ipiv, rhs, mat_size, err);

        if (err)
        {
            TBOX_ERROR("DirectMobilitySolver::computeSolution(). Solution failed using "
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
