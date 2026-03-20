// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2025 by the IBAMR developers
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

#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PETScLevelSolver.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Database.h>
#include <tbox/PIO.h>
#include <tbox/Pointer.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

#include <petscis.h>
#include <petscistypes.h>
#include <petsclog.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscpctypes.h>
#include <petscsys.h>
#include <petscvec.h>
#include <petscversion.h>

#include <Eigen/SVD>

#include <CoarseFineBoundary.h>
#include <IntVector.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SAMRAIVectorReal.h>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_solve_system;
static Timer* t_initialize_solver_state;
static Timer* t_deallocate_solver_state;

void
destroy_index_sets(std::vector<IS>& index_sets)
{
    int ierr;
    for (auto& index_set : index_sets)
    {
        ierr = ISDestroy(&index_set);
        IBTK_CHKERRQ(ierr);
    }
    index_sets.clear();
    return;
} // destroy_index_sets

void
build_petsc_subdomain_index_sets(std::vector<IS>& subdomain_is,
                                 std::vector<IS>& nonoverlap_subdomain_is,
                                 const std::vector<std::vector<int>>& subdomain_dofs,
                                 const std::vector<std::vector<int>>& nonoverlap_subdomain_dofs)
{
    int ierr;
    destroy_index_sets(subdomain_is);
    destroy_index_sets(nonoverlap_subdomain_is);

    subdomain_is.resize(subdomain_dofs.size());
    for (std::size_t subdomain_num = 0; subdomain_num < subdomain_dofs.size(); ++subdomain_num)
    {
        PetscInt* dof_arr = nullptr;
        const PetscInt n_dofs = static_cast<PetscInt>(subdomain_dofs[subdomain_num].size());
        ierr = PetscMalloc1(n_dofs, &dof_arr);
        IBTK_CHKERRQ(ierr);
        std::copy(subdomain_dofs[subdomain_num].begin(), subdomain_dofs[subdomain_num].end(), dof_arr);
        ierr = ISCreateGeneral(PETSC_COMM_SELF, n_dofs, dof_arr, PETSC_OWN_POINTER, &subdomain_is[subdomain_num]);
        IBTK_CHKERRQ(ierr);
    }

    nonoverlap_subdomain_is.resize(nonoverlap_subdomain_dofs.size());
    for (std::size_t subdomain_num = 0; subdomain_num < nonoverlap_subdomain_dofs.size(); ++subdomain_num)
    {
        PetscInt* dof_arr = nullptr;
        const PetscInt n_dofs = static_cast<PetscInt>(nonoverlap_subdomain_dofs[subdomain_num].size());
        ierr = PetscMalloc1(n_dofs, &dof_arr);
        IBTK_CHKERRQ(ierr);
        std::copy(
            nonoverlap_subdomain_dofs[subdomain_num].begin(), nonoverlap_subdomain_dofs[subdomain_num].end(), dof_arr);
        ierr = ISCreateGeneral(
            PETSC_COMM_SELF, n_dofs, dof_arr, PETSC_OWN_POINTER, &nonoverlap_subdomain_is[subdomain_num]);
        IBTK_CHKERRQ(ierr);
    }
    return;
} // build_petsc_subdomain_index_sets

template <class SolverType>
Eigen::MatrixXd
build_qr_solve_matrix(const Eigen::MatrixXd& local_operator, const double threshold)
{
    SolverType solver(local_operator);
    if (threshold >= 0.0) solver.setThreshold(threshold);
    return solver.solve(Eigen::MatrixXd::Identity(local_operator.rows(), local_operator.cols()));
} // build_qr_solve_matrix

Eigen::MatrixXd
build_llt_solve_matrix(const Eigen::MatrixXd& local_operator)
{
    Eigen::LLT<Eigen::MatrixXd> solver(local_operator);
    return solver.solve(Eigen::MatrixXd::Identity(local_operator.rows(), local_operator.cols()));
} // build_llt_solve_matrix

Eigen::MatrixXd
build_ldlt_solve_matrix(const Eigen::MatrixXd& local_operator)
{
    Eigen::LDLT<Eigen::MatrixXd> solver(local_operator);
    return solver.solve(Eigen::MatrixXd::Identity(local_operator.rows(), local_operator.cols()));
} // build_ldlt_solve_matrix

Eigen::MatrixXd
build_householder_qr_solve_matrix(const Eigen::MatrixXd& local_operator)
{
    Eigen::HouseholderQR<Eigen::MatrixXd> solver(local_operator);
    return solver.solve(Eigen::MatrixXd::Identity(local_operator.rows(), local_operator.cols()));
} // build_householder_qr_solve_matrix

Eigen::MatrixXd
build_partial_piv_lu_solve_matrix(const Eigen::MatrixXd& local_operator)
{
    Eigen::PartialPivLU<Eigen::MatrixXd> solver(local_operator);
    return solver.solve(Eigen::MatrixXd::Identity(local_operator.rows(), local_operator.cols()));
} // build_partial_piv_lu_solve_matrix

Eigen::MatrixXd
build_full_piv_lu_solve_matrix(const Eigen::MatrixXd& local_operator, const double threshold)
{
    Eigen::FullPivLU<Eigen::MatrixXd> solver(local_operator);
    if (threshold >= 0.0) solver.setThreshold(threshold);
    return solver.solve(Eigen::MatrixXd::Identity(local_operator.rows(), local_operator.cols()));
} // build_full_piv_lu_solve_matrix

template <class SolverType>
struct eigen_solver_supports_threshold : std::false_type
{
};

template <>
struct eigen_solver_supports_threshold<Eigen::FullPivLU<Eigen::MatrixXd>> : std::true_type
{
};

template <>
struct eigen_solver_supports_threshold<Eigen::ColPivHouseholderQR<Eigen::MatrixXd>> : std::true_type
{
};

template <>
struct eigen_solver_supports_threshold<Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd>> : std::true_type
{
};

template <>
struct eigen_solver_supports_threshold<Eigen::FullPivHouseholderQR<Eigen::MatrixXd>> : std::true_type
{
};

template <>
struct eigen_solver_supports_threshold<Eigen::JacobiSVD<Eigen::MatrixXd>> : std::true_type
{
};

template <>
struct eigen_solver_supports_threshold<Eigen::BDCSVD<Eigen::MatrixXd>> : std::true_type
{
};

template <class SolverType>
struct eigen_solver_uses_svd_compute : std::false_type
{
};

template <>
struct eigen_solver_uses_svd_compute<Eigen::JacobiSVD<Eigen::MatrixXd>> : std::true_type
{
};

template <>
struct eigen_solver_uses_svd_compute<Eigen::BDCSVD<Eigen::MatrixXd>> : std::true_type
{
};

template <class SolverType>
void
initialize_eigen_solver(SolverType& solver, const Eigen::MatrixXd& local_operator, const double threshold)
{
    if constexpr (eigen_solver_uses_svd_compute<SolverType>::value)
    {
        solver.compute(local_operator, Eigen::ComputeThinU | Eigen::ComputeThinV);
    }
    else
    {
        solver.compute(local_operator);
    }

    if constexpr (eigen_solver_supports_threshold<SolverType>::value)
    {
        if (threshold >= 0.0) solver.setThreshold(threshold);
    }
}

template <class SolverType, class RhsType>
auto
solve_eigen_system(const SolverType& solver, const RhsType& rhs) -> decltype(solver.solve(rhs))
{
    return solver.solve(rhs);
}

template <class SVDType>
Eigen::MatrixXd
build_svd_pseudoinverse(const Eigen::MatrixXd& local_operator, const double threshold)
{
    SVDType svd(local_operator, Eigen::ComputeThinU | Eigen::ComputeThinV);
    if (threshold >= 0.0) svd.setThreshold(threshold);

    const auto& singular_values = svd.singularValues();
    Eigen::VectorXd inv_singular_values = Eigen::VectorXd::Zero(singular_values.size());
    if (singular_values.size() > 0)
    {
        const double cutoff = singular_values[0] * svd.threshold();
        for (Eigen::Index k = 0; k < singular_values.size(); ++k)
        {
            if (singular_values[k] > cutoff) inv_singular_values[k] = 1.0 / singular_values[k];
        }
    }

    return svd.matrixV() * inv_singular_values.asDiagonal() * svd.matrixU().adjoint();
} // build_svd_pseudoinverse

Eigen::MatrixXd
build_complete_orthogonal_decomposition_pseudoinverse(const Eigen::MatrixXd& local_operator, const double threshold)
{
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> solver(local_operator);
    if (threshold >= 0.0) solver.setThreshold(threshold);
    Eigen::MatrixXd pseudoinverse = solver.pseudoInverse();
    return pseudoinverse;
} // build_complete_orthogonal_decomposition_pseudoinverse
} // namespace

Eigen::SparseMatrix<double, Eigen::RowMajor>
PETScLevelSolver::copyPETScMatToEigenSparse(Mat mat)
{
    PetscInt m = 0, n = 0;
    int ierr = MatGetSize(mat, &m, &n);
    IBTK_CHKERRQ(ierr);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<std::size_t>(m) * 16);
    for (PetscInt row = 0; row < m; ++row)
    {
        PetscInt row_nnz = 0;
        const PetscInt* cols = nullptr;
        const PetscScalar* vals = nullptr;
        ierr = MatGetRow(mat, row, &row_nnz, &cols, &vals);
        IBTK_CHKERRQ(ierr);
        for (PetscInt k = 0; k < row_nnz; ++k)
        {
            triplets.emplace_back(static_cast<int>(row), static_cast<int>(cols[k]), PetscRealPart(vals[k]));
        }
        ierr = MatRestoreRow(mat, row, &row_nnz, &cols, &vals);
        IBTK_CHKERRQ(ierr);
    }
    Eigen::SparseMatrix<double, Eigen::RowMajor> A(static_cast<int>(m), static_cast<int>(n));
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();
    return A;
} // copyPETScMatToEigenSparse

/////////////////////////////// PUBLIC ///////////////////////////////////////

PETScLevelSolver::PETScLevelSolver()
{
    // Setup default options.
    d_max_iterations = 10000;
    d_abs_residual_tol = 1.0e-50;
    d_rel_residual_tol = 1.0e-5;
    d_initial_guess_nonzero = true;
    d_enable_logging = false;

    // Setup Timers.
    IBTK_DO_ONCE(t_solve_system = TimerManager::getManager()->getTimer("IBTK::PETScLevelSolver::solveSystem()");
                 t_initialize_solver_state =
                     TimerManager::getManager()->getTimer("IBTK::PETScLevelSolver::initializeSolverState()");
                 t_deallocate_solver_state =
                     TimerManager::getManager()->getTimer("IBTK::PETScLevelSolver::deallocateSolverState()"););
    return;
} // PETScLevelSolver

PETScLevelSolver::~PETScLevelSolver()
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::~PETScLevelSolver()\n"
                                 << "  subclass must call deallocateSolverState in subclass destructor" << std::endl);
    }
    destroy_index_sets(d_petsc_nonoverlap_subdomain_is);
    destroy_index_sets(d_petsc_subdomain_is);
    destroy_index_sets(d_field_is);
    return;
} // ~PETScLevelSolver

void
PETScLevelSolver::setKSPType(const std::string& ksp_type)
{
    d_ksp_type = ksp_type;
    return;
} // setKSPType

void
PETScLevelSolver::setOptionsPrefix(const std::string& options_prefix)
{
    d_options_prefix = options_prefix;
    return;
} // setOptionsPrefix

const KSP&
PETScLevelSolver::getPETScKSP() const
{
    return d_petsc_ksp;
} // getPETScKSP

void
PETScLevelSolver::getASMSubdomains(std::vector<std::vector<int>>** nonoverlap_subdomain_dofs,
                                   std::vector<std::vector<int>>** subdomain_dofs)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_is_initialized);
#endif
    if (nonoverlap_subdomain_dofs) *nonoverlap_subdomain_dofs = &d_nonoverlap_subdomain_dofs;
    if (subdomain_dofs) *subdomain_dofs = &d_subdomain_dofs;
    return;
} // getASMSubdomains

void
PETScLevelSolver::setNullSpace(bool contains_constant_vec,
                               const std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>>& nullspace_basis_vecs)
{
    LinearSolver::setNullSpace(contains_constant_vec, nullspace_basis_vecs);
    if (d_is_initialized) setupNullSpace();
    return;
} // setNullSpace

bool
PETScLevelSolver::solveSystem(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& b)
{
    IBTK_TIMER_START(t_solve_system);

    int ierr;

    if (d_enable_logging) plog << d_object_name << "::solveSystem():" << std::endl;

    // Initialize the solver, when necessary.
    const bool deallocate_after_solve = !d_is_initialized;
    if (deallocate_after_solve) initializeSolverState(x, b);

    // Configure solver.
    ierr = KSPSetTolerances(d_petsc_ksp, d_rel_residual_tol, d_abs_residual_tol, PETSC_DEFAULT, d_max_iterations);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(d_petsc_ksp, d_initial_guess_nonzero ? PETSC_TRUE : PETSC_FALSE);
    IBTK_CHKERRQ(ierr);

    // Solve the system.
    setupKSPVecs(d_petsc_x, d_petsc_b, x, b);
    ierr = KSPSolve(d_petsc_ksp, d_petsc_b, d_petsc_x);
    IBTK_CHKERRQ(ierr);
    copyFromPETScVec(d_petsc_x, x);

    // Log solver info.
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(d_petsc_ksp, &reason);
    IBTK_CHKERRQ(ierr);
    const bool converged = reason > 0;
    if (d_enable_logging)
    {
        plog << d_object_name << "::solveSystem(): solver " << (converged ? "converged" : "diverged") << "\n"
             << "iterations = " << d_current_iterations << "\n"
             << "residual norm = " << d_current_residual_norm << std::endl;
    }

    // Deallocate the solver, when necessary.
    if (deallocate_after_solve) deallocateSolverState();

    IBTK_TIMER_STOP(t_solve_system);
    return converged;
} // solveSystem

void
PETScLevelSolver::initializeSolverState(const SAMRAIVectorReal<NDIM, double>& x,
                                        const SAMRAIVectorReal<NDIM, double>& b)
{
    IBTK_TIMER_START(t_initialize_solver_state);

// Rudimentary error checking.
#if !defined(NDEBUG)
    if (x.getNumberOfComponents() != b.getNumberOfComponents())
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  vectors must have the same number of components" << std::endl);
    }

    const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy = x.getPatchHierarchy();
    if (patch_hierarchy != b.getPatchHierarchy())
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  vectors must have the same hierarchy" << std::endl);
    }

    const int coarsest_ln = x.getCoarsestLevelNumber();
    if (coarsest_ln < 0)
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  coarsest level number must not be negative" << std::endl);
    }
    if (coarsest_ln != b.getCoarsestLevelNumber())
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  vectors must have same coarsest level number" << std::endl);
    }

    const int finest_ln = x.getFinestLevelNumber();
    if (finest_ln < coarsest_ln)
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  finest level number must be >= coarsest level number" << std::endl);
    }
    if (finest_ln != b.getFinestLevelNumber())
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  vectors must have same finest level number" << std::endl);
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!patch_hierarchy->getPatchLevel(ln))
        {
            TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                     << "  hierarchy level " << ln << " does not exist" << std::endl);
        }
    }

    if (coarsest_ln != finest_ln)
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState()\n"
                                 << "  coarsest_ln != finest_ln in PETScLevelSolver" << std::endl);
    }
#endif
    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized) deallocateSolverState();

    // Get the hierarchy information.
    d_hierarchy = x.getPatchHierarchy();
    d_level_num = x.getCoarsestLevelNumber();
    TBOX_ASSERT(d_level_num == x.getFinestLevelNumber());
    d_level = d_hierarchy->getPatchLevel(d_level_num);
    if (d_level_num > 0)
    {
        d_cf_boundary = new CoarseFineBoundary<NDIM>(*d_hierarchy, d_level_num, IntVector<NDIM>(1));
    }

    // Setup data cache.
    d_cached_eulerian_data.setPatchHierarchy(d_hierarchy);
    d_cached_eulerian_data.resetLevels(d_level_num, d_level_num);

    // Perform specialized operations to initialize solver state();
    initializeSolverStateSpecialized(x, b);

    // Setup PETSc objects.
    int ierr;
    ierr = KSPCreate(PETSC_COMM_WORLD, &d_petsc_ksp);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetOperators(d_petsc_ksp, d_petsc_mat, d_petsc_pc);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetReusePreconditioner(d_petsc_ksp, PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetType(d_petsc_ksp, d_ksp_type.c_str());
    IBTK_CHKERRQ(ierr);
    PetscBool initial_guess_nonzero = d_initial_guess_nonzero ? PETSC_TRUE : PETSC_FALSE;
    ierr = KSPSetInitialGuessNonzero(d_petsc_ksp, initial_guess_nonzero);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetTolerances(d_petsc_ksp, d_rel_residual_tol, d_abs_residual_tol, PETSC_DEFAULT, d_max_iterations);
    IBTK_CHKERRQ(ierr);

    // Setup KSP PC.
    PC ksp_pc;
    ierr = KSPGetPC(d_petsc_ksp, &ksp_pc);
    IBTK_CHKERRQ(ierr);
    PCType pc_type = d_pc_type.c_str();
    ierr = PCSetType(ksp_pc, pc_type);
    IBTK_CHKERRQ(ierr);
    if (d_options_prefix != "")
    {
        ierr = KSPSetOptionsPrefix(d_petsc_ksp, d_options_prefix.c_str());
        IBTK_CHKERRQ(ierr);
    }
    ierr = KSPSetFromOptions(d_petsc_ksp);
    IBTK_CHKERRQ(ierr);

    // Reset class data structure to correspond to command-line options.
    ierr = KSPGetTolerances(d_petsc_ksp, &d_rel_residual_tol, &d_abs_residual_tol, nullptr, &d_max_iterations);
    IBTK_CHKERRQ(ierr);
    ierr = PCGetType(ksp_pc, &pc_type);
    IBTK_CHKERRQ(ierr);
    d_pc_type = pc_type;
    d_preconditioner_type = parsePreconditionerType(d_pc_type);

    // Set the nullspace.
    if (d_nullspace_contains_constant_vec || !d_nullspace_basis_vecs.empty()) setupNullSpace();

    // Setup the preconditioner.
    switch (d_preconditioner_type)
    {
    case PreconditionerType::ASM:
        configureASMPreconditioner(ksp_pc);
        break;
    case PreconditionerType::FIELDSPLIT:
        configureFieldSplitPreconditioner(ksp_pc);
        break;
    case PreconditionerType::SHELL:
        configureShellPreconditioner(ksp_pc);
        break;
    case PreconditionerType::OTHER:
        break;
    }

    // Indicate that the solver is initialized.
    d_is_initialized = true;

    IBTK_TIMER_STOP(t_initialize_solver_state);
    return;
} // initializeSolverState

void
PETScLevelSolver::deallocateSolverState()
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_deallocate_solver_state);

    // Perform specialized operations to deallocate solver state.
    deallocateSolverStateSpecialized();

    // Deallocate PETSc objects.
    int ierr;
    ierr = KSPDestroy(&d_petsc_ksp);
    IBTK_CHKERRQ(ierr);
    if (d_petsc_mat != d_petsc_pc)
    {
        ierr = MatDestroy(&d_petsc_pc);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatDestroy(&d_petsc_mat);
    IBTK_CHKERRQ(ierr);
    if (d_nullspace_contains_constant_vec || !d_nullspace_basis_vecs.empty())
    {
        ierr = MatNullSpaceDestroy(&d_petsc_nullsp);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecDestroy(&d_petsc_x);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&d_petsc_b);
    IBTK_CHKERRQ(ierr);

    // Deallocate PETSc objects for shell preconditioner.
    deallocateShellData();
    destroy_index_sets(d_field_is);
    d_field_name.clear();
    // Indicate that the solver is NOT initialized.
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_solver_state);
    return;
} // deallocateSolverState

/////////////////////////////// PROTECTED ////////////////////////////////////

void
PETScLevelSolver::init(Pointer<Database> input_db, const std::string& default_options_prefix)
{
    d_options_prefix = default_options_prefix;
    if (input_db)
    {
        if (input_db->keyExists("options_prefix")) d_options_prefix = input_db->getString("options_prefix");
        if (input_db->keyExists("enable_logging")) d_enable_logging = input_db->getBool("enable_logging");
        if (input_db->keyExists("max_iterations")) d_max_iterations = input_db->getInteger("max_iterations");
        if (input_db->keyExists("abs_residual_tol")) d_abs_residual_tol = input_db->getDouble("abs_residual_tol");
        if (input_db->keyExists("rel_residual_tol")) d_rel_residual_tol = input_db->getDouble("rel_residual_tol");
        if (input_db->keyExists("ksp_type")) d_ksp_type = input_db->getString("ksp_type");
        if (input_db->keyExists("pc_type")) d_pc_type = input_db->getString("pc_type");
        if (input_db->keyExists("shell_pc_type")) d_shell_pc_type = input_db->getString("shell_pc_type");
        if (input_db->keyExists("eigen_subdomain_solver_type"))
        {
            d_eigen_subdomain_solver_type =
                parseEigenSubdomainSolverType(input_db->getString("eigen_subdomain_solver_type"));
        }
        if (input_db->keyExists("eigen_subdomain_solver_threshold"))
        {
            d_eigen_subdomain_solver_threshold = input_db->getDouble("eigen_subdomain_solver_threshold");
        }
        if (input_db->keyExists("eigen_subdomain_pseudoinverse_type"))
        {
            d_eigen_subdomain_pseudoinverse_type =
                parseEigenSubdomainSolverType(input_db->getString("eigen_subdomain_pseudoinverse_type"));
        }
        if (input_db->keyExists("eigen_subdomain_pseudoinverse_threshold"))
        {
            d_eigen_subdomain_pseudoinverse_threshold = input_db->getDouble("eigen_subdomain_pseudoinverse_threshold");
        }
        if (input_db->keyExists("initial_guess_nonzero"))
            d_initial_guess_nonzero = input_db->getBool("initial_guess_nonzero");
    }
    d_preconditioner_type = parsePreconditionerType(d_pc_type);
    configureShellSmootherType();
    return;
} // init

void
PETScLevelSolver::configureShellSmootherType()
{
    if (d_preconditioner_type != PreconditionerType::SHELL)
    {
        d_shell_pc_type.clear();
        d_shell_smoother_backend = ShellSmootherBackend::PETSC;
        d_shell_smoother_composition = ShellSmootherComposition::MULTIPLICATIVE;
        d_shell_smoother_partition = ShellSmootherPartition::BASIC;
        return;
    }

    d_shell_pc_type = normalizeShellSmootherType(d_shell_pc_type);
    if (d_shell_pc_type.empty()) return;

    d_shell_smoother_backend = parseShellSmootherBackend(d_shell_pc_type);
    d_shell_smoother_composition = parseShellSmootherComposition(d_shell_pc_type);
    d_shell_smoother_partition = parseShellSmootherPartition(d_shell_pc_type, d_shell_smoother_composition);
    return;
} // configureShellSmootherType

std::string
PETScLevelSolver::normalizeShellSmootherType(const std::string& type) const
{
    if (type == "additive-petsc") return "additive";
    if (type == "multiplicative-petsc") return "multiplicative";
    if (type == "additive-eigen-pinv") return "additive-eigen-pseudoinverse";
    if (type == "multiplicative-eigen-pinv") return "multiplicative-eigen-pseudoinverse";
    return type;
} // normalizeShellSmootherType

PETScLevelSolver::PreconditionerType
PETScLevelSolver::parsePreconditionerType(const std::string& type) const
{
    if (type == "asm") return PreconditionerType::ASM;
    if (type == "fieldsplit") return PreconditionerType::FIELDSPLIT;
    if (type == "shell") return PreconditionerType::SHELL;
    return PreconditionerType::OTHER;
} // parsePreconditionerType

PETScLevelSolver::ShellSmootherBackend
PETScLevelSolver::parseShellSmootherBackend(const std::string& type) const
{
    if (type.find("-eigen-pseudoinverse") != std::string::npos) return ShellSmootherBackend::EIGEN_PSEUDOINVERSE;
    return type.find("-eigen") != std::string::npos ? ShellSmootherBackend::EIGEN : ShellSmootherBackend::PETSC;
} // parseShellSmootherBackend

PETScLevelSolver::ShellSmootherComposition
PETScLevelSolver::parseShellSmootherComposition(const std::string& type) const
{
    return type.find("additive") == 0 ? ShellSmootherComposition::ADDITIVE : ShellSmootherComposition::MULTIPLICATIVE;
} // parseShellSmootherComposition

PETScLevelSolver::ShellSmootherPartition
PETScLevelSolver::parseShellSmootherPartition(const std::string& type, ShellSmootherComposition composition) const
{
    if (type.find("-restrict") != std::string::npos) return ShellSmootherPartition::RESTRICT;
    if (type.find("-basic") != std::string::npos) return ShellSmootherPartition::BASIC;
    return composition == ShellSmootherComposition::ADDITIVE ? ShellSmootherPartition::RESTRICT :
                                                               ShellSmootherPartition::BASIC;
} // parseShellSmootherPartition

PETScLevelSolver::EigenSubdomainSolverType
PETScLevelSolver::parseEigenSubdomainSolverType(const std::string& type) const
{
    std::string normalized = type;
    std::transform(normalized.begin(),
                   normalized.end(),
                   normalized.begin(),
                   [](const unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (normalized == "custom" || normalized == "staggered_stokes_custom" || normalized == "staggered-stokes-custom")
    {
        return EigenSubdomainSolverType::CUSTOM;
    }
    if (normalized == "llt")
    {
        return EigenSubdomainSolverType::LLT;
    }
    if (normalized == "ldlt")
    {
        return EigenSubdomainSolverType::LDLT;
    }
    if (normalized == "partialpivlu" || normalized == "partial-piv-lu" || normalized == "partial_piv_lu")
    {
        return EigenSubdomainSolverType::PARTIAL_PIV_LU;
    }
    if (normalized == "fullpivlu" || normalized == "full-piv-lu" || normalized == "full_piv_lu")
    {
        return EigenSubdomainSolverType::FULL_PIV_LU;
    }
    if (normalized == "householderqr" || normalized == "householder-qr" || normalized == "householder_qr")
    {
        return EigenSubdomainSolverType::HOUSEHOLDER_QR;
    }
    if (normalized == "colpivhouseholderqr" || normalized == "col-piv-householder-qr" ||
        normalized == "col_piv_householder_qr")
    {
        return EigenSubdomainSolverType::COL_PIV_HOUSEHOLDER_QR;
    }
    if (normalized == "completeorthogonaldecomposition" || normalized == "complete-orthogonal-decomposition" ||
        normalized == "complete_orthogonal_decomposition" || normalized == "cod")
    {
        return EigenSubdomainSolverType::COMPLETE_ORTHOGONAL_DECOMPOSITION;
    }
    if (normalized == "fullpivhouseholderqr" || normalized == "full-piv-householder-qr" ||
        normalized == "full_piv_householder_qr")
    {
        return EigenSubdomainSolverType::FULL_PIV_HOUSEHOLDER_QR;
    }
    if (normalized == "jacobisvd" || normalized == "jacobi-svd" || normalized == "jacobi_svd")
    {
        return EigenSubdomainSolverType::JACOBI_SVD;
    }
    if (normalized == "bdcsvd" || normalized == "bdc-svd" || normalized == "bdc_svd")
    {
        return EigenSubdomainSolverType::BDC_SVD;
    }

    TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::parseEigenSubdomainSolverType()\n"
                             << "Unknown Eigen subdomain solver type: " << type << std::endl);
    return EigenSubdomainSolverType::JACOBI_SVD;
} // parseEigenSubdomainSolverType

PETScLevelSolver::EigenSubdomainSolverType
PETScLevelSolver::getEigenSolveSolverType() const
{
    switch (d_shell_smoother_backend)
    {
    case ShellSmootherBackend::EIGEN:
        return d_eigen_subdomain_solver_type;
    case ShellSmootherBackend::EIGEN_PSEUDOINVERSE:
    case ShellSmootherBackend::PETSC:
        break;
    }

    TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::getEigenSolveSolverType()\n"
                             << "Eigen solve-based solver type requested for a non-solve shell backend." << std::endl);
    return EigenSubdomainSolverType::PARTIAL_PIV_LU;
} // getEigenSolveSolverType

PETScLevelSolver::EigenSubdomainSolverType
PETScLevelSolver::getEigenPseudoinverseSolverType() const
{
    switch (d_shell_smoother_backend)
    {
    case ShellSmootherBackend::EIGEN_PSEUDOINVERSE:
        return d_eigen_subdomain_pseudoinverse_type;
    case ShellSmootherBackend::EIGEN:
    case ShellSmootherBackend::PETSC:
        break;
    }

    TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::getEigenPseudoinverseSolverType()\n"
                             << "Eigen matrix-based solver type requested for a non-pseudoinverse shell backend."
                             << std::endl);
    return EigenSubdomainSolverType::COL_PIV_HOUSEHOLDER_QR;
} // getEigenPseudoinverseSolverType

bool
PETScLevelSolver::usingCustomEigenSolveSubdomainSolver() const
{
    return (d_shell_smoother_backend == ShellSmootherBackend::EIGEN &&
            d_eigen_subdomain_solver_type == EigenSubdomainSolverType::CUSTOM);
} // usingCustomEigenSolveSubdomainSolver

bool
PETScLevelSolver::usingCustomEigenPseudoinverseSubdomainSolver() const
{
    return (d_shell_smoother_backend == ShellSmootherBackend::EIGEN_PSEUDOINVERSE &&
            d_eigen_subdomain_pseudoinverse_type == EigenSubdomainSolverType::CUSTOM);
} // usingCustomEigenPseudoinverseSubdomainSolver

template <class SolverType>
PETScLevelSolver::EigenTypedSolveStorage<SolverType>&
PETScLevelSolver::getEigenSolveStorage()
{
    TBOX_ASSERT(d_shell_data.eigen_data);
    auto* storage = static_cast<EigenTypedSolveStorage<SolverType>*>(d_shell_data.eigen_data->solve_storage.get());
    TBOX_ASSERT(storage);
    return *storage;
}

template <class SolverType>
const PETScLevelSolver::EigenTypedSolveStorage<SolverType>&
PETScLevelSolver::getEigenSolveStorage() const
{
    TBOX_ASSERT(d_shell_data.eigen_data);
    const auto* storage =
        static_cast<const EigenTypedSolveStorage<SolverType>*>(d_shell_data.eigen_data->solve_storage.get());
    TBOX_ASSERT(storage);
    return *storage;
}

template <class SolverType>
void
PETScLevelSolver::initializeBuiltinEigenSolveStorage()
{
    auto& shell = d_shell_data;
    shell.eigen_data = std::make_unique<EigenShellSmootherData>();
    auto storage = std::make_unique<EigenTypedSolveStorage<SolverType>>();
    storage->solvers.resize(static_cast<std::size_t>(shell.n_local_subdomains));
    shell.eigen_data->solve_storage = std::move(storage);
    auto& typed_storage = getEigenSolveStorage<SolverType>();
    initializeEigenShellDataWithLocalOperatorHook(
        [this, &typed_storage](const Eigen::MatrixXd& local_operator, const std::size_t subdomain_num)
        {
            initialize_eigen_solver(
                typed_storage.solvers[subdomain_num], local_operator, d_eigen_subdomain_solver_threshold);
        });
}

template <class SolverType>
void
PETScLevelSolver::applyAdditiveEigenSolveImpl(Vec x, Vec y)
{
    TBOX_ASSERT(IBTK_MPI::getNodes() == 1);
    auto& eigen = *d_shell_data.eigen_data;
    auto& typed_storage = getEigenSolveStorage<SolverType>();
    const Eigen::Index n = eigen.n_dofs;
    TBOX_ASSERT(n > 0);

    {
        ConstPetscVecArrayMap x_array(x, n);
        PetscVecArrayMap y_array(y, n);
        const auto x_map = x_array.getMap();
        auto y_map = y_array.getMap();
        y_map.setZero();
        for (std::size_t subdomain_num = 0; subdomain_num < eigen.common_subdomains.size(); ++subdomain_num)
        {
            auto& cache = eigen.common_subdomains[subdomain_num];
            auto& solver = typed_storage.solvers[subdomain_num];
            std::size_t rhs_idx = 0;
            for (const int dof : cache.overlap_dofs)
            {
                cache.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = x_map[dof];
            }
            cache.delta_workspace = solve_eigen_system(solver, cache.rhs_workspace);
            std::size_t update_idx = 0;
            for (const int dof : cache.update_dofs)
            {
                y_map[dof] +=
                    cache.delta_workspace[static_cast<Eigen::Index>(cache.update_local_positions[update_idx++])];
            }
        }
    }
    postprocessShellResult(y);
}

template <class SolverType>
void
PETScLevelSolver::applyMultiplicativeEigenSolveImpl(Vec x, Vec y)
{
    TBOX_ASSERT(IBTK_MPI::getNodes() == 1);
    auto& shell = d_shell_data;
    auto& eigen = *shell.eigen_data;
    auto& typed_storage = getEigenSolveStorage<SolverType>();
    const Eigen::Index n = eigen.n_dofs;
    TBOX_ASSERT(n > 0);

    {
        ConstPetscVecArrayMap x_array(x, n);
        PetscVecArrayMap y_array(y, n);
        const auto x_map = x_array.getMap();
        auto y_map = y_array.getMap();
        Eigen::VectorXd residual(n);
        y_map.setZero();
        residual = x_map;
        const std::size_t n_subdomains = static_cast<std::size_t>(shell.n_local_subdomains);
        for (std::size_t subdomain_num = 0; subdomain_num + 1 < n_subdomains; ++subdomain_num)
        {
            auto& cache = eigen.common_subdomains[subdomain_num];
            auto& solver = typed_storage.solvers[subdomain_num];
            std::size_t rhs_idx = 0;
            for (const int dof : cache.overlap_dofs)
            {
                cache.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = residual[dof];
            }
            cache.delta_workspace = solve_eigen_system(solver, cache.rhs_workspace);

            std::size_t update_idx = 0;
            for (const int dof : cache.update_dofs)
            {
                y_map[dof] +=
                    cache.delta_workspace[static_cast<Eigen::Index>(cache.update_local_positions[update_idx++])];
            }
            if (cache.active_residual_update_mat.rows() > 0)
            {
                std::size_t residual_input_idx = 0;
                for (const int local_pos : cache.update_local_positions)
                {
                    cache.residual_input_workspace[static_cast<Eigen::Index>(residual_input_idx++)] =
                        cache.delta_workspace[static_cast<Eigen::Index>(local_pos)];
                }
                cache.residual_delta_workspace.noalias() =
                    cache.active_residual_update_mat * cache.residual_input_workspace;
                std::size_t row_idx = 0;
                for (const int row : cache.active_residual_update_rows)
                {
                    residual[row] -= cache.residual_delta_workspace[static_cast<Eigen::Index>(row_idx++)];
                }
            }
        }
        if (n_subdomains > 0)
        {
            auto& cache = eigen.common_subdomains.back();
            auto& solver = typed_storage.solvers.back();
            std::size_t rhs_idx = 0;
            for (const int dof : cache.overlap_dofs)
            {
                cache.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = residual[dof];
            }
            cache.delta_workspace = solve_eigen_system(solver, cache.rhs_workspace);

            std::size_t update_idx = 0;
            for (const int dof : cache.update_dofs)
            {
                y_map[dof] +=
                    cache.delta_workspace[static_cast<Eigen::Index>(cache.update_local_positions[update_idx++])];
            }
        }
    }
    postprocessShellResult(y);
}

const std::vector<int>&
PETScLevelSolver::getEigenOverlapDofs(const std::size_t subdomain_num) const
{
    TBOX_ASSERT(d_shell_data.eigen_data);
    TBOX_ASSERT(subdomain_num < d_shell_data.eigen_data->common_subdomains.size());
    return d_shell_data.eigen_data->common_subdomains[subdomain_num].overlap_dofs;
} // getEigenOverlapDofs

const std::vector<int>&
PETScLevelSolver::getEigenUpdateDofs(const std::size_t subdomain_num) const
{
    TBOX_ASSERT(d_shell_data.eigen_data);
    TBOX_ASSERT(subdomain_num < d_shell_data.eigen_data->common_subdomains.size());
    return d_shell_data.eigen_data->common_subdomains[subdomain_num].update_dofs;
} // getEigenUpdateDofs

const std::vector<int>&
PETScLevelSolver::getEigenUpdateLocalPositions(const std::size_t subdomain_num) const
{
    TBOX_ASSERT(d_shell_data.eigen_data);
    TBOX_ASSERT(subdomain_num < d_shell_data.eigen_data->common_subdomains.size());
    return d_shell_data.eigen_data->common_subdomains[subdomain_num].update_local_positions;
} // getEigenUpdateLocalPositions

const std::vector<int>&
PETScLevelSolver::getEigenActiveResidualUpdateRows(const std::size_t subdomain_num) const
{
    TBOX_ASSERT(d_shell_data.eigen_data);
    TBOX_ASSERT(subdomain_num < d_shell_data.eigen_data->common_subdomains.size());
    return d_shell_data.eigen_data->common_subdomains[subdomain_num].active_residual_update_rows;
} // getEigenActiveResidualUpdateRows

const Eigen::SparseMatrix<double, Eigen::RowMajor>&
PETScLevelSolver::getEigenActiveResidualUpdateMat(const std::size_t subdomain_num) const
{
    TBOX_ASSERT(d_shell_data.eigen_data);
    TBOX_ASSERT(subdomain_num < d_shell_data.eigen_data->common_subdomains.size());
    return d_shell_data.eigen_data->common_subdomains[subdomain_num].active_residual_update_mat;
} // getEigenActiveResidualUpdateMat

void
PETScLevelSolver::initializeBuiltinEigenSubdomainSolver(const EigenSubdomainSolverType solver_type,
                                                        const Eigen::MatrixXd& local_operator,
                                                        const std::size_t subdomain_num)
{
    dispatchEigenSolverType(solver_type,
                            [&](auto solver_tag)
                            {
                                using SolverType = typename decltype(solver_tag)::type;
                                initialize_eigen_solver(getEigenSolveStorage<SolverType>().solvers[subdomain_num],
                                                        local_operator,
                                                        d_eigen_subdomain_solver_threshold);
                            });
    return;
} // initializeBuiltinEigenSubdomainSolver

Eigen::VectorXd
PETScLevelSolver::solveBuiltinEigenSubdomainSystem(const EigenSubdomainSolverType solver_type,
                                                   const Eigen::VectorXd& rhs,
                                                   const std::size_t subdomain_num) const
{
    Eigen::VectorXd solution;
    dispatchEigenSolverType(solver_type,
                            [&](auto solver_tag)
                            {
                                using SolverType = typename decltype(solver_tag)::type;
                                solution =
                                    solve_eigen_system(getEigenSolveStorage<SolverType>().solvers[subdomain_num], rhs);
                            });
    return solution;
} // solveBuiltinEigenSubdomainSystem

void
PETScLevelSolver::cacheASMSubdomains(const std::vector<std::set<int>>& subdomain_dofs,
                                     const std::vector<std::set<int>>& nonoverlap_subdomain_dofs)
{
    d_subdomain_dofs.resize(subdomain_dofs.size());
    for (std::size_t k = 0; k < subdomain_dofs.size(); ++k)
    {
        d_subdomain_dofs[k].assign(subdomain_dofs[k].begin(), subdomain_dofs[k].end());
    }

    d_nonoverlap_subdomain_dofs.resize(nonoverlap_subdomain_dofs.size());
    for (std::size_t k = 0; k < nonoverlap_subdomain_dofs.size(); ++k)
    {
        d_nonoverlap_subdomain_dofs[k].assign(nonoverlap_subdomain_dofs[k].begin(), nonoverlap_subdomain_dofs[k].end());
    }
    return;
} // cacheASMSubdomains

void
PETScLevelSolver::cacheGeneratedASMSubdomains()
{
    std::vector<std::set<int>> subdomain_dofs, nonoverlap_subdomain_dofs;
    generateASMSubdomains(subdomain_dofs, nonoverlap_subdomain_dofs);
    cacheASMSubdomains(subdomain_dofs, nonoverlap_subdomain_dofs);
    return;
} // cacheGeneratedASMSubdomains

void
PETScLevelSolver::configureASMPreconditioner(PC ksp_pc)
{
    int ierr;
    cacheGeneratedASMSubdomains();
    build_petsc_subdomain_index_sets(
        d_petsc_subdomain_is, d_petsc_nonoverlap_subdomain_is, d_subdomain_dofs, d_nonoverlap_subdomain_dofs);
    const int num_subdomains = static_cast<int>(d_petsc_subdomain_is.size());
    if (num_subdomains == 0)
    {
        IS index_set = nullptr;
        ierr = ISCreateGeneral(PETSC_COMM_SELF, 0, nullptr, PETSC_OWN_POINTER, &index_set);
        IBTK_CHKERRQ(ierr);
        ierr = PCASMSetLocalSubdomains(ksp_pc, 1, &index_set, &index_set);
        IBTK_CHKERRQ(ierr);
        ierr = ISDestroy(&index_set);
        IBTK_CHKERRQ(ierr);
    }
    else
    {
        ierr = PCASMSetLocalSubdomains(
            ksp_pc, num_subdomains, d_petsc_subdomain_is.data(), d_petsc_nonoverlap_subdomain_is.data());
        IBTK_CHKERRQ(ierr);
    }
    return;
} // configureASMPreconditioner

void
PETScLevelSolver::configureFieldSplitPreconditioner(PC ksp_pc)
{
    std::vector<std::set<int>> field_is;
    std::vector<std::string> field_name;
    generateFieldSplitSubdomains(field_name, field_is);
    d_field_name = field_name;
    const int n_fields = static_cast<int>(field_is.size());

    destroy_index_sets(d_field_is);
    d_field_is.resize(n_fields);
    int ierr;
    for (int field_num = 0; field_num < n_fields; ++field_num)
    {
        PetscInt* field_dof_arr = nullptr;
        const int n_field_dofs = static_cast<int>(field_is[field_num].size());
        ierr = PetscMalloc1(n_field_dofs, &field_dof_arr);
        IBTK_CHKERRQ(ierr);
        std::copy(field_is[field_num].begin(), field_is[field_num].end(), field_dof_arr);
        ierr =
            ISCreateGeneral(PETSC_COMM_WORLD, n_field_dofs, field_dof_arr, PETSC_OWN_POINTER, &d_field_is[field_num]);
        IBTK_CHKERRQ(ierr);
        ierr = PCFieldSplitSetIS(ksp_pc, d_field_name[field_num].c_str(), d_field_is[field_num]);
        IBTK_CHKERRQ(ierr);
    }
    return;
} // configureFieldSplitPreconditioner

void
PETScLevelSolver::initializePetscShellData()
{
    auto& shell = d_shell_data;
    shell.petsc_data = std::make_unique<PetscShellSmootherData>();
    auto& petsc = *shell.petsc_data;
    const bool use_restrict_partition = d_shell_smoother_partition == ShellSmootherPartition::RESTRICT;
    const bool use_multiplicative = d_shell_smoother_composition == ShellSmootherComposition::MULTIPLICATIVE;
    int ierr;
    if (use_multiplicative)
    {
        ierr = VecDuplicate(d_petsc_x, &petsc.shell_r);
        IBTK_CHKERRQ(ierr);
    }
    petsc.prolongation_insert_mode = use_restrict_partition ? INSERT_VALUES : ADD_VALUES;

    ierr = MatCreateSubMatrices(d_petsc_mat,
                                shell.n_local_subdomains,
                                d_petsc_subdomain_is.data(),
                                d_petsc_subdomain_is.data(),
                                MAT_INITIAL_MATRIX,
                                &petsc.sub_mat);
    IBTK_CHKERRQ(ierr);

    petsc.local_overlap_is.resize(shell.n_local_subdomains);
    petsc.restriction.resize(shell.n_local_subdomains);
    petsc.prolongation.resize(shell.n_local_subdomains);
    petsc.sub_x.resize(shell.n_local_subdomains);
    petsc.sub_y.resize(shell.n_local_subdomains);
    if (use_restrict_partition) petsc.local_nonoverlap_is.resize(shell.n_local_subdomains);
    if (use_multiplicative)
    {
        petsc.active_update_local_positions.resize(shell.n_local_subdomains);
        petsc.active_residual_update_x.resize(shell.n_local_subdomains);
        petsc.active_residual_update_y.resize(shell.n_local_subdomains);
    }
#if !defined(NDEBUG)
    std::set<int> idxs;
#endif
    for (int subdomain_num = 0; subdomain_num < shell.n_local_subdomains; ++subdomain_num)
    {
        const auto& overlap_dofs = d_subdomain_dofs[static_cast<std::size_t>(subdomain_num)];
        const int overlap_is_size = static_cast<int>(overlap_dofs.size());
        PetscInt* overlap_indices = nullptr;
        ierr = PetscMalloc1(overlap_is_size, &overlap_indices);
        IBTK_CHKERRQ(ierr);
        for (int overlap_local_idx = 0; overlap_local_idx < overlap_is_size; ++overlap_local_idx)
        {
            overlap_indices[overlap_local_idx] = overlap_local_idx;
#if !defined(NDEBUG)
            TBOX_ASSERT(idxs.find(overlap_dofs[static_cast<std::size_t>(overlap_local_idx)]) == idxs.end());
            idxs.insert(overlap_dofs[static_cast<std::size_t>(overlap_local_idx)]);
#endif
        }

        const auto& nonoverlap_dofs = d_nonoverlap_subdomain_dofs[static_cast<std::size_t>(subdomain_num)];
        const int nonoverlap_is_size = static_cast<int>(nonoverlap_dofs.size());
        PetscInt* nonoverlap_indices = nullptr;
        if (use_restrict_partition)
        {
            ierr = PetscMalloc1(nonoverlap_is_size, &nonoverlap_indices);
            IBTK_CHKERRQ(ierr);
            int nonoverlap_local_idx = 0;
            for (int overlap_local_idx = 0; overlap_local_idx < overlap_is_size; ++overlap_local_idx)
            {
                if (nonoverlap_local_idx < nonoverlap_is_size &&
                    overlap_dofs[static_cast<std::size_t>(overlap_local_idx)] ==
                        nonoverlap_dofs[static_cast<std::size_t>(nonoverlap_local_idx)])
                {
                    nonoverlap_indices[nonoverlap_local_idx] = overlap_local_idx;
                    ++nonoverlap_local_idx;
                }
            }
            TBOX_ASSERT(nonoverlap_local_idx == nonoverlap_is_size);
        }

        if (use_multiplicative)
        {
            auto& active_update_local_positions = petsc.active_update_local_positions[subdomain_num];
            switch (d_shell_smoother_partition)
            {
            case ShellSmootherPartition::BASIC:
                active_update_local_positions.resize(static_cast<std::size_t>(overlap_is_size));
                for (int local_pos = 0; local_pos < overlap_is_size; ++local_pos)
                {
                    active_update_local_positions[static_cast<std::size_t>(local_pos)] = local_pos;
                }
                break;
            case ShellSmootherPartition::RESTRICT:
                active_update_local_positions.resize(static_cast<std::size_t>(nonoverlap_is_size));
                for (int local_pos = 0; local_pos < nonoverlap_is_size; ++local_pos)
                {
                    active_update_local_positions[static_cast<std::size_t>(local_pos)] =
                        static_cast<int>(nonoverlap_indices[local_pos]);
                }
                break;
            }
        }
        ierr = MatCreateVecs(petsc.sub_mat[subdomain_num], &petsc.sub_x[subdomain_num], &petsc.sub_y[subdomain_num]);
        IBTK_CHKERRQ(ierr);

        ierr = ISCreateGeneral(PETSC_COMM_SELF,
                               overlap_is_size,
                               overlap_indices,
                               PETSC_OWN_POINTER,
                               &petsc.local_overlap_is[subdomain_num]);
        IBTK_CHKERRQ(ierr);
        if (use_restrict_partition)
        {
            ierr = ISCreateGeneral(PETSC_COMM_SELF,
                                   nonoverlap_is_size,
                                   nonoverlap_indices,
                                   PETSC_OWN_POINTER,
                                   &petsc.local_nonoverlap_is[subdomain_num]);
            IBTK_CHKERRQ(ierr);
        }

        ierr = VecScatterCreate(d_petsc_x,
                                d_petsc_subdomain_is[subdomain_num],
                                petsc.sub_x[subdomain_num],
                                petsc.local_overlap_is[subdomain_num],
                                &petsc.restriction[subdomain_num]);
        IBTK_CHKERRQ(ierr);
        switch (d_shell_smoother_partition)
        {
        case ShellSmootherPartition::BASIC:
            ierr = VecScatterCreate(petsc.sub_y[subdomain_num],
                                    petsc.local_overlap_is[subdomain_num],
                                    d_petsc_b,
                                    d_petsc_subdomain_is[subdomain_num],
                                    &petsc.prolongation[subdomain_num]);
            IBTK_CHKERRQ(ierr);
            break;
        case ShellSmootherPartition::RESTRICT:
            ierr = VecScatterCreate(petsc.sub_y[subdomain_num],
                                    petsc.local_nonoverlap_is[subdomain_num],
                                    d_petsc_b,
                                    d_petsc_nonoverlap_subdomain_is[subdomain_num],
                                    &petsc.prolongation[subdomain_num]);
            IBTK_CHKERRQ(ierr);
            break;
        }
    }
#if !defined(NDEBUG)
    int n_local_dofs = 0;
    ierr = VecGetLocalSize(d_petsc_x, &n_local_dofs);
    IBTK_CHKERRQ(ierr);
    TBOX_ASSERT(n_local_dofs == static_cast<int>(idxs.size()));
#endif

    if (use_multiplicative)
    {
        PetscInt row_begin = 0, row_end = 0;
        ierr = MatGetOwnershipRange(d_petsc_mat, &row_begin, &row_end);
        IBTK_CHKERRQ(ierr);
        const PetscInt n_owned_rows = row_end - row_begin;
        std::vector<IS> active_update_global_is(shell.n_local_subdomains);
        std::vector<IS> owned_residual_update_rows_is(shell.n_local_subdomains, nullptr);
        ierr = ISCreateStride(PETSC_COMM_SELF, n_owned_rows, row_begin, 1, &petsc.owned_residual_update_rows_is);
        IBTK_CHKERRQ(ierr);
        for (int subdomain_num = 0; subdomain_num < shell.n_local_subdomains; ++subdomain_num)
        {
            switch (d_shell_smoother_partition)
            {
            case ShellSmootherPartition::BASIC:
                active_update_global_is[subdomain_num] = d_petsc_subdomain_is[subdomain_num];
                break;
            case ShellSmootherPartition::RESTRICT:
                active_update_global_is[subdomain_num] = d_petsc_nonoverlap_subdomain_is[subdomain_num];
                break;
            }
            owned_residual_update_rows_is[subdomain_num] = petsc.owned_residual_update_rows_is;
        }
        ierr = MatCreateSubMatrices(d_petsc_mat,
                                    shell.n_local_subdomains,
                                    owned_residual_update_rows_is.data(),
                                    active_update_global_is.data(),
                                    MAT_INITIAL_MATRIX,
                                    &petsc.active_residual_update_mat);
        IBTK_CHKERRQ(ierr);

        for (int subdomain_num = 0; subdomain_num < shell.n_local_subdomains; ++subdomain_num)
        {
            ierr = MatCreateVecs(petsc.active_residual_update_mat[subdomain_num],
                                 &petsc.active_residual_update_x[subdomain_num],
                                 &petsc.active_residual_update_y[subdomain_num]);
            IBTK_CHKERRQ(ierr);
        }
    }

    petsc.sub_ksp.resize(shell.n_local_subdomains);
    for (int subdomain_num = 0; subdomain_num < shell.n_local_subdomains; ++subdomain_num)
    {
        KSP& sub_ksp = petsc.sub_ksp[subdomain_num];
        Mat& sub_mat = petsc.sub_mat[subdomain_num];
        ierr = KSPCreate(PETSC_COMM_SELF, &sub_ksp);
        IBTK_CHKERRQ(ierr);
        const std::string sub_prefix = d_options_prefix + "_sub";
        ierr = KSPSetOptionsPrefix(sub_ksp, sub_prefix.c_str());
        IBTK_CHKERRQ(ierr);
        ierr = KSPSetOperators(sub_ksp, sub_mat, sub_mat);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSetReusePreconditioner(sub_ksp, PETSC_TRUE);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSetType(sub_ksp, KSPPREONLY);
        IBTK_CHKERRQ(ierr);
        PC sub_pc = nullptr;
        ierr = KSPGetPC(sub_ksp, &sub_pc);
        IBTK_CHKERRQ(ierr);
        ierr = PCSetType(sub_pc, PCSVD);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSetFromOptions(sub_ksp);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSetInitialGuessNonzero(sub_ksp, PETSC_FALSE);
        IBTK_CHKERRQ(ierr);
    }
    return;
} // initializePetscShellData

void
PETScLevelSolver::initializeEigenShellData()
{
    switch (d_shell_smoother_backend)
    {
    case ShellSmootherBackend::EIGEN:
        if (usingCustomEigenSolveSubdomainSolver())
        {
            initializeCustomEigenShellData();
            break;
        }
        dispatchEigenSolverType(getEigenSolveSolverType(),
                                [&](auto solver_tag)
                                {
                                    using SolverType = typename decltype(solver_tag)::type;
                                    initializeBuiltinEigenSolveStorage<SolverType>();
                                });
        break;
    case ShellSmootherBackend::EIGEN_PSEUDOINVERSE:
        if (usingCustomEigenPseudoinverseSubdomainSolver())
        {
            TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::initializeEigenShellData()\n"
                                     << "Custom Eigen matrix-based shell data must be implemented by a subclass."
                                     << std::endl);
        }
        switch (getEigenPseudoinverseSolverType())
        {
        case EigenSubdomainSolverType::CUSTOM:
            TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::initializeEigenShellData()\n"
                                     << "CUSTOM is not a built-in Eigen matrix-based subdomain solver type."
                                     << std::endl);
            break;
        case EigenSubdomainSolverType::LLT:
        case EigenSubdomainSolverType::LDLT:
        case EigenSubdomainSolverType::PARTIAL_PIV_LU:
        case EigenSubdomainSolverType::FULL_PIV_LU:
        case EigenSubdomainSolverType::HOUSEHOLDER_QR:
        case EigenSubdomainSolverType::COL_PIV_HOUSEHOLDER_QR:
        case EigenSubdomainSolverType::COMPLETE_ORTHOGONAL_DECOMPOSITION:
        case EigenSubdomainSolverType::FULL_PIV_HOUSEHOLDER_QR:
        case EigenSubdomainSolverType::JACOBI_SVD:
        case EigenSubdomainSolverType::BDC_SVD:
            initializeEigenShellDataWithLocalOperatorHook(
                [this](const Eigen::MatrixXd& local_operator, const std::size_t subdomain_num)
                {
                    auto& cache = d_shell_data.eigen_data->common_subdomains[subdomain_num];
                    cache.local_pseudoinverse = buildEigenSubdomainPseudoinverse(local_operator);
                });
            break;
        }
        break;
    case ShellSmootherBackend::PETSC:
        TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::initializeEigenShellData()\n"
                                 << "PETSc shell backend should not initialize Eigen shell data." << std::endl);
        break;
    }
    return;
} // initializeEigenShellData

void
PETScLevelSolver::initializeEigenSubdomainSolver(const Eigen::MatrixXd& local_operator, std::size_t subdomain_num)
{
    switch (d_shell_smoother_backend)
    {
    case ShellSmootherBackend::EIGEN:
        if (usingCustomEigenSolveSubdomainSolver())
        {
            TBOX_ERROR(d_object_name << " " << d_options_prefix
                                     << " PETScLevelSolver::initializeEigenSubdomainSolver()\n"
                                     << "Custom Eigen solve-based subdomain solvers must be implemented by a "
                                     << "subclass." << std::endl);
        }
        initializeBuiltinEigenSubdomainSolver(getEigenSolveSolverType(), local_operator, subdomain_num);
        break;
    case ShellSmootherBackend::EIGEN_PSEUDOINVERSE:
        if (usingCustomEigenPseudoinverseSubdomainSolver())
        {
            TBOX_ERROR(d_object_name << " " << d_options_prefix
                                     << " PETScLevelSolver::initializeEigenSubdomainSolver()\n"
                                     << "Custom Eigen matrix-based subdomain solvers must be implemented by a "
                                     << "subclass." << std::endl);
        }
        d_shell_data.eigen_data->common_subdomains[subdomain_num].local_pseudoinverse =
            buildEigenSubdomainPseudoinverse(local_operator);
        break;
    case ShellSmootherBackend::PETSC:
        TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::initializeEigenSubdomainSolver()\n"
                                 << "PETSc shell backend should not initialize Eigen subdomain solvers." << std::endl);
        break;
    }
    return;
} // initializeEigenSubdomainSolver

Eigen::VectorXd
PETScLevelSolver::solveEigenSubdomainSystem(const Eigen::VectorXd& rhs, std::size_t subdomain_num) const
{
    switch (d_shell_smoother_backend)
    {
    case ShellSmootherBackend::EIGEN:
        if (usingCustomEigenSolveSubdomainSolver())
        {
            TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::solveEigenSubdomainSystem()\n"
                                     << "Custom Eigen solve-based subdomain solvers must be implemented by a "
                                     << "subclass." << std::endl);
        }
        return solveBuiltinEigenSubdomainSystem(getEigenSolveSolverType(), rhs, subdomain_num);
        break;
    case ShellSmootherBackend::EIGEN_PSEUDOINVERSE:
        if (usingCustomEigenPseudoinverseSubdomainSolver())
        {
            TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::solveEigenSubdomainSystem()\n"
                                     << "Custom Eigen matrix-based subdomain solvers must be implemented by a "
                                     << "subclass." << std::endl);
        }
        return d_shell_data.eigen_data->common_subdomains[subdomain_num].local_pseudoinverse * rhs;
    case ShellSmootherBackend::PETSC:
        TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::solveEigenSubdomainSystem()\n"
                                 << "PETSc shell backend should not use the Eigen solve path." << std::endl);
        break;
    }

    return Eigen::VectorXd();
} // solveEigenSubdomainSystem

void
PETScLevelSolver::initializeCustomEigenShellData()
{
    TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::initializeCustomEigenShellData()\n"
                             << "Custom Eigen shell data must be implemented by a subclass." << std::endl);
} // initializeCustomEigenShellData

void
PETScLevelSolver::applyAdditiveCustomEigen(Vec /*x*/, Vec /*y*/)
{
    TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::applyAdditiveCustomEigen()\n"
                             << "Custom additive Eigen shell application must be implemented by a subclass."
                             << std::endl);
} // applyAdditiveCustomEigen

void
PETScLevelSolver::applyMultiplicativeCustomEigen(Vec /*x*/, Vec /*y*/)
{
    TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::applyMultiplicativeCustomEigen()\n"
                             << "Custom multiplicative Eigen shell application must be implemented by a subclass."
                             << std::endl);
} // applyMultiplicativeCustomEigen

Eigen::MatrixXd
PETScLevelSolver::buildEigenSubdomainPseudoinverse(const Eigen::MatrixXd& local_operator) const
{
    switch (getEigenPseudoinverseSolverType())
    {
    case EigenSubdomainSolverType::CUSTOM:
        TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::buildEigenSubdomainPseudoinverse()\n"
                                 << "Custom Eigen matrix-based subdomain solvers must be implemented by a subclass."
                                 << std::endl);
        break;
    case EigenSubdomainSolverType::LLT:
        return build_llt_solve_matrix(local_operator);
    case EigenSubdomainSolverType::LDLT:
        return build_ldlt_solve_matrix(local_operator);
    case EigenSubdomainSolverType::PARTIAL_PIV_LU:
        return build_partial_piv_lu_solve_matrix(local_operator);
    case EigenSubdomainSolverType::FULL_PIV_LU:
        return build_full_piv_lu_solve_matrix(local_operator, d_eigen_subdomain_pseudoinverse_threshold);
    case EigenSubdomainSolverType::HOUSEHOLDER_QR:
        return build_householder_qr_solve_matrix(local_operator);
    case EigenSubdomainSolverType::COL_PIV_HOUSEHOLDER_QR:
        return build_qr_solve_matrix<Eigen::ColPivHouseholderQR<Eigen::MatrixXd>>(
            local_operator, d_eigen_subdomain_pseudoinverse_threshold);
    case EigenSubdomainSolverType::COMPLETE_ORTHOGONAL_DECOMPOSITION:
        return build_complete_orthogonal_decomposition_pseudoinverse(local_operator,
                                                                     d_eigen_subdomain_pseudoinverse_threshold);
    case EigenSubdomainSolverType::FULL_PIV_HOUSEHOLDER_QR:
        return build_qr_solve_matrix<Eigen::FullPivHouseholderQR<Eigen::MatrixXd>>(
            local_operator, d_eigen_subdomain_pseudoinverse_threshold);
    case EigenSubdomainSolverType::JACOBI_SVD:
        return build_svd_pseudoinverse<Eigen::JacobiSVD<Eigen::MatrixXd>>(local_operator,
                                                                          d_eigen_subdomain_pseudoinverse_threshold);
    case EigenSubdomainSolverType::BDC_SVD:
        return build_svd_pseudoinverse<Eigen::BDCSVD<Eigen::MatrixXd>>(local_operator,
                                                                       d_eigen_subdomain_pseudoinverse_threshold);
    }

    TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::buildEigenSubdomainPseudoinverse()\n"
                             << "Unknown Eigen subdomain pseudoinverse type." << std::endl);
    return Eigen::MatrixXd();
} // buildEigenSubdomainPseudoinverse

void
PETScLevelSolver::deallocatePetscShellData()
{
    auto& shell = d_shell_data;
    if (!shell.petsc_data) return;

    auto& petsc = *shell.petsc_data;
    int ierr;
    for (int subdomain_num = 0; subdomain_num < shell.n_local_subdomains; ++subdomain_num)
    {
        ierr = KSPDestroy(&petsc.sub_ksp[subdomain_num]);
        IBTK_CHKERRQ(ierr);
    }
    for (int subdomain_num = 0; subdomain_num < shell.n_local_subdomains; ++subdomain_num)
    {
        ierr = ISDestroy(&petsc.local_overlap_is[subdomain_num]);
        IBTK_CHKERRQ(ierr);
        if (subdomain_num < static_cast<int>(petsc.local_nonoverlap_is.size()))
        {
            ierr = ISDestroy(&petsc.local_nonoverlap_is[subdomain_num]);
            IBTK_CHKERRQ(ierr);
        }
        ierr = VecScatterDestroy(&petsc.prolongation[subdomain_num]);
        IBTK_CHKERRQ(ierr);
        ierr = VecScatterDestroy(&petsc.restriction[subdomain_num]);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&petsc.sub_x[subdomain_num]);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&petsc.sub_y[subdomain_num]);
        IBTK_CHKERRQ(ierr);
        if (subdomain_num < static_cast<int>(petsc.active_residual_update_x.size()))
        {
            ierr = VecDestroy(&petsc.active_residual_update_x[subdomain_num]);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&petsc.active_residual_update_y[subdomain_num]);
            IBTK_CHKERRQ(ierr);
        }
    }
    ierr = MatDestroyMatrices(shell.n_local_subdomains, &petsc.sub_mat);
    IBTK_CHKERRQ(ierr);
    if (petsc.active_residual_update_mat)
    {
        ierr = MatDestroyMatrices(shell.n_local_subdomains, &petsc.active_residual_update_mat);
        IBTK_CHKERRQ(ierr);
    }
    if (petsc.owned_residual_update_rows_is)
    {
        ierr = ISDestroy(&petsc.owned_residual_update_rows_is);
        IBTK_CHKERRQ(ierr);
    }
    if (petsc.shell_r)
    {
        ierr = VecDestroy(&petsc.shell_r);
        IBTK_CHKERRQ(ierr);
    }
    shell.petsc_data.reset();
    return;
} // deallocatePetscShellData

void
PETScLevelSolver::deallocateEigenShellData()
{
    auto& shell = d_shell_data;
    if (!shell.eigen_data) return;

    shell.eigen_data->n_dofs = 0;
    shell.eigen_data->common_subdomains.clear();
    shell.eigen_data->solve_storage.reset();
    shell.eigen_data.reset();
    return;
} // deallocateEigenShellData

void
PETScLevelSolver::deallocateShellData()
{
    deallocatePetscShellData();
    deallocateEigenShellData();
    d_shell_data.n_local_subdomains = 0;
    d_subdomain_dofs.clear();
    d_nonoverlap_subdomain_dofs.clear();
    destroy_index_sets(d_petsc_nonoverlap_subdomain_is);
    destroy_index_sets(d_petsc_subdomain_is);
    return;
} // deallocateShellData

void
PETScLevelSolver::configureShellApply(PC ksp_pc)
{
    int ierr;
    ierr = PCSetType(ksp_pc, PCSHELL);
    IBTK_CHKERRQ(ierr);
    ierr = PCShellSetContext(ksp_pc, static_cast<void*>(this));
    IBTK_CHKERRQ(ierr);

    PetscErrorCode (*apply_op)(PC, Vec, Vec) = nullptr;
    const char* pc_name_suffix = nullptr;
    switch (d_shell_smoother_backend)
    {
    case ShellSmootherBackend::PETSC:
        switch (d_shell_smoother_composition)
        {
        case ShellSmootherComposition::ADDITIVE:
            apply_op = PETScLevelSolver::PCApply_AdditivePetsc;
            pc_name_suffix = "PC_AdditivePetsc";
            break;
        case ShellSmootherComposition::MULTIPLICATIVE:
            apply_op = PETScLevelSolver::PCApply_MultiplicativePetsc;
            pc_name_suffix = "PC_MultiplicativePetsc";
            break;
        }
        break;
    case ShellSmootherBackend::EIGEN:
    case ShellSmootherBackend::EIGEN_PSEUDOINVERSE:
        switch (d_shell_smoother_composition)
        {
        case ShellSmootherComposition::ADDITIVE:
            apply_op = PETScLevelSolver::PCApply_AdditiveEigen;
            switch (d_shell_smoother_backend)
            {
            case ShellSmootherBackend::EIGEN:
                pc_name_suffix = "PC_AdditiveEigen";
                break;
            case ShellSmootherBackend::EIGEN_PSEUDOINVERSE:
                pc_name_suffix = "PC_AdditiveEigenPseudoinverse";
                break;
            case ShellSmootherBackend::PETSC:
                break;
            }
            break;
        case ShellSmootherComposition::MULTIPLICATIVE:
            apply_op = PETScLevelSolver::PCApply_MultiplicativeEigen;
            switch (d_shell_smoother_backend)
            {
            case ShellSmootherBackend::EIGEN:
                pc_name_suffix = "PC_MultiplicativeEigen";
                break;
            case ShellSmootherBackend::EIGEN_PSEUDOINVERSE:
                pc_name_suffix = "PC_MultiplicativeEigenPseudoinverse";
                break;
            case ShellSmootherBackend::PETSC:
                break;
            }
            break;
        }
        break;
    }
    if (!apply_op || !pc_name_suffix)
    {
        TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::configureShellApply()\n"
                                 << "Unknown PCSHELL specified. Supported values are additive, additive-basic, "
                                    "additive-restrict, additive-eigen, additive-eigen-basic, "
                                    "additive-eigen-restrict, "
                                    "additive-eigen-pseudoinverse, "
                                    "additive-eigen-pseudoinverse-basic, additive-eigen-pseudoinverse-restrict, "
                                    "multiplicative, multiplicative-basic, multiplicative-restrict, "
                                    "multiplicative-eigen, multiplicative-eigen-basic, "
                                    "multiplicative-eigen-restrict, "
                                    "multiplicative-eigen-pseudoinverse, "
                                    "multiplicative-eigen-pseudoinverse-basic, and "
                                    "multiplicative-eigen-pseudoinverse-restrict."
                                 << std::endl);
    }
    ierr = PCShellSetApply(ksp_pc, apply_op);
    IBTK_CHKERRQ(ierr);
    const std::string pc_name = d_options_prefix + pc_name_suffix;
    ierr = PCShellSetName(ksp_pc, pc_name.c_str());
    IBTK_CHKERRQ(ierr);
    return;
} // configureShellApply

void
PETScLevelSolver::configureShellPreconditioner(PC ksp_pc)
{
    auto& shell = d_shell_data;
    cacheGeneratedASMSubdomains();
    shell.n_local_subdomains = static_cast<int>(d_subdomain_dofs.size());

    switch (d_shell_smoother_backend)
    {
    case ShellSmootherBackend::PETSC:
        build_petsc_subdomain_index_sets(
            d_petsc_subdomain_is, d_petsc_nonoverlap_subdomain_is, d_subdomain_dofs, d_nonoverlap_subdomain_dofs);
        initializePetscShellData();
        break;
    case ShellSmootherBackend::EIGEN:
    case ShellSmootherBackend::EIGEN_PSEUDOINVERSE:
        initializeEigenShellData();
        break;
    }

    configureShellApply(ksp_pc);
    return;
} // configureShellPreconditioner

void
PETScLevelSolver::generateASMSubdomains(std::vector<std::set<int>>& /*subdomain_dofs*/,
                                        std::vector<std::set<int>>& /*nonoverlap_subdomain_dofs*/)
{
    TBOX_ERROR("PETScLevelSolver::generateASMSubdomains(): Subclasses need to generate ASM subdomains. \n");

    return;
} // generateASMSubdomains

void
PETScLevelSolver::generateFieldSplitSubdomains(std::vector<std::string>& /*field_names*/,
                                               std::vector<std::set<int>>& /*field_is*/)
{
    TBOX_ERROR(
        "PETScLevelSolver::generateFieldSplitSubdomains(): Subclasses need to generate FieldSplit subdomains. \n");

    return;
} // generateFieldSplitSubdomains

void
PETScLevelSolver::postprocessShellResult(Vec& /*y*/)
{
    return;
} // postprocessShellResult

void
PETScLevelSolver::beginAccumulateCorrectionPetsc(const int subdomain_num, Vec sub_y, Vec y)
{
    auto& petsc = *d_shell_data.petsc_data;
    int ierr;
    ierr = VecScatterBegin(
        petsc.prolongation[subdomain_num], sub_y, y, petsc.prolongation_insert_mode, SCATTER_FORWARD_LOCAL);
    IBTK_CHKERRQ(ierr);
    return;
} // beginAccumulateCorrectionPetsc

void
PETScLevelSolver::endAccumulateCorrectionPetsc(const int subdomain_num, Vec sub_y, Vec y)
{
    auto& petsc = *d_shell_data.petsc_data;
    int ierr;
    ierr = VecScatterEnd(
        petsc.prolongation[subdomain_num], sub_y, y, petsc.prolongation_insert_mode, SCATTER_FORWARD_LOCAL);
    IBTK_CHKERRQ(ierr);
    return;
} // endAccumulateCorrectionPetsc

void
PETScLevelSolver::accumulateCorrectionPetsc(const int subdomain_num, Vec sub_y, Vec y)
{
    beginAccumulateCorrectionPetsc(subdomain_num, sub_y, y);
    endAccumulateCorrectionPetsc(subdomain_num, sub_y, y);
    return;
} // accumulateCorrectionPetsc

void
PETScLevelSolver::applyAdditivePetsc(Vec x, Vec y)
{
    auto& shell = d_shell_data;
    auto& petsc = *shell.petsc_data;
    int ierr = VecZeroEntries(y);
    IBTK_CHKERRQ(ierr);
    for (int i = 0; i < shell.n_local_subdomains; ++i)
    {
        ierr = VecScatterBegin(petsc.restriction[i], x, petsc.sub_x[i], INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
    }
    for (int i = 0; i < shell.n_local_subdomains; ++i)
    {
        ierr = VecScatterEnd(petsc.restriction[i], x, petsc.sub_x[i], INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSolve(petsc.sub_ksp[i], petsc.sub_x[i], petsc.sub_y[i]);
        IBTK_CHKERRQ(ierr);
        accumulateCorrectionPetsc(i, petsc.sub_y[i], y);
    }
    postprocessShellResult(y);
    return;
} // applyAdditivePetsc

void
PETScLevelSolver::applyMultiplicativePetsc(Vec x, Vec y)
{
    auto& shell = d_shell_data;
    auto& petsc = *shell.petsc_data;
    int ierr = VecZeroEntries(y);
    IBTK_CHKERRQ(ierr);
    ierr = VecCopy(x, petsc.shell_r);
    IBTK_CHKERRQ(ierr);
    for (int i = 0; i < shell.n_local_subdomains; ++i)
    {
        ierr = VecScatterBegin(petsc.restriction[i], petsc.shell_r, petsc.sub_x[i], INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
        ierr = VecScatterEnd(petsc.restriction[i], petsc.shell_r, petsc.sub_x[i], INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSolve(petsc.sub_ksp[i], petsc.sub_x[i], petsc.sub_y[i]);
        IBTK_CHKERRQ(ierr);
        accumulateCorrectionPetsc(i, petsc.sub_y[i], y);
        if (i + 1 < shell.n_local_subdomains)
        {
            updateResidualPetsc(i, petsc.sub_y[i], petsc.shell_r);
        }
    }
    postprocessShellResult(y);
    return;
} // applyMultiplicativePetsc

void
PETScLevelSolver::updateResidualPetsc(const int subdomain_num, Vec sub_y, Vec residual)
{
    auto& shell = d_shell_data;
    auto& petsc = *shell.petsc_data;
    TBOX_ASSERT(petsc.active_residual_update_mat);

    const auto& active_update_local_positions = petsc.active_update_local_positions[subdomain_num];
    if (active_update_local_positions.empty()) return;

    int ierr;
    const PetscScalar* sub_y_arr = nullptr;
    PetscScalar* update_x_arr = nullptr;
    ierr = VecGetArrayRead(sub_y, &sub_y_arr);
    IBTK_CHKERRQ(ierr);
    ierr = VecGetArray(petsc.active_residual_update_x[subdomain_num], &update_x_arr);
    IBTK_CHKERRQ(ierr);
    for (std::size_t local_pos = 0; local_pos < active_update_local_positions.size(); ++local_pos)
    {
        update_x_arr[local_pos] = sub_y_arr[active_update_local_positions[local_pos]];
    }
    ierr = VecRestoreArray(petsc.active_residual_update_x[subdomain_num], &update_x_arr);
    IBTK_CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(sub_y, &sub_y_arr);
    IBTK_CHKERRQ(ierr);

    ierr = MatMult(petsc.active_residual_update_mat[subdomain_num],
                   petsc.active_residual_update_x[subdomain_num],
                   petsc.active_residual_update_y[subdomain_num]);
    IBTK_CHKERRQ(ierr);

    PetscScalar* residual_arr = nullptr;
    const PetscScalar* update_y_arr = nullptr;
    PetscInt n_local_entries = 0;
    PetscInt n_update_entries = 0;
    ierr = VecGetLocalSize(residual, &n_local_entries);
    IBTK_CHKERRQ(ierr);
    ierr = VecGetLocalSize(petsc.active_residual_update_y[subdomain_num], &n_update_entries);
    IBTK_CHKERRQ(ierr);
    TBOX_ASSERT(n_local_entries == n_update_entries);
    ierr = VecGetArray(residual, &residual_arr);
    IBTK_CHKERRQ(ierr);
    ierr = VecGetArrayRead(petsc.active_residual_update_y[subdomain_num], &update_y_arr);
    IBTK_CHKERRQ(ierr);
    for (PetscInt local_idx = 0; local_idx < n_local_entries; ++local_idx)
    {
        residual_arr[local_idx] -= update_y_arr[local_idx];
    }
    ierr = VecRestoreArrayRead(petsc.active_residual_update_y[subdomain_num], &update_y_arr);
    IBTK_CHKERRQ(ierr);
    ierr = VecRestoreArray(residual, &residual_arr);
    IBTK_CHKERRQ(ierr);
    return;
} // updateResidualPetsc

void
PETScLevelSolver::applyAdditiveEigen(Vec x, Vec y)
{
    if (usingCustomEigenSolveSubdomainSolver())
    {
        applyAdditiveCustomEigen(x, y);
        return;
    }

    switch (d_shell_smoother_backend)
    {
    case ShellSmootherBackend::EIGEN:
        dispatchEigenSolverType(getEigenSolveSolverType(),
                                [&](auto solver_tag)
                                {
                                    using SolverType = typename decltype(solver_tag)::type;
                                    applyAdditiveEigenSolveImpl<SolverType>(x, y);
                                });
        break;
    case ShellSmootherBackend::EIGEN_PSEUDOINVERSE:
    {
        TBOX_ASSERT(IBTK_MPI::getNodes() == 1);
        auto& eigen = *d_shell_data.eigen_data;
        const Eigen::Index n = eigen.n_dofs;
        TBOX_ASSERT(n > 0);

        {
            ConstPetscVecArrayMap x_array(x, n);
            PetscVecArrayMap y_array(y, n);
            const auto x_map = x_array.getMap();
            auto y_map = y_array.getMap();
            y_map.setZero();
            for (std::size_t subdomain_num = 0; subdomain_num < eigen.common_subdomains.size(); ++subdomain_num)
            {
                auto& cache = eigen.common_subdomains[subdomain_num];
                std::size_t rhs_idx = 0;
                for (const int dof : cache.overlap_dofs)
                {
                    cache.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = x_map[dof];
                }
                cache.delta_workspace.noalias() = cache.local_pseudoinverse * cache.rhs_workspace;
                std::size_t update_idx = 0;
                for (const int dof : cache.update_dofs)
                {
                    y_map[dof] +=
                        cache.delta_workspace[static_cast<Eigen::Index>(cache.update_local_positions[update_idx++])];
                }
            }
        }
        postprocessShellResult(y);
        break;
    }
    case ShellSmootherBackend::PETSC:
        TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::applyAdditiveEigen()\n"
                                 << "PETSc shell backend should not use the Eigen additive path." << std::endl);
        break;
    }
    return;
} // applyAdditiveEigen

void
PETScLevelSolver::applyMultiplicativeEigen(Vec x, Vec y)
{
    if (usingCustomEigenSolveSubdomainSolver())
    {
        applyMultiplicativeCustomEigen(x, y);
        return;
    }

    switch (d_shell_smoother_backend)
    {
    case ShellSmootherBackend::EIGEN:
        dispatchEigenSolverType(getEigenSolveSolverType(),
                                [&](auto solver_tag)
                                {
                                    using SolverType = typename decltype(solver_tag)::type;
                                    applyMultiplicativeEigenSolveImpl<SolverType>(x, y);
                                });
        break;
    case ShellSmootherBackend::EIGEN_PSEUDOINVERSE:
    {
        TBOX_ASSERT(IBTK_MPI::getNodes() == 1);
        auto& shell = d_shell_data;
        auto& eigen = *shell.eigen_data;
        const Eigen::Index n = eigen.n_dofs;
        TBOX_ASSERT(n > 0);

        {
            ConstPetscVecArrayMap x_array(x, n);
            PetscVecArrayMap y_array(y, n);
            const auto x_map = x_array.getMap();
            auto y_map = y_array.getMap();
            Eigen::VectorXd residual(n);
            y_map.setZero();
            residual = x_map;
            const std::size_t n_subdomains = static_cast<std::size_t>(shell.n_local_subdomains);
            for (std::size_t subdomain_num = 0; subdomain_num + 1 < n_subdomains; ++subdomain_num)
            {
                auto& cache = eigen.common_subdomains[subdomain_num];
                std::size_t rhs_idx = 0;
                for (const int dof : cache.overlap_dofs)
                {
                    cache.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = residual[dof];
                }
                cache.delta_workspace.noalias() = cache.local_pseudoinverse * cache.rhs_workspace;

                std::size_t update_idx = 0;
                for (const int dof : cache.update_dofs)
                {
                    y_map[dof] +=
                        cache.delta_workspace[static_cast<Eigen::Index>(cache.update_local_positions[update_idx++])];
                }
                if (cache.active_residual_update_mat.rows() > 0)
                {
                    std::size_t residual_input_idx = 0;
                    for (const int local_pos : cache.update_local_positions)
                    {
                        cache.residual_input_workspace[static_cast<Eigen::Index>(residual_input_idx++)] =
                            cache.delta_workspace[static_cast<Eigen::Index>(local_pos)];
                    }
                    cache.residual_delta_workspace.noalias() =
                        cache.active_residual_update_mat * cache.residual_input_workspace;
                    std::size_t row_idx = 0;
                    for (const int row : cache.active_residual_update_rows)
                    {
                        residual[row] -= cache.residual_delta_workspace[static_cast<Eigen::Index>(row_idx++)];
                    }
                }
            }
            if (n_subdomains > 0)
            {
                auto& cache = eigen.common_subdomains.back();
                std::size_t rhs_idx = 0;
                for (const int dof : cache.overlap_dofs)
                {
                    cache.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = residual[dof];
                }
                cache.delta_workspace.noalias() = cache.local_pseudoinverse * cache.rhs_workspace;

                std::size_t update_idx = 0;
                for (const int dof : cache.update_dofs)
                {
                    y_map[dof] +=
                        cache.delta_workspace[static_cast<Eigen::Index>(cache.update_local_positions[update_idx++])];
                }
            }
        }
        postprocessShellResult(y);
        break;
    }
    case ShellSmootherBackend::PETSC:
        TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::applyMultiplicativeEigen()\n"
                                 << "PETSc shell backend should not use the Eigen multiplicative path." << std::endl);
        break;
    }
    return;
} // applyMultiplicativeEigen

void
PETScLevelSolver::setupNullSpace()
{
    int ierr;
    std::vector<Vec> petsc_nullspace_basis_vecs(d_nullspace_basis_vecs.size());
    for (unsigned k = 0; k < d_nullspace_basis_vecs.size(); ++k)
    {
        Vec& petsc_nullspace_vec = petsc_nullspace_basis_vecs[k];
        ierr = MatCreateVecs(d_petsc_mat, nullptr, &petsc_nullspace_vec);
        IBTK_CHKERRQ(ierr);
        copyToPETScVec(petsc_nullspace_vec, *d_nullspace_basis_vecs[k]);
        double norm;
        ierr = VecNorm(petsc_nullspace_vec, NORM_2, &norm);
        IBTK_CHKERRQ(ierr);
        ierr = VecScale(petsc_nullspace_vec, 1.0 / norm);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,
                              d_nullspace_contains_constant_vec ? PETSC_TRUE : PETSC_FALSE,
                              static_cast<int>(petsc_nullspace_basis_vecs.size()),
                              petsc_nullspace_basis_vecs.data(),
                              &d_petsc_nullsp);
    IBTK_CHKERRQ(ierr);
    ierr = MatSetNullSpace(d_petsc_mat, d_petsc_nullsp);
    IBTK_CHKERRQ(ierr);
    for (unsigned k = 0; k < d_nullspace_basis_vecs.size(); ++k)
    {
        ierr = VecDestroy(&petsc_nullspace_basis_vecs[k]);
        IBTK_CHKERRQ(ierr);
    }
    return;
} // setupNullSpace

/////////////////////////////// PRIVATE //////////////////////////////////////

PetscErrorCode
PETScLevelSolver::PCApply_AdditivePetsc(PC pc, Vec x, Vec y)
{
    PetscFunctionBeginUser;
    int ierr;
    void* ctx;
    ierr = PCShellGetContext(pc, &ctx);
    CHKERRQ(ierr);
    auto solver = static_cast<PETScLevelSolver*>(ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(solver);
#endif
    solver->applyAdditivePetsc(x, y);
    PetscFunctionReturn(0);
} // PCApply_AdditivePetsc

PetscErrorCode
PETScLevelSolver::PCApply_MultiplicativePetsc(PC pc, Vec x, Vec y)
{
    PetscFunctionBeginUser;
    int ierr;
    void* ctx;
    ierr = PCShellGetContext(pc, &ctx);
    CHKERRQ(ierr);
    auto solver = static_cast<PETScLevelSolver*>(ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(solver);
#endif
    solver->applyMultiplicativePetsc(x, y);
    PetscFunctionReturn(0);
} // PCApply_MultiplicativePetsc

PetscErrorCode
PETScLevelSolver::PCApply_AdditiveEigen(PC pc, Vec x, Vec y)
{
    PetscFunctionBeginUser;
    int ierr;
    void* ctx;
    ierr = PCShellGetContext(pc, &ctx);
    CHKERRQ(ierr);
    auto solver = static_cast<PETScLevelSolver*>(ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(solver);
#endif
    if (IBTK_MPI::getNodes() != 1)
    {
        TBOX_ERROR("PETScLevelSolver::PCApply_AdditiveEigen(): only supports 1 MPI rank" << std::endl);
    }
    solver->applyAdditiveEigen(x, y);
    PetscFunctionReturn(0);
} // PCApply_AdditiveEigen

PetscErrorCode
PETScLevelSolver::PCApply_MultiplicativeEigen(PC pc, Vec x, Vec y)
{
    PetscFunctionBeginUser;
    int ierr;
    void* ctx;
    ierr = PCShellGetContext(pc, &ctx);
    CHKERRQ(ierr);
    auto solver = static_cast<PETScLevelSolver*>(ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(solver);
#endif
    if (IBTK_MPI::getNodes() != 1)
    {
        TBOX_ERROR("PETScLevelSolver::PCApply_MultiplicativeEigen(): only supports 1 MPI rank" << std::endl);
    }
    solver->applyMultiplicativeEigen(x, y);
    PetscFunctionReturn(0);
} // PCApply_MultiplicativeEigen

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
