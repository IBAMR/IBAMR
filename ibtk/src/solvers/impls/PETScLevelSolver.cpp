// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
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

#include <ibtk/DOFCoverage.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PETScLevelSolver.h>
#include <ibtk/PETScLevelSolverSubdomainSolver.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/string_utilities.h>

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

#include <CoarseFineBoundary.h>
#include <IntVector.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SAMRAIVectorReal.h>

#include <algorithm>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <utility>
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

// Return the built-in subdomain solver with the given name, or nothing if the name is not that of one.
std::optional<PETScLevelSolverSubdomainSolver>
make_built_in_subdomain_solver(const std::string& subdomain_solver_type, Pointer<Database> input_db)
{
    if (equals_ignore_case(subdomain_solver_type, "petsc"))
    {
        return make_petsc_subdomain_solver();
    }
    if (equals_ignore_case(subdomain_solver_type, "blas-lapack"))
    {
        return make_blas_lapack_subdomain_solver(input_db);
    }
    if (equals_ignore_case(subdomain_solver_type, "eigen"))
    {
        return make_eigen_subdomain_solver(input_db);
    }
    if (equals_ignore_case(subdomain_solver_type, "eigen-pseudoinverse"))
    {
        return make_eigen_pseudoinverse_subdomain_solver(input_db);
    }
    return std::nullopt;
}

// The columns of the sequential matrix rows that are in the list of global column indices, in the order of the
// list, scaled by scale. Entries of rows in other columns are dropped. The list has no repeated indices. This does
// O(nnz log n_columns) work and, unlike MatCreateSubMatrix() with an index set of the columns, allocates nothing
// of the size of the columns of rows.
Mat
extract_columns(Mat rows, const PetscInt* columns, const PetscInt n_columns, const PetscScalar scale)
{
    std::vector<std::pair<PetscInt, PetscInt>> positions(n_columns);
    for (PetscInt k = 0; k < n_columns; ++k) positions[k] = std::make_pair(columns[k], k);
    std::sort(positions.begin(), positions.end());
    const auto position_of = [&](const PetscInt column) -> PetscInt
    {
        const auto found = std::lower_bound(positions.begin(), positions.end(), std::make_pair(column, PetscInt(0)));
        return found != positions.end() && found->first == column ? found->second : -1;
    };
    PetscInt n_rows = 0;
    int ierr = MatGetSize(rows, &n_rows, nullptr);
    IBTK_CHKERRQ(ierr);
    std::vector<PetscInt> row_counts(n_rows, 0);
    for (PetscInt row = 0; row < n_rows; ++row)
    {
        PetscInt count = 0;
        const PetscInt* indices = nullptr;
        ierr = MatGetRow(rows, row, &count, &indices, nullptr);
        IBTK_CHKERRQ(ierr);
        for (PetscInt k = 0; k < count; ++k) row_counts[row] += position_of(indices[k]) >= 0;
        ierr = MatRestoreRow(rows, row, &count, &indices, nullptr);
        IBTK_CHKERRQ(ierr);
    }
    Mat extracted = nullptr;
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, n_rows, n_columns, 0, row_counts.data(), &extracted);
    IBTK_CHKERRQ(ierr);
    std::vector<PetscInt> extracted_columns;
    std::vector<PetscScalar> extracted_values;
    for (PetscInt row = 0; row < n_rows; ++row)
    {
        PetscInt count = 0;
        const PetscInt* indices = nullptr;
        const PetscScalar* values = nullptr;
        ierr = MatGetRow(rows, row, &count, &indices, &values);
        IBTK_CHKERRQ(ierr);
        extracted_columns.clear();
        extracted_values.clear();
        for (PetscInt k = 0; k < count; ++k)
        {
            const PetscInt position = position_of(indices[k]);
            if (position >= 0)
            {
                extracted_columns.push_back(position);
                extracted_values.push_back(scale * values[k]);
            }
        }
        ierr = MatRestoreRow(rows, row, &count, &indices, &values);
        IBTK_CHKERRQ(ierr);
        ierr = MatSetValues(extracted,
                            1,
                            &row,
                            static_cast<PetscInt>(extracted_columns.size()),
                            extracted_columns.data(),
                            extracted_values.data(),
                            INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(extracted, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(extracted, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    return extracted;
}

void
generate_petsc_is_from_std_is(std::vector<std::set<int>>& overlap_std,
                              std::vector<std::set<int>>& nonoverlap_std,
                              std::vector<IS>& overlap_petsc,
                              std::vector<IS>& nonoverlap_petsc)
{
    // Destroy old IS'es and generate new ones.
    int ierr;
    for (auto& is : overlap_petsc)
    {
        ierr = ISDestroy(&is);
        IBTK_CHKERRQ(ierr);
    }
    overlap_petsc.clear();
    for (auto& is : nonoverlap_petsc)
    {
        ierr = ISDestroy(&is);
        IBTK_CHKERRQ(ierr);
    }
    nonoverlap_petsc.clear();

    const int n_overlap_subdomains = static_cast<int>(overlap_std.size());
    overlap_petsc.resize(n_overlap_subdomains);
    for (int k = 0; k < n_overlap_subdomains; ++k)
    {
        PetscInt* overlap_dof_arr;
        const int n_overlap_dofs = static_cast<int>(overlap_std[k].size());
        ierr = PetscMalloc1(n_overlap_dofs, &overlap_dof_arr);
        IBTK_CHKERRQ(ierr);
        std::copy(overlap_std[k].begin(), overlap_std[k].end(), overlap_dof_arr);
        ierr = ISCreateGeneral(PETSC_COMM_SELF, n_overlap_dofs, overlap_dof_arr, PETSC_OWN_POINTER, &overlap_petsc[k]);
        IBTK_CHKERRQ(ierr);
    }

    const int n_nonoverlap_subdomains = static_cast<int>(nonoverlap_std.size());
    nonoverlap_petsc.resize(n_nonoverlap_subdomains);
    for (int k = 0; k < n_nonoverlap_subdomains; ++k)
    {
        PetscInt* nonoverlap_dof_arr;
        const int n_nonoverlap_dofs = static_cast<int>(nonoverlap_std[k].size());
        ierr = PetscMalloc1(n_nonoverlap_dofs, &nonoverlap_dof_arr);
        IBTK_CHKERRQ(ierr);
        std::copy(nonoverlap_std[k].begin(), nonoverlap_std[k].end(), nonoverlap_dof_arr);
        ierr = ISCreateGeneral(
            PETSC_COMM_SELF, n_nonoverlap_dofs, nonoverlap_dof_arr, PETSC_OWN_POINTER, &nonoverlap_petsc[k]);
        IBTK_CHKERRQ(ierr);
    }

    return;
} // generate_petsc_is_from_std_is
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

PETScLevelSolver::PETScLevelSolver()
{
    // Setup default options.
    d_max_iterations = 10000;
    d_abs_residual_tol = 1.0e-50;
    d_rel_residual_tol = 1.0e-5;
    d_initial_guess_nonzero = true;
    d_enable_logging = false;
    d_box_size = 2;
    d_overlap_size = 1;

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

    int ierr;
    for (auto& is : d_nonoverlap_is)
    {
        ierr = ISDestroy(&is);
        IBTK_CHKERRQ(ierr);
    }
    for (auto& is : d_overlap_is)
    {
        ierr = ISDestroy(&is);
        IBTK_CHKERRQ(ierr);
    }
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

void
PETScLevelSolver::setSubdomainSolver(PETScLevelSolverSubdomainSolver subdomain_solver)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setSubdomainSolver():\n"
                                 << "  cannot be called while solver state is initialized.\n");
    }
    if (!subdomain_solver)
    {
        TBOX_ERROR(d_object_name << "::setSubdomainSolver():\n"
                                 << "  the subdomain solver must not be empty.\n");
    }
    d_subdomain_solver = std::move(subdomain_solver);
    return;
} // setSubdomainSolver

const KSP&
PETScLevelSolver::getPETScKSP() const
{
    return d_petsc_ksp;
} // getPETScKSP

void
PETScLevelSolver::getASMSubdomains(std::vector<IS>** nonoverlapping_subdomains,
                                   std::vector<IS>** overlapping_subdomains)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_is_initialized);
#endif
    *nonoverlapping_subdomains = &d_nonoverlap_is;
    *overlapping_subdomains = &d_overlap_is;
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
    validatePreconditionerType();

    // Set the nullspace.
    if (d_nullspace_contains_constant_vec || !d_nullspace_basis_vecs.empty()) setupNullSpace();

    // Setup the preconditioner.
    if (d_pc_type == "asm")
    {
        // Generate user-defined subdomains.
        std::vector<std::set<int>> overlap_is, nonoverlap_is;
        generateASMSubdomains(overlap_is, nonoverlap_is);

        // Generate PETSc IS in cases where they have not been generated directly.
        if (!d_overlap_is.size())
        {
            generate_petsc_is_from_std_is(overlap_is, nonoverlap_is, d_overlap_is, d_nonoverlap_is);
            d_generated_subdomain_is = true;
        }

        int num_subdomains = static_cast<int>(d_overlap_is.size());
        if (num_subdomains == 0)
        {
            IS is;
            ierr = ISCreateGeneral(PETSC_COMM_SELF, 0, nullptr, PETSC_OWN_POINTER, &is);
            IBTK_CHKERRQ(ierr);
            ierr = PCASMSetLocalSubdomains(ksp_pc, 1, &is, &is);
            IBTK_CHKERRQ(ierr);
            ierr = ISDestroy(&is);
            IBTK_CHKERRQ(ierr);
        }
        else
        {
            ierr = PCASMSetLocalSubdomains(ksp_pc, num_subdomains, d_overlap_is.data(), d_nonoverlap_is.data());
            IBTK_CHKERRQ(ierr);
        }
    }

    if (d_pc_type == "fieldsplit")
    {
        std::vector<std::set<int>> field_is;
        std::vector<std::string> field_name;
        generateFieldSplitSubdomains(field_name, field_is);
        d_field_name = field_name;
        const int n_fields = static_cast<int>(field_is.size());

        // Destroy old IS'es and generate new ones.
        for (auto& is : d_field_is)
        {
            ierr = ISDestroy(&is);
            IBTK_CHKERRQ(ierr);
        }
        d_field_is.clear();

        d_field_is.resize(n_fields);
        for (int k = 0; k < n_fields; ++k)
        {
            PetscInt* field_dof_arr;
            const int n_field_dofs = static_cast<int>(field_is[k].size());
            ierr = PetscMalloc1(n_field_dofs, &field_dof_arr);
            IBTK_CHKERRQ(ierr);
            std::copy(field_is[k].begin(), field_is[k].end(), field_dof_arr);
            ierr = ISCreateGeneral(PETSC_COMM_WORLD, n_field_dofs, field_dof_arr, PETSC_OWN_POINTER, &d_field_is[k]);
            IBTK_CHKERRQ(ierr);
            ierr = PCFieldSplitSetIS(ksp_pc, d_field_name[k].c_str(), d_field_is[k]);
            IBTK_CHKERRQ(ierr);
        }
    }

    if (d_pc_type == "shell")
    {
        if (!d_shell_composition)
        {
            TBOX_ERROR(d_object_name << "::initializeSolverState():\n"
                                     << "  shell_pc_type must be \"additive\" or \"multiplicative\" when the "
                                        "preconditioner is a shell.\n");
        }

        const bool multiplicative = d_shell_composition == ShellComposition::MULTIPLICATIVE;
        if (!multiplicative && d_shell_traversal != ShellTraversal::FORWARD)
        {
            TBOX_ERROR(d_object_name << "::initializeSolverState():\n"
                                     << "  shell_pc_subdomain_traversal applies only to multiplicative composition.\n");
        }

        // Generate user-defined subdomains.
        std::vector<std::set<int>> overlap_is, nonoverlap_is;
        generateASMSubdomains(overlap_is, nonoverlap_is);

        // Generate PETSc IS in cases where they have not been generated directly.
        if (!d_overlap_is.size())
        {
            generate_petsc_is_from_std_is(overlap_is, nonoverlap_is, d_overlap_is, d_nonoverlap_is);
            d_generated_subdomain_is = true;
        }
        d_n_local_subdomains = static_cast<int>(d_overlap_is.size());

        // Get the local submatrices. The multiplicative shell also needs the rows of each subdomain
        // with every column, from which the submatrices follow without another collective extraction.
        Mat* subdomain_rows = nullptr;
        if (multiplicative)
        {
            PetscInt n_columns = 0;
            ierr = MatGetSize(d_petsc_mat, nullptr, &n_columns);
            IBTK_CHKERRQ(ierr);
            IS all_columns = nullptr;
            ierr = ISCreateStride(PETSC_COMM_SELF, n_columns, 0, 1, &all_columns);
            IBTK_CHKERRQ(ierr);
            std::vector<IS> all_columns_of_subdomains(d_n_local_subdomains, all_columns);
            ierr = MatCreateSubMatrices(d_petsc_mat,
                                        d_n_local_subdomains,
                                        get_data_or_null(d_overlap_is),
                                        get_data_or_null(all_columns_of_subdomains),
                                        MAT_INITIAL_MATRIX,
                                        &subdomain_rows);
            IBTK_CHKERRQ(ierr);
            ierr = ISDestroy(&all_columns);
            IBTK_CHKERRQ(ierr);
            // d_sub_mat is populated below one matrix at a time, not by MatCreateSubMatrices(), but
            // MatDestroySubMatrices() (used for both branches) still expects an array of exactly this
            // shape: one extra trailing slot, which it reads to decide whether a type-specific reuse
            // context is present. PetscCalloc1(), not PetscMalloc1(), zero-initializes that slot so
            // MatDestroySubMatrices() finds no such context and falls back to destroying each matrix
            // individually, instead of reading uninitialized memory there.
            ierr = PetscCalloc1(d_n_local_subdomains + 1, &d_sub_mat);
            IBTK_CHKERRQ(ierr);
            for (int i = 0; i < d_n_local_subdomains; ++i)
            {
                PetscInt overlap_size = 0;
                const PetscInt* overlap_indices = nullptr;
                ierr = ISGetLocalSize(d_overlap_is[i], &overlap_size);
                IBTK_CHKERRQ(ierr);
                ierr = ISGetIndices(d_overlap_is[i], &overlap_indices);
                IBTK_CHKERRQ(ierr);
                d_sub_mat[i] = extract_columns(subdomain_rows[i], overlap_indices, overlap_size, 1.0);
                ierr = ISRestoreIndices(d_overlap_is[i], &overlap_indices);
                IBTK_CHKERRQ(ierr);
            }
        }
        else
        {
            ierr = MatCreateSubMatrices(d_petsc_mat,
                                        d_n_local_subdomains,
                                        d_overlap_is.data(),
                                        d_overlap_is.data(),
                                        MAT_INITIAL_MATRIX,
                                        &d_sub_mat);
            IBTK_CHKERRQ(ierr);
        }

        // The right-hand sides and solutions of the local problems of all subdomains are packed,
        // in order, into two sequential vectors, so that one scatter gathers every right-hand
        // side and subdomain solvers can work on the vectors in place.
        PetscInt n_lo, n_hi;
        ierr = VecGetOwnershipRange(d_petsc_x, &n_lo, &n_hi);
        IBTK_CHKERRQ(ierr);
        d_subdomain_offsets.assign(d_n_local_subdomains + 1, 0);
        d_write_offsets.assign(d_n_local_subdomains + 1, 0);
        d_write_sources.clear();
        d_write_targets.clear();
        std::vector<PetscInt> gathered_indices;
        for (int i = 0; i < d_n_local_subdomains; ++i)
        {
            PetscInt overlap_size = 0;
            ierr = ISGetLocalSize(d_overlap_is[i], &overlap_size);
            IBTK_CHKERRQ(ierr);
            const PetscInt* overlap_indices = nullptr;
            ierr = ISGetIndices(d_overlap_is[i], &overlap_indices);
            IBTK_CHKERRQ(ierr);
            gathered_indices.insert(gathered_indices.end(), overlap_indices, overlap_indices + overlap_size);

            if (!multiplicative)
            {
                // The subdomain's nonoverlapping DOFs are owned by this rank. Since the IS are sorted, they
                // are found in one pass over the overlapping DOFs.
                PetscInt nonoverlap_size = 0;
                ierr = ISGetLocalSize(d_nonoverlap_is[i], &nonoverlap_size);
                IBTK_CHKERRQ(ierr);
                const PetscInt* nonoverlap_indices = nullptr;
                ierr = ISGetIndices(d_nonoverlap_is[i], &nonoverlap_indices);
                IBTK_CHKERRQ(ierr);
                PetscInt jj = 0;
                for (PetscInt ii = 0; ii < overlap_size && jj < nonoverlap_size; ++ii)
                {
                    if (overlap_indices[ii] == nonoverlap_indices[jj])
                    {
                        if (nonoverlap_indices[jj] < n_lo || nonoverlap_indices[jj] >= n_hi)
                        {
                            TBOX_ERROR(d_object_name << "::initializeSolverState():\n"
                                                     << "  nonoverlapping DOF " << nonoverlap_indices[jj]
                                                     << " of subdomain " << i << " is not owned by this rank.\n");
                        }
                        d_write_sources.push_back(d_subdomain_offsets[i] + ii);
                        d_write_targets.push_back(nonoverlap_indices[jj] - n_lo);
                        ++jj;
                    }
                }
                if (jj != nonoverlap_size)
                {
                    TBOX_ERROR(d_object_name << "::initializeSolverState():\n"
                                             << "  the nonoverlapping set of subdomain " << i
                                             << " is not a subset of its overlapping set.\n");
                }
                ierr = ISRestoreIndices(d_nonoverlap_is[i], &nonoverlap_indices);
                IBTK_CHKERRQ(ierr);
            }
            d_subdomain_offsets[i + 1] = d_subdomain_offsets[i] + overlap_size;
            d_write_offsets[i + 1] = static_cast<PetscInt>(d_write_sources.size());
            ierr = ISRestoreIndices(d_overlap_is[i], &overlap_indices);
            IBTK_CHKERRQ(ierr);
        }
        // The additive preconditioner needs the nonoverlapping subsets to partition the DOFs.
        if (!multiplicative && d_check_subdomain_coverage)
        {
            check_dof_coverage(
                d_object_name + "::initializeSolverState()", d_nonoverlap_is, n_hi - n_lo, DOFCoverage::EXACTLY_ONCE);
        }
        const PetscInt n_packed = d_subdomain_offsets[d_n_local_subdomains];
        ierr = VecCreateSeq(PETSC_COMM_SELF, n_packed, &d_subdomain_rhs);
        IBTK_CHKERRQ(ierr);
        ierr = VecDuplicate(d_subdomain_rhs, &d_subdomain_solution);
        IBTK_CHKERRQ(ierr);
        // The gather needs communication only if a subdomain of some rank includes DOFs that another rank
        // owns. Every rank makes the same choice, so the collective scatter is created and applied by all
        // ranks or by none.
        int communicates = 0;
        for (const PetscInt dof : gathered_indices)
        {
            communicates = communicates || dof < n_lo || dof >= n_hi;
        }
        IBTK_MPI::maxReduction(&communicates, 1);
        d_restriction_communicates = communicates != 0;
        if (d_restriction_communicates)
        {
            IS gathered_is, packed_is;
            ierr = ISCreateGeneral(PETSC_COMM_SELF, n_packed, gathered_indices.data(), PETSC_COPY_VALUES, &gathered_is);
            IBTK_CHKERRQ(ierr);
            ierr = ISCreateStride(PETSC_COMM_SELF, n_packed, 0, 1, &packed_is);
            IBTK_CHKERRQ(ierr);
            ierr = VecScatterCreate(d_petsc_x, gathered_is, d_subdomain_rhs, packed_is, &d_restriction);
            IBTK_CHKERRQ(ierr);
            ierr = ISDestroy(&gathered_is);
            IBTK_CHKERRQ(ierr);
            ierr = ISDestroy(&packed_is);
            IBTK_CHKERRQ(ierr);
        }
        else
        {
            d_gather_indices = gathered_indices;
            for (PetscInt& dof : d_gather_indices)
            {
                dof -= n_lo;
            }
        }
        if (multiplicative)
        {
            initializeMultiplicativeShell(subdomain_rows);
            ierr = MatDestroySubMatrices(d_n_local_subdomains, &subdomain_rows);
            IBTK_CHKERRQ(ierr);
        }

        // Set up the subdomain solvers.
        if (!d_subdomain_solver)
        {
            TBOX_ERROR(
                d_object_name << "::initializeSolverState():\n"
                              << "  the subdomain solver " << d_subdomain_solver_type
                              << " has not been set; a derived class that names it must call setSubdomainSolver().\n");
        }
        d_subdomain_solver->initializeSolverState(
            std::vector<Mat>(d_sub_mat, d_sub_mat + d_n_local_subdomains), d_overlap_is, d_options_prefix);
        d_subdomain_solver_initialized = true;
        ierr = PCSetType(ksp_pc, PCSHELL);
        IBTK_CHKERRQ(ierr);
        ierr = PCShellSetContext(ksp_pc, static_cast<void*>(this));
        IBTK_CHKERRQ(ierr);
        if (*d_shell_composition == ShellComposition::ADDITIVE)
        {
            ierr = PCShellSetApply(ksp_pc, PETScLevelSolver::PCApply_Additive);
            IBTK_CHKERRQ(ierr);
            std::string pc_name = d_options_prefix + "PC_Additive";
            ierr = PCShellSetName(ksp_pc, pc_name.c_str());
            IBTK_CHKERRQ(ierr);
        }
        else
        {
            ierr = PCShellSetApply(ksp_pc, PETScLevelSolver::PCApply_Multiplicative);
            IBTK_CHKERRQ(ierr);
            std::string pc_name = d_options_prefix + "PC_Multiplicative";
            ierr = PCShellSetName(ksp_pc, pc_name.c_str());
            IBTK_CHKERRQ(ierr);
        }
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

    if (d_subdomain_solver_initialized)
    {
        d_subdomain_solver->deallocateSolverState();
        d_subdomain_solver_initialized = false;
    }

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
    if (d_pc_type == "shell")
    {
        for (std::vector<Vec>* vectors :
             { &d_subdomain_rhs_views, &d_subdomain_residual_views, &d_subdomain_solution_views, &d_halo_vectors })
        {
            for (Vec& vector : *vectors)
            {
                ierr = VecDestroy(&vector);
                IBTK_CHKERRQ(ierr);
            }
            vectors->clear();
        }
        for (std::vector<VecScatter>* scatters : { &d_halo_scatters, &d_correction_scatters })
        {
            for (VecScatter& scatter : *scatters)
            {
                ierr = VecScatterDestroy(&scatter);
                IBTK_CHKERRQ(ierr);
            }
            scatters->clear();
        }
        for (Mat& matrix : d_residual_matrices)
        {
            ierr = MatDestroy(&matrix);
            IBTK_CHKERRQ(ierr);
        }
        d_residual_matrices.clear();
        d_n_stages = 0;
        d_halo_communicates.clear();
        d_correction_communicates.clear();
        d_halo_local_indices.clear();
        d_correction_local_indices.clear();
        d_stage_subdomains.clear();
        ierr = VecDestroy(&d_empty_vector);
        IBTK_CHKERRQ(ierr);
        ierr = VecScatterDestroy(&d_restriction);
        IBTK_CHKERRQ(ierr);
        d_gather_indices.clear();
        ierr = VecDestroy(&d_subdomain_rhs);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&d_subdomain_residual);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&d_subdomain_solution);
        IBTK_CHKERRQ(ierr);
        d_subdomain_offsets.clear();
        d_write_offsets.clear();
        d_write_sources.clear();
        d_write_targets.clear();
        ierr = MatDestroySubMatrices(d_n_local_subdomains, &d_sub_mat);
        IBTK_CHKERRQ(ierr);
        d_sub_mat = nullptr;
        d_n_local_subdomains = 0;
    }

    // Discard PETSc index sets converted from the subclass's std::set<int> lists, so that they
    // are rebuilt for the next layout. Subclasses that construct PETSc index sets directly
    // manage their regeneration.
    if (d_generated_subdomain_is)
    {
        for (IS& is : d_nonoverlap_is)
        {
            ierr = ISDestroy(&is);
            IBTK_CHKERRQ(ierr);
        }
        for (IS& is : d_overlap_is)
        {
            ierr = ISDestroy(&is);
            IBTK_CHKERRQ(ierr);
        }
        d_nonoverlap_is.clear();
        d_overlap_is.clear();
        d_generated_subdomain_is = false;
    }

    d_petsc_ksp = nullptr;
    d_petsc_mat = nullptr;
    d_petsc_x = nullptr;
    d_petsc_b = nullptr;

    // Indicate that the solver is NOT initialized.
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_solver_state);
    return;
} // deallocateSolverState

/////////////////////////////// PROTECTED ////////////////////////////////////

void
PETScLevelSolver::init(Pointer<Database> input_db,
                       const std::string& default_options_prefix,
                       const std::vector<std::string>& additional_subdomain_solver_names)
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
        if (input_db->keyExists("shell_pc_type"))
        {
            const std::string shell_pc_type = input_db->getString("shell_pc_type");
            if (equals_ignore_case(shell_pc_type, "additive"))
            {
                d_shell_composition = ShellComposition::ADDITIVE;
            }
            else if (equals_ignore_case(shell_pc_type, "multiplicative"))
            {
                d_shell_composition = ShellComposition::MULTIPLICATIVE;
            }
            else
            {
                TBOX_ERROR(d_object_name << "::init():\n"
                                         << "  unknown shell_pc_type = " << shell_pc_type
                                         << "; valid values are \"additive\" and \"multiplicative\".\n");
            }
        }
        d_subdomain_solver_type = input_db->getStringWithDefault("subdomain_solver", d_subdomain_solver_type);
        if (input_db->keyExists("shell_pc_subdomain_traversal"))
        {
            const std::string traversal = input_db->getString("shell_pc_subdomain_traversal");
            if (equals_ignore_case(traversal, "FORWARD"))
            {
                d_shell_traversal = ShellTraversal::FORWARD;
            }
            else if (equals_ignore_case(traversal, "REVERSE"))
            {
                d_shell_traversal = ShellTraversal::REVERSE;
            }
            else if (equals_ignore_case(traversal, "SYMMETRIC"))
            {
                d_shell_traversal = ShellTraversal::SYMMETRIC;
            }
            else
            {
                TBOX_ERROR(d_object_name << "::init():\n"
                                         << "  unknown shell_pc_subdomain_traversal = " << traversal
                                         << "; valid values are FORWARD, REVERSE, and SYMMETRIC.\n");
            }
        }
        if (input_db->keyExists("initial_guess_nonzero"))
            d_initial_guess_nonzero = input_db->getBool("initial_guess_nonzero");
        d_check_subdomain_coverage =
            input_db->getBoolWithDefault("check_subdomain_coverage", d_check_subdomain_coverage);
        if (input_db->keyExists("subdomain_box_size"))
            input_db->getIntegerArray("subdomain_box_size", d_box_size, NDIM);
        if (input_db->keyExists("subdomain_overlap_size"))
            input_db->getIntegerArray("subdomain_overlap_size", d_overlap_size, NDIM);
    }
    // A subclass that calls setSubdomainSolver() before PETScLevelSolver::init() (none does today) takes
    // precedence over subdomain_solver; this guard is what makes that ordering safe.
    if (!d_subdomain_solver)
    {
        d_subdomain_solver = make_built_in_subdomain_solver(d_subdomain_solver_type, input_db);
    }
    if (!d_subdomain_solver &&
        std::none_of(additional_subdomain_solver_names.begin(),
                     additional_subdomain_solver_names.end(),
                     [&](const std::string& name) { return equals_ignore_case(name, d_subdomain_solver_type); }))
    {
        std::vector<std::string> supported_names = { "petsc", "blas-lapack", "eigen", "eigen-pseudoinverse" };
        supported_names.insert(
            supported_names.end(), additional_subdomain_solver_names.begin(), additional_subdomain_solver_names.end());
        std::string supported;
        for (std::size_t k = 0; k < supported_names.size(); ++k)
        {
            supported += (k == 0                          ? "" :
                          k + 1 == supported_names.size() ? ", and " :
                                                            ", ") +
                         std::string("\"") + supported_names[k] + "\"";
        }
        TBOX_ERROR(d_object_name << "::init():\n"
                                 << "  unsupported subdomain_solver = " << d_subdomain_solver_type
                                 << "; supported values are " << supported << ".\n");
    }
    return;
} // init

void
PETScLevelSolver::generateASMSubdomains(std::vector<std::set<int>>& /*overlap_is*/,
                                        std::vector<std::set<int>>& /*nonoverlap_is*/)
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
                              get_data_or_null(petsc_nullspace_basis_vecs),
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

std::vector<int>
PETScLevelSolver::subdomainVisitOrder(const ShellTraversal traversal, const int n)
{
    std::vector<int> order;
    if (traversal == ShellTraversal::REVERSE)
    {
        for (int i = n - 1; i >= 0; --i) order.push_back(i);
    }
    else
    {
        for (int i = 0; i < n; ++i) order.push_back(i);
        if (traversal == ShellTraversal::SYMMETRIC)
        {
            for (int i = n - 2; i >= 0; --i) order.push_back(i);
        }
    }
    return order;
}

void
PETScLevelSolver::initializeMultiplicativeShell(Mat* rows)
{
    d_stage_subdomains = subdomainVisitOrder(d_shell_traversal, d_n_local_subdomains);
    d_n_stages = IBTK_MPI::maxReduction(static_cast<int>(d_stage_subdomains.size()));

    // The packed right-hand sides, residuals and solutions hold the subdomains one after another, so a
    // view of each subdomain shares their storage. The residual of a subdomain is computed from its
    // right-hand side at each visit, so that a subdomain may be visited more than once.
    int ierr = VecDuplicate(d_subdomain_rhs, &d_subdomain_residual);
    IBTK_CHKERRQ(ierr);
    Vec packed[3] = { d_subdomain_rhs, d_subdomain_residual, d_subdomain_solution };
    std::vector<Vec>* views[3] = { &d_subdomain_rhs_views, &d_subdomain_residual_views, &d_subdomain_solution_views };
    for (int k = 0; k < 3; ++k)
    {
        PetscScalar* storage = nullptr;
        ierr = VecGetArray(packed[k], &storage);
        IBTK_CHKERRQ(ierr);
        views[k]->resize(d_n_local_subdomains);
        for (int i = 0; i < d_n_local_subdomains; ++i)
        {
            ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,
                                         1,
                                         d_subdomain_offsets[i + 1] - d_subdomain_offsets[i],
                                         storage + d_subdomain_offsets[i],
                                         &(*views[k])[i]);
            IBTK_CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(packed[k], &storage);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecCreateSeq(PETSC_COMM_SELF, 0, &d_empty_vector);
    IBTK_CHKERRQ(ierr);

    // The rows of the operator for each subdomain contain the columns that the residual on the subdomain
    // depends on. The columns are numbered globally, so that the rows of subdomains that include DOFs of
    // other ranks are available. The matrix of each subdomain keeps only these columns, and is negated so
    // that the residual is a sum.
    d_residual_matrices.assign(d_n_local_subdomains, nullptr);
    std::vector<std::vector<PetscInt>> halos(d_n_local_subdomains);
    for (int i = 0; i < d_n_local_subdomains; ++i)
    {
        PetscInt n_rows = 0;
        ierr = MatGetSize(rows[i], &n_rows, nullptr);
        IBTK_CHKERRQ(ierr);
        for (PetscInt row = 0; row < n_rows; ++row)
        {
            PetscInt count = 0;
            const PetscInt* indices = nullptr;
            ierr = MatGetRow(rows[i], row, &count, &indices, nullptr);
            IBTK_CHKERRQ(ierr);
            halos[i].insert(halos[i].end(), indices, indices + count);
            ierr = MatRestoreRow(rows[i], row, &count, &indices, nullptr);
            IBTK_CHKERRQ(ierr);
        }
        std::sort(halos[i].begin(), halos[i].end());
        halos[i].erase(std::unique(halos[i].begin(), halos[i].end()), halos[i].end());
        d_residual_matrices[i] =
            extract_columns(rows[i], halos[i].data(), static_cast<PetscInt>(halos[i].size()), -1.0);
    }

    // A stage needs communication to gather the halo, or to add the correction, only if some rank visits a
    // subdomain at that stage with a DOF that another rank owns. Every rank learns this from one reduction at
    // setup, so the collective scatters of a stage are created and applied by all ranks or by none, and a stage
    // without communication works on the local array of the output.
    PetscInt n_lo = 0, n_hi = 0;
    ierr = VecGetOwnershipRange(d_petsc_x, &n_lo, &n_hi);
    IBTK_CHKERRQ(ierr);
    const auto has_remote_dof = [&](const PetscInt* dofs, const PetscInt n)
    {
        for (PetscInt k = 0; k < n; ++k)
        {
            if (dofs[k] < n_lo || dofs[k] >= n_hi) return true;
        }
        return false;
    };
    std::vector<int> communicates(2 * d_n_stages, 0);
    for (std::size_t stage = 0; stage < d_stage_subdomains.size(); ++stage)
    {
        const int i = d_stage_subdomains[stage];
        // The output is zero at the first stage, which needs no halo.
        communicates[2 * stage] = stage > 0 && has_remote_dof(halos[i].data(), static_cast<PetscInt>(halos[i].size()));
        PetscInt overlap_size = 0;
        const PetscInt* overlap_dofs = nullptr;
        ierr = ISGetLocalSize(d_overlap_is[i], &overlap_size);
        IBTK_CHKERRQ(ierr);
        ierr = ISGetIndices(d_overlap_is[i], &overlap_dofs);
        IBTK_CHKERRQ(ierr);
        communicates[2 * stage + 1] = has_remote_dof(overlap_dofs, overlap_size);
        ierr = ISRestoreIndices(d_overlap_is[i], &overlap_dofs);
        IBTK_CHKERRQ(ierr);
    }
    IBTK_MPI::maxReduction(communicates.data(), static_cast<int>(communicates.size()));

    d_halo_vectors.assign(d_n_stages, nullptr);
    d_halo_scatters.assign(d_n_stages, nullptr);
    d_correction_scatters.assign(d_n_stages, nullptr);
    d_halo_communicates.assign(d_n_stages, false);
    d_correction_communicates.assign(d_n_stages, false);
    d_halo_local_indices.assign(d_n_stages, {});
    d_correction_local_indices.assign(d_n_stages, {});
    for (int stage = 0; stage < d_n_stages; ++stage)
    {
        const bool visits = stage < static_cast<int>(d_stage_subdomains.size());
        const int i = visits ? d_stage_subdomains[stage] : -1;
        d_halo_communicates[stage] = communicates[2 * stage] != 0;
        d_correction_communicates[stage] = communicates[2 * stage + 1] != 0;
        const std::vector<PetscInt> no_dofs;
        const std::vector<PetscInt>& halo = visits && stage > 0 ? halos[i] : no_dofs;
        const PetscInt halo_size = static_cast<PetscInt>(halo.size());
        ierr = VecCreateSeq(PETSC_COMM_SELF, halo_size, &d_halo_vectors[stage]);
        IBTK_CHKERRQ(ierr);
        if (d_halo_communicates[stage])
        {
            // A rank that has visited all of its subdomains uses empty index sets.
            IS halo_is = nullptr, halo_positions = nullptr;
            ierr = ISCreateGeneral(PETSC_COMM_SELF, halo_size, halo.data(), PETSC_COPY_VALUES, &halo_is);
            IBTK_CHKERRQ(ierr);
            ierr = ISCreateStride(PETSC_COMM_SELF, halo_size, 0, 1, &halo_positions);
            IBTK_CHKERRQ(ierr);
            ierr = VecScatterCreate(d_petsc_x, halo_is, d_halo_vectors[stage], halo_positions, &d_halo_scatters[stage]);
            IBTK_CHKERRQ(ierr);
            ierr = ISDestroy(&halo_is);
            IBTK_CHKERRQ(ierr);
            ierr = ISDestroy(&halo_positions);
            IBTK_CHKERRQ(ierr);
        }
        else if (visits)
        {
            d_halo_local_indices[stage] = halo;
            for (PetscInt& dof : d_halo_local_indices[stage])
            {
                dof -= n_lo;
            }
        }

        PetscInt overlap_size = 0;
        const PetscInt* overlap_dofs = nullptr;
        if (visits)
        {
            ierr = ISGetLocalSize(d_overlap_is[i], &overlap_size);
            IBTK_CHKERRQ(ierr);
            ierr = ISGetIndices(d_overlap_is[i], &overlap_dofs);
            IBTK_CHKERRQ(ierr);
        }
        if (d_correction_communicates[stage])
        {
            IS overlap_is = nullptr, overlap_positions = nullptr;
            ierr = ISCreateGeneral(PETSC_COMM_SELF, overlap_size, overlap_dofs, PETSC_COPY_VALUES, &overlap_is);
            IBTK_CHKERRQ(ierr);
            ierr = ISCreateStride(PETSC_COMM_SELF, overlap_size, 0, 1, &overlap_positions);
            IBTK_CHKERRQ(ierr);
            ierr = VecScatterCreate(d_petsc_x,
                                    overlap_is,
                                    visits ? d_subdomain_solution_views[i] : d_empty_vector,
                                    overlap_positions,
                                    &d_correction_scatters[stage]);
            IBTK_CHKERRQ(ierr);
            ierr = ISDestroy(&overlap_is);
            IBTK_CHKERRQ(ierr);
            ierr = ISDestroy(&overlap_positions);
            IBTK_CHKERRQ(ierr);
        }
        else if (visits)
        {
            d_correction_local_indices[stage].assign(overlap_dofs, overlap_dofs + overlap_size);
            for (PetscInt& dof : d_correction_local_indices[stage])
            {
                dof -= n_lo;
            }
        }
        if (visits)
        {
            ierr = ISRestoreIndices(d_overlap_is[i], &overlap_dofs);
            IBTK_CHKERRQ(ierr);
        }
    }
}

PetscErrorCode
PETScLevelSolver::gatherSubdomainRhs(Vec x) const
{
    PetscFunctionBeginUser;
    int ierr = 0;
    if (d_restriction_communicates)
    {
        ierr = VecScatterBegin(d_restriction, x, d_subdomain_rhs, INSERT_VALUES, SCATTER_FORWARD);
        CHKERRQ(ierr);
        ierr = VecScatterEnd(d_restriction, x, d_subdomain_rhs, INSERT_VALUES, SCATTER_FORWARD);
        CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }
    // All of the DOFs of the subdomains are owned by their ranks, so gather them from the local array.
    const PetscScalar* x_values = nullptr;
    PetscScalar* rhs_values = nullptr;
    ierr = VecGetArrayRead(x, &x_values);
    CHKERRQ(ierr);
    ierr = VecGetArray(d_subdomain_rhs, &rhs_values);
    CHKERRQ(ierr);
    for (std::size_t k = 0; k < d_gather_indices.size(); ++k)
    {
        rhs_values[k] = x_values[d_gather_indices[k]];
    }
    ierr = VecRestoreArray(d_subdomain_rhs, &rhs_values);
    CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(x, &x_values);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // gatherSubdomainRhs

PetscErrorCode
PETScLevelSolver::writeSubdomainSolutions(const int first, const int last, Vec y) const
{
    PetscFunctionBeginUser;
    const PetscScalar* solution = nullptr;
    PetscScalar* y_values = nullptr;
    int ierr = VecGetArrayRead(d_subdomain_solution, &solution);
    CHKERRQ(ierr);
    ierr = VecGetArray(y, &y_values);
    CHKERRQ(ierr);
    for (PetscInt k = d_write_offsets[first]; k < d_write_offsets[last]; ++k)
    {
        y_values[d_write_targets[k]] = solution[d_write_sources[k]];
    }
    ierr = VecRestoreArray(y, &y_values);
    CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(d_subdomain_solution, &solution);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // writeSubdomainSolutions

PetscErrorCode
PETScLevelSolver::PCApply_Additive(PC pc, Vec x, Vec y)
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
    // writeSubdomainSolutions() only assigns the entries owned by a subdomain. A DOF that no subdomain
    // covers is otherwise left holding whatever y already contained instead of the zero correction an
    // uncovered DOF should get; check_dof_coverage() only catches this when check_subdomain_coverage is
    // on, and not every subclass requires a full covering.
    ierr = VecZeroEntries(y);
    CHKERRQ(ierr);
    ierr = solver->gatherSubdomainRhs(x);
    CHKERRQ(ierr);
    solver->d_subdomain_solver->solve(
        0, solver->d_n_local_subdomains, solver->d_subdomain_rhs, solver->d_subdomain_solution);
    ierr = solver->writeSubdomainSolutions(0, solver->d_n_local_subdomains, y);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // PCApply_Additive

PetscErrorCode
PETScLevelSolver::PCApply_Multiplicative(PC pc, Vec x, Vec y)
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
    ierr = VecZeroEntries(y);
    CHKERRQ(ierr);
    ierr = solver->gatherSubdomainRhs(x);
    CHKERRQ(ierr);
    for (int stage = 0; stage < solver->d_n_stages; ++stage)
    {
        // A rank that has visited all of its subdomains takes part in the scatters with empty ones.
        const bool visits = stage < static_cast<int>(solver->d_stage_subdomains.size());
        const int i = visits ? solver->d_stage_subdomains[stage] : -1;
        Vec rhs = solver->d_subdomain_rhs;
        // The output is zero at the first stage, so the residual there is the right-hand side.
        if (stage > 0)
        {
            if (solver->d_halo_communicates[stage])
            {
                ierr = VecScatterBegin(
                    solver->d_halo_scatters[stage], y, solver->d_halo_vectors[stage], INSERT_VALUES, SCATTER_FORWARD);
                CHKERRQ(ierr);
                ierr = VecScatterEnd(
                    solver->d_halo_scatters[stage], y, solver->d_halo_vectors[stage], INSERT_VALUES, SCATTER_FORWARD);
                CHKERRQ(ierr);
            }
            else if (visits)
            {
                const PetscScalar* y_values = nullptr;
                PetscScalar* halo_values = nullptr;
                ierr = VecGetArrayRead(y, &y_values);
                CHKERRQ(ierr);
                ierr = VecGetArray(solver->d_halo_vectors[stage], &halo_values);
                CHKERRQ(ierr);
                const std::vector<PetscInt>& indices = solver->d_halo_local_indices[stage];
                for (std::size_t k = 0; k < indices.size(); ++k)
                {
                    halo_values[k] = y_values[indices[k]];
                }
                ierr = VecRestoreArray(solver->d_halo_vectors[stage], &halo_values);
                CHKERRQ(ierr);
                ierr = VecRestoreArrayRead(y, &y_values);
                CHKERRQ(ierr);
            }
            if (visits)
            {
                ierr = MatMultAdd(solver->d_residual_matrices[i],
                                  solver->d_halo_vectors[stage],
                                  solver->d_subdomain_rhs_views[i],
                                  solver->d_subdomain_residual_views[i]);
                CHKERRQ(ierr);
                rhs = solver->d_subdomain_residual;
            }
        }
        if (visits)
        {
            solver->d_subdomain_solver->solve(i, i + 1, rhs, solver->d_subdomain_solution);
        }
        // Add the whole correction of the subdomain, including the DOFs of other ranks.
        Vec correction = visits ? solver->d_subdomain_solution_views[i] : solver->d_empty_vector;
        if (solver->d_correction_communicates[stage])
        {
            ierr = VecScatterBegin(solver->d_correction_scatters[stage], correction, y, ADD_VALUES, SCATTER_REVERSE);
            CHKERRQ(ierr);
            ierr = VecScatterEnd(solver->d_correction_scatters[stage], correction, y, ADD_VALUES, SCATTER_REVERSE);
            CHKERRQ(ierr);
        }
        else if (visits)
        {
            const PetscScalar* correction_values = nullptr;
            PetscScalar* y_values = nullptr;
            ierr = VecGetArrayRead(correction, &correction_values);
            CHKERRQ(ierr);
            ierr = VecGetArray(y, &y_values);
            CHKERRQ(ierr);
            const std::vector<PetscInt>& indices = solver->d_correction_local_indices[stage];
            for (std::size_t k = 0; k < indices.size(); ++k)
            {
                y_values[indices[k]] += correction_values[k];
            }
            ierr = VecRestoreArray(y, &y_values);
            CHKERRQ(ierr);
            ierr = VecRestoreArrayRead(correction, &correction_values);
            CHKERRQ(ierr);
        }
    }
    PetscFunctionReturn(0);
} // PCApply_Multiplicative

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
