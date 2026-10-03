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
#include <cstddef>
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

// The columns of the sequential matrix rows that are in the list of global column indices, in the order of the
// list, scaled by scale. Entries of rows in other columns are dropped. The list has no repeated indices. This does
// O(nnz log n_columns) work and, unlike MatCreateSubMatrix() with an index set of the columns, allocates nothing
// of the size of the columns of rows.
Mat
extract_columns(Mat rows, const PetscInt* columns, const PetscInt n_columns, const PetscScalar scale)
{
    std::vector<std::pair<PetscInt, PetscInt>> positions(n_columns);
    for (PetscInt k = 0; k < n_columns; ++k)
    {
        positions[k] = std::make_pair(columns[k], k);
    }
    std::sort(positions.begin(), positions.end());
    const auto position_of = [&](const PetscInt column) -> PetscInt
    {
        const std::vector<std::pair<PetscInt, PetscInt>>::const_iterator found =
            std::lower_bound(positions.begin(), positions.end(), std::make_pair(column, PetscInt(0)));
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
        for (PetscInt k = 0; k < count; ++k)
        {
            row_counts[row] += position_of(indices[k]) >= 0;
        }
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

// The position of dof in the sorted list dofs, or -1 if it is not in the list.
PetscInt
position_in(const std::vector<PetscInt>& dofs, const PetscInt dof)
{
    const std::vector<PetscInt>::const_iterator found = std::lower_bound(dofs.begin(), dofs.end(), dof);
    return found != dofs.end() && *found == dof ? static_cast<PetscInt>(found - dofs.begin()) : -1;
}

// Report a DOF of this rank, from n_lo up to n_hi, that OWNED output does not write exactly once, given the local
// indices of the entries that it writes.
void
check_owned_output(const std::string& context,
                   const std::vector<PetscInt>& targets,
                   const PetscInt n_lo,
                   const PetscInt n_hi)
{
    std::vector<int> writes(n_hi - n_lo, 0);
    for (const PetscInt target : targets)
    {
        ++writes[target];
    }
    const auto wrong = [](const int count) { return count != 1; };
    const std::vector<int>::const_iterator first_wrong = std::find_if(writes.cbegin(), writes.cend(), wrong);
    if (first_wrong != writes.cend())
    {
        const std::ptrdiff_t n_wrong = std::count_if(writes.cbegin(), writes.cend(), wrong);
        TBOX_ERROR(context << ":\n"
                           << "  the owned DOFs do not partition the DOFs of this rank: " << n_wrong
                           << (n_wrong == 1 ? " DOF is" : " DOFs are") << " not owned exactly once; the lowest is "
                           << n_lo + (first_wrong - writes.cbegin()) << ".\n");
    }
}

// The names of the built-in subdomain solvers.
const std::vector<std::string>&
built_in_subdomain_solver_names()
{
    static const std::vector<std::string> names = { "petsc", "blas-lapack" };
    return names;
}

// The factory of the built-in subdomain solver with the given name, or an empty function if the name is not that
// of one.
PETScLevelSolver::SubdomainSolverFactory
built_in_subdomain_solver_factory(const std::string& subdomain_solver_type)
{
    if (equals_ignore_case(subdomain_solver_type, "petsc"))
    {
        return make_petsc_subdomain_solver;
    }
    if (equals_ignore_case(subdomain_solver_type, "blas-lapack"))
    {
        return make_blas_lapack_subdomain_solver;
    }
    return {};
}

// Report a factory that cannot be selected by name: one with an empty name or function, or with the name of a
// built-in subdomain solver or of an earlier factory.
void
check_subdomain_solver_factories(const std::string& object_name,
                                 const PETScLevelSolver::SubdomainSolverFactories& factories)
{
    for (std::size_t i = 0; i < factories.size(); ++i)
    {
        const std::string& name = factories[i].first;
        if (name.empty() || !factories[i].second)
        {
            TBOX_ERROR(object_name << "::init():\n"
                                   << "  a subdomain solver factory needs a nonempty name and a nonempty function.\n");
        }
        const auto same_name = [&](const std::string& other) { return equals_ignore_case(other, name); };
        if (std::any_of(built_in_subdomain_solver_names().begin(), built_in_subdomain_solver_names().end(), same_name))
        {
            TBOX_ERROR(object_name << "::init():\n"
                                   << "  the subdomain solver factory " << name
                                   << " has the name of a built-in subdomain solver.\n");
        }
        for (std::size_t j = 0; j < i; ++j)
        {
            if (same_name(factories[j].first))
            {
                TBOX_ERROR(object_name << "::init():\n"
                                       << "  the subdomain solver factory " << name
                                       << " is supplied more than once.\n");
            }
        }
    }
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
    PCType pc_type = d_selected_pc_type.c_str();
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
        const bool multiplicative = d_subdomain_composition == SubdomainComposition::MULTIPLICATIVE;
        if (!d_subdomain_composition || !d_subdomain_output || (multiplicative && !d_subdomain_grouping))
        {
            TBOX_ERROR(d_object_name << "::initializeSolverState():\n"
                                     << "  subdomain_relaxation needs composition and output, and grouping for "
                                        "MULTIPLICATIVE composition, when the preconditioner is a shell.\n");
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

        // Get the local submatrices. Multiplicative relaxation also needs the rows of each subdomain with every
        // column, from which the submatrices follow without another collective extraction.
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
            // d_sub_mat is populated one matrix at a time, but MatDestroySubMatrices() still expects an array
            // with one extra trailing slot, which it reads to decide whether a type-specific reuse context is
            // present. PetscCalloc1() zero-initializes that slot, so that MatDestroySubMatrices() destroys each
            // matrix individually instead of reading uninitialized memory.
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
            ierr = ISRestoreIndices(d_overlap_is[i], &overlap_indices);
            IBTK_CHKERRQ(ierr);
            d_subdomain_offsets[i + 1] = d_subdomain_offsets[i] + overlap_size;
        }
        if (d_check_subdomain_coverage)
        {
            check_dof_coverage(
                d_object_name + "::initializeSolverState()", d_overlap_is, n_hi - n_lo, DOFCoverage::AT_LEAST_ONCE);
        }
        const PetscInt n_packed = d_subdomain_offsets[d_n_local_subdomains];
        ierr = VecCreateSeq(PETSC_COMM_SELF, n_packed, &d_subdomain_rhs);
        IBTK_CHKERRQ(ierr);
        ierr = VecDuplicate(d_subdomain_rhs, &d_subdomain_solution);
        IBTK_CHKERRQ(ierr);
        // The gather needs communication only if a subdomain of some rank includes DOFs that another rank
        // owns. Every rank makes the same choice, so the collective scatter is created and applied by all
        // ranks or by none. FULL output adds the solutions back through the same scatter.
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
            // The corrections are stored at the positions of their DOFs in the sorted list of the distinct DOFs of
            // the subdomains of this rank.
            std::vector<PetscInt> dofs = gathered_indices;
            std::sort(dofs.begin(), dofs.end());
            dofs.erase(std::unique(dofs.begin(), dofs.end()), dofs.end());
            initializeSubdomainResiduals(subdomain_rows, dofs);
            ierr = MatDestroySubMatrices(d_n_local_subdomains, &subdomain_rows);
            IBTK_CHKERRQ(ierr);
            initializeSubdomainGroups(gathered_indices, dofs, n_lo, n_hi);
        }
        else if (d_subdomain_output == SubdomainOutput::OWNED)
        {
            initializeOwnedSubdomainOutput(n_lo, n_hi);
        }

        // Set up the subdomain solvers.
        if (!d_subdomain_solver)
        {
            d_subdomain_solver = d_subdomain_solver_factory(d_subdomain_solver_db);
            if (!*d_subdomain_solver)
            {
                TBOX_ERROR(d_object_name << "::initializeSolverState():\n"
                                         << "  the factory of the subdomain solver " << d_subdomain_solver_type
                                         << " returned an empty subdomain solver.\n");
            }
        }
        d_subdomain_solver->initializeSolverState(
            std::vector<Mat>(d_sub_mat, d_sub_mat + d_n_local_subdomains), d_overlap_is, d_options_prefix);
        d_subdomain_solver_initialized = true;
        ierr = PCSetType(ksp_pc, PCSHELL);
        IBTK_CHKERRQ(ierr);
        ierr = PCShellSetContext(ksp_pc, static_cast<void*>(this));
        IBTK_CHKERRQ(ierr);
        ierr = PCShellSetApply(
            ksp_pc, multiplicative ? PETScLevelSolver::PCApply_Multiplicative : PETScLevelSolver::PCApply_Additive);
        IBTK_CHKERRQ(ierr);
        ierr = PCShellSetName(ksp_pc, "subdomain_relaxation");
        IBTK_CHKERRQ(ierr);
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

    // Deallocate PETSc objects for subdomain relaxation. The views share the arrays of the packed vectors, so they
    // are destroyed first.
    if (d_pc_type == "shell")
    {
        for (std::vector<Vec>* vectors : { &d_subdomain_rhs_views, &d_subdomain_residual_views, &d_halo_vectors })
        {
            for (Vec& v : *vectors)
            {
                ierr = VecDestroy(&v);
                IBTK_CHKERRQ(ierr);
            }
            vectors->clear();
        }
        for (Mat& residual_matrix : d_residual_matrices)
        {
            ierr = MatDestroy(&residual_matrix);
            IBTK_CHKERRQ(ierr);
        }
        d_residual_matrices.clear();
        ierr = VecDestroy(&d_subdomain_residual);
        IBTK_CHKERRQ(ierr);
        ierr = VecScatterDestroy(&d_restriction);
        IBTK_CHKERRQ(ierr);
        d_gather_indices.clear();
        ierr = VecDestroy(&d_subdomain_rhs);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&d_subdomain_solution);
        IBTK_CHKERRQ(ierr);
        d_subdomain_offsets.clear();
        d_write_offsets.clear();
        d_write_sources.clear();
        d_write_targets.clear();
        for (std::vector<PetscInt>* indices : { &d_packed_positions,
                                                &d_position_slots,
                                                &d_halo_offsets,
                                                &d_halo_positions,
                                                &d_group_support_offsets,
                                                &d_group_support,
                                                &d_group_output_offsets,
                                                &d_group_output_sources,
                                                &d_group_output_targets })
        {
            indices->clear();
        }
        d_group_visit_offsets.clear();
        d_group_visits.clear();
        d_group_correction.clear();
        d_full_correction.clear();
        ierr = MatDestroySubMatrices(d_n_local_subdomains, &d_sub_mat);
        IBTK_CHKERRQ(ierr);
        d_sub_mat = nullptr;
        d_n_local_subdomains = 0;
    }

    // Destroy the index sets that initializeSolverState() created.
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
                       const SubdomainSolverFactories& subdomain_solver_factories)
{
    check_subdomain_solver_factories(d_object_name, subdomain_solver_factories);
    d_options_prefix = default_options_prefix;
    if (input_db)
    {
        if (input_db->keyExists("options_prefix")) d_options_prefix = input_db->getString("options_prefix");
        if (input_db->keyExists("enable_logging")) d_enable_logging = input_db->getBool("enable_logging");
        if (input_db->keyExists("max_iterations")) d_max_iterations = input_db->getInteger("max_iterations");
        if (input_db->keyExists("abs_residual_tol")) d_abs_residual_tol = input_db->getDouble("abs_residual_tol");
        if (input_db->keyExists("rel_residual_tol")) d_rel_residual_tol = input_db->getDouble("rel_residual_tol");
        if (input_db->keyExists("ksp_type")) d_ksp_type = input_db->getString("ksp_type");
        if (input_db->keyExists("pc_type"))
        {
            d_selected_pc_type = input_db->getString("pc_type");
            d_pc_type = d_selected_pc_type;
        }
        if (input_db->keyExists("shell_pc_type"))
        {
            TBOX_ERROR(d_object_name << "::init():\n"
                                     << "  shell_pc_type has been removed; the subdomain_relaxation database sets "
                                        "the composition and output of the shell preconditioner.\n");
        }
        if (input_db->keyExists("subdomain_relaxation"))
        {
            Pointer<Database> relaxation_db = input_db->getDatabase("subdomain_relaxation");
            check_database_keys(
                d_object_name + "::init()", relaxation_db, { "composition", "grouping", "output", "subdomain_solver" });
            if (relaxation_db->keyExists("composition"))
            {
                const std::string composition = relaxation_db->getString("composition");
                if (equals_ignore_case(composition, "ADDITIVE"))
                {
                    d_subdomain_composition = SubdomainComposition::ADDITIVE;
                }
                else if (equals_ignore_case(composition, "MULTIPLICATIVE"))
                {
                    d_subdomain_composition = SubdomainComposition::MULTIPLICATIVE;
                }
                else
                {
                    TBOX_ERROR(d_object_name << "::init():\n"
                                             << "  unsupported subdomain_relaxation composition = " << composition
                                             << "; supported values are \"ADDITIVE\" and \"MULTIPLICATIVE\".\n");
                }
            }
            if (relaxation_db->keyExists("grouping"))
            {
                const std::string grouping = relaxation_db->getString("grouping");
                if (!equals_ignore_case(grouping, "RANK"))
                {
                    TBOX_ERROR(d_object_name << "::init():\n"
                                             << "  unsupported subdomain_relaxation grouping = " << grouping
                                             << "; the supported value is \"RANK\".\n");
                }
                if (d_subdomain_composition == SubdomainComposition::ADDITIVE)
                {
                    TBOX_ERROR(d_object_name << "::init():\n"
                                             << "  subdomain_relaxation grouping applies only to MULTIPLICATIVE "
                                                "composition.\n");
                }
                d_subdomain_grouping = SubdomainGrouping::RANK;
            }
            if (relaxation_db->keyExists("output"))
            {
                const std::string output = relaxation_db->getString("output");
                if (equals_ignore_case(output, "FULL"))
                {
                    d_subdomain_output = SubdomainOutput::FULL;
                }
                else if (equals_ignore_case(output, "OWNED"))
                {
                    d_subdomain_output = SubdomainOutput::OWNED;
                }
                else
                {
                    TBOX_ERROR(d_object_name << "::init():\n"
                                             << "  unsupported subdomain_relaxation output = " << output
                                             << "; supported values are \"FULL\" and \"OWNED\".\n");
                }
            }
            if (relaxation_db->keyExists("subdomain_solver"))
            {
                d_subdomain_solver_db = relaxation_db->getDatabase("subdomain_solver");
                d_subdomain_solver_type = d_subdomain_solver_db->getStringWithDefault("type", d_subdomain_solver_type);
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
    // The selected subdomain solver is created, and its settings are validated, only when a shell preconditioner
    // needs it.
    d_subdomain_solver_factory = built_in_subdomain_solver_factory(d_subdomain_solver_type);
    if (!d_subdomain_solver_factory)
    {
        const SubdomainSolverFactories::const_iterator factory =
            std::find_if(subdomain_solver_factories.begin(),
                         subdomain_solver_factories.end(),
                         [&](const auto& entry) { return equals_ignore_case(entry.first, d_subdomain_solver_type); });
        if (factory == subdomain_solver_factories.end())
        {
            std::vector<std::string> supported_names = built_in_subdomain_solver_names();
            for (const auto& entry : subdomain_solver_factories)
            {
                supported_names.push_back(entry.first);
            }
            std::string supported;
            for (std::size_t k = 0; k < supported_names.size(); ++k)
            {
                supported += (k == 0                          ? "" :
                              k + 1 == supported_names.size() ? ", and " :
                                                                ", ") +
                             std::string("\"") + supported_names[k] + "\"";
            }
            TBOX_ERROR(d_object_name << "::init():\n"
                                     << "  unsupported subdomain_solver type = " << d_subdomain_solver_type
                                     << "; supported values are " << supported << ".\n");
        }
        d_subdomain_solver_factory = factory->second;
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

void
PETScLevelSolver::initializeOwnedSubdomainOutput(const PetscInt n_lo, const PetscInt n_hi)
{
    if (static_cast<int>(d_nonoverlap_is.size()) != d_n_local_subdomains)
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState():\n"
                                 << "  OWNED output needs a nonoverlapping set for each subdomain.\n");
    }
    d_write_offsets.assign(d_n_local_subdomains + 1, 0);
    d_write_sources.clear();
    d_write_targets.clear();
    for (int i = 0; i < d_n_local_subdomains; ++i)
    {
        PetscInt overlap_size = 0, nonoverlap_size = 0;
        int ierr = ISGetLocalSize(d_overlap_is[i], &overlap_size);
        IBTK_CHKERRQ(ierr);
        ierr = ISGetLocalSize(d_nonoverlap_is[i], &nonoverlap_size);
        IBTK_CHKERRQ(ierr);
        const PetscInt *overlap_indices = nullptr, *nonoverlap_indices = nullptr;
        ierr = ISGetIndices(d_overlap_is[i], &overlap_indices);
        IBTK_CHKERRQ(ierr);
        ierr = ISGetIndices(d_nonoverlap_is[i], &nonoverlap_indices);
        IBTK_CHKERRQ(ierr);

        // The subdomain's nonoverlapping DOFs are owned by this rank. Since the IS are sorted, they
        // are found in one pass over the overlapping DOFs.
        PetscInt jj = 0;
        for (PetscInt ii = 0; ii < overlap_size && jj < nonoverlap_size; ++ii)
        {
            if (overlap_indices[ii] == nonoverlap_indices[jj])
            {
                if (nonoverlap_indices[jj] < n_lo || nonoverlap_indices[jj] >= n_hi)
                {
                    TBOX_ERROR(d_object_name << "::initializeSolverState():\n"
                                             << "  nonoverlapping DOF " << nonoverlap_indices[jj] << " of subdomain "
                                             << i << " is not owned by this rank.\n");
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
        d_write_offsets[i + 1] = static_cast<PetscInt>(d_write_sources.size());
        ierr = ISRestoreIndices(d_overlap_is[i], &overlap_indices);
        IBTK_CHKERRQ(ierr);
        ierr = ISRestoreIndices(d_nonoverlap_is[i], &nonoverlap_indices);
        IBTK_CHKERRQ(ierr);
    }
    check_owned_output(d_object_name + "::initializeSolverState()", d_write_targets, n_lo, n_hi);
} // initializeOwnedSubdomainOutput

void
PETScLevelSolver::initializeSubdomainResiduals(Mat* rows, const std::vector<PetscInt>& dofs)
{
    // The views of the packed right-hand sides and residuals of each subdomain share the storage of the packed
    // vectors.
    int ierr = VecDuplicate(d_subdomain_rhs, &d_subdomain_residual);
    IBTK_CHKERRQ(ierr);
    const Vec packed[2] = { d_subdomain_rhs, d_subdomain_residual };
    std::vector<Vec>* const views[2] = { &d_subdomain_rhs_views, &d_subdomain_residual_views };
    for (int k = 0; k < 2; ++k)
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

    // The residual matrix of each subdomain keeps the columns of its rows that have entries. The columns are
    // numbered globally, so that the rows of subdomains that include DOFs of other ranks are available.
    d_residual_matrices.assign(d_n_local_subdomains, nullptr);
    d_halo_vectors.assign(d_n_local_subdomains, nullptr);
    d_halo_offsets.assign(1, 0);
    d_halo_positions.clear();
    std::vector<PetscInt> columns;
    for (int i = 0; i < d_n_local_subdomains; ++i)
    {
        PetscInt n_rows = 0;
        ierr = MatGetSize(rows[i], &n_rows, nullptr);
        IBTK_CHKERRQ(ierr);
        columns.clear();
        for (PetscInt row = 0; row < n_rows; ++row)
        {
            PetscInt count = 0;
            const PetscInt* indices = nullptr;
            ierr = MatGetRow(rows[i], row, &count, &indices, nullptr);
            IBTK_CHKERRQ(ierr);
            columns.insert(columns.end(), indices, indices + count);
            ierr = MatRestoreRow(rows[i], row, &count, &indices, nullptr);
            IBTK_CHKERRQ(ierr);
        }
        std::sort(columns.begin(), columns.end());
        columns.erase(std::unique(columns.begin(), columns.end()), columns.end());
        const PetscInt n_columns = static_cast<PetscInt>(columns.size());
        d_residual_matrices[i] = extract_columns(rows[i], columns.data(), n_columns, -1.0);
        for (const PetscInt column : columns)
        {
            d_halo_positions.push_back(position_in(dofs, column));
        }
        d_halo_offsets.push_back(static_cast<PetscInt>(d_halo_positions.size()));
        ierr = VecCreateSeq(PETSC_COMM_SELF, n_columns, &d_halo_vectors[i]);
        IBTK_CHKERRQ(ierr);
    }
} // initializeSubdomainResiduals

void
PETScLevelSolver::initializeSubdomainGroups(const std::vector<PetscInt>& gathered_indices,
                                            const std::vector<PetscInt>& dofs,
                                            const PetscInt n_lo,
                                            const PetscInt n_hi)
{
    d_packed_positions.resize(gathered_indices.size());
    d_position_slots.assign(dofs.size(), -1);
    for (std::size_t k = 0; k < gathered_indices.size(); ++k)
    {
        const PetscInt position = position_in(dofs, gathered_indices[k]);
        d_packed_positions[k] = position;
        if (d_position_slots[position] < 0)
        {
            d_position_slots[position] = static_cast<PetscInt>(k);
        }
    }

    // Each rank forms one group of its subdomains, which it visits in order.
    std::vector<std::vector<int>> group_subdomains(1);
    for (int i = 0; i < d_n_local_subdomains; ++i)
    {
        group_subdomains[0].push_back(i);
    }

    const bool owned_output = d_subdomain_output == SubdomainOutput::OWNED;
    if (owned_output && static_cast<int>(d_nonoverlap_is.size()) != d_n_local_subdomains)
    {
        TBOX_ERROR(d_object_name << "::initializeSolverState():\n"
                                 << "  OWNED output needs a nonoverlapping set for each subdomain.\n");
    }
    d_group_visit_offsets.assign(1, 0);
    d_group_visits.clear();
    d_group_support_offsets.assign(1, 0);
    d_group_support.clear();
    d_group_output_offsets.assign(1, 0);
    d_group_output_sources.clear();
    d_group_output_targets.clear();
    // The positions in the support of the group being set up are marked with the number of the group plus one.
    std::vector<std::size_t> support_group(dofs.size(), 0);
    for (std::size_t g = 0; g < group_subdomains.size(); ++g)
    {
        for (const int i : group_subdomains[g])
        {
            d_group_visits.push_back(i);
            for (PetscInt k = d_subdomain_offsets[i]; k < d_subdomain_offsets[i + 1]; ++k)
            {
                const PetscInt position = d_packed_positions[k];
                if (support_group[position] != g + 1)
                {
                    support_group[position] = g + 1;
                    d_group_support.push_back(position);
                }
            }
        }
        d_group_visit_offsets.push_back(static_cast<int>(d_group_visits.size()));
        d_group_support_offsets.push_back(static_cast<PetscInt>(d_group_support.size()));

        // A group owns the nonoverlapping DOFs of its subdomains, which must be DOFs of this rank in its support.
        if (owned_output)
        {
            for (const int i : group_subdomains[g])
            {
                PetscInt nonoverlap_size = 0;
                const PetscInt* nonoverlap_indices = nullptr;
                int ierr = ISGetLocalSize(d_nonoverlap_is[i], &nonoverlap_size);
                IBTK_CHKERRQ(ierr);
                ierr = ISGetIndices(d_nonoverlap_is[i], &nonoverlap_indices);
                IBTK_CHKERRQ(ierr);
                for (PetscInt k = 0; k < nonoverlap_size; ++k)
                {
                    const PetscInt dof = nonoverlap_indices[k];
                    const PetscInt position = position_in(dofs, dof);
                    if (dof < n_lo || dof >= n_hi || position < 0 || support_group[position] != g + 1)
                    {
                        TBOX_ERROR(d_object_name << "::initializeSolverState():\n"
                                                 << "  DOF " << dof << " owned by group " << g
                                                 << " is not a DOF of this rank or is not in the subdomains of the "
                                                    "group.\n");
                    }
                    d_group_output_sources.push_back(position);
                    d_group_output_targets.push_back(dof - n_lo);
                }
                ierr = ISRestoreIndices(d_nonoverlap_is[i], &nonoverlap_indices);
                IBTK_CHKERRQ(ierr);
            }
        }
        d_group_output_offsets.push_back(static_cast<PetscInt>(d_group_output_sources.size()));
    }
    if (owned_output)
    {
        check_owned_output(d_object_name + "::initializeSolverState()", d_group_output_targets, n_lo, n_hi);
    }
    d_group_correction.assign(dofs.size(), 0.0);
    d_full_correction.assign(owned_output ? 0 : dofs.size(), 0.0);
} // initializeSubdomainGroups

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
PETScLevelSolver::addPackedSolutions(Vec y) const
{
    PetscFunctionBeginUser;
    int ierr = 0;
    if (d_restriction_communicates)
    {
        ierr = VecScatterBegin(d_restriction, d_subdomain_solution, y, ADD_VALUES, SCATTER_REVERSE);
        CHKERRQ(ierr);
        ierr = VecScatterEnd(d_restriction, d_subdomain_solution, y, ADD_VALUES, SCATTER_REVERSE);
        CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }
    // All of the DOFs of the subdomains are owned by their ranks, so add to the local array.
    const PetscScalar* solution = nullptr;
    PetscScalar* y_values = nullptr;
    ierr = VecGetArrayRead(d_subdomain_solution, &solution);
    CHKERRQ(ierr);
    ierr = VecGetArray(y, &y_values);
    CHKERRQ(ierr);
    for (std::size_t k = 0; k < d_gather_indices.size(); ++k)
    {
        y_values[d_gather_indices[k]] += solution[k];
    }
    ierr = VecRestoreArray(y, &y_values);
    CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(d_subdomain_solution, &solution);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
} // addPackedSolutions

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
    // FULL output adds the corrections to y, and a DOF that no subdomain covers gets a zero correction.
    ierr = VecZeroEntries(y);
    CHKERRQ(ierr);
    ierr = solver->gatherSubdomainRhs(x);
    CHKERRQ(ierr);
    solver->d_subdomain_solver->solve(
        0, solver->d_n_local_subdomains, solver->d_subdomain_rhs, solver->d_subdomain_solution);
    if (solver->d_subdomain_output == SubdomainOutput::OWNED)
    {
        ierr = solver->writeSubdomainSolutions(0, solver->d_n_local_subdomains, y);
    }
    else
    {
        ierr = solver->addPackedSolutions(y);
    }
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
    std::vector<PetscScalar>& correction = solver->d_group_correction;
    std::vector<PetscScalar>& full_correction = solver->d_full_correction;
    std::fill(correction.begin(), correction.end(), 0.0);
    std::fill(full_correction.begin(), full_correction.end(), 0.0);
    const bool owned_output = solver->d_subdomain_output == SubdomainOutput::OWNED;
    PetscScalar* y_values = nullptr;
    if (owned_output)
    {
        ierr = VecGetArray(y, &y_values);
        CHKERRQ(ierr);
    }
    for (std::size_t g = 0; g + 1 < solver->d_group_visit_offsets.size(); ++g)
    {
        const int first_visit = solver->d_group_visit_offsets[g];
        for (int visit = first_visit; visit < solver->d_group_visit_offsets[g + 1]; ++visit)
        {
            const int i = solver->d_group_visits[visit];
            Vec rhs = solver->d_subdomain_rhs;
            // The correction of the group is zero at its first visit, so the residual there is the right-hand side.
            if (visit > first_visit)
            {
                PetscScalar* halo_values = nullptr;
                ierr = VecGetArray(solver->d_halo_vectors[i], &halo_values);
                CHKERRQ(ierr);
                for (PetscInt k = solver->d_halo_offsets[i]; k < solver->d_halo_offsets[i + 1]; ++k)
                {
                    const PetscInt position = solver->d_halo_positions[k];
                    halo_values[k - solver->d_halo_offsets[i]] = position >= 0 ? correction[position] : 0.0;
                }
                ierr = VecRestoreArray(solver->d_halo_vectors[i], &halo_values);
                CHKERRQ(ierr);
                ierr = MatMultAdd(solver->d_residual_matrices[i],
                                  solver->d_halo_vectors[i],
                                  solver->d_subdomain_rhs_views[i],
                                  solver->d_subdomain_residual_views[i]);
                CHKERRQ(ierr);
                rhs = solver->d_subdomain_residual;
            }
            solver->d_subdomain_solver->solve(i, i + 1, rhs, solver->d_subdomain_solution);
            const PetscScalar* solution = nullptr;
            ierr = VecGetArrayRead(solver->d_subdomain_solution, &solution);
            CHKERRQ(ierr);
            for (PetscInt k = solver->d_subdomain_offsets[i]; k < solver->d_subdomain_offsets[i + 1]; ++k)
            {
                correction[solver->d_packed_positions[k]] += solution[k];
            }
            ierr = VecRestoreArrayRead(solver->d_subdomain_solution, &solution);
            CHKERRQ(ierr);
        }

        // Only the completed correction of the group contributes to the output. It is then reset for the next
        // group.
        if (owned_output)
        {
            for (PetscInt k = solver->d_group_output_offsets[g]; k < solver->d_group_output_offsets[g + 1]; ++k)
            {
                y_values[solver->d_group_output_targets[k]] = correction[solver->d_group_output_sources[k]];
            }
        }
        for (PetscInt k = solver->d_group_support_offsets[g]; k < solver->d_group_support_offsets[g + 1]; ++k)
        {
            const PetscInt position = solver->d_group_support[k];
            if (!owned_output)
            {
                full_correction[position] += correction[position];
            }
            correction[position] = 0.0;
        }
    }
    if (owned_output)
    {
        ierr = VecRestoreArray(y, &y_values);
        CHKERRQ(ierr);
    }
    else
    {
        // Add the sum of the corrections to y through the scatter of the packed vectors, from the first entry of
        // each DOF.
        PetscScalar* solution = nullptr;
        ierr = VecGetArray(solver->d_subdomain_solution, &solution);
        CHKERRQ(ierr);
        std::fill(solution, solution + solver->d_subdomain_offsets.back(), 0.0);
        for (std::size_t position = 0; position < full_correction.size(); ++position)
        {
            solution[solver->d_position_slots[position]] = full_correction[position];
        }
        ierr = VecRestoreArray(solver->d_subdomain_solution, &solution);
        CHKERRQ(ierr);
        ierr = solver->addPackedSolutions(y);
        CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
} // PCApply_Multiplicative

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
