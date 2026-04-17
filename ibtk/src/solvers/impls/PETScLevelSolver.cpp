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
#include <ibtk/PETScLevelSolver.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/private/PETScLevelSolverBlasLapackShellBackend.h>
#include <ibtk/private/PETScLevelSolverEigenPseudoinverseShellBackend.h>
#include <ibtk/private/PETScLevelSolverEigenReferenceShellBackend.h>
#include <ibtk/private/PETScLevelSolverEigenShellBackend.h>
#include <ibtk/private/PETScLevelSolverPetscShellBackend.h>

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
#include <sstream>
#include <set>
#include <string>
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

std::unique_ptr<PETScLevelSolverShellBackend>
allocate_petsc_shell_backend(PETScLevelSolver& solver, Pointer<Database> input_db)
{
    auto backend = std::make_unique<PETScLevelSolverPetscShellBackend>(solver);
    backend->configure(input_db);
    return backend;
}

std::unique_ptr<PETScLevelSolverShellBackend>
allocate_blas_lapack_shell_backend(PETScLevelSolver& solver, Pointer<Database> input_db)
{
    auto backend = std::make_unique<PETScLevelSolverBlasLapackShellBackend>(solver);
    backend->configure(input_db);
    return backend;
}

std::unique_ptr<PETScLevelSolverShellBackend>
allocate_eigen_shell_backend(PETScLevelSolver& solver, Pointer<Database> input_db)
{
    auto backend = std::make_unique<PETScLevelSolverEigenShellBackend>(solver);
    backend->configure(input_db);
    return backend;
}

std::unique_ptr<PETScLevelSolverShellBackend>
allocate_eigen_pseudoinverse_shell_backend(PETScLevelSolver& solver, Pointer<Database> input_db)
{
    auto backend = std::make_unique<PETScLevelSolverEigenPseudoinverseShellBackend>(solver);
    backend->configure(input_db);
    return backend;
}

std::unique_ptr<PETScLevelSolverShellBackend>
allocate_eigen_reference_shell_backend(PETScLevelSolver& solver, Pointer<Database> input_db)
{
    auto backend = std::make_unique<PETScLevelSolverEigenReferenceShellBackend>(solver);
    backend->configure(input_db);
    return backend;
}
} // namespace

PETScLevelSolverShellBackendManager* PETScLevelSolverShellBackendManager::s_shell_backend_manager_instance = nullptr;
bool PETScLevelSolverShellBackendManager::s_registered_callback = false;
unsigned char PETScLevelSolverShellBackendManager::s_shutdown_priority = 200;

PETScLevelSolverShellBackendManager*
PETScLevelSolverShellBackendManager::getManager()
{
    if (!s_shell_backend_manager_instance)
    {
        s_shell_backend_manager_instance = new PETScLevelSolverShellBackendManager();
    }
    if (!s_registered_callback)
    {
        ShutdownRegistry::registerShutdownRoutine(freeManager, s_shutdown_priority);
        s_registered_callback = true;
    }
    return s_shell_backend_manager_instance;
}

void
PETScLevelSolverShellBackendManager::freeManager()
{
    delete s_shell_backend_manager_instance;
    s_shell_backend_manager_instance = nullptr;
}

std::unique_ptr<PETScLevelSolverShellBackend>
PETScLevelSolverShellBackendManager::allocateShellBackend(const std::string& type_key,
                                                          PETScLevelSolver& solver,
                                                          Pointer<Database> input_db) const
{
    const auto it = d_shell_backend_maker_map.find(type_key);
    if (it == d_shell_backend_maker_map.end())
    {
        TBOX_ERROR("PETScLevelSolverShellBackendManager::allocateShellBackend():\n"
                   << "  unrecognized shell backend type: " << type_key << "\n");
    }
    return (it->second)(solver, input_db);
}

void
PETScLevelSolverShellBackendManager::registerShellBackendFactoryFunction(const std::string& type_key,
                                                                         ShellBackendMaker backend_maker)
{
    if (d_shell_backend_maker_map.find(type_key) != d_shell_backend_maker_map.end())
    {
        pout << "PETScLevelSolverShellBackendManager::registerShellBackendFactoryFunction():\n"
             << "  NOTICE: overriding initialization function for shell backend type = " << type_key << "\n";
    }
    d_shell_backend_maker_map[type_key] = backend_maker;
}

std::vector<std::string>
PETScLevelSolverShellBackendManager::getRegisteredShellBackendTypes() const
{
    std::vector<std::string> type_keys;
    type_keys.reserve(d_shell_backend_maker_map.size());
    for (const auto& entry : d_shell_backend_maker_map)
    {
        type_keys.push_back(entry.first);
    }
    return type_keys;
}

PETScLevelSolverShellBackendManager::PETScLevelSolverShellBackendManager() : d_shell_backend_maker_map()
{
    registerShellBackendFactoryFunction("petsc", allocate_petsc_shell_backend);
    registerShellBackendFactoryFunction("blas-lapack", allocate_blas_lapack_shell_backend);
    registerShellBackendFactoryFunction("eigen", allocate_eigen_shell_backend);
    registerShellBackendFactoryFunction("eigen-pseudoinverse", allocate_eigen_pseudoinverse_shell_backend);
    registerShellBackendFactoryFunction("eigen-reference", allocate_eigen_reference_shell_backend);
}

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

const std::string&
PETScLevelSolver::getObjectNameForBackend() const
{
    return d_object_name;
}

const std::string&
PETScLevelSolver::getOptionsPrefixForBackend() const
{
    return d_options_prefix;
}

bool
PETScLevelSolver::isShellMultiplicativeForBackend() const
{
    return d_shell_smoother_composition == ShellSmootherComposition::MULTIPLICATIVE;
}

bool
PETScLevelSolver::useRestrictPartitionForBackend() const
{
    return d_shell_smoother_partition == ShellSmootherPartition::RESTRICT;
}

const std::vector<std::vector<int>>&
PETScLevelSolver::getSubdomainDOFsForBackend() const
{
    return d_subdomain_dofs;
}

const std::vector<std::vector<int>>&
PETScLevelSolver::getNonoverlapSubdomainDOFsForBackend() const
{
    return d_nonoverlap_subdomain_dofs;
}

Mat
PETScLevelSolver::getPETScMatForBackend() const
{
    return d_petsc_mat;
}

Vec
PETScLevelSolver::getPETScXForBackend() const
{
    return d_petsc_x;
}

Vec
PETScLevelSolver::getPETScBForBackend() const
{
    return d_petsc_b;
}

void
PETScLevelSolver::postprocessShellResultForBackend(Vec y)
{
    postprocessShellResult(y);
}

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
        if (input_db->keyExists("initial_guess_nonzero"))
            d_initial_guess_nonzero = input_db->getBool("initial_guess_nonzero");
    }
    loadShellBackends(input_db);
    d_preconditioner_type = parsePreconditionerType(d_pc_type);
    configureShellSmootherType();
    return;
} // init

void
PETScLevelSolver::loadShellBackends(Pointer<Database> input_db)
{
    d_shell_backends.clear();
    auto* manager = PETScLevelSolverShellBackendManager::getManager();
    for (const auto& type_key : manager->getRegisteredShellBackendTypes())
    {
        auto backend = manager->allocateShellBackend(type_key, *this, input_db);
        if (!backend) continue;
        const std::string key = backend->getTypeKey();
        d_shell_backends[key] = std::move(backend);
    }
}

void
PETScLevelSolver::configureShellSmootherType()
{
    if (d_preconditioner_type != PreconditionerType::SHELL)
    {
        d_shell_pc_type.clear();
        d_shell_smoother_backend_key = "petsc";
        d_shell_smoother_composition = ShellSmootherComposition::MULTIPLICATIVE;
        d_shell_smoother_partition = ShellSmootherPartition::BASIC;
        d_active_shell_backend = nullptr;
        return;
    }

    d_shell_pc_type = normalizeShellSmootherType(d_shell_pc_type);
    if (d_shell_pc_type.empty())
    {
        d_shell_smoother_backend_key = "petsc";
        d_shell_smoother_composition = ShellSmootherComposition::MULTIPLICATIVE;
        d_shell_smoother_partition = ShellSmootherPartition::BASIC;
        d_active_shell_backend = getShellBackend(d_shell_smoother_backend_key);
        return;
    }

    const std::string type_key = extractShellSmootherTypeKey(d_shell_pc_type);
    d_shell_smoother_backend_key = parseShellSmootherBackendKey(type_key);
    d_shell_smoother_composition = parseShellSmootherComposition(d_shell_pc_type);
    d_shell_smoother_partition = parseShellSmootherPartition(d_shell_pc_type, d_shell_smoother_composition);
    d_active_shell_backend = getShellBackend(d_shell_smoother_backend_key);
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

std::string
PETScLevelSolver::extractShellSmootherTypeKey(const std::string& type) const
{
    std::string key = type;
    if (key.rfind("additive", 0) == 0)
    {
        key.erase(0, std::strlen("additive"));
    }
    else if (key.rfind("multiplicative", 0) == 0)
    {
        key.erase(0, std::strlen("multiplicative"));
    }
    if (!key.empty() && key.front() == '-') key.erase(0, 1);
    if (key.size() >= std::strlen("-restrict") &&
        key.compare(key.size() - std::strlen("-restrict"), std::strlen("-restrict"), "-restrict") == 0)
    {
        key.erase(key.size() - std::strlen("-restrict"));
    }
    else if (key.size() >= std::strlen("-basic") &&
             key.compare(key.size() - std::strlen("-basic"), std::strlen("-basic"), "-basic") == 0)
    {
        key.erase(key.size() - std::strlen("-basic"));
    }
    return key;
}

PETScLevelSolver::PreconditionerType
PETScLevelSolver::parsePreconditionerType(const std::string& type) const
{
    if (type == "asm") return PreconditionerType::ASM;
    if (type == "fieldsplit") return PreconditionerType::FIELDSPLIT;
    if (type == "shell") return PreconditionerType::SHELL;
    return PreconditionerType::OTHER;
} // parsePreconditionerType

std::string
PETScLevelSolver::parseShellSmootherBackendKey(const std::string& type_key) const
{
    const std::string normalized_key = type_key.empty() ? "petsc" : type_key;
    if (getShellBackend(normalized_key)) return normalized_key;

    std::ostringstream supported_types;
    bool first = true;
    for (const auto& entry : d_shell_backends)
    {
        if (!first) supported_types << ", ";
        supported_types << entry.first;
        first = false;
    }
    TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::parseShellSmootherBackendKey()\n"
                             << "Unknown shell smoother backend key: " << normalized_key << "\n"
                             << "Available shell backend keys: " << supported_types.str() << std::endl);
    return "petsc";
} // parseShellSmootherBackendKey

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
    std::vector<IS> subdomain_is, nonoverlap_subdomain_is;
    cacheGeneratedASMSubdomains();
    build_petsc_subdomain_index_sets(
        subdomain_is, nonoverlap_subdomain_is, d_subdomain_dofs, d_nonoverlap_subdomain_dofs);
    const int num_subdomains = static_cast<int>(subdomain_is.size());
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
        ierr = PCASMSetLocalSubdomains(ksp_pc, num_subdomains, subdomain_is.data(), nonoverlap_subdomain_is.data());
        IBTK_CHKERRQ(ierr);
    }
    destroy_index_sets(nonoverlap_subdomain_is);
    destroy_index_sets(subdomain_is);
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
PETScLevelSolver::deallocateShellData()
{
    for (auto& entry : d_shell_backends)
    {
        entry.second->deallocate();
    }
    d_active_shell_backend = nullptr;
    d_subdomain_dofs.clear();
    d_nonoverlap_subdomain_dofs.clear();
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
    switch (d_shell_smoother_composition)
    {
    case ShellSmootherComposition::ADDITIVE:
        apply_op = PETScLevelSolver::PCApply_AdditiveShell;
        pc_name_suffix = d_active_shell_backend ? d_active_shell_backend->getPCNameSuffixAdditive() : nullptr;
        break;
    case ShellSmootherComposition::MULTIPLICATIVE:
        apply_op = PETScLevelSolver::PCApply_MultiplicativeShell;
        pc_name_suffix = d_active_shell_backend ? d_active_shell_backend->getPCNameSuffixMultiplicative() : nullptr;
        break;
    }
    if (!d_active_shell_backend || !apply_op || !pc_name_suffix)
    {
        TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::configureShellApply()\n"
                                 << "Invalid shell smoother configuration for shell_pc_type = " << d_shell_pc_type
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
    cacheGeneratedASMSubdomains();
    d_active_shell_backend = getShellBackend(d_shell_smoother_backend_key);
    if (!d_active_shell_backend)
    {
        TBOX_ERROR(d_object_name << " " << d_options_prefix << " PETScLevelSolver::configureShellPreconditioner()\n"
                                 << "Selected shell backend key has no allocated backend: "
                                 << d_shell_smoother_backend_key << std::endl);
    }
    d_active_shell_backend->initialize();

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

PETScLevelSolverShellBackend*
PETScLevelSolver::getShellBackend(const std::string& type_key)
{
    auto it = d_shell_backends.find(type_key);
    return it == d_shell_backends.end() ? nullptr : it->second.get();
}

const PETScLevelSolverShellBackend*
PETScLevelSolver::getShellBackend(const std::string& type_key) const
{
    auto it = d_shell_backends.find(type_key);
    return it == d_shell_backends.end() ? nullptr : it->second.get();
}

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
PETScLevelSolver::PCApply_AdditiveShell(PC pc, Vec x, Vec y)
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
    if (!solver->d_active_shell_backend)
    {
        TBOX_ERROR("PETScLevelSolver::PCApply_AdditiveShell(): no active shell backend." << std::endl);
    }
    solver->d_active_shell_backend->applyAdditive(x, y);
    PetscFunctionReturn(0);
} // PCApply_AdditiveShell

PetscErrorCode
PETScLevelSolver::PCApply_MultiplicativeShell(PC pc, Vec x, Vec y)
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
    if (!solver->d_active_shell_backend)
    {
        TBOX_ERROR("PETScLevelSolver::PCApply_MultiplicativeShell(): no active shell backend." << std::endl);
    }
    solver->d_active_shell_backend->applyMultiplicative(x, y);
    PetscFunctionReturn(0);
} // PCApply_MultiplicativeShell

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
