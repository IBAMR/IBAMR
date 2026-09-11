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

#include <ibamr/StaggeredStokesPETScLevelSolver.h>
#include <ibamr/StaggeredStokesPETScMatUtilities.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>
#include <ibamr/private/CouplingAwareASMSubdomains.h>
#include <ibamr/private/StaggeredStokesEigenSchurComplementShellBackend.h>

#include <ibtk/GeneralSolver.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LinearSolver.h>
#include <ibtk/PETScLevelSolver.h>
#include <ibtk/PoissonUtilities.h>

#include <tbox/Array.h>
#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <petscao.h>
#include <petscis.h>
#include <petsclog.h>
#include <petscvec.h>

#include <BoundaryBox.h>
#include <BoxList.h>
#include <CellData.h>
#include <CellVariable.h>
#include <CoarseFineBoundary.h>
#include <IntVector.h>
#include <MultiblockDataTranslator.h>
#include <Patch.h>
#include <PatchGeometry.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RefineSchedule.h>
#include <SAMRAIVectorReal.h>
#include <SideData.h>
#include <SideVariable.h>
#include <Variable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;
static const int SIDEG = 1;
static const int NOGHOST = 0;

/*! \brief Determine whether the level's boxes cover its physical domain. */
bool
level_covers_entire_physical_domain(const Pointer<PatchLevel<NDIM>>& level)
{
    BoxList<NDIM> uncovered_domain(level->getPhysicalDomain());
    uncovered_domain.removeIntersections(BoxList<NDIM>(level->getBoxes()));
    return uncovered_domain.isEmpty();
}

/*! \brief Construct the named field's local coupled DOF indices. */
void
construct_field_is(const std::vector<std::set<int>>& field_is,
                   const std::vector<std::string>& field_names,
                   const std::string& field_name,
                   IS& local_is)
{
    const auto field_name_it = std::find(field_names.begin(), field_names.end(), field_name);
    if (field_name_it == field_names.end())
    {
        TBOX_ERROR("construct_field_is():\n"
                   << "  unable to locate " << field_name << " field DOFs.\n");
    }

    const std::size_t field_idx = static_cast<std::size_t>(std::distance(field_names.begin(), field_name_it));
    std::vector<PetscInt> field_dofs(field_is[field_idx].begin(), field_is[field_idx].end());
    int ierr = ISCreateGeneral(PETSC_COMM_WORLD,
                               static_cast<PetscInt>(field_dofs.size()),
                               field_dofs.empty() ? nullptr : field_dofs.data(),
                               PETSC_COPY_VALUES,
                               &local_is);
    IBTK_CHKERRQ(ierr);

    return;
}

/*! \brief Construct the mapping from compact velocity indices to coupled indices. */
void
construct_velocity_field_ao(IS velocity_field_is_local, Mat velocity_block_mat, AO& velocity_field_ao)
{
    PetscInt n_velocity_local = 0;
    PetscInt row_start = 0;
    PetscInt row_end = 0;
    int ierr = ISGetLocalSize(velocity_field_is_local, &n_velocity_local);
    IBTK_CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(velocity_block_mat, &row_start, &row_end);
    IBTK_CHKERRQ(ierr);
    if (n_velocity_local != row_end - row_start)
    {
        TBOX_ERROR("construct_velocity_field_ao():\n"
                   << "  local velocity-field DOF count (" << n_velocity_local
                   << ") does not match local velocity-block row count (" << row_end - row_start << ").\n");
    }

    std::vector<PetscInt> velocity_field_ids(static_cast<std::size_t>(n_velocity_local));
    for (PetscInt k = 0; k < n_velocity_local; ++k)
    {
        velocity_field_ids[static_cast<std::size_t>(k)] = row_start + k;
    }
    const PetscInt* velocity_full_ids = nullptr;
    ierr = ISGetIndices(velocity_field_is_local, &velocity_full_ids);
    IBTK_CHKERRQ(ierr);
    ierr = AOCreateBasic(
        PETSC_COMM_WORLD, n_velocity_local, velocity_field_ids.data(), velocity_full_ids, &velocity_field_ao);
    IBTK_CHKERRQ(ierr);
    ierr = ISRestoreIndices(velocity_field_is_local, &velocity_full_ids);
    IBTK_CHKERRQ(ierr);

    return;
}

/*! \brief Insert compact velocity rows into coupled numbering. */
void
insert_velocity_block_rows(Mat source, AO mapping, Mat destination)
{
    PetscInt row_start = 0, row_end = 0;
    int ierr = MatGetOwnershipRange(source, &row_start, &row_end);
    IBTK_CHKERRQ(ierr);
    std::vector<PetscInt> mapped_cols;
    for (PetscInt row = row_start; row < row_end; ++row)
    {
        PetscInt ncols = 0;
        const PetscInt* cols = nullptr;
        const PetscScalar* vals = nullptr;
        ierr = MatGetRow(source, row, &ncols, &cols, &vals);
        IBTK_CHKERRQ(ierr);

        mapped_cols.resize(static_cast<std::size_t>(ncols));
        for (PetscInt k = 0; k < ncols; ++k)
        {
            mapped_cols[static_cast<std::size_t>(k)] = cols[k];
        }
        ierr = AOApplicationToPetsc(mapping, ncols, mapped_cols.data());
        IBTK_CHKERRQ(ierr);

        PetscInt full_row = row;
        ierr = AOApplicationToPetsc(mapping, 1, &full_row);
        IBTK_CHKERRQ(ierr);

        ierr = MatSetValues(destination, 1, &full_row, ncols, mapped_cols.data(), vals, INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
        ierr = MatRestoreRow(source, row, &ncols, &cols, &vals);
        IBTK_CHKERRQ(ierr);
    }
}

std::unique_ptr<IBTK::PETScLevelSolverShellBackend>
allocate_eigen_schur_backend(Pointer<Database> input_db)
{
    return std::make_unique<StaggeredStokesEigenSchurComplementShellBackend>(input_db);
}
const bool registered_eigen_schur = []()
{
    IBTK::PETScLevelSolverShellBackendManager::get_manager().registerFactory("eigen-schur-complement",
                                                                             allocate_eigen_schur_backend);
    return true;
}();
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesPETScLevelSolver::StaggeredStokesPETScLevelSolver(const std::string& object_name,
                                                                 Pointer<Database> input_db,
                                                                 const std::string& default_options_prefix)
{
    GeneralSolver::init(object_name, /*homogeneous_bc*/ false);
    PETScLevelSolver::init(input_db, default_options_prefix);
    if (input_db)
    {
        d_asm_mode = string_to_enum<ASMSubdomainConstructionMode>(
            input_db->getStringWithDefault("asm_subdomain_construction_mode", "GEOMETRICAL"));
        d_ca_seed_axis = input_db->getIntegerWithDefault("coupling_aware_asm_seed_axis", d_ca_seed_axis);
        d_ca_seed_stride = input_db->getIntegerWithDefault("coupling_aware_asm_seed_stride", d_ca_seed_stride);
        d_ca_order = string_to_enum<CouplingAwareASMSeedTraversalOrder>(
            input_db->getStringWithDefault("coupling_aware_asm_seed_traversal_order", enum_to_string(d_ca_order)));
        d_ca_policy = string_to_enum<CouplingAwareASMClosurePolicy>(
            input_db->getStringWithDefault("coupling_aware_asm_closure_policy", "RELAXED"));
        d_ca_relative_zero_tol =
            input_db->getDoubleWithDefault("coupling_aware_asm_relative_zero_tol", d_ca_relative_zero_tol);
    }
#if (NDIM == 2)
    const bool valid_order =
        d_ca_order == CouplingAwareASMSeedTraversalOrder::I_J || d_ca_order == CouplingAwareASMSeedTraversalOrder::J_I;
#else
    const bool valid_order = d_ca_order == CouplingAwareASMSeedTraversalOrder::I_J_K ||
                             d_ca_order == CouplingAwareASMSeedTraversalOrder::J_K_I ||
                             d_ca_order == CouplingAwareASMSeedTraversalOrder::K_I_J;
#endif
    if (d_asm_mode == ASMSubdomainConstructionMode::UNKNOWN || d_ca_policy == CouplingAwareASMClosurePolicy::UNKNOWN ||
        !valid_order || d_ca_seed_axis < 0 || d_ca_seed_axis >= NDIM || d_ca_seed_stride < 1 ||
        !std::isfinite(d_ca_relative_zero_tol) || d_ca_relative_zero_tol < 0.0)
    {
        TBOX_ERROR(d_object_name << ": invalid coupling-aware ASM construction settings.\n");
    }
    // Construct the DOF index variable/context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(object_name + "::CONTEXT");
    d_u_dof_index_var = new SideVariable<NDIM, int>(object_name + "::u_dof_index");
    if (var_db->checkVariableExists(d_u_dof_index_var->getName()))
    {
        d_u_dof_index_var = var_db->getVariable(d_u_dof_index_var->getName());
        d_u_dof_index_idx = var_db->mapVariableAndContextToIndex(d_u_dof_index_var, d_context);
        var_db->removePatchDataIndex(d_u_dof_index_idx);
    }
    const int u_gcw = std::max(d_overlap_size.max(), SIDEG);
    d_u_dof_index_idx = var_db->registerVariableAndContext(d_u_dof_index_var, d_context, u_gcw);
    d_p_dof_index_var = new CellVariable<NDIM, int>(object_name + "::p_dof_index");
    if (var_db->checkVariableExists(d_p_dof_index_var->getName()))
    {
        d_p_dof_index_var = var_db->getVariable(d_p_dof_index_var->getName());
        d_p_dof_index_idx = var_db->mapVariableAndContextToIndex(d_p_dof_index_var, d_context);
        var_db->removePatchDataIndex(d_p_dof_index_idx);
    }
    const int p_gcw = std::max(d_overlap_size.max(), CELLG);
    d_p_dof_index_idx = var_db->registerVariableAndContext(d_p_dof_index_var, d_context, p_gcw);

    // Construct the nullspace variable/index.
    d_u_nullspace_var = new SideVariable<NDIM, double>(object_name + "::u_nullspace_var");
    if (var_db->checkVariableExists(d_u_nullspace_var->getName()))
    {
        d_u_nullspace_var = var_db->getVariable(d_u_nullspace_var->getName());
        d_u_nullspace_idx = var_db->mapVariableAndContextToIndex(d_u_nullspace_var, d_context);
        var_db->removePatchDataIndex(d_u_nullspace_idx);
    }
    d_u_nullspace_idx = var_db->registerVariableAndContext(d_u_nullspace_var, d_context, NOGHOST);
    d_p_nullspace_var = new CellVariable<NDIM, double>(object_name + "::p_nullspace_var");
    if (var_db->checkVariableExists(d_p_nullspace_var->getName()))
    {
        d_p_nullspace_var = var_db->getVariable(d_p_nullspace_var->getName());
        d_p_nullspace_idx = var_db->mapVariableAndContextToIndex(d_p_nullspace_var, d_context);
        var_db->removePatchDataIndex(d_p_nullspace_idx);
    }
    d_p_nullspace_idx = var_db->registerVariableAndContext(d_p_nullspace_var, d_context, NOGHOST);

    return;
} // StaggeredStokesPETScLevelSolver

StaggeredStokesPETScLevelSolver::~StaggeredStokesPETScLevelSolver()
{
    if (d_is_initialized) deallocateSolverState();
    int ierr = MatDestroy(&d_operator_mat);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&d_augmented_operator_mat);
    IBTK_CHKERRQ(ierr);
    return;
} // ~StaggeredStokesPETScLevelSolver

void
StaggeredStokesPETScLevelSolver::setOperatorMat(Mat operator_mat)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setOperatorMat(): deallocate solver state before changing the matrix.");
    }
    // Retain before releasing so that same-handle replacement is safe.
    int ierr = PetscObjectReference(reinterpret_cast<PetscObject>(operator_mat));
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&d_operator_mat);
    IBTK_CHKERRQ(ierr);
    d_operator_mat = operator_mat;
    return;
} // setOperatorMat

void
StaggeredStokesPETScLevelSolver::setAugmentedOperatorMat(Mat augmented_operator_mat)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setAugmentedOperatorMat(): deallocate solver state before changing the matrix.");
    }
    int ierr = PetscObjectReference(reinterpret_cast<PetscObject>(augmented_operator_mat));
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&d_augmented_operator_mat);
    IBTK_CHKERRQ(ierr);
    d_augmented_operator_mat = augmented_operator_mat;
    return;
} // setAugmentedOperatorMat

/////////////////////////////// PROTECTED ////////////////////////////////////

void
StaggeredStokesPETScLevelSolver::initializeShellBackend(IBTK::PETScLevelSolverShellBackend& backend,
                                                        const bool use_multiplicative,
                                                        const IBTK::PETScLevelSolverShellTraversal traversal)
{
    StaggeredStokesEigenSchurComplementShellBackend* schur =
        dynamic_cast<StaggeredStokesEigenSchurComplementShellBackend*>(&backend);
    if (!schur)
    {
        PETScLevelSolver::initializeShellBackend(backend, use_multiplicative, traversal);
        return;
    }
    if (IBTK_MPI::getNodes() != 1)
    {
        TBOX_ERROR("Eigen shell backends require one MPI rank.\n");
    }
    std::vector<std::string> names;
    std::vector<std::set<int>> fields;
    generateFieldSplitSubdomains(names, fields);
    const std::vector<std::string>::const_iterator velocity = std::find(names.cbegin(), names.cend(), "velocity");
    const std::vector<std::string>::const_iterator pressure = std::find(names.cbegin(), names.cend(), "pressure");
    if (velocity == names.cend() || pressure == names.cend())
    {
        TBOX_ERROR("Eigen Schur shell backend requires named velocity and pressure fields.\n");
    }
    schur->initializeSolverState(d_petsc_mat,
                                 d_petsc_x,
                                 d_petsc_b,
                                 d_overlap_is,
                                 d_nonoverlap_is,
                                 fields[velocity - names.cbegin()],
                                 fields[pressure - names.cbegin()],
                                 d_options_prefix,
                                 use_multiplicative,
                                 traversal);
}

void
StaggeredStokesPETScLevelSolver::generateASMSubdomains(std::vector<std::set<int>>& overlap_is,
                                                       std::vector<std::set<int>>& nonoverlap_is)
{
    if (d_asm_mode == ASMSubdomainConstructionMode::COUPLING_AWARE)
    {
        if (!d_ca_subdomains)
        {
            d_ca_subdomains =
                std::make_unique<CouplingAwareASMSubdomains>(d_u_dof_index_idx, d_p_dof_index_idx, d_level);
        }
        d_ca_subdomains->constructSubdomains(overlap_is,
                                             nonoverlap_is,
                                             d_num_dofs_per_proc,
                                             d_petsc_mat,
                                             d_ca_seed_axis,
                                             d_ca_seed_stride,
                                             d_ca_order,
                                             d_ca_policy,
                                             d_ca_relative_zero_tol);
        return;
    }
    // Construct subdomains for ASM and MSM preconditioner.
    StaggeredStokesPETScMatUtilities::constructPatchLevelASMSubdomains(overlap_is,
                                                                       nonoverlap_is,
                                                                       d_box_size,
                                                                       d_overlap_size,
                                                                       d_num_dofs_per_proc,
                                                                       d_u_dof_index_idx,
                                                                       d_p_dof_index_idx,
                                                                       d_level,
                                                                       d_cf_boundary);

    return;
} // generateASMSubdomains

void
StaggeredStokesPETScLevelSolver::generateFieldSplitSubdomains(std::vector<std::string>& field_names,
                                                              std::vector<std::set<int>>& field_is)
{
    // Set IS'es for field split preconditioner.
    StaggeredStokesPETScMatUtilities::constructPatchLevelFields(
        field_is, field_names, d_num_dofs_per_proc, d_u_dof_index_idx, d_p_dof_index_idx, d_level);

    return;
} // generateFieldSplitSubdomains

void
StaggeredStokesPETScLevelSolver::initializeSolverStateSpecialized(const SAMRAIVectorReal<NDIM, double>& x,
                                                                  const SAMRAIVectorReal<NDIM, double>& /*b*/)
{
    // Allocate DOF index data.
    if (!d_level->checkAllocated(d_u_dof_index_idx)) d_level->allocatePatchData(d_u_dof_index_idx);
    if (!d_level->checkAllocated(d_p_dof_index_idx)) d_level->allocatePatchData(d_p_dof_index_idx);
    StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
        d_num_dofs_per_proc, d_u_dof_index_idx, d_p_dof_index_idx, d_level);

    // Setup PETSc objects.
    int ierr;
    const int mpi_rank = IBTK_MPI::getRank();
    ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[mpi_rank], PETSC_DETERMINE, &d_petsc_x);
    IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[mpi_rank], PETSC_DETERMINE, &d_petsc_b);
    IBTK_CHKERRQ(ierr);
    if (d_operator_mat)
    {
        if (d_augmented_operator_mat)
        {
            // Augmentation changes entries: preserve the installed input.
            ierr = MatDuplicate(d_operator_mat, MAT_COPY_VALUES, &d_petsc_mat);
            IBTK_CHKERRQ(ierr);
        }
        else
        {
            // Solver state aliases the separately retained input reference.
            d_petsc_mat = d_operator_mat;
        }
    }
    else
    {
        StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(d_petsc_mat,
                                                                         d_U_problem_coefs,
                                                                         d_U_bc_coefs,
                                                                         d_new_time,
                                                                         d_num_dofs_per_proc,
                                                                         d_u_dof_index_idx,
                                                                         d_p_dof_index_idx,
                                                                         d_level);
    }

    if (d_augmented_operator_mat)
    {
        PetscInt full_m = 0, full_n = 0, aug_m = 0, aug_n = 0;
        ierr = MatGetSize(d_petsc_mat, &full_m, &full_n);
        IBTK_CHKERRQ(ierr);
        ierr = MatGetSize(d_augmented_operator_mat, &aug_m, &aug_n);
        IBTK_CHKERRQ(ierr);

        if (aug_m == full_m && aug_n == full_n)
        {
            ierr = MatAXPY(d_petsc_mat, 1.0, d_augmented_operator_mat, DIFFERENT_NONZERO_PATTERN);
            IBTK_CHKERRQ(ierr);
        }
        else
        {
            std::vector<std::set<int>> field_is;
            std::vector<std::string> field_names;
            StaggeredStokesPETScMatUtilities::constructPatchLevelFields(
                field_is, field_names, d_num_dofs_per_proc, d_u_dof_index_idx, d_p_dof_index_idx, d_level);
            IS velocity_field_is_local = nullptr;
            AO velocity_field_ao = nullptr;
            construct_field_is(field_is, field_names, "velocity", velocity_field_is_local);
            PetscInt n_velocity_global = 0;
            ierr = ISGetSize(velocity_field_is_local, &n_velocity_global);
            IBTK_CHKERRQ(ierr);
            if (aug_m != n_velocity_global || aug_n != n_velocity_global)
            {
                TBOX_ERROR("StaggeredStokesPETScLevelSolver::initializeSolverStateSpecialized():\n"
                           << "  augmented operator has incompatible size: (" << aug_m << " x " << aug_n << ").\n"
                           << "  expected either full operator size (" << full_m << " x " << full_n
                           << ") or velocity block size (" << n_velocity_global << " x " << n_velocity_global
                           << ").\n");
            }

            construct_velocity_field_ao(velocity_field_is_local, d_augmented_operator_mat, velocity_field_ao);

            PetscInt full_m_local = 0, full_n_local = 0;
            ierr = MatGetLocalSize(d_petsc_mat, &full_m_local, &full_n_local);
            IBTK_CHKERRQ(ierr);

            Mat preallocator = nullptr;
            ierr = MatCreate(PETSC_COMM_WORLD, &preallocator);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetSizes(preallocator, full_m_local, full_n_local, full_m, full_n);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetType(preallocator, MATPREALLOCATOR);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetUp(preallocator);
            IBTK_CHKERRQ(ierr);

            insert_velocity_block_rows(d_augmented_operator_mat, velocity_field_ao, preallocator);
            ierr = MatAssemblyBegin(preallocator, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);
            ierr = MatAssemblyEnd(preallocator, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);

            Mat embedded_augmented_operator_mat = nullptr;
            ierr = MatCreate(PETSC_COMM_WORLD, &embedded_augmented_operator_mat);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetSizes(embedded_augmented_operator_mat, full_m_local, full_n_local, full_m, full_n);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetType(embedded_augmented_operator_mat, MATAIJ);
            IBTK_CHKERRQ(ierr);
            ierr = MatPreallocatorPreallocate(preallocator, PETSC_TRUE, embedded_augmented_operator_mat);
            IBTK_CHKERRQ(ierr);
            ierr = MatDestroy(&preallocator);
            IBTK_CHKERRQ(ierr);

            insert_velocity_block_rows(d_augmented_operator_mat, velocity_field_ao, embedded_augmented_operator_mat);
            ierr = MatAssemblyBegin(embedded_augmented_operator_mat, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);
            ierr = MatAssemblyEnd(embedded_augmented_operator_mat, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);

            ierr = MatAXPY(d_petsc_mat, 1.0, embedded_augmented_operator_mat, DIFFERENT_NONZERO_PATTERN);
            IBTK_CHKERRQ(ierr);
            ierr = MatDestroy(&embedded_augmented_operator_mat);
            IBTK_CHKERRQ(ierr);
            ierr = AODestroy(&velocity_field_ao);
            IBTK_CHKERRQ(ierr);
            ierr = ISDestroy(&velocity_field_is_local);
            IBTK_CHKERRQ(ierr);
        }
    }
    d_petsc_pc = d_petsc_mat;

    // Set pressure nullspace if the level covers the entire domain.
    if (d_has_pressure_nullspace)
    {
        const bool level_covers_entire_domain = level_covers_entire_physical_domain(d_level);

        if (level_covers_entire_domain)
        {
            // Allocate pressure nullspace data.
            if (!d_level->checkAllocated(d_u_nullspace_idx)) d_level->allocatePatchData(d_u_nullspace_idx);
            if (!d_level->checkAllocated(d_p_nullspace_idx)) d_level->allocatePatchData(d_p_nullspace_idx);

            Pointer<SAMRAIVectorReal<NDIM, double>> nullspace_vec = new SAMRAIVectorReal<NDIM, double>(
                d_object_name + "nullspace_vec", d_hierarchy, d_level_num, d_level_num);
            nullspace_vec->addComponent(d_u_nullspace_var, d_u_nullspace_idx);
            nullspace_vec->addComponent(d_p_nullspace_var, d_p_nullspace_idx);
            for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = d_level->getPatch(p());
                Pointer<SideData<NDIM, double>> u_patch_data = nullspace_vec->getComponentPatchData(0, *patch);
                u_patch_data->fill(0.0);
                Pointer<CellData<NDIM, double>> p_patch_data = nullspace_vec->getComponentPatchData(1, *patch);
                p_patch_data->fill(1.0);
            }

            LinearSolver::setNullSpace(
                /*const vec*/ false, std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>>(1, nullspace_vec));
        }
    }

    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    d_data_synch_sched = StaggeredStokesPETScVecUtilities::constructDataSynchSchedule(u_idx, p_idx, d_level);
    d_ghost_fill_sched = StaggeredStokesPETScVecUtilities::constructGhostFillSchedule(u_idx, p_idx, d_level);
    return;
} // initializeSolverStateSpecialized

void
StaggeredStokesPETScLevelSolver::deallocateSolverStateSpecialized()
{
    if (d_operator_mat && d_petsc_mat == d_operator_mat)
    {
        // The installed reference survives solver-state teardown. The base
        // class releases the KSP references, but must not release this alias.
        TBOX_ASSERT(d_petsc_pc == d_petsc_mat);
        d_petsc_mat = nullptr;
        d_petsc_pc = nullptr;
    }

    // Couplings and DOF numbering can change between solver lifetimes.
    if (d_asm_mode == ASMSubdomainConstructionMode::COUPLING_AWARE)
    {
        for (IS& is : d_overlap_is)
        {
            const int ierr = ISDestroy(&is);
            IBTK_CHKERRQ(ierr);
        }
        d_overlap_is.clear();
        for (IS& is : d_nonoverlap_is)
        {
            const int ierr = ISDestroy(&is);
            IBTK_CHKERRQ(ierr);
        }
        d_nonoverlap_is.clear();
    }
    d_ca_subdomains.reset();

    // Deallocate DOF index data.
    if (d_level->checkAllocated(d_u_dof_index_idx)) d_level->deallocatePatchData(d_u_dof_index_idx);
    if (d_level->checkAllocated(d_p_dof_index_idx)) d_level->deallocatePatchData(d_p_dof_index_idx);
    return;
} // deallocateSolverStateSpecialized

void
StaggeredStokesPETScLevelSolver::copyToPETScVec(Vec& petsc_x, SAMRAIVectorReal<NDIM, double>& x)
{
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        petsc_x, u_idx, d_u_dof_index_idx, p_idx, d_p_dof_index_idx, d_level);
    return;
} // copyToPETScVec

void
StaggeredStokesPETScLevelSolver::copyFromPETScVec(Vec& petsc_x, SAMRAIVectorReal<NDIM, double>& x)
{
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(
        petsc_x, u_idx, d_u_dof_index_idx, p_idx, d_p_dof_index_idx, d_level, d_data_synch_sched, d_ghost_fill_sched);
    return;
} // copyFromPETScVec

void
StaggeredStokesPETScLevelSolver::setupKSPVecs(Vec& petsc_x,
                                              Vec& petsc_b,
                                              SAMRAIVectorReal<NDIM, double>& x,
                                              SAMRAIVectorReal<NDIM, double>& b)
{
    if (d_initial_guess_nonzero) copyToPETScVec(petsc_x, x);
    const bool level_zero = (d_level_num == 0);
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    const int f_idx = b.getComponentDescriptorIndex(0);
    const int h_idx = b.getComponentDescriptorIndex(1);
    const auto f_adj_idx = d_cached_eulerian_data.getCachedPatchDataIndex(f_idx);
    const auto h_adj_idx = d_cached_eulerian_data.getCachedPatchDataIndex(h_idx);
    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = d_level->getPatch(p());
        Pointer<PatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        Pointer<SideData<NDIM, double>> f_data = patch->getPatchData(f_idx);
        Pointer<CellData<NDIM, double>> h_data = patch->getPatchData(h_idx);
        Pointer<SideData<NDIM, double>> f_adj_data = patch->getPatchData(f_adj_idx);
        Pointer<CellData<NDIM, double>> h_adj_data = patch->getPatchData(h_adj_idx);
        f_adj_data->copy(*f_data);
        h_adj_data->copy(*h_data);
        const bool at_physical_bdry = pgeom->intersectsPhysicalBoundary();
        // TODO: should we be using target data idx's here?
        StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
            d_U_bc_coefs, d_P_bc_coef, u_idx, p_idx, d_homogeneous_bc);
        if (at_physical_bdry)
        {
            PoissonUtilities::adjustRHSAtPhysicalBoundary(
                *f_adj_data, patch, d_U_problem_coefs, d_U_bc_coefs, d_solution_time, d_homogeneous_bc);
            d_bc_helper->enforceNormalVelocityBoundaryConditions(
                f_adj_idx, h_adj_idx, d_U_bc_coefs, d_solution_time, d_homogeneous_bc, d_level_num, d_level_num);
        }
        const Array<BoundaryBox<NDIM>>& type_1_cf_bdry = level_zero ?
                                                             Array<BoundaryBox<NDIM>>() :
                                                             d_cf_boundary->getBoundaries(patch->getPatchNumber(),
                                                                                          /* boundary type */ 1);
        const bool at_cf_bdry = type_1_cf_bdry.size() > 0;
        if (at_cf_bdry)
        {
            PoissonUtilities::adjustRHSAtCoarseFineBoundary(
                *f_adj_data, *u_data, patch, d_U_problem_coefs, type_1_cf_bdry);
        }
    }

    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        petsc_b, f_adj_idx, d_u_dof_index_idx, h_adj_idx, d_p_dof_index_idx, d_level);

    return;
} // setupKSPVecs

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
