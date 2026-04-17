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

#include <ibamr/StaggeredStokesPETScLevelSolver.h>
#include <ibamr/StaggeredStokesPETScMatUtilities.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/GeneralSolver.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LinearSolver.h>
#include <ibtk/PETScLevelSolver.h>
#include <ibtk/PoissonUtilities.h>

#include <tbox/Array.h>
#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <petsclog.h>
#include <petscvec.h>

#include <BoundaryBox.h>
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
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesPETScLevelSolver::StaggeredStokesPETScLevelSolver(const std::string& object_name,
                                                                 Pointer<Database> input_db,
                                                                 const std::string& default_options_prefix)
{
    GeneralSolver::init(object_name, /*homogeneous_bc*/ false);
    PETScLevelSolver::init(input_db, default_options_prefix);
    if (input_db && input_db->keyExists("subdomain_box_size"))
    {
        input_db->getIntegerArray("subdomain_box_size", d_box_size, NDIM);
    }
    if (input_db && input_db->keyExists("subdomain_overlap_size"))
    {
        input_db->getIntegerArray("subdomain_overlap_size", d_overlap_size, NDIM);
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
    return;
} // ~StaggeredStokesPETScLevelSolver

/////////////////////////////// PROTECTED ////////////////////////////////////

void
StaggeredStokesPETScLevelSolver::generateASMSubdomains(std::vector<std::set<int>>& overlap_is,
                                                       std::vector<std::set<int>>& nonoverlap_is)
{
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
        ierr = MatDuplicate(d_operator_mat, MAT_COPY_VALUES, &d_petsc_mat);
        IBTK_CHKERRQ(ierr);
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

    std::vector<std::set<int>> field_is;
    std::vector<std::string> field_names;
    StaggeredStokesPETScMatUtilities::constructPatchLevelFields(
        field_is, field_names, d_num_dofs_per_proc, d_u_dof_index_idx, d_p_dof_index_idx, d_level);
    const auto velocity_name_it = std::find(field_names.begin(), field_names.end(), "velocity");
    if (velocity_name_it == field_names.end())
    {
        TBOX_ERROR("StaggeredStokesPETScLevelSolver::initializeSolverStateSpecialized():\n"
                   << "  unable to locate velocity field DOFs.\n");
    }
    const std::size_t velocity_idx = static_cast<std::size_t>(std::distance(field_names.begin(), velocity_name_it));
    std::vector<PetscInt> velocity_dofs(field_is[velocity_idx].begin(), field_is[velocity_idx].end());

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
            IS velocity_is = nullptr, velocity_is_all = nullptr;
            ierr = ISCreateGeneral(PETSC_COMM_WORLD,
                                   static_cast<PetscInt>(velocity_dofs.size()),
                                   velocity_dofs.empty() ? nullptr : velocity_dofs.data(),
                                   PETSC_COPY_VALUES,
                                   &velocity_is);
            IBTK_CHKERRQ(ierr);
            ierr = ISAllGather(velocity_is, &velocity_is_all);
            IBTK_CHKERRQ(ierr);

            PetscInt n_velocity_global = 0;
            ierr = ISGetSize(velocity_is_all, &n_velocity_global);
            IBTK_CHKERRQ(ierr);
            if (aug_m != n_velocity_global || aug_n != n_velocity_global)
            {
                ierr = ISDestroy(&velocity_is_all);
                IBTK_CHKERRQ(ierr);
                ierr = ISDestroy(&velocity_is);
                IBTK_CHKERRQ(ierr);
                TBOX_ERROR("StaggeredStokesPETScLevelSolver::initializeSolverStateSpecialized():\n"
                           << "  augmented operator has incompatible size: (" << aug_m << " x " << aug_n << ").\n"
                           << "  expected either full operator size (" << full_m << " x " << full_n
                           << ") or velocity block size (" << n_velocity_global << " x " << n_velocity_global
                           << ").\n");
            }

            const PetscInt* velocity_global_ids = nullptr;
            ierr = ISGetIndices(velocity_is_all, &velocity_global_ids);
            IBTK_CHKERRQ(ierr);

            PetscInt row_start = 0, row_end = 0;
            ierr = MatGetOwnershipRange(d_augmented_operator_mat, &row_start, &row_end);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetOption(d_petsc_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
            IBTK_CHKERRQ(ierr);
            std::vector<PetscInt> mapped_cols;
            for (PetscInt row = row_start; row < row_end; ++row)
            {
                PetscInt ncols = 0;
                const PetscInt* cols = nullptr;
                const PetscScalar* vals = nullptr;
                ierr = MatGetRow(d_augmented_operator_mat, row, &ncols, &cols, &vals);
                IBTK_CHKERRQ(ierr);
                mapped_cols.resize(static_cast<std::size_t>(ncols));
                for (PetscInt k = 0; k < ncols; ++k)
                {
                    mapped_cols[static_cast<std::size_t>(k)] = velocity_global_ids[cols[k]];
                }
                const PetscInt full_row = velocity_global_ids[row];
                ierr = MatSetValues(d_petsc_mat, 1, &full_row, ncols, mapped_cols.data(), vals, ADD_VALUES);
                IBTK_CHKERRQ(ierr);
                ierr = MatRestoreRow(d_augmented_operator_mat, row, &ncols, &cols, &vals);
                IBTK_CHKERRQ(ierr);
            }
            ierr = MatAssemblyBegin(d_petsc_mat, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);
            ierr = MatAssemblyEnd(d_petsc_mat, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);

            ierr = ISRestoreIndices(velocity_is_all, &velocity_global_ids);
            IBTK_CHKERRQ(ierr);
            ierr = ISDestroy(&velocity_is_all);
            IBTK_CHKERRQ(ierr);
            ierr = ISDestroy(&velocity_is);
            IBTK_CHKERRQ(ierr);
        }
    }
    d_petsc_pc = d_petsc_mat;
    // Set pressure nullspace if the level covers the entire domain.
    if (d_has_pressure_nullspace)
    {
        bool level_covers_entire_domain = d_level_num == 0;
        if (d_level_num > 0)
        {
            int local_cf_bdry_box_size = 0;
            for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = d_level->getPatch(p());
                const Array<BoundaryBox<NDIM>>& type_1_cf_bdry = d_cf_boundary->getBoundaries(patch->getPatchNumber(),
                                                                                              /* boundary type */ 1);
                local_cf_bdry_box_size += type_1_cf_bdry.size();
            }
            level_covers_entire_domain = IBTK_MPI::sumReduction(local_cf_bdry_box_size) == 0;
        }

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

    copyToPETScVec(petsc_b, b);
    return;
} // setupKSPVecs

void
StaggeredStokesPETScLevelSolver::setOperatorMat(Mat operator_mat)
{
    d_operator_mat = operator_mat;
    return;
} // setOperatorMat

void
StaggeredStokesPETScLevelSolver::setAugmentedOperatorMat(Mat augmented_operator_mat)
{
    d_augmented_operator_mat = augmented_operator_mat;
    return;
} // setAugmentedOperatorMat

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
