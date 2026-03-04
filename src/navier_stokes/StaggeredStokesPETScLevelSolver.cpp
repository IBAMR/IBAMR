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
    if (input_db && input_db->keyExists("asm_subdomain_construction_mode"))
    {
        const std::string mode = input_db->getString("asm_subdomain_construction_mode");
        d_asm_subdomain_construction_mode = string_to_enum<ASMSubdomainConstructionMode>(mode);
        if (d_asm_subdomain_construction_mode == ASMSubdomainConstructionMode::UNKNOWN)
        {
            TBOX_ERROR("StaggeredStokesPETScLevelSolver::StaggeredStokesPETScLevelSolver():\n"
                       << "  invalid asm_subdomain_construction_mode = " << mode << "\n"
                       << "  expected values are " << enum_to_string(ASMSubdomainConstructionMode::GEOMETRICAL)
                       << " and " << enum_to_string(ASMSubdomainConstructionMode::COUPLING_AWARE) << ".\n");
        }
    }
    if (input_db && input_db->keyExists("coupling_aware_asm_seed_axis"))
    {
        d_coupling_aware_asm_seed_axis = input_db->getInteger("coupling_aware_asm_seed_axis");
        if (d_coupling_aware_asm_seed_axis < 0 || d_coupling_aware_asm_seed_axis >= NDIM)
        {
            TBOX_ERROR("StaggeredStokesPETScLevelSolver::StaggeredStokesPETScLevelSolver():\n"
                       << "  invalid coupling_aware_asm_seed_axis = " << d_coupling_aware_asm_seed_axis << "\n"
                       << "  expected value in [0, " << NDIM - 1 << "].\n");
        }
    }
    if (input_db && input_db->keyExists("coupling_aware_asm_seed_stride"))
    {
        d_coupling_aware_asm_seed_stride = input_db->getInteger("coupling_aware_asm_seed_stride");
        if (d_coupling_aware_asm_seed_stride < 1)
        {
            TBOX_ERROR("StaggeredStokesPETScLevelSolver::StaggeredStokesPETScLevelSolver():\n"
                       << "  invalid coupling_aware_asm_seed_stride = " << d_coupling_aware_asm_seed_stride << "\n"
                       << "  expected value >= 1.\n");
        }
    }
    if (input_db && input_db->keyExists("coupling_aware_asm_closure_policy"))
    {
        const std::string closure_policy = input_db->getString("coupling_aware_asm_closure_policy");
        d_coupling_aware_asm_closure_policy = string_to_enum<CouplingAwareASMClosurePolicy>(closure_policy);
        if (d_coupling_aware_asm_closure_policy == CouplingAwareASMClosurePolicy::UNKNOWN)
        {
            TBOX_ERROR("StaggeredStokesPETScLevelSolver::StaggeredStokesPETScLevelSolver():\n"
                       << "  invalid coupling_aware_asm_closure_policy = " << closure_policy << "\n"
                       << "  expected values are " << enum_to_string(CouplingAwareASMClosurePolicy::RELAXED) << " and "
                       << enum_to_string(CouplingAwareASMClosurePolicy::STRICT) << ".\n");
        }
    }
    if (input_db && input_db->keyExists("log_ASM_subdomains"))
    {
        d_log_asm_subdomains = input_db->getBool("log_ASM_subdomains");
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
    switch (d_asm_subdomain_construction_mode)
    {
    case ASMSubdomainConstructionMode::GEOMETRICAL:
        StaggeredStokesPETScMatUtilities::constructPatchLevelGeometricalASMSubdomains(overlap_is,
                                                                                      nonoverlap_is,
                                                                                      d_box_size,
                                                                                      d_overlap_size,
                                                                                      d_num_dofs_per_proc,
                                                                                      d_u_dof_index_idx,
                                                                                      d_p_dof_index_idx,
                                                                                      d_level,
                                                                                      d_cf_boundary);
        break;
    case ASMSubdomainConstructionMode::COUPLING_AWARE:
    {
        if (!d_petsc_mat)
        {
            TBOX_ERROR("StaggeredStokesPETScLevelSolver::generateASMSubdomains():\n"
                       << "  level matrix is not initialized for coupling-aware ASM subdomains.\n");
        }
        if (!d_coupling_aware_asm_map_data_is_initialized)
        {
            TBOX_ERROR("StaggeredStokesPETScLevelSolver::generateASMSubdomains():\n"
                       << "  coupling-aware ASM map data is not initialized.\n");
        }
        Mat A00_velocity_mat = nullptr;
        StaggeredStokesPETScMatUtilities::constructA00VelocitySubmatrix(
            A00_velocity_mat, d_petsc_mat, d_num_dofs_per_proc, d_u_dof_index_idx, d_p_dof_index_idx, d_level);
        StaggeredStokesPETScMatUtilities::constructPatchLevelCouplingAwareASMSubdomains(
            overlap_is,
            nonoverlap_is,
            d_box_size,
            d_overlap_size,
            d_num_dofs_per_proc,
            d_u_dof_index_idx,
            d_level,
            d_cf_boundary,
            A00_velocity_mat,
            d_coupling_aware_asm_map_data,
            d_coupling_aware_asm_seed_axis,
            d_coupling_aware_asm_seed_stride,
            d_coupling_aware_asm_closure_policy);
        int ierr = MatDestroy(&A00_velocity_mat);
        IBTK_CHKERRQ(ierr);
        break;
    }
    case ASMSubdomainConstructionMode::UNKNOWN:
    default:
        TBOX_ERROR("StaggeredStokesPETScLevelSolver::generateASMSubdomains():\n"
                   << "  unsupported asm_subdomain_construction_mode = "
                   << enum_to_string(d_asm_subdomain_construction_mode) << ".\n");
    }

    if (d_log_asm_subdomains)
    {
        plog << d_object_name
             << "::generateASMSubdomains(): mode = " << enum_to_string(d_asm_subdomain_construction_mode) << "\n";
        std::size_t nonoverlap_sum = 0;
        std::size_t overlap_sum = 0;
        for (std::size_t k = 0; k < overlap_is.size(); ++k)
        {
            const std::size_t n_nonoverlap = k < nonoverlap_is.size() ? nonoverlap_is[k].size() : 0;
            const std::size_t n_overlap = overlap_is[k].size();
            nonoverlap_sum += n_nonoverlap;
            overlap_sum += n_overlap;
            plog << "  subdomain " << k << ": nonoverlap = " << n_nonoverlap << ", overlap = " << n_overlap
                 << ", delta = " << static_cast<long long>(n_overlap) - static_cast<long long>(n_nonoverlap) << "\n";
        }
        plog << "  totals: nonoverlap = " << nonoverlap_sum << ", overlap = " << overlap_sum
             << ", delta = " << static_cast<long long>(overlap_sum) - static_cast<long long>(nonoverlap_sum) << "\n";
    }

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
    StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(d_petsc_mat,
                                                                     d_U_problem_coefs,
                                                                     d_U_bc_coefs,
                                                                     d_new_time,
                                                                     d_num_dofs_per_proc,
                                                                     d_u_dof_index_idx,
                                                                     d_p_dof_index_idx,
                                                                     d_level);
    d_petsc_pc = d_petsc_mat;
    d_coupling_aware_asm_map_data = StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData();
    d_coupling_aware_asm_map_data_is_initialized = false;
    if (d_asm_subdomain_construction_mode == ASMSubdomainConstructionMode::COUPLING_AWARE)
    {
        StaggeredStokesPETScMatUtilities::buildPatchLevelCellClosureMaps(
            d_coupling_aware_asm_map_data, d_u_dof_index_idx, d_p_dof_index_idx, d_level);
        d_coupling_aware_asm_map_data_is_initialized = true;
    }

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
    d_coupling_aware_asm_map_data = StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData();
    d_coupling_aware_asm_map_data_is_initialized = false;
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

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
