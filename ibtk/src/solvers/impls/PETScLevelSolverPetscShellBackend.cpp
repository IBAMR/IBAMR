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

#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScLevelSolver.h>
#include <ibtk/private/PETScLevelSolverPetscShellBackend.h>

namespace IBTK
{
namespace
{
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
}

void
build_subdomain_index_sets(std::vector<IS>& subdomain_is,
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
}
} // namespace

PETScLevelSolverPetscShellBackend::PETScLevelSolverPetscShellBackend(PETScLevelSolver& solver) : d_context(solver)
{
}

const std::string&
PETScLevelSolverPetscShellBackend::getTypeKey() const
{
    return d_type_key;
}

void
PETScLevelSolverPetscShellBackend::configure(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> /*input_db*/)
{
}

const char*
PETScLevelSolverPetscShellBackend::getPCNameSuffixAdditive() const
{
    return "PC_AdditivePetsc";
}

const char*
PETScLevelSolverPetscShellBackend::getPCNameSuffixMultiplicative() const
{
    return "PC_MultiplicativePetsc";
}

void
PETScLevelSolverPetscShellBackend::initialize()
{
    d_data = std::make_unique<Data>();
    auto& petsc = *d_data;
    const int n_local_subdomains = static_cast<int>(d_context.getSubdomainDOFsForBackend().size());
    const bool use_restrict_partition =
        d_context.useRestrictPartitionForBackend();
    const bool use_multiplicative =
        d_context.isShellMultiplicativeForBackend();
    int ierr;
    if (use_multiplicative)
    {
        ierr = VecDuplicate(d_context.getPETScXForBackend(), &petsc.shell_r);
        IBTK_CHKERRQ(ierr);
    }
    petsc.prolongation_insert_mode = use_restrict_partition ? INSERT_VALUES : ADD_VALUES;
    build_subdomain_index_sets(petsc.global_overlap_is,
                               petsc.global_nonoverlap_is,
                               d_context.getSubdomainDOFsForBackend(),
                               d_context.getNonoverlapSubdomainDOFsForBackend());

    ierr = MatCreateSubMatrices(d_context.getPETScMatForBackend(),
                                n_local_subdomains,
                                petsc.global_overlap_is.data(),
                                petsc.global_overlap_is.data(),
                                MAT_INITIAL_MATRIX,
                                &petsc.sub_mat);
    IBTK_CHKERRQ(ierr);

    petsc.local_overlap_is.resize(n_local_subdomains);
    petsc.restriction.resize(n_local_subdomains);
    petsc.prolongation.resize(n_local_subdomains);
    petsc.sub_x.resize(n_local_subdomains);
    petsc.sub_y.resize(n_local_subdomains);
    if (use_restrict_partition) petsc.local_nonoverlap_is.resize(n_local_subdomains);
    if (use_multiplicative)
    {
        petsc.active_update_local_positions.resize(n_local_subdomains);
        petsc.active_residual_update_x.resize(n_local_subdomains);
        petsc.active_residual_update_y.resize(n_local_subdomains);
    }

#if !defined(NDEBUG)
    std::set<int> idxs;
#endif
    for (int subdomain_num = 0; subdomain_num < n_local_subdomains; ++subdomain_num)
    {
        const auto& overlap_dofs = d_context.getSubdomainDOFsForBackend()[static_cast<std::size_t>(subdomain_num)];
        const int overlap_is_size = static_cast<int>(overlap_dofs.size());
        PetscInt* overlap_indices = nullptr;
        ierr = PetscMalloc1(overlap_is_size, &overlap_indices);
        IBTK_CHKERRQ(ierr);
        for (int overlap_local_idx = 0; overlap_local_idx < overlap_is_size; ++overlap_local_idx)
        {
            overlap_indices[overlap_local_idx] = overlap_local_idx;
        }

        const auto& nonoverlap_dofs = d_context.getNonoverlapSubdomainDOFsForBackend()[static_cast<std::size_t>(subdomain_num)];
        const int nonoverlap_is_size = static_cast<int>(nonoverlap_dofs.size());
#if !defined(NDEBUG)
        for (const int dof : nonoverlap_dofs)
        {
            TBOX_ASSERT(idxs.find(dof) == idxs.end());
            idxs.insert(dof);
        }
#endif
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
            if (!use_restrict_partition)
            {
                active_update_local_positions.resize(static_cast<std::size_t>(overlap_is_size));
                for (int local_pos = 0; local_pos < overlap_is_size; ++local_pos)
                {
                    active_update_local_positions[static_cast<std::size_t>(local_pos)] = local_pos;
                }
            }
            else
            {
                active_update_local_positions.resize(static_cast<std::size_t>(nonoverlap_is_size));
                for (int local_pos = 0; local_pos < nonoverlap_is_size; ++local_pos)
                {
                    active_update_local_positions[static_cast<std::size_t>(local_pos)] =
                        static_cast<int>(nonoverlap_indices[local_pos]);
                }
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

        ierr = VecScatterCreate(d_context.getPETScXForBackend(),
                                petsc.global_overlap_is[subdomain_num],
                                petsc.sub_x[subdomain_num],
                                petsc.local_overlap_is[subdomain_num],
                                &petsc.restriction[subdomain_num]);
        IBTK_CHKERRQ(ierr);
        if (!use_restrict_partition)
        {
            ierr = VecScatterCreate(petsc.sub_y[subdomain_num],
                                    petsc.local_overlap_is[subdomain_num],
                                    d_context.getPETScBForBackend(),
                                    petsc.global_overlap_is[subdomain_num],
                                    &petsc.prolongation[subdomain_num]);
            IBTK_CHKERRQ(ierr);
        }
        else
        {
            ierr = VecScatterCreate(petsc.sub_y[subdomain_num],
                                    petsc.local_nonoverlap_is[subdomain_num],
                                    d_context.getPETScBForBackend(),
                                    petsc.global_nonoverlap_is[subdomain_num],
                                    &petsc.prolongation[subdomain_num]);
            IBTK_CHKERRQ(ierr);
        }
    }

#if !defined(NDEBUG)
    int n_local_dofs = 0;
    ierr = VecGetLocalSize(d_context.getPETScXForBackend(), &n_local_dofs);
    IBTK_CHKERRQ(ierr);
    TBOX_ASSERT(n_local_dofs == static_cast<int>(idxs.size()));
#endif

    if (use_multiplicative)
    {
        PetscInt row_begin = 0, row_end = 0;
        ierr = MatGetOwnershipRange(d_context.getPETScMatForBackend(), &row_begin, &row_end);
        IBTK_CHKERRQ(ierr);
        const PetscInt n_owned_rows = row_end - row_begin;
        std::vector<IS> active_update_global_is(static_cast<std::size_t>(n_local_subdomains));
        std::vector<IS> owned_residual_update_rows_is(static_cast<std::size_t>(n_local_subdomains), nullptr);
        ierr = ISCreateStride(PETSC_COMM_SELF, n_owned_rows, row_begin, 1, &petsc.owned_residual_update_rows_is);
        IBTK_CHKERRQ(ierr);
        for (int subdomain_num = 0; subdomain_num < n_local_subdomains; ++subdomain_num)
        {
            if (use_restrict_partition)
            {
                active_update_global_is[subdomain_num] = petsc.global_nonoverlap_is[subdomain_num];
            }
            else
            {
                active_update_global_is[subdomain_num] = petsc.global_overlap_is[subdomain_num];
            }
            owned_residual_update_rows_is[subdomain_num] = petsc.owned_residual_update_rows_is;
        }
        ierr = MatCreateSubMatrices(d_context.getPETScMatForBackend(),
                                    n_local_subdomains,
                                    owned_residual_update_rows_is.data(),
                                    active_update_global_is.data(),
                                    MAT_INITIAL_MATRIX,
                                    &petsc.active_residual_update_mat);
        IBTK_CHKERRQ(ierr);

        for (int subdomain_num = 0; subdomain_num < n_local_subdomains; ++subdomain_num)
        {
            ierr = MatCreateVecs(petsc.active_residual_update_mat[subdomain_num],
                                 &petsc.active_residual_update_x[subdomain_num],
                                 &petsc.active_residual_update_y[subdomain_num]);
            IBTK_CHKERRQ(ierr);
        }
    }

    petsc.sub_ksp.resize(n_local_subdomains);
    for (int subdomain_num = 0; subdomain_num < n_local_subdomains; ++subdomain_num)
    {
        KSP& sub_ksp = petsc.sub_ksp[subdomain_num];
        Mat& sub_mat = petsc.sub_mat[subdomain_num];
        ierr = KSPCreate(PETSC_COMM_SELF, &sub_ksp);
        IBTK_CHKERRQ(ierr);
        const std::string sub_prefix = d_context.getOptionsPrefixForBackend() + "_sub";
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
}

void
PETScLevelSolverPetscShellBackend::deallocate()
{
    if (!d_data) return;

    auto& petsc = *d_data;
    const int n_local_subdomains = static_cast<int>(petsc.sub_ksp.size());
    int ierr;
    for (int subdomain_num = 0; subdomain_num < n_local_subdomains; ++subdomain_num)
    {
        ierr = KSPDestroy(&petsc.sub_ksp[subdomain_num]);
        IBTK_CHKERRQ(ierr);
    }
    for (int subdomain_num = 0; subdomain_num < n_local_subdomains; ++subdomain_num)
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
    ierr = MatDestroyMatrices(n_local_subdomains, &petsc.sub_mat);
    IBTK_CHKERRQ(ierr);
    if (petsc.active_residual_update_mat)
    {
        ierr = MatDestroyMatrices(n_local_subdomains, &petsc.active_residual_update_mat);
        IBTK_CHKERRQ(ierr);
    }
    if (petsc.owned_residual_update_rows_is)
    {
        ierr = ISDestroy(&petsc.owned_residual_update_rows_is);
        IBTK_CHKERRQ(ierr);
    }
    destroy_index_sets(petsc.global_nonoverlap_is);
    destroy_index_sets(petsc.global_overlap_is);
    if (petsc.shell_r)
    {
        ierr = VecDestroy(&petsc.shell_r);
        IBTK_CHKERRQ(ierr);
    }
    d_data.reset();
}

void
PETScLevelSolverPetscShellBackend::beginAccumulateCorrection(const int subdomain_num, Vec sub_y, Vec y)
{
    auto& petsc = *d_data;
    const int ierr = VecScatterBegin(
        petsc.prolongation[subdomain_num], sub_y, y, petsc.prolongation_insert_mode, SCATTER_FORWARD_LOCAL);
    IBTK_CHKERRQ(ierr);
}

void
PETScLevelSolverPetscShellBackend::endAccumulateCorrection(const int subdomain_num, Vec sub_y, Vec y)
{
    auto& petsc = *d_data;
    const int ierr = VecScatterEnd(
        petsc.prolongation[subdomain_num], sub_y, y, petsc.prolongation_insert_mode, SCATTER_FORWARD_LOCAL);
    IBTK_CHKERRQ(ierr);
}

void
PETScLevelSolverPetscShellBackend::accumulateCorrection(const int subdomain_num, Vec sub_y, Vec y)
{
    beginAccumulateCorrection(subdomain_num, sub_y, y);
    endAccumulateCorrection(subdomain_num, sub_y, y);
}

void
PETScLevelSolverPetscShellBackend::applyAdditive(Vec x, Vec y)
{
    auto& petsc = *d_data;
    const int n_local_subdomains = static_cast<int>(petsc.sub_ksp.size());
    int ierr = VecZeroEntries(y);
    IBTK_CHKERRQ(ierr);
    for (int i = 0; i < n_local_subdomains; ++i)
    {
        ierr = VecScatterBegin(petsc.restriction[i], x, petsc.sub_x[i], INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
    }
    for (int i = 0; i < n_local_subdomains; ++i)
    {
        ierr = VecScatterEnd(petsc.restriction[i], x, petsc.sub_x[i], INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSolve(petsc.sub_ksp[i], petsc.sub_x[i], petsc.sub_y[i]);
        IBTK_CHKERRQ(ierr);
        accumulateCorrection(i, petsc.sub_y[i], y);
    }
    d_context.postprocessShellResultForBackend(y);
}

void
PETScLevelSolverPetscShellBackend::applyMultiplicative(Vec x, Vec y)
{
    auto& petsc = *d_data;
    const int n_local_subdomains = static_cast<int>(petsc.sub_ksp.size());
    int ierr = VecZeroEntries(y);
    IBTK_CHKERRQ(ierr);
    ierr = VecCopy(x, petsc.shell_r);
    IBTK_CHKERRQ(ierr);
    for (int i = 0; i < n_local_subdomains; ++i)
    {
        ierr = VecScatterBegin(petsc.restriction[i], petsc.shell_r, petsc.sub_x[i], INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
        ierr = VecScatterEnd(petsc.restriction[i], petsc.shell_r, petsc.sub_x[i], INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
        ierr = KSPSolve(petsc.sub_ksp[i], petsc.sub_x[i], petsc.sub_y[i]);
        IBTK_CHKERRQ(ierr);
        accumulateCorrection(i, petsc.sub_y[i], y);
        if (i + 1 < n_local_subdomains)
        {
            updateResidual(i, petsc.sub_y[i], petsc.shell_r);
        }
    }
    d_context.postprocessShellResultForBackend(y);
}

void
PETScLevelSolverPetscShellBackend::updateResidual(const int subdomain_num, Vec sub_y, Vec residual)
{
    auto& petsc = *d_data;
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
}
} // namespace IBTK
