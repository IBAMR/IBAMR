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

#include <set>

namespace IBTK
{
namespace
{
std::unique_ptr<PETScLevelSolverShellBackend>
allocate_petsc_shell_backend(PETScLevelSolver& /*solver*/)
{
    return std::make_unique<PETScLevelSolverPetscShellBackend>();
}

class PETScLevelSolverPetscShellBackendRegistrar
{
public:
    PETScLevelSolverPetscShellBackendRegistrar()
    {
        PETScLevelSolverShellBackendManager::getManager()->registerShellBackendFactoryFunction(
            "petsc", allocate_petsc_shell_backend);
    }
};

static PETScLevelSolverPetscShellBackendRegistrar petsc_shell_backend_registrar;

} // namespace

const std::string PETScLevelSolverPetscShellBackend::s_backend_name = "Petsc";

const std::string&
PETScLevelSolverPetscShellBackend::getName() const
{
    return s_backend_name;
}

void
PETScLevelSolverPetscShellBackend::initializeSolverState(const PETScLevelSolverShellBackendState& solver_state)
{
    d_solver_state = solver_state;
    d_data = std::make_unique<Data>();
    auto& petsc = *d_data;
    const int n_local_subdomains = static_cast<int>(d_solver_state.subdomain_dofs->size());
    const bool use_restrict_partition = d_solver_state.use_restrict_partition;
    const bool use_multiplicative = d_solver_state.use_multiplicative;
    int ierr;
    if (use_multiplicative)
    {
        ierr = VecDuplicate(d_solver_state.petsc_x, &petsc.shell_r);
        IBTK_CHKERRQ(ierr);
    }
    petsc.prolongation_insert_mode = use_restrict_partition ? INSERT_VALUES : ADD_VALUES;
    build_petsc_subdomain_index_sets(petsc.global_overlap_is,
                                     petsc.global_nonoverlap_is,
                                     *d_solver_state.subdomain_dofs,
                                     *d_solver_state.nonoverlap_subdomain_dofs);

    ierr = MatCreateSubMatrices(d_solver_state.petsc_mat,
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

#if !defined(NDEBUG)
    std::set<int> idxs;
#endif
    for (int subdomain_num = 0; subdomain_num < n_local_subdomains; ++subdomain_num)
    {
        const auto& overlap_dofs = (*d_solver_state.subdomain_dofs)[static_cast<std::size_t>(subdomain_num)];
        const int overlap_is_size = static_cast<int>(overlap_dofs.size());
        PetscInt* overlap_indices = nullptr;
        ierr = PetscMalloc1(overlap_is_size, &overlap_indices);
        IBTK_CHKERRQ(ierr);
        for (int overlap_local_idx = 0; overlap_local_idx < overlap_is_size; ++overlap_local_idx)
        {
            overlap_indices[overlap_local_idx] = overlap_local_idx;
        }

        const auto& nonoverlap_dofs =
            (*d_solver_state.nonoverlap_subdomain_dofs)[static_cast<std::size_t>(subdomain_num)];
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

        ierr = VecScatterCreate(d_solver_state.petsc_x,
                                petsc.global_overlap_is[subdomain_num],
                                petsc.sub_x[subdomain_num],
                                petsc.local_overlap_is[subdomain_num],
                                &petsc.restriction[subdomain_num]);
        IBTK_CHKERRQ(ierr);
        if (!use_restrict_partition)
        {
            ierr = VecScatterCreate(petsc.sub_y[subdomain_num],
                                    petsc.local_overlap_is[subdomain_num],
                                    d_solver_state.petsc_b,
                                    petsc.global_overlap_is[subdomain_num],
                                    &petsc.prolongation[subdomain_num]);
            IBTK_CHKERRQ(ierr);
        }
        else
        {
            ierr = VecScatterCreate(petsc.sub_y[subdomain_num],
                                    petsc.local_nonoverlap_is[subdomain_num],
                                    d_solver_state.petsc_b,
                                    petsc.global_nonoverlap_is[subdomain_num],
                                    &petsc.prolongation[subdomain_num]);
            IBTK_CHKERRQ(ierr);
        }
    }

#if !defined(NDEBUG)
    int n_local_dofs = 0;
    ierr = VecGetLocalSize(d_solver_state.petsc_x, &n_local_dofs);
    IBTK_CHKERRQ(ierr);
    if (use_restrict_partition) TBOX_ASSERT(n_local_dofs == static_cast<int>(idxs.size()));
#endif

    petsc.sub_ksp.resize(n_local_subdomains);
    for (int subdomain_num = 0; subdomain_num < n_local_subdomains; ++subdomain_num)
    {
        KSP& sub_ksp = petsc.sub_ksp[subdomain_num];
        Mat& sub_mat = petsc.sub_mat[subdomain_num];
        ierr = KSPCreate(PETSC_COMM_SELF, &sub_ksp);
        IBTK_CHKERRQ(ierr);
        const std::string sub_prefix = d_solver_state.options_prefix + "_sub";
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
PETScLevelSolverPetscShellBackend::deallocateSolverState()
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
    }
    ierr = MatDestroyMatrices(n_local_subdomains, &petsc.sub_mat);
    IBTK_CHKERRQ(ierr);
    destroy_petsc_index_sets(petsc.global_nonoverlap_is);
    destroy_petsc_index_sets(petsc.global_overlap_is);
    if (petsc.shell_r)
    {
        ierr = VecDestroy(&petsc.shell_r);
        IBTK_CHKERRQ(ierr);
    }
    d_data.reset();
    d_solver_state = PETScLevelSolverShellBackendState();
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

std::size_t
PETScLevelSolverPetscShellBackend::getNumberOfSubdomains() const
{
    return d_data->sub_ksp.size();
}

void
PETScLevelSolverPetscShellBackend::initializeSubdomainSweep(Vec /*x*/, Vec /*y*/)
{
}

void
PETScLevelSolverPetscShellBackend::beginSubdomainRhs(const std::size_t subdomain_num, Vec x, Vec y)
{
    auto& petsc = *d_data;
    Vec source = x;
    int ierr;
    if (d_solver_state.use_multiplicative)
    {
        ierr = MatMult(d_solver_state.petsc_mat, y, petsc.shell_r);
        IBTK_CHKERRQ(ierr);
        ierr = VecAYPX(petsc.shell_r, -1.0, x);
        IBTK_CHKERRQ(ierr);
        source = petsc.shell_r;
    }
    ierr = VecScatterBegin(
        petsc.restriction[subdomain_num], source, petsc.sub_x[subdomain_num], INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);
}

void
PETScLevelSolverPetscShellBackend::endSubdomainRhs(const std::size_t subdomain_num, Vec x, Vec /*y*/)
{
    auto& petsc = *d_data;
    Vec source = d_solver_state.use_multiplicative ? petsc.shell_r : x;
    const int ierr = VecScatterEnd(
        petsc.restriction[subdomain_num], source, petsc.sub_x[subdomain_num], INSERT_VALUES, SCATTER_FORWARD);
    IBTK_CHKERRQ(ierr);
}

void
PETScLevelSolverPetscShellBackend::solveSubdomain(const std::size_t subdomain_num)
{
    auto& petsc = *d_data;
    const int ierr = KSPSolve(petsc.sub_ksp[subdomain_num], petsc.sub_x[subdomain_num], petsc.sub_y[subdomain_num]);
    IBTK_CHKERRQ(ierr);
}

void
PETScLevelSolverPetscShellBackend::accumulateSubdomainCorrection(const std::size_t subdomain_num, Vec y)
{
    auto& petsc = *d_data;
    if (!d_solver_state.use_multiplicative)
    {
        beginAccumulateCorrection(static_cast<int>(subdomain_num), petsc.sub_y[subdomain_num], y);
        endAccumulateCorrection(static_cast<int>(subdomain_num), petsc.sub_y[subdomain_num], y);
    }
    else
    {
        Vec y_sub = nullptr;
        int ierr = VecGetSubVector(y, petsc.global_overlap_is[subdomain_num], &y_sub);
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPY(y_sub, 1.0, petsc.sub_y[subdomain_num]);
        IBTK_CHKERRQ(ierr);
        ierr = VecRestoreSubVector(y, petsc.global_overlap_is[subdomain_num], &y_sub);
        IBTK_CHKERRQ(ierr);
    }
}

void
PETScLevelSolverPetscShellBackend::updateSubdomainResidual(std::size_t /*subdomain_num*/, Vec /*x*/, Vec /*y*/)
{
    // This reference path recomputes the original residual before each patch.
}
} // namespace IBTK
