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

#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/private/PETScLevelSolverBlasLapackShellBackend.h>
#include <ibtk/private/PETScLevelSolverPetscShellBackend.h>

#include <tbox/Utilities.h>

#include <algorithm>
#include <limits>

namespace IBTK
{
namespace
{
std::unique_ptr<PETScLevelSolverShellBackend>
allocate_petsc_backend(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> /*input_db*/)
{
    return std::make_unique<PETScLevelSolverPetscShellBackend>();
}

std::unique_ptr<PETScLevelSolverShellBackend>
allocate_blas_lapack_backend(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
{
    return std::make_unique<PETScLevelSolverBlasLapackShellBackend>(input_db);
}
} // namespace

PETScLevelSolverShellBackendManager&
PETScLevelSolverShellBackendManager::get_manager()
{
    static PETScLevelSolverShellBackendManager manager;
    return manager;
}

void
PETScLevelSolverShellBackendManager::registerFactory(const std::string& key, Factory factory)
{
    if (key.empty() || !factory)
    {
        TBOX_ERROR("PETScLevelSolverShellBackendManager::registerFactory():\n"
                   << "  require a nonempty key and a nonnull factory.\n");
    }
    d_factories[key] = factory;
}

std::unique_ptr<PETScLevelSolverShellBackend>
PETScLevelSolverShellBackendManager::allocateBackend(const std::string& key,
                                                     SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db) const
{
    const std::map<std::string, Factory>::const_iterator factory = d_factories.find(key);
    if (factory == d_factories.end())
    {
        TBOX_ERROR("PETScLevelSolverShellBackendManager::allocateBackend():\n"
                   << "  unknown additive shell backend: " << key << "\n");
    }
    std::unique_ptr<PETScLevelSolverShellBackend> backend = factory->second(input_db);
    if (!backend)
    {
        TBOX_ERROR("PETScLevelSolverShellBackendManager::allocateBackend():\n"
                   << "  factory returned a null backend: " << key << "\n");
    }
    return backend;
}

PETScLevelSolverShellBackendManager::PETScLevelSolverShellBackendManager()
    : d_factories{ { "petsc", allocate_petsc_backend }, { "blas-lapack", allocate_blas_lapack_backend } }
{
}

PETScLevelSolverPetscShellBackend::~PETScLevelSolverPetscShellBackend()
{
    deallocateSolverState();
}

void
PETScLevelSolverPetscShellBackend::initializeSolverState(Mat mat,
                                                         Vec x,
                                                         Vec b,
                                                         const std::vector<IS>& overlap,
                                                         const std::vector<IS>& nonoverlap,
                                                         const std::string& options_prefix)
{
    deallocateSolverState();
    TBOX_ASSERT(overlap.size() == nonoverlap.size());
    const int n_local = static_cast<int>(overlap.size());
    const int n_max = IBTK_MPI::maxReduction(n_local);
    Mat* sub_mat = nullptr;
    int ierr = MatCreateSubMatrices(mat, n_local, overlap.data(), overlap.data(), MAT_INITIAL_MATRIX, &sub_mat);
    IBTK_CHKERRQ(ierr);
    d_sub_ksp.resize(n_local);
    d_sub_x.resize(n_max);
    d_sub_y.resize(n_max);
    d_restriction.resize(n_max);
    d_prolongation.resize(n_max);
    for (int i = 0; i < n_max; ++i)
    {
        PetscInt n_overlap = 0, n_nonoverlap = 0;
        std::vector<PetscInt> local_nonoverlap;
        if (i < n_local)
        {
            ierr = ISGetLocalSize(overlap[i], &n_overlap);
            IBTK_CHKERRQ(ierr);
            ierr = ISGetLocalSize(nonoverlap[i], &n_nonoverlap);
            IBTK_CHKERRQ(ierr);
            const PetscInt* overlap_indices = nullptr;
            const PetscInt* nonoverlap_indices = nullptr;
            ierr = ISGetIndices(overlap[i], &overlap_indices);
            IBTK_CHKERRQ(ierr);
            ierr = ISGetIndices(nonoverlap[i], &nonoverlap_indices);
            IBTK_CHKERRQ(ierr);
            for (PetscInt j = 0; j < n_nonoverlap; ++j)
            {
                const PetscInt* position =
                    std::lower_bound(overlap_indices, overlap_indices + n_overlap, nonoverlap_indices[j]);
                TBOX_ASSERT(position != overlap_indices + n_overlap && *position == nonoverlap_indices[j]);
                local_nonoverlap.push_back(static_cast<PetscInt>(position - overlap_indices));
            }
            ierr = ISRestoreIndices(nonoverlap[i], &nonoverlap_indices);
            IBTK_CHKERRQ(ierr);
            ierr = ISRestoreIndices(overlap[i], &overlap_indices);
            IBTK_CHKERRQ(ierr);
            ierr = MatCreateVecs(sub_mat[i], &d_sub_x[i], &d_sub_y[i]);
            IBTK_CHKERRQ(ierr);
            ierr = KSPCreate(PETSC_COMM_SELF, &d_sub_ksp[i]);
            IBTK_CHKERRQ(ierr);
            const std::string prefix = options_prefix + "_sub";
            ierr = KSPSetOptionsPrefix(d_sub_ksp[i], prefix.c_str());
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetOperators(d_sub_ksp[i], sub_mat[i], sub_mat[i]);
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetReusePreconditioner(d_sub_ksp[i], PETSC_TRUE);
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetType(d_sub_ksp[i], KSPPREONLY);
            IBTK_CHKERRQ(ierr);
            PC pc = nullptr;
            ierr = KSPGetPC(d_sub_ksp[i], &pc);
            IBTK_CHKERRQ(ierr);
            ierr = PCSetType(pc, PCLU);
            IBTK_CHKERRQ(ierr);
            ierr = PCFactorReorderForNonzeroDiagonal(pc, std::numeric_limits<double>::epsilon());
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetFromOptions(d_sub_ksp[i]);
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetInitialGuessNonzero(d_sub_ksp[i], PETSC_FALSE);
            IBTK_CHKERRQ(ierr);
        }
        else
        {
            ierr = VecCreateSeq(PETSC_COMM_SELF, 0, &d_sub_x[i]);
            IBTK_CHKERRQ(ierr);
            ierr = VecDuplicate(d_sub_x[i], &d_sub_y[i]);
            IBTK_CHKERRQ(ierr);
        }
        // All ranks create the same number of scatters, including empty subdomains.
        IS local_overlap = nullptr, local_partition = nullptr;
        ierr = ISCreateStride(PETSC_COMM_WORLD, n_overlap, 0, 1, &local_overlap);
        IBTK_CHKERRQ(ierr);
        ierr = ISCreateGeneral(
            PETSC_COMM_WORLD, n_nonoverlap, local_nonoverlap.data(), PETSC_COPY_VALUES, &local_partition);
        IBTK_CHKERRQ(ierr);
        IS global_overlap = i < n_local ? overlap[i] : local_overlap;
        IS global_partition = i < n_local ? nonoverlap[i] : local_partition;
        ierr = VecScatterCreate(x, global_overlap, d_sub_x[i], local_overlap, &d_restriction[i]);
        IBTK_CHKERRQ(ierr);
        ierr = VecScatterCreate(d_sub_y[i], local_partition, b, global_partition, &d_prolongation[i]);
        IBTK_CHKERRQ(ierr);
        ierr = ISDestroy(&local_overlap);
        IBTK_CHKERRQ(ierr);
        ierr = ISDestroy(&local_partition);
        IBTK_CHKERRQ(ierr);
    }
    // Each sub-KSP retains its submatrix.
    ierr = MatDestroyMatrices(n_local, &sub_mat);
    IBTK_CHKERRQ(ierr);
}

void
PETScLevelSolverPetscShellBackend::deallocateSolverState()
{
    for (KSP& ksp : d_sub_ksp)
    {
        const int ierr = KSPDestroy(&ksp);
        IBTK_CHKERRQ(ierr);
    }
    d_sub_ksp.clear();
    for (std::size_t i = 0; i < d_sub_x.size(); ++i)
    {
        int ierr = VecScatterDestroy(&d_restriction[i]);
        IBTK_CHKERRQ(ierr);
        ierr = VecScatterDestroy(&d_prolongation[i]);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&d_sub_x[i]);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&d_sub_y[i]);
        IBTK_CHKERRQ(ierr);
    }
    d_restriction.clear();
    d_prolongation.clear();
    d_sub_x.clear();
    d_sub_y.clear();
}

void
PETScLevelSolverPetscShellBackend::apply(Vec x, Vec y)
{
    int ierr = VecZeroEntries(y);
    IBTK_CHKERRQ(ierr);
    for (std::size_t i = 0; i < d_sub_x.size(); ++i)
    {
        ierr = VecScatterBegin(d_restriction[i], x, d_sub_x[i], INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
    }
    for (std::size_t i = 0; i < d_sub_x.size(); ++i)
    {
        ierr = VecScatterEnd(d_restriction[i], x, d_sub_x[i], INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
        if (i < d_sub_ksp.size())
        {
            ierr = KSPSolve(d_sub_ksp[i], d_sub_x[i], d_sub_y[i]);
            IBTK_CHKERRQ(ierr);
        }
        ierr = VecScatterBegin(d_prolongation[i], d_sub_y[i], y, INSERT_VALUES, SCATTER_FORWARD_LOCAL);
        IBTK_CHKERRQ(ierr);
        ierr = VecScatterEnd(d_prolongation[i], d_sub_y[i], y, INSERT_VALUES, SCATTER_FORWARD_LOCAL);
        IBTK_CHKERRQ(ierr);
    }
}
} // namespace IBTK
