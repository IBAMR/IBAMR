// Filename: PETScMatUtilities.cpp
// Created on 24 Aug 2010 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <algorithm>
#include <numeric>
#include <ostream>
#include <vector>

#include "Box.h"
#include "BoxArray.h"
#include "BoxTree.h"
#include "CartesianGridGeometry.h"
#include "CellData.h"
#include "CellGeometry.h"
#include "CellIndex.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchLevel.h"
#include "PoissonSpecifications.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/PETScMatUtilities.h"
#include "ibtk/PoissonUtilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "petscmat.h"
#include "petscsys.h"
#include "petscvec.h"
#include "tbox/Array.h"
#include "tbox/Pointer.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI
// IWYU pragma: no_include "petsc-private/vecimpl.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

void PETScMatUtilities::constructPatchLevelCCLaplaceOp(Mat& mat,
                                                       const PoissonSpecifications& poisson_spec,
                                                       RobinBcCoefStrategy<NDIM>* bc_coef,
                                                       double data_time,
                                                       const std::vector<int>& num_dofs_per_proc,
                                                       const int dof_index_idx,
                                                       Pointer<PatchLevel<NDIM> > patch_level)
{
    constructPatchLevelCCLaplaceOp(mat,
                                   poisson_spec,
                                   std::vector<RobinBcCoefStrategy<NDIM>*>(1, bc_coef),
                                   data_time,
                                   num_dofs_per_proc,
                                   dof_index_idx,
                                   patch_level);
    return;
} // constructPatchLevelCCLaplaceOp

void PETScMatUtilities::constructPatchLevelCCComplexLaplaceOp(Mat& mat,
                                                              const PoissonSpecifications& poisson_spec_real,
                                                              const PoissonSpecifications& poisson_spec_imag,
                                                              RobinBcCoefStrategy<NDIM>* bc_coef,
                                                              double data_time,
                                                              const std::vector<int>& num_dofs_per_proc,
                                                              const int dof_index_idx,
                                                              Pointer<PatchLevel<NDIM> > patch_level)
{
    constructPatchLevelCCComplexLaplaceOp(mat,
                                          poisson_spec_real,
                                          poisson_spec_imag,
                                          std::vector<RobinBcCoefStrategy<NDIM>*>(2, bc_coef),
                                          data_time,
                                          num_dofs_per_proc,
                                          dof_index_idx,
                                          patch_level);
    return;
} // constructPatchLevelCCComplexLaplaceOp

void PETScMatUtilities::constructPatchLevelCCLaplaceOp(Mat& mat,
                                                       const PoissonSpecifications& poisson_spec,
                                                       const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                                       double data_time,
                                                       const std::vector<int>& num_dofs_per_proc,
                                                       const int dof_index_idx,
                                                       Pointer<PatchLevel<NDIM> > patch_level)
{
    int ierr;
    if (mat)
    {
        ierr = MatDestroy(&mat);
        IBTK_CHKERRQ(ierr);
    }

    const int depth = static_cast<int>(bc_coefs.size());

    // Setup the finite difference stencil.
    static const int stencil_sz = 2 * NDIM + 1;
    std::vector<Index<NDIM> > stencil(stencil_sz, Index<NDIM>(0));
    for (unsigned int axis = 0, stencil_index = 1; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side, ++stencil_index)
        {
            stencil[stencil_index](axis) = (side == 0 ? -1 : +1);
        }
    }

    // Determine the index ranges.
    const int mpi_rank = SAMRAI_MPI::getRank();
    const int n_local = num_dofs_per_proc[mpi_rank];
    const int i_lower = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int i_upper = i_lower + n_local;
    const int n_total = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.end(), 0);

    // Determine the non-zero structure of the matrix.
    std::vector<int> d_nnz(n_local, 0), o_nnz(n_local, 0);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, int> > dof_index_data = patch->getPatchData(dof_index_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(depth == dof_index_data->getDepth());
#endif
        for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
        {
            const CellIndex<NDIM>& i = b();
            for (int d = 0; d < depth; ++d)
            {
                const int dof_index = (*dof_index_data)(i, d);
                if (i_lower <= dof_index && dof_index < i_upper)
                {
                    // Stencil for finite difference operator.
                    const int local_idx = dof_index - i_lower;
                    d_nnz[local_idx] += 1;
                    for (unsigned int axis = 0, stencil_index = 1; axis < NDIM; ++axis)
                    {
                        for (int side = 0; side <= 1; ++side, ++stencil_index)
                        {
                            const int dof_index = (*dof_index_data)(i + stencil[stencil_index], d);
                            if (dof_index >= i_lower && dof_index < i_upper)
                            {
                                d_nnz[local_idx] += 1;
                            }
                            else
                            {
                                o_nnz[local_idx] += 1;
                            }
                        }
                    }
                    d_nnz[local_idx] = std::min(n_local, d_nnz[local_idx]);
                    o_nnz[local_idx] = std::min(n_total - n_local, o_nnz[local_idx]);
                }
            }
        }
    }

    // Create an empty matrix.
    ierr = MatCreateAIJ(PETSC_COMM_WORLD,
                        n_local,
                        n_local,
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        PETSC_DEFAULT,
                        &d_nnz[0],
                        PETSC_DEFAULT,
                        &o_nnz[0],
                        &mat);
    IBTK_CHKERRQ(ierr);

    // Set some general matrix options.
    ierr = MatSetBlockSize(mat, depth);
    IBTK_CHKERRQ(ierr);
#if !defined(NDEBUG)
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
#endif

    // Set the matrix coefficients to correspond to the standard finite
    // difference approximation to the Laplacian.
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();

        // Compute matrix coefficients.
        const IntVector<NDIM> no_ghosts(0);
        CellData<NDIM, double> matrix_coefs(patch_box, stencil_sz * depth, no_ghosts);
        PoissonUtilities::computeCCMatrixCoefficients(patch, matrix_coefs, stencil, poisson_spec, bc_coefs, data_time);

        // Copy matrix entries to the PETSc matrix structure.
        Pointer<CellData<NDIM, int> > dof_index_data = patch->getPatchData(dof_index_idx);
        std::vector<double> mat_vals(stencil_sz);
        std::vector<int> mat_cols(stencil_sz);
        for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
        {
            const CellIndex<NDIM>& i = b();
            for (int d = 0; d < depth; ++d)
            {
                const int dof_index = (*dof_index_data)(i, d);
                if (i_lower <= dof_index && dof_index < i_upper)
                {
                    // Notice that the order in which values are set corresponds
                    // to that of the stencil defined above.
                    const int offset = d * stencil_sz;
                    mat_vals[0] = matrix_coefs(i, offset);
                    mat_cols[0] = dof_index;
                    for (unsigned int axis = 0, stencil_index = 1; axis < NDIM; ++axis)
                    {
                        for (int side = 0; side <= 1; ++side, ++stencil_index)
                        {
                            mat_vals[stencil_index] = matrix_coefs(i, offset + stencil_index);
                            mat_cols[stencil_index] = (*dof_index_data)(i + stencil[stencil_index], d);
                        }
                    }
                    ierr = MatSetValues(mat, 1, &dof_index, stencil_sz, &mat_cols[0], &mat_vals[0], INSERT_VALUES);
                    IBTK_CHKERRQ(ierr);
                }
            }
        }
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    return;
} // constructPatchLevelCCLaplaceOp

void PETScMatUtilities::constructPatchLevelCCComplexLaplaceOp(Mat& mat,
                                                              const PoissonSpecifications& poisson_spec_real,
                                                              const PoissonSpecifications& poisson_spec_imag,
                                                              const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                                              double data_time,
                                                              const std::vector<int>& num_dofs_per_proc,
                                                              const int dof_index_idx,
                                                              Pointer<PatchLevel<NDIM> > patch_level)
{
    int ierr;
    if (mat)
    {
        ierr = MatDestroy(&mat);
        IBTK_CHKERRQ(ierr);
    }

    const int depth = static_cast<int>(bc_coefs.size());

    // Setup the finite difference stencil.
    static const int stencil_sz = 2 * NDIM + 1;
    std::vector<Index<NDIM> > stencil(stencil_sz, Index<NDIM>(0));
    for (unsigned int axis = 0, stencil_index = 1; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side, ++stencil_index)
        {
            stencil[stencil_index](axis) = (side == 0 ? -1 : +1);
        }
    }

    // Determine the index ranges.
    const int mpi_rank = SAMRAI_MPI::getRank();
    const int n_local = num_dofs_per_proc[mpi_rank];
    const int i_lower = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int i_upper = i_lower + n_local;
    const int n_total = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.end(), 0);

    // Determine the non-zero structure of the matrix.
    std::vector<int> d_nnz(n_local, 0), o_nnz(n_local, 0);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, int> > dof_index_data = patch->getPatchData(dof_index_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(depth == dof_index_data->getDepth());
#endif
        for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
        {
            const CellIndex<NDIM>& i = b();
            for (int d = 0; d < depth; ++d)
            {
                const int dof_index = (*dof_index_data)(i, d);
                if (i_lower <= dof_index && dof_index < i_upper)
                {
                    // Stencil for finite difference operator.
                    const int local_idx = dof_index - i_lower;
                    d_nnz[local_idx] += 2;
                    for (unsigned int axis = 0, stencil_index = 1; axis < NDIM; ++axis)
                    {
                        for (int side = 0; side <= 1; ++side, ++stencil_index)
                        {
                            const int dof_index = (*dof_index_data)(i + stencil[stencil_index], d);
                            if (dof_index >= i_lower && dof_index < i_upper)
                            {
                                d_nnz[local_idx] += 2;
                            }
                            else
                            {
                                o_nnz[local_idx] += 2;
                            }
                        }
                    }
                    d_nnz[local_idx] = std::min(n_local, d_nnz[local_idx]);
                    o_nnz[local_idx] = std::min(n_total - n_local, o_nnz[local_idx]);
                }
            }
        }
    }

    // Create an empty matrix.
    ierr = MatCreateAIJ(PETSC_COMM_WORLD,
                        n_local,
                        n_local,
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        PETSC_DEFAULT,
                        &d_nnz[0],
                        PETSC_DEFAULT,
                        &o_nnz[0],
                        &mat);
    IBTK_CHKERRQ(ierr);

    // Set some general matrix options.
    ierr = MatSetBlockSize(mat, depth);
    IBTK_CHKERRQ(ierr);
#if !defined(NDEBUG)
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
#endif

    // Set the matrix coefficients to correspond to the standard finite
    // difference approximation to the Laplacian.
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();

        // Compute matrix coefficients.
        const IntVector<NDIM> no_ghosts(0);
        CellData<NDIM, double> matrix_coefs(patch_box, 2 * stencil_sz * depth, no_ghosts);
        PoissonUtilities::computeCCComplexMatrixCoefficients(
            patch, matrix_coefs, stencil, poisson_spec_real, poisson_spec_imag, bc_coefs, data_time);

        // Copy matrix entries to the PETSc matrix structure.
        Pointer<CellData<NDIM, int> > dof_index_data = patch->getPatchData(dof_index_idx);
        std::vector<double> mat_vals_real(2 * stencil_sz), mat_vals_imag(2 * stencil_sz);
        std::vector<int> mat_cols_real(2 * stencil_sz), mat_cols_imag(2 * stencil_sz);
        for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
        {
            const CellIndex<NDIM>& i = b();
            for (int d = 0; d < depth; d = d + 2)
            {
                const int dof_index_real = (*dof_index_data)(i, d);
                const int dof_index_imag = (*dof_index_data)(i, d + 1);
                if ((i_lower <= dof_index_real && dof_index_real < i_upper) ||
                    (i_lower <= dof_index_imag && dof_index_imag < i_upper))
                {
#if !defined(NDEBUG)
                    TBOX_ASSERT(i_lower <= dof_index_real && dof_index_real < i_upper);
                    TBOX_ASSERT(i_lower <= dof_index_imag && dof_index_imag < i_upper);
#endif
                    // Notice that the order in which values are set corresponds
                    // to that of the stencil defined above.
                    const int offset = d * stencil_sz * 2;

                    mat_vals_real[0] = matrix_coefs(i, offset);
                    mat_vals_real[stencil_sz] = matrix_coefs(i, offset + stencil_sz);
                    mat_cols_real[0] = dof_index_real;
                    mat_cols_real[stencil_sz] = dof_index_imag;

                    mat_vals_imag[0] = matrix_coefs(i, offset + 2 * stencil_sz);
                    mat_vals_imag[stencil_sz] = matrix_coefs(i, offset + 3 * stencil_sz);
                    mat_cols_imag[0] = dof_index_real;
                    mat_cols_imag[stencil_sz] = dof_index_imag;

                    for (unsigned int axis = 0, stencil_index = 1; axis < NDIM; ++axis)
                    {
                        for (int side = 0; side <= 1; ++side, ++stencil_index)
                        {
                            mat_vals_real[stencil_index] = matrix_coefs(i, offset + stencil_index);
                            mat_vals_real[stencil_index + stencil_sz] =
                                matrix_coefs(i, offset + stencil_index + stencil_sz);
                            mat_cols_real[stencil_index] = (*dof_index_data)(i + stencil[stencil_index], d);
                            mat_cols_real[stencil_index + stencil_sz] =
                                (*dof_index_data)(i + stencil[stencil_index], d + 1);

                            mat_vals_imag[stencil_index] = matrix_coefs(i, offset + 2 * stencil_sz + stencil_index);
                            mat_vals_imag[stencil_index + stencil_sz] =
                                matrix_coefs(i, offset + stencil_index + 3 * stencil_sz);
                            mat_cols_imag[stencil_index] = (*dof_index_data)(i + stencil[stencil_index], d);
                            mat_cols_imag[stencil_index + stencil_sz] =
                                (*dof_index_data)(i + stencil[stencil_index], d + 1);
                        }
                    }
                    ierr = MatSetValues(
                        mat, 1, &dof_index_real, 2 * stencil_sz, &mat_cols_real[0], &mat_vals_real[0], INSERT_VALUES);
                    IBTK_CHKERRQ(ierr);
                    ierr = MatSetValues(
                        mat, 1, &dof_index_imag, 2 * stencil_sz, &mat_cols_imag[0], &mat_vals_imag[0], INSERT_VALUES);
                    IBTK_CHKERRQ(ierr);
                }
            }
        }
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    return;
} // constructPatchLevelCCComplexLaplaceOp

void PETScMatUtilities::constructPatchLevelSCLaplaceOp(Mat& mat,
                                                       const PoissonSpecifications& poisson_spec,
                                                       const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                                       double data_time,
                                                       const std::vector<int>& num_dofs_per_proc,
                                                       const int dof_index_idx,
                                                       Pointer<PatchLevel<NDIM> > patch_level)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_coefs.size() == NDIM);
#endif

    int ierr;
    if (mat)
    {
        ierr = MatDestroy(&mat);
        IBTK_CHKERRQ(ierr);
    }

    // Setup the finite difference stencil.
    static const int stencil_sz = 2 * NDIM + 1;
    std::vector<Index<NDIM> > stencil(stencil_sz, Index<NDIM>(0));
    for (unsigned int axis = 0, stencil_index = 1; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side, ++stencil_index)
        {
            stencil[stencil_index](axis) = (side == 0 ? -1 : +1);
        }
    }

    // Determine the index ranges.
    const int mpi_rank = SAMRAI_MPI::getRank();
    const int n_local = num_dofs_per_proc[mpi_rank];
    const int i_lower = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int i_upper = i_lower + n_local;
    const int n_total = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.end(), 0);

    // Determine the non-zero structure of the matrix.
    std::vector<int> d_nnz(n_local, 0), o_nnz(n_local, 0);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM, int> > dof_index_data = patch->getPatchData(dof_index_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(dof_index_data->getDepth() == 1);
#endif
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
            {
                const SideIndex<NDIM> i(b(), axis, SideIndex<NDIM>::Lower);
                const int dof_index = (*dof_index_data)(i);
                if (i_lower <= dof_index && dof_index < i_upper)
                {
                    // Stencil for finite difference operator.
                    const int local_idx = dof_index - i_lower;
                    d_nnz[local_idx] += 1;
                    for (unsigned int axis = 0, stencil_index = 1; axis < NDIM; ++axis)
                    {
                        for (int side = 0; side <= 1; ++side, ++stencil_index)
                        {
                            const int dof_index = (*dof_index_data)(i + stencil[stencil_index]);
                            if (dof_index >= i_lower && dof_index < i_upper)
                            {
                                d_nnz[local_idx] += 1;
                            }
                            else
                            {
                                o_nnz[local_idx] += 1;
                            }
                        }
                    }
                    d_nnz[local_idx] = std::min(n_local, d_nnz[local_idx]);
                    o_nnz[local_idx] = std::min(n_total - n_local, o_nnz[local_idx]);
                }
            }
        }
    }

    // Create an empty matrix.
    ierr = MatCreateAIJ(PETSC_COMM_WORLD,
                        n_local,
                        n_local,
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        PETSC_DEFAULT,
                        &d_nnz[0],
                        PETSC_DEFAULT,
                        &o_nnz[0],
                        &mat);
    IBTK_CHKERRQ(ierr);

// Set some general matrix options.
#if !defined(NDEBUG)
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
#endif

    // Set the matrix coefficients to correspond to the standard finite
    // difference approximation to the Laplacian.
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();

        // Compute matrix coefficients.
        const IntVector<NDIM> no_ghosts(0);
        SideData<NDIM, double> matrix_coefs(patch_box, stencil_sz, no_ghosts);
        PoissonUtilities::computeSCMatrixCoefficients(patch, matrix_coefs, stencil, poisson_spec, bc_coefs, data_time);

        // Copy matrix entries to the PETSc matrix structure.
        Pointer<SideData<NDIM, int> > dof_index_data = patch->getPatchData(dof_index_idx);
        std::vector<double> mat_vals(stencil_sz);
        std::vector<int> mat_cols(stencil_sz);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
            {
                const SideIndex<NDIM> i(b(), axis, SideIndex<NDIM>::Lower);
                const int dof_index = (*dof_index_data)(i);
                if (i_lower <= dof_index && dof_index < i_upper)
                {
                    // Notice that the order in which values are set corresponds
                    // to that of the stencil defined above.
                    mat_vals[0] = matrix_coefs(i, 0);
                    mat_cols[0] = dof_index;
                    for (unsigned int axis = 0, stencil_index = 1; axis < NDIM; ++axis)
                    {
                        for (int side = 0; side <= 1; ++side, ++stencil_index)
                        {
                            mat_vals[stencil_index] = matrix_coefs(i, stencil_index);
                            mat_cols[stencil_index] = (*dof_index_data)(i + stencil[stencil_index]);
                        }
                    }
                    ierr = MatSetValues(mat, 1, &dof_index, stencil_sz, &mat_cols[0], &mat_vals[0], INSERT_VALUES);
                    IBTK_CHKERRQ(ierr);
                }
            }
        }
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    return;
} // constructPatchLevelSCLaplaceOp

void PETScMatUtilities::constructPatchLevelSCInterpOp(Mat& mat,
                                                      void (*interp_fcn)(double r_lower, double* w),
                                                      int interp_stencil,
                                                      Vec& X_vec,
                                                      const std::vector<int>& num_dofs_per_proc,
                                                      const int dof_index_idx,
                                                      Pointer<PatchLevel<NDIM> > patch_level)
{
    // \todo Properly support odd stencil sizes.
    if (interp_stencil % 2 != 0) interp_stencil += 1;

    int ierr;
    if (mat)
    {
        ierr = MatDestroy(&mat);
        IBTK_CHKERRQ(ierr);
    }

    // Determine the grid extents.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = patch_level->getGridGeometry();
    const double* const x_lower = grid_geom->getXLower();
    const double* const x_upper = grid_geom->getXUpper();
    const double* const dx0 = grid_geom->getDx();
    const IntVector<NDIM>& ratio = patch_level->getRatio();
    double dx[NDIM];
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        dx[d] = dx0[d] / static_cast<double>(ratio(d));
    }
    const BoxArray<NDIM>& domain_boxes = patch_level->getPhysicalDomain();
#if !defined(NDEBUG)
    TBOX_ASSERT(domain_boxes.size() == 1);
#endif
    const Index<NDIM>& domain_lower = domain_boxes[0].lower();
    const Index<NDIM>& domain_upper = domain_boxes[0].upper();

    // Determine the matrix dimensions and index ranges.
    int m_local;
    ierr = VecGetLocalSize(X_vec, &m_local);
    IBTK_CHKERRQ(ierr);
    int i_lower, i_upper;
    ierr = VecGetOwnershipRange(X_vec, &i_lower, &i_upper);
    IBTK_CHKERRQ(ierr);

    const int mpi_rank = SAMRAI_MPI::getRank();
    const int n_local = num_dofs_per_proc[mpi_rank];
    const int j_lower = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int j_upper = j_lower + n_local;
    const int n_total = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.end(), 0);

    // Determine the index of the Cartesian grid cell containing each local IB
    // point; find that index in a local patch or in the ghost cell region of a
    // local patch; compute the stencil boxes for each local IB point; and
    // compute the nonzero structure of the matrix.
    const int n_local_points = m_local / NDIM;
    double* X_arr;
    ierr = VecGetArray(X_vec, &X_arr);
    IBTK_CHKERRQ(ierr);
    std::vector<int> patch_num(n_local_points);
    std::vector<std::vector<Box<NDIM> > > stencil_box(n_local_points, std::vector<Box<NDIM> >(NDIM));
    std::vector<int> d_nnz(m_local, 0), o_nnz(m_local, 0);
    for (int k = 0; k < n_local_points; ++k)
    {
        const double* const X = &X_arr[NDIM * k];
        const Index<NDIM>& X_idx = IndexUtilities::getCellIndex(X, x_lower, x_upper, dx, domain_lower, domain_upper);

// Determine the position of the center of the Cartesian grid cell
// containing the IB point.
#if (NDIM == 2)
        const double X_cell[NDIM] = { (static_cast<double>(X_idx(0) - domain_lower(0)) + 0.5) * dx[0] + x_lower[0],
                                      (static_cast<double>(X_idx(1) - domain_lower(1)) + 0.5) * dx[1] + x_lower[1] };
#endif
#if (NDIM == 3)
        const double X_cell[NDIM] = { (static_cast<double>(X_idx(0) - domain_lower(0)) + 0.5) * dx[0] + x_lower[0],
                                      (static_cast<double>(X_idx(1) - domain_lower(1)) + 0.5) * dx[1] + x_lower[1],
                                      (static_cast<double>(X_idx(2) - domain_lower(2)) + 0.5) * dx[2] + x_lower[2] };
#endif
        // Find a local patch that contains the IB point in either its patch
        // interior or ghost cell region.
        Box<NDIM> box(X_idx, X_idx);
        Array<int> patch_num_arr;
        patch_level->getBoxTree()->findOverlapIndices(patch_num_arr, box);
        if (patch_num_arr.size() == 0)
        {
            box.grow(IntVector<NDIM>(1));
            patch_level->getBoxTree()->findOverlapIndices(patch_num_arr, box);
#if !defined(NDEBUG)
            TBOX_ASSERT(patch_num_arr.size() != 0);
#endif
        }
        patch_num[k] = patch_num_arr[0];
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(patch_num[k]);
        Pointer<SideData<NDIM, int> > dof_index_data = patch->getPatchData(dof_index_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(dof_index_data->getDepth() == 1);
#endif

        // Compute the stencil box and setup the nonzero structure.
        for (int axis = 0; axis < NDIM; ++axis)
        {
            // Determine the stencil box.
            if (interp_stencil % 2 != 0)
            {
                TBOX_ERROR(
                    "PETScMatUtilities::constructPatchLevelSCInterpOp(): support for odd "
                    "stencil "
                    "sizes not currently implemented\n");
            }
            Box<NDIM>& stencil_box_axis = stencil_box[k][axis];
            Index<NDIM>& stencil_box_lower = stencil_box_axis.lower();
            Index<NDIM>& stencil_box_upper = stencil_box_axis.upper();
            for (int d = 0; d < NDIM; ++d)
            {
                if (d == axis)
                {
                    stencil_box_lower(d) = X_idx(d) - interp_stencil / 2 + 1;
                    stencil_box_upper(d) = X_idx(d) + interp_stencil / 2;
                }
                else if (X[d] <= X_cell[d])
                {
                    stencil_box_lower(d) = X_idx(d) - interp_stencil / 2;
                    stencil_box_upper(d) = X_idx(d) + interp_stencil / 2 - 1;
                }
                else
                {
                    stencil_box_lower(d) = X_idx(d) - interp_stencil / 2 + 1;
                    stencil_box_upper(d) = X_idx(d) + interp_stencil / 2;
                }
            }
            const int local_idx = NDIM * k + axis;
#if !defined(NDEBUG)
            TBOX_ASSERT(SideGeometry<NDIM>::toSideBox(dof_index_data->getGhostBox(), axis).contains(stencil_box_axis));
#endif
            for (Box<NDIM>::Iterator b(stencil_box_axis); b; b++)
            {
                const int dof_index = (*dof_index_data)(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Lower));
                if (dof_index >= j_lower && dof_index < j_upper)
                {
                    d_nnz[local_idx] += 1;
                }
                else
                {
                    o_nnz[local_idx] += 1;
                }
            }
            d_nnz[local_idx] = std::min(n_local, d_nnz[local_idx]);
            o_nnz[local_idx] = std::min(n_total - n_local, o_nnz[local_idx]);
        }
    }

    // Create an empty matrix.
    ierr = MatCreateAIJ(PETSC_COMM_WORLD,
                        m_local,
                        n_local,
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        PETSC_DEFAULT,
                        &d_nnz[0],
                        PETSC_DEFAULT,
                        &o_nnz[0],
                        &mat);
    IBTK_CHKERRQ(ierr);

// Set some general matrix options.
#if !defined(NDEBUG)
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
#endif

    // Set the matrix coefficients.
    for (int k = 0; k < m_local / NDIM; ++k)
    {
        const double* const X = &X_arr[NDIM * k];

        // Look-up the local patch that we have associated with this IB point.
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(patch_num[k]);
        Pointer<SideData<NDIM, int> > dof_index_data = patch->getPatchData(dof_index_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(dof_index_data->getDepth() == 1);
#endif

        // Construct the interpolation weights for this IB point.
        std::vector<double> w[NDIM];
        for (int d = 0; d < NDIM; ++d) w[d].resize(interp_stencil);
        int stencil_box_nvals = 1;
        for (unsigned int d = 0; d < NDIM; ++d) stencil_box_nvals *= interp_stencil;
        std::vector<double> stencil_box_vals(stencil_box_nvals);
        std::vector<int> stencil_box_cols(stencil_box_nvals);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            // Look-up the stencil box.
            const Box<NDIM>& stencil_box_axis = stencil_box[k][axis];
            const Index<NDIM>& stencil_box_lower = stencil_box_axis.lower();

            // Compute the weights of the 1-dimensional delta functions.
            for (int d = 0; d < NDIM; ++d)
            {
                const int i = stencil_box_lower(d);
                const double X_stencil_lower =
                    (static_cast<double>(i + 1 - domain_lower(d)) + (d == axis ? 0.0 : 0.5)) * dx[d] + x_lower[d];
                interp_fcn((X[d] - X_stencil_lower) / dx[d], &w[d][0]);
            }

            // Compute the weights of the d-dimensional delta function as the
            // tensor products of the 1-dimensional delta functions.
            int stencil_box_row = i_lower + NDIM * k + axis;
            int stencil_idx = 0;
            std::fill(stencil_box_vals.begin(), stencil_box_vals.end(), 1.0);
            for (Box<NDIM>::Iterator b(stencil_box_axis); b; b++, ++stencil_idx)
            {
                const SideIndex<NDIM> i(b(), axis, SideIndex<NDIM>::Lower);
                for (int d = 0; d < NDIM; ++d)
                {
                    stencil_box_vals[stencil_idx] *= w[d][i(d) - stencil_box_lower(d)];
                }
                stencil_box_cols[stencil_idx] = (*dof_index_data)(i);
            }

            // Set the values for this IB point.
            ierr = MatSetValues(
                mat, 1, &stencil_box_row, stencil_box_nvals, &stencil_box_cols[0], &stencil_box_vals[0], INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
    }
    ierr = VecRestoreArray(X_vec, &X_arr);
    IBTK_CHKERRQ(ierr);

    // Assemble the matrix.
    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    return;
} // constructPatchLevelSCInterpOp

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
