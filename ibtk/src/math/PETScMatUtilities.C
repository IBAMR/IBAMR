// Filename: PETScMatUtilities.C
// Created on 24 Aug 2010 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#include "PETScMatUtilities.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PoissonUtilities.h>
#include <ibtk/namespaces.h>

// C++ STDLIB INCLUDES
#include <numeric>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
PETScMatUtilities::constructPatchLevelCCLaplaceOp(
    Mat& mat,
    const PoissonSpecifications& poisson_spec,
    RobinBcCoefStrategy<NDIM>* bc_coef,
    double data_time,
    const std::vector<int>& num_dofs_per_proc,
    const int dof_index_idx,
    Pointer<PatchLevel<NDIM> > patch_level)
{
    constructPatchLevelCCLaplaceOp(mat, poisson_spec, std::vector<RobinBcCoefStrategy<NDIM>*>(1,bc_coef), data_time, num_dofs_per_proc, dof_index_idx, patch_level);
    return;
}// constructPatchLevelCCLaplaceOp

void
PETScMatUtilities::constructPatchLevelCCLaplaceOp(
    Mat& mat,
    const PoissonSpecifications& poisson_spec,
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& bc_coefs,
    double data_time,
    const std::vector<int>& num_dofs_per_proc,
    const int dof_index_idx,
    Pointer<PatchLevel<NDIM> > patch_level)
{
    constructPatchLevelCCLaplaceOp(mat, poisson_spec, std::vector<RobinBcCoefStrategy<NDIM>*>(&bc_coefs[0],&bc_coefs[0]+NDIM), data_time, num_dofs_per_proc, dof_index_idx, patch_level);
    return;
}// constructPatchLevelCCLaplaceOp

void
PETScMatUtilities::constructPatchLevelCCLaplaceOp(
    Mat& mat,
    const PoissonSpecifications& poisson_spec,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
    double data_time,
    const std::vector<int>& num_dofs_per_proc,
    const int dof_index_idx,
    Pointer<PatchLevel<NDIM> > patch_level)
{
    int ierr;
    if (mat != PETSC_NULL)
    {
        ierr = MatDestroy(&mat); IBTK_CHKERRQ(ierr);
    }

    const int depth = bc_coefs.size();
    static const int stencil_sz = 2*NDIM+1;

    // Determine the index ranges.
    const int mpi_rank = SAMRAI_MPI::getRank();
    const int nlocal = num_dofs_per_proc[mpi_rank];
    const int ilower = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin()+mpi_rank, 0);
    const int iupper = ilower+nlocal;

    // Determine the non-zero structure of the matrix.
    std::vector<int> d_nnz(nlocal,0), o_nnz(nlocal,0);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM,int> > dof_index_data = patch->getPatchData(dof_index_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(depth == dof_index_data->getDepth());
#endif
        const Box<NDIM>& data_box = CellGeometry<NDIM>::toCellBox(patch_box);
        for (Box<NDIM>::Iterator b(data_box); b; b++)
        {
            const CellIndex<NDIM>& i = b();
            for (int d = 0; d < depth; ++d)
            {
                const int dof_index = (*dof_index_data)(i,d);
                if (UNLIKELY(ilower > dof_index || dof_index >= iupper)) continue;

                // Stencil for finite difference operator.
                const int local_idx = dof_index-ilower;
                d_nnz[local_idx] += 1;
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    for (int shift = -1; shift <= 1; shift += 2)
                    {
                        CellIndex<NDIM> i_shift(i);
                        i_shift(axis) += shift;
                        const int idx_shift = (*dof_index_data)(i_shift,d);
                        if (idx_shift >= ilower && idx_shift < iupper)
                        {
                            d_nnz[local_idx] += 1;
                        }
                        else
                        {
                            o_nnz[local_idx] += 1;
                        }
                    }
                }
            }
        }
    }

    // Create an empty matrix.
    ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,
                           nlocal, nlocal,
                           PETSC_DETERMINE, PETSC_DETERMINE,
                           PETSC_DEFAULT, &d_nnz[0],
                           PETSC_DEFAULT, &o_nnz[0],
                           &mat); IBTK_CHKERRQ(ierr);

    // Set some general matrix options.
    ierr = MatSetBlockSize(mat, depth); IBTK_CHKERRQ(ierr);
#ifdef DEBUG_CHECK_ASSERTIONS
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_LOCATION_ERR  , PETSC_TRUE); IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE); IBTK_CHKERRQ(ierr);
#endif

    // Set the matrix coefficients to correspond to the standard finite
    // difference Laplacian discretization.
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const IntVector<NDIM> no_ghosts(0);
        CellData<NDIM,double>     diagonal(patch_box, 1, no_ghosts);
        SideData<NDIM,double> off_diagonal(patch_box, 1, no_ghosts);
        PoissonUtilities::computeCCMatrixCoefficients(patch, diagonal, off_diagonal, poisson_spec, bc_coefs, data_time);
        Pointer<CellData<NDIM,int> > dof_index_data = patch->getPatchData(dof_index_idx);
        const Box<NDIM>& data_box = CellGeometry<NDIM>::toCellBox(patch_box);

        // Copy matrix entries to the PETSc matrix structure.
        std::vector<double> mat_vals(stencil_sz);
        std::vector<int>    mat_cols(stencil_sz);
        for (Box<NDIM>::Iterator b(data_box); b; b++)
        {
            const CellIndex<NDIM>& i = b();
            for (int d = 0; d < depth; ++d)
            {
                const int dof_index = (*dof_index_data)(i,d);
                if (UNLIKELY(ilower > dof_index || dof_index >= iupper)) continue;

                // Matrix coefficients for finite difference operator.
                mat_vals[0] = diagonal(i);
                mat_cols[0] = dof_index;
                for (unsigned int axis = 0, k = 1; axis < NDIM; ++axis)
                {
                    for (int side = 0; side <= 1; ++side, ++k)
                    {
                        CellIndex<NDIM> i_shift(i);
                        i_shift(axis) += (side == 0 ? -1 : +1);
                        mat_vals[k] = off_diagonal(SideIndex<NDIM>(i, axis, side));
                        mat_cols[k] = (*dof_index_data)(i_shift,d);
                    }
                }
                ierr = MatSetValues(mat, 1, &mat_cols[0], stencil_sz, &mat_cols[0], &mat_vals[0], INSERT_VALUES); IBTK_CHKERRQ(ierr);
            }
        }
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY); IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(  mat, MAT_FINAL_ASSEMBLY); IBTK_CHKERRQ(ierr);
    return;
}// constructPatchLevelLaplaceOp

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
