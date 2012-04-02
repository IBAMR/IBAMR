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

    // Setup the finite difference stencil.
    static const int stencil_sz = 2*NDIM+1;
    std::vector<Index<NDIM> > stencil(stencil_sz,Index<NDIM>(0));
    for (unsigned int axis = 0, stencil_index = 1; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side, ++stencil_index)
        {
            stencil[stencil_index](axis) = (side == 0 ? -1 : +1);
        }
    }

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
                for (unsigned int axis = 0, stencil_index = 1; axis < NDIM; ++axis)
                {
                    for (int side = 0; side <= 1; ++side, ++stencil_index)
                    {
                        const int dof_index = (*dof_index_data)(i+stencil[stencil_index],d);
                        if (dof_index >= ilower && dof_index < iupper)
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
        CellData<NDIM,double> matrix_coefs(patch_box, stencil_sz*depth, no_ghosts);
        PoissonUtilities::computeCCMatrixCoefficients(patch, matrix_coefs, stencil, poisson_spec, bc_coefs, data_time);
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

                // Matrix coefficients for finite difference operator.  Notice
                // that the order in which values are set corresponds to that of
                // the stencil defined above.
                const int offset = d*stencil_sz;
                mat_vals[0] = matrix_coefs(i,offset);
                mat_cols[0] = dof_index;
                for (unsigned int axis = 0, stencil_index = 1; axis < NDIM; ++axis)
                {
                    for (int side = 0; side <= 1; ++side, ++stencil_index)
                    {
                        mat_vals[stencil_index] = matrix_coefs(i,offset+stencil_index);
                        mat_cols[stencil_index] = (*dof_index_data)(i+stencil[stencil_index],d);
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

void
PETScMatUtilities::constructPatchLevelSCLaplaceOp(
    Mat& mat,
    const PoissonSpecifications& poisson_spec,
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& bc_coefs,
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

    // Setup the finite difference stencil.
    static const int stencil_sz = 2*NDIM+1;
    std::vector<Index<NDIM> > stencil(stencil_sz,Index<NDIM>(0));
    for (unsigned int axis = 0, stencil_index = 1; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side, ++stencil_index)
        {
            stencil[stencil_index](axis) = (side == 0 ? -1 : +1);
        }
    }

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
        Pointer<SideData<NDIM,int> > dof_index_data = patch->getPatchData(dof_index_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(dof_index_data->getDepth() == 1);
#endif
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            const Box<NDIM>& data_box = SideGeometry<NDIM>::toSideBox(patch_box,axis);
            for (Box<NDIM>::Iterator b(data_box); b; b++)
            {
                const SideIndex<NDIM> i(b(),axis,SideIndex<NDIM>::Lower);
                const int dof_index = (*dof_index_data)(i);
                if (UNLIKELY(ilower > dof_index || dof_index >= iupper)) continue;

                // Stencil for finite difference operator.
                const int local_idx = dof_index-ilower;
                d_nnz[local_idx] += 1;
                for (unsigned int axis = 0, stencil_index = 1; axis < NDIM; ++axis)
                {
                    for (int side = 0; side <= 1; ++side, ++stencil_index)
                    {
                        const int dof_index = (*dof_index_data)(i+stencil[stencil_index]);
                        if (dof_index >= ilower && dof_index < iupper)
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
        SideData<NDIM,double> matrix_coefs(patch_box, stencil_sz, no_ghosts);
        PoissonUtilities::computeSCMatrixCoefficients(patch, matrix_coefs, stencil, poisson_spec, bc_coefs, data_time);
        Pointer<SideData<NDIM,int> > dof_index_data = patch->getPatchData(dof_index_idx);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            const Box<NDIM>& data_box = SideGeometry<NDIM>::toSideBox(patch_box,axis);

            // Copy matrix entries to the PETSc matrix structure.
            std::vector<double> mat_vals(stencil_sz);
            std::vector<int>    mat_cols(stencil_sz);
            for (Box<NDIM>::Iterator b(data_box); b; b++)
            {
                const SideIndex<NDIM> i(b(),axis,SideIndex<NDIM>::Lower);
                const int dof_index = (*dof_index_data)(i);
                if (UNLIKELY(ilower > dof_index || dof_index >= iupper)) continue;

                // Matrix coefficients for finite difference operator.  Notice
                // that the order in which values are set corresponds to that of
                // the stencil defined above.
                mat_vals[0] = matrix_coefs(i,0);
                mat_cols[0] = dof_index;
                for (unsigned int axis = 0, stencil_index = 1; axis < NDIM; ++axis)
                {
                    for (int side = 0; side <= 1; ++side, ++stencil_index)
                    {
                        mat_vals[stencil_index] = matrix_coefs(i,stencil_index);
                        mat_cols[stencil_index] = (*dof_index_data)(i+stencil[stencil_index]);
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
