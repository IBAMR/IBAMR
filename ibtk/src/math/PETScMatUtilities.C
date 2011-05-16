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
#include <ibtk/IndexUtilities.h>
#include <ibtk/SideSynchCopyFillPattern.h>
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <RefineAlgorithm.h>

// C++ STDLIB INCLUDES
#include <algorithm>
#include <map>
#include <numeric>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline void
getIndexRange(
    int& n_local,
    int& i_lower,
    int& i_upper,
    const int data_idx,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > data_var,
    Pointer<PatchLevel<NDIM> > patch_level)
{
    n_local = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<CellData<NDIM,double> > data = patch->getPatchData(data_idx);
        n_local += data->getDepth()*CellGeometry<NDIM>::toCellBox(data->getGhostBox()).size();
    }
    const int mpi_size = SAMRAI_MPI::getNodes();
    const int mpi_rank = SAMRAI_MPI::getRank();
    std::vector<int> n_proc(mpi_size,0);
    SAMRAI_MPI::allGather(n_local, &n_proc[0]);
    i_lower = std::accumulate(n_proc.begin(), n_proc.begin()+mpi_rank, 0);
    i_upper = i_lower+n_local;
    return;
}// getIndexRange

inline void
getIndexRange(
    int& n_local,
    int& i_lower,
    int& i_upper,
    const int data_idx,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > data_var,
    Pointer<PatchLevel<NDIM> > patch_level)
{
    n_local = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<SideData<NDIM,double> > data = patch->getPatchData(data_idx);
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            n_local += data->getDepth()*SideGeometry<NDIM>::toSideBox(data->getGhostBox(), component_axis).size();
        }
    }
    const int mpi_size = SAMRAI_MPI::getNodes();
    const int mpi_rank = SAMRAI_MPI::getRank();
    std::vector<int> n_proc(mpi_size,0);
    SAMRAI_MPI::allGather(n_local, &n_proc[0]);
    i_lower = std::accumulate(n_proc.begin(), n_proc.begin()+mpi_rank, 0);
    i_upper = i_lower+n_local;
    return;
}// getIndexRange
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
PETScMatUtilities::constructPatchLaplaceOp(
    Mat& mat,
    const double C,
    const double D,
    CellData<NDIM,double>& src_data,
    CellData<NDIM,double>& dst_data,
    Patch<NDIM>& patch)
{
    const Box<NDIM>& patch_box = patch.getBox();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_box == src_data.getBox());
    TBOX_ASSERT(patch_box == dst_data.getBox());
#endif
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    if (mat != static_cast<Mat>(NULL))
    {
        int ierr = MatDestroy(mat); IBTK_CHKERRQ(ierr);
    }
    const Box<NDIM>& dst_ghost_box = dst_data.getGhostBox();
    const Box<NDIM>& src_ghost_box = src_data.getGhostBox();
    const int data_depth = dst_data.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_depth == src_data.getDepth());
#endif
    constructBoxLaplaceOp(mat, C, D, src_ghost_box, dst_ghost_box, patch_box, data_depth, dx);
    return;
}// constructPatchLaplaceOp

void
PETScMatUtilities::constructPatchLaplaceOps(
    blitz::TinyVector<Mat,NDIM>& mats,
    const double C,
    const double D,
    SideData<NDIM,double>& src_data,
    SideData<NDIM,double>& dst_data,
    Patch<NDIM>& patch)
{
    const Box<NDIM>& patch_box = patch.getBox();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_box == src_data.getBox());
    TBOX_ASSERT(patch_box == dst_data.getBox());
#endif
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
    {
        if (mats[component_axis] != static_cast<Mat>(NULL))
        {
            int ierr = MatDestroy(mats[component_axis]); IBTK_CHKERRQ(ierr);
        }
        const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, component_axis);
        const Box<NDIM> dst_ghost_box = SideGeometry<NDIM>::toSideBox(dst_data.getGhostBox(), component_axis);
        const Box<NDIM> src_ghost_box = SideGeometry<NDIM>::toSideBox(src_data.getGhostBox(), component_axis);
        const int data_depth = dst_data.getDepth();
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_depth == src_data.getDepth());
#endif
        constructBoxLaplaceOp(mats[component_axis], C, D, src_ghost_box, dst_ghost_box, side_box, data_depth, dx);
    }
    return;
}// constructPatchLaplaceOps

void
PETScMatUtilities::constructPatchLevelLaplaceOp(
    Mat& mat,
    const double C,
    const double D,
    const int data_idx,
    Pointer<CellVariable<NDIM,double> > data_var,
    const int dof_index_idx,
    Pointer<CellVariable<NDIM,int> > dof_index_var,
    Pointer<PatchLevel<NDIM> > patch_level,
    Pointer<RefineSchedule<NDIM> > dof_index_fill)
{
    int ierr;
    if (mat != static_cast<Mat>(NULL))
    {
        ierr = MatDestroy(mat); IBTK_CHKERRQ(ierr);
    }

    static const int stencil_sz = 2*NDIM+1;

    // Create a clone of the DOF index data and fill ghost cell values using the
    // original DOF index data.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int dof_master_index_idx = var_db->registerClonedPatchDataIndex(dof_index_var, dof_index_idx);
    patch_level->allocatePatchData(dof_master_index_idx);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<CellData<NDIM,int> > dof_index = patch->getPatchData(dof_index_idx);
        Pointer<CellData<NDIM,int> > dof_master_index = patch->getPatchData(dof_master_index_idx);
        dof_master_index->copy(*dof_index);
    }
    RefineAlgorithm<NDIM> ref_algorithm;
    ref_algorithm.registerRefine(dof_master_index_idx, dof_index_idx, dof_master_index_idx, NULL);
    if (dof_index_fill.isNull())
    {
        dof_index_fill = ref_algorithm.createSchedule(patch_level);
    }
    else
    {
        ref_algorithm.resetSchedule(dof_index_fill);
    }
    dof_index_fill->fillData(0.0);

    // Determine the index ranges.
    int n_local, i_lower, i_upper;
    getIndexRange(n_local, i_lower, i_upper, data_idx, data_var, patch_level);

    // Determine the non-zero structure of the matrix.
    std::vector<int> d_nnz(n_local,0), o_nnz(n_local,0);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<CellData<NDIM,double> > data = patch->getPatchData(data_idx);
        Pointer<CellData<NDIM,int> > dof_index = patch->getPatchData(dof_index_idx);
        Pointer<CellData<NDIM,int> > dof_master_index = patch->getPatchData(dof_master_index_idx);
        const Box<NDIM>& patch_box = CellGeometry<NDIM>::toCellBox(patch->getBox());
        const Box<NDIM>& ghost_box = CellGeometry<NDIM>::toCellBox(data->getGhostBox());
        const int data_depth = data->getDepth();
        for (int d = 0; d < data_depth; ++d)
        {
            for (Box<NDIM>::Iterator b(ghost_box); b; b++)
            {
                const Index<NDIM>& i = b();
                const int idx = (*dof_index)(i,d);
                const int mastr_idx = (*dof_master_index)(i,d);
                const int local_idx = idx-i_lower;
                if (patch_box.contains(i) && idx == mastr_idx)
                {
                    // Stencil for finite difference operator.
                    d_nnz[local_idx] += 1;
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        for (int shift = -1; shift <= 1; shift += 2)
                        {
                            CellIndex<NDIM> i_shift(i);
                            i_shift(axis) += shift;
                            const int idx_shift = (*dof_master_index)(i_shift,d);
                            if (idx_shift >= i_lower && idx_shift < i_upper)
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
                else if (idx == mastr_idx)
                {
                    // Local boundary value.
                    d_nnz[local_idx] += 1;
                }
                else if (mastr_idx >= i_lower && mastr_idx < i_upper)
                {
                    // Local DOF constraint.
                    d_nnz[local_idx] += 2;
                }
                else
                {
                    // Non-local DOF constraint.
                    d_nnz[local_idx] += 1;
                    o_nnz[local_idx] += 1;
                }
            }
        }
    }

    // Create an empty matrix.
    ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,
                           n_local, n_local,
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
    //
    // Note that boundary conditions both at physical boundaries and at
    // coarse-fine interfaces are assumed to be treated implicitly by setting
    // ghost cell values.  This allows the matrix coefficients to be independent
    // of the boundary conditions.
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<CellData<NDIM,double> > data = patch->getPatchData(data_idx);
        Pointer<CellData<NDIM,int> > dof_index = patch->getPatchData(dof_index_idx);
        Pointer<CellData<NDIM,int> > dof_master_index = patch->getPatchData(dof_master_index_idx);
        const Box<NDIM>& patch_box = CellGeometry<NDIM>::toCellBox(patch->getBox());
        const Box<NDIM>& ghost_box = CellGeometry<NDIM>::toCellBox(data->getGhostBox());
        const int data_depth = data->getDepth();

        const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();
        std::vector<double> mat_vals(stencil_sz,0.0);
        mat_vals[NDIM] = C;
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            const double& h = dx[axis];
            mat_vals[NDIM-axis-1] +=     D/(h*h);
            mat_vals[NDIM       ] -= 2.0*D/(h*h);
            mat_vals[NDIM+axis+1] +=     D/(h*h);
        }

        // Set matrix coefficients.
        for (int d = 0; d < data_depth; ++d)
        {
            for (Box<NDIM>::Iterator b(ghost_box); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                static const int m = 1;
                const int idxm = (*dof_index)(i,d);
                const int mastr_idxm = (*dof_master_index)(i,d);
                if (patch_box.contains(i) && idxm == mastr_idxm)
                {
                    // Finite difference stencil.
                    const int n = stencil_sz;
                    std::vector<int> idxn(stencil_sz);
                    idxn[NDIM] = (*dof_master_index)(i,d);
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        for (int shift = -1; shift <= 1; shift += 2)
                        {
                            CellIndex<NDIM> i_shift(i);
                            i_shift(axis) += shift;
                            idxn[NDIM+shift*(axis+1)] = (*dof_master_index)(i_shift,d);
                        }
                    }
                    ierr = MatSetValues(mat, m, &idxm, n, &idxn[0], &mat_vals[0], INSERT_VALUES); IBTK_CHKERRQ(ierr);
                }
                else if (idxm == mastr_idxm)
                {
                    // Boundary value.
                    ierr = MatSetValue(mat, idxm, idxm, +1.0, INSERT_VALUES); IBTK_CHKERRQ(ierr);
                }
                else
                {
                    // DOF constraint.
                    ierr = MatSetValue(mat, idxm, idxm, +1.0, INSERT_VALUES); IBTK_CHKERRQ(ierr);
                    const int idxn = (*dof_master_index)(i,d);
                    ierr = MatSetValue(mat, idxm, idxn, -1.0, INSERT_VALUES); IBTK_CHKERRQ(ierr);
                }
            }
        }
    }

    // Free the temporary data and assemble the matrix.
    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY); IBTK_CHKERRQ(ierr);
    patch_level->deallocatePatchData(dof_master_index_idx);
    var_db->removePatchDataIndex(dof_master_index_idx);
    ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY); IBTK_CHKERRQ(ierr);
    return;
}// constructPatchLevelLaplaceOp

void
PETScMatUtilities::constructPatchLevelLaplaceOp(
    Mat& mat,
    const double C,
    const double D,
    const int data_idx,
    Pointer<SideVariable<NDIM,double> > data_var,
    const int dof_index_idx,
    Pointer<SideVariable<NDIM,int> > dof_index_var,
    Pointer<PatchLevel<NDIM> > patch_level,
    Pointer<RefineSchedule<NDIM> > dof_index_fill)
{
    int ierr;
    if (mat != static_cast<Mat>(NULL))
    {
        ierr = MatDestroy(mat); IBTK_CHKERRQ(ierr);
    }

    static const int stencil_sz = 2*NDIM+1;

    // Create a clone of the DOF index data and fill ghost cell values using the
    // original DOF index data.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int dof_master_index_idx = var_db->registerClonedPatchDataIndex(dof_index_var, dof_index_idx);
    patch_level->allocatePatchData(dof_master_index_idx);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<SideData<NDIM,int> > dof_index = patch->getPatchData(dof_index_idx);
        Pointer<SideData<NDIM,int> > dof_master_index = patch->getPatchData(dof_master_index_idx);
        dof_master_index->copy(*dof_index);
    }
    RefineAlgorithm<NDIM> ref_algorithm;
    ref_algorithm.registerRefine(dof_master_index_idx, dof_index_idx, dof_master_index_idx, NULL);
    ref_algorithm.registerRefine(dof_master_index_idx, dof_master_index_idx, dof_master_index_idx, NULL, new SideSynchCopyFillPattern());
    if (dof_index_fill.isNull())
    {
        dof_index_fill = ref_algorithm.createSchedule(patch_level);
    }
    else
    {
        ref_algorithm.resetSchedule(dof_index_fill);
    }
    dof_index_fill->fillData(0.0);

    // Determine the index ranges.
    int n_local, i_lower, i_upper;
    getIndexRange(n_local, i_lower, i_upper, data_idx, data_var, patch_level);

    // Determine the non-zero structure of the matrix.
    std::vector<int> d_nnz(n_local,0), o_nnz(n_local,0);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<SideData<NDIM,double> > data = patch->getPatchData(data_idx);
        Pointer<SideData<NDIM,int> > dof_index = patch->getPatchData(dof_index_idx);
        Pointer<SideData<NDIM,int> > dof_master_index = patch->getPatchData(dof_master_index_idx);
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            const Box<NDIM>& patch_box = SideGeometry<NDIM>::toSideBox(patch->getBox(), component_axis);
            const Box<NDIM>& ghost_box = SideGeometry<NDIM>::toSideBox(data->getGhostBox(), component_axis);
            const int data_depth = data->getDepth();
            for (int d = 0; d < data_depth; ++d)
            {
                for (Box<NDIM>::Iterator b(ghost_box); b; b++)
                {
                    const SideIndex<NDIM> i(b(), component_axis, SideIndex<NDIM>::Lower);
                    const int idx = (*dof_index)(i,d);
                    const int mastr_idx = (*dof_master_index)(i,d);
                    const int local_idx = idx-i_lower;
                    if (patch_box.contains(i) && idx == mastr_idx)
                    {
                        // Stencil for finite difference operator.
                        d_nnz[local_idx] += 1;
                        for (unsigned int axis = 0; axis < NDIM; ++axis)
                        {
                            for (int shift = -1; shift <= 1; shift += 2)
                            {
                                SideIndex<NDIM> i_shift(i);
                                i_shift(axis) += shift;
                                const int idx_shift = (*dof_master_index)(i_shift,d);
                                if (idx_shift >= i_lower && idx_shift < i_upper)
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
                    else if (idx == mastr_idx)
                    {
                        // Local boundary value.
                        d_nnz[local_idx] += 1;
                    }
                    else if (mastr_idx >= i_lower && mastr_idx < i_upper)
                    {
                        // Local DOF constraint.
                        d_nnz[local_idx] += 2;
                    }
                    else
                    {
                        // Non-local DOF constraint.
                        d_nnz[local_idx] += 1;
                        o_nnz[local_idx] += 1;
                    }
                }
            }
        }
    }

    // Create an empty matrix.
    ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,
                           n_local, n_local,
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
    //
    // Note that boundary conditions both at physical boundaries and at
    // coarse-fine interfaces are assumed to be treated implicitly by setting
    // ghost cell values.  This allows the matrix coefficients to be independent
    // of the boundary conditions.
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<SideData<NDIM,double> > data = patch->getPatchData(data_idx);
        Pointer<SideData<NDIM,int> > dof_index = patch->getPatchData(dof_index_idx);
        Pointer<SideData<NDIM,int> > dof_master_index = patch->getPatchData(dof_master_index_idx);
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            const Box<NDIM>& patch_box = SideGeometry<NDIM>::toSideBox(patch->getBox(), component_axis);
            const Box<NDIM>& ghost_box = SideGeometry<NDIM>::toSideBox(data->getGhostBox(), component_axis);
            const int data_depth = data->getDepth();

            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            std::vector<double> mat_vals(stencil_sz,0.0);
            mat_vals[NDIM] = C;
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const double& h = dx[axis];
                mat_vals[NDIM-axis-1] +=     D/(h*h);
                mat_vals[NDIM       ] -= 2.0*D/(h*h);
                mat_vals[NDIM+axis+1] +=     D/(h*h);
            }

            // Set matrix coefficients.
            for (int d = 0; d < data_depth; ++d)
            {
                for (Box<NDIM>::Iterator b(ghost_box); b; b++)
                {
                    const SideIndex<NDIM> i(b(), component_axis, SideIndex<NDIM>::Lower);
                    static const int m = 1;
                    const int idxm = (*dof_index)(i,d);
                    const int mastr_idxm = (*dof_master_index)(i,d);
                    if (patch_box.contains(i) && idxm == mastr_idxm)
                    {
                        // Finite difference stencil.
                        const int n = stencil_sz;
                        std::vector<int> idxn(stencil_sz);
                        idxn[NDIM] = (*dof_master_index)(i,d);
                        for (unsigned int axis = 0; axis < NDIM; ++axis)
                        {
                            for (int shift = -1; shift <= 1; shift += 2)
                            {
                                SideIndex<NDIM> i_shift(i);
                                i_shift(axis) += shift;
                                idxn[NDIM+shift*(axis+1)] = (*dof_master_index)(i_shift,d);
                            }
                        }
                        ierr = MatSetValues(mat, m, &idxm, n, &idxn[0], &mat_vals[0], INSERT_VALUES); IBTK_CHKERRQ(ierr);
                    }
                    else if (idxm == mastr_idxm)
                    {
                        // Boundary value.
                        ierr = MatSetValue(mat, idxm, idxm, +1.0, INSERT_VALUES); IBTK_CHKERRQ(ierr);
                    }
                    else
                    {
                        // DOF constraint.
                        ierr = MatSetValue(mat, idxm, idxm, +1.0, INSERT_VALUES); IBTK_CHKERRQ(ierr);
                        const int idxn = (*dof_master_index)(i,d);
                        ierr = MatSetValue(mat, idxm, idxn, -1.0, INSERT_VALUES); IBTK_CHKERRQ(ierr);
                    }
                }
            }
        }
    }

    // Free the temporary data and assemble the matrix.
    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY); IBTK_CHKERRQ(ierr);
    patch_level->deallocatePatchData(dof_master_index_idx);
    var_db->removePatchDataIndex(dof_master_index_idx);
    ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY); IBTK_CHKERRQ(ierr);
    return;
}// constructPatchLevelLaplaceOp

namespace
{
// Basic 4-point delta function kernel.
static const int stencil_sz = 4;
inline double
ib4_delta_fcn(
    double r)
{
    r = std::abs(r);
    if (r < 1.0)
    {
        const double t2 = r * r;
        const double t6 = sqrt(-0.4e1 * t2 + 0.4e1 * r + 0.1e1);
        return -r / 0.4e1 + 0.3e1 / 0.8e1 + t6 / 0.8e1;
    }
    else if (r < 2.0)
    {
        const double t2 = r * r;
        const double t6 = sqrt(0.12e2 * r - 0.7e1 - 0.4e1 * t2);
        return -r / 0.4e1 + 0.5e1 / 0.8e1 - t6 / 0.8e1;
    }
    else
    {
        return 0.0;
    }
}// ib4_delta_fcn
}

void
PETScMatUtilities::constructPatchLevelInterpOp(
    Mat& mat,
    Vec& X_vec,
    const int data_idx,
    Pointer<SideVariable<NDIM,double> > data_var,
    const int dof_index_idx,
    Pointer<SideVariable<NDIM,int> > dof_index_var,
    Pointer<PatchLevel<NDIM> > patch_level,
    Pointer<RefineSchedule<NDIM> > dof_index_fill)
{
    int ierr;
    if (mat != static_cast<Mat>(NULL))
    {
        ierr = MatDestroy(mat); IBTK_CHKERRQ(ierr);
    }

    // Create a clone of the DOF index data and fill ghost cell values using the
    // original DOF index data.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int dof_master_index_idx = var_db->registerClonedPatchDataIndex(dof_index_var, dof_index_idx);
    patch_level->allocatePatchData(dof_master_index_idx);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<SideData<NDIM,int> > dof_index = patch->getPatchData(dof_index_idx);
        Pointer<SideData<NDIM,int> > dof_master_index = patch->getPatchData(dof_master_index_idx);
        dof_master_index->copy(*dof_index);
    }
    RefineAlgorithm<NDIM> ref_algorithm;
    ref_algorithm.registerRefine(dof_master_index_idx, dof_index_idx, dof_master_index_idx, NULL);
    ref_algorithm.registerRefine(dof_master_index_idx, dof_master_index_idx, dof_master_index_idx, NULL, new SideSynchCopyFillPattern());
    if (dof_index_fill.isNull())
    {
        dof_index_fill = ref_algorithm.createSchedule(patch_level);
    }
    else
    {
        ref_algorithm.resetSchedule(dof_index_fill);
    }
    dof_index_fill->fillData(0.0);

    // Determine the matrix dimensions and index ranges.
    int m_local;
    ierr = VecGetLocalSize(X_vec, &m_local); IBTK_CHKERRQ(ierr);
    int i_lower, i_upper;
    ierr = VecGetOwnershipRange(X_vec, &i_lower, &i_upper); IBTK_CHKERRQ(ierr);

    int n_local, j_lower, j_upper;
    getIndexRange(n_local, j_lower, j_upper, data_idx, data_var, patch_level);

    // Crudely approximate the non-zero structure of the matrix.
    std::vector<int> d_nnz(m_local,static_cast<int>(pow(stencil_sz,NDIM)));
    std::vector<int> o_nnz(m_local,static_cast<int>(pow(stencil_sz,NDIM)));

    // Create an empty matrix.
    ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,
                           m_local, n_local,
                           PETSC_DETERMINE, PETSC_DETERMINE,
                           PETSC_DEFAULT, &d_nnz[0],
                           PETSC_DEFAULT, &o_nnz[0],
                           &mat); IBTK_CHKERRQ(ierr);

    // Set some general matrix options.
#ifdef DEBUG_CHECK_ASSERTIONS
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_LOCATION_ERR  , PETSC_TRUE); IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE); IBTK_CHKERRQ(ierr);
#endif

    // Set the matrix coefficients to correspond to the IB interpolation
    // operator.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = patch_level->getGridGeometry();
    const double* const XLower = grid_geom->getXLower();
    const double* const XUpper = grid_geom->getXUpper();
    const double* const dx0 = grid_geom->getDx();
    const IntVector<NDIM>& ratio = patch_level->getRatio();
    double dx[NDIM];
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        dx[d] = dx0[d] / double(ratio(d));
    }

    const BoxArray<NDIM>& domain_boxes = patch_level->getPhysicalDomain();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(domain_boxes.size() == 1);
#endif
    const Index<NDIM>& domain_lower = domain_boxes[0].lower();
    const Index<NDIM>& domain_upper = domain_boxes[0].upper();

    double* X_arr;
    ierr = VecGetArray(X_vec, &X_arr); IBTK_CHKERRQ(ierr);
    for (int k = 0; k < m_local/NDIM; ++k)
    {
        // Determine the index of the Cartesian grid cell containing the IB
        // point.
        const double* const X = &X_arr[NDIM*k];
        const Index<NDIM> idx = IndexUtilities::getCellIndex(X, XLower, XUpper, dx, domain_lower, domain_upper);

        // Determine the position of the center of the Cartesian grid cell
        // containing the IB point.
#if (NDIM == 2)
        const double X_cell[NDIM] = { (double(idx(0) - domain_lower(0))+0.5)*dx[0] + XLower[0] ,
                                      (double(idx(1) - domain_lower(1))+0.5)*dx[1] + XLower[1] };
#endif
#if (NDIM == 3)
        const double X_cell[NDIM] = { (double(idx(0) - domain_lower(0))+0.5)*dx[0] + XLower[0] ,
                                      (double(idx(1) - domain_lower(1))+0.5)*dx[1] + XLower[1] ,
                                      (double(idx(2) - domain_lower(2))+0.5)*dx[2] + XLower[2] };
#endif

        // Find the local patch that contains the IB point.
        const Box<NDIM> box(idx,idx);
        Array<int> patch_idx_arr;
        patch_level->getBoxTree()->findOverlapIndices(patch_idx_arr, box);
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(patch_idx_arr.size() == 1);
#endif
        const int patch_num = patch_idx_arr[0];
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(patch_num);
        Pointer<SideData<NDIM,int> > dof_master_index = patch->getPatchData(dof_master_index_idx);

        // Construct the interpolation for the current IB point.
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            // Compute the stencil box.
            Box<NDIM> stencil_box(idx,idx);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d == component_axis)
                {
                    stencil_box.lower()(d) -= (stencil_sz/2 - 1);
                    stencil_box.upper()(d) += (stencil_sz/2    );
                }
                else
                {
                    if (X[d] <= X_cell[d])
                    {
                        stencil_box.lower()(d) -= (stencil_sz/2    );
                        stencil_box.upper()(d) += (stencil_sz/2 - 1);
                    }
                    else
                    {
                        stencil_box.lower()(d) -= (stencil_sz/2 - 1);
                        stencil_box.upper()(d) += (stencil_sz/2    );
                    }
                }
            }
            const Index<NDIM>& stencil_box_lower = stencil_box.lower();

            // Compute the weights of the 1-dimensional delta functions.
            blitz::TinyVector<std::vector<double>,NDIM> w(std::vector<double>(stencil_sz,0.0));
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                for (int i = stencil_box_lower(d), j = 0; i <= stencil_box.upper()(d); ++i, ++j)
                {
                    const double X_grid = (double(i - domain_lower(d)) + (d == component_axis ? 0.0 : 0.5))*dx[d] + XLower[d];
                    w[d][j] = ib4_delta_fcn((X_grid - X[d])/dx[d]);
                }
            }

            // Compute the weights of the d-dimensional delta functions as the
            // tensor products of the 1-dimensional delta functions.
            const int stencil_box_nvals = stencil_box.size();
            std::vector<double> stencil_box_vals(stencil_box_nvals,1.0);
            std::vector<int> stencil_box_n_idxs(stencil_box_nvals);
            int stencil_idx = 0;
            for (Box<NDIM>::Iterator b(stencil_box); b; b++, stencil_idx++)
            {
                const SideIndex<NDIM> s_i(b(), component_axis, SideIndex<NDIM>::Lower);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    stencil_box_vals[stencil_idx] *= w[d][s_i(d) - stencil_box_lower(d)];
                }
                stencil_box_n_idxs[stencil_idx] = (*dof_master_index)(s_i);
            }

            int stencil_box_m_idx = i_lower + NDIM*k + component_axis;

            // Set the values.
            ierr = MatSetValues(mat, 1, &stencil_box_m_idx, stencil_box_nvals, &stencil_box_n_idxs[0], &stencil_box_vals[0], INSERT_VALUES); IBTK_CHKERRQ(ierr);
        }
    }
    ierr = VecRestoreArray(X_vec, &X_arr); IBTK_CHKERRQ(ierr);

    // Free the temporary data and assemble the matrix.
    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY); IBTK_CHKERRQ(ierr);
    patch_level->deallocatePatchData(dof_master_index_idx);
    var_db->removePatchDataIndex(dof_master_index_idx);
    ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY); IBTK_CHKERRQ(ierr);
    return;
}// constructPatchLevelInterpOp

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
PETScMatUtilities::constructBoxLaplaceOp(
    Mat& mat,
    const double C,
    const double D,
    const Box<NDIM>& src_ghost_box,
    const Box<NDIM>& dst_ghost_box,
    const Box<NDIM>& interior_box,
    const int data_depth,
    const double* const dx)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    static const int minimum_ghost_width = 1;
    TBOX_ASSERT(Box<NDIM>::grow(src_ghost_box, -minimum_ghost_width).contains(interior_box));
#endif

    int ierr;

    // Determine the shape and non-zero structure of the matrix.
    const int M = data_depth*dst_ghost_box.size();
    const int N = data_depth*src_ghost_box.size();

    BoxList<NDIM> dst_ghost_boxes(dst_ghost_box);
    dst_ghost_boxes.removeIntersections(interior_box);
    static const int stencil_sz = 2*NDIM+1;
    std::vector<int> nnz(M, stencil_sz);
    for (int d = 0; d < data_depth; ++d)
    {
        const int dst_offset = d*dst_ghost_box.size();
        for (BoxList<NDIM>::Iterator bl(dst_ghost_boxes); bl; bl++)
        {
            for (Box<NDIM>::Iterator b(bl()); b; b++)
            {
                nnz[dst_ghost_box.offset(b())+dst_offset] = 1;
            }
        }
    }

    // Create an empty matrix.
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, M, N, PETSC_DEFAULT, &nnz[0], &mat); IBTK_CHKERRQ(ierr);

    // Set some general matrix options.
#ifdef DEBUG_CHECK_ASSERTIONS
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_LOCATION_ERR  , PETSC_TRUE); IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE); IBTK_CHKERRQ(ierr);
#endif

    // Setup the finite difference stencil.  The stencil order is chosen to try
    // to optimize performance when setting the matrix coefficients.
    static const int x_axis = 0; (void) x_axis;
    static const int y_axis = 1; (void) y_axis;
    static const int z_axis = 2; (void) z_axis;
    blitz::TinyVector<int,NDIM> num_cells;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        num_cells[d] = src_ghost_box.numberCells(d);
    }
    std::vector<int> mat_stencil(stencil_sz);
#if (NDIM == 2)
    mat_stencil[0] = -num_cells[x_axis]; // ylower
    mat_stencil[1] = -1;                 // xlower
    mat_stencil[2] = 0;
    mat_stencil[3] = +1;                 // xupper
    mat_stencil[4] = +num_cells[x_axis]; // yupper
#endif
#if (NDIM == 3)
    mat_stencil[0] = -num_cells[x_axis]*num_cells[y_axis]; // zlower
    mat_stencil[1] = -num_cells[x_axis];                   // ylower
    mat_stencil[2] = -1;                                   // xlower
    mat_stencil[3] = 0;
    mat_stencil[4] = +1;                                   // xupper
    mat_stencil[5] = +num_cells[x_axis];                   // yupper
    mat_stencil[6] = +num_cells[x_axis]*num_cells[y_axis]; // zupper
#endif

    // Set the matrix coefficients to correspond to the standard finite
    // difference Laplacian discretization.
    //
    // Note that boundary conditions both at physical boundaries and at
    // coarse-fine interfaces are assumed to be treated implicitly by setting
    // ghost cell values.  This allows the matrix coefficients to be independent
    // of the boundary conditions.
    std::vector<double> mat_vals(stencil_sz,0.0);
    mat_vals[NDIM] = C;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const double& h = dx[axis];
        mat_vals[NDIM-axis-1] +=     D/(h*h);
        mat_vals[NDIM       ] -= 2.0*D/(h*h);
        mat_vals[NDIM+axis+1] +=     D/(h*h);
    }

    for (int d = 0; d < data_depth; ++d)
    {
        const int dst_offset = d*dst_ghost_box.size();
        const int src_offset = d*src_ghost_box.size();
        for (Box<NDIM>::Iterator b(interior_box); b; b++)
        {
            const Index<NDIM>& i = b();

            static const int m = 1;
            const int idxm = dst_ghost_box.offset(i)+dst_offset;

            static const int n = stencil_sz;
            std::vector<int> idxn(stencil_sz);
            const int stencil_shift = src_ghost_box.offset(i)+src_offset;
            std::transform(mat_stencil.begin(), mat_stencil.end(),
                           idxn.begin(), std::bind2nd(std::plus<int>(), stencil_shift));

            ierr = MatSetValues(mat, m, &idxm, n, &idxn[0], &mat_vals[0],
                                INSERT_VALUES); IBTK_CHKERRQ(ierr);
        }
    }

    // Set the entries in the ghost cell region of the dst data so that ghost
    // cell values are not modified by the application of the operator.
    for (int d = 0; d < data_depth; ++d)
    {
        const int dst_offset = d*dst_ghost_box.size();
        const int src_offset = d*src_ghost_box.size();
        for (BoxList<NDIM>::Iterator bl(dst_ghost_boxes); bl; bl++)
        {
            for (Box<NDIM>::Iterator b(bl()); b; b++)
            {
                const Index<NDIM>& i = b();
                const int idxm = dst_ghost_box.offset(i)+dst_offset;
                const int idxn = src_ghost_box.offset(i)+src_offset;
                ierr = MatSetValue(mat, idxm, idxn, 1.0, INSERT_VALUES); IBTK_CHKERRQ(ierr);
            }
        }
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY); IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd  (mat, MAT_FINAL_ASSEMBLY); IBTK_CHKERRQ(ierr);
    return;
}// constructBoxLaplaceOp

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
