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
#include <set>
#include <vector>

#include "Box.h"
#include "BoxArray.h"
#include "BoxTree.h"
#include "CartesianGridGeometry.h"
#include "CellData.h"
#include "CellGeometry.h"
#include "CellIndex.h"
#include "CoarseFineBoundary.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchLevel.h"
#include "PoissonSpecifications.h"
#include "RefineAlgorithm.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "boost/array.hpp"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/PETScMatUtilities.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "ibtk/PoissonUtilities.h"
#include "ibtk/SideSynchCopyFillPattern.h"
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

namespace
{
bool inline is_cf_bdry_idx(const Index<NDIM>& idx, const std::vector<Box<NDIM> >& cf_bdry_boxes)
{
    bool contains_idx = false;
    int n_cf_bdry_boxes = static_cast<int>(cf_bdry_boxes.size());
    for (int k = 0; !contains_idx || k < n_cf_bdry_boxes; ++k)
    {
        contains_idx = contains_idx || cf_bdry_boxes[k].contains(idx);
    }
    return contains_idx;
} // is_cf_bdry_idx

static const int LOWER = 0;
static const int UPPER = 1;
static const std::string CONSERVATIVE = "CONSERVATIVE";
static const std::string RT0 = "RT0";
static const std::string LINEAR = "LINEAR";

#define SCD(a) static_cast<double>(a)
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
PETScMatUtilities::constructPatchLevelCCLaplaceOp(Mat& mat,
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

void
PETScMatUtilities::constructPatchLevelCCLaplaceOp(Mat& mat,
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
                        0,
                        n_local ? &d_nnz[0] : NULL,
                        0,
                        n_local ? &o_nnz[0] : NULL,
                        &mat);
    IBTK_CHKERRQ(ierr);

    // Set block size.
    ierr = MatSetBlockSize(mat, depth);
    IBTK_CHKERRQ(ierr);

    // Set the matrix coefficients to correspond to the standard finite
    // difference approximation to the Laplacian.
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();

        // Compute matrix coefficients.
        const IntVector<NDIM> no_ghosts(0);
        CellData<NDIM, double> matrix_coefs(patch_box, stencil_sz * depth, no_ghosts);
        PoissonUtilities::computeMatrixCoefficients(matrix_coefs, patch, stencil, poisson_spec, bc_coefs, data_time);

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

void
PETScMatUtilities::constructPatchLevelSCLaplaceOp(Mat& mat,
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
                        0,
                        n_local ? &d_nnz[0] : NULL,
                        0,
                        n_local ? &o_nnz[0] : NULL,
                        &mat);
    IBTK_CHKERRQ(ierr);

    // Set the matrix coefficients to correspond to the standard finite
    // difference approximation to the Laplacian.
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();

        // Compute matrix coefficients.
        const IntVector<NDIM> no_ghosts(0);
        SideData<NDIM, double> matrix_coefs(patch_box, stencil_sz, no_ghosts);
        PoissonUtilities::computeMatrixCoefficients(matrix_coefs, patch, stencil, poisson_spec, bc_coefs, data_time);

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

void
PETScMatUtilities::constructPatchLevelSCInterpOp(Mat& mat,
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

    // The processor mapping determines which patches are assigned to which processors.
    const ProcessorMapping& proc_mapping = patch_level->getProcessorMapping();

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
        const Index<NDIM> X_idx = IndexUtilities::getCellIndex(X, grid_geom, ratio);

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
        bool found_local_patch = false;
        for (int growth_size = 0; growth_size <= 1; ++growth_size)
        {
            Box<NDIM> box(X_idx, X_idx);
            box.grow(IntVector<NDIM>(growth_size));
            Array<int> patch_num_arr;
            patch_level->getBoxTree()->findOverlapIndices(patch_num_arr, box);
            for (int j = 0; j < patch_num_arr.size() && !found_local_patch; ++j)
            {
                const int n = patch_num_arr[j];
                if (proc_mapping.isMappingLocal(n))
                {
                    patch_num[k] = n;
                    found_local_patch = true;
                }
            }
        }
#if !defined(NDEBUG)
        TBOX_ASSERT(found_local_patch);
#endif
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
                        0,
                        m_local ? &d_nnz[0] : NULL,
                        0,
                        m_local ? &o_nnz[0] : NULL,
                        &mat);
    IBTK_CHKERRQ(ierr);

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
                    (static_cast<double>(i - domain_lower(d)) + (d == axis ? 0.0 : 0.5)) * dx[d] + x_lower[d];
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

void
PETScMatUtilities::constructProlongationOp(Mat& mat,
                                           const std::string& op_type,
                                           int dof_index_idx,
                                           const std::vector<int>& num_fine_dofs_per_proc,
                                           const std::vector<int>& num_coarse_dofs_per_proc,
                                           Pointer<PatchLevel<NDIM> > fine_patch_level,
                                           Pointer<PatchLevel<NDIM> > coarse_patch_level,
                                           const AO& coarse_level_ao,
                                           const int coarse_ao_offset)
{
    // Determine the data-centering type.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > dof_index_var;
    var_db->mapIndexToVariable(dof_index_idx, dof_index_var);
    Pointer<CellVariable<NDIM, int> > dof_index_cc_var = dof_index_var;
    Pointer<SideVariable<NDIM, int> > dof_index_sc_var = dof_index_var;
    if (dof_index_cc_var)
    {
        if (op_type == CONSERVATIVE)
        {
            constructConservativeProlongationOp_cell(mat,
                                                     dof_index_idx,
                                                     num_fine_dofs_per_proc,
                                                     num_coarse_dofs_per_proc,
                                                     fine_patch_level,
                                                     coarse_patch_level,
                                                     coarse_level_ao,
                                                     coarse_ao_offset);
        }
        else
        {
            TBOX_ERROR(
                "PETScMatUtilities::constructProlongationOp(): Unsupported prolongation operator for cc-variable. "
                "Given operator is "
                << op_type
                << ". Supported ops are: "
                << CONSERVATIVE
                << std::endl);
        }
    }
    else if (dof_index_sc_var)
    {
        if (op_type == RT0)
        {
            constructRT0ProlongationOp_side(mat,
                                            dof_index_idx,
                                            num_fine_dofs_per_proc,
                                            num_coarse_dofs_per_proc,
                                            fine_patch_level,
                                            coarse_patch_level,
                                            coarse_level_ao,
                                            coarse_ao_offset);
        }
        else if (op_type == LINEAR)
        {
            constructLinearProlongationOp_side(mat,
                                               dof_index_idx,
                                               num_fine_dofs_per_proc,
                                               num_coarse_dofs_per_proc,
                                               fine_patch_level,
                                               coarse_patch_level,
                                               coarse_level_ao,
                                               coarse_ao_offset);
        }
        else
        {
            TBOX_ERROR(
                "PETScMatUtilities::constructProlongationOp(): Unsupported prolongation operator for sc-variable. "
                "Given operator is "
                << op_type
                << ". Supported ops are: "
                << RT0
                << " and "
                << LINEAR
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR("PETScVecUtilities::constructPatchLevelProlongationOp():\n"
                   << "  unsupported data centering type for variable "
                   << dof_index_var->getName()
                   << "\n");
    }

} // constructPatchLevelProlongationOp

void
PETScMatUtilities::constructRestrictionScalingOp(Mat& P, Vec& L)
{
    // Note that enteries of P are positive, so we will use column norm-1 function
    // of PETSc which appears to be faster than row sum call from the documentation.
    // We might have a column of all zeros for some DOFs (say pressure for combined
    // u-p vec, if we are only prolongating u). Therefore, care should be taken for
    // this case.

    int ierr;
    PetscInt M, N;
    ierr = MatGetSize(P, &M, &N);
    IBTK_CHKERRQ(ierr);
    std::vector<double> column_sum_inv(N);
    ierr = MatGetColumnNorms(P, NORM_1, &column_sum_inv[0]);
    IBTK_CHKERRQ(ierr);

    for (int k = 0; k < N; ++k)
    {
        const double sum = column_sum_inv[k];
        if (!MathUtilities<double>::equalEps(sum, 0.0))
        {
            column_sum_inv[k] = 1.0 / sum;
        }
        else
        {
            column_sum_inv[k] = 0.0;
        }
    }

    // Get the right vector of P, which becomes the left-vector of R.
    if (L)
    {
        ierr = VecDestroy(&L);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatCreateVecs(P, &L, PETSC_NULL);
    IBTK_CHKERRQ(ierr);
    PetscInt ilower, iupper, num_elems;
    ierr = VecGetOwnershipRange(L, &ilower, &iupper);
    IBTK_CHKERRQ(ierr);
    num_elems = iupper - ilower;
    if (num_elems != 0)
    {
        std::vector<PetscScalar> L_vals(num_elems);
        std::vector<PetscInt> L_idxs(num_elems);
        for (int k = ilower; k < iupper; ++k)
        {
            L_idxs[k - ilower] = k;
            L_vals[k - ilower] = column_sum_inv[k];
        }
        ierr = VecSetValues(L, num_elems, &L_idxs[0], &L_vals[0], INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
    }

    // Assemble the diagonal matrix.
    ierr = VecAssemblyBegin(L);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(L);
    IBTK_CHKERRQ(ierr);
    return;
} // constructRestrictionScalingOp

void
PETScMatUtilities::constructPatchLevelASMSubdomains(std::vector<IS>& is_overlap,
                                                    std::vector<IS>& is_nonoverlap,
                                                    const IntVector<NDIM>& box_size,
                                                    const IntVector<NDIM>& overlap_size,
                                                    const std::vector<int>& num_dofs_per_proc,
                                                    int dof_index_idx,
                                                    Pointer<PatchLevel<NDIM> > patch_level,
                                                    Pointer<CoarseFineBoundary<NDIM> > cf_boundary)
{
    int ierr;
    for (unsigned int k = 0; k < is_overlap.size(); ++k)
    {
        ierr = ISDestroy(&is_overlap[k]);
        IBTK_CHKERRQ(ierr);
    }
    is_overlap.clear();
    for (unsigned int k = 0; k < is_nonoverlap.size(); ++k)
    {
        ierr = ISDestroy(&is_nonoverlap[k]);
        IBTK_CHKERRQ(ierr);
    }
    is_nonoverlap.clear();

    // Determine the data-centering type.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > dof_index_var;
    var_db->mapIndexToVariable(dof_index_idx, dof_index_var);
    Pointer<CellVariable<NDIM, int> > dof_index_cc_var = dof_index_var;
    Pointer<SideVariable<NDIM, int> > dof_index_sc_var = dof_index_var;
    if (dof_index_cc_var)
    {
        constructPatchLevelASMSubdomains_cell(is_overlap,
                                              is_nonoverlap,
                                              box_size,
                                              overlap_size,
                                              num_dofs_per_proc,
                                              dof_index_idx,
                                              patch_level,
                                              cf_boundary);
    }
    else if (dof_index_sc_var)
    {
        constructPatchLevelASMSubdomains_side(is_overlap,
                                              is_nonoverlap,
                                              box_size,
                                              overlap_size,
                                              num_dofs_per_proc,
                                              dof_index_idx,
                                              patch_level,
                                              cf_boundary);
    }
    else
    {
        TBOX_ERROR("PETScVecUtilities::constructPatchLevelASMSubdomains():\n"
                   << "  unsupported data centering type for variable "
                   << dof_index_var->getName()
                   << "\n");
    }
    return;
} // constructPatchLevelASMSubdomains

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
PETScMatUtilities::constructConservativeProlongationOp_cell(Mat& mat,
                                                            int dof_index_idx,
                                                            const std::vector<int>& num_fine_dofs_per_proc,
                                                            const std::vector<int>& num_coarse_dofs_per_proc,
                                                            Pointer<PatchLevel<NDIM> > fine_patch_level,
                                                            Pointer<PatchLevel<NDIM> > coarse_patch_level,
                                                            const AO& coarse_level_ao,
                                                            const int coarse_ao_offset)
{
    int ierr;
    if (mat)
    {
        ierr = MatDestroy(&mat);
        IBTK_CHKERRQ(ierr);
    }

    // Determine the grid and data extents for the coarse level.
    const BoxArray<NDIM>& coarse_domain_boxes = coarse_patch_level->getPhysicalDomain();
#if !defined(NDEBUG)
    TBOX_ASSERT(coarse_domain_boxes.size() == 1);
#endif
    const Index<NDIM>& coarse_domain_lower = coarse_domain_boxes[0].lower();
    const Index<NDIM>& coarse_domain_upper = coarse_domain_boxes[0].upper();
    Index<NDIM> coarse_num_cells = 1;
    coarse_num_cells += coarse_domain_upper - coarse_domain_lower;

    // Ratio between fine and coarse levels.
    const IntVector<NDIM>& coarse_ratio = coarse_patch_level->getRatio();
    const IntVector<NDIM>& fine_ratio = fine_patch_level->getRatio();
    const IntVector<NDIM> fine_coarse_ratio = fine_ratio / coarse_ratio;

    // Determine the matrix dimensions and index ranges.
    const int mpi_rank = SAMRAI_MPI::getRank();
    const int m_local = num_fine_dofs_per_proc[mpi_rank];
    const int n_local = num_coarse_dofs_per_proc[mpi_rank];
    const int i_fine_lower =
        std::accumulate(num_fine_dofs_per_proc.begin(), num_fine_dofs_per_proc.begin() + mpi_rank, 0);
    const int i_fine_upper = i_fine_lower + m_local;
    const int j_coarse_lower =
        std::accumulate(num_coarse_dofs_per_proc.begin(), num_coarse_dofs_per_proc.begin() + mpi_rank, 0);
    const int j_coarse_upper = j_coarse_lower + n_local;

    // Determine the non-zero matrix structure for constant refine.
    std::vector<int> d_nnz(m_local, 0), o_nnz(m_local, 0);
    for (PatchLevel<NDIM>::Iterator p(fine_patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > fine_patch = fine_patch_level->getPatch(p());
        const Box<NDIM>& fine_patch_box = fine_patch->getBox();
        Pointer<CellData<NDIM, int> > dof_fine_data = fine_patch->getPatchData(dof_index_idx);
        const unsigned depth = dof_fine_data->getDepth();
        std::vector<int> samrai_petsc_map(depth), local_row(depth);

        for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(fine_patch_box)); b; b++)
        {
            const CellIndex<NDIM>& i_fine = b();
            for (unsigned d = 0; d < depth; ++d)
            {
                local_row[d] = (*dof_fine_data)(i_fine, d);
#if !defined(NDEBUG)
                TBOX_ASSERT(local_row[d] >= i_fine_lower && local_row[d] < i_fine_upper);
#endif
                local_row[d] -= i_fine_lower;
            }

            const CellIndex<NDIM> i_coarse = IndexUtilities::coarsen(i_fine, fine_coarse_ratio);
            for (unsigned d = 0; d < depth; ++d)
            {
                samrai_petsc_map[d] = IndexUtilities::mapIndexToInteger(
                    i_coarse, coarse_domain_lower, coarse_num_cells, d, coarse_ao_offset);
            }
            AOApplicationToPetsc(coarse_level_ao, depth, &samrai_petsc_map[0]);

            for (unsigned d = 0; d < depth; ++d)
            {
                if (samrai_petsc_map[d] >= j_coarse_lower && samrai_petsc_map[d] < j_coarse_upper)
                {
                    d_nnz[local_row[d]] = 1;
                }
                else
                {
                    o_nnz[local_row[d]] = 1;
                }
            }
        }
    }

    // Create an empty matrix.
    ierr = MatCreateAIJ(PETSC_COMM_WORLD,
                        m_local,
                        n_local,
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        0,
                        m_local ? &d_nnz[0] : NULL,
                        0,
                        m_local ? &o_nnz[0] : NULL,
                        &mat);
    IBTK_CHKERRQ(ierr);

    // Determine matrix-coefficients
    for (PatchLevel<NDIM>::Iterator p(fine_patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > fine_patch = fine_patch_level->getPatch(p());
        const Box<NDIM>& fine_patch_box = fine_patch->getBox();
        Pointer<CellData<NDIM, int> > dof_fine_data = fine_patch->getPatchData(dof_index_idx);
        const unsigned depth = dof_fine_data->getDepth();
        std::vector<int> samrai_petsc_map(depth);

        for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(fine_patch_box)); b; b++)
        {
            const CellIndex<NDIM>& i_fine = b();
            const CellIndex<NDIM> i_coarse = IndexUtilities::coarsen(i_fine, fine_coarse_ratio);

            for (unsigned d = 0; d < depth; ++d)
            {
                samrai_petsc_map[d] = IndexUtilities::mapIndexToInteger(
                    i_coarse, coarse_domain_lower, coarse_num_cells, d, coarse_ao_offset);
            }
            AOApplicationToPetsc(coarse_level_ao, depth, &samrai_petsc_map[0]);
            for (unsigned d = 0; d < depth; ++d)
            {
                int row = (*dof_fine_data)(i_fine, d);
                int col = samrai_petsc_map[d];
                PetscScalar val = 1.0;
                ierr = MatSetValues(mat, 1, &row, 1, &col, &val, INSERT_VALUES);
                IBTK_CHKERRQ(ierr);
            }
        }
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    return;

} // constructConservativeProlongationOp_cell

void
PETScMatUtilities::constructRT0ProlongationOp_side(Mat& mat,
                                                   int dof_index_idx,
                                                   const std::vector<int>& num_fine_dofs_per_proc,
                                                   const std::vector<int>& num_coarse_dofs_per_proc,
                                                   Pointer<PatchLevel<NDIM> > fine_patch_level,
                                                   Pointer<PatchLevel<NDIM> > coarse_patch_level,
                                                   const AO& coarse_level_ao,
                                                   const int coarse_ao_offset)
{
    int ierr;
    if (mat)
    {
        ierr = MatDestroy(&mat);
        IBTK_CHKERRQ(ierr);
    }

    // Determine the grid and data extents for the coarse level and fine levels.
    const BoxArray<NDIM>& coarse_domain_boxes = coarse_patch_level->getPhysicalDomain();
#if !defined(NDEBUG)
    TBOX_ASSERT(coarse_domain_boxes.size() == 1);
#endif
    const Index<NDIM>& coarse_domain_lower = coarse_domain_boxes[0].lower();
    const Index<NDIM>& coarse_domain_upper = coarse_domain_boxes[0].upper();
    Box<NDIM> coarse_domain_side_boxes[NDIM];
    for (int axis = 0; axis < NDIM; ++axis)
    {
        coarse_domain_side_boxes[axis] = SideGeometry<NDIM>::toSideBox(coarse_domain_boxes[0], axis);
    }
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = coarse_patch_level->getGridGeometry();
    IntVector<NDIM> coarse_periodic_shift = grid_geom->getPeriodicShift(coarse_patch_level->getRatio());
    boost::array<Index<NDIM>, NDIM> coarse_num_cells;
    for (unsigned d = 0; d < NDIM; ++d)
    {
        Index<NDIM> offset = 1;
        offset(d) = coarse_periodic_shift(d) ? 1 : 2;
        coarse_num_cells[d] = coarse_domain_upper - coarse_domain_lower + offset;
    }

    // Ratio between fine and coarse levels.
    const IntVector<NDIM>& coarse_ratio = coarse_patch_level->getRatio();
    const IntVector<NDIM>& fine_ratio = fine_patch_level->getRatio();
    const IntVector<NDIM> fine_coarse_ratio = fine_ratio / coarse_ratio;

    // Determine the matrix dimensions and index ranges.
    const int mpi_rank = SAMRAI_MPI::getRank();
    const int m_local = num_fine_dofs_per_proc[mpi_rank];
    const int n_local = num_coarse_dofs_per_proc[mpi_rank];
    const int i_fine_lower =
        std::accumulate(num_fine_dofs_per_proc.begin(), num_fine_dofs_per_proc.begin() + mpi_rank, 0);
    const int i_fine_upper = i_fine_lower + m_local;
    const int j_coarse_lower =
        std::accumulate(num_coarse_dofs_per_proc.begin(), num_coarse_dofs_per_proc.begin() + mpi_rank, 0);
    const int j_coarse_upper = j_coarse_lower + n_local;

    // Determine the non-zero matrix structure for the refine operator.
    std::vector<int> d_nnz(m_local, 0), o_nnz(m_local, 0);
    for (PatchLevel<NDIM>::Iterator p(fine_patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > fine_patch = fine_patch_level->getPatch(p());
        const Box<NDIM>& fine_patch_box = fine_patch->getBox();
        Pointer<SideData<NDIM, int> > fine_dof_data = fine_patch->getPatchData(dof_index_idx);
        const unsigned depth = fine_dof_data->getDepth();
        const int n_interpolants = 2;
        std::vector<int> samrai_petsc_map(n_interpolants * depth), local_row(depth);

        for (int axis = 0; axis < NDIM; ++axis)
        {
            IntVector<NDIM> offset = 0;
            offset(axis) = 1;

            int data_offset = 0;
            for (int side = 0; side < axis; ++side)
            {
                int side_offset = depth;
                for (unsigned d = 0; d < NDIM; ++d) side_offset *= coarse_num_cells[side](d);
                data_offset += side_offset;
            }

            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(fine_patch_box, axis)); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                bool on_proc_fine_loc = true;
                for (unsigned d = 0; d < depth; ++d)
                {
                    local_row[d] = (*fine_dof_data)(i_s, d);

                    on_proc_fine_loc = on_proc_fine_loc && local_row[d] >= i_fine_lower && local_row[d] < i_fine_upper;

                    local_row[d] -= i_fine_lower;
                }
                if (!on_proc_fine_loc) continue;

                const CellIndex<NDIM> I = IndexUtilities::coarsen(i, fine_coarse_ratio);
                const CellIndex<NDIM>& I_L = I;
                const CellIndex<NDIM> I_U = I_L + offset;
                for (unsigned d = 0; d < depth; ++d)
                {
                    samrai_petsc_map[d * n_interpolants] =
                        IndexUtilities::mapIndexToInteger(I_L,
                                                          coarse_domain_lower,
                                                          coarse_num_cells[axis],
                                                          d,
                                                          coarse_ao_offset + data_offset,
                                                          coarse_periodic_shift);

                    samrai_petsc_map[d * n_interpolants + 1] =
                        IndexUtilities::mapIndexToInteger(I_U,
                                                          coarse_domain_lower,
                                                          coarse_num_cells[axis],
                                                          d,
                                                          coarse_ao_offset + data_offset,
                                                          coarse_periodic_shift);
                }
                AOApplicationToPetsc(coarse_level_ao, n_interpolants * depth, &samrai_petsc_map[0]);
#if !defined(NDEBUG)
                for (unsigned d = 0; d < depth; ++d)
                {
                    if (samrai_petsc_map[d * n_interpolants] < 0)
                    {
                        int domain = IndexUtilities::mapIndexToInteger(I_L,
                                                                       coarse_domain_lower,
                                                                       coarse_num_cells[axis],
                                                                       d,
                                                                       coarse_ao_offset + data_offset,
                                                                       coarse_periodic_shift);
                        TBOX_ERROR("Component axis = " << axis << " with coarse grid index " << I_L
                                                       << " and SAMRAI mapping "
                                                       << domain
                                                       << " is mapped to "
                                                       << samrai_petsc_map[d * n_interpolants]
                                                       << " by AO object constructed for level "
                                                       << coarse_patch_level->getLevelNumber()
                                                       << ". Index "
                                                       << I
                                                       << " contained in coarse domain = "
                                                       << coarse_domain_side_boxes[axis].contains(I_L)
                                                       << std::endl);
                    }
                }
#endif

                for (unsigned d = 0; d < depth; ++d)
                {
                    if (samrai_petsc_map[d * n_interpolants] >= 0 &&
                        samrai_petsc_map[d * n_interpolants] >= j_coarse_lower && samrai_petsc_map[d] < j_coarse_upper)
                    {
                        d_nnz[local_row[d]] = 1;
                    }
                    else if (samrai_petsc_map[d * n_interpolants] >= 0)
                    {
                        o_nnz[local_row[d]] = 1;
                    }

                    if (samrai_petsc_map[d * n_interpolants + 1] >= 0 &&
                        samrai_petsc_map[d * n_interpolants + 1] >= j_coarse_lower &&
                        samrai_petsc_map[d * n_interpolants + 1] < j_coarse_upper)
                    {
                        d_nnz[local_row[d]] += 1;
                    }
                    else if (samrai_petsc_map[d * n_interpolants + 1] >= 0)
                    {
                        o_nnz[local_row[d]] += 1;
                    }
                }
            }
        }
    }

    // Create an empty matrix.
    ierr = MatCreateAIJ(PETSC_COMM_WORLD,
                        m_local,
                        n_local,
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        0,
                        m_local ? &d_nnz[0] : NULL,
                        0,
                        m_local ? &o_nnz[0] : NULL,
                        &mat);
    IBTK_CHKERRQ(ierr);

    // Determine the matrix-coefficients
    for (PatchLevel<NDIM>::Iterator p(fine_patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > fine_patch = fine_patch_level->getPatch(p());
        const Box<NDIM>& fine_patch_box = fine_patch->getBox();
        Pointer<SideData<NDIM, int> > fine_dof_data = fine_patch->getPatchData(dof_index_idx);
        const unsigned depth = fine_dof_data->getDepth();
        const int n_interpolants = 2;
        std::vector<int> samrai_petsc_map(n_interpolants * depth);

        for (int axis = 0; axis < NDIM; ++axis)
        {
            IntVector<NDIM> offset = 0;
            offset(axis) = 1;

            int data_offset = 0;
            for (int side = 0; side < axis; ++side)
            {
                int side_offset = depth;
                for (unsigned d = 0; d < NDIM; ++d) side_offset *= coarse_num_cells[side](d);
                data_offset += side_offset;
            }
            IntVector<NDIM> coarse_periodic_shift_axis = 0;
            coarse_periodic_shift_axis(axis) = coarse_periodic_shift(axis);

            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(fine_patch_box, axis)); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                bool on_proc_fine_loc = true;
                for (unsigned d = 0; d < depth; ++d)
                {
                    const int fine_dof_idx = (*fine_dof_data)(i_s, d);
                    on_proc_fine_loc = on_proc_fine_loc && fine_dof_idx >= i_fine_lower && fine_dof_idx < i_fine_upper;
                }
                if (!on_proc_fine_loc) continue;

                const CellIndex<NDIM> I = IndexUtilities::coarsen(i, fine_coarse_ratio);
                const CellIndex<NDIM>& I_L = I;
                const CellIndex<NDIM> I_U = I_L + offset;
                for (unsigned d = 0; d < depth; ++d)
                {
                    samrai_petsc_map[d * n_interpolants] =
                        IndexUtilities::mapIndexToInteger(I_L,
                                                          coarse_domain_lower,
                                                          coarse_num_cells[axis],
                                                          d,
                                                          coarse_ao_offset + data_offset,
                                                          coarse_periodic_shift);

                    samrai_petsc_map[d * n_interpolants + 1] =
                        IndexUtilities::mapIndexToInteger(I_U,
                                                          coarse_domain_lower,
                                                          coarse_num_cells[axis],
                                                          d,
                                                          coarse_ao_offset + data_offset,
                                                          coarse_periodic_shift);
                }
                AOApplicationToPetsc(coarse_level_ao, n_interpolants * depth, &samrai_petsc_map[0]);

                for (unsigned d = 0; d < depth; ++d)
                {
                    int row = (*fine_dof_data)(i_s, d);
                    std::vector<int> col;
                    std::vector<double> col_val(2);
                    for (int ip = 0; ip < n_interpolants; ++ip)
                    {
                        if (samrai_petsc_map[d * n_interpolants + ip] < 0) continue;
                        col.push_back(samrai_petsc_map[d * n_interpolants + ip]);
                    }
                    int col_size = static_cast<int>(col.size());
                    TBOX_ASSERT(col_size >= 1);

                    // w_L = 1 - [i(axis) - refine(I_L,ratio)(axis)]/ratio(axis)
                    double w_L = 1.0 -
                                 (static_cast<double>(i(axis)) -
                                  static_cast<double>(IndexUtilities::refine(I_L, fine_coarse_ratio)(axis))) /
                                     static_cast<double>(fine_coarse_ratio(axis));

                    col_val[0] = w_L;
                    col_val[1] = 1.0 - w_L;
                    if (col_size == 1)
                    {
                        col_val[0] += col_val[1];
                    }

                    ierr = MatSetValues(mat, 1, &row, col_size, &col[0], &col_val[0], INSERT_VALUES);
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
} // constructRT0ProlongationOp_side

void
PETScMatUtilities::constructLinearProlongationOp_side(Mat& mat,
                                                      int dof_index_idx,
                                                      const std::vector<int>& num_fine_dofs_per_proc,
                                                      const std::vector<int>& num_coarse_dofs_per_proc,
                                                      Pointer<PatchLevel<NDIM> > fine_patch_level,
                                                      Pointer<PatchLevel<NDIM> > coarse_patch_level,
                                                      const AO& coarse_level_ao,
                                                      const int coarse_ao_offset)
{
    int ierr;
    if (mat)
    {
        ierr = MatDestroy(&mat);
        IBTK_CHKERRQ(ierr);
    }

    // Determine the grid and data extents for the coarse level and fine levels.
    const BoxArray<NDIM>& coarse_domain_boxes = coarse_patch_level->getPhysicalDomain();
#if !defined(NDEBUG)
    TBOX_ASSERT(coarse_domain_boxes.size() == 1);
#endif
    const Index<NDIM>& coarse_domain_lower = coarse_domain_boxes[0].lower();
    const Index<NDIM>& coarse_domain_upper = coarse_domain_boxes[0].upper();
    Box<NDIM> coarse_domain_side_boxes[NDIM];
    for (int axis = 0; axis < NDIM; ++axis)
    {
        coarse_domain_side_boxes[axis] = SideGeometry<NDIM>::toSideBox(coarse_domain_boxes[0], axis);
    }
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = coarse_patch_level->getGridGeometry();
    IntVector<NDIM> coarse_periodic_shift = grid_geom->getPeriodicShift(coarse_patch_level->getRatio());
    boost::array<Index<NDIM>, NDIM> coarse_num_cells;
    for (unsigned d = 0; d < NDIM; ++d)
    {
        Index<NDIM> offset = 1;
        offset(d) = coarse_periodic_shift(d) ? 1 : 2;
        coarse_num_cells[d] = coarse_domain_upper - coarse_domain_lower + offset;
    }

    // Ratio between fine and coarse levels.
    const IntVector<NDIM>& coarse_ratio = coarse_patch_level->getRatio();
    const IntVector<NDIM>& fine_ratio = fine_patch_level->getRatio();
    const IntVector<NDIM> fine_coarse_ratio = fine_ratio / coarse_ratio;

    // Determine the matrix dimensions and index ranges.
    const int mpi_rank = SAMRAI_MPI::getRank();
    const int m_local = num_fine_dofs_per_proc[mpi_rank];
    const int n_local = num_coarse_dofs_per_proc[mpi_rank];
    const int i_fine_lower =
        std::accumulate(num_fine_dofs_per_proc.begin(), num_fine_dofs_per_proc.begin() + mpi_rank, 0);
    const int i_fine_upper = i_fine_lower + m_local;
    const int j_coarse_lower =
        std::accumulate(num_coarse_dofs_per_proc.begin(), num_coarse_dofs_per_proc.begin() + mpi_rank, 0);
    const int j_coarse_upper = j_coarse_lower + n_local;

    // Determine the non-zero matrix structure for the refine operator.
    std::vector<int> d_nnz(m_local, 0), o_nnz(m_local, 0);
    for (PatchLevel<NDIM>::Iterator p(fine_patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > fine_patch = fine_patch_level->getPatch(p());
        const Box<NDIM>& fine_patch_box = fine_patch->getBox();
        Pointer<SideData<NDIM, int> > fine_dof_data = fine_patch->getPatchData(dof_index_idx);
        const unsigned depth = fine_dof_data->getDepth();
        const int n_interpolants = 4 * (NDIM - 1);
        std::vector<int> samrai_petsc_map(n_interpolants * depth), local_row(depth);

        for (int axis = 0; axis < NDIM; ++axis)
        {
            int data_offset = 0;
            for (int side = 0; side < axis; ++side)
            {
                int side_offset = depth;
                for (unsigned d = 0; d < NDIM; ++d) side_offset *= coarse_num_cells[side](d);
                data_offset += side_offset;
            }

            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(fine_patch_box, axis)); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                bool on_proc_fine_loc = true;
                for (unsigned d = 0; d < depth; ++d)
                {
                    local_row[d] = (*fine_dof_data)(i_s, d);

                    on_proc_fine_loc = on_proc_fine_loc && local_row[d] >= i_fine_lower && local_row[d] < i_fine_upper;

                    local_row[d] -= i_fine_lower;
                }
                if (!on_proc_fine_loc) continue;

                const CellIndex<NDIM> I = IndexUtilities::coarsen(i, fine_coarse_ratio);
                const CellIndex<NDIM> i_lower = IndexUtilities::refine(I, fine_coarse_ratio);

                std::vector<CellIndex<NDIM> > interpolants(n_interpolants);
                std::vector<IntVector<NDIM> > offsets(n_interpolants, 0);
                int upperlower[NDIM];
                for (int side = 0; side < NDIM; ++side)
                {
                    if (side == axis)
                    {
                        upperlower[side] = 1;
                    }
                    else
                    {
                        upperlower[side] = i(side) - i_lower(side) >= fine_coarse_ratio(side) / 2 ? 1 : -1;
                    }
                }

                if (axis == 0)
                {
                    /***************************************
                     *        2 -- 3          6 -- 7
                     *     y  |    |  ---> z  |    |
                     *        0 -- 1          4 -- 5
                     *           x
                     ***************************************/
                    offsets[1](0) = upperlower[0];

                    offsets[2](1) = upperlower[1];

                    offsets[3](0) = upperlower[0];
                    offsets[3](1) = upperlower[1];
#if (NDIM == 3)
                    offsets[4](2) = upperlower[2];

                    offsets[5](0) = upperlower[0];
                    offsets[5](2) = upperlower[2];

                    offsets[6](1) = upperlower[1];
                    offsets[6](2) = upperlower[2];

                    offsets[7](0) = upperlower[0];
                    offsets[7](1) = upperlower[1];
                    offsets[7](2) = upperlower[2];
#endif
                }
                else if (axis == 1)
                {
                    /************************************
                     *        1 -- 3           5 -- 7
                     *    y   |    |  ---> z   |    |
                     *        0 -- 2           4 -- 6
                     *           x
                     ************************************/

                    offsets[1](1) = upperlower[1];

                    offsets[2](0) = upperlower[0];

                    offsets[3](0) = upperlower[0];
                    offsets[3](1) = upperlower[1];
#if (NDIM == 3)
                    offsets[4](2) = upperlower[2];

                    offsets[5](1) = upperlower[1];
                    offsets[5](2) = upperlower[2];

                    offsets[6](0) = upperlower[0];
                    offsets[6](2) = upperlower[2];

                    offsets[7](0) = upperlower[0];
                    offsets[7](1) = upperlower[1];
                    offsets[7](2) = upperlower[2];
#endif
                }
#if (NDIM == 3)
                else if (axis == 2)
                {
                    /***********************************
                     *        1 -- 3          5 -- 7
                     *     z  |    |  ---> y  |    |
                     *        0 -- 2          4 -- 6
                     *           x
                     **********************************/

                    offsets[1](2) = upperlower[2];

                    offsets[2](0) = upperlower[0];

                    offsets[3](0) = upperlower[0];
                    offsets[3](2) = upperlower[2];

                    offsets[4](1) = upperlower[1];

                    offsets[5](1) = upperlower[1];
                    offsets[5](2) = upperlower[2];

                    offsets[6](0) = upperlower[0];
                    offsets[6](1) = upperlower[1];

                    offsets[7](0) = upperlower[0];
                    offsets[7](1) = upperlower[1];
                    offsets[7](2) = upperlower[2];
                }
#endif

                for (int ip = 0; ip < n_interpolants; ++ip)
                {
                    interpolants[ip] = I + offsets[ip];
                }

                for (unsigned d = 0; d < depth; ++d)
                {
                    for (int ip = 0; ip < n_interpolants; ++ip)
                    {
                        samrai_petsc_map[d * n_interpolants + ip] =
                            IndexUtilities::mapIndexToInteger(interpolants[ip],
                                                              coarse_domain_lower,
                                                              coarse_num_cells[axis],
                                                              d,
                                                              coarse_ao_offset + data_offset,
                                                              coarse_periodic_shift);
                    }
                }
                AOApplicationToPetsc(coarse_level_ao, n_interpolants * depth, &samrai_petsc_map[0]);
#if !defined(NDEBUG)
                for (unsigned d = 0; d < depth; ++d)
                {
                    if (samrai_petsc_map[d * n_interpolants] < 0)
                    {
                        int domain = IndexUtilities::mapIndexToInteger(interpolants[0],
                                                                       coarse_domain_lower,
                                                                       coarse_num_cells[axis],
                                                                       d,
                                                                       coarse_ao_offset + data_offset,
                                                                       coarse_periodic_shift);
                        TBOX_ERROR("Component axis = " << axis << " with coarse grid index " << I
                                                       << " and SAMRAI mapping "
                                                       << domain
                                                       << " is mapped to "
                                                       << samrai_petsc_map[d * n_interpolants]
                                                       << " by AO object constructed for level "
                                                       << coarse_patch_level->getLevelNumber()
                                                       << ". Index "
                                                       << I
                                                       << " contained in coarse domain = "
                                                       << coarse_domain_side_boxes[axis].contains(interpolants[0])
                                                       << std::endl);
                    }
                }
#endif

                for (unsigned d = 0; d < depth; ++d)
                {
                    for (int ip = 0; ip < n_interpolants; ++ip)
                    {
                        const int idx = d * n_interpolants + ip;
                        if (samrai_petsc_map[idx] < 0) continue;

                        if (samrai_petsc_map[idx] >= j_coarse_lower && samrai_petsc_map[idx] < j_coarse_upper)
                            d_nnz[local_row[d]] += 1;
                        else
                            o_nnz[local_row[d]] += 1;
                    }
                }
            }
        }
    }

    // Create an empty matrix.
    ierr = MatCreateAIJ(PETSC_COMM_WORLD,
                        m_local,
                        n_local,
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        0,
                        m_local ? &d_nnz[0] : NULL,
                        0,
                        m_local ? &o_nnz[0] : NULL,
                        &mat);
    IBTK_CHKERRQ(ierr);

    // Determine the matrix-coefficients
    for (PatchLevel<NDIM>::Iterator p(fine_patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > fine_patch = fine_patch_level->getPatch(p());
        const Box<NDIM>& fine_patch_box = fine_patch->getBox();
        Pointer<SideData<NDIM, int> > fine_dof_data = fine_patch->getPatchData(dof_index_idx);
        const unsigned depth = fine_dof_data->getDepth();
        const int n_interpolants = 4 * (NDIM - 1);
        std::vector<int> samrai_petsc_map(n_interpolants * depth);

        for (int axis = 0; axis < NDIM; ++axis)
        {
            int data_offset = 0;
            for (int side = 0; side < axis; ++side)
            {
                int side_offset = depth;
                for (unsigned d = 0; d < NDIM; ++d) side_offset *= coarse_num_cells[side](d);
                data_offset += side_offset;
            }

            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(fine_patch_box, axis)); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                bool on_proc_fine_loc = true;
                for (unsigned d = 0; d < depth; ++d)
                {
                    const int fine_dof_idx = (*fine_dof_data)(i_s, d);
                    on_proc_fine_loc = on_proc_fine_loc && fine_dof_idx >= i_fine_lower && fine_dof_idx < i_fine_upper;
                }
                if (!on_proc_fine_loc) continue;

                const CellIndex<NDIM> I = IndexUtilities::coarsen(i, fine_coarse_ratio);
                const CellIndex<NDIM> i_lower = IndexUtilities::refine(I, fine_coarse_ratio);

                std::vector<CellIndex<NDIM> > interpolants(n_interpolants);
                std::vector<IntVector<NDIM> > offsets(n_interpolants, 0);
                int upperlower[NDIM];
                for (int side = 0; side < NDIM; ++side)
                {
                    if (side == axis)
                    {
                        upperlower[side] = 1;
                    }
                    else
                    {
                        upperlower[side] = i(side) - i_lower(side) >= fine_coarse_ratio(side) / 2 ? 1 : -1;
                    }
                }

                if (axis == 0)
                {
                    /***************************************
                     *        2 -- 3          6 -- 7
                     *     y  |    |  ---> z  |    |
                     *        0 -- 1          4 -- 5
                     *           x
                     ***************************************/
                    offsets[1](0) = upperlower[0];

                    offsets[2](1) = upperlower[1];

                    offsets[3](0) = upperlower[0];
                    offsets[3](1) = upperlower[1];
#if (NDIM == 3)
                    offsets[4](2) = upperlower[2];

                    offsets[5](0) = upperlower[0];
                    offsets[5](2) = upperlower[2];

                    offsets[6](1) = upperlower[1];
                    offsets[6](2) = upperlower[2];

                    offsets[7](0) = upperlower[0];
                    offsets[7](1) = upperlower[1];
                    offsets[7](2) = upperlower[2];
#endif
                }
                else if (axis == 1)
                {
                    /************************************
                     *        1 -- 3           5 -- 7
                     *    y   |    |  ---> z   |    |
                     *        0 -- 2           4 -- 6
                     *           x
                     ************************************/

                    offsets[1](1) = upperlower[1];

                    offsets[2](0) = upperlower[0];

                    offsets[3](0) = upperlower[0];
                    offsets[3](1) = upperlower[1];
#if (NDIM == 3)
                    offsets[4](2) = upperlower[2];

                    offsets[5](1) = upperlower[1];
                    offsets[5](2) = upperlower[2];

                    offsets[6](0) = upperlower[0];
                    offsets[6](2) = upperlower[2];

                    offsets[7](0) = upperlower[0];
                    offsets[7](1) = upperlower[1];
                    offsets[7](2) = upperlower[2];
#endif
                }
#if (NDIM == 3)
                else if (axis == 2)
                {
                    /***********************************
                     *        1 -- 3          5 -- 7
                     *     z  |    |  ---> y  |    |
                     *        0 -- 2          4 -- 6
                     *           x
                     **********************************/

                    offsets[1](2) = upperlower[2];

                    offsets[2](0) = upperlower[0];

                    offsets[3](0) = upperlower[0];
                    offsets[3](2) = upperlower[2];

                    offsets[4](1) = upperlower[1];

                    offsets[5](1) = upperlower[1];
                    offsets[5](2) = upperlower[2];

                    offsets[6](0) = upperlower[0];
                    offsets[6](1) = upperlower[1];

                    offsets[7](0) = upperlower[0];
                    offsets[7](1) = upperlower[1];
                    offsets[7](2) = upperlower[2];
                }
#endif

                for (int ip = 0; ip < n_interpolants; ++ip)
                {
                    interpolants[ip] = I + offsets[ip];
                }

                for (unsigned d = 0; d < depth; ++d)
                {
                    for (int ip = 0; ip < n_interpolants; ++ip)
                    {
                        samrai_petsc_map[d * n_interpolants + ip] =
                            IndexUtilities::mapIndexToInteger(interpolants[ip],
                                                              coarse_domain_lower,
                                                              coarse_num_cells[axis],
                                                              d,
                                                              coarse_ao_offset + data_offset,
                                                              coarse_periodic_shift);
                    }
                }
                AOApplicationToPetsc(coarse_level_ao, n_interpolants * depth, &samrai_petsc_map[0]);

                for (unsigned d = 0; d < depth; ++d)
                {
                    // Interpolation weights in Cartesian axis.
                    double w[NDIM];
                    if (axis == 0)
                    {
                        w[0] = 1.0 -
                               (SCD(i(0)) - SCD(IndexUtilities::refine(interpolants[0], fine_coarse_ratio)(0))) /
                                   SCD(fine_coarse_ratio(0));

                        w[1] = 1.0 -
                               (0.5 + SCD(i(1)) -
                                SCD(IndexUtilities::refine(IntVector<NDIM>::min(interpolants[0], interpolants[2]),
                                                           fine_coarse_ratio)(1)) -
                                SCD(fine_coarse_ratio(1) / 2.0)) /
                                   SCD(fine_coarse_ratio(1));
#if (NDIM == 3)

                        w[2] = 1.0 -
                               (0.5 + SCD(i(2)) -
                                SCD(IndexUtilities::refine(IntVector<NDIM>::min(interpolants[0], interpolants[4]),
                                                           fine_coarse_ratio)(2)) -
                                SCD(fine_coarse_ratio(2) / 2.0)) /
                                   SCD(fine_coarse_ratio(2));
#endif
                    }
                    else if (axis == 1)
                    {
                        w[1] = 1.0 -
                               (SCD(i(1)) - SCD(IndexUtilities::refine(interpolants[0], fine_coarse_ratio)(1))) /
                                   SCD(fine_coarse_ratio(1));

                        w[0] = 1.0 -
                               (0.5 + SCD(i(0)) -
                                SCD(IndexUtilities::refine(IntVector<NDIM>::min(interpolants[0], interpolants[2]),
                                                           fine_coarse_ratio)(0)) -
                                SCD(fine_coarse_ratio(0) / 2.0)) /
                                   SCD(fine_coarse_ratio(0));
#if (NDIM == 3)

                        w[2] = 1.0 -
                               (0.5 + SCD(i(2)) -
                                SCD(IndexUtilities::refine(IntVector<NDIM>::min(interpolants[0], interpolants[4]),
                                                           fine_coarse_ratio)(2)) -
                                SCD(fine_coarse_ratio(2) / 2.0)) /
                                   SCD(fine_coarse_ratio(2));
#endif
                    }
#if (NDIM == 3)
                    else if (axis == 2)
                    {
                        w[2] = 1.0 -
                               (SCD(i(2)) - SCD(IndexUtilities::refine(interpolants[0], fine_coarse_ratio)(2))) /
                                   SCD(fine_coarse_ratio(2));

                        w[0] = 1.0 -
                               (0.5 + SCD(i(0)) -
                                SCD(IndexUtilities::refine(IntVector<NDIM>::min(interpolants[0], interpolants[2]),
                                                           fine_coarse_ratio)(0)) -
                                SCD(fine_coarse_ratio(0) / 2.0)) /
                                   SCD(fine_coarse_ratio(0));

                        w[1] = 1.0 -
                               (0.5 + SCD(i(1)) -
                                SCD(IndexUtilities::refine(IntVector<NDIM>::min(interpolants[0], interpolants[4]),
                                                           fine_coarse_ratio)(1)) -
                                SCD(fine_coarse_ratio(1) / 2.0)) /
                                   SCD(fine_coarse_ratio(1));
                    }
#endif

                    int row = (*fine_dof_data)(i_s, d);
                    std::vector<int> col;
                    std::vector<double> col_val;
                    int n_cols = 0;
                    for (int ip = 0; ip < n_interpolants; ++ip)
                    {
                        if (samrai_petsc_map[d * n_interpolants + ip] >= 0)
                        {
                            ++n_cols;
                        }
                    }
                    col.resize(n_cols);
                    col_val.resize(n_cols);

                    // Find weights.
                    double w0, w1, w2 = 1.0;
                    double values[8];
                    if (axis == 0)
                    {
                        w0 = w[0];
                        w1 = interpolants[0](1) < interpolants[2](1) ? w[1] : 1.0 - w[1];
#if (NDIM == 3)
                        w2 = interpolants[0](2) < interpolants[4](2) ? w[2] : 1.0 - w[2];
#endif

                        values[0] = w0 * w1 * w2;
                        values[1] = (1.0 - w0) * w1 * w2;
                        values[2] = w0 * (1.0 - w1) * w2;
                        if (samrai_petsc_map[d * n_interpolants + 2] < 0)
                        {
                            values[0] += values[2];
                        }
                        values[3] = (1.0 - w0) * (1.0 - w1) * w2;
                        if (samrai_petsc_map[d * n_interpolants + 3] < 0)
                        {
                            values[1] += values[3];
                        }

#if (NDIM == 3)
                        values[4] = w0 * w1 * (1.0 - w2);
                        values[5] = (1.0 - w0) * w1 * (1.0 - w2);
                        values[6] = w0 * (1.0 - w1) * (1.0 - w2);
                        if (samrai_petsc_map[d * n_interpolants + 6] < 0)
                        {
                            values[4] += values[6];
                        }
                        values[7] = (1.0 - w0) * (1.0 - w1) * (1.0 - w2);
                        if (samrai_petsc_map[d * n_interpolants + 7] < 0)
                        {
                            values[5] += values[7];
                        }
                        if (samrai_petsc_map[d * n_interpolants + 4] < 0)
                        {
                            values[0] += values[4];
                        }
                        if (samrai_petsc_map[d * n_interpolants + 5] < 0)
                        {
                            values[1] += values[5];
                        }
#endif
                        if (samrai_petsc_map[d * n_interpolants + 1] < 0)
                        {
                            values[0] += values[1];
                        }

                        int k = 0;
                        for (int ip = 0; ip < n_interpolants; ++ip)
                        {
                            if (samrai_petsc_map[d * n_interpolants + ip] >= 0)
                            {
                                col[k] = samrai_petsc_map[d * n_interpolants + ip];
                                col_val[k] = values[ip];
                                k = k + 1;
                            }
                        }
                        TBOX_ASSERT(k == n_cols);
                    }
                    else if (axis == 1)
                    {
                        w0 = interpolants[0](0) < interpolants[2](0) ? w[0] : 1.0 - w[0];
                        w1 = w[1];
#if (NDIM == 3)
                        w2 = interpolants[0](2) < interpolants[4](2) ? w[2] : 1.0 - w[2];
#endif

                        values[0] = w0 * w1 * w2;
                        values[1] = w0 * (1.0 - w1) * w2;
                        values[2] = (1.0 - w0) * w1 * w2;
                        if (samrai_petsc_map[d * n_interpolants + 2] < 0)
                        {
                            values[0] += values[2];
                        }
                        values[3] = (1.0 - w0) * (1.0 - w1) * w2;
                        if (samrai_petsc_map[d * n_interpolants + 3] < 0)
                        {
                            values[1] += values[3];
                        }

#if (NDIM == 3)
                        values[4] = w0 * w1 * (1.0 - w2);
                        values[5] = w0 * (1.0 - w1) * (1.0 - w2);
                        values[6] = (1.0 - w0) * w1 * (1.0 - w2);
                        if (samrai_petsc_map[d * n_interpolants + 6] < 0)
                        {
                            values[4] += values[6];
                        }
                        values[7] = (1.0 - w0) * (1.0 - w1) * (1.0 - w2);
                        if (samrai_petsc_map[d * n_interpolants + 7] < 0)
                        {
                            values[5] += values[7];
                        }
                        if (samrai_petsc_map[d * n_interpolants + 4] < 0)
                        {
                            values[0] += values[4];
                        }
                        if (samrai_petsc_map[d * n_interpolants + 5] < 0)
                        {
                            values[1] += values[5];
                        }
#endif
                        if (samrai_petsc_map[d * n_interpolants + 1] < 0)
                        {
                            values[0] += values[1];
                        }

                        int k = 0;
                        for (int ip = 0; ip < n_interpolants; ++ip)
                        {
                            if (samrai_petsc_map[d * n_interpolants + ip] >= 0)
                            {
                                col[k] = samrai_petsc_map[d * n_interpolants + ip];
                                col_val[k] = values[ip];
                                k = k + 1;
                            }
                        }
                        TBOX_ASSERT(k == n_cols);
                    }
#if (NDIM == 3)
                    else if (axis == 2)
                    {
                        w0 = interpolants[0](0) < interpolants[2](0) ? w[0] : 1.0 - w[0];
                        w1 = interpolants[0](1) < interpolants[4](1) ? w[1] : 1.0 - w[1];
                        w2 = w[2];

                        values[0] = w0 * w1 * w2;
                        values[1] = w0 * w1 * (1.0 - w2);
                        values[2] = (1.0 - w0) * w1 * w2;
                        if (samrai_petsc_map[d * n_interpolants + 2] < 0)
                        {
                            values[0] += values[2];
                        }
                        values[3] = (1.0 - w0) * w1 * (1.0 - w2);
                        if (samrai_petsc_map[d * n_interpolants + 3] < 0)
                        {
                            values[1] += values[3];
                        }

                        values[4] = w0 * (1.0 - w1) * w2;
                        values[5] = w0 * (1.0 - w1) * (1.0 - w2);
                        values[6] = (1.0 - w0) * (1.0 - w1) * w2;
                        if (samrai_petsc_map[d * n_interpolants + 6] < 0)
                        {
                            values[4] += values[6];
                        }
                        values[7] = (1.0 - w0) * (1.0 - w1) * (1.0 - w2);
                        if (samrai_petsc_map[d * n_interpolants + 7] < 0)
                        {
                            values[5] += values[7];
                        }
                        if (samrai_petsc_map[d * n_interpolants + 4] < 0)
                        {
                            values[0] += values[4];
                        }
                        if (samrai_petsc_map[d * n_interpolants + 5] < 0)
                        {
                            values[1] += values[5];
                        }
                        if (samrai_petsc_map[d * n_interpolants + 1] < 0)
                        {
                            values[0] += values[1];
                        }

                        int k = 0;
                        for (int ip = 0; ip < n_interpolants; ++ip)
                        {
                            if (samrai_petsc_map[d * n_interpolants + ip] >= 0)
                            {
                                col[k] = samrai_petsc_map[d * n_interpolants + ip];
                                col_val[k] = values[ip];
                                k = k + 1;
                            }
                        }
                        TBOX_ASSERT(k == n_cols);
                    }
#endif
                    ierr = MatSetValues(mat, 1, &row, n_cols, &col[0], &col_val[0], INSERT_VALUES);
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
} // constructLinearProlongationOp_side

void
PETScMatUtilities::constructPatchLevelASMSubdomains_cell(std::vector<IS>& is_overlap,
                                                         std::vector<IS>& is_nonoverlap,
                                                         const IntVector<NDIM>& box_size,
                                                         const IntVector<NDIM>& overlap_size,
                                                         const std::vector<int>& /*num_dofs_per_proc*/,
                                                         int dof_index_idx,
                                                         Pointer<PatchLevel<NDIM> > patch_level,
                                                         Pointer<CoarseFineBoundary<NDIM> > /*cf_boundary*/)
{
    // Check if there is an overlap.
    const bool there_is_overlap = overlap_size.max();

    // Determine the subdomains associated with this processor.
    const int n_local_patches = patch_level->getProcessorMapping().getNumberOfLocalIndices();
    std::vector<std::vector<Box<NDIM> > > overlap_boxes(n_local_patches), nonoverlap_boxes(n_local_patches);
    int patch_counter = 0, subdomain_counter = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++, ++patch_counter)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        IndexUtilities::partitionPatchBox(
            overlap_boxes[patch_counter], nonoverlap_boxes[patch_counter], patch_box, box_size, overlap_size);
        subdomain_counter += overlap_boxes[patch_counter].size();
    }
    is_overlap.resize(subdomain_counter);
    is_nonoverlap.resize(subdomain_counter);

    // Fill in the IS'es.
    patch_counter = 0, subdomain_counter = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++, ++patch_counter)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<CellData<NDIM, int> > dof_data = patch->getPatchData(dof_index_idx);
        const int data_depth = dof_data->getDepth();
#if !defined(NDEBUG)
        TBOX_ASSERT(dof_data->getGhostCellWidth().min() >= overlap_size.max());
#endif
        size_t n_patch_subdomains = overlap_boxes[patch_counter].size();
        for (size_t i = 0; i < n_patch_subdomains; ++i, ++subdomain_counter)
        {
            // The nonoverlapping subdomains.
            const Box<NDIM>& box_local = nonoverlap_boxes[patch_counter][i];
            std::set<int> box_local_dofs;
            for (Box<NDIM>::Iterator b(box_local); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                for (int d = 0; d < data_depth; ++d)
                {
                    box_local_dofs.insert((*dof_data)(i, d));
                }
            }
            const int n_local = static_cast<int>(box_local_dofs.size());
            PetscInt* box_local_dof_arr;
            PetscMalloc1(n_local, &box_local_dof_arr);
            std::copy(box_local_dofs.begin(), box_local_dofs.end(), box_local_dof_arr);
            ISCreateGeneral(
                PETSC_COMM_SELF, n_local, box_local_dof_arr, PETSC_OWN_POINTER, &is_nonoverlap[subdomain_counter]);

            // The overlapping subdomains.
            if (!there_is_overlap)
            {
                PetscObjectReference(reinterpret_cast<PetscObject>(is_nonoverlap[subdomain_counter]));
                is_overlap[subdomain_counter] = is_nonoverlap[subdomain_counter];
            }
            else
            {
                const Box<NDIM>& box_overlap = overlap_boxes[patch_counter][i];
                std::set<int> box_overlap_dofs;
                for (Box<NDIM>::Iterator b(box_overlap); b; b++)
                {
                    const CellIndex<NDIM>& i = b();
                    for (int d = 0; d < data_depth; ++d)
                    {
                        // We keep only those DOFs that are inside the physical
                        // domain and away from c-f interface.  Some of the DOFs
                        // may be on other processors.  Cell-centered DOFs can
                        // never lie on boundaries.
                        const int dof_idx = (*dof_data)(i, d);
                        if (dof_idx >= 0)
                        {
                            box_overlap_dofs.insert(dof_idx);
                        }
                    }
                }
                const int n_overlap = static_cast<int>(box_overlap_dofs.size());
                PetscInt* box_overlap_dof_arr;
                PetscMalloc1(n_overlap, &box_overlap_dof_arr);
                std::copy(box_overlap_dofs.begin(), box_overlap_dofs.end(), box_overlap_dof_arr);
                ISCreateGeneral(
                    PETSC_COMM_SELF, n_overlap, box_overlap_dof_arr, PETSC_OWN_POINTER, &is_overlap[subdomain_counter]);
            }
        }
    }
    return;
} // constructPatchLevelASMSubdomains_cell

void
PETScMatUtilities::constructPatchLevelASMSubdomains_side(std::vector<IS>& is_overlap,
                                                         std::vector<IS>& is_nonoverlap,
                                                         const IntVector<NDIM>& box_size,
                                                         const IntVector<NDIM>& overlap_size,
                                                         const std::vector<int>& /*num_dofs_per_proc*/,
                                                         int dof_index_idx,
                                                         Pointer<PatchLevel<NDIM> > patch_level,
                                                         Pointer<CoarseFineBoundary<NDIM> > cf_boundary)
{
    // Determine the subdomains associated with this processor.
    const int n_local_patches = patch_level->getProcessorMapping().getNumberOfLocalIndices();
    std::vector<std::vector<Box<NDIM> > > overlap_boxes(n_local_patches), nonoverlap_boxes(n_local_patches);
    int patch_counter = 0, subdomain_counter = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++, ++patch_counter)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        IndexUtilities::partitionPatchBox(
            overlap_boxes[patch_counter], nonoverlap_boxes[patch_counter], patch_box, box_size, overlap_size);
        subdomain_counter += overlap_boxes[patch_counter].size();
    }
    is_overlap.resize(subdomain_counter);
    is_nonoverlap.resize(subdomain_counter);

    // Fill in the IS'es
    const int level_num = patch_level->getLevelNumber();
    subdomain_counter = 0, patch_counter = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++, ++patch_counter)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Box<NDIM> side_patch_box[NDIM];
        for (int axis = 0; axis < NDIM; ++axis)
        {
            side_patch_box[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
        }
        Pointer<SideData<NDIM, int> > dof_data = patch->getPatchData(dof_index_idx);
        const int data_depth = dof_data->getDepth();
#if !defined(NDEBUG)
        TBOX_ASSERT(data_depth == 1);
        TBOX_ASSERT(dof_data->getGhostCellWidth().min() >= overlap_size.max());
#endif

        // Check if the patch touches physical boundary.
        Array<Array<bool> > touches_physical_bdry(NDIM);
        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const bool patch_touches_physical_bdry = pgeom->intersectsPhysicalBoundary();
        if (patch_touches_physical_bdry)
        {
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                touches_physical_bdry[axis].resizeArray(2);
                for (int upperlower = 0; upperlower < 2; ++upperlower)
                {
                    touches_physical_bdry[axis][upperlower] = pgeom->getTouchesRegularBoundary(axis, upperlower);
                }
            }
        }

        // Check if the patch touches the c-f interface on the upper side of the patch.
        Array<Array<bool> > touches_cf_bdry(NDIM);
        Array<std::vector<Box<NDIM> > > upper_side_cf_bdry_box(NDIM);
        const Array<BoundaryBox<NDIM> >& cf_codim1_boxes =
            (level_num == 0) ? Array<BoundaryBox<NDIM> >() :
                               cf_boundary->getBoundaries(patch->getPatchNumber(), /* boundary type */ 1);
        const int n_cf_codim1_boxes = cf_codim1_boxes.size();
        const bool patch_touches_cf_bdry = n_cf_codim1_boxes;
        if (patch_touches_cf_bdry)
        {
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                touches_cf_bdry[axis].resizeArray(2);
                touches_cf_bdry[axis][LOWER] = false;
                touches_cf_bdry[axis][UPPER] = false;
            }
            for (int k = 0; k < n_cf_codim1_boxes; ++k)
            {
                const BoundaryBox<NDIM>& cf_bdry_box = cf_codim1_boxes[k];
                const Box<NDIM>& bdry_box = cf_bdry_box.getBox();
                const unsigned int location_index = cf_bdry_box.getLocationIndex();
                const unsigned int bdry_normal_axis = location_index / 2;
                touches_cf_bdry[bdry_normal_axis][location_index % 2] = true;
                if (location_index % 2 == UPPER)
                {
                    upper_side_cf_bdry_box[bdry_normal_axis].push_back(bdry_box);
                }
            }
        }

        int n_patch_subdomains = static_cast<int>(nonoverlap_boxes[patch_counter].size());
        for (int k = 0; k < n_patch_subdomains; ++k, ++subdomain_counter)
        {
            // The nonoverlapping subdomains.
            const Box<NDIM>& box_local = nonoverlap_boxes[patch_counter][k];
            Box<NDIM> side_box_local[NDIM];
            for (int axis = 0; axis < NDIM; ++axis)
            {
                side_box_local[axis] = SideGeometry<NDIM>::toSideBox(box_local, axis);
            }
            std::set<int> box_local_dofs;

            // Get the local DOFs.
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator b(side_box_local[axis]); b; b++)
                {
                    const CellIndex<NDIM>& i = b();
                    const bool at_upper_subdomain_bdry = (i(axis) == side_box_local[axis].upper(axis));
                    const bool at_upper_patch_bdry = (i(axis) == side_patch_box[axis].upper(axis));
                    const bool at_upper_physical_bdry =
                        at_upper_patch_bdry && patch_touches_physical_bdry && touches_physical_bdry[axis][UPPER];
                    const bool at_upper_cf_bdry = at_upper_patch_bdry && patch_touches_cf_bdry &&
                                                  touches_cf_bdry[axis][UPPER] &&
                                                  is_cf_bdry_idx(i, upper_side_cf_bdry_box[axis]);
                    if (!at_upper_subdomain_bdry || at_upper_physical_bdry || at_upper_cf_bdry)
                    {
                        const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                        box_local_dofs.insert((*dof_data)(i_s));
                    }
                }
            }
            const int n_local = static_cast<int>(box_local_dofs.size());
            PetscInt* box_local_dof_arr;
            PetscMalloc1(n_local, &box_local_dof_arr);
            std::copy(box_local_dofs.begin(), box_local_dofs.end(), box_local_dof_arr);
            ISCreateGeneral(
                PETSC_COMM_SELF, n_local, box_local_dof_arr, PETSC_OWN_POINTER, &is_nonoverlap[subdomain_counter]);

            // The overlapping subdomains.
            const Box<NDIM>& box_overlap = overlap_boxes[patch_counter][k];
            Box<NDIM> side_box_overlap[NDIM];
            for (int axis = 0; axis < NDIM; ++axis)
            {
                side_box_overlap[axis] = SideGeometry<NDIM>::toSideBox(box_overlap, axis);
            }
            std::set<int> box_overlap_dofs;

            // Get the overlap DOFs.
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator b(side_box_overlap[axis]); b; b++)
                {
                    const CellIndex<NDIM>& i = b();
                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                    const int dof_idx = (*dof_data)(i_s);

                    // We keep only those DOFs that are inside the
                    // physical domain and on physical and c-f boundaries.
                    // Some of the DOFs may be on other processors.
                    if (dof_idx >= 0)
                    {
                        box_overlap_dofs.insert(dof_idx);
                    }
                }
            }
            const int n_overlap = static_cast<int>(box_overlap_dofs.size());
            PetscInt* box_overlap_dof_arr;
            PetscMalloc1(n_overlap, &box_overlap_dof_arr);
            std::copy(box_overlap_dofs.begin(), box_overlap_dofs.end(), box_overlap_dof_arr);
            ISCreateGeneral(
                PETSC_COMM_SELF, n_overlap, box_overlap_dof_arr, PETSC_OWN_POINTER, &is_overlap[subdomain_counter]);
        }
    }
    return;
} // constructPatchLevelASMSubdomains_side

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
