// Filename: StaggeredStokesPETScMatUtilities.cpp
// Created on 03 Apr 2012 by Boyce Griffith
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

#include <stddef.h>
#include <algorithm>
#include <numeric>
#include <ostream>
#include <vector>

#include "ArrayData.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellGeometry.h"
#include "CellIndex.h"
#include "CoarseFineBoundary.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "PatchLevel.h"
#include "PoissonSpecifications.h"
#include "RobinBcCoefStrategy.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "Variable.h"
#include "boost/array.hpp"
#include "ibamr/StaggeredStokesPETScMatUtilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "ibtk/compiler_hints.h"
#include "ibtk/ibtk_utilities.h"
#include "petscmat.h"
#include "petscsys.h"
#include "tbox/Array.h"
#include "tbox/MathUtilities.h"
#include "tbox/Pointer.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline Box<NDIM>
compute_tangential_extension(const Box<NDIM>& box, const int data_axis)
{
    Box<NDIM> extended_box = box;
    extended_box.upper()(data_axis) += 1;
    return extended_box;
} // compute_tangential_extension

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
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(
    Mat& mat,
    const PoissonSpecifications& u_problem_coefs,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    double data_time,
    const std::vector<int>& num_dofs_per_proc,
    int u_dof_index_idx,
    int p_dof_index_idx,
    Pointer<PatchLevel<NDIM> > patch_level)
{
    int ierr;
    if (mat)
    {
        ierr = MatDestroy(&mat);
        IBTK_CHKERRQ(ierr);
    }

    // Setup the finite difference stencils.
    static const int uu_stencil_sz = 2 * NDIM + 1;
    boost::array<Index<NDIM>, uu_stencil_sz> uu_stencil(array_constant<Index<NDIM>, uu_stencil_sz>(Index<NDIM>(0)));
    for (unsigned int axis = 0, uu_stencil_index = 1; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side, ++uu_stencil_index)
        {
            uu_stencil[uu_stencil_index](axis) = (side == 0 ? -1 : +1);
        }
    }
    static const int up_stencil_sz = 2;
    boost::array<boost::array<Index<NDIM>, up_stencil_sz>, NDIM> up_stencil(
        array_constant<boost::array<Index<NDIM>, up_stencil_sz>, NDIM>(
            array_constant<Index<NDIM>, up_stencil_sz>(Index<NDIM>(0))));
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side)
        {
            up_stencil[axis][side](axis) = (side == 0 ? -1 : 0);
        }
    }
    static const int pu_stencil_sz = 2 * NDIM;
    boost::array<Index<NDIM>, pu_stencil_sz> pu_stencil(array_constant<Index<NDIM>, pu_stencil_sz>(Index<NDIM>(0)));
    for (unsigned int axis = 0, pu_stencil_index = 0; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side, ++pu_stencil_index)
        {
            pu_stencil[pu_stencil_index](axis) = (side == 0 ? 0 : +1);
        }
    }

    // Determine the index ranges.
    const int mpi_rank = SAMRAI_MPI::getRank();
    const int nlocal = num_dofs_per_proc[mpi_rank];
    const int ilower = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int iupper = ilower + nlocal;
    const int ntotal = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.end(), 0);

    // Determine the non-zero structure of the matrix.
    std::vector<int> d_nnz(nlocal, 0), o_nnz(nlocal, 0);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM, int> > u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int> > p_dof_index_data = patch->getPatchData(p_dof_index_idx);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
            {
                const CellIndex<NDIM>& ic = b();
                const SideIndex<NDIM> is(ic, axis, SideIndex<NDIM>::Lower);
                const int u_dof_index = (*u_dof_index_data)(is);
                if (UNLIKELY(ilower > u_dof_index || u_dof_index >= iupper)) continue;
                const int u_local_idx = u_dof_index - ilower;
                d_nnz[u_local_idx] += 1;
                for (unsigned int d = 0, uu_stencil_index = 1; d < NDIM; ++d)
                {
                    for (int side = 0; side <= 1; ++side, ++uu_stencil_index)
                    {
                        const int uu_dof_index = (*u_dof_index_data)(is + uu_stencil[uu_stencil_index]);
                        if (LIKELY(uu_dof_index >= ilower && uu_dof_index < iupper))
                        {
                            d_nnz[u_local_idx] += 1;
                        }
                        else
                        {
                            o_nnz[u_local_idx] += 1;
                        }
                    }
                }
                for (int side = 0, up_stencil_index = 0; side <= 1; ++side, ++up_stencil_index)
                {
                    const int up_dof_index = (*p_dof_index_data)(ic + up_stencil[axis][up_stencil_index]);
                    if (LIKELY(up_dof_index >= ilower && up_dof_index < iupper))
                    {
                        d_nnz[u_local_idx] += 1;
                    }
                    else
                    {
                        o_nnz[u_local_idx] += 1;
                    }
                }
                d_nnz[u_local_idx] = std::min(nlocal, d_nnz[u_local_idx]);
                o_nnz[u_local_idx] = std::min(ntotal - nlocal, o_nnz[u_local_idx]);
            }
        }
        for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
        {
            const CellIndex<NDIM>& ic = b();
            const int p_dof_index = (*p_dof_index_data)(ic);
            if (UNLIKELY(ilower > p_dof_index || p_dof_index >= iupper)) continue;
            const int p_local_idx = p_dof_index - ilower;
            d_nnz[p_local_idx] += 1;
            for (unsigned int axis = 0, pu_stencil_index = 0; axis < NDIM; ++axis)
            {
                for (int side = 0; side <= 1; ++side, ++pu_stencil_index)
                {
                    const int pu_dof_index = (*u_dof_index_data)(
                        SideIndex<NDIM>(ic + pu_stencil[pu_stencil_index], axis, SideIndex<NDIM>::Lower));
                    if (LIKELY(pu_dof_index >= ilower && pu_dof_index < iupper))
                    {
                        d_nnz[p_local_idx] += 1;
                    }
                    else
                    {
                        o_nnz[p_local_idx] += 1;
                    }
                }
            }
            d_nnz[p_local_idx] = std::min(nlocal, d_nnz[p_local_idx]);
            o_nnz[p_local_idx] = std::min(ntotal - nlocal, o_nnz[p_local_idx]);
        }
    }

    // Create an empty matrix.
    ierr = MatCreateAIJ(PETSC_COMM_WORLD,
                        nlocal,
                        nlocal,
                        PETSC_DETERMINE,
                        PETSC_DETERMINE,
                        0,
                        nlocal ? &d_nnz[0] : NULL,
                        0,
                        nlocal ? &o_nnz[0] : NULL,
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
    const double C = u_problem_coefs.getCConstant();
    const double D = u_problem_coefs.getDConstant();
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();

        const IntVector<NDIM> no_ghosts(0);
        SideData<NDIM, double> uu_matrix_coefs(patch_box, uu_stencil_sz, no_ghosts);
        SideData<NDIM, double> up_matrix_coefs(patch_box, up_stencil_sz, no_ghosts);
        CellData<NDIM, double> pu_matrix_coefs(patch_box, pu_stencil_sz, no_ghosts);

        // Compute all matrix coefficients, including those on the physical
        // boundary; however, do not yet take physical boundary conditions into
        // account.  Boundary conditions are handled subsequently.
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            std::vector<double> uu_mat_vals(uu_stencil_sz, 0.0);
            uu_mat_vals[0] = C; // diagonal
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const double dx_sq = dx[d] * dx[d];
                uu_mat_vals[0] -= 2 * D / dx_sq;    // diagonal
                uu_mat_vals[2 * d + 1] = D / dx_sq; // lower off-diagonal
                uu_mat_vals[2 * d + 2] = D / dx_sq; // upper off-diagonal
            }
            for (int uu_stencil_index = 0; uu_stencil_index < uu_stencil_sz; ++uu_stencil_index)
            {
                uu_matrix_coefs.fill(uu_mat_vals[uu_stencil_index], uu_stencil_index);
            }

            // grad p
            for (int d = 0; d < NDIM; ++d)
            {
                up_matrix_coefs.getArrayData(d).fill(-1.0 / dx[d], 0);
                up_matrix_coefs.getArrayData(d).fill(+1.0 / dx[d], 1);
            }

            // -div u
            std::vector<double> pu_mat_vals(pu_stencil_sz, 0.0);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                pu_matrix_coefs.fill(+1.0 / dx[d], 2 * d);
                pu_matrix_coefs.fill(-1.0 / dx[d], 2 * d + 1);
            }
        }

        // Data structures required to set physical boundary conditions.
        const Array<BoundaryBox<NDIM> > physical_codim1_boxes =
            PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
        const int n_physical_codim1_boxes = physical_codim1_boxes.size();
        const double* const patch_x_lower = pgeom->getXLower();
        const double* const patch_x_upper = pgeom->getXUpper();
        const IntVector<NDIM>& ratio_to_level_zero = pgeom->getRatio();
        Array<Array<bool> > touches_regular_bdry(NDIM), touches_periodic_bdry(NDIM);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            touches_regular_bdry[axis].resizeArray(2);
            touches_periodic_bdry[axis].resizeArray(2);
            for (int upperlower = 0; upperlower < 2; ++upperlower)
            {
                touches_regular_bdry[axis][upperlower] = pgeom->getTouchesRegularBoundary(axis, upperlower);
                touches_periodic_bdry[axis][upperlower] = pgeom->getTouchesPeriodicBoundary(axis, upperlower);
            }
        }

        // Modify matrix coefficients to account for physical boundary
        // conditions along boundaries which ARE NOT aligned with the data axis.
        //
        // NOTE: It important to set these values first to avoid problems at
        // corners in the physical domain.  In particular, since Dirichlet
        // boundary conditions for values located on the physical boundary
        // override all other boundary conditions, we set those values last.
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                const unsigned int location_index = bdry_box.getLocationIndex();
                const unsigned int bdry_normal_axis = location_index / 2;
                const bool is_lower = location_index % 2 == 0;

                if (bdry_normal_axis == axis) continue;

                const Box<NDIM> bc_fill_box =
                    pgeom->getBoundaryFillBox(bdry_box, patch_box, /* ghost_width_to_fill */ IntVector<NDIM>(1));
                const BoundaryBox<NDIM> trimmed_bdry_box =
                    PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
                const Box<NDIM> bc_coef_box = compute_tangential_extension(
                    PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box), axis);

                Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double> > gcoef_data;

                // Temporarily reset the patch geometry object associated with
                // the patch so that boundary conditions are set at the correct
                // spatial locations.
                boost::array<double, NDIM> shifted_patch_x_lower, shifted_patch_x_upper;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    shifted_patch_x_lower[d] = patch_x_lower[d];
                    shifted_patch_x_upper[d] = patch_x_upper[d];
                }
                shifted_patch_x_lower[axis] -= 0.5 * dx[axis];
                shifted_patch_x_upper[axis] -= 0.5 * dx[axis];
                patch->setPatchGeometry(new CartesianPatchGeometry<NDIM>(ratio_to_level_zero,
                                                                         touches_regular_bdry,
                                                                         touches_periodic_bdry,
                                                                         dx,
                                                                         shifted_patch_x_lower.data(),
                                                                         shifted_patch_x_upper.data()));

                // Set the boundary condition coefficients.
                static const bool homogeneous_bc = true;
                ExtendedRobinBcCoefStrategy* extended_bc_coef =
                    dynamic_cast<ExtendedRobinBcCoefStrategy*>(u_bc_coefs[axis]);
                if (extended_bc_coef)
                {
                    extended_bc_coef->clearTargetPatchDataIndex();
                    extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                }
                u_bc_coefs[axis]->setBcCoefs(
                    acoef_data, bcoef_data, gcoef_data, NULL, *patch, trimmed_bdry_box, data_time);
                if (gcoef_data && homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

                // Restore the original patch geometry object.
                patch->setPatchGeometry(pgeom);

                // Modify the matrix coefficients to account for homogeneous
                // boundary conditions.
                for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
                {
                    const Index<NDIM>& i = bc();
                    const double& a = (*acoef_data)(i, 0);
                    const double& b = (*bcoef_data)(i, 0);
                    const bool velocity_bc = (a == 1.0 || MathUtilities<double>::equalEps(a, 1.0));
                    const bool traction_bc = (b == 1.0 || MathUtilities<double>::equalEps(b, 1.0));
#if !defined(NDEBUG)
                    TBOX_ASSERT((velocity_bc || traction_bc) && !(velocity_bc && traction_bc));
#endif
                    Index<NDIM> i_intr = i;
                    if (is_lower)
                    {
                        i_intr(bdry_normal_axis) += 0;
                    }
                    else
                    {
                        i_intr(bdry_normal_axis) -= 1;
                    }
                    const SideIndex<NDIM> i_s(i_intr, axis, SideIndex<NDIM>::Lower);

                    if (velocity_bc)
                    {
                        if (is_lower)
                        {
                            uu_matrix_coefs(i_s, 0) -= uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 1);
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 1) = 0.0;
                        }
                        else
                        {
                            uu_matrix_coefs(i_s, 0) -= uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 2);
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 2) = 0.0;
                        }
                    }
                    else if (traction_bc)
                    {
                        if (is_lower)
                        {
                            uu_matrix_coefs(i_s, 0) += uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 1);
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 1) = 0.0;
                        }
                        else
                        {
                            uu_matrix_coefs(i_s, 0) -= uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 2);
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 2) = 0.0;
                        }
                    }
                    else
                    {
                        TBOX_ERROR(
                            "StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(): Unknown BC type for "
                            "tangential velocity specified.");
                    }
                }
            }
        }

        // Modify matrix coefficients to account for physical boundary
        // conditions along boundaries which ARE aligned with the data axis.
        //
        // NOTE: It important to set these values last to avoid problems at corners
        // in the physical domain.  In particular, since Dirichlet boundary
        // conditions for values located on the physical boundary override all other
        // boundary conditions, we set those values last.
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                const unsigned int location_index = bdry_box.getLocationIndex();
                const unsigned int bdry_normal_axis = location_index / 2;
                const bool is_lower = location_index % 2 == 0;

                if (bdry_normal_axis != axis) continue;

                const Box<NDIM> bc_fill_box =
                    pgeom->getBoundaryFillBox(bdry_box, patch_box, /* ghost_width_to_fill */ IntVector<NDIM>(1));
                const BoundaryBox<NDIM> trimmed_bdry_box =
                    PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
                const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

                Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double> > gcoef_data;

                // Set the boundary condition coefficients.
                static const bool homogeneous_bc = true;
                ExtendedRobinBcCoefStrategy* extended_bc_coef =
                    dynamic_cast<ExtendedRobinBcCoefStrategy*>(u_bc_coefs[axis]);
                if (extended_bc_coef)
                {
                    extended_bc_coef->clearTargetPatchDataIndex();
                    extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                }
                u_bc_coefs[axis]->setBcCoefs(
                    acoef_data, bcoef_data, gcoef_data, NULL, *patch, trimmed_bdry_box, data_time);
                if (gcoef_data && homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

                // Modify the matrix coefficients to account for homogeneous
                // boundary conditions.
                for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
                {
                    const Index<NDIM>& i = bc();
                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                    const double& a = (*acoef_data)(i, 0);
                    const double& b = (*bcoef_data)(i, 0);
                    const bool velocity_bc = (a == 1.0 || MathUtilities<double>::equalEps(a, 1.0));
                    const bool traction_bc = (b == 1.0 || MathUtilities<double>::equalEps(b, 1.0));
#if !defined(NDEBUG)
                    TBOX_ASSERT((velocity_bc || traction_bc) && !(velocity_bc && traction_bc));
#endif
                    if (velocity_bc)
                    {
                        uu_matrix_coefs(i_s, 0) = 1.0;
                        for (int k = 1; k < uu_stencil_sz; ++k)
                        {
                            uu_matrix_coefs(i_s, k) = 0.0;
                        }
                        for (int k = 0; k < up_stencil_sz; ++k)
                        {
                            up_matrix_coefs(i_s, k) = 0.0;
                        }
                    }
                    else if (traction_bc)
                    {
                        if (is_lower)
                        {
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 2) +=
                                uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 1);
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 1) = 0.0;
                        }
                        else
                        {
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 1) +=
                                uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 2);
                            uu_matrix_coefs(i_s, 2 * bdry_normal_axis + 2) = 0.0;
                        }
                    }
                    else
                    {
                        TBOX_ERROR(
                            "StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(): Unknown BC type for "
                            "normal velocity specified.");
                    }
                }
            }
        }

        // Set matrix coefficients.
        Pointer<SideData<NDIM, int> > u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int> > p_dof_index_data = patch->getPatchData(p_dof_index_idx);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
            {
                const CellIndex<NDIM>& ic = b();
                const SideIndex<NDIM> is(ic, axis, SideIndex<NDIM>::Lower);
                const int u_dof_index = (*u_dof_index_data)(is);
                if (UNLIKELY(ilower > u_dof_index || u_dof_index >= iupper)) continue;

                const int u_stencil_sz = uu_stencil_sz + up_stencil_sz;
                std::vector<double> u_mat_vals(u_stencil_sz);
                std::vector<int> u_mat_cols(u_stencil_sz);

                u_mat_vals[0] = uu_matrix_coefs(is, 0);
                u_mat_cols[0] = u_dof_index;
                for (unsigned int d = 0, uu_stencil_index = 1; d < NDIM; ++d)
                {
                    for (int side = 0; side <= 1; ++side, ++uu_stencil_index)
                    {
                        u_mat_vals[uu_stencil_index] = uu_matrix_coefs(is, uu_stencil_index);
                        u_mat_cols[uu_stencil_index] = (*u_dof_index_data)(is + uu_stencil[uu_stencil_index]);
                    }
                }
                for (int side = 0, up_stencil_index = 0; side <= 1; ++side, ++up_stencil_index)
                {
                    u_mat_vals[uu_stencil_sz + side] = up_matrix_coefs(is, up_stencil_index);
                    u_mat_cols[uu_stencil_sz + side] = (*p_dof_index_data)(ic + up_stencil[axis][up_stencil_index]);
                }

                ierr = MatSetValues(mat, 1, &u_dof_index, u_stencil_sz, &u_mat_cols[0], &u_mat_vals[0], INSERT_VALUES);
                IBTK_CHKERRQ(ierr);
            }
        }

        for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
        {
            const CellIndex<NDIM>& ic = b();
            const int p_dof_index = (*p_dof_index_data)(ic);
            if (UNLIKELY(ilower > p_dof_index || p_dof_index >= iupper)) continue;

            const int p_stencil_sz = pu_stencil_sz + 1;
            std::vector<double> p_mat_vals(p_stencil_sz);
            std::vector<int> p_mat_cols(p_stencil_sz);

            for (unsigned int axis = 0, pu_stencil_index = 0; axis < NDIM; ++axis)
            {
                for (int side = 0; side <= 1; ++side, ++pu_stencil_index)
                {
                    p_mat_vals[pu_stencil_index] = pu_matrix_coefs(ic, pu_stencil_index);
                    p_mat_cols[pu_stencil_index] = (*u_dof_index_data)(
                        SideIndex<NDIM>(ic + pu_stencil[pu_stencil_index], axis, SideIndex<NDIM>::Lower));
                }
            }
            p_mat_vals[pu_stencil_sz] = 0.0;
            p_mat_cols[pu_stencil_sz] = p_dof_index;

            ierr = MatSetValues(mat, 1, &p_dof_index, p_stencil_sz, &p_mat_cols[0], &p_mat_vals[0], INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    return;
} // constructPatchLevelMACStokesOp

void
StaggeredStokesPETScMatUtilities::constructPatchLevelASMSubdomains(std::vector<IS>& is_overlap,
                                                                   std::vector<IS>& is_nonoverlap,
                                                                   const IntVector<NDIM>& box_size,
                                                                   const IntVector<NDIM>& overlap_size,
                                                                   const std::vector<int>& /*num_dofs_per_proc*/,
                                                                   int u_dof_index_idx,
                                                                   int p_dof_index_idx,
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
        Pointer<SideData<NDIM, int> > u_dof_data = patch->getPatchData(u_dof_index_idx);
        const int u_data_depth = u_dof_data->getDepth();
        Pointer<CellData<NDIM, int> > p_dof_data = patch->getPatchData(p_dof_index_idx);
        const int p_data_depth = p_dof_data->getDepth();
#if !defined(NDEBUG)
        TBOX_ASSERT(u_data_depth == 1);
        TBOX_ASSERT(p_data_depth == 1);
        TBOX_ASSERT(u_dof_data->getGhostCellWidth().min() >= overlap_size.max());
        TBOX_ASSERT(p_dof_data->getGhostCellWidth().min() >= overlap_size.max());
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
            int box_local_dofs_size = 0;
            for (int axis = 0; axis < NDIM; ++axis)
            {
                side_box_local[axis] = SideGeometry<NDIM>::toSideBox(box_local, axis);
                box_local_dofs_size += side_box_local[axis].size();
            }
            box_local_dofs_size += box_local.size();
            std::vector<int> box_local_dofs;
            box_local_dofs.reserve(box_local_dofs_size);

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
                        box_local_dofs.push_back((*u_dof_data)(i_s));
                    }
                }
            }
            for (Box<NDIM>::Iterator b(box_local); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                box_local_dofs.push_back((*p_dof_data)(i));
            }
            std::sort(box_local_dofs.begin(), box_local_dofs.end());
            const int n_local = static_cast<int>(box_local_dofs.size());
            ISCreateGeneral(
                PETSC_COMM_SELF, n_local, &box_local_dofs[0], PETSC_COPY_VALUES, &is_nonoverlap[subdomain_counter]);

            // The overlapping subdomains.
            const Box<NDIM>& box_overlap = overlap_boxes[patch_counter][k];
            Box<NDIM> side_box_overlap[NDIM];
            int box_overlap_dofs_size = 0;
            for (int axis = 0; axis < NDIM; ++axis)
            {
                side_box_overlap[axis] = SideGeometry<NDIM>::toSideBox(box_overlap, axis);
                box_overlap_dofs_size += side_box_overlap[axis].size();
            }
            box_overlap_dofs_size += box_overlap.size();
            std::vector<int> box_overlap_dofs;
            box_overlap_dofs.reserve(box_overlap_dofs_size);

            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator b(side_box_overlap[axis]); b; b++)
                {
                    const CellIndex<NDIM>& i = b();
                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                    const int dof_idx = (*u_dof_data)(i_s);

                    // We keep only those DOFs that are inside the
                    // physical domain and on physical and c-f boundaries.
                    // Some of the DOFs may be on other processors.
                    if (dof_idx >= 0)
                    {
                        box_overlap_dofs.push_back(dof_idx);
                    }
                }
            }
            for (Box<NDIM>::Iterator b(box_overlap); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                const int dof_idx = (*p_dof_data)(i);

                // We keep only those DOFs that are inside the
                // physical domain and on physical and c-f boundaries.
                // Some of the DOFs may be on other processors.
                if (dof_idx >= 0)
                {
                    box_overlap_dofs.push_back(dof_idx);
                }
            }

            std::sort(box_overlap_dofs.begin(), box_overlap_dofs.end());
            box_overlap_dofs.erase(std::unique(box_overlap_dofs.begin(), box_overlap_dofs.end()),
                                   box_overlap_dofs.end());
            const int n_overlap = static_cast<int>(box_overlap_dofs.size());
            ISCreateGeneral(
                PETSC_COMM_SELF, n_overlap, &box_overlap_dofs[0], PETSC_COPY_VALUES, &is_overlap[subdomain_counter]);
        }
    }
    return;
} // constructPatchLevelASMSubdomains

void
StaggeredStokesPETScMatUtilities::constructPatchLevelMSMSubdomains(std::vector<IS>& is_row,
                                                                   std::vector<IS>& is_col,
                                                                   const IntVector<NDIM>& box_size,
                                                                   const IntVector<NDIM>& overlap_size,
                                                                   const std::vector<int>& num_dofs_per_proc,
                                                                   int u_dof_index_idx,
                                                                   int p_dof_index_idx,
                                                                   Pointer<PatchLevel<NDIM> > patch_level)
{
    // Destroy the previously stored IS'es
    int ierr;
    for (unsigned int k = 0; k < is_row.size(); ++k)
    {
        ierr = ISDestroy(&is_row[k]);
        IBTK_CHKERRQ(ierr);
    }
    is_row.clear();
    for (unsigned int k = 0; k < is_col.size(); ++k)
    {
        ierr = ISDestroy(&is_col[k]);
        IBTK_CHKERRQ(ierr);
    }
    is_col.clear();

    // Determine the subdomains associated with this processor.
    const int n_local_patches = patch_level->getProcessorMapping().getNumberOfLocalIndices();
    std::vector<std::vector<Box<NDIM> > > overlap_boxes(n_local_patches), nonoverlap_boxes(n_local_patches);
    int patch_counter = 0, n_local_subdomains = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++, ++patch_counter)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        IndexUtilities::partitionPatchBox(
            overlap_boxes[patch_counter], nonoverlap_boxes[patch_counter], patch_box, box_size, overlap_size);

        // Keep only local part of overlap boxes.
        const int patch_subdomains = static_cast<int>(overlap_boxes[patch_counter].size());
        for (int k = 0; k < patch_subdomains; ++k)
        {
            Box<NDIM> box_intersection = overlap_boxes[patch_counter][k] * patch_box;
            overlap_boxes[patch_counter][k] = box_intersection;
        }
        n_local_subdomains += patch_subdomains;
    }
    is_row.resize(n_local_subdomains);
    is_col.resize(n_local_subdomains);

    // DOFs on this processor.
    const int mpi_rank = SAMRAI_MPI::getRank();
    const int n_local_dofs = num_dofs_per_proc[mpi_rank];
    const int first_local_dof = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int last_local_dof = first_local_dof + n_local_dofs;
    std::vector<int> all_local_dofs(n_local_dofs);
    for (int k = 0; k < n_local_dofs; ++k)
    {
        all_local_dofs[k] = first_local_dof + k;
    }

    // Fill in the IS'es
    int subdomain_counter = 0;
    patch_counter = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++, ++patch_counter)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<SideData<NDIM, int> > u_dof_data = patch->getPatchData(u_dof_index_idx);
        const int u_data_depth = u_dof_data->getDepth();
        Pointer<CellData<NDIM, int> > p_dof_data = patch->getPatchData(p_dof_index_idx);
        const int p_data_depth = p_dof_data->getDepth();
#if !defined(NDEBUG)
        TBOX_ASSERT(u_data_depth == 1);
        TBOX_ASSERT(p_data_depth == 1);
#endif

        int n_patch_subdomains = static_cast<int>(overlap_boxes[patch_counter].size());
        for (int s = 0; s < n_patch_subdomains; ++s, ++subdomain_counter)
        {
            // Get the overlapping subdomain.
            const Box<NDIM>& box_overlap = overlap_boxes[patch_counter][s];
            Box<NDIM> side_box_overlap[NDIM];
            int box_overlap_dofs_size = 0;
            for (int axis = 0; axis < NDIM; ++axis)
            {
                side_box_overlap[axis] = SideGeometry<NDIM>::toSideBox(box_overlap, axis);
                box_overlap_dofs_size += side_box_overlap[axis].size();
            }
            box_overlap_dofs_size += box_overlap.size();
            std::vector<int> box_row_dofs;
            box_row_dofs.reserve(box_overlap_dofs_size);

            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator b(side_box_overlap[axis]); b; b++)
                {
                    const CellIndex<NDIM>& i = b();
                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                    const int dof_idx = (*u_dof_data)(i_s);
                    if (dof_idx >= first_local_dof && dof_idx < last_local_dof)
                    {
                        box_row_dofs.push_back(dof_idx);
                    }
                }
            }
            for (Box<NDIM>::Iterator b(box_overlap); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                const int dof_idx = (*p_dof_data)(i);
                if (dof_idx >= first_local_dof && dof_idx < last_local_dof)
                {
                    box_row_dofs.push_back(dof_idx);
                }
            }

            // Create the list of rows
            std::sort(box_row_dofs.begin(), box_row_dofs.end());
            const int n_rows = static_cast<int>(box_row_dofs.size());
            ISCreateGeneral(PETSC_COMM_SELF, n_rows, &box_row_dofs[0], PETSC_COPY_VALUES, &is_row[subdomain_counter]);

            // Create the list of cols.
            const int n_cols = n_local_dofs - n_rows;
            std::vector<int> box_col_dofs(n_cols);
            std::set_difference(all_local_dofs.begin(),
                                all_local_dofs.end(),
                                box_row_dofs.begin(),
                                box_row_dofs.end(),
                                box_col_dofs.begin());
            ISCreateGeneral(PETSC_COMM_SELF, n_cols, &box_col_dofs[0], PETSC_COPY_VALUES, &is_col[subdomain_counter]);
        }
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(subdomain_counter == n_local_subdomains);
#endif

    return;
} // constructPatchLevelMSMSubdomains

void
StaggeredStokesPETScMatUtilities::constructPatchLevelMSMSubdomains(std::vector<IS>& is_red_row,
                                                                   std::vector<IS>& is_red_col,
                                                                   std::vector<IS>& is_black_row,
                                                                   std::vector<IS>& is_black_col,
                                                                   const IntVector<NDIM>& box_size,
                                                                   const IntVector<NDIM>& overlap_size,
                                                                   const std::vector<int>& num_dofs_per_proc,
                                                                   int u_dof_index_idx,
                                                                   int p_dof_index_idx,
                                                                   Pointer<PatchLevel<NDIM> > patch_level)
{
    // Destroy the previously stored IS'es
    int ierr;
    for (unsigned int k = 0; k < is_red_row.size(); ++k)
    {
        ierr = ISDestroy(&is_red_row[k]);
        IBTK_CHKERRQ(ierr);
    }
    is_red_row.clear();
    for (unsigned int k = 0; k < is_black_row.size(); ++k)
    {
        ierr = ISDestroy(&is_black_row[k]);
        IBTK_CHKERRQ(ierr);
    }
    is_black_row.clear();
    for (unsigned int k = 0; k < is_red_col.size(); ++k)
    {
        ierr = ISDestroy(&is_red_col[k]);
        IBTK_CHKERRQ(ierr);
    }
    is_red_col.clear();
    for (unsigned int k = 0; k < is_black_col.size(); ++k)
    {
        ierr = ISDestroy(&is_black_col[k]);
        IBTK_CHKERRQ(ierr);
    }
    is_black_col.clear();

    // Determine the subdomains associated with this processor.
    const int n_local_patches = patch_level->getProcessorMapping().getNumberOfLocalIndices();
    std::vector<std::vector<Box<NDIM> > > overlap_boxes(n_local_patches), nonoverlap_boxes(n_local_patches);
    std::vector<SAMRAI::hier::IntVector<NDIM> > n_patch_partitions(n_local_patches);
    int patch_counter = 0, n_local_subdomains = 0, n_red_subdomains = 0, n_black_subdomains = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++, ++patch_counter)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        n_patch_partitions[patch_counter] = IndexUtilities::partitionPatchBox(
            overlap_boxes[patch_counter], nonoverlap_boxes[patch_counter], patch_box, box_size, overlap_size);

        const int Nx = n_patch_partitions[patch_counter][0];
        const int Ny = n_patch_partitions[patch_counter][1];
        const int Nz = NDIM == 2 ? 1 : n_patch_partitions[patch_counter][2];

        int n_patch_red = 0, n_patch_black = 0;
        for (int K = 0; K < Nz; ++K)
        {
            for (int J = 0; J < Ny; ++J)
            {
                for (int I = 0; I < Nx; ++I)
                {
                    if ((K % 2 == 0 && (I + J) % 2 == 0) || (K % 2 == 1 && (I + J) % 2 == 1))
                    {
                        n_patch_red += 1;
                    }

                    if ((K % 2 == 0 && (I + J) % 2 == 1) || (K % 2 == 1 && (I + J) % 2 == 0))
                    {
                        n_patch_black += 1;
                    }
                }
            }
        }

        n_red_subdomains += n_patch_red;
        n_black_subdomains += n_patch_black;

        // Keep only local part of overlap boxes.
        const int patch_subdomains = static_cast<int>(overlap_boxes[patch_counter].size());
        for (int k = 0; k < patch_subdomains; ++k)
        {
            Box<NDIM> box_intersection = overlap_boxes[patch_counter][k] * patch_box;
            overlap_boxes[patch_counter][k] = box_intersection;
        }
        n_local_subdomains += patch_subdomains;
    }
    is_red_row.resize(n_red_subdomains);
    is_red_col.resize(n_red_subdomains);
    is_black_row.resize(n_black_subdomains);
    is_black_col.resize(n_black_subdomains);

    // DOFs on this processor.
    const int mpi_rank = SAMRAI_MPI::getRank();
    const int n_local_dofs = num_dofs_per_proc[mpi_rank];
    const int first_local_dof = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int last_local_dof = first_local_dof + n_local_dofs;
    std::vector<int> all_local_dofs(n_local_dofs);
    for (int k = 0; k < n_local_dofs; ++k)
    {
        all_local_dofs[k] = first_local_dof + k;
    }

    // Fill in the IS'es
    int subdomain_counter = 0, red_subdomain_counter = 0, black_subdomain_counter = 0;
    patch_counter = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++, ++patch_counter)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        Pointer<SideData<NDIM, int> > u_dof_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int> > p_dof_data = patch->getPatchData(p_dof_index_idx);
#if !defined(NDEBUG)
        const int u_data_depth = u_dof_data->getDepth();
        const int p_data_depth = p_dof_data->getDepth();
        TBOX_ASSERT(u_data_depth == 1);
        TBOX_ASSERT(p_data_depth == 1);
#endif

        const int Nx = n_patch_partitions[patch_counter][0];
        const int Ny = n_patch_partitions[patch_counter][1];
        const int Nz = NDIM == 2 ? 1 : n_patch_partitions[patch_counter][2];

        for (int K = 0; K < Nz; ++K)
        {
            for (int J = 0; J < Ny; ++J)
            {
                for (int I = 0; I < Nx; ++I, ++subdomain_counter)
                {
                    // Get the overlapping subdomain.
                    int subdomain_id = K * Nx * Ny + J * Nx + I;
                    const Box<NDIM>& box_overlap = overlap_boxes[patch_counter][subdomain_id];
                    Box<NDIM> side_box_overlap[NDIM];
                    int box_overlap_dofs_size = 0;
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        side_box_overlap[axis] = SideGeometry<NDIM>::toSideBox(box_overlap, axis);
                        box_overlap_dofs_size += side_box_overlap[axis].size();
                    }
                    box_overlap_dofs_size += box_overlap.size();
                    std::vector<int> box_row_dofs;
                    box_row_dofs.reserve(box_overlap_dofs_size);

                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for (Box<NDIM>::Iterator b(side_box_overlap[axis]); b; b++)
                        {
                            const CellIndex<NDIM>& i = b();
                            const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                            const int dof_idx = (*u_dof_data)(i_s);
                            if (dof_idx >= first_local_dof && dof_idx < last_local_dof)
                            {
                                box_row_dofs.push_back(dof_idx);
                            }
                        }
                    }
                    for (Box<NDIM>::Iterator b(box_overlap); b; b++)
                    {
                        const CellIndex<NDIM>& i = b();
                        const int dof_idx = (*p_dof_data)(i);
                        if (dof_idx >= first_local_dof && dof_idx < last_local_dof)
                        {
                            box_row_dofs.push_back(dof_idx);
                        }
                    }

                    // Create the list of rows and cols.
                    std::sort(box_row_dofs.begin(), box_row_dofs.end());
                    const int n_rows = static_cast<int>(box_row_dofs.size());
                    const int n_cols = n_local_dofs - n_rows;
                    std::vector<int> box_col_dofs(n_cols);
                    std::set_difference(all_local_dofs.begin(),
                                        all_local_dofs.end(),
                                        box_row_dofs.begin(),
                                        box_row_dofs.end(),
                                        box_col_dofs.begin());

                    if ((K % 2 == 0 && (I + J) % 2 == 0) || (K % 2 == 1 && (I + J) % 2 == 1))
                    {
                        ISCreateGeneral(PETSC_COMM_SELF,
                                        n_rows,
                                        &box_row_dofs[0],
                                        PETSC_COPY_VALUES,
                                        &is_red_row[red_subdomain_counter]);
                        ISCreateGeneral(PETSC_COMM_SELF,
                                        n_cols,
                                        &box_col_dofs[0],
                                        PETSC_COPY_VALUES,
                                        &is_red_col[red_subdomain_counter]);

                        red_subdomain_counter += 1;
                    }

                    if ((K % 2 == 0 && (I + J) % 2 == 1) || (K % 2 == 1 && (I + J) % 2 == 0))
                    {
                        ISCreateGeneral(PETSC_COMM_SELF,
                                        n_rows,
                                        &box_row_dofs[0],
                                        PETSC_COPY_VALUES,
                                        &is_black_row[black_subdomain_counter]);
                        ISCreateGeneral(PETSC_COMM_SELF,
                                        n_cols,
                                        &box_col_dofs[0],
                                        PETSC_COPY_VALUES,
                                        &is_black_col[black_subdomain_counter]);

                        black_subdomain_counter += 1;
                    }
                }
            }
        }
    }

// Account for all the domains.
#if !defined(NDEBUG)
    TBOX_ASSERT(subdomain_counter == n_local_subdomains);
    TBOX_ASSERT(black_subdomain_counter == n_black_subdomains);
    TBOX_ASSERT(red_subdomain_counter == n_red_subdomains);
    TBOX_ASSERT((n_red_subdomains + n_black_subdomains) == n_local_subdomains);
#endif

    return;
} // constructPatchLevelMSMSubdomains

void
StaggeredStokesPETScMatUtilities::constructPatchLevelFields(
    std::vector<IS>& is_field,
    std::vector<std::string>& is_field_name,
    const std::vector<int>& num_dofs_per_proc,
    int u_dof_index_idx,
    int p_dof_index_idx,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level)
{
    // Destroy the previously stored IS'es
    int ierr;
    for (unsigned int k = 0; k < is_field.size(); ++k)
    {
        ierr = ISDestroy(&is_field[k]);
        IBTK_CHKERRQ(ierr);
    }
    is_field.resize(2);
    is_field_name.resize(2);

    // Name of the fields.
    is_field_name[0] = "velocity";
    is_field_name[1] = "pressure";

    // DOFs on this processor.
    const int mpi_rank = SAMRAI_MPI::getRank();
    const int n_local_dofs = num_dofs_per_proc[mpi_rank];
    std::vector<int> local_vel_dofs, local_pressure_dofs;
    local_vel_dofs.reserve(n_local_dofs);
    local_pressure_dofs.reserve(n_local_dofs);

    const int first_local_dof = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int last_local_dof = first_local_dof + n_local_dofs;

    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Box<NDIM> side_patch_box[NDIM];
        for (int axis = 0; axis < NDIM; ++axis)
        {
            side_patch_box[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
        }

        Pointer<SideData<NDIM, int> > u_dof_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int> > p_dof_data = patch->getPatchData(p_dof_index_idx);
#if !defined(NDEBUG)
        const int u_data_depth = u_dof_data->getDepth();
        const int p_data_depth = p_dof_data->getDepth();
        TBOX_ASSERT(u_data_depth == 1);
        TBOX_ASSERT(p_data_depth == 1);
#endif

        // Get the local velocity DOFs.
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(side_patch_box[axis]); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                const int dof_idx = (*u_dof_data)(i_s);
                if (dof_idx >= first_local_dof && dof_idx < last_local_dof)
                {
                    local_vel_dofs.push_back(dof_idx);
                }
            }
        }

        // Get the local pressure DOFs.
        for (Box<NDIM>::Iterator b(patch_box); b; b++)
        {
            const CellIndex<NDIM>& i = b();
            const int dof_idx = (*p_dof_data)(i);
            if (dof_idx >= first_local_dof && dof_idx < last_local_dof)
            {
                local_pressure_dofs.push_back(dof_idx);
            }
        }
    }

    std::sort(local_vel_dofs.begin(), local_vel_dofs.end());
    const int n_local_vel_dofs = static_cast<int>(local_vel_dofs.size());
    ISCreateGeneral(PETSC_COMM_WORLD, n_local_vel_dofs, &local_vel_dofs[0], PETSC_COPY_VALUES, &is_field[0]);

    std::sort(local_pressure_dofs.begin(), local_pressure_dofs.end());
    const int n_local_pressure_dofs = static_cast<int>(local_pressure_dofs.size());
    ISCreateGeneral(PETSC_COMM_WORLD, n_local_pressure_dofs, &local_pressure_dofs[0], PETSC_COPY_VALUES, &is_field[1]);

    return;
} // constructPatchLevelFields

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
