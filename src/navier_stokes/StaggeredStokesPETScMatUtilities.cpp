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

#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/solv/PoissonSpecifications.h"
#include "SAMRAI/solv/RobinBcCoefStrategy.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideIndex.h"
#include "SAMRAI/hier/Variable.h"
#include "boost/array.hpp"
#include "ibamr/StaggeredStokesPETScMatUtilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "ibtk/compiler_hints.h"
#include "ibtk/ibtk_utilities.h"
#include "petscmat.h"
#include "petscsys.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline Box compute_tangential_extension(const Box& box, const int data_axis)
{
    Box extended_box = box;
    extended_box.upper(data_axis) += 1;
    return extended_box;
}
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(Mat& mat,
                                                                 const PoissonSpecifications& u_problem_coefs,
                                                                 const std::vector<boost::shared_ptr<RobinBcCoefStrategy>>& u_bc_coefs,
                                                                 double data_time,
                                                                 const std::vector<int>& num_dofs_per_proc,
                                                                 int u_dof_index_idx,
                                                                 int p_dof_index_idx,
                                                                 boost::shared_ptr<PatchLevel> patch_level)
{
    int ierr;
    if (mat)
    {
        ierr = MatDestroy(&mat);
        IBTK_CHKERRQ(ierr);
    }

    // Setup the finite difference stencils.
    static const int uu_stencil_sz = 2 * NDIM + 1;
    std::vector<Index> uu_stencil(uu_stencil_sz, Index::getZeroIndex(DIM));
    for (unsigned int axis = 0, uu_stencil_index = 1; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side, ++uu_stencil_index)
        {
            uu_stencil[uu_stencil_index](axis) = (side == 0 ? -1 : +1);
        }
    }
    static const int up_stencil_sz = 2;
    std::vector<std::vector<Index> > up_stencil(NDIM, std::vector<Index>(up_stencil_sz, Index::getZeroIndex(DIM)));
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side)
        {
            up_stencil[axis][side](axis) = (side == 0 ? -1 : 0);
        }
    }
    static const int pu_stencil_sz = 4;
    std::vector<Index> pu_stencil(pu_stencil_sz, Index::getZeroIndex(DIM));
    for (unsigned int axis = 0, pu_stencil_index = 0; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side, ++pu_stencil_index)
        {
            pu_stencil[pu_stencil_index](axis) = (side == 0 ? 0 : +1);
        }
    }

    // Determine the index ranges.
    tbox::SAMRAI_MPI comm(MPI_COMM_WORLD);
    const int mpi_rank = comm.getRank();
    const int nlocal = num_dofs_per_proc[mpi_rank];
    const int ilower = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int iupper = ilower + nlocal;
    const int ntotal = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.end(), 0);

    // Determine the non-zero structure of the matrix.
    std::vector<int> d_nnz(nlocal, 0), o_nnz(nlocal, 0);
    for (auto p = patch_level->begin(); p != patch_level->end(); ++p)
    {
        auto patch = *p;
        const Box& patch_box = patch->getBox();
        auto u_dof_index_data = BOOST_CAST<SideData<int> >(patch->getPatchData(u_dof_index_idx));
        auto p_dof_index_data = BOOST_CAST<CellData<int> >(patch->getPatchData(p_dof_index_idx));
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (auto b = SideGeometry::begin(patch_box, axis), e = SideGeometry::end(patch_box, axis); b != e; ++b)
            {
                const SideIndex& is = b();
                const CellIndex ic(is.toCell(SideIndex::Upper));
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
        for (auto b = CellGeometry::begin(patch_box), e = CellGeometry::end(patch_box); b != e; ++b)
        {
            const CellIndex& ic = b();
            const int p_dof_index = (*p_dof_index_data)(ic);
            if (UNLIKELY(ilower > p_dof_index || p_dof_index >= iupper)) continue;
            const int p_local_idx = p_dof_index - ilower;
            d_nnz[p_local_idx] += 1;
            for (unsigned int axis = 0, pu_stencil_index = 0; axis < NDIM; ++axis)
            {
                for (int side = 0; side <= 1; ++side, ++pu_stencil_index)
                {
                    const int pu_dof_index =
                        (*u_dof_index_data)(SideIndex(ic + pu_stencil[pu_stencil_index], axis, SideIndex::Lower));
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
    ierr = MatCreateAIJ(PETSC_COMM_WORLD, nlocal, nlocal, PETSC_DETERMINE, PETSC_DETERMINE, PETSC_DEFAULT, &d_nnz[0],
                        PETSC_DEFAULT, &o_nnz[0], &mat);
    IBTK_CHKERRQ(ierr);

    // Set the matrix coefficients.
    const double C = u_problem_coefs.getCConstant();
    const double D = u_problem_coefs.getDConstant();
    for (auto p = patch_level->begin(); p != patch_level->end(); ++p)
    {
        auto patch = *p;
        const Box& patch_box = patch->getBox();
        auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch->getPatchGeometry());
        const double* const dx = pgeom->getDx();

        const IntVector no_ghosts = IntVector::getZero(DIM);
        SideData<double> uu_matrix_coefs(patch_box, uu_stencil_sz, no_ghosts);
        SideData<double> up_matrix_coefs(patch_box, up_stencil_sz, no_ghosts);
        CellData<double> pu_matrix_coefs(patch_box, pu_stencil_sz, no_ghosts);

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
                uu_mat_vals[2 * d + 1] -= 0.5 * D / dx_sq; // lower off-diagonal
                uu_mat_vals[2 * d + 2] -= 0.5 * D / dx_sq; // upper off-diagonal
                uu_mat_vals[0] += 0.5 * D / dx_sq;         // diagonal
                uu_mat_vals[0] += 0.5 * D / dx_sq;         // diagonal
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
        const std::vector<BoundaryBox> physical_codim1_boxes =
            PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
        const int n_physical_codim1_boxes = physical_codim1_boxes.size();
        const double* const patch_x_lower = pgeom->getXLower();
        const double* const patch_x_upper = pgeom->getXUpper();
        const IntVector& ratio_to_level_zero = pgeom->getRatio();
        PatchGeometry::TwoDimBool touches_regular_bdry(DIM), touches_periodic_bdry(DIM);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (int side = 0; side < 2; ++side)
            {
                touches_regular_bdry(axis, side) = pgeom->getTouchesRegularBoundary(axis, side);
                touches_periodic_bdry(axis, side) = pgeom->getTouchesPeriodicBoundary(axis, side);
            }
        }

        // Modify matrix coefficients to account for physical boundary
        // conditions along boundaries which ARE NOT aligned with the data axis.
        //
        // NOTE: We only handle Dirichlet boundary conditions at this point.
        //
        // NOTE: It important to set these values first to avoid problems at
        // corners in the physical domain.  In particular, since Dirichlet
        // boundary conditions for values located on the physical boundary
        // override all other boundary conditions, we set those values last.
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                const BoundaryBox& bdry_box = physical_codim1_boxes[n];
                const unsigned int location_index = bdry_box.getLocationIndex();
                const unsigned int bdry_normal_axis = location_index / 2;
                const bool is_lower = location_index % 2 == 0;

                if (bdry_normal_axis == axis) continue;

                const Box bc_fill_box =
                    pgeom->getBoundaryFillBox(bdry_box, patch_box, /* ghost_width_to_fill */ IntVector::getOne(DIM));
                const BoundaryBox trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
                const Box bc_coef_box = compute_tangential_extension(
                    PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box), axis);

                auto acoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
                auto bcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
                boost::shared_ptr<ArrayData<double> > gcoef_data;

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
                patch->setPatchGeometry(boost::make_shared<CartesianPatchGeometry>(
                    ratio_to_level_zero, touches_regular_bdry, touches_periodic_bdry, dx, shifted_patch_x_lower.data(),
                    shifted_patch_x_upper.data()));

                // Set the boundary condition coefficients.
                static const bool homogeneous_bc = true;
                auto extended_bc_coef = boost::dynamic_pointer_cast<ExtendedRobinBcCoefStrategy>(u_bc_coefs[axis]);
                if (extended_bc_coef)
                {
                    extended_bc_coef->clearTargetPatchDataIndex();
                    extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                }
                u_bc_coefs[axis]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, NULL, *patch, trimmed_bdry_box,
                                             data_time);
                if (gcoef_data && homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

                // Restore the original patch geometry object.
                patch->setPatchGeometry(pgeom);

                // Modify the matrix coefficients to account for homogeneous
                // boundary conditions.
                for (auto bc = CellGeometry::begin(bc_coef_box), e = CellGeometry::end(bc_coef_box); bc != e; ++bc)
                {
                    const CellIndex& i_c = bc();
                    const double& a = (*acoef_data)(i_c, 0);
                    const double& b = (*bcoef_data)(i_c, 0);
                    TBOX_ASSERT(a == 1.0 || MathUtilities<double>::equalEps(a, 1.0));
                    TBOX_ASSERT(b == 0.0 || MathUtilities<double>::equalEps(b, 0.0));

                    Index i_intr = i_c;
                    if (is_lower)
                    {
                        i_intr(bdry_normal_axis) += 0;
                    }
                    else
                    {
                        i_intr(bdry_normal_axis) -= 1;
                    }
                    const SideIndex i_s(i_intr, axis, SideIndex::Lower);
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
            }
        }

        // Modify matrix coefficients to account for physical boundary
        // conditions along boundaries which ARE aligned with the data axis.
        //
        // NOTE: We only handle Dirichlet boundary conditions at this point.
        //
        // NOTE: It important to set these values last to avoid problems at corners
        // in the physical domain.  In particular, since Dirichlet boundary
        // conditions for values located on the physical boundary override all other
        // boundary conditions, we set those values last.
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (int n = 0; n < n_physical_codim1_boxes; ++n)
            {
                const BoundaryBox& bdry_box = physical_codim1_boxes[n];
                const unsigned int location_index = bdry_box.getLocationIndex();
                const unsigned int bdry_normal_axis = location_index / 2;

                if (bdry_normal_axis != axis) continue;

                const Box bc_fill_box =
                    pgeom->getBoundaryFillBox(bdry_box, patch_box, /* ghost_width_to_fill */ IntVector::getOne(DIM));
                const BoundaryBox trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
                const Box bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

                auto acoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
                auto bcoef_data = boost::make_shared<ArrayData<double> >(bc_coef_box, 1);
                boost::shared_ptr<ArrayData<double> > gcoef_data;

                // Set the boundary condition coefficients.
                static const bool homogeneous_bc = true;
                auto extended_bc_coef = boost::dynamic_pointer_cast<ExtendedRobinBcCoefStrategy>(u_bc_coefs[axis]);
                if (extended_bc_coef)
                {
                    extended_bc_coef->clearTargetPatchDataIndex();
                    extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                }
                u_bc_coefs[axis]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, NULL, *patch, trimmed_bdry_box,
                                             data_time);
                if (gcoef_data && homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

                // Modify the matrix coefficients to account for homogeneous
                // boundary conditions.
                for (auto bc = CellGeometry::begin(bc_coef_box), e = CellGeometry::end(bc_coef_box); bc != e; ++bc)
                {
                    const CellIndex& i_c = bc();
                    const SideIndex i_s(i_c, axis, SideIndex::Lower);
                    const double& a = (*acoef_data)(i_c, 0);
                    const double& b = (*bcoef_data)(i_c, 0);
                    TBOX_ASSERT(a == 1.0 || !MathUtilities<double>::equalEps(a, 1.0));
                    TBOX_ASSERT(b == 0.0 || !MathUtilities<double>::equalEps(b, 0.0));
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
            }
        }

        // Set matrix coefficients.
        auto u_dof_index_data = BOOST_CAST<SideData<int> >(patch->getPatchData(u_dof_index_idx));
        auto p_dof_index_data = BOOST_CAST<CellData<int> >(patch->getPatchData(p_dof_index_idx));
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (auto b = SideGeometry::begin(patch_box, axis), e = SideGeometry::end(patch_box, axis); b != e; ++b)
            {
                const SideIndex& is = b();
                const CellIndex ic(is.toCell(SideIndex::Upper));
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

        for (auto b = CellGeometry::begin(patch_box), e = CellGeometry::end(patch_box); b != e; ++b)
        {
            const CellIndex& ic = b();
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
                    p_mat_cols[pu_stencil_index] =
                        (*u_dof_index_data)(SideIndex(ic + pu_stencil[pu_stencil_index], axis, SideIndex::Lower));
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
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
