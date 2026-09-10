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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/StaggeredStokesPETScMatUtilities.h>

#include <ibtk/ExtendedRobinBcCoefStrategy.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/SideSynchCopyFillPattern.h>
#include <ibtk/compiler_hints.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Array.h>
#include <tbox/MathUtilities.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <petsclog.h>
#include <petscmat.h>

#include <ArrayData.h>
#include <BoundaryBox.h>
#include <Box.h>
#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellGeometry.h>
#include <CellIndex.h>
#include <Index.h>
#include <IntVector.h>
#include <MultiblockDataTranslator.h>
#include <Patch.h>
#include <PatchGeometry.h>
#include <PatchLevel.h>
#include <PoissonSpecifications.h>
#include <ProcessorMapping.h>
#include <RefineAlgorithm.h>
#include <RefineOperator.h>
#include <RefineSchedule.h>
#include <RobinBcCoefStrategy.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>
#include <SideVariable.h>
#include <Variable.h>
#include <VariableDatabase.h>
#include <VariableFillPattern.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <limits>
#include <memory>
#include <numeric>
#include <ostream>
#include <set>
#include <utility>
#include <vector>

#include <ibamr/namespaces.h> // IWYU pragma: keep

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
} // namespace

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
    Pointer<PatchLevel<NDIM>> patch_level)
{
    int ierr;
    if (mat)
    {
        ierr = MatDestroy(&mat);
        IBTK_CHKERRQ(ierr);
    }

    // Setup the finite difference stencils.
    static const int uu_stencil_sz = 2 * NDIM + 1;
    std::array<hier::Index<NDIM>, uu_stencil_sz> uu_stencil(
        array_constant<hier::Index<NDIM>, uu_stencil_sz>(hier::Index<NDIM>(0)));
    for (unsigned int axis = 0, uu_stencil_index = 1; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side, ++uu_stencil_index)
        {
            uu_stencil[uu_stencil_index](axis) = (side == 0 ? -1 : +1);
        }
    }
    static const int up_stencil_sz = 2;
    std::array<std::array<hier::Index<NDIM>, up_stencil_sz>, NDIM> up_stencil(
        array_constant<std::array<hier::Index<NDIM>, up_stencil_sz>, NDIM>(
            array_constant<hier::Index<NDIM>, up_stencil_sz>(hier::Index<NDIM>(0))));
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side)
        {
            up_stencil[axis][side](axis) = (side == 0 ? -1 : 0);
        }
    }
    static const int pu_stencil_sz = 2 * NDIM;
    std::array<hier::Index<NDIM>, pu_stencil_sz> pu_stencil(
        array_constant<hier::Index<NDIM>, pu_stencil_sz>(hier::Index<NDIM>(0)));
    for (unsigned int axis = 0, pu_stencil_index = 0; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side, ++pu_stencil_index)
        {
            pu_stencil[pu_stencil_index](axis) = (side == 0 ? 0 : +1);
        }
    }

    // Determine the index ranges.
    const int mpi_rank = IBTK_MPI::getRank();
    const int nlocal = num_dofs_per_proc[mpi_rank];
    const int ilower = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int iupper = ilower + nlocal;
    const int ntotal = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.end(), 0);

    // Determine the non-zero structure of the matrix.
    std::vector<int> d_nnz(nlocal, 0), o_nnz(nlocal, 0);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM, int>> u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int>> p_dof_index_data = patch->getPatchData(p_dof_index_idx);
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
                        get_data_or_null(d_nnz),
                        0,
                        get_data_or_null(o_nnz),
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
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
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
        const Array<BoundaryBox<NDIM>> physical_codim1_boxes =
            PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
        const int n_physical_codim1_boxes = physical_codim1_boxes.size();
        const double* const patch_x_lower = pgeom->getXLower();
        const double* const patch_x_upper = pgeom->getXUpper();
        const IntVector<NDIM>& ratio_to_level_zero = pgeom->getRatio();
        Array<Array<bool>> touches_regular_bdry(NDIM), touches_periodic_bdry(NDIM);
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

                Pointer<ArrayData<NDIM, double>> acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double>> bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double>> gcoef_data;

                // Temporarily reset the patch geometry object associated with
                // the patch so that boundary conditions are set at the correct
                // spatial locations.
                std::array<double, NDIM> shifted_patch_x_lower, shifted_patch_x_upper;
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
                auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(u_bc_coefs[axis]);
                if (extended_bc_coef)
                {
                    extended_bc_coef->clearTargetPatchDataIndex();
                    extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                }
                u_bc_coefs[axis]->setBcCoefs(
                    acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
                if (gcoef_data && homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

                // Restore the original patch geometry object.
                patch->setPatchGeometry(pgeom);

                // Modify the matrix coefficients to account for homogeneous
                // boundary conditions.
                for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
                {
                    const hier::Index<NDIM>& i = bc();
                    const double& a = (*acoef_data)(i, 0);
                    const double& b = (*bcoef_data)(i, 0);
                    const bool velocity_bc = (a == 1.0 || IBTK::rel_equal_eps(a, 1.0));
                    const bool traction_bc = (b == 1.0 || IBTK::rel_equal_eps(b, 1.0));
#if !defined(NDEBUG)
                    TBOX_ASSERT((velocity_bc || traction_bc) && !(velocity_bc && traction_bc));
#endif
                    hier::Index<NDIM> i_intr = i;
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
                            "StaggeredStokesPETScMatUtilities::"
                            "constructPatchLevelMACStokesOp(): Unknown BC type for "
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

                Pointer<ArrayData<NDIM, double>> acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double>> bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double>> gcoef_data;

                // Set the boundary condition coefficients.
                static const bool homogeneous_bc = true;
                auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(u_bc_coefs[axis]);
                if (extended_bc_coef)
                {
                    extended_bc_coef->clearTargetPatchDataIndex();
                    extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                }
                u_bc_coefs[axis]->setBcCoefs(
                    acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
                if (gcoef_data && homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

                // Modify the matrix coefficients to account for homogeneous
                // boundary conditions.
                for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
                {
                    const hier::Index<NDIM>& i = bc();
                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                    const double& a = (*acoef_data)(i, 0);
                    const double& b = (*bcoef_data)(i, 0);
                    const bool velocity_bc = (a == 1.0 || IBTK::rel_equal_eps(a, 1.0));
                    const bool traction_bc = (b == 1.0 || IBTK::rel_equal_eps(b, 1.0));
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
                            "StaggeredStokesPETScMatUtilities::"
                            "constructPatchLevelMACStokesOp(): Unknown BC type for "
                            "normal velocity specified.");
                    }
                }
            }
        }

        // Set matrix coefficients.
        Pointer<SideData<NDIM, int>> u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int>> p_dof_index_data = patch->getPatchData(p_dof_index_idx);
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

                ierr = MatSetValues(
                    mat, 1, &u_dof_index, u_stencil_sz, u_mat_cols.data(), u_mat_vals.data(), INSERT_VALUES);
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

            ierr =
                MatSetValues(mat, 1, &p_dof_index, p_stencil_sz, p_mat_cols.data(), p_mat_vals.data(), INSERT_VALUES);
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
StaggeredStokesPETScMatUtilities::constructPatchLevelASMSubdomains(std::vector<std::set<int>>& is_overlap,
                                                                   std::vector<std::set<int>>& is_nonoverlap,
                                                                   const IntVector<NDIM>& box_size,
                                                                   const IntVector<NDIM>& overlap_size,
                                                                   const std::vector<int>& /*num_dofs_per_proc*/,
                                                                   int u_dof_index_idx,
                                                                   int p_dof_index_idx,
                                                                   Pointer<PatchLevel<NDIM>> patch_level,
                                                                   Pointer<CoarseFineBoundary<NDIM>> /*cf_boundary*/)
{
    // Clear previously stored index sets.
    for (auto& k : is_overlap)
    {
        k.clear();
    }
    is_overlap.clear();
    for (auto& k : is_nonoverlap)
    {
        k.clear();
    }
    is_nonoverlap.clear();

    // Create variables to keep track of whether a particular velocity location
    // is the "master" location.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<SideVariable<NDIM, int>> patch_num_var = new SideVariable<NDIM, int>(
        "StaggeredStokesPETScMatUtilities::constructPatchLevelASMSubdomains()::"
        "patch_num_var");
    static const int patch_num_idx = var_db->registerPatchDataIndex(patch_num_var);
    patch_level->allocatePatchData(patch_num_idx);
    Pointer<SideVariable<NDIM, bool>> u_mastr_loc_var = new SideVariable<NDIM, bool>(
        "StaggeredStokesPETScMatUtilities::"
        "constructPatchLevelASMSubdomains()::u_"
        "mastr_loc_var");
    static const int u_mastr_loc_idx = var_db->registerPatchDataIndex(u_mastr_loc_var);
    patch_level->allocatePatchData(u_mastr_loc_idx);
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const int patch_num = patch->getPatchNumber();
        Pointer<SideData<NDIM, int>> patch_num_data = patch->getPatchData(patch_num_idx);
        Pointer<SideData<NDIM, bool>> u_mastr_loc_data = patch->getPatchData(u_mastr_loc_idx);
        patch_num_data->fillAll(patch_num);
        u_mastr_loc_data->fillAll(false);
    }

    // Synchronize the patch number at patch boundaries to determine which patch
    // owns a given DOF along patch boundaries.
    RefineAlgorithm<NDIM> bdry_synch_alg;
    bdry_synch_alg.registerRefine(patch_num_idx, patch_num_idx, patch_num_idx, nullptr, new SideSynchCopyFillPattern());
    bdry_synch_alg.createSchedule(patch_level)->fillData(0.0);

    // For a single patch in a periodic domain, the far side DOFs are not master.
    Pointer<CartesianGridGeometry<NDIM>> grid_geom = patch_level->getGridGeometry();
    IntVector<NDIM> periodic_shift = grid_geom->getPeriodicShift(patch_level->getRatio());
    const BoxArray<NDIM>& domain_boxes = patch_level->getPhysicalDomain();
#if !defined(NDEBUG)
    TBOX_ASSERT(domain_boxes.size() == 1);
#endif
    const hier::Index<NDIM>& domain_upper = domain_boxes[0].upper();

    // Determine the number of local DOFs.
    int local_dof_count = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const int patch_num = patch->getPatchNumber();
        const Box<NDIM>& patch_box = patch->getBox();
        const IntVector<NDIM> patch_size = patch_box.numberCells();

        Pointer<SideData<NDIM, int>> u_dof_index_data = patch->getPatchData(u_dof_index_idx);
        Pointer<SideData<NDIM, int>> patch_num_data = patch->getPatchData(patch_num_idx);
        Pointer<SideData<NDIM, bool>> u_mastr_loc_data = patch->getPatchData(u_mastr_loc_idx);
        for (unsigned int component_axis = 0; component_axis < NDIM; ++component_axis)
        {
            const int upper_domain_side_idx = domain_upper(component_axis) + 1;
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, component_axis)); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                const SideIndex<NDIM> is(i, component_axis, SideIndex<NDIM>::Lower);
                bool fully_periodic_patch_in_axis =
                    periodic_shift(component_axis) && (patch_size(component_axis) == periodic_shift(component_axis));
                bool periodic_image = fully_periodic_patch_in_axis && (i(component_axis) == upper_domain_side_idx);
                if ((*patch_num_data)(is) == patch_num && !periodic_image)
                {
                    (*u_mastr_loc_data)(is) = true;
                    ++local_dof_count;
                }
            }
        }
        local_dof_count += CellGeometry<NDIM>::toCellBox(patch_box).size();
    }

    // Determine the subdomains associated with this processor.
    const int n_local_patches = patch_level->getProcessorMapping().getNumberOfLocalIndices();
    std::vector<std::vector<Box<NDIM>>> overlap_boxes(n_local_patches), nonoverlap_boxes(n_local_patches);
    int patch_counter = 0, subdomain_counter = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++, ++patch_counter)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        IndexUtilities::partitionPatchBox(
            overlap_boxes[patch_counter], nonoverlap_boxes[patch_counter], patch_box, box_size, overlap_size);
        subdomain_counter += overlap_boxes[patch_counter].size();
    }
    is_overlap.resize(subdomain_counter);
    is_nonoverlap.resize(subdomain_counter);

    // Fill in the IS'es.
    int nonoverlap_dof_counter = 0;
    subdomain_counter = 0, patch_counter = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++, ++patch_counter)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Box<NDIM> side_patch_box[NDIM];
        for (int axis = 0; axis < NDIM; ++axis)
        {
            side_patch_box[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
        }
        Pointer<SideData<NDIM, bool>> u_mastr_loc_data = patch->getPatchData(u_mastr_loc_idx);
        Pointer<SideData<NDIM, int>> u_dof_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int>> p_dof_data = patch->getPatchData(p_dof_index_idx);
#if !defined(NDEBUG)
        {
            const int u_data_depth = u_dof_data->getDepth();
            const int p_data_depth = p_dof_data->getDepth();
            TBOX_ASSERT(u_data_depth == 1);
            TBOX_ASSERT(p_data_depth == 1);
            TBOX_ASSERT(u_dof_data->getGhostCellWidth().min() >= overlap_size.max());
            TBOX_ASSERT(p_dof_data->getGhostCellWidth().min() >= overlap_size.max());
        }
#endif
        int n_patch_subdomains = static_cast<int>(nonoverlap_boxes[patch_counter].size());
        for (int k = 0; k < n_patch_subdomains; ++k, ++subdomain_counter)
        {
            // The nonoverlapping subdomains.
            const Box<NDIM>& sub_box = nonoverlap_boxes[patch_counter][k];
            Box<NDIM> side_sub_box[NDIM];
            for (int axis = 0; axis < NDIM; ++axis)
            {
                side_sub_box[axis] = SideGeometry<NDIM>::toSideBox(sub_box, axis);
            }

            // Get the local DOFs.
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator b(side_sub_box[axis]); b; b++)
                {
                    const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
                    const bool at_upper_subdomain_bdry = (i_s(axis) == side_sub_box[axis].upper(axis));
                    const bool at_upper_patch_bdry = (i_s(axis) == side_patch_box[axis].upper(axis));
                    if (!at_upper_subdomain_bdry || (at_upper_patch_bdry && (*u_mastr_loc_data)(i_s)))
                    {
                        const int dof_idx = (*u_dof_data)(i_s);
                        if (dof_idx >= 0) is_nonoverlap[subdomain_counter].insert(dof_idx);
                    }
                }
            }
            for (Box<NDIM>::Iterator b(sub_box); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                const int dof_idx = (*p_dof_data)(i);
                if (dof_idx >= 0) is_nonoverlap[subdomain_counter].insert(dof_idx);
            }
            const int n_nonoverlap = static_cast<int>(is_nonoverlap[subdomain_counter].size());
            nonoverlap_dof_counter += n_nonoverlap;

            // The overlapping subdomains.
            const Box<NDIM>& overlap_sub_box = overlap_boxes[patch_counter][k];
            Box<NDIM> side_overlap_sub_box[NDIM];
            for (int axis = 0; axis < NDIM; ++axis)
            {
                side_overlap_sub_box[axis] = SideGeometry<NDIM>::toSideBox(overlap_sub_box, axis);
            }

            // Get the overlap DOFs.
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator b(side_overlap_sub_box[axis]); b; b++)
                {
                    const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
                    const int dof_idx = (*u_dof_data)(i_s);
                    if (dof_idx >= 0) is_overlap[subdomain_counter].insert(dof_idx);
                }
            }
            for (Box<NDIM>::Iterator b(overlap_sub_box); b; b++)
            {
                const CellIndex<NDIM>& i = b();
                const int dof_idx = (*p_dof_data)(i);
                if (dof_idx >= 0) is_overlap[subdomain_counter].insert(dof_idx);
            }
        }
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(local_dof_count == nonoverlap_dof_counter);
#else
    NULL_USE(local_dof_count);
    NULL_USE(nonoverlap_dof_counter);
#endif

    // Deallocate patch_num variable data.
    patch_level->deallocatePatchData(patch_num_idx);
    patch_level->deallocatePatchData(u_mastr_loc_idx);
    return;
} // constructPatchLevelASMSubdomains

namespace
{
using CouplingAwareASMSeedRecord = std::pair<std::array<int, NDIM>, int>;

// Rows of sorted, distinct integers in two flat arrays, numbered by DOF. On one rank the DOFs are the integers
// from 0 up, so a row is found by indexing, without hashing or pointer chasing. The entries of a row are the range
// begin(row) to end(row), which is empty for a row that was never given an entry.
class DofRows
{
public:
    // Build the rows from (row, entry) pairs. A pair that occurs more than once gives one entry.
    void build(const std::vector<std::pair<int, int>>& pairs, const int n_rows)
    {
        d_offsets.assign(n_rows + 1, 0);
        for (const auto& pair : pairs) ++d_offsets[pair.first + 1];
        std::partial_sum(d_offsets.begin(), d_offsets.end(), d_offsets.begin());
        d_entries.resize(pairs.size());
        std::vector<int> cursor(d_offsets.begin(), d_offsets.end() - 1);
        for (const auto& pair : pairs) d_entries[cursor[pair.first]++] = pair.second;
        int write = 0;
        for (int row = 0; row < n_rows; ++row)
        {
            const auto first = d_entries.begin() + d_offsets[row];
            const auto last = d_entries.begin() + d_offsets[row + 1];
            std::sort(first, last);
            const int start = write;
            for (auto entry = first; entry != last; ++entry)
            {
                if (write == start || d_entries[write - 1] != *entry) d_entries[write++] = *entry;
            }
            d_offsets[row] = start;
        }
        d_offsets[n_rows] = write;
        d_entries.resize(write);
    }

    const int* begin(const int row) const
    {
        return d_entries.data() + (row < static_cast<int>(d_offsets.size()) - 1 ? d_offsets[row] : d_entries.size());
    }

    const int* end(const int row) const
    {
        return d_entries.data() +
               (row < static_cast<int>(d_offsets.size()) - 1 ? d_offsets[row + 1] : d_entries.size());
    }

private:
    std::vector<int> d_offsets, d_entries;
};

// Sort a vector and remove repeated entries.
void
sort_unique(std::vector<int>& values)
{
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
}

// Sort the records geometrically, drop repeats of a DOF (a periodic or shared side is found more than once), and keep
// every seed_stride-th of the rest.
std::vector<int>
select_coupling_aware_seeds(std::vector<CouplingAwareASMSeedRecord> records, const int n_dofs, const int seed_stride)
{
    std::sort(records.begin(), records.end());
    std::vector<char> seen(n_dofs, 0);
    int n_seen = 0;
    std::vector<int> seeds;
    for (const CouplingAwareASMSeedRecord& record : records)
    {
        if (seen[record.second])
        {
            continue;
        }
        seen[record.second] = 1;
        if (n_seen++ % seed_stride == 0)
        {
            seeds.push_back(record.second);
        }
    }
    return seeds;
}

// Matrix entries of a row at or below this size are not couplings: the larger of a roundoff bound for the number of
// entries that were compared and the relative tolerance, times the largest entry of the row.
double
coupling_threshold(const PetscInt entries, const double row_max, const double relative_zero_tol)
{
    return std::max(entries * std::numeric_limits<double>::epsilon(), relative_zero_tol) * row_max;
}

// Require a square matrix that has the full coupled numbering and the ownership of this rank.
void
require_full_coupled_numbering(Mat matrix, const std::vector<int>& num_dofs_per_proc, const char* const name)
{
    PetscInt nrows = 0, ncols = 0, first = 0, last = 0;
    int ierr = MatGetSize(matrix, &nrows, &ncols);
    IBTK_CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(matrix, &first, &last);
    IBTK_CHKERRQ(ierr);
    const int rank = IBTK_MPI::getRank();
    const int local_begin = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + rank, 0);
    const int local_end = local_begin + num_dofs_per_proc[rank];
    if (nrows != ncols || nrows != std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.end(), 0) ||
        first != local_begin || last != local_end)
    {
        TBOX_ERROR("require_full_coupled_numbering():\n"
                   << "  " << name << " must use full coupled numbering and ownership.\n");
    }
}

/*! \brief Level geometry for coupling-aware patch construction, built from live DOF data.
 *
 * The supplied patch data must outlive this object. Only one MPI rank is supported. Construction semantics are
 * defined by
 * StaggeredStokesPETScMatUtilities::construct_patch_level_coupling_aware_asm_subdomains()
 * and StaggeredStokesPETScMatUtilities::construct_patch_level_pressure_cell_seeded_cav_patches().
 */
class CouplingAwareASMSubdomains
{
public:
    /*! \brief Build velocity adjacency and cell closures from live DOF data. */
    CouplingAwareASMSubdomains(int u_idx, int p_idx, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level);

    /*! \brief Construct ordered overlap sets and their first-owner partition. */
    void constructSubdomains(std::vector<std::set<int>>& overlap,
                             std::vector<std::set<int>>& nonoverlap,
                             const std::vector<int>& num_dofs_per_proc,
                             Mat matrix,
                             int seed_axis,
                             int seed_stride,
                             CouplingAwareASMSeedTraversalOrder order,
                             CouplingAwareASMClosurePolicy policy,
                             double relative_zero_tol);

    /*! \brief Construct pressure patches as defined by
     * StaggeredStokesPETScMatUtilities::construct_patch_level_pressure_cell_seeded_cav_patches().
     */
    void constructPressureCellPatches(std::vector<std::set<int>>& patches,
                                      std::vector<int>& pressure_seeds,
                                      const std::vector<int>& num_dofs_per_proc,
                                      Mat elasticity,
                                      int seed_stride,
                                      CouplingAwareASMSeedTraversalOrder order,
                                      CouplingAwareASMClosurePolicy policy,
                                      double relative_zero_tol);

private:
    /*! \brief Add lower-face component pairing when STRICT first needs it. */
    void buildSeedPairs();

    /*! \brief Join the standard Vanka patches of the cells that touch the expanded velocity DOFs and of extra_cells.
     *
     * The STRICT policy skips a cell that has a velocity DOF outside expanded. The RELAXED policy also adds expanded.
     */
    std::vector<int> closeExpandedDOFs(const std::vector<int>& expanded,
                                       std::vector<int> extra_cells,
                                       CouplingAwareASMClosurePolicy policy) const;

    bool isVelocity(const int dof) const
    {
        return dof >= 0 && dof < d_n_dofs && d_is_velocity[dof];
    }

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> d_level;
    int d_u_idx, d_p_idx;
    // The velocity DOFs are flagged by DOF, and the rows are numbered by DOF, up to d_n_dofs.
    int d_n_dofs = 0;
    std::vector<char> d_is_velocity;
    DofRows d_adjacent_cells, d_cell_closures, d_seed_pairs;
    bool d_pairs_built = false;
};

std::array<int, NDIM>
get_coupling_aware_axis_order(const CouplingAwareASMSeedTraversalOrder order)
{
    std::array<int, NDIM> axis_order{};
#if (NDIM == 2)
    if (order == CouplingAwareASMSeedTraversalOrder::I_J)
    {
        axis_order = { 0, 1 };
    }
    else if (order == CouplingAwareASMSeedTraversalOrder::J_I)
    {
        axis_order = { 1, 0 };
    }
#else
    if (order == CouplingAwareASMSeedTraversalOrder::I_J_K)
    {
        axis_order = { 0, 1, 2 };
    }
    else if (order == CouplingAwareASMSeedTraversalOrder::J_K_I)
    {
        axis_order = { 1, 2, 0 };
    }
    else if (order == CouplingAwareASMSeedTraversalOrder::K_I_J)
    {
        axis_order = { 2, 0, 1 };
    }
#endif
    else
    {
        TBOX_ERROR("get_coupling_aware_axis_order():\n"
                   << "  invalid logical traversal order.\n");
    }
    return axis_order;
}
} // namespace

CouplingAwareASMSubdomains::CouplingAwareASMSubdomains(const int u_idx,
                                                       const int p_idx,
                                                       Pointer<PatchLevel<NDIM>> level)
    : d_level(level), d_u_idx(u_idx), d_p_idx(p_idx)
{
    if (!level || u_idx < 0 || p_idx < 0 || !level->checkAllocated(u_idx) || !level->checkAllocated(p_idx))
    {
        TBOX_ERROR("CouplingAwareASMSubdomains::CouplingAwareASMSubdomains():\n"
                   << "  allocated level DOF data are required.\n");
    }
    std::vector<std::pair<int, int>> adjacent_cells, cell_closures;
    std::vector<int> velocities;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, int>> u = patch->getPatchData(u_idx);
        Pointer<CellData<NDIM, int>> pressure = patch->getPatchData(p_idx);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            const Box<NDIM> sides = SideGeometry<NDIM>::toSideBox(patch->getBox(), axis);
            for (Box<NDIM>::Iterator b(sides); b; b++)
            {
                const SideIndex<NDIM> side(b(), axis, SideIndex<NDIM>::Lower);
                const int velocity = (*u)(side);
                if (velocity < 0)
                {
                    continue;
                }
                d_n_dofs = std::max(d_n_dofs, velocity + 1);
                velocities.push_back(velocity);
                for (int face = 0; face < 2; ++face)
                {
                    const int cell = (*pressure)(side.toCell(face));
                    if (cell >= 0)
                    {
                        adjacent_cells.emplace_back(velocity, cell);
                    }
                }
            }
        }
        for (Box<NDIM>::Iterator b(patch->getBox()); b; b++)
        {
            const int cell = (*pressure)(b());
            if (cell < 0)
            {
                continue;
            }
            d_n_dofs = std::max(d_n_dofs, cell + 1);
            cell_closures.emplace_back(cell, cell);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (int face = 0; face < 2; ++face)
                {
                    const int velocity = (*u)(SideIndex<NDIM>(b(), axis, face));
                    if (velocity >= 0)
                    {
                        d_n_dofs = std::max(d_n_dofs, velocity + 1);
                        cell_closures.emplace_back(cell, velocity);
                    }
                }
            }
        }
    }
    d_is_velocity.assign(d_n_dofs, 0);
    for (const int velocity : velocities) d_is_velocity[velocity] = 1;
    d_adjacent_cells.build(adjacent_cells, d_n_dofs);
    d_cell_closures.build(cell_closures, d_n_dofs);
}

void
CouplingAwareASMSubdomains::buildSeedPairs()
{
    std::vector<std::pair<int, int>> seed_pairs;
    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = d_level->getPatch(p());
        Pointer<SideData<NDIM, int>> u = patch->getPatchData(d_u_idx);
        for (Box<NDIM>::Iterator b(patch->getBox()); b; b++)
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                const int seed = (*u)(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Lower));
                if (seed < 0)
                {
                    continue;
                }
                for (int other = 0; other < NDIM; ++other)
                {
                    const int paired = (*u)(SideIndex<NDIM>(b(), other, SideIndex<NDIM>::Lower));
                    if (other != axis && paired >= 0)
                    {
                        seed_pairs.emplace_back(seed, paired);
                    }
                }
            }
        }
    }
    d_seed_pairs.build(seed_pairs, d_n_dofs);
    d_pairs_built = true;
}

std::vector<int>
CouplingAwareASMSubdomains::closeExpandedDOFs(const std::vector<int>& expanded,
                                              std::vector<int> cells,
                                              const CouplingAwareASMClosurePolicy policy) const
{
    for (const int velocity : expanded)
    {
        cells.insert(cells.end(), d_adjacent_cells.begin(velocity), d_adjacent_cells.end(velocity));
    }
    sort_unique(cells);
    std::vector<int> closure;
    for (const int cell : cells)
    {
        const int* const first = d_cell_closures.begin(cell);
        const int* const last = d_cell_closures.end(cell);
        if (policy == CouplingAwareASMClosurePolicy::STRICT &&
            std::any_of(first,
                        last,
                        [&](const int dof)
                        { return isVelocity(dof) && !std::binary_search(expanded.begin(), expanded.end(), dof); }))
        {
            continue;
        }
        closure.insert(closure.end(), first, last);
    }
    if (policy == CouplingAwareASMClosurePolicy::RELAXED)
    {
        closure.insert(closure.end(), expanded.begin(), expanded.end());
    }
    sort_unique(closure);
    return closure;
}

void
CouplingAwareASMSubdomains::constructSubdomains(std::vector<std::set<int>>& overlap,
                                                std::vector<std::set<int>>& nonoverlap,
                                                const std::vector<int>& num_dofs_per_proc,
                                                Mat matrix,
                                                const int seed_axis,
                                                const int seed_stride,
                                                const CouplingAwareASMSeedTraversalOrder order,
                                                const CouplingAwareASMClosurePolicy policy,
                                                const double relative_zero_tol)
{
    if (IBTK_MPI::getNodes() != 1)
    {
        TBOX_ERROR("CouplingAwareASMSubdomains::constructSubdomains():\n"
                   << "  velocity-seeded construction requires one MPI rank.\n");
    }
    const std::array<int, NDIM> axis_order = get_coupling_aware_axis_order(order);
    if (seed_axis < 0 || seed_axis >= NDIM || seed_stride < 1 ||
        (policy != CouplingAwareASMClosurePolicy::RELAXED && policy != CouplingAwareASMClosurePolicy::STRICT) ||
        !std::isfinite(relative_zero_tol) || relative_zero_tol < 0.0 || !matrix ||
        num_dofs_per_proc.size() != static_cast<std::size_t>(IBTK_MPI::getNodes()))
    {
        TBOX_ERROR("CouplingAwareASMSubdomains::constructSubdomains():\n"
                   << "  invalid construction arguments.\n");
    }
    require_full_coupled_numbering(matrix, num_dofs_per_proc, "matrix");
    if (d_n_dofs > num_dofs_per_proc.front())
    {
        TBOX_ERROR("CouplingAwareASMSubdomains::constructSubdomains():\n"
                   << "  the DOF data numbers DOFs up to " << d_n_dofs - 1 << ", but the matrix has only "
                   << num_dofs_per_proc.front() << " DOFs.\n");
    }
    int ierr;
    if (policy == CouplingAwareASMClosurePolicy::STRICT && !d_pairs_built)
    {
        buildSeedPairs();
    }
    // Sort geometric records before removing periodic/shared-side duplicates.
    std::vector<CouplingAwareASMSeedRecord> records;
    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = d_level->getPatch(p());
        Pointer<SideData<NDIM, int>> u = patch->getPatchData(d_u_idx);
        const Box<NDIM> sides = SideGeometry<NDIM>::toSideBox(patch->getBox(), seed_axis);
        for (Box<NDIM>::Iterator b(sides); b; b++)
        {
            const int seed = (*u)(SideIndex<NDIM>(b(), seed_axis, SideIndex<NDIM>::Lower));
            if (seed < 0)
            {
                continue;
            }
            std::array<int, NDIM> index{};
            for (int d = 0; d < NDIM; ++d)
            {
                index[d] = b()(axis_order[d]);
            }
            records.emplace_back(index, seed);
        }
    }
    overlap.clear();
    for (const int seed_dof : select_coupling_aware_seeds(std::move(records), d_n_dofs, seed_stride))
    {
        std::vector<int> seeds{ seed_dof };
        if (policy == CouplingAwareASMClosurePolicy::STRICT)
        {
            seeds.insert(seeds.end(), d_seed_pairs.begin(seed_dof), d_seed_pairs.end(seed_dof));
        }
        std::vector<int> expanded = seeds;
        for (const int seed : seeds)
        {
            PetscInt count = 0;
            const PetscInt* columns = nullptr;
            const PetscScalar* values = nullptr;
            ierr = MatGetRow(matrix, seed, &count, &columns, &values);
            IBTK_CHKERRQ(ierr);
            PetscInt velocity_count = 0;
            double row_max = 0.0;
            for (PetscInt k = 0; k < count; ++k)
            {
                if (isVelocity(columns[k]))
                {
                    ++velocity_count; // Stored velocity zeros contribute to the roundoff threshold.
                    row_max = std::max(row_max, static_cast<double>(PetscAbsScalar(values[k])));
                }
            }
            const double threshold = coupling_threshold(velocity_count, row_max, relative_zero_tol);
            for (PetscInt k = 0; k < count; ++k)
            {
                if (isVelocity(columns[k]) && PetscAbsScalar(values[k]) > threshold)
                {
                    expanded.push_back(columns[k]);
                }
            }
            ierr = MatRestoreRow(matrix, seed, &count, &columns, &values);
            IBTK_CHKERRQ(ierr);
        }
        sort_unique(expanded);
        std::vector<int> closure = closeExpandedDOFs(expanded, {}, policy);
        if (policy == CouplingAwareASMClosurePolicy::RELAXED && closure.size() < 2 * NDIM + 1)
        {
            TBOX_ERROR("CouplingAwareASMSubdomains::constructSubdomains():\n"
                       << "  incomplete relaxed cell closure.\n");
        }
        overlap.emplace_back(closure.begin(), closure.end());
    }
    nonoverlap.assign(overlap.size(), {});
    std::vector<char> assigned(d_n_dofs, 0);
    for (std::size_t k = 0; k < overlap.size(); ++k)
    {
        for (const int dof : overlap[k])
        {
            if (!assigned[dof])
            {
                assigned[dof] = 1;
                nonoverlap[k].insert(nonoverlap[k].end(), dof);
            }
        }
    }
}

void
CouplingAwareASMSubdomains::constructPressureCellPatches(std::vector<std::set<int>>& patches,
                                                         std::vector<int>& pressure_seeds,
                                                         const std::vector<int>& num_dofs_per_proc,
                                                         Mat elasticity,
                                                         const int seed_stride,
                                                         const CouplingAwareASMSeedTraversalOrder order,
                                                         const CouplingAwareASMClosurePolicy policy,
                                                         const double relative_zero_tol)
{
    if (IBTK_MPI::getNodes() != 1)
    {
        TBOX_ERROR("CouplingAwareASMSubdomains::constructPressureCellPatches():\n"
                   << "  pressure-cell CAV construction requires one MPI rank.\n");
    }
    const std::array<int, NDIM> axis_order = get_coupling_aware_axis_order(order);
    if (!elasticity || num_dofs_per_proc.size() != 1 || seed_stride < 1 ||
        (policy != CouplingAwareASMClosurePolicy::RELAXED && policy != CouplingAwareASMClosurePolicy::STRICT) ||
        !std::isfinite(relative_zero_tol) || relative_zero_tol < 0.0)
    {
        TBOX_ERROR("CouplingAwareASMSubdomains::constructPressureCellPatches():\n"
                   << "  invalid pressure-cell CAV construction arguments.\n");
    }
    require_full_coupled_numbering(elasticity, num_dofs_per_proc, "elasticity");
    const PetscInt rows = num_dofs_per_proc.front();
    if (d_n_dofs > rows)
    {
        TBOX_ERROR("CouplingAwareASMSubdomains::constructPressureCellPatches():\n"
                   << "  the DOF data numbers DOFs up to " << d_n_dofs - 1 << ", but the elasticity matrix has only "
                   << rows << " DOFs.\n");
    }
    int ierr;

    // Each retained entry adds both directions without materializing a transpose.
    std::vector<std::pair<int, int>> adjacency_pairs;
    for (PetscInt row = 0; row < rows; ++row)
    {
        PetscInt count = 0;
        const PetscInt* cols = nullptr;
        const PetscScalar* values = nullptr;
        ierr = MatGetRow(elasticity, row, &count, &cols, &values);
        IBTK_CHKERRQ(ierr);
        double row_max = 0.0;
        for (PetscInt k = 0; k < count; ++k)
        {
            row_max = std::max(row_max, static_cast<double>(PetscAbsScalar(values[k])));
        }
        const double threshold = coupling_threshold(count, row_max, relative_zero_tol);
        bool invalid_pressure_entry = false;
        for (PetscInt k = 0; k < count; ++k)
        {
            if (PetscAbsScalar(values[k]) <= threshold)
            {
                continue;
            }
            if (!isVelocity(row) || !isVelocity(cols[k]))
            {
                invalid_pressure_entry = true;
                break;
            }
            adjacency_pairs.emplace_back(row, cols[k]);
            adjacency_pairs.emplace_back(cols[k], row);
        }
        ierr = MatRestoreRow(elasticity, row, &count, &cols, &values);
        IBTK_CHKERRQ(ierr);
        if (invalid_pressure_entry)
        {
            TBOX_ERROR("CouplingAwareASMSubdomains::constructPressureCellPatches():\n"
                       << "  elasticity pressure rows and columns must be numerically zero.\n");
        }
    }
    DofRows adjacency;
    adjacency.build(adjacency_pairs, rows);

    std::vector<CouplingAwareASMSeedRecord> records;
    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = d_level->getPatch(p());
        Pointer<CellData<NDIM, int>> pressure = patch->getPatchData(d_p_idx);
        for (Box<NDIM>::Iterator b(patch->getBox()); b; b++)
        {
            const int seed = (*pressure)(b());
            if (seed < 0)
            {
                continue;
            }
            std::array<int, NDIM> index{};
            for (int d = 0; d < NDIM; ++d)
            {
                index[d] = b()(axis_order[d]);
            }
            records.emplace_back(index, seed);
        }
    }
    patches.clear();
    pressure_seeds.clear();
    for (const int seed : select_coupling_aware_seeds(std::move(records), d_n_dofs, seed_stride))
    {
        const int* const standard_first = d_cell_closures.begin(seed);
        const int* const standard_last = d_cell_closures.end(seed);
        std::vector<int> velocities;
        std::copy_if(standard_first,
                     standard_last,
                     std::back_inserter(velocities),
                     [&](const int dof) { return isVelocity(dof); });
        if (velocities.size() != 2 * NDIM)
        {
            TBOX_ERROR("CouplingAwareASMSubdomains::constructPressureCellPatches():\n"
                       << "  a pressure seed requires a complete MAC velocity stencil.\n");
        }
        std::vector<int> expanded = velocities;
        for (const int velocity : velocities)
        {
            expanded.insert(expanded.end(), adjacency.begin(velocity), adjacency.end(velocity));
        }
        sort_unique(expanded);
        pressure_seeds.push_back(seed);
        // Even RELAXED must not close neighboring cells when elasticity adds nothing.
        if (expanded == velocities)
        {
            patches.emplace_back(standard_first, standard_last);
            continue;
        }
        const std::vector<int> closure = closeExpandedDOFs(expanded, { seed }, policy);
        patches.emplace_back(closure.begin(), closure.end());
    }
}

void
StaggeredStokesPETScMatUtilities::construct_patch_level_coupling_aware_asm_subdomains(
    std::vector<std::set<int>>& overlap,
    std::vector<std::set<int>>& nonoverlap,
    const std::vector<int>& num_dofs_per_proc,
    const int u_idx,
    const int p_idx,
    Pointer<PatchLevel<NDIM>> level,
    Mat matrix,
    const int seed_axis,
    const int seed_stride,
    const CouplingAwareASMSeedTraversalOrder order,
    const CouplingAwareASMClosurePolicy policy,
    const double relative_zero_tol)
{
    CouplingAwareASMSubdomains maps(u_idx, p_idx, level);
    maps.constructSubdomains(
        overlap, nonoverlap, num_dofs_per_proc, matrix, seed_axis, seed_stride, order, policy, relative_zero_tol);
}

void
StaggeredStokesPETScMatUtilities::construct_patch_level_pressure_cell_seeded_cav_patches(
    std::vector<std::set<int>>& patches,
    std::vector<int>& pressure_seeds,
    const std::vector<int>& num_dofs_per_proc,
    const int u_idx,
    const int p_idx,
    Pointer<PatchLevel<NDIM>> level,
    Mat elasticity,
    const int seed_stride,
    const CouplingAwareASMSeedTraversalOrder order,
    const CouplingAwareASMClosurePolicy policy,
    const double relative_zero_tol)
{
    CouplingAwareASMSubdomains maps(u_idx, p_idx, level);
    maps.constructPressureCellPatches(
        patches, pressure_seeds, num_dofs_per_proc, elasticity, seed_stride, order, policy, relative_zero_tol);
}

void
StaggeredStokesPETScMatUtilities::constructPatchLevelFields(
    std::vector<std::set<int>>& is_field,
    std::vector<std::string>& is_field_name,
    const std::vector<int>& num_dofs_per_proc,
    int u_dof_index_idx,
    int p_dof_index_idx,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level)
{
    // Destroy the previously stored IS'es
    for (auto& k : is_field)
    {
        k.clear();
    }
    is_field.clear();
    is_field_name.clear();

    // Resize vectors
    is_field.resize(2);
    is_field_name.resize(2);

    // Name of the fields.
    static const int U_FIELD_IDX = 0;
    static const int P_FIELD_IDX = 1;
    is_field_name[U_FIELD_IDX] = "velocity";
    is_field_name[P_FIELD_IDX] = "pressure";

    // DOFs on this processor.
    const int mpi_rank = IBTK_MPI::getRank();
    const int n_local_dofs = num_dofs_per_proc[mpi_rank];

    const int first_local_dof = std::accumulate(num_dofs_per_proc.begin(), num_dofs_per_proc.begin() + mpi_rank, 0);
    const int last_local_dof = first_local_dof + n_local_dofs;

    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Box<NDIM> side_patch_box[NDIM];
        for (int axis = 0; axis < NDIM; ++axis)
        {
            side_patch_box[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
        }

        Pointer<SideData<NDIM, int>> u_dof_data = patch->getPatchData(u_dof_index_idx);
        Pointer<CellData<NDIM, int>> p_dof_data = patch->getPatchData(p_dof_index_idx);
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
                    is_field[0].insert(dof_idx);
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
                is_field[1].insert(dof_idx);
            }
        }
    }

    return;
} // constructPatchLevelFields

void
StaggeredStokesPETScMatUtilities::constructProlongationOp(Mat& mat,
                                                          const std::string& u_op_type,
                                                          const std::string& p_op_type,
                                                          int u_dof_index_idx,
                                                          int p_dof_index_idx,
                                                          const std::vector<int>& num_fine_dofs_per_proc,
                                                          const std::vector<int>& num_coarse_dofs_per_proc,
                                                          Pointer<PatchLevel<NDIM>> fine_patch_level,
                                                          Pointer<PatchLevel<NDIM>> coarse_patch_level,
                                                          const AO& coarse_level_ao,
                                                          const int u_coarse_ao_offset,
                                                          const int p_coarse_ao_offset)
{
    int ierr;
    Mat p_prolong_mat = nullptr;
    PETScMatUtilities::constructProlongationOp(mat,
                                               u_op_type,
                                               u_dof_index_idx,
                                               num_fine_dofs_per_proc,
                                               num_coarse_dofs_per_proc,
                                               fine_patch_level,
                                               coarse_patch_level,
                                               coarse_level_ao,
                                               u_coarse_ao_offset);

    PETScMatUtilities::constructProlongationOp(p_prolong_mat,
                                               p_op_type,
                                               p_dof_index_idx,
                                               num_fine_dofs_per_proc,
                                               num_coarse_dofs_per_proc,
                                               fine_patch_level,
                                               coarse_patch_level,
                                               coarse_level_ao,
                                               p_coarse_ao_offset);

    // P{u,p} = (P_u + P_p){u,p}
    ierr = MatAXPY(mat, 1.0, p_prolong_mat, DIFFERENT_NONZERO_PATTERN);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&p_prolong_mat);
    IBTK_CHKERRQ(ierr);

} // constructPatchLevelProlongationOp

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
