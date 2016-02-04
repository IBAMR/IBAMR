// Filename: PoissonUtilities.cpp
// Created on 31 Mar 2012 by Boyce Griffith
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
#include <functional>
#include <map>
#include <ostream>
#include <vector>

#include "ArrayData.h"
#include "ArrayDataBasicOps.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "Index.h"
#include "IntVector.h"
#include "OutersideData.h"
#include "Patch.h"
#include "PatchData.h"
#include "PatchGeometry.h"
#include "PoissonSpecifications.h"
#include "RobinBcCoefStrategy.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "Variable.h"
#include "boost/array.hpp"
#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "ibtk/PoissonUtilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Array.h"
#include "tbox/MathUtilities.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
struct IndexComp : std::binary_function<Index<NDIM>, Index<NDIM>, bool>
{
    inline bool operator()(const Index<NDIM>& lhs, const Index<NDIM>& rhs) const
    {
        return ((lhs(0) < rhs(0))
#if (NDIM > 1)
                ||
                (lhs(0) == rhs(0) && lhs(1) < rhs(1))
#if (NDIM > 2)
                ||
                (lhs(0) == rhs(0) && lhs(1) == rhs(1) && lhs(2) < rhs(2))
#endif
#endif
                    );
    }
};

inline Box<NDIM>
compute_tangential_extension(const Box<NDIM>& box, const int data_axis)
{
    Box<NDIM> extended_box = box;
    extended_box.upper()(data_axis) += 1;
    return extended_box;
}
}

void
PoissonUtilities::computeMatrixCoefficients(CellData<NDIM, double>& matrix_coefficients,
                                            Pointer<Patch<NDIM> > patch,
                                            const std::vector<Index<NDIM> >& stencil,
                                            const PoissonSpecifications& poisson_spec,
                                            RobinBcCoefStrategy<NDIM>* bc_coef,
                                            double data_time)
{
    std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs(1, bc_coef);
    computeMatrixCoefficients(matrix_coefficients, patch, stencil, poisson_spec, bc_coefs, data_time);
    return;
}

void
PoissonUtilities::computeMatrixCoefficients(CellData<NDIM, double>& matrix_coefficients,
                                            Pointer<Patch<NDIM> > patch,
                                            const std::vector<Index<NDIM> >& stencil,
                                            const PoissonSpecifications& poisson_spec,
                                            const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                            double data_time)
{
    const int stencil_sz = static_cast<int>(stencil.size());
#if !defined(NDEBUG)
    TBOX_ASSERT(stencil_sz == 2 * NDIM + 1);
#endif
    std::map<Index<NDIM>, int, IndexComp> stencil_map;
    for (int k = 0; k < stencil_sz; ++k)
    {
        stencil_map[stencil[k]] = k;
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(stencil_map.find(Index<NDIM>(0)) != stencil_map.end());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        Index<NDIM> ilower(0), iupper(0);
        ilower(axis) = -1;
        iupper(axis) = +1;
        TBOX_ASSERT(stencil_map.find(ilower) != stencil_map.end());
        TBOX_ASSERT(stencil_map.find(iupper) != stencil_map.end());
    }
#endif
    const int stencil_index_diag = stencil_map[Index<NDIM>(0)];
    boost::array<int, NDIM> stencil_index_lower, stencil_index_upper;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        Index<NDIM> ilower(0), iupper(0);
        ilower(axis) = -1;
        iupper(axis) = +1;
        stencil_index_lower[axis] = stencil_map[ilower];
        stencil_index_upper[axis] = stencil_map[iupper];
    }
    const int depth = static_cast<int>(bc_coefs.size());
#if !defined(NDEBUG)
    TBOX_ASSERT(matrix_coefficients.getDepth() == depth * stencil_sz);
#endif
    const Box<NDIM>& patch_box = patch->getBox();
    CellData<NDIM, double> diagonal(patch_box, depth, IntVector<NDIM>(0));
    SideData<NDIM, double> off_diagonal(patch_box, depth, IntVector<NDIM>(0));

    ArrayDataBasicOps<NDIM, double> array_ops;
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    // Compute all off-diagonal matrix coefficients for all cell sides,
    // including those that touch the physical boundary; however, do not yet
    // take physical boundary conditions into account.  Boundary conditions are
    // handled subsequently.
    if (!poisson_spec.dIsConstant())
    {
        off_diagonal.copy(*patch->getPatchData(poisson_spec.getDPatchDataId()));
    }
    else
    {
        off_diagonal.fillAll(poisson_spec.getDConstant());
    }

    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
        array_ops.scale(
            off_diagonal.getArrayData(axis), 1.0 / (dx[axis] * dx[axis]), off_diagonal.getArrayData(axis), side_box);
    }

    // Compute all diagonal matrix coefficients for all cells, including those
    // that are adjacent to the physical boundary; however, do not yet take
    // physical boundary conditions into account.  Boundary conditions are
    // handled subsequently.
    if (!poisson_spec.cIsZero() && !poisson_spec.cIsConstant())
    {
        diagonal.copy(*patch->getPatchData(poisson_spec.getCPatchDataId()));
    }
    else
    {
        if (poisson_spec.cIsZero())
            diagonal.fillAll(0.0);
        else
            diagonal.fillAll(poisson_spec.getCConstant());
    }

    for (int d = 0; d < depth; ++d)
    {
        for (Box<NDIM>::Iterator b(patch_box); b; b++)
        {
            const Index<NDIM>& i = b();
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const SideIndex<NDIM> ilower(i, axis, SideIndex<NDIM>::Lower);
                diagonal(i, d) -= off_diagonal(ilower, d);
                const SideIndex<NDIM> iupper(i, axis, SideIndex<NDIM>::Upper);
                diagonal(i, d) -= off_diagonal(iupper, d);
            }
        }
    }

    // Modify matrix coefficients to account for physical boundary conditions.
    const Array<BoundaryBox<NDIM> > physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    for (int n = 0; n < n_physical_codim1_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
        const BoundaryBox<NDIM> trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
        const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

        Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

        for (int d = 0; d < depth; ++d)
        {
            static const bool homogeneous_bc = true;
            ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[d]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[d]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, NULL, *patch, trimmed_bdry_box, data_time);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = (location_index % 2) == 0;
            const bool is_upper = (location_index % 2) != 0;

            // Modify the diagonal and off-diagonal entries to account for
            // homogeneous boundary conditions.
            //
            // Here, we follow the same linear extrapolation approach
            // implemented in class CartesianRobinBcHelper.  Namely, with u_i
            // denoting the interior cell, u_o denoting the ghost cell, and u_b
            // and u_n denoting the value and normal derivative of u at the
            // boundary,
            //
            //     u_b = (u_i + u_o)/2   and   u_n = (u_o - u_i)/h
            //
            // Now, if
            //
            //     a*u_b + b*u_n = 0
            //
            // then
            //
            //     u_o = -((a*h - 2*b)/(a*h + 2*b))*u_i
            for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
            {
                const Index<NDIM>& i_s_bdry = bc();
                const double& a = (*acoef_data)(i_s_bdry, 0);
                const double& b = (*bcoef_data)(i_s_bdry, 0);
                const double& h = dx[bdry_normal_axis];

                // i_s_bdry: side index located on physical boundary
                //
                // i_c_intr: cell index located adjacent to physical boundary
                // in the patch interior
                Index<NDIM> i_c_intr = i_s_bdry;
                if (is_upper)
                {
                    i_c_intr(bdry_normal_axis) -= 1;
                }

                if (is_lower)
                {
                    const SideIndex<NDIM> ilower(i_c_intr, bdry_normal_axis, SideIndex<NDIM>::Lower);
                    diagonal(i_c_intr, d) += off_diagonal(ilower, d) * (-(a * h - 2.0 * b) / (a * h + 2.0 * b));
                    off_diagonal(ilower, d) = 0.0;
                }

                if (is_upper)
                {
                    const SideIndex<NDIM> iupper(i_c_intr, bdry_normal_axis, SideIndex<NDIM>::Upper);
                    diagonal(i_c_intr, d) += off_diagonal(iupper, d) * (-(a * h - 2.0 * b) / (a * h + 2.0 * b));
                    off_diagonal(iupper, d) = 0.0;
                }
            }
        }
    }

    // Setup the matrix coefficients.
    for (int d = 0; d < depth; ++d)
    {
        const unsigned int offset = static_cast<unsigned int>(d * stencil_sz);
        for (Box<NDIM>::Iterator b(patch_box); b; b++)
        {
            const Index<NDIM>& i = b();
            matrix_coefficients(i, offset + stencil_index_diag) = diagonal(i, d);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const SideIndex<NDIM> ilower(i, axis, SideIndex<NDIM>::Lower);
                matrix_coefficients(i, offset + stencil_index_lower[axis]) = off_diagonal(ilower, d);
                const SideIndex<NDIM> iupper(i, axis, SideIndex<NDIM>::Upper);
                matrix_coefficients(i, offset + stencil_index_upper[axis]) = off_diagonal(iupper, d);
            }
        }
    }
    return;
}

void
PoissonUtilities::computeMatrixCoefficients(SideData<NDIM, double>& matrix_coefficients,
                                            Pointer<Patch<NDIM> > patch,
                                            const std::vector<Index<NDIM> >& stencil,
                                            const PoissonSpecifications& poisson_spec,
                                            const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                            double data_time)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_coefs.size() == NDIM);
#endif
    const int stencil_sz = static_cast<int>(stencil.size());
#if !defined(NDEBUG)
    TBOX_ASSERT(stencil_sz == 2 * NDIM + 1);
#endif
    std::map<Index<NDIM>, int, IndexComp> stencil_map;
    for (int k = 0; k < stencil_sz; ++k)
    {
        stencil_map[stencil[k]] = k;
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(stencil_map.find(Index<NDIM>(0)) != stencil_map.end());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        Index<NDIM> ilower(0), iupper(0);
        ilower(axis) = -1;
        iupper(axis) = +1;
        TBOX_ASSERT(stencil_map.find(ilower) != stencil_map.end());
        TBOX_ASSERT(stencil_map.find(iupper) != stencil_map.end());
    }
#endif
    const int stencil_index_diag = stencil_map[Index<NDIM>(0)];
    boost::array<int, NDIM> stencil_index_lower, stencil_index_upper;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        Index<NDIM> ilower(0), iupper(0);
        ilower(axis) = -1;
        iupper(axis) = +1;
        stencil_index_lower[axis] = stencil_map[ilower];
        stencil_index_upper[axis] = stencil_map[iupper];
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(matrix_coefficients.getDepth() == stencil_sz);
#endif
    if (!(poisson_spec.cIsZero() || poisson_spec.cIsConstant()) || !poisson_spec.dIsConstant())
    {
        TBOX_ERROR(
            "PoissonUtilities::computeSCMatrixCoefficients() does not support non-constant "
            "coefficient problems\n");
    }
    const double C = (poisson_spec.cIsZero() ? 0.0 : poisson_spec.getCConstant());
    const double D = poisson_spec.getDConstant();

    const Box<NDIM>& patch_box = patch->getBox();
    const Array<BoundaryBox<NDIM> > physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
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

    // Compute all matrix coefficients, including those on the physical
    // boundary; however, do not yet take physical boundary conditions into
    // account.  Boundary conditions are handled subsequently.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        std::vector<double> mat_vals(stencil_sz, 0.0);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            const double dx_sq = dx[d] * dx[d];
            mat_vals[stencil_index_lower[d]] += D / dx_sq; // lower off-diagonal
            mat_vals[stencil_index_upper[d]] += D / dx_sq; // upper off-diagonal
            mat_vals[stencil_index_diag] -= D / dx_sq;     // diagonal
            mat_vals[stencil_index_diag] -= D / dx_sq;     // diagonal
        }
        mat_vals[stencil_index_diag] += C; // diagonal

        for (int stencil_index = 0; stencil_index < stencil_sz; ++stencil_index)
        {
            matrix_coefficients.fill(mat_vals[stencil_index], stencil_index);
        }
    }

    // Modify matrix coefficients to account for physical boundary conditions
    // along boundaries which ARE NOT aligned with the data axis.
    //
    // NOTE: It important to set these values first to avoid problems at corners
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

            if (bdry_normal_axis == axis) continue;

            const Box<NDIM> bc_fill_box =
                pgeom->getBoundaryFillBox(bdry_box, patch_box, /* ghost_width_to_fill */ IntVector<NDIM>(1));
            const BoundaryBox<NDIM> trimmed_bdry_box =
                PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
            const Box<NDIM> bc_coef_box = compute_tangential_extension(
                PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box), axis);

            Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

            // Temporarily reset the patch geometry object associated with the
            // patch so that boundary conditions are set at the correct spatial
            // locations.
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
            ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, NULL, *patch, trimmed_bdry_box, data_time);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

            // Restore the original patch geometry object.
            patch->setPatchGeometry(pgeom);

            // Modify the matrix coefficients to account for homogeneous
            // boundary conditions.
            //
            // Here, we follow the same linear extrapolation approach
            // implemented in class CartesianRobinBcHelper.  Namely, with u_i
            // denoting the interior cell, u_o denoting the ghost cell, and u_b
            // and u_n denoting the value and normal derivative of u at the
            // boundary,
            //
            //     u_b = (u_i + u_o)/2   and   u_n = (u_o - u_i)/h
            //
            // Now, if
            //
            //     a*u_b + b*u_n = 0
            //
            // then
            //
            //     u_o = -((a*h - 2*b)/(a*h + 2*b))*u_i
            for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
            {
                const Index<NDIM>& i = bc();
                const double& a = (*acoef_data)(i, 0);
                const double& b = (*bcoef_data)(i, 0);
                const double& h = dx[bdry_normal_axis];

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

                if (is_lower)
                {
                    matrix_coefficients(i_s, stencil_index_diag) +=
                        matrix_coefficients(i_s, stencil_index_lower[bdry_normal_axis]) *
                        (-(a * h - 2.0 * b) / (a * h + 2.0 * b));
                    matrix_coefficients(i_s, stencil_index_lower[bdry_normal_axis]) = 0.0;
                }
                else
                {
                    matrix_coefficients(i_s, stencil_index_diag) +=
                        matrix_coefficients(i_s, stencil_index_upper[bdry_normal_axis]) *
                        (-(a * h - 2.0 * b) / (a * h + 2.0 * b));
                    matrix_coefficients(i_s, stencil_index_upper[bdry_normal_axis]) = 0.0;
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
            Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

            // Set the boundary condition coefficients.
            static const bool homogeneous_bc = true;
            ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, NULL, *patch, trimmed_bdry_box, data_time);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

            // Modify the matrix coefficients to account for homogeneous
            // boundary conditions.
            //
            // With u_i denoting the interior cell, u_o denoting the ghost cell,
            // and u_b and u_n denoting the value and normal derivative of u at
            // the boundary,
            //
            //     u_n = (u_o - u_i)/(2*h)
            //
            // Now, if
            //
            //     a*u_b + b*u_n = 0
            //
            // then if b != 0,
            //
            //     u_o = u_i - 2*h*(a/b)*u_b
            //
            // and
            //
            //     -(D/h^2)*u_o = (D*2*(a/b)/h)*u_b - (D/h^2)*u_i
            //
            // If b == 0, then u_b = 0, which we enforce directly.
            for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
            {
                const Index<NDIM>& i = bc();
                const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                const double& a = (*acoef_data)(i, 0);
                const double& b = (*bcoef_data)(i, 0);
                if (b == 0.0)
                {
                    for (int k = 0; k < stencil_sz; ++k)
                    {
                        matrix_coefficients(i_s, k) = 0.0;
                    }
                    matrix_coefficients(i_s, stencil_index_diag) = 1.0;
                }
                else
                {
#if !defined(NDEBUG)
                    TBOX_ASSERT(!MathUtilities<double>::equalEps(b, 0.0));
#endif
                    if (is_lower)
                    {
                        matrix_coefficients(i_s, stencil_index_diag) -=
                            matrix_coefficients(i_s, stencil_index_lower[bdry_normal_axis]) * 2 * dx[bdry_normal_axis] *
                            a / b;
                        matrix_coefficients(i_s, stencil_index_upper[bdry_normal_axis]) +=
                            matrix_coefficients(i_s, stencil_index_lower[bdry_normal_axis]);
                        matrix_coefficients(i_s, stencil_index_lower[bdry_normal_axis]) = 0.0;
                    }
                    else
                    {
                        matrix_coefficients(i_s, stencil_index_diag) -=
                            matrix_coefficients(i_s, stencil_index_upper[bdry_normal_axis]) * 2 * dx[bdry_normal_axis] *
                            a / b;
                        matrix_coefficients(i_s, stencil_index_lower[bdry_normal_axis]) +=
                            matrix_coefficients(i_s, stencil_index_upper[bdry_normal_axis]);
                        matrix_coefficients(i_s, stencil_index_upper[bdry_normal_axis]) = 0.0;
                    }
                }
            }
        }
    }
    return;
}

void
PoissonUtilities::adjustRHSAtPhysicalBoundary(CellData<NDIM, double>& rhs_data,
                                              Pointer<Patch<NDIM> > patch,
                                              const PoissonSpecifications& poisson_spec,
                                              RobinBcCoefStrategy<NDIM>* bc_coef,
                                              double data_time,
                                              bool homogeneous_bc)
{
    std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs(1, bc_coef);
    adjustRHSAtPhysicalBoundary(rhs_data, patch, poisson_spec, bc_coefs, data_time, homogeneous_bc);
    return;
}

void
PoissonUtilities::adjustRHSAtPhysicalBoundary(CellData<NDIM, double>& rhs_data,
                                              Pointer<Patch<NDIM> > patch,
                                              const PoissonSpecifications& poisson_spec,
                                              const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                              double data_time,
                                              bool homogeneous_bc)
{
    const int depth = rhs_data.getDepth();
#if !defined(NDEBUG)
    TBOX_ASSERT(static_cast<int>(bc_coefs.size()) == depth);
#endif
    const Box<NDIM>& patch_box = patch->getBox();
    OutersideData<NDIM, double> D_data(patch_box, depth);
    if (!poisson_spec.dIsConstant())
    {
        D_data.copy(*patch->getPatchData(poisson_spec.getDPatchDataId()));
    }
    else
    {
        D_data.fillAll(poisson_spec.getDConstant());
    }
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    // Modify the rhs entries to account for inhomogeneous boundary conditions.
    const Array<BoundaryBox<NDIM> > codim1_boxes = PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
    const int n_bdry_boxes = codim1_boxes.size();
    for (int n = 0; n < n_bdry_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = codim1_boxes[n];
        const BoundaryBox<NDIM> trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
        const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

        Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
        Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

        for (int d = 0; d < depth; ++d)
        {
            ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[d]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[d]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, NULL, *patch, trimmed_bdry_box, data_time);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_upper = (location_index % 2) != 0;
            const int bdry_side = (is_upper ? 1 : 0);

            // Here, we follow the same linear extrapolation approach
            // implemented in class CartesianRobinBcHelper.  Namely, with u_i
            // denoting the interior cell, u_o denoting the ghost cell, and u_b
            // and u_n denoting the value and normal derivative of u at the
            // boundary,
            //
            //     u_b = (u_i + u_o)/2   and   u_n = (u_o - u_i)/h
            //
            // Now, if
            //
            //     a*u_b + b*u_n = g
            //
            // then, with u_i = 0,
            //
            //     u_o = 2*h*g/(2*b + a*h)
            //
            // so that the boundary flux is
            //
            //     (u_i - u_o)/h = -2*g/(2*b + h*a)
            //
            // In this loop, we modify the rhs entries appropriately.
            //
            // NOTE: i_s_bdry: side index located on physical boundary
            //       i_c_intr: cell index located adjacent to physical boundary
            //                 in the patch interior
            for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
            {
                const Index<NDIM>& i_s_bdry = bc();
                const double& a = (*acoef_data)(i_s_bdry, 0);
                const double& b = (*bcoef_data)(i_s_bdry, 0);
                const double& g = (*gcoef_data)(i_s_bdry, 0);
                const double& h = dx[bdry_normal_axis];
                Index<NDIM> i_c_intr = i_s_bdry;
                if (is_upper)
                {
                    i_c_intr(bdry_normal_axis) -= 1;
                }
                const double& D = D_data.getArrayData(bdry_normal_axis, bdry_side)(i_s_bdry, d);
                rhs_data(i_c_intr, d) += (D / h) * (-2.0 * g) / (2.0 * b + h * a);
            }
        }
    }
    return;
}

void
PoissonUtilities::adjustRHSAtPhysicalBoundary(SideData<NDIM, double>& rhs_data,
                                              Pointer<Patch<NDIM> > patch,
                                              const PoissonSpecifications& poisson_spec,
                                              const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                              double data_time,
                                              bool homogeneous_bc)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(static_cast<int>(bc_coefs.size()) == NDIM);
#endif
    if (!(poisson_spec.cIsZero() || poisson_spec.cIsConstant()) || !poisson_spec.dIsConstant())
    {
        TBOX_ERROR(
            "PoissonUtilities::adjustRHSAtPhysicalBoundary() does not support non-constant "
            "coefficient problems\n");
    }
    const Box<NDIM>& patch_box = patch->getBox();
    const double D = poisson_spec.getDConstant();
    const Array<BoundaryBox<NDIM> > physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
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

    // Modify the rhs entries to account for inhomogeneous boundary conditions
    // along boundaries which ARE NOT aligned with the data axis.
    //
    // NOTE: It important to set these values first to avoid problems at corners
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

            if (bdry_normal_axis == axis) continue;

            const Box<NDIM> bc_fill_box =
                pgeom->getBoundaryFillBox(bdry_box, patch_box, /* ghost_width_to_fill */ IntVector<NDIM>(1));
            const BoundaryBox<NDIM> trimmed_bdry_box =
                PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
            const Box<NDIM> bc_coef_box = compute_tangential_extension(
                PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box), axis);

            Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

            // Temporarily reset the patch geometry object associated with the
            // patch so that boundary conditions are set at the correct spatial
            // locations.
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
            ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, NULL, *patch, trimmed_bdry_box, data_time);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

            // Restore the original patch geometry object.
            patch->setPatchGeometry(pgeom);

            // Here, we follow the same linear extrapolation approach
            // implemented in class CartesianRobinBcHelper.  Namely, with u_i
            // denoting the interior cell, u_o denoting the ghost cell, and u_b
            // and u_n denoting the value and normal derivative of u at the
            // boundary,
            //
            //     u_b = (u_i + u_o)/2   and   u_n = (u_o - u_i)/h
            //
            // Now, if
            //
            //     a*u_b + b*u_n = g
            //
            // then, with u_i = 0,
            //
            //     u_o = 2*h*g/(2*b + a*h)
            //
            // so that the boundary flux is
            //
            //     (u_i - u_o)/h = -2*g/(2*b + h*a)
            //
            // In this loop, we modify the rhs entries appropriately.
            for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
            {
                const Index<NDIM>& i = bc();
                const double& a = (*acoef_data)(i, 0);
                const double& b = (*bcoef_data)(i, 0);
                const double& g = (*gcoef_data)(i, 0);
                const double& h = dx[bdry_normal_axis];
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
                rhs_data(i_s) += (D / h) * (-2.0 * g) / (2.0 * b + h * a);
            }
        }
    }

    // Modify the rhs entries to account for inhomogeneous boundary conditions
    // along boundaries which ARE aligned with the data axis.
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

            if (bdry_normal_axis != axis) continue;

            const Box<NDIM> bc_fill_box =
                pgeom->getBoundaryFillBox(bdry_box, patch_box, /* ghost_width_to_fill */ IntVector<NDIM>(1));
            const BoundaryBox<NDIM> trimmed_bdry_box =
                PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
            const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

            // Set the boundary condition coefficients.
            ExtendedRobinBcCoefStrategy* extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, NULL, *patch, trimmed_bdry_box, data_time);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

            // For the non-symmetric boundary treatment,
            //
            //     u_n = (u_o - u_i)/(2*h)
            //
            // Now, if
            //
            //     a*u_b + b*u_n = g
            //
            // then, with u_b = u_i = 0,
            //
            //     u_o = 2*h*g/b
            //
            // so that the boundary flux is
            //
            //     (u_b - u_o)/h = -2*g/b
            //
            // In this loop, we modify the rhs entries appropriately.
            //
            // NOTE: At Dirichlet boundaries, boundary values are provided by
            // the right-hand side vector.
            for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
            {
                const Index<NDIM>& i = bc();
                const double& b = (*bcoef_data)(i, 0);
                const double& g = (*gcoef_data)(i, 0);
                const double& h = dx[bdry_normal_axis];
                const SideIndex<NDIM> i_s_bdry(i, bdry_normal_axis, SideIndex<NDIM>::Lower);
                if (b != 0.0)
                {
#if !defined(NDEBUG)
                    TBOX_ASSERT(!MathUtilities<double>::equalEps(b, 0.0));
#endif
                    rhs_data(i_s_bdry) += (D / h) * (-2.0 * g) / b;
                }
            }
        }
    }
    return;
}

void
PoissonUtilities::adjustRHSAtCoarseFineBoundary(CellData<NDIM, double>& rhs_data,
                                                const CellData<NDIM, double>& sol_data,
                                                Pointer<Patch<NDIM> > patch,
                                                const PoissonSpecifications& poisson_spec,
                                                const Array<BoundaryBox<NDIM> >& type1_cf_bdry)
{
    const int depth = rhs_data.getDepth();
    TBOX_ASSERT(depth == sol_data.getDepth());
    const Box<NDIM>& patch_box = patch->getBox();
    OutersideData<NDIM, double> D_data(patch_box, depth);
    if (!poisson_spec.dIsConstant())
    {
        D_data.copy(*patch->getPatchData(poisson_spec.getDPatchDataId()));
    }
    else
    {
        D_data.fillAll(poisson_spec.getDConstant());
    }
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    // Modify the rhs entries to account for coarse-fine interface boundary conditions.
    const int n_bdry_boxes = type1_cf_bdry.size();
    const IntVector<NDIM> ghost_width_to_fill(1);
    for (int n = 0; n < n_bdry_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = type1_cf_bdry[n];
        const BoundaryBox<NDIM> trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
        const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(trimmed_bdry_box, patch_box, ghost_width_to_fill);
        const unsigned int location_index = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;
        const double h = dx[bdry_normal_axis];
        const int bdry_side = location_index % 2;
        const bool is_lower = bdry_side == 0;
        const bool is_upper = bdry_side == 1;
        for (int d = 0; d < depth; ++d)
        {
            for (Box<NDIM>::Iterator bc(bc_fill_box); bc; bc++)
            {
                Index<NDIM> i_s_bdry = bc();
                CellIndex<NDIM> i_c_intr(i_s_bdry);
                CellIndex<NDIM> i_c_bdry(i_s_bdry);
                if (is_lower)
                {
                    i_c_intr(bdry_normal_axis) += 1;
                    i_s_bdry(bdry_normal_axis) += 1;
                }
                if (is_upper) i_c_intr(bdry_normal_axis) -= 1;
                const double& D = D_data.getArrayData(bdry_normal_axis, bdry_side)(i_s_bdry, d);
                rhs_data(i_c_intr, d) -= (D / h) * sol_data(i_c_bdry, d) / h;
            }
        }
    }
    return;
} // adjustRHSAtCoarseFineBoundary

void
PoissonUtilities::adjustRHSAtCoarseFineBoundary(SideData<NDIM, double>& rhs_data,
                                                const SideData<NDIM, double>& sol_data,
                                                Pointer<Patch<NDIM> > patch,
                                                const PoissonSpecifications& poisson_spec,
                                                const Array<BoundaryBox<NDIM> >& type1_cf_bdry)
{
    const int depth = rhs_data.getDepth();
    TBOX_ASSERT(depth == sol_data.getDepth());
    if (!(poisson_spec.cIsZero() || poisson_spec.cIsConstant()) || !poisson_spec.dIsConstant())
    {
        TBOX_ERROR(
            "PoissonUtilities::adjustRHSAtCoarseFineBoundary() does not support non-constant "
            "coefficient problems\n");
    }
    const Box<NDIM>& patch_box = patch->getBox();
    const double D = poisson_spec.getDConstant();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    // Modify the rhs entries to account for coarse-fine interface boundary conditions.
    const int n_bdry_boxes = type1_cf_bdry.size();
    const IntVector<NDIM> ghost_width_to_fill(1);
    for (int n = 0; n < n_bdry_boxes; ++n)
    {
        const BoundaryBox<NDIM>& bdry_box = type1_cf_bdry[n];
        const BoundaryBox<NDIM> trimmed_bdry_box = PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
        const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(trimmed_bdry_box, patch_box, ghost_width_to_fill);
        const unsigned int location_index = bdry_box.getLocationIndex();
        const unsigned int bdry_normal_axis = location_index / 2;
        const double h = dx[bdry_normal_axis];
        const int bdry_side = location_index % 2;
        const bool is_lower = bdry_side == 0;
        const bool is_upper = bdry_side == 1;
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            const Box<NDIM> bc_fill_box_axis =
                bdry_normal_axis == axis ? bc_fill_box : compute_tangential_extension(bc_fill_box, axis);
            for (int d = 0; d < depth; ++d)
            {
                for (Box<NDIM>::Iterator bc(bc_fill_box); bc; bc++)
                {
                    SideIndex<NDIM> i_s_intr(bc(), axis, SideIndex<NDIM>::Lower);
                    SideIndex<NDIM> i_s_bdry(bc(), axis, SideIndex<NDIM>::Lower);
                    if (is_lower)
                    {
                        i_s_intr(bdry_normal_axis) += 1;
                    }
                    if (is_upper)
                    {
                        i_s_intr(bdry_normal_axis) -= axis == bdry_normal_axis ? 0 : 1;
                        i_s_bdry(axis) += axis == bdry_normal_axis ? 1 : 0;
                    }
                    rhs_data(i_s_intr, d) -= (D / h) * sol_data(i_s_bdry, d) / h;
                }
            }
        }
    }
    return;
} // adjustRHSAtCoarseFineBoundary

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
