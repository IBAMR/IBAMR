// Filename: PoissonUtilities.cpp
// Created on 31 Mar 2012 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
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

#include <functional>
#include <limits>
#include <map>
#include <ostream>
#include <stddef.h>
#include <vector>

#include "ArrayData.h"
#include "ArrayDataBasicOps.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "EdgeData.h"
#include "Index.h"
#include "IntVector.h"
#include "NodeData.h"
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
inline Box<NDIM>
compute_tangential_extension(const Box<NDIM>& box, const int data_axis)
{
    Box<NDIM> extended_box = box;
    extended_box.upper()(data_axis) += 1;
    return extended_box;
}

#if (NDIM == 2)
inline double
compute_mu_avg(const Index<NDIM>& i, const NodeData<NDIM, double>& mu_data)
{
    Box<NDIM> node_box(i, i);
    const int n_nodes = pow(2, NDIM);

    double avg_mu = 0.0;
    int total_nodes = 0;
    for (NodeIterator<NDIM> n(node_box); n; n++, total_nodes++)
    {
        avg_mu += mu_data(n(), /*depth*/ 0);
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(total_nodes == n_nodes);
#endif
    return avg_mu / n_nodes;
} // compute_mu_avg

inline double
compute_mu_harmonic_avg(const Index<NDIM>& i, const NodeData<NDIM, double>& mu_data)
{
    Box<NDIM> node_box(i, i);
    const int n_nodes = pow(2, NDIM);

    double avg_mu = 0.0;
    int total_nodes = 0;
    for (NodeIterator<NDIM> n(node_box); n; n++, total_nodes++)
    {
        avg_mu += 1.0 / mu_data(n(), /*depth*/ 0);
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(total_nodes == n_nodes);
#endif
    return n_nodes / avg_mu;
} // compute_mu_harmonic_avg
#endif

#if (NDIM == 3)
inline double
compute_mu_avg(const Index<NDIM>& i, const EdgeData<NDIM, double>& mu_data)
{
    Box<NDIM> edge_box(i, i);
    const int n_edges = 12;

    double avg_mu = 0.0;
    int total_edges = 0;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        for(EdgeIterator<NDIM> e(edge_box, axis); e; e++, total_edges++)
        {
            avg_mu += mu_data(e(), /*depth*/ 0);
        }
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(total_edges == n_edges);
#endif
    return avg_mu / n_edges;
} // compute_mu_avg

inline double
compute_mu_harmonic_avg(const Index<NDIM>& i, const EdgeData<NDIM, double>& mu_data)
{
    Box<NDIM> edge_box(i, i);
    const int n_edges = 12;

    double avg_mu = 0.0;
    int total_edges = 0;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (EdgeIterator<NDIM> e(edge_box, axis); e; e++, total_edges++)
        {
            avg_mu += 1.0 / mu_data(e(), /*depth*/ 0);
        }
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(total_edges == n_edges);
#endif
    return n_edges / avg_mu;
} // compute_mu_harmonic_avg

inline double
get_mu_edge(const Index<NDIM>& i, const int perp, const Pointer<EdgeData<NDIM, double> > mu_data)
{
    const ArrayData<NDIM, double>& mu_array_data = mu_data->getArrayData(perp);
    return mu_array_data(i, /*depth*/ 0);
}
#endif

inline Index<NDIM>
get_shift(int dir, int shift)
{
    SAMRAI::hier::Index<NDIM> iv(0);
    iv(dir) = shift;
    return iv;
} // get_shift
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
    std::map<Index<NDIM>, int, IndexFortranOrder> stencil_map;
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
    std::map<Index<NDIM>, int, IndexFortranOrder> stencil_map;
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
PoissonUtilities::computeVCSCViscousOpMatrixCoefficients(
    SAMRAI::pdat::SideData<NDIM, double>& matrix_coefficients,
    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    const std::vector<std::map<Index<NDIM>, int, IndexFortranOrder> >& stencil_map_vec,
    const SAMRAI::solv::PoissonSpecifications& poisson_spec,
    double alpha,
    double beta,
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
    double data_time,
    VCInterpType mu_interp_type)
{

#if !defined(NDEBUG)
    TBOX_ASSERT(bc_coefs.size() == NDIM);
#endif
    const int stencil_sz = (2 * NDIM + 1) + 4 * (NDIM - 1);

#if !defined(NDEBUG)
    TBOX_ASSERT(matrix_coefficients.getDepth() == stencil_sz);
#endif
    matrix_coefficients.fillAll(0.0);

    const bool C_is_varying = poisson_spec.cIsVariable();
    Pointer<SideData<NDIM, double> > C_data = NULL;
    if (C_is_varying) C_data = patch->getPatchData(poisson_spec.getCPatchDataId());

#if (NDIM == 2)
    Pointer<NodeData<NDIM, double> > mu_data = patch->getPatchData(poisson_spec.getDPatchDataId());
#elif (NDIM == 3)
    Pointer<EdgeData<NDIM, double> > mu_data = patch->getPatchData(poisson_spec.getDPatchDataId());
#endif

#if !defined(NDEBUG)
    if (C_is_varying) TBOX_ASSERT(!C_data.isNull());
    TBOX_ASSERT(!mu_data.isNull());
#endif

#if (NDIM == 2)
    const ArrayData<NDIM, double>& mu_array_data = mu_data->getArrayData();
#endif

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
typedef std::map<Index<NDIM>, int, IndexFortranOrder> StencilMapType;
#if (NDIM == 2)
    StencilMapType stencil_map = stencil_map_vec[0];
#endif
    static const Index<NDIM> ORIGIN(0);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
#if (NDIM == 3)
        StencilMapType stencil_map = stencil_map_vec[axis];
#endif
        for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
        {
            const Index<NDIM>& cc = b();
            const SideIndex<NDIM> i(cc, axis, SideIndex<NDIM>::Lower);

            if (C_is_varying)
            {
                matrix_coefficients(i, stencil_map[ORIGIN]) = beta * (*C_data)(i, 0);
            }
            else
            {
                matrix_coefficients(i, stencil_map[ORIGIN]) =
                    poisson_spec.cIsZero() ? 0.0 : beta * poisson_spec.getCConstant();
            }
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d == axis)
                {
                    const Index<NDIM> shift_axis_plus = get_shift(axis, 1);
                    const Index<NDIM> shift_axis_minus = get_shift(axis, -1);

                    double mu_upper = std::numeric_limits<double>::quiet_NaN();
                    double mu_lower = std::numeric_limits<double>::quiet_NaN();
                    if (mu_interp_type == VC_AVERAGE_INTERP)
                    {
                        mu_upper = compute_mu_avg(cc, *mu_data);
                        mu_lower = compute_mu_avg(cc + shift_axis_minus, *mu_data);
                    }
                    else if (mu_interp_type == VC_HARMONIC_INTERP)
                    {
                        mu_upper = compute_mu_harmonic_avg(cc, *mu_data);
                        mu_lower = compute_mu_harmonic_avg(cc + shift_axis_minus, *mu_data);
                    }
                    else
                    {
                        TBOX_ERROR("this statement should not be reached");
                    }

                    const double coef_plus = (2.0 * alpha * mu_upper) / (dx[axis] * dx[axis]);
                    const double coef_minus = (2.0 * alpha * mu_lower) / (dx[axis] * dx[axis]);
                    matrix_coefficients(i, stencil_map[shift_axis_plus]) = coef_plus;
                    matrix_coefficients(i, stencil_map[shift_axis_minus]) = coef_minus;
                    matrix_coefficients(i, stencil_map[ORIGIN]) -= coef_plus + coef_minus;
                }
                else
                {
                    const Index<NDIM> shift_d_plus = get_shift(d, 1);
                    const Index<NDIM> shift_d_minus = get_shift(d, -1);
                    const Index<NDIM> shift_axis_plus = get_shift(axis, 1);
                    const Index<NDIM> shift_axis_minus = get_shift(axis, -1);

#if (NDIM == 2)
                    const double mu_upper = mu_array_data(cc + shift_d_plus, 0);
                    const double mu_lower = mu_array_data(cc, 0);
#elif (NDIM == 3)   
                    // Get edge data aligned with perp dir. (perpendicular to d and axis) and shifted in the d dir.
                    const int perp = 2*(d + axis) % 3; // 2 if {0,1}, 1 if {0,2} and 0 if {1,2} 
                    const double mu_upper = get_mu_edge(cc + shift_d_plus, perp, mu_data);
                    const double mu_lower = get_mu_edge(cc, perp, mu_data);
#endif

                    matrix_coefficients(i, stencil_map[shift_d_plus]) = (alpha * mu_upper) / (dx[d] * dx[d]);
                    matrix_coefficients(i, stencil_map[shift_d_minus]) = (alpha * mu_lower) / (dx[d] * dx[d]);
                    matrix_coefficients(i, stencil_map[ORIGIN]) -= matrix_coefficients(i, stencil_map[shift_d_plus]) +
                                                                   matrix_coefficients(i, stencil_map[shift_d_minus]);

                    matrix_coefficients(i, stencil_map[shift_d_plus + shift_axis_plus]) =
                        (alpha * mu_upper) / (dx[axis] * dx[d]);
                    matrix_coefficients(i, stencil_map[shift_d_plus + shift_axis_minus]) =
                        -(alpha * mu_upper) / (dx[axis] * dx[d]);
                    matrix_coefficients(i, stencil_map[shift_d_minus + shift_axis_plus]) =
                        -(alpha * mu_lower) / (dx[axis] * dx[d]);
                    matrix_coefficients(i, stencil_map[shift_d_minus + shift_axis_minus]) =
                        (alpha * mu_lower) / (dx[axis] * dx[d]);
                }
            }
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
#if (NDIM == 3)
        StencilMapType stencil_map = stencil_map_vec[axis];
#endif
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
                    Index<NDIM> shift = get_shift(bdry_normal_axis, -1);
                    matrix_coefficients(i_s, stencil_map[ORIGIN]) +=
                        matrix_coefficients(i_s, stencil_map[shift]) * (-(a * h - 2.0 * b) / (a * h + 2.0 * b));
                    matrix_coefficients(i_s, stencil_map[shift]) = 0.0;
                }
                else
                {
                    Index<NDIM> shift = get_shift(bdry_normal_axis, 1);
                    matrix_coefficients(i_s, stencil_map[ORIGIN]) +=
                        matrix_coefficients(i_s, stencil_map[shift]) * (-(a * h - 2.0 * b) / (a * h + 2.0 * b));
                    matrix_coefficients(i_s, stencil_map[shift]) = 0.0;
                }
            }
        }
    }

    // Modify matrix coefficients to account for physical boundary
    // conditions for other components of velocity along boundaries which ARE
    // aligned with the data axis.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
#if (NDIM == 3)
        StencilMapType stencil_map = stencil_map_vec[axis];
#endif
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
            const Box<NDIM> side_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            for (unsigned int comp = 0; comp < NDIM; ++comp)
            {
                if (comp == axis) continue;

                const Box<NDIM> bc_coef_box = compute_tangential_extension(side_box, comp);

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
                shifted_patch_x_lower[comp] -= 0.5 * dx[comp];
                shifted_patch_x_upper[comp] -= 0.5 * dx[comp];
                patch->setPatchGeometry(new CartesianPatchGeometry<NDIM>(ratio_to_level_zero,
                                                                         touches_regular_bdry,
                                                                         touches_periodic_bdry,
                                                                         dx,
                                                                         shifted_patch_x_lower.data(),
                                                                         shifted_patch_x_upper.data()));

                // Set the boundary condition coefficients.
                static const bool homogeneous_bc = true;
                ExtendedRobinBcCoefStrategy* extended_bc_coef =
                    dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[comp]);
                if (extended_bc_coef)
                {
                    extended_bc_coef->clearTargetPatchDataIndex();
                    extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                }
                bc_coefs[comp]->setBcCoefs(
                    acoef_data, bcoef_data, gcoef_data, NULL, *patch, trimmed_bdry_box, data_time);
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
                for (Box<NDIM>::Iterator bc(side_box); bc; bc++)
                {
                    const Index<NDIM>& i = bc();
                    const Index<NDIM> i_upper = i + get_shift(comp, 1);
                    const double& a_lower = (*acoef_data)(i, 0);
                    const double& b_lower = (*bcoef_data)(i, 0);
                    const double& a_upper = (*acoef_data)(i_upper, 0);
                    const double& b_upper = (*bcoef_data)(i_upper, 0);
                    const double& h = dx[bdry_normal_axis];

                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);

                    if (is_lower)
                    {
                        Index<NDIM> shift_outer_lower = get_shift(bdry_normal_axis, -1) + get_shift(comp, -1);
                        Index<NDIM> shift_inner_lower = get_shift(bdry_normal_axis, 1) + get_shift(comp, -1);
                        matrix_coefficients(i_s, stencil_map[shift_inner_lower]) +=
                            matrix_coefficients(i_s, stencil_map[shift_outer_lower]) *
                            (-(a_lower * h - 2.0 * b_lower) / (a_lower * h + 2.0 * b_lower));
                        matrix_coefficients(i_s, stencil_map[shift_outer_lower]) = 0.0;

                        Index<NDIM> shift_outer_upper = get_shift(bdry_normal_axis, -1) + get_shift(comp, 1);
                        Index<NDIM> shift_inner_upper = get_shift(bdry_normal_axis, 1) + get_shift(comp, 1);
                        matrix_coefficients(i_s, stencil_map[shift_inner_upper]) +=
                            matrix_coefficients(i_s, stencil_map[shift_outer_upper]) *
                            (-(a_upper * h - 2.0 * b_upper) / (a_upper * h + 2.0 * b_upper));
                        matrix_coefficients(i_s, stencil_map[shift_outer_upper]) = 0.0;
                    }
                    else
                    {
                        Index<NDIM> shift_outer_lower = get_shift(bdry_normal_axis, 1) + get_shift(comp, -1);
                        Index<NDIM> shift_inner_lower = get_shift(bdry_normal_axis, -1) + get_shift(comp, -1);
                        matrix_coefficients(i_s, stencil_map[shift_inner_lower]) +=
                            matrix_coefficients(i_s, stencil_map[shift_outer_lower]) *
                            (-(a_lower * h - 2.0 * b_lower) / (a_lower * h + 2.0 * b_lower));
                        matrix_coefficients(i_s, stencil_map[shift_outer_lower]) = 0.0;

                        Index<NDIM> shift_outer_upper = get_shift(bdry_normal_axis, 1) + get_shift(comp, 1);
                        Index<NDIM> shift_inner_upper = get_shift(bdry_normal_axis, -1) + get_shift(comp, 1);
                        matrix_coefficients(i_s, stencil_map[shift_inner_upper]) +=
                            matrix_coefficients(i_s, stencil_map[shift_outer_upper]) *
                            (-(a_upper * h - 2.0 * b_upper) / (a_upper * h + 2.0 * b_upper));
                        matrix_coefficients(i_s, stencil_map[shift_outer_upper]) = 0.0;
                    }
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
#if (NDIM == 3)
        StencilMapType stencil_map = stencil_map_vec[axis];
#endif
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
                    matrix_coefficients(i_s, stencil_map[ORIGIN]) = 1.0;
                }
                else
                {
#if !defined(NDEBUG)
                    TBOX_ASSERT(!MathUtilities<double>::equalEps(b, 0.0));
#endif
                    if (is_lower)
                    {
                        const Index<NDIM> shift_outer = get_shift(bdry_normal_axis, -1);
                        const Index<NDIM> shift_inner = get_shift(bdry_normal_axis, 1);
                        matrix_coefficients(i_s, stencil_map[ORIGIN]) -=
                            matrix_coefficients(i_s, stencil_map[shift_outer]) * 2 * dx[bdry_normal_axis] * a / b;
                        matrix_coefficients(i_s, stencil_map[shift_inner]) +=
                            matrix_coefficients(i_s, stencil_map[shift_outer]);
                        matrix_coefficients(i_s, stencil_map[shift_outer]) = 0.0;
                    }
                    else
                    {
                        const Index<NDIM> shift_outer = get_shift(bdry_normal_axis, 1);
                        const Index<NDIM> shift_inner = get_shift(bdry_normal_axis, -1);
                        matrix_coefficients(i_s, stencil_map[ORIGIN]) -=
                            matrix_coefficients(i_s, stencil_map[shift_outer]) * 2 * dx[bdry_normal_axis] * a / b;
                        matrix_coefficients(i_s, stencil_map[shift_inner]) +=
                            matrix_coefficients(i_s, stencil_map[shift_outer]);
                        matrix_coefficients(i_s, stencil_map[shift_outer]) = 0.0;
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
} // adjustRHSAtPhysicalBoundary

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
} // adjustRHSAtPhysicalBoundary

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
} // adjustRHSAtPhysicalBoundary

void
PoissonUtilities::adjustVCSCViscousOpRHSAtPhysicalBoundary(SideData<NDIM, double>& rhs_data,
                                                           Pointer<Patch<NDIM> > patch,
                                                           const PoissonSpecifications& poisson_spec,
                                                           double alpha,
                                                           const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                                           double data_time,
                                                           bool homogeneous_bc,
                                                           VCInterpType mu_interp_type)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(static_cast<int>(bc_coefs.size()) == NDIM);
#endif

#if (NDIM == 2)
    Pointer<NodeData<NDIM, double> > mu_data = patch->getPatchData(poisson_spec.getDPatchDataId());
#elif (NDIM == 3)
    Pointer<EdgeData<NDIM, double> > mu_data = patch->getPatchData(poisson_spec.getDPatchDataId());
#endif

#if !defined(NDEBUG)
    TBOX_ASSERT(!mu_data.isNull());
#endif

#if (NDIM == 2)
    const ArrayData<NDIM, double>& mu_array_data = mu_data->getArrayData();
#endif

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

                const Index<NDIM> shift_bdry_normal = get_shift(bdry_normal_axis, 1);
#if (NDIM == 2)
                const double mu_upper = mu_array_data(i_intr + shift_bdry_normal, 0);
                const double mu_lower = mu_array_data(i_intr, 0);
#elif (NDIM == 3)
                const int perp = 2*(bdry_normal_axis + axis) % 3; // 2 if {0,1}, 1 if {0,2} and 0 if {1,2} 
                const double mu_upper = get_mu_edge(i_intr + shift_bdry_normal, perp, mu_data);
                const double mu_lower = get_mu_edge(i_intr, perp, mu_data);
#endif

                const double D = is_lower ? mu_lower : mu_upper;
                rhs_data(i_s) += alpha * (D / h) * (-2.0 * g) / (2.0 * b + h * a);
            }
        }
    }

    // Modify the rhs entries to account for inhomogeneous boundary conditions
    // for other components of velocity along boundaries which ARE aligned with
    // the data axis.
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
            const Box<NDIM> side_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            for (unsigned int comp = 0; comp < NDIM; ++comp)
            {
                if (comp == axis) continue;

                const Box<NDIM> bc_coef_box = compute_tangential_extension(side_box, comp);

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
                shifted_patch_x_lower[comp] -= 0.5 * dx[comp];
                shifted_patch_x_upper[comp] -= 0.5 * dx[comp];
                patch->setPatchGeometry(new CartesianPatchGeometry<NDIM>(ratio_to_level_zero,
                                                                         touches_regular_bdry,
                                                                         touches_periodic_bdry,
                                                                         dx,
                                                                         shifted_patch_x_lower.data(),
                                                                         shifted_patch_x_upper.data()));

                // Set the boundary condition coefficients.
                ExtendedRobinBcCoefStrategy* extended_bc_coef =
                    dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[comp]);
                if (extended_bc_coef)
                {
                    extended_bc_coef->clearTargetPatchDataIndex();
                    extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                }
                bc_coefs[comp]->setBcCoefs(
                    acoef_data, bcoef_data, gcoef_data, NULL, *patch, trimmed_bdry_box, data_time);
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
                //     a*u_b + b*u_n = g
                //
                // then, with u_i = 0,
                //
                //     u_o = 2*h*g/(2*b + a*h)
                for (Box<NDIM>::Iterator bc(side_box); bc; bc++)
                {
                    const Index<NDIM>& i = bc();
                    const Index<NDIM>& i_upper = i + get_shift(comp, 1);
                    const double& a_lower = (*acoef_data)(i, 0);
                    const double& b_lower = (*bcoef_data)(i, 0);
                    const double& g_lower = (*gcoef_data)(i, 0);
                    const double& a_upper = (*acoef_data)(i_upper, 0);
                    const double& b_upper = (*bcoef_data)(i_upper, 0);
                    const double& g_upper = (*gcoef_data)(i_upper, 0);
                    const double& h = dx[bdry_normal_axis];
                    const double& hd = dx[comp];

                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);

                    const Index<NDIM> shift_d = get_shift(comp, 1);
#if (NDIM == 2)
                    const double mu_upper = mu_array_data(i + shift_d, 0);
                    const double mu_lower = mu_array_data(i, 0);
#elif (NDIM == 3)
                    const int perp = 2*(comp + axis) % 3; // 2 if {0,1}, 1 if {0,2} and 0 if {1,2} 
                    const double mu_upper = get_mu_edge(i + shift_d, perp, mu_data);
                    const double mu_lower = get_mu_edge(i, perp, mu_data);
#endif
                    if (is_lower)
                    {
                        rhs_data(i_s) -= alpha * (mu_lower / hd) * (2.0 * g_lower) / (2.0 * b_lower + a_lower * h);
                        rhs_data(i_s) += alpha * (mu_upper / hd) * (2.0 * g_upper) / (2.0 * b_upper + a_upper * h);
                    }
                    else
                    {
                        rhs_data(i_s) += alpha * (mu_lower / hd) * (2.0 * g_lower) / (2.0 * b_lower + a_lower * h);
                        rhs_data(i_s) -= alpha * (mu_upper / hd) * (2.0 * g_upper) / (2.0 * b_upper + a_upper * h);
                    }
                }
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

                const Index<NDIM> shift_axis_minus = get_shift(bdry_normal_axis, -1);
                double mu_upper = std::numeric_limits<double>::quiet_NaN();
                double mu_lower = std::numeric_limits<double>::quiet_NaN();
                if (mu_interp_type == VC_AVERAGE_INTERP)
                {
                    mu_upper = compute_mu_avg(i, *mu_data);
                    mu_lower = compute_mu_avg(i + shift_axis_minus, *mu_data);
                }
                else if (mu_interp_type == VC_HARMONIC_INTERP)
                {
                    mu_upper = compute_mu_harmonic_avg(i, *mu_data);
                    mu_lower = compute_mu_harmonic_avg(i + shift_axis_minus, *mu_data);
                }
                else
                {
                    TBOX_ERROR("this statement should not be reached");
                }
                const double D = is_lower ? mu_lower : mu_upper;

                if (b != 0.0)
                {
#if !defined(NDEBUG)
                    TBOX_ASSERT(!MathUtilities<double>::equalEps(b, 0.0));
#endif
                    rhs_data(i_s_bdry) += (2.0 * alpha) * (D / h) * (-2.0 * g) / b;
                }
            }
        }
    }
    return;
} // adjustVCSCViscousOpRHSAtPhysicalBoundary

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

void
PoissonUtilities::adjustVCSCViscousOpRHSAtCoarseFineBoundary(SideData<NDIM, double>& rhs_data,
                                                             const SideData<NDIM, double>& sol_data,
                                                             Pointer<Patch<NDIM> > patch,
                                                             const PoissonSpecifications& poisson_spec,
                                                             double alpha,
                                                             const Array<BoundaryBox<NDIM> >& type1_cf_bdry,
                                                             VCInterpType mu_interp_type)
{
#if (NDIM == 2)
    Pointer<NodeData<NDIM, double> > mu_data = patch->getPatchData(poisson_spec.getDPatchDataId());
#elif (NDIM == 3)
    Pointer<EdgeData<NDIM, double> > mu_data = patch->getPatchData(poisson_spec.getDPatchDataId());
#endif

#if !defined(NDEBUG)
    TBOX_ASSERT(!mu_data.isNull());
#endif

#if (NDIM == 2)
    const ArrayData<NDIM, double>& mu_array_data = mu_data->getArrayData();
#endif

    // Modify the rhs entries to account for coarse-fine interface boundary conditions.
    const int n_cf_bdry_boxes = type1_cf_bdry.size();
    const IntVector<NDIM> ghost_width_to_fill(1);

    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    // Modify the rhs entries to account for inhomogeneous Dirichelt boundary
    // conditions along boundaries which ARE NOT aligned with the data axis.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (int n = 0; n < n_cf_bdry_boxes; ++n)
        {
            const BoundaryBox<NDIM>& bdry_box = type1_cf_bdry[n];
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

            for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
            {
                const Index<NDIM>& i = bc();
                const double& h = dx[bdry_normal_axis];
                Index<NDIM> i_intr = i, i_bdry = i;
                if (is_lower)
                {
                    i_intr(bdry_normal_axis) += 0;
                    i_bdry(bdry_normal_axis) -= 1;
                }
                else
                {
                    i_intr(bdry_normal_axis) -= 1;
                    i_bdry(bdry_normal_axis) += 0;
                }
                const SideIndex<NDIM> i_s_intr(i_intr, axis, SideIndex<NDIM>::Lower);
                const SideIndex<NDIM> i_s_bdry(i_bdry, axis, SideIndex<NDIM>::Lower);

                const Index<NDIM> shift_bdry_normal = get_shift(bdry_normal_axis, 1);
#if (NDIM == 2)
                const double mu_upper = mu_array_data(i_intr + shift_bdry_normal, 0);
                const double mu_lower = mu_array_data(i_intr, 0);
#elif (NDIM == 3)
                const int perp = 2*(bdry_normal_axis + axis) % 3; // 2 if {0,1}, 1 if {0,2} and 0 if {1,2} 
                const double mu_upper = get_mu_edge(i_intr + shift_bdry_normal, perp, mu_data);
                const double mu_lower = get_mu_edge(i_intr, perp, mu_data);
#endif

                const double D = is_lower ? mu_lower : mu_upper;
                rhs_data(i_s_intr) -= alpha * (D / h) * sol_data(i_s_bdry) / h;
            }
        }
    }

    // Modify the rhs entries to account for inhomogeneous boundary conditions
    // for other components of velocity along boundaries which ARE aligned with
    // the data axis.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (int n = 0; n < n_cf_bdry_boxes; ++n)
        {
            const BoundaryBox<NDIM>& bdry_box = type1_cf_bdry[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = location_index % 2 == 0;

            if (bdry_normal_axis != axis) continue;

            const Box<NDIM> bc_fill_box =
                pgeom->getBoundaryFillBox(bdry_box, patch_box, /* ghost_width_to_fill */ IntVector<NDIM>(1));
            const BoundaryBox<NDIM> trimmed_bdry_box =
                PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
            const Box<NDIM> side_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            for (unsigned int comp = 0; comp < NDIM; ++comp)
            {
                if (comp == axis) continue;

                for (Box<NDIM>::Iterator bc(side_box); bc; bc++)
                {
                    const Index<NDIM>& i = bc();
                    const double& h = dx[bdry_normal_axis];
                    const double& hd = dx[comp];

                    const Index<NDIM> shift_d = get_shift(comp, 1);
#if (NDIM == 2)
                    const double mu_upper = mu_array_data(i + shift_d, 0);
                    const double mu_lower = mu_array_data(i, 0);
#elif (NDIM == 3)
                    const int perp = 2*(comp + axis) % 3; // 2 if {0,1}, 1 if {0,2} and 0 if {1,2} 
                    const double mu_upper = get_mu_edge(i + shift_d, perp, mu_data);
                    const double mu_lower = get_mu_edge(i, perp, mu_data);
#endif

                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                    if (is_lower)
                    {
                        const Index<NDIM> shift_axis_minus = get_shift(axis, -1);
                        const SideIndex<NDIM> sw(i + shift_axis_minus, comp, SideIndex<NDIM>::Lower);
                        const SideIndex<NDIM> nw(i + shift_axis_minus, comp, SideIndex<NDIM>::Upper);

                        rhs_data(i_s) -= alpha * (mu_lower / (h * hd)) * sol_data(sw);
                        rhs_data(i_s) += alpha * (mu_upper / (h * hd)) * sol_data(nw);
                    }
                    else
                    {
                        const Index<NDIM> shift_axis_plus = get_shift(axis, 1);
                        const SideIndex<NDIM> se(i, comp, SideIndex<NDIM>::Lower);
                        const SideIndex<NDIM> ne(i, comp, SideIndex<NDIM>::Upper);

                        rhs_data(i_s) += alpha * (mu_lower / (h * hd)) * sol_data(se);
                        rhs_data(i_s) -= alpha * (mu_upper / (h * hd)) * sol_data(ne);
                    }
                }
            }
        }
    }

    // Modify the rhs entries to account for inhomogeneous boundary conditions
    // along boundaries which ARE aligned with the data axis.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (int n = 0; n < n_cf_bdry_boxes; ++n)
        {
            const BoundaryBox<NDIM>& bdry_box = type1_cf_bdry[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = location_index % 2 == 0;

            if (bdry_normal_axis != axis) continue;

            const Box<NDIM> bc_fill_box =
                pgeom->getBoundaryFillBox(bdry_box, patch_box, /* ghost_width_to_fill */ IntVector<NDIM>(1));
            const BoundaryBox<NDIM> trimmed_bdry_box =
                PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
            const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
            {
                const Index<NDIM>& i = bc();
                const double& h = dx[bdry_normal_axis];

                const Index<NDIM> shift_axis_minus = get_shift(bdry_normal_axis, -1);
                const Index<NDIM> shift_axis_plus = get_shift(bdry_normal_axis, 1);

                const SideIndex<NDIM> i_s(i, bdry_normal_axis, SideIndex<NDIM>::Lower);
                const SideIndex<NDIM> i_s_bdry(
                    i + (is_lower ? shift_axis_minus : shift_axis_plus), bdry_normal_axis, SideIndex<NDIM>::Lower);

                double mu_upper = std::numeric_limits<double>::quiet_NaN();
                double mu_lower = std::numeric_limits<double>::quiet_NaN();
                if (mu_interp_type == VC_AVERAGE_INTERP)
                {
                    mu_upper = compute_mu_avg(i, *mu_data);
                    mu_lower = compute_mu_avg(i + shift_axis_minus, *mu_data);
                }
                else if (mu_interp_type == VC_HARMONIC_INTERP)
                {
                    mu_upper = compute_mu_harmonic_avg(i, *mu_data);
                    mu_lower = compute_mu_harmonic_avg(i + shift_axis_minus, *mu_data);
                }
                else
                {
                    TBOX_ERROR("this statement should not be reached");
                }
                const double D = is_lower ? mu_lower : mu_upper;

                rhs_data(i_s) -= (2.0 * alpha) * (D / h) * sol_data(i_s_bdry) / h;
            }
        }
    }
    return;
} // adjustVCSCViscousOpRHSAtCoarseFineBoundary

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
