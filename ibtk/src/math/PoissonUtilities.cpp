// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2023 by the IBAMR developers
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

#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "ibtk/PoissonUtilities.h"
#include "ibtk/ProblemSpecification.h"
#include "ibtk/ibtk_enums.h"

#include "ArrayData.h"
#include "ArrayDataBasicOps.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "EdgeData.h"
#include "EdgeIterator.h"
#include "IntVector.h"
#include "NodeData.h"
#include "NodeIterator.h"
#include "OutersideData.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "PoissonSpecifications.h"
#include "RobinBcCoefStrategy.h"
#include "SideData.h"
#include "SideIndex.h"
#include "Variable.h"
#include "tbox/Array.h"
#include "tbox/Pointer.h"

#include <array>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

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
compute_mu_avg(const hier::Index<NDIM>& i, const NodeData<NDIM, double>& mu_data)
{
    Box<NDIM> node_box(i, i);
    const int n_nodes = std::pow(2, NDIM);

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
compute_mu_harmonic_avg(const hier::Index<NDIM>& i, const NodeData<NDIM, double>& mu_data)
{
    Box<NDIM> node_box(i, i);
    const int n_nodes = std::pow(2, NDIM);

    double avg_mu = 0.0;
    int total_nodes = 0;
    for (NodeIterator<NDIM> n(node_box); n; n++, total_nodes++)
    {
        const double mu = mu_data(n(), /*depth*/ 0);
        if (IBTK::abs_equal_eps(mu, 0.0))
        {
            return 0.0;
        }
        avg_mu += 1.0 / mu;
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(total_nodes == n_nodes);
#endif
    return n_nodes / avg_mu;
} // compute_mu_harmonic_avg
#endif

#if (NDIM == 3)
inline double
compute_mu_avg(const hier::Index<NDIM>& i, const EdgeData<NDIM, double>& mu_data)
{
    Box<NDIM> edge_box(i, i);
    const int n_edges = 12;

    double avg_mu = 0.0;
    int total_edges = 0;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (EdgeIterator<NDIM> e(edge_box, axis); e; e++, total_edges++)
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
compute_mu_harmonic_avg(const hier::Index<NDIM>& i, const EdgeData<NDIM, double>& mu_data)
{
    Box<NDIM> edge_box(i, i);
    const int n_edges = 12;

    double avg_mu = 0.0;
    int total_edges = 0;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (EdgeIterator<NDIM> e(edge_box, axis); e; e++, total_edges++)
        {
            const double mu = mu_data(e(), /*depth*/ 0);
            if (IBTK::abs_equal_eps(mu, 0.0))
            {
                return 0.0;
            }
            avg_mu += 1.0 / mu;
        }
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(total_edges == n_edges);
#endif
    return n_edges / avg_mu;
} // compute_mu_harmonic_avg

inline double
get_mu_edge(const hier::Index<NDIM>& i, const int perp, const Pointer<EdgeData<NDIM, double> > mu_data)
{
    const ArrayData<NDIM, double>& mu_array_data = mu_data->getArrayData(perp);
    return mu_array_data(i, /*depth*/ 0);
}
#endif

inline hier::Index<NDIM>
get_shift(int dir, int shift)
{
    SAMRAI::hier::Index<NDIM> iv(0);
    iv(dir) = shift;
    return iv;
} // get_shift
} // namespace

void
PoissonUtilities::computeMatrixCoefficients(CellData<NDIM, double>& matrix_coefficients,
                                            Pointer<Patch<NDIM> > patch,
                                            const std::vector<hier::Index<NDIM> >& stencil,
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
                                            const std::vector<hier::Index<NDIM> >& stencil,
                                            const PoissonSpecifications& poisson_spec,
                                            const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                            double data_time)
{
    const int stencil_sz = static_cast<int>(stencil.size());
#if !defined(NDEBUG)
    TBOX_ASSERT(stencil_sz == 2 * NDIM + 1);
#endif
    std::map<hier::Index<NDIM>, int, IndexFortranOrder> stencil_map;
    for (int k = 0; k < stencil_sz; ++k)
    {
        stencil_map[stencil[k]] = k;
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(stencil_map.find(hier::Index<NDIM>(0)) != stencil_map.end());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        hier::Index<NDIM> ilower(0), iupper(0);
        ilower(axis) = -1;
        iupper(axis) = +1;
        TBOX_ASSERT(stencil_map.find(ilower) != stencil_map.end());
        TBOX_ASSERT(stencil_map.find(iupper) != stencil_map.end());
    }
#endif
    const int stencil_index_diag = stencil_map[hier::Index<NDIM>(0)];
    std::array<int, NDIM> stencil_index_lower, stencil_index_upper;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        hier::Index<NDIM> ilower(0), iupper(0);
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
            const hier::Index<NDIM>& i = b();
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
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[d]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[d]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
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
                const hier::Index<NDIM>& i_s_bdry = bc();
                const double& a = (*acoef_data)(i_s_bdry, 0);
                const double& b = (*bcoef_data)(i_s_bdry, 0);
                const double& h = dx[bdry_normal_axis];

                // i_s_bdry: side index located on physical boundary
                //
                // i_c_intr: cell index located adjacent to physical boundary
                // in the patch interior
                hier::Index<NDIM> i_c_intr = i_s_bdry;
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
        const auto offset = static_cast<unsigned int>(d * stencil_sz);
        for (Box<NDIM>::Iterator b(patch_box); b; b++)
        {
            const hier::Index<NDIM>& i = b();
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
                                            const std::vector<hier::Index<NDIM> >& stencil,
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
    std::map<hier::Index<NDIM>, int, IndexFortranOrder> stencil_map;
    for (int k = 0; k < stencil_sz; ++k)
    {
        stencil_map[stencil[k]] = k;
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(stencil_map.find(hier::Index<NDIM>(0)) != stencil_map.end());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        hier::Index<NDIM> ilower(0), iupper(0);
        ilower(axis) = -1;
        iupper(axis) = +1;
        TBOX_ASSERT(stencil_map.find(ilower) != stencil_map.end());
        TBOX_ASSERT(stencil_map.find(iupper) != stencil_map.end());
    }
#endif
    const int stencil_index_diag = stencil_map[hier::Index<NDIM>(0)];
    std::array<int, NDIM> stencil_index_lower, stencil_index_upper;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        hier::Index<NDIM> ilower(0), iupper(0);
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
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(
                acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
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
                const hier::Index<NDIM>& i = bc();
                const double& a = (*acoef_data)(i, 0);
                const double& b = (*bcoef_data)(i, 0);
                const double& h = dx[bdry_normal_axis];

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
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(
                acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
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
                const hier::Index<NDIM>& i = bc();
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
                    TBOX_ASSERT(!IBTK::abs_equal_eps(b, 0.0));
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
PoissonUtilities::computeVCCCViscousOpMatrixCoefficients(
    SAMRAI::pdat::CellData<NDIM, double>& matrix_coefficients,
    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    const std::vector<CCVCStencilMap>& stencil_map_vec,
    const SAMRAI::solv::PoissonSpecifications& poisson_spec,
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
    double data_time,
    VCInterpType mu_interp_type)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(static_cast<int>(bc_coefs.size()) == NDIM);
#endif

    using StencilMapType = CCVCStencilMap;
    const auto make_key = makeCCVCStencilKey;

    static const hier::Index<NDIM> ORIGIN(0);
    static const int stencil_sz = 1 + 2 * NDIM + 9 * (NDIM - 1);

#if !defined(NDEBUG)
    TBOX_ASSERT(matrix_coefficients.getDepth() == stencil_sz * NDIM);
#endif

    matrix_coefficients.fillAll(0.0);

    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    const bool C_is_variable = poisson_spec.cIsVariable();
    Pointer<CellData<NDIM, double> > C_data = nullptr;
    if (C_is_variable) C_data = patch->getPatchData(poisson_spec.getCPatchDataId());

    const bool D_is_variable = poisson_spec.dIsVariable();
    Pointer<CellData<NDIM, double> > mu_cc_data = nullptr;
    if (D_is_variable) mu_cc_data = patch->getPatchData(poisson_spec.getDPatchDataId());

#if !defined(NDEBUG)
    if (C_is_variable) TBOX_ASSERT(!C_data.isNull());
    if (D_is_variable) TBOX_ASSERT(!mu_cc_data.isNull());
#endif

    auto mu_cc = [&](const CellIndex<NDIM>& i) -> double
    { return D_is_variable ? (*mu_cc_data)(i, 0) : poisson_spec.getDConstant(); };

    auto interp_mu = [&](const CellIndex<NDIM>& i, int dir, int shift) -> double
    {
        const CellIndex<NDIM> j = i + get_shift(dir, shift);
        const double mu0 = mu_cc(i);
        const double mu1 = mu_cc(j);

        switch (mu_interp_type)
        {
        case VC_AVERAGE_INTERP:
            return 0.5 * (mu0 + mu1);

        case VC_HARMONIC_INTERP:
            return IBTK::abs_equal_eps(mu0 + mu1, 0.0) ? 0.0 : 2.0 * mu0 * mu1 / (mu0 + mu1);

        default:
            TBOX_ERROR("PoissonUtilities::computeCCVCViscousOpMatrixCoefficients():\n"
                       << "  unsupported mu_interp_type\n");
            return 0.0;
        }
    };

    auto robin_elim_factor = [&](double a, double b, double h) -> double
    { return -(a * h - 2.0 * b) / (a * h + 2.0 * b); };

    // ---------------------------------------------------------------------
    // Interior operator.
    // ---------------------------------------------------------------------
    for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
    {
        const CellIndex<NDIM>& i = b();

        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            const StencilMapType& sm = stencil_map_vec[axis];
            const int offset = axis * stencil_sz;

            if (C_is_variable)
            {
                matrix_coefficients(i, offset + sm.at(make_key(ORIGIN, axis))) = (*C_data)(i, 0);
            }
            else
            {
                matrix_coefficients(i, offset + sm.at(make_key(ORIGIN, axis))) =
                    poisson_spec.cIsZero() ? 0.0 : poisson_spec.getCConstant();
            }

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const hier::Index<NDIM> ed_p = get_shift(d, +1);
                const hier::Index<NDIM> ed_m = get_shift(d, -1);

                if (d == axis)
                {
                    const double mu_p = interp_mu(i, axis, +1);
                    const double mu_m = interp_mu(i, axis, -1);

                    const double cp = 2.0 * mu_p / (dx[axis] * dx[axis]);
                    const double cm = 2.0 * mu_m / (dx[axis] * dx[axis]);

                    matrix_coefficients(i, offset + sm.at(make_key(ed_p, axis))) += cp;
                    matrix_coefficients(i, offset + sm.at(make_key(ed_m, axis))) += cm;
                    matrix_coefficients(i, offset + sm.at(make_key(ORIGIN, axis))) -= cp + cm;
                }
                else
                {
                    // d/dx_d [ mu du_axis/dx_d ]
                    const double mu_p = interp_mu(i, d, +1);
                    const double mu_m = interp_mu(i, d, -1);

                    const double same_p = mu_p / (dx[d] * dx[d]);
                    const double same_m = mu_m / (dx[d] * dx[d]);

                    matrix_coefficients(i, offset + sm.at(make_key(ed_p, axis))) += same_p;
                    matrix_coefficients(i, offset + sm.at(make_key(ed_m, axis))) += same_m;
                    matrix_coefficients(i, offset + sm.at(make_key(ORIGIN, axis))) -= same_p + same_m;

                    // d/dx_d [ mu du_d/dx_axis ]
                    const hier::Index<NDIM> ea_p = get_shift(axis, +1);
                    const hier::Index<NDIM> ea_m = get_shift(axis, -1);

                    const double cp = mu_p / (4.0 * dx[d] * dx[axis]);
                    const double cm = mu_m / (4.0 * dx[d] * dx[axis]);

                    // +d face
                    matrix_coefficients(i, offset + sm.at(make_key(ea_p + ed_p, d))) += cp;
                    matrix_coefficients(i, offset + sm.at(make_key(ea_p, d))) += cp;
                    matrix_coefficients(i, offset + sm.at(make_key(ea_m + ed_p, d))) += -cp;
                    matrix_coefficients(i, offset + sm.at(make_key(ea_m, d))) += -cp;

                    // -d face
                    matrix_coefficients(i, offset + sm.at(make_key(ea_p, d))) += -cm;
                    matrix_coefficients(i, offset + sm.at(make_key(ea_p + ed_m, d))) += -cm;
                    matrix_coefficients(i, offset + sm.at(make_key(ea_m, d))) += cm;
                    matrix_coefficients(i, offset + sm.at(make_key(ea_m + ed_m, d))) += cm;
                }
            }
        }
    }

    const Array<BoundaryBox<NDIM> > physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();

    // ---------------------------------------------------------------------
    // Pass 1:
    // eliminate same-component ghost values on boundaries not aligned
    // with the row axis.
    // ---------------------------------------------------------------------
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const StencilMapType& sm = stencil_map_vec[axis];
        const int offset = axis * stencil_sz;

        for (int n = 0; n < n_physical_codim1_boxes; ++n)
        {
            const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = (location_index % 2 == 0);

            if (bdry_normal_axis == axis) continue;

            const BoundaryBox<NDIM> trimmed_bdry_box =
                PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
            const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

            static const bool homogeneous_bc = true;
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }

            bc_coefs[axis]->setBcCoefs(
                acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);

            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

            for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
            {
                const CellIndex<NDIM>& ib = bc();
                const double acoef = (*acoef_data)(ib, 0);
                const double bcoef = (*bcoef_data)(ib, 0);
                const double h = dx[bdry_normal_axis];

                CellIndex<NDIM> i_intr = ib;
                if (!is_lower) i_intr(bdry_normal_axis) -= 1;

                const hier::Index<NDIM> outer =
                    is_lower ? get_shift(bdry_normal_axis, -1) : get_shift(bdry_normal_axis, +1);

                const double r = robin_elim_factor(acoef, bcoef, h);

                matrix_coefficients(i_intr, offset + sm.at(make_key(ORIGIN, axis))) +=
                    matrix_coefficients(i_intr, offset + sm.at(make_key(outer, axis))) * r;
                matrix_coefficients(i_intr, offset + sm.at(make_key(outer, axis))) = 0.0;
            }
        }
    }

    // ---------------------------------------------------------------------
    // Pass 2:
    // eliminate cross-component ghost strip entries on aligned boundaries
    // using the correct interior partner across the boundary face:
    //
    // lower boundary:
    //   (-e_axis + shift, comp) -> (shift, comp)
    //
    // upper boundary:
    //   (+e_axis + shift, comp) -> (shift, comp)
    //
    // where shift in { +e_comp, 0, -e_comp }.
    // ---------------------------------------------------------------------
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const StencilMapType& sm = stencil_map_vec[axis];
        const int offset = axis * stencil_sz;

        for (int n = 0; n < n_physical_codim1_boxes; ++n)
        {
            const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = (location_index % 2 == 0);

            if (bdry_normal_axis != axis) continue;

            const BoundaryBox<NDIM> trimmed_bdry_box =
                PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
            const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            for (unsigned int comp = 0; comp < NDIM; ++comp)
            {
                if (comp == axis) continue;

                Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

                static const bool homogeneous_bc = true;
                auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[comp]);
                if (extended_bc_coef)
                {
                    extended_bc_coef->clearTargetPatchDataIndex();
                    extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                }

                bc_coefs[comp]->setBcCoefs(
                    acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);

                if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

                const hier::Index<NDIM> ea_p = get_shift(axis, +1);
                const hier::Index<NDIM> ea_m = get_shift(axis, -1);
                const hier::Index<NDIM> ec_p = get_shift(comp, +1);
                const hier::Index<NDIM> ec_m = get_shift(comp, -1);

                for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
                {
                    const CellIndex<NDIM>& i = bc();
                    const double acoef = (*acoef_data)(i, 0);
                    const double bcoef = (*bcoef_data)(i, 0);
                    const double h = dx[axis];
                    const double r = robin_elim_factor(acoef, bcoef, h);

                    CellIndex<NDIM> i_intr = i;
                    if (!is_lower) i_intr(axis) -= 1;

                    if (is_lower)
                    {
                        // (-e_axis + e_comp, comp) -> (+e_comp, comp)
                        matrix_coefficients(i_intr, offset + sm.at(make_key(ec_p, comp))) +=
                            matrix_coefficients(i_intr, offset + sm.at(make_key(ea_m + ec_p, comp))) * r;
                        matrix_coefficients(i_intr, offset + sm.at(make_key(ea_m + ec_p, comp))) = 0.0;

                        // (-e_axis, comp) -> (0, comp)
                        matrix_coefficients(i_intr, offset + sm.at(make_key(ORIGIN, comp))) +=
                            matrix_coefficients(i_intr, offset + sm.at(make_key(ea_m, comp))) * r;
                        matrix_coefficients(i_intr, offset + sm.at(make_key(ea_m, comp))) = 0.0;

                        // (-e_axis - e_comp, comp) -> (-e_comp, comp)
                        matrix_coefficients(i_intr, offset + sm.at(make_key(ec_m, comp))) +=
                            matrix_coefficients(i_intr, offset + sm.at(make_key(ea_m + ec_m, comp))) * r;
                        matrix_coefficients(i_intr, offset + sm.at(make_key(ea_m + ec_m, comp))) = 0.0;
                    }
                    else
                    {
                        // (+e_axis + e_comp, comp) -> (+e_comp, comp)
                        matrix_coefficients(i_intr, offset + sm.at(make_key(ec_p, comp))) +=
                            matrix_coefficients(i_intr, offset + sm.at(make_key(ea_p + ec_p, comp))) * r;
                        matrix_coefficients(i_intr, offset + sm.at(make_key(ea_p + ec_p, comp))) = 0.0;

                        // (+e_axis, comp) -> (0, comp)
                        matrix_coefficients(i_intr, offset + sm.at(make_key(ORIGIN, comp))) +=
                            matrix_coefficients(i_intr, offset + sm.at(make_key(ea_p, comp))) * r;
                        matrix_coefficients(i_intr, offset + sm.at(make_key(ea_p, comp))) = 0.0;

                        // (+e_axis - e_comp, comp) -> (-e_comp, comp)
                        matrix_coefficients(i_intr, offset + sm.at(make_key(ec_m, comp))) +=
                            matrix_coefficients(i_intr, offset + sm.at(make_key(ea_p + ec_m, comp))) * r;
                        matrix_coefficients(i_intr, offset + sm.at(make_key(ea_p + ec_m, comp))) = 0.0;
                    }
                }
            }
        }
    }

    // ---------------------------------------------------------------------
    // Pass 3:
    // enforce BC for row component on aligned boundaries.
    // ---------------------------------------------------------------------
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const StencilMapType& sm = stencil_map_vec[axis];
        const int offset = axis * stencil_sz;

        for (int n = 0; n < n_physical_codim1_boxes; ++n)
        {
            const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = (location_index % 2 == 0);

            if (bdry_normal_axis != axis) continue;

            const BoundaryBox<NDIM> trimmed_bdry_box =
                PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
            const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

            static const bool homogeneous_bc = true;
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }

            bc_coefs[axis]->setBcCoefs(
                acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);

            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

            for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
            {
                const CellIndex<NDIM>& ib = bc();
                const double acoef = (*acoef_data)(ib, 0);
                const double bcoef = (*bcoef_data)(ib, 0);
                const double h = dx[axis];

                CellIndex<NDIM> i_intr = ib;
                if (!is_lower) i_intr(axis) -= 1;

                const hier::Index<NDIM> outer = is_lower ? get_shift(axis, -1) : get_shift(axis, +1);

                const double r = robin_elim_factor(acoef, bcoef, h);

                matrix_coefficients(i_intr, offset + sm.at(make_key(ORIGIN, axis))) +=
                    matrix_coefficients(i_intr, offset + sm.at(make_key(outer, axis))) * r;
                matrix_coefficients(i_intr, offset + sm.at(make_key(outer, axis))) = 0.0;
            }
        }
    }

    return;
} // computeCCVCViscousOpMatrixCoefficients

void
PoissonUtilities::computeVCCCViscousDilatationalOpMatrixCoefficients(
    CellData<NDIM, double>& matrix_coefficients,
    Pointer<Patch<NDIM> > patch,
    const std::vector<CCVCStencilMap>& stencil_map_vec,
    const ProblemSpecification* problem_spec,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
    double data_time,
    VCInterpType mu_interp_type)
{
    const auto& vc_op_spec = static_cast<const VCViscousDilatationalOpSpec&>(*problem_spec);

#if !defined(NDEBUG)
    TBOX_ASSERT(static_cast<int>(bc_coefs.size()) == NDIM);
#endif

    using StencilMapType = CCVCStencilMap;
    const auto make_key = makeCCVCStencilKey;

    static const hier::Index<NDIM> ORIGIN(0);
    static const int stencil_sz = 1 + 2 * NDIM + 9 * (NDIM - 1);

#if !defined(NDEBUG)
    TBOX_ASSERT(matrix_coefficients.getDepth() == stencil_sz * NDIM);
#endif

    matrix_coefficients.fillAll(0.0);

    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    const bool C_is_const = vc_op_spec.d_C_is_const;
    const double C_const = vc_op_spec.d_C_const;
    const int C_idx = vc_op_spec.d_C_idx;
    Pointer<CellData<NDIM, double> > C_data = nullptr;
    if (!C_is_const)
    {
        C_data = patch->getPatchData(C_idx);
    }

    Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(vc_op_spec.d_D_idx);
    Pointer<CellData<NDIM, double> > lambda_data = patch->getPatchData(vc_op_spec.d_L_idx);

#if !defined(NDEBUG)
    TBOX_ASSERT(!mu_data.isNull());
    TBOX_ASSERT(!lambda_data.isNull());
    if (!C_is_const) TBOX_ASSERT(!C_data.isNull());
#endif

    auto mu_interp_cc = [&](const CellIndex<NDIM>& i, int dir, int shift) -> double
    {
        const CellIndex<NDIM> j = i + get_shift(dir, shift);
        const double q0 = (*mu_data)(i, 0);
        const double q1 = (*mu_data)(j, 0);

        switch (mu_interp_type)
        {
        case VC_AVERAGE_INTERP:
            return 0.5 * (q0 + q1);
        case VC_HARMONIC_INTERP:
            return IBTK::abs_equal_eps(q0 + q1, 0.0) ? 0.0 : 2.0 * q0 * q1 / (q0 + q1);
        default:
            TBOX_ERROR("PoissonUtilities::computeCCVCViscousDilatationalOpMatrixCoefficients():\n"
                       << "  unsupported interpolation type for shear viscosity.\n");
            return 0.0;
        }
    };

    auto lambda_interp_cc = [&](const CellIndex<NDIM>& i, int dir, int shift) -> double
    {
        const CellIndex<NDIM> j = i + get_shift(dir, shift);
        const double q0 = (*lambda_data)(i, 0);
        const double q1 = (*lambda_data)(j, 0);

        return 0.5 * (q0 + q1);
    };

    auto robin_elim_factor = [&](double a, double b, double h) -> double
    { return -(a * h - 2.0 * b) / (a * h + 2.0 * b); };

    // Interior operator:
    //   C u + div(mu (grad u + grad u^T)) + grad(lambda div u)
    for (Box<NDIM>::Iterator b(CellGeometry<NDIM>::toCellBox(patch_box)); b; b++)
    {
        const CellIndex<NDIM>& i = b();

        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            const StencilMapType& sm = stencil_map_vec[axis];
            const int offset = axis * stencil_sz;

            matrix_coefficients(i, offset + sm.at(make_key(ORIGIN, axis))) = C_is_const ? C_const : (*C_data)(i, 0);

            // Contributions from d/dx_axis of the normal stress:
            // (2 mu + lambda) du_axis/dx_axis + lambda sum_{comp != axis} du_comp/dx_comp
            {
                const hier::Index<NDIM> ea_p = get_shift(axis, +1);
                const hier::Index<NDIM> ea_m = get_shift(axis, -1);

                const double mu_p = mu_interp_cc(i, axis, +1);
                const double mu_m = mu_interp_cc(i, axis, -1);
                const double la_p = lambda_interp_cc(i, axis, +1);
                const double la_m = lambda_interp_cc(i, axis, -1);

                const double cp = (2.0 * mu_p + la_p) / (dx[axis] * dx[axis]);
                const double cm = (2.0 * mu_m + la_m) / (dx[axis] * dx[axis]);

                matrix_coefficients(i, offset + sm.at(make_key(ea_p, axis))) += cp;
                matrix_coefficients(i, offset + sm.at(make_key(ea_m, axis))) += cm;
                matrix_coefficients(i, offset + sm.at(make_key(ORIGIN, axis))) -= cp + cm;

                for (unsigned int comp = 0; comp < NDIM; ++comp)
                {
                    if (comp == axis) continue;

                    const hier::Index<NDIM> ea_p = get_shift(axis, +1);
                    const hier::Index<NDIM> ea_m = get_shift(axis, -1);
                    const hier::Index<NDIM> ec_p = get_shift(comp, +1);
                    const hier::Index<NDIM> ec_m = get_shift(comp, -1);

                    const double lcp = la_p / (4.0 * dx[axis] * dx[comp]);
                    const double lcm = la_m / (4.0 * dx[axis] * dx[comp]);

                    // +axis face contribution:
                    // +(lambda_+ / dx_axis) * (du_comp/dx_comp)|_{+axis face}
                    matrix_coefficients(i, offset + sm.at(make_key(ea_p + ec_p, comp))) += lcp;
                    matrix_coefficients(i, offset + sm.at(make_key(ec_p, comp))) += lcp;
                    matrix_coefficients(i, offset + sm.at(make_key(ea_p + ec_m, comp))) += -lcp;
                    matrix_coefficients(i, offset + sm.at(make_key(ec_m, comp))) += -lcp;

                    // -axis face contribution enters with opposite sign:
                    // -(lambda_- / dx_axis) * (du_comp/dx_comp)|_{-axis face}
                    matrix_coefficients(i, offset + sm.at(make_key(ec_p, comp))) += -lcm;
                    matrix_coefficients(i, offset + sm.at(make_key(ea_m + ec_p, comp))) += -lcm;
                    matrix_coefficients(i, offset + sm.at(make_key(ec_m, comp))) += lcm;
                    matrix_coefficients(i, offset + sm.at(make_key(ea_m + ec_m, comp))) += lcm;
                }
            }

            // Contributions from d/dx_d of transverse stresses, d != axis:
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d == axis) continue;

                const hier::Index<NDIM> ed_p = get_shift(d, +1);
                const hier::Index<NDIM> ed_m = get_shift(d, -1);
                const hier::Index<NDIM> ea_p = get_shift(axis, +1);
                const hier::Index<NDIM> ea_m = get_shift(axis, -1);

                const double mu_p = mu_interp_cc(i, d, +1);
                const double mu_m = mu_interp_cc(i, d, -1);

                // d/dx_d [ mu du_axis/dx_d ]
                const double same_p = mu_p / (dx[d] * dx[d]);
                const double same_m = mu_m / (dx[d] * dx[d]);

                matrix_coefficients(i, offset + sm.at(make_key(ed_p, axis))) += same_p;
                matrix_coefficients(i, offset + sm.at(make_key(ed_m, axis))) += same_m;
                matrix_coefficients(i, offset + sm.at(make_key(ORIGIN, axis))) -= same_p + same_m;

                // d/dx_d [ mu du_d/dx_axis ]
                const double cp = mu_p / (4.0 * dx[d] * dx[axis]);
                const double cm = mu_m / (4.0 * dx[d] * dx[axis]);

                // +d face
                matrix_coefficients(i, offset + sm.at(make_key(ea_p + ed_p, d))) += cp;
                matrix_coefficients(i, offset + sm.at(make_key(ea_p, d))) += cp;
                matrix_coefficients(i, offset + sm.at(make_key(ea_m + ed_p, d))) += -cp;
                matrix_coefficients(i, offset + sm.at(make_key(ea_m, d))) += -cp;

                // -d face
                matrix_coefficients(i, offset + sm.at(make_key(ea_p, d))) += -cm;
                matrix_coefficients(i, offset + sm.at(make_key(ea_p + ed_m, d))) += -cm;
                matrix_coefficients(i, offset + sm.at(make_key(ea_m, d))) += cm;
                matrix_coefficients(i, offset + sm.at(make_key(ea_m + ed_m, d))) += cm;
            }
        }
    }

    const Array<BoundaryBox<NDIM> > physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();

    // Pass 1: boundaries not aligned with row axis.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const StencilMapType& sm = stencil_map_vec[axis];
        const int offset = axis * stencil_sz;

        for (int n = 0; n < n_physical_codim1_boxes; ++n)
        {
            const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = (location_index % 2 == 0);

            if (bdry_normal_axis == axis) continue;

            const BoundaryBox<NDIM> trimmed_bdry_box =
                PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
            const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

            static const bool homogeneous_bc = true;
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }

            bc_coefs[axis]->setBcCoefs(
                acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

            for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
            {
                const CellIndex<NDIM>& i = bc();

                CellIndex<NDIM> i_intr = i;
                if (!is_lower) i_intr(bdry_normal_axis) -= 1;

                const double acoef = (*acoef_data)(i, 0);
                const double bcoef = (*bcoef_data)(i, 0);
                const double r = robin_elim_factor(acoef, bcoef, dx[bdry_normal_axis]);

                const hier::Index<NDIM> outer =
                    is_lower ? get_shift(bdry_normal_axis, -1) : get_shift(bdry_normal_axis, +1);

                matrix_coefficients(i_intr, offset + sm.at(make_key(ORIGIN, axis))) +=
                    r * matrix_coefficients(i_intr, offset + sm.at(make_key(outer, axis)));
                matrix_coefficients(i_intr, offset + sm.at(make_key(outer, axis))) = 0.0;
            }
        }
    }

    // Pass 2: eliminate cross-component ghost strip entries on aligned boundaries.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const StencilMapType& sm = stencil_map_vec[axis];
        const int offset = axis * stencil_sz;

        for (int n = 0; n < n_physical_codim1_boxes; ++n)
        {
            const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = (location_index % 2 == 0);

            if (bdry_normal_axis != axis) continue;

            const BoundaryBox<NDIM> trimmed_bdry_box =
                PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
            const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            for (unsigned int comp = 0; comp < NDIM; ++comp)
            {
                if (comp == axis) continue;

                Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

                static const bool homogeneous_bc = true;
                auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[comp]);
                if (extended_bc_coef)
                {
                    extended_bc_coef->clearTargetPatchDataIndex();
                    extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                }

                bc_coefs[comp]->setBcCoefs(
                    acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
                if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

                const hier::Index<NDIM> ea_p = get_shift(axis, +1);
                const hier::Index<NDIM> ea_m = get_shift(axis, -1);
                const hier::Index<NDIM> ec_p = get_shift(comp, +1);
                const hier::Index<NDIM> ec_m = get_shift(comp, -1);

                for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
                {
                    const CellIndex<NDIM>& i = bc();

                    CellIndex<NDIM> i_intr = i;
                    if (!is_lower) i_intr(axis) -= 1;

                    const double acoef = (*acoef_data)(i, 0);
                    const double bcoef = (*bcoef_data)(i, 0);
                    const double r = robin_elim_factor(acoef, bcoef, dx[axis]);

                    if (is_lower)
                    {
                        matrix_coefficients(i_intr, offset + sm.at(make_key(ec_p, comp))) +=
                            r * matrix_coefficients(i_intr, offset + sm.at(make_key(ea_m + ec_p, comp)));
                        matrix_coefficients(i_intr, offset + sm.at(make_key(ea_m + ec_p, comp))) = 0.0;

                        matrix_coefficients(i_intr, offset + sm.at(make_key(ORIGIN, comp))) +=
                            r * matrix_coefficients(i_intr, offset + sm.at(make_key(ea_m, comp)));
                        matrix_coefficients(i_intr, offset + sm.at(make_key(ea_m, comp))) = 0.0;

                        matrix_coefficients(i_intr, offset + sm.at(make_key(ec_m, comp))) +=
                            r * matrix_coefficients(i_intr, offset + sm.at(make_key(ea_m + ec_m, comp)));
                        matrix_coefficients(i_intr, offset + sm.at(make_key(ea_m + ec_m, comp))) = 0.0;
                    }
                    else
                    {
                        matrix_coefficients(i_intr, offset + sm.at(make_key(ec_p, comp))) +=
                            r * matrix_coefficients(i_intr, offset + sm.at(make_key(ea_p + ec_p, comp)));
                        matrix_coefficients(i_intr, offset + sm.at(make_key(ea_p + ec_p, comp))) = 0.0;

                        matrix_coefficients(i_intr, offset + sm.at(make_key(ORIGIN, comp))) +=
                            r * matrix_coefficients(i_intr, offset + sm.at(make_key(ea_p, comp)));
                        matrix_coefficients(i_intr, offset + sm.at(make_key(ea_p, comp))) = 0.0;

                        matrix_coefficients(i_intr, offset + sm.at(make_key(ec_m, comp))) +=
                            r * matrix_coefficients(i_intr, offset + sm.at(make_key(ea_p + ec_m, comp)));
                        matrix_coefficients(i_intr, offset + sm.at(make_key(ea_p + ec_m, comp))) = 0.0;
                    }
                }
            }
        }
    }

    // Pass 3: eliminate row-component normal ghost into the center.
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const StencilMapType& sm = stencil_map_vec[axis];
        const int offset = axis * stencil_sz;

        for (int n = 0; n < n_physical_codim1_boxes; ++n)
        {
            const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = (location_index % 2 == 0);

            if (bdry_normal_axis != axis) continue;

            const BoundaryBox<NDIM> trimmed_bdry_box =
                PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
            const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

            static const bool homogeneous_bc = true;
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }

            bc_coefs[axis]->setBcCoefs(
                acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

            for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
            {
                const CellIndex<NDIM>& i = bc();

                CellIndex<NDIM> i_intr = i;
                if (!is_lower) i_intr(axis) -= 1;

                const double acoef = (*acoef_data)(i, 0);
                const double bcoef = (*bcoef_data)(i, 0);
                const double h = dx[axis];
                const double r = robin_elim_factor(acoef, bcoef, h);

                const hier::Index<NDIM> outer = is_lower ? get_shift(axis, -1) : get_shift(axis, +1);

                matrix_coefficients(i_intr, offset + sm.at(make_key(ORIGIN, axis))) +=
                    r * matrix_coefficients(i_intr, offset + sm.at(make_key(outer, axis)));
                matrix_coefficients(i_intr, offset + sm.at(make_key(outer, axis))) = 0.0;
            }
        }
    }
    return;
} // computeVCCCViscousDilatationalOpMatrixCoefficients

void
PoissonUtilities::computeVCSCViscousOpMatrixCoefficients(
    SAMRAI::pdat::SideData<NDIM, double>& matrix_coefficients,
    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    const std::vector<std::map<hier::Index<NDIM>, int, IndexFortranOrder> >& stencil_map_vec,
    const SAMRAI::solv::PoissonSpecifications& poisson_spec,
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
    double data_time,
    VCInterpType mu_interp_type,
    bool impose_physical_bcs_normal_comp)
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
    Pointer<SideData<NDIM, double> > C_data = nullptr;
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
    using StencilMapType = std::map<hier::Index<NDIM>, int, IndexFortranOrder>;
#if (NDIM == 2)
    StencilMapType stencil_map = stencil_map_vec[0];
#endif
    static const hier::Index<NDIM> ORIGIN(0);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
#if (NDIM == 3)
        StencilMapType stencil_map = stencil_map_vec[axis];
#endif
        for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
        {
            const hier::Index<NDIM>& cc = b();
            const SideIndex<NDIM> i(cc, axis, SideIndex<NDIM>::Lower);

            if (C_is_varying)
            {
                matrix_coefficients(i, stencil_map[ORIGIN]) = (*C_data)(i, 0);
            }
            else
            {
                matrix_coefficients(i, stencil_map[ORIGIN]) =
                    poisson_spec.cIsZero() ? 0.0 : poisson_spec.getCConstant();
            }
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d == axis)
                {
                    const hier::Index<NDIM> shift_axis_plus = get_shift(axis, 1);
                    const hier::Index<NDIM> shift_axis_minus = get_shift(axis, -1);

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

                    const double coef_plus = (2.0 * mu_upper) / (dx[axis] * dx[axis]);
                    const double coef_minus = (2.0 * mu_lower) / (dx[axis] * dx[axis]);
                    matrix_coefficients(i, stencil_map[shift_axis_plus]) = coef_plus;
                    matrix_coefficients(i, stencil_map[shift_axis_minus]) = coef_minus;
                    matrix_coefficients(i, stencil_map[ORIGIN]) -= coef_plus + coef_minus;
                }
                else
                {
                    const hier::Index<NDIM> shift_d_plus = get_shift(d, 1);
                    const hier::Index<NDIM> shift_d_minus = get_shift(d, -1);
                    const hier::Index<NDIM> shift_axis_plus = get_shift(axis, 1);
                    const hier::Index<NDIM> shift_axis_minus = get_shift(axis, -1);

#if (NDIM == 2)
                    const double mu_upper = mu_array_data(cc + shift_d_plus, 0);
                    const double mu_lower = mu_array_data(cc, 0);
#elif (NDIM == 3)
                    // Get edge data aligned with perp dir. (perpendicular to d and axis) and shifted in the d dir.
                    const int perp = 2 * (d + axis) % 3; // 2 if {0,1}, 1 if {0,2} and 0 if {1,2}
                    const double mu_upper = get_mu_edge(cc + shift_d_plus, perp, mu_data);
                    const double mu_lower = get_mu_edge(cc, perp, mu_data);
#endif

                    matrix_coefficients(i, stencil_map[shift_d_plus]) = (mu_upper) / (dx[d] * dx[d]);
                    matrix_coefficients(i, stencil_map[shift_d_minus]) = (mu_lower) / (dx[d] * dx[d]);
                    matrix_coefficients(i, stencil_map[ORIGIN]) -= matrix_coefficients(i, stencil_map[shift_d_plus]) +
                                                                   matrix_coefficients(i, stencil_map[shift_d_minus]);

                    matrix_coefficients(i, stencil_map[shift_d_plus + shift_axis_plus]) =
                        (mu_upper) / (dx[axis] * dx[d]);
                    matrix_coefficients(i, stencil_map[shift_d_plus + shift_axis_minus]) =
                        -(mu_upper) / (dx[axis] * dx[d]);
                    matrix_coefficients(i, stencil_map[shift_d_minus + shift_axis_plus]) =
                        -(mu_lower) / (dx[axis] * dx[d]);
                    matrix_coefficients(i, stencil_map[shift_d_minus + shift_axis_minus]) =
                        (mu_lower) / (dx[axis] * dx[d]);
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
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(
                acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
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
                const hier::Index<NDIM>& i = bc();
                const double& a = (*acoef_data)(i, 0);
                const double& b = (*bcoef_data)(i, 0);
                const double& h = dx[bdry_normal_axis];

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

                if (is_lower)
                {
                    hier::Index<NDIM> shift = get_shift(bdry_normal_axis, -1);
                    matrix_coefficients(i_s, stencil_map[ORIGIN]) +=
                        matrix_coefficients(i_s, stencil_map[shift]) * (-(a * h - 2.0 * b) / (a * h + 2.0 * b));
                    matrix_coefficients(i_s, stencil_map[shift]) = 0.0;
                }
                else
                {
                    hier::Index<NDIM> shift = get_shift(bdry_normal_axis, 1);
                    matrix_coefficients(i_s, stencil_map[ORIGIN]) +=
                        matrix_coefficients(i_s, stencil_map[shift]) * (-(a * h - 2.0 * b) / (a * h + 2.0 * b));
                    matrix_coefficients(i_s, stencil_map[shift]) = 0.0;
                }
            }
        }
    }

    if (!impose_physical_bcs_normal_comp) return;

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
                std::array<double, NDIM> shifted_patch_x_lower, shifted_patch_x_upper;
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
                auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[comp]);
                if (extended_bc_coef)
                {
                    extended_bc_coef->clearTargetPatchDataIndex();
                    extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                }
                bc_coefs[comp]->setBcCoefs(
                    acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
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
                    const hier::Index<NDIM>& i = bc();
                    const hier::Index<NDIM> i_upper = i + get_shift(comp, 1);
                    const double& a_lower = (*acoef_data)(i, 0);
                    const double& b_lower = (*bcoef_data)(i, 0);
                    const double& a_upper = (*acoef_data)(i_upper, 0);
                    const double& b_upper = (*bcoef_data)(i_upper, 0);
                    const double& h = dx[bdry_normal_axis];

                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);

                    if (is_lower)
                    {
                        hier::Index<NDIM> shift_outer_lower = get_shift(bdry_normal_axis, -1) + get_shift(comp, -1);
                        hier::Index<NDIM> shift_inner_lower = get_shift(bdry_normal_axis, 1) + get_shift(comp, -1);
                        matrix_coefficients(i_s, stencil_map[shift_inner_lower]) +=
                            matrix_coefficients(i_s, stencil_map[shift_outer_lower]) *
                            (-(a_lower * h - 2.0 * b_lower) / (a_lower * h + 2.0 * b_lower));
                        matrix_coefficients(i_s, stencil_map[shift_outer_lower]) = 0.0;

                        hier::Index<NDIM> shift_outer_upper = get_shift(bdry_normal_axis, -1) + get_shift(comp, 1);
                        hier::Index<NDIM> shift_inner_upper = get_shift(bdry_normal_axis, 1) + get_shift(comp, 1);
                        matrix_coefficients(i_s, stencil_map[shift_inner_upper]) +=
                            matrix_coefficients(i_s, stencil_map[shift_outer_upper]) *
                            (-(a_upper * h - 2.0 * b_upper) / (a_upper * h + 2.0 * b_upper));
                        matrix_coefficients(i_s, stencil_map[shift_outer_upper]) = 0.0;
                    }
                    else
                    {
                        hier::Index<NDIM> shift_outer_lower = get_shift(bdry_normal_axis, 1) + get_shift(comp, -1);
                        hier::Index<NDIM> shift_inner_lower = get_shift(bdry_normal_axis, -1) + get_shift(comp, -1);
                        matrix_coefficients(i_s, stencil_map[shift_inner_lower]) +=
                            matrix_coefficients(i_s, stencil_map[shift_outer_lower]) *
                            (-(a_lower * h - 2.0 * b_lower) / (a_lower * h + 2.0 * b_lower));
                        matrix_coefficients(i_s, stencil_map[shift_outer_lower]) = 0.0;

                        hier::Index<NDIM> shift_outer_upper = get_shift(bdry_normal_axis, 1) + get_shift(comp, 1);
                        hier::Index<NDIM> shift_inner_upper = get_shift(bdry_normal_axis, -1) + get_shift(comp, 1);
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
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(
                acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
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
                const hier::Index<NDIM>& i = bc();
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
                    TBOX_ASSERT(!IBTK::abs_equal_eps(b, 0.0));
#endif
                    if (is_lower)
                    {
                        const hier::Index<NDIM> shift_outer = get_shift(bdry_normal_axis, -1);
                        const hier::Index<NDIM> shift_inner = get_shift(bdry_normal_axis, 1);
                        matrix_coefficients(i_s, stencil_map[ORIGIN]) -=
                            matrix_coefficients(i_s, stencil_map[shift_outer]) * 2 * dx[bdry_normal_axis] * a / b;
                        matrix_coefficients(i_s, stencil_map[shift_inner]) +=
                            matrix_coefficients(i_s, stencil_map[shift_outer]);
                        matrix_coefficients(i_s, stencil_map[shift_outer]) = 0.0;
                    }
                    else
                    {
                        const hier::Index<NDIM> shift_outer = get_shift(bdry_normal_axis, 1);
                        const hier::Index<NDIM> shift_inner = get_shift(bdry_normal_axis, -1);
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
} // computeVCSCViscousOpMatrixCoefficients

void
PoissonUtilities::computeVCSCViscousDilatationalOpMatrixCoefficients(
    SideData<NDIM, double>& matrix_coefficients,
    Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    const std::vector<std::map<Index<NDIM>, int, IndexFortranOrder> >& stencil_map_vec,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
    double data_time,
    int mu_idx,
    int lambda_idx,
    VCInterpType mu_interp_type,
    bool impose_physical_bcs_normal_comp)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_coefs.size() == NDIM);
#endif
    const int stencil_sz = (2 * NDIM + 1) + 4 * (NDIM - 1);

#if !defined(NDEBUG)
    TBOX_ASSERT(matrix_coefficients.getDepth() == stencil_sz);
#endif
    matrix_coefficients.fillAll(0.0);

    Pointer<CellData<NDIM, double> > lambda_data = nullptr;
    if (lambda_idx != IBTK::invalid_index) lambda_data = patch->getPatchData(lambda_idx);

#if (NDIM == 2)
    Pointer<NodeData<NDIM, double> > mu_data = patch->getPatchData(mu_idx);
#elif (NDIM == 3)
    Pointer<EdgeData<NDIM, double> > mu_data = patch->getPatchData(mu_idx);
#endif

#if !defined(NDEBUG)
    if (lambda_idx != IBTK::invalid_index) TBOX_ASSERT(!lambda_data.isNull());
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
    using StencilMapType = std::map<hier::Index<NDIM>, int, IndexFortranOrder>;
#if (NDIM == 2)
    StencilMapType stencil_map = stencil_map_vec[0];
#endif
    static const hier::Index<NDIM> ORIGIN(0);
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
#if (NDIM == 3)
        StencilMapType stencil_map = stencil_map_vec[axis];
#endif
        for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
        {
            const hier::Index<NDIM>& cc = b();
            const SideIndex<NDIM> i(cc, axis, SideIndex<NDIM>::Lower);

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d == axis)
                {
                    const hier::Index<NDIM> shift_axis_plus = get_shift(axis, 1);
                    const hier::Index<NDIM> shift_axis_minus = get_shift(axis, -1);

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
                    const double lambda_upper = lambda_data.isNull() ? 0.0 : (*lambda_data)(cc);
                    const double lambda_lower = lambda_data.isNull() ? 0.0 : (*lambda_data)(cc + shift_axis_minus);

                    const double coef_plus = (2.0 * mu_upper + lambda_upper) / (dx[axis] * dx[axis]);
                    const double coef_minus = (2.0 * mu_lower + lambda_lower) / (dx[axis] * dx[axis]);
                    matrix_coefficients(i, stencil_map[shift_axis_plus]) = coef_plus;
                    matrix_coefficients(i, stencil_map[shift_axis_minus]) = coef_minus;
                    matrix_coefficients(i, stencil_map[ORIGIN]) -= coef_plus + coef_minus;
                }
                else
                {
                    const hier::Index<NDIM> shift_d_plus = get_shift(d, 1);
                    const hier::Index<NDIM> shift_d_minus = get_shift(d, -1);
                    const hier::Index<NDIM> shift_axis_plus = get_shift(axis, 1);
                    const hier::Index<NDIM> shift_axis_minus = get_shift(axis, -1);

#if (NDIM == 2)
                    const double mu_upper = mu_array_data(cc + shift_d_plus, 0);
                    const double mu_lower = mu_array_data(cc, 0);
#elif (NDIM == 3)
                    // Get edge data aligned with perp dir. (perpendicular to d and axis) and shifted in the d dir.
                    const int perp = 2 * (d + axis) % 3; // 2 if {0,1}, 1 if {0,2} and 0 if {1,2}
                    const double mu_upper = get_mu_edge(cc + shift_d_plus, perp, mu_data);
                    const double mu_lower = get_mu_edge(cc, perp, mu_data);
#endif

                    const double lambda_upper = lambda_data.isNull() ? 0.0 : (*lambda_data)(cc);
                    const double lambda_lower = lambda_data.isNull() ? 0.0 : (*lambda_data)(cc + shift_axis_minus);

                    matrix_coefficients(i, stencil_map[shift_d_plus]) = (mu_upper) / (dx[d] * dx[d]);
                    matrix_coefficients(i, stencil_map[shift_d_minus]) = (mu_lower) / (dx[d] * dx[d]);
                    matrix_coefficients(i, stencil_map[ORIGIN]) -= matrix_coefficients(i, stencil_map[shift_d_plus]) +
                                                                   matrix_coefficients(i, stencil_map[shift_d_minus]);

                    matrix_coefficients(i, stencil_map[shift_d_plus + shift_axis_plus]) =
                        (mu_upper + lambda_upper) / (dx[axis] * dx[d]);
                    matrix_coefficients(i, stencil_map[shift_d_plus + shift_axis_minus]) =
                        -(mu_upper + lambda_lower) / (dx[axis] * dx[d]);
                    matrix_coefficients(i, stencil_map[shift_d_minus + shift_axis_plus]) =
                        -(mu_lower + lambda_upper) / (dx[axis] * dx[d]);
                    matrix_coefficients(i, stencil_map[shift_d_minus + shift_axis_minus]) =
                        (mu_lower + lambda_lower) / (dx[axis] * dx[d]);
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
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(
                acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
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
            //
            for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
            {
                const hier::Index<NDIM>& i = bc();
                const double& a = (*acoef_data)(i, 0);
                const double& b = (*bcoef_data)(i, 0);
                const double& h = dx[bdry_normal_axis];

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

                if (is_lower)
                {
                    hier::Index<NDIM> shift = get_shift(bdry_normal_axis, -1);
                    matrix_coefficients(i_s, stencil_map[ORIGIN]) +=
                        matrix_coefficients(i_s, stencil_map[shift]) * (-(a * h - 2.0 * b) / (a * h + 2.0 * b));
                    matrix_coefficients(i_s, stencil_map[shift]) = 0.0;
                }
                else
                {
                    hier::Index<NDIM> shift = get_shift(bdry_normal_axis, 1);
                    matrix_coefficients(i_s, stencil_map[ORIGIN]) +=
                        matrix_coefficients(i_s, stencil_map[shift]) * (-(a * h - 2.0 * b) / (a * h + 2.0 * b));
                    matrix_coefficients(i_s, stencil_map[shift]) = 0.0;
                }
            }
        }
    }

    if (!impose_physical_bcs_normal_comp) return;

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
                std::array<double, NDIM> shifted_patch_x_lower, shifted_patch_x_upper;
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
                auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[comp]);
                if (extended_bc_coef)
                {
                    extended_bc_coef->clearTargetPatchDataIndex();
                    extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                }
                bc_coefs[comp]->setBcCoefs(
                    acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
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
                    const hier::Index<NDIM>& i = bc();
                    const hier::Index<NDIM> i_upper = i + get_shift(comp, 1);
                    const double& a_lower = (*acoef_data)(i, 0);
                    const double& b_lower = (*bcoef_data)(i, 0);
                    const double& a_upper = (*acoef_data)(i_upper, 0);
                    const double& b_upper = (*bcoef_data)(i_upper, 0);
                    const double& h = dx[bdry_normal_axis];

                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);

                    if (is_lower)
                    {
                        hier::Index<NDIM> shift_outer_lower = get_shift(bdry_normal_axis, -1) + get_shift(comp, -1);
                        hier::Index<NDIM> shift_inner_lower = get_shift(bdry_normal_axis, 1) + get_shift(comp, -1);
                        matrix_coefficients(i_s, stencil_map[shift_inner_lower]) +=
                            matrix_coefficients(i_s, stencil_map[shift_outer_lower]) *
                            (-(a_lower * h - 2.0 * b_lower) / (a_lower * h + 2.0 * b_lower));
                        matrix_coefficients(i_s, stencil_map[shift_outer_lower]) = 0.0;

                        hier::Index<NDIM> shift_outer_upper = get_shift(bdry_normal_axis, -1) + get_shift(comp, 1);
                        hier::Index<NDIM> shift_inner_upper = get_shift(bdry_normal_axis, 1) + get_shift(comp, 1);
                        matrix_coefficients(i_s, stencil_map[shift_inner_upper]) +=
                            matrix_coefficients(i_s, stencil_map[shift_outer_upper]) *
                            (-(a_upper * h - 2.0 * b_upper) / (a_upper * h + 2.0 * b_upper));
                        matrix_coefficients(i_s, stencil_map[shift_outer_upper]) = 0.0;
                    }
                    else
                    {
                        hier::Index<NDIM> shift_outer_lower = get_shift(bdry_normal_axis, 1) + get_shift(comp, -1);
                        hier::Index<NDIM> shift_inner_lower = get_shift(bdry_normal_axis, -1) + get_shift(comp, -1);
                        matrix_coefficients(i_s, stencil_map[shift_inner_lower]) +=
                            matrix_coefficients(i_s, stencil_map[shift_outer_lower]) *
                            (-(a_lower * h - 2.0 * b_lower) / (a_lower * h + 2.0 * b_lower));
                        matrix_coefficients(i_s, stencil_map[shift_outer_lower]) = 0.0;

                        hier::Index<NDIM> shift_outer_upper = get_shift(bdry_normal_axis, 1) + get_shift(comp, 1);
                        hier::Index<NDIM> shift_inner_upper = get_shift(bdry_normal_axis, -1) + get_shift(comp, 1);
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
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(
                acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
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
                const hier::Index<NDIM>& i = bc();
                const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                const double& a = (*acoef_data)(i, 0);
                const double& b = (*bcoef_data)(i, 0);
                if (IBTK::abs_equal_eps(b, 0.0))
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
                    TBOX_ASSERT(!IBTK::abs_equal_eps(b, 0.0));
#endif
                    if (is_lower)
                    {
                        const hier::Index<NDIM> shift_outer = get_shift(bdry_normal_axis, -1);
                        const hier::Index<NDIM> shift_inner = get_shift(bdry_normal_axis, 1);
                        matrix_coefficients(i_s, stencil_map[ORIGIN]) -=
                            matrix_coefficients(i_s, stencil_map[shift_outer]) * 2 * dx[bdry_normal_axis] * a / b;
                        matrix_coefficients(i_s, stencil_map[shift_inner]) +=
                            matrix_coefficients(i_s, stencil_map[shift_outer]);
                        matrix_coefficients(i_s, stencil_map[shift_outer]) = 0.0;
                    }
                    else
                    {
                        const hier::Index<NDIM> shift_outer = get_shift(bdry_normal_axis, 1);
                        const hier::Index<NDIM> shift_inner = get_shift(bdry_normal_axis, -1);
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
} // computeVCSCViscousDilatationalOpMatrixCoefficients

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
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[d]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[d]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
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
                const hier::Index<NDIM>& i_s_bdry = bc();
                const double& a = (*acoef_data)(i_s_bdry, 0);
                const double& b = (*bcoef_data)(i_s_bdry, 0);
                const double& g = (*gcoef_data)(i_s_bdry, 0);
                const double& h = dx[bdry_normal_axis];
                hier::Index<NDIM> i_c_intr = i_s_bdry;
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
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(
                acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
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
                const hier::Index<NDIM>& i = bc();
                const double& a = (*acoef_data)(i, 0);
                const double& b = (*bcoef_data)(i, 0);
                const double& g = (*gcoef_data)(i, 0);
                const double& h = dx[bdry_normal_axis];
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
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(
                acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
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
                const hier::Index<NDIM>& i = bc();
                const double& b = (*bcoef_data)(i, 0);
                const double& g = (*gcoef_data)(i, 0);
                const double& h = dx[bdry_normal_axis];
                const SideIndex<NDIM> i_s_bdry(i, bdry_normal_axis, SideIndex<NDIM>::Lower);
                if (b != 0.0)
                {
#if !defined(NDEBUG)
                    TBOX_ASSERT(!IBTK::abs_equal_eps(b, 0.0));
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
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(
                acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
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
                const hier::Index<NDIM>& i = bc();
                const double& a = (*acoef_data)(i, 0);
                const double& b = (*bcoef_data)(i, 0);
                const double& g = (*gcoef_data)(i, 0);
                const double& h = dx[bdry_normal_axis];
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

                const hier::Index<NDIM> shift_bdry_normal = get_shift(bdry_normal_axis, 1);
#if (NDIM == 2)
                const double mu_upper = mu_array_data(i_intr + shift_bdry_normal, 0);
                const double mu_lower = mu_array_data(i_intr, 0);
#elif (NDIM == 3)
                const int perp = 2 * (bdry_normal_axis + axis) % 3; // 2 if {0,1}, 1 if {0,2} and 0 if {1,2}
                const double mu_upper = get_mu_edge(i_intr + shift_bdry_normal, perp, mu_data);
                const double mu_lower = get_mu_edge(i_intr, perp, mu_data);
#endif

                const double D = is_lower ? mu_lower : mu_upper;
                rhs_data(i_s) += (D / h) * (-2.0 * g) / (2.0 * b + h * a);
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
                std::array<double, NDIM> shifted_patch_x_lower, shifted_patch_x_upper;
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
                auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[comp]);
                if (extended_bc_coef)
                {
                    extended_bc_coef->clearTargetPatchDataIndex();
                    extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                }
                bc_coefs[comp]->setBcCoefs(
                    acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
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
                    const hier::Index<NDIM>& i = bc();
                    const hier::Index<NDIM>& i_upper = i + get_shift(comp, 1);
                    const double& a_lower = (*acoef_data)(i, 0);
                    const double& b_lower = (*bcoef_data)(i, 0);
                    const double& g_lower = (*gcoef_data)(i, 0);
                    const double& a_upper = (*acoef_data)(i_upper, 0);
                    const double& b_upper = (*bcoef_data)(i_upper, 0);
                    const double& g_upper = (*gcoef_data)(i_upper, 0);
                    const double& h = dx[bdry_normal_axis];
                    const double& hd = dx[comp];

                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);

                    const hier::Index<NDIM> shift_d = get_shift(comp, 1);
#if (NDIM == 2)
                    const double mu_upper = mu_array_data(i + shift_d, 0);
                    const double mu_lower = mu_array_data(i, 0);
#elif (NDIM == 3)
                    const int perp = 2 * (comp + axis) % 3; // 2 if {0,1}, 1 if {0,2} and 0 if {1,2}
                    const double mu_upper = get_mu_edge(i + shift_d, perp, mu_data);
                    const double mu_lower = get_mu_edge(i, perp, mu_data);
#endif
                    if (is_lower)
                    {
                        rhs_data(i_s) -= (mu_lower / hd) * (2.0 * g_lower) / (2.0 * b_lower + a_lower * h);
                        rhs_data(i_s) += (mu_upper / hd) * (2.0 * g_upper) / (2.0 * b_upper + a_upper * h);
                    }
                    else
                    {
                        rhs_data(i_s) += (mu_lower / hd) * (2.0 * g_lower) / (2.0 * b_lower + a_lower * h);
                        rhs_data(i_s) -= (mu_upper / hd) * (2.0 * g_upper) / (2.0 * b_upper + a_upper * h);
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
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(
                acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
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
                const hier::Index<NDIM>& i = bc();
                const double& b = (*bcoef_data)(i, 0);
                const double& g = (*gcoef_data)(i, 0);
                const double& h = dx[bdry_normal_axis];
                const SideIndex<NDIM> i_s_bdry(i, bdry_normal_axis, SideIndex<NDIM>::Lower);

                const hier::Index<NDIM> shift_axis_minus = get_shift(bdry_normal_axis, -1);
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
                    TBOX_ASSERT(!IBTK::abs_equal_eps(b, 0.0));
#endif
                    rhs_data(i_s_bdry) += 2.0 * (D / h) * (-2.0 * g) / b;
                }
            }
        }
    }
    return;
} // adjustVCSCViscousOpRHSAtPhysicalBoundary

void
PoissonUtilities::adjustVCCCViscousDilatationalOpRHSAtPhysicalBoundary(
    CellData<NDIM, double>& rhs_data,
    Pointer<Patch<NDIM> > patch,
    int mu_idx,
    int lambda_idx,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
    double data_time,
    bool homogeneous_bc,
    VCInterpType mu_interp_type)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(static_cast<int>(bc_coefs.size()) == NDIM);
    TBOX_ASSERT(rhs_data.getDepth() == NDIM);
#endif

    Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_idx);
    Pointer<CellData<NDIM, double> > lambda_data = nullptr;
    if (lambda_idx != IBTK::invalid_index) lambda_data = patch->getPatchData(lambda_idx);

#if !defined(NDEBUG)
    TBOX_ASSERT(!mu_data.isNull());
    if (lambda_idx != IBTK::invalid_index) TBOX_ASSERT(!lambda_data.isNull());
#endif

    const Box<NDIM>& patch_box = patch->getBox();
    const Array<BoundaryBox<NDIM> > physical_codim1_boxes =
        PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
    const int n_physical_codim1_boxes = physical_codim1_boxes.size();

    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    auto mu_interp_cc = [&](const CellIndex<NDIM>& i, int dir, int shift) -> double
    {
        const CellIndex<NDIM> j = i + get_shift(dir, shift);
        const double q0 = (*mu_data)(i, 0);
        const double q1 = (*mu_data)(j, 0);

        switch (mu_interp_type)
        {
        case VC_AVERAGE_INTERP:
            return 0.5 * (q0 + q1);

        case VC_HARMONIC_INTERP:
            return IBTK::abs_equal_eps(q0 + q1, 0.0) ? 0.0 : 2.0 * q0 * q1 / (q0 + q1);

        default:
            TBOX_ERROR("PoissonUtilities::adjustVCCCViscousDilatationalOpRHSAtPhysicalBoundary():\n"
                       << "  unsupported shear viscosity interpolation type.\n");
            return 0.0;
        }
    };

    auto lambda_interp_cc = [&](const CellIndex<NDIM>& i, int dir, int shift) -> double
    {
        if (lambda_data.isNull()) return 0.0;

        const CellIndex<NDIM> j = i + get_shift(dir, shift);
        const double q0 = (*lambda_data)(i, 0);
        const double q1 = (*lambda_data)(j, 0);

        return 0.5 * (q0 + q1);
    };

    auto robin_ghost_value = [&](double a, double b, double g, double h) -> double
    {
        // Face-based CC closure:
        //   u_b = (u_i + u_o)/2, u_n = (u_o - u_i)/h
        // With u_i = 0 for RHS correction:
        //   u_o = 2 h g / (2 b + a h)
        return (2.0 * h * g) / (2.0 * b + a * h);
    };

    auto tangentially_extend_both = [](const Box<NDIM>& box, int axis) -> Box<NDIM>
    {
        Box<NDIM> extended_box = box;
        extended_box.lower()(axis) -= 1;
        extended_box.upper()(axis) += 1;
        return extended_box;
    };

    // ---------------------------------------------------------------------
    // Pass 1:
    // boundaries NOT aligned with the row component axis.
    //
    // Adjust rhs for the same-component transverse term
    //   d/dx_n [ mu du_axis/dx_n ]
    // where n = bdry_normal_axis != axis.
    // ---------------------------------------------------------------------
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (int n = 0; n < n_physical_codim1_boxes; ++n)
        {
            const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = (location_index % 2 == 0);

            if (bdry_normal_axis == axis) continue;

            const BoundaryBox<NDIM> trimmed_bdry_box =
                PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
            const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(
                acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

            for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
            {
                const CellIndex<NDIM>& i_bdry = bc();

                CellIndex<NDIM> i_intr = i_bdry;
                if (!is_lower) i_intr(bdry_normal_axis) -= 1;

                const double a = (*acoef_data)(i_bdry, 0);
                const double b = (*bcoef_data)(i_bdry, 0);
                const double g = (*gcoef_data)(i_bdry, 0);
                const double h = dx[bdry_normal_axis];

                const double u_o = robin_ghost_value(a, b, g, h);
                const double mu_face = mu_interp_cc(i_intr, bdry_normal_axis, is_lower ? -1 : +1);
                const double A_outer = mu_face / (h * h);

                rhs_data(i_intr, axis) -= A_outer * u_o;
            }
        }
    }

    // ---------------------------------------------------------------------
    // Pass 2:
    // boundaries aligned with row axis; adjust rhs for OTHER components
    // entering the cross terms.
    //
    // For row "axis" and comp != axis, the aligned-boundary ghost entries are:
    //   lower: (-e_axis + e_comp), (-e_axis), (-e_axis - e_comp)
    //   upper: (+e_axis + e_comp), (+e_axis), (+e_axis - e_comp)
    //
    // Correct CC elimination maps those into (+e_comp), (0), (-e_comp),
    // so the RHS must absorb the inhomogeneous ghost values at those three
    // tangential positions.
    // ---------------------------------------------------------------------
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (int n = 0; n < n_physical_codim1_boxes; ++n)
        {
            const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = (location_index % 2 == 0);

            if (bdry_normal_axis != axis) continue;

            const BoundaryBox<NDIM> trimmed_bdry_box =
                PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
            const Box<NDIM> side_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            for (unsigned int comp = 0; comp < NDIM; ++comp)
            {
                if (comp == axis) continue;

                // Need BC data at tangential shifts 0, +/- e_comp.
                const Box<NDIM> bc_coef_box = tangentially_extend_both(side_box, comp);

                Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

                auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[comp]);
                if (extended_bc_coef)
                {
                    extended_bc_coef->clearTargetPatchDataIndex();
                    extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                }
                bc_coefs[comp]->setBcCoefs(
                    acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
                if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

                const double ha = dx[axis];
                const double hd = dx[comp];

                for (Box<NDIM>::Iterator bc(side_box); bc; bc++)
                {
                    const CellIndex<NDIM>& i_bdry = bc();

                    CellIndex<NDIM> i_intr = i_bdry;
                    if (!is_lower) i_intr(axis) -= 1;

                    const CellIndex<NDIM> i_bdry_p = i_bdry + get_shift(comp, +1);
                    const CellIndex<NDIM> i_bdry_0 = i_bdry;
                    const CellIndex<NDIM> i_bdry_m = i_bdry + get_shift(comp, -1);

                    const double a_p = (*acoef_data)(i_bdry_p, 0);
                    const double b_p = (*bcoef_data)(i_bdry_p, 0);
                    const double g_p = (*gcoef_data)(i_bdry_p, 0);

                    const double a_0 = (*acoef_data)(i_bdry_0, 0);
                    const double b_0 = (*bcoef_data)(i_bdry_0, 0);
                    const double g_0 = (*gcoef_data)(i_bdry_0, 0);

                    const double a_m = (*acoef_data)(i_bdry_m, 0);
                    const double b_m = (*bcoef_data)(i_bdry_m, 0);
                    const double g_m = (*gcoef_data)(i_bdry_m, 0);

                    const double uo_p = robin_ghost_value(a_p, b_p, g_p, ha);
                    const double uo_0 = robin_ghost_value(a_0, b_0, g_0, ha);
                    const double uo_m = robin_ghost_value(a_m, b_m, g_m, ha);

                    const double mu_p = mu_interp_cc(i_intr, comp, +1);
                    const double mu_m = mu_interp_cc(i_intr, comp, -1);
                    const double lambda_n =
                        lambda_data.isNull() ? 0.0 : lambda_interp_cc(i_intr, axis, is_lower ? -1 : +1);

                    double A_p = 0.0, A_0 = 0.0, A_m = 0.0;

                    if (is_lower)
                    {
                        // Ghost coefficients in the row for:
                        //   (-e_axis + e_comp), (-e_axis), (-e_axis - e_comp)
                        A_p = -(mu_p + lambda_n) / (4.0 * ha * hd);
                        A_0 = (-mu_p + mu_m) / (4.0 * ha * hd);
                        A_m = +(mu_m + lambda_n) / (4.0 * ha * hd);
                    }
                    else
                    {
                        // Ghost coefficients in the row for:
                        //   (+e_axis + e_comp), (+e_axis), (+e_axis - e_comp)
                        A_p = +(mu_p + lambda_n) / (4.0 * ha * hd);
                        A_0 = (+mu_p - mu_m) / (4.0 * ha * hd);
                        A_m = -(mu_m + lambda_n) / (4.0 * ha * hd);
                    }

                    rhs_data(i_intr, axis) -= A_p * uo_p;
                    rhs_data(i_intr, axis) -= A_0 * uo_0;
                    rhs_data(i_intr, axis) -= A_m * uo_m;
                }
            }
        }
    }

    // ---------------------------------------------------------------------
    // Pass 3:
    // boundaries aligned with row axis; adjust rhs for the ROW component.
    //
    // Since the unknown is cell-centered, use the same face-based CC Robin
    // closure as in the matrix construction:
    //   u_b = (u_i + u_o)/2,  u_n = (u_o - u_i)/h
    // so the inhomogeneous ghost value is again u_o = 2 h g / (2 b + a h).
    // ---------------------------------------------------------------------
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (int n = 0; n < n_physical_codim1_boxes; ++n)
        {
            const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = (location_index % 2 == 0);

            if (bdry_normal_axis != axis) continue;

            const BoundaryBox<NDIM> trimmed_bdry_box =
                PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, *patch);
            const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

            Pointer<ArrayData<NDIM, double> > acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
            Pointer<ArrayData<NDIM, double> > gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(
                acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
            if (homogeneous_bc && !extended_bc_coef) gcoef_data->fillAll(0.0);

            for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
            {
                const CellIndex<NDIM>& i_bdry = bc();

                CellIndex<NDIM> i_intr = i_bdry;
                if (!is_lower) i_intr(axis) -= 1;

                const double a = (*acoef_data)(i_bdry, 0);
                const double b = (*bcoef_data)(i_bdry, 0);
                const double g = (*gcoef_data)(i_bdry, 0);
                const double h = dx[axis];

                const double u_o = robin_ghost_value(a, b, g, h);
                const double mu_face = mu_interp_cc(i_intr, axis, is_lower ? -1 : +1);
                const double lambda_face =
                    lambda_data.isNull() ? 0.0 : lambda_interp_cc(i_intr, axis, is_lower ? -1 : +1);

                const double A_outer = (2.0 * mu_face + lambda_face) / (h * h);
                rhs_data(i_intr, axis) -= A_outer * u_o;
            }
        }
    }

    return;
} // adjustVCCCViscousDilatationalOpRHSAtPhysicalBoundary

void
PoissonUtilities::adjustVCSCViscousDilatationalOpRHSAtPhysicalBoundary(
    SideData<NDIM, double>& rhs_data,
    const int rhs_depth,
    Pointer<Patch<NDIM> > patch,
    int mu_idx,
    int lambda_idx,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
    double data_time,
    bool homogeneous_bc,
    VCInterpType mu_interp_type)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(static_cast<int>(bc_coefs.size()) == NDIM);
#endif

#if (NDIM == 2)
    Pointer<NodeData<NDIM, double> > mu_data = patch->getPatchData(mu_idx);
#elif (NDIM == 3)
    Pointer<EdgeData<NDIM, double> > mu_data = patch->getPatchData(mu_idx);
#endif

    Pointer<CellData<NDIM, double> > lambda_data = nullptr;
    if (lambda_idx != IBTK::invalid_index) lambda_data = patch->getPatchData(lambda_idx);

#if !defined(NDEBUG)
    if (lambda_idx != IBTK::invalid_index) TBOX_ASSERT(!lambda_data.isNull());
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
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(
                acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
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
                const hier::Index<NDIM>& i = bc();
                const double& a = (*acoef_data)(i, 0);
                const double& b = (*bcoef_data)(i, 0);
                const double& g = (*gcoef_data)(i, 0);
                const double& h = dx[bdry_normal_axis];
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

                const hier::Index<NDIM> shift_bdry_normal = get_shift(bdry_normal_axis, 1);
#if (NDIM == 2)
                const double mu_upper = mu_array_data(i_intr + shift_bdry_normal, 0);
                const double mu_lower = mu_array_data(i_intr, 0);
#elif (NDIM == 3)
                const int perp = 2 * (bdry_normal_axis + axis) % 3; // 2 if {0,1}, 1 if {0,2} and 0 if {1,2}
                const double mu_upper = get_mu_edge(i_intr + shift_bdry_normal, perp, mu_data);
                const double mu_lower = get_mu_edge(i_intr, perp, mu_data);
#endif

                const double D = is_lower ? mu_lower : mu_upper;
                rhs_data(i_s, rhs_depth) += (D / h) * (-2.0 * g) / (2.0 * b + h * a);
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
                std::array<double, NDIM> shifted_patch_x_lower, shifted_patch_x_upper;
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
                auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[comp]);
                if (extended_bc_coef)
                {
                    extended_bc_coef->clearTargetPatchDataIndex();
                    extended_bc_coef->setHomogeneousBc(homogeneous_bc);
                }
                bc_coefs[comp]->setBcCoefs(
                    acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
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
                    const hier::Index<NDIM>& i = bc();
                    const hier::Index<NDIM>& i_upper = i + get_shift(comp, 1);
                    const double& a_lower = (*acoef_data)(i, 0);
                    const double& b_lower = (*bcoef_data)(i, 0);
                    const double& g_lower = (*gcoef_data)(i, 0);
                    const double& a_upper = (*acoef_data)(i_upper, 0);
                    const double& b_upper = (*bcoef_data)(i_upper, 0);
                    const double& g_upper = (*gcoef_data)(i_upper, 0);
                    const double& h = dx[bdry_normal_axis];
                    const double& hd = dx[comp];

                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);

                    const hier::Index<NDIM> shift_d = get_shift(comp, 1);
#if (NDIM == 2)
                    const double mu_upper = mu_array_data(i + shift_d, 0);
                    const double mu_lower = mu_array_data(i, 0);
#elif (NDIM == 3)
                    const int perp = 2 * (comp + axis) % 3; // 2 if {0,1}, 1 if {0,2} and 0 if {1,2}
                    const double mu_upper = get_mu_edge(i + shift_d, perp, mu_data);
                    const double mu_lower = get_mu_edge(i, perp, mu_data);
#endif

                    const double lambda_upper = lambda_data.isNull() ? 0.0 : (*lambda_data)(i);
                    const double lambda_lower =
                        lambda_data.isNull() ? 0.0 : (*lambda_data)(i + get_shift(bdry_normal_axis, -1));

                    if (is_lower)
                    {
                        rhs_data(i_s, rhs_depth) -=
                            ((mu_lower + lambda_lower) / hd) * (2.0 * g_lower) / (2.0 * b_lower + a_lower * h);
                        rhs_data(i_s, rhs_depth) +=
                            ((mu_upper + lambda_lower) / hd) * (2.0 * g_upper) / (2.0 * b_upper + a_upper * h);
                    }
                    else
                    {
                        rhs_data(i_s, rhs_depth) +=
                            ((mu_lower + lambda_upper) / hd) * (2.0 * g_lower) / (2.0 * b_lower + a_lower * h);
                        rhs_data(i_s, rhs_depth) -=
                            ((mu_upper + lambda_upper) / hd) * (2.0 * g_upper) / (2.0 * b_upper + a_upper * h);
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
            auto extended_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coefs[axis]);
            if (extended_bc_coef)
            {
                extended_bc_coef->clearTargetPatchDataIndex();
                extended_bc_coef->setHomogeneousBc(homogeneous_bc);
            }
            bc_coefs[axis]->setBcCoefs(
                acoef_data, bcoef_data, gcoef_data, nullptr, *patch, trimmed_bdry_box, data_time);
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
                const hier::Index<NDIM>& i = bc();
                const double& b = (*bcoef_data)(i, 0);
                const double& g = (*gcoef_data)(i, 0);
                const double& h = dx[bdry_normal_axis];
                const SideIndex<NDIM> i_s_bdry(i, bdry_normal_axis, SideIndex<NDIM>::Lower);

                const hier::Index<NDIM> shift_axis_minus = get_shift(bdry_normal_axis, -1);
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

                const double lambda_upper = lambda_data.isNull() ? 0.0 : (*lambda_data)(i);
                const double lambda_lower = lambda_data.isNull() ? 0.0 : (*lambda_data)(i + shift_axis_minus);

                const double D = is_lower ? 2.0 * mu_lower + lambda_lower : 2.0 * mu_upper + lambda_upper;

                if (b != 0.0)
                {
#if !defined(NDEBUG)
                    TBOX_ASSERT(!IBTK::abs_equal_eps(b, 0.0));
#endif
                    rhs_data(i_s_bdry, rhs_depth) += (D / h) * (-2.0 * g) / b;
                }
            }
        }
    }
    return;
} // adjustVCSCViscousDilatationalOpRHSAtPhysicalBoundary

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
                hier::Index<NDIM> i_s_bdry = bc();
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
                const hier::Index<NDIM>& i = bc();
                const double& h = dx[bdry_normal_axis];
                hier::Index<NDIM> i_intr = i, i_bdry = i;
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

                const hier::Index<NDIM> shift_bdry_normal = get_shift(bdry_normal_axis, 1);
#if (NDIM == 2)
                const double mu_upper = mu_array_data(i_intr + shift_bdry_normal, 0);
                const double mu_lower = mu_array_data(i_intr, 0);
#elif (NDIM == 3)
                const int perp = 2 * (bdry_normal_axis + axis) % 3; // 2 if {0,1}, 1 if {0,2} and 0 if {1,2}
                const double mu_upper = get_mu_edge(i_intr + shift_bdry_normal, perp, mu_data);
                const double mu_lower = get_mu_edge(i_intr, perp, mu_data);
#endif

                const double D = is_lower ? mu_lower : mu_upper;
                rhs_data(i_s_intr) -= (D / h) * sol_data(i_s_bdry) / h;
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
                    const hier::Index<NDIM>& i = bc();
                    const double& h = dx[bdry_normal_axis];
                    const double& hd = dx[comp];

                    const hier::Index<NDIM> shift_d = get_shift(comp, 1);
#if (NDIM == 2)
                    const double mu_upper = mu_array_data(i + shift_d, 0);
                    const double mu_lower = mu_array_data(i, 0);
#elif (NDIM == 3)
                    const int perp = 2 * (comp + axis) % 3; // 2 if {0,1}, 1 if {0,2} and 0 if {1,2}
                    const double mu_upper = get_mu_edge(i + shift_d, perp, mu_data);
                    const double mu_lower = get_mu_edge(i, perp, mu_data);
#endif

                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                    if (is_lower)
                    {
                        const hier::Index<NDIM> shift_axis_minus = get_shift(axis, -1);
                        const SideIndex<NDIM> sw(i + shift_axis_minus, comp, SideIndex<NDIM>::Lower);
                        const SideIndex<NDIM> nw(i + shift_axis_minus, comp, SideIndex<NDIM>::Upper);

                        rhs_data(i_s) -= (mu_lower / (h * hd)) * sol_data(sw);
                        rhs_data(i_s) += (mu_upper / (h * hd)) * sol_data(nw);
                    }
                    else
                    {
                        const SideIndex<NDIM> se(i, comp, SideIndex<NDIM>::Lower);
                        const SideIndex<NDIM> ne(i, comp, SideIndex<NDIM>::Upper);

                        rhs_data(i_s) += (mu_lower / (h * hd)) * sol_data(se);
                        rhs_data(i_s) -= (mu_upper / (h * hd)) * sol_data(ne);
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
                const hier::Index<NDIM>& i = bc();
                const double& h = dx[bdry_normal_axis];

                const hier::Index<NDIM> shift_axis_minus = get_shift(bdry_normal_axis, -1);
                const hier::Index<NDIM> shift_axis_plus = get_shift(bdry_normal_axis, 1);

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

                rhs_data(i_s) -= 2.0 * (D / h) * sol_data(i_s_bdry) / h;
            }
        }
    }
    return;
} // adjustVCSCViscousOpRHSAtCoarseFineBoundary

void
PoissonUtilities::adjustVCCCViscousDilatationalOpRHSAtCoarseFineBoundary(CellData<NDIM, double>& rhs_data,
                                                                         const CellData<NDIM, double>& sol_data,
                                                                         Pointer<Patch<NDIM> > patch,
                                                                         const int mu_idx,
                                                                         const int lambda_idx,
                                                                         const Array<BoundaryBox<NDIM> >& type1_cf_bdry,
                                                                         VCInterpType mu_interp_type)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(rhs_data.getDepth() == NDIM);
    TBOX_ASSERT(sol_data.getDepth() == NDIM);
#endif

    Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_idx);
    Pointer<CellData<NDIM, double> > lambda_data = nullptr;
    if (lambda_idx != IBTK::invalid_index) lambda_data = patch->getPatchData(lambda_idx);

#if !defined(NDEBUG)
    TBOX_ASSERT(!mu_data.isNull());
    if (lambda_idx != IBTK::invalid_index) TBOX_ASSERT(!lambda_data.isNull());
#endif

    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    const int n_cf_bdry_boxes = type1_cf_bdry.size();
    const IntVector<NDIM> ghost_width_to_fill(1);

    auto mu_interp_cc = [&](const CellIndex<NDIM>& i, const int dir, const int shift) -> double
    {
        const CellIndex<NDIM> j = i + get_shift(dir, shift);
        const double q0 = (*mu_data)(i, 0);
        const double q1 = (*mu_data)(j, 0);

        switch (mu_interp_type)
        {
        case VC_AVERAGE_INTERP:
            return 0.5 * (q0 + q1);

        case VC_HARMONIC_INTERP:
            return IBTK::abs_equal_eps(q0 + q1, 0.0) ? 0.0 : 2.0 * q0 * q1 / (q0 + q1);

        default:
            TBOX_ERROR("PoissonUtilities::adjustVCCCViscousDilatationalOpRHSAtCoarseFineBoundary():\n"
                       << "  unsupported shear viscosity interpolation type\n");
            return 0.0;
        }
    };

    auto lambda_interp_cc = [&](const CellIndex<NDIM>& i, const int dir, const int shift) -> double
    {
        if (lambda_data.isNull()) return 0.0;

        const CellIndex<NDIM> j = i + get_shift(dir, shift);
        const double q0 = (*lambda_data)(i, 0);
        const double q1 = (*lambda_data)(j, 0);

        return 0.5 * (q0 + q1);
    };

    auto tangentially_extend_both = [](const Box<NDIM>& box, int axis) -> Box<NDIM>
    {
        Box<NDIM> extended_box = box;
        extended_box.lower()(axis) -= 1;
        extended_box.upper()(axis) += 1;
        return extended_box;
    };

    // ---------------------------------------------------------------------
    // Pass 1:
    // boundaries NOT aligned with the row component axis.
    //
    // This accounts for the same-component transverse viscous term
    //   d/dx_n [ mu du_axis/dx_n ]
    // with n = bdry_normal_axis != axis.
    // ---------------------------------------------------------------------
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (int n = 0; n < n_cf_bdry_boxes; ++n)
        {
            const BoundaryBox<NDIM>& bdry_box = type1_cf_bdry[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = (location_index % 2 == 0);

            if (bdry_normal_axis == axis) continue;

            const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);

            const double h = dx[bdry_normal_axis];

            for (Box<NDIM>::Iterator bc(bc_fill_box); bc; bc++)
            {
                const CellIndex<NDIM>& i_bdry = bc();

                CellIndex<NDIM> i_intr = i_bdry;
                if (is_lower)
                {
                    i_intr(bdry_normal_axis) += 1;
                }
                else
                {
                    i_intr(bdry_normal_axis) -= 1;
                }

                // Face coefficient at the coarse-fine interface.
                const double mu_face = mu_interp_cc(i_intr, bdry_normal_axis, is_lower ? -1 : +1);
                const double A_outer = mu_face / (h * h);

                rhs_data(i_intr, axis) -= A_outer * sol_data(i_bdry, axis);
            }
        }
    }

    // ---------------------------------------------------------------------
    // Pass 2:
    // boundaries aligned with the row axis; account for ghost values of the
    // OTHER velocity components appearing in the cross terms.
    //
    // For row component "axis" and comp != axis, the ghost entries are:
    //
    // lower boundary:
    //   (-e_axis + e_comp), (-e_axis), (-e_axis - e_comp)
    //
    // upper boundary:
    //   (+e_axis + e_comp), (+e_axis), (+e_axis - e_comp)
    //
    // with coefficients:
    //
    // lower:
    //   A_p = -(mu_p + lambda_n)/(4 h hd)
    //   A_0 = (-mu_p + mu_m)/(4 h hd)
    //   A_m = +(mu_m + lambda_n)/(4 h hd)
    //
    // upper:
    //   A_p = +(mu_p + lambda_n)/(4 h hd)
    //   A_0 = (+mu_p - mu_m)/(4 h hd)
    //   A_m = -(mu_m + lambda_n)/(4 h hd)
    // ---------------------------------------------------------------------
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (int n = 0; n < n_cf_bdry_boxes; ++n)
        {
            const BoundaryBox<NDIM>& bdry_box = type1_cf_bdry[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = (location_index % 2 == 0);

            if (bdry_normal_axis != axis) continue;

            const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);

            const double h = dx[axis];

            for (unsigned int comp = 0; comp < NDIM; ++comp)
            {
                if (comp == axis) continue;

                // Need BC data at tangential shifts 0, +/- e_comp.
                const Box<NDIM> bc_coef_box = tangentially_extend_both(bc_fill_box, comp);

                const double hd = dx[comp];
                const hier::Index<NDIM> shift_axis = get_shift(axis, is_lower ? -1 : +1);
                const hier::Index<NDIM> shift_comp_p = get_shift(comp, +1);
                const hier::Index<NDIM> shift_comp_m = get_shift(comp, -1);

                for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
                {
                    const CellIndex<NDIM>& i_bdry_0 = bc();

                    CellIndex<NDIM> i_intr = i_bdry_0;
                    if (is_lower)
                    {
                        i_intr(axis) += 1;
                    }
                    else
                    {
                        i_intr(axis) -= 1;
                    }

                    const CellIndex<NDIM> i_bdry_p = i_bdry_0 + shift_comp_p;
                    const CellIndex<NDIM> i_bdry_m = i_bdry_0 + shift_comp_m;

                    const double mu_p = mu_interp_cc(i_intr, comp, +1);
                    const double mu_m = mu_interp_cc(i_intr, comp, -1);
                    const double lambda_n =
                        lambda_data.isNull() ? 0.0 : lambda_interp_cc(i_intr, axis, is_lower ? -1 : +1);

                    double A_p = 0.0, A_0 = 0.0, A_m = 0.0;

                    if (is_lower)
                    {
                        A_p = -(mu_p + lambda_n) / (4.0 * h * hd);
                        A_0 = (-mu_p + mu_m) / (4.0 * h * hd);
                        A_m = +(mu_m + lambda_n) / (4.0 * h * hd);
                    }
                    else
                    {
                        A_p = +(mu_p + lambda_n) / (4.0 * h * hd);
                        A_0 = (+mu_p - mu_m) / (4.0 * h * hd);
                        A_m = -(mu_m + lambda_n) / (4.0 * h * hd);
                    }

                    rhs_data(i_intr, axis) -= A_p * sol_data(i_bdry_p, comp);
                    rhs_data(i_intr, axis) -= A_0 * sol_data(i_bdry_0, comp);
                    rhs_data(i_intr, axis) -= A_m * sol_data(i_bdry_m, comp);
                }
            }
        }
    }

    // ---------------------------------------------------------------------
    // Pass 3:
    // boundaries aligned with the row axis; account for the row component's
    // normal ghost value in
    //   d/dx_axis [ (2 mu + lambda) du_axis/dx_axis ].
    // ---------------------------------------------------------------------
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (int n = 0; n < n_cf_bdry_boxes; ++n)
        {
            const BoundaryBox<NDIM>& bdry_box = type1_cf_bdry[n];
            const unsigned int location_index = bdry_box.getLocationIndex();
            const unsigned int bdry_normal_axis = location_index / 2;
            const bool is_lower = (location_index % 2 == 0);

            if (bdry_normal_axis != axis) continue;

            const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, ghost_width_to_fill);

            const double h = dx[axis];

            for (Box<NDIM>::Iterator bc(bc_fill_box); bc; bc++)
            {
                const CellIndex<NDIM>& i_bdry = bc();

                CellIndex<NDIM> i_intr = i_bdry;
                if (is_lower)
                {
                    i_intr(axis) += 1;
                }
                else
                {
                    i_intr(axis) -= 1;
                }

                const double mu_face = mu_interp_cc(i_intr, axis, is_lower ? -1 : +1);
                const double lambda_face =
                    lambda_data.isNull() ? 0.0 : lambda_interp_cc(i_intr, axis, is_lower ? -1 : +1);

                const double A_outer = (2.0 * mu_face + lambda_face) / (h * h);
                rhs_data(i_intr, axis) -= A_outer * sol_data(i_bdry, axis);
            }
        }
    }

    return;
} // adjustVCCCViscousDilatationalOpRHSAtCoarseFineBoundary

void
PoissonUtilities::adjustVCSCViscousDilatationalOpRHSAtCoarseFineBoundary(SideData<NDIM, double>& rhs_data,
                                                                         const SideData<NDIM, double>& sol_data,
                                                                         const int data_depth,
                                                                         Pointer<Patch<NDIM> > patch,
                                                                         const int mu_idx,
                                                                         const int lambda_idx,
                                                                         const Array<BoundaryBox<NDIM> >& type1_cf_bdry,
                                                                         VCInterpType mu_interp_type)
{
#if (NDIM == 2)
    Pointer<NodeData<NDIM, double> > mu_data = patch->getPatchData(mu_idx);
#elif (NDIM == 3)
    Pointer<EdgeData<NDIM, double> > mu_data = patch->getPatchData(mu_idx);
#endif

    Pointer<CellData<NDIM, double> > lambda_data = nullptr;
    if (lambda_idx != IBTK::invalid_index) lambda_data = patch->getPatchData(lambda_idx);

#if !defined(NDEBUG)
    if (lambda_idx != IBTK::invalid_index) TBOX_ASSERT(!lambda_data.isNull());
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
                const hier::Index<NDIM>& i = bc();
                const double& h = dx[bdry_normal_axis];
                hier::Index<NDIM> i_intr = i, i_bdry = i;
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

                const hier::Index<NDIM> shift_bdry_normal = get_shift(bdry_normal_axis, 1);
#if (NDIM == 2)
                const double mu_upper = mu_array_data(i_intr + shift_bdry_normal, 0);
                const double mu_lower = mu_array_data(i_intr, 0);
#elif (NDIM == 3)
                const int perp = 2 * (bdry_normal_axis + axis) % 3; // 2 if {0,1}, 1 if {0,2} and 0 if {1,2}
                const double mu_upper = get_mu_edge(i_intr + shift_bdry_normal, perp, mu_data);
                const double mu_lower = get_mu_edge(i_intr, perp, mu_data);
#endif

                const double D = is_lower ? mu_lower : mu_upper;
                rhs_data(i_s_intr, data_depth) -= (D / h) * sol_data(i_s_bdry, data_depth) / h;
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
                    const hier::Index<NDIM>& i = bc();
                    const double& h = dx[bdry_normal_axis];
                    const double& hd = dx[comp];

                    const hier::Index<NDIM> shift_d = get_shift(comp, 1);
#if (NDIM == 2)
                    const double mu_upper = mu_array_data(i + shift_d, 0);
                    const double mu_lower = mu_array_data(i, 0);
#elif (NDIM == 3)
                    const int perp = 2 * (comp + axis) % 3; // 2 if {0,1}, 1 if {0,2} and 0 if {1,2}
                    const double mu_upper = get_mu_edge(i + shift_d, perp, mu_data);
                    const double mu_lower = get_mu_edge(i, perp, mu_data);
#endif

                    const double lambda_upper = lambda_data.isNull() ? 0.0 : (*lambda_data)(i);
                    const double lambda_lower =
                        lambda_data.isNull() ? 0.0 : (*lambda_data)(i + get_shift(bdry_normal_axis, -1));

                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                    if (is_lower)
                    {
                        const hier::Index<NDIM> shift_axis_minus = get_shift(axis, -1);
                        const SideIndex<NDIM> sw(i + shift_axis_minus, comp, SideIndex<NDIM>::Lower);
                        const SideIndex<NDIM> nw(i + shift_axis_minus, comp, SideIndex<NDIM>::Upper);

                        rhs_data(i_s, data_depth) -= ((mu_lower + lambda_lower) / (h * hd)) * sol_data(sw, data_depth);
                        rhs_data(i_s, data_depth) += ((mu_upper + lambda_lower) / (h * hd)) * sol_data(nw, data_depth);
                    }
                    else
                    {
                        const SideIndex<NDIM> se(i, comp, SideIndex<NDIM>::Lower);
                        const SideIndex<NDIM> ne(i, comp, SideIndex<NDIM>::Upper);

                        rhs_data(i_s, data_depth) += ((mu_lower + lambda_upper) / (h * hd)) * sol_data(se, data_depth);
                        rhs_data(i_s, data_depth) -= ((mu_upper + lambda_upper) / (h * hd)) * sol_data(ne, data_depth);
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
                const hier::Index<NDIM>& i = bc();
                const double& h = dx[bdry_normal_axis];

                const hier::Index<NDIM> shift_axis_minus = get_shift(bdry_normal_axis, -1);
                const hier::Index<NDIM> shift_axis_plus = get_shift(bdry_normal_axis, 1);

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

                const double lambda_upper = lambda_data.isNull() ? 0.0 : (*lambda_data)(i);
                const double lambda_lower = lambda_data.isNull() ? 0.0 : (*lambda_data)(i + shift_axis_minus);

                const double D = is_lower ? 2.0 * mu_lower + lambda_lower : 2.0 * mu_upper + lambda_upper;

                rhs_data(i_s, data_depth) -= (D / h) * sol_data(i_s_bdry, data_depth) / h;
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
