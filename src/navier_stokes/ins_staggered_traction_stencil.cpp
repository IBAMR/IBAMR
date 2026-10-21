// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
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

#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/ShiftedPatchGeometry.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <ArrayData.h>
#include <BoundaryBox.h>
#include <Box.h>
#include <BoxArray.h>
#include <Index.h>
#include <Patch.h>
#include <RobinBcCoefStrategy.h>
#include <SideIndex.h>
#include <Variable.h>

#include <algorithm>
#include <array>
#include <optional>

#include "./ins_staggered_traction_stencil.h"

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
namespace traction_stencil
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Return the normal velocity that normal_bc_coef prescribes on the adjacent boundary at the extension of the face
// i_face, if it prescribes one. The geometry normal_geometry is that of patch centered on the normal component.
std::optional<double>
get_prescribed_normal_velocity(const hier::Index<NDIM>& i_face,
                               const unsigned int bdry_normal_axis,
                               const unsigned int adj_location_index,
                               RobinBcCoefStrategy<NDIM>* const normal_bc_coef,
                               const Patch<NDIM>& patch,
                               const ShiftedPatchGeometry& normal_geometry,
                               const double fill_time)
{
    const Box<NDIM>& patch_box = patch.getBox();
    hier::Index<NDIM> i_cell = i_face;
    i_cell(bdry_normal_axis) = std::min(i_face(bdry_normal_axis), patch_box.upper(bdry_normal_axis));
    const BoundaryBox<NDIM> adj_bdry_box(Box<NDIM>(i_cell, i_cell), 1, adj_location_index);
    Box<NDIM> adj_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(adj_bdry_box);
    adj_coef_box.lower()(bdry_normal_axis) = i_face(bdry_normal_axis);
    adj_coef_box.upper()(bdry_normal_axis) = i_face(bdry_normal_axis);

    Pointer<ArrayData<NDIM, double>> acoef_data = new ArrayData<NDIM, double>(adj_coef_box, 1);
    Pointer<ArrayData<NDIM, double>> bcoef_data = new ArrayData<NDIM, double>(adj_coef_box, 1);
    Pointer<ArrayData<NDIM, double>> gcoef_data = new ArrayData<NDIM, double>(adj_coef_box, 1);
    {
        const ShiftedPatchGeometry::Scope scope(patch, normal_geometry);
        normal_bc_coef->setBcCoefs(
            acoef_data, bcoef_data, gcoef_data, Pointer<Variable<NDIM>>(), patch, adj_bdry_box, fill_time);
    }
    const hier::Index<NDIM>& i_coef = adj_coef_box.lower();
    const double alpha = (*acoef_data)(i_coef, 0);
    if (!IBTK::rel_equal_eps(alpha, 1.0))
    {
        return std::nullopt;
    }
    return (*gcoef_data)(i_coef, 0) / alpha;
} // get_prescribed_normal_velocity
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

// Return the stencil of the normal velocity that the TRACTION condition differences along tangential_axis at the face
// i_face. Beyond a corner it is 2*u_b - u if the adjacent boundary prescribes u_b, else the linear extrapolation along
// the boundary; outside ghost_box it extrapolates from ghost_box.
NormalVelocityStencil
get_normal_velocity_stencil(hier::Index<NDIM> i_face,
                            const unsigned int bdry_normal_axis,
                            const bool bdry_is_lower,
                            const unsigned int tangential_axis,
                            RobinBcCoefStrategy<NDIM>* const normal_bc_coef,
                            const bool homogeneous_bc,
                            const Patch<NDIM>& patch,
                            const ShiftedPatchGeometry& normal_geometry,
                            const BoxArray<NDIM>& domain,
                            const Box<NDIM>& ghost_box,
                            const double fill_time)
{
    const int j = i_face(tangential_axis);
    if (PhysicalBoundaryUtilities::isBoundaryFace(i_face, bdry_normal_axis, bdry_is_lower, domain))
    {
        const int j_lower = ghost_box.lower(tangential_axis);
        const int j_upper = ghost_box.upper(tangential_axis);
        i_face(tangential_axis) = std::min(std::max(j, j_lower), j_upper);
        const SideIndex<NDIM> idx(i_face, bdry_normal_axis, SideIndex<NDIM>::Lower);
        if ((j < j_lower || j > j_upper) && j_lower < j_upper)
        {
            i_face(tangential_axis) += (j < j_lower ? +1 : -1);
            const SideIndex<NDIM> idx_inside(i_face, bdry_normal_axis, SideIndex<NDIM>::Lower);
            return { { idx, idx_inside }, { 2.0, -1.0 }, 0.0 };
        }
        return { { idx, idx }, { 1.0, 0.0 }, 0.0 };
    }

    hier::Index<NDIM> i_below = i_face, i_above = i_face;
    i_below(tangential_axis) = j - 1;
    i_above(tangential_axis) = j + 1;
    int step = 0;
    if (PhysicalBoundaryUtilities::isBoundaryFace(i_below, bdry_normal_axis, bdry_is_lower, domain))
    {
        step = -1;
    }
    else if (PhysicalBoundaryUtilities::isBoundaryFace(i_above, bdry_normal_axis, bdry_is_lower, domain))
    {
        step = +1;
    }
    else
    {
        TBOX_ERROR("INSStaggeredVelocityBcCoef::get_normal_velocity_stencil():\n"
                   << "  the normal velocity face " << i_face << " is beyond a corner of the physical boundary\n"
                   << "  but neither adjacent face along axis " << tangential_axis << " lies on the boundary\n");
    }
    hier::Index<NDIM> i_near = i_face, i_next = i_face;
    i_near(tangential_axis) = j + step;
    i_next(tangential_axis) = j + 2 * step;
    TBOX_ASSERT(ghost_box.contains(i_near));
    const SideIndex<NDIM> near_idx(i_near, bdry_normal_axis, SideIndex<NDIM>::Lower);
    const SideIndex<NDIM> next_idx(i_next, bdry_normal_axis, SideIndex<NDIM>::Lower);

    // The adjacent boundary is the one at which the boundary ends.
    const unsigned int adj_location_index = 2 * tangential_axis + (step < 0 ? 1 : 0);
    const std::optional<double> u_prescribed = get_prescribed_normal_velocity(
        i_face, bdry_normal_axis, adj_location_index, normal_bc_coef, patch, normal_geometry, fill_time);
    if (u_prescribed)
    {
        return { { near_idx, near_idx }, { -1.0, 0.0 }, homogeneous_bc ? 0.0 : 2.0 * *u_prescribed };
    }

    // A boundary segment of one face has no second face to extrapolate from.
    if (PhysicalBoundaryUtilities::isBoundaryFace(i_next, bdry_normal_axis, bdry_is_lower, domain))
    {
        TBOX_ASSERT(ghost_box.contains(i_next));
        return { { near_idx, next_idx }, { 2.0, -1.0 }, 0.0 };
    }
    return { { near_idx, near_idx }, { 1.0, 0.0 }, 0.0 };
} // get_normal_velocity_stencil

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace traction_stencil
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
