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
#include <CartesianPatchGeometry.h>
#include <Index.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <RobinBcCoefStrategy.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>
#include <Variable.h>

#include <algorithm>
#include <array>
#include <optional>

#include "./ins_staggered_traction_stencil.h"
#include "./ins_velocity_bc.h"

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
namespace traction_stencil
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Return the stencil of the normal velocity at the face i_face on the boundary normal to bdry_normal_axis, for a
// condition that differences the normal velocity along tangential_axis. A face beyond a corner is found from domain,
// one outside ghost_box is extrapolated from the two nearest faces in ghost_box.
NormalVelocityStencil
get_normal_velocity_stencil(hier::Index<NDIM> i_face,
                            const unsigned int bdry_normal_axis,
                            const bool bdry_is_lower,
                            const unsigned int tangential_axis,
                            const BoxArray<NDIM>& domain,
                            const Box<NDIM>& ghost_box)
{
    const int j = i_face(tangential_axis);
    if (PhysicalBoundaryUtilities::isBoundaryFace(i_face, bdry_normal_axis, bdry_is_lower, domain))
    {
        NormalVelocityStencil stencil;
        const int j_lower = ghost_box.lower(tangential_axis);
        const int j_upper = ghost_box.upper(tangential_axis);
        i_face(tangential_axis) = std::min(std::max(j, j_lower), j_upper);
        const SideIndex<NDIM> idx(i_face, bdry_normal_axis, SideIndex<NDIM>::Lower);
        if ((j < j_lower || j > j_upper) && j_lower < j_upper)
        {
            i_face(tangential_axis) += (j < j_lower ? +1 : -1);
            const SideIndex<NDIM> idx_inside(i_face, bdry_normal_axis, SideIndex<NDIM>::Lower);
            stencil.idx = { idx, idx_inside };
            stencil.weight = { 2.0, -1.0 };
        }
        else
        {
            stencil.idx = { idx, idx };
            stencil.weight = { 1.0, 0.0 };
        }
        return stencil;
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
        TBOX_ERROR("traction_stencil::get_normal_velocity_stencil():\n"
                   << "  the normal velocity face " << i_face << " is beyond a corner of the physical boundary\n"
                   << "  but neither adjacent face along axis " << tangential_axis << " lies on the boundary\n");
    }
    hier::Index<NDIM> i_near = i_face, i_next = i_face;
    i_near(tangential_axis) = j + step;
    i_next(tangential_axis) = j + 2 * step;
    TBOX_ASSERT(ghost_box.contains(i_near));

    // A boundary segment of one face has no second face to extrapolate from.
    const bool can_extrapolate =
        PhysicalBoundaryUtilities::isBoundaryFace(i_next, bdry_normal_axis, bdry_is_lower, domain);
    TBOX_ASSERT(!can_extrapolate || ghost_box.contains(i_next));
    return get_corner_normal_velocity_stencil(i_face, bdry_normal_axis, tangential_axis, step, can_extrapolate);
} // get_normal_velocity_stencil

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
    const double beta = (*bcoef_data)(i_coef, 0);
    if (velocity_bc::get_type(alpha, beta, "traction_stencil::get_prescribed_normal_velocity()") !=
        velocity_bc::Type::VELOCITY)
    {
        return std::nullopt;
    }
    return (*gcoef_data)(i_coef, 0) / alpha;
} // get_prescribed_normal_velocity

// The normal velocity used by the TRACTION condition at a boundary face: the stencil plus offset.
struct NormalVelocity
{
    NormalVelocityStencil stencil;
    double offset;
};

// Return the normal velocity that the TRACTION condition takes at the boundary face i_face: beyond a corner whose
// adjacent boundary prescribes the normal velocity u_b (zero if homogeneous_bc), 2*u_b - u, with u the nearest face on
// the boundary.
NormalVelocity
get_normal_velocity(const hier::Index<NDIM>& i_face,
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
    NormalVelocity u = {
        get_normal_velocity_stencil(i_face, bdry_normal_axis, bdry_is_lower, tangential_axis, domain, ghost_box), 0.0
    };
    if (u.stencil.beyond_corner)
    {
        const std::optional<double> u_prescribed = get_prescribed_normal_velocity(i_face,
                                                                                  bdry_normal_axis,
                                                                                  u.stencil.adjacent_location_index,
                                                                                  normal_bc_coef,
                                                                                  patch,
                                                                                  normal_geometry,
                                                                                  fill_time);
        if (u_prescribed)
        {
            u.stencil.useAdjacentBoundaryValue();
            u.offset = homogeneous_bc ? 0.0 : 2.0 * *u_prescribed;
        }
    }
    return u;
} // get_normal_velocity
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
NormalVelocityStencil::useAdjacentBoundaryValue()
{
    TBOX_ASSERT(beyond_corner);
    idx[1] = idx[0];
    weight = { -1.0, 0.0 };
    uses_adjacent_boundary_value = true;
} // useAdjacentBoundaryValue

NormalVelocityStencil
get_corner_normal_velocity_stencil(hier::Index<NDIM> i_face,
                                   const unsigned int bdry_normal_axis,
                                   const unsigned int tangential_axis,
                                   const int step,
                                   const bool can_extrapolate)
{
    const int j = i_face(tangential_axis);
    hier::Index<NDIM> i_near = i_face, i_next = i_face;
    i_near(tangential_axis) = j + step;
    i_next(tangential_axis) = j + 2 * step;
    const SideIndex<NDIM> near_idx(i_near, bdry_normal_axis, SideIndex<NDIM>::Lower);
    const SideIndex<NDIM> next_idx(i_next, bdry_normal_axis, SideIndex<NDIM>::Lower);
    NormalVelocityStencil stencil;
    stencil.idx = { near_idx, can_extrapolate ? next_idx : near_idx };
    stencil.weight = can_extrapolate ? std::array<double, 2>{ 2.0, -1.0 } : std::array<double, 2>{ 1.0, 0.0 };
    stencil.beyond_corner = true;
    stencil.adjacent_location_index = 2 * tangential_axis + (step < 0 ? 1 : 0);
    return stencil;
} // get_corner_normal_velocity_stencil

BoxArray<NDIM>
get_physical_domain(const Pointer<PatchHierarchy<NDIM>>& hierarchy, const IntVector<NDIM>& ratio)
{
    if (!hierarchy)
    {
        TBOX_ERROR("traction_stencil::get_physical_domain():\n"
                   << "  the fluid solver has no patch hierarchy\n"
                   << "  call initializePatchHierarchy() on the fluid solver before applying these boundary "
                      "conditions\n");
    }
    return PhysicalBoundaryUtilities::getPeriodicPhysicalDomain(hierarchy->getGridGeometry(), ratio);
} // get_physical_domain

ShiftedPatchGeometry
get_normal_velocity_geometry(const Patch<NDIM>& patch,
                             const unsigned int bdry_normal_axis,
                             const unsigned int tangential_axis)
{
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    std::array<double, NDIM> shift;
    shift.fill(0.0);
    shift[tangential_axis] = 0.5 * dx[tangential_axis];
    shift[bdry_normal_axis] = -0.5 * dx[bdry_normal_axis];
    return ShiftedPatchGeometry(patch, shift);
} // get_normal_velocity_geometry

double
get_normal_velocity_difference(const SideData<NDIM, double>& u_data,
                               const hier::Index<NDIM>& i,
                               const unsigned int bdry_normal_axis,
                               const bool bdry_is_lower,
                               const unsigned int tangential_axis,
                               RobinBcCoefStrategy<NDIM>* const normal_bc_coef,
                               const bool homogeneous_bc,
                               const Patch<NDIM>& patch,
                               const ShiftedPatchGeometry& normal_geometry,
                               const BoxArray<NDIM>& domain,
                               const double fill_time)
{
    const Box<NDIM>& ghost_box = u_data.getGhostBox();
    hier::Index<NDIM> i_lower(i);
    i_lower(tangential_axis) -= 1;
    const auto get_u = [&](const hier::Index<NDIM>& i_face)
    {
        return get_normal_velocity(i_face,
                                   bdry_normal_axis,
                                   bdry_is_lower,
                                   tangential_axis,
                                   normal_bc_coef,
                                   homogeneous_bc,
                                   patch,
                                   normal_geometry,
                                   domain,
                                   ghost_box,
                                   fill_time);
    };
    const NormalVelocity u_lower = get_u(i_lower);
    const NormalVelocity u_upper = get_u(i);
    double du = u_upper.offset - u_lower.offset;
    for (int k = 0; k < 2; ++k)
    {
        du += u_upper.stencil.weight[k] * u_data(u_upper.stencil.idx[k]) -
              u_lower.stencil.weight[k] * u_data(u_lower.stencil.idx[k]);
    }
    return du;
} // get_normal_velocity_difference

void
accumulate_from_traction_bc_coefs(SideData<NDIM, double>& u_data,
                                  const ArrayData<NDIM, double>& gcoef_data,
                                  const unsigned int tangential_axis,
                                  RobinBcCoefStrategy<NDIM>* const tangential_bc_coef,
                                  RobinBcCoefStrategy<NDIM>* const normal_bc_coef,
                                  const bool homogeneous_bc,
                                  const Patch<NDIM>& patch,
                                  const BoundaryBox<NDIM>& bdry_box,
                                  const Pointer<PatchHierarchy<NDIM>>& hierarchy,
                                  const double fill_time)
{
    const Box<NDIM>& bc_coef_box = gcoef_data.getBox();
    Pointer<ArrayData<NDIM, double>> acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
    Pointer<ArrayData<NDIM, double>> bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
    Pointer<ArrayData<NDIM, double>> gcoef_unused;
    tangential_bc_coef->setBcCoefs(
        acoef_data, bcoef_data, gcoef_unused, Pointer<Variable<NDIM>>(), patch, bdry_box, fill_time);

    const unsigned int location_index = bdry_box.getLocationIndex();
    const unsigned int bdry_normal_axis = location_index / 2;
    const bool bdry_is_lower = location_index % 2 == 0;
    const double sgn = bdry_is_lower ? -1.0 : +1.0;
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const Box<NDIM>& ghost_box = u_data.getGhostBox();
    const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(ghost_box, bdry_normal_axis);
    std::optional<BoxArray<NDIM>> domain;
    std::optional<ShiftedPatchGeometry> normal_geometry;
    for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
    {
        const hier::Index<NDIM>& i = it();
        if (velocity_bc::get_type(
                (*acoef_data)(i, 0), (*bcoef_data)(i, 0), "traction_stencil::accumulate_from_traction_bc_coefs()") !=
            velocity_bc::Type::TRACTION)
        {
            continue;
        }
        if (!domain)
        {
            domain = get_physical_domain(hierarchy, pgeom->getRatio());
            normal_geometry.emplace(get_normal_velocity_geometry(patch, bdry_normal_axis, tangential_axis));
        }
        hier::Index<NDIM> i_lower(i);
        i_lower(tangential_axis) -= 1;
        const auto get_stencil = [&](const hier::Index<NDIM>& i_face)
        {
            return get_normal_velocity(i_face,
                                       bdry_normal_axis,
                                       bdry_is_lower,
                                       tangential_axis,
                                       normal_bc_coef,
                                       homogeneous_bc,
                                       patch,
                                       *normal_geometry,
                                       *domain,
                                       ghost_box,
                                       fill_time)
                .stencil;
        };
        const NormalVelocityStencil u_lower = get_stencil(i_lower);
        const NormalVelocityStencil u_upper = get_stencil(i);
        const double du_transpose = sgn * gcoef_data(i, 0) / dx[tangential_axis];
        for (int k = 0; k < 2; ++k)
        {
            TBOX_ASSERT(side_box.contains(u_upper.idx[k]));
            TBOX_ASSERT(side_box.contains(u_lower.idx[k]));
            u_data(u_upper.idx[k]) -= u_upper.weight[k] * du_transpose;
            u_data(u_lower.idx[k]) += u_lower.weight[k] * du_transpose;
        }
    }
    return;
} // accumulate_from_traction_bc_coefs

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace traction_stencil
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
