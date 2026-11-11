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

#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <ArrayData.h>
#include <BoundaryBox.h>
#include <Box.h>
#include <BoxArray.h>
#include <CartesianPatchGeometry.h>
#include <ComponentSelector.h>
#include <Index.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchDescriptor.h>
#include <PatchHierarchy.h>
#include <RobinBcCoefStrategy.h>
#include <SideData.h>
#include <SideIndex.h>
#include <Variable.h>

#include <algorithm>
#include <array>
#include <limits>
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
// Return a copy of patch that shares its patch data. Replacing the geometry of the copy leaves patch unchanged.
Pointer<Patch<NDIM>>
copy_patch(const Patch<NDIM>& patch)
{
    Pointer<Patch<NDIM>> copy = new Patch<NDIM>(patch.getBox(), patch.getPatchDescriptor());
    copy->setPatchNumber(patch.getPatchNumber());
    copy->setPatchLevelNumber(patch.getPatchLevelNumber());
    copy->setPatchInHierarchy(patch.inHierarchy());
    copy->setPatchGeometry(patch.getPatchGeometry());
    // allocatePatchData() with an empty selector sizes the patch data array of the copy without allocating data.
    copy->allocatePatchData(ComponentSelector());
    const int num_components = patch.getPatchDescriptor()->getMaxNumberRegisteredComponents();
    for (int id = 0; id < num_components; ++id)
    {
        if (patch.checkAllocated(id))
        {
            copy->setPatchData(id, patch.getPatchData(id));
        }
    }
    return copy;
} // copy_patch

// Return the stencil of the normal velocity at the face i_face on the boundary normal to bdry_normal_axis, for the
// difference along tangential_axis. The stencil faces must lie in ghost_box.
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
        TBOX_ASSERT(ghost_box.contains(i_face));
        NormalVelocityStencil stencil;
        const SideIndex<NDIM> idx(i_face, bdry_normal_axis, SideIndex<NDIM>::Lower);
        stencil.idx = { idx, idx };
        stencil.weight = { 1.0, 0.0 };
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
                   << "  the normal velocity face " << i_face
                   << " is beyond a corner of the physical boundary with no adjacent face on the boundary along axis "
                   << tangential_axis << "\n");
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
// i_face, or std::nullopt if it prescribes the traction.
std::optional<double>
get_prescribed_normal_velocity(const hier::Index<NDIM>& i_face,
                               const unsigned int bdry_normal_axis,
                               const unsigned int adj_location_index,
                               RobinBcCoefStrategy<NDIM>* const normal_bc_coef,
                               NormalVelocityPatch& normal_velocity_patch,
                               const double fill_time)
{
    const Box<NDIM>& patch_box = normal_velocity_patch.getBox();
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
        Patch<NDIM>& shifted_patch = normal_velocity_patch.getPatch();
        const ShiftedPatchGeometry::ScopedGeometryReplacement geometry_replacement(
            shifted_patch, normal_velocity_patch.getShiftedGeometry());
        normal_bc_coef->setBcCoefs(
            acoef_data, bcoef_data, gcoef_data, Pointer<Variable<NDIM>>(), shifted_patch, adj_bdry_box, fill_time);
    }
    const hier::Index<NDIM>& i_coef = adj_coef_box.lower();
    const double alpha = (*acoef_data)(i_coef, 0);
    const double beta = (*bcoef_data)(i_coef, 0);
    if (velocity_bc::get_type(alpha, beta, "traction_stencil::get_prescribed_normal_velocity()") !=
        velocity_bc::Type::VELOCITY)
    {
        return std::nullopt;
    }
    return (*gcoef_data)(i_coef, 0);
} // get_prescribed_normal_velocity

// The normal velocity in the TRACTION condition at a boundary face: the stencil plus offset.
struct NormalVelocity
{
    NormalVelocityStencil stencil;
    double offset;
};

// Return the normal velocity in the TRACTION condition at the boundary face i_face. Beyond a corner with u_b prescribed
// on the adjacent boundary (zero if homogeneous_bc), it is 2*u_b - u, in which u is the value at the nearest boundary
// face.
NormalVelocity
get_normal_velocity(const hier::Index<NDIM>& i_face,
                    const unsigned int bdry_normal_axis,
                    const bool bdry_is_lower,
                    const unsigned int tangential_axis,
                    RobinBcCoefStrategy<NDIM>* const normal_bc_coef,
                    const bool homogeneous_bc,
                    NormalVelocityPatch& normal_velocity_patch,
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
                                                                                  normal_velocity_patch,
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
get_periodic_physical_domain(const Pointer<PatchHierarchy<NDIM>>& hierarchy, const IntVector<NDIM>& ratio)
{
    if (!hierarchy)
    {
        TBOX_ERROR("traction_stencil::get_periodic_physical_domain():\n"
                   << "  the fluid solver has no patch hierarchy\n"
                   << "  call initializePatchHierarchy() on the fluid solver before applying these boundary "
                      "conditions\n");
    }
    return PhysicalBoundaryUtilities::getPeriodicPhysicalDomain(hierarchy->getGridGeometry(), ratio);
} // get_periodic_physical_domain

NormalVelocityPatch::NormalVelocityPatch(const Patch<NDIM>& patch,
                                         const unsigned int bdry_normal_axis,
                                         const unsigned int tangential_axis)
    : d_patch(patch)
{
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    d_shift.fill(0.0);
    d_shift[tangential_axis] = 0.5 * dx[tangential_axis];
    d_shift[bdry_normal_axis] = -0.5 * dx[bdry_normal_axis];
} // NormalVelocityPatch

const Box<NDIM>&
NormalVelocityPatch::getBox() const
{
    return d_patch.getBox();
} // getBox

Patch<NDIM>&
NormalVelocityPatch::getPatch()
{
    buildShiftedPatch();
    return *d_shifted_patch;
} // getPatch

const ShiftedPatchGeometry&
NormalVelocityPatch::getShiftedGeometry()
{
    buildShiftedPatch();
    return *d_shifted_geometry;
} // getShiftedGeometry

double
get_normal_velocity_difference(const SideData<NDIM, double>& u_data,
                               const hier::Index<NDIM>& i,
                               const unsigned int bdry_normal_axis,
                               const bool bdry_is_lower,
                               const unsigned int tangential_axis,
                               RobinBcCoefStrategy<NDIM>* const normal_bc_coef,
                               const bool homogeneous_bc,
                               NormalVelocityPatch& normal_velocity_patch,
                               const BoxArray<NDIM>& domain,
                               const double fill_time)
{
    const Box<NDIM>& ghost_box = u_data.getGhostBox();
    if (i(tangential_axis) - 1 < ghost_box.lower(tangential_axis) ||
        i(tangential_axis) > ghost_box.upper(tangential_axis))
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
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
                                   normal_velocity_patch,
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

/////////////////////////////// PRIVATE //////////////////////////////////////

void
NormalVelocityPatch::buildShiftedPatch()
{
    if (!d_shifted_patch)
    {
        d_shifted_patch = copy_patch(d_patch);
        d_shifted_geometry.emplace(*d_shifted_patch, d_shift);
    }
} // buildShiftedPatch

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace traction_stencil
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
