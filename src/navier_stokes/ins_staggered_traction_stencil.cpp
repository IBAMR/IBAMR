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
#include <ibtk/ibtk_utilities.h>

#include <tbox/Array.h>
#include <tbox/MathUtilities.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <ArrayData.h>
#include <BoundaryBox.h>
#include <Box.h>
#include <BoxArray.h>
#include <BoxList.h>
#include <CartesianPatchGeometry.h>
#include <GridGeometry.h>
#include <Index.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <RobinBcCoefStrategy.h>
#include <SideIndex.h>

#include <algorithm>
#include <array>
#include <functional>

#include "./ins_staggered_traction_stencil.h"

#include <ibamr/namespaces.h> // IWYU pragma: keep

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Variable;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
namespace traction_stencil
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Make a patch with the box and patch descriptor of patch whose geometry is
// that of patch with x_lower and x_upper shifted by shift. Boundary condition
// objects locate their coefficients using the geometry of the patch, so a
// shifted patch makes them evaluate the coefficients at shifted locations.
Pointer<Patch<NDIM>>
make_shifted_patch(const Patch<NDIM>& patch, const std::array<double, NDIM>& shift)
{
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    std::array<double, NDIM> x_lower, x_upper;
    Array<Array<bool>> touches_regular_bdry(NDIM), touches_periodic_bdry(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        x_lower[d] = pgeom->getXLower()[d] + shift[d];
        x_upper[d] = pgeom->getXUpper()[d] + shift[d];
        touches_regular_bdry[d].resizeArray(2);
        touches_periodic_bdry[d].resizeArray(2);
        for (int upperlower = 0; upperlower < 2; ++upperlower)
        {
            touches_regular_bdry[d][upperlower] = pgeom->getTouchesRegularBoundary(d, upperlower);
            touches_periodic_bdry[d][upperlower] = pgeom->getTouchesPeriodicBoundary(d, upperlower);
        }
    }
    Pointer<Patch<NDIM>> shifted_patch = new Patch<NDIM>(patch.getBox(), patch.getPatchDescriptor());
    shifted_patch->setPatchGeometry(new CartesianPatchGeometry<NDIM>(pgeom->getRatio(),
                                                                     touches_regular_bdry,
                                                                     touches_periodic_bdry,
                                                                     pgeom->getDx(),
                                                                     x_lower.data(),
                                                                     x_upper.data()));
    return shifted_patch;
}

// Return true and set u_prescribed if normal_bc_coef prescribes the normal
// velocity on the adjacent boundary (location index adj_location_index) at the
// extension of the boundary face i_face. The geometry of patch is centered on
// the tangential component; the coefficients are evaluated with a geometry
// centered on the normal component.
bool
get_prescribed_normal_velocity(double& u_prescribed,
                               const hier::Index<NDIM>& i_face,
                               const unsigned int bdry_normal_axis,
                               const unsigned int tangential_axis,
                               const unsigned int adj_location_index,
                               RobinBcCoefStrategy<NDIM>* const normal_bc_coef,
                               const Patch<NDIM>& patch,
                               const double fill_time)
{
    const Box<NDIM>& patch_box = patch.getBox();
    hier::Index<NDIM> i_cell = i_face;
    i_cell(bdry_normal_axis) = std::min(i_face(bdry_normal_axis), patch_box.upper(bdry_normal_axis));
    const BoundaryBox<NDIM> adj_bdry_box(Box<NDIM>(i_cell, i_cell), 1, adj_location_index);
    Box<NDIM> adj_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(adj_bdry_box);
    adj_coef_box.lower()(bdry_normal_axis) = i_face(bdry_normal_axis);
    adj_coef_box.upper()(bdry_normal_axis) = i_face(bdry_normal_axis);

    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    std::array<double, NDIM> shift;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        shift[d] = (d == tangential_axis ? 0.5 * dx[d] : 0.0) - (d == bdry_normal_axis ? 0.5 * dx[d] : 0.0);
    }
    Pointer<Patch<NDIM>> adj_patch = make_shifted_patch(patch, shift);
    Pointer<ArrayData<NDIM, double>> acoef_data = new ArrayData<NDIM, double>(adj_coef_box, 1);
    Pointer<ArrayData<NDIM, double>> bcoef_data = new ArrayData<NDIM, double>(adj_coef_box, 1);
    Pointer<ArrayData<NDIM, double>> gcoef_data = new ArrayData<NDIM, double>(adj_coef_box, 1);
    normal_bc_coef->setBcCoefs(
        acoef_data, bcoef_data, gcoef_data, Pointer<Variable<NDIM>>(), *adj_patch, adj_bdry_box, fill_time);
    const hier::Index<NDIM>& i_coef = adj_coef_box.lower();
    const double alpha = (*acoef_data)(i_coef, 0);
    if (!IBTK::rel_equal_eps(alpha, 1.0))
    {
        return false;
    }
    u_prescribed = (*gcoef_data)(i_coef, 0) / alpha;
    return true;
}

// Return true if the normal velocity face i_face lies on the physical boundary
// normal to bdry_normal_axis, i.e., if the cell adjacent to the face on the
// interior side of that boundary lies in domain.
bool
is_boundary_face(hier::Index<NDIM> i_face,
                 const unsigned int bdry_normal_axis,
                 const bool bdry_is_lower,
                 const BoxArray<NDIM>& domain)
{
    if (!bdry_is_lower)
    {
        i_face(bdry_normal_axis) -= 1;
    }
    return domain.contains(i_face);
}

} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

BoxArray<NDIM>
get_physical_domain(const Pointer<PatchHierarchy<NDIM>>& hierarchy, const Patch<NDIM>& patch, const char* const method)
{
    if (!hierarchy)
    {
        TBOX_ERROR(method << ": the fluid solver has no patch hierarchy; a TRACTION boundary condition requires the "
                             "hierarchy to locate the corners of the physical boundary.\n");
    }
    Pointer<GridGeometry<NDIM>> grid_geom = hierarchy->getGridGeometry();
    const IntVector<NDIM>& ratio = patch.getPatchGeometry()->getRatio();
    BoxArray<NDIM> domain = grid_geom->getPhysicalDomain();
    domain.refine(ratio);
    const IntVector<NDIM> shift = grid_geom->getPeriodicShift(ratio);
    BoxList<NDIM> domain_list(domain);
    int num_offsets = 1;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        num_offsets *= 3;
    }
    for (int offset = 0; offset < num_offsets; ++offset)
    {
        IntVector<NDIM> image_shift(0);
        bool is_image = false;
        bool is_valid = true;
        for (unsigned int d = 0, k = offset; d < NDIM; ++d, k /= 3)
        {
            const int n = static_cast<int>(k % 3) - 1;
            image_shift(d) = n * shift(d);
            is_image = is_image || n != 0;
            is_valid = is_valid && (n == 0 || shift(d) != 0);
        }
        if (is_image && is_valid)
        {
            BoxArray<NDIM> image_domain = domain;
            for (int b = 0; b < image_domain.size(); ++b)
            {
                image_domain[b].shift(image_shift);
            }
            domain_list.unionBoxes(BoxList<NDIM>(image_domain));
        }
    }
    return BoxArray<NDIM>(domain_list);
} // get_physical_domain

NormalVelocityStencil
get_end_normal_velocity_stencil(hier::Index<NDIM> i_face,
                                const unsigned int bdry_normal_axis,
                                const unsigned int tangential_axis,
                                const int step,
                                const TractionBcCornerType corner_type,
                                const std::function<bool(double&, unsigned int)>& get_adjacent_wall_velocity)
{
    const int j = i_face(tangential_axis);
    hier::Index<NDIM> i_near = i_face, i_next = i_face;
    i_near(tangential_axis) = j + step;
    i_next(tangential_axis) = j + 2 * step;
    const SideIndex<NDIM> near_idx(i_near, bdry_normal_axis, SideIndex<NDIM>::Lower);
    const SideIndex<NDIM> next_idx(i_next, bdry_normal_axis, SideIndex<NDIM>::Lower);
    const NormalVelocityStencil linear_extrapolation = { { near_idx, next_idx }, { 2.0, -1.0 }, 0.0 };
    switch (corner_type)
    {
    case TractionBcCornerType::ZERO_DIFFERENCE:
        return { { near_idx, near_idx }, { 1.0, 0.0 }, 0.0 };
    case TractionBcCornerType::LINEAR_EXTRAPOLATION:
        return linear_extrapolation;
    case TractionBcCornerType::WALL_AWARE:
    {
        const unsigned int adj_location_index = 2 * tangential_axis + (step < 0 ? 1 : 0);
        double u_prescribed = 0.0;
        if (get_adjacent_wall_velocity(u_prescribed, adj_location_index))
        {
            return { { near_idx, near_idx }, { -1.0, 0.0 }, 2.0 * u_prescribed, true, adj_location_index };
        }
        return linear_extrapolation;
    }
    }
    return { { near_idx, near_idx }, { 1.0, 0.0 }, 0.0 };
} // get_end_normal_velocity_stencil

NormalVelocityStencil
get_normal_velocity_stencil(hier::Index<NDIM> i_face,
                            const unsigned int bdry_normal_axis,
                            const bool bdry_is_lower,
                            const unsigned int tangential_axis,
                            const TractionBcCornerType corner_type,
                            RobinBcCoefStrategy<NDIM>* const normal_bc_coef,
                            const bool homogeneous_bc,
                            const Patch<NDIM>& patch,
                            const BoxArray<NDIM>& domain,
                            const Box<NDIM>& ghost_box,
                            const double fill_time)
{
    const int j = i_face(tangential_axis);
    if (is_boundary_face(i_face, bdry_normal_axis, bdry_is_lower, domain))
    {
        i_face(tangential_axis) =
            std::min(std::max(j, ghost_box.lower(tangential_axis)), ghost_box.upper(tangential_axis));
        const SideIndex<NDIM> idx(i_face, bdry_normal_axis, SideIndex<NDIM>::Lower);
        return { { idx, idx }, { 1.0, 0.0 }, 0.0 };
    }

    // The nearest two faces on the boundary: below the face if the boundary ends
    // above it, above the face if the boundary ends below it.
    hier::Index<NDIM> i_below = i_face, i_above = i_face;
    i_below(tangential_axis) = j - 1;
    i_above(tangential_axis) = j + 1;
    int step = 0;
    if (is_boundary_face(i_below, bdry_normal_axis, bdry_is_lower, domain))
    {
        step = -1;
    }
    else if (is_boundary_face(i_above, bdry_normal_axis, bdry_is_lower, domain))
    {
        step = +1;
    }
    else
    {
        TBOX_ERROR("traction_stencil::get_normal_velocity_stencil(): the normal velocity face "
                   << i_face << " is beyond a corner of the physical boundary, but neither adjacent face along axis "
                   << tangential_axis << " lies on the boundary.\n");
    }
#if !defined(NDEBUG)
    hier::Index<NDIM> i_near = i_face, i_next = i_face;
    i_near(tangential_axis) = j + step;
    i_next(tangential_axis) = j + 2 * step;
    TBOX_ASSERT(ghost_box.contains(i_near) && ghost_box.contains(i_next));
#endif
    return get_end_normal_velocity_stencil(i_face,
                                           bdry_normal_axis,
                                           tangential_axis,
                                           step,
                                           corner_type,
                                           [&](double& u_prescribed, const unsigned int adj_location_index)
                                           {
                                               const bool prescribed =
                                                   get_prescribed_normal_velocity(u_prescribed,
                                                                                  i_face,
                                                                                  bdry_normal_axis,
                                                                                  tangential_axis,
                                                                                  adj_location_index,
                                                                                  normal_bc_coef,
                                                                                  patch,
                                                                                  fill_time);
                                               if (prescribed && homogeneous_bc)
                                               {
                                                   u_prescribed = 0.0;
                                               }
                                               return prescribed;
                                           });
} // get_normal_velocity_stencil

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace traction_stencil
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
