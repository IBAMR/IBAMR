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
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>
#include <Variable.h>

#include <algorithm>
#include <array>
#include <memory>

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
// Return the physical domain refined by ratio, together with its periodic
// images.
BoxArray<NDIM>
compute_physical_domain(const Pointer<GridGeometry<NDIM>>& grid_geom, const IntVector<NDIM>& ratio)
{
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
} // compute_physical_domain

// Return true if the face i_face, in the plane of the boundary normal to
// bdry_normal_axis, carries a normal velocity value, that is, if the cell
// adjacent to it on the interior side of that boundary lies in the domain.
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
} // is_boundary_face
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
NormalVelocityStencil::useAdjacentBoundaryValue()
{
#if !defined(NDEBUG)
    TBOX_ASSERT(beyond_corner);
#endif
    idx[1] = idx[0];
    weight = { -1.0, 0.0 };
    uses_adjacent_boundary_value = true;
} // useAdjacentBoundaryValue

ShiftedPatchGeometry::ShiftedPatchGeometry(const Patch<NDIM>& patch, const std::array<double, NDIM>& shift)
{
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
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
    d_geometry = new CartesianPatchGeometry<NDIM>(
        pgeom->getRatio(), touches_regular_bdry, touches_periodic_bdry, dx, x_lower.data(), x_upper.data());
} // ShiftedPatchGeometry

ShiftedPatchGeometry::Scope::Scope(const Patch<NDIM>& patch, const ShiftedPatchGeometry& shifted)
    : d_patch(const_cast<Patch<NDIM>&>(patch)), d_original(patch.getPatchGeometry())
{
    d_patch.setPatchGeometry(shifted.d_geometry);
} // Scope

ShiftedPatchGeometry::Scope::~Scope()
{
    d_patch.setPatchGeometry(d_original);
} // ~Scope

const BoxArray<NDIM>&
get_physical_domain(PhysicalDomainCache& cache,
                    const Pointer<PatchHierarchy<NDIM>>& hierarchy,
                    const IntVector<NDIM>& ratio)
{
    if (!hierarchy)
    {
        TBOX_ERROR("the fluid solver has no patch hierarchy, which is needed to locate the physical boundary.\n");
    }
    std::array<int, NDIM> key;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        key[d] = ratio(d);
    }
    auto it = cache.find(key);
    if (it == cache.end())
    {
        // The physical domain and its periodic shifts do not change during a run, so the refined domain is computed
        // once for each refinement ratio.
        it = cache.emplace(key, compute_physical_domain(hierarchy->getGridGeometry(), ratio)).first;
    }
    return it->second;
} // get_physical_domain

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

NormalVelocityStencil
get_normal_velocity_stencil(hier::Index<NDIM> i_face,
                            const unsigned int bdry_normal_axis,
                            const bool bdry_is_lower,
                            const unsigned int tangential_axis,
                            const BoxArray<NDIM>& domain,
                            const Box<NDIM>& ghost_box)
{
    const int j = i_face(tangential_axis);
    if (is_boundary_face(i_face, bdry_normal_axis, bdry_is_lower, domain))
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
        TBOX_ERROR("the normal velocity face "
                   << i_face << " is beyond a corner of the physical boundary, but neither adjacent face along axis "
                   << tangential_axis << " lies on the boundary.\n");
    }
    hier::Index<NDIM> i_near = i_face, i_next = i_face;
    i_near(tangential_axis) = j + step;
    i_next(tangential_axis) = j + 2 * step;
#if !defined(NDEBUG)
    TBOX_ASSERT(ghost_box.contains(i_near));
#endif

    // Linear extrapolation requires a second face on the boundary; a boundary segment that is one cell long has none.
    const bool can_extrapolate = is_boundary_face(i_next, bdry_normal_axis, bdry_is_lower, domain);
#if !defined(NDEBUG)
    TBOX_ASSERT(!can_extrapolate || ghost_box.contains(i_next));
#endif
    return get_corner_normal_velocity_stencil(i_face, bdry_normal_axis, tangential_axis, step, can_extrapolate);
} // get_normal_velocity_stencil

namespace
{
// Return true and set u_prescribed if normal_bc_coef prescribes the normal
// velocity on the adjacent boundary (location index adj_location_index) at the
// extension of the boundary face i_face. The geometry of patch is centered on
// the tangential component; the coefficients are evaluated with a geometry
// centered on the normal component, which is built when first needed and kept
// in normal_geometry.
bool
get_prescribed_normal_velocity(double& u_prescribed,
                               const hier::Index<NDIM>& i_face,
                               const unsigned int bdry_normal_axis,
                               const unsigned int tangential_axis,
                               const unsigned int adj_location_index,
                               RobinBcCoefStrategy<NDIM>* const normal_bc_coef,
                               const Patch<NDIM>& patch,
                               std::unique_ptr<ShiftedPatchGeometry>& normal_geometry,
                               const double fill_time)
{
    const Box<NDIM>& patch_box = patch.getBox();
    hier::Index<NDIM> i_cell = i_face;
    i_cell(bdry_normal_axis) = std::min(i_face(bdry_normal_axis), patch_box.upper(bdry_normal_axis));
    const BoundaryBox<NDIM> adj_bdry_box(Box<NDIM>(i_cell, i_cell), 1, adj_location_index);
    Box<NDIM> adj_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(adj_bdry_box);
    adj_coef_box.lower()(bdry_normal_axis) = i_face(bdry_normal_axis);
    adj_coef_box.upper()(bdry_normal_axis) = i_face(bdry_normal_axis);

    if (!normal_geometry)
    {
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
        const double* const dx = pgeom->getDx();
        std::array<double, NDIM> shift;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            shift[d] = (d == tangential_axis ? 0.5 * dx[d] : 0.0) - (d == bdry_normal_axis ? 0.5 * dx[d] : 0.0);
        }
        normal_geometry = std::make_unique<ShiftedPatchGeometry>(patch, shift);
    }
    Pointer<ArrayData<NDIM, double>> acoef_data = new ArrayData<NDIM, double>(adj_coef_box, 1);
    Pointer<ArrayData<NDIM, double>> bcoef_data = new ArrayData<NDIM, double>(adj_coef_box, 1);
    Pointer<ArrayData<NDIM, double>> gcoef_data = new ArrayData<NDIM, double>(adj_coef_box, 1);
    {
        ShiftedPatchGeometry::Scope scope(patch, *normal_geometry);
        normal_bc_coef->setBcCoefs(
            acoef_data, bcoef_data, gcoef_data, Pointer<Variable<NDIM>>(), patch, adj_bdry_box, fill_time);
    }
    const hier::Index<NDIM>& i_coef = adj_coef_box.lower();
    const double alpha = (*acoef_data)(i_coef, 0);
    const double beta = (*bcoef_data)(i_coef, 0);
    const bool velocity_bc = (alpha == 1.0 && beta == 0.0);
    const bool traction_bc = (alpha == 0.0 && beta == 1.0);
    if (!velocity_bc && !traction_bc)
    {
        TBOX_ERROR("traction_stencil::get_prescribed_normal_velocity():\n"
                   << "  unsupported boundary condition coefficients (a, b) = (" << alpha << ", " << beta
                   << ") for the normal velocity.\n"
                   << "  Only a prescribed velocity, (a, b) = (1, 0), or a prescribed traction, (a, b) = (0, 1),\n"
                   << "  is supported.\n");
    }
    if (!velocity_bc)
    {
        return false;
    }
    u_prescribed = (*gcoef_data)(i_coef, 0) / alpha;
    return true;
} // get_prescribed_normal_velocity

// The normal velocity used by the TRACTION condition at a boundary face, as
// the stencil plus offset.
struct NormalVelocity
{
    NormalVelocityStencil stencil;
    double offset;
};

// Return the normal velocity that the TRACTION condition takes at the boundary
// face with index i_face, where the boundary is normal to bdry_normal_axis and
// the condition differences the normal velocity along tangential_axis. If the
// face is beyond a corner and the adjacent boundary prescribes the normal
// velocity, the value is 2*u_b - u, where u is the nearest face on the boundary
// and u_b is the prescribed velocity, which is zero in a homogeneous condition.
NormalVelocity
get_normal_velocity(const hier::Index<NDIM>& i_face,
                    const unsigned int bdry_normal_axis,
                    const bool bdry_is_lower,
                    const unsigned int tangential_axis,
                    RobinBcCoefStrategy<NDIM>* const normal_bc_coef,
                    const bool homogeneous_bc,
                    const Patch<NDIM>& patch,
                    std::unique_ptr<ShiftedPatchGeometry>& normal_geometry,
                    const BoxArray<NDIM>& domain,
                    const Box<NDIM>& ghost_box,
                    const double fill_time)
{
    NormalVelocity u = {
        get_normal_velocity_stencil(i_face, bdry_normal_axis, bdry_is_lower, tangential_axis, domain, ghost_box), 0.0
    };
    double u_prescribed = 0.0;
    if (u.stencil.beyond_corner && get_prescribed_normal_velocity(u_prescribed,
                                                                  i_face,
                                                                  bdry_normal_axis,
                                                                  tangential_axis,
                                                                  u.stencil.adjacent_location_index,
                                                                  normal_bc_coef,
                                                                  patch,
                                                                  normal_geometry,
                                                                  fill_time))
    {
        u.stencil.useAdjacentBoundaryValue();
        u.offset = homogeneous_bc ? 0.0 : 2.0 * u_prescribed;
    }
    return u;
} // get_normal_velocity
} // namespace

double
get_normal_velocity_difference(const SideData<NDIM, double>& u_data,
                               const hier::Index<NDIM>& i,
                               const unsigned int bdry_normal_axis,
                               const bool bdry_is_lower,
                               const unsigned int tangential_axis,
                               RobinBcCoefStrategy<NDIM>* const normal_bc_coef,
                               const bool homogeneous_bc,
                               const Patch<NDIM>& patch,
                               std::unique_ptr<ShiftedPatchGeometry>& normal_geometry,
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

std::vector<std::pair<SideIndex<NDIM>, double>>
get_traction_stencil(const hier::Index<NDIM>& i,
                     const unsigned int bdry_normal_axis,
                     const bool bdry_is_lower,
                     const unsigned int tangential_axis,
                     RobinBcCoefStrategy<NDIM>* const normal_bc_coef,
                     const Patch<NDIM>& patch,
                     std::unique_ptr<ShiftedPatchGeometry>& normal_geometry,
                     const BoxArray<NDIM>& domain,
                     const Box<NDIM>& ghost_box,
                     const double fill_time)
{
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    const double sgn = bdry_is_lower ? -1.0 : +1.0;
    const double derivative_scale = sgn / pgeom->getDx()[tangential_axis];
    hier::Index<NDIM> i_lower(i);
    i_lower(tangential_axis) -= 1;
    const auto get_stencil = [&](const hier::Index<NDIM>& i_face)
    {
        return get_normal_velocity(i_face,
                                   bdry_normal_axis,
                                   bdry_is_lower,
                                   tangential_axis,
                                   normal_bc_coef,
                                   /*homogeneous_bc*/ true,
                                   patch,
                                   normal_geometry,
                                   domain,
                                   ghost_box,
                                   fill_time)
            .stencil;
    };
    const NormalVelocityStencil u_lower = get_stencil(i_lower);
    const NormalVelocityStencil u_upper = get_stencil(i);
    std::vector<std::pair<SideIndex<NDIM>, double>> stencil;
    for (int k = 0; k < 2; ++k)
    {
        if (u_upper.weight[k] != 0.0)
        {
            stencil.emplace_back(u_upper.idx[k], -u_upper.weight[k] * derivative_scale);
        }
        if (u_lower.weight[k] != 0.0)
        {
            stencil.emplace_back(u_lower.idx[k], +u_lower.weight[k] * derivative_scale);
        }
    }
    return stencil;
} // get_traction_stencil

void
accumulate_from_traction_bc_coefs(SideData<NDIM, double>& u_data,
                                  const ArrayData<NDIM, double>& gcoef_data,
                                  const unsigned int tangential_axis,
                                  RobinBcCoefStrategy<NDIM>* const tangential_bc_coef,
                                  RobinBcCoefStrategy<NDIM>* const normal_bc_coef,
                                  const bool homogeneous_bc,
                                  const Patch<NDIM>& patch,
                                  const BoundaryBox<NDIM>& bdry_box,
                                  PhysicalDomainCache& domain_cache,
                                  const Pointer<PatchHierarchy<NDIM>>& hierarchy,
                                  const double fill_time)
{
    // Determine where the physical boundary conditions prescribe the traction.
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
#if !defined(NDEBUG)
    const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(ghost_box, bdry_normal_axis);
#endif
    const BoxArray<NDIM>* domain = nullptr;
    std::unique_ptr<ShiftedPatchGeometry> normal_geometry;
    for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
    {
        const hier::Index<NDIM>& i = it();
        const double alpha = (*acoef_data)(i, 0);
        const double beta = (*bcoef_data)(i, 0);
        const bool velocity_bc = (alpha == 1.0 && beta == 0.0);
        const bool traction_bc = (alpha == 0.0 && beta == 1.0);
        if (!velocity_bc && !traction_bc)
        {
            TBOX_ERROR("traction_stencil::accumulate_from_traction_bc_coefs():\n"
                       << "  unsupported boundary condition coefficients (a, b) = (" << alpha << ", " << beta
                       << ") for the tangential velocity.\n"
                       << "  Only a prescribed velocity, (a, b) = (1, 0), or a prescribed traction, (a, b) = (0, 1),\n"
                       << "  is supported.\n");
        }
        if (!traction_bc)
        {
            continue;
        }
        if (!domain)
        {
            domain = &get_physical_domain(domain_cache, hierarchy, patch.getPatchGeometry()->getRatio());
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
                                       normal_geometry,
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
#if !defined(NDEBUG)
            TBOX_ASSERT(side_box.contains(u_upper.idx[k]));
            TBOX_ASSERT(side_box.contains(u_lower.idx[k]));
#endif
            u_data(u_upper.idx[k]) -= u_upper.weight[k] * du_transpose;
            u_data(u_lower.idx[k]) += u_lower.weight[k] * du_transpose;
        }
    }
    return;
} // accumulate_from_traction_bc_coefs

double
get_divergence_free_ghost_value(const SideData<NDIM, double>& u_data,
                                const hier::Index<NDIM>& i_g,
                                const unsigned int normal_axis,
                                const bool is_lower,
                                const double* const dx)
{
    double div_u_g = 0.0;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const SideIndex<NDIM> i_g_s_upper(i_g, axis, SideIndex<NDIM>::Upper);
        const SideIndex<NDIM> i_g_s_lower(i_g, axis, SideIndex<NDIM>::Lower);
        double u_upper = u_data(i_g_s_upper);
        double u_lower = u_data(i_g_s_lower);
        if (axis == normal_axis)
        {
            (is_lower ? u_lower : u_upper) = 0.0;
        }
        div_u_g += (u_upper - u_lower) * dx[normal_axis] / dx[axis];
    }
    return (is_lower ? +1.0 : -1.0) * div_u_g;
} // get_divergence_free_ghost_value

std::vector<std::pair<SideIndex<NDIM>, double>>
get_divergence_free_ghost_value_stencil(const hier::Index<NDIM>& i_g,
                                        const unsigned int normal_axis,
                                        const bool is_lower,
                                        const double* const dx)
{
    const double sign = is_lower ? +1.0 : -1.0;
    std::vector<std::pair<SideIndex<NDIM>, double>> stencil;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const SideIndex<NDIM> i_g_s_upper(i_g, axis, SideIndex<NDIM>::Upper);
        const SideIndex<NDIM> i_g_s_lower(i_g, axis, SideIndex<NDIM>::Lower);
        const double weight = sign * dx[normal_axis] / dx[axis];
        if (axis == normal_axis)
        {
            // The face of i_g that is farther from the boundary is not read.
            if (is_lower)
            {
                stencil.emplace_back(i_g_s_upper, weight);
            }
            else
            {
                stencil.emplace_back(i_g_s_lower, -weight);
            }
        }
        else
        {
            stencil.emplace_back(i_g_s_upper, weight);
            stencil.emplace_back(i_g_s_lower, -weight);
        }
    }
    return stencil;
} // get_divergence_free_ghost_value_stencil

std::vector<std::pair<SideIndex<NDIM>, double>>
get_normal_stress_stencil(const hier::Index<NDIM>& i,
                          const unsigned int bdry_normal_axis,
                          const bool bdry_is_lower,
                          const double* const dx)
{
    hier::Index<NDIM> i_g = i;
    if (bdry_is_lower)
    {
        i_g(bdry_normal_axis) -= 1;
    }
    const SideIndex<NDIM> i_s_boundary(i, bdry_normal_axis, SideIndex<NDIM>::Lower);
    SideIndex<NDIM> i_s_inside = i_s_boundary;
    i_s_inside(bdry_normal_axis) += bdry_is_lower ? 1 : -1;
    std::vector<std::pair<SideIndex<NDIM>, double>> stencil;
    stencil.emplace_back(i_s_inside, 1.0);
    for (const auto& entry : get_divergence_free_ghost_value_stencil(i_g, bdry_normal_axis, bdry_is_lower, dx))
    {
        stencil.emplace_back(entry.first, -entry.second);
    }
    return stencil;
} // get_normal_stress_stencil

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace traction_stencil
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
