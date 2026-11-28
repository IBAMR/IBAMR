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

#include <ibamr/StokesBcCoefStrategy.h>
#include <ibamr/ibamr_enums.h>
#include <ibamr/private/ins_staggered_traction_stencil.h>
#include <ibamr/private/ins_velocity_bc.h>

#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/private/ScopedShiftedPatchGeometry.h>

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
#include <functional>
#include <limits>
#include <optional>
#include <vector>

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
namespace detail
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Return a copy of patch that shares its patch data. Replacing the geometry of the copy leaves patch unchanged.
Pointer<Patch<NDIM>>
copy_patch(const Patch<NDIM>& patch)
{
    Pointer<Patch<NDIM>> copy = new Patch<NDIM>(patch.getBox(), patch.getPatchDescriptor());
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

// Return the displacement of the geometry of patch from the faces of velocity component tangential_axis to the faces of
// component bdry_normal_axis: half a cell along tangential_axis and minus half a cell along bdry_normal_axis.
std::array<double, NDIM>
get_half_cell_shift(const Patch<NDIM>& patch, const unsigned int bdry_normal_axis, const unsigned int tangential_axis)
{
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    std::array<double, NDIM> shift{};
    shift[tangential_axis] = 0.5 * dx[tangential_axis];
    shift[bdry_normal_axis] = -0.5 * dx[bdry_normal_axis];
    return shift;
} // get_half_cell_shift

// One term (u(upper) - u(lower)) * dx_normal / dx_axis of the discrete divergence in a ghost cell. A face that is not
// used contributes zero.
struct DivergenceFreeGhostTerm
{
    SideIndex<NDIM> upper;
    SideIndex<NDIM> lower;
    bool use_upper = true;
    bool use_lower = true;
    double dx_normal = 0.0;
    double dx_axis = 0.0;
};

// The divergence-free ghost value is sign times the sum of the terms.
struct DivergenceFreeGhostFormula
{
    double sign = 1.0;
    std::array<DivergenceFreeGhostTerm, NDIM> terms;
};

// Return the faces, signs, and scales that define get_divergence_free_ghost_value() and its stencil.
DivergenceFreeGhostFormula
get_divergence_free_ghost_formula(const hier::Index<NDIM>& i_g,
                                  const unsigned int bdry_normal_axis,
                                  const bool bdry_is_lower,
                                  const double* const dx)
{
    DivergenceFreeGhostFormula formula;
    formula.sign = bdry_is_lower ? +1.0 : -1.0;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        DivergenceFreeGhostTerm& term = formula.terms[axis];
        term.upper = SideIndex<NDIM>(i_g, axis, SideIndex<NDIM>::Upper);
        term.lower = SideIndex<NDIM>(i_g, axis, SideIndex<NDIM>::Lower);
        // The face of i_g that is farther from the boundary is not used.
        term.use_upper = axis != bdry_normal_axis || bdry_is_lower;
        term.use_lower = axis != bdry_normal_axis || !bdry_is_lower;
        term.dx_normal = dx[bdry_normal_axis];
        term.dx_axis = dx[axis];
    }
    return formula;
} // get_divergence_free_ghost_formula

} // namespace

void
NormalVelocityStencil::useAdjacentBoundaryValue()
{
    TBOX_ASSERT(beyond_corner);
    idx[1] = idx[0];
    weight = { -1.0, 0.0 };
    uses_adjacent_boundary_value = true;
} // useAdjacentBoundaryValue

// Return the stencil of the normal velocity at the face i_face beyond a corner of the boundary normal to
// bdry_normal_axis, where the boundary segment ends along tangential_axis. The nearest face on the boundary is i_face +
// step along tangential_axis. The stencil also uses the second nearest face if can_extrapolate is true.
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

namespace
{
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
        TBOX_ERROR("get_normal_velocity_stencil():\n"
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

// Return the normal velocity that the normal velocity patch prescribes on the adjacent boundary at the extension of the
// face i_face, or std::nullopt if it prescribes the traction.
std::optional<double>
get_prescribed_normal_velocity(const hier::Index<NDIM>& i_face,
                               const unsigned int adj_location_index,
                               RobinBcCoefStrategy<NDIM>* const normal_bc_coef,
                               NormalVelocityPatch& normal_velocity_patch,
                               const double fill_time)
{
    const unsigned int bdry_normal_axis = normal_velocity_patch.getBdryNormalAxis();
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
    normal_velocity_patch.setBcCoefs(normal_bc_coef, acoef_data, bcoef_data, gcoef_data, adj_bdry_box, fill_time);
    const hier::Index<NDIM>& i_coef = adj_coef_box.lower();
    const double alpha = (*acoef_data)(i_coef, 0);
    const double beta = (*bcoef_data)(i_coef, 0);
    if (classify_velocity_bc(alpha, beta) != VelocityBcType::VELOCITY)
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
                    const bool bdry_is_lower,
                    RobinBcCoefStrategy<NDIM>* const normal_bc_coef,
                    const bool homogeneous_bc,
                    NormalVelocityPatch& normal_velocity_patch,
                    const BoxArray<NDIM>& domain,
                    const Box<NDIM>& ghost_box,
                    const double fill_time)
{
    NormalVelocity u = { get_normal_velocity_stencil(i_face,
                                                     normal_velocity_patch.getBdryNormalAxis(),
                                                     bdry_is_lower,
                                                     normal_velocity_patch.getTangentialAxis(),
                                                     domain,
                                                     ghost_box),
                         0.0 };
    if (u.stencil.beyond_corner)
    {
        const std::optional<double> u_prescribed = get_prescribed_normal_velocity(
            i_face, u.stencil.adjacent_location_index, normal_bc_coef, normal_velocity_patch, fill_time);
        if (u_prescribed)
        {
            u.stencil.useAdjacentBoundaryValue();
            u.offset = homogeneous_bc ? 0.0 : 2.0 * *u_prescribed;
        }
    }
    return u;
} // get_normal_velocity

// Return the difference along the tangential axis of normal_velocity_patch of the normal velocity in the TRACTION
// condition for the tangential velocity component at the face i on the boundary normal to its normal axis, or NaN if
// the face i - 1 or i lies outside the ghost box of u_data. The patch of normal_velocity_patch must be the patch of
// u_data.
double
get_normal_velocity_difference(const SideData<NDIM, double>& u_data,
                               const hier::Index<NDIM>& i,
                               const bool bdry_is_lower,
                               RobinBcCoefStrategy<NDIM>* const normal_bc_coef,
                               const bool homogeneous_bc,
                               NormalVelocityPatch& normal_velocity_patch,
                               const BoxArray<NDIM>& domain,
                               const double fill_time)
{
    const unsigned int tangential_axis = normal_velocity_patch.getTangentialAxis();
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
        return get_normal_velocity(
            i_face, bdry_is_lower, normal_bc_coef, homogeneous_bc, normal_velocity_patch, domain, ghost_box, fill_time);
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
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

NormalVelocityPatch::NormalVelocityPatch(const Patch<NDIM>& patch,
                                         const unsigned int bdry_normal_axis,
                                         const unsigned int tangential_axis)
    : d_patch(copy_patch(patch)), d_bdry_normal_axis(bdry_normal_axis), d_tangential_axis(tangential_axis)
{
} // NormalVelocityPatch

const Box<NDIM>&
NormalVelocityPatch::getBox() const
{
    return d_patch->getBox();
} // getBox

unsigned int
NormalVelocityPatch::getBdryNormalAxis() const
{
    return d_bdry_normal_axis;
} // getBdryNormalAxis

unsigned int
NormalVelocityPatch::getTangentialAxis() const
{
    return d_tangential_axis;
} // getTangentialAxis

double
NormalVelocityPatch::getTangentialDx() const
{
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = d_patch->getPatchGeometry();
    return pgeom->getDx()[d_tangential_axis];
} // getTangentialDx

void
NormalVelocityPatch::setBcCoefs(RobinBcCoefStrategy<NDIM>* const normal_bc_coef,
                                Pointer<ArrayData<NDIM, double>>& acoef_data,
                                Pointer<ArrayData<NDIM, double>>& bcoef_data,
                                Pointer<ArrayData<NDIM, double>>& gcoef_data,
                                const BoundaryBox<NDIM>& bdry_box,
                                const double fill_time)
{
    const IBTK::detail::ScopedShiftedPatchGeometry shifted_geometry(
        *d_patch, get_half_cell_shift(*d_patch, d_bdry_normal_axis, d_tangential_axis));
    normal_bc_coef->setBcCoefs(
        acoef_data, bcoef_data, gcoef_data, Pointer<Variable<NDIM>>(), *d_patch, bdry_box, fill_time);
} // setBcCoefs

double
get_divergence_free_ghost_value(const SideData<NDIM, double>& u_data,
                                const hier::Index<NDIM>& i_g,
                                const unsigned int bdry_normal_axis,
                                const bool bdry_is_lower,
                                const double* const dx)
{
    const DivergenceFreeGhostFormula formula =
        get_divergence_free_ghost_formula(i_g, bdry_normal_axis, bdry_is_lower, dx);
    double div_u_g = 0.0;
    for (const DivergenceFreeGhostTerm& term : formula.terms)
    {
        const double u_upper = term.use_upper ? u_data(term.upper) : 0.0;
        const double u_lower = term.use_lower ? u_data(term.lower) : 0.0;
        div_u_g += (u_upper - u_lower) * term.dx_normal / term.dx_axis;
    }
    return formula.sign * div_u_g;
} // get_divergence_free_ghost_value

StokesBcCoefStrategy::SideStencil
get_divergence_free_ghost_value_stencil(const hier::Index<NDIM>& i_g,
                                        const unsigned int bdry_normal_axis,
                                        const bool bdry_is_lower,
                                        const double* const dx)
{
    const DivergenceFreeGhostFormula formula =
        get_divergence_free_ghost_formula(i_g, bdry_normal_axis, bdry_is_lower, dx);
    StokesBcCoefStrategy::SideStencil stencil;
    for (const DivergenceFreeGhostTerm& term : formula.terms)
    {
        const double weight = formula.sign * term.dx_normal / term.dx_axis;
        if (term.use_upper)
        {
            stencil.emplace_back(term.upper, weight);
        }
        if (term.use_lower)
        {
            stencil.emplace_back(term.lower, -weight);
        }
    }
    return stencil;
} // get_divergence_free_ghost_value_stencil

bool
has_normal_stress_term(const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
{
    for (const auto& u_bc_coef : u_bc_coefs)
    {
        const auto stokes_u_bc_coef = dynamic_cast<const StokesBcCoefStrategy*>(u_bc_coef);
        if (stokes_u_bc_coef && stokes_u_bc_coef->getTractionBcType() == TRACTION)
        {
            return true;
        }
    }
    return false;
} // has_normal_stress_term

bool
has_traction_face(const ArrayData<NDIM, double>& acoef_data, const ArrayData<NDIM, double>& bcoef_data)
{
    for (Box<NDIM>::Iterator it(acoef_data.getBox()); it; it++)
    {
        if (classify_velocity_bc(acoef_data(it(), 0), bcoef_data(it(), 0)) == VelocityBcType::TRACTION)
        {
            return true;
        }
    }
    return false;
} // has_traction_face

void
correct_traction_coefficients(ArrayData<NDIM, double>& acoef_data,
                              ArrayData<NDIM, double>& bcoef_data,
                              ArrayData<NDIM, double>& gcoef_data,
                              const SideData<NDIM, double>& u_data,
                              const Patch<NDIM>& patch,
                              const unsigned int bdry_normal_axis,
                              const bool bdry_is_lower,
                              const unsigned int tangential_axis,
                              RobinBcCoefStrategy<NDIM>* const normal_bc_coef,
                              const bool homogeneous_bc,
                              const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                              const std::function<double(const hier::Index<NDIM>&, double)>& get_mu,
                              const double fill_time)
{
    TBOX_ASSERT(patch_hierarchy);
    const Box<NDIM>& bc_coef_box = acoef_data.getBox();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const BoxArray<NDIM> domain =
        PhysicalBoundaryUtilities::getPeriodicPhysicalDomain(patch_hierarchy->getGridGeometry(), pgeom->getRatio());
    NormalVelocityPatch normal_velocity_patch(patch, bdry_normal_axis, tangential_axis);
    for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
    {
        const hier::Index<NDIM>& i = it();
        double& alpha = acoef_data(i, 0);
        double& beta = bcoef_data(i, 0);
        double& gamma = gcoef_data(i, 0);
        if (classify_velocity_bc(alpha, beta) == VelocityBcType::TRACTION)
        {
            // mu*(du_tan/dx_norm + du_norm/dx_tan) = g. Compute the tangential derivative of the normal component of
            // the velocity at the boundary and correct the boundary condition value.
            const double du_norm_dx_tan = get_normal_velocity_difference(u_data,
                                                                         i,
                                                                         bdry_is_lower,
                                                                         normal_bc_coef,
                                                                         homogeneous_bc,
                                                                         normal_velocity_patch,
                                                                         domain,
                                                                         fill_time) /
                                          dx[tangential_axis];
            alpha = 0.0;
            beta = 1.0;
            gamma = (bdry_is_lower ? -1.0 : +1.0) * (gamma / get_mu(i, gamma) - du_norm_dx_tan);
        }
    }
} // correct_traction_coefficients

void
set_component_bc_coefs(RobinBcCoefStrategy<NDIM>& bc_coef,
                       const unsigned int axis,
                       Patch<NDIM>& patch,
                       const BoundaryBox<NDIM>& bdry_box,
                       Pointer<ArrayData<NDIM, double>>& acoef_data,
                       Pointer<ArrayData<NDIM, double>>& bcoef_data,
                       Pointer<ArrayData<NDIM, double>>& gcoef_data,
                       const double fill_time)
{
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    std::array<double, NDIM> shift{};
    shift[axis] = -0.5 * pgeom->getDx()[axis];
    const IBTK::detail::ScopedShiftedPatchGeometry shifted_geometry(patch, shift);
    bc_coef.setBcCoefs(acoef_data, bcoef_data, gcoef_data, nullptr, patch, bdry_box, fill_time);
    return;
} // set_component_bc_coefs

bool
traction_data_depend_on_velocity(RobinBcCoefStrategy<NDIM>* const comp_bc_coef,
                                 const unsigned int comp_idx,
                                 const Patch<NDIM>& patch,
                                 const BoundaryBox<NDIM>& bdry_box,
                                 const double fill_time)
{
    if (static_cast<int>(comp_idx) == bdry_box.getLocationIndex() / 2)
    {
        return false;
    }

    // Evaluate the coefficients at the faces where the fill sets the component: the cells of the boundary box and one
    // more face along the component axis, with the patch geometry shifted by half a cell along that axis.
    Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(bdry_box);
    bc_coef_box.upper()(comp_idx) += 1;
    Pointer<ArrayData<NDIM, double>> acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
    Pointer<ArrayData<NDIM, double>> bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
    Pointer<ArrayData<NDIM, double>> gcoef_data;
    std::array<double, NDIM> shift{};
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    shift[comp_idx] = -0.5 * pgeom->getDx()[comp_idx];
    // Shifting the geometry of a copy of the patch leaves the patch geometry unchanged.
    Pointer<Patch<NDIM>> shifted_patch = copy_patch(patch);
    const IBTK::detail::ScopedShiftedPatchGeometry shifted_geometry(*shifted_patch, shift);
    comp_bc_coef->setBcCoefs(
        acoef_data, bcoef_data, gcoef_data, Pointer<Variable<NDIM>>(), *shifted_patch, bdry_box, fill_time);
    for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
    {
        if (classify_velocity_bc((*acoef_data)(it(), 0), (*bcoef_data)(it(), 0)) == VelocityBcType::TRACTION)
        {
            return true;
        }
    }
    return false;
} // traction_data_depend_on_velocity

StokesBcCoefStrategy::SideStencil
get_traction_stencil(const hier::Index<NDIM>& i,
                     const bool bdry_is_lower,
                     RobinBcCoefStrategy<NDIM>* const normal_bc_coef,
                     NormalVelocityPatch& normal_velocity_patch,
                     const BoxArray<NDIM>& domain,
                     const Box<NDIM>& ghost_box,
                     const double fill_time)
{
    const unsigned int tangential_axis = normal_velocity_patch.getTangentialAxis();
    const double derivative_scale = (bdry_is_lower ? -1.0 : +1.0) / normal_velocity_patch.getTangentialDx();
    hier::Index<NDIM> i_lower(i);
    i_lower(tangential_axis) -= 1;
    const auto get_stencil = [&](const hier::Index<NDIM>& i_face)
    {
        return get_normal_velocity(i_face,
                                   bdry_is_lower,
                                   normal_bc_coef,
                                   /*homogeneous_bc*/ true,
                                   normal_velocity_patch,
                                   domain,
                                   ghost_box,
                                   fill_time)
            .stencil;
    };
    const NormalVelocityStencil lower_stencil = get_stencil(i_lower);
    const NormalVelocityStencil upper_stencil = get_stencil(i);
    StokesBcCoefStrategy::SideStencil stencil;
    for (int k = 0; k < 2; ++k)
    {
        if (upper_stencil.weight[k] != 0.0)
        {
            stencil.emplace_back(upper_stencil.idx[k], -upper_stencil.weight[k] * derivative_scale);
        }
        if (lower_stencil.weight[k] != 0.0)
        {
            stencil.emplace_back(lower_stencil.idx[k], +lower_stencil.weight[k] * derivative_scale);
        }
    }
    return stencil;
} // get_traction_stencil

StokesBcCoefStrategy::SideStencil
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
    SideIndex<NDIM> i_s_inside(i, bdry_normal_axis, SideIndex<NDIM>::Lower);
    i_s_inside(bdry_normal_axis) += bdry_is_lower ? 1 : -1;
    StokesBcCoefStrategy::SideStencil stencil;
    stencil.emplace_back(i_s_inside, 1.0);
    for (const auto& entry : get_divergence_free_ghost_value_stencil(i_g, bdry_normal_axis, bdry_is_lower, dx))
    {
        stencil.emplace_back(entry.first, -entry.second);
    }
    return stencil;
} // get_normal_stress_stencil

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace detail
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
