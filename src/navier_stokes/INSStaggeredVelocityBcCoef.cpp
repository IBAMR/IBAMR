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

#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/INSStaggeredVelocityBcCoef.h>
#include <ibamr/StokesBcCoefStrategy.h>
#include <ibamr/StokesSpecifications.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/ExtendedRobinBcCoefStrategy.h>

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
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>

#include <algorithm>
#include <array>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

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
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Return the physical domain refined to the level of patch, together with its
// periodic images, so that a cell across a periodic boundary is not mistaken
// for a cell outside the domain.
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

// The normal velocity used by the TRACTION condition at a boundary face, as
// weight[0]*u(idx[0]) + weight[1]*u(idx[1]).
struct NormalVelocityStencil
{
    std::array<SideIndex<NDIM>, 2> idx;
    std::array<double, 2> weight;
};

// Return the stencil of the normal velocity from which the TRACTION condition
// takes the normal velocity at the boundary face with index i_face, where the
// boundary is normal to bdry_normal_axis and the condition differences the
// normal velocity along tangential_axis. A face is beyond a corner if the cell
// adjacent to it on the interior side of the boundary lies outside the physical
// domain; this is determined from the domain, not from the patch. The value at
// such a face is not available, so the nearest face on the boundary is used,
// which makes the difference zero across the corner. A face on the boundary
// that lies outside ghost_box along tangential_axis is extrapolated linearly
// from the two nearest faces in ghost_box, which makes the difference there the
// one-sided difference inside ghost_box; if ghost_box holds a single face along
// tangential_axis, the nearest face is used.
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
        const int j_lower = ghost_box.lower(tangential_axis);
        const int j_upper = ghost_box.upper(tangential_axis);
        NormalVelocityStencil stencil;
        i_face(tangential_axis) = std::min(std::max(j, j_lower), j_upper);
        stencil.idx[0] = SideIndex<NDIM>(i_face, bdry_normal_axis, SideIndex<NDIM>::Lower);
        stencil.idx[1] = stencil.idx[0];
        stencil.weight = { 1.0, 0.0 };
        if ((j < j_lower || j > j_upper) && j_lower < j_upper)
        {
            i_face(tangential_axis) += (j < j_lower ? +1 : -1);
            stencil.idx[1] = SideIndex<NDIM>(i_face, bdry_normal_axis, SideIndex<NDIM>::Lower);
            stencil.weight = { 2.0, -1.0 };
        }
        return stencil;
    }

    // The nearest face on the boundary: below the face if the boundary ends above
    // it, above the face if the boundary ends below it.
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
        TBOX_ERROR("INSStaggeredVelocityBcCoef: the normal velocity face "
                   << i_face << " is beyond a corner of the physical boundary, but neither adjacent face along axis "
                   << tangential_axis << " lies on the boundary.\n");
    }
    hier::Index<NDIM> i_near = i_face;
    i_near(tangential_axis) = j + step;
#if !defined(NDEBUG)
    TBOX_ASSERT(ghost_box.contains(i_near));
#endif
    NormalVelocityStencil stencil;
    stencil.idx[0] = SideIndex<NDIM>(i_near, bdry_normal_axis, SideIndex<NDIM>::Lower);
    stencil.idx[1] = stencil.idx[0];
    stencil.weight = { 1.0, 0.0 };
    return stencil;
} // get_normal_velocity_stencil
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredVelocityBcCoef::INSStaggeredVelocityBcCoef(const unsigned int comp_idx,
                                                       const INSStaggeredHierarchyIntegrator* fluid_solver,
                                                       const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                                       const TractionBcType traction_bc_type,
                                                       const bool homogeneous_bc)
    : d_comp_idx(comp_idx), d_fluid_solver(fluid_solver), d_bc_coefs(NDIM, nullptr)
{
    setStokesSpecifications(d_fluid_solver->getStokesSpecifications());
    setPhysicalBcCoefs(bc_coefs);
    setTractionBcType(traction_bc_type);
    setHomogeneousBc(homogeneous_bc);
    return;
} // INSStaggeredVelocityBcCoef

void
INSStaggeredVelocityBcCoef::setStokesSpecifications(const StokesSpecifications* problem_coefs)
{
    StokesBcCoefStrategy::setStokesSpecifications(problem_coefs);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setStokesSpecifications(problem_coefs);
    }
    return;
} // setStokesSpecifications

void
INSStaggeredVelocityBcCoef::setTargetVelocityPatchDataIndex(int u_target_data_idx)
{
    StokesBcCoefStrategy::setTargetVelocityPatchDataIndex(u_target_data_idx);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setTargetVelocityPatchDataIndex(u_target_data_idx);
    }
    return;
} // setTargetVelocityPatchDataIndex

void
INSStaggeredVelocityBcCoef::clearTargetVelocityPatchDataIndex()
{
    StokesBcCoefStrategy::clearTargetVelocityPatchDataIndex();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->clearTargetVelocityPatchDataIndex();
    }
    return;
} // clearTargetVelocityPatchDataIndex

void
INSStaggeredVelocityBcCoef::setTargetPressurePatchDataIndex(int p_target_data_idx)
{
    StokesBcCoefStrategy::setTargetPressurePatchDataIndex(p_target_data_idx);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setTargetPressurePatchDataIndex(p_target_data_idx);
    }
    return;
} // setTargetPressurePatchDataIndex

void
INSStaggeredVelocityBcCoef::clearTargetPressurePatchDataIndex()
{
    StokesBcCoefStrategy::clearTargetPressurePatchDataIndex();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->clearTargetPressurePatchDataIndex();
    }
    return;
} // clearTargetPressurePatchDataIndex

void
INSStaggeredVelocityBcCoef::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_coefs.size() == NDIM);
#endif
    d_bc_coefs = bc_coefs;
    return;
} // setPhysicalBcCoefs

void
INSStaggeredVelocityBcCoef::setSolutionTime(const double /*solution_time*/)
{
    // intentionally blank
    return;
} // setSolutionTime

void
INSStaggeredVelocityBcCoef::setTimeInterval(const double /*current_time*/, const double /*new_time*/)
{
    // intentionally blank
    return;
} // setTimeInterval

void
INSStaggeredVelocityBcCoef::setTargetPatchDataIndex(int target_idx)
{
    StokesBcCoefStrategy::setTargetPatchDataIndex(target_idx);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setTargetPatchDataIndex(target_idx);
    }
    return;
} // setTargetPatchDataIndex

void
INSStaggeredVelocityBcCoef::clearTargetPatchDataIndex()
{
    StokesBcCoefStrategy::clearTargetPatchDataIndex();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->clearTargetPatchDataIndex();
    }
    return;
} // clearTargetPatchDataIndex

void
INSStaggeredVelocityBcCoef::setHomogeneousBc(bool homogeneous_bc)
{
    ExtendedRobinBcCoefStrategy::setHomogeneousBc(homogeneous_bc);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setHomogeneousBc(homogeneous_bc);
    }
    return;
} // setHomogeneousBc

void
INSStaggeredVelocityBcCoef::accumulateFromBcCoefs(const ArrayData<NDIM, double>& gcoef_data,
                                                  const Patch<NDIM>& patch,
                                                  const BoundaryBox<NDIM>& bdry_box,
                                                  const double fill_time) const
{
    // Only a TRACTION condition on a tangential component depends on the
    // velocity.
    const unsigned int location_index = bdry_box.getLocationIndex();
    const unsigned int bdry_normal_axis = location_index / 2;
    if (d_traction_bc_type != TRACTION || d_comp_idx == bdry_normal_axis)
    {
        return;
    }

    // Determine where the physical boundary conditions prescribe the traction.
    const Box<NDIM>& bc_coef_box = gcoef_data.getBox();
    Pointer<ArrayData<NDIM, double>> acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
    Pointer<ArrayData<NDIM, double>> bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
    Pointer<ArrayData<NDIM, double>> gcoef_unused;
    d_bc_coefs[d_comp_idx]->setBcCoefs(
        acoef_data, bcoef_data, gcoef_unused, Pointer<Variable<NDIM>>(), patch, bdry_box, fill_time);

    // setBcCoefs() sets gamma = sgn*(g/mu - (u_upper - u_lower)/dx_tan) using
    // the normal velocity at the boundary. Accumulate the transpose of the
    // velocity term into the target velocity, including its ghost cells.
    Pointer<SideData<NDIM, double>> u_target_data;
    if (d_u_target_data_idx >= 0)
    {
        u_target_data = patch.getPatchData(d_u_target_data_idx);
    }
    else if (d_target_data_idx >= 0)
    {
        u_target_data = patch.getPatchData(d_target_data_idx);
    }
#if !defined(NDEBUG)
    TBOX_ASSERT(u_target_data);
#endif
    const Box<NDIM>& ghost_box = u_target_data->getGhostBox();
    const Box<NDIM> target_side_box = SideGeometry<NDIM>::toSideBox(ghost_box, bdry_normal_axis);
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const bool is_lower = location_index % 2 == 0;
    const double sgn = is_lower ? -1.0 : +1.0;
    BoxArray<NDIM> domain;
    bool have_domain = false;
    for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
    {
        const hier::Index<NDIM>& i = it();
        if (!IBTK::rel_equal_eps((*bcoef_data)(i, 0), 1.0))
        {
            continue;
        }
        if (!have_domain)
        {
            domain = get_physical_domain(
                d_fluid_solver->getPatchHierarchy(), patch, "INSStaggeredVelocityBcCoef::accumulateFromBcCoefs()");
            have_domain = true;
        }
        hier::Index<NDIM> i_lower(i);
        i_lower(d_comp_idx) -= 1;
        const NormalVelocityStencil lower_stencil =
            get_normal_velocity_stencil(i_lower, bdry_normal_axis, is_lower, d_comp_idx, domain, ghost_box);
        const NormalVelocityStencil upper_stencil =
            get_normal_velocity_stencil(i, bdry_normal_axis, is_lower, d_comp_idx, domain, ghost_box);
        const double du_transpose = sgn * gcoef_data(i, 0) / dx[d_comp_idx];
        for (unsigned int k = 0; k < 2; ++k)
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(target_side_box.contains(upper_stencil.idx[k]));
            TBOX_ASSERT(target_side_box.contains(lower_stencil.idx[k]));
#endif
            (*u_target_data)(upper_stencil.idx[k]) -= upper_stencil.weight[k] * du_transpose;
            (*u_target_data)(lower_stencil.idx[k]) += lower_stencil.weight[k] * du_transpose;
        }
    }
    return;
} // accumulateFromBcCoefs

void
INSStaggeredVelocityBcCoef::setBcCoefs(Pointer<ArrayData<NDIM, double>>& acoef_data,
                                       Pointer<ArrayData<NDIM, double>>& bcoef_data,
                                       Pointer<ArrayData<NDIM, double>>& gcoef_data,
                                       const Pointer<Variable<NDIM>>& variable,
                                       const Patch<NDIM>& patch,
                                       const BoundaryBox<NDIM>& bdry_box,
                                       double fill_time) const
{
#if !defined(NDEBUG)
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d_bc_coefs[d]);
    }
#endif
    // Set the unmodified velocity bc coefs.
    d_bc_coefs[d_comp_idx]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);

    // We do not make any further modifications to the values of acoef_data and
    // bcoef_data beyond this point.
    if (!gcoef_data) return;
#if !defined(NDEBUG)
    TBOX_ASSERT(acoef_data);
    TBOX_ASSERT(bcoef_data);
#endif

    // Ensure homogeneous boundary conditions are enforced.
    if (d_homogeneous_bc) gcoef_data->fillAll(0.0);

    // Get the target velocity data.
    Pointer<SideData<NDIM, double>> u_target_data;
    if (d_u_target_data_idx >= 0)
        u_target_data = patch.getPatchData(d_u_target_data_idx);
    else if (d_target_data_idx >= 0)
        u_target_data = patch.getPatchData(d_target_data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(u_target_data);
#endif

    // Where appropriate, update boundary condition coefficients.
    //
    // Dirichlet boundary conditions are not modified.
    //
    // Neumann boundary conditions on the normal component of the velocity are
    // interpreted as "open" boundary conditions, and we set du/dn = 0.
    //
    // Neumann boundary conditions on the tangential component of the velocity
    // are interpreted as traction (stress) boundary conditions, and we update
    // the boundary condition coefficients accordingly.
    const unsigned int location_index = bdry_box.getLocationIndex();
    const unsigned int bdry_normal_axis = location_index / 2;
    const bool is_lower = location_index % 2 == 0;
    const Box<NDIM>& bc_coef_box = acoef_data->getBox();
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_coef_box == acoef_data->getBox());
    TBOX_ASSERT(bc_coef_box == bcoef_data->getBox());
    TBOX_ASSERT(bc_coef_box == gcoef_data->getBox());
#endif
    const Box<NDIM>& ghost_box = u_target_data->getGhostBox();
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double mu = d_problem_coefs->getMu();
    BoxArray<NDIM> domain;
    bool have_domain = false;
    for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
    {
        const hier::Index<NDIM>& i = it();
        double& alpha = (*acoef_data)(i, 0);
        double& beta = (*bcoef_data)(i, 0);
        double& gamma = (*gcoef_data)(i, 0);
        const bool velocity_bc = IBTK::rel_equal_eps(alpha, 1.0);
        const bool traction_bc = IBTK::rel_equal_eps(beta, 1.0);
#if !defined(NDEBUG)
        TBOX_ASSERT((velocity_bc || traction_bc) && !(velocity_bc && traction_bc));
#endif
        if (velocity_bc)
        {
            alpha = 1.0;
            beta = 0.0;
        }
        else if (traction_bc)
        {
            if (d_comp_idx == bdry_normal_axis)
            {
                // Set du/dn = 0.
                //
                // NOTE: We would prefer to determine the ghost cell value of
                // the normal velocity so that div u = 0 in the ghost cell.
                // This could be done here, but it is more convenient to do so
                // as a post-processing step after the tangential velocity ghost
                // cell values have all been set.
                alpha = 0.0;
                beta = 1.0;
                gamma = 0.0;
            }
            else
            {
                switch (d_traction_bc_type)
                {
                case TRACTION: // mu*(du_tan/dx_norm + du_norm/dx_tan) = g.
                {
                    // Compute the tangential derivative of the normal
                    // component of the velocity at the boundary.
                    if (!have_domain)
                    {
                        domain = get_physical_domain(
                            d_fluid_solver->getPatchHierarchy(), patch, "INSStaggeredVelocityBcCoef::setBcCoefs()");
                        have_domain = true;
                    }
                    hier::Index<NDIM> i_lower(i);
                    i_lower(d_comp_idx) -= 1;
                    const NormalVelocityStencil lower_stencil =
                        get_normal_velocity_stencil(i_lower, bdry_normal_axis, is_lower, d_comp_idx, domain, ghost_box);
                    const NormalVelocityStencil upper_stencil =
                        get_normal_velocity_stencil(i, bdry_normal_axis, is_lower, d_comp_idx, domain, ghost_box);
                    double du_norm = 0.0;
                    for (unsigned int k = 0; k < 2; ++k)
                    {
                        du_norm += upper_stencil.weight[k] * (*u_target_data)(upper_stencil.idx[k]) -
                                   lower_stencil.weight[k] * (*u_target_data)(lower_stencil.idx[k]);
                    }
                    const double du_norm_dx_tan = du_norm / dx[d_comp_idx];

                    // Correct the boundary condition value.
                    alpha = 0.0;
                    beta = 1.0;
                    gamma = (is_lower ? -1.0 : +1.0) * (gamma / mu - du_norm_dx_tan);
                    break;
                }
                case PSEUDO_TRACTION: // mu*du_tan/dx_norm = g.
                {
                    alpha = 0.0;
                    beta = 1.0;
                    gamma = (is_lower ? -1.0 : +1.0) * (gamma / mu);
                    break;
                }
                default:
                {
                    TBOX_ERROR(
                        "INSStaggeredVelocityBcCoef::setBcCoefs(): unrecognized or "
                        "unsupported "
                        "traction boundary condition type: "
                        << enum_to_string<TractionBcType>(d_traction_bc_type) << "\n");
                }
                }
            }
        }
        else
        {
            TBOX_ERROR("this statement should not be reached!\n");
        }
    }
    return;
} // setBcCoefs

IntVector<NDIM>
INSStaggeredVelocityBcCoef::numberOfExtensionsFillable() const
{
#if !defined(NDEBUG)
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d_bc_coefs[d]);
    }
#endif
    IntVector<NDIM> ret_val(std::numeric_limits<int>::max());
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        ret_val = IntVector<NDIM>::min(ret_val, d_bc_coefs[d]->numberOfExtensionsFillable());
    }
    return ret_val;
} // numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
