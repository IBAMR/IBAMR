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
#include <ibtk/PhysicalBoundaryUtilities.h>

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
#include <SideData.h>
#include <SideIndex.h>

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
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
// A patch geometry shifted by a given displacement, so that boundary coefficients are evaluated at the locations of
// another component of the velocity. The geometry is built once, from the geometry that the patch has at that time.
// Constructing a ShiftedPatchGeometry::Scope replaces the geometry of the patch by the shifted geometry, and destroying
// the scope restores the geometry that the patch had, on every exit path. The patch is not otherwise modified.
class ShiftedPatchGeometry
{
public:
    ShiftedPatchGeometry(const Patch<NDIM>& patch, const std::array<double, NDIM>& shift)
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
    }

    class Scope
    {
    public:
        // The patch passed to a boundary coefficient object is a reference to a patch that the caller owns and can
        // modify, and its geometry is replaced only until the scope ends.
        Scope(const Patch<NDIM>& patch, const ShiftedPatchGeometry& shifted)
            : d_patch(const_cast<Patch<NDIM>&>(patch)), d_original(patch.getPatchGeometry())
        {
            d_patch.setPatchGeometry(shifted.d_geometry);
        }

        ~Scope()
        {
            d_patch.setPatchGeometry(d_original);
        }

        Scope(const Scope&) = delete;
        Scope& operator=(const Scope&) = delete;

    private:
        Patch<NDIM>& d_patch;
        Pointer<PatchGeometry<NDIM>> d_original;
    };

private:
    Pointer<PatchGeometry<NDIM>> d_geometry;
};

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
    if (!IBTK::rel_equal_eps(alpha, 1.0))
    {
        return false;
    }
    u_prescribed = (*gcoef_data)(i_coef, 0) / alpha;
    return true;
} // get_prescribed_normal_velocity

// Return the physical domain refined by ratio, together with its periodic
// images, so that a cell across a periodic boundary is not mistaken for a cell
// outside the domain.
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

// The normal velocity used by the TRACTION condition at a boundary face, as
// weight[0]*u(idx[0]) + weight[1]*u(idx[1]) + offset.
struct NormalVelocityStencil
{
    std::array<SideIndex<NDIM>, 2> idx;
    std::array<double, 2> weight;
    double offset;
};

// Return the stencil of the normal velocity from which the TRACTION condition
// takes the normal velocity at the boundary face with index i_face, where the
// boundary is normal to bdry_normal_axis and the condition differences the
// normal velocity along tangential_axis. A face is beyond a corner if the cell
// adjacent to it on the interior side of the boundary lies outside the physical
// domain; this is determined from the domain, not from the patch. The value at
// such a face is not available. If the adjacent boundary prescribes the normal
// velocity, the value is 2*u_b - u, where u is the nearest face on the boundary
// and u_b is the prescribed velocity, which is zero in a homogeneous condition;
// otherwise it is extrapolated linearly from the two nearest faces on the
// boundary, or is the value of the nearest face if the boundary segment has
// only one. A face on the boundary that lies outside ghost_box along
// tangential_axis is extrapolated linearly from the two nearest faces in
// ghost_box, which makes the difference there the one-sided difference inside
// ghost_box; if ghost_box holds a single face along tangential_axis, the
// nearest face is used.
NormalVelocityStencil
get_normal_velocity_stencil(hier::Index<NDIM> i_face,
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
    const int j = i_face(tangential_axis);
    if (is_boundary_face(i_face, bdry_normal_axis, bdry_is_lower, domain))
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
        TBOX_ERROR("INSStaggeredVelocityBcCoef: the normal velocity face "
                   << i_face << " is beyond a corner of the physical boundary, but neither adjacent face along axis "
                   << tangential_axis << " lies on the boundary.\n");
    }
    hier::Index<NDIM> i_near = i_face, i_next = i_face;
    i_near(tangential_axis) = j + step;
    i_next(tangential_axis) = j + 2 * step;
#if !defined(NDEBUG)
    TBOX_ASSERT(ghost_box.contains(i_near));
#endif
    const SideIndex<NDIM> near_idx(i_near, bdry_normal_axis, SideIndex<NDIM>::Lower);
    const SideIndex<NDIM> next_idx(i_next, bdry_normal_axis, SideIndex<NDIM>::Lower);

    // The adjacent boundary is the one at which the boundary ends.
    const unsigned int adj_location_index = 2 * tangential_axis + (step < 0 ? 1 : 0);
    double u_prescribed = 0.0;
    if (get_prescribed_normal_velocity(u_prescribed,
                                       i_face,
                                       bdry_normal_axis,
                                       tangential_axis,
                                       adj_location_index,
                                       normal_bc_coef,
                                       patch,
                                       normal_geometry,
                                       fill_time))
    {
        return { { near_idx, near_idx }, { -1.0, 0.0 }, homogeneous_bc ? 0.0 : 2.0 * u_prescribed };
    }

    // Linear extrapolation requires a second face on the boundary; a boundary segment that is one cell long has none.
    if (is_boundary_face(i_next, bdry_normal_axis, bdry_is_lower, domain))
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(ghost_box.contains(i_next));
#endif
        return { { near_idx, next_idx }, { 2.0, -1.0 }, 0.0 };
    }
    return { { near_idx, near_idx }, { 1.0, 0.0 }, 0.0 };
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
    const BoxArray<NDIM>* domain = nullptr;
    std::unique_ptr<ShiftedPatchGeometry> normal_geometry;
    const auto get_stencil = [&](const hier::Index<NDIM>& i_face)
    {
        return get_normal_velocity_stencil(i_face,
                                           bdry_normal_axis,
                                           is_lower,
                                           d_comp_idx,
                                           d_bc_coefs[bdry_normal_axis],
                                           d_homogeneous_bc,
                                           patch,
                                           normal_geometry,
                                           *domain,
                                           ghost_box,
                                           fill_time);
    };
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
                    if (!domain)
                    {
                        domain = &getPhysicalDomain(patch);
                    }
                    hier::Index<NDIM> i_lower(i);
                    i_lower(d_comp_idx) -= 1;
                    const NormalVelocityStencil u_lower = get_stencil(i_lower);
                    const NormalVelocityStencil u_upper = get_stencil(i);
                    double du_norm = u_upper.offset - u_lower.offset;
                    for (int k = 0; k < 2; ++k)
                    {
                        du_norm += u_upper.weight[k] * (*u_target_data)(u_upper.idx[k]) -
                                   u_lower.weight[k] * (*u_target_data)(u_lower.idx[k]);
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

const BoxArray<NDIM>&
INSStaggeredVelocityBcCoef::getPhysicalDomain(const Patch<NDIM>& patch) const
{
    const Pointer<PatchHierarchy<NDIM>> hierarchy = d_fluid_solver->getPatchHierarchy();
    if (!hierarchy)
    {
        TBOX_ERROR(
            "INSStaggeredVelocityBcCoef::getPhysicalDomain(): the fluid solver has no patch hierarchy; a "
            "TRACTION boundary condition requires the hierarchy to locate the corners of the physical "
            "boundary.\n");
    }
    const IntVector<NDIM>& ratio = patch.getPatchGeometry()->getRatio();
    std::array<int, NDIM> key;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        key[d] = ratio(d);
    }
    auto it = d_physical_domain.find(key);
    if (it == d_physical_domain.end())
    {
        // The physical domain and its periodic shifts do not change during a run, so the refined domain is computed
        // once for each refinement ratio.
        it = d_physical_domain.emplace(key, compute_physical_domain(hierarchy->getGridGeometry(), ratio)).first;
    }
    return it->second;
} // getPhysicalDomain

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
