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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_ins_staggered_traction_stencil
#define included_IBAMR_ins_staggered_traction_stencil

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <tbox/Pointer.h>

#include <ArrayData.h>
#include <BoundaryBox.h>
#include <Box.h>
#include <BoxArray.h>
#include <Index.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchGeometry.h>
#include <PatchHierarchy.h>
#include <RobinBcCoefStrategy.h>
#include <SideData.h>
#include <SideIndex.h>

#include <array>
#include <map>
#include <memory>
#include <utility>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

// Helpers for the TRACTION velocity boundary condition on the staggered grid. This header is not installed.
namespace IBAMR
{
namespace traction_stencil
{
/*!
 * The normal velocity used by the TRACTION condition at a face on a boundary: the weighted sum
 * weight[0]*u(idx[0]) + weight[1]*u(idx[1]) of the normal velocity at faces on that boundary.
 *
 * At a face beyond a corner (beyond_corner is true), the stencil extrapolates linearly from the nearest faces. If the
 * adjacent boundary, with location index adjacent_location_index, prescribes the normal velocity u_b, the value is
 * 2*u_b minus the value at the nearest face: the caller calls useAdjacentBoundaryValue(), which sets
 * uses_adjacent_boundary_value, and adds 2*u_b.
 */
struct NormalVelocityStencil
{
    std::array<SAMRAI::pdat::SideIndex<NDIM>, 2> idx;
    std::array<double, 2> weight;

    bool beyond_corner = false;
    unsigned int adjacent_location_index = 0;
    bool uses_adjacent_boundary_value = false;

    /*!
     * Replace the weighted sum by minus the value at the nearest face. The face must be beyond a corner.
     */
    void useAdjacentBoundaryValue();
};

/*!
 * A patch geometry shifted by a displacement, so that boundary condition objects evaluate their coefficients at the
 * locations of another velocity component.
 */
class ShiftedPatchGeometry
{
public:
    ShiftedPatchGeometry(const SAMRAI::hier::Patch<NDIM>& patch, const std::array<double, NDIM>& shift);

    /*!
     * Replace the geometry of a patch by a ShiftedPatchGeometry until the scope ends.
     */
    class Scope
    {
    public:
        Scope(const SAMRAI::hier::Patch<NDIM>& patch, const ShiftedPatchGeometry& shifted);

        ~Scope();

        Scope(const Scope&) = delete;
        Scope& operator=(const Scope&) = delete;

    private:
        SAMRAI::hier::Patch<NDIM>& d_patch;
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchGeometry<NDIM>> d_original;
    };

private:
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchGeometry<NDIM>> d_geometry;
};

/*!
 * The physical domain with its periodic images, by refinement ratio, as filled by get_physical_domain().
 */
using PhysicalDomainCache = std::map<std::array<int, NDIM>, SAMRAI::hier::BoxArray<NDIM>>;

/*!
 * Return the physical domain of hierarchy refined by ratio, with its periodic images, computed once for each ratio and
 * stored in cache.
 */
const SAMRAI::hier::BoxArray<NDIM>&
get_physical_domain(PhysicalDomainCache& cache,
                    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
                    const SAMRAI::hier::IntVector<NDIM>& ratio);

/*!
 * Return the stencil of the normal velocity at the face i_face beyond a corner of the boundary normal to
 * bdry_normal_axis, where the boundary segment ends along tangential_axis. The nearest face on the boundary is at
 * i_face + step along tangential_axis; the second nearest is used if can_extrapolate is true.
 */
NormalVelocityStencil get_corner_normal_velocity_stencil(SAMRAI::hier::Index<NDIM> i_face,
                                                         unsigned int bdry_normal_axis,
                                                         unsigned int tangential_axis,
                                                         int step,
                                                         bool can_extrapolate);

/*!
 * Return the difference along tangential_axis of the normal velocity that the TRACTION condition for the velocity
 * component tangential_axis uses at the location i on the boundary normal to bdry_normal_axis: the value at the face i
 * minus the value at the next face below it, taken from u_data with the stencils for domain. At a face beyond a
 * corner, normal_bc_coef decides whether the adjacent boundary prescribes the normal velocity (zero if homogeneous_bc).
 * The geometry of patch must be centered on the tangential component; normal_geometry caches the geometry centered on
 * the normal component.
 */
double get_normal_velocity_difference(const SAMRAI::pdat::SideData<NDIM, double>& u_data,
                                      const SAMRAI::hier::Index<NDIM>& i,
                                      unsigned int bdry_normal_axis,
                                      bool bdry_is_lower,
                                      unsigned int tangential_axis,
                                      SAMRAI::solv::RobinBcCoefStrategy<NDIM>* normal_bc_coef,
                                      bool homogeneous_bc,
                                      const SAMRAI::hier::Patch<NDIM>& patch,
                                      std::unique_ptr<ShiftedPatchGeometry>& normal_geometry,
                                      const SAMRAI::hier::BoxArray<NDIM>& domain,
                                      double fill_time);

/*!
 * Return the dependence on the normal velocity of the inhomogeneous Robin coefficient (+/-)(g/mu - d/h) that the
 * TRACTION condition for the velocity component tangential_axis sets at the location i on the boundary normal to
 * bdry_normal_axis, with d = get_normal_velocity_difference() at i. The result lists the normal velocity values that
 * the coefficient depends on, each with its derivative; a value may appear more than once, and the derivatives add.
 * The arguments are those of get_normal_velocity_difference(), except that ghost_box replaces u_data and there is no
 * homogeneous_bc.
 */
std::vector<std::pair<SAMRAI::pdat::SideIndex<NDIM>, double>>
get_traction_stencil(const SAMRAI::hier::Index<NDIM>& i,
                     unsigned int bdry_normal_axis,
                     bool bdry_is_lower,
                     unsigned int tangential_axis,
                     SAMRAI::solv::RobinBcCoefStrategy<NDIM>* normal_bc_coef,
                     const SAMRAI::hier::Patch<NDIM>& patch,
                     std::unique_ptr<ShiftedPatchGeometry>& normal_geometry,
                     const SAMRAI::hier::BoxArray<NDIM>& domain,
                     const SAMRAI::hier::Box<NDIM>& ghost_box,
                     double fill_time);

/*!
 * Apply the transpose of the dependence on the normal velocity of the TRACTION condition for the velocity component
 * tangential_axis on the boundary box bdry_box. The condition sets the inhomogeneous Robin coefficient at a location i
 * to (+/-)(g/mu - d/h), with d = get_normal_velocity_difference() at i, h the grid spacing along tangential_axis, and
 * the sign negative on a lower boundary. For each location i of gcoef_data at which tangential_bc_coef prescribes the
 * traction, the function adds gcoef_data(i) times the derivative of the coefficient with respect to each normal
 * velocity value to that value of u_data.
 */
void accumulate_from_traction_bc_coefs(SAMRAI::pdat::SideData<NDIM, double>& u_data,
                                       const SAMRAI::pdat::ArrayData<NDIM, double>& gcoef_data,
                                       unsigned int tangential_axis,
                                       SAMRAI::solv::RobinBcCoefStrategy<NDIM>* tangential_bc_coef,
                                       SAMRAI::solv::RobinBcCoefStrategy<NDIM>* normal_bc_coef,
                                       bool homogeneous_bc,
                                       const SAMRAI::hier::Patch<NDIM>& patch,
                                       const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
                                       PhysicalDomainCache& domain_cache,
                                       const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
                                       double fill_time);
} // namespace traction_stencil
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_ins_staggered_traction_stencil
