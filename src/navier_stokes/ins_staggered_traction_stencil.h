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

#include <ibtk/ShiftedPatchGeometry.h>

#include <tbox/Pointer.h>

#include <ArrayData.h>
#include <BoundaryBox.h>
#include <BoxArray.h>
#include <Index.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <RobinBcCoefStrategy.h>
#include <SideData.h>
#include <SideIndex.h>

#include <array>
#include <utility>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

// Helpers for the TRACTION velocity boundary condition on the staggered grid. This header is not installed.
namespace IBAMR
{
namespace traction_stencil
{
/*!
 * The normal velocity that the TRACTION condition uses at a face on a boundary: weight[0]*u(idx[0]) +
 * weight[1]*u(idx[1]). At a face beyond a corner, the stencil extrapolates linearly from the nearest faces; if the
 * adjacent boundary (adjacent_location_index) prescribes the normal velocity u_b, the caller calls
 * useAdjacentBoundaryValue() and adds 2*u_b.
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
 * Return the physical domain of hierarchy refined by ratio, with its periodic images.
 */
SAMRAI::hier::BoxArray<NDIM>
get_physical_domain(const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
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
 * Return the geometry of patch centered on the velocity component bdry_normal_axis, given that the geometry of patch is
 * centered on the velocity component tangential_axis.
 */
IBTK::ShiftedPatchGeometry get_normal_velocity_geometry(const SAMRAI::hier::Patch<NDIM>& patch,
                                                        unsigned int bdry_normal_axis,
                                                        unsigned int tangential_axis);

/*!
 * Return the difference along tangential_axis of the normal velocity that the TRACTION condition for velocity component
 * tangential_axis uses at i on the boundary normal to bdry_normal_axis. Beyond a corner, normal_bc_coef decides whether
 * the adjacent boundary prescribes the normal velocity (zero if homogeneous_bc). The geometry of patch must be centered
 * on the tangential component; normal_geometry is get_normal_velocity_geometry() of patch.
 */
double get_normal_velocity_difference(const SAMRAI::pdat::SideData<NDIM, double>& u_data,
                                      const SAMRAI::hier::Index<NDIM>& i,
                                      unsigned int bdry_normal_axis,
                                      bool bdry_is_lower,
                                      unsigned int tangential_axis,
                                      SAMRAI::solv::RobinBcCoefStrategy<NDIM>* normal_bc_coef,
                                      bool homogeneous_bc,
                                      const SAMRAI::hier::Patch<NDIM>& patch,
                                      const IBTK::ShiftedPatchGeometry& normal_geometry,
                                      const SAMRAI::hier::BoxArray<NDIM>& domain,
                                      double fill_time);

/*!
 * Return the normal velocity values on which the inhomogeneous Robin coefficient (+/-)(g/mu - d/h) of the TRACTION
 * condition for tangential_axis at the location i depends, each with its derivative; a value may appear more than once
 * and the derivatives add. The arguments are those of get_normal_velocity_difference(), except that ghost_box replaces
 * u_data and there is no homogeneous_bc.
 */
std::vector<std::pair<SAMRAI::pdat::SideIndex<NDIM>, double>>
get_traction_stencil(const SAMRAI::hier::Index<NDIM>& i,
                     unsigned int bdry_normal_axis,
                     bool bdry_is_lower,
                     unsigned int tangential_axis,
                     SAMRAI::solv::RobinBcCoefStrategy<NDIM>* normal_bc_coef,
                     const SAMRAI::hier::Patch<NDIM>& patch,
                     const IBTK::ShiftedPatchGeometry& normal_geometry,
                     const SAMRAI::hier::BoxArray<NDIM>& domain,
                     const SAMRAI::hier::Box<NDIM>& ghost_box,
                     double fill_time);

/*!
 * Add the transpose of the dependence of the TRACTION condition for velocity component tangential_axis on the normal
 * velocity to u_data: for each location of gcoef_data where tangential_bc_coef prescribes the traction, gcoef_data
 * times the derivative of the inhomogeneous Robin coefficient with respect to each normal velocity value.
 */
void accumulate_from_traction_bc_coefs(SAMRAI::pdat::SideData<NDIM, double>& u_data,
                                       const SAMRAI::pdat::ArrayData<NDIM, double>& gcoef_data,
                                       unsigned int tangential_axis,
                                       SAMRAI::solv::RobinBcCoefStrategy<NDIM>* tangential_bc_coef,
                                       SAMRAI::solv::RobinBcCoefStrategy<NDIM>* normal_bc_coef,
                                       bool homogeneous_bc,
                                       const SAMRAI::hier::Patch<NDIM>& patch,
                                       const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
                                       const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
                                       double fill_time);
} // namespace traction_stencil
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_ins_staggered_traction_stencil
