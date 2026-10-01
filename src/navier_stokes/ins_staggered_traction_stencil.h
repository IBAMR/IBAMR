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

#include <ibamr/ibamr_enums.h>

#include <tbox/Pointer.h>

#include <BoxArray.h>
#include <Index.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <RobinBcCoefStrategy.h>
#include <SideIndex.h>

#include <array>
#include <functional>

/////////////////////////////// NAMESPACE ////////////////////////////////////

// This header is private to the library sources. It is shared by the classes
// that interpret the TRACTION velocity boundary condition on the staggered grid
// and is not installed.
namespace IBAMR
{
namespace traction_stencil
{
/*!
 * The normal velocity used by the TRACTION condition at a boundary face, as
 * weight[0]*u(idx[0]) + weight[1]*u(idx[1]) + offset.
 */
struct NormalVelocityStencil
{
    std::array<SAMRAI::pdat::SideIndex<NDIM>, 2> idx;
    std::array<double, 2> weight;
    double offset;

    /*!
     * Whether offset is twice the normal velocity that the condition on an
     * adjacent physical boundary prescribes. That boundary has location index
     * adjacent_location_index.
     */
    bool uses_adjacent_wall_velocity = false;
    unsigned int adjacent_location_index = 0;
};

/*!
 * Return the physical domain refined to the level of patch, together with its
 * periodic images, so that a cell across a periodic boundary is not mistaken
 * for a cell outside the domain. The name of the calling method is used in the
 * error message.
 */
SAMRAI::hier::BoxArray<NDIM>
get_physical_domain(const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
                    const SAMRAI::hier::Patch<NDIM>& patch,
                    const char* method);

/*!
 * Return the stencil of the normal velocity at a face that is one cell beyond
 * a corner of the boundary normal to bdry_normal_axis, across the boundary
 * normal to tangential_axis. The face lies beyond that boundary, at index
 * i_face(tangential_axis), and step is +1 if it lies beyond the lower boundary
 * and -1 if it lies beyond the upper boundary. The value at such a face is not
 * available, so it is determined by corner_type from the two nearest faces on
 * the boundary and, for WALL_AWARE, the condition on the adjacent boundary. For
 * WALL_AWARE, get_adjacent_wall_velocity is called with the location index of
 * the adjacent boundary. It returns true and sets the normal velocity
 * u_prescribed if the condition on the adjacent boundary prescribes the
 * velocity, and the stencil then includes the offset 2*u_prescribed.
 */
NormalVelocityStencil get_end_normal_velocity_stencil(
    SAMRAI::hier::Index<NDIM> i_face,
    unsigned int bdry_normal_axis,
    unsigned int tangential_axis,
    int step,
    TractionBcCornerType corner_type,
    const std::function<bool(double& u_prescribed, unsigned int adj_location_index)>& get_adjacent_wall_velocity);

/*!
 * Return the normal velocity on the boundary face with index i_face, where the
 * boundary is normal to bdry_normal_axis and the TRACTION condition differences
 * the normal velocity along tangential_axis. A face is beyond a corner if the
 * cell adjacent to it on the interior side of the boundary lies outside the
 * physical domain; this is determined from the domain, not from the patch. The
 * stencil of such a face is determined by get_end_normal_velocity_stencil(). A
 * face on the boundary is clamped to ghost_box along tangential_axis.
 */
NormalVelocityStencil get_normal_velocity_stencil(SAMRAI::hier::Index<NDIM> i_face,
                                                  unsigned int bdry_normal_axis,
                                                  bool bdry_is_lower,
                                                  unsigned int tangential_axis,
                                                  TractionBcCornerType corner_type,
                                                  SAMRAI::solv::RobinBcCoefStrategy<NDIM>* normal_bc_coef,
                                                  bool homogeneous_bc,
                                                  const SAMRAI::hier::Patch<NDIM>& patch,
                                                  const SAMRAI::hier::BoxArray<NDIM>& domain,
                                                  const SAMRAI::hier::Box<NDIM>& ghost_box,
                                                  double fill_time);
} // namespace traction_stencil
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_ins_staggered_traction_stencil
