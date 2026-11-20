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

#include <Box.h>
#include <BoxArray.h>
#include <Index.h>
#include <IntVector.h>
#include <Patch.h>
#include <RobinBcCoefStrategy.h>
#include <SideData.h>
#include <SideIndex.h>

#include <array>

/////////////////////////////// NAMESPACE ////////////////////////////////////

// Helpers for the TRACTION velocity boundary condition on the staggered grid. This header is not installed.
namespace IBAMR
{
/*!
 * \brief Stencil weight[0]*u(idx[0]) + weight[1]*u(idx[1]) for the normal velocity in the TRACTION condition at a face
 * on a boundary. Beyond a corner, the weights are those of linear extrapolation, or of constant extrapolation from a
 * boundary segment of one face. If the normal velocity u_b is prescribed on the adjacent boundary
 * (adjacent_location_index), the caller calls useAdjacentBoundaryValue() and adds 2*u_b.
 */
struct NormalVelocityStencil
{
    std::array<SAMRAI::pdat::SideIndex<NDIM>, 2> idx;
    std::array<double, 2> weight;

    bool beyond_corner = false;
    unsigned int adjacent_location_index = 0;

    /*!
     * Replace the weighted sum by minus the value at the nearest face. The member beyond_corner must be true.
     */
    void useAdjacentBoundaryValue();
};

/*!
 * Return the stencil of the normal velocity at the face i_face beyond a corner of the boundary normal to
 * bdry_normal_axis, where the boundary segment ends along tangential_axis. The nearest face on the boundary is i_face +
 * step along tangential_axis. The stencil also uses the second nearest face if can_extrapolate is true.
 */
NormalVelocityStencil get_corner_normal_velocity_stencil(SAMRAI::hier::Index<NDIM> i_face,
                                                         unsigned int bdry_normal_axis,
                                                         unsigned int tangential_axis,
                                                         int step,
                                                         bool can_extrapolate);

/*!
 * \brief A copy of a patch whose geometry is shifted from the faces of velocity component tangential_axis to the faces
 * of component bdry_normal_axis, for evaluating the boundary condition coefficients of the normal component. The
 * original patch is unchanged.
 */
class NormalVelocityPatch
{
public:
    /*!
     * Constructor. Build the copy and its shifted geometry.
     */
    NormalVelocityPatch(const SAMRAI::hier::Patch<NDIM>& patch,
                        unsigned int bdry_normal_axis,
                        unsigned int tangential_axis);

    /*!
     * Return the box of the original patch.
     */
    const SAMRAI::hier::Box<NDIM>& getBox() const;

    /*!
     * Return the copy of the patch.
     */
    SAMRAI::hier::Patch<NDIM>& getPatch() const;

    /*!
     * Return the shifted geometry of the copy.
     */
    const IBTK::ShiftedPatchGeometry& getShiftedGeometry() const;

private:
    /*!
     * The original patch, the copy, and the shifted geometry of the copy.
     */
    const SAMRAI::hier::Patch<NDIM>& d_patch;
    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> d_shifted_patch;
    IBTK::ShiftedPatchGeometry d_shifted_geometry;
};

/*!
 * Return the difference along tangential_axis of the normal velocity in the TRACTION condition for velocity component
 * tangential_axis at i on the boundary normal to bdry_normal_axis, or NaN if the face i - 1 or i lies outside the ghost
 * box of u_data. The patch of normal_velocity_patch must be the patch of u_data.
 */
double get_normal_velocity_difference(const SAMRAI::pdat::SideData<NDIM, double>& u_data,
                                      const SAMRAI::hier::Index<NDIM>& i,
                                      unsigned int bdry_normal_axis,
                                      bool bdry_is_lower,
                                      unsigned int tangential_axis,
                                      SAMRAI::solv::RobinBcCoefStrategy<NDIM>* normal_bc_coef,
                                      bool homogeneous_bc,
                                      const NormalVelocityPatch& normal_velocity_patch,
                                      const SAMRAI::hier::BoxArray<NDIM>& domain,
                                      double fill_time);

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_ins_staggered_traction_stencil
