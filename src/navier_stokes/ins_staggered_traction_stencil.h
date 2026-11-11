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
#include <PatchHierarchy.h>
#include <RobinBcCoefStrategy.h>
#include <SideData.h>
#include <SideIndex.h>

#include <array>
#include <optional>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

// Helpers for the TRACTION velocity boundary condition on the staggered grid. This header is not installed.
namespace IBAMR
{
namespace traction_stencil
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
 * Return the physical domain of hierarchy refined by ratio, with its periodic images.
 */
SAMRAI::hier::BoxArray<NDIM>
get_periodic_physical_domain(const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
                             const SAMRAI::hier::IntVector<NDIM>& ratio);

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
     * Constructor. The copy and its shifted geometry are built on the first call of getPatch() or getShiftedGeometry().
     */
    NormalVelocityPatch(const SAMRAI::hier::Patch<NDIM>& patch,
                        unsigned int bdry_normal_axis,
                        unsigned int tangential_axis);

    /*!
     * Return the box of the original patch.
     */
    const SAMRAI::hier::Box<NDIM>& getBox() const;

    /*!
     * Return the copy of the patch, built on first call.
     */
    SAMRAI::hier::Patch<NDIM>& getPatch();

    /*!
     * Return the shifted geometry of the copy, built on first call.
     */
    const IBTK::ShiftedPatchGeometry& getShiftedGeometry();

private:
    /*!
     * Build the copy and its shifted geometry if they do not exist.
     */
    void buildShiftedPatch();

    /*!
     * The original patch, the shift, and the copy and its shifted geometry.
     */
    const SAMRAI::hier::Patch<NDIM>& d_patch;
    std::array<double, NDIM> d_shift;
    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> d_shifted_patch;
    std::optional<IBTK::ShiftedPatchGeometry> d_shifted_geometry;
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
                                      NormalVelocityPatch& normal_velocity_patch,
                                      const SAMRAI::hier::BoxArray<NDIM>& domain,
                                      double fill_time);

/*!
 * Return the normal velocity on the face, farther from the boundary, of the ghost cell i_g outside the physical
 * boundary normal to bdry_normal_axis (the lower boundary if bdry_is_lower), such that the discrete divergence of
 * u_data vanishes in i_g. The function ignores u_data on that face.
 */
double get_divergence_free_ghost_value(const SAMRAI::pdat::SideData<NDIM, double>& u_data,
                                       const SAMRAI::hier::Index<NDIM>& i_g,
                                       unsigned int bdry_normal_axis,
                                       bool bdry_is_lower,
                                       const double* dx);

/*!
 * Return whether any of u_bc_coefs is a StokesBcCoefStrategy with TRACTION conditions, for which the operator adds a
 * normal stress term.
 */
bool has_normal_stress_term(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs);

} // namespace traction_stencil
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_ins_staggered_traction_stencil
