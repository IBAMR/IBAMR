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
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

// Helpers for the TRACTION velocity boundary condition on the staggered grid. This header is not installed.
namespace IBAMR
{
namespace traction_stencil
{
/*!
 * \brief Stencil sum_k weight[k]*u(idx[k]) for the normal velocity in the TRACTION condition at a face on a boundary.
 * An unused entry has weight 0 and the index idx[0]. Beyond a corner, the stencil is an extrapolation along the
 * boundary segment: quadratic from three faces with "QUADRATIC" interpolation, otherwise linear from two faces, or the
 * value of the nearest face for a segment of one face. If the normal velocity u_b is prescribed on the adjacent
 * boundary (adjacent_location_index), the caller calls useAdjacentBoundaryValue() and adds prescribed_weight*u_b.
 */
struct NormalVelocityStencil
{
    static constexpr std::size_t NUM_FACES = 3;

    std::array<SAMRAI::pdat::SideIndex<NDIM>, NUM_FACES> idx;
    std::array<double, NUM_FACES> weight;

    bool beyond_corner = false;
    unsigned int adjacent_location_index = 0;
    bool uses_adjacent_boundary_value = false;
    double prescribed_weight = 2.0;

    /*!
     * Replace the weighted sum by -q_0, or by -2*q_0 + q_1/3 if quadratic and the boundary segment has two faces, in
     * which q_0 and q_1 are the values at the two nearest faces on the boundary, and set prescribed_weight to 2 or 8/3.
     * The member beyond_corner must be true.
     */
    void useAdjacentBoundaryValue(bool quadratic);
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
 * step along tangential_axis. The stencil is the extrapolation through the num_boundary_faces faces of the segment from
 * there, with 1 <= num_boundary_faces <= NormalVelocityStencil::NUM_FACES.
 */
NormalVelocityStencil get_corner_normal_velocity_stencil(SAMRAI::hier::Index<NDIM> i_face,
                                                         unsigned int bdry_normal_axis,
                                                         unsigned int tangential_axis,
                                                         int step,
                                                         int num_boundary_faces);

/*!
 * Return a patch with the box, numbers, geometry, and patch data of patch. Replacing the geometry of the copy leaves
 * patch unchanged.
 */
SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> copy_patch(const SAMRAI::hier::Patch<NDIM>& patch);

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
 * box of u_data. The patch of normal_velocity_patch must be the patch of u_data. Beyond a corner with a prescribed
 * normal velocity, the function passes quadratic to NormalVelocityStencil::useAdjacentBoundaryValue().
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
                                      bool quadratic,
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
 * Return whether any of u_bc_coefs is a StokesBcCoefStrategy with TRACTION conditions or, if quadratic, PSEUDO_TRACTION
 * conditions: the conditions for which StaggeredStokesPhysicalBoundaryHelper::addNormalTractionViscousTerm() adds a
 * term.
 */
bool has_normal_stress_term(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs, bool quadratic);

/*!
 * Return the side indices of the normal velocities on which the inhomogeneous Robin coefficient of the TRACTION
 * condition depends at the location i, each with the derivative of the coefficient with respect to it. An index may
 * appear more than once: its derivative is the sum of its entries. The arguments are those of
 * get_normal_velocity_difference(), except that ghost_box replaces u_data.
 */
std::vector<std::pair<SAMRAI::pdat::SideIndex<NDIM>, double>>
get_traction_stencil(const SAMRAI::hier::Index<NDIM>& i,
                     unsigned int bdry_normal_axis,
                     bool bdry_is_lower,
                     unsigned int tangential_axis,
                     SAMRAI::solv::RobinBcCoefStrategy<NDIM>* normal_bc_coef,
                     NormalVelocityPatch& normal_velocity_patch,
                     const SAMRAI::hier::BoxArray<NDIM>& domain,
                     const SAMRAI::hier::Box<NDIM>& ghost_box,
                     bool quadratic,
                     double fill_time);

/*!
 * Return the velocity faces, with coefficients that add for a repeated face, in the normal-stress term at the boundary
 * face i normal to bdry_normal_axis. The term is (11/3)*u_I - 2*u_B - (2/3)*u_I2 - u_div, less (4/3)*(u_I - u_div) at
 * PSEUDO_TRACTION, if quadratic, and u_I - u_div (none at PSEUDO_TRACTION) if linear. u_B, u_I, and u_I2 are the normal
 * velocities at i and the next two faces inside. u_div is get_divergence_free_ghost_value() at the ghost cell beyond i.
 */
std::vector<std::pair<SAMRAI::pdat::SideIndex<NDIM>, double>>
get_normal_stress_stencil(const SAMRAI::hier::Index<NDIM>& i,
                          unsigned int bdry_normal_axis,
                          bool bdry_is_lower,
                          const double* dx,
                          bool quadratic,
                          bool pseudo_traction);
} // namespace traction_stencil
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_ins_staggered_traction_stencil
