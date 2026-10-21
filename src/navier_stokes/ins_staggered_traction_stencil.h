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
#include <cstddef>
#include <utility>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

// Helpers for the TRACTION velocity boundary condition on the staggered grid. This header is not installed.
namespace IBAMR
{
namespace traction_stencil
{
/*!
 * The normal velocity that the TRACTION condition uses at a face on a boundary: the sum over k of weight[k]*u(idx[k]).
 * An entry that is not used has weight 0 and the index idx[0]. At a face beyond a corner, the stencil extrapolates from
 * the nearest faces on the boundary: quadratically from three faces, linearly from two, and with the value of the
 * nearest face if the boundary has one. If the adjacent boundary (adjacent_location_index) prescribes the normal
 * velocity u_b, the value is 2*u_b minus the value at the nearest face, or (8/3)*u_b - 2*q_0 + q_1/3 with q_0 and q_1
 * the two nearest faces if the caller asks for the quadratic value and the boundary has two faces: the caller calls
 * useAdjacentBoundaryValue() and adds prescribed_weight*u_b.
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
     * Replace the weighted sum by the part of the value that depends on the faces of the boundary, with the weight
     * prescribed_weight of the prescribed adjacent boundary value. The face must be beyond a corner.
     */
    void useAdjacentBoundaryValue(bool quadratic_tangential);
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
 * i_face + step along tangential_axis, and num_boundary_faces is the number of faces of the boundary from there to its
 * end: the stencil extrapolates quadratically from the nearest three faces if num_boundary_faces is at least 3,
 * linearly if it is 2, and takes the value of the nearest face if it is 1.
 */
NormalVelocityStencil get_corner_normal_velocity_stencil(SAMRAI::hier::Index<NDIM> i_face,
                                                         unsigned int bdry_normal_axis,
                                                         unsigned int tangential_axis,
                                                         int step,
                                                         int num_boundary_faces);

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
 * the adjacent boundary prescribes the normal velocity (zero if homogeneous_bc), and quadratic_tangential selects the
 * value described for NormalVelocityStencil. The geometry of patch must be centered on the tangential component;
 * normal_geometry is get_normal_velocity_geometry() of patch.
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
                                      bool quadratic_tangential,
                                      double fill_time);

/*!
 * Return the normal velocity values on which the inhomogeneous Robin coefficient (+/-)(g/mu - d/h) of the TRACTION
 * condition for tangential_axis at the location i depends, each with its derivative; a value may appear more than once
 * and the derivatives add. The arguments are those of get_normal_velocity_difference(), except that ghost_box replaces
 * u_data and there is no homogeneous_bc. The stencil takes the reflection value beyond a corner, not the quadratic
 * value.
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
 * quadratic_tangential is that of get_normal_velocity_difference().
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
                                       bool quadratic_tangential,
                                       double fill_time);

/*!
 * Return the normal velocity on the face of the ghost cell i_g outside a physical boundary that is farther from the
 * boundary, such that the discrete divergence of u_data vanishes in the ghost cell. The boundary is normal to
 * normal_axis and is the lower boundary if is_lower. The value of u_data on the face that is returned is not read.
 */
double get_divergence_free_ghost_value(const SAMRAI::pdat::SideData<NDIM, double>& u_data,
                                       const SAMRAI::hier::Index<NDIM>& i_g,
                                       unsigned int normal_axis,
                                       bool is_lower,
                                       const double* dx);

/*!
 * Return the faces of u_data that u_I - u_div reads at the face i of a physical boundary normal to bdry_normal_axis
 * (the lower boundary if bdry_is_lower), each with its coefficient; a face may appear more than once and the
 * coefficients add. u_I is the next normal velocity inside the domain and u_div is get_divergence_free_ghost_value() in
 * the ghost cell at i; the tangential faces are ghost faces.
 */
std::vector<std::pair<SAMRAI::pdat::SideIndex<NDIM>, double>>
get_normal_stress_stencil(const SAMRAI::hier::Index<NDIM>& i,
                          unsigned int bdry_normal_axis,
                          bool bdry_is_lower,
                          const double* dx);

/*!
 * Return whether any of u_bc_coefs is a StokesBcCoefStrategy with TRACTION conditions.
 */
bool has_traction_conditions(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs);
} // namespace traction_stencil
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_ins_staggered_traction_stencil
