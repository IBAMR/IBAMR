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

// Helpers for the TRACTION velocity boundary condition on the staggered grid: the stencil of the normal velocity
// that its tangential difference uses, including beyond corners of the physical boundary, the geometry of the
// physical domain, and the tangential difference of the normal velocity together with its transpose. This header is
// not installed.
namespace IBAMR
{
namespace traction_stencil
{
/*!
 * The normal velocity used by the TRACTION condition at a face on a boundary normal to bdry_normal_axis, as the
 * weighted sum weight[0]*u(idx[0]) + weight[1]*u(idx[1]) of the normal velocity at faces on that boundary.
 *
 * A face is beyond a corner if the cell adjacent to it on the interior side of the boundary lies outside the physical
 * domain. The normal velocity is not available at such a face, and the stencil extrapolates linearly from the two
 * nearest faces on the boundary, or uses the nearest face if the boundary segment has only one. The boundary at which
 * the segment ends, which is adjacent to the face, is the one with location index adjacent_location_index. If that
 * boundary prescribes the normal velocity, the value is instead twice the prescribed velocity at the extension of the
 * face minus the value at the nearest face. Whether it does, and the prescribed velocity, depend on the boundary
 * conditions, which the caller evaluates; useAdjacentBoundaryValue() converts the stencil, and the caller adds twice
 * the prescribed velocity.
 */
struct NormalVelocityStencil
{
    std::array<SAMRAI::pdat::SideIndex<NDIM>, 2> idx;
    std::array<double, 2> weight;

    /*!
     * Whether the face is beyond a corner.
     */
    bool beyond_corner = false;

    /*!
     * The location index of the boundary adjacent to the face, if it is beyond a corner.
     */
    unsigned int adjacent_location_index = 0;

    /*!
     * Whether the value is twice the velocity that the adjacent boundary prescribes, plus the weighted sum.
     */
    bool uses_adjacent_boundary_value = false;

    /*!
     * Replace the weighted sum by minus the value at the nearest face on the boundary, and record that the adjacent
     * boundary's value is used. The face must be beyond a corner.
     */
    void useAdjacentBoundaryValue();
};

/*!
 * A patch geometry shifted by a given displacement, so that boundary condition objects evaluate their coefficients at
 * the locations of another component of the velocity. The geometry is built once, from the geometry that the patch has
 * at that time. Constructing a ShiftedPatchGeometry::Scope replaces the geometry of the patch by the shifted geometry,
 * and destroying the scope restores the geometry that the patch had, on every exit path. The patch is not otherwise
 * modified.
 */
class ShiftedPatchGeometry
{
public:
    /*!
     * Build the geometry of patch with its coordinates shifted by shift.
     */
    ShiftedPatchGeometry(const SAMRAI::hier::Patch<NDIM>& patch, const std::array<double, NDIM>& shift);

    /*!
     * Replace the geometry of a patch by a ShiftedPatchGeometry for the lifetime of the scope.
     *
     * \note The patch is passed to boundary condition objects as a reference to a patch that its owner can modify.
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
 * The physical domain refined by each refinement ratio, with its periodic images, by the components of the ratio. It
 * is filled by get_physical_domain() and owned by the object that uses it.
 */
using PhysicalDomainCache = std::map<std::array<int, NDIM>, SAMRAI::hier::BoxArray<NDIM>>;

/*!
 * Return the physical domain of hierarchy refined by ratio, together with its periodic images, so that a cell across
 * a periodic boundary is not mistaken for a cell outside the domain. The domain is computed once for each ratio and
 * stored in cache; the returned reference is valid as long as cache is.
 */
const SAMRAI::hier::BoxArray<NDIM>&
get_physical_domain(PhysicalDomainCache& cache,
                    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
                    const SAMRAI::hier::IntVector<NDIM>& ratio);

/*!
 * Return the stencil of the normal velocity at the face i_face one cell beyond a corner of the boundary normal to
 * bdry_normal_axis, at which the boundary segment ends along tangential_axis. The nearest face on the boundary is one
 * cell from i_face along tangential_axis in the direction step, which is +1 (-1) if the segment ends below (above)
 * i_face. The linear extrapolation uses the second nearest face if can_extrapolate is true.
 */
NormalVelocityStencil get_corner_normal_velocity_stencil(SAMRAI::hier::Index<NDIM> i_face,
                                                         unsigned int bdry_normal_axis,
                                                         unsigned int tangential_axis,
                                                         int step,
                                                         bool can_extrapolate);

/*!
 * Return the stencil of the normal velocity at the face i_face on the boundary normal to bdry_normal_axis, which is the
 * lower boundary if bdry_is_lower, for a condition that differences the normal velocity along tangential_axis.
 * Whether the face is beyond a corner is determined from domain, which is the physical domain with its periodic
 * images at the level of the face, and not from the patch. A face on the boundary that lies outside ghost_box along
 * tangential_axis is extrapolated linearly from the two nearest faces in ghost_box, which makes the difference there
 * the one-sided difference inside ghost_box; if ghost_box holds a single face along tangential_axis, the nearest face
 * is used.
 */
NormalVelocityStencil get_normal_velocity_stencil(SAMRAI::hier::Index<NDIM> i_face,
                                                  unsigned int bdry_normal_axis,
                                                  bool bdry_is_lower,
                                                  unsigned int tangential_axis,
                                                  const SAMRAI::hier::BoxArray<NDIM>& domain,
                                                  const SAMRAI::hier::Box<NDIM>& ghost_box);

/*!
 * Return the difference along tangential_axis of the normal velocity that the TRACTION condition for the velocity
 * component tangential_axis uses at the location i on the boundary normal to bdry_normal_axis, which is the lower
 * boundary if bdry_is_lower. It is the value at the face i minus the value at the next face below it along
 * tangential_axis, each taken from u_data with the stencil of get_normal_velocity_stencil() for the ghost box of
 * u_data and for domain.
 *
 * For a face beyond a corner, normal_bc_coef, the boundary condition object of the normal velocity component,
 * determines whether the adjacent boundary prescribes the normal velocity. If it does, the value at that face is
 * twice the prescribed velocity, which is zero if homogeneous_bc, minus the value at the nearest face on the boundary.
 * The geometry of patch must be centered on the tangential component; the coefficients of normal_bc_coef are
 * evaluated at fill_time with a geometry centered on the normal component, which is built when first needed and kept
 * in normal_geometry.
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
 * Return the dependence on the normal velocity of the inhomogeneous Robin coefficient that the TRACTION condition for
 * the velocity component tangential_axis sets at the location i on the boundary normal to bdry_normal_axis, which is
 * the lower boundary if bdry_is_lower. That coefficient is (+/-)(g/mu - d/h), in which d is the difference
 * get_normal_velocity_difference() computes, h is the grid spacing along tangential_axis, and the sign is negative on
 * a lower boundary. The result lists the normal velocity values that the coefficient depends on, each with the
 * derivative of the coefficient with respect to it. The part of the coefficient that does not depend on the normal
 * velocity, which includes twice the velocity prescribed by an adjacent boundary at a corner, is not part of the
 * result. A value may be listed more than once, and the derivatives add.
 *
 * The arguments are those of get_normal_velocity_difference(), except that ghost_box is the ghost box of the data
 * that the stencil indexes, which replaces u_data, and that there is no homogeneous_bc.
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
 * tangential_axis on the boundary box bdry_box. That condition sets the inhomogeneous Robin coefficient at a location
 * i to (+/-)(g/mu - d/h), in which d is get_normal_velocity_difference() at i, h is the grid spacing along
 * tangential_axis, and the sign is negative on a lower boundary. For each location i of gcoef_data at which
 * tangential_bc_coef prescribes the traction, this function adds gcoef_data(i) times the derivative of that
 * coefficient with respect to each normal velocity value to that value of u_data, in the patch interior or in its
 * ghost cells.
 *
 * The physical domain is taken from domain_cache and hierarchy when a location prescribes the traction. The remaining
 * arguments are those of get_normal_velocity_difference().
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

/*!
 * Return the normal velocity on the face of the ghost cell i_g outside a physical boundary that is farther from the
 * boundary, such that the discrete divergence of u_data vanishes in the ghost cell. The boundary is normal to
 * normal_axis and is the lower boundary if is_lower. The face of i_g on the boundary and the tangential faces of i_g
 * are read, and the value of u_data on the face that is returned is not.
 */
double get_divergence_free_ghost_value(const SAMRAI::pdat::SideData<NDIM, double>& u_data,
                                       const SAMRAI::hier::Index<NDIM>& i_g,
                                       unsigned int normal_axis,
                                       bool is_lower,
                                       const double* dx);

/*!
 * Return the dependence of get_divergence_free_ghost_value() on u_data as the list of the faces that it reads, each
 * with its coefficient. The faces are those of i_g, which are the face on the boundary and the faces tangential to the
 * boundary, which are ghost values.
 */
std::vector<std::pair<SAMRAI::pdat::SideIndex<NDIM>, double>>
get_divergence_free_ghost_value_stencil(const SAMRAI::hier::Index<NDIM>& i_g,
                                        unsigned int normal_axis,
                                        bool is_lower,
                                        const double* dx);

/*!
 * Return the dependence on u_data of u_I - u_div at the face i of a physical boundary normal to bdry_normal_axis,
 * which is the lower boundary if bdry_is_lower, as the list of the faces that it reads, each with its coefficient.
 * Here u_I is the normal velocity on the next face inside the domain and u_div is the value of
 * get_divergence_free_ghost_value() in the ghost cell abutting the boundary at i. The viscous normal-stress term that
 * StaggeredStokesPhysicalBoundaryHelper::addNormalTractionViscousTerm() adds at the boundary face is D/h^2 times this
 * difference. A face may be listed more than once, and the coefficients add. The faces tangential to the boundary
 * are ghost faces; their values are defined by the boundary conditions for the tangential velocity.
 */
std::vector<std::pair<SAMRAI::pdat::SideIndex<NDIM>, double>>
get_normal_stress_stencil(const SAMRAI::hier::Index<NDIM>& i,
                          unsigned int bdry_normal_axis,
                          bool bdry_is_lower,
                          const double* dx);
} // namespace traction_stencil
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_ins_staggered_traction_stencil
