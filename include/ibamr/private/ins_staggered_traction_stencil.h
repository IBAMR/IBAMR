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
#include <Index.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <RobinBcCoefStrategy.h>
#include <SideData.h>

#include <functional>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
namespace detail
{
/*!
 * \brief Return whether any face of \a acoef_data and \a bcoef_data has a TRACTION condition.
 */
bool has_traction_face(const SAMRAI::pdat::ArrayData<NDIM, double>& acoef_data,
                       const SAMRAI::pdat::ArrayData<NDIM, double>& bcoef_data);

/*!
 * \brief Correct the coefficients of the TRACTION condition mu*(du_tan/dx_norm + du_norm/dx_tan) = g on the tangential
 * velocity component \a tangential_axis at the boundary normal to \a bdry_normal_axis. The corrected coefficients are
 * (0, 1) and a value that includes the tangential derivative du_norm/dx_tan of the normal velocity. The coefficients
 * of faces with a velocity condition are unchanged.
 *
 * \a get_mu returns the viscosity at a face, given the traction datum g there. The function evaluates it only at faces
 * with a TRACTION condition.
 * \a patch_hierarchy must be set.
 */
void correct_traction_coefficients(SAMRAI::pdat::ArrayData<NDIM, double>& acoef_data,
                                   SAMRAI::pdat::ArrayData<NDIM, double>& bcoef_data,
                                   SAMRAI::pdat::ArrayData<NDIM, double>& gcoef_data,
                                   const SAMRAI::pdat::SideData<NDIM, double>& u_data,
                                   const SAMRAI::hier::Patch<NDIM>& patch,
                                   unsigned int bdry_normal_axis,
                                   bool bdry_is_lower,
                                   unsigned int tangential_axis,
                                   SAMRAI::solv::RobinBcCoefStrategy<NDIM>* normal_bc_coef,
                                   bool homogeneous_bc,
                                   const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& patch_hierarchy,
                                   const std::function<double(const SAMRAI::hier::Index<NDIM>&, double)>& get_mu,
                                   double fill_time);

/*!
 * \brief Return the normal velocity on the face, farther from the boundary, of the ghost cell \a i_g outside the
 * physical boundary normal to \a bdry_normal_axis (the lower boundary if \a bdry_is_lower), such that the discrete
 * divergence of \a u_data vanishes in \a i_g. The function ignores \a u_data on that face.
 */
double get_divergence_free_ghost_value(const SAMRAI::pdat::SideData<NDIM, double>& u_data,
                                       const SAMRAI::hier::Index<NDIM>& i_g,
                                       unsigned int bdry_normal_axis,
                                       bool bdry_is_lower,
                                       const double* dx);

/*!
 * \brief Return whether any of \a u_bc_coefs is a StokesBcCoefStrategy with TRACTION conditions.
 */
bool has_normal_stress_term(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs);
} // namespace detail
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_ins_staggered_traction_stencil
