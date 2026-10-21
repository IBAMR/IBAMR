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

#include <Box.h>
#include <BoxArray.h>
#include <Index.h>
#include <Patch.h>
#include <RobinBcCoefStrategy.h>
#include <SideIndex.h>

#include <array>

/////////////////////////////// NAMESPACE ////////////////////////////////////

// Helpers for the TRACTION velocity boundary condition on the staggered grid. This header is not installed.
namespace IBAMR
{
namespace traction_stencil
{
// The normal velocity used by the TRACTION condition: weight[0]*u(idx[0]) + weight[1]*u(idx[1]) + offset.
struct NormalVelocityStencil
{
    std::array<SAMRAI::pdat::SideIndex<NDIM>, 2> idx;
    std::array<double, 2> weight;
    double offset;
};

NormalVelocityStencil get_normal_velocity_stencil(SAMRAI::hier::Index<NDIM> i_face,
                                                  unsigned int bdry_normal_axis,
                                                  bool bdry_is_lower,
                                                  unsigned int tangential_axis,
                                                  SAMRAI::solv::RobinBcCoefStrategy<NDIM>* normal_bc_coef,
                                                  bool homogeneous_bc,
                                                  const SAMRAI::hier::Patch<NDIM>& patch,
                                                  const IBTK::ShiftedPatchGeometry& normal_geometry,
                                                  const SAMRAI::hier::BoxArray<NDIM>& domain,
                                                  const SAMRAI::hier::Box<NDIM>& ghost_box,
                                                  double fill_time);
} // namespace traction_stencil
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_ins_staggered_traction_stencil
