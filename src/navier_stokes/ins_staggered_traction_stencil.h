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

#include <Box.h>
#include <BoxArray.h>
#include <GridGeometry.h>
#include <Index.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchGeometry.h>
#include <RobinBcCoefStrategy.h>
#include <SideIndex.h>

#include <array>
#include <memory>

/////////////////////////////// NAMESPACE ////////////////////////////////////

// Helpers for the TRACTION velocity boundary condition on the staggered grid. This header is not installed.
namespace IBAMR
{
namespace traction_stencil
{
// A patch geometry shifted by a displacement, so that boundary coefficients are evaluated at the locations of another
// velocity component. A Scope replaces the geometry of a patch by it until the scope ends.
class ShiftedPatchGeometry
{
public:
    ShiftedPatchGeometry(const SAMRAI::hier::Patch<NDIM>& patch, const std::array<double, NDIM>& shift);

    class Scope
    {
    public:
        // The scope modifies the const patch only until it ends.
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

// The normal velocity used by the TRACTION condition: weight[0]*u(idx[0]) + weight[1]*u(idx[1]) + offset.
struct NormalVelocityStencil
{
    std::array<SAMRAI::pdat::SideIndex<NDIM>, 2> idx;
    std::array<double, 2> weight;
    double offset;
};

SAMRAI::hier::BoxArray<NDIM>
compute_physical_domain(const SAMRAI::tbox::Pointer<SAMRAI::hier::GridGeometry<NDIM>>& grid_geom,
                        const SAMRAI::hier::IntVector<NDIM>& ratio);

NormalVelocityStencil get_normal_velocity_stencil(SAMRAI::hier::Index<NDIM> i_face,
                                                  unsigned int bdry_normal_axis,
                                                  bool bdry_is_lower,
                                                  unsigned int tangential_axis,
                                                  SAMRAI::solv::RobinBcCoefStrategy<NDIM>* normal_bc_coef,
                                                  bool homogeneous_bc,
                                                  const SAMRAI::hier::Patch<NDIM>& patch,
                                                  std::unique_ptr<ShiftedPatchGeometry>& normal_geometry,
                                                  const SAMRAI::hier::BoxArray<NDIM>& domain,
                                                  const SAMRAI::hier::Box<NDIM>& ghost_box,
                                                  double fill_time);
} // namespace traction_stencil
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_ins_staggered_traction_stencil
