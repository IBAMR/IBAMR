// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/StokesBcCoefStrategy.h>
#include <ibamr/ibamr_enums.h>

#include <tbox/Utilities.h>

#include <BoundaryBox.h>
#include <Box.h>
#include <BoxArray.h>
#include <Index.h>
#include <Patch.h>
#include <SideIndex.h>

#include <utility>
#include <vector>

#include <ibamr/namespaces.h> // IWYU pragma: keep

namespace IBAMR
{
class StokesSpecifications;
} // namespace IBAMR

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
StokesBcCoefStrategy::setStokesSpecifications(const StokesSpecifications* const problem_coefs)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(problem_coefs);
#endif
    d_problem_coefs = problem_coefs;
    return;
} // setStokesSpecifications

void
StokesBcCoefStrategy::setTargetVelocityPatchDataIndex(int u_target_data_idx)
{
    d_u_target_data_idx = u_target_data_idx;
    return;
} // setTargetVelocityPatchDataIndex

void
StokesBcCoefStrategy::clearTargetVelocityPatchDataIndex()
{
    d_u_target_data_idx = invalid_index;
    return;
} // clearTargetVelocityPatchDataIndex

void
StokesBcCoefStrategy::setTargetPressurePatchDataIndex(int p_target_data_idx)
{
    d_p_target_data_idx = p_target_data_idx;
    return;
} // setPressurePatchDataIndex

void
StokesBcCoefStrategy::clearTargetPressurePatchDataIndex()
{
    d_p_target_data_idx = invalid_index;
    return;
} // clearPressurePatchDataIndex

void
StokesBcCoefStrategy::setTractionBcType(TractionBcType bc_type)
{
    d_traction_bc_type = bc_type;
    return;
} // setTractionBcType

TractionBcType
StokesBcCoefStrategy::getTractionBcType() const
{
    return d_traction_bc_type;
} // getTractionBcType

std::vector<std::pair<SideIndex<NDIM>, double>>
StokesBcCoefStrategy::getTractionStencil(const hier::Index<NDIM>& /*i*/,
                                         const Patch<NDIM>& /*patch*/,
                                         const BoundaryBox<NDIM>& /*bdry_box*/,
                                         const BoxArray<NDIM>& /*domain*/,
                                         const Box<NDIM>& /*ghost_box*/,
                                         const double /*fill_time*/) const
{
    if (d_traction_bc_type == TRACTION)
    {
        TBOX_ERROR("StokesBcCoefStrategy::getTractionStencil():\n"
                   << "  this boundary condition object does not provide the stencil of its TRACTION conditions\n"
                   << "  override getTractionStencil() or use the PSEUDO_TRACTION type\n");
    }
    return {};
} // getTractionStencil

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
