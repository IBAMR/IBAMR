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

#ifndef included_IBAMR_ib_velocity_ghost_fill
#define included_IBAMR_ib_velocity_ghost_fill

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <tbox/Pointer.h>

#include <PatchHierarchy.h>

namespace IBTK
{
class RobinPhysBdryPatchStrategy;
} // namespace IBTK

/////////////////////////////// NAMESPACE ////////////////////////////////////

// Helper that sets the ghost values of an Eulerian velocity at the physical boundaries with a velocity boundary
// operator. This header is not installed.
namespace IBAMR
{
/*!
 * Set the ghost values outside the physical boundaries of the velocity with patch data index velocity_idx with the
 * boundary operator u_phys_bdry_op, on every level on which the velocity is allocated. The ghost values inside the
 * domain must already be filled. The operator is left with homogeneous_bc false.
 */
void fill_velocity_phys_bdry_ghost_values(IBTK::RobinPhysBdryPatchStrategy& u_phys_bdry_op,
                                          int velocity_idx,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                          bool homogeneous_bc,
                                          double fill_time);
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_ib_velocity_ghost_fill
