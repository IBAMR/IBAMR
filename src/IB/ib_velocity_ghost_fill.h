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
 * Use the boundary operator u_phys_bdry_op to set the ghost values of the velocity with patch data index u_data_idx
 * outside the physical boundaries, on every level where it is allocated. The ghost values inside the domain must be
 * filled before the call. The homogeneous flag of the operator is false after the call.
 */
void fill_velocity_phys_bdry_ghost_values(IBTK::RobinPhysBdryPatchStrategy& u_phys_bdry_op,
                                          int u_data_idx,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                          bool homogeneous_bc,
                                          double fill_time);
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_ib_velocity_ghost_fill
