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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/RobinPhysBdryPatchStrategy.h>

#include <tbox/Pointer.h>

#include <Patch.h>
#include <PatchGeometry.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>

#include "./ib_velocity_ghost_fill.h"

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

void
fill_velocity_phys_bdry_ghost_values(RobinPhysBdryPatchStrategy& u_phys_bdry_op,
                                     const int velocity_idx,
                                     Pointer<PatchHierarchy<NDIM>> hierarchy,
                                     const bool homogeneous_bc,
                                     const double fill_time)
{
    u_phys_bdry_op.setPatchDataIndex(velocity_idx);
    u_phys_bdry_op.setHomogeneousBc(homogeneous_bc);
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(velocity_idx))
        {
            continue;
        }
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            if (patch->getPatchGeometry()->getTouchesRegularBoundary())
            {
                u_phys_bdry_op.setPhysicalBoundaryConditions(
                    *patch, fill_time, patch->getPatchData(velocity_idx)->getGhostCellWidth());
            }
        }
    }
    u_phys_bdry_op.setHomogeneousBc(false);
    return;
} // fill_velocity_phys_bdry_ghost_values

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR
