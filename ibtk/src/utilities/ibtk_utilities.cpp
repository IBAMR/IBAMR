// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2026 by the IBAMR developers
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

#include <ibtk/IBTK_MPI.h>
#include <ibtk/ibtk_utilities.h>

#include <BoxGeometry.h>
#include <BoxOverlap.h>
#include <CartesianPatchGeometry.h>
#include <Patch.h>
#include <PatchData.h>
#include <PatchDataFactory.h>
#include <PatchDescriptor.h>
#include <PatchLevel.h>

#include <algorithm>
#include <cmath>

#include <ibtk/app_namespaces.h>

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////
double
get_min_patch_dx(const PatchLevel<NDIM>& patch_level)
{
    double result = std::numeric_limits<double>::max();

    // Some processors might not have any patches so its easier to just quit
    // after one loop operation than to check
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level.getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
        const double* const patch_dx = patch_geom->getDx();
        const double patch_dx_min = *std::min_element(patch_dx, patch_dx + NDIM);
        result = std::min(result, patch_dx_min);
        break; // all patches on the same level have the same dx values
    }

    result = IBTK_MPI::minReduction(result);

    return result;
} // get_min_patch_dx

void
copy_ghost_region(const PatchLevel<NDIM>& patch_level, const int dst_idx, const int src_idx)
{
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level.getPatch(p());
        Pointer<PatchData<NDIM>> dst_data = patch->getPatchData(dst_idx);
        Pointer<PatchData<NDIM>> src_data = patch->getPatchData(src_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(dst_data);
        TBOX_ASSERT(src_data);
#endif
        Pointer<PatchDescriptor<NDIM>> descriptor = patch->getPatchDescriptor();
        Pointer<BoxGeometry<NDIM>> dst_geometry =
            descriptor->getPatchDataFactory(dst_idx)->getBoxGeometry(patch->getBox());
        Pointer<BoxGeometry<NDIM>> src_geometry =
            descriptor->getPatchDataFactory(src_idx)->getBoxGeometry(patch->getBox());
        // Not overwriting the interior makes the overlap the part of dst's ghost box (restricted to src's ghost box)
        // that lies outside dst's interior, in the data's own centering.
        Pointer<BoxOverlap<NDIM>> overlap = dst_geometry->calculateOverlap(
            *src_geometry, src_data->getGhostBox(), /*overwrite_interior*/ false, IntVector<NDIM>(0));
        dst_data->copy(*src_data, *overlap);
    }
    return;
} // copy_ghost_region
} // namespace IBTK
