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

#include <ibtk/IBTK_MPI.h>
#include <ibtk/IndexUtilities.h>

#include <tbox/Array.h>

#include <BoundaryBox.h>
#include <CoarseFineBoundary.h>
#include <Patch.h>
#include <PatchLevel.h>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

bool
IndexUtilities::patchLevelCoversEntireDomain(const int level_num,
                                             Pointer<PatchLevel<NDIM>> patch_level,
                                             Pointer<CoarseFineBoundary<NDIM>> cf_boundary)
{
    if (level_num == 0) return true;

    int local_cf_bdry_box_size = 0;
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
        const Array<BoundaryBox<NDIM>>& type_1_cf_bdry = cf_boundary->getBoundaries(patch->getPatchNumber(),
                                                                                    /* boundary type */ 1);
        local_cf_bdry_box_size += type_1_cf_bdry.size();
    }
    return IBTK_MPI::sumReduction(local_cf_bdry_box_size) == 0;
} // patchLevelCoversEntireDomain

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
