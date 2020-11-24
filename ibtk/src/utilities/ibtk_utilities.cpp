// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2020 by the IBAMR developers
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

#include <CartesianPatchGeometry.h>
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
        Pointer<Patch<NDIM> > patch = patch_level.getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const patch_dx = patch_geom->getDx();
        const double patch_dx_min = *std::min_element(patch_dx, patch_dx + NDIM);
        result = std::min(result, patch_dx_min);
        break; // all patches on the same level have the same dx values
    }

    result = IBTK_MPI::minReduction(result);

    return result;
} // get_min_patch_dx
} // namespace IBTK
