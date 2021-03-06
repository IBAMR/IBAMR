// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#include "ibtk/PhysicalBoundaryUtilities.h"

#include "BoundaryBox.h"
#include "Box.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "tbox/Array.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <algorithm>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

/*
 * Conventions:
 * ============
 *
 * 1d
 * --
 * node (co-dimension 1):
 * x_lo : 0
 * x_hi : 1
 *
 * 2d
 * --
 * edge (co-dimension 1):
 * x_lo: 0
 * x_hi: 1
 * y_lo: 2
 * y_hi: 3
 *
 * node (co-dimension 2):
 * x_lo, y_lo: 0
 * x_hi, y_lo: 1
 * x_lo, y_hi: 2
 * x_hi, y_hi: 3
 *
 * 3d
 * --
 *
 * face (co-dimension 1):
 * x_lo: 0
 * x_hi: 1
 * y_lo: 2
 * y_hi: 3
 * z_lo: 4
 * z_hi: 5
 *
 * edge (co-dimension 2):
 * y_lo, z_lo: 0
 * y_hi, z_lo: 1
 * y_lo, z_hi: 2
 * y_hi, z_hi: 3
 * x_lo, z_lo: 4
 * x_lo, z_hi: 5
 * x_hi, z_lo: 6
 * x_hi, z_hi: 7
 * x_lo, y_lo: 8
 * x_hi, y_lo: 9
 * x_lo, y_hi: 10
 * x_hi, y_hi: 11
 *
 * node (co-dimension 3):
 * x_lo, y_lo, z_lo: 0
 * x_hi, y_lo, z_lo: 1
 * x_lo, y_hi, z_lo: 2
 * x_hi, y_hi, z_lo: 3
 * x_lo, y_lo, z_hi: 4
 * x_hi, y_lo, z_hi: 5
 * x_lo, y_hi, z_hi: 6
 * x_hi, y_hi, z_hi: 7
 */

bool
PhysicalBoundaryUtilities::isLower(int loc, int codim, int direction)
{
    const BoundaryLookupTable<NDIM>* const bdry_lookup_table = BoundaryLookupTable<NDIM>::getLookupTable();

    if (codim == NDIM) return bdry_lookup_table->isLower(loc, codim, direction);

    if (codim == 1 && loc == 2 * direction)
    {
        return true;
    }
    else if (codim == 1)
    {
        return false;
    }

#if (NDIM == 3)
    if (codim == 2)
    {
        if (direction == 0) // x direction
        {
            return (loc == 4 || loc == 5 || loc == 8 || loc == 10);
        }
        if (direction == 1) // y direction
        {
            return (loc == 0 || loc == 2 || loc == 8 || loc == 9);
        }
        if (direction == 2) // z direction
        {
            return (loc == 0 || loc == 1 || loc == 4 || loc == 6);
        }
    }
#endif

    TBOX_ERROR("this statement should not be reached!\n");

    return false;
} // isLower

bool
PhysicalBoundaryUtilities::isUpper(int loc, int codim, int direction)
{
    const BoundaryLookupTable<NDIM>* const bdry_lookup_table = BoundaryLookupTable<NDIM>::getLookupTable();

    if (codim == NDIM) return bdry_lookup_table->isUpper(loc, codim, direction);

    if (codim == 1 && loc == 2 * direction + 1)
    {
        return true;
    }
    else if (codim == 1)
    {
        return false;
    }

#if (NDIM == 3)
    if (codim == 2)
    {
        if (direction == 0) // x direction
        {
            return (loc == 6 || loc == 7 || loc == 9 || loc == 11);
        }
        if (direction == 1) // y direction
        {
            return (loc == 1 || loc == 3 || loc == 10 || loc == 11);
        }
        if (direction == 2) // z direction
        {
            return (loc == 2 || loc == 3 || loc == 5 || loc == 7);
        }
    }
#endif

    TBOX_ERROR("this statement should not be reached!\n");

    return false;
} // isUpper

Array<BoundaryBox<NDIM> >
PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(const Patch<NDIM>& patch)
{
    return patch.getPatchGeometry()->getCodimensionBoundaries(1);
} // getPhysicalBoundaryCodim1Boxes

Array<BoundaryBox<NDIM> >
PhysicalBoundaryUtilities::getPhysicalBoundaryCodim2Boxes(const Patch<NDIM>& patch)
{
    return patch.getPatchGeometry()->getCodimensionBoundaries(2);
} // getPhysicalBoundaryCodim2Boxes

Array<BoundaryBox<NDIM> >
PhysicalBoundaryUtilities::getPhysicalBoundaryCodim3Boxes(const Patch<NDIM>& patch)
{
    return patch.getPatchGeometry()->getCodimensionBoundaries(3);
} // getPhysicalBoundaryCodim3Boxes

BoundaryBox<NDIM>
PhysicalBoundaryUtilities::trimBoundaryCodim1Box(const BoundaryBox<NDIM>& bdry_box, const Patch<NDIM>& patch)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bdry_box.getBoundaryType() == 1);
#endif
    // Trim a boundary box so it does not stick out past the corners of a patch.
    const Box<NDIM>& b_box = bdry_box.getBox();
    const Box<NDIM>& patch_box = patch.getBox();
    const unsigned int bdry_normal_axis = bdry_box.getLocationIndex() / 2;

    Box<NDIM> trimmed_b_box = b_box;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (d != bdry_normal_axis)
        {
            trimmed_b_box.lower()[d] = std::max(b_box.lower()[d], patch_box.lower()[d]);
            trimmed_b_box.upper()[d] = std::min(b_box.upper()[d], patch_box.upper()[d]);
        }
    }
    const BoundaryBox<NDIM> trimmed_bdry_box(trimmed_b_box, bdry_box.getBoundaryType(), bdry_box.getLocationIndex());
    return trimmed_bdry_box;
} // trimBoundaryCodim1Box

Box<NDIM>
PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(const BoundaryBox<NDIM>& bdry_box)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bdry_box.getBoundaryType() == 1);
#endif
    // Make surface box on boundary.
    Box<NDIM> side_bdry_box = bdry_box.getBox();
    const unsigned int location_index = bdry_box.getLocationIndex();
    const unsigned int bdry_normal_axis = location_index / 2;
    const bool bdry_lower_side = (location_index % 2) == 0;
    if (bdry_lower_side)
    {
        // On the lower side of a patch, the side indices are one higher than
        // the boundary cell indices in the direction normal to the boundary.
        side_bdry_box.shift(bdry_normal_axis, 1);
    }
    return side_bdry_box;
} // makeSideBoundaryCodim1Box

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
