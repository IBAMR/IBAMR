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

#include <ibtk/PhysicalBoundaryUtilities.h>

#include <tbox/Array.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <BoundaryBox.h>
#include <Box.h>
#include <BoxArray.h>
#include <BoxList.h>
#include <GridGeometry.h>
#include <Index.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchGeometry.h>

#include <algorithm>
#include <vector>

#include <ibtk/namespaces.h> // IWYU pragma: keep

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

Array<BoundaryBox<NDIM>>
PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(const Patch<NDIM>& patch)
{
    return patch.getPatchGeometry()->getCodimensionBoundaries(1);
} // getPhysicalBoundaryCodim1Boxes

bool
PhysicalBoundaryUtilities::touchesOppositeBoundary(const Array<BoundaryBox<NDIM>>& physical_codim1_boxes,
                                                   const unsigned int location_index)
{
    // The boundaries of an axis have the location indices 2 * axis and 2 * axis + 1.
    const int opposite_index = location_index % 2 == 0 ? location_index + 1 : location_index - 1;
    for (int n = 0; n < physical_codim1_boxes.size(); ++n)
    {
        if (physical_codim1_boxes[n].getLocationIndex() == opposite_index) return true;
    }
    return false;
} // touchesOppositeBoundary

Array<BoundaryBox<NDIM>>
PhysicalBoundaryUtilities::getPhysicalBoundaryCodim2Boxes(const Patch<NDIM>& patch)
{
    return patch.getPatchGeometry()->getCodimensionBoundaries(2);
} // getPhysicalBoundaryCodim2Boxes

Array<BoundaryBox<NDIM>>
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

std::vector<IntVector<NDIM>>
PhysicalBoundaryUtilities::getPeriodicImageShifts(const IntVector<NDIM>& periodic_shift)
{
    std::vector<IntVector<NDIM>> shifts;
    int num_offsets = 1;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        num_offsets *= 3;
    }
    for (int offset = 0; offset < num_offsets; ++offset)
    {
        IntVector<NDIM> image_shift(0);
        bool is_image = false;
        bool is_valid = true;
        for (unsigned int d = 0, k = offset; d < NDIM; ++d, k /= 3)
        {
            const int n = static_cast<int>(k % 3) - 1;
            image_shift(d) = n * periodic_shift(d);
            is_image = is_image || n != 0;
            is_valid = is_valid && (n == 0 || periodic_shift(d) != 0);
        }
        if (is_image && is_valid)
        {
            shifts.push_back(image_shift);
        }
    }
    return shifts;
} // getPeriodicImageShifts

BoxArray<NDIM>
PhysicalBoundaryUtilities::getPeriodicPhysicalDomain(const Pointer<GridGeometry<NDIM>>& grid_geom,
                                                     const IntVector<NDIM>& ratio)
{
    BoxArray<NDIM> domain = grid_geom->getPhysicalDomain();
    domain.refine(ratio);
    BoxList<NDIM> domain_list(domain);
    for (const IntVector<NDIM>& image_shift : getPeriodicImageShifts(grid_geom->getPeriodicShift(ratio)))
    {
        BoxArray<NDIM> image_domain = domain;
        for (int b = 0; b < image_domain.size(); ++b)
        {
            image_domain[b].shift(image_shift);
        }
        domain_list.unionBoxes(BoxList<NDIM>(image_domain));
    }
    return BoxArray<NDIM>(domain_list);
} // getPeriodicPhysicalDomain

bool
PhysicalBoundaryUtilities::isBoundaryFace(hier::Index<NDIM> i_face,
                                          const unsigned int bdry_normal_axis,
                                          const bool bdry_is_lower,
                                          const BoxArray<NDIM>& domain)
{
    if (!bdry_is_lower)
    {
        i_face(bdry_normal_axis) -= 1;
    }
    return domain.contains(i_face);
} // isBoundaryFace

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
