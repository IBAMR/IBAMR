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

#include <ibtk/private/CoarseFineBoundaryCache.h>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
namespace detail
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

Pointer<CoarseFineBoundary<NDIM>>
CoarseFineBoundaryCache::getCoarseFineBoundary(const PatchHierarchy<NDIM>& hierarchy,
                                               const int ln,
                                               const IntVector<NDIM>& ghost_width)
{
    // Update both keys, even when the first differs.
    const bool level_changed = d_level_key.changed(hierarchy, ln, ln);
    const bool level0_changed = d_level0_key.changed(hierarchy, 0, 0);
    d_ln = ln;
    if (level_changed)
    {
        d_level_box_index_is_current = false;
    }
    if (d_coarse_fine_boundary && !level_changed && !level0_changed && ghost_width == d_ghost_width)
    {
        return d_coarse_fine_boundary;
    }
    d_coarse_fine_boundary = new CoarseFineBoundary<NDIM>(hierarchy, ln, ghost_width);
    d_ghost_width = ghost_width;
    return d_coarse_fine_boundary;
} // getCoarseFineBoundary

const BoxArrayIndex&
CoarseFineBoundaryCache::getLevelBoxIndex(const PatchHierarchy<NDIM>& hierarchy)
{
    if (d_ln >= 0 && !d_level_box_index_is_current)
    {
        const Pointer<PatchLevel<NDIM>> level = hierarchy.getPatchLevel(d_ln);
        d_level_box_index = periodic_level_box_index(*level);
        d_level_box_index_is_current = true;
    }
    return d_level_box_index;
} // getLevelBoxIndex

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
