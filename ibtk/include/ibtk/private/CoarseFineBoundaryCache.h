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

#ifndef included_IBTK_CoarseFineBoundaryCache
#define included_IBTK_CoarseFineBoundaryCache

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/private/HierarchyLayoutKey.h>

#include <tbox/Pointer.h>

#include <CoarseFineBoundary.h>
#include <IntVector.h>
#include <PatchHierarchy.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
namespace detail
{
/*!
 * \brief Hold a CoarseFineBoundary and construct it again only when the layout that determines it changes.
 *
 * Constructing a CoarseFineBoundary costs O(P^2) in the number of patches P on the level.
 */
class CoarseFineBoundaryCache
{
public:
    /*!
     * \brief Return the coarse-fine boundary of level \a ln, constructing it only if the layout of level \a ln or of
     * level 0 (boxes, processor mapping, ratio, and periodic shift) or the ghost width differs from that of the
     * previous call. The boundary does not depend on level ln - 1.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::CoarseFineBoundary<NDIM>>
    getCoarseFineBoundary(const SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy,
                          int ln,
                          const SAMRAI::hier::IntVector<NDIM>& ghost_width);

private:
    /*!
     * \name The data that determine the boundary, and the boundary. The layout of level 0 stands in for the physical
     * domain.
     */
    //\{
    HierarchyLayoutKey d_level_key;
    HierarchyLayoutKey d_level0_key;
    SAMRAI::hier::IntVector<NDIM> d_ghost_width;
    SAMRAI::tbox::Pointer<SAMRAI::hier::CoarseFineBoundary<NDIM>> d_coarse_fine_boundary;
    //\}
};
} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_CoarseFineBoundaryCache
