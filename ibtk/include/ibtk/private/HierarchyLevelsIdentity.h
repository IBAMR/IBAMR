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

#ifndef included_IBTK_HierarchyLevelsIdentity
#define included_IBTK_HierarchyLevelsIdentity

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <tbox/Pointer.h>

#include <PatchHierarchy.h>
#include <PatchLevel.h>

#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
namespace detail
{
/*!
 * \brief Identity of the PatchLevel objects of a range of levels, which detects when any of them is replaced.
 *
 * A class that caches an object that holds PatchLevel objects or their patch data, such as a communication schedule
 * or a fill operator, and that receives no notice when the hierarchy regrids can keep a HierarchyLevelsIdentity and
 * rebuild the object whenever changed() returns true. A regrid can replace the PatchLevel object of a level even when
 * the boxes and processor mapping are unchanged, and the identity reports that change.
 *
 * The identity stores the PatchLevel objects, so the levels that a regrid replaces, and their patch data, remain in
 * memory until the next call to changed(). A stored level cannot be freed, so a level that a later regrid creates never
 * has its address. The identity needs no record of the hierarchy: each PatchLevel object belongs to one hierarchy, and
 * the range of levels is never empty.
 *
 * \see HierarchyLayoutKey, for caches of objects built from index-space data alone.
 */
class HierarchyLevelsIdentity
{
public:
    /*!
     * \brief Construct an empty identity, which differs from the identity of every level range.
     */
    HierarchyLevelsIdentity() = default;

    /*!
     * \brief Replace this identity by the identity of levels \a coarsest_ln through \a finest_ln of \a hierarchy.
     *
     * \return true if the level range or any PatchLevel object in the range differs from the previous call, which
     * includes the first call. The range must contain at least one level.
     */
    bool changed(const SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy, int coarsest_ln, int finest_ln);

private:
    /*!
     * \name Data that determine the identity: the number of the first level in the range and the level objects in
     * the range.
     */
    //\{
    int d_coarsest_ln = -1;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>> d_levels;
    //\}
};
} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_HierarchyLevelsIdentity
