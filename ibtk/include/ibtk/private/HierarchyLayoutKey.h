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

#ifndef included_IBTK_HierarchyLayoutKey
#define included_IBTK_HierarchyLayoutKey

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <BoxArray.h>
#include <IntVector.h>
#include <PatchHierarchy.h>

#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
namespace detail
{
/*!
 * \brief Value copy of the layout of a range of hierarchy levels, which detects when the layout changes.
 *
 * A class that caches an object built only from the index-space layout, such as a coarse-fine boundary, and that
 * receives no notice when the hierarchy regrids can keep a HierarchyLayoutKey and rebuild the object only when
 * changed() returns true. For each level in the range, the key holds the boxes, the processor mapping, the
 * refinement ratio to level 0, and the periodic shift of the grid geometry. The key holds no PatchLevel or hierarchy,
 * so it does not extend a level's lifetime.
 *
 * A regrid can replace a PatchLevel with one that has the same layout, and the key does not detect that. An object
 * that holds a PatchLevel or its patch data, such as a communication schedule or a fill operator, is therefore not
 * safe to reuse when only the key is unchanged.
 */
class HierarchyLayoutKey
{
public:
    /*!
     * \brief Construct an empty key, which differs from the key of every level range.
     */
    HierarchyLayoutKey() = default;

    /*!
     * \brief Replace this key by the key of levels \a coarsest_ln through \a finest_ln of \a hierarchy.
     *
     * \return true if the new key differs from the previous one, which includes the first call.
     */
    bool changed(const SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy, int coarsest_ln, int finest_ln);

private:
    /*!
     * \brief Construct the key of levels \a coarsest_ln through \a finest_ln of \a hierarchy.
     */
    HierarchyLayoutKey(const SAMRAI::hier::PatchHierarchy<NDIM>& hierarchy, int coarsest_ln, int finest_ln);

    /*!
     * \brief Return true if both keys cover the same level range with the same layouts.
     */
    bool operator==(const HierarchyLayoutKey& other) const;

    /*!
     * \name Data that determine the key: the layout of one level, the number of the first level in the range, and the
     * layout of each level in the range.
     */
    //\{
    struct Level
    {
        SAMRAI::hier::BoxArray<NDIM> boxes;
        std::vector<int> mapping;
        SAMRAI::hier::IntVector<NDIM> ratio;
        SAMRAI::hier::IntVector<NDIM> periodic_shift;

        bool operator==(const Level& other) const;
    };

    int d_coarsest_ln = -1;
    std::vector<Level> d_levels;
    //\}
};
} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_HierarchyLayoutKey
