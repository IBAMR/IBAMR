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

#ifndef included_IBTK_CachedGhostCellFill
#define included_IBTK_CachedGhostCellFill

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/private/HierarchyLevelsIdentity.h>
#include <ibtk/private/data_kind.h>

#include <tbox/Pointer.h>

#include <PatchHierarchy.h>

#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
namespace detail
{
/*!
 * \brief Ghost-cell fill operator for a range of levels, rebuilt after a level is replaced or when a component changes
 * in more than its indices, Robin coefficient objects, and physical boundary operator, and reset otherwise.
 *
 * Replaced levels stay alive until the next call to get() or clear(). The indices must stay registered meanwhile. The
 * operator checks levels 0 through finest_ln, because the schedules of a fill can read every coarser level.
 */
class CachedGhostCellFill
{
public:
    /*!
     * \brief Construct an empty cache.
     */
    CachedGhostCellFill() = default;

    /*!
     * \brief Return the fill operator for \a comps on levels \a coarsest_ln through \a finest_ln of \a hierarchy, with
     * \a homogeneous_bc. IBTK::invalid_level_number means level 0 or the finest level.
     */
    HierarchyGhostCellInterpolation&
    get(const std::vector<HierarchyGhostCellInterpolation::InterpolationTransactionComponent>& comps,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
        int coarsest_ln = invalid_level_number,
        int finest_ln = invalid_level_number,
        bool homogeneous_bc = false);

    /*!
     * \brief Release the operator and the levels, so that the next call to get() builds a new operator.
     */
    void clear();

private:
    /*!
     * Properties of the patch data indices of a component that must stay equal for a reset to follow a change of the
     * indices.
     */
    struct IndexProperties
    {
        DataKind dst_kind{}, src_kind{};
        SAMRAI::hier::IntVector<NDIM> dst_ghosts = SAMRAI::hier::IntVector<NDIM>(0);
        SAMRAI::hier::IntVector<NDIM> src_ghosts = SAMRAI::hier::IntVector<NDIM>(0);
        bool in_place = false;
        std::vector<bool> null_robin_bc_coefs;

        bool operator==(const IndexProperties& that) const
        {
            return dst_kind == that.dst_kind && src_kind == that.src_kind && dst_ghosts == that.dst_ghosts &&
                   src_ghosts == that.src_ghosts && in_place == that.in_place &&
                   null_robin_bc_coefs == that.null_robin_bc_coefs;
        }
    };

    SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> d_fill;
    std::vector<HierarchyGhostCellInterpolation::InterpolationTransactionComponent> d_comps;
    std::vector<IndexProperties> d_index_properties;
    HierarchyLevelsIdentity d_levels;
    int d_coarsest_ln = invalid_level_number;
};
} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_CachedGhostCellFill
