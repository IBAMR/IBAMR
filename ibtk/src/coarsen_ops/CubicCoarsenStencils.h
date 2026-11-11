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

#ifndef included_IBTK_CubicCoarsenStencils
#define included_IBTK_CubicCoarsenStencils

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <tbox/Pointer.h>

#include <ArrayData.h>
#include <BoxTree.h>
#include <Index.h>
#include <IntVector.h>
#include <PatchLevel.h>

#include <array>
#include <functional>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * Determines whether the boxes of a patch level, including their periodic images, cover a cell.
 */
class LevelCoverage
{
public:
    explicit LevelCoverage(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level);

    bool covers(SAMRAI::hier::Index<NDIM> cell) const;

private:
    SAMRAI::tbox::Pointer<SAMRAI::hier::BoxTree<NDIM>> d_box_tree;
    SAMRAI::hier::IntVector<NDIM> d_periodic_shift;
    SAMRAI::hier::Index<NDIM> d_domain_lower;
};

/*!
 * How the cubic coarsening of a coarse position treats an axis.
 */
enum class CubicCoarsenAxis
{
    NORMAL,  // the normal axis of side data: the coarse side is the fine side on it
    RATIO_2, // refinement ratio 2: the stencil depends on which coarse neighbors the fine level covers
    INTERIOR // an even refinement ratio of 4 or more: the cubic stencil inside the coarse cell
};

/*!
 * The stencil for a coarse position whose cubic stencil involves a fine cell outside the fine level. Along each RATIO_2
 * axis, it is the centered cubic, a one-sided quadratic, or the two-point average. Among the products that use only
 * covered coarse neighbors, the stencil is the mean of those with the fewest averaged axes, then the fewest one-sided
 * axes.
 */
class CubicCoarsenStencil
{
public:
    /*!
     * Choose the stencil. covered(offset) returns whether the fine level covers the coarse neighbor at the given
     * offset, whose components are -1, 0, or 1 along the RATIO_2 axes and 0 along the others.
     */
    CubicCoarsenStencil(const std::array<CubicCoarsenAxis, NDIM>& axes,
                        const std::function<bool(const SAMRAI::hier::IntVector<NDIM>&)>& covered);

    /*!
     * Return the coarse value from the fine data. Along a RATIO_2 or INTERIOR axis, fine_base_index is the lower of the
     * two fine indices nearest to the coarse center. Along the NORMAL axis, fine_base_index is the index of the fine
     * side.
     */
    double apply(const SAMRAI::pdat::ArrayData<NDIM, double>& fine_data,
                 const SAMRAI::hier::Index<NDIM>& fine_base_index,
                 int depth) const;

private:
    // The one-dimensional stencil of each axis in each of the selected products.
    int d_num_products = 0;
    std::array<std::array<int, NDIM>, (NDIM == 2) ? 16 : 64> d_products;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_CubicCoarsenStencils
