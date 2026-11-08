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

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
namespace detail
{
/*!
 * \brief Return true if every component of the refinement ratio is even.
 */
inline bool
ratios_are_even(const SAMRAI::hier::IntVector<NDIM>& ratio)
{
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (ratio(d) % 2 != 0)
        {
            return false;
        }
    }
    return true;
} // ratios_are_even

/*!
 * \brief Class LevelCoverage determines whether the boxes of a patch level, including their periodic images, cover a
 * cell.
 */
class LevelCoverage
{
public:
    /*!
     * \brief Constructor. Use the boxes of the level and their periodic images.
     */
    explicit LevelCoverage(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level);

    /*!
     * \brief Return true if the boxes of the level or their periodic images contain the cell.
     */
    bool covers(SAMRAI::hier::Index<NDIM> cell) const;

private:
    SAMRAI::tbox::Pointer<SAMRAI::hier::BoxTree<NDIM>> d_box_tree;
    SAMRAI::hier::IntVector<NDIM> d_periodic_shift;
    SAMRAI::hier::Index<NDIM> d_domain_lower;
};

/*!
 * \brief How the cubic coarsening of a coarse position treats an axis.
 */
enum class CubicCoarsenAxis
{
    NORMAL,  // the normal axis of side data: the coarse side is the fine side on it
    RATIO_2, // refinement ratio 2: the stencil depends on which coarse neighbors the fine level covers
    INTERIOR // an even refinement ratio of 4 or more: the cubic stencil inside the coarse cell
};

/*!
 * \brief Class CubicCoarsenStencil is the stencil for a coarse position whose cubic stencil involves a fine cell
 * outside the fine level.
 *
 * Along each RATIO_2 axis, the stencil is the centered cubic, a one-sided quadratic, or the two-point average. Among
 * the products that use only covered coarse neighbors, the stencil is the mean of those with the fewest averaged axes,
 * then the fewest one-sided axes.
 */
class CubicCoarsenStencil
{
public:
    /*!
     * \brief Constructor. Choose the stencil. covered(offset) returns whether the fine level covers the coarse
     * neighbor at the given offset, whose components are -1, 0, or 1 along the RATIO_2 axes and 0 along the others.
     */
    template <class CoveredPredicate>
    CubicCoarsenStencil(const std::array<CubicCoarsenAxis, NDIM>& axes, const CoveredPredicate& covered)
        : CubicCoarsenStencil(axes, get_neighbor_coverage(axes, covered))
    {
    }

    /*!
     * \brief Return the coarse value from the fine data. Along a RATIO_2 or INTERIOR axis, fine_base_index is the
     * lower of the two fine indices nearest to the coarse center. Along the NORMAL axis, fine_base_index is the index
     * of the fine side.
     */
    double apply(const SAMRAI::pdat::ArrayData<NDIM, double>& fine_data,
                 const SAMRAI::hier::Index<NDIM>& fine_base_index,
                 int depth) const;

private:
    // The number of coarse neighbors of a position, itself included.
    static constexpr int NUM_NEIGHBORS = (NDIM == 2) ? 9 : 27;

    // Whether the fine level covers each coarse neighbor, indexed as in get_neighbor_coverage().
    using NeighborCoverage = std::array<bool, NUM_NEIGHBORS>;

    // Choose the stencil from the coverage of the coarse neighbors.
    CubicCoarsenStencil(const std::array<CubicCoarsenAxis, NDIM>& axes, const NeighborCoverage& neighbor_covered);

    // Return whether the fine level covers each coarse neighbor along the RATIO_2 axes. The neighbor with index
    // k = sum_d (offset(d) + 1) 3^d has the offset with components -1, 0, or 1; the others count as covered.
    template <class CoveredPredicate>
    static NeighborCoverage get_neighbor_coverage(const std::array<CubicCoarsenAxis, NDIM>& axes,
                                                  const CoveredPredicate& covered)
    {
        NeighborCoverage neighbor_covered;
        neighbor_covered.fill(true);
        for (int k = 0; k < NUM_NEIGHBORS; ++k)
        {
            SAMRAI::hier::IntVector<NDIM> offset(0);
            bool is_neighbor = false;
            bool along_ratio_2 = true;
            int rest = k;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                offset(d) = rest % 3 - 1;
                rest /= 3;
                if (offset(d) != 0)
                {
                    is_neighbor = true;
                    along_ratio_2 = along_ratio_2 && axes[d] == CubicCoarsenAxis::RATIO_2;
                }
            }
            if (is_neighbor && along_ratio_2)
            {
                neighbor_covered[k] = covered(offset);
            }
        }
        return neighbor_covered;
    } // get_neighbor_coverage

    // The largest number of products: one choice among four one-dimensional stencils on each axis.
    static constexpr int MAX_PRODUCTS = (NDIM == 2) ? 16 : 64;

    // The one-dimensional stencil of each axis in each of the selected products.
    int d_num_products = 0;
    std::array<std::array<int, NDIM>, MAX_PRODUCTS> d_products{};
};
} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_CubicCoarsenStencils
