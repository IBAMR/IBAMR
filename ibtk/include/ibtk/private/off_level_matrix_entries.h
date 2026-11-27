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

#ifndef included_IBTK_off_level_matrix_entries
#define included_IBTK_off_level_matrix_entries

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/private/BoxArrayIndex.h>

#include <CellData.h>
#include <Index.h>
#include <PatchLevel.h>
#include <SideData.h>

#include <vector>

/////////////////////////////// FUNCTION DECLARATIONS ////////////////////////

namespace IBTK
{
namespace detail
{
/*!
 * \brief Set to zero the entries of a cell-centered matrix that couple a cell of a patch to a cell of the physical
 * domain outside the patch level.
 *
 * Such a cell lies on the coarse side of a coarse-fine boundary. hypre's multigrid solvers form the matrices of the
 * coarser grids from such entries and converge only if they are zero. Cells of other patches of the level, including
 * across a periodic boundary, are cells of the level. Entries that couple to cells outside the physical domain are left
 * to the boundary condition treatment. The diagonal entry is unchanged.
 * PoissonUtilities::adjustRHSAtCoarseFineBoundary() accounts for the dropped terms in the right-hand side. This
 * function does nothing on the coarsest level.
 *
 * \param[in,out] matrix_coefficients Matrix coefficients on the box of a patch of \p level with one depth per entry of
 * \p stencil and no ghost cells.
 * \param[in] level Patch level that contains the patch.
 * \param[in] level_boxes Index of the boxes of \p level and their periodic images, from periodic_level_box_index().
 * \param[in] stencil Offsets of the cells that the entries couple to.
 */
void clearOffLevelMatrixEntries(SAMRAI::pdat::CellData<NDIM, double>& matrix_coefficients,
                                const SAMRAI::hier::PatchLevel<NDIM>& level,
                                const BoxArrayIndex& level_boxes,
                                const std::vector<SAMRAI::hier::Index<NDIM>>& stencil);

/*!
 * \brief Set to zero the entries of a side-centered matrix that couple a side of a patch to a side of the same
 * component in the physical domain outside the patch level.
 *
 * The sides of the level are the sides of the cells of its patches, including the sides that patches share. See the
 * cell-centered version for the reason.
 *
 * \param[in,out] matrix_coefficients Matrix coefficients on the sides of the box of a patch of \p level with one depth
 * per entry of \p stencil and no ghost sides.
 * \param[in] level Patch level that contains the patch.
 * \param[in] level_boxes Index of the boxes of \p level and their periodic images, from periodic_level_box_index().
 * \param[in] stencil Offsets of the sides that the entries couple to, along the axes.
 */
void clearOffLevelMatrixEntries(SAMRAI::pdat::SideData<NDIM, double>& matrix_coefficients,
                                const SAMRAI::hier::PatchLevel<NDIM>& level,
                                const BoxArrayIndex& level_boxes,
                                const std::vector<SAMRAI::hier::Index<NDIM>>& stencil);
} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_off_level_matrix_entries
