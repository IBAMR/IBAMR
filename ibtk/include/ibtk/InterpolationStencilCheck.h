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

#ifndef included_IBTK_InterpolationStencilCheck
#define included_IBTK_InterpolationStencilCheck

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <tbox/Pointer.h>

#include <Box.h>
#include <BoxList.h>
#include <IntVector.h>
#include <PatchLevel.h>

#include <string>
#include <vector>

namespace IBTK
{
template <class T>
class LIndexSetData;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
template <int DIM>
class PatchData;
template <int DIM>
class PatchHierarchy;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class InterpolationStencilCheck stops the interpolation of Eulerian data to a point whose kernel stencil
 * includes data in a flagged cell.
 *
 * A cell is flagged if the physical domain contains it and the patch level does not cover it, or if it lies outside the
 * domain beside a boundary that the patch does not touch. The caller must call reset() whenever the hierarchy changes.
 * The ghost cell width of the data must not exceed the maximum ghost width of the patch data descriptor at reset().
 */
class InterpolationStencilCheck
{
public:
    /*!
     * \brief Find the flagged cells near each local patch of the hierarchy.
     * \param fill_sets_uncovered True if the ghost fill sets the cells of the domain that the level does not cover. The
     * check then flags only cells outside the domain.
     */
    void reset(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy, bool fill_sets_uncovered = false);

    /*!
     * \brief Emit an error, naming \a caller, if the \a kernel_fcn stencil of a point selected from \a idx_data in
     * \a interp_box of \a patch includes data of \a q_data in a flagged cell. The local PETSc indices of \a idx_data
     * index the positions in \a X_data.
     */
    template <class T>
    void check(const std::string& caller,
               SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
               const SAMRAI::hier::PatchData<NDIM>& q_data,
               const SAMRAI::hier::Box<NDIM>& interp_box,
               const double* X_data,
               SAMRAI::tbox::Pointer<LIndexSetData<T>> idx_data,
               const std::string& kernel_fcn) const;

    /*!
     * \brief Emit an error, naming \a caller, if the \a kernel_fcn stencil of a point in \a interp_box of \a patch
     * includes data of \a q_data in a flagged cell. \a X_data contains \a X_size values: NDIM coordinates for each
     * point.
     */
    void check(const std::string& caller,
               SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
               const SAMRAI::hier::PatchData<NDIM>& q_data,
               const SAMRAI::hier::Box<NDIM>& interp_box,
               const double* X_data,
               int X_size,
               const std::string& kernel_fcn) const;

private:
    /*!
     * Return whether the check flags cells of the patch. Emit an error if the last reset() did not find the patch.
     */
    bool hasFlaggedCells(const SAMRAI::hier::Patch<NDIM>& patch) const;

    /*!
     * Emit an error if the stencils of the points \a local_indices, shifted by \a periodic_shifts, include flagged
     * data.
     */
    void checkPoints(const std::string& caller,
                     const SAMRAI::hier::Patch<NDIM>& patch,
                     const SAMRAI::hier::PatchData<NDIM>& q_data,
                     const std::string& kernel_fcn,
                     const double* X_data,
                     const std::vector<int>& local_indices,
                     const std::vector<double>& periodic_shifts) const;

    /*!
     * Flagged data indices in the ghost region of each local patch, by level and patch number. The first lists, one per
     * staggering, contain the flagged data. The next lists, one per staggering, contain the flagged data that touch a
     * ghost cell outside the domain.
     */
    std::vector<std::vector<std::vector<SAMRAI::hier::BoxList<NDIM>>>> d_flagged_boxes;

    /*!
     * The levels of the hierarchy at the time of the last reset, to detect patches of another hierarchy configuration.
     */
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>> d_levels;

    /*!
     * Maximum ghost cell width of the patch data descriptor.
     */
    SAMRAI::hier::IntVector<NDIM> d_max_ghost_width = SAMRAI::hier::IntVector<NDIM>(0);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_InterpolationStencilCheck
