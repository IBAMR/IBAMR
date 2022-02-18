// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2021 by the IBAMR developers
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

#ifndef included_IBTK_IndexUtilities_inl_h
#define included_IBTK_IndexUtilities_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/IndexUtilities.h"

#include <cmath>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

inline SAMRAI::hier::Index<NDIM>
IndexUtilities::coarsen(const SAMRAI::hier::Index<NDIM>& i_fine, const SAMRAI::hier::Index<NDIM>& ratio)
{
    SAMRAI::hier::Index<NDIM> i_coarse;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        i_coarse(d) = i_fine(d) < 0 ? (i_fine(d) + 1) / ratio(d) - 1 : i_fine(d) / ratio(d);
    }
    return i_coarse;
} // coarsen

inline SAMRAI::hier::Index<NDIM>
IndexUtilities::refine(const SAMRAI::hier::Index<NDIM>& i_coarse, const SAMRAI::hier::Index<NDIM>& ratio)
{
    return i_coarse * ratio;
} // refine

template <class DoubleArray>
inline SAMRAI::hier::Index<NDIM>
IndexUtilities::getCellIndex(const DoubleArray& X,
                             const double* const x_lower,
                             const double* const x_upper,
                             const double* const dx,
                             const SAMRAI::hier::Index<NDIM>& ilower,
                             const SAMRAI::hier::Index<NDIM>& iupper)
{
    // NOTE: This expression guarantees consistency between neighboring patches, but it is still possible to get
    // inconsitent mappings on disjoint patches.
    SAMRAI::hier::Index<NDIM> idx;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        double dX_lower = X[d] - x_lower[d], dX_upper = X[d] - x_upper[d];
        if (std::abs(dX_lower) <= std::abs(dX_upper))
        {
            idx(d) = ilower(d) + floor(dX_lower / dx[d]);
        }
        else
        {
            idx(d) = iupper(d) + floor(dX_upper / dx[d]) + 1;
        }
    }
    return idx;
} // getCellIndex

template <typename Vector>
inline Vector
IndexUtilities::getSideCenter(const SAMRAI::hier::Patch<NDIM>& patch, const SAMRAI::pdat::SideIndex<NDIM>& side_idx)
{
    const int axis = side_idx.getAxis();
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
    const double* const patch_x_lower = patch_geom->getXLower();
    const SAMRAI::hier::Box<NDIM>& patch_box = patch.getBox();
    const SAMRAI::hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
    const double* const patch_dx = patch_geom->getDx();

    Vector side_coord;
    for (int d = 0; d < NDIM; ++d)
    {
        if (d == axis)
        {
            side_coord[d] = patch_x_lower[d] + patch_dx[d] * (static_cast<double>(side_idx(d) - patch_lower_idx(d)));
        }
        else
        {
            side_coord[d] =
                patch_x_lower[d] + patch_dx[d] * (static_cast<double>(side_idx(d) - patch_lower_idx(d)) + 0.5);
        }
    }
    for (int d = NDIM; d < side_coord.size(); ++d)
    {
        side_coord[d] = 0.0;
    }
    return side_coord;
} // getSideCenter

inline VectorNd
IndexUtilities::getSideCenter(const SAMRAI::hier::Patch<NDIM>& patch, const SAMRAI::pdat::SideIndex<NDIM>& side_idx)
{
    return getSideCenter<VectorNd>(patch, side_idx);
} // getSideCenter

template <typename Vector>
inline Vector
IndexUtilities::getSideCenter(const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> >& grid_geom,
                              const SAMRAI::hier::IntVector<NDIM>& ratio,
                              const SAMRAI::pdat::SideIndex<NDIM>& side_idx)
{
    const int axis = side_idx.getAxis();
    const double* const dx0 = grid_geom->getDx();
    double patch_dx[NDIM];
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        patch_dx[d] = dx0[d] / static_cast<double>(ratio(d));
    }
    const SAMRAI::hier::Box<NDIM> domain_box =
        SAMRAI::hier::Box<NDIM>::refine(grid_geom->getPhysicalDomain()[0], ratio);
    const double* const grid_x_lower = grid_geom->getXLower();
    const SAMRAI::hier::Index<NDIM>& grid_lower_idx = domain_box.lower();

    Vector side_coord;
    for (int d = 0; d < NDIM; ++d)
    {
        if (d == axis)
        {
            side_coord[d] = grid_x_lower[d] + patch_dx[d] * (static_cast<double>(side_idx(d) - grid_lower_idx(d)));
        }
        else
        {
            side_coord[d] =
                grid_x_lower[d] + patch_dx[d] * (static_cast<double>(side_idx(d) - grid_lower_idx(d)) + 0.5);
        }
    }
    for (int d = NDIM; d < side_coord.size(); ++d)
    {
        side_coord[d] = 0.0;
    }
    return side_coord;
} // getSideCenter

inline VectorNd
IndexUtilities::getSideCenter(const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> >& grid_geom,
                              const SAMRAI::hier::IntVector<NDIM>& ratio,
                              const SAMRAI::pdat::SideIndex<NDIM>& side_idx)
{
    return getSideCenter<VectorNd>(grid_geom, ratio, side_idx);
} // getSideCenter

template <class DoubleArray>
inline SAMRAI::hier::Index<NDIM>
IndexUtilities::getCellIndex(const DoubleArray& X,
                             const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> >& patch_geom,
                             const SAMRAI::hier::Box<NDIM>& patch_box)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(patch_geom);
#endif
    return getCellIndex(
        X, patch_geom->getXLower(), patch_geom->getXUpper(), patch_geom->getDx(), patch_box.lower(), patch_box.upper());
} // getCellIndex

template <class DoubleArray>
inline SAMRAI::hier::Index<NDIM>
IndexUtilities::getCellIndex(const DoubleArray& X,
                             const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> >& grid_geom,
                             const SAMRAI::hier::IntVector<NDIM>& ratio)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(grid_geom);
#endif
    const double* const dx0 = grid_geom->getDx();
    double dx[NDIM];
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        dx[d] = dx0[d] / static_cast<double>(ratio(d));
    }
    const SAMRAI::hier::Box<NDIM> domain_box =
        SAMRAI::hier::Box<NDIM>::refine(grid_geom->getPhysicalDomain()[0], ratio);
    return getCellIndex(X, grid_geom->getXLower(), grid_geom->getXUpper(), dx, domain_box.lower(), domain_box.upper());
} // getCellIndex

template <class DoubleArray>
inline SAMRAI::hier::Index<NDIM>
IndexUtilities::getAssignedCellIndex(const DoubleArray& X,
                                     const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> >& grid_geom,
                                     const SAMRAI::hier::IntVector<NDIM>& ratio)
{
    auto idx = getCellIndex(X, grid_geom, ratio);
    const auto periodic_shift = grid_geom->getPeriodicShift();
    for (int d = 0; d < NDIM; ++d)
    {
        // We don't shift things on periodic boundaries - the caller should have already shifted them
        if (!periodic_shift(d) && X[d] == grid_geom->getXUpper()[d])
        {
            idx(d) -= 1;
        }
    }
    return idx;
}

inline int
IndexUtilities::mapIndexToInteger(const SAMRAI::hier::Index<NDIM>& i,
                                  const SAMRAI::hier::Index<NDIM>& domain_lower,
                                  const SAMRAI::hier::Index<NDIM>& n_cells,
                                  const int depth,
                                  const int offset,
                                  const SAMRAI::hier::IntVector<NDIM>& periodic_shift)
{
    SAMRAI::hier::Index<NDIM> idx = i;
    for (int d = 0; d < NDIM; ++d)
    {
        if (periodic_shift(d) && idx(d) >= domain_lower(d) + periodic_shift(d))
        {
            idx(d) -= periodic_shift(d);
        }

        if (periodic_shift(d) && idx(d) < domain_lower(d))
        {
            idx(d) += periodic_shift(d);
        }
    }

#if (NDIM == 1)
    return (idx(0) - domain_lower(0) + depth * n_cells(0) + offset);
#elif (NDIM == 2)
    return (idx(0) - domain_lower(0) + (idx(1) - domain_lower(1)) * n_cells(0) + depth * n_cells(0) * n_cells(1) +
            offset);
#elif (NDIM == 3)
    return (idx(0) - domain_lower(0) + (idx(1) - domain_lower(1)) * n_cells(0) +
            (idx(2) - domain_lower(2)) * n_cells(0) * n_cells(1) + depth * n_cells(0) * n_cells(1) * n_cells(2) +
            offset);

#else
    TBOX_ASSERT(false);
    return -1;
#endif

} // mapIndexToInteger

inline SAMRAI::hier::IntVector<NDIM>
IndexUtilities::partitionPatchBox(std::vector<SAMRAI::hier::Box<NDIM> >& overlap_boxes,
                                  std::vector<SAMRAI::hier::Box<NDIM> >& nonoverlap_boxes,
                                  const SAMRAI::hier::Box<NDIM>& patch_box,
                                  const SAMRAI::hier::IntVector<NDIM>& box_size,
                                  const SAMRAI::hier::IntVector<NDIM>& overlap_size)
{
    // Compute number of nonoverlapping subdomains.
    const SAMRAI::hier::IntVector<NDIM>& patch_lower = patch_box.lower();
    const SAMRAI::hier::IntVector<NDIM>& patch_upper = patch_box.upper();
    SAMRAI::hier::IntVector<NDIM> cells = 1;
    cells += patch_upper - patch_lower;
    const SAMRAI::hier::IntVector<NDIM> subdomains = SAMRAI::hier::IntVector<NDIM>::max(cells / box_size, 1);
    const int n_subdomains = subdomains.getProduct();

    // Compute the boxes.
    nonoverlap_boxes.resize(n_subdomains);
    overlap_boxes.resize(n_subdomains);
    int counter = 0;
    SAMRAI::hier::IntVector<NDIM> I, box_lower, box_upper;
#if (NDIM == 3)
    for (I(2) = 0; I(2) < subdomains(2); ++I(2))
    {
#endif
        for (I(1) = 0; I(1) < subdomains(1); ++I(1))
        {
            for (I(0) = 0; I(0) < subdomains(0); ++I(0))
            {
                for (int d = 0; d < NDIM; ++d)
                {
                    box_lower(d) = patch_lower(d) + I(d) * box_size(d);
                    box_upper(d) = box_lower(d) + box_size(d) - 1;
                    if (I(d) + 1 == subdomains(d)) box_upper(d) += cells(d) % box_size(d);
                    box_upper(d) = std::min(box_upper(d), patch_upper(d));
                }
                nonoverlap_boxes[counter] = SAMRAI::hier::Box<NDIM>(box_lower, box_upper);
                overlap_boxes[counter] = SAMRAI::hier::Box<NDIM>::grow(nonoverlap_boxes[counter], overlap_size);
                ++counter;
            }
        }
#if (NDIM == 3)
    }
#endif
#if !defined(NDEBUG)
    TBOX_ASSERT(counter == n_subdomains);
#endif
    return subdomains;
} // partitionPatchBox

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_IndexUtilities_inl_h
