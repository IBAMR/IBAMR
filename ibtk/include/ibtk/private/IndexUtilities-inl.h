// Filename: IndexUtilities-inl.h
// Created on 18 Jun 2005 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_IndexUtilities_inl_h
#define included_IndexUtilities_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <cmath>

#include "boost/math/special_functions/round.hpp"
#include "ibtk/IndexUtilities.h"

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

inline int
IndexUtilities::mapIndexToInteger(const SAMRAI::hier::Index<NDIM>& i,
                                  const SAMRAI::hier::Index<NDIM>& domain_lower,
                                  const SAMRAI::hier::Index<NDIM>& num_cells,
                                  const int depth,
                                  const int offset)
{
#if (NDIM == 1)
    return (i(0) - domain_lower(0) + depth * num_cells(0) + offset);
#elif(NDIM == 2)
    return (i(0) - domain_lower(0) + (i(1) - domain_lower(1)) * num_cells(0) + depth * num_cells(0) * num_cells(1) +
            offset);
#elif(NDIM == 3)
    return (i(0) - domain_lower(0) + (i(1) - domain_lower(1)) * num_cells(0) +
            (i(2) - domain_lower(2)) * num_cells(0) * num_cells(1) +
            depth * num_cells(0) * num_cells(1) * num_cells(2) + offset);

#else
    return -1;
#endif

} // getIntegerMapping

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

    // Resize vectors
    nonoverlap_boxes.resize(n_subdomains);
    overlap_boxes.resize(n_subdomains);

    int counter = 0;
#if (NDIM == 2)

    const int Nx = subdomains(0);
    const int Ny = subdomains(1);

    int j_lower = patch_lower(1);
    for (int J = 0; J < Ny; ++J)
    {
        const int height = (cells(1) / box_size(1) ? box_size(1) : 0) + (J == (Ny - 1) ? cells(1) % box_size(1) : 0);
        int i_lower = patch_lower(0);
        for (int I = 0; I < Nx; ++I)
        {
            const int width = (cells(0) / box_size(0) ? box_size(0) : 0) + (I == (Nx - 1) ? cells(0) % box_size(0) : 0);

            SAMRAI::hier::IntVector<NDIM> box_lower(i_lower, j_lower);
            SAMRAI::hier::IntVector<NDIM> box_upper(i_lower + width - 1, j_lower + height - 1);
            nonoverlap_boxes[counter] = SAMRAI::hier::Box<NDIM>(box_lower, box_upper);
            overlap_boxes[counter] = SAMRAI::hier::Box<NDIM>::grow(nonoverlap_boxes[counter], overlap_size);

            ++counter;
            i_lower += width;
        }
        j_lower += height;
    }
#elif(NDIM == 3)

    const int Nx = subdomains(0);
    const int Ny = subdomains(1);
    const int Nz = subdomains(2);

    int k_lower = patch_lower(2);
    for (int K = 0; K < Nz; ++K)
    {
        const int depth = (cells(2) / box_size(2) ? box_size(2) : 0) + (K == (Nz - 1) ? cells(2) % box_size(2) : 0);

        int j_lower = patch_lower(1);
        for (int J = 0; J < Ny; ++J)
        {
            const int height =
                (cells(1) / box_size(1) ? box_size(1) : 0) + (J == (Ny - 1) ? cells(1) % box_size(1) : 0);

            int i_lower = patch_lower(0);
            for (int I = 0; I < Nx; ++I)
            {
                const int width =
                    (cells(0) / box_size(0) ? box_size(0) : 0) + (I == (Nx - 1) ? cells(0) % box_size(0) : 0);

                SAMRAI::hier::IntVector<NDIM> box_lower(i_lower, j_lower, k_lower);
                SAMRAI::hier::IntVector<NDIM> box_upper(i_lower + width - 1, j_lower + height - 1, k_lower + depth - 1);
                nonoverlap_boxes[counter] = SAMRAI::hier::Box<NDIM>(box_lower, box_upper);
                overlap_boxes[counter] = SAMRAI::hier::Box<NDIM>::grow(nonoverlap_boxes[counter], overlap_size);

                ++counter;
                i_lower += width;
            }
            j_lower += height;
        }
        k_lower += depth;
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

#endif //#ifndef included_IndexUtilities_inl_h
