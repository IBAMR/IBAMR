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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/private/ScopedShiftedPatchGeometry.h>

#include <tbox/Array.h>
#include <tbox/Pointer.h>

#include <CartesianPatchGeometry.h>
#include <Patch.h>
#include <PatchGeometry.h>

#include <array>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
namespace detail
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

ScopedShiftedPatchGeometry::ScopedShiftedPatchGeometry(Patch<NDIM>& patch, const std::array<double, NDIM>& shift)
    : d_patch(patch), d_original_geometry(patch.getPatchGeometry())
{
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = d_original_geometry;
    const double* const dx = pgeom->getDx();
    std::array<double, NDIM> x_lower, x_upper;
    Array<Array<bool>> touches_regular_bdry(NDIM), touches_periodic_bdry(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        x_lower[d] = pgeom->getXLower()[d] + shift[d];
        x_upper[d] = pgeom->getXUpper()[d] + shift[d];
        touches_regular_bdry[d].resizeArray(2);
        touches_periodic_bdry[d].resizeArray(2);
        for (int upperlower = 0; upperlower < 2; ++upperlower)
        {
            touches_regular_bdry[d][upperlower] = pgeom->getTouchesRegularBoundary(d, upperlower);
            touches_periodic_bdry[d][upperlower] = pgeom->getTouchesPeriodicBoundary(d, upperlower);
        }
    }
    d_patch.setPatchGeometry(new CartesianPatchGeometry<NDIM>(
        pgeom->getRatio(), touches_regular_bdry, touches_periodic_bdry, dx, x_lower.data(), x_upper.data()));
} // ScopedShiftedPatchGeometry

ScopedShiftedPatchGeometry::~ScopedShiftedPatchGeometry()
{
    d_patch.setPatchGeometry(d_original_geometry);
} // ~ScopedShiftedPatchGeometry

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
