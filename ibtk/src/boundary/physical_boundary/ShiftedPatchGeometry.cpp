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

#include <ibtk/ShiftedPatchGeometry.h>

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
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

ShiftedPatchGeometry::ShiftedPatchGeometry(const Patch<NDIM>& patch, const std::array<double, NDIM>& shift)
{
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
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
    d_geometry = new CartesianPatchGeometry<NDIM>(
        pgeom->getRatio(), touches_regular_bdry, touches_periodic_bdry, dx, x_lower.data(), x_upper.data());
} // ShiftedPatchGeometry

ShiftedPatchGeometry::Scope::Scope(const Patch<NDIM>& patch, const ShiftedPatchGeometry& shifted)
    : d_patch(const_cast<Patch<NDIM>&>(patch)), d_original(patch.getPatchGeometry())
{
    d_patch.setPatchGeometry(shifted.d_geometry);
} // Scope

ShiftedPatchGeometry::Scope::~Scope()
{
    d_patch.setPatchGeometry(d_original);
} // ~Scope

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
