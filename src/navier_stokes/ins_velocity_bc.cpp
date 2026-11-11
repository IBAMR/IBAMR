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

#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/ShiftedPatchGeometry.h>

#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <ArrayData.h>
#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <Patch.h>
#include <Variable.h>

#include <array>

#include "./ins_staggered_traction_stencil.h"
#include "./ins_velocity_bc.h"

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
namespace velocity_bc
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

Type
get_type(const double a, const double b, const char* const caller)
{
    const bool is_velocity = (a == 1.0 && b == 0.0);
    const bool is_traction = (a == 0.0 && b == 1.0);
    if (!is_velocity && !is_traction)
    {
        TBOX_ERROR(caller << ":\n"
                          << "  unsupported velocity boundary condition coefficients (a, b) = (" << a << ", " << b
                          << ")\n"
                          << "  use (a, b) = (1, 0) to prescribe the velocity or (0, 1) to prescribe the traction\n");
    }
    return is_velocity ? Type::VELOCITY : Type::TRACTION;
} // get_type

bool
traction_data_depend_on_velocity(RobinBcCoefStrategy<NDIM>* const comp_bc_coef,
                                 const unsigned int comp_idx,
                                 const Patch<NDIM>& patch,
                                 const BoundaryBox<NDIM>& bdry_box,
                                 const double fill_time)
{
    if (static_cast<int>(comp_idx) == bdry_box.getLocationIndex() / 2) return false;

    // Evaluate the coefficients at the faces where the fill sets the component: the cells of the boundary box and one
    // more face along the component axis, with the patch geometry shifted by half a cell along that axis.
    Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(bdry_box);
    bc_coef_box.upper()(comp_idx) += 1;
    Pointer<ArrayData<NDIM, double>> acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
    Pointer<ArrayData<NDIM, double>> bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
    Pointer<ArrayData<NDIM, double>> gcoef_data;
    std::array<double, NDIM> shift;
    shift.fill(0.0);
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    shift[comp_idx] = -0.5 * pgeom->getDx()[comp_idx];
    // Shifting the geometry of a copy of the patch leaves the patch geometry unchanged.
    Pointer<Patch<NDIM>> shifted_patch = traction_stencil::copy_patch(patch);
    const ShiftedPatchGeometry shifted_geometry(*shifted_patch, shift);
    const ShiftedPatchGeometry::ScopedGeometryReplacement geometry_replacement(*shifted_patch, shifted_geometry);
    comp_bc_coef->setBcCoefs(
        acoef_data, bcoef_data, gcoef_data, Pointer<Variable<NDIM>>(), *shifted_patch, bdry_box, fill_time);
    for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
    {
        if (get_type((*acoef_data)(it(), 0),
                     (*bcoef_data)(it(), 0),
                     "velocity_bc::traction_data_depend_on_velocity()") == Type::TRACTION)
        {
            return true;
        }
    }
    return false;
} // traction_data_depend_on_velocity

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace velocity_bc
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
