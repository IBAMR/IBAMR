// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
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

#include <ibamr/EnthalpyCorrectionBcCoefs.h>

#include <tbox/Utilities.h>

#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIndex.h>

#include <cmath>

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

EnthalpyCorrectionBcCoefs::EnthalpyCorrectionBcCoefs(const std::string& object_name) : d_object_name(object_name)
{
    return;
} // EnthalpyCorrectionBcCoefs

void
EnthalpyCorrectionBcCoefs::setTemperatureAndDerivativePatchDataIndices(const int T_idx, const int dT_dh_idx)
{
    d_T_idx = T_idx;
    d_dT_dh_idx = dT_dh_idx;

    return;
} // setTemperatureAndDerivativePatchDataIndices

void
EnthalpyCorrectionBcCoefs::setTemperatureBcCoef(RobinBcCoefStrategy<NDIM>* temperature_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(temperature_bc_coef);
#endif

    d_temperature_bc_coef = temperature_bc_coef;

    return;
} // setTemperatureBcCoef

void
EnthalpyCorrectionBcCoefs::setBcCoefs(ArrayPointer& acoef_data,
                                      ArrayPointer& bcoef_data,
                                      ArrayPointer& gcoef_data,
                                      const VariablePointer& variable,
                                      const Patch<NDIM>& patch,
                                      const BoundaryBox<NDIM>& bdry_box,
                                      const double fill_time) const
{
    if (!acoef_data && !bcoef_data && !gcoef_data) return;

#if !defined(NDEBUG)
    if (d_T_idx == IBTK::invalid_index || d_dT_dh_idx == IBTK::invalid_index)
    {
        TBOX_ERROR("EnthalpyCorrectionBcCoefs::setBcCoefs():\n"
                   << "Set the temperature and derivative patch-data indices before using the boundary condition.\n");
    }
#endif

    // Convert the temperature Robin coefficients to enthalpy correction coefficients.
    fillConvertedCoefficient(acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);

    return;
} // setBcCoefs

IntVector<NDIM>
EnthalpyCorrectionBcCoefs::numberOfExtensionsFillable() const
{
    // The enthalpy correction coefficients are defined on the boundary, so no extensions are fillable.
    return SAMRAI::hier::IntVector<NDIM>(0);
} // numberOfExtensionsFillable

/////////////////////////////// PRIVATE //////////////////////////////////////

void
EnthalpyCorrectionBcCoefs::fillConvertedCoefficient(ArrayPointer& acoef_data,
                                                    ArrayPointer& bcoef_data,
                                                    ArrayPointer& gcoef_data,
                                                    const VariablePointer& variable,
                                                    const Patch<NDIM>& patch,
                                                    const BoundaryBox<NDIM>& bdry_box,
                                                    const double fill_time) const
{
    if (!acoef_data && !bcoef_data && !gcoef_data) return;

    const Box<NDIM>& coef_box = acoef_data ? acoef_data->getBox() :
                                bcoef_data ? bcoef_data->getBox() :
                                             gcoef_data->getBox();

    // Initially fill the arrays with the temperature Robin coefficients.
    d_temperature_bc_coef->setBcCoefs(acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);

    Pointer<CellData<NDIM, double>> derivative_data = patch.getPatchData(d_dT_dh_idx);
    Pointer<CellData<NDIM, double>> T_data = patch.getPatchData(d_T_idx);

    if (!(derivative_data && T_data))
    {
        TBOX_ERROR("EnthalpyCorrectionBcCoefs::fillConvertedCoefficient():\n"
                   << "Required temperature and derivative patch data are not allocated.\n");
    }

    Pointer<CartesianPatchGeometry<NDIM>> patch_geometry = patch.getPatchGeometry();

    const int location_index = bdry_box.getLocationIndex();
    const int axis = location_index / 2;
    const bool upper_boundary = location_index % 2;
    const double dx = patch_geometry->getDx()[axis];

    const Box<NDIM>& patch_box = patch.getBox();

    for (Box<NDIM>::Iterator it(coef_box); it; it++)
    {
        const Index<NDIM>& boundary_index = it();

        CellIndex<NDIM> interior_index(boundary_index);

        // The normal index is set to the adjacent interior cell.
        // Clamp the tangential indices to not go outside the patch box.
        for (int d = 0; d < NDIM; ++d)
        {
            if (d == axis)
            {
                interior_index(d) = upper_boundary ? patch_box.upper()(d) : patch_box.lower()(d);
            }
            else
            {
                interior_index(d) = std::max(patch_box.lower()(d), std::min(interior_index(d), patch_box.upper()(d)));
            }
        }

        // Lag beta using only the adjacent interior-cell value.
        const double beta = (*derivative_data)(interior_index);

        // Compute g_h before overwriting the temperature coefficients a_T and b_T.
        if (gcoef_data)
        {
            CellIndex<NDIM> ghost_index(interior_index);
            ghost_index(axis) += upper_boundary ? 1 : -1;

#if !defined(NDEBUG)
            if (!T_data->getGhostBox().contains(ghost_index))
            {
                TBOX_ERROR("EnthalpyCorrectionBcCoefs::fillConvertedCoefficient():\n"
                           << "Temperature ghost data are required for the boundary residual.\n");
            }
#endif

            const double T_i = (*T_data)(interior_index);
            const double T_g = (*T_data)(ghost_index);
            const double T_boundary = 0.5 * (T_i + T_g);
            const double dT_dn = (T_g - T_i) / dx;

            (*gcoef_data)(boundary_index, 0) -=
                (*acoef_data)(boundary_index, 0) * T_boundary + (*bcoef_data)(boundary_index, 0) * dT_dn;
        }

        if (acoef_data) (*acoef_data)(boundary_index, 0) *= beta;
        if (bcoef_data) (*bcoef_data)(boundary_index, 0) *= beta;
    }

    return;
} // fillConvertedCoefficient

} // namespace IBAMR
