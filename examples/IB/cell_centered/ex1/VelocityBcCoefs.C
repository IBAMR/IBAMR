// Filename: VelocityBcCoefs.C
// Created on 18 Dec 2007 by Boyce Griffith
//
// Copyright (c) 2002-2010 Boyce Griffith
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "VelocityBcCoefs.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/PhysicalBoundaryUtilities.h>

// SAMRAI INCLUDES
#include <CartesianPatchGeometry.h>
#include <tbox/Utilities.h>

// NAMESPACE
using namespace IBTK;

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

VelocityBcCoefs::VelocityBcCoefs(
    const string& object_name,
    const tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry)
    : d_object_name(object_name),
      d_grid_geometry(grid_geometry)
{
    // intentionally blank
    return;
}// VelocityBcCoefs

VelocityBcCoefs::~VelocityBcCoefs()
{
    // intentionally blank
    return;
}// ~VelocityBcCoefs

void
VelocityBcCoefs::setBcCoefs(
    tbox::Pointer<pdat::ArrayData<NDIM,double> >& acoef_data,
    tbox::Pointer<pdat::ArrayData<NDIM,double> >& bcoef_data,
    tbox::Pointer<pdat::ArrayData<NDIM,double> >& gcoef_data,
    const tbox::Pointer<hier::Variable<NDIM> >& variable,
    const hier::Patch<NDIM>& patch,
    const hier::BoundaryBox<NDIM>& bdry_box,
    double fill_time) const
{
    const hier::Box<NDIM>& patch_box = patch.getBox();
    const hier::Index<NDIM>& patch_lower = patch_box.lower();
    tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const XLower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();
    double X[NDIM];

    const double* const grid_XUpper = d_grid_geometry->getXUpper();
    const double* const grid_XLower = d_grid_geometry->getXLower();
    double L[NDIM];
    for (int d = 0; d < NDIM; ++d)
    {
        L[d] = grid_XUpper[d] - grid_XLower[d];
    }

    const int location_index = bdry_box.getLocationIndex();
    const int bdry_normal_axis =  location_index / 2;
    const hier::Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(bdry_box);

    if (location_index == 0) // x lower
    {
        if (!acoef_data.isNull()) acoef_data->fillAll(1.0);
        if (!bcoef_data.isNull()) bcoef_data->fillAll(0.0);

        if (!gcoef_data.isNull())
        {
            static const double T_ramp1 = 20.0;
            static const double T_ramp2 = 25.0;
            const double T_wgt1 =       (1.0/(1.0+tanh(2.0)))*(tanh(4.0*fill_time/T_ramp1-2.0)+tanh(2.0));
            const double T_wgt2 = 1.0 - (1.0/(1.0+tanh(2.0)))*(tanh(4.0*fill_time/T_ramp2-2.0)+tanh(2.0));

            for (hier::Box<NDIM>::Iterator b(bc_coef_box); b; b++)
            {
                const hier::Index<NDIM>& i = b();
                for (int d = 0; d < NDIM; ++d)
                {
                    if (d != bdry_normal_axis)
                    {
                        X[d] = XLower[d] + dx[d]*(double(i(d)-patch_lower(d))+0.5);
                    }
                    else
                    {
                        X[d] = XLower[d] + dx[d]*(double(i(d)-patch_lower(d)));
                    }
                }
                (*gcoef_data)(i,0) = T_wgt1*(1.0 + T_wgt2*0.5*sin(2.0*M_PI*X[1]/L[1]));
            }
        }
    }
    else // x upper, y lower, y upper, z lower, z upper
    {
        if (!acoef_data.isNull()) acoef_data->fillAll(0.0);
        if (!bcoef_data.isNull()) bcoef_data->fillAll(1.0);
        if (!gcoef_data.isNull()) gcoef_data->fillAll(0.0);
    }
    return;
}// setBcCoefs

hier::IntVector<NDIM>
VelocityBcCoefs::numberOfExtensionsFillable() const
{
    return 128;
}// numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
