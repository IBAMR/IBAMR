// Filename: PressureBcCoefs.C
// Created on 04 May 2007 by Boyce Griffith
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

#include "PressureBcCoefs.h"

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

PressureBcCoefs::PressureBcCoefs(
    const string& object_name)
    : d_object_name(object_name)
{
    // intentionally blank.
    return;
}// PressureBcCoefs

PressureBcCoefs::~PressureBcCoefs()
{
    return;
}// ~PressureBcCoefs

void
PressureBcCoefs::setBcCoefs(
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

    const hier::BoundaryBox<NDIM> trimmed_bdry_box =
        PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, patch);
    const hier::Box<NDIM> bc_coef_box =
        PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

    const int location_index = bdry_box.getLocationIndex();
    const int bdry_normal_axis =  location_index / 2;

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

        double dummy;
        double& a = (!acoef_data.isNull() ? (*acoef_data)(i,0) : dummy);
        double& b = (!bcoef_data.isNull() ? (*bcoef_data)(i,0) : dummy);
        double& g = (!gcoef_data.isNull() ? (*gcoef_data)(i,0) : dummy);

        switch (location_index)
        {
            case 0:   // lower x
                a = 0.0;
                b = 1.0;
                g = 0.0;
                break;
            case 1:   // upper x
                a = 0.0;
                b = 1.0;
                g = 0.0;
                break;
            case 2:   // lower y
                a = 1.0;
                b = 0.0;
                g = 0.0;
                break;
            case 3:   // upper y
                a = 1.0;
                b = 0.0;
                g = 50.0*0.25*(tanh(1000.0*fill_time-3.5)+1.0)*(tanh(1000.0*(0.03-fill_time)-3.5)+1.0);
                break;
            default:
                TBOX_ASSERT(false);
        }
    }
    return;
}// setBcCoefs

hier::IntVector<NDIM>
PressureBcCoefs::numberOfExtensionsFillable() const
{
    return 1;
}// numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
