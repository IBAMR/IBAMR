// Filename: PressureBcCoefs.C
// Created on 04 May 2007 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

// BLITZ++ INCLUDES
#include <blitz/array.h>

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

    blitz::TinyVector<double,NDIM> X;

    const hier::BoundaryBox<NDIM> trimmed_bdry_box =
        PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, patch);
    const hier::Box<NDIM> bc_coef_box =
        PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

    const unsigned int location_index = bdry_box.getLocationIndex();
    const unsigned int bdry_normal_axis =  location_index / 2;

    for (hier::Box<NDIM>::Iterator b(bc_coef_box); b; b++)
    {
        const hier::Index<NDIM>& i = b();
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (d != bdry_normal_axis)
            {
                X[d] = XLower[d] + dx[d]*(static_cast<double>(i(d)-patch_lower(d))+0.5);
            }
            else
            {
                X[d] = XLower[d] + dx[d]*(static_cast<double>(i(d)-patch_lower(d)));
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
