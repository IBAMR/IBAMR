// Filename: VelocityBcCoefs.C
// Created on 18 Dec 2007 by Boyce Griffith
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
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        L[d] = grid_XUpper[d] - grid_XLower[d];
    }

    const unsigned int location_index = bdry_box.getLocationIndex();
    const unsigned int bdry_normal_axis =  location_index / 2;
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
                for (unsigned int d = 0; d < NDIM; ++d)
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
