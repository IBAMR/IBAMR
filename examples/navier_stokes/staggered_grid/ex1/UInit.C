// Filename: UInit.C
// Created on 19 Mar 2004 by Boyce Griffith
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

#include "UInit.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// SAMRAI INCLUDES
#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <SideData.h>
#include <SideIndex.h>
#include <SideIterator.h>
#include <Index.h>

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

UInit::UInit(
    const string& object_name,
    tbox::Pointer<tbox::Database> input_db)
    : CartGridFunction(object_name),
      d_object_name(object_name),
      d_rho(80.0),
      d_delta(0.05)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
#endif

    if (!input_db.isNull())
    {
        d_rho   = input_db->getDoubleWithDefault("rho"  , d_rho  );
        d_delta = input_db->getDoubleWithDefault("delta", d_delta);
    }

    return;
}// UInit

UInit::~UInit()
{
    // intentionally blank
    return;
}// ~UInit

void
UInit::setDataOnPatch(
    const int data_idx,
    tbox::Pointer<hier::Variable<NDIM> > var,
    tbox::Pointer<hier::Patch<NDIM> > patch,
    const double data_time,
    const bool initial_time,
    tbox::Pointer<hier::PatchLevel<NDIM> > level)
{
    tbox::Pointer< pdat::SideData<NDIM,double> > U_data = patch->getPatchData(data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!U_data.isNull());
#endif
    const hier::Box<NDIM>& patch_box = patch->getBox();
    const hier::Index<NDIM>& patch_lower = patch_box.lower();
    tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

    const double* const XLower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();

    double X[NDIM];
    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (pdat::SideIterator<NDIM> ic(patch_box,axis); ic; ic++)
        {
            const pdat::SideIndex<NDIM>& i = ic();

            for (int d = 0; d < NDIM; ++d)
            {
                if (d == axis)
                {
                    X[d] = XLower[d] + dx[d]*(double(i(d)-patch_lower(d))+0.0);
                }
                else
                {
                    X[d] = XLower[d] + dx[d]*(double(i(d)-patch_lower(d))+0.5);
                }
            }

            if (axis == 0)
            {
                if (X[1] < 0.5)
                {
                    (*U_data)(i) = tanh(d_rho*(X[1]-0.25));
                }
                else
                {
                    (*U_data)(i) = tanh(d_rho*(0.75-X[1]));
                }
            }
            else
            {
                (*U_data)(i) = d_delta*sin(2.0*M_PI*(X[0]+0.25));
            }
        }
    }
    return;
}// setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
