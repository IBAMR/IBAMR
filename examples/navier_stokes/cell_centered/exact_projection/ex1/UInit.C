// Filename: UInit.C
// Created on 19 Mar 2004 by Boyce Griffith
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
#include <CellData.h>
#include <CellIterator.h>
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

    if (NDIM != 2) TBOX_ERROR("only NDIM=2 is presently implemented!\n");

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
    tbox::Pointer< pdat::CellData<NDIM,double> > U_data = patch->getPatchData(data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!U_data.isNull());
#endif
    const hier::Box<NDIM>& patch_box = patch->getBox();
    const hier::Index<NDIM>& patch_lower = patch_box.lower();
    tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

    const double* const XLower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();

    double X[NDIM];
    for (pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
    {
        const hier::Index<NDIM>& i = ic();
        for (int d = 0; d < NDIM; ++d)
        {
            X[d] = XLower[d] + dx[d]*(double(i(d)-patch_lower(d))+0.5);
        }

        if (X[1] < 0.5)
        {
            (*U_data)(i,0) = tanh(d_rho*(X[1]-0.25));
        }
        else
        {
            (*U_data)(i,0) = tanh(d_rho*(0.75-X[1]));
        }

        (*U_data)(i,1) = d_delta*sin(2.0*M_PI*(X[0]+0.25));
    }
    return;
}// setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
