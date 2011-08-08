// Filename: FeedbackForcer.C
// Created on 19 Oct 2007 by Boyce Griffith
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

#include "FeedbackForcer.h"

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
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <tbox/Utilities.h>

// BLITZ++ INCLUDES
#include <blitz/array.h>

// C++ STDLIB INCLUDES
#include <numeric>

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline double
delta_4(
    const double& r)
{
    const double r1 = abs(r);
    const double r2 = r*r;
    if (r1 < 1.0)
    {
        return 0.125*(3.0 - 2.0*r1 + sqrt( 1.0+ 4.0*r1-4.0*r2));
    }
    else if (r1 < 2.0)
    {
        return 0.125*(5.0 - 2.0*r1 - sqrt(-7.0+12.0*r1-4.0*r2));
    }
    else
    {
        return 0.0;
    }
}// delta_4
}

////////////////////////////// PUBLIC ///////////////////////////////////////

FeedbackForcer::FeedbackForcer(
    const string& object_name,
    const tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry)
    : d_U_data_idx(-1),
      d_kappa(0.0),
      d_object_name(object_name),
      d_grid_geometry(grid_geometry)
{
    return;
}// FeedbackForcer

FeedbackForcer::~FeedbackForcer()
{
    // intentionally blank
    return;
}// ~FeedbackForcer

void
FeedbackForcer::setDataOnPatch(
    const int data_idx,
    tbox::Pointer<hier::Variable<NDIM> > /*var*/,
    tbox::Pointer<hier::Patch<NDIM> > patch,
    const double data_time,
    const bool initial_time,
    tbox::Pointer<hier::PatchLevel<NDIM> > /*level*/)
{
    tbox::Pointer<pdat::CellData<NDIM,double> > U_data = patch->getPatchData(d_U_data_idx);
    tbox::Pointer<pdat::CellData<NDIM,double> > F_data = patch->getPatchData(data_idx);
    F_data->fillAll(0.0);

    if (initial_time) return;

    static const double T_ramp1 = 20.0;
    static const double T_ramp2 = 25.0;
    const double T_wgt1 =       (1.0/(1.0+tanh(2.0)))*(tanh(4.0*data_time/T_ramp1-2.0)+tanh(2.0));
    const double T_wgt2 = 1.0 - (1.0/(1.0+tanh(2.0)))*(tanh(4.0*data_time/T_ramp2-2.0)+tanh(2.0));

    const hier::Box<NDIM>& patch_box = patch->getBox();
    const hier::Index<NDIM>& patch_lower = patch_box.lower();
    tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const XLower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();

    const double* const grid_XUpper = d_grid_geometry->getXUpper();
    const double* const grid_XLower = d_grid_geometry->getXLower();
    blitz::TinyVector<double,NDIM> L;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        L[d] = grid_XUpper[d] - grid_XLower[d];
    }

    const double* const dx_coarsest = d_grid_geometry->getDx();
    const double H = dx_coarsest[0];
    const double R = 4.0*H;
    const int offset = static_cast<int>(R/dx[0])-1;

    const hier::IntVector<NDIM>& ratio = pgeom->getRatio();
    const hier::Box<NDIM> domain_box = hier::Box<NDIM>::refine(d_grid_geometry->getPhysicalDomain()[0],ratio);
    const hier::Index<NDIM>& domain_lower = domain_box.lower();

    hier::Box<NDIM> bdry_box = domain_box;
    bdry_box.upper()(0) = domain_lower(0)+offset;
    for (hier::Box<NDIM>::Iterator b(bdry_box*patch_box); b; b++)
    {
        const hier::Index<NDIM>& i = b();
        blitz::TinyVector<double,NDIM> X;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X[d] = XLower[d] + dx[d]*(static_cast<double>(i(d)-patch_lower(d))+0.5);
        }

        // Penalize deviations from the specified normal and
        // tangential velocity boundary conditions.
        const double edge_grader = delta_4((X[0]-grid_XLower[0])/(0.5*R))/delta_4(0.0);
        const double U_oo = T_wgt1*(1.0 + T_wgt2*0.5*sin(2.0*M_PI*X[1]/L[1]));
        (*F_data)(i,0) = d_kappa*edge_grader*(U_oo - (*U_data)(i,0));
    }
    return;
}// setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
