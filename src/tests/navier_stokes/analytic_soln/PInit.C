// Filename: PInit.C
// Last modified: <20.Feb.2007 02:35:52 boyce@bigboy.nyconnect.com>
// Created on 19 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

#include "PInit.h"

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

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

PInit::PInit(
    const string& object_name,
    const double nu)
    : SetDataStrategy(object_name),
      d_nu(nu)
{
    // intentionally blank
    return;
}// PInit

PInit::~PInit()
{
    // intentionally blank
    return;
}// ~PInit

void
PInit::setDataOnPatch(
    const int data_idx,
    tbox::Pointer<hier::Variable<NDIM> > var,
    hier::Patch<NDIM>& patch,
    const double data_time,
    const bool initial_time)
{
    tbox::Pointer< pdat::CellData<NDIM,double> > P_data = patch.getPatchData(data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!P_data.isNull());
#endif
    const hier::Box<NDIM>& patch_box = patch.getBox();
    const hier::Index<NDIM>& patch_lower = patch_box.lower();
    tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();

    const double* const XLower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();

    double X[NDIM];
    const double t = data_time;

    P_data->fillAll(0.0);
    for (pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
    {
        const hier::Index<NDIM>& i = ic();
        for (int d = 0; d < NDIM; ++d)
        {
            X[d] = XLower[d] +
                dx[d]*(static_cast<double>(i(d)-patch_lower(d))+0.5);
        }
#if (NDIM == 2)
        (*P_data)(i) = -(cos(4.0*M_PI*(X[0]-t)) + cos(4.0*M_PI*(X[1]-t)))*
            exp(-16.0*M_PI*M_PI*d_nu*t);
#endif
#if (NDIM == 3)
        const double A = 1.0;
        const double B = 1.0;
        const double C = 1.0;

        (*P_data)(i) = -exp(-8.0*M_PI*M_PI*d_nu*t)*
            (A*C*cos(2*M_PI*(X[1]-t))*sin(2*M_PI*(X[2]-t)) +
             A*B*sin(2*M_PI*(X[0]-t))*cos(2*M_PI*(X[2]-t)) +
             B*C*cos(2*M_PI*(X[0]-t))*sin(2*M_PI*(X[1]-t)));
#endif
    }
    return;
}// setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
