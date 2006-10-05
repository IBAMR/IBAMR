// Filename: UInit.C
// Last modified: <03.Oct.2006 13:55:24 boyce@boyce-griffiths-powerbook-g4-15.local>
// Created on 19 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "UInit.h"

// IBAMR INCLUDES
#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#endif

// SAMRAI INCLUDES
#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#endif

#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <Index.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

UInit::UInit(
    const string& object_name,
    tbox::Pointer<hier::GridGeometry<NDIM> > grid_geom,
    tbox::Pointer<tbox::Database> input_db,
    const double nu)
    : SetDataStrategy(object_name),
      d_object_name(object_name),
      d_grid_geom(grid_geom),
      d_nu(nu)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
    assert(!grid_geom.isNull());
#endif
    d_object_name = object_name;
    d_grid_geom = grid_geom;
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!d_grid_geom.isNull());
#endif

    // Initialize object with data read from the input database.
    getFromInput(input_db);

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
    hier::Patch<NDIM>& patch,
    const double data_time,
    const bool initial_time)
{
    tbox::Pointer< pdat::CellData<NDIM,double> > U_data = patch.getPatchData(data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!U_data.isNull());
#endif
    const hier::Box<NDIM>& patch_box = patch.getBox();
    const hier::Index<NDIM>& patch_lower = patch_box.lower();
    tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();

    const double* const XLower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();

    double X[NDIM];
    const double t = data_time;

    U_data->fillAll(0.0);
    for (pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
    {
        const hier::Index<NDIM>& i = ic();
        for (int d = 0; d < NDIM; ++d)
        {
            X[d] = XLower[d] +
                dx[d]*(static_cast<double>(i(d)-patch_lower(d))+0.5);
        }

///////////////////////////////////////////////////////////////////////////
#if 0
        (void) t;

        const double rho = 80.0;
        const double delta = 0.05;

        if (X[1] < 0.5)
        {
            (*U_data)(i,0) = tanh(rho*(X[1]-0.25));
        }
        else
        {
            (*U_data)(i,0) = tanh(rho*(0.75-X[1]));
        }

        (*U_data)(i,1) = delta*sin(2.0*M_PI*(X[0]+0.25));
#endif
///////////////////////////////////////////////////////////////////////////
#if (NDIM == 2)
        (*U_data)(i,0) = 1.0 -
            2.0*(cos(2.0*M_PI*(X[0]-t))*sin(2.0*M_PI*(X[1]-t)))*
            exp(-8.0*M_PI*M_PI*d_nu*t);

        (*U_data)(i,1) = 1.0 +
            2.0*(sin(2.0*M_PI*(X[0]-t))*cos(2.0*M_PI*(X[1]-t)))*
            exp(-8.0*M_PI*M_PI*d_nu*t);
#endif
#if (NDIM == 3)
        const double A = 1.0;
        const double B = 1.0;
        const double C = 1.0;

        (*U_data)(i,0) = 1.0 + exp(-4.0*M_PI*M_PI*d_nu*t)*
            (C*cos(2*M_PI*(X[1]-t)) + A*sin(2*M_PI*(X[2]-t)));
        (*U_data)(i,1) = 1.0 + exp(-4.0*M_PI*M_PI*d_nu*t)*
            (B*sin(2*M_PI*(X[0]-t)) + A*cos(2*M_PI*(X[2]-t)));
        (*U_data)(i,2) = 1.0 + exp(-4.0*M_PI*M_PI*d_nu*t)*
            (B*cos(2*M_PI*(X[0]-t)) + C*sin(2*M_PI*(X[1]-t)));
#endif
    }
    return;
}// setDataOnPatch

/////////////////////////////// PRIVATE //////////////////////////////////////

void
UInit::getFromInput(
    tbox::Pointer<tbox::Database> db)
{
    if (!db.isNull())
    {
        // intentionally blank
    }
    return;
}// getFromInput

//////////////////////////////////////////////////////////////////////////////
