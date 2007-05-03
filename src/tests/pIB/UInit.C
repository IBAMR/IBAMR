// Filename: UInit.C
// Last modified: <03.May.2007 15:47:56 griffith@box221.cims.nyu.edu>
// Created on 03 May 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

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
#include <ArrayData.h>
#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <Index.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

UInit::UInit(
    const string& object_name,
    tbox::Pointer<hier::GridGeometry<NDIM> > grid_geom,
    tbox::Pointer<tbox::Database> input_db)
    : SetDataStrategy(object_name),
      d_object_name(object_name),
      d_grid_geom(grid_geom),
      d_V0(0.0)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
    assert(!grid_geom.isNull());
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
    tbox::Pointer< pdat::CellData<NDIM,double> > u_data = patch.getPatchData(data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!u_data.isNull());
#endif
    const hier::Box<NDIM>& patch_box = patch.getBox();
    const hier::Index<NDIM>& patch_lower = patch_box.lower();
    tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();

    const double* const grid_XLower = d_grid_geom->getXLower();
    const double* const grid_XUpper = d_grid_geom->getXUpper();

    const double L = grid_XUpper[0]-grid_XLower[0];

    const double* const XLower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();

    double X[NDIM];

    for (hier::Box<NDIM>::Iterator it(patch_box); it; it++)
    {
        const hier::Index<NDIM>& i = it();
        for (int d = 0; d < NDIM; ++d)
        {
            X[d] = XLower[d] + dx[d]*(static_cast<double>(i(d)-patch_lower(d))+0.5);
        }

        for (int d = 0; d < NDIM; ++d)
        {
            if (d != NDIM-1)
            {
                (*u_data)(i,d) = 0.0;
            }
            else
            {
                const double& x = X[0];
                (*u_data)(i,d) = (d_V0/(-0.25*L*L))*x*(x-L);
            }
        }
    }
    return;
}// setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
UInit::getFromInput(
    tbox::Pointer<tbox::Database> db)
{
    if (!db.isNull())
    {
        if (db->keyExists("V0"))
        {
            d_V0 = db->getDouble("V0");
        }
    }
    return;
}// getFromInput

//////////////////////////////////////////////////////////////////////////////
