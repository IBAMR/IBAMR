// Filename: UInit.C
// Last modified: <01.Apr.2008 17:12:02 griffith@box221.cims.nyu.edu>
// Created on 19 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

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
    : SetDataStrategy(object_name),
      d_object_name(object_name),
      d_rho(80.0),
      d_delta(0.05)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
#endif

    if (NDIM != 2) TBOX_ERROR("only NDIM=2 is presently implemented!\n");

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
    TBOX_ASSERT(!U_data.isNull());
#endif
    const hier::Box<NDIM>& patch_box = patch.getBox();
    const hier::Index<NDIM>& patch_lower = patch_box.lower();
    tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();

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
