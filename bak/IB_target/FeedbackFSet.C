// Filename: FeedbackFSet.C
// Last modified: <22.Jan.2007 23:24:19 boyce@bigboy.nyconnect.com>
// Created on 20 Nov 2006 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "FeedbackFSet.h"

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

FeedbackFSet::FeedbackFSet(
    const string& object_name,
    tbox::Pointer<hier::GridGeometry<NDIM> > grid_geom,
    tbox::Pointer<tbox::Database> input_db)
    : SetDataStrategy(object_name),
      d_U_var(NULL),
      d_U_context(NULL),
      d_object_name(object_name),
      d_grid_geom(grid_geom),
      d_kappa(1000.0),
      d_width0(1.0/32.0),
      d_width1(4.0/32.0)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
    assert(!grid_geom.isNull());
#endif
    if (!input_db.isNull())
    {
        d_kappa = input_db->getDoubleWithDefault("kappa",d_kappa);
        d_width0 = input_db->getDoubleWithDefault("width0",d_width0);
        d_width1 = input_db->getDoubleWithDefault("width1",d_width1);
    }
    return;
}// FeedbackFSet

FeedbackFSet::~FeedbackFSet()
{
    // intentionally blank
    return;
}// ~FeedbackFSet

void
FeedbackFSet::setDataOnPatch(
    const int data_idx,
    tbox::Pointer<hier::Variable<NDIM> > var,
    hier::Patch<NDIM>& patch,
    const double data_time,
    const bool initial_time)
{
    tbox::Pointer<pdat::CellData<NDIM,double> > F_data = patch.getPatchData(data_idx);
    tbox::Pointer<pdat::CellData<NDIM,double> > U_data = patch.getPatchData(d_U_var, d_U_context);
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!F_data.isNull());
    assert(!U_data.isNull());
#endif

    // Supply a feedback force so that the fluid velocity is U = (1,0)
    // at the x-periodic boundary.
    const hier::Box<NDIM>& patch_box = patch.getBox();
    const hier::Index<NDIM>& patch_lower = patch_box.lower();
    tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();

    const double* const XLower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();

    const double* const grid_XLower = d_grid_geom->getXLower();
    const double* const grid_XUpper = d_grid_geom->getXUpper();

    double X[NDIM];
    F_data->fillAll(0.0);
    for (pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
    {
        const hier::Index<NDIM>& i = ic();
        for (int d = 0; d < NDIM; ++d)
        {
            X[d] = XLower[d] +
                dx[d]*(static_cast<double>(i(d)-patch_lower(d))+0.5);
        }

        if ((tbox::Utilities::dabs(X[0]-grid_XLower[0]) <= d_width0) ||
            (tbox::Utilities::dabs(X[0]-grid_XUpper[0]) <= d_width0))
        {
            (*F_data)(i,0) = d_kappa*(1.0-(*U_data)(i,0));
            (*F_data)(i,1) = d_kappa*(0.0-(*U_data)(i,1));
        }
        else
        {
            if (X[0]-grid_XLower[0] <= d_width1)
            {
                const double fac = 0.5*(cos(M_PI*(X[0]-grid_XLower[0]-d_width0)/(d_width1-d_width0))+1.0);
                (*F_data)(i,0) = d_kappa*fac*(1.0-(*U_data)(i,0));
                for (int d = 1; d < NDIM; ++d)
                {
                    (*F_data)(i,d) = d_kappa*fac*(0.0-(*U_data)(i,d));
                }
            }
            else if (grid_XUpper[0]-X[0] <= d_width1)
            {
                const double fac = 0.5*(cos(M_PI*(grid_XUpper[0]-X[0]-d_width0)/(d_width1-d_width0))+1.0);
                (*F_data)(i,0) = d_kappa*fac*(1.0-(*U_data)(i,0));
                for (int d = 1; d < NDIM; ++d)
                {
                    (*F_data)(i,d) = d_kappa*fac*(0.0-(*U_data)(i,d));
                }
            }
        }
    }

    return;
}// setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
