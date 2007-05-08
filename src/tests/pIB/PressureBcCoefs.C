// Filename: PressureBcCoefs.C
// Last modified: <05.May.2007 21:16:50 griffith@box221.cims.nyu.edu>
// Created on 04 May 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

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

// STOOLS INCLUDES
#include <stools/PhysicalBoundaryUtilities.h>

// SAMRAI INCLUDES
#include <CartesianPatchGeometry.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <cassert>

// NAMESPACE
using namespace STOOLS;

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
#if USING_OLD_ROBIN_BC_INTERFACE
    TBOX_ERROR("PressureBcCoefs::setBcCoefs():\n"
               << "  using incorrect solv::RobinBcCoefStrategy interface." << endl);
#endif

    const hier::Box<NDIM>& patch_box = patch.getBox();
    const hier::Index<NDIM>& patch_lower = patch_box.lower();
    tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();

    const double* const XLower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();

    double X[NDIM];

    const SAMRAI::hier::BoundaryBox<NDIM> trimmed_bdry_box =
        PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, patch);
    const SAMRAI::hier::Box<NDIM> bc_coef_box =
        PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

    const int location_index = bdry_box.getLocationIndex();
    const int bdry_normal_axis =  location_index / 2;

    for (SAMRAI::hier::Box<NDIM>::Iterator b(bc_coef_box); b; b++)
    {
        const SAMRAI::hier::Index<NDIM>& i = b();
        for (int d = 0; d < NDIM; ++d)
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
                assert(false);
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
