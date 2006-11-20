// Filename: FeedbackFSet.C
// Last modified: <20.Nov.2006 17:18:01 griffith@box221.cims.nyu.edu>
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
      d_kappa(100.0)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
    assert(!grid_geom.isNull());
#endif
    if (!input_db.isNull())
    {
        input_db->getDoubleWithDefault("kappa",d_kappa);
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
    //
    // NOTE: The following code DOES NOT assume the x-periodic
    // boundary is uniformly refined.
    const hier::Box<NDIM>& patch_box = patch.getBox();
    const hier::IntVector<NDIM>& ratio = patch.getPatchGeometry()->getRatio();
    const hier::Box<NDIM>& domain_box = d_grid_geom->getPhysicalDomain()(0);
    hier::Box<NDIM> upper_box = domain_box;
    upper_box.lower(0) = upper_box.upper(0)-1;
    upper_box = hier::Box<NDIM>::refine(upper_box,ratio);

    F_data->fillAll(0.0);
    for (pdat::CellIterator<NDIM> ic(patch_box*upper_box); ic; ic++)
    {
        const hier::Index<NDIM>& i = ic();
        (*F_data)(i,0) = d_kappa*(1.0-(*U_data)(i,0));
        (*F_data)(i,1) = d_kappa*(0.0-(*U_data)(i,1));
    }

    return;
}// setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
